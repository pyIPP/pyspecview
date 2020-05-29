from .loader import * 
import os


def check(shot):
    #fastest check if the shotfile exist

    status = False
    path = shot_path+'/%d/L0/CON/%d'
    status |= os.path.isfile(path%(shot//10,shot))
    
    path = shot_path+'/%d/XX/FHC/%d'
    status |= os.path.isfile(path%(shot//10,shot))&(shot>30150)&(shot<31700)
    
    return status


class loader_DCN(loader):

    radial_profile=True

    def __init__(self,*args, **kargs):

        super(loader_DCN,self).__init__(*args, **kargs)
        
        self.groups = ['raw','unwrap']
        self.diag = 'FHC' if self.shot < 31700 else 'CON' 
        
        if not self.dd.Open(self.diag, self.shot):
            raise Exception('Fast density data are not avalible!')

        self.dd.Open('DCN', self.shot)

        self.signals = [s for s in self.dd.GetNames() if s[0] in ('H',)] 

        self.signals.sort()
 
        if self.dd.GetParameter('DCNgeo','H5LFSphi') is None:
            raise Exception('DCN geometry is not avalible')
        
        self.Phi     = {s:self.dd.GetParameter('DCNgeo','H%sLFSphi'%s[2])[0] for s in self.signals}
        self.R_start = {s:self.dd.GetParameter('DCNgeo','H%sLFS_R'%s[2])[0]/1e3 for s in self.signals}
        self.z_start = {s:self.dd.GetParameter('DCNgeo','H%sLFS_z'%s[2])[0]/1e3 for s in self.signals}
        self.R_end   = {s:self.dd.GetParameter('DCNgeo','H%sHFS_R'%s[2])[0]/1e3 for s in self.signals}
        self.z_end   = {s:self.dd.GetParameter('DCNgeo','H%sHFS_z'%s[2])[0]/1e3 for s in self.signals}
        self.dd.Close()
        
        #values from diaggeom
        self.R_start['V1'] = 1.785
        self.R_end['V1'] =  1.785
        self.z_start['V1'] =  1.200
        self.z_end['V1'] = -1.200
        self.Phi['V1'] = np.nan
        
        self.R_start['V2'] = 1.200
        self.R_end['V2'] =  1.200
        self.z_start['V2'] =  1.200
        self.z_end['V2'] = -1.200
        self.Phi['V2'] = np.nan
        
        if self.diag == 'CON': self.signals+= ['V1','V2'] 


    def unwrap_ne(self,tvec, sig,ref):
        #remove referece signal and get unwraped phase
        #BUG slow!!
        
        #not working for 32427, H1, why? 
        
        from scipy.signal import argrelmax,argrelmin,argrelextrema
        from numpy import random

        print( 'unwraping DCN....')
        #preprocess signals
        sig = np.single(sig)
        sig -= np.mean(sig)
        sig/= np.std(sig)*np.sqrt(2)*1.05
        ref = np.single(ref)
        ref -= np.mean(ref)
        ref/= np.std(ref)*np.sqrt(2)*1.05    

        density_constant = 0.572e19/(2*np.pi)#m^-2

        #troubles with dicretizations steps
        ref+= np.random.rand(len(tvec))*1e-5
        sig+= np.random.rand(len(tvec))*1e-5
        dt = (tvec[-1]-tvec[0])/(len(tvec)-1)
        
        #idenitify an envelope
        N = int(round(1/dt/5e4))

        def loc_max(sig,N):
            #find local maxima (much fastter than argrelmin)
            imax = np.argmax(sig[:len(sig)//(2*N)*(2*N)].reshape(-1,2*N),axis=1)
            ind = ~((imax==0)|(imax==2*N-1))
            ind |= (imax==0) & (np.r_[0,imax[:-1]]==2*N-1)
            imax = imax[ind]+np.arange(len(ind))[ind]*2*N
    
            iimaxL =  ( np.r_[np.diff(imax)] < N )  
            iimaxR =  ( np.r_[N,np.diff(imax)[:-1]] < N )  
            comp_ind = sig[imax][:-1][iimaxL]<sig[imax][:-1][iimaxR]
            iimaxL[iimaxL] &= comp_ind
            iimaxR[iimaxR] &= ~comp_ind
            iimax = iimaxL|iimaxR
            imax = imax[:-1][~iimax]
            return imax
        
        pos = loc_max(sig,N)
        max_sig = np.interp(np.arange(len(tvec)),pos, sig[pos])

        pos = loc_max(-sig,N)
        min_sig = np.interp(np.arange(len(tvec)),pos, sig[pos])
        

        pos = loc_max(ref,N)

        max_ref = np.interp(np.arange(len(tvec)),pos, ref[pos])

        pos = loc_max(-ref,N)

        min_ref = np.interp(np.arange(len(tvec)),pos, ref[pos])

        #normalize for unitary envelope (BUG can introduce fake higher harminics!!!!)
        sig = (sig-(max_sig+min_sig)/2)/(max_sig-min_sig)*2
        ref = (ref-(max_ref+min_ref)/2)/(max_ref-min_ref)*2

        ref[ref>1] = 1-1e-6
        ref[ref<-1] = -1+1e-6
        sig[sig>1] = 1-1e-6
        sig[sig<-1] = -1+1e-6


        #unwrap reference phase
        ref = np.arcsin(ref,ref)

        ref_sign = np.ones_like(ref)
        ref_sign[np.abs(ref) == np.pi/2] = -1

        ref = np.cumprod(ref_sign)*ref

        ref = np.unwrap(ref*2)/2
        ref*= np.sign(np.mean(ref))

        #unwrap signal phase
        sig = np.arcsin(sig,sig)

        sig_sign = np.ones_like(sig)
        sig_sign[np.abs(sig) == np.pi/2] = -1

        sig = np.cumprod(sig_sign)*sig

        sig = np.unwrap(sig*2, np.pi/2)/2
        sign_sig = np.sign(np.mean(sig))
        sig*= sign_sig

        #remove phase jumps and some low frequency information
        knots,vals = [],[]
        ind_knots = np.r_[0:len(tvec):N*1000,-1]
        peaks = np.infty
        ind = slice(None,None)
        dsig,tvec_ = ref-sig,tvec
        for i in range(7):
            if not np.any(peaks>11):
                break
            knots.extend(tvec_[ind_knots])
            vals.extend(dsig[ind_knots])
            sort_ind = np.argsort(knots)
            retrofit = np.interp(tvec_, np.array(knots)[sort_ind], np.array(vals)[sort_ind])
            peaks = np.abs(dsig-retrofit)
            ind = peaks>10
            ind_knots = argrelmax(peaks[ind],0,100)[0]
            dsig = dsig[ind]
            tvec_ = tvec_[ind]

        sort_ind = np.argsort(knots)
        retrofit = np.interp(tvec, np.array(knots)[sort_ind], np.array(vals)[sort_ind])

        dn = (ref-sig-retrofit)*density_constant

        return np.single(dn)

    def get_signal(self,group, name,calib=False,tmin=None,tmax=None):
        
        if tmin is None:    tmin = self.tmin
        if tmax is None:    tmax = self.tmax
            
        if np.size(name) > 1:
            return [self.get_signal(group, n, calib=calib, tmin=tmin,tmax=tmax) for n in name]
        
        self.dd.Open(self.diag , self.shot, experiment=self.exp, edition=self.ed)

        if name[0] == 'H':  #DCN lasers
            tvec = self.dd.GetTimebase('DCN_'+name,cal=True)
            nbeg, nend = tvec.searchsorted((tmin,tmax))

            sig = self.dd.GetSignal('DCN_'+name, nbeg= nbeg, nend = nend-1 )

            tvec = tvec[nbeg:nend]
            
            if group == 'unwrap':  #it will also increase noise and small modes can be lost
                ref = self.dd.GetSignal('DCN_ref', nbeg= nbeg, nend = nend-1)
                sig = self.unwrap_ne(tvec, sig,ref)

        elif name[0] == 'V':  #CO2 lasers
            n = int(name[-1])
            tvec  = self.dd.GetTimebase('V%dHePh2'%n)
            
            #the He-Ne signal is not used because it was reducing the SNR at high frequencies!
            #cos_1 = self.dd.GetSignalCalibrated('V%dHePh2'%n)-1
            cos_2 = self.dd.GetSignalCalibrated('V%dCOPh2'%n)-1
            #sin_1 = self.dd.GetSignalCalibrated('V%dHePh1'%n)-1
            sin_2 = self.dd.GetSignalCalibrated('V%dCOPh1'%n)-1
            #phi1  = unwrap(np.arctan2(cos_1,sin_1))
            phi2  = np.unwrap(np.arctan2(cos_2,sin_2)) #this signal has stronger mode signal
            sig = phi2[:min(len(phi2), len(tvec))]
            #sig2 = phi1[:min(len(phi1), len(tvec))]

            tvec = tvec[:len(sig)]

        return tvec,sig

    
    def get_names(self,group):
      
        return self.signals


    def get_rho(self,group,names,time,dR=0,dZ=0):

        R_start = np.array([self.R_start[name] for name in names])
        z_start = np.array([self.z_start[name] for name in names])
        R_end = np.array([self.R_end[name] for name in names])
        z_end = np.array([self.z_end[name] for name in names])
        Phi = np.array([self.Phi[name] for name in names])
        

        rho_tg,theta_tg, R,Z = super(loader_DCN,self).get_rho(time,R_start,z_start, \
                                        Phi,R_end,z_end,Phi,dR=dR,dZ=dZ)
        
        return rho_tg,theta_tg,R,Z

        
    def signal_info(self,group,name,time):
        
        rho_tg = self.get_rho(group,[name,],time)[0]
        
        phi = self.Phi[name]
        
        info = str(name)+' Phi: %f deg, rho_tg: %.2f'%(phi,rho_tg)
        return info
