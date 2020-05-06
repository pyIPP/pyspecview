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
        
        
        self.Phi     = {s:self.dd.GetParameter('DCNgeo','H%sLFSphi'%s[2]) for s in self.signals}
        self.R_start = {s:self.dd.GetParameter('DCNgeo','H%sLFS_R'%s[2])/1e3 for s in self.signals}
        self.z_start = {s:self.dd.GetParameter('DCNgeo','H%sLFS_z'%s[2])/1e3 for s in self.signals}
        self.R_end   = {s:self.dd.GetParameter('DCNgeo','H%sHFS_R'%s[2])/1e3 for s in self.signals}
        self.z_end   = {s:self.dd.GetParameter('DCNgeo','H%sHFS_z'%s[2])/1e3 for s in self.signals}
        self.dd.Close()

        
        #values from diaggeom
        self.R_start['V1'] = 1.785
        self.R_end['V1'] =  1.785
        self.z_start['V1'] =  1.200
        self.z_end['V1'] = -1.200
        self.Phi['V1'] = nan
        
        self.R_start['V2'] = 1.200
        self.R_end['V2'] =  1.200
        self.z_start['V2'] =  1.200
        self.z_end['V2'] = -1.200
        self.Phi['V2'] = nan
        
        if self.diag == 'CON': self.signals+= ['V1','V2'] 


    def unwrap_ne(self,tvec, sig,ref):
        #remove referece signal and get unwraped phase
        #BUG slow!!
        
        #not working for 32427, H1, why? 
        
        from scipy.signal import argrelmax,argrelmin,argrelextrema
        from numpy import random

        print( 'unwraping DCN....')
        #preprocess signals
        sig = single(sig)
        sig -= mean(sig)
        sig/= std(sig)*sqrt(2)*1.05
        ref = single(ref)
        ref -= mean(ref)
        ref/= std(ref)*sqrt(2)*1.05    

        density_constant = 0.572e19/(2*pi)#m^-2

        #troubles with dicretizations steps
        ref+= random.rand(len(tvec))*1e-5
        sig+= random.rand(len(tvec))*1e-5
        dt = (tvec[-1]-tvec[0])/len(tvec)
        
        #idenitify an envelope
        N = int(round(1/dt/5e4))

        def loc_max(sig,N):
            #find local maxima (much fastter than argrelmin)
            imax = argmax(sig[:len(sig)//(2*N)*(2*N)].reshape(-1,2*N),axis=1)
            ind = ~((imax==0)|(imax==2*N-1))
            ind |= (imax==0) & (r_[0,imax[:-1]]==2*N-1)
            imax = imax[ind]+arange(len(ind))[ind]*2*N
    
            iimaxL =  ( r_[diff(imax)] < N )  
            iimaxR =  ( r_[N,diff(imax)[:-1]] < N )  
            comp_ind = sig[imax][iimaxL]<sig[imax][iimaxR]
            iimaxL[iimaxL] &= comp_ind
            iimaxR[iimaxR] &= ~comp_ind
            iimax = iimaxL|iimaxR
            imax = imax[~iimax]
            return imax
        
        pos = loc_max(sig,N)
        max_sig = interp(arange(len(tvec)),pos, sig[pos])

        pos = loc_max(-sig,N)
        min_sig = interp(arange(len(tvec)),pos, sig[pos])
        

        pos = loc_max(ref,N)

        max_ref = interp(arange(len(tvec)),pos, ref[pos])

        pos = loc_max(-ref,N)

        min_ref = interp(arange(len(tvec)),pos, ref[pos])

        #normalize for unitary envelope (BUG can introduce fake higher harminics!!!!)
        sig = (sig-(max_sig+min_sig)/2)/(max_sig-min_sig)*2
        ref = (ref-(max_ref+min_ref)/2)/(max_ref-min_ref)*2

        ref[ref>1] = 1-1e-6
        ref[ref<-1] = -1+1e-6
        sig[sig>1] = 1-1e-6
        sig[sig<-1] = -1+1e-6


        #unwrap reference phase
        ref = arcsin(ref,ref)

        ref_sign = ones_like(ref)
        ref_sign[abs(ref)==pi/2] = -1

        ref = cumprod(ref_sign)*ref

        ref = unwrap(ref*2)/2
        ref*= sign(mean(ref))

        
        #unwrap signal phase
        sig = arcsin(sig,sig)

        sig_sign = ones_like(sig)
        sig_sign[abs(sig)==pi/2] = -1

        sig = cumprod(sig_sign)*sig

        sig = unwrap(sig*2,pi/2)/2
        sign_sig = sign(mean(sig))
        sig*= sign_sig


        #remove phase jumps and some low frequency information
        knots,vals = [],[]
        ind_knots = r_[list(range(0,len(tvec),N*1000)),-1]
        peaks = infty
        ind = slice(None,None)
        dsig,tvec_ = ref-sig,tvec
        for i in range(7):
            if not any(peaks>11):
                break
            knots.extend(tvec_[ind_knots])
            vals.extend(dsig[ind_knots])
            sort_ind = argsort(knots)
            retrofit = interp(tvec_, array(knots)[sort_ind], array(vals)[sort_ind])
            peaks = abs(dsig-retrofit)
            ind = peaks>10
            ind_knots = argrelmax(peaks[ind],0,100)[0]
            dsig = dsig[ind]
            tvec_ = tvec_[ind]

        sort_ind = argsort(knots)
        retrofit = interp(tvec, array(knots)[sort_ind], array(vals)[sort_ind])

        dn = (ref-sig-retrofit)*density_constant

        return dn



    def get_signal(self,group, name,calib=False,tmin=None,tmax=None):
        
        if tmin is None:    tmin = self.tmin
        if tmax is None:    tmax = self.tmax
            
        if size(name) > 1:
            return [self.get_signal(group, n, calib=calib, tmin=tmin,tmax=tmax) for n in name]
        
        
        self.dd.Open(self.diag , self.shot, experiment=self.exp, edition=self.ed)
      
        if name[0] == 'H':  #DCN lasers
            tvec = self.dd.GetTimebase('DCN_'+name)
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
            #phi1  = unwrap(arctan2(cos_1,sin_1))
            phi2  = unwrap(arctan2(cos_2,sin_2)) #this signal has stronger mode signal
            sig = phi2[:min(len(phi2), len(tvec))]
            #sig2 = phi1[:min(len(phi1), len(tvec))]

            tvec = tvec[:len(sig)]






            #plot(sig[::1000]*10.6e3)
            #plot(sig2[::1000]*633)
            
            
            #A1 = self.dd.GetSignalCalibrated('V1HeAmp')
            #A2 = self.dd.GetSignalCalibrated('V1COAmp')

            #ne_co2 = -(sig2*633-sig*10.6e3)/(2.82e-15*(633**2-10.6e3**2)*1e-9)
            #ne_co2-= ne_co2[tvec<0].mean()
            ##plot(tvec[::10],  ne_co2[::10]-ne_co2[tvec<0].mean())
            ##plot(self.dd.GetTimebase('V-1'),self.dd.GetSignalCalibrated('V-1'))
            #co2 = fft.irfft(fft.rfft(ne_co2)[:1601])/len(ne_co2)*1601*2
            #tvecco2 = linspace(tvec[0],tvec[-1], len(co2))
            #plot(abs(fft.rfft(ne_co2)));show()
          
            #plot(tvec,ne_co2)
            #corr = 
   
            #import dd   
            #dd = dd.shotfile() 
            #dd.Open('DCN',self.shot)
            #H1 = dd.GetSignalCalibrated('H-1')
            #tvecH1 = dd.GetTimebase('H-1')
            #H1-= double(H1[tvecH1<0.05]).mean()
            #dd.Open('GQH',self.shot )
            #GQH_tvec = dd.GetTimebase('TIMEF')
            #lenH1 = interp(tvecH1,GQH_tvec,dd.GetSignal('lenH-1'))
            #lenV1 = interp(tvec,GQH_tvec,dd.GetSignal('lenV-1'))
            #co2-= mean(co2[(tvecco2 > 7)&(tvecco2<8)])
            #from scipy.signal import savgol_filter,decimate,medfilt
            #import IPython
            #IPython.embed()
            
            ##ne_co2_smooth = savgol_filter(ne_co2,1001,5)
            ##c = mean(ne_co2)
            ##ne_co2_smooth = medfilt(ne_co2,1001)
            #n = 2000
            #ne_co2_smooth = reshape(ne_co2[:(ne_co2.size/n)*n],(-1,n)).mean(1)
            #tvec_smooth = reshape(tvec[:(tvec.size/n)*n],(-1,n)).mean(1)
            #lenV1_smooth = reshape(lenV1[:(lenV1.size/n)*n],(-1,n)).mean(1)
            #ne_co2_smooth-= mean( ne_co2_smooth[tvec_smooth>7] )

            #plot(tvecH1,H1/lenH1)
            #plot(tvec_smooth, ne_co2_smooth/lenV1_smooth*1.08)
            #xlabel('time [s]')
            #ylabel('ne [m$^{-3}$]')
            #show()

            
            #figure()
            #plot(tvec[::10], A1[::10]);plot(tvec[::10], A2[::10])
            
            
        #import IPython
        #IPython.embed()
   
        #self.dd.Close()

        return tvec,sig

    
    def get_names(self,group):
  
    
        return self.signals
    

    
    def get_rho(self,group,names,time,dR=0,dZ=0):

    
        
        R_start = array([self.R_start[name] for name in names])
        z_start = array([self.z_start[name] for name in names])
        R_end = array([self.R_end[name] for name in names])
        z_end = array([self.z_end[name] for name in names])
        Phi = array([self.Phi[name] for name in names])

        rho_tg,theta_tg, R,Z = super(loader_DCN,self).get_rho(time,R_start,z_start,
                                                Phi,R_end,z_end,Phi,dR=dR,dZ=dZ)
        
        return rho_tg,theta_tg,R,Z

    
        
    def signal_info(self,group,name,time):
        
        rho_tg = self.get_rho(group,[name,],time)[0]
        
        phi = self.Phi[name]
        
        info = str(name)+' Phi: %f deg, rho_tg: %.2f'%(phi,rho_tg)
        return info
    
  
 

def main():
    import os,sys
    sys.path.append('/afs/ipp/home/t/todstrci/TRANSP/')
    from . import dd   



    loader = loader_DCN(32961)
    print((  loader.get_names('DCN')))

    
    t,x = loader.get_signal('DCN','H-0')
    #plot(t,x)       ;show()


if __name__ == "__main__":
    main()


