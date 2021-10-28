from .loader import * 
from scipy.signal import argrelmin

logger = logging.getLogger('pyspecview.density')
logger.setLevel(logging.INFO)


def check(shot):
    """
    Check if any shotfile exists
    """

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
        
        self.groups = ['raw', 'unwrap']
        self.diag = 'FHC' if self.shot < 31700 else 'CON' 
        
        sfo = sf.SFREAD(self.diag, self.shot)
        if not sfo.status:
            raise Exception('Fast density data are not avalible!')

        dcn = sf.SFREAD('DCN', self.shot)

        self.signals = [s for s in dcn.objects if s[0] in ('H',)] 

        self.signals.sort()
 
        dcngeo = dcn('DCNgeo')
        if 'H5LFSphi' not in dcngeo.keys():
            raise Exception('DCN geometry is not avalible')

        self.Phi     = {s: dcngeo['H%sLFSphi' %s[2]]     for s in self.signals}
        self.R_start = {s: dcngeo['H%sLFS_R'  %s[2]]/1e3 for s in self.signals}
        self.z_start = {s: dcngeo['H%sLFS_z'  %s[2]]/1e3 for s in self.signals}
        self.R_end   = {s: dcngeo['H%sHFS_R'  %s[2]]/1e3 for s in self.signals}
        self.z_end   = {s: dcngeo['H%sHFS_z'  %s[2]]/1e3 for s in self.signals}

        #values from diaggeom
        self.R_start['V1'] =  1.785
        self.R_end  ['V1'] =  1.785
        self.z_start['V1'] =  1.200
        self.z_end  ['V1'] = -1.200
        self.Phi    ['V1'] = np.nan
        
        self.R_start['V2'] =  1.200
        self.R_end  ['V2'] =  1.200
        self.z_start['V2'] =  1.200
        self.z_end  ['V2'] = -1.200
        self.Phi    ['V2'] = np.nan
        
        if self.diag == 'CON': self.signals+= ['V1','V2'] 


    def unwrap_ne(self,tvec, sig, ref):
        #remove referece signal and get unwraped phase
        #BUG slow!!
        
        #not working for 32427, H1, why? 
        
        from scipy.signal import argrelmax,argrelmin,argrelextrema
        from numpy import random

        logger.info('Unwrapping DCN...')
        logger.debug('Unwrap: %d %d %d', len(ref), len(sig), len(tvec))
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

        pos = argrelmax(sig, order=N)[0]
        max_sig = np.interp(np.arange(len(tvec)), pos, sig[pos])

        pos = argrelmax(-sig, order=N)[0]
        min_sig = np.interp(np.arange(len(tvec)), pos, sig[pos])

        pos = argrelmax(ref, order=N)[0]
        max_ref = np.interp(np.arange(len(tvec)), pos, ref[pos])

        pos = argrelmax(-ref, order=N)[0]
        min_ref = np.interp(np.arange(len(tvec)), pos, ref[pos])

        #normalize for unitary envelope (BUG can introduce fake higher harminics!!!!)
        sig = (sig-(max_sig+min_sig)/2)/(max_sig-min_sig)*2
        ref = (ref-(max_ref+min_ref)/2)/(max_ref-min_ref)*2

        ref[ref> 1] =  1-1e-6
        ref[ref<-1] = -1+1e-6
        sig[sig> 1] =  1-1e-6
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


    def get_signal(self, group, name,calib=False,tmin=None,tmax=None):

        if tmin is None:    tmin = self.tmin
        if tmax is None:    tmax = self.tmax
            
        if np.size(name) > 1:
            return [self.get_signal(group, n, calib=calib, tmin=tmin,tmax=tmax) for n in name]
        
        sfo = sf.SFREAD(self.diag , self.shot, experiment=self.exp, edition=self.ed)

        if name[0] == 'H':  #DCN lasers
            tvec = sfo.gettimebase('DCN_'+name)
            nbeg, nend = tvec.searchsorted((tmin, tmax))
            logger.debug('%s %.2f %.2f %.2f %.2f %d %d', self.diag, tvec[0], tvec[-1], tmin, tmax, nbeg, nend)
            sig = sfo.getobject('DCN_'+name, nbeg=nbeg, nend=nend)

            tvec = tvec[nbeg:nend]
            
            if group == 'unwrap':  #it will also increase noise and small modes can be lost
                ref = sfo.getobject('DCN_ref', nbeg=nbeg, nend=nend)
                sig = self.unwrap_ne(tvec, sig, ref)

        elif name[0] == 'V':  #CO2 lasers
            n = int(name[-1])
            tvec  = sfo.gettimebase('V%dHePh2'%n)
            
            #the He-Ne signal is not used because it was reducing the SNR at high frequencies!
            #cos_1 = self.dd.GetSignalCalibrated('V%dHePh2'%n)-1
            cos_2 = sfo.getobject('V%dCOPh2' %n, cal=True) - 1
            #sin_1 = self.dd.GetSignalCalibrated('V%dHePh1'%n)-1
            sin_2 = sfo.getobject('V%dCOPh1' %n, cal=True) - 1
            #phi1  = unwrap(np.arctan2(cos_1,sin_1))
            phi2  = np.unwrap(np.arctan2(cos_2,sin_2)) #this signal has stronger mode signal
            sig = phi2[:min(len(phi2), len(tvec))]
            #sig2 = phi1[:min(len(phi1), len(tvec))]

            tvec = tvec[:len(sig)]

        return tvec,sig

    
    def get_names(self,group):
      
        return self.signals


    def get_rho(self, group, names, time, dR=0, dZ=0):

        R_start = np.array([self.R_start[name] for name in names])
        z_start = np.array([self.z_start[name] for name in names])
        R_end   = np.array([self.R_end  [name] for name in names])
        z_end   = np.array([self.z_end  [name] for name in names])
        Phi     = np.array([self.Phi    [name] for name in names])

        rho_tg,theta_tg, R,Z = super(loader_DCN,self).get_rho(time,R_start,z_start, \
                                        Phi,R_end,z_end,Phi,dR=dR,dZ=dZ)

        return rho_tg,theta_tg,R,Z

        
    def signal_info(self,group,name,time):
        
        rho_tg = self.get_rho(group,[name,],time)[0]
        
        phi = self.Phi[name]
        
        info = str(name)+' Phi: %f deg, rho_tg: %.2f'%(phi,rho_tg)
        return info
