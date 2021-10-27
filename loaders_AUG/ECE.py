from .loader import * 
from scipy.interpolate import interp1d

logger = logging.getLogger('pyspecview.ece')
logger.setLevel(logging.INFO)


def check(shot):
    """
    Check if any shotfile exists
    """

    path = shot_path+'/%d/XX/%s/%d'

    for diag in ('RMC', 'RMA','RMB'):
        if os.path.isfile(path%(shot//10,diag, shot)):
            return True
        
    path = shot_path+'/%d/L0/%s/%d'

    if os.path.isfile(path%(shot//10,'RAD', shot)):
        return True

    return False

 
class loader_ECE(loader):

    radial_profile = True
    units = 'eV'

    def __init__(self,*args, **kargs):

        super(loader_ECE,self).__init__(*args, **kargs)

        self.mixers = None
        self.corrupted_calibration=False
        self.names = []
        new_rmd_shotfile = False

        rmc = sf.SFREAD('RMC', self.shot, experiment=self.exp, edition=self.ed)
        rma = sf.SFREAD('RMA', self.shot, experiment=self.exp, edition=self.ed)
        rad = sf.SFREAD('RAD', self.shot, experiment=self.exp, edition=self.ed)
        if rmc.status:
            self.shotfile = 'RMC' # new diag
            self.mixers = rmc('parms-A')['IFGROUP']
        elif rma.status:
            self.shotfile = 'RMA' # old diag
        elif rad.status:
            self.shotfile = 'RAD' # old+slow diagnostic
        else:
            raise Exception('Te shotfile is missing')

        if len(self.names) == 0: 
            self.names = np.arange(60)+1

        if self.mixers is None: # it was missing in old shotfiles 
            self.mixers = np.int_(np.r_[(1,)*36,(2,)*12, (3,)*12])
            
        self.groups = ['mixer %d'%d for d in np.unique(self.mixers)]

        #get the final clalibration factors
        if not new_rmd_shotfile:
            sfo = sf.SFREAD('CEC', self.shot)
        else:
            sfo = sf.SFREAD('RMD', self.shot)

        self.RZtime = sfo('rztime')
        self.R = sfo('R-A')
        self.z = sfo('z-A')
        if 'Trad-B' in sfo.objects:
            self.R = np.hstack((self.R, sfo('R-A')))
            self.z = np.hstack((self.z, sfo('z-A')))

        self.calfact = sfo('parms-A')['calfact']
        if self.calfact is None or all(self.calfact == 0):
            #BUG extract calibration directly from CEC!
    
            tvec = sfo.gettimebase('Trad-A')
            i0, imin, imax = tvec.searchsorted((0, self.tmin, self.tmax))
            TeCEC = sfo.getobject('Trad-A',nbeg=imin, nend=imax-1).mean(0)

            rad = sf.SFREAD('RAD',self.shot)
            tvec = rad('TIME-AD')
            i0 = tvec.searchsorted(0)
            Te1 = rad('SI-MI-A', cal=True, nend=imax+1)
            Te2 = rad('SI-MI-B', cal=True, nend=imax+1)
            Te3 = rad('SI-MI-C', cal=True, nend=imax+1)
            Te4 = rad('SI-MI-D', cal=True, nend=imax+1)
            
            TeRAD0 = np.hstack((Te1,Te2,Te3,Te4))[:i0].mean(0)
            TeRAD  = np.hstack((Te1,Te2,Te3,Te4))[i0:].mean(0)
            TeRAD -= TeRAD0
            self.calfact = TeCEC/abs(TeRAD)

        #sometimes missing calibration
        if hasattr(self, 'CALA2_M1') and all(self.CALA2_M2 == 1):
            self.CALA1_M2 = -self.calfact[:30]
            self.CALA2_M2 = -self.calfact[30:]


    def get_signal_groups(self):
        return  self.groups

    def get_names(self,group):
        mixer = int(group[-1])

        return self.names[self.mixers[self.names-1] == mixer]

    def get_signal(self, group, name, calib=False, tmin=None, tmax=None):

        if tmin is None:    tmin = self.tmin
        if tmax is None:    tmax = self.tmax
            
        if np.size(name) > 1:
            data = [self.get_signal(group, n, calib=calib, tmin=tmin, tmax=tmax) for n in name]
            return data

        n = int(name)-1
        sig_name = 'Trad-A1' if n < 30 else 'Trad-A2'
        
        #cache loaded Te data
        if n < 30:
            if not hasattr(self, 'TradA1'):
                if self.shotfile == 'RMC':
                    rmc = sf.SFREAD(self.shotfile, self.shot, experiment=self.exp, edition=self.ed)
                    self.tvecA1 = rmc.gettimebase(sig_name)
                    self.TradA1 = rmc.getobject(sig_name, cal=True)
                    
                if self.shotfile in ('RMA', 'RMB', 'RAD') :
                    if self.shotfile == 'RMA':
                        self.shotfile = 'RMB'
                    sfo = sf.SFREAD(self.shotfile, self.shot, experiment=self.exp, edition=self.ed)
                    self.tvecA1 = sfo.gettimebase('SI-MI-A')
                    self.TradA1 = np.hstack((sfo.getobject('SI-MI-A', cal=True), sfo.getobject('SI-MI-B', cal=True)))
            tvec = self.tvecA1
            sig = self.TradA1[:, n]

        else:
            if not hasattr(self, 'TradA2'):
                if self.shotfile == 'RMC':
                    rmc = sf.SFREAD(self.shotfile, self.shot, experiment=self.exp, edition=self.ed)
                    self.tvecA2 = rmc.gettimebase(sig_name)
                    self.TradA2 = rmc.getobject(sig_name, cal=True)
                    
                if self.shotfile in ('RMA', 'RMB', 'RAD') :
                    if self.shotfile == 'RMA':
                        self.shotfile = 'RMB'
                    sfo = sf.SFREAD(self.shotfile, self.shot, experiment=self.exp, edition=self.ed)
                    self.tvecA2 = sfo.gettimebase('SI-MI-C')
                    self.TradA2 = np.hstack((sfo.getobject('SI-MI-C', cal=True), sfo.getobject('SI-MI-D', cal=True)))
            tvec = self.tvecA2
            sig = self.TradA2[:, n-30]

        ioff, imin, imax = tvec.searchsorted((0, tmin, tmax))
        tvec = tvec[imin: imax].astype(np.float32)
        sig = sig[imin: imax].astype(np.float32)
        logger.debug('tvec-size: %d, sig-size: %d', len(tvec), len(sig))
        return [tvec, sig]


    def get_Te0(self, tmin, tmax, IDA=True, dR=0, dZ=0):
        

        time = (tmin + tmax)/2
        rho = self.get_rho('', self.names, time, dR=dR, dZ=dZ)[0]
        if IDA:
            idasf = sf.SFREAD('IDA', self.shot)
            if idasf.status:
                ida_tvec = idasf('time')
                
                if tmin > ida_tvec[-1]:
                    return
                    
                imin, imax = ida_tvec.searchsorted((tmin, tmax))

                Te = idasf('Te')[:, imin-1:imax].mean(-1)
                rho_ida = idasf('rhop')[:, imin]
                rho_ida = sf.rho2rho(self.eqm, rho_ida, time)[0]
                TeEce2 = np.interp(abs(rho), rho_ida, Te)  #dR, dZ are not accounted 
   
                return rho, TeEce2

# ugly robust fit of ECE data
        cec = sf.SFREAD('CEC', self.shot)
        if not cec.status:
            return

        tvec = cec.gettimebase('Trad-A')
        imin, imax = tvec[:-1].searchsorted((tmin, tmax))

        Te = cec.getobject('Trad-A', nbeg=imin, nend=imax+1)
        Te = np.nanmedian(Te[:,self.names-1], 0)
        from scipy.signal import medfilt
  
        from scipy.interpolate import LSQUnivariateSpline

        ind = np.argsort(rho)
        ind = ind[(np.abs(rho[ind]) < 1) & (Te[ind] > 0)]

        x = np.r_[rho[ind],1]
        y = np.log(np.r_[medfilt(Te[ind],5),min(Te[ind].min(),100)])
        ind_ = np.argsort(np.r_[-x,x])
        
        try:
            S   = LSQUnivariateSpline(np.r_[-x,x][ind_],np.r_[y,y][ind_],np.linspace(-0.999,0.999,21),ext=3, bbox=[-1,1],k=2)
            c = (S(x)-y)/S(x) < 0.005
            x,y = x[c], y[c]
            ind_ = np.argsort(np.r_[-x,x])
            S   = LSQUnivariateSpline(np.r_[-x,x][ind_], np.r_[y,y][ind_],np.linspace(-0.999,0.999,21),ext=3, bbox=[-1,1],k=2)
            Te_ = np.exp(S(rho))

            if not all(np.isfinite(Te_)):
                raise Exception()

        except Exception as e:
            logger.warning('spline fit has failed: %s', str(e))
            Te_ = np.exp(np.interp(rho, np.r_[-x,x][ind_], np.r_[y,y][ind_]))
        
        return rho, Te_

    def get_names_phase(self):
        pass
    def get_signal_phase(self,name,calib=False):
        pass
    def get_phi_tor(self,name=None):
        return np.deg2rad(300+40.5) # position of SXR and relative distance between SXR anbd ECE #BUG what abiut these 45deg in SXR loader? 
    def get_phase_corrections(self,name):
        pass
        
    
    def get_RZ_theta(self, time,names,dR=0,dZ=0,cold=False):

        if not hasattr(self,'RZtime'):
            err = np.ones(len(names))*np.nan
            return err,err,err

        time = np.clip(time,*self.RZtime[[0,-1]])
        ch_ind = np.in1d(self.names, np.int_(names))
        
        R = interp1d(self.RZtime, self.R[:,np.array(names)-1],axis=0)(time)
        z = interp1d(self.RZtime, self.z[:,np.array(names)-1],axis=0)(time)

        r0 = np.interp(time, self.eqm.time, self.eqm.Rmag)+dR
        z0 = np.interp(time, self.eqm.time, self.eqm.Zmag)+dZ
   
        return R,z, np.arctan2(z-z0, R-r0)
        
    
    def get_rho(self,group,names,time,dR=0,dZ=0):

        if hasattr(self,'RZtime'):
            time = np.clip(time, self.tmin,self.tmax)
            R,z,theta = self.get_RZ_theta(time,names,dR=dR,dZ=dZ)
            
            rho = super(loader_ECE,self).get_rho(time,R,z,dR=dR,dZ=dZ)

            r0 = np.interp(time, self.eqm.time, self.eqm.Rmag)+dR
            R = np.atleast_1d( R)
            rho[R < r0] *= -1
        else:
            err = np.ones(len(names))*np.nan
            rho,theta,R,z = err,err,err,err
        
        return rho, theta, R, z
    
        
    def signal_info(self,group,name,time):
        name = int(name)
        rho,theta,R,z = self.get_rho(group,[name,],time)

        info = 'ch: '+str(name)+'  R:%.3f m   '%R+self.rho_lbl+': %.3f'%rho
        return info
    
    def get_description(self,group,name):
        return 'AUG %d diag: %s sig: %s'%(self.shot,self.shotfile,name)



