from .loader import * 

logger = logging.getLogger('pyspecview.sxr')
logger.setLevel(logging.INFO)


def check(shot):
    #fastest check if the shotfile exist
    status = False

    path = shot_path+'/%d/SX/SXA/%d'
    status |= os.path.isfile(path%(shot/10,shot))

    return status


#not implemented for the old SXR system 
class loader_SXR(loader):

    tor_mode_num = True
    pol_mode_num = False  #not implemented yet!! 

    mode_range = (-3,4)
    radial_profile=True
    units = 'W/m$^2$'


    def __init__(self, *args, **kargs):

        import aug_sfutils as sf
        super(loader_SXR, self).__init__(*args, **kargs)

        self.calib_shot = sf.previousshot('CSX', self.shot, exp='AUGD')

        csx = sf.SFREAD('CSX', self.calib_shot)
        names = csx.objects + csx.parsets

        signals = [b[1:] for b in names if b[0] == 'C']

        shotfiles = []
        self.Phi      = {}
        self.R_start  = {}
        self.z_start  = {}
        self.R_end    = {}
        self.z_end    = {}
        self.status   = {}
        different_det = {}
        self.ADCrange = {}

        for s in signals:
            par = csx('C' + s)
            shotfiles.append(sf.str_byt.to_str(  b''.join(par['SX_DIAG']) ) )
            self.Phi[s] = par['Tor_Pos'] + 45
        #BUG missing phi_end for T camera!!
            self.R_start[s] = par['RPINHOLE']
            self.z_start[s] = par['ZPINHOLE']
            self.R_end[s]   = par['REND'    ]
            self.z_end[s]   = par['ZEND'    ]
            self.status[s]  = par['ADDRESS' ] !=256
            thickness    = par['THICKNES']
            filt_mat     = sf.str_byt.to_str(par['FILT-MAT'].item())
            different_det[s] = (abs(thickness - 75e-6) > 1e-5) | (filt_mat != 'Be')
            logger.debug('%s %12.4e %s %s', s, thickness, filt_mat, different_det[s])
            self.ADCrange[s]  = par['ADCrange']

        self.ADCmin = 0

        self.SXR_diods = {}
        for sf in np.unique(shotfiles):
            self.SXR_diods[sf] = []
            
        for d, s in zip(shotfiles, signals):
            self.SXR_diods[d].append(s)
            
        try:
            self.SXR_diods.pop('OOO')  #null detector
        except:
            pass
        
        self.all_signals = []
        for signals in self.SXR_diods.values():
            self.all_signals.extend(signals)
        self.groups = np.unique([i[0] for i in self.all_signals])
        self.openshotfile = ''

    def get_names(self, group):

        logger.debug('get_names')
        names = [s for s in self.all_signals if s[0] == group]
        names.sort()
        return names

            
    def get_signal(self, group, name, calib=False, tmin=None, tmax=None):
        
        logger.debug('get_signal')
        if tmin is None:    tmin = self.tmin
        if tmax is None:    tmax = self.tmax
        
        if isinstance(name,str):
            names = (name, )
        else:  names = name
        
        output = np.empty(len(names), dtype=object)
        for shotfile,signals in self.SXR_diods.items():
            
            for j, name in enumerate(names):
                if name in signals :
                    if self.openshotfile != shotfile:
                        sfo = sf.SFREAD(shotfile, self.shot, experiment=self.exp, edition=self.ed)

                    self.ed = sfo.ed
                    self.openshotfile = shotfile
                    
                    tvec = sfo('Time').astype(np.float32)
                    sig = sfo.getobject(name, copy=True) #*1. # git create a copy of read-only array

# remove corrupted points
                    wrong =  (sig < self.ADCmin) |  (sig > self.ADCrange[name]) 

                    if any(wrong):
                        sig[wrong]=np.interp(np.where(wrong)[0], np.where(~wrong)[0],sig[~wrong])
                    if calib:
                        sig = sfo.raw2calib(sig)

                    output[j] = [tvec, sig.astype(np.float32)]

        if len(names) == 1:
            return output[0]

        return output


    def get_names_phase(self):

        logger.debug('get_names_phase')
        Fsig = self.get_names('F')
        Gsig = self.get_names('G')
        
        FGnames = []
        for f in Fsig:
            for g in Gsig:
                if f[1:]==g[1:]:  FGnames.append('FG'+f[1:])

        return FGnames
        
    
    def get_signal_phase(self,name,calib=False,tmin=None,tmax=None):
        
        logger.debug('get_signal_phase')
        if name[:2]!= 'FG':    raise Exception()
      
        tvec1, sig1 = self.get_signal('F','F'+name[2:],calib=calib,tmin=tmin,tmax=tmax)
        tvec2, sig2 = self.get_signal('G','G'+name[2:],calib=calib,tmin=tmin,tmax=tmax)

        #downsample  to the slower DAS
        def reduce(x,y,x_out):
            r = int(round(float(len(x))/len(x_out)))
            x = np.mean(x[:(len(x)//r)*r].reshape(len(x)//r, r),1)
            y = np.mean(y[:(len(y)//r)*r].reshape(len(y)//r, r),1)
            return np.interp( x_out, x, y)

        if len(sig1)!= len(sig2):
            if len(sig1) > len(sig2):
                sig1,tvec1 = reduce(tvec1,sig1,tvec2),tvec2
            else:
                sig2,tvec2 = reduce(tvec2,sig2,tvec1),tvec1

        return tvec1, np.single(np.vstack((sig1,sig2)).T)
    
    
    def get_phi_tor(self,name=None):

        logger.debug('get_phi_tor')
        if name in self.get_names_phase():
            
            phi1 = self.Phi['F'+name[2:]]
            phi2 = self.Phi['G'+name[2:]]
            
            return np.deg2rad( np.r_[phi1,phi2])
        elif name in self.Phi:
            return np.deg2rad(self.Phi[name])
        
        else: 
            try:
                return np.deg2rad(np.median([v for k,v in self.Phi.items()]))
            except:
                logger.error(extra=self.Phi)
                raise
            
    
    def get_rho(self, group, names, time, dR=0, dZ=0):

        logger.debug('get_rho')
        #BUG not working for a tangential camera!!!
        R_start = np.array([self.R_start[name] for name in names])
        z_start = np.array([self.z_start[name] for name in names])
        R_end = np.array([self.R_end[name] for name in names])
        z_end = np.array([self.z_end[name] for name in names])
        Phi = np.array([self.Phi[name] for name in names])
        rho_tg,theta_tg,R,Z = super(loader_SXR,self).get_rho(time,R_start,
                                    z_start,Phi,R_end,z_end,Phi,dR=dR, dZ=dZ)

        return rho_tg, theta_tg,R,Z

    
    def signal_info(self,group,name,time):
        
        logger.debug('signal_info')
        if name[:2] == 'FG': 
            rho_tg1 = self.get_rho(group,[ 'F'+name[2:],],time)[0]
            rho_tg2 = self.get_rho(group,[ 'G'+name[2:],],time)[0]

            phi1 = self.Phi['F'+name[2:]]
            phi2 = self.Phi['G'+name[2:]]
            
            info = str(name)+'   Phi: %d and %d deg, rho_tg: %.2f'%(phi1,phi2,rho_tg1)

        else:
            rho_tg = self.get_rho(group,[ name,],time)[0]
            phi = self.Phi[name]

            info = str(name)+' Phi: %d deg, rho_tg: %.2f'%(phi,rho_tg)

        return info
