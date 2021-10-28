from .loader import * 

logger = logging.getLogger('pyspecview.mirnov')
logger.setLevel(logging.INFO)

shotfiles = 'MHI', 'MHA', 'MHB', 'MHC', 'MHD', 'MHE', 'MHF', 'MHG', 'MHH', 'MTR'

def check(shot):
    #fastest check if the shotfile exist

    for shotfile in shotfiles:
        path = shot_path+str(shot//10)+'/MH/'+shotfile+'/'+str(shot)
        if os.path.isfile(path):
            return True 
                
    path = shot_path+str(shot//10)+'/L0/'+'MTR'+'/'+str(shot)
    if os.path.isfile(path):
        return True 

    return False


class loader_MirnovCoils(loader):


    mode_range = (-6, 5)
    pol_mode_num = True

    def __init__(self, *args, **kargs):

        super(loader_MirnovCoils, self).__init__(*args, **kargs)
        
        self.old_shotfile = self.shot < 19101

        if self.old_shotfile:
            mtr = sf.SFREAD('MTR', self.shot)
            if mtr.status:
                names = mtr.objects + mtr.parsets
                names = [n[1:] for n in names if n[:2] == 'PC' and n[-3] == '-']
        elif self.shot < 21496:
            mha = sf.SFREAD('MHA', self.shot)
            if mha.status:
                names = mha.objects + mha.parsets
                names = [n[1:] for n in names if n[:2] == 'CC' and n[-3] == '-']
                self.old_shotfile = True
        else:
            calib_shot = sf.previousshot('CMH', self.shot, exp='AUGD')
            cmh = sf.SFREAD('CMH', calib_shot)
            names = cmh.parsets
            names = [n[1:] for n in names if n[:2] == 'CC' and n[-3] == '-']

        names.sort()

        self.groups = np.unique([n[:3] for n in names])
        self.names = {g:[] for g in self.groups  }
        
        for shotfile in shotfiles:
            path = shot_path+str(self.shot//10)+'/MH/'+shotfile+'/'+str(self.shot)
            if shotfile == 'MTR':
                path = shot_path+str(self.shot//10)+'/L0/MTR/'+str(self.shot)
    
            if not os.path.isfile(path):
                continue
            
            sfo = sf.SFREAD( shotfile, self.shot, experiment=self.exp, edition=self.ed)

            sf_names = sfo.objects
            sf_names.sort()

            for n in sf_names:
                if n in names:
                    self.names[n[:3]].append((shotfile, n))
                    names.remove(n)  #do not include the same coil twice, if it is in more different shotfiles
    
    def get_names(self, group):
        names = [n for sf, n in self.names[group]]
        names.sort()
        return names  
    
    def get_signal(self, group, name, calib=False, tmin=None, tmax=None):
        
        if tmin is None:    tmin = self.tmin
        if tmax is None:    tmax = self.tmax

        sig_names = self.names[group]
        for self.shotfile, sig in sig_names:
            if sig == name:
                break
        
        sfo = sf.SFREAD(self.shotfile, self.shot, experiment=self.exp, edition=self.ed)
        tvec = sfo.gettimebase(name) 
        nbeg, nend = tvec.searchsorted((tmin, tmax))

        sig = sfo.getobject(name, cal=calib, nbeg=nbeg, nend=nend)
        return tvec[nbeg:nend+1], sig
        
    def get_names_phase(self):
        return ['Mcoils', ]
    
    def get_signal_phase(self, name, calib=False, tmin=None, tmax=None, rem_elms=True):
        if tmin is None:    tmin = self.tmin
        if tmax is None:    tmax = self.tmax

        phase_signals = self.names['C09']

        sigs = []
        length = np.infty
        self.shotfile = None
        for sfile, sig_p in phase_signals:
            if self.old_shotfile  and sfile == 'MHA': continue
            if self.shotfile != sfile:
                sfo = sf.SFREAD(sfile, self.shot, experiment=self.exp, edition=self.ed)
                self.shotfile = sfile
                tvec = sfo.gettimebase(sig_p)
                nbeg, nend = tvec.searchsorted((tmin, tmax))

            sig = sfo.getobject(sig_p, cal=calib, nbeg=nbeg, nend=nend)
            length = min(length, len(sig))
            sigs.append(sig)

        data = np.vstack([s[:length] for s in sigs]).T
        tvec = tvec[nbeg:nend+1]
        
        #remove elms from the data 
        logger.info('Removing ELMs')
        elm = sf.SFREAD('ELM', self.shot)
        if rem_elms and elm.status:
            nbeg = tvec.searchsorted(elm('t_begELM') - 0.001)
            nend = tvec.searchsorted(elm('t_endELM'))
            
            for n1, n2 in zip(nbeg, nend):
                n1, n2 = max(1, n1), min(n2, len(tvec) - 2)
                data[n1: n2] = np.linspace(0, 1, n2-n1)[:, None] * \
                    (np.float_(data[n2+1]) - data[n1-1]) + data[n1-1]
   
        logger.info('Done')
        return tvec, data


    def get_theta_pol(self, name, tshot = 4 , rhop=.2 ):

        import aug_sfutils as sf # for call cross-phaseogram

        try:
            magr, magz, theta0, theta_star = self.mag_theta_star( tshot, rhop )
        except:
            logger.error( 'Error: get_theta_pol: mag_theta_star do not work!')
            logger.error( traceback.format_exc())

            theta0     = np.linspace(-np.pi, np.pi, 10)
            theta_star = np.linspace(-np.pi, np.pi, 10)

        if name in  ['Mcoils', ]:
            theta = []
            
            if self.shot > 19090:
                calib_shot = sf.previousshot('CMH', self.shot, exp='AUGD')
                sfo = sf.SFREAD('CMH', calib_shot)
                pre = 'C'
            else:
                sfo = sf.SFREAD('MTR', self.shot)
                pre = 'P'

            for sf, sig in self.names['C09']:
                theta.append(sfo(pre+sig)['theta'])

            theta = np.interp(theta,  theta0, theta_star)

        return -np.squeeze(theta) #BUG minus is there to have a positive m numbers for ordinary MHD modes with coinjected NBI 
            

    def signal_info(self, group, name, time):
        if name == 'Mcoils' or self.old_shotfile: return ' '

        calib_shot = sf.previousshot('CMH', self.shot, exp='AUGD')
        cmh = sf.SFREAD('CMH', calib_shot)
            
        theta = cmh('C'+name)['theta']

        try:
            info = str(name) + ' theta: %.2fdeg' %np.rad2deg(theta)
        except:
            print(group, name, time, theta)
            raise
        
        return info
