from .loader import * 

logger = logging.getLogger('specview.elm_c')
logger.setLevel(logging.INFO)

shotfiles = 'MHI', 'MHA', 'MHB', 'MHC', 'MHD', 'MHE', 'MHF', 'MHG', 'MHH'
 
def check(shot):
    """
    Check if any shotfile exists
    """

    for shotfile in shotfiles:
        path = shot_path+str(shot//10)+'/MH/'+shotfile+'/'+str(shot)
        if os.path.isfile(path):
            return True 

    return False


class loader_BalooningCoils(loader):
    tor_mode_num = True
    mode_range = (-6, 5)
    
    #BUG 13 is broken in new discharges!!!!!!!!!!!!!
    phase_balooning_coils = [1, 2, 3, 12, 13, 14]


    def __init__(self, *args, **kargs):

        super(loader_BalooningCoils, self).__init__(*args, **kargs)

        if self.shot > 21496:
            calib_shot = sf.previousshot('CMH', self.shot, exp='AUGD')
            sfo = sf.SFREAD('CMH', calib_shot)
        else:
            sfo = sf.SFREAD('MHA', self.shot)
    
        sig_names = ['CB31-%.2d'%c for c in self.phase_balooning_coils]
        self.phi = [sfo(n)['phi'] for n in sig_names]
        names = sfo.objects + sfo.parsets
        names =  [n[1:] for n in names if n[:2] == 'CB' and n[-3] == '-']

        self.groups = np.unique([n[:3] for n in names])
        self.names = {g:[] for g in self.groups  }

        for shotfile in shotfiles:

            path = shot_path+str(self.shot//10)+'/MH/'+shotfile+'/'+str(self.shot)
    
            if not os.path.isfile(path):
                continue
            
            sfo = sf.SFREAD(shotfile, self.shot, experiment=self.exp, edition=self.ed)
            sf_names = sfo.objects + sfo.parsets
            sf_names.sort()
            
            for n in sf_names:
                if n in names:
                    self.names[n[:3]].append((shotfile, n))

   
    def get_signal(self, group, name, calib=False, tmin=None, tmax=None):
        
        if tmin is None:    tmin = self.tmin
        if tmax is None:    tmax = self.tmax
           
        sig_names = self.names[group]
        for self.shotfile, sig in sig_names:
            if sig == name: break
         
        sfo = sf.SFREAD(self.shotfile, self.shot, experiment=self.exp, edition=self.ed)
        tvec = sfo.gettimebase(name)
        nbeg, nend = tvec.searchsorted((tmin, tmax))

        sig = sfo.getobject(name, cal=calib, nbeg=nbeg, nend=nend)

        return tvec[nbeg:nend+1], sig

    
    def get_names(self, group):
        names = [n for sf, n in self.names[group]]
        names.sort()
        return names  

        
    def get_names_phase(self):
        return ['BCoils', ]


    def get_signal_phase(self, name, calib=False, tmin=None, tmax=None):
        
        shotfiles = [sf for sf, n in self.names['B31']]
        
        if 'MHI' in shotfiles:
            self.shotfile = 'MHI'
        else:
            self.shotfile = 'MHA'

        if tmin is None:    tmin = self.tmin
        if tmax is None:    tmax = self.tmax

        sfo = sf.SFREAD(self.shotfile, self.shot, experiment=self.exp, edition=self.ed)
        tvec = sfo.gettimebase('B31-%.2d'%self.phase_balooning_coils[0])
        nbeg, nend = tvec.searchsorted((tmin, tmax))

        if name in  ['BCoils', ]:
            sigs = [sfo.getobject('B31-%.2d' %c, cal=calib, nbeg=nbeg, nend=nend) for c in self.phase_balooning_coils]

        return tvec[nbeg:nend+1], np.vstack(sigs).T


    def get_phi_tor(self, name):
        #works only for balooning coils yet

        return np.squeeze(self.phi)
            

    def get_phase_corrections(self, name):
        path = os.path.abspath(__file__)
        path = path[:path.rfind('/')]
        return np.loadtxt(path+'/balooning_coils_correction.txt')


    def signal_info(self, group, name, time):
        
        if name == 'BCoils': return ' '
        
        if self.shot > 19101:
            calib_shot = sf.previousshot('CMH', self.shot, exp='AUGD')
            sfo = sf.SFREAD('CMH', calib_shot)
        else:
            sfo = sf.SFREAD('MHA', self.shot)

        parset = sfo('C' + name)
        theta = parset['theta']
        phi   = parset['phi']

        try:
            info = str(name)+' theta: %.2fdeg, phi: %.2fdeg'%(np.rad2deg(theta), np.rad2deg(phi))
        except:
            
            print(theta, phi, name)
            raise
        return info
