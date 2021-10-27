from .loader import * 

logger = logging.getLogger('pyspecview.heb')
logger.setLevel(logging.INFO)


def check(shot):
    #fastest check if the shotfile exist

    status = False
    
    path = shot_path+'/%d/YY/HEB/%d'
    status = os.path.isfile(path%(shot//10,shot))
    
    return status


class loader_HEB(loader):

    def __init__(self,*args, **kargs):
        
        super(loader_HEB,self).__init__(*args, **kargs)

        self.shotfile = 'HEB'
        heb = sf.SFREAD(self.shotfile, self.shot)
        if heb.status:
            LOS = heb('INFO')['LOS']
            self.groups = ['PMT%d'%i for i in range(1,5)]
            self.names = {g:LOS.tolist() for g in self.groups}
  

    def get_names(self,group):
        return self.names[group]
    
    def get_signal(self,group, name,calib=False,tmin=None,tmax=None):

        if tmin is None:    tmin = self.tmin
        if tmax is None:    tmax = self.tmax

        ich = self.names[group].index(name)
        
        #load data only once and than keep then in the memory
        if not hasattr(self, group):
            sfo = sf.SFREAD(self.shotfile , self.shot, experiment=self.exp, edition=self.ed)
        
            if not hasattr(self, 'tvec'):
                print('HEB ', group)
                self.tvec = sfo.gettimebase(group)

            setattr(self, group, sfo(group))
        
        nbeg, nend = self.tvec.searchsorted((tmin,tmax))
        signal = getattr(self,group)

        return self.tvec[nbeg:nend+1],signal[nbeg:nend+1,ich]


    def signal_info(self,group,name,time):
        info = sf.str_byt.to_str(group) + ': ' + sf.str_byt.to_str(name)
        return info
