from .loader import * 
import os

shotfiles = 'FHA', 'FHB', 'FHC'

def check(shot):
    #fastest check if the shotfile exist

    status = False
    
    for s in shotfiles:
        path = shot_path+'/%d/XX/%s/%d'
        status |= os.path.isfile(path%(shot//10,s,shot))
    
    
    return status





class loader_FILD(loader):

    def __init__(self,*args, **kargs):
        
        super(loader_FILD,self).__init__(*args, **kargs)


        self.names = {}
        for shotfile in shotfiles:
            if self.dd.Open(shotfile,self.shot):
                names = self.dd.GetNames()  
                self.dd.Close()   
                names = [n for n in names if n[:2] == 'FI']
                self.names[shotfile] = names
                
        self.groups = list(self.names.keys())

    def get_names(self,group):
        return self.names[group]
    
    def get_signal(self,group, name,calib=False,tmin=None,tmax=None):
        

        if tmin is None:    tmin = self.tmin
        if tmax is None:    tmax = self.tmax
            
        self.shotfile =  group
        self.dd.Open(group, self.shot, experiment=self.exp, edition=self.ed)
        tvec = self.dd.GetTimebase(name)
        nbeg, nend = tvec.searchsorted((tmin,tmax))

        sig = self.dd.GetSignal(name,cal=calib, nbeg= nbeg,nend = nend)
        #self.dd.Close()

        return tvec[nbeg:nend+1],sig

            

    def signal_info(self,group,name,time):
        info = group+': '+name
        return info
    
 
