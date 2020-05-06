from .loader import * 
import os


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
        if self.dd.Open(self.shotfile,self.shot):
            LOS = self.dd.GetParameter('INFO', 'LOS')  
            self.dd.Close()   
            self.groups = ['PMT%d'%i for i in range(1,5)]
            self.names = {g:LOS.tolist() for g in self.groups}
  

    def get_names(self,group):
        return self.names[group]
    
    def get_signal(self,group, name,calib=False,tmin=None,tmax=None):
        

        if tmin is None:    tmin = self.tmin
        if tmax is None:    tmax = self.tmax
        
        ich = self.names[group].index(name)
        
        #load data only once and than keep then in the memory
        if not hasattr(self,group):
            self.dd.Open(self.shotfile , self.shot, experiment=self.exp, edition=self.ed)
        
            if not hasattr(self,'tvec'):
                self.tvec = self.dd.GetTimebase(group)

            setattr(self,group,self.dd.GetSignal(group))
        self.dd.Close()
        
        nbeg, nend = self.tvec.searchsorted((tmin,tmax))
        signal = getattr(self,group)

        return self.tvec[nbeg:nend+1],signal[nbeg:nend+1,ich]

            

    def signal_info(self,group,name,time):
        info = group+': '+name
        return info
    
 
 
def main():
    import os,sys
    sys.path.append('/afs/ipp/home/t/todstrci/TRANSP/')
    from . import dd   
    #from matplotlib.pylab import *
    from . import map_equ
    #import IPython
    #IPython.embed()
    
    shot =  34874


    eqm = map_equ.equ_map(debug=True)
    eqm.Open(shot, diag='EQI')
    #eqm.read_ssq()
    heb = loader_HEB(shot,eqm= eqm,rho_lbl='rho_tor')
    ##sxr =  loader_SXR(33871)
    #G = sxr.get_signal_groups()
    S = heb.get_signal('PMT1','HEB-4+0 ',tmin=3,tmax=3.3)
    
    #t,x = heb.get_signal_phase('FG_020')
    
  
if __name__ == "__main__":
    main()

