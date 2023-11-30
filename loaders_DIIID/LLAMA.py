try:
    from loaders_DIIID.loader import * 
except:
    from loader import * 
 

def check(shot):
    return True



class loader_LAMMA(loader):

    def __init__(self,*args, **kargs):
        
        super(loader_LAMMA,self).__init__(*args, **kargs)
        
        n_los = 20
        name = 'LYA1%s%.2d'
        self.groups =  ['LFS','HFS']
        self.names = {}
        
        for g in self.groups:
            self.names[g]  = [name%(g[0],n+1) for n in range(n_los)]
      
        self.tvec = None
        self.catch = { }
 
    def get_names(self,group):
        return self.names[group]
    
    def get_signal(self,group, name,calib=False,tmin=None,tmax=None):
        

        PTNAME = 'PTDATA2("%sRAW",%d,1)'%(name, self.shot) 
        
        if name not in self.catch:
            #from IPython import embed
            #embed()
            self.catch[name] = self.MDSconn.get(PTNAME)
     
            if len(self.catch[name]) == 1:
                raise Exception('No LLAMA data')
            
            if self.tvec is None:
                self.tvec = self.MDSconn.get(f'dim_of({PTNAME})')/1e3
    
    
        if tmin is None:    tmin = self.tmin
        if tmax is None:    tmax = self.tmax
        
        
        imin,imax = self.tvec.searchsorted([tmin,tmax])
        
        
        return self.tvec[imin:imax], self.catch[name][imin:imax]
    
      
  


    def signal_info(self,group,name, time):
        return 'sig:%s %.2d   %s:%.2f  '%(group,name )
    
 
 

     
if __name__ == "__main__":
    main()
    
    
    
    
    
