from loaders_DIIID.loader import * 
import os
from multiprocessing import  Pool

from time import time as T

import MDSplus as mds

def check(shot):
    #fastest check if the shotfile exist
    #BUG 
    return True



def mds_load(tmp):
    (mds_server,  TDI) = tmp
    MDSconn = mds.Connection(mds_server)
    data = []
    for tdi in TDI:
        try:
            data.append(MDSconn.get(tdi).data())
            if data[-1].dtype == 'int32':
                #print data[-1].max(),
                data[-1] = data[-1].astype('int16',copy=False)
                #print data[-1].max()

            #print tdi
        except:
            data.append([])
   
    return data


class loader_ICE(loader):

    def __init__(self,*args, **kargs):
        
        super(loader_ICE,self).__init__(*args, **kargs)

        self.names = {'ICE': list(range(1,6))}
        self.groups = ('ICE',)


        self.tvec = None
        self.catch = {}


    
    def get_names(self,group):
        return self.names[group]
    
    def get_signal(self,group, names,calib=False,tmin=None,tmax=None):
        #WARNING data are really huge!!
        
        if tmin is None:    tmin = self.tmin
        if tmax is None:    tmax = self.tmax
        

        if size(names) >1:
            return [self.get_signal(group, name,calib=calib,tmin=tmin,tmax=tmax) for name in names]
        name = names
        
   
        if name in self.catch:
            ind = self.tvec.searchsorted(tmin,tmax)
            return self.tvec[ind], self.catch[name][ind]

        
    
        print( '\n fast parallel fetch... be patient!! it can take few minutes')
        numTasks = 8
        n_parts = 62
        #n_parts = 7

        t = T()
        server = self.MDSconn.hostspec
        TDIcall = 'ptdata2("ice%.2d_%s",%d,0)'
        
        TDI = [TDIcall%(name, n, self.shot) for n in range(1,n_parts+1)]
        if self.tvec is None:  
            TDI.append('dim_of('+TDI[0]+')')
        
        TDI = array_split(TDI, numTasks)


        args = [(server, tdi) for tdi in TDI]
             

        pool = Pool()
        out = pool.map(mds_load,args)
      
        pool.close()
        pool.join()
        
        out_ = []
        for o in out:
            out_.extend(o)
        
        

        
        
        print(( 'data loaded in %.2f'%( T()-t)))
        
        
        if self.tvec is None:
            tvec = out_[ -1]
            out_  = out_[:-1]
            if len(tvec) == 0: raise Exception('No ICE data!')
            dt = (tvec[-1]-tvec[0])/(len(tvec)-1)
            N = sum([o.size for o in out_])
            #fast object emulating equally spaced vector
            self.tvec = spaced_vector(tvec[0]/1e3,(tvec[0]+dt*N)/1e3,dt/1e3)
                  #sum(out)

        ind = slice(*self.tvec.searchsorted([tmin,tmax]))
        sig = hstack(out_) 
        del out_
        
        self.catch[name] = sig

          
        return self.tvec[ind], sig[ind]

    


    def signal_info(self,group,name,time):
        return ''
    
 
def main():
    
    
    mds_server = "localhost"
    mds_server = "atlas.gat.com"
    import MDSplus as mds

    MDSconn = mds.Connection(mds_server )
 
    ice = loader_ICE(169132,exp='DIII-D',MDSconn=MDSconn)
    #cd 
    ice.get_signal('',1, )
    

if __name__ == "__main__":
    main()
    
