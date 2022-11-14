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
    MDSconn = mds.Connection(mds_server )
    data = [MDSconn.get(tdi).data() for tdi in TDI]
    
    #print TDI
    return data



class loader_BES(loader):
    radial_profile=True

    def __init__(self,*args, **kargs):
        
        super(loader_BES,self).__init__(*args, **kargs)

        self.names = {'BESFU': list(range(1,65)),'BESSU': list(range(1,65))}
        self.groups =  list(self.names.keys())


        self.tvec =  {'BESFU':None, 'BESSU':None }
        self.catch = {'BESFU':{}, 'BESSU':{} }
        self.MDSconn.openTree('BES', self.shot)
        try:
            self.R = self.MDSconn.get('\\bes::bes_r').data()/100.
            self.Z = self.MDSconn.get('\\bes::bes_z').data()/100.
        except:
            print('BES R,Z coordinates were not found')
            self.R = 0
            self.Z = 0
        self.MDSconn.closeTree('BES', self.shot)

    
    def get_names(self,group):
        return self.names[group]
    
    def get_signal(self,group, names,calib=False,tmin=None,tmax=None):
        
        
        if tmin is None:    tmin = self.tmin
        if tmax is None:    tmax = self.tmax
        
        names = atleast_1d(names)
        ch2load = []
        for name in names:
            if not name in self.catch[group]:
                ch2load.append(name)
            
        TDIcall ='_x=PTDATA2("%s%.2d",%d,0)'    
        
        
        if  size(ch2load) == 1:
        
            TDIcall = TDIcall%(group,ch2load[0],self.shot)
            #print TDIcall
            sig = self.MDSconn.get(TDIcall).data()
            if len(sig) < 2:
                raise Exception(group+' channel %d do not exist'%name)

            self.catch[group][name] = sig
            if self.tvec[group] is None:
                self.tvec[group] = self.MDSconn.get('dim_of(_x)').data()
                self.tvec[group] /= 1e3 #s
                
                        
            #imin,imax = self.tvec[group].searchsorted([tmin,tmax])
            #ind = slice(imin,imax+1)
            #return self.tvec[group][ind], sig[ind]
                    
        elif size(ch2load) > 1:
            print( '\n fast parallel fetch...' )
            numTasks = 8
                            
            t = T()
            server = self.MDSconn.hostspec
            TDI = [TDIcall%(group,n,self.shot) for n in ch2load]
                        
            TDI = array_split(TDI, numTasks)
            ch2load = array_split(ch2load, numTasks)

            if self.tvec[group] is None:  
                TDI[-1] = append(TDI[-1],'dim_of(_x)' )
            
            args = [(server, tdi) for tdi in TDI]
            
      
            pool = Pool()
            out = pool.map(mds_load,args)
            pool.close()
            pool.join()
            print(( 'in %.1fs'%(T()-t)))
    
            if self.tvec[group] is None:  
                self.tvec[group] = out[-1][-1]/1000#s
        
            for chnls, data in zip(ch2load, out):
                for ch,sig in zip(chnls, data):
                    self.catch[group][ch] = sig
         
        imin,imax = self.tvec[group].searchsorted([tmin,tmax])
        ind = slice(imin,imax+1)
        
        
        if len(names) == 1:
            return self.tvec[group][ind], self.catch[group][names[0]][ind]

        
        return [[self.tvec[group][ind], self.catch[group][n][ind]] for n in names]

        

    
    def get_rho(self,group,names,time,dR=0,dZ=0):
        ind = asarray(names)-1
        return  super(loader_BES,self).get_rho(time,self.R[ind],self.Z[ind],dR=dR, dZ=dZ)
    


    def signal_info(self,group,name, time):
        
        rho_tg = self.get_rho(group,[name,],time)[0]
        
        R,Z = self.R[name-1],self.Z[name-1]
        return 'sig:%s %.2d   %s:%.2f R:%.3fm z:%.3fm'%(group,name, self.rho_lbl,rho_tg,R,Z)
    
 
 
def main():
    

    mds_server = "localhost"
    #mds_server = "atlas.gat.com"

    import MDSplus as mds
    MDSconn = mds.Connection(mds_server )
    from .map_equ import equ_map
    eqm = equ_map(MDSconn,debug=False)
    
    shots = r_[175849:175869, 175882:175890, 175898:175905]
    for shot in shots:
        #try:
        eqm.Close( )
        eqm.Open(shot,diag='EFIT01' )
        bes = loader_BES(shot,exp='DIII-D',eqm=eqm,rho_lbl='rho_tor',MDSconn=MDSconn)
        rho_tg,theta_tg,R_tg,Z_tg = bes.get_rho('BESFU',bes.get_names('BESFU'),3 )

        print(( '%d %.2f %.2f'%(shot, min(rho_tg), max(rho_tg))))
        #except Exception as e:
            #print e
            #continue
    

     
if __name__ == "__main__":
    main()
    
    
    
    
    
