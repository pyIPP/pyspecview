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
    #print TDI 
    data = [MDSconn.get(tdi).data() for tdi in TDI]
   
    return data



class loader_DBS(loader):
    radial_profile=True

    def __init__(self,*args, **kargs):
        
        super(loader_DBS,self).__init__(*args, **kargs)

        self.names =  list(range(1,9))
        self.groups =  ['D','DV','CPS']


        self.tvec =  None
        self.catch = {}
  
    def get_names(self,group):
        return self.names
    
    def get_signal(self,group, names,calib=False,tmin=None,tmax=None):
        
        
        if tmin is None:    tmin = self.tmin
        if tmax is None:    tmax = self.tmax
        
        names = atleast_1d(names)
        ch2load = []
        for name in names:
            if not group+str(name) in self.catch:
                ch2load.append(group+str(name))
                
             
        TDIcall ='_x=PTDATA2("%s%s",%d,0)'    
        

        if size(ch2load) > 0:
            print( '\n fast parallel fetch...')
            numTasks = min(8,len(ch2load)*2)
                            
            t = T()
            server = self.MDSconn.hostspec
            TDI  = [TDIcall%(n,'a',self.shot) for n in ch2load] #sin component
            TDI += [TDIcall%(n,'b',self.shot) for n in ch2load] #cos component
    
            TDI = array_split(TDI, numTasks)
            ch2load = array_split(ch2load, numTasks)

            if self.tvec is None:  
                TDI[-1] = append(TDI[-1],'dim_of(_x)' )
            
            args = [(server, tdi) for tdi in TDI]
            
      
            pool = Pool()
            out = pool.map(mds_load,args)
            pool.close()
            pool.join()
            print(( 'in %.1fs'%(T()-t)))
    
            if self.tvec is None:  
                self.tvec = out[-1][-1]/1000#s
                
        
            for chnls, data1,data2 in zip(ch2load, out[:len(ch2load)//2],out[len(ch2load)//2:]):
                for ch,siga,sigb in zip(chnls, data1,data2):
                    siga-= int(siga.mean()) #remove offset
                    sigb-= int(sigb.mean())
                    #convert to complex signal 
                    self.catch[ch] = siga.astype('csingle')
                    self.catch[ch]*= 1j
                    self.catch[ch]+= sigb.astype('csingle')
         
        imin,imax = self.tvec.searchsorted([tmin,tmax])
        ind = slice(imin,imax+1)

        
        if len(names) == 1:
            return self.tvec[ind], self.catch[group+str(names[0])][ind]

        
        return [[self.tvec[ind], self.catch[group+str(n)][ind]] for n in names]

        

    
    def get_rho(self,group,names,time,dR=0,dZ=0):
        ind = asarray(names)-1
        return  super(loader_DBS,self).get_rho(time,self.R[ind],self.Z[ind],dR=dR, dZ=dZ)
    


    def signal_info(self,group,name, time):
        
        ###rho_tg = self.get_rho(group,[name,],time)[0]
        
        ###R,Z = self.R[name-1],self.Z[name-1]
        return 'sig:%s %.2d '%(group,name) 
    
 
 
def main():
    

    mds_server = "localhost"
    mds_server = "atlas.gat.com"

    import MDSplus as mds
    MDSconn = mds.Connection(mds_server )
    from .map_equ import equ_map
    eqm = equ_map(MDSconn,debug=False)
    
    shots = r_[175849:175869, 175882:175890, 175898:175905]
    for shot in shots:
        #try:
        #eqm.Close( )
        #eqm.Open(shot,diag='EFIT01' )
        bes = loader_DBS(175862,exp='DIII-D',rho_lbl='rho_tor',MDSconn=MDSconn)
        #rho_tg,theta_tg,R_tg,Z_tg = bes.get_rho('DBSFU',bes.get_names('DBSFU'),3 )
        A=bes.get_signal('DBS',1)
        #B=bes.get_signal('DBSb',1)
        
        #a = A[1]
        #b = B[1]
        #a-= int(mean(a))
        #b-= int(mean(b))
        #c=hypot(a,b)
        
        #plot(a[::10])
        #plot(b[::10])
        #plot(c[::10])
        #show()
        
        

 
        
        print(( '%d %.2f %.2f'%(shot, min(rho_tg), max(rho_tg))))
        #except Exception as e:
            #print e
            #continue
    

     
if __name__ == "__main__":
    main()
    
    
    
    
    
