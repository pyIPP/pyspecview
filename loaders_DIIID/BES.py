try:
    from loaders_DIIID.loader import * 
except:
    from loader import * 
 


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
    tvec_fast = None

    def __init__(self,*args, **kargs):
        
        super(loader_BES,self).__init__(*args, **kargs)

        self.names = {'BESFU': list(range(1,65)) }
        self.groups =  ['BESFU']

 
        self.catch = { }
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
        ch2load = [n for n in names if not n in self.catch] 
        
        #NOTE this PTsignal cannot be splitted in segments
        TDIcall ='PTDATA2("%s%.2d",%d)' 

 
        if self.tvec_fast is None:
            header = self.MDSconn.get(f'PTHEAD2("BESFU{ch2load[0]:02d}",{self.shot}); __real64').data()
            tbeg, tend, dt =  header[2:5]
            self.tvec_fast = np.arange(tbeg, tend, dt)
            header = self.MDSconn.get(f'PTHEAD2("BESSU{ch2load[0]:02d}",{self.shot}); __real64').data()
            tbeg, tend, dt =  header[2:5]
            self.tvec_slow = np.arange(tbeg, tend, dt)
           
         
        
        if  size(ch2load) == 1:
        
            fTDI = TDIcall%('BESFU',ch2load[0],self.shot)
            fsig = self.MDSconn.get(fTDI).data()
            if len(fsig) < 2:
                raise Exception(group+' channel %d do not exist'%name)
            assert len(fsig) == len(self.tvec_fast)

            sTDI = TDIcall%('BESSU',ch2load[0],self.shot)
            ssig = self.MDSconn.get(sTDI).data()
 
            #relative perturbations
            #self.catch[ch2load[0]] = np.single(fsig / np.interp(self.tvec_fast, self.tvec_slow, ssig))
            ssig = maximum(ssig, ssig.max()/1e3)
            islow = np.single(np.interp(self.tvec_fast, self.tvec_slow, ssig))
            self.catch[ch2load[0]] = (fsig)/islow
                    
        elif size(ch2load) > 1:
            print( '\n fast parallel fetch...' )
            numTasks = 8
                            



            t = T()
            server = self.MDSconn.hostspec
            TDI = [TDIcall%('BESFU',n,self.shot) for n in ch2load]
                        
            TDI = array_split(TDI, numTasks)
 
            pool = Pool()
            fast_out = pool.map(mds_load,[(server, tdi) for tdi in TDI])
     
            TDI = [TDIcall%('BESSU',n,self.shot) for n in ch2load]    
            TDI = array_split(TDI, numTasks)
       
            slow_out = pool.map(mds_load,[(server, tdi) for tdi in TDI])
            pool.close()
            pool.join()
            ch2load = array_split(ch2load, numTasks)

            print(( 'in %.1fs'%(T()-t)))
     
        
            for chnls, fdata, sdata in zip(ch2load, fast_out, slow_out):
                for ch,fsig, ssig in zip(chnls, fdata, sdata):
                    ssig = maximum(ssig, ssig.max()/1e3)
                    islow = np.single(np.interp(self.tvec_fast, self.tvec_slow, ssig))
                    self.catch[ch] = (fsig)/islow
         

        imin,imax = self.tvec_fast.searchsorted([tmin,tmax])
        ind = slice(imin,imax+1)
        
        
        if len(names) == 1:
            return self.tvec_fast[ind], self.catch[names[0]][ind]

        
        return [[self.tvec_fast[ind], self.catch[n][ind]] for n in names]

        

    
    def get_rho(self,group,names,time,dR=0,dZ=0):
        ind = asarray(names)-1
        return  super(loader_BES,self).get_rho(time,self.R[ind],self.Z[ind],dR=dR, dZ=dZ)
    


    def signal_info(self,group,name, time):
        
        rho_tg = self.get_rho(group,[name,],time)[0]
        
        R,Z = self.R[name-1],self.Z[name-1]
        return 'sig:%s %.2d   %s:%.2f R:%.3fm z:%.3fm'%(group,name, self.rho_lbl,rho_tg,R,Z)
    
 
 
def main():
    

    mds_server = "localhost"
    mds_server = "atlas.gat.com"

    import MDSplus as mds
    MDSconn = mds.Connection(mds_server )
    from map_equ import equ_map
    shot = 175860
    eqm = equ_map(MDSconn,debug=False)

    eqm.Open(shot,diag='EFIT01' )
    bes = loader_BES(shot,exp='DIII-D',eqm=eqm,rho_lbl='rho_tor',MDSconn=MDSconn)
    bes.get_signal('BESFU',bes.get_names('BESFU'),3, 3.2 )
    exit()
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
    
    
    
    
    
