from loaders_DIIID.loader import *
#from loader import * 
import os
from multiprocessing import  Pool

import MDSplus as mds

def check(shot):
    #fastest check if the shotfile exist
    #BUG 
    return True





class loader_VB(loader):
    radial_profile=True

    def __init__(self,*args, **kargs):
        
        super(loader_VB,self).__init__(*args, **kargs)
        nchans = 16

        self.names = {'VB': list(range(1,nchans+1))}
        self.groups = ('VB',)


        self.tvec = None
        self.catch = {}


        path = '\\TOP.VB.CHORD%.2d:'
        tdi = []
        nodes = ['PHI_START','PHI_END','R_START','R_END','Z_START','Z_END']
        for ch in range(1,nchans+1):
            for n in nodes:
                tdi.append(path%ch+n)
            
        tdi = '['+ ','.join(tdi)+']'
        
 
    
        self.MDSconn.openTree('SPECTROSCOPY', self.shot)            

        out = self.MDSconn.get(tdi).data().reshape(-1,len(nodes)).T.astype('single')
        
        self.phi_start,self.phi_end = out[:2] #deg
        self.R_start,self.R_end,self.z_start,self.z_end = out[2:]/1000 #m

 
    
    def get_names(self,group):
        return self.names[group]
    
    def get_signal(self,group, names,calib=False,tmin=None,tmax=None):       
        
        if tmin is None:    tmin = self.tmin
        if tmax is None:    tmax = self.tmax
        

        if size(names) >1:
            return [self.get_signal(group, name,calib=calib,tmin=tmin,tmax=tmax) for name in names]
        name = names

        if name in self.catch:
            sig = self.catch[name]
        else:
            
            tree = 'SPECTROSCOPY'
            self.MDSconn.openTree(tree,self.shot)
            TDIcall = "_x=\\%s::%s%.2d"%(tree,group,name )
  
            sig = np.single(self.MDSconn.get(TDIcall).data())
            if self.tvec is None:
                self.tvec = self.MDSconn.get('dim_of(_x)').data()/1e3
                
            self.MDSconn.closeTree(tree,self.shot)
            self.catch[name] = sig

        imin,imax = self.tvec.searchsorted([tmin,tmax])
        ind = slice(imin,imax+1) 
        
        return self.tvec[ind], sig[ind]

    
    def get_rho(self,group,names,time,dR=0,dZ=0):
 
        R_start = array([self.R_start[int(n)-1] for n in names])
        z_start = array([self.z_start[int(n)-1] for n in names])
        R_end = array([self.R_end[int(n)-1] for n in names])
        z_end = array([self.z_end[int(n)-1] for n in names])
        Phi_start = array([self.phi_start[int(n)-1] for n in names])
        Phi_end = array([self.phi_end[int(n)-1] for n in names])

        
        
        rho_tg,theta_tg,R,Z = super(loader_VB,self).get_rho(time,R_start,
                                    z_start,Phi_start,R_end,z_end,Phi_end,dR=dR, dZ=dZ, tangential=True)
     
      
                
        return rho_tg, theta_tg,R,Z
    
    
    def signal_info(self,group,name,time):

        rho_tg = self.get_rho(group,[ name,],time)[0]

        info = 'ch:'+str(name)+'  rho_tg: %.2f'%(rho_tg)

        return info
    

    def signal_info(self,group,name,time):
        
        return ''
    
 
 
    
 

from matplotlib.pylab import *
def main():
    
    mds_server = "localhost"
    #mds_server = "atlas.gat.com"

    import MDSplus as mds
    MDSconn = mds.Connection(mds_server )
    
    from map_equ import equ_map
    eqm = equ_map(MDSconn,debug=False)
    eqm.Open(176770,diag='EFIT01' )
    eqm._read_pfm()
    eqm._read_scalars()
    eqm._read_profiles()

    vb = loader_VB(176770,exp='DIII-D',eqm=eqm,rho_lbl='rho_pol',MDSconn=MDSconn)
 
  		
    rho_tg, theta_tg,R,Z = vb.get_rho('', arange(1,17),3)
    
    print(rho_tg)
    #import IPython
    #IPython.embed()
    
    
if __name__ == "__main__":
    main()
    
