from loaders_DIIID.loader import * 
import MDSplus as mds

def check(shot):
    return True





class loader_RIP(loader):

    def __init__(self,*args, **kargs):
        
        super(loader_RIP,self).__init__(*args, **kargs)

        self.names = {'n_e':list(range(1,4)),'B_r':list(range(1,4))}
        self.groups =  list(self.names.keys())
        self.z = 0,0.135, -0.135 #m vertical position

        self.tvec = None
        self.catch = {}


    
    def get_names(self,group):
        return self.names[group]
    
    def get_signal(self,group, names,calib=False,tmin=None,tmax=None):
         
        
        if tmin is None:    tmin = self.tmin
        if tmax is None:    tmax = self.tmax
        

        if size(names) >1:
            return [self.get_signal(group, name,calib=calib,tmin=tmin,tmax=tmax) for name in names]
        name = names
        
        
        if group == 'B_r':
            tagname = 'rpich%dphi'%(name*2-1)
            
        if group == 'n_e':
            tagname = 'rpich%dphi'%(name*2  )
            
        tree = 'RPI'
        if tagname in self.catch:
            sig = self.catch[tagname]
        else:
            
            self.MDSconn.openTree(tree,self.shot)

            TDIcall ='_x= \\%s::%s'%(tree,tagname)
            sig = self.MDSconn.get(TDIcall).data()
 
            self.catch[tagname] = sig

            if self.tvec is None:
                self.tvec = self.MDSconn.get('dim_of(_x)').data()
                self.tvec /= 1e3 #s
            self.MDSconn.closeTree(tree,self.shot)

                    
        imin,imax = self.tvec.searchsorted([tmin,tmax])
        ind = slice(imin,imax+1) 
        
        return self.tvec[ind], sig[ind]

                
    def get_rho(self,group,names,time,dR=0,dZ=0):

        R_start = array(2.6,ndmin=1)
        R_end = array(1,ndmin=1)
        name = names[0]
        z = array(self.z[name-1],ndmin=1)

        rho_tg,theta_tg,R,Z = super(loader_RIP,self).get_rho(time,R_start,
                                    z,0,R_end,z,0,dR=dR, dZ=dZ)
        
  
        return rho_tg, theta_tg,R,Z
    
    
    
    def signal_info(self,group,name,time):

        rho_tg = self.get_rho(group,[ name,],time)[0]
 
        info = 'ch:'+str(name)+' z:%.1fcm  rho_tg: %.2f'%(self.z[name-1]*100,rho_tg)

        return info
    
 
from matplotlib.pylab import *
def main():
         
    mds_server = "localhost"
    #mds_server = "atlas.gat.com"

    import MDSplus as mds
    MDSconn = mds.Connection(mds_server )
    
    from .map_equ import equ_map
    eqm = equ_map(MDSconn,debug=False)
    eqm.Open(175847,diag='EFIT01' )
    rip = loader_RIP(175847,exp='DIII-D',eqm=eqm,rho_lbl='rho_tor',MDSconn=MDSconn)
    
    f,ax = subplots(2)
    t,n = rip.get_signal('n_e',1,tmin=-1)
    n-= n[t<0].mean()
    ax[0].plot(t*1000,-n,'r')
    ax[0].axvline(1800,c='k')
    ax[0].set_title('RIP on-axis LOS')
    ax[0].set_xlim(1700,2000)
    ax[0].set_ylim(1950,2150)
    ax[0].set_ylabel('ne (a.u.)')

 
    rip = loader_RIP(175894,exp='DIII-D',eqm=eqm,rho_lbl='rho_tor',MDSconn=MDSconn)
 
    t,n = rip.get_signal('n_e',1,tmin=-1)
    n-= n[t<0].mean()

    ax[1].plot(t*1000,-n)
    ax[1].axvline(2000,c='k')
    
    
    ax[1].set_xlim(1900,2200)
    ax[1].set_xlabel('Time (ms)')
    ax[1].set_ylabel('ne (a.u.)')
    ax[1].set_ylim(1600,1800)

    show()

if __name__ == "__main__":
    main()
    
