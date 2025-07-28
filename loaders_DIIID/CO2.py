try:
    from loaders_DIIID.loader import * 
except:
    from loader import * 
 
import os
from multiprocessing import  Pool, cpu_count
import MDSplus as mds

#T

def check(shot):
    #fastest check if the shotfile exist
    #BUG 
    status = True

    return status

verbose = False


def mds_load(tmp):
    (mds_server,tree, shot,  TDI) = tmp
    #print TDI
    MDSconn = mds.Connection(mds_server )
    MDSconn.openTree(tree, shot)
    data = MDSconn.get(TDI).data()
    MDSconn.closeTree(tree, shot)

    #print mds_server,tree, shot,  TDI

    return data


def mds_par_load(mds_server,tree, shot,  TDI, indexes):

    #load a junks of a single vector

    args = [(mds_server,tree, shot, TDI%i) for i in indexes]
  
    nconn = len(args)  #brutal force
    pool = Pool()
    out = pool.map(mds_load,args)
    pool.close()
    pool.join()
    
    return out
    



class loader_CO2(loader):
    

    radial_profile=True
    units = ''


    def __init__(self,*args, **kargs):
        
        super(loader_CO2,self).__init__(*args, **kargs)

        from time import time 
        
        t = time()
        
        self.tree = 'ELECTRONS'
        #self.MDSconn.openTree(self.tree,self.shot)
        
        #geometry 
        #Set the z top and bottom for a purely vertical LOS
        Z_top = 1.24
        Z_bottom = -1.375

        #Set the R left and right for a purely horizontal LOS
        R_lfs = 2.36
        R_hfs = 1.01

        #phi locations of the vertical and radial chords
        phi_vertical = 240
        phi_radial = 225

        #Format: (name, stat_error_threshold, LOS_pt1, LOS_pt2)
        self.geom = {'V1':([1.48,Z_top,phi_vertical], [1.48,Z_bottom,phi_vertical]),
                     'V2':([1.94,Z_top,phi_vertical], [1.94,Z_bottom,phi_vertical]),
                     'V3':([2.10,Z_top,phi_vertical], [2.10,Z_bottom,phi_vertical]),
                     'R0':([R_lfs,0,phi_radial], [R_hfs,0,phi_radial])}
        
        
        self.names = 'V1', 'V2', 'V3', 'R0'
        self.groups = 'DEN', 'PHASE'
        
        self.n_chunks = 8 #BUG is it the same in all shots? 
        
        self.cache = {'DEN':  {n:empty(self.n_chunks,dtype=object) for n in self.names},
                      'PHASE':{n:empty(self.n_chunks,dtype=object) for n in self.names}}
        
         
        
    def get_names(self,group):
        return self.names
        
    def get_signal(self,group, names, calib=False, tmin=None, tmax=None):
        
        #TODO raw data ? fast?  
        if tmin is None:    tmin = self.tmin
        if tmax is None:    tmax = self.tmax
        
        if not isinstance(names, str):
            data = [self.get_signal(group, n, calib=calib, tmin=tmin,tmax=tmax) for n in names]
            return data
        
        name = names
        
        node = 'PL1' if group == "PHASE" else 'DEN'
        
        if self.shot < 198000: #before 2024 campaign
            TDI = '\\ELECTRONS::TOP.BCI.DPD.%s.%s:%s_UF_'%(name, group,node)+'%d'
        
            
            if not hasattr(self,'tvec'):
                index = list(range(self.n_chunks))
                print(TDI)
                self.tvec = mds_par_load(self.MDSconn.hostspec,
                        self.tree,self.shot,'dim_of('+TDI+')',index)
            
                for i in range(self.n_chunks):
                    self.tvec[i] /= 1e3 #s

      
            indmin = where([t[-1] > tmin for t in self.tvec])[0][0]
            indmax = where([t[ 0] < tmax for t in self.tvec])[0][-1]+1
          
            index = list(range(indmin,indmax))
            index_toload = []
            for i in index:
                if self.cache[group][name][i] is None:
                    index_toload.append(i)
            
            if len(index_toload):
                data = mds_par_load(self.MDSconn.hostspec,self.tree,self.shot,TDI,index_toload)
                for i, d in zip(index_toload, data):
                    d*=1e6/1e19 #10**19*m^-2
                    self.cache[group][name][i] = d 
            
            
            tvec = hstack(self.tvec[indmin:indmax])
            sig  = hstack(self.cache[group][name][indmin:indmax])
            
        else:
        
            MDSconn = mds.Connection(self.MDSconn.hostspec)
            MDSconn.openTree(self.tree, self.shot)
            TDI = '_x=\\ELECTRONS::TOP.BCI.DPD.%s:%sUF'%(name,node)
            if len(self.cache[group][name]) == self.n_chunks:
                self.cache[group][name] = MDSconn.get(TDI).data().astype('single')
                
            if not hasattr(self,'tvec'):
                self.tvec =  MDSconn.get('dim_of(_x)').data()
                self.tvec /= 1e3 #s
                
             
            tvec = self.tvec
            sig = self.cache[group][name]


    
        

        sig = self.remove_elms(tvec, sig)        
        
 
        imin,imax = tvec.searchsorted([tmin,tmax])
        ind = slice(imin,imax+1)
        
        
        
        return tvec[ind], sig[ind]
            



            
    def get_rho(self,group,names,time,dR=0,dZ=0):
        
        R_start= zeros(size(names))
        z_start= zeros(size(names))
        Phi1= zeros(size(names))
        R_end= zeros(size(names))
        z_end= zeros(size(names))
        Phi2 = zeros(size(names))
        
        for i,n in enumerate(names):
            (R_start[i],z_start[i],Phi1[i]),(R_end[i],z_end[i],Phi2[i]) = self.geom[n]
        
        rho_tg,theta_tg,R,Z = super(loader_CO2,self).get_rho(time,R_start,
                                    z_start,Phi1,R_end,z_end,Phi2,dR=dR, dZ=dZ)

        return abs(rho_tg), theta_tg,R,Z
    

    
            
    def signal_info(self,group,name,time):
        
        rho_tg = self.get_rho(group,[name,],time)[0]
        
        #phi = self.Phi[name]
        
        info = str(name)+':'+group+', '+self.rho_lbl+': %.2f'%rho_tg
        return info
    
  
 
    

    #def signal_info(self,group,name,time):
        #info = group+': '+str(name)
        #return info
    
    
from matplotlib.pylab import *
def main():

    shot = 199035
    mds_server = "localhost"
    mds_server = "atlas.gat.com"

    import MDSplus as mds
    MDSconn = mds.Connection(mds_server )
    from map_equ import equ_map
    eqm = equ_map(MDSconn,debug=False)
    eqm.Open(shot,diag='EFIT01' )
    sxr = loader_CO2(shot,exp='DIII-D',eqm=eqm,rho_lbl='rho_pol',MDSconn=MDSconn)
    data1 = sxr.get_signal( 'DEN',('V1',),tmin=-infty, tmax=infty,calib=True)


    import IPython
    IPython.embed()


    
    
    
    
if __name__ == "__main__":
    main()
    
