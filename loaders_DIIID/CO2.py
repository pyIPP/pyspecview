from loaders_DIIID.loader import * 
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
        
    def get_signal(self,group, names,calib=False,tmin=None,tmax=None):
        
        #TODO raw data ? fast?  
        if tmin is None:    tmin = self.tmin
        if tmax is None:    tmax = self.tmax
        
        if size(names) > 1:
            data = [self.get_signal(group, n, calib=calib, tmin=tmin,tmax=tmax) for n in names]
            return data
        
        name = names
        
        node = 'PL1' if group == "PHASE" else 'DEN'
        TDI = '\\ELECTRONS::TOP.BCI.DPD.%s.%s:%s_UF_'%(name, group,node)+'%d'
        
            
        if not hasattr(self,'tvec'):
            index = list(range(self.n_chunks))
            self.tvec = mds_par_load(self.MDSconn.hostspec,
                    self.tree,self.shot,'dim_of('+TDI+')',index)
        
            for i in range(self.n_chunks):
                self.tvec[i] /= 1e3 #s

        #print self.tvec
        #print tmin, tmax
        #print [t[-1] > tmin for t in self.tvec]
        #print [t[ 0] < tmax for t in self.tvec]
        indmin = where([t[-1] > tmin for t in self.tvec])[0][0]
        indmax = where([t[ 0] < tmax for t in self.tvec])[0][-1]+1
        #print indmin,indmax
        index = list(range(indmin,indmax))
        index_toload = []
        for i in index:
            if self.cache[group][name][i] is None:
                index_toload.append(i)
        #print index_toload
        if len(index_toload):
            data = mds_par_load(self.MDSconn.hostspec,self.tree,self.shot,TDI,index_toload)
            for i, d in zip(index_toload, data):
                d*=1e6 #m^-2
                self.cache[group][name][i] = d 
        
        
        tvec = hstack(self.tvec[indmin:indmax])
        sig  = hstack(self.cache[group][name][indmin:indmax])
        

        sig = self.remove_elms(tvec, sig)
        print('elms removed')
        
        
        ##tvec=np.linspace(0,5,100000000*5)
        ##sig = np.random.randint(0,2**16,100000000*5)
        
        ##return tvec,sig
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

    
    mds_server = "localhost"
    mds_server = "atlas.gat.com"

    import MDSplus as mds
    MDSconn = mds.Connection(mds_server )
    from .map_equ import equ_map
    eqm = equ_map(MDSconn,debug=False)
    eqm.Open(153291,diag='EFIT01' )
    sxr = loader_CO2(153291,exp='DIII-D',eqm=eqm,rho_lbl='rho_pol',MDSconn=MDSconn)
    data1 = sxr.get_signal( 'DEN',('V1',),tmin=-infty, tmax=infty,calib=True)


    #import IPython
    #IPython.embed()


    g = sxr.groups[0]
    n = sxr.get_names(g)
    data1 = sxr.get_signal( '45R1',list(range(1,13)),tmin=-infty, tmax=infty,calib=True)
    data2 = sxr.get_signal('165R1',list(range(1,13)),tmin=-infty, tmax=infty,calib=True)
    data3 = sxr.get_signal('195R1',list(range(1,13)),tmin=-infty, tmax=infty,calib=True)


    data = array([d for t,d in data1]+[d for t,d in data2]+[d for t,d in data3])
    tvec = data1[0][0]
    
    
    #\\ELECTRONS::TOP.BCI.DPD.V2.DEN:DEN_UF_8
    #\\ELECTRONS::TOP.BCI.DPD.V1.DEN:PL1_UF_8

    offset = tvec < .1
    data_ = data[:,offset]-data[:,offset].mean(1)[:,None]
    u1,s1,v1 = linalg.svd(data_[:12], full_matrices=False)
    u2,s2,v2 = linalg.svd(data_[12:24], full_matrices=False)
    u3,s3,v3 = linalg.svd(data_[24:], full_matrices=False)


    from scipy import signal
    fnq = len(tvec)/(tvec[-1]-tvec[0])/2
    fmax = 50
    b, a = signal.butter(4, fmax/fnq, 'low')
    noise1 = inner(u1[:,0], data[:12].T)
    noise1 -= signal.filtfilt(b,a,noise1)
    noise2 = inner(u2[:,0], data[12:24].T)
    noise2 -= signal.filtfilt(b,a,noise2)
    noise3 = inner(u3[:,0], data[24:].T)
    noise3 -= signal.filtfilt(b,a,noise3)
    
    
    
    filtered_data = copy(data)
    filtered_data[:12]   -= outer( u1[:,0],noise1)
    filtered_data[12:24] -= outer( u2[:,0],noise2)
    filtered_data[24:]   -= outer( u3[:,0],noise3)

    fmax = 3000
    b, a = signal.butter(6, fmax/fnq, 'low')
    filtered_data = signal.filtfilt(b,a,filtered_data,axis=1)
    
    #fmax = 2000
    #b, a = signal.butter(6, fmax/fnq, 'low')
    #filtered_data = signal.filtfilt(b,a,filtered_data,axis=1)
    i1 = tvec<.1
    i2 = tvec> tvec[-1]-.5
    b1,b2 = filtered_data[:,i1].mean(1), filtered_data[:,i2].mean(1)
    a1,a2 = tvec[i1].mean(), tvec[i2].mean()

    filtered_data -= ((b2-b1)/(a2-a1)*(tvec[:,None]-a1)+b1).T
    offset_err = sqrt(mean(filtered_data[:,tvec<.3]**2,1))
    error =  offset_err[:,None]+filtered_data*0.05
    cov_mat = corrcoef(filtered_data[:,tvec<.3])

    
    f,ax=subplots(2,1,sharex=True, sharey=True)
    ax[0].plot(tvec, data[:12].T)    
    #ax[1].plot(tvec, filtered_data.T)
    ax[1].plot(tvec, filtered_data[:12].T)

    ax[1].set_xlabel('time [s]')
    ax[1].set_ylabel('filtered SXR')
    ax[0].set_ylabel('raw SXR')

    show()
    
    
    [errorbar(list(range(36)), d,e) for d,e in zip(filtered_data[:,::10000].T,error[:,::10000].T)];show()

    

    imshow(filtered_data,interpolation='nearest',aspect='auto',vmax=-0.1,vmin=0.1);colorbar();show()

    
    import IPython
    IPython.embed()
    
    exit()
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    #plot(range(32),data[:32,35000:35100])
    #plot(range(32,64),data[32:,35000:35100])
    #show()
    
    

    data[2] = 0
    offset = tvec < .1
    data_ = data[:,offset]-data[:,offset].mean(1)[:,None]
    u1,s1,v1 = linalg.svd(data_[:32], full_matrices=False)
    u2,s2,v2 = linalg.svd(data_[32:], full_matrices=False)

    #plot(range(28),u1[:,0] )
    #plot(range(28,28*2), u2[:,0])
    from scipy import signal
    fnq = len(tvec)/(tvec[-1]-tvec[0])/2
    fmax = 50
    b, a = signal.butter(4, fmax/fnq, 'low')
    noise1 = inner(u1[:,0], data[:32].T)
    noise1 -= signal.filtfilt(b,a,noise1)
    noise2 = inner(u2[:,0], data[32:].T)
    noise2 -= signal.filtfilt(b,a,noise2)
    filtered_data = copy(data)
    filtered_data[:32]  -= outer( u1[:,0],noise1)
    filtered_data[32:]  -= outer( u2[:,0],noise2)
    
    
    fmax = 2000
    b, a = signal.butter(6, fmax/fnq, 'low')
    filtered_data = signal.filtfilt(b,a,filtered_data,axis=1)
    
    
    i1 = tvec<.1
    i2 = tvec> tvec[-1]-.5
    b1,b2 = filtered_data[:,i1].mean(1), filtered_data[:,i2].mean(1)
    a1,a2 = tvec[i1].mean(), tvec[i2].mean()

    filtered_data -= ((b2-b1)/(a2-a1)*(tvec[:,None]-a1)+b1).T
    offset_err = sqrt(mean(filtered_data[:,tvec<.3]**2,1))
    error =  offset_err[:,None]+filtered_data*0.05
    cov_mat = corrcoef(filtered_data[:,tvec<.3])

    #plot(v1[0])
    #plot(v2[0])

    #print data.shape, tvec.shape
    import IPython
    IPython.embed()
    #plot(abs(fft.rfft(data_[0])))

    #NOTE odecist pozadi pred aplikaci IIR  filtru!
    
    f,ax=subplots(2,1,sharex=True, sharey=True)
    ax[0].plot(tvec, data.T)    
    #ax[1].plot(tvec, filtered_data.T)
    ax[1].plot(tvec, filtered_data2.T)

    ax[1].set_xlabel('time [s]')
    ax[1].set_ylabel('filtered SXR')
    ax[0].set_ylabel('raw SXR')

    show()
    
    
    
    offset = filtered_data[:,80000:].mean(1)[:,None]

    contourf(data,20)
    
    #filtered_data-= filtered_data[:,80000:].mean(1)[:,None]
    
    imshow(filtered_data,interpolation='nearest',aspect='auto',vmin=-1000, vmax=1000, extent=(tvec[0], tvec[-1], 0,1));colorbar();show()
    
    
    [errorbar(list(range(64)), d,e) for d,e in zip(filtered_data[:,::10000].T,error[:,::10000].T)];show()

    imshow(filtered_data, interpolation='nearest',aspect='auto');colorbar();show()
    plot((filtered_data-offset)[:,40000]);show()

    import matplotlib.pylab as plt
    plt.plot(tvec, sig)
    plt.show()
    
if __name__ == "__main__":
    main()
    
