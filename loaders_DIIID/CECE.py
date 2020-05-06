from loaders_DIIID.loader import * 
import os
from multiprocessing import  Pool

#from time import time 

import MDSplus as mds

def check(shot):
    #fastest check if the shotfile exist
    #BUG 
    return True





class loader_CECE_(loader):

    def __init__(self,*args, **kargs):
        
        super(loader_CECE_,self).__init__(*args, **kargs)

        self.names = {'A': ['cece%d'%i for i in range(1,9)],'B':['cece%d'%i for i in range(1,9)]}
        self.groups = list(self.names.keys())


        self.tvec = None
        self.catch = {}


    
    def get_names(self,group):
        return self.names[group]
    
    def get_signal(self,group, names,calib=False,tmin=None,tmax=None):
        #BUG load each chunk separatelly? 
        
        
        if tmin is None:    tmin = self.tmin
        if tmax is None:    tmax = self.tmax
        

        if size(names) >1:
            return [self.get_signal(group, name,calib=calib,tmin=tmin,tmax=tmax) for name in names]
        name = names+group

        if name in self.catch:
            sig = self.catch[name]
        else:
            TDIcall ='_x=PTDATA2("%s",%d)'%(name,self.shot)
            #print TDIcall
            sig = self.MDSconn.get(TDIcall).data()
            if len(sig) < 2:
                raise Exception('CECE channel %d do not exist'%name)
            self.catch[name] = sig

        if self.tvec is None:
                self.tvec = self.MDSconn.get('dim_of(_x)').data()
                self.tvec /= 1e3 #s
                
        imin,imax = self.tvec.searchsorted([tmin,tmax])
        ind = slice(imin,imax+1) 
        
        return self.tvec[ind], sig[ind]

    


    def signal_info(self,group,name,time):
        return ''
    
 



def main():
    
    import matplotlib.pylab as plt    
    mds_server = "localhost"
    #mds_server = "atlas.gat.com"

    import MDSplus as mds
    MDSconn = mds.Connection(mds_server )
    from .map_equ import equ_map
    eqm = equ_map(MDSconn,debug=False)
    
    shots = r_[175849:175869, 175882:175890, 175898:175905]
    shots = r_[175860, 175861]
    shots = r_[175859, 175861]

    for shot in shots:
        #try:
        #print shot
        eqm.Close( )
        eqm.Open(shot,diag='EFIT01' )
        cece = loader_CECE_(shot,exp='DIII-D',rho_lbl='rho_tor',MDSconn=MDSconn)
        #rho_tg,theta_tg,R_tg,Z_tg = bes.get_rho('DBSFU',bes.get_names('DBSFU'),3 )
        names = cece.get_names('A')
        tvec,sigA = cece.get_signal('A',names[0], tmin = 3.5, tmax = 4)
        tvec,sigB = cece.get_signal('B',names[0], tmin = 3.5, tmax = 4)

        tvec2,sigA2 = cece.get_signal('A',names[0], tmin = 4.1, tmax = 4.6)
        _,sigB2 = cece.get_signal('B',names[0], tmin = 4.1, tmax = 4.6)

        from scipy.signal import coherence
        
        f,C = coherence(sigA, sigB, fs=1/mean(diff(tvec)), window='hann', nperseg=1000, noverlap=None, nfft=None,detrend='constant')

        plt.plot(f/1000,C)
        plt.xlim(0,1000) 
        
         
        f2,C2 = coherence(sigA2, sigB2, fs=1/mean(diff(tvec2)), window='hann', nperseg=1000, noverlap=None, nfft=None,detrend='constant')

        plt.plot(f2/1000,C2)
        plt.xlim(0,1000);plt.show()
        
        
        plt.plot( sigA2),  plt.plot( sigB2);plt.show()
        
        
        
        
        
        #tvec,off = cece.get_signal('A',names[0], tmin = -5,tmax = .1)
        #tvec,sig2= cece.get_signal('A',names[0], tmin = 3, tmax = 3.5)
        #tvec,sig3= cece.get_signal('A',names[0], tmin = 4.5, tmax = 4.5)
        
        #print shot, off.mean(), sig.mean(),  sig2.mean(), sig3.mean()
        #B=bes.get_signal('DBSb',1)
        import IPython
        
        IPython.embed()   
        tvec,sig=cece.get_signal('A',names[-1], tmin = -2)
        sig = sig[:len(sig)/1000*1000].reshape(-1,1000).mean(1)
        tvec=tvec[:len(tvec)/1000*1000].reshape(-1,1000).mean(1)
        off = sig[tvec<0].mean()
        plt.plot(tvec, off-sig)
        plt.xlim(0,5)
        plt.ylim(0,5)
        plt.savefig('%d.png'%shot)
        plt.clf()
        
        
        
        
        #a = A[1]
        #b = B[1]
        #a-= int(mean(a))
        #b-= int(mean(b))
        #c=hypot(a,b)
        
        #plot(a[::10])
        #plot(b[::10])
        #plot(c[::10])
        #show()
        
        

 
        
        #print '%d %.2f %.2f'%(shot, min(rho_tg), max(rho_tg))
        #except Exception as e:
            #print e
            #continue
    

     
if __name__ == "__main__":
    main()
    
    
    
    
    
