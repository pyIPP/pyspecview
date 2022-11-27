try:
    from loaders_DIIID.loader import * 
except:
    from loader import * 

import os
from multiprocessing import  Pool
import MDSplus as mds

#TODO lod only a data slice
#TODO calibration for all diagnostics 
#TODO 



def check(shot):
    #fastest check if the shotfile exist
    #BUG 
    status = True

    return status

verbose = False




from time import time as T




def mds_load(tmp):
    mds_server,  TDI = tmp
    MDSconn = mds.Connection(mds_server )
    data = []
    for tdi in TDI:
        try:
            data.append(MDSconn.get(tdi).data())
        except:
            data.append([])
    return data


def mds_load_par( MDSconn, TDI, numTasks=8):
    
    if len(TDI) == 1:
        return [MDSconn.get(TDI[0]).data()]
    
    numTasks = 8
    server = MDSconn.hostspec
    TDI = array_split(TDI, min(numTasks,len(TDI)))

    args = [(server, tdi) for tdi in TDI]
    
    out = []
    pool = Pool()
    for o in pool.map(mds_load,args):
        out.extend(o)
    pool.close()
    pool.join()
    
    return out 



from IPython import embed
class loader_ECE(loader):
    

    radial_profile=True
    units = 'eV'
    
    #can improve performace over very slow networks

    def __init__(self,*args, **kargs):
        
        super(loader_ECE,self).__init__(*args, **kargs)

   
        
        self.data_dict = {}

        self.tree = 'ece'
        self.MDSconn.openTree(self.tree,self.shot)
        
        
        # set up strings to use in mds loading
        nodes = '\\ECE::TOP.SETUP.'
        nodec = '\\ECE::TOP.CALF:'
        noded = '\\ECE::TOP.TECEF'
        subnodes_phi = nodes + 'ECEPHI'
        subnodes_nch = nodec +  'NUMCHF'
        subnodes_theta = nodes + 'ECETHETA'
        subnodes_z = nodes + 'ECEZH'
        subnodes_freq = nodes + 'FREQ'
        subnodes_lofreq = nodes + 'LOFREQ' #Also not yet used - but may be something folks will want later.
        subnodes_valid = nodec +  'VALIDF'
        subnodes_channels = nodec +  'ECPTNAMF'
        calid_cc1f = nodec +  'CC1F'
        calid_cc2f = nodec +  'CC2F'
        calid_cc3f = nodec +  'CC3F'
        calid_adosf = nodec +  'ADOSF'

        
        tag='\TECEF'
        TDIcall= tag+"01"

        self.tvec = None
            
        if self.shot < 109400:
            self.nchs = 32
        else:
            self.nchs = self.MDSconn.get(subnodes_nch).data()
        #BUG!!!1   which channels are turned off? 
        #self.nchs = 48
        self.names = arange(1,self.nchs+1)
        self.groups = ('ECE', ) 

        # get frequency settings
        self.freq = self.MDSconn.get(subnodes_freq).data()*1e9
        
        self.valid = bool_(self.MDSconn.get(subnodes_valid))
        self.channels = self.MDSconn.get( subnodes_channels  )
        
        #calibration factors
        self.cc1f = self.MDSconn.get(calid_cc1f)
        self.cc2f = self.MDSconn.get(calid_cc2f)
        self.cc3f = self.MDSconn.get(calid_cc3f)
        self.adosf= self.MDSconn.get(calid_adosf)
 
        # get position settings
        self.z = self.MDSconn.get(subnodes_z).data()
        #print('ECE PHI', self.MDSconn.get(subnodes_phi).data())
        #self.Phi  =  21.#81.60 #from diagnostic webpage #self.MDSconn.get(subnodes_phi).data()#deg
        self.Phi = 81
        self.Phi = 60        
        self.MDSconn.closeTree(self.tree, self.shot)

        if self.freq is None or len(self.freq)<2:
            raise Exception("Failed to get ECE channels")
       
        #BUG!!!1   which channels are turned off? 
        self.freq = self.freq[:self.nchs]
        self.channels = self.channels[:self.nchs]
        #raw_ch_name = self.channels[nch-1]
        try:
            self.channels = [ch.decode for ch in self.channels]
        except:
            pass
    
    def get_signal(self,group, names,calib=False,tmin=None,tmax=None):    
    
        if tmin is None:    tmin = self.tmin
        if tmax is None:    tmax = self.tmax
        
        
        nch = atleast_1d(int_(names))
        
         
        load_ch = [n for n in nch if n not in self.data_dict]


        #BUG https://diii-d.gat.com/diii-d/ECE#pointnames
        #Ask if ther are any issues with timing in the new discharges 
        
        if len(load_ch) >= 1:
            #NOTE ece data can be splitted in halve and fetch separatelly
            
            TDI = [f'PTDATA2("{self.channels[n-1]}", {self.shot}, 4)' for n in load_ch]
            tt = T()
            out = mds_load_par( self.MDSconn, TDI)
            for n, Te in zip(load_ch, out ):
                self.data_dict[n] = Te
     
        
        if self.tvec is None:
            time_header = self.MDSconn.get(f'PTHEAD2("{self.channels[0]}",{self.shot}); __real64').data()
            tbeg, tend = time_header[2], time_header[-1]
            self.tvec = np.linspace(tbeg, tend, len(Te))
        
        imin,imax = self.tvec.searchsorted([tmin,tmax])
        imax += 1
            
      
        #calibrate data locally, faster then fetching of calibrated data
        for n in load_ch:
            Te  = self.data_dict[n]
            if len(Te) in [0,1]: 
                self.data_dict[n] = np.zeros_like(self.tvec, dtype='single')
                continue
         
            C1 = float(self.cc1f[n-1])
            C2 = float(self.cc2f[n-1])
            C3 = float(self.cc3f[n-1])
            adosf = float(self.adosf[n-1])
                            
            if C2 == 0 and C3 == 0 and adosf == 0:
                Te *= C1*1e3 #eV
            else:
                Te_ = np.copy(Te)
                Te *= C1*1e3 #eV
                if C2 != 0:
                    Te /= np.sqrt(1.-C2*Te_**2)
                if C3 != 0:
                    Te += (C3*1e3)*Te_**4
                    
            if adosf != 0:
                Te += adosf*1e3
                
            self.data_dict[n] = Te
            
     
        output = [[self.tvec[imin:imax],  self.data_dict[n][imin:imax]] for n in nch]
        
        if len(nch) == 1:
            return output[0]
        else:
            return output 
 
        
    def get_names(self,group):
        return self.names
    
    
    def get_RZ_theta(self, time,names,dR=0,dZ=0):

        if hasattr(self,'tvec'):
            time = clip(time, *self.tvec[[0,-1]])
  
        B = self.eqm.rz2brzt(r_in=self.eqm.Rmesh, z_in=self.z, t_in=time)
        Btot = squeeze(linalg.norm(B,axis=0))
        
        ch_ind = atleast_1d(in1d(self.names, int_(names)))
        
        from scipy.constants import m_e, e, c

        #Accounting for relativistic mass downshift
        try:
            Te =  self.get_Te0(time,time,R=self.eqm.Rmesh,
                                Z=self.z*ones_like(self.eqm.Rmesh))[1][0]

            v=sqrt(2*Te*e/m_e)
            gamma = 1/sqrt(1-(v/c)**2)
        except:
            print( 'relativistic mass downshift could not be done')
            gamma = 1
            
        wce = e*Btot/(m_e*gamma)

        nharm = 2
        R = interp(-2*pi*self.freq[ch_ind],-wce*nharm,self.eqm.Rmesh)
        z = self.z*ones_like(R)+.04

        r0 = interp(time, self.eqm.t_eq, self.eqm.ssq['Rmag'])+dR
        z0 = interp(time, self.eqm.t_eq, self.eqm.ssq['Zmag'])+dZ

        return R,z, arctan2(z-z0, R-r0)
        
    def get_Te0(self,tmin,tmax,IDA=True,dR=0,dZ=0,R=None,Z=None):
        #load ziprofiles
        try:
            if not hasattr(self,'zipcache'):
                self.zipcache = []
                TDIcall = '\ELECTRONS::TOP.PROFILE_FITS.ZIPFIT:'

                self.MDSconn.openTree('ELECTRONS' ,self.shot)
                ZipTe = self.MDSconn.get('_x='+TDIcall+'ETEMPFIT').data()
                ZipNe = self.MDSconn.get('_x='+TDIcall+'EDENSFIT').data()
                ZipTe*= 1e3 #eV
                ZipNe*= 1e19 #m^-3
                zip_tvec = self.MDSconn.get('dim_of(_x,1)').data()
                zip_tvec/= 1e3#s
                zip_rho = self.MDSconn.get('dim_of(_x,0)').data()
                zip_rho = self.eqm.rho2rho(zip_rho,zip_tvec,coord_in='rho_tor', coord_out=self.rho_lbl)
                self.MDSconn.closeTree('ELECTRONS',self.shot)
                self.zipcache = zip_tvec,zip_rho, ZipTe, ZipNe
            
            if self.zipcache is []:
                return
            
            zip_tvec,zip_rho,ZipTe, ZipNe = self.zipcache
            
            imin,imax = zip_tvec.searchsorted((tmin,tmax))
            Te = median(ZipTe[imin:imax+1], 0)
            zip_rho = mean(zip_rho[imin:imax+1], 0)
            time = (tmin+tmax)/2
            
            if R is None and Z is None:
                rho = self.get_rho('',self.names,time,dR=dR,dZ=dZ)[0]
            else:
                rho = self.eqm.rz2rho(R,Z,time,self.rho_lbl)
    
            
            Te_ece = interp(abs(rho), zip_rho,Te )
        except Exception as e:
            print('Zipfit error', e)
            return 
        
        return rho,Te_ece
        
        
    
    def ne_crit(self):
        
        from scipy.constants import m_e, e, c, epsilon_0
        self.get_Te0()
        tvec,rho,_,n_e = self.zipcache

        #everything in GHz
        f_CE=ece['freq']/2e9
        #f_RHC=(0.5*f_CE)+sqrt((0.5*f_CE)**2 + n_e*e**2/m_e/epsilon_0/(2*pi*1E9)**2)
        C2=e**2/m_e/epsilon_0/(2*pi*1E9)**2
        n_e_cutoff = 8*(f_CE/2)**2/C2

        plot(rho,n_e )
        plot(rho,n_e_cutoff)
        show()
        
        
    
    def get_rho(self,group,names,time,dR=0,dZ=0):
        time = clip(time, self.tmin, self.tmax)
        R,z,theta = self.get_RZ_theta(time,names,dR=dR,dZ=dZ)
        
        rho = super(loader_ECE,self).get_rho(time,R,z,dR=dR,dZ=dZ)[0]

        r0 = interp(time, self.eqm.t_eq, self.eqm.ssq['Rmag'])+dR
        R = atleast_1d(R)
        rho[R < r0] *= -1
 
        return rho,theta,R,z
    
    def get_phi_tor(self,name=None):
        return deg2rad(self.Phi)
        
    def signal_info(self,group,name,time):
        
        rho,theta,R,z = self.get_rho(group,[name,],time)
        try:
            info = 'ch: '+str(name)+'  R:%.3f m   '%R+self.rho_lbl+': %.3f'%rho
        except:
            print( 'signal_info err')
            info = ''
        return info
    
    def get_description(self,group,name):
        return 'DIII-D %d diag: %s sig: %s'%(self.shot,'ECE',name)

  

    
from matplotlib.pylab import *



def main():

    
    mds_server = "localhost"
    #mds_server = "atlas.gat.com"
    import MDSplus as mds

    MDSconn = mds.Connection(mds_server )
    from map_equ import equ_map
    eqm = equ_map(MDSconn,debug=False)
    eqm.Open(178924,diag='EFIT01' )
    MDSconn2 = mds.Connection(mds_server )


    ece = loader_ECE(178924,exp='DIII-D',eqm=eqm,rho_lbl='rho_pol',MDSconn=MDSconn2)
    #cd 
    #ece.get_RZ_theta( 3,range(1,40),dR=0,dZ=0)
    #ece.get_Te0(2,3)
    
    #print ece.get_rho('ece', ece.names,0.1)
    
    
    #MDSconn.openTree('ECE',163303)

    #tag = '\TECEF'
    #pool = Pool(5)
    #pool.map(mds_load, [ (MDSconn, "_x="+tag+"%02d"%nch) for nch in range(1,11) ])

    #import IPython
    #IPython.embed()
    #exit()

    #f = file('pokus', 'wb')
    #pickle.dump(MDSconn , f, 2)
    #f.close

    
    #TDIcall="_x="+tag+"%02d"%nch
        
    #self.MDSconn.openTree(self.tree,self.shot)
    #Te = self.MDSconn.get(TDIcall).data()
    
    
    #data_ = ece.get_signal("", ece.get_names(ece.groups[0]), tmin=2, tmax = 2.1)
    TT =T()
    data = ece.get_signal("",ece.get_names(ece.groups[0]), tmin=1, tmax = 5)
    print(T()-TT)
    exit()

    tvec, sig = data
    plot(tvec,sig)
    xlim(tvec[0],tvec[-1])
    show()
      
      
      
    
    import IPython
    IPython.embed()
    

    
    data = array([d for t,d in data_])
    tvec = ece.tvec
    


    
    
    
    
    
    
    
    
    
    
    exit()
    
    g = sxr.groups[0]
    #print sxr.groups
    n = sxr.get_names(g)
    data1 = sxr.get_signal( '45R1',list(range(1,13)),tmin=-infty, tmax=infty,calib=True)
    data2 = sxr.get_signal('165R1',list(range(1,13)),tmin=-infty, tmax=infty,calib=True)
    data3 = sxr.get_signal('195R1',list(range(1,13)),tmin=-infty, tmax=infty,calib=True)


    data = array([d for t,d in data1]+[d for t,d in data2]+[d for t,d in data3])
    tvec = data1[0][0]
    
    
    
    

    offset = tvec < .1
    data_ = data[:,offset]-data[:,offset].mean(1)[:,None]
    u1,s1,v1 = linalg.svd(data_[:12], full_matrices=False)
    u2,s2,v2 = linalg.svd(data_[12:24], full_matrices=False)
    u3,s3,v3 = linalg.svd(data_[24:], full_matrices=False)

    #plot(range(28),u1[:,0] )
    #plot(range(28,28*2), u2[:,0])
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
    #print tvec.shape, sig.shape
    plt.plot(tvec, sig)
    plt.show()
    
if __name__ == "__main__":
    main()
    
