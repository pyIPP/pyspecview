try:
    from loaders_DIIID.loader import * 
except:
    from loader import * 

import os
from multiprocessing import  Pool
import MDSplus as mds
import traceback
 


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
    TT = T()
    data = []
    for tdi in TDI:
        try:
            data.append(MDSconn.get(tdi).data())
        except:
            data.append([])
    return data


def mds_load_par( MDSconn, TDI, numTasks=8):
    
    if len(TDI) <= 4:
        return [MDSconn.get(tdi).data() for tdi in TDI]
    
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
class loader_ECEI(loader):
    

    radial_profile=True
    units = 'eV'
 

    def __init__(self,*args, **kargs):
        
        super(loader_ECEI,self).__init__(*args, **kargs)
   
        self.tree = 'ELECTRONS'
        self.MDSconn.openTree(self.tree,self.shot)
        self.time_header = None
        self.data_dict = {}
        self.groups = []
        self.tvec = {}
        self.goodchannels = []
        self.calib = {}

        
        #MDS values are probably not right 
        try:
            self.LFSGOOD = self.MDSconn.get('\\ECEI::TOP.LFSGOOD')
            self.groups.append('LFS')
        except:
            self.LFSGOOD = False
            
        try:
            self.HFSGOOD = self.MDSconn.get('\\ECEI::TOP.HFSGOOD')
            self.groups.append('HFS')
        except:
            self.HFSGOOD = False
    

        self.Phi = 270
        
        self.MDSconn.closeTree(self.tree, self.shot)

   
        
        def WN(data):  # wide or narrow zoom?
            return {'W':1, 'N': 0}.get(data, -1)


        def readline(line):  # digitize the settings of each day
            date = int(line[0])  # date
            ss = int(line[1])  # start shot
            es = int(line[2])  # end shot
            Op = WN(line[3])  # optics
            Lz = WN(line[4])  # LFS zoom
            Hz = WN(line[5])  # HFS zoom
            LL = float(line[6])  # LFS LO
            HL = float(line[7])  # HFS LO
            return date, ss, es, Op, Lz, Hz, LL, HL


        ###Load hardward setting record from a local file
        filename = '/fusion/projects/results/ecei/ECEI_record_uptodate'
        if not os.path.exists(filename):
            path = os.path.dirname(os.path.realpath(__file__))
            filename = path+'/'+filename.split('/')[-1]
        
        #print('Read Setting from Folder')
        valid = False
        with open(filename) as f:
            for line in f.read().splitlines():
                data = line.split()
                if len(data) == 8:  # valid records
                    data = readline(data)
                    date, ss, es, Op, Lz, Hz, LL, HL = data
                    if ss <= self.shot <= es:
                        valid = True
                        break

        if valid:
            self.LFSLO = LL
            self.HFSLO = HL
            self.Optics = 'Wide' if Op else 'Narrow'
            

            self.LFSGOOD = True
            if Lz == 0:
                self.LFSIF = 'Narrow'
            elif Lz == 1:
                self.LFSIF = 'Wide'
            else:
                self.LFSIF =  None
                self.LFSGOOD = False

            self.HFSGOOD = True
            if Hz == 0:
                self.HFSIF = 'Narrow'
            elif Hz == 1:
                self.HFSIF = 'Wide'
            else:
                self.HFSIF =  None
                self.HFSGOOD = False
                
        if not valid or not (self.HFSGOOD or self.LFSGOOD):
            print('ECEI is not available for this shot')
            self.HFSGOOD = self.LFSGOOD = False
            return

        #print('Hardware settings of LFS array is loaded')

        ###Settings of good modules
        filename = '/fusion/projects/results/ecei/ECEI_goodchannels_uptodate'
        if not os.path.exists(filename):
            path = os.path.dirname(os.path.realpath(__file__))
            filename = path+'/'+filename.split('/')[-1]
        
        with open(filename) as f:
            for line in f.read().splitlines():
                data = line.split()
                if len(data) != 0:  # valid records
                    ss = int(data[0])
                    es = int(data[1])
                    if ss <= self.shot <=  es:
                        OpN = float(data[2])
                        OpW = float(data[3])
                        self.goodchannels = [int(x) for x in data[4::]]
                        #print('Good Modules are ', self.goodchannels)
                        self.OpticsWide = OpW
                        self.OpticsNarrow = OpN
                        break

                
    def get_RZ_theta(self, time, system, names,dR=0,dZ=0):
        
        if isinstance(names, str):
            names = [names]
            
        if system == 'LFS':
            LO = self.LFSLO
            IF = self.LFSIF
            
        elif system == 'HFS':
            LO = self.HFSLO
            IF = self.HFSIF            
 
                
        if IF == 'Narrow':
            df = np.array([7.5, 6.9, 6.3, 5.7, 5.1, 4.5, 3.9, 3.3])
        else:
            df = np.array([8.9, 7.8, 6.9, 6.0, 5.1, 4.2, 3.3, 2.4])
            
        #BUG this is in source ode of ECEI module ECEIlayout.py, moreover they reversed names order :( 
 
        #deterctor frequency
        fc = LO + df
     
        if self.Optics == 'Wide':
            OpSpan = self.OpticsWide / 100 #m
        else:
            OpSpan = self.OpticsNarrow / 100 #m

        Zmod = np.linspace(-OpSpan, OpSpan, 20) 
                
        if self.time_header is not None:
            time = np.clip(time,self.time_header[0,0], self.time_header[1,-1])
            

        Rmesh = np.linspace(self.eqm.Rmesh[0], self.eqm.Rmesh[-1], 200)
        Zmesh = np.linspace(self.eqm.Zmesh[0], self.eqm.Zmesh[-1], 300)

        
        #position including relativistic shift and refraction
        try:
            B = self.eqm.rz2brzt(r_in=Rmesh, z_in=Zmesh, t_in=time)
            Btot = squeeze(linalg.norm(B,axis=0)).T
            try:
                from ._ECE_chord import ece_los
            except:
                from ._ECE_chord import ece_los

            Rcm,Zcm = np.meshgrid(Rmesh, Zmesh)
            _, Te, Ne = self.get_Te0(time,time,R=Rcm,Z=Zcm, returnTeNe=True)
            #pcolor(Rcm, Zcm, Ne)
            Rc,Zc = [],[]
            for n in names:
                z_los = Zmod[len(Zmod)-int(n[:2])+3-1]
                f_los = fc[len(fc)-int(n[2:])+1-1]*1e9
                r,z, R_res, Z_res = ece_los(f_los, z_los, Rmesh, Zmesh, Btot, Ne, Te)
                Rc.append(R_res)
                Zc.append(Z_res)
                #plot(r,z,lw=.3)
            #plot(Rc,Zc,'o')
            #ylim(Zmesh[0], Zmesh[-1])
            #oo
        except Exception as e:
            print(e)
            from scipy.constants import c, e,m_e, epsilon_0

            B = self.eqm.rz2brzt(r_in=self.eqm.Rmesh, z_in=Zmod, t_in=time)
            Btot = squeeze(linalg.norm(B,axis=0)).T
            wce = e*Btot/m_e/1e9 #GHz

            nharm = 2
    
            Rc = np.zeros((len(Zmod), len(fc)))
            for i, wce_row in enumerate(wce):
                Rc[i] = np.interp(-2*pi*fc,-wce_row*nharm,self.eqm.Rmesh)
            

            Rc = [Rc[len(Zmod)-int(n[:2])+3-1,len(fc)-int(n[2:])+1-1] for n in names]
            Zc = [Zmod[len(Zmod)-int(n[:2])+3-1] for n in names]
             
            #plot(Rc,Zc,'x')
            #show()
        
        Rc,Zc = array(Rc),array(Zc)
  
        r0 = interp(time, self.eqm.t_eq, self.eqm.ssq['Rmag'])+dR
        z0 = interp(time, self.eqm.t_eq, self.eqm.ssq['Zmag'])+dZ

        
        return Rc,Zc, arctan2(Zc-z0, Rc-r0)
    
    

    
    def get_signal(self,group, names,calib=False,tmin=None,tmax=None):
        if len(names) == 0:
            return [],[]

        if isinstance(names, str):
            names = [names]
 
      
        if tmin is None:    tmin = self.tmin
        if tmax is None:    tmax = self.tmax
 
        #it takes ~2s to fetch the header
        if self.time_header is None: 
            channel = f'{group}{names[0]}'
            #https://diii-d.gat.com/diii-d/Ptdatacall
            time_header = self.MDSconn.get(f'PTHEAD2("{channel}",{self.shot}); __real64')[2:]
            self.time_header = time_header.reshape(-1,2).T
        
        imin = max(0,self.time_header[0].searchsorted(tmin)-1)
        imax = min(len(self.time_header.T)-1,self.time_header[1].searchsorted(tmax))
        time_intervals = range(imin, imax+1)
        
        #use catch if possible 
        load_seg = []
        load_name = []
        for n in atleast_1d(names):
            for it in time_intervals:
                seg = f'{group+n}_{it}'
                if len(self.data_dict.get(seg, [])) <= 1:
                    load_seg.append(seg)
                    load_name.append(n)
        
 
        TT = T()

        #fast paraell fetch 
        if len(load_seg) > 0:
            #fetch raw data as integeres and make calibration locally
            #it allows to identify saturated timepoints
            for n in load_name:
                if n not in self.calib:
                    #https://diii-d.gat.com/diii-d/Ptdatacall
                    self.calib[n] = self.MDSconn.get(f'PTHEAD2("{group+n}",{self.shot}); __rarray').data()[4:6] 
    
            #NOTE it is faster to fetch the whole time range than each segment one by one.
            TDI = [f'PTDATA2("{s}", {self.shot}, 0)' for s in load_seg]
            out = mds_load_par(self.MDSconn, TDI)
           
            #identify saturated data and calibrate signals
            for s,n, o in zip(load_seg,load_name, out):
                if len(o) > 1:
                    saturated  = np.abs(o) == 2**15
                    o -= int32(self.calib[n][1])
                    o = float32(o)
                    o*= -float32(self.calib[n][0])
                    self.data_dict[s] = ma.array(o, mask=saturated)
                else:
                    self.data_dict[s] = o
    
        for it in time_intervals:
            if it not in self.tvec:
                for n in names:
                    nt = len(self.data_dict[ f'{group+n}_{it}'])
                    if nt > 1: break
                self.tvec[it] = np.linspace(self.time_header[0,it],self.time_header[1,it], nt)
        
        tvec = np.hstack([self.tvec[it] for it in time_intervals])
        imin,imax = tvec.searchsorted([tmin,tmax]) 
        imax += 1
        
        #collect all data
        output = []
        for n in names:
            #it could be done a bit more efficiently
            if len(time_intervals) > 1:
            	Te = np.ma.hstack([self.data_dict[f'{group+n}_{it}'] for it in time_intervals])	
            else:
                Te = self.data_dict[f'{group+n}_{time_intervals[0]}']
            
            if len(Te) != len(tvec):
                print(f'data from  channel {n} are not available')
                Te = np.zeros_like(tvec, dtype='single')
            output.append([tvec[imin:imax], Te[imin:imax]])
 
        #calibrate using zipfit data, offset is already set to be zero
        if calib:
            R,Z,Theta = self.get_RZ_theta((tmin+tmax)/2, group, names)
            rho,Te0 = self.get_Te0(tmin,tmax,R=R,Z=Z)

            if Te0 is not None:
                for out, Te, n in zip(output, Te0, names):
                    #don't use inplace operation
                    out[1] = out[1]*(Te/asarray(out[1]).mean())
            else:
                print('ECEI data are not calibrated')

        if len(names) == 1:
            return output[0]
        else:
            return output
                
        
    def get_names(self,group):
        return ['%.2d%.2d'%(gc,R) for gc in self.goodchannels for R in range(1,9)]
    
  
        
    def get_Te0(self,tmin,tmax,dR=0,dZ=0,R=None,Z=None, returnTeNe=False):
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
                
                zip_rho = self.eqm.rho2rho(zip_rho,zip_tvec,coord_in='rho_tor', coord_out=self.rho_lbl, extrapolate=True)
                self.MDSconn.closeTree('ELECTRONS',self.shot)
                self.zipcache = zip_tvec,zip_rho, ZipTe, ZipNe
            
            if self.zipcache is []:
                return
            
            zip_tvec,zip_rho,ZipTe, ZipNe = self.zipcache
            
            imin,imax = zip_tvec.searchsorted((tmin,tmax))
            Te = median(ZipTe[imin:imax+1], 0)
            Ne = median(ZipNe[imin:imax+1], 0)

            zip_rho = mean(zip_rho[imin:imax+1], 0)
            time = (tmin+tmax)/2
            
            if R is None and Z is None:
                rho = self.get_rho('',self.names,time,dR=dR,dZ=dZ)[0]
            else:
                rho = self.eqm.rz2rho(R[None],Z[None],time,self.rho_lbl)[0]

       
            Te = interp(abs(rho), zip_rho,Te )
            Ne = interp(abs(rho), zip_rho,Ne )
      

        except Exception as e:
            print('Zipfit error', e)
            return None,None
        
        if returnTeNe:
            return rho,Te, Ne
        
        return rho,Te
        
        
    
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
        
        if self.time_header is not None:
            time = np.clip(time,self.time_header[0,0], self.time_header[1,-1])
        
        R,z,theta = self.get_RZ_theta(time,group, names,dR=dR,dZ=dZ)
        rho = super(loader_ECEI,self).get_rho(time,R,z,dR=dR,dZ=dZ)[0]

        r0 = interp(time, self.eqm.t_eq, self.eqm.ssq['Rmag'])+dR
        R = atleast_1d(R)
        rho[R < r0] *= -1
  
        
        return rho,theta,R,z
    
    def get_phi_tor(self,name=None):
        return deg2rad(self.Phi)
        
    def signal_info(self,group,name,time):
        
        rho,theta,R,z = self.get_rho(group,[name,],time)
        try:
            info = 'sys: '+group+ ' ch: '+str(name)+'  R:%.3f m   '%R+self.rho_lbl+': %.3f'%rho
        except:
            print( traceback.format_exc())
            info = ''
        return info
    
    def get_description(self,group,name):
        return 'DIII-D %d diag: %s sig: %s'%(self.shot,'ECEI',name)

  

    
from matplotlib.pylab import *



def main():
   
 
    from time import time as T
    mds_server = "localhost"
    #mds_server = "atlas.gat.com"
    import MDSplus as mds

    MDSconn = mds.Connection(mds_server )
    from map_equ import equ_map
    eqm = equ_map(MDSconn,debug=False)
    eqm.Open(191823,diag='EFIT01' )
    MDSconn2 = mds.Connection(mds_server )


    ecei = loader_ECEI(191823 ,exp='DIII-D',eqm=eqm,rho_lbl='rho_pol',MDSconn=MDSconn2)
    names = ecei.get_names('LFS')
    R,Z,T = ecei.get_RZ_theta(2.1, 'LFS', names)
        #def get_RZ_theta(self, time, system, names,dR=0,dZ=0):
    plot(R,Z,'o')
    #TT =  T()
    #signals = ecei.get_signal('LFS',names, calib=True, tmin=3, tmax=3.6)
    #print(T()-TT)
    exit()
 
     
    signals = ecei.get_signal('LFS',names, calib=True, tmin=-1.1, tmax=8)
    #embed()
    for t,s in signals:
        plot(t[718:].reshape(-1,1000).mean(1),s[718:].reshape(-1,1000).mean(1))
    
    
    Te = [s.mean() for t,s in signals]
    R,Z,T = ecei.get_RZ_theta(2.1, 'LFS', names)
    
    
    #pcolor(R.reshape(-1,8),Z.reshape(-1,8), array(Te).reshape(-1,8) )
    #show()
    
    contourf(R.reshape(-1,8),Z.reshape(-1,8), array(Te).reshape(-1,8) )
    names2D = array(names).reshape(-1,8)
    for i in range(10):
        for j in range(8):
            plot(R.reshape(-1,8)[i,j],Z.reshape(-1,8)[i,j], '.')
            text(R.reshape(-1,8)[i,j],Z.reshape(-1,8)[i,j], names2D[i,j])
    
    show()
    
    
    
    print(T()-t)
    embed()
    
    
    exit()

    
    #######################################
    from time import time
    t = time()
    #mds_server = "localhost"
    mds_server = "atlas.gat.com"
    import MDSplus as mds

    MDSconn = mds.Connection(mds_server )
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
    group = 'LFS'
    names = ['0301', '0302', '0303', '0304', '0305', '0306', '0307', '0308', '0501', '0502', '0503', '0504', '0505', '0506', '0507', '0508', '0701', '0702', '0703', '0704', '0705', '0706', '0707', '0708', '1101', '1102', '1103', '1104', '1105', '1106', '1107', '1108', '1201', '1202', '1203', '1204', '1205', '1206', '1207', '1208', '1301', '1302', '1303', '1304', '1305', '1306', '1307', '1308', '1501', '1502', '1503', '1504', '1505', '1506', '1507', '1508', '1701', '1702', '1703', '1704', '1705', '1706', '1707', '1708', '1901', '1902', '1903', '1904', '1905', '1906', '1907', '1908', '2101', '2102', '2103', '2104', '2105', '2106', '2107', '2108']
    numTasks = 8
    TDI = [f'PTDATA2("LFS{n}", 191823)' for n in names]
    TDI = array_split(TDI, numTasks)
    args = [(mds_server, tdi) for tdi in TDI]
    from multiprocessing import  Pool

    out = []
    pool = Pool()
    for o in map(mds_load,args):
        out.extend(o)
    pool.close()
    pool.join()
    print('-----------',time()-t)
   
   
    #######################################

    
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
    #T =T()
    #embed()
    #data = ecei.get_signal("LFS",names, tmin=1, tmax = 5)
    #data = ecei.get_signal("LFS",names[:5], tmin=1, tmax = 5)
    ecei.get_RZ_theta(2.3, "LFS",names)
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
    
