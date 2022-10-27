from loaders_DIIID.loader import * 
#from loader import * 


import os
from multiprocessing import  Pool

from time import time as T

import MDSplus as mds

def check(shot):
    #fastest check if the shotfile exist
    #BUG 
    return True


def mds_load(par):
    (mds_server, shot,  TDI) = par
    MDSconn = mds.Connection(mds_server )
    #print TDI
    try:
        data = MDSconn.get('_x='+TDI).data()
        
        if len(data) > 1:
            tvec = MDSconn.get('dim_of(_x)').data()
        else:
            tvec = data
    except:
        return [] ,[]
    
    return tvec,data



#			-----------------
#			 TAG DEFINITIONS
#			-----------------
#
# For magnetic sensors, the elements of the structure are
#
#{ name: R:0.0, Z:0.0, Phi:0.0, Tilt:0.0, L:0.0, W:0.0}
#
# where	
#	type		= Bp, Br, or Bt
#	name		= pointname for data acquisition
# 	R (m) 		= radius of center of probe
# 	Z (m) 		= vertical position of center of probe
# 	Phi (deg)	= toroidal angle of center of probe
# 	Tilt (deg)	= attitude angle of probe axis, 
# 			    relative to radial direction (Bp and Br probes)
# 			    relative to vertical 	 (Bt probes)
# 	L (m)		= length of probe,
# 			    in the poloidal plane	 (Bp and Br probes)
# 			    in the toroidal direction	 (Bt probes)
# 	W (deg)		= width of probe in the toroidal direction
# 

BpDot_probes_322 = { 
'MPI11M322D':( 0.973, -0.002, 322.5,  89.9, 0.115, 0.0), 
 'MPI1A322D':(  0.974,  0.182, 322.5,  90.0, 0.140, 0.0), 
 'MPI2A322D':(  0.974,  0.512, 322.5,  89.7, 0.140, 0.0), 
 'MPI3A322D':(  0.975,  0.850, 322.5,  90.2, 0.140, 0.0), 
 'MPI4A322D':(  0.972,  1.161, 322.5,  90.5, 0.140, 0.0), 
 'MPI5A322D':(  1.051,  1.330, 322.5,  44.6, 0.141, 0.0), 
 'MPI8A322D':(  1.219,  1.406, 322.5,   0.3, 0.139, 0.0), 
 'MPI89A322D':( 1.402,  1.407, 322.5,   0.7, 0.115, 0.0), 
 'MPI9A322D':(  1.584,  1.408, 322.5,   0.6, 0.140, 0.0), 
 'MPI79F322D':(1.783,  1.323, 322.5, -39.3, 0.141, 0.0), 
 'MPI79N322D':(1.924,  1.206, 322.5, -39.2, 0.114, 0.0), 
 'MPI7FA322D':( 2.067,  1.090, 322.5, -39.4, 0.140, 0.0), 
 'MPI7NA322D':( 2.219,  0.870, 323.3, -68.1, 0.140, 0.0), 
 'MPI67A322D':( 2.270,  0.746, 321.7, -67.9, 0.141, 0.0), 
 'MPI6FA322D':( 2.319,  0.623, 323.2, -68.0, 0.140, 0.0), 
 'MPI6NA322D':( 2.416,  0.249, 317.4, -89.3, 0.141, 0.0), 
 'MPI66M322D':( 2.418, -0.001, 317.4, -89.9, 0.140, 0.0), 
 'MPI1B322D':(  0.974, -0.187, 322.5,  90.1, 0.140, 0.0), 
 'MPI2B322D':(  0.975, -0.512, 322.5,  90.3, 0.140, 0.0), 
 'MPI3B322D':(  0.974, -0.854, 322.5,  89.8, 0.140, 0.0), 
 'MPI4B322D':(  0.972, -1.159, 322.5,  89.6, 0.141, 0.0), 
 'MPI5B322D':(  1.048, -1.330, 322.5, 136.4, 0.140, 0.0), 
 'MPI8B322D':(  1.254, -1.405, 322.5,-180.0, 0.140, 0.0), 
 'MPI89B322D':( 1.477, -1.406, 322.5,-179.9, 0.140, 0.0), 
 'MPI9B322D':(  1.699, -1.406, 322.5, 179.9, 0.142, 0.0), 
 'MPI79B322D':( 1.894, -1.333, 322.5,-129.5, 0.141, 0.0), 
 'MPI7FB322D':( 2.085, -1.102, 322.5,-129.2, 0.140, 0.0), 
 'MPI7NB322D':( 2.212, -0.873, 323.4,-113.5, 0.141, 0.0), 
 'MPI67B322D':( 2.263, -0.749, 321.7,-113.0, 0.141, 0.0), 
 'MPI6FB322D':( 2.315, -0.624, 323.3,-113.4, 0.140, 0.0), 
 'MPI6NB322D':( 2.416, -0.244, 317.4, -90.8, 0.140, 0.0)}


BpDot_probes_R0 = {
'MPI66M020D':(  2.410,  0.002,  19.5, -89.9, 0.139, 0.0), 
'MPI66M067D':(  2.413,  0.003,  67.5, -89.9, 0.141, 0.0), 
'MPI66M097D':(  2.413, -0.005,  97.4, -89.8, 0.137, 0.0), 
'MPI66M127D':(  2.413,  0.002, 127.9, -90.0, 0.140, 0.0), 
'MPI66M132D':(  2.409,  0.007, 132.5, -90.3, 0.055, 0.0), 
'MPI66M137D':(  2.416,  0.006, 137.4, -90.4, 0.054, 0.0), 
'MPI66M157D':(  2.413, -0.001, 157.6, -89.8, 0.137, 0.0), 
'MPI66M200D':(  2.412,  0.003, 199.7, -89.9, 0.139, 0.0), 
'MPI66M247D':(  2.413, -0.003, 246.4, -90.5, 0.140, 0.0), 
'MPI66M277D':(  2.413, -0.009, 277.5, -89.7, 0.137, 0.0), 
'MPI66M307D':(  2.413,  0.001, 307.0, -90.2, 0.140, 0.0), 
'MPI66M312D':(  2.411, -0.001, 312.4, -90.0, 0.054, 0.0), 
'MPI66M322D':(  2.418, -0.001, 317.4, -89.9, 0.140, 0.0), 
'MPI66M340D':(  2.413, -0.002, 339.7, -89.8, 0.140, 0.0)} 



class loader_Magnetics(loader):
    mode_range = (-5,4)

    tor_mode_num = True
    pol_mode_num = True

    def __init__(self,*args, **kargs):
        
        super(loader_Magnetics,self).__init__(*args, **kargs)

        self.groups = ('hi-freq', 'Midplane','Below mdpl', 'Lower div.','Poloidal')
        
        
        self.names = {\
            'hi-freq':['B%d'%i for i in range(1,9)]+['DXD'],\
            'Midplane': ['MPI66M'+s for s in ('067D','097D','137D','157D','277D','307D','322B','322C','322D','322E','340D')],\
            'Below mdpl': ['MPI6FB'+s for s in ('322D','327D','333D','342d','003D','012d','067D','157d')],\
            'Lower div.':['MPI9B'+s for s in ('007D','037D','052D','067D','157D','x187d','Y187D','z187d','217D','232D','322D')],\
            'Poloidal':['MPI11M322D','MPI1A322D','MPI2A322D','MPI3A322D','MPI4A322D','MPI5A322D','MPI67A322D',
                        'MPI6FA322D','MPI6NA322D','MPI79F322D','MPI79N322D','MPI7FA322D','MPI7NA322D',
                        'MPI8A322D','MPI89A322D','MPI9A322D','MPI66M322D','MPI1B322D','MPI2B322D',
                        'MPI3B322D','MPI4B322D','MPI5B322D','MPI67B322D','MPI6FB322D','MPI6NB322D',
                        'MPI79B322D','MPI7FB322D','MPI7NB322D','MPI8B322D','MPI89B322D','MPI9B322D']}

        self.active = {}
        
       
        self.tor_num_phi = {\
            'Toroidal':  (67,97,137, 157,277,307,322,340),\
            'Lower div.':(7,37,52,67,157,217,232,322)}
        
        #poloidal coils with the best SNR 
        self.pol_num_names = {'Poloidal':('MPI11M322D','MPI1A322D','MPI2A322D','MPI3A322D','MPI67A322D',
                 'MPI6FA322D','MPI6NA322D','MPI7FA322D','MPI7NA322D','MPI66M322D','MPI1B322D','MPI2B322D',
                 'MPI67B322D','MPI6FB322D','MPI6NB322D','MPI7FB322D','MPI7NB322D')}
        
        
        self.tor_num_names = {\
            'Toroidal':  ['MPI66M%.3dD'%n for n in self.tor_num_phi['Toroidal']],\
            'Lower div.':['MPI9B%.3dD'%n  for n in self.tor_num_phi['Lower div.']]}
   
        
        self.Phi = {
             'hi-freq': [135,]*2+[150]*6+[135,],
             'Midplane': (67,97,137,157,277,307,322,322,322,322,340),\
             'Below mdpl':(322,327,333,342,3,12,67,157),\
             'Lower div.':(7,37,52,67,157,187,187,187,217,232,322),\
             'Poloidal' : [322,]*31,\
            }
        
  
        self.tvec = {}
  
        self.catch = {}


    
    def get_names(self,group):
        return self.names[group]
    
    #def get_signal(self,group, names,calib=False,tmin=None,tmax=None):
        
        #if tmin is None:    tmin = self.tmin
        #if tmax is None:    tmax = self.tmax
        
        
        
        #if not group in self.tvec:
            #self.tvec = {}
        
        #if size(names) == 1: names = (names, )
        #names_to_load = []
        #for n in names:
            #if not n in self.catch:
                #names_to_load.append(n)
                

        
        #TDIcall = 'PTDATA2("%s",%d)'
        
        ##load data in paraell 
        #t = time()
        #server = self.MDSconn.hostspec
        
        #args = [(server, self.shot, TDIcall%(s,self.shot)) for s in names_to_load]
        
        #if not group in self.tvec:        #load tvec if necessary
           #args+= [(server, self.shot,'dim_of('+args[0][2]+')'),] 
        
        #nconn = len(args)  #brutal force
        #pool = Pool()
        #out = pool.map(mds_load,args)
        #pool.close()
        #pool.join()
    
        #print 'data loaded in %.2f'%( time()-t)
        

        #if not group in self.tvec: 
            #self.tvec[group] = out[-1]
            #self.tvec[group]/= 1e3  #
        ##print names_to_load
        #for n,sig in zip(names_to_load, out):
            #print n,not all(sig == sig[0]), shape(sig)
            #self.active[n] = not all(sig == sig[0])
            ##if len(sig) < 2 : raise Exception('Signal %s was not found'%n) 
            #self.catch[n] = sig
            
        #tvec = self.tvec[group]
        
        
        #imin,imax = tvec.searchsorted([tmin,tmax])
        #ind = slice(imin,imax+1)
        
        #data = []
        #for n in names:
            #if self.active[n]:
                #data.append(self.catch[n][ind])
        #if len(data) == 0:
            #raise Exception('Broken coil, use another one')
            
        #data = squeeze(vstack(data).T)

        #return tvec[ind] ,data

    
    def get_signal(self,group, names,calib=False,tmin=None,tmax=None):
        
        if tmin is None:    tmin = self.tmin
        if tmax is None:    tmax = self.tmax

        if size(names) == 1: names = (names, )
        names_to_load = []

        for n in names:
            if not n in self.catch:
                names_to_load.append(n)
                
        if len(names_to_load):
            
            TDIcall = 'PTDATA2("%s",%d)'
            
            #load data in paraell 
            t = T()
            server = self.MDSconn.hostspec
            args = []
            for s in names_to_load:
                args.append((server, self.shot, TDIcall%(s,self.shot)))
                            
            print( '\nParallel fetch...')
            pool = Pool()
            out = pool.map(mds_load,args)
            pool.close()
            pool.join()
            
        else:
            out = []

        for n, (tvec, sig) in zip(names_to_load, out):
            if len(tvec) == 0:
                print('Signal %s was not found'%n)
                self.active[n] = False
                continue
            
            self.active[n] = not all(sig == sig[0])
             
            #------------------------------ SPECIAL FIX FOR SIGN ERROR DURING 2011 STARTUP
            name_inv=['MPI66M067E','MPI66M097E','MPI66M127E','MPI66M157E', 
                'MPI66M247E','MPI66M277E','MPI66M307E','MPI66M340E'  ]
            if n.upper() in name_inv and 143400 <= self.shot <=  144100:
                sig*= -1
 
            if self.active[n]:
                sig = self.remove_elms(tvec, sig)
                tvec/= 1e3 #s
                self.catch[n] = tvec,sig
        
        data = []
        for n in names:
            if self.active[n]:
                tvec, sig = self.catch[n]
                imin,imax = tvec.searchsorted([tmin,tmax])
                ind = slice(imin,imax+1)
                tvec, sig = tvec[ind], sig[ind] 
                data.append((tvec, sig))

        if len(data) == 1:
            return data[0]
        
        
        if len(data) == 0:
            raise Exception('Broken coils, use another one')
            
       
        return data

    

    def remove_elms(self, tvec, signal):
        if not hasattr(self,'elm_start'):
            try:
                #load ELMs from Tom Osborns MDS+ tree
                self.MDSconn.openTree('PEDESTAL', self.shot)
                
                #TDI = r'\PEDESTAL::TOP.ELM:ELMDANAME'
                #fname = self.MDSconn.get(TDI).data()
                
                TDI = r'\PEDESTAL::TOP.ELM:ELMSTART'
                elm_val = self.MDSconn.get(r'_x=\PEDESTAL::TOP.ELM:ELMSTART').data()
                elm_start = self.MDSconn.get(r'dim_of(_x)').data()
                TDI = r'\PEDESTAL::TOP.ELM:ELMEND'
                elm_end = self.MDSconn.get('dim_of('+TDI+')').data()
                TDI = r'\PEDESTAL::TOP.ELM:ELMPEAK'
                elm_peak = self.MDSconn.get('dim_of('+TDI+')').data()
                
                
        
                #estimate size of the ELM
                area = (elm_end-elm_start)*elm_val
                ind = area > 1e15  #BUG hardcoded value
                
                self.elm_start = elm_start[ind]
                self.elm_peak = elm_peak[ind]
                self.elm_end = elm_end[ind]
                self.elm_val = elm_val[ind]
                
                

                #self.MDSconn.openTree('SPECTROSCOPY', self.shot)
                #TDI = '_x=\\SPECTROSCOPY::'+ fname
                #filterscope = self.MDSconn.get(TDI).data()
                #filterscope_t = self.MDSconn.get('dim_of(_x)').data()
 
            except Exception as e:
                print('ELM detection issue: ', e)
                return signal
     
        
        valid = (self.elm_start > tvec[0])&(self.elm_start < tvec[-1])
        
        #delete also 0.1ms before and after due to uncertainties in elm time
        ind_start = tvec.searchsorted(self.elm_start[valid]-0.1)
        ind_end = tvec.searchsorted(self.elm_peak[valid]+0.1)
        
 
        corrected_signal = copy(signal)
        for i1,i2 in zip(ind_start, ind_end):
            corrected_signal[i1:i2] = signal[i1:i2].mean(0)
            
            
        #plot(tvec, signal)
        #plot(tvec, corrected_signal)
        #plot(filterscope_t,filterscope/1e15 )
        #plot(c_[self.elm_start, self.elm_peak, self.elm_end].flatten(), c_[0*self.elm_val,self.elm_val/1e15,self.elm_val*0].flatten())        
        #show()

        return corrected_signal
            
    
        
        
    def get_names_phase(self):
        return ('Toroidal','Poloidal')
        
 
    
        
    def get_signal_phase(self,name,calib=False,tmin=None,tmax=None,rem_elms=True):
  
        if name in self.tor_num_names:
            coils = self.tor_num_names[name]
        elif name in self.pol_num_names:
            coils = self.pol_num_names[name]
  
        data = self.get_signal(name, coils,calib=calib,tmin=tmin,tmax=tmax)
        tvec = []
        for t,d in data:
            if len(t) > len(tvec): tvec = t
        output = empty((len(tvec), len(data)),dtype='single')
        
        for i,(t,d) in enumerate(data):
            if len(tvec)!= len(t):
                output[:,i] = interp(tvec, t, d, left=0, right=0)
            else:
                output[:,i] = d

        
            
        #print data.shape

        return tvec, output

    def get_phi_tor(self,name):

        phi_tor = []
        for ic,c in enumerate(self.tor_num_names[name]):
            if c in self.active and self.active[c]:
                phi_tor.append(self.tor_num_phi[name][ic])
        
        
        return deg2rad(phi_tor)
            
    def get_theta_pol(self,name,tshot = 4 ,rhop=.2 ):

        
        try:
            magr,magz, theta0, theta_star = self.mag_theta_star( tshot,rhop )
        except:
            print( 'Error: get_theta_pol: mag_theta_star do not work!')
            print( traceback.format_exc())
            theta0     = linspace(-pi,pi,10)
            theta_star = linspace(-pi,pi,10)
            
            
        self.eqm.read_ssq()
        R0 = interp(tshot, self.eqm.t_eq, self.eqm.ssq['Rmag'])
        Z0 = interp(tshot, self.eqm.t_eq, self.eqm.ssq['Zmag'])
        theta = [] 

        for n in self.pol_num_names[name]:
            R,Z, Phi, Tilt, L, W = BpDot_probes_322[n]
            #print(( n, name, self.active[n]))
            if self.active[n]:
                theta.append(arctan2(Z-Z0, R-R0 ))
    
            
        theta = interp(theta,  theta0, theta_star)

        return -array(theta) #BUG minus is there to have a positive m numbers for ordinary MHD modes with coinjected NBI 
            
        
        
            
            
    def signal_info(self,group,name,time):
        #if name == 'Mcoils' or self.old_shotfile: return ' '

        #calib_shot = self.dd.cShotNr('AUGD', 'CMH',self.shot)
        #self.dd.Open('CMH',calib_shot)
            
        #theta = self.dd.GetParameter('C'+name,'theta')
        #self.dd.Close()
        #try:
            #info = str(name)+' theta: %.2fdeg'%(theta*180/pi)
        #except:
            #print group,name,time,theta
            #raise
        
        return ''
    
 
 
def main():\
    

 
    mds_server = "localhost"
    #mds_server = "atlas.gat.com"

    import MDSplus as mds
    MDSconn = mds.Connection(mds_server )
    
    #MDSconn.openTree('d3d',-1)

    #diags = MDSconn.get('getnci(".SX*:*","node_name")')
    
    
    from map_equ import equ_map
    eqm = equ_map(MDSconn,debug=False)
    eqm.Open(175860,diag='EFIT01' )
    coils = loader_Magnetics(175860,exp='DIII-D',eqm=eqm,rho_lbl='rho_pol',MDSconn=MDSconn)
    
    g = coils.groups[0]
    print(( coils.groups))
    n = coils.get_names(g)
    data1 = coils.get_signal(  g,n,tmin=-infty, tmax=infty,calib=True)
    data = array([d for t,d in data1])

    tvec = data1[0][0]
    
    
    
    import IPython
    IPython.embed()
    
        

if __name__ == "__main__":
    main()
