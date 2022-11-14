from loaders_NSTX.loader import *
#from loader import *
import os
from multiprocessing import  Pool

from time import time as T
import time
import MDSplus as mds

def check(shot):
    #fastest check if the shotfile exist
    #BUG 
    return True


def mds_load(par):
    (mds_server, shot,tree, TDI) = par
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



class loader_Magnetics(loader):
    mode_range = (-5,4)

    tor_mode_num = True
    pol_mode_num = True

    def __init__(self,*args, **kargs):
        
        super(loader_Magnetics,self).__init__(*args, **kargs)

        self.groups = ('HF','HFP','HN')
        self.names = {}
        self.coil_info = {}

        #self.tor_num_names = ['HN']
        
        self.tvec = {}
  
        self.cache = {}

        
        #self.names = {\
            #'hi-freq':['B%d'%i for i in range(1,9)]+['DXD'],\
            #'Midplane': ['MPI66M'+s for s in ('067D','097D','137D','157D','277D','307D','322B','322C','322D','322E','340D')],\
            #'Below mdpl': ['MPI6FB'+s for s in ('322D','327D','333D','342d','003D','012d','067D','157d')],\
            #'Lower div.':['MPI9B'+s for s in ('007D','037D','052D','067D','157D','x187d','Y187D','z187d','217D','232D','322D')],\
            #'Poloidal':['MPI11M322D','MPI1A322D','MPI2A322D','MPI3A322D','MPI4A322D','MPI5A322D','MPI67A322D',
                        #'MPI6FA322D','MPI6NA322D','MPI79F322D','MPI79N322D','MPI7FA322D','MPI7NA322D',
                        #'MPI8A322D','MPI89A322D','MPI9A322D','MPI66M322D','MPI1B322D','MPI2B322D',
                        #'MPI3B322D','MPI4B322D','MPI5B322D','MPI67B322D','MPI6FB322D','MPI6NB322D',
                        #'MPI79B322D','MPI7FB322D','MPI7NB322D','MPI8B322D','MPI89B322D','MPI9B322D']}

        #self.active = {}
        
       
        #self.tor_num_phi = {\
            #'Toroidal':  (67,97,137, 157,277,307,322,340),\
            #'Lower div.':(7,37,52,67,157,217,232,322)}
        
        ##poloidal coils with the best SNR 
        #self.pol_num_names = {'Poloidal':('MPI11M322D','MPI1A322D','MPI2A322D','MPI3A322D','MPI67A322D',
                 #'MPI6FA322D','MPI6NA322D','MPI7FA322D','MPI7NA322D','MPI66M322D','MPI1B322D','MPI2B322D',
                 #'MPI67B322D','MPI6FB322D','MPI6NB322D','MPI7FB322D','MPI7NB322D')}
        
        
        self.tor_num_names = {'toroidal':'HN'}
            #'Toroidal':  ['MPI66M%.3dD'%n for n in self.tor_num_phi['Toroidal']],\
            #'Lower div.':['MPI9B%.3dD'%n  for n in self.tor_num_phi['Lower div.']]}
   
        
        #self.Phi = {
             #'hi-freq': [135,]*2+[150]*6+[135,],
             #'Midplane': (67,97,137,157,277,307,322,322,322,322,340),\
             #'Below mdpl':(322,327,333,342,3,12,67,157),\
             #'Lower div.':(7,37,52,67,157,187,187,187,217,232,322),\
             #'Poloidal' : [322,]*31,\
            #}
        
  


    
    def get_names(self,group):
        group = group.lower()
        path = os.path.dirname(os.path.realpath(__file__))
        
        if group not in self.names:
            
    
            shots, files = loadtxt(path+'/mm/config_'+group+'.mm', skiprows=1, comments= '!',
                                dtype={'names': ('shot', 'file'),'formats': ( 'i4', 'U15')}, unpack=True)
            
            
            for ind, shot in enumerate(shots):
                if self.shot >= shot:
                    ifi = ind
                    
            try:
                with open(path+'/mm/'+files[ifi]) as f:
                    lines = f.readlines()
            except:
                print('File with magnetic configuration was not found '+files[ifi])
                return []
            
            if group == 'hf':
                ntor, nfp, npol = int_(lines[0].split())
            else:
                ntor,nfp, npol = int(lines[0]),0,0
    
                
            na = float(lines[1])
            tree = lines[2]
                
            nch = ntor + nfp + npol
            
            chn = zeros(nch, dtype=int)
            tor = zeros(nch, dtype=float)
            pol = zeros(nch, dtype=float)
            sig = zeros(nch, dtype=int)
            tdl  = zeros(nch, dtype=float)
            cnam = empty(nch, dtype=object)

            for i in range(nch):
                index = i + 3
                if group == 'hf':
                    if self.shot < 200000:
                        to, po, sign, chan, point = lines[i+3].split()
                        td = 0.0
                    else:
                        to, po, sign, td, chan, point = lines[i+3].split()
                else:
                    to, sign, chan, point = lines[i+3].split()
                    po = -33.47
                    td = 0

                chn[i] = int(chan)
                cnam[i] = point
                sig[i] = float(sign)
                tor[i] = float(to)
                pol[i] = float(po)
                tdl[i] = int(td)

            oe = int_(lines[-1].split())
            
            self.coil_info[group] = {'chn':chn,'cnam':cnam,'sig': sig, 'tor':tor, 
                                     'pol':pol, 'tdl':tdl, 'tree':tree, 'na':na, 'oe': oe }
            
            
            self.names[group] = arange(nch)
            
        
        return self.names[group]
 
    
    def get_signal(self,group, names,calib=False,tmin=None,tmax=None):
        
        if tmin is None:    tmin = self.tmin
        if tmax is None:    tmax = self.tmax
        group = group.lower()
        
        if size(names) == 1: names = (names, )
        names_to_load = []
        try:
            cache = load('cache_%d.npz'%self.shot, allow_pickle=True)
            cache = {'tvec':cache['tvec'].item(), 'cache':cache['cache'].item()}
        except:
            cache = {'tvec':{}, 'cache':{}}
            
        #print('cached',list(cache['cache'].keys()))
        
        self.cache.update(cache['cache'])
        self.tvec.update(cache['tvec'])

            
        for n in names:
            if not (group,n) in self.cache:
                names_to_load.append(n)
        #print(names_to_load)
        TDIcall = [self.coil_info[group]['cnam'][i] for i in names_to_load]
        #TDIcall = 'PTDATA2("%s",%d)'
        #exit()
        #load data in paraell 
        #server = self.MDSconn.hostspec
        #args = []
        if len(names_to_load):
            t = T()

            self.MDSconn.openTree(self.coil_info[group]['tree'], self.shot)
            #print('connected')
            if group not in self.tvec:
                
                self.tvec[group] = self.MDSconn.get('dim_of('+TDIcall[0]+')').data()
                #print('tvec', len(elf.tvec[group]))

            for n, tdi in zip(names_to_load,  TDIcall):
                #print(n)
                try:
                    data = self.MDSconn.get(tdi, timeout=1000000).data()
                except:
                    continue
                    #try:
                        #time.sleep(1)
                        #data = self.MDSconn.get(tdi, timeout=1000000).data()
                    #except:
                        #print('Not loaded', group, n, tdi)
                        #continue
                    #import IPython
                    #IPython.embed()
                    #if calib:
                data *= self.coil_info[group]['sig'][n]/self.coil_info[group]['na']
                self.cache[(group, n)] = data
                #print(self.cache[(group, n)].dtype, self.tvec[group].dtype)
            try:
                self.MDSconn.closeTree(self.coil_info[group]['tree'], self.shot)
            except:
                pass
            print(( 'data loaded in %.2f'%( T()-t)))
            #savez_compressed('cache_%d'%self.shot,tvec = self.tvec, cache = self.cache )

            #args.append((server, self.shot, TDIcall%(s,self.shot)))
            #args.append((server, self.shot, 'dim_of(_x)'))  #timevectors can be different for each coil

                        
        #args = [(server, self.shot, TDIcall%(s,self.shot)) for s in names_to_load]
        #print( '\n fast parallel fetch...')

        #if not group in self.tvec:        #load tvec if necessary
        #args+= [(server, self.shot,'dim_of(_x)'),] [(server, self.shot, TDIcall%(s,self.shot)) for s in names_to_load]
        #args = array_split(args,  len(names_to_load))
        #nconn = len(args)  #brutal force
        #pool = Pool()
        #out = pool.map(mds_load,args)
        #pool.close()
        #pool.join()
    

        #if not group in self.tvec: 
            #self.tvec[group][] = out[-1]
            #self.tvec[group]/= 1e3  #
            
        #tvec = self.tvec[group]
        
        #for n,sig in zip(names_to_load, out):
            #print n,not all(sig == sig[0]), shape(sig)
            #self.active[n] = not all(sig == sig[0])
            ##if len(sig) < 2 : raise Exception('Signal %s was not found'%n) 
            #self.cache[n] = sig
        #from IPython import embed
        #embed()
        #for n, (tvec, sig) in zip(names_to_load, out):
            #if len(tvec) == 0:
                #print('Signal %s was not found'%n)
                #self.active[n] = False
                #continue
            
            #self.active[n] = not all(sig == sig[0])
             
            ##------------------------------ SPECIAL FIX FOR SIGN ERROR DURING 2011 STARTUP
            #name_inv=['MPI66M067E','MPI66M097E','MPI66M127E','MPI66M157E', 
                #'MPI66M247E','MPI66M277E','MPI66M307E','MPI66M340E'  ]
            #if n.upper() in name_inv and 143400 <= self.shot <=  144100:
                #sig*= -1
 
            #if self.active[n]:

                ##if len(sig) < 2 : raise Exception('Signal %s was not found'%n)
                ##if len(sig) != len(tvec) : raise Exception('signal length do not match length of group timevec'%n)
                #tvec/= 1e3 #s
                
                #self.cache[n] = tvec,sig
        
        data = []
        for n in names:
            #if self.active[n]:
            sig = self.cache[(group, n)]
            tvec = self.tvec[group]
            #print tvec.shape, sig.shape
            imin,imax = tvec.searchsorted([tmin,tmax])
            #print imin,imax,tvec[0],  tvec[-1], std(sig)
            ind = slice(imin,imax+1)
            tvec, sig = tvec[ind], sig[ind] 
            data.append((tvec, sig))
        
        #print tvec.shape, [self.cache[n].shape for n in names]
        
        
        #data = vstack([self.cache['group'][n][ind] for n in names]).T
        #tvec = tvec[ind]
        
        if len(data) == 1:
            #print(( data[0][0].shape, data[0][1].shape ))
            #if not self.active[names[0]]:
                #raise Exception('Broken coil, use another one')
            return data[0]
        
        
        if len(data) == 0:
            raise Exception('Broken coils, use another one')
            
        #print tvec.shape, shape(data)
        return data

    

        
        
    def get_names_phase(self):
        return ('toroidal',)

        
    def get_signal_phase(self,name,calib=False,tmin=None,tmax=None,rem_elms=True):
  
        #if name in self.tor_num_names:
        group = self.tor_num_names[name]
        coils = self.get_names(group)
        #else:
            #raise Exception('Not in tor_num_names')
  
        data = self.get_signal(group, coils,calib=calib,tmin=tmin,tmax=tmax)
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
        group = self.tor_num_names[name]
        group = group.lower()

        
        return deg2rad(self.coil_info[group]['tor'])
            
    #def get_theta_pol(self,name,tshot = 4 ,rhop=.2 ):

        
        #try:
            #magr,magz, theta0, theta_star = self.mag_theta_star( tshot,rhop )
        #except:
            #print( 'Error: get_theta_pol: mag_theta_star do not work!')
            #print( traceback.format_exc())
            #theta0     = linspace(-pi,pi,10)
            #theta_star = linspace(-pi,pi,10)

            
        #R0 = interp(tshot, self.eqm.t_eq, self.eqm.ssq['Rmag'])
        #Z0 = interp(tshot, self.eqm.t_eq, self.eqm.ssq['Zmag'])
        #theta = [] 

        #for n in self.pol_num_names[name]:
            #R,Z, Phi, Tilt, L, W = BpDot_probes_322[n]
            #print(( n, name, self.active[n]))
            #if self.active[n]:
                #theta.append(arctan2(Z-Z0, R-R0 ))
    
            
        #theta = interp(theta,  theta0, theta_star)

        #return -array(theta) #BUG minus is there to have a positive m numbers for ordinary MHD modes with coinjected NBI 
            
        
        
            
            
    def signal_info(self,group,name,time):
        print(group,name,time)
        if group is None or name == 'toroidal':
            return ''
        #if name == 'Mcoils' or self.old_shotfile: return ' '
        group = group.lower()

        #calib_shot = self.dd.cShotNr('AUGD', 'CMH',self.shot)
        #self.dd.Open('CMH',calib_shot)
            
        #theta = self.dd.GetParameter('C'+name,'theta')
        #self.dd.Close()
        #try:
            #info = str(name)+' theta: %.2fdeg'%(theta*180/pi)
        #except:
            #print group,name,time,theta
            #raise
        #return deg2rad(self.coil_info[name]['tor'])

        return 'Tor: %.1f, Pol: %.1f'%(self.coil_info[group]['tor'][name], self.coil_info[group]['pol'][name])
    
 
 
def main():
    

 
    mds_server = "localhost"
    mds_server = "skylark.pppl.gov"

    import MDSplus as mds
    MDSconn = mds.Connection(mds_server )
    
    #MDSconn.openTree('d3d',-1)

    #diags = MDSconn.get('getnci(".SX*:*","node_name")')
    
    
    from map_equ import equ_map
    eqm = equ_map(MDSconn,debug=False)
    eqm.Open(121009,diag='EFIT01' )
    sxr = loader_Magnetics(121009,exp='NSTX',eqm=eqm,rho_lbl='rho_pol',MDSconn=MDSconn)
    
    g = sxr.groups[2]
    #print(( sxr.groups))
    n = sxr.get_names(g)
    data1 = sxr.get_signal(  g,n,tmin=-infty, tmax=infty,calib=True)
    data = array([d for t,d in data1])

    tvec = data1[0][0]
    
    
    
    #import IPython
    #IPython.embed()
    
        

if __name__ == "__main__":
    main()
