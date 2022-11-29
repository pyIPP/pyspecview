
from numpy import *
import sys, os, random
from scipy.interpolate import interp1d
import traceback

#shot_path = '/afs/ipp-garching.mpg.de/augd/shots/'

def min_fine(x,y): #find extrem with subpixel precision
    i = argmin(y)
    if i==0 or i == len(y)-1: return y[i],x[i]    #not working at the edge! 
    A = x[i-1:i+2]**2, x[i-1:i+2], ones(3)
    a,b,c = linalg.solve(array(A).T,y[i-1:i+2])
    
    ymin =  c-b**2/(4*a)
    xmin = -b/(2*a) 
    return  ymin,xmin

class loader(object):
    pol_mode_num = False
    tor_mode_num = False
    n_mode_range = None
    radial_profile=False
    units = 'a.u.'
    #plasma_freq=None

    def __init__(self,shot,exp=None,ed=0,eqm=None,rho_lbl='rho_pol',MDSconn=None):
        
        #use DD library on AUG
        if MDSconn is None:
            import dd  
            self.dd = dd.shotfile() 
        else:
            self.MDSconn = MDSconn
            
            
        self.eqm = eqm

        self.shot = int(shot)
        self.exp = str(exp)
        self.ed = int(ed)
        self.names = {}
        self.groups = {}
        self.shotfile = 'XXX'
        self.rho_lbl= rho_lbl
        
        
        if eqm is not None:
            self.eqm.read_ssq()

        if eqm is not None and hasattr(self.eqm,'t_eq'):
                time = self.eqm.t_eq
                self.tmin = time[0]
                self.tmax = time[-1]+ (time[-1]-time[0])/len(time)
        else:
            print ('No magnetic equlibrium is avalible!')
            self.tmin = 0
            self.tmax = 10

        #print 'tmin tmax',self.tmin,self.tmax 
        #print self.eqm.t_eq, self.eqm.diag
        
             
    def get_signal_groups(self):
        return  self.groups
            
    def get_names(self,group):
        pass 
    def get_signal(self,group, name,calib=False, tmin=None, tmax=None):
        pass
    def get_names_phase(self):
        pass
    def get_signal_phase(self,name,calib=False):
        pass
    def get_phi_tor(self,name):
        pass
    def get_phase_corrections(self,name):
        pass

            
    def mag_theta_star(self,tshot = 4 ,rho=0,dR=0,dZ=0,n_theta = 100):

        rhop = self.eqm.rho2rho(rho,tshot, coord_in=self.rho_lbl, coord_out='rho_pol' )[0]
   
        single_radius = False
        
        
        if size(rhop) == 1:
            rhop = r_[rhop,rhop+1e-2]
            single_radius = True

        angle = linspace(0, 2*pi,n_theta)

        
        magr,magz = self.eqm.rhoTheta2rz(rhop,angle, tshot)
      
        magr = magr[0].T + dR
        magz = magz[0].T + dZ
 
        
        r0 = interp(tshot, self.eqm.t_eq, self.eqm.ssq['Rmag'])+dR
        z0 = interp(tshot, self.eqm.t_eq, self.eqm.ssq['Zmag'])+dZ

        n_rho, n_theta = magr.shape

        drdrho,drtheta = gradient(magr)
        dzdrho,dztheta = gradient(magz)
        dpsidrho,dpsitheta = gradient(tile(rhop**2, (n_theta,1)).T )

        grad_rho = dstack((drdrho,dzdrho,dpsidrho ))
        grad_theta = dstack((drtheta,dztheta,dpsitheta))
        normal = cross(grad_rho,grad_theta,axis=-1)

                

        dpsi_dr = -normal[:,:,0]/(normal[:,:,2]+1e-8) #Bz
        dpsi_dz = -normal[:,:,1]/(normal[:,:,2]+1e-8) #Br
 
        dtheta_star = ((magr-r0)**2+(magz-z0)**2)/(dpsi_dz*(magz-z0)+dpsi_dr*(magr-r0))/magr
    
        if not all(isfinite(dtheta_star)):
            print(( 'dtheta_star is corrupted',tshot))
            print(( sum(isfinite(dtheta_star))/float(size(dtheta_star))))

        theta = -arctan2(magz-z0, -magr+r0 )

        theta = unwrap(theta-theta[:,(0,)],axis=1)
        from scipy.integrate import cumtrapz

        #definition of the thetat star by integral
        theta_star = cumtrapz(dtheta_star,theta,axis=1,initial=0)
        correction = (n_theta-1.)/n_theta
        if sum(abs(magr[:,0]-magr[:,-1]))<1e-5 and sum(abs(magz[:,0]-magz[:,-1]))<1e-5:
            correction = 1
        theta_star/= theta_star[:,(-1,)]/(2*pi)/correction     #normalize to 2piz
        
           
        if single_radius:
            return magr[0],magz[0],angle,theta_star.mean(0)
            
        return magr,magz, angle, theta_star

    
    def get_rho(self,time,R_start,z_start,Phi_start=None, R_end=None,z_end=None,Phi_end=None, dR=0,dZ=0):
        #get rho_pol/rho_tor for points, R2, and Z2 must be specified for LOSs
        
        

        
        n_points = size(R_start)
        time = max(min(time, self.tmax), self.tmin)

        r0 = interp(time, self.eqm.t_eq, self.eqm.ssq['Rmag'])+dR
        z0 = interp(time, self.eqm.t_eq, self.eqm.ssq['Zmag'])+dZ
        
        if R_end is None or z_end is None:
            rho = self.eqm.rz2rho(R_start-dR,z_start-dZ,time, coord_out=self.rho_lbl)[0]

            theta = arctan2(z_start-z0,R_start-r0)
            return rho,theta, R_start,z_start
        

        
        Phi_start = deg2rad(Phi_start)
        Phi_end = deg2rad(Phi_end)

        Start = array((R_start*cos(Phi_start),R_start*sin(Phi_start), z_start))
        End   = array((R_end  *cos(Phi_end  ),R_end  *sin(Phi_end  ), z_end  ))
        
        
        #import IPython
        #IPython.embed()
                
        t = linspace(0,1,200)
        X,Y,Z = (End-Start)[:,:,None]*t[None,None,:]+Start[:,:,None]
        R = hypot(X,Y)

        self.los_length = linalg.norm(End-Start)

        #X,Y,Z = array(n)[:,:,None]*t[None,None,:]+c_[R_start,z_start].T[:,:,None]
        
        rho = self.eqm.rz2rho(R[None]-dR,Z[None]-dZ,time,coord_out=self.rho_lbl)[0]


        rho_tg = zeros(n_points)
        xlim = zeros((3,n_points))
        ylim = zeros((3,n_points))
        zlim = zeros((3,n_points))
        


        for i in range(n_points):
            rho_tg[i],tmin = min_fine(t,rho[i])  #not very accurate close to the core :(

            if (len(t[t>tmin]) > 1)&(len(t[t<tmin]) > 1): 
                t1 = interp(1,rho[i,t>tmin]      ,t[t>tmin]      )
                t2 = interp(1,rho[i,t<tmin][::-1],t[t<tmin][::-1])
            else:
                t1,t2 = 0,1
                
            xlim[:,i] = interp([t1,tmin,t2], t, X[i])
            ylim[:,i] = interp([t1,tmin,t2], t, Y[i])
            zlim[:,i] = interp([t1,tmin,t2], t, Z[i])

        #get sign of rho
        n =  (R_end-R_start, z_end-z_start)
        n_perp = c_[-n[1],n[0]]/hypot(n[0],n[1])[:,None]
        ddist = -sum(n_perp*(c_[R_start,z_start]-r_[r0,z0]),1)
        rho_tg*= sign(ddist)
        
        

        
        R_tg,Z_tg = hypot(xlim[1],ylim[1]), zlim[1] 
        theta_tg = arctan2(Z_tg-z0, R_tg-r0)#not very accurate close to the core :(
        
       
        
        
            
        return rho_tg,theta_tg,R_tg,Z_tg
    
        
    def get_q_surfs(self,time,qvalues):
        
        rho = linspace(0,0.98,100)
        q = self.eqm.getQuantity(rho, 'Qpsi', time,coord_in=self.rho_lbl)[0]
        #qvalues = asarray(qvalues)
        
        #print 'q_min = %.3f'%q.min()
        
        ind =  argmax((diff(abs(q))>0))
        #find outermost position of the resonance surface
        rho_surf = interp( qvalues, abs(q)[ind:], rho[ind:],left=nan)
        
        return rho_surf
        
        

        
        
    def signal_info(self,group,name,time):
        return ''
    
    def get_description(self,group,name):
        return 'DIII-D '+str(self.shot)+' diag: '+str(group)+' sig: '+str(name)

        
    def get_plasma_freq(self,rho):
        
        #fetch rotation data from CER
        if not hasattr(self,'plasma_freq'):
            analysis_type = 'CERAUTO'
            system = 'tangential'
            path = 'CER.%s.%s.CHANNEL*'%(analysis_type,system)
            self.MDSconn.openTree('IONS',self.shot)

            nodes = self.MDSconn.get('getnci("'+path+'","path")')
            lengths = self.MDSconn.get('getnci("'+path+':ROT","LENGTH")').data()
            TDI = []
            for node,length in zip(nodes,lengths):
                if length == 0: continue
                node = node.item()
                try:
                    node = node.decode()
                except:
                    pass
                TDI.append(node+':R')
                TDI.append(node+':Z')
                TDI.append(node+':TIME')
                TDI.append(node+':ROTC')
                
            TDI = array(TDI, dtype=object).reshape(-1,4)
            TDI_list = '['+ ','.join(list(TDI.T.flatten())) +']'
            
            R,Z,T,Vtor = self.MDSconn.get(TDI_list).reshape(4,-1)
            T /= 1e3
            Vtor *= 1e3
            self.plasma_freq_rho = self.eqm.rz2rho(R[:,None],Z[:,None],T, coord_out=self.rho_lbl)[:,0]
    
            self.plasma_freq = Vtor/(R*2*pi)
            r0 = interp(T,  self.eqm.t_eq, self.eqm.ssq['Rmag'] )
            self.plasma_freq_rho[R<r0]*= -1
            self.plasma_freq_tvec =  T

        
        from scipy.interpolate import griddata
        tvec = unique(self.plasma_freq_tvec)
        freq=griddata((self.plasma_freq_tvec, self.plasma_freq_rho), self.plasma_freq,(tvec, rho*ones_like(tvec)))



        #plasma_freq_tvec = copy(self.plasma_freq_tvec)
        
        #breaks = where(diff(self.plasma_freq_tvec) > 2* median(diff(self.plasma_freq_tvec)))[0]+1
        #plasma_freq_tvec = insert(plasma_freq_tvec, breaks,  plasma_freq_tvec[breaks])
        #F = insert(F, breaks, nan)
 
        return tvec, np.abs(freq)
        
        
         
        

                
        #import time
        #t = time.time()
        if not  hasattr(self,'plasma_freq'):
            shotfiles = 'CEZ', 'COZ', 'CHZ'
            self.openshotfile = ''
            for s in shotfiles:
                if self.dd.Open(s, self.shot):
                    self.openshotfile = s
                    break
            #not CXRS shotfile avalible
            if self.openshotfile == '':
                self.plasma_freq = None
                return [],[]
            
            #if not self.dd.Open('CEZ', self.shot):
                #if not self.dd.Open('COZ', self.shot):
                    #return [],[]
            #self.openshotfile = 'CEZ'
            names = self.dd.GetNames()
            if 'R_time' in names:
                R = self.dd.GetAreabase('R_time')
                z = self.dd.GetAreabase('z_time') 
            else:
                R = self.dd.GetAreabase('R')[None]
                z = self.dd.GetAreabase('z')[None]

            ind = R.mean(0)!= 0
            R = R[:,ind]
            z = z[:,ind]

            V = self.dd.GetSignalGroup('vrot')[:,ind]
            T = self.dd.GetTimebase('vrot')

            self.plasma_freq = V/(R*2*pi)
            self.plasma_freq_rho = self.eqm.rz2rho(R,z,T, coord_out=self.rho_lbl)
            r0 = interp(T,  self.eqm.t_eq, self.eqm.ssq['Rmag'] )
            self.plasma_freq_rho[R<r0[:,None]]*= -1
            self.plasma_freq_tvec =  T


        if self.plasma_freq is None:
            #not CXRS shotfile avalible
            return [],[]
            
            
        F = zeros_like(self.plasma_freq_tvec)
        for it in range(len(F)):
            sind = argsort(self.plasma_freq_rho[it])
            F[it] = interp(abs(rho),self.plasma_freq_rho[it,sind], self.plasma_freq[it,sind]) 
            #[interp(rho,f,r) for f,r in zip(self.plasma_freq,self.plasma_freq_rho)]
        
        
        #print 'get_plasma_freq',time.time()-t
        plasma_freq_tvec = copy(self.plasma_freq_tvec)
        
        breaks = where(diff(self.plasma_freq_tvec) > 2* median(diff(self.plasma_freq_tvec)))[0]+1
        plasma_freq_tvec = insert(plasma_freq_tvec, breaks,  plasma_freq_tvec[breaks])
        F = insert(F, breaks, nan)
 
        return plasma_freq_tvec,F
        
    
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
                
                self.elm_start = elm_start[ind]/1e3 #[s]
                self.elm_peak = elm_peak[ind]/1e3#[s]
                self.elm_end = elm_end[ind]/1e3#[s]
                self.elm_val = elm_val[ind]
                
                

                #self.MDSconn.openTree('SPECTROSCOPY', self.shot)
                #TDI = '_x=\\SPECTROSCOPY::'+ fname
                #filterscope = self.MDSconn.get(TDI).data()
                #filterscope_t = self.MDSconn.get('dim_of(_x)').data()/1e3
 
            except Exception as e:
                self.elm_start = None
                print('ELM detection issue: ', e)
                return signal
     
        if self.elm_start is None:
            return signal
        
        valid = (self.elm_start > tvec[0])&(self.elm_start < tvec[-1])
        
        #delete also 0.1ms before and after due to uncertainties in elm time
        ind_start = tvec.searchsorted(self.elm_start[valid]-1e-4)
        ind_end = np.minimum(tvec.searchsorted(self.elm_peak[valid]+1e-4), len(tvec)-1)
        
 
        corrected_signal = copy(signal)
        for i1,i2 in zip(ind_start, ind_end):
            corrected_signal[i1:i2] = np.linspace(signal[i1], signal[i2], i2-i1)
            
            
        #plot(tvec, signal)
        #plot(tvec, corrected_signal)
        #plot(filterscope_t,filterscope/1e15 )
        #plot(c_[self.elm_start, self.elm_peak, self.elm_end].flatten(), c_[0*self.elm_val,self.elm_val/1e15,self.elm_val*0].flatten())        
        #show()

        return corrected_signal
            

    #self.dd.Open('COZ', )


class spaced_vector:
    #fast equaly spaced vector object
    def __init__(self,start,stop,step=1):
        self.start = start
        self.stop = stop
        self.step = step
        self.size = int((stop-start)/step)
        self.shape = (self.size,) 
    
    def searchsorted(self,vals):
        vals = asarray(vals)
        vals = minimum(maximum(self.start,vals),self.stop)
        return int_((vals-self.start)/self.step)
    
    def __getitem__(self, key):
        #print 'key',key
        if isinstance(key,slice):
            #print 'get slide'
            return self.__getslice__(key)
        
        key = array(key,copy=True)
        key[key<0] += self.size
        #print ',',key, self.size

        return self.start+key*self.step

    def __array__(self,data=None):
        return arange(self.start, self.stop,self.step)
    
    def __getslice__(self, x1,x2=None):
 
        if isinstance(x1,slice):

            start = 0 if x1.start is None else  x1.start
            stop = self.size if x1.stop is None else x1.stop
            step = 1 if x1.step is None else x1.step
        else:
            start = x1
            stop = x2
            step = 1

        start,stop = self.__getitem__((start,stop))
        step = step*self.step
        
        return spaced_vector(start,stop,step)


from matplotlib.pylab import *


def main():
    
 

    shot = 175860
    
    
    mds_server = "localhost"
    #mds_server = "atlas.gat.com"
    import MDSplus as mds

    MDSconn = mds.Connection(mds_server )
    from map_equ import equ_map
    eqm = equ_map(MDSconn,debug=False)
    eqm.Open(shot,diag='EFIT01' )
    MDSconn2 = mds.Connection(mds_server)

    
    ll = loader(shot,exp='DIII-D',eqm=eqm,rho_lbl='rho_pol',MDSconn=MDSconn2 )
    t,f = ll.get_plasma_freq(0.5)
    plot(t,f)
    show()
    
   
    
if __name__ == "__main__":
    main()
    
