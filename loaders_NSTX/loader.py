
from numpy import *
import sys, os, random
from scipy.interpolate import interp1d
import traceback, time

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
    qprofile_profile=False    
    units = 'a.u.'
    #plasma_freq=None

    def __init__(self,shot,exp=None,ed=0,eqm=None,rho_lbl='rho_pol',MDSconn=None):
        
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
        
        #get_rho(self.diag_group, (self.signal, ), 
                                    #np.mean(xlims), dR=self.dR_corr, dZ=self.dZ_corr)[0]
        

        print(time, self.tmax, self.tmin)
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
    
        
    def get_q_surfs(self,time,qvalues, rho_lbl=None):
        
        if rho_lbl is None:
            rho_lbl = self.rho_lbl
        rho = linspace(0,0.98,100)
        q = self.eqm.getQuantity(rho, 'Qpsi', time,coord_in=rho_lbl)
        #qvalues = asarray(qvalues)
        
        #print 'q_min = %.3f'%q.min()
        #embed()
        
        #ind = argmax((diff(abs(q))>0))
        #find outermost position of the resonance surface
        #rho_surf = interp(qvalues, abs(q)[ind:], rho[ind:],left=nan)
        #import IPython
        #IPython.embed()
        rho_surf = zeros((len(time), len(qvalues)))

        for it,t in enumerate(time):
            imin = argmin(q[it])
            rho_surf[it] = interp(qvalues, q[it,imin:], rho[imin:],left=nan)
               
       
        return rho_surf
        
        

        
        
    def signal_info(self,group,name,time):
        return ''
    
    def get_description(self,group,name):
        return 'NSTX '+str(self.shot)+' diag: '+str(group)+' sig: '+str(name)

    def get_plasma_freq_q(self,Qvalues, rho_lbl=None):
        
        
        if rho_lbl is None:
            rho_lbl = self.rho_lbl
        #tt = time.time()
        #tvec = []
        #rvec = []

        #rhop = self.eqm.PSIN.T**.5
        #q = abs(self.eqm.q)
        t_eq = self.eqm.t_eq
                
        #rvec = zeros((len(t_eq), len(Qvalues)))
        ##tvec = zeros((len(t_eq), len(Qvalues)))
        
        #for iq,qval in enumerate(Qvalues):
            #for it,t in enumerate(t_eq):
               #imin = argmin(q[:,it])
               #rvec[it, iq] = interp(qval, q[imin:,it], rhop[imin:,it],left=nan)
               
        rvec = self.get_q_surfs(t_eq, Qvalues, rho_lbl=rho_lbl)
     
        return self.get_plasma_freq(rvec,tvec=t_eq, rho_lbl=rho_lbl)
            
        

        

        

    def get_plasma_freq(self,rho,tvec=None,rho_lbl=None):
        if rho_lbl is None:
            rho_lbl = self.rho_lbl


    
        if not  hasattr(self,'plasma_freq'):
   
            tree = 'ACTIVESPEC'
            chers_edition = 'CT1'
            chers_signals = [  'VTS', 'RS', 'TIME']
            TDI = [f'\\{tree}::TOP.CHERS.ANALYSIS.{chers_edition}:{sig}' for sig in chers_signals]
            self.MDSconn.openTree(tree, self.shot)
            V, R, T = [self.MDSconn.get(tdi).data() for tdi in TDI]
            self.MDSconn.closeTree(tree, self.shot )

 
            #use SI units!!
            V *= 1e3 #m/s
            R /= 100 #m
            z = R*0
    

            self.plasma_freq = V/(R*2*pi)
            self.plasma_freq_rho = self.eqm.rz2rho(R,z,T, coord_out=rho_lbl)
            r0 = interp(T,  self.eqm.t_eq, self.eqm.ssq['Rmag'] )
            self.plasma_freq_rho[R<r0[:,None]]*= -1
            self.plasma_freq_tvec =  T


        if self.plasma_freq is None:
            #not CXRS shotfile avalible
            return [],[]
                    

        
        if tvec is None:
            out_tvec = copy(self.plasma_freq_tvec)
            rho = ones_like(out_tvec)*rho
            tvec = out_tvec
        else:
            #out_tvec = linspace(tvec[0], tvec[-1], 1000)
            out_tvec = unique(hstack((tvec.flatten(), self.plasma_freq_tvec)))

 
        F = zeros((len(out_tvec),rho.shape[1]))*nan
        for it,t in enumerate(out_tvec):
            it_f = argmin(abs(t-self.plasma_freq_tvec))
            it_e = argmin(abs(t-tvec))
            if abs(t-tvec[it_e]) > 0.05 or abs(t-self.plasma_freq_tvec[it_f]) > 0.05:
                continue            
            sind = argsort(self.plasma_freq_rho[it_f])
            F[it] = interp(abs(rho[it_e]),self.plasma_freq_rho[it_f,sind], self.plasma_freq[it_f,sind]) 
        #import IPython
        #IPython.embed()
            #embed()
        out_tvec = tile(out_tvec, rho.shape[1])
        F = F.T.flatten()
        
        breaks = where(abs(diff(out_tvec)) > 2* median(abs(diff(out_tvec))))[0]+1
        out_tvec = insert(out_tvec, breaks,  out_tvec[breaks])
        F = insert(F, breaks, nan)

                 
        return out_tvec,F

#import MDSplus
##mdsserver='skylark.pppl.gov:8501'
#mdsserver='localhost'

#MDSconn = MDSplus.Connection(mdsserver)
#print('connected')
##TT = time()
#shot =  121174
#rho_coord = 'rho_tor'
#print(shot)
##print_line( '  * Fetching EFIT01 data ...')
#from map_equ import equ_map
#eqm = equ_map(MDSconn)
#from matplotlib.pylab import *
#eqm.Open(shot, 'LRDFIT09', exp='NSTX')
##eqm._read_pfm()
#eqm.read_ssq()
#eqm._read_scalars()
#eqm._read_profiles()
#LOAD = loader(shot,eqm=eqm,rho_lbl=rho_coord,MDSconn=MDSconn)
#plot(LOAD.get_plasma_freq_q([1,2,3,4], rho_lbl='rho_tor'))


