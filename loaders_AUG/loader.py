import numpy as np
import os, logging
import traceback
import aug_sfutils as sf

logger = logging.getLogger('pyspecview.loader')
logger.setLevel(logging.DEBUG)

shot_path = '/afs/ipp-garching.mpg.de/augd/shots/'


def min_fine(x,y): #find extrem with subpixel precision

    i = np.argmin(y)
    if i==0 or i == len(y)-1: return y[i],x[i]    #not working at the edge! 
    A = x[i-1:i+2]**2, x[i-1:i+2], np.ones(3)
    a,b,c = np.linalg.solve(np.array(A).T, y[i-1:i+2])

    ymin =  c-b**2/(4*a)
    xmin = -b/(2*a) 
    return  ymin,xmin


class loader(object):

    pol_mode_num = False
    tor_mode_num = False
    n_mode_range = None
    radial_profile=False
    units = 'a.u.'


    def __init__(self, shot, exp='AUGD', ed=0, eqm=None, rho_lbl='rho_pol', MDSconn=None):
        

        self.eqm = eqm

        self.shot = int(shot)
        self.exp = str(exp)
        self.ed = int(ed)
        self.names = {}
        self.groups = {}
        self.shotfile = 'XXX'
        self.rho_lbl= rho_lbl
        if hasattr(self.eqm, 'time'):
            time = self.eqm.time
            self.tmin = time[0]
            self.tmax = time[-1] + (time[-1] - time[0])/len(time)
        else:
            logger.error('No magnetic equlibrium is avalible!')
            self.tmin = 0
            self.tmax = 10
             
    def get_signal_groups(self):
        return self.groups

    def get_names(self,group):
        pass 
    def get_signal(self,group, name,calib=False):
        pass
    def get_names_phase(self):
        pass
    def get_signal_phase(self,name,calib=False):
        pass
    def get_phi_tor(self,name=None):
        pass
    def get_phase_corrections(self,name):
        pass

            
    def mag_theta_star(self, tshot=4, rho=0, dR=0, dZ=0, n_theta=100):

        rhop = sf.rho2rho(self.eqm, rho, tshot, coord_in=self.rho_lbl, coord_out='rho_pol' )[0]
   
        single_radius = False
        
        if np.size(rhop) == 1:
            rhop = np.r_[rhop,rhop+1e-2]
            single_radius = True

        angle = np.linspace(0, 2*np.pi,n_theta)

        magr,magz = sf.rhoTheta2rz(self.eqm, rhop, angle, tshot)
      
        magr = magr[0].T + dR
        magz = magz[0].T + dZ

        r0 = np.interp(tshot, self.eqm.time, self.eqm.Rmag) + dR
        z0 = np.interp(tshot, self.eqm.time, self.eqm.Zmag) + dZ

        n_rho, n_theta = magr.shape

        drdrho,drtheta = np.gradient(magr)
        dzdrho,dztheta = np.gradient(magz)
        dpsidrho,dpsitheta = np.gradient(np.tile(rhop**2, (n_theta,1)).T )

        grad_rho = np.dstack((drdrho,dzdrho,dpsidrho ))
        grad_theta = np.dstack((drtheta,dztheta,dpsitheta))
        normal = np.cross(grad_rho,grad_theta,axis=-1)
        dpsi_dr = -normal[:,:,0]/(normal[:,:,2]+1e-8) #Bz
        dpsi_dz = -normal[:,:,1]/(normal[:,:,2]+1e-8) #Br
 
        dtheta_star = ((magr-r0)**2+(magz-z0)**2)/(dpsi_dz*(magz-z0)+dpsi_dr*(magr-r0))/magr
    
        if not np.all(np.isfinite(dtheta_star)):
            print(( 'dtheta_star is corrupted', tshot))
            print(( np.sum(np.isfinite(dtheta_star))/float(np.size(dtheta_star))))

        theta = -np.arctan2(magz-z0, -magr+r0 )

        theta = np.unwrap(theta-theta[:,(0,)],axis=1)
        from scipy.integrate import cumtrapz

        #definition of the thetat star by integral
        theta_star = cumtrapz(dtheta_star,theta,axis=1,initial=0)
        correction = (n_theta-1.)/n_theta
        if np.sum(np.abs(magr[:,0]-magr[:,-1]))<1e-5 and np.sum(np.abs(magz[:,0]-magz[:,-1]))<1e-5:
            correction = 1
        theta_star/= theta_star[:,(-1,)]/(2*np.pi)/correction     #normalize to 2piz

        if single_radius:
            return magr[0],magz[0],angle,theta_star.mean(0)
            
        return magr,magz, angle, theta_star

    
    def get_rho(self,time,R_start,z_start,Phi_start=None, R_end=None,z_end=None,Phi_end=None, dR=0,dZ=0):
        #get rho_pol/rho_tor for points, R2, and Z2 must be specified for LOSs
        #BUG phi coordinate is ignored!!! it will not work for toroidal LOS! 

        n_points = np.size(R_start)
        time = np.clip(time, self.tmin, self.tmax)

        if R_end is None or z_end is None:
            return sf.rz2rho(self.eqm, R_start-dR,z_start-dZ,time, coord_out=self.rho_lbl)[0]

        r0 = np.interp(time, self.eqm.time, self.eqm.Rmag)+dR
        z0 = np.interp(time, self.eqm.time, self.eqm.Zmag)+dZ

        Phi_start = np.deg2rad(Phi_start)
        Phi_end   = np.deg2rad(Phi_end)

        Start = np.array((R_start*np.cos(Phi_start),R_start*np.sin(Phi_start), z_start))
        End   = np.array((R_end  *np.cos(Phi_end  ),R_end  *np.sin(Phi_end  ), z_end  ))

        t = np.linspace(0,1,200)
        X,Y,Z = (End-Start)[:,:,None]*t[None,None,:]+Start[:,:,None]
        R = np.hypot(X,Y)
        rho = sf.rz2rho(self.eqm, R[None]-dR,Z[None]-dZ,time,coord_out=self.rho_lbl)[0]
        
        rho_tg = np.zeros(n_points)
        xlim = np.zeros((3,n_points))
        ylim = np.zeros((3,n_points))
        zlim = np.zeros((3,n_points))

        for i in range(n_points):
            rho_tg[i], tmin = min_fine(t,rho[i])  #not very accurate close to the core :(

            if (len(t[t>tmin]) > 1)&(len(t[t<tmin]) > 1): 
                t1 = np.interp(1,rho[i,t>tmin]      , t[t>tmin]      )
                t2 = np.interp(1,rho[i,t<tmin][::-1], t[t<tmin][::-1])
            else:
                t1,t2 = 0,1
                
            xlim[:,i] = np.interp([t1,tmin,t2], t, X[i])
            ylim[:,i] = np.interp([t1,tmin,t2], t, Y[i])
            zlim[:,i] = np.interp([t1,tmin,t2], t, Z[i])

        #get sign of rho
        n =  (R_end-R_start, z_end-z_start)
        n_perp = np.c_[-n[1], n[0]]/np.hypot(n[0], n[1])[:,None]
        ddist = -np.sum(n_perp*(np.c_[R_start, z_start] - np.r_[r0, z0]), 1)
        rho_tg*= np.sign(ddist)
        
        R_tg, Z_tg = np.hypot(xlim[1], ylim[1]), zlim[1] 
        theta_tg = np.arctan2(Z_tg-z0, R_tg-r0)#not very accurate close to the core :(
            
        return rho_tg,theta_tg,R_tg,Z_tg
    
        
    def get_q_surfs(self, time, qvalues):

        rho_surf = []
        for qval in qvalues:
            rho_surf.append(sf.get_q_surf(self.eqm, qvalue=qval, t_in=time))

        return rho_surf


    def signal_info(self,group,name,time):
        return ''
    
    def get_description(self,group,name):
        return 'AUG '+str(self.shot)+' diag: '+str(group)+' sig: '+str(name)

    def get_plasma_freq(self,rho):
        
        if not  hasattr(self,'plasma_freq'):
            shotfiles = 'CEZ', 'COZ', 'CHZ'
            self.openshotfile = ''
            for s in shotfiles:
                sfo = sf.SFREAD(s, self.shot)
                if sfo.status:
                    self.openshotfile = s
                    break
            #not CXRS shotfile avalible
            if self.openshotfile == '':
                self.plasma_freq = None
                return [],[]

            if 'R_time' in sfo.objects:
                R = sfo('R_time').T
                z = sfo('z_time').T
            else:
                R = sfo('R')[None]
                z = sfo('z')[None]

            ind = R.mean(0) != 0
            R = R[:, ind]
            z = z[:, ind]

            V = sfo('vrot')[:, ind]
            T = sfo.gettimebase('vrot')

            self.plasma_freq = V/(R*2*np.pi)
            self.plasma_freq_rho = sf.rz2rho(self.eqm, R, z, T, coord_out=self.rho_lbl)
            r0 = np.interp(T,  self.eqm.time, self.eqm.Rmag )
            self.plasma_freq_rho[R<r0[:, None]] *= -1
            self.plasma_freq_tvec =  T

        if self.plasma_freq is None:
            #not CXRS shotfile avalible
            return [], []
            
        F = np.zeros_like(self.plasma_freq_tvec)
        for it in range(len(F)):
            sind = np.argsort(self.plasma_freq_rho[it])
            F[it] = np.interp(np.abs(rho),self.plasma_freq_rho[it,sind], self.plasma_freq[it,sind]) 

        plasma_freq_tvec = np.copy(self.plasma_freq_tvec)
        
        breaks = np.where(np.diff(self.plasma_freq_tvec) > 2* np.median(np.diff(self.plasma_freq_tvec)))[0]+1
        plasma_freq_tvec = np.insert(plasma_freq_tvec, breaks,  plasma_freq_tvec[breaks])
        F = np.insert(F, breaks, np.nan)
 
        return plasma_freq_tvec,F
