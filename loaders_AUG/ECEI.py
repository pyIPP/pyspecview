from scipy.interpolate import interp1d
from .loader import * 
import os
import aug_sfutils as sf

def check(shot):
    #fastest check if the shotfile exist

    status = False

    path = shot_path+'%d/L1/ECI/%d.1' 
    status |= os.path.isfile(path%(shot//10,shot)) and os.stat(path%(shot//10,shot)).st_size

    path = shot_path+'%d/XX/TDI/%d' 
    status |= os.path.isfile(path%(shot//10,shot))

    return status


class loader_ECEI(loader):

    radial_profile=False

    def __init__(self,*args, **kargs):

        super(loader_ECEI,self).__init__(*args, **kargs)

        self.groups = []
        
        path = shot_path+'%d/L1/ECI/%d.1' 
        if os.path.isfile(path%(self.shot/10,self.shot)): 
            self.groups.append('ECI')

        path = shot_path+'%d/XX/TDI/%d' 
        if os.path.isfile(path%(self.shot/10,self.shot)): 
            self.groups.append('ECO')
            self.groups.append('ECN')

        rzo = sf.SFREAD('RZO', self.shot, experiment="ECEI", edition=0)
        rzn = sf.SFREAD('RZN', self.shot, experiment="ECEI", edition=0)
        if rzo.status:
            self.RZtime_old = rzo('time')
            self.R_old = rzo('R')
            self.z_old = rzo('z')  
            self.rho_old = rzo('rho')  

        if rzn.status:
            self.RZtime_new = rzn('time')
            self.R_new = rzn('R')
            self.z_new = rzn('z')  
            self.rho_new = rzn('rho')

        eci_names = ['%.2d:%.2d'%(i+1,j+1) for i in range(16) for j in range(8)]
        eco_names = ['%.2d:%.2d'%(i+1,j+1) for i in range(16) for j in range(8)]
        ecn_names = ['%.2d:%.2d'%(i+1,j+1) for i in range(20) for j in range(8)]
        self.names = {'ECI':eci_names, 'ECO':eco_names, 'ECN':ecn_names  }
        
    def get_signal_groups(self):
        return  self.groups
            
    def get_names(self,group):
        return self.names[group]

    def get_signal(self,group, name,calib=False, tmin=None,tmax=None):
                        
        if tmin is None:    tmin = self.tmin
        if tmax is None:    tmax = self.tmax
        
        #cache loaded Te data
        if group == 'ECI':
            LOS, ind = np.int_(name.split(':'))
            LOS = "LOS" + str(LOS)
            if not hasattr(self, 'ECI_'+LOS):
                eci = sf.SFREAD('ECI', self.shot, experiment="AUGD")
                self.ECI_time = eci("time")
                data = eci(LOS)
                setattr(self, 'ECI_'+LOS, data)

            imin,imax = self.ECI_time.searchsorted([tmin,tmax])
            sig = getattr(self, 'ECI_'+LOS)[imin:imax+1,ind-1]
            tvec = self.ECI_time[imin:imax+1]
            
        elif group in ('ECO','ECN'):
            los,ind = np.int_(name.split(':'))
            ch = (los-1)*(ind-1)
            
            if group == 'ECN': ch += 128

            sig = ch//72+1
            ind = ch%72

            if not hasattr(self, 'TDI_%d'%sig):
                tdi = sf.SFREAD('TDI', self.shot, experiment="AUGD")
                self.TDI_time = tdi.gettimebase("Sig1")
                data = tdi("Sig%d" %sig)
                setattr(self, 'TDI_%d'%sig, data)

            imin,imax = self.TDI_time.searchsorted([tmin,tmax])
            sig = getattr(self, 'TDI_%d'%sig)[imin:imax+1,ind]
            tvec = self.TDI_time[imin:imax+1]

        else:
            raise Exception('ECEI group %s do not exists!'%group)
  
        return tvec,sig 

    def get_names_phase(self):
        pass

    def get_signal_phase(self,name,calib=False):
        pass
    def get_phi_tor(self,name):
        pass
    def get_phase_corrections(self,name):
        pass
    
    def get_rho(self,group,names,time,dR=0,dZ=0):

        if group == 'ECN' and hasattr(self, 'RZtime_new'):
            RZtime = self.RZtime_new
            R = self.R_new
            z = self.z_new
            rho = self.rho_new
        elif hasattr(self, 'RZtime_old'):
            RZtime = self.RZtime_old
            R = self.R_old
            z = self.z_old
            rho = self.rho_old
            
        else:
            print( 'Warning: shotfile with resonance  positions was not found ')
            return np.nan,np.nan,np.nan,np.nan
            
        time = np.clip(time, *RZtime[[0,-1]])
        
        i,j =  np.int_(names[0].split(':'))-1
        R = interp1d(RZtime, R[:,i,j],axis=0)(time)
        z = interp1d(RZtime, z[:,i,j],axis=0)(time)
        rho_p = interp1d(RZtime, rho[:,i,j],axis=0)(time) #dR,dZ are ignored!

        r0 = np.interp(time, self.eqm.time, self.eqm.ssq['Rmag'])+dR
        z0 = np.interp(time, self.eqm.time, self.eqm.ssq['Zmag'])+dZ

        return rho_p, np.arctan2(z-z0, R-r0), R,z, #,R_, z_
        

    def signal_info(self,group,name,time):

        rho, theta, R,z = self.get_rho(group,(name,),time)

        info =  'ch: '+str(name)+'  R: %.3fm  z: %.3fm rho_p: %.2f'%(R,z,rho)
  
        return info
    
    def get_description(self,group,name):
        diag = 'ECI' if group=='ECI' else 'TDI'
        return 'AUG %d diag: %s sig: %s'%(self.shot,diag,name)
