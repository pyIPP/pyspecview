import numpy as np
from .loader import * 
import os, sys
from scipy.interpolate import interp1d
 
 
def check(shot):
    #fastest check if the shotfile exist

    path = shot_path+'/%d/XX/%s/%d'

    for diag in ('RMC', 'RMA','RMB'):
        if os.path.isfile(path%(shot//10,diag, shot)):
            return True
        
    path = shot_path+'/%d/L0/%s/%d'

    if os.path.isfile(path%(shot//10,'RAD', shot)):
        return True

    return False

 
class loader_ECE(loader):
    radial_profile=True
    units = 'eV'

    def __init__(self,*args, **kargs):
        super(loader_ECE,self).__init__(*args, **kargs)

        self.mixers = None
        self.corrupted_calibration=False
        self.names = []
        new_rmd_shotfile = False
   
        if self.dd.Open('RMC',self.shot ,experiment=self.exp, edition=self.ed):
            #new diagnostic
    
            self.shotfile = 'RMC'

            self.mixers = self.dd.GetParameter('parms-A', 'IFGROUP')

            self.CALA1_M0 = self.dd.GetParameter('CAL-A1', 'MULTIA00')
            self.CALA1_M1 = self.dd.GetParameter('CAL-A1', 'MULTIA01')
            self.CALA1_S0 = self.dd.GetParameter('CAL-A1', 'SHIFTB00')
            self.CALA1_S1 = self.dd.GetParameter('CAL-A1', 'SHIFTB01')
            
            self.CALA2_M0 = self.dd.GetParameter('CAL-A2', 'MULTIA00')
            self.CALA2_M1 = self.dd.GetParameter('CAL-A2', 'MULTIA01')
            self.CALA2_S0 = self.dd.GetParameter('CAL-A2', 'SHIFTB00')
            self.CALA2_S1 = self.dd.GetParameter('CAL-A2', 'SHIFTB01')
            self.dd.Close()
                    
            calib_shot = dd.PreviousShot('RMD', self.shot, experiment=self.exp)
            if self.dd.Open('RMD',calib_shot,experiment=self.exp, edition=self.ed):
                prefix = '' if np.any(np.in1d(self.dd.GetNames(),'CAL-A1')) else 'e'
                        
                self.CALA1_M2 = self.dd.GetParameter(prefix+'CAL-A1', 'MULTIA00')
                self.CALA1_S2 = self.dd.GetParameter(prefix+'CAL-A1', 'SHIFTB00')
                self.CALA2_M2 = self.dd.GetParameter(prefix+'CAL-A2', 'MULTIA00')
                self.CALA2_S2 = self.dd.GetParameter(prefix+'CAL-A2', 'SHIFTB00')
                self.names = np.where(self.dd.GetParameter('parms-A','AVAILABL'))[0]+1
                if 'rztime' in self.dd.GetNames() and calib_shot == self.shot:
                    new_rmd_shotfile = True
                    
            
            else:
                print( 'calibratin shotfile was not found!')
                self.CALA1_M2 = np.ones_like(self.CALA1_M0)
                self.CALA1_S2 = np.zeros_like(self.CALA1_M0)
                self.CALA2_M2 = np.ones_like(self.CALA1_M0)
                self.CALA2_S2 = np.zeros_like(self.CALA1_M0)
         
            
        elif self.dd.Open('RMA',self.shot ,experiment=self.exp, edition=self.ed):
            #old diagnostic
            self.shotfile = 'RMA'
            self.CALA1_M0 = self.dd.GetParameter('CALIB', 'MULTIA00')
            self.CALA1_S0 = self.dd.GetParameter('CALIB', 'SHIFTB00')

            self.dd.Open('RMB',self.shot ,experiment=self.exp, edition=self.ed)
            self.CALA2_M0 = self.dd.GetParameter('CALIB', 'MULTIA00')
            self.CALA2_S0 = self.dd.GetParameter('CALIB', 'SHIFTB00')
            
        elif self.dd.Open('RAD',self.shot ,experiment=self.exp, edition=self.ed):
            #old+slow diagnostic
            self.shotfile = 'RAD'
            self.CALA1_M0 = self.dd.GetParameter('CALIB', 'MULTIA00')
            self.CALA1_S0 = self.dd.GetParameter('CALIB', 'SHIFTB00')
            self.CALA2_M0 = self.dd.GetParameter('CALIB', 'MULTIA00')
            self.CALA2_S0 = self.dd.GetParameter('CALIB', 'SHIFTB00')
        else:
            raise Exception('Te shotfile is missing')
        
        self.dd.Close()

        if len(self.names) == 0: 
            self.names = arange(60)+1

        if self.mixers is None: # it was missing in old shotfiles 
            self.mixers = int_(r_[(1,)*36,(2,)*12, (3,)*12])
            
        self.groups = ['mixer %d'%d for d in np.unique(self.mixers)]

        #get the final clalibration factors
        if not new_rmd_shotfile:
            open_sf = self.dd.Open('CEC',self.shot)
        else:
            open_sf = self.dd.Open('RMD',self.shot)
        
        if not open_sf:
            print( 'error: the calibration shotfile is not avalible!!')
            self.calfact = np.ones(len(self.names))
            return 
            
        self.RZtime = self.dd.GetTimebase('rztime')
        self.R = self.dd.GetSignal('R-A')
        self.z = self.dd.GetSignal('z-A')
        if 'Trad-B' in self.dd.GetNames():
            self.R = np.hstack((self.R,self.dd.GetSignal('R-A')))
            self.z = np.hstack((self.z,self.dd.GetSignal('z-A')))

        self.calfact = self.dd.GetParameter('parms-A', 'calfact')
        if self.calfact is None or all(self.calfact == 0):
            #BUG extract calibration directly from CEC!
    
            tvec = self.dd.GetTimebase('Trad-A')
            i0,imin,imax = tvec.searchsorted((0,self.tmin,self.tmax))
            TeCEC = self.dd.GetSignal('Trad-A',nbeg=imin, nend=imax-1).mean(0)

            self.dd.Open('RAD',self.shot)
            tvec = self.dd.GetTimebase('TIME-AD', cal=True)
            i0 = tvec.searchsorted(0)
            Te1 = self.dd.GetSignalCalibrated('SI-MI-A', nend=imax+1)
            Te2 = self.dd.GetSignalCalibrated('SI-MI-B', nend=imax+1)
            Te3 = self.dd.GetSignalCalibrated('SI-MI-C', nend=imax+1)
            Te4 = self.dd.GetSignalCalibrated('SI-MI-D', nend=imax+1)
            
            TeRAD0 = np.hstack((Te1,Te2,Te3,Te4))[:i0].mean(0)
            TeRAD  = np.hstack((Te1,Te2,Te3,Te4))[i0:].mean(0)
            TeRAD -= TeRAD0
            self.calfact = TeCEC/abs(TeRAD)

        self.dd.Close()
        #sometimes missing calibration
        if hasattr(self, 'CALA2_M1') and all(self.CALA2_M2 == 1):
            self.CALA1_M2 = -self.calfact[:30]
            self.CALA2_M2 = -self.calfact[30:]

    def get_signal_groups(self):
        return  self.groups
            
    def get_names(self,group):
        mixer = int(group[-1])

        return self.names[self.mixers[self.names-1]==mixer]

    def get_signal(self,group, name,calib=False, tmin=None,tmax=None):
        #print 'get_signal', tmin,tmax 
        if tmin is None:    tmin = self.tmin
        if tmax is None:    tmax = self.tmax
            
        if np.size(name) > 1:
            data = [self.get_signal(group, n, calib=calib, tmin=tmin,tmax=tmax) for n in name]
            return data

        n = int(name)-1
        sig_name = 'Trad-A1' if n < 30 else 'Trad-A2'
        
        #cache loaded Te data
        if n < 30:
            if not hasattr(self,'TradA1'):
                if self.shotfile == 'RMC':
                    self.dd.Open(  self.shotfile, self.shot, experiment=self.exp, edition=self.ed)
                    self.tvecA1 = self.dd.GetTimebase(sig_name, cal=True)
                    imax = self.tvecA1.searchsorted(self.tmax)+1
                    self.TradA1 = self.dd.GetSignal(sig_name,nend=imax)
                    
                if self.shotfile in ('RMA','RMB','RAD') :
                    if self.shotfile == 'RMA': self.shotfile = 'RMB'
                    self.dd.Open(  self.shotfile, self.shot, experiment=self.exp, edition=self.ed)
                    self.tvecA1 = self.dd.GetTimebase('SI-MI-A', cal=True)
                    imax = self.tvecA1.searchsorted(self.tmax)+1
                    self.TradA1 = self.dd.GetSignal('SI-MI-A',nend=imax)
                    self.TradA1 = np.hstack((self.TradA1, self.dd.GetSignal('SI-MI-B',nend=imax)))
        
                self.dd.Close()
                self.tvecA1 = self.tvecA1[:imax+1]

            tvec,sig = self.tvecA1,self.TradA1[:,n]

        else:
            if not hasattr(self,'TradA2'):
                if self.shotfile == 'RMC':
                    self.dd.Open( self.shotfile, self.shot, experiment=self.exp, edition=self.ed)
                    self.tvecA2 = self.dd.GetTimebase(sig_name, cal=True)
                    imax = self.tvecA2.searchsorted(self.tmax)+1
                    self.TradA2 = self.dd.GetSignal(sig_name,nend=imax)
                    
                if self.shotfile in ('RMA','RMB','RAD'):
                    if self.shotfile == 'RMB': self.shotfile = 'RMA'
                    self.dd.Open(self.shotfile, self.shot, experiment=self.exp, edition=self.ed)
                    self.tvecA2 = self.dd.GetTimebase('SI-MI-C', cal=True)
                    imax = self.tvecA2.searchsorted(self.tmax)+1
                    self.TradA2 = self.dd.GetSignal('SI-MI-C',nend=imax)
                    self.TradA2 = np.hstack((self.TradA2, self.dd.GetSignal('SI-MI-D',nend=imax)))
                    
                self.tvecA2 = self.tvecA2[:imax+1]
                self.dd.Close()
                                
            tvec,sig = self.tvecA2,self.TradA2[:,n-30]


        ioff,imin,imax = tvec.searchsorted((0,tmin,tmax))
        tvec = tvec[imin:imax]
    
        
   
        #data are loaded as uncalibrated because it is much faster!!
        if calib:
            offset = sig[:ioff].mean()
            sig = np.single(sig[imin:imax])

            if n < 30:
                if self.shotfile == 'RMC':
                    #print(( self.CALA2_M0[n],self.CALA2_M1[n],self.CALA2_M2[n]))

                    offset*=  self.CALA1_M0[n]*self.CALA1_M1[n]*self.CALA1_M2[n]
                    offset+= (self.CALA1_S0[n]*self.CALA1_M1[n]+self.CALA1_S1[n])*self.CALA1_M2[n]+self.CALA1_S2[n]
                    
                    sig   *=  self.CALA1_M0[n]*self.CALA1_M1[n]*self.CALA1_M2[n]
                    sig   += (self.CALA1_S0[n]*self.CALA1_M1[n]+self.CALA1_S1[n])*self.CALA1_M2[n]+self.CALA1_S2[n]-offset
                
                else:
                    offset*=  self.CALA1_M0[n%16]
                    offset+= self.CALA1_S0[n%16]
                    sig*=  self.CALA1_M0[n%16]
                    sig+=  self.CALA1_S0[n%16]-offset
                    sig*= sign(mean(sig[::10]))*self.calfact[n]
                
            else:
                if self.shotfile == 'RMC':
                    n-= 30
                    
                    offset*=  self.CALA2_M0[n]*self.CALA2_M1[n]*self.CALA2_M2[n]
                    offset+= (self.CALA2_S0[n]*self.CALA2_M1[n]+self.CALA2_S1[n])*self.CALA2_M2[n]+self.CALA2_S2[n]
                    
                    sig*= self.CALA2_M0[n]*self.CALA2_M1[n]*self.CALA2_M2[n]
                    sig+= (self.CALA2_S0[n]*self.CALA2_M1[n]+self.CALA2_S1[n])*self.CALA2_M2[n]+self.CALA2_S2[n]-offset
                    
                else:
                    offset*=  self.CALA1_M0[n%16]
                    offset+= self.CALA1_S0[n%16]
                    sig*=  self.CALA1_M0[n%16]
                    sig+=  self.CALA1_S0[n%16]-offset
                    sig*= sign(mean(sig[::10]))*self.calfact[n]
        else:
            #change sign 
            if self.shotfile == 'RMC':
                sig = (2**16-1)-sig[imin:imax]
            else:
                sig = sig[imin:imax]
            
        return [tvec,sig]

    def get_Te0(self,tmin,tmax,IDA=True,dR=0,dZ=0):
        
        
        time = (tmin+tmax)/2
        rho = self.get_rho('',self.names,time,dR=dR,dZ=dZ)[0]
        status = False
        if IDA:
            status = False
            if not self.dd.diag is 'IDA':
                status = self.dd.Open('IDA', self.shot)
            else:
                status = True
            if status:
                ida_tvec = self.dd.GetTimebase('time')
                
                if tmin > ida_tvec[-1]:
                    return 
                    
                imin,imax = ida_tvec[:-1].searchsorted((tmin,tmax))

                Te = self.dd.GetSignal('Te')[:,imin-1:imax].mean(-1)
                rho_ida = self.dd.GetAreabase('rhop',imin,imin)
                rho_ida = self.eqm.rho2rho(rho_ida, time)[0]
                TeEce2 = np.interp(abs(rho), rho_ida, Te)  #dR, dZ are not accounted 
   
                self.dd.Close()
                return rho, TeEce2
        
        
        #ugly robust fit of ECE data
        if not self.dd.diag is 'CEC':
            status = self.dd.Open('CEC',self.shot)
        
        if not status:
            return

        tvec = self.dd.GetTimebase('Trad-A', cal=True)
        imin,imax = tvec[:-1].searchsorted((tmin, tmax))

        Te = self.dd.GetSignal('Trad-A',nbeg=imin, nend=imax+1)
        Te = np.nanmedian(Te[:,self.names-1], 0)
        from scipy.signal import medfilt
  

        from scipy.interpolate import LSQUnivariateSpline

        ind = np.argsort(rho)
        ind = ind[np.abs(rho[ind]) < 1 & Te[ind] > 0]
        x = np.r_[rho[ind],1]
        y = np.log(np.r_[medfilt(Te[ind],5),min(Te[ind].min(),100)])
        ind_ = np.argsort(r_[-x,x])
        
        try:
            S   = LSQUnivariateSpline(np.r_[-x,x][ind_],np.r_[y,y][ind_],np.linspace(-0.999,0.999,21),ext=3, bbox=[-1,1],k=2)
            c = (S(x)-y)/S(x) < 0.005
            x,y = x[c], y[c]
            ind_ = np.argsort(r_[-x,x])
            S   = LSQUnivariateSpline(np.r_[-x,x][ind_], np.r_[y,y][ind_],np.linspace(-0.999,0.999,21),ext=3, bbox=[-1,1],k=2)
            Te_ = np.exp(S(rho))
            #if 
            if not all(np.isfinite(Te_)):
                raise Exception()

        except Exception as e:
            print(( 'spline fit has failured: ', e))
            Te_ = np.exp(np.interp(rho, np.r_[-x,x][ind_], np.r_[y,y][ind_]))
        
        return   rho, Te_

    def get_names_phase(self):
        pass
    def get_signal_phase(self,name,calib=False):
        pass
    def get_phi_tor(self,name=None):
        return np.deg2rad(300+40.5) # position of SXR and relative distance between SXR anbd ECE #BUG what abiut these 45deg in SXR loader? 
    def get_phase_corrections(self,name):
        pass
        
    
    def get_RZ_theta(self, time,names,dR=0,dZ=0,cold=False):

        if not hasattr(self,'RZtime'):
            err = np.ones(len(names))*np.nan
            return err,err,err

        time = np.clip(time,*self.RZtime[[0,-1]])
        ch_ind = np.in1d(self.names, np.int_(names))
        
        R = interp1d(self.RZtime, self.R[:,np.array(names)-1],axis=0)(time)
        z = interp1d(self.RZtime, self.z[:,np.array(names)-1],axis=0)(time)

        r0 = np.interp(time, self.eqm.t_eq, self.eqm.ssq['Rmag'])+dR
        z0 = np.interp(time, self.eqm.t_eq, self.eqm.ssq['Zmag'])+dZ
   
        return R,z, np.arctan2(z-z0, R-r0)
        
    
    def get_rho(self,group,names,time,dR=0,dZ=0):

        if hasattr(self,'RZtime'):
            time = np.clip(time, self.tmin,self.tmax)
            R,z,theta = self.get_RZ_theta(time,names,dR=dR,dZ=dZ)
            
            rho = super(loader_ECE,self).get_rho(time,R,z,dR=dR,dZ=dZ)

            r0 = np.interp(time, self.eqm.t_eq, self.eqm.ssq['Rmag'])+dR
            R = np.atleast_1d( R)
            rho[R < r0] *= -1
        else:
            err = np.ones(len(names))*np.nan
            rho,theta,R,z = err,err,err,err
        
        return rho,theta,R,z
    
        
    def signal_info(self,group,name,time):
        name = int(name)
        rho,theta,R,z = self.get_rho(group,[name,],time)

        info = 'ch: '+str(name)+'  R:%.3f m   '%R+self.rho_lbl+': %.3f'%rho
        return info
    
    def get_description(self,group,name):
        return 'AUG %d diag: %s sig: %s'%(self.shot,self.shotfile,name)



