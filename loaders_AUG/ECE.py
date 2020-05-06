from .loader import * 
import os

 
 
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
        #import IPython
        #IPython.embed()
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
                    
            calib_shot = self.dd.cShotNr(self.exp, 'RMD',self.shot)
            if self.dd.Open('RMD',calib_shot,experiment=self.exp, edition=self.ed):
                prefix = '' if any(in1d(self.dd.GetNames(),'CAL-A1')) else 'e'
                        
                self.CALA1_M2 = self.dd.GetParameter(prefix+'CAL-A1', 'MULTIA00')
                self.CALA1_S2 = self.dd.GetParameter(prefix+'CAL-A1', 'SHIFTB00')
                self.CALA2_M2 = self.dd.GetParameter(prefix+'CAL-A2', 'MULTIA00')
                self.CALA2_S2 = self.dd.GetParameter(prefix+'CAL-A2', 'SHIFTB00')
                self.names = where(self.dd.GetParameter('parms-A','AVAILABL'))[0]+1
                if 'rztime' in self.dd.GetNames() and calib_shot == self.shot:
                    new_rmd_shotfile = True
                    
            
            else:
                print( 'calibratin shotfile was not found!')

                #self.dd.Open('CEC',self.shot)
                #calfact = self.dd.GetParameter('parms-A', 'calfact')
                self.CALA1_M2 = ones_like(self.CALA1_M0)
                self.CALA1_S2 = zeros_like(self.CALA1_M0)
                self.CALA2_M2 = ones_like(self.CALA1_M0)
                self.CALA2_S2 = zeros_like(self.CALA1_M0)
                #calib_shot = None
            #import IPython
            #IPython.embed()
            
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
            
            
        self.groups = ['mixer %d'%d for d in unique(self.mixers)]


        #get the final clalibration factors
        if not new_rmd_shotfile:
            open_sf = self.dd.Open('CEC',self.shot)
        else:
            open_sf = self.dd.Open('RMD',self.shot)
        
        
        if not open_sf:
            print( 'error: the calibration shotfile is not avalible!!')
            self.calfact = ones(len(self.names))
            return 
            
            
        self.RZtime = self.dd.GetTimebase('rztime')
        self.R = self.dd.GetSignalGroup('R-A')
        self.z = self.dd.GetSignalGroup('z-A')
        if 'Trad-B' in self.dd.GetNames():
            self.R = hstack((self.R,self.dd.GetSignalGroup('R-A')))
            self.z = hstack((self.z,self.dd.GetSignalGroup('z-A')))

        self.calfact = self.dd.GetParameter('parms-A', 'calfact')
        if self.calfact is None or all(self.calfact == 0):
            #BUG extract calibration directly from CEC!
    
            tvec = self.dd.GetTimebase('Trad-A')
            i0,imin,imax = tvec.searchsorted((0,self.tmin,self.tmax))
            TeCEC = self.dd.GetSignalGroup('Trad-A',nbeg=imin, nend=imax-1).mean(0)

            self.dd.Open('RAD',self.shot)
            tvec = self.dd.GetTimebase('TIME-AD')
            i0 = tvec.searchsorted(0)
            Te1 = self.dd.GetSignalGroupCalibrated('SI-MI-A', nend=imax+1)
            Te2 = self.dd.GetSignalGroupCalibrated('SI-MI-B', nend=imax+1)
            Te3 = self.dd.GetSignalGroupCalibrated('SI-MI-C', nend=imax+1)
            Te4 = self.dd.GetSignalGroupCalibrated('SI-MI-D', nend=imax+1)
            
            TeRAD0 = hstack((Te1,Te2,Te3,Te4))[:i0].mean(0)
            TeRAD  = hstack((Te1,Te2,Te3,Te4))[i0:].mean(0)
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
            
        if size(name) > 1:
            data = [self.get_signal(group, n, calib=calib, tmin=tmin,tmax=tmax) for n in name]
            return data
        
        
        n = int(name)-1
        sig_name = 'Trad-A1' if n < 30 else 'Trad-A2'
     
        
        #cache loaded Te data
        if n < 30:
            if not hasattr(self,'TradA1'):
                if self.shotfile == 'RMC':
                    self.dd.Open(  self.shotfile, self.shot, experiment=self.exp, edition=self.ed)
                    self.tvecA1 = self.dd.GetTimebase(sig_name)
                    imax = self.tvecA1.searchsorted(self.tmax)+1
                    self.TradA1 = self.dd.GetSignalGroup(sig_name,nend=imax)
                    
                if self.shotfile in ('RMA','RMB','RAD') :
                    if self.shotfile == 'RMA': self.shotfile = 'RMB'
                    self.dd.Open(  self.shotfile, self.shot, experiment=self.exp, edition=self.ed)
                    self.tvecA1 = self.dd.GetTimebase('SI-MI-A')
                    imax = self.tvecA1.searchsorted(self.tmax)+1
                    self.TradA1 = self.dd.GetSignalGroup('SI-MI-A',nend=imax)
                    self.TradA1 = hstack((self.TradA1, self.dd.GetSignalGroup('SI-MI-B',nend=imax)))
            
                
                self.dd.Close()
                self.tvecA1 = self.tvecA1[:imax+1]

            tvec,sig = self.tvecA1,self.TradA1[:,n]


        else:
            if not hasattr(self,'TradA2'):
                if self.shotfile == 'RMC':
                    self.dd.Open( self.shotfile, self.shot, experiment=self.exp, edition=self.ed)
                    self.tvecA2 = self.dd.GetTimebase(sig_name)
                    imax = self.tvecA2.searchsorted(self.tmax)+1
                    self.TradA2 = self.dd.GetSignalGroup(sig_name,nend=imax)
                    
                if self.shotfile in ('RMA','RMB','RAD'):
                    if self.shotfile == 'RMB': self.shotfile = 'RMA'
                    self.dd.Open(self.shotfile, self.shot, experiment=self.exp, edition=self.ed)
                    self.tvecA2 = self.dd.GetTimebase('SI-MI-C')
                    imax = self.tvecA2.searchsorted(self.tmax)+1
                    self.TradA2 = self.dd.GetSignalGroup('SI-MI-C',nend=imax)
                    self.TradA2 = hstack((self.TradA2, self.dd.GetSignalGroup('SI-MI-D',nend=imax)))
              
                    
                self.tvecA2 = self.tvecA2[:imax+1]
                self.dd.Close()
                                
            tvec,sig = self.tvecA2,self.TradA2[:,n-30]


        ioff,imin,imax = tvec.searchsorted((0,tmin,tmax))
        tvec = tvec[imin:imax]
    
        
   
        #data are loaded as uncalibrated because it is much faster!!
        if calib:
            offset = sig[:ioff].mean()
            sig = single(sig[imin:imax])


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
                    #print(( self.CALA2_M0[n],self.CALA2_M1[n],self.CALA2_M2[n]))
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
            if not self.dd.diagname is 'IDA':
                status = self.dd.Open('IDA', self.shot)
            else:
                status = True
            if status:
                ida_tvec = self.dd.GetTimebase('time')
                
                if tmin > ida_tvec[-1]:
                    return 
                    
                imin,imax = ida_tvec[:-1].searchsorted((tmin,tmax))

                Te = self.dd.GetSignalGroup('Te')[:,imin-1:imax].mean(-1)
                rho_ida = self.dd.GetAreabase('rhop',imin,imin)
                rho_ida = self.eqm.rho2rho(rho_ida, time)[0]
                TeEce2 = interp(abs(rho), rho_ida, Te)  #dR, dZ are not accounted 
   
                self.dd.Close()
                return rho, TeEce2
        
        
        #ugly robust fit of ECE data
        if not self.dd.diagname is 'CEC':
            status = self.dd.Open('CEC',self.shot)
        
        if not status:
            return
            
        tvec = self.dd.GetTimebase('Trad-A')
        imin,imax = tvec[:-1].searchsorted((tmin, tmax))

        Te = self.dd.GetSignalGroup('Trad-A',nbeg=imin, nend=imax+1)
        Te = nanmedian(Te[:,self.names-1], 0)
        from scipy.signal import medfilt
  

        from scipy.interpolate import LSQUnivariateSpline

        ind = argsort(rho)
        ind = ind[abs(rho[ind])<1&(Te[ind]>0)]
        x = r_[rho[ind],1]
        y = log(r_[medfilt(Te[ind],5),min(Te[ind].min(),100)])
        ind_ = argsort(r_[-x,x])
        
        try:
            S   = LSQUnivariateSpline(r_[-x,x][ind_], r_[y,y][ind_],linspace(-0.999,0.999,21),ext=3, bbox=[-1,1],k=2)
            c = (S(x)-y)/S(x) < 0.005
            x,y = x[c], y[c]
            ind_ = argsort(r_[-x,x])
            S   = LSQUnivariateSpline(r_[-x,x][ind_], r_[y,y][ind_],linspace(-0.999,0.999,21),ext=3, bbox=[-1,1],k=2)
            Te_ = exp(S(rho))
            #if 
            if not all(isfinite(Te_)):
                raise Exception()

        except Exception as e:
            print(( 'spline fit has failured: ', e))
            Te_ = exp(interp(rho, r_[-x,x][ind_], r_[y,y][ind_]))
        
        
        
        return   rho, Te_

 

    def get_names_phase(self):
        pass
    def get_signal_phase(self,name,calib=False):
        pass
    def get_phi_tor(self,name=None):
        return deg2rad(300+40.5) # position of SXR and relative distance between SXR anbd ECE #BUG what abiut these 45deg in SXR loader? 
    def get_phase_corrections(self,name):
        pass
        
    
    def get_RZ_theta(self, time,names,dR=0,dZ=0,cold=False):
        
        
        if not hasattr(self,'RZtime'):
            err = ones(len(names))*nan
            return err,err,err

        
        time = max(min(time,self.RZtime[-1] ), self.RZtime[0] )
        ch_ind = in1d(self.names, int_(names))

              
        try:
            #severin's file!
            assert not cold
            from scipy.io import loadmat 
            ECFM = loadmat('ECFM_%d_ed 2.mat'%self.shot)
            R = ECFM['R_warm'][0,ch_ind]
            z = ECFM['z_warm'][0,ch_ind]
            
            print(( 'warm R,z used!!! at %.3fs'%ECFM['time']))
        except:
            R = interp1d(self.RZtime, self.R[:,array(names)-1],axis=0)(time)
            z = interp1d(self.RZtime, self.z[:,array(names)-1],axis=0)(time)
            

        r0 = interp(time, self.eqm.t_eq, self.eqm.ssq['Rmag'])+dR
        z0 = interp(time, self.eqm.t_eq, self.eqm.ssq['Zmag'])+dZ
   
        return R,z, arctan2(z-z0, R-r0)#,R_, z_
        
    
    def get_rho(self,group,names,time,dR=0,dZ=0):

        if hasattr(self,'RZtime'):
            time = max(min(time, self.tmax), self.tmin)
            R,z,theta = self.get_RZ_theta(time,names,dR=dR,dZ=dZ)
            
            rho = super(loader_ECE,self).get_rho(time,R,z,dR=dR,dZ=dZ)

            r0 = interp(time, self.eqm.t_eq, self.eqm.ssq['Rmag'])+dR
            R = atleast_1d( R)
            rho[R < r0] *= -1
        else:
            err = ones(len(names))*nan
            rho,theta,R,z = err,err,err,err
        
        return rho,theta,R,z
    
        
    def signal_info(self,group,name,time):
        name = int(name)
        rho,theta,R,z = self.get_rho(group,[name,],time)

        info = 'ch: '+str(name)+'  R:%.3f m   '%R+self.rho_lbl+': %.3f'%rho
        return info
    
    def get_description(self,group,name):
        return 'AUG %d diag: %s sig: %s'%(self.shot,self.shotfile,name)

  








from matplotlib.pylab import *


def main():
    import os
    import os,sys
    sys.path.append('/afs/ipp-garching.mpg.de/home/t/todstrci/pyspecviewer/')
    from . import dd   
    dd = dd.shotfile()
    sys.path.append('/afs/ipp/home/t/todstrci/TRANSP/')


    from . import map_equ
    
    
    
    shot = 	34797
    
    #shot = 	30516

    
    T = 2.5

    eqm = map_equ.equ_map(debug=True)
    eqm.Open(shot, diag='EQI')
    eqm.read_ssq()
    ece =  loader_ECE(shot,eqm= eqm,rho_lbl='rho_tor')
    
    #import IPython
    #IPython.embed()
    dd.Open('CEC', shot)
    tvec = dd.GetTimebase('time-A')
    Te = dd.GetSignalGroup('Trad-A')

    
    
    
    data = ece.get_signal('',arange(60)+1,calib=True, tmin=5, tmax =5.1)
    
    #exit()
    #plot([d[1].mean() for d in data])
    #show()
    
    
    import IPython
    IPython.embed()
    
    plot(tvec, Te[:,0])
    plot(data[0], data[1])
    show()
    
    

    
    import IPython
    IPython.embed()

    
    
    TeRMC = array(data)[:,1].mean(-1)
    rho,theta_tg,R,Z = ece.get_rho('',arange(60)+1,T,dR=-0.00,dZ=-0.00)
    
    
    r2,t2 = ece.get_Te0(T-1e-2,T+1e-2,False,dR=-0.0,dZ=-0.0)
    
    ece.dd.Open('CEC',ece.shot)
    ind = ece.dd.GetTimebase('time-A').searchsorted([T-1e-2,T+1e-2])
    TeCEC = ece.dd.GetSignalGroupCalibrated('Trad-A',ind[0],ind[1]).mean(0)
        
    plot(r2,t2,'o')
    plot(rho,TeRMC  ,'x')
    plot(rho,-TeRMC  ,'x')

    plot(rho,TeCEC  ,'+')
    show()
    
    
    import IPython
    IPython.embed()
    
    

    #eqm = map_equ.equ_map(debug=False)
    
    #eqm.Open(27274,diag='EQI')
    #eqm.read_ssq()
    #print eqm.ssq 
    #exit()
    #R = linspace(0,2.5,100)
    #rho = eqm.rz2rho( R,zeros(100),3.23 , coord_out = 'rho_tor')[0]
    #plot(R,rho);show()

    
    shot = 25890
    shot = 25299
    shot = 25447
    shot = 30382
    
    shot = 33120
    shot = 27274
    T = 2.1


    eqm = map_equ.equ_map(debug=True)
    eqm.Open(shot, diag='EQI')
    eqm.read_ssq()
    ece =  loader_ECE(shot,eqm= eqm,rho_lbl='rho_tor')
    
    
    

    
    #calfact = copy(ece.calfact)
    #ece.calfact[:] = 1
    
    data = ece.get_signal('',arange(60)+1,calib=True, tmin=T-1e-2, tmax = T+1e-2)
    TeRMC = array(data)[:,1].mean(-1)
    rho,theta_tg,R,Z = ece.get_rho('',arange(60)+1,T,dR=-0.06,dZ=-0.02)
    
    
    r2,t2 = ece.get_Te0(2,False,dR=-0.06,dZ=-0.02)
    
    ece.dd.Open('CEC',ece.shot)
    ind = ece.dd.GetTimebase('time-A').searchsorted([T-1e-2,T+1e-2])
    TeCEC = ece.dd.GetSignalGroupCalibrated('Trad-A',ind[0],ind[1]).mean(0)
        
    plot(r2,t2,'o')
    plot(rho,TeRMC  ,'x')
    plot(rho,TeCEC  ,'+')
    show()
    
    

    

    ece.dd.Open('CEC',ece.shot)
    R_ = ece.dd.GetSignalGroup('R-A')
    RZtime = ece.dd.GetTimebase('rztime')
    R_ = R_[RZtime.searchsorted(T)]


    plot(R, array(data)[:,1].mean(-1))
    plot(R[ece.names-1], array(data)[ece.names-1,1].mean(-1),'o')
    show()
    
    
    
    
    
    
       
    plot(R,Z,'o-')
    show()
    
    
    plot(rho, theta_tg/pi,'o-')
    show()


    
 
    
    #plot(R, theta_tg/pi,'o-')
    #show()
    
    #plot(R, rho,'o-')
    #show()
    ece.tmin = -2
    #data = ece.get_signal('',30,calib=True)
    #plot(data[40][1])

    

        
    r,t = ece.get_Te0(2,True)
    r2,t2 = ece.get_Te0(2,False)

    data = ece.get_signal('',arange(60)+1,calib=True, tmin=1.99, tmax = 2.01)
    #Te = array([mean(d[1]) for d in data] ) 
    
    #data = ece.get_signal('', arange(60)+1,calib=False, tmin=1.99, tmax = 2.01)
    #Te_ = array([mean(d[1]) for d in data]  )
    
        
    #dd.Open('CEC', shot)
 
    #tvec_ = dd.GetTimebase('Trad-A')
    #i,j = tvec_.searchsorted((1.99,2.01))
    #TeCEC = dd.GetSignalGroup('Trad-A')[i:j].mean(0)

    #calfact = dd.GetParameter('parms-A', 'calfact')
        

        
        
    #dd.Open('RAD', shot)
    
    #tvec_rad = dd.GetTimebase('SI-MI-A')
    #teRAD = dd.GetSignalGroup('SI-MI-A')[:,14]
    
    #dd.Open('CEC', shot)
    
    #tvec_cec = dd.GetTimebase('Trad-A')
    #teCEC = dd.GetSignalGroup('Trad-A')[:,14]
    
    #dd.Open('RMB', shot)

    #tvec_RMB = dd.GetTimebase('SI-MI-A')
    #teRMB = dd.GetSignalGroup('SI-MI-A')[:,14]
    T = 2.25
    data = ece.get_signal('',arange(60)+1,calib=True)
    
            #get a smooth calibration for the stationary profile 
    _, Te0 = ece.get_Te0(T,False)
    for t0,d in zip(Te0,data):
        i,j = d[0].searchsorted((T-0.01,T+0.01))
        m = mean(d[1][i:j])
        if m == 0: continue
        d[1] = d[1]*(t0/m)
        #print t0/m, mean(d[1]), t0

    dd.Open('CEC', shot)
    tvec_ = dd.GetTimebase('Trad-A')
    TeCEC = dd.GetSignalGroup('Trad-A')
 
    from scipy.signal import medfilt
    
    
    def MovingAveradge(sig, n,axis=-1):
        #Fast algorithm for calculation of the moving averadge
        #can failure due to rounding errors!!
        #print axis
        sig = sig.swapaxes(axis, -1)
        n = int(n)


        sig.cumsum(axis=-1,out=sig)
        right_pad  = ones(n//2, dtype=sig.dtype)*sig[...,-1][...,None]
        left_pad = zeros(shape(sig[...,0])+((n+1)//2,), dtype=sig.dtype)
        

        cs = concatenate((sig[...,n//2:],right_pad), axis=-1)
        cs-= concatenate((left_pad,sig[...,:-n//2]), axis=-1)
        cs *=1./n
        edge = 1-arange(n//2+1.)/n
        cs[...,:(n+1)//2] /= edge[-2+n%2::-1]
        cs[...,-(n+1)//2:]/= edge
        cs.swapaxes( -1,axis)
        

        return cs

    #import IPython
    #IPython.embed()

    for i in range(60):
        title(i)
        plot(data[i][0], MovingAveradge(double(data[i][1]),50))
        plot(tvec_,TeCEC[:,i])
        show()
    

    
    Te = array([mean(d[1]) for d in data] ) 

    #plot(tvec_RMB,medfilt((teRMB-mean(teRMB[:tvec_rad.searchsorted(0)]))*2.5,5),'-')
    #plot(data[0][0],data[14][1])
    
    #plot(tvec_rad,(teRAD-mean(teRAD[:tvec_rad.searchsorted(0)]))*10)
    #plot(tvec_cec,teCEC)
    #show()
    
    


    
    
    
    dd.Open('CEC', shot)
    tvec_ = dd.GetTimebase('Trad-A')
    i,j = tvec_.searchsorted((T-0.01,T+0.01))
    TeCEC = dd.GetSignalGroup('Trad-A')[i:j].mean(0)
 
        

    plot(r2,t2,'.',label='get_Te0')
    plot(r,t,'o',label='get_Te0')
    plot(r,Te[ece.names-1],'x',label='get_signal')
    #plot(r,TeCEC[ece.names-1],'s',label='CEC')

    show()

    


    #ece.CALA1_M0*ece.CALA1_M1*ece.CALA1_M2
    

        
    dd.Open('RAD', shot)
 
    tvec_2 = dd.GetTimebase('TIME-AD')
    Te1 = dd.GetSignalGroupCalibrated('SI-MI-A')
    Te2 = dd.GetSignalGroupCalibrated('SI-MI-B')
    Te3 = dd.GetSignalGroupCalibrated('SI-MI-C')
    Te4 = dd.GetSignalGroupCalibrated('SI-MI-D')
    #Te5 = dd.GetSignalGroupCalibrated('SI-MI-E')
    #Te6 = dd.GetSignalGroupCalibrated('SI-MI-F')

    
    k,i,j = tvec_2.searchsorted((0,1.99,2.01))
    
    TeRAD = hstack((Te1,Te2,Te3,Te4))[i:j].mean(0)
    TeRAD0 = hstack((Te1,Te2,Te3,Te4))[:k].mean(0)
    
    TeRAD_ = abs(TeRAD-TeRAD0)*ece.calfact
    
    
    Te1 = dd.GetSignalGroup('SI-MI-A')
    Te2 = dd.GetSignalGroup('SI-MI-B')
    Te3 = dd.GetSignalGroup('SI-MI-C')
    Te4 = dd.GetSignalGroup('SI-MI-D')

    TeRAD__ = hstack((Te1,Te2,Te3,Te4))[i:j].mean(0)
    TeRAD0__ = hstack((Te1,Te2,Te3,Te4))[:k].mean(0)
    TeRAD__-= TeRAD0__
    import IPython
    IPython.embed()
    
    
    #data = ece.get_signal('',ece.names,calib=True)
    #Te_ = array([d[1] for d in data]  )
    
    
    #i = 50
    #plot(data[0][0],Te_[i] )
    ##plot(tvec_2 ,-TeRAD_[:,30]/300 )
    #plot(tvec_ ,TeCEC[:,i] )
    #show()
    
    #plot(Te_.mean(1))
    #plot(TeCEC.mean(0))
    #show()
    
    
    

    


    ind = sign(TeRAD__)
    plot((TeRAD__)*ind)
    plot( array(Te_)-65536)

    show()
    
    
    
    Te = array(Te)
    Te[Te==0] = nan
    plot(TeCEC)
    plot((TeRAD_))
    plot( (Te),'o')
    show()
    
    


    #plot(ece.CALA1_M0*ece.CALA1_M1*ece.CALA1_M2*-1)
    #plot(1/calfact*2**16/2)
    #show()




    #plot(tvec, data );
    #plot(tvec_,TeCEC[:,n])
    #show()
    
    

    
    dd.Open('CEC', shot)
    
    

    
    
    
    
    tvec, data = ece.get_signal('',10,calib=False)

    
    
    tvec = dd.GetTimebase('Trad-A')
    TeCEC = dd.GetSignalGroupCalibrated('Trad-A')
    
    TeCEC = TeCEC[tvec.searchsorted(4)]
    
    
    
    ece =  loader_ECE(shot)
    
    data = ece.get_signal('',arange(60),calib=False, tmin=4,tmax=4.001)
    
    IDAte,IDAte2, rho = ece.get_TeIDA(4)
    
    rho2 = ece.get_rho('',arange(60),4)[0]
    
    ecete = array([mean(te) for t, te in data])
    

    plot(rho2, ecete,'o')
    plot(rho2, TeCEC,'x')

    plot(rho, IDAte,'-')
    
    show()
    
            

    plot(ecete[ece.R[0]!= 0]);plot(IDAte);show()
    
    
    ind = argmin(abs((rho2[None,:]-rho[:,None])),axis=1)
    
    
    plot(rho2[ind],'o-')
    plot(rho,'x-')
    show()
    
    plot(ecete[ind],'o-')
    plot(IDAte,'.-')
    plot(IDAte2,'x-')
    show()
   
    


    

    
    
    dd = dd.shotfile()
    dd.Open('IDA', 33000)
    
    t_ida =  2.
    tvec = dd.GetTimebase('time')
    i_tvec = tvec.searchsorted(t_ida)
    Te = dd.GetSignalGroup('Te')[:,i_tvec]
    rhop = dd.GetAreabase('rhop',i_tvec,i_tvec)
    ece_rhop = dd.GetSignalGroup('ece_rhop')[:,i_tvec]
    x_ece = dd.GetAreabase('x_ece')
    TeEce = dd.GetSignalGroup('ece_mod')[:,i_tvec]
    
    TeEce2 = interp(ece_rhop, rhop, Te)
    
    
    plot(ece_rhop,TeEce,'x')
    plot(ece_rhop,TeEce2,'o')
    show()
    
    

    
    

    plot(ece_rhop, TeEce,'o')
    show()
    
    
    t_ida =  2.
    tvec = dd.GetTimebase('time')
    i_tvec = tvec.searchsorted(t_ida)
    #Te = dd.GetSignalGroup('Te')[:,i_tvec]
    #rhop = dd.GetAreabase('rhop',i_tvec,i_tvec)
    #ece_rhop = dd.GetSignalGroup('ece_rhop')[:,i_tvec]
    x_ece = dd.GetAreabase('x_ece')
    TeEce = dd.GetSignalGroup('ece_mod')[:,i_tvec]
    TeIDA = zeros(60)
    TeIDA[x_ece] = TeEce
    
    
    
    #plot(ece_rhop, TeEce,'o')
    show()
    
    
    
    
    

    
    
    
    
    
    
    
    
    
    
    dd.Close()
    
    ece =  loader_ECE(27068)
    theta, theta_star = ece.mag_theta_star(tshot = 4 ,rhop=linspace(0.01,1,60) )
    
    out = ece.get_theta()
    

    
    G = ece.get_signal_groups()
    S = ece.get_names(G[0])
   

    out = ece.get_signal('',arange(60) ,True)
    #print s.size
    #print t.size
    import IPython
    IPython.embed()
    exit()
    
    ece2 =  loader_ECE(27068)
    


if __name__ == "__main__":
    main()

