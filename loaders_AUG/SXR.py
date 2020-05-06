from .loader import * 
from scipy.interpolate import interp1d
import os
from functools import reduce



def check(shot):
    #fastest check if the shotfile exist
    status = False

    path = shot_path+'/%d/SX/SXA/%d'
    status |= os.path.isfile(path%(shot/10,shot))

    
    return status




#not implemented for the old SXR system 
class loader_SXR(loader):
    tor_mode_num = True
    pol_mode_num = False  #not implemented yet!! 

    mode_range = (-3,4)
    radial_profile=True
    units = 'W/m$^2$'

    def __init__(self,*args, **kargs):
        super(loader_SXR,self).__init__(*args, **kargs)

        self.calib_shot = self.dd.cShotNr('AUGD', 'CSX',self.shot)
        
        self.dd.Open( 'CSX', self.calib_shot)
        names = self.dd.GetNames()  

        
        signals = [b[1:] for b in names if b[0]== 'C']
        
        shotfiles    = [  self.dd.GetParameter('C'+s, 'SX_DIAG' ).tostring().decode('utf-8') for s in signals]
        self.Phi     = {s:self.dd.GetParameter('C'+s, 'Tor_Pos' )+45 for s in signals}
        #BUG missing phi_end for T camera!!
        self.R_start = {s:self.dd.GetParameter('C'+s, 'RPINHOLE') for s in signals}
        self.z_start = {s:self.dd.GetParameter('C'+s, 'ZPINHOLE') for s in signals}
        self.R_end   = {s:self.dd.GetParameter('C'+s, 'REND'    ) for s in signals}
        self.z_end   = {s:self.dd.GetParameter('C'+s, 'ZEND'    ) for s in signals}
        self.status  = {s:self.dd.GetParameter('C'+s, 'ADDRESS' )!=256 for s in signals}
        thickness    = {s:self.dd.GetParameter('C'+s, 'THICKNES') for s in signals}
        filt_mat     = {s:self.dd.GetParameter('C'+s, 'FILT-MAT').item().decode('utf-8') for s in signals}
        different_det= {s:(abs(thickness[s]-75e-6)>1e-5) | (filt_mat[s]!= 'Be') for s in signals}
        self.ADCrange  = {s:self.dd.GetParameter('C'+s, 'ADCrange') for s in signals}
        self.ADCmin = 0

        self.MULTIA   = {}
        self.SHIFTB   = {}

        for k,s in enumerate(signals):
            n = int(self.dd.GetParameter('C'+s, 'NCALSTEP'))
            self.MULTIA[s] = [self.dd.GetParameter('C'+s,'MULTIA%.2d'%i) for i in range(n)]
            self.SHIFTB[s] = [self.dd.GetParameter('C'+s,'SHIFTB%.2d'%i) for i in range(n)]
            if self.SHIFTB[s][-1] == 0 or different_det[s]:  #corrupted signal
                shotfiles[k] = 'OOO'
        

        self.SXR_diods = {}

        for sf in unique(shotfiles):
            self.SXR_diods[sf] = []
            
        for d,s in zip(shotfiles,signals):    
            self.SXR_diods[d].append(s)
            
        try:
            self.SXR_diods.pop('OOO')  #null detector
            #self.SXR_diods.pop('SXO')  #null detector
        except:
            pass
        
        self.all_signals = []
        for signals in list(self.SXR_diods.values()):
            self.all_signals.extend(signals)
        self.groups = unique([i[0] for i in self.all_signals])
        self.openshotfile = ''
        
        #print(self.all_signals)
        #print(self.groups )
        #import IPython
        #IPython.embed()
    
    def get_names(self,group):
        names = [s for s in self.all_signals if s[0] == group]
        names.sort()
        return  names

        
        

            
    def get_signal(self,group, name,calib=False,tmin=None,tmax=None):
        
        if tmin is None:    tmin = self.tmin
        if tmax is None:    tmax = self.tmax
            
        
        if isinstance(name,str):
            names = (name, )
        else:  names = name
        
        output = empty(len(names),dtype=object)
        for shotfile,signals in self.SXR_diods.items():
            
            for j, name in enumerate(names):
                if name in signals :
                    if self.openshotfile != shotfile:
                        self.dd.Open(shotfile,self.shot, experiment=self.exp, edition=self.ed)
                        #self.tvec  = self.dd.GetTimebase(name)
                    
                    self.ed = self.dd.edition
                    self.openshotfile = shotfile
                    
                    info = self.dd.GetInfo('Time')
                    tlen = info.tlen
                    tbeg = self.dd.GetTimebase('Time',1,1)[0]
                    tend = self.dd.GetTimebase('Time',tlen,tlen)[0]
                    #BUG assume equally spaced time vector
                    #print(tmin,tmax,tbeg,tlen)
                    nbeg,nend = int_((r_[tmin,tmax]-tbeg)/(tend-tbeg)*tlen)
                    #tvec = self.dd.GetTimebase('Time',nbeg+1,nend)
                    tvec = linspace(tmin,tmax, nend-nbeg)
    
                    sig = self.dd.GetSignal(name, nbeg= nbeg+1,nend = nend)
                    
                    #remove corrupted points!! SLOW!
                    wrong =  (sig < self.ADCmin) |  (sig > self.ADCrange[name]) 

                    if any(wrong):
                        print(( 'correcting wrong points!', sum(wrong)))
                        sig[wrong]=interp(where(wrong)[0], where(~wrong)[0],sig[~wrong])
                    #print sig.max(), name
                    output[j] = [tvec, sig]
                        
                    #if name == 'I_060':
                        #import IPython
                        #IPython.embed()
                        
                        #from matplotlib.pylab import *
                        #self.ADCrange[name]
                        #plot(tvec, sig)
                        #axhline(self.ADCrange[name])
                        #show()
                        
                    
    
        
        #the calibration could be done also by method dd.GetSignalCalibrated, 
        #but in this case a corrupted points could not be identified
        if calib:
            for j, name in enumerate(names):
                sig = single(output[j][1])
                M = self.MULTIA[name]
                S = self.SHIFTB[name]

                sig*= prod(M)
                sig+= ((S[0]*M[1]+S[1])*M[2]+S[2])*M[3]+S[3]
                output[j][1] = sig
 
        
        if len(names) == 1: return  output[0]
        
        #print 'done',name
        #from matplotlib.pylab import *
        #for x,y in output:
            #title(group)
            #plot(x,y)
        #show()

        return output 
        
        
        
    def get_names_phase(self):
        Fsig = self.get_names('F')
        Gsig = self.get_names('G')
        
        FGnames = []
        for f in Fsig:
            for g in Gsig:
                if f[1:]==g[1:]:  FGnames.append('FG'+f[1:])

        return FGnames
        
    
    def get_signal_phase(self,name,calib=False,tmin=None,tmax=None):
        
        if name[:2]!= 'FG':    raise Exception()
      
        tvec1,sig1 = self.get_signal('F','F'+name[2:],calib=calib,tmin=tmin,tmax=tmax)
        tvec2,sig2 = self.get_signal('G','G'+name[2:],calib=calib,tmin=tmin,tmax=tmax)

        #downsample  to the slower DAS
        def reduce(x,y,x_out):
            r = int(round(float(len(x))/len(x_out)))
            x = mean(x[:(len(x)/r)*r].reshape(len(x)/r, r),1)
            y = mean(y[:(len(y)/r)*r].reshape(len(y)/r, r),1)
            return interp( x_out, x,y)
        

        if len(sig1)!= len(sig2):
            if len(sig1) > len(sig2):
                sig1,tvec1 = reduce(tvec1,sig1,tvec2),tvec2
            else:
                sig2,tvec2 = reduce(tvec2,sig2,tvec1),tvec1

        return tvec1, vstack((sig1,sig2)).T
    
    
    def get_phi_tor(self,name=None):
        if name in self.get_names_phase():
            
            phi1 = self.Phi['F'+name[2:]]
            phi2 = self.Phi['G'+name[2:]]
            
            return deg2rad( r_[phi1,phi2])
        elif name in self.Phi:
            return deg2rad(self.Phi[name])
        
        else: 
            try:
                return deg2rad(median([v for k,v in self.Phi.items()]))
            except:
                print(self.Phi.values())
                raise
            
            


    
    def get_rho(self,group,names,time,dR=0,dZ=0):

        #BUG not working for a tangential camera!!!
        R_start = array([self.R_start[name] for name in names])
        z_start = array([self.z_start[name] for name in names])
        R_end = array([self.R_end[name] for name in names])
        z_end = array([self.z_end[name] for name in names])
        Phi = array([self.Phi[name] for name in names])
        rho_tg,theta_tg,R,Z = super(loader_SXR,self).get_rho(time,R_start,
                                    z_start,Phi,R_end,z_end,Phi,dR=dR, dZ=dZ)

        return rho_tg, theta_tg,R,Z
    
    def signal_info(self,group,name,time):
        
        if name[:2] == 'FG': 
            rho_tg1 = self.get_rho(group,[ 'F'+name[2:],],time)[0]
            rho_tg2 = self.get_rho(group,[ 'G'+name[2:],],time)[0]

            phi1 = self.Phi['F'+name[2:]]
            phi2 = self.Phi['G'+name[2:]]
            
            info = str(name)+'   Phi: %d and %d deg, rho_tg: %.2f'%(phi1,phi2,rho_tg1)

        else:
            rho_tg = self.get_rho(group,[ name,],time)[0]
            phi = self.Phi[name]

            info = str(name)+' Phi: %d deg, rho_tg: %.2f'%(phi,rho_tg)
            
            
            
            
        return info
    
    
 

 

def main():
    import os
    import os,sys
    sys.path.append('/afs/ipp/home/t/todstrci/TRANSP/')
    from . import dd   
    #from matplotlib.pylab import *
    from . import map_equ

    shot = 33871


    eqm = map_equ.equ_map(debug=True)
    eqm.Open(shot, diag='EQI')
    eqm.read_ssq()
    sxr =  loader_SXR(shot,eqm= eqm,rho_lbl='rho_tor')
    ##sxr =  loader_SXR(33871)
    #G = sxr.get_signal_groups()
    #S = sxr.get_names('I')
    
    t,x = sxr.get_signal_phase('FG_020')
    
    import IPython
    IPython.embed()
    
    
    #dd = dd.shotfile()

    #dd.Open('SXG',33871 )
    
    #sig2 = dd.GetSignalCalibrated('I_045')
    
    #tvec, sig = sxr.get_signal('J','I_045',calib=True)
    
    #plot(sig2,'-')
    #plot(sig,'-')

    #show()

    
    
    
    
    
    #dd = dd.shotfile()

    
            
 

if __name__ == "__main__":
    main()
