from .loader import * 
import os



def check(shot):
    #fastest check if the shotfile exist

    status = False
    path = shot_path+'/%d/XX/%s/%d'
    for diag in ('XVR', 'XVS','XVT'):
        status |= os.path.isfile(path%(shot/10,diag,shot))
    #print status
    return status


class loader_diod_bolometers(loader):
    tor_mode_num = False  #it was useless
    n_mode_range = (1,2)
    radial_profile=True
    units = 'W/m$^2$'

    def __init__(self,*args, **kargs):
        super(loader_diod_bolometers,self).__init__(*args, **kargs)
        

        bolo_shotfiles = 'XVR', 'XVS'#, 'XVT'
       
        calib_shot = self.dd.cShotNr('AUGD', 'BLC',self.shot)
       
        if not self.dd.Open( 'BLC', calib_shot):
            #second try
            calib_shot = self.dd.cShotNr('AUGD', 'BLC', calib_shot-1)
            if not self.dd.Open( 'BLC', calib_shot):
                raise Exception('No BLC shotfile!')
       
        names = [n.strip() for n in self.dd.GetNames()]

        cams = [b for b in names if b[0] == 'D' and len(b) == 3]
       
        active_los = {s:bool_(self.dd.GetParameter(s, 'active' )) for s in cams}
        activ_cam = [c for c in cams if any(active_los[c])]
        self.Phi_start= {s:self.dd.GetParameter(s, 'P_Blende' )[active_los[s]] for s in activ_cam}
        self.R_start = {s:self.dd.GetParameter(s, 'R_Blende' )[active_los[s]] for s in activ_cam}
        self.z_start = {s:self.dd.GetParameter(s, 'z_Blende' )[active_los[s]] for s in activ_cam}
        self.R_end   = {s:self.dd.GetParameter(s, 'R_end' )[active_los[s]] for s in activ_cam}
        self.z_end   = {s:self.dd.GetParameter(s, 'z_end' )[active_los[s]] for s in activ_cam}
        self.Phi_end = {s:self.dd.GetParameter(s, 'P_end' )[active_los[s]] for s in activ_cam}
      
        self.theta   = {s:arctan2(self.z_end[s]-self.z_start[s],self.R_end[s]-self.R_start[s])  for s in activ_cam}
        self.delta   = {s:self.dd.GetParameter(s, 'delta' )[active_los[s]] for s in activ_cam}
        
        self.all_names = [n for n in self.dd.GetNames() if n[:2] == 'CS']
        self.calibration= {}
        
        for name in self.all_names:
            ncalib = int(self.dd.GetParameter(name, 'NCALSTEP' ))
            multi,shift = [],[]

            for i in range(ncalib):
                multi.append(self.dd.GetParameter(name, 'MULTIA%.2d'%i))
                shift.append(self.dd.GetParameter(name, 'SHIFTB%.2d'%i))
            self.calibration[name[1:]] = multi,shift


        sig_dict = {s:array([c.decode('utf-8') for c in self.dd.GetParameter(s, 'RAW' )]) for s in activ_cam}
        
        self.subcam_ind = {c:[self.delta[c]+self.R_start[c] == r for r in unique(self.delta[c]+self.R_start[c])] for c in activ_cam}
        

        #geometric corrections  estimated by hand!
        n =  self.R_end['DHC']-self.R_start['DHC'], self.z_end['DHC']-self.z_start['DHC']
        corr = 0.04
        self.R_end['DHC'][self.subcam_ind['DHC'][1]]+= -n[1][self.subcam_ind['DHC'][1]]*tan(corr)
        self.z_end['DHC'][self.subcam_ind['DHC'][1]]+= +n[0][self.subcam_ind['DHC'][1]]*tan(corr)
            
        n =  self.R_end['DVC']-self.R_start['DVC'], self.z_end['DVC']-self.z_start['DVC']
        corr = 0.02
        self.R_end['DVC'][self.subcam_ind['DVC'][0]]+= -n[1][self.subcam_ind['DVC'][0]]*tan(corr)
        self.z_end['DVC'][self.subcam_ind['DVC'][0]]+= +n[0][self.subcam_ind['DVC'][0]]*tan(corr)
        
            
        self.groups = activ_cam
       
        xv_names = []
        bolo_los = {}
        for bs in bolo_shotfiles:
            bolo_los[bs] = []
            self.dd.Open(bs, self.shot)
            names = self.dd.GetNames()
            #print( bs, names)
            xv_names.append([n.strip() for n in names])
            #self.dd.Close()
           
        self.signals = {}
        for c in activ_cam:
            sig = sig_dict[c]
            activ = active_los[c]
            ch = 1+arange(len(activ))
            sf = []
                       
            for s in sig[activ]:
                s = s.strip()
                if s in xv_names[0]:
                    sf.append(bolo_shotfiles[0])
                elif s in xv_names[1]:
                    sf.append(bolo_shotfiles[1])
                else:
                    print('missing signal %s'%s)
                    #raise Exception('missing camera %s'%s)
            #print(c, sig[activ],ch[activ],array(sf))
            if len(sf) == 0:
                print('Missing data for camera '+c)
                self.groups.pop(self.groups.index(c))
                continue
            self.signals[c] = list(zip(sig[activ],ch[activ],array(sf)))
           



    def get_names(self,group):
        return sorted([int(ch) for sig,ch,sf in  self.signals[group]])

    def get_signal(self,group, names,calib=False,tmin=None,tmax=None):
        if isinstance(names,str) or not  hasattr(names, '__iter__'):
            names = (names,)
        
        sf_open = None
        data = []
        
        if tmin is None:    tmin = self.tmin
        if tmax is None:    tmax = self.tmax
            
            
        for name in names:
                
            for sig,ch,sf in self.signals[group]:
                if ch == int(name):
                    break
            
            if sf_open!= sf:
                self.dd.Open(sf,self.shot, experiment=self.exp, edition=self.ed)
                sf_open = sf
                tvec = self.dd.GetTimebase('Dio-Time')
                offset_start = tvec.searchsorted(-.1)
                offset_end = tvec.searchsorted(0)
                ind = slice(tvec.searchsorted(tmin), tvec.searchsorted(tmax))
                tvec = tvec[ind]


            offset = self.dd.GetSignal(sig, nbeg=offset_start+1,nend=offset_end).mean().astype('int16')
            signal = self.dd.GetSignal(sig,nbeg=ind.start+1,nend=ind.stop)
            signal-= offset
            
            #the calibration was not sometimes avalible, backup option:
            if  calib:
                signal = float_(signal)

                multi, shift = self.calibration[sig.strip()]
                signal*= prod(multi)
                offset = shift[0]
                for m,s in zip(multi[1:],shift[1:]):
                    offset = offset*m+s
                signal+= offset
            
            data.append(signal)
            
        calib = True #BUG 

        #self.dd.Close()
        #corrections of the diod bolometers estimated by hand!
        if calib:
            if self.shot > 34000:
                pass
            elif self.shot > 31591:
                pass
            elif self.shot > 30135:
                pass
            elif self.shot > 28523:
                pass      
            elif self.shot > 27352:
                if group == 'DHC':
                    ind = where(in1d(names, where(self.subcam_ind['DHC'][1])[0]))[0]
                    for i in ind:  data[i]/= 1.05
                    ind = where(in1d(names, where(self.subcam_ind['DHC'][2])[0]))[0]
                    for i in ind:  data[i]*= 1.3
                if group == 'DVC':
                    ind = where(in1d(names, where(self.subcam_ind['DVC'][0])[0]))[0]
                    for i in ind:  data[i]*= 1.1
       
        if self.shot > 34000:
            if group == 'DHC' and 29 in names:
                ind = where(array(names) == 29)[0][0]
                print(( ind, data[0].shape, len(data)))
                data[ind][:] *= 0  #corrupted LOS
                #data[ind][:] += 1  #corrupted LOS

    
        if len(data) == 1: 
            return tvec, data[0]
        
        
        
        return [(tvec, d) for d in data]
    
        
    def get_names_phase(self):

        
        D = []
        for r1,z1,r2,z2 in zip(self.R_start['DVC'],self.z_start['DVC'],self.R_end['DVC'],self.z_end['DVC']):
            X1 = r1-(r2-r1)/(z2-z1)*z1
            D.append([])
            for r1,z1,r2,z2 in zip(self.R_start['D13'],self.z_start['D13'],self.R_end['D13'],self.z_end['D13']):
                X2 = r1-(r2-r1)/(z2-z1)*z1
                D[-1].append(abs(X1-X2))  #distance of cross-setions with midplane 

        indDVC,indD13  = arange(len(self.signals['DVC'])),argmin(D,1)
        
        
        return ['DVC-%.2d:D13-%.2d'%(i,j) for i,j in zip(indDVC,indD13)]
       
    def get_signal_phase(self,name,calib=False,tmin=None,tmax=None):
        i = int(name[4:6])+1
        j = int(name[-2:])+1

        tvec , sig1 = self.get_signal('DVC', i,calib=calib,tmin=tmin,tmax=tmax)
        tvec2, sig2 = self.get_signal('D13', j,calib=calib,tmin=tmin,tmax=tmax)
        min_len = min(len(sig1), len(sig2))

        return tvec[:min_len], vstack((sig1[:min_len],sig2[:min_len])).T
    

    
    
    def get_phi_tor(self,name):
        i = int(name[4:6])
        j = int(name[-2:])
        
        return r_[self.Phi_['DVC'][i],self.Phi_start['D13'][j]]/180*pi

    
    def get_rho(self,group,names,time,dR=0,dZ=0):
        
        rho_tg,theta_tg, R,Z = super(loader_diod_bolometers,self).get_rho(time, 
                            self.R_start[group],self.z_start[group],self.Phi_start[group],
                            self.R_end[group],self.z_end[group],self.Phi_end[group],dR=dR,dZ=dZ)
        
        
        all_names = self.get_names(group)
       
        ind = where(in1d(  all_names, [int(n) for n in names]))

        return rho_tg[ind],theta_tg[ind],R[ind],Z[ind]

    
    def signal_info(self,group,name,time):
        name = int(name)
        rho_tg = self.get_rho(group,[name,],time)[0]
        
        phi = self.Phi_start[group][int(name)] 
        info = group+' '+str(name)+' Phi: %d deg, rho_tg: %.3f'%(phi,rho_tg)
        return info
    



def main():
    import os
    import os,sys
    sys.path.append('/afs/ipp/home/t/todstrci/TRANSP/')
    #from . import dd  
    import dd
    import matplotlib.pylab  as plt
    #from matplotlib.pylab import *
    shot = 37289
    import map_equ

    eqm = map_equ.equ_map(debug=True)
    eqm.Open(shot, diag='EQI')
    eqm.read_ssq()
    #ece =  loader_ECE(shot,eqm= eqm,rho_lbl='rho_tor')
    
    
    bolo =  loader_diod_bolometers(shot,eqm= eqm,rho_lbl='rho_tor')
    
    names = bolo.get_names('DT1')

    rho_tg,theta_tg,R_tg,Z_tg = bolo.get_rho('DT1',names,4.5)
    
    plt.plot(rho_tg.T)
    plt.show()
    
    
    tvec, sig = bolo.get_signal('DVC', names[0],calib=False)
    plot(tvec, sig)
    tvec, sig = bolo.get_signal('DVC', names[0],calib=True)
    plot(tvec, sig)
    show()

    
    
    
    #plot(R,Z,'o')
    #plot(R2,Z2,'x')
    #show()
    
    
if __name__ == "__main__":
    main()
