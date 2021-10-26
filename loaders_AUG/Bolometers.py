from .loader import * 
import os, sys
import aug_sfutils as sf


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

        calib_shot = sf.previousshot('BLC', self.shot, exp='AUGD')
        if calib_shot < 0:
            raise Exception('No BLC shotfile!')
        blc = sf.SFREAD('BLC', calib_shot, exp='AUGD')
        names = [n.strip() for n in blc.objects + blc.parsets]

        cams = [b for b in names if b[0] == 'D' and len(b) == 3]
       
        active_los = {s: np.bool_( blc(s)['active'] ) for s in cams}
        activ_cam = [c for c in cams if any(active_los[c])]
        self.Phi_start = {s: blc(s)['P_Blende'][active_los[s]] for s in activ_cam}
        self.R_start   = {s: blc(s)['R_Blende'][active_los[s]] for s in activ_cam}
        self.z_start   = {s: blc(s)['z_Blende'][active_los[s]] for s in activ_cam}
        self.R_end     = {s: blc(s)['R_end'   ][active_los[s]] for s in activ_cam}
        self.z_end     = {s: blc(s)['z_end'   ][active_los[s]] for s in activ_cam}
        self.Phi_end   = {s: blc(s)['P_end'   ][active_los[s]] for s in activ_cam}

        self.theta   = {s: np.arctan2(self.z_end[s]-self.z_start[s],self.R_end[s]-self.R_start[s])  for s in activ_cam}
        self.delta   = {s: blc(s)['delta'][active_los[s]] for s in activ_cam}
        
        self.all_names = [n for n in blc.objects + blc.parsets if n[:2] == 'CS']
        self.calibration= {}

        for name in self.all_names:
            ncalib = int(blc(name)['NCALSTEP'])
            multi,shift = [],[]

            for i in range(ncalib):
                multi.append(blc(name)['MULTIA%.2d' %i])
                shift.append(blc(name)['SHIFTB%.2d' %i])
            self.calibration[name[1:]] = multi,shift

        sig_dict = {s:np.array([ c for c in blc(s)['RAW'] ]) for s in activ_cam}

        self.subcam_ind = {c:[self.delta[c]+self.R_start[c] == r for r in np.unique(self.delta[c]+self.R_start[c])] for c in activ_cam}
        

        #geometric corrections  estimated by hand!
        n =  self.R_end['DHC']-self.R_start['DHC'], self.z_end['DHC']-self.z_start['DHC']
        corr = 0.04
        self.R_end['DHC'][self.subcam_ind['DHC'][1]]+= -n[1][self.subcam_ind['DHC'][1]]*np.tan(corr)
        self.z_end['DHC'][self.subcam_ind['DHC'][1]]+= +n[0][self.subcam_ind['DHC'][1]]*np.tan(corr)
            
        n =  self.R_end['DVC']-self.R_start['DVC'], self.z_end['DVC']-self.z_start['DVC']
        corr = 0.02
        self.R_end['DVC'][self.subcam_ind['DVC'][0]]+= -n[1][self.subcam_ind['DVC'][0]]*np.tan(corr)
        self.z_end['DVC'][self.subcam_ind['DVC'][0]]+= +n[0][self.subcam_ind['DVC'][0]]*np.tan(corr)
 
        self.groups = activ_cam
       
        xv_names = []
        bolo_los = {}
        for bs in bolo_shotfiles:
            bolo_los[bs] = []
            sfo = sf.SFREAD(bs, self.shot)
            names = sfo.objects + sfo.parsets
            xv_names.append([n.strip() for n in names])

        self.signals = {}
        for c in activ_cam:
            sig = sig_dict[c]
            activ = active_los[c]
            ch = 1 + np.arange(len(activ))
            sfiles = []
                       
            for s in sig[activ]:
                s = s.strip()
               # print(s)
                if s == '': continue
                if s in xv_names[0]:
                    sfiles.append(bolo_shotfiles[0])
                elif s in xv_names[1]:
                    sfiles.append(bolo_shotfiles[1])
                else:
                    print('missing signal %s'%s)

            if len(sfiles) == 0:
                print('Missing data for camera '+c)
                self.groups.pop(self.groups.index(c))
                continue
            self.signals[c] = list(zip(sig[activ], ch[activ], np.array(sfiles)))
           

    def get_names(self,group):
        if group in self.signals:
            return sorted([int(ch) for sig, ch, sfile in self.signals[group]])
        else:
            print('Bolometry: no activa data for shotfile %s' %group)

    def get_signal(self,group, names, calib=False, tmin=None, tmax=None):
        if isinstance(names,str) or not  hasattr(names, '__iter__'):
            names = (names,)
        
        sf_open = None
        data = []
        
        if tmin is None:    tmin = self.tmin
        if tmax is None:    tmax = self.tmax
     
        for name in names:
                
            for sig,ch,sfile in self.signals[group]:
                if ch == int(name):
                    break
            
            if sf_open != sfile:
                sfo = sf.SFREAD(sfile, self.shot, experiment=self.exp, edition=self.ed)
                sf_open = sfile
                tvec = sfo('Dio-Time')
                offset_start = tvec.searchsorted(-.1)
                offset_end = tvec.searchsorted(0)
                ind = slice(tvec.searchsorted(tmin), tvec.searchsorted(tmax))
                tvec = tvec[ind]

            offset = sfo.getobject(sig, nbeg=offset_start+1, nend=offset_end).mean().astype('int16')
            signal = sfo.getobject(sig, nbeg=ind.start+1, nend=ind.stop)
            signal-= offset
            
            #the calibration was not sometimes avalible, backup option:
            if  calib:
                signal = np.float_(signal)

                multi, shift = self.calibration[sig.strip()]
                signal*= np.prod(multi)
                offset = shift[0]
                for m,s in zip(multi[1:],shift[1:]):
                    offset = offset*m+s
                signal+= offset
            
            data.append(signal)
            
        calib = True #BUG 

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
                    ind = np.where(np.in1d(names, np.where(self.subcam_ind['DHC'][1])[0]))[0]
                    for i in ind:  data[i]/= 1.05
                    ind = np.where(np.in1d(names, np.where(self.subcam_ind['DHC'][2])[0]))[0]
                    for i in ind:  data[i]*= 1.3
                if group == 'DVC':
                    ind = np.where(np.in1d(names, np.where(self.subcam_ind['DVC'][0])[0]))[0]
                    for i in ind:  data[i]*= 1.1
       
        if self.shot > 34000:
            if group == 'DHC' and 29 in names:
                ind = np.where(np.array(names) == 29)[0][0]
                print(( ind, data[0].shape, len(data)))
                data[ind][:] *= 0  #corrupted LOS
    
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

        indDVC,indD13  = np.arange(len(self.signals['DVC'])),np.argmin(D,1)

        return ['DVC-%.2d:D13-%.2d'%(i,j) for i,j in zip(indDVC,indD13)]
       
    def get_signal_phase(self,name,calib=False,tmin=None,tmax=None):
        i = int(name[4:6])+1
        j = int(name[-2:])+1

        tvec , sig1 = self.get_signal('DVC', i,calib=calib,tmin=tmin,tmax=tmax)
        tvec2, sig2 = self.get_signal('D13', j,calib=calib,tmin=tmin,tmax=tmax)
        min_len = min(len(sig1), len(sig2))

        return tvec[:min_len], np.vstack((sig1[:min_len],sig2[:min_len])).T

    def get_phi_tor(self,name):
        i = int(name[4:6])
        j = int(name[-2:])

        return np.r_[self.Phi_['DVC'][i],self.Phi_start['D13'][j]]/180*np.pi
    
    def get_rho(self,group,names,time,dR=0,dZ=0):
        
        rho_tg,theta_tg, R,Z = super(loader_diod_bolometers,self).get_rho(time, 
                            self.R_start[group],self.z_start[group],self.Phi_start[group],
                            self.R_end[group],self.z_end[group],self.Phi_end[group],dR=dR,dZ=dZ)
        
        all_names = self.get_names(group)
       
        ind = np.where(np.in1d(  all_names, [int(n) for n in names]))

        return rho_tg[ind],theta_tg[ind],R[ind],Z[ind]

    
    def signal_info(self,group,name,time):

        rho_tg = self.get_rho(group,[name,],time)[0]
                
        all_names = self.get_names(group)
       
        ind = all_names.index(name)

        phi = self.Phi_start[group][ind] 
        info = group+' '+str(name)+' Phi: %d deg, rho_tg: %.3f'%(phi,rho_tg)
        return info
