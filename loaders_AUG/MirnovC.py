from .loader import * 
import os, sys


shotfiles = 'MHI', 'MHA','MHB','MHC','MHD','MHE','MHF','MHG','MHH','MTR'

def check(shot):
    #fastest check if the shotfile exist

    for shotfile in shotfiles:
        path = shot_path+str(shot//10)+'/MH/'+shotfile+'/'+str(shot)
        if os.path.isfile(path):
            return True 
                
    path = shot_path+str(shot//10)+'/L0/'+'MTR'+'/'+str(shot)
    if os.path.isfile(path):
        return True 
                
    return False


class loader_MirnovCoils(loader):

    mode_range = (-6,5)
    pol_mode_num = True

    def __init__(self,*args, **kargs):

        super(loader_MirnovCoils,self).__init__(*args, **kargs)
        
        self.old_shotfile = self.shot < 19101
    
        if self.old_shotfile:
            self.dd.Open('MTR',self.shot)
            names = self.dd.GetNames()  
            names =  [n[1:] for n in names if n[:2] == 'PC' and n[-3] == '-']
        elif self.shot < 21496:
            self.dd.Open('MHA',self.shot)
            names = self.dd.GetNames()  
            names =  [n[1:] for n in names if n[:2] == 'CC' and n[-3] == '-']
            self.old_shotfile = True
        else:
            calib_shot = dd.PreviousShot('CMH', self.shot, experiment='AUGD')
            self.dd.Open('CMH',calib_shot)
            names = self.dd.GetNames()  
            names =  [n[1:] for n in names if n[:2] == 'CC' and n[-3] == '-']
            
        self.dd.Close()

        names.sort()
        self.groups = np.unique([n[:3] for n in names])
        self.names = {g:[] for g in self.groups  }
        
        for shotfile in shotfiles:

            path = shot_path+str(self.shot//10)+'/MH/'+shotfile+'/'+str(self.shot)
            if shotfile == 'MTR':
                path = shot_path+str(self.shot//10)+'/L0/MTR/'+str(self.shot)
    
            if not os.path.isfile(path):    continue
            
            self.dd.Open( shotfile,self.shot, experiment=self.exp, edition=self.ed)

            sf_names = self.dd.GetNames()
            self.dd.Close()

            sf_names.sort()
            
            for n in sf_names:
                if n in names:
                    self.names[n[:3]].append((shotfile, n))
                    names.remove(n)  #do not include the same coil twice, if it is in more different shotfiles
    
    def get_names(self,group):
        names = [n for sf,n in self.names[group]]
        names.sort()
        return names  
    
    def get_signal(self,group, name,calib=False,tmin=None,tmax=None):
        
        if tmin is None:    tmin = self.tmin
        if tmax is None:    tmax = self.tmax
            
                
        sig_names = self.names[group]
        for self.shotfile,sig in sig_names:
            if sig == name: break
        
        self.dd.Open(self.shotfile, self.shot, experiment=self.exp, edition=self.ed)
        tvec = self.dd.GetTimebase(name, cal=True)  #BUG make it fastrer!! do not reed the whole time vector!!
        nbeg, nend = tvec.searchsorted((tmin,tmax))

        sig = self.dd.GetSignal(name,cal=calib, nbeg= nbeg,nend = nend)

        return tvec[nbeg:nend+1],sig
        
    def get_names_phase(self):
        return ['Mcoils',]
    
    def get_signal_phase(self,name,calib=False,tmin=None,tmax=None,rem_elms=True):
        if tmin is None:    tmin = self.tmin
        if tmax is None:    tmax = self.tmax

        phase_signals = self.names['C09']
        
        sigs = []
        length = infty
        self.shotfile = None
        for sf, sig_p in phase_signals:
            if self.old_shotfile  and sf == 'MHA': continue
            if self.shotfile != sf :
                self.dd.Open(sf, self.shot, experiment=self.exp, edition=self.ed)
                self.shotfile = sf
                tvec = self.dd.GetTimebase(sig_p)
                nbeg, nend = tvec.searchsorted((tmin,tmax))

            sig = self.dd.GetSignal(sig_p,cal=calib,nbeg=nbeg,nend=nend)
            length = min(length, len(sig))
            sigs.append(sig)

        self.dd.Close()
        data = np.vstack([s[:length] for s in sigs]).T
        tvec = tvec[nbeg:nend+1]
        
        #remove elms from the data 
        print( 'removing elms')
        if rem_elms and self.dd.Open('ELM', self.shot):
            nbeg = tvec.searchsorted(self.dd.GetTimebase('t_begELM')-0.001)
            nend = tvec.searchsorted(self.dd.GetSignal('t_endELM'))
            
            for n1,n2 in zip(nbeg,nend):
                n1,n2 = max(1,n1), min(n2, len(tvec)-2)
                data[n1:n2] =  (np.linspace(0,1,n2-n1)[:,None]*(np.float_(data[n2+1])-data[n1-1])+data[n1-1])
   
        print( 'done')
        return tvec,data
    
    def get_theta_pol(self,name,tshot = 4 ,rhop=.2 ):
        
        try:
            magr,magz, theta0, theta_star = self.mag_theta_star( tshot,rhop )
        except:
            print( 'Error: get_theta_pol: mag_theta_star do not work!')
            print( traceback.format_exc())

            theta0     = np.linspace(-np.pi, np.pi, 10)
            theta_star = np.linspace(-np.pi, np.pi, 10)

        if name in  ['Mcoils',]:
            theta = []
            
            if self.shot > 19090:
                calib_shot = dd.PreviousShot('CMH', self.shot, experiment='AUGD')
                self.dd.Open('CMH', calib_shot)
                pre = 'C'
            else:
                self.dd.Open('MTR', self.shot)
                pre = 'P'

            for sf, sig in self.names['C09']:
                theta.append(self.dd.GetParameter(pre+sig,'theta'))
                
            self.dd.Close()

            theta = np.interp(theta,  theta0, theta_star)

        return -np.squeeze(theta) #BUG minus is there to have a positive m numbers for ordinary MHD modes with coinjected NBI 
            

    def signal_info(self,group,name,time):
        if name == 'Mcoils' or self.old_shotfile: return ' '

        calib_shot = dd.PreviousShot('CMH', self.shot, experiment='AUGD')
        self.dd.Open('CMH',calib_shot)
            
        theta = self.dd.GetParameter('C'+name,'theta')
        self.dd.Close()
        try:
            info = str(name)+' theta: %.2fdeg'%np.rad2deg(theta)
        except:
            print(( group,name,time,theta))
            raise
        
        return info
