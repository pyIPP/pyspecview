from .loader import * 
import os
import aug_sfutils as sf

#simple interface to load genetic signal 

def check(shot):
    return True


class loader_gen(loader):
    

    def __init__(self,*args, **kargs):
        
        super(loader_gen,self).__init__(*args, **kargs)

        self.names = {}
        
        self.groups = []
        shot_folder = shot_path+'%d'%(self.shot//10)
        for folder in ['L0', 'XX']:
            path = os.path.join(shot_folder,folder)
            if os.path.isdir(path):
                names = os.listdir(os.path.join(shot_folder,folder))
                for name in names:
                    if len(os.listdir(os.path.join(shot_folder,folder,name)))>0:
                        self.groups.append(name)
                
        self.groups.sort()
  

    def get_names(self,group):

        names_list = []
        sfo = sf.SFREAD(group, self.shot, experiment=self.exp, edition=self.ed)
        if sfo.status:
            names = sfo.objects + sfo.parsets

            for n in names:
                typ = sfo.sfh[n].obj_type
                if typ == 7:
                    names_list.append(n)
                if typ == 6:
                    ndim = min(info.ind[:2])  #BUG assume that 2D array, and ntime > nch
                    for i in range(1,ndim+1):
                        names_list.append(n+':%d'%i)
            
        else:
            raise Exception('Shotfile could not be opened')
        
        return names_list

    
    def get_signal(self, group, name, calib=False, tmin=None, tmax=None):


        if tmin is None:    tmin = self.tmin
        if tmax is None:    tmax = self.tmax
            
        sfo = sf.SFREAD(group, self.shot, experiment=self.exp, edition=self.ed)
        
        name+= ' '
        i_split = name.find(':')

        tvec = sfo.gettimebase(name[:i_split])

        if tvec is None:
            raise Exception('No timebase has been found!!')
        
        dt = np.diff(tvec)
        if np.std(dt) > np.mean(dt)/5:
            print('NWarning: non equally spaced timebase!!')

        nbeg, nend = tvec.searchsorted((tmin,tmax))
        
        sig = sfo(name[:i_split])

        if i_split != -1:
            ch = int(name[i_split+1:])-1
            sig = sig if sig.shape[0]  == len(tvec) else sig.T
            sig = sig[nbeg:nend+1,ch]

        return tvec[nbeg:nend+1],sig


    def signal_info(self,group,name,time):
        info = group+': '+name
        return info
    
 
