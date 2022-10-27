from loaders_DIIID.loader import * 
import os
from multiprocessing import  Pool

#from time import time 

import MDSplus as mds

def check(shot):
    #fastest check if the shotfile exist
    #BUG 
    return True


lang_data = [\
('ulang18001',      1.20660,     -1.36550),
('ulang18003',      1.22150,     -1.36550),
('ulang18005',      1.23650,     -1.36550),
('ulang18007',      1.25150,     -1.36550),
('ulang18009',      1.26650,     -1.36550),
('ulang18011',      1.29150,     -1.36550),
('ulang18013',      1.30660,     -1.36550),
('ulang18015',      1.32150,     -1.36550),
('ulang18017',      1.33650,     -1.36550),
('ulang18019',      1.35150,     -1.36550),
('ulang18021',      1.36650,     -1.36550),
('ulang18023',      1.50030,     -1.25290),
('ulang18025',      1.52820,     -1.25290),
('ulang18027',      1.55620,     -1.25290),
('ulang18029',      1.58410,     -1.25290),
('ulang18031',      1.61210,     -1.25290),
('ulang18033',      1.64000,     -1.25290),
('ulang18035',      1.66790,     -1.25290),
('ulang18037',      1.69590,     -1.25290),
('ulang18039',      1.72380,     -1.25290),
('ulang18057',      1.01422,      1.03759),
('ulang18059',      1.01422,      1.08204),
('ulang18061',      1.01422,      1.12243),
('ulang18063',      1.01422,      1.14732),
('ulang18065',      1.01168,      1.16408),
('ulang18067',      1.06299,      1.17424),
('ulang18069',      1.12116,      1.17754),
('ulang18071',      1.23800,      1.25984),
('ulang18073',      1.26771,      1.31877),
('ulang18075',      1.28549,      1.35153),
('ulang18077',      1.30175 ,     1.35153),
('ulang18079',      1.31801,      1.35153),
('ulang18081',      1.33426,      1.35153),
('ulang18083',      1.35077,      1.35153),
('ulang18085',      1.36703,      1.35153),
('ulang18087',      1.38328,      1.35153),
('ulang18089',      1.28549,      1.35153),
('ulang18091',      1.30175,      1.35153),
('ulang18093',      1.31801,      1.35153),
('ulang18095',      1.33426,      1.35153),
('ulang18097',      1.35077,      1.35153),
('ulang18099',      1.36703,      1.35153),
('ulang18101',      1.38328,      1.35153),
('ulang18103',      1.37185,      1.29743),
('ulang18105',      1.42316,      1.24943),
('ulang18107',      1.48438,      1.19812),
('ulang18109',      1.56032,      1.13462),
('ulang18111',      1.68631,      1.07696),
('ulang18113',      1.019 ,      -1.026),
('ulang18115',      1.019,       -1.054),
('ulang18117',      1.019,       -1.081),
('ulang18119',      1.019,       -1.115),
('ulang18121',      1.019,       -1.143),
('ulang18123',      1.019,       -1.170),
('ulang18125',      1.019,       -1.198),
('ulang18127',      1.039,       -1.239),
('ulang18129',      1.058,       -1.259),
('ulang18131',      1.078,       -1.278),
('ulang18133',      1.107,       -1.308),
('ulang18135',      1.137,       -1.338)]

class loader_TPLANG(loader):

    def __init__(self,*args, **kargs):
        
        super(loader_TPLANG,self).__init__(*args, **kargs)
        #volatges, currents are even numbers 
        self.names = {'floor': r_[1:22:2],
                      'shelf': r_[23:40:2],
                      'up. div.': r_[57:113:2],
                      'low cnt. post.': r_[113:136:2]}
        
        self.groups = list(self.names.keys())


        self.tvec = None
        self.catch = {}


    
    def get_names(self,group):
        return self.names[group]
    
    def get_signal(self,group, name,calib=False,tmin=None,tmax=None):
        #BUG load each chunk separatelly? 
        
        
        if tmin is None:    tmin = self.tmin
        if tmax is None:    tmax = self.tmax
        
 
 
        if name in self.catch:
            sig = self.catch[name]
        else:
            TDIcall ='_x=PTDATA2("%s%.2d",%d)'%('tplang',name,self.shot)
           
            sig = self.MDSconn.get(TDIcall).data()
            if len(sig) < 2:
                raise Exception('LANGMUIR channel %d do not exist'%name)
            self.catch[name] = sig

        if self.tvec is None:
            self.tvec = self.MDSconn.get('dim_of(_x)').data()
            self.tvec /= 1e3 #s
                
        
        imin,imax = self.tvec.searchsorted([tmin,tmax])
        ind = slice(imin,imax+1) 
        
        return self.tvec[ind], sig[ind]

    


    def signal_info(self,group,name,time):
        for lname, r, z in lang_data:
            if lname == 'lang18%.3d'%name:
                return 'R: %.3fm  Z:%.3fm  rho: %.3f'(r,z,nan)
    
 
 
