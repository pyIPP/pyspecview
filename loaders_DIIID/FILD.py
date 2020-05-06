from loaders_DIIID.loader import * 
import os
from multiprocessing import  Pool

#from time import time 

import MDSplus as mds

def check(shot):
    #fastest check if the shotfile exist
    #BUG 
    return True





class loader_FILD(loader):

    def __init__(self,*args, **kargs):
        
        super(loader_FILD,self).__init__(*args, **kargs)

        self.names = {'FILD': list(range(1,25))}
        self.groups = ('FILD',)


        self.tvec = None
        self.catch = {}


    
    def get_names(self,group):
        return self.names[group]
    
    def get_signal(self,group, names,calib=False,tmin=None,tmax=None):
        #BUG load each chunk separatelly? 
        
        
        if tmin is None:    tmin = self.tmin
        if tmax is None:    tmax = self.tmax
        

        if size(names) >1:
            return [self.get_signal(group, name,calib=calib,tmin=tmin,tmax=tmax) for name in names]
        name = names

        if name in self.catch:
            sig = self.catch[name]
        else:
            TDIcall ='_x=PTDATA2("%s%.2d",%d)'%(group,name,self.shot)
            #print TDIcall
            sig = self.MDSconn.get(TDIcall).data()
            if len(sig) < 2:
                raise Exception('FILD channel %d do not exist'%name)
            #print sig
            self.catch[name] = sig

        if self.tvec is None:
                self.tvec = self.MDSconn.get('dim_of(_x)').data()
                self.tvec /= 1e3 #s
                
        imin,imax = self.tvec.searchsorted([tmin,tmax])
        ind = slice(imin,imax+1) 
        
        return self.tvec[ind], sig[ind]

    


    def signal_info(self,group,name,time):
        return ''
    
 
 
