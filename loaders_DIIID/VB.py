from loaders_DIIID.loader import * 
import os
from multiprocessing import  Pool

import MDSplus as mds

def check(shot):
    #fastest check if the shotfile exist
    #BUG 
    return True





class loader_VB(loader):

    def __init__(self,*args, **kargs):
        
        super(loader_VB,self).__init__(*args, **kargs)

        self.names = {'PHDMIDVB': list(range(1,17))}
        self.groups = ('PHDMIDVB',)


        self.tvec = None
        self.catch = {}


    
    def get_names(self,group):
        return self.names[group]
    
    def get_signal(self,group, names,calib=False,tmin=None,tmax=None):       
        
        if tmin is None:    tmin = self.tmin
        if tmax is None:    tmax = self.tmax
        

        if size(names) >1:
            return [self.get_signal(group, name,calib=calib,tmin=tmin,tmax=tmax) for name in names]
        name = names

        if name in self.catch:
            sig = self.catch[name]
        else:
            
            tree = 'SPECTROSCOPY'
            self.MDSconn.openTree(tree,self.shot)
            TDIcall = "_x=\\%s::%s%.2d"%(tree,group,name )
            #"PTDATA(\PHDMIDVB%.2d,0)'%ch
            #print TDIcall

            sig = self.MDSconn.get(TDIcall).data()
            if self.tvec is None:
                self.tvec = self.MDSconn.get('dim_of(_x)').data()/1e3
                
            self.MDSconn.closeTree(tree,self.shot)
            self.catch[name] = single(sig)

        imin,imax = self.tvec.searchsorted([tmin,tmax])
        ind = slice(imin,imax+1) 
        
        return self.tvec[ind], sig[ind]

    


    def signal_info(self,group,name,time):
        return ''
    
 
 
