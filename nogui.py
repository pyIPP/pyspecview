#!/usr/bin/env python
# -*- coding: utf-8 -*-

try:
    from PyQt5.QtCore import *
    from PyQt5.QtGui import *
    from PyQt5.QtWidgets import *
    import configparser
except:
    from PyQt4.QtCore import *
    from PyQt4.QtGui import *
    import ConfigParser as configparser

import numpy as np
import sys, os
import time
from pyspecview import MainGUI


def debug():
    
    shot = 167542
    app = QApplication(sys.argv)

    form = MainGUI(shot )
    time.sleep(5)
    for i in range(10): 
        print(i, form.cross_diagnostics.itemText(i))
        if form.cross_diagnostics.itemText(i) == 'SXR': break
    form.change_diagnostics(i)

    form.change_sig_group(4)
    form.change_signal(1)
    form.SpecWin.t_range = 3.3,3.4
    form.SpecWin.f_range = 4e3,5e3

    for i in range(10):  
        if form.cross_diagnostics.itemText(i) == 'SXR': break
    form.change_diagnostics_radial(i)
    for i in range(10):  
        if form.cross_sig_group.itemText(i) == '90RM1': break
    form.change_sig_group_radial(i)

    form.radial_view.f_range = form.SpecWin.f_range

    print(i)
    form.curr_tab = 2

    form.radial_view.t_range =  form.SpecWin.t_range
    
    form.change_signal_radial( 7,update=True)
    form.updatePanels(3)


def make_movie():
    
    shot = 175900
    app = QApplication(sys.argv)

    form = MainGUI(shot )
    time.sleep(5)

    for idiag in range(10):  
        if form.cross_diagnostics.itemText(idiag) == 'ECE': break
    form.change_diagnostics(idiag)
    form.change_sig_group(0)
    form.change_signal(0)
    form.dR_corr = .02
    form.dZ_corr = .01

    for i in range(10):  
        if form.cross_diagnostics.itemText(i) == 'ECE': break
    form.change_diagnostics_radial(i)
    for i in range(10):  
        if form.cross_sig_group.itemText(i) == 'ECE': break
    form.change_sig_group_radial(i)
    
    # Full 2D profile 0, substract background 1
    form.Te2Dmap.UpdatePlotType(1)
    out_folder = '~/pyspecviewer/movie'
    out_folder = os.path.expanduser(out_folder)
    
    if not os.path.exists(out_folder):
        os.mkdir(out_folder)

    N = 100
    T = np.linspace(2.03,2.08 ,N)
    F = np.linspace(6, 5.5,N)*1e3
    for isig in range(30):  
        if form.cross_signal.itemText(isig) == '1': break
    isig = 0
    form.Te2Dmap.m = 1
    form.Te2Dmap.substract = True
    form.Te2Dmap.rmax = .2
    
    for i in range(N):

        form.radial_view.f_range = F[i]*0.7, F[i]*1.3
        form.frange = F[i]*0.7, F[i]*1.3

        print(i)
        form.curr_tab = 2
        #dT = mean(diff(T))
        dT = max(0.001, mean(diff(T))*2)
        
        form.radial_view.t_range = T[i]-dT, T[i]+dT
        
        form.change_signal_radial( isig,update=True)
        form.updatePanels(3)
        form.canvas_2Dmap.print_figure('./movie/%.3d.png'%i)

    os.system('mencoder "mf://'+out_folder+'/*.png" -o '+out_folder\
        +'/movie_'+str(shot)+'.avi -ovc lavc -lavcopts vcodec=mpeg4:mbd=1:vbitrate=3500')


def trace_mode():
    
    shot = 32148
    app = QApplication(sys.argv)

    form = MainGUI(shot )
    time.sleep(5)
    
    #use ECE
    for i in range(10): 
        if form.cross_diagnostics.itemText(i) == 'ECE': break
    form.change_diagnostics(i)#te
    form.change_sig_group(2) #anything
    form.change_signal(1)#anything
    form.dR_corr = -0.0
    form.dZ_corr = -0.0
 
    #set scrossphase diagnostic
    for i in range(10):  
        if form.cross_diagnostics.itemText(i) == 'SXR': break
    form.change_diagnostics_radial(i)
    for i in range(10):  
        if form.cross_sig_group.itemText(i) == 'I': break
    form.change_sig_group_radial(i)
    for isig in range(30):  
        if form.cross_signal.itemText(isig) == 'I_049': break

    N = 200
    T = np.linspace(2., 8.5 ,N)
    form.radial_view.f_range = 15e3, 20e3
    
    print(form.radial_view.f_range)

    R = []
    Z = []
    phi = []
    Amp = []

    for i in range(N-1):
        print(i)
        form.curr_tab = 2
        dT = np.mean(np.diff(T))
        form.radial_view.t_range = T[i]-dT, T[i]+dT
        form.change_signal_radial( isig,update=True)
        phi.append(form.radial_view.phi2)
        R.append(form.radial_view.R)
        Z.append(form.radial_view.Z)
        Amp.append(form.radial_view.amplitude2)
        
    import IPython
    IPython.embed()


if __name__ == "__main__":

    make_movie()
