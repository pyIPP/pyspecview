#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os, sys, logging
from scipy import sparse
import matplotlib.pylab as plt
import numpy as np
import time 
from scipy.ndimage.interpolation import map_coordinates
from scipy.interpolate import interp1d
from scipy.linalg import eigh
from scipy.stats.mstats import mquantiles
from multiprocessing import Process,Queue
from matplotlib.ticker import MaxNLocator
from matplotlib.widgets import MultiCursor
#from .pyspecview import extract_harmonics

from IPython import embed
from copy import deepcopy
try: #only for AUG
    import aug_sfutils as sf
except:
    pass


try:
    from PyQt5.QtCore import *
    from PyQt5.QtGui import *
    from PyQt5.QtWidgets import *
    from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
    from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT
except:
    from PyQt4.QtCore import *
    from PyQt4.QtGui import *
    from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
    from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT
 
from . import fconf

logger = logging.getLogger('pyspecview.roto_tomo')
logger.setLevel(logging.DEBUG)

from .SVDfilter import SVDFilter
from sksparse.cholmod import  cholesky,CholmodWarning,CholmodError
from matplotlib import colors,cbook,cm

cdict = {
'red': ((0, 0.1, 0.1),(0.04, 0.15, 0.15),(0.07, 0.25, 0.25), (0.12, 0.4, 0.4),\
    (0.2, 0.6, 0.6),(0.35, 0.9, 0.9),(0.4, 0.9, 0.9),(0.6, 1, 1),(0.8, 1, 1),(0.9, 1, 1),(1, 1, 1)),
'green': ((0, 0.1, 0.1),(0.04, 0.15, 0.15),(0.07, 0.2, 0.2), (0.12, 0.15,0.15 ),\
    (0.2, 0.1, 0.1),(0.35, 0.1, 0.1),(0.4, 0.1, 0.1),(0.6, 0.6, 0.6),(0.8, 1, 1),(0.9, 1, 1),(1, 1, 1)),
'blue' : ((0, 0.3, 0.3),(0.04, 0.5, 0.5),(0.07, 0.6, 0.6), (0.12, 0.6,0.6),\
    (0.2, 0.6, 0.6),(0.35, 0.6, 0.6),(0.4, 0.1, 0.1), (0.6, 0, 0),(0.8, 0, 0),(0.9, 0.8, 0.8),(1, 1, 1))
}
my_cmap = colors.LinearSegmentedColormap('my_colormap', cdict, 1024)

tomo_local_path = '~/tomography/'
loc_dir = os.path.dirname(os.path.abspath(__file__))

if os.path.exists('$PYTOMO'):  #AUG
    tomo_code_path = '$PYTOMO/'
#afs/ipp-garching.mpg.de/home/t/todstrci/pytomo_orig/'
elif os.path.exists('/fusion/projects/codes/pytomo/'): #DIII-D
     tomo_code_path = '/fusion/projects/codes/pytomo/'
else:  #local instalation
     tomo_code_path = tomo_local_path

 


tomo_code_path  = os.path.expanduser(os.path.expandvars(tomo_code_path))
tomo_local_path = os.path.expanduser(os.path.expandvars(tomo_local_path))

sys.path.append(tomo_code_path)

#load some modules from pytomo 
from mat_deriv_B import mat_deriv_B
from shared_modules import  read_config
from graph_derivation import graph_derivation
from geom_mat_setting import geom_mat_setting,prune_geo_matrix
from annulus import get_bd_mat, get_rho_field_mat
import config
 
class NavigationToolbar2(NavigationToolbar2QT):
    # only display the buttons we need
    toolitems = [t for t in NavigationToolbar2QT.toolitems if t[0] in ['Pan', 'Zoom']]


class DataSettingWindow(QMainWindow):
    dpi = 100
    fontsize = 10
    def __init__(self, parent, tomo):
        QMainWindow.__init__(self, parent)
        self.tomo = tomo
        self.show_harm = False
        if self.tomo.sxr_harmonics is not None:
            self.show_harm = True
        elif not hasattr(self.tomo,'tvec'):
            QMessageBox.warning(self.parent,"No data",
                    "Data wre not loaded! Missing fast SXR?",QMessageBox.Ok)
            return
        
            
        
        self.parent = parent
        self.det_num = []
        cWidget = QWidget(self)
        self.bg_col = cWidget.palette().color(QPalette.Base)

        self.setWindowTitle('Data preprocessing tool')
        self.create_status_bar()
        self.create_main_frame()

        if self.show_harm:
            self.create_Harm_table()
        else:
            #self.tomo.SVDF = SVDF
            self.create_Data_table()
            self.create_SVD_table()
            self.create_residuum_table()
        
        self.main_tab.setCurrentIndex(0)
        self.main_tab.setTabEnabled(0,True)

    def create_main_frame(self):

        self.cWidget = QWidget(self)
        self.gridLayout = QGridLayout(self.cWidget)
        self.setWindowTitle('Data selection / Preprocessing')
        self.resize(600, 650)

        self.Expand = QSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        self.Fixed  = QSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed)

        self.main_tab = QTabWidget(self)
        self.main_tab.setTabPosition(QTabWidget.North)

        self.gridLayout.addWidget(self.main_tab, 0, 0, 1,2)
        self.setCentralWidget(self.cWidget)        
        self.setCenter()
        self.main_tab.currentChanged.connect(self.switch_tab)
        self.switch_tab(0)
        
    def setCenter(self):
        frect = QDesktopWidget.frameGeometry(self)
        frect.moveCenter(QDesktopWidget().availableGeometry(self).center());
        self.move(frect.topLeft())
            

   
    def create_Harm_table(self):
        #plot 0th and 1-n th harmonics of the input data
        self.tab_widget_harm = QWidget(self)
        
        self.main_tab.insertTab(0,self.tab_widget_harm, 'Data Harmonics')

        self.tab_widget_harm.setSizePolicy(self.Expand)
        self.verticalLayout_harm = QVBoxLayout(self.tab_widget_harm)
        self.horizontalLayout_harm = QHBoxLayout()


        self.fig_harm = plt.Figure((7.0, 5.0), dpi=self.dpi)
        self.tab_widget_harm.setToolTip('Use left/right mouse button to remove/return point')

        self.canvas_harm = FigureCanvas(self.fig_harm)
        self.canvas_harm.setParent(self.cWidget)
        
        self.fig_harm.patch.set_facecolor((self.bg_col.red()/255., 
                                           self.bg_col.green()/255., 
                                           self.bg_col.blue()/255.))

        self.canvas_harm.setSizePolicy(self.Expand)

        self.verticalLayout_harm.addWidget(self.canvas_harm)

        self.link = self.fig_harm.canvas.mpl_connect('pick_event', self.onpick_data)
    
        #prepare first data
        self.groupBox = QGroupBox("")
        self.gridLayout = QGridLayout(self.groupBox)
        self.labelDetectors = QLabel( "Wrong detectors:")
        self.Edit_wrongDetectors = QLineEdit()
        self.gridLayout.addWidget(self.labelDetectors, 0, 0)
        self.gridLayout.addWidget(self.Edit_wrongDetectors, 0, 1)

        self.verticalLayout_harm.addWidget(self.groupBox)

        self.OkButton=QPushButton("Close")
    
        self.gridLayout.addWidget(self.OkButton, 0, 2)


        #self.mpl_toolbar_harm = NavigationToolbar2QT(self.canvas_harm, self.cWidget)
        #self.verticalLayout_harm.addWidget(self.mpl_toolbar_harm)
        
        
        #vbox = QVBoxLayout()
        #self.verticalLayout_harm.addWidget(self.canvas_svd)
        

        label_slider_nharm = QLabel('N HARM:')
        self.nharm_slider = QSlider(Qt.Horizontal)
        self.nharm_slider.setToolTip('Set number of filtered harmononics')
        
        nharm_max, self.nch = np.shape(self.tomo.all_bb)

        self.nharm_slider.setRange(1, nharm_max)
        self.nharm_slider.setTracking(False)
        self.nharm_slider.setTickPosition(QSlider.TicksBelow)
        self.nharm_slider.setTickInterval(1)
        self.nharm_slider.setValue(nharm_max)
        self.nharm_slider.setMinimumSize(50,20)
         

       # spacer = QSpacerItem(0,0,)

        hbox = QHBoxLayout()
        hbox.addWidget(label_slider_nharm)
        hbox.setAlignment(label_slider_nharm, Qt.AlignRight)
        hbox.addWidget(self.nharm_slider)
        hbox.setAlignment(self.nharm_slider, Qt.AlignRight)
         
        self.verticalLayout_harm.addLayout(hbox) 
        def apply_nharm_slider(self):
            self.nharm = self.nharm_slider.value()
            self.RefreshEvent()

            
        self.nharm.valueChanged.connect(self.apply_nharm_slider)


        #tooltips
        self.Edit_wrongDetectors.setToolTip('Insert wrong detectors, use ":" for the intervals - 1,,3,5:40,...')
        self.Edit_wrongDetectors.editingFinished.connect(self.dets_list_edited)

  
        plt.rcParams['xtick.direction'] = 'out'
        plt.rcParams['ytick.direction'] = 'out'
        self.fig_harm.clf()
        self.fig_harm.subplots_adjust(left=.15 ,top=.90,right=.97,bottom=.05, hspace=0.03, wspace = 0.03)

        self.ax_harm0 = self.fig_harm.add_subplot(211)
        self.ax_harm1 = self.fig_harm.add_subplot(212, sharex=self.ax_harm0)

        self.ax_harm0.point_label = self.ax_harm0.text(0,0,'', fontsize=8 ,zorder=99)
        self.ax_harm1.point_label = self.ax_harm1.text(0,0,'', fontsize=8 ,zorder=99)

        
 
        self.err_plots = []
        self.err_plots.append( self.ax_harm0.errorbar(0, 0, 0, fmt='.-',picker=True,
                                capsize=0, label='0. harmonic'))


        self.invalid_plot, = self.ax_harm0.plot([], [], 'ro',picker=True)

        c = 'b', 'g', 'r', 'y', 'm'
        for iharm in range(1,self.nharm):
            self.err_plots.append(self.ax_harm1.errorbar(0, 0, 0, fmt='.-'+c[iharm%len(c)], 
                            capsize=0, zorder=self.nharm-iharm, picker=True,
                            label=f'{iharm}. harmonic'))

        self.ax_harm0.legend(loc='upper right', fontsize  =  self.fontsize) 
        self.ax_harm1.legend(loc='upper right', fontsize  =  self.fontsize) 
        
        self.ax_harm1.set_xlabel('Channel',fontsize= self.fontsize)

        self.ax_harm0.xaxis.set_pickradius(2)
        self.ax_harm0.yaxis.set_pickradius(2)
        self.ax_harm1.xaxis.set_pickradius(2)
        self.ax_harm1.yaxis.set_pickradius(2)


            
        for ind in self.tomo.tok.dets_index[:-1]:
            self.ax_harm0.axvline(x=1.5+np.amax(ind), linestyle='--' ,color='k')
            self.ax_harm1.axvline(x=1.5+np.amax(ind), linestyle='--',color='k')

        for label in (self.ax_harm0.get_yticklabels() +self.ax_harm0.get_xticklabels()\
                     + self.ax_harm1.get_yticklabels() + self.ax_harm1.get_xticklabels()):
            label.set_fontsize( self.fontsize) # Size here overrides font_prop
    
        self.ax_harm1.xaxis.set_major_locator(MaxNLocator(12))

        xax = self.ax_harm0.axes.get_xaxis()
        xax = xax.set_visible(False)

        #add upper axis with camera names 
        xlabels = [np.median(ind) for ind in self.tomo.tok.dets_index]
        labels = list(self.tomo.tok.detectors_dict.keys())
        ax2 = self.ax_harm0.twiny()
        ax2.set_xticks(xlabels)
        ax2.set_xticklabels(labels,fontsize= self.fontsize)
        ax2.set_xlim(0,self.nch)
        ax2.tick_params(axis='x', which='major', pad=0) 
        # critical for picking!!!! 
        ax2.set_zorder(0) 
        ax2.set_autoscaley_on(True) 
        self.ax_harm0.set_zorder(1) 
        self.fig_harm.sca(self.ax_harm0) 
        self.multi_harm = MultiCursor(self.canvas_harm, [self.ax_harm0,self.ax_harm1], horizOn=False, color='k', lw=1)

        self.ax_harm0.set_xlim(0,self.nch)

        wrong =  self.ReadWrongDets(init=True)
        self.SetWrongDets(wrong)
        
        #fill plots with data
        self.RefreshEvent()
        self.OkButton.clicked.connect(self.closeEvent)
        
    def create_Data_table(self):
        
        self.tab_widget_data = QWidget(self)
        
        self.main_tab.insertTab(0,self.tab_widget_data, 'Data')

        self.tab_widget_data.setSizePolicy(self.Expand)
        self.verticalLayout_data = QVBoxLayout(self.tab_widget_data)
        self.horizontalLayout_data = QHBoxLayout()


        self.fig_data = plt.Figure((7.0, 5.0), dpi=self.dpi)
        self.tab_widget_data.setToolTip('Use left/right mouse button to remove/return point')

        self.canvas_data = FigureCanvas(self.fig_data)
        self.canvas_data.setParent(self.cWidget)
        
        self.fig_data.patch.set_facecolor((self.bg_col.red()/255., 
                                           self.bg_col.green()/255., 
                                           self.bg_col.blue()/255.))

        self.canvas_data.setSizePolicy(self.Expand)

        self.verticalLayout_data.addWidget(self.canvas_data )

        self.link = self.fig_data.canvas.mpl_connect('pick_event', self.onpick_data)
      

        #prepare first data
        self.groupBox = QGroupBox("")
        self.gridLayout = QGridLayout(self.groupBox)
        self.labelDetectors = QLabel( "Wrong detectors:")
        self.Edit_wrongDetectors = QLineEdit()
        self.gridLayout.addWidget(self.labelDetectors, 0, 0)
        self.gridLayout.addWidget(self.Edit_wrongDetectors, 0, 1)

        self.verticalLayout_data.addWidget(self.groupBox)

        self.OkButton=QPushButton("Close")
    
        self.gridLayout.addWidget(self.OkButton, 0, 2)

        #tooltips
        self.Edit_wrongDetectors.setToolTip('Insert wrong detectors, use ":" for the intervals - 1,,3,5:40,...')
        self.Edit_wrongDetectors.editingFinished.connect(self.dets_list_edited)

        if not hasattr(self.tomo,'tvec'):
            QMessageBox.warning(self.parent,"No data",
                    "Data wre not loaded! Missing fast SXR?",QMessageBox.Ok)
            return
            
        tvec = self.tomo.tvec
        data = self.tomo.signals
        self.nch = data.shape[1]
        
 
        plt.rcParams['xtick.direction'] = 'out'
        plt.rcParams['ytick.direction'] = 'out'
        self.fig_data.clf()
        self.fig_data.subplots_adjust(left=.15-.05,top=.95,right=.95+.05,bottom=.12)

        self.ax_data = self.fig_data.add_subplot(111)
        
 


        self.ax_data.point_label = self.ax_data.text(0,0,'', fontsize=8 ,zorder=99)
        extent = (0.5, self.nch+.5,tvec[0],tvec[-1])
        self.data_im = self.ax_data.imshow(np.ones((2,2)),cmap = 'jet',extent=extent
                    ,interpolation='nearest', origin='lower', picker=True)
        self.ax_data.axis('tight')
        self.ax_data.set_ylabel('Time [s]',fontsize= self.fontsize)
        self.ax_data.set_xlabel('Channel',fontsize= self.fontsize)
        self.ax_data.xaxis.set_pickradius(2)
        self.ax_data.yaxis.set_pickradius(2)
        from make_graphs import LogFormatterTeXExponent

        cbar = self.fig_data.colorbar(self.data_im,format=LogFormatterTeXExponent('%.2e'))
        self.ax_data.axis((0.5, self.nch+.5,tvec[0],tvec[-1]))
        cbar.ax.tick_params(labelsize= self.fontsize ) 
        cbar.locator = MaxNLocator(nbins=5)
        cbar.update_ticks()

 
            
        for ind in self.tomo.tok.dets_index[:-1]:
            self.ax_data.axvline(x=1.5+np.amax(ind), linestyle='-' ,color='k')
            self.ax_data.axvline(x=1.5+np.amax(ind), linestyle='--',color='w')

        for label in (self.ax_data.get_xticklabels() + self.ax_data.get_yticklabels()):
            label.set_fontsize( self.fontsize) # Size here overrides font_prop
    
        self.ax_data.xaxis.set_major_locator(MaxNLocator(12))


        xlabels = [np.median(ind) for ind in self.tomo.tok.dets_index]
        labels = list(self.tomo.tok.detectors_dict.keys())
        ax2 = self.ax_data.twiny()
        ax2.set_xticks(xlabels)
        ax2.set_xticklabels(labels,fontsize= self.fontsize)
        ax2.set_xlim(0, self.nch)
        ax2.tick_params(axis='x', which='major', pad=0) 
        # critical for picking!!!! 
        ax2.set_zorder(0) 
        ax2.set_autoscaley_on(True) 
        self.ax_data.set_zorder(1) 
        self.fig_data.sca(self.ax_data) 


        wrong = self.ReadWrongDets(init=True)
        self.SetWrongDets(wrong)

        self.RefreshEvent()
        self.OkButton.clicked.connect(self.closeEvent)

    def RefreshEvent(self): 

        if self.show_harm:
           bb = np.copy(self.tomo.all_bb)
           bb_err = np.copy(self.tomo.all_bb_err)
           x = np.arange(self.nch)+1
           
           ind_correct = self.tomo.tok.get_correct_dets()
 

           data_max = np.nanmax(np.abs(bb[0][ind_correct]))
           if data_max > 2e5:
               fact,pre = 1e6, 'MW'
           elif  data_max > 2e2:
               fact,pre = 1e3, 'kW'
           else:
               fact,pre = 1e0, 'W'


           self.ax_harm0.set_ylabel('Mean signal  ['+pre+'/m$^2$]', fontsize=self.fontsize)
           self.ax_harm1.set_ylabel('Amplitudes ['+pre+'/m$^2$]', fontsize=self.fontsize)
           self.invalid_plot.set_data(x[~ind_correct], np.abs(bb[0,~ind_correct])/fact)
           bb[:,~ind_correct] = np.nan
           #TODO plot also phase?
           for i, err_plot in enumerate(self.err_plots[:self.nharm]):
               y,yerr = np.abs(bb[i])/fact, bb_err[i]/fact
               plotline, caplines, barlinecols = err_plot
               plotline.set_data(x, y)
               # Find the ending points of the errorbars
               error_positions = (x, y-yerr), (x, y+yerr)
               # Update the error bars
               barlinecols[0].set_segments(zip(zip(x, y-yerr), zip(x, y+yerr)))
     
           self.ax_harm0.set_ylim(0, data_max/fact)
           self.ax_harm1.set_ylim(0, np.nanmax(np.abs(bb[1:]))/fact)
           self.canvas_harm.draw()
           #load data with new set of enabled detectors by self.dets
           #ind_correct = self.tomo.tok.get_correct_dets()

           #self.tok.dets = self.dets = np.array(dets)
           #self.tomo.bb = [bb[ind_correct] for bb in self.tomo.all_bb]
           #self.tomo.bb_err = [be[ind_correct] for be in self.tomo.all_bb_err]


        else:
            tvec = self.tomo.tvec
            data = np.copy(self.tomo.signals)
            if np.nanmax(data) > 2e5:
                fact,pre = 1e6, 'MW'
            elif np.nanmax(data) > 2e2:
                fact,pre = 1e3, 'kW'
            else:
                fact,pre = 1e0, 'W'

            ind_correct = self.tomo.tok.get_correct_dets(data.T)
            data[:,~ind_correct] = np.nan
            #ind = slice(None,None, max(1, len(tvec)//2000)) #downsample if necessary
            ndown = int(np.ceil(len(tvec)/1000))
            nt = len(tvec)//ndown*ndown
            #BUG test it!!!
            data = data[:nt].reshape(nt//ndown, -1, len(data.T)).mean(1)
            tvec = tvec[:nt].reshape(nt//ndown, -1).mean(1)
            vmax = mquantiles(data[np.isfinite(data)], 0.99)[0]

            self.data_im.set_array(data/fact)
            self.data_im.set_clim(0,vmax/fact)
            self.canvas_data.draw()
 
            if hasattr(self.tomo,'SVDF'):
                self.tomo.SVDF.set_corrupted_ch(~ind_correct)
        
        #write back to roto_tomo
        self.tomo.dets = self.tomo.tok.dets[ind_correct]
        
     

    def dets_list_edited(self,event=None):    
        wrong = self.ReadWrongDets()
        self.SetWrongDets(wrong)
        self.RefreshEvent()
        
        
    def onpick_data(self,event):
        #mouse interaction, remove corrupted channels
        ax = event.mouseevent.inaxes
        ch = int(round(event.mouseevent.xdata)) 

        if event.mouseevent.dblclick:
            wrong = list(self.ReadWrongDets())
            
            if event.mouseevent.button == 1: #delete
                wrong.append(ch-1)
            
            if event.mouseevent.button == 3 and ch-1 in wrong: #delete
                 wrong.remove(ch-1)

            self.SetWrongDets(wrong)
            self.ReadWrongDets()
            self.RefreshEvent()


    def ReadWrongDets(self,init=False):
        
        wrong_dets_pref = []
        wrong_dets_str = str(self.Edit_wrongDetectors.text())

        if wrong_dets_str.strip() != '':
            
            try:
                dets_list = eval('np.r_['+str(wrong_dets_str)+']')
                wrong_dets_pref = np.array(np.int_(dets_list)-1, ndmin=1)
            except:
                QMessageBox.warning(self.parent,"Input problem",
                    "Wrong detectors are in bad format, use ...,10,11,13:15,...",QMessageBox.Ok)
                #raise
        elif init:     # do not remove on init
            wrong_dets_pref = np.int_(config.wrong_dets_pref)
        
        wrong = self.tomo.tok.dets[~self.tomo.tok.get_correct_dets(include_pref=False)]
      
        config.wrong_dets_pref  = np.unique(np.setdiff1d( wrong_dets_pref,wrong))
        
        return  wrong_dets_pref
    
    def SetWrongDets(self, wrong):
        wrong = np.unique(wrong)

        wrong_dets_str = ''
        i = 0
        while i < len(wrong):
            wrong_dets_str = wrong_dets_str+str(wrong[i]+1)
            if i+1 < len(wrong) and wrong[i+1] == wrong[i]+1:
                while i+1 < len(wrong) and wrong[i+1] == wrong[i]+1:
                    i += 1
                wrong_dets_str = wrong_dets_str+':'+str(wrong[i]+2)
            i += 1
            wrong_dets_str = wrong_dets_str+','
        wrong_dets_str = wrong_dets_str[:-1]
            
        self.Edit_wrongDetectors.setText(wrong_dets_str)


        
    def create_SVD_table(self):
    
        self.tab_widget_svd = QWidget(self)
        self.main_tab.insertTab(2,self.tab_widget_svd, 'SVD filter')
        self.tab_widget_svd.setSizePolicy(self.Expand)
        self.verticalLayout_svd = QVBoxLayout(self.tab_widget_svd)
        self.horizontalLayout_svd = QHBoxLayout()
 

        self.verticalLayout_svd.addLayout(self.horizontalLayout_svd)
 

        self.fig_svd = plt.Figure((6.0, 5.0), dpi=self.dpi)
        self.canvas_svd = FigureCanvas(self.fig_svd)
        self.canvas_svd.setParent(self.cWidget)
          
        self.fig_svd.patch.set_facecolor((self.bg_col.red()/255., 
                                          self.bg_col.green()/255., 
                                          self.bg_col.blue()/255.))
         
        vbox = QVBoxLayout()
        vbox.addWidget(self.canvas_svd)
        
        hbox = QHBoxLayout()

        slider_label0 = QLabel('TAU:')
        slider_label1 = QLabel('N SVD:')
        slider_label2 = QLabel('N HARM:')
        
        self.tau  = QSlider(Qt.Horizontal)
        self.nsvd  = QSlider(Qt.Horizontal)
        self.nharm = QSlider(Qt.Horizontal)
        self.tau.setToolTip('Set time resolution')
        self.nsvd.setToolTip('Set number of used SVD components')
        self.nharm.setToolTip('Set number of filtered harmononics')
        
        
        self.tau.setRange(0,9)
        self.tau.setTracking(False)
        self.tau.setTickPosition(QSlider.TicksBelow)
        self.tau.setTickInterval(1)
        self.tau.setValue(self.tomo.tau)
        self.tau.setMinimumSize(50,20)
        
        self.nharm.setRange(2, 6)
        self.nharm.setTracking(False)
        self.nharm.setTickPosition(QSlider.TicksBelow)
        self.nharm.setTickInterval(1)
        self.nharm.setValue(self.tomo.n_harm)
        self.nharm.setMinimumSize(50,20)
        

        self.nsvd.setRange(2, 6)
        self.nsvd.setTracking(False)
        self.nsvd.setTickPosition(QSlider.TicksBelow)
        self.nsvd.setTickInterval(1)
        self.nsvd.setValue( self.tomo.n_svd)
        self.nsvd.setMinimumSize(50,20)

        spacer = QSpacerItem(0,0,)
 
        for i,w in enumerate([  slider_label0, self.tau, slider_label1, self.nsvd,slider_label2,self.nharm]):
            hbox.addWidget(w)
            hbox.setAlignment(w, Qt.AlignRight)
         
        vbox.addLayout(hbox)

        self.horizontalLayout_svd.addLayout(vbox)
    

        self.nsvd.valueChanged.connect(self.apply_slider)
        self.nharm.valueChanged.connect(self.apply_slider)
        self.tau.valueChanged.connect(self.apply_slider)

                      
        tvec = self.tomo.tvec
        data = np.copy(self.tomo.signals)
        err  = np.copy(self.tomo.error_sig)

        self.tomo.SVDF.actual = False
        self.tomo.SVDF.fig_svd = self.fig_svd
    
        
    def run_svd_filt(self):
        self.tomo.SVDF.run_filter()
        self.canvas_svd.draw()
        self.canvas_res.draw()
     
    def apply_slider(self):
        if self.tomo.SVDF.n_harm != self.nharm.value():
            self.tomo.SVDF.n_harm = self.nharm.value()
            self.statusBar().showMessage('N HARM = %d' % self.nharm.value(), 2000)

        if self.tomo.SVDF.n_svd != self.nsvd.value():
            self.tomo.SVDF.n_svd = self.nsvd.value()
            self.statusBar().showMessage('N SVD = %d' %  self.nsvd.value(), 2000)

        if self.tomo.SVDF.tau != self.tau.value():
            self.tomo.SVDF.tau = self.tau.value()
            self.statusBar().showMessage('TAU = %d' %  self.tomo.SVDF.tau, 2000)


        
        self.tomo.n_harm = self.nharm.value()
        self.tomo.n_svd = self.nsvd.value()
        self.tomo.tau = self.tau.value()
        

        self.tomo.SVDF.actual = False
        self.run_svd_filt()


    def create_residuum_table(self):
        
        self.tab_widget_res = QWidget(self)
        self.main_tab.insertTab(3,self.tab_widget_res, 'Residuum')
        self.tab_widget_res.setSizePolicy(self.Expand)
        self.verticalLayout_res = QVBoxLayout(self.tab_widget_res)
        self.horizontalLayout_res = QHBoxLayout()
 

        self.verticalLayout_res.addLayout(self.horizontalLayout_res)

        self.fig_res = plt.Figure((6.0, 5.0), dpi=self.dpi)
        self.canvas_res = FigureCanvas(self.fig_res)
        self.canvas_res.setParent(self.cWidget)
                
        self.fig_res.patch.set_facecolor((self.bg_col.red()/255., 
                                          self.bg_col.green()/255.,
                                          self.bg_col.blue()/255.))

        vbox = QVBoxLayout()
        vbox.addWidget(self.canvas_res)
        
        self.tomo.SVDF.fig_retro = self.fig_res
        self.tomo.SVDF.actual = False

        self.mpl_toolbar_res = NavigationToolbar2QT(self.canvas_res, self.cWidget)
        vbox.addWidget(self.mpl_toolbar_res)

        self.horizontalLayout_res.addLayout(vbox)


    
    def create_status_bar(self):
        self.status_text = QLabel("")
        self.statusBar().addWidget(self.status_text, 1)
     
    def switch_tab(self, tab):
    
        if tab == 0:
            self.status_text.setText("Remove corrupted channels by a double-click")
            
        if tab == 1:
            self.run_svd_filt()
            self.status_text.setText("f0: %.4gkHz,\t RSS/TSS: %.3g%%"%(self.tomo.SVDF.f0/1000,
                                             self.tomo.SVDF.RSS/self.tomo.SVDF.TSS*100))
        if tab == 2:
            self.run_svd_filt()
            self.status_text.setText("Mouse: center button - select channel, right - full view")  
            
            
    def closeEvent(self, event=None):
        #set the wrong dets back to GUI and recalculate tomography
        self.tomo.calculate_tomo()
        self.close() 
 





class tokamak:
    #primitive version of the tokamak class from the tomography

    norm = 1
    allow_negative = False
    allow_divertor = False
    transform_index = 0
    #rgmin = 1e-8
    camera = False
    beta = 0
    boundary = 0
    weight_power = 1
    use_pfm = False
    smooth_boundary = True
    def __init__(self, rhop,magr, magz,xgrid,ygrid):
        self.dx = np.mean(np.diff(xgrid))
        self.dy = np.mean(np.diff(ygrid))
        self.xgrid = xgrid
        self.ygrid = ygrid
        self.nx = len(xgrid)
        self.ny = len(ygrid)
        self.npix = self.nx*self.ny
        self.xmin = xgrid[0]
        self.ymin = ygrid[0]
        self.xmax = self.dx+xgrid[-1]
        self.ymax = self.dy+ygrid[-1]
        self.rhop = rhop
        self.magr = magr
        self.magz = magz

        
    def get_boundary(self,*var, **kwarg):
        #expand boundary by 50% to improve fitting 
        R0 = self.magr[:,0].mean()
        Z0 = self.magz[:,0].mean()
        extrapol = 1. if self.boundary else 1.2
        R = (self.magr[::-1,-1]-R0)*extrapol+R0
        Z = (self.magz[::-1,-1]-Z0)*extrapol+Z0
        
        #R = (self.magr[::-1,-1]-R0)*1.+R0
        #Z = (self.magz[::-1,-1]-Z0)*1.+Z0
        
        

        return np.array((R,Z)).T
    
    def mag_equilibrium(self, tvec, return_mean):
        if return_mean:
            return self.rhop,self.magr,self.magz
        return self.rhop,self.magr[:,:,None],self.magz[:,:,None]


def flux_mean(G, Tok,magr,magz):
    #flux surface averadge G profile

    n_mag = np.size(magr,1)
    
    scaling = np.array([Tok.dx,Tok.dy])
    offset = np.array([Tok.xmin,Tok.ymin])    
    coords = np.c_[magr.ravel(),magz.ravel()].T
    idx = (coords-offset[:,None])/ scaling[:,None]
    map_prof = map_coordinates(G.T,idx,order=1)
    map_prof = map_prof.reshape(magz.shape)
    return np.mean(map_prof,0)



def build_reg_mat_time(rho_mat,theta_star_rz,theta_star,rhop,magr,magz,xgrid,ygrid,BdMat, dtheta):
    #prepare regularization matrix for the solid body rotation
    
    BdMat = (BdMat| (theta_star_rz == 0)).flatten('F')
    theta_star_grid = np.linspace(0,2*np.pi, len(theta_star))
    magr = [np.interp(theta_star_grid, t,r) for r,t in zip(magr.T,theta_star.T)]
    magz = [np.interp(theta_star_grid, t,z) for z,t in zip(magz.T,theta_star.T)]

    RHO = rho_mat.flatten('F')[~BdMat]
    RHO = np.minimum(RHO, rhop[-1])
    THETA = theta_star_rz.flatten('F')[~BdMat]

    i_theta_r = ((THETA-dtheta)/(2*np.pi))%1*(len(theta_star_grid)-1)
    i_rho = RHO/rhop[-1]*(len(rhop)-1)

    xi_new_r = map_coordinates(magr,(i_rho, i_theta_r),order=3)
    yi_new_r = map_coordinates(magz,(i_rho, i_theta_r),order=3)

    ny,nx = rho_mat.shape
    dx = np.mean(np.diff(xgrid))
    dy = np.mean(np.diff(ygrid))
    npix = nx*ny
    xgrid+= dx/nx*2 #BUG??
    npoints = np.sum(~BdMat)

    x_int_r = (xi_new_r-xgrid[0])/dx
    y_int_r = (yi_new_r-ygrid[0])/dy

    x_int_c,y_int_c = np.meshgrid(np.arange(nx),np.arange(ny))
    x_int_c = x_int_c.flatten('F')[~BdMat]
    y_int_c = y_int_c.flatten('F')[~BdMat]

    xcoords = np.zeros((5,npoints),dtype=int)
    ycoords = np.zeros((5,npoints),dtype=int)
    weights = np.zeros((5,npoints))


    #forward
    xcoords[0] = x_int_r
    ycoords[0] = y_int_r
    weights[0] = (1-x_int_r+xcoords[0])*(1-y_int_r+ycoords[0])

    xcoords[1] = xcoords[0]+1
    ycoords[1] = ycoords[0]
    weights[1] = (  x_int_r-xcoords[0])*(1-y_int_r+ycoords[0])

    xcoords[2] = xcoords[0]
    ycoords[2] = ycoords[0]+1
    weights[2] = (1-x_int_r+xcoords[0])*(  y_int_r-ycoords[0])

    xcoords[3] = xcoords[1]
    ycoords[3] = ycoords[2]
    weights[3] = (  x_int_r-xcoords[0])*(  y_int_r-ycoords[0])

    #diagonal 
    xcoords[4] = x_int_c
    ycoords[4] = y_int_c
    weights[4] = 1

    col_coords = ycoords+ny*xcoords  
    row_coord = np.tile(np.arange(npix)[~BdMat],(9,1))

    #make it sparser
    weights[weights<0.1] = 0
    weights[:4]/= np.sum(weights[:4],0)

    #higher harmonics should vanish exactly in the core!!
    weights[4] += np.maximum(0.05,1-RHO)**30
    
    
    #left position
    L = sparse.coo_matrix(( weights[:4].flatten('C'), 
                         (row_coord[:4].flatten('C'),
                         col_coords[:4].flatten('C'))), 
                    shape=( npix, npix))
    #central 
    C = sparse.coo_matrix(( weights[4], (row_coord[4], col_coords[4])), 
                    shape=( npix, npix))
    
    #forward difference
    return C-L
    



def w_i(g,D):
    w = 1./(1.+np.exp(g)/D**2)
    w[~np.isfinite(w)] = 0
    return w


def GCV(g, prod,D,U,resid=0):
    #generalized crossvalidation
    w = w_i(g,D)
    return (np.sum(np.abs(((w-1)*prod))**2)+resid)/len(w)/(1-np.mean(w))**2
    
def Press(g, prod,D,U,resid=0):
    w = w_i(g,D)
    u = np.array(U,copy=False)
    ind = np.einsum('ij,ij->i',u,np.conj(u)).real < np.single(u.shape[1]*.1/u.shape[0])  #  diag(U*U.H) - smaller are wrong, linearly dependent? 
    if ind.any(): u = u[~ind]
    return np.sum(np.abs(np.dot(u, (1-w)*prod)/np.einsum('ij,ij,j->i', u,np.conj(u), 1-w))**2)/len(prod)




def FindMin(F, x0,dx0,prod,resid,D,U,tol=0.01):
    #stupid but robust minimum searching algorithm.

    fg = F(x0, prod,D,U,resid)

    while abs(dx0) > tol:
        fg2 = F(x0+dx0, prod,D,U,resid)
                            
        if fg2 < fg:
            fg = fg2
            x0 += dx0
            continue
        else:
            dx0/=-2.
            
    return x0, fg2





####          separate tomographic inversion of each harmonics               ###  
class HarmSolver(Process):

    """fast method for the inversion of each harmonics"""
    def __init__(self,qin,qout,n, E,W, Ht,Hper,Hpar,T,BdMat,b,berr,rank):
        self.qin=qin
        self.qout=qout

        self.n = n
        self.E = E
        self.Ht = Ht
        self.Hper = Hper
        self.Hpar = Hpar
        self.BdMat = BdMat.flatten('F')
        self.b = b/berr
        
        iErr = sparse.spdiags(1/berr, 0, len(b), len(b))

        self.T = iErr*T
        self.W = W
        self.rank = rank
        
        super(HarmSolver, self).__init__()


    def guess_lam(self,U,S,d):
        g0 = 2*np.log(np.median(S))
        prod = np.asarray(np.dot(U.T,self.b[~d]))
        g, f = FindMin(Press,g0,1,prod,0,S,U)
        return g
   
        
    def run(self):

        H = self.Hper +  self.Hpar  
        if self.n > 0:
            H = H + self.Ht.H*self.Ht#TODO set ration between them?
            H = self.W*H*self.W.T  #compansate for asymmetry of 0. order profile

        wrong_dets = np.squeeze(np.array(self.T.sum(1)==0))
        K = self.T[np.where(~wrong_dets)[0]]  

        npix = K.shape[1]

        ind = np.where(~self.BdMat)[0] 
        B = H[:,ind][ind,:]  #real regularisation matrix
        K = K[:,ind]   #complex contribution matrix

        if self.n > 0:
            En = self.E.tocsc()[:,ind][ind,:]**self.n
            K *= En  #include a phaase shift in the contribution matrix

        try:
            F = cholesky(B, ordering_method='best') #167ms
        except Exception as e:
            print( 'Cholmod error: ', e)
            print('Check the input data, remove corrupted channels!!!')
            #try to add regularization
            F = cholesky(B, ordering_method='best',beta=1) 
        
        #first inversion of the geometry matrix - slowest step?
        LPK = F.solve_L(F.apply_P(deepcopy(K.H.real)))
        if self.n > 0:#add complex parts, cholmod do not allow to applu solve_L on complex matrix if L is real!
            LPK += 1j*F.solve_L(F.apply_P(deepcopy(K.H.imag)))


        LPK = LPK.toarray().T
        sqrtD = np.sqrt(F.D().real)
        LPK/= sqrtD

        LL = np.dot(LPK, np.conj(LPK).T)  #Gram matrix
        s,u = eigh(LL,overwrite_a=True, check_finite=False,lower=True)  
        s,u = np.sqrt(np.maximum(s,0))[::-1], u[:,::-1] 

        try:
            rank = np.where(np.cumprod((np.diff(np.log(s[s!=0]))>-5)|(s[s!=0][1:]>np.median(s))))[0][-1]+2
        except:
            print('rank s', s)
        rank_ = min(self.rank,rank)
        S = s[:rank_]
        U = u[:,:rank_]
        

        while True:
            #wait for a value of lambda and return results

            lam = self.qin.get()
            
            if lam is None:
                return
            elif lam <= 0.05:
                g = self.guess_lam(U,S,wrong_dets)
            else:
                g = mquantiles(2*np.log(S),lam)
            t2 = time.time()

            #evaluate solution
            w = w_i(g,S)
            prod = np.dot(U.T,self.b[~wrong_dets])
            f = prod*(w/S)
            #G_ solution before the last inversion 
            G_ = np.dot(np.dot(f,np.conj(U.T/S[:,None])), LPK)
            G_ /= sqrtD
            
            #final inversion
            G = np.zeros(npix,dtype=LPK.dtype)
            if self.n > 0:
                G[~self.BdMat] = En*(F.apply_Pt(F.solve_Lt(G_.real))+1j*F.apply_Pt(F.solve_Lt(G_.imag)))
            else:
                G[~self.BdMat] = F.apply_Pt(F.solve_Lt(G_))

            retro = np.zeros_like(self.b)
            retro[~wrong_dets] = np.dot(np.conj(U),f*S)
            chi2 = (np.linalg.norm(((w-1)*prod))**2)/np.sum(~wrong_dets)

            lam = np.interp(g, 2*np.log(S)[::-1], np.linspace(0,1,len(S)))

            self.qout.put((G,retro,lam,chi2 ))


def OptimizeF0( tvec, sig, f0,df0 = 200,n_steps=400):
    #simple algorithm to fine tune basic frequency
    
    from matplotlib.mlab import detrend_linear
    sig = detrend_linear(sig  )
    difference = np.ones(n_steps)*np.infty
    
    test_fun = np.exp(1j*2*np.pi*(f0-df0)*tvec)
    test_fun/= np.linalg.norm(test_fun)
    dtest_fun = np.exp(1j*2*np.pi*(2*df0/float(n_steps))*tvec)
    
    for i in range(n_steps):
        
        retro = np.outer(test_fun, np.dot(np.conj(test_fun), sig))[:,0]*2
        difference[i] = np.linalg.norm(retro.real-sig)
        test_fun*= dtest_fun

    return f0+(np.argmin(difference)*2.-n_steps+1)/(n_steps)*df0

 

 
def create_derivation_matrix(g, Bmat, danis,rgmin=1e-8):
    """
    Prepare matrix of derivation according to choosed regularization method
    """
 
    g_tmp = np.mean(np.asarray(g),1)
    
    npix = len(g_tmp)
    max_g = np.nanmax(abs(g_tmp))

    g_tmp /= max_g
    

    g_tmp[g_tmp < rgmin] = rgmin
    

    g_tmp **=-1

 
    W=sparse.spdiags(g_tmp,0, npix,npix,format='csr')

    #Bmat = [Bper1_frwd,Bpar1_frwd,Bper2_bckwrd,Bpar2_bckwrd]
    def sigmoid(x):    return 1 / (1 + np.exp(-x))

    Hper = sigmoid(-danis)*(Bmat[0].T*(W*Bmat[0]) + Bmat[2].T*(W*Bmat[2]))
    Hpar = sigmoid( danis)*(Bmat[1].T*(W*Bmat[1]) + Bmat[3].T*(W*Bmat[3]))

    
        
    return Hper.tocsc(), Hpar.tocsc()



class Roto_tomo:
    
    keyCtrl = False
    keyShift = False
    initialized = False
    loaded = False

    def __init__(self, parent, fig,m_num,n_num,add_back, TeOver,show_flux,
                 plot_limit, reg, n_harm,n_svd,map_equ,rho_lbl, show_contours ):

        
        self.m_num = m_num
        self.n_num = n_num
        self.TeOver = TeOver
        self.show_flux = show_flux
        self.slider_lim = plot_limit
        self.slider_reg = reg
        self.add_back = add_back
        self.eqm = map_equ
        self.rho_lbl = rho_lbl
        self.show_contours = show_contours
        self.t_range = 0, np.inf
        self.f_range = 0, np.inf

        
        self.parent = parent
        self.fig = fig
        self.fig.subplots_adjust(right=0.8,top=.95)
        self.ax = self.fig.add_subplot(111)
        self.substract = not self.add_back.isChecked()
        self.dPhi = 0

        self.n_harm = n_harm
        self.n_svd = n_svd
        self.tau = 4 #time resolution for signal extraction
        self.nr = 100 #radial resulution for calculation of magnetic equalibrium in radial/ poloidal coordinates
        self.ntheta = 120 #angular resulution for calculation of magnetic equalibrium in radial/ poloidal coordinates
        self.dtheta =  2*np.pi/30 #poloidal step for regularisation matrix
        self.rgmin = 1e-8 #zero zupression parameter of the tomography
        self.nfisher = 3 #number of fisher iterations
        self.plot_lim = self.slider_lim.value()/100.
        self.danis = 4  #ratio between poloidal and radialcorelation in tomgraphy
        self.nx = 80    #horisontal resolution of the tomography
        self.ny = 120   #vertical resolution of the tomography
        self.m = int(self.m_num.currentText())
        self.n = int(self.n_num.currentText())
        self.n_rho_plot = 10  #number of plotted radial magnetic contours
        self.n_theta_plot = 16#number of plotted poloidal magnetic contours
        self.n_contour = 15 #numberof Te contours 
        self.fontsize = 10
        self.hsolvers = []
        self.t0 = np.nan
        self.showTe = False
        self.cmap = my_cmap
 
   
        if self.show_contours:
            self.levels = np.linspace(0, 1, self.n_contour)            
            self.tomo_img = self.ax.contourf([0,0], [0,0], np.ones((2,2)),levels=self.levels,
                                             extend='both', cmap= self.cmap, vmin=0, vmax = 1) 

        else:
            self.tomo_img = self.ax.imshow(np.zeros((self.nx,self.ny)), animated=True,origin='lower'
                                        ,vmin = 0,vmax=1,cmap=self.cmap, interpolation='quadric')
        #magnetic flux surfaces
        X = np.zeros((1,self.n_rho_plot))
        self.plot_mag_rho   = self.ax.plot(X,X,'--',c='0.5',lw=.5, visible=self.show_flux.isChecked())
        X = np.zeros((1,self.n_theta_plot))
        self.plot_mag_theta = self.ax.plot(X,X,'--',c='0.5',lw=.5, visible=self.show_flux.isChecked())

        #tokamak chamber
        logger.debug('tok lbl: %s', self.parent.tokamak)
        if self.parent.tokamak == 'AUG':
            gc_d = sf.getgc()
            for gcc in gc_d.values():
                self.ax.plot(gcc.r, gcc.z, c='0.5', lw=.5)
        elif self.parent.tokamak == 'DIIID':
            from loaders_DIIID import map_equ
            gc_r, gc_z = map_equ.get_gc()
            for key in gc_r.keys():
               self.ax.plot(gc_r[key], gc_z[key], '0.5', lw=.5)

        for label in (self.ax.get_xticklabels() + self.ax.get_yticklabels()):
            label.set_fontsize( self.fontsize) # Size here overrides font_prop

        self.cbar_ax = self.fig.add_axes([0.85, 0.1, 0.04, 0.85])
        self.cbar_ax.xaxis.set_major_formatter(plt.NullFormatter())
        self.cbar_ax.yaxis.set_major_formatter(plt.NullFormatter())
        self.cbar_ax.tick_params(labelsize=  self.fontsize) 
        
        self.update_colorbar()
        #self.cbar = self.fig.colorbar(self.tomo_img, cax=self.cbar_ax)
        self.ax.set_title('SXR emissivity [kW/m$^3$]',fontsize= self.fontsize)
        self.plot_description = self.ax.text(1.008,.05,'',rotation='vertical', 
                transform=self.ax.transAxes,verticalalignment='bottom',
                size='xx-small',backgroundcolor='none')
        
        self.ax.axis('equal')
        self.ax.set_xlabel('R [m]',fontsize= self.fontsize)
        self.ax.set_ylabel('z [m]',fontsize= self.fontsize)
        self.ax.yaxis.labelpad = -5
        
        self.ax.patch.set_facecolor(self.cmap(0))


        self.ax.axis([1.4,2.2,-0.25,0.35]) 
        self.shift_phi = 0
        self.dR = 0
        self.dZ = 0

        self.cid_scroll  = self.fig.canvas.mpl_connect('scroll_event',self.MouseWheelInteraction)
        self.cid_press   = self.fig.canvas.mpl_connect('key_press_event',   self.onKeyPress)
        self.cid_release = self.fig.canvas.mpl_connect('key_release_event', self.onKeyRelease)

  


    def prepare_tok_object(self, tok_lbl, shot):

       

        #Prepare tokamak object from original tomography code
        input_parameters = read_config(tomo_code_path+"tomography_D3D.cfg")
        input_parameters['shot'] = shot
        input_parameters['local_path'] = tomo_local_path
        input_parameters['program_path'] = tomo_code_path
        input_parameters['nx'] = self.nx
        input_parameters['ny'] = self.ny

        if not hasattr(config, 'wrong_dets_pref'):
            config.wrong_dets_pref = input_parameters['wrong_dets']
        

        if tok_lbl == 'DIIID':
            import geometry.DIIID as Tok
            #if self.sxr_harmonics is None:
            diag = 'SXR fast'
 
            diag_path = tomo_local_path+ 'geometry/DIIID/SXR'
            
            #$use only two poloidal cameras 
            config.wrong_dets_pref = np.unique(list(config.wrong_dets_pref)+list(range(64,88)))
            
        elif tok_lbl == 'AUG':
            import geometry.ASDEX as Tok
            diag = 'SXR_fast'
            diag_path = tomo_local_path+ 'geometry/ASDEX/SXR'
        else:
            raise Exception('Support of the tokamak %s was not implemented' %tok_lbl)

        logger.debug('diag_path %s', diag_path)
        if not os.path.exists(diag_path):
            os.makedirs(diag_path)
        
         
        try:
            self.tok = Tok(diag, input_parameters, load_data_only=True, only_prepare=True)#BUG 
        except:
            import traceback
            print( traceback.format_exc())
            QMessageBox.warning(self.parent,"Load SXR data issue", "Fast SXR data probably do not exist",QMessageBox.Ok)
            return 

        #calculate complex geometry matrix including all informations about diagnostic            
        self.T_full = geom_mat_setting( self.tok,self.nx, self.ny, 
                             self.tok.virt_chord, path=None)[0]
       

        
  
        
    def prepare_tomo(self, tok_lbl, shot ,t_range, f_range, t0, f0,  cross_sig_name = None, 
                              sxr_harmonics=None, cmplx_cross_sig=None):

        QApplication.setOverrideCursor(Qt.WaitCursor)

        t = time.time()

      
        self.t_range = t_range
        self.f_range = f_range
        self.t0 = t0  #time where exactly sxr_harmonics are evaluated
        self.F0 = self.F1 = f0  #median frequency in the t_range
        self.sxr_harmonics = sxr_harmonics
        self.cmplx_cross_sig = cmplx_cross_sig
        self.cross_sig_name  = cross_sig_name
        self.shot = shot
        self.tok_lbl = tok_lbl
  

        self.rhop = np.linspace(0,1,self.nr)
        theta_in  = np.linspace(0,2*np.pi,self.ntheta)

        if not hasattr(self, 'tok'):
            self.prepare_tok_object(tok_lbl, shot)


        if self.tok_lbl == 'AUG':
            magr, magz = sf.rhoTheta2rz(self.eqm, self.rhop,theta_in,t0,coord_in=self.rho_lbl)
            self.rho = sf.rho2rho(self.eqm, self.rhop, t0, coord_in='rho_pol', coord_out=self.rho_lbl)[0]
        elif self.tok_lbl == 'DIIID':
            magr, magz = self.eqm.rhoTheta2rz(self.rhop,theta_in,t0,coord_in=self.rho_lbl)
            self.rho = self.eqm.rho2rho(self.rhop, t0, coord_in='rho_pol', coord_out=self.rho_lbl)[0]
        self.magr = magr[0]
        self.magz = magz[0]
        
        #set view to contain rho = 0.7 contour 
        ir = np.argmin(np.abs(self.rhop - 0.5))
        self.ax.axis([self.magr[:,ir].min(),self.magr[:,ir].max(),self.magz[:,ir].min(),self.magz[:,ir].max()])

        if hasattr(self.eqm, 'getQuantity'): #different eqm_map version on DIII-D
            self.q_prof = self.eqm.getQuantity(self.rhop, 'Qpsi', t_in=t0)	
        else:#  and AUG
            jtq = np.argmin(abs(eqm.time - t0))
            nrho = np.max(self.eqm.lpfp) + 1
            psi = self.eqm.pfl[jtq, :nrho]
            q   = self.eqm.q[jtq, :nrho]
            rhop_q = sf.rho2rho(self.eqm, psi, t_in=t0, coord_in='Psi', coord_out='rho_pol')[0]
            self.q_prof = np.interp(self.rhop, rhop_q, q)
 
        
        
        
        self.xgridc = (self.tok.xgrid+self.tok.dx/2)/self.tok.norm  #centers of the pixels
        self.ygridc = (self.tok.ygrid+self.tok.dy/2)/self.tok.norm  #centers of the pixels
        self.nl = self.tok.nl
        
        if not self.show_contours:
            self.tomo_img.set_extent([self.tok.xgrid[0],self.tok.xgrid[-1]+self.tok.dx,
                                      self.tok.ygrid[0],self.tok.ygrid[-1]+self.tok.dy])

        
        if self.tok_lbl == 'DIIID':
            config.useCache = False  #force not using cached SXR geometry
 
        


        #Auxiliarly tokamak to avoid loading of equilibrium from tomography
        self.aux_tok = tokamak(self.rhop,self.magr, self.magz,self.tok.xgrid,self.tok.ygrid)
        self.aux_tok.vessel_boundary = self.tok.vessel_boundary

        #Load diagnostic infomation
        n_harm_max = 10
        self.A   = np.ones((self.nl,n_harm_max-1)) #assume n_harm_max harmonic at most!
        self.phi = np.zeros((self.nl,n_harm_max-1))
        self.Phi = np.asarray(self.tok.Phi)  #BUG co ten posuv o 45 stupnu? 
        self.Phi0 = np.median(self.Phi) #toroidal position of SXR cameras in radians!

        #add correction for finite DAS IIR
        if self.tok_lbl == 'AUG':

            slow_das = np.loadtxt(loc_dir+'/slow_sxr_diag.txt' )
            fast_das = np.loadtxt(loc_dir+'/fast_sxr_das.txt'  )
            SampFreq =  self.tok.SampFreq 
            #amplitude reduction of DAS IIR
            self.A[  SampFreq == 5e5] = np.interp(self.F0*np.arange(1,n_harm_max),slow_das[:,0],slow_das[:,1])
            self.A[  SampFreq == 2e6] = np.interp(self.F0*np.arange(1,n_harm_max),fast_das[:,0],fast_das[:,1])
            #phase shift of DAS IIR
            self.phi[SampFreq == 5e5] = np.interp(self.F0*np.arange(1,n_harm_max),slow_das[:,0],slow_das[:,2])
            self.phi[SampFreq == 2e6] = np.interp(self.F0*np.arange(1,n_harm_max),fast_das[:,0],fast_das[:,2])
        if self.tok_lbl == 'DIIID':
            pass
            #response of unknown... :( 

        #calculate equilibrium related properties
        theta_star_rz, theta_star = self.tok.mag_theta_star(t0,self.rhop,self.magr,self.magz,rz_grid=True)
        self.theta_star_rz, theta_star = theta_star_rz.T, theta_star.T
        self.rho_mat = get_rho_field_mat(self.aux_tok, t0)
        self.BdMat = get_bd_mat(self.aux_tok, time=t0).reshape(self.ny,self.nx, order='F')

        #self.BdMat|= self.theta_star_rz==0 
        #remove points outside of boundary
        self.T_full_prunned = prune_geo_matrix(self.T_full,self.BdMat)
 
        #regularization matrix
        self.Bmat, diag_mat = mat_deriv_B(self.aux_tok, 0, 1,None)
        self.Ht = build_reg_mat_time(self.rho_mat,self.theta_star_rz,theta_star,self.rhop,self.magr,
                                self.magz,self.xgridc,self.ygridc,self.BdMat,self.dtheta)

        #contours of constant theta star
        t = np.linspace(0,2*np.pi,self.n_theta_plot,endpoint=False)
        self.isotheta_R = np.array([np.interp(t,ts,r) for ts,r in zip(theta_star.T,self.magr.T)]).T
        self.isotheta_Z = np.array([np.interp(t,ts,z) for ts,z in zip(theta_star.T,self.magz.T)]).T

        #contours of constant coordinate rho_lbl
        r = np.linspace(0,self.rho[-1],self.n_rho_plot)
        self.isoflux_R = interp1d(self.rho,self.magr)(r).T
        self.isoflux_Z = interp1d(self.rho,self.magz)(r).T
        
              
        #update plots of magnetics map_coordinates 
        for ip,p in enumerate(self.plot_mag_theta):
            p.set_data(self.isotheta_R[ip],self.isotheta_Z[ip])
        
        for ip,p in enumerate(self.plot_mag_rho):
            p.set_data(self.isoflux_R[ip],self.isoflux_Z[ip])
        
      
        
        self.load_data()

        self.initialized = True


        self.calculate_tomo()
        
        #update also Te if already shown
        self.TeOverplot()

        QApplication.restoreOverrideCursor()


    def calculate_tomo(self):
        tcalc = time.time()
        ##close threats from previous calculation
        if not self.initialized:
            return
            
        QApplication.setOverrideCursor(Qt.WaitCursor)

        for h in self.hsolvers:
            h.qin.put(None) #send kill signal
            h.terminate()
            h.join(1e-2)  #join process
            h.qin.close()
            h.qout.close()
            del h  #release memory 
            
        self.hsolvers = []
        QApplication.restoreOverrideCursor()
         
        #convert time dependent signals in complex harmonics 
        self.prepare_harms()
        #prepare inversion
        self.precalc_tomo()
        #use m and n number from the GUI and evaluate tomography 
        self.update_mode_number()
        
        print( 'Tomographic inversion calculated in %.1fs'%(time.time()-tcalc))

        
    def prepare_harms(self):
        #calculate complex harmonics out of the measured signals

        QApplication.setOverrideCursor(Qt.WaitCursor)
        if self.sxr_harmonics is not None:
 
            ind_correct = self.tok.get_correct_dets()

            self.dets = np.where(ind_correct)[0]
            self.bb = [bb[ind_correct] for bb in self.all_bb[:self.nharm]]
            self.bb_err = [be[ind_correct] for be in self.all_bb_err[:self.nharm]]
            
        else:
            #apply SVD filter to get complex harmonic 
            self.SVDF.run_filter(update_plots=False)
            
            self.F1 = self.SVDF.f0
            self.t0 = self.SVDF.t0
            self.bb = list(self.SVDF.harm)
            self.bb[0] = self.bb[0].real
            self.bb_err = [np.copy(self.SVDF.harm_err)+1e-6 for h in self.bb]
            
            #NOTE errorbars found this ways are too small for zeroth and first harmonics
            self.bb_err[0] += self.SVDF.err[~self.SVDF.invalid]

            QApplication.restoreOverrideCursor()
 
            #add at least 5% noise to zero'th harmonics    
            self.bb_err[0] += np.abs(self.bb[0])*0.05+0.03*np.abs(self.bb[0]).mean()


        QApplication.restoreOverrideCursor()


    def load_data(self):


        if self.sxr_harmonics is not None:
            
            #SXR data were loaded directly by PYSPECVIEW, works only for DIII-D now

            bb = []
            bb_err = []
            dets = []
            idet = 0
            #calibration factors for each camera estimated by PYTOMO 
            camera_calib = self.tok.get_calb()
            #sort signals n the order expected by the geometry matrix 
            for c, (cam, channels) in zip(camera_calib, self.tok.detectors_dict.items()):
                for ch in channels:
                    ch = int(ch.split('_')[1])
                    hdata = self.sxr_harmonics[(cam[:5], ch)]
                    bb.append([h*c for h in hdata['harm']])
                    bb_err.append([c*hdata['error']]*len(bb[-1]))
                    if hdata['valid'] and idet in self.tok.dets:
                        dets.append(idet)
                    idet += 1 
            self.all_bb = [np.array(bb) for bb in zip(*bb)]
            self.all_bb_err = [np.array(be)+1e-6 for be in zip(*bb_err)]

            #add at least 5% noise to 0th an 1th harmonics    
            self.all_bb_err[0] += np.abs(self.all_bb[0])*0.05+0.03*np.abs(self.all_bb[0]).mean()
            self.all_bb_err[1] += np.abs(self.all_bb[1])*0.05+0.03*np.abs(self.all_bb[1]).mean()
            self.tok.dets = np.array(dets)

        else:
            print( 'Loading SXR data ... ')

            QApplication.setOverrideCursor(Qt.WaitCursor)
            config.useCache = True  #force storing of the data

            #load SXR data 
            t =  time.time()
            signals,error_sig, self.tvec, dets,_,_ = self.tok.prepare_data(self.t_range[0],self.t_range[1],1,1,1,detsCutOff=False)

            
            self.signals,self.error_sig = signals.T,error_sig.T

        
            self.dets = dets[np.all(np.isfinite(error_sig),1)[dets]]
    
            #initialise it only once 
            dets_index = self.tok.dets_index[:-1]  
            dF = np.diff(self.f_range)[0]/2
            self.SVDF = SVDFilter(self.tvec, self.signals,np.nanmean(self.error_sig ,0),dets_index,
                                self.F0, dF, self.n_harm,self.n_svd, self.tau, self.cmplx_cross_sig)
            invalid = ~np.in1d(np.arange(self.signals.shape[1]), self.dets)
            self.SVDF.set_corrupted_ch(invalid)       
               

   

    def precalc_tomo(self):
       
        #first iterations of the tomography algorithm
        QApplication.setOverrideCursor(Qt.WaitCursor)

        #initial guess of the SXR profile
        P = (self.rhop**2+4**-2)**(-2)*np.tanh((1-self.rhop)*10)
        P = np.exp(-((self.rhop)/0.7)**2)
        G0 = np.interp(self.rho_mat, self.rhop, P/P[0])
        G0[self.BdMat] = 0
   
        #use only valid detectors
        self.T = self.T_full_prunned[self.dets]
     
        #solution of tomography for stationary profile 
        t = time.time()

        for ifisher in range(self.nfisher):
            Hper,Hpar = create_derivation_matrix(G0.reshape(-1,1,order='F'), self.Bmat,self.danis)
            qin, qout = Queue(),Queue() 

            if ifisher < self.nfisher-1:  #the last step will be evaluated with the others  #BUG slow!!! 
                h = HarmSolver(qin, qout,0,None,None,None,Hper,Hpar,self.T,
                               self.BdMat,self.bb[0],self.bb_err[0],len(self.dets))
                h.start()
                qin.put(0.7)  #BUG Regularzation is fixed to 0.7 to avoid overfitting
                G0,retro,lam,chi2 = qout.get()
                qin.put(None) #kill process
                qin.close()
                qout.close()
                h.terminate()
                h.join(1e-2)
                
                
    
        G0 = G0.reshape(self.ny,self.nx,order='F')
        #use asymmetry of . 0. harmonic in the regularization operator
        lim = G0.max()/50
        G0 = np.hypot(G0, lim)
        fm = flux_mean(G0, self.aux_tok,self.magr,self.magz)
        G0/= np.interp(self.rho_mat,self.rhop, np.hypot(fm, lim))
    
        #weight matrix compensating poloidal asymmetry 
        npix = self.nx*self.ny
        self.W = sparse.spdiags( 1/G0.flatten('F'),0, npix, npix, format='csc')
        self.Hper = Hper
        self.Hpar = Hpar

        QApplication.restoreOverrideCursor()

         
        
    def update_node_m_number(self,ind):
        M = self.parent.m_numbers[ind]
        if M is self.m: return 
        self.m = M
        self.m_num.setCurrentIndex(ind)
        self.update_mode_number()
        #update also ECE mapping if availible
        if self.showTe:
            self.Te2Dmap.clear_contours(ax=self.ax)
            self.Te2Dmap.UpdateModeM(ind)
        
        
    def update_node_n_number(self,ind):
        N = self.parent.n_numbers[ind]
        if N is self.n: return 
        self.n = N
        self.n_num.setCurrentIndex(ind)
        self.update_mode_number()

    def update_mode_number(self):
        self.set_mode_numbes()
        self.lam0 = -1 #guess regularization in the first step 
        self.eval_tomo()
        self.lam0 = int(100*self.gamma[1])/100.
        self.slider_reg.setValue(int(100*self.lam0))
        self.update(update_cax=True)
        
    def set_plot_lim(self, val):
        if self.plot_lim == val/100.: return 
        self.plot_lim = val/100.
        self.update(update_cax=True)

    def set_reg(self, val):
        if self.lam0*100 == val: return

        self.lam0 = val/100.
        self.eval_tomo()
        self.update(update_cax=True)

    def set_substract(self,ind):
        self.substract = not self.add_back.isChecked()
        if self.substract:
            self.cmap = cm.get_cmap('seismic')
            c = 0.5
        else:
            self.cmap = my_cmap
            c = 0
        
        self.tomo_img.set_cmap(self.cmap)
        #set background color of ax
        self.ax.patch.set_facecolor(self.cmap(c))   
        self.update(update_cax=True)
        
    def update_colorbar(self):
        #BUG update colorbar by creating a new one :( 
        self.cbar_ax.cla()
        cb = self.fig.colorbar(self.tomo_img, cax=self.cbar_ax)
        tick_locator = MaxNLocator(nbins=7)
        cb.locator = tick_locator
        cb.update_ticks()     
                    
    def TeOverplot(self):
        #show contours of Te
        QApplication.setOverrideCursor(Qt.WaitCursor)

        self.showTe = self.TeOver.isChecked()
        if not self.showTe:
            if hasattr(self, 'Te2Dmap'):
                self.Te2Dmap.clear_contours(self.ax)
                self.fig.canvas.draw_idle()
            QApplication.restoreOverrideCursor()
            return

        
        #load ECE data, prepare 2D profile
        if not self.parent.Te2Dmap.initialized:
            self.parent.radial_view.t_range = self.t_range
            self.parent.radial_view.f_range = self.f_range

            id_2D = self.parent.tables_names.index('2D Te')
            self.parent.updatePanels(id_2D)
        
        self.Te2Dmap = self.parent.Te2Dmap
        self.Te2Dmap.m = 0
        self.Te2Dmap.lim = 0
        self.Te2Dmap.f0 = self.F1
        self.Te2Dmap.phase_locked = False

        self.Te2Dmap.UpdateModeM(np.searchsorted(self.parent.m_numbers,self.m), animate=True)
        self.fig.canvas.draw_idle()
        self.update( update_mag=False)
        QApplication.restoreOverrideCursor()

        
    def show_magflx(self):
        #show stationary magnetic equilibrium
        v = self.show_flux.isChecked()
        for p in self.plot_mag_theta+self.plot_mag_rho:
            p.set_visible(v)
        self.fig.canvas.draw_idle()
        
    def show_retrofit(self):
        #show how well were fitted experimental data 
        
        
        dets_dict = self.tok.detectors_dict
        det_ind = [0,]
        for k, i in  dets_dict.items():
            det_ind.append(det_ind[-1]+len(i))
        
        fig_exists = False
        if plt.fignum_exists('Retrofit'):
            fig_exists = True
            plt.figure('Retrofit').clf()
            
        f,axes = plt.subplots(max((self.n_harm+1)//2,1),2,sharex=True,num='Retrofit',figsize=(9,8))
        f.suptitle('t = %.3fs'%self.t0)
        for ax in np.atleast_2d(axes)[:,0]:
            ax.set_ylabel('SXR [kW/m$^2$]')
            
        axes = axes.flatten()
        axes[0].set_xlim(0, det_ind[-1])
        f.subplots_adjust(hspace=0.07, wspace = 0.07)

        for n,ax in enumerate(axes):
            if n >= len(self.bb):continue 
        
            phase_data = np.angle(self.bb[n])
            aplitude_data = np.abs( self.bb[n])
            phase_retro = np.angle(self.retro[n])
            aplitude_retro = np.abs(self.retro[n])
            aplitude_err =  self.bb_err[n]

            for i in range(len(det_ind)-1):
                ind = (self.dets>= det_ind[i])&(self.dets <  det_ind[i+1])
                if np.sum(ind) < 3:  continue
                phase_retro[ind] = np.unwrap(phase_retro[ind],discont=1.5*np.pi)
                phase_retro[ind] -= np.round(np.average(phase_retro[ind], weights=aplitude_retro[ind])/(2*np.pi))*2*np.pi
                phase_data[ind] = np.unwrap(phase_data[ind]-phase_retro[ind])+phase_retro[ind]
                phase_data[ind] -=  np.round(np.average(phase_data[ind]-phase_retro[ind], weights=aplitude_data[ind])/(2*np.pi))*2*np.pi

            phase_shift =  np.average(phase_retro, weights=aplitude_retro )
            phase_data  -= phase_shift
            phase_retro -= phase_shift
            
            weak = aplitude_data < np.nanmax(aplitude_data)*0.05
            phase_data[weak] = np.nan
            phase_retro[weak] = np.nan
   
            ax.plot(self.dets, aplitude_retro/1e3, 'b',label='retro')
            ax.errorbar(self.dets,aplitude_data/1e3, aplitude_err/1e3,c='r',label='data')
            ymax = max(aplitude_data.max(),aplitude_retro.max())/1e3
            
            ax_phase = ax.twinx()

            ax.set_ylim(0,ymax)
       
            w = aplitude_data/np.nanmax(aplitude_data)
            w[weak] = 0
            if any(~np.isreal(self.bb[n])):
                ax_phase.scatter(self.dets, phase_data,facecolor=(np.outer((1,0,0), w)+1-w).T,      
                           edgecolors=(1,0,0),s=50,linewidths=.5)
                ax_phase.plot( self.dets,  phase_retro ,'.--b')
            ax.set_xlim(self.dets[0], self.dets[-1]+1)
            ax_phase.set_ylim(-2*np.pi, 2*np.pi)
          
            y_tick = np.linspace(-np.pi*2, np.pi*2, 9)
            y_label = [r"$-2\pi$", r"$-3\pi/2$", r"$-\pi$", r"$-\pi/2$", "$0$", 
                   r"$\pi/2$", r"$\pi$", r"$3\pi/2$", r"$2\pi$"]
            ax_phase.set_yticks(y_tick)
            ax_phase.set_yticklabels(y_label ) 
            ax.legend()
            [ax.axvline(i-0.5,c='k') for i in det_ind]
            ax.axhline(y=0,c='k')
            ax.set_title('%d. harmonic: $\gamma$ = %.2f  $\chi^2$=%.2f'%(n,self.gamma[n],self.chi2[n]))
            
        plt.tight_layout()
        self.multi_retrofit = MultiCursor(f.canvas, axes ,horizOn=False, color='k', lw=1)

        self.AxZoom = fconf.AxZoom()
        self.zoom_cid = f.canvas.mpl_connect('button_press_event', self.AxZoom.on_click)
        
        if fig_exists:
            f.canvas.draw()
        else:
            f.show()







        
    def set_mode_numbes(self):
   
        #do reconstruction assuming certain M and N mode number
        QApplication.setOverrideCursor(Qt.WaitCursor)

        t = time.time()
        npix = self.nx*self.ny
        ndet = len(self.dets)
        A = self.A[self.dets,:self.n_harm-1]#amplitude correction for detector response
        phi = self.phi[self.dets,:self.n_harm-1] #phase correction for detector response
        Phi = self.Phi[self.dets] #toroidal position of SXR cameras
        #build complex contribution matrix 
        orientation = np.sign(self.q_prof.mean())
        #phase and amplitude correction for different toroidal positions and response of detectors 
        cmplxA = A*np.exp(-1j*(phi+(Phi[:,None]-self.Phi0)*self.n*orientation*np.arange(1,self.n_harm)))
        #complex geometry matrix
        cmplxT = [self.T,]+[sparse.spdiags(cmplxA[:,i],0,ndet,ndet)*self.T for i in range(self.n_harm-1)]
        E = sparse.spdiags(np.exp(1j*self.theta_star_rz.flatten('F')*self.m),0, npix, npix)
        
        #close threats from previous calculation
        for h in self.hsolvers:
            h.qin.put(None) #send kill signal
            h.terminate()
            h.join(1e-2)  #join process
            h.qin.close()
            h.qout.close()
            del h  #release memory 
            
        self.hsolvers = []

        #precalculate decompositions - slow!!
        for n in range(self.n_harm):
            qin,qout = Queue(),Queue()
            h = HarmSolver(qin,qout,n,E,self.W,self.Ht,self.Hper, self.Hpar,
                           cmplxT[n],self.BdMat,self.bb[n],self.bb_err[n], ndet)
            h.start()
            self.hsolvers.append(h)


        self.shift_phase(0)
        
        QApplication.restoreOverrideCursor()


    def eval_tomo(self):
        #evaluate tomographic reconstructed profiles for specific value of the regularization 
        
        #special choise of the lambda for 0.th harmonic
        self.hsolvers[0].qin.put(self.lam0)

        #manual choise
        for h in self.hsolvers[1:]:
            h.qin.put(self.lam0) 

        #get output, usually fast
        out = [h.qout.get() for h in self.hsolvers]
        self.G     = [o[0].reshape(self.ny,self.nx,order='F') for o in out]
        self.retro = [o[1]*be for o,be in zip(out,self.bb_err)]
        self.gamma = [o[2] for o in out]
        self.chi2  = [o[3] for o in out]
        
        #estimate vmax for the tomography image
        G_t = [np.abs(g) for g in self.G]
        self.vmax = np.sum(G_t,0).max()
        self.vmax_bcg = np.sum(G_t[1:],0).max()
        
        return
    
        #rest is just for debugging

        cmplx_phase = np.exp(2*np.pi*self.F1*1j*(self.tvec-self.t0))
        

        #retrofit in time domain
        self.retro_t    = sum([np.outer(r,cmplx_phase**i).real for i,r in enumerate(self.retro)],0).T
        #filtered data in time domain
        self.filtered_t = sum([np.outer(r,cmplx_phase**i).real for i,r in enumerate(self.bb)],0).T
 
        #this is exactly equal to self.retro_t
        self.retro_gt = np.sum([np.real(self.T*g.flatten(order='F')*cmplx_phase**i) for i,g in enumerate(self.G)],0).T

        
        

    def shift_phase(self,shift_phi):
        
        if not self.initialized:
            return
        
        
        self.shift_phi %= 2*np.pi #periodicity
        self.time = self.t0 + self.shift_phi/(2*np.pi*self.F1)*np.abs(self.m)
        description = '#%d  at %.6fs, $\psi_0$ = %d, f$_0$ = %.4fkHz and m/n=%d/%d'%(self.shot,
                                    self.time,np.rad2deg(self.Phi0), self.F1/1e3, self.m,self.n)
        self.plot_description.set_text(description)


    def update(self,update_mag=True,animate=False,update_cax=False,):

        w = 2*np.pi*self.F1
        G_t = [np.real(g*np.exp(i*1j*w*(self.time-self.t0))) for i,g in enumerate(self.G)]
        G_t = np.sum(G_t[1:] if self.substract else G_t,0)
        
        if self.substract:
            vmin,vmax = -self.plot_lim*self.vmax_bcg/1e3, self.vmax_bcg*self.plot_lim/1e3 
            self.levels = np.linspace(vmin,vmax, self.n_contour*2)            
        else:
            vmin,vmax = (1-self.plot_lim)*self.vmax/1e3, self.vmax/1e3 
            self.levels = np.linspace(vmin,vmax, self.n_contour)            
       
        
        if self.show_contours:
            #remove previous contours
            if hasattr(self.tomo_img, 'collections'):
                for coll in self.tomo_img.collections:
                    try:    self.ax.collections.remove(coll)
                    except ValueError: pass#Everything is not removed for some reason!    
                
            self.tomo_img = self.ax.contourf( self.xgridc, self.ygridc, G_t/1e3, extend='both', 
                                    vmin=vmin, vmax=vmax, levels=self.levels, cmap=self.cmap)
            if update_cax:
                self.update_colorbar()
            anim_obj = (self.tomo_img.collections,)

        else:
            self.tomo_img.set_array(G_t/1e3)
            self.tomo_img.set_clim(vmin,vmax)
            anim_obj = (self.tomo_img,)
 
 
        
        if self.showTe:
            #relative shift in the mode phase between SXR and ECE
            Phi_ece = self.Te2Dmap.Phi0
            self.Te2Dmap.f0 = self.F1
            
            dPhi = self.n*(Phi_ece-self.Phi0)/self.m

            self.Te2Dmap.shift_phi = self.shift_phi +dPhi+self.dPhi/self.n
            anim_obj += self.Te2Dmap.update(update_cax=False,update_mag=False,animate=animate,
                                ax = self.ax, filled_contours=False)
            
        if animate:
            return anim_obj
        
        self.fig.canvas.draw_idle()
    
 
 
    def MouseWheelInteraction(self,event):
        if event.inaxes is self.ax:
            if self.keyCtrl:
                self.dPhi += event.step*2*np.pi/360 
                print('ECE is shifted by %d deg'%(np.rad2deg(self.dPhi)%360))

            else:
                steps = 360 if self.keyCtrl else 36
                self.shift_phi +=  event.step*2*np.pi/steps
             
            self.shift_phase(self.shift_phi)
            self.update( update_mag=False)

        
    def onKeyPress(self,event):
        step = 36
        if 'control' == event.key:
            self.keyCtrl=True
            step = 360
        if 'shift' == event.key:
            self.keyShift=True

        sign = None
        if 'left' == event.key:
            sign = -1

        elif 'right' == event.key:
            sign = 1

        if sign is not None:
            self.shift_phi+= sign*2*np.pi/step
            self.shift_phase(self.shift_phi)
            self.update(update_mag=False)
            

    def onKeyRelease(self,event):
        if 'control' == event.key:
            self.keyCtrl=False
        if 'shift' == event.key:
            self.keyShift=False
  
    def __del__(self, event=None):
        try:
            for h in self.hsolvers:
                h.qin.put(None) #send kill signal
                h.terminate()
                h.join(1e-2)  #join process
                h.qin.close()
                h.qout.close()
                del h  #release memory 
            #print('solvers are closed')
            self.hsolvers = []
            self.initialized = False
            if plt.fignum_exists('Retrofit'):
                plt.close(plt.figure('Retrofit'))
            
            import gc
            gc.collect()
        except Exception as e:
            print('Error __del__', e)
            

