#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys, os, random, time, argparse, logging
import traceback
import matplotlib  
from copy import deepcopy

from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
import configparser
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT
matplotlib.rcParams['backend'] = 'Qt5Agg' 

fmt = logging.Formatter('%(asctime)s | %(name)s | %(levelname)s: %(message)s', '%H:%M:%S')
hnd = logging.StreamHandler()
hnd.setFormatter(fmt)
logger = logging.getLogger('pyspecview')
logger.addHandler(hnd)
logger.setLevel(logging.DEBUG)
#logger.setLevel(logging.INFO)


 
# http://matplotlib.org/users/dflt_style_changes.html
params = {'legend.fontsize': 'large', 
            'axes.labelsize': 'large', 
            'axes.titlesize': 'large', 
            'xtick.labelsize' :'medium', 
            'ytick.labelsize': 'medium', 
            'font.size':12, 
            'mathtext.fontset': 'cm', 
            'mathtext.rm': 'serif', 
            'grid.color': 'k', 
            'grid.linestyle': ':', 
            'grid.linewidth': 0.5, 
            'lines.linewidth'   : 1.0, 
            'lines.dashed_pattern' : (6, 6), 
            'lines.dashdot_pattern' : (3, 5, 1, 5), 
            'lines.dotted_pattern' : (1, 3), 
            'lines.scale_dashes': False, 
            'errorbar.capsize':3, 
            'mathtext.fontset': 'cm', 
            'mathtext.rm' : 'serif' }
matplotlib.rcParams.update(params)

from matplotlib.figure import Figure
from matplotlib.widgets import Slider, RadioButtons, RectangleSelector
from matplotlib.ticker import NullFormatter, ScalarFormatter, MaxNLocator, AutoMinorLocator, FixedLocator
from matplotlib.widgets import MultiCursor
import matplotlib.pylab as plt

import numpy as np
from scipy.stats.mstats import mquantiles
from scipy.interpolate import interp1d
from scipy.signal import butter, filtfilt, hilbert, detrend
        
try:
    from scipy.signal.signaltools import _next_regular as next_fast_len
except:
    try:
        from scipy.fftpack.helper import next_fast_len
    except:
        from scipy.fftpack import next_fast_len 

from stft import stft
from sfft import sfft
from sstft import sstft

random_seed = 10
font_size = 10

from IPython import embed

 

 
#prevents killing the GUI whenexception is raised
def excepthook(exc_type, exc_value, exc_tb):
    # Note: sys.__excepthook__(...) would not work here.
    # We need to use print_exception(...):
    print('Something went wrong. Sent error message odstrcilt@fusion.gat.com')
    traceback.print_exception(exc_type, exc_value, exc_tb)

 

 
sys.excepthook = excepthook



class colorize:
    def __init__( self, hue, invert=True):

        # hue is a list of hue values for which the colors will be evaluated
        from colorsys import hsv_to_rgb
        self.colors = np.array([hsv_to_rgb(h, 1, 1) for h in hue])
        self.invert = invert

        
    def __call__(self, index, r):
        
        value = r/np.amax(r)
        
        #get fully saturated colors
        color = self.colors[index, :]
        
        #reduce the lightness 
        color *= np.asarray(value)[..., None]
        
        #invert colors 
        if self.invert:   color = 1-color

        #remove a rounding errors
        color[color<0] = 0 
        color[color*256>255] = 255./256

        return color


def min_fine(x, y):
 
    #find extrem with a subpixel precision
    i = np.argmin(y)
    if i==0 or i == len(y)-1: return y[i], x[i]    #not working at the edge! 
    A = x[i-1:i+2]**2, x[i-1:i+2], np.ones(3)
    a, b, c = np.linalg.solve(np.array(A).T, y[i-1:i+2])
    
    ymin =  c-b**2/(4*a)
    xmin = -b/(2*a) 
    return  ymin, xmin


class MyFormatter(ScalarFormatter):   # improved format of axis
    def __call__(self, x, pos=None):
        self.set_scientific(True)
        self.set_useOffset(True)
        self.set_powerlimits((-3, 3))
        return ScalarFormatter.__call__(self, x, pos)

class NavigationToolbar(NavigationToolbar2QT):
    # only display the buttons we need
    toolitems = [t for t in NavigationToolbar2QT.toolitems if
                 t[0] in ('Home', 'Pan', 'Zoom')]


class SpectraViewer(object):
    
    #class for the spectrum and data plot handling 
    
    def __init__(self, sgamma, stau, method, message_out, show_raw=True, show_colorbar=False, 
                 allow_selector=True, fig_name=None, fig=None, phase_analysis=False, cmap= 'gnuplot2'):
        """
        
        
        """
        matplotlib.rcParams['xtick.direction'] = 'out'
        matplotlib.rcParams['ytick.direction'] = 'out'
        matplotlib.rcParams['xtick.major.size'] = 4
        matplotlib.rcParams['xtick.minor.size'] = 2
        matplotlib.rcParams['ytick.major.size'] = 4
        matplotlib.rcParams['ytick.minor.size'] = 2

        self.initialized=False
        self.mode_traced = False

        axcolor = 'lightgoldenrodyellow'
        self.show_raw = show_raw
        self.show_colorbar = show_colorbar

        self.phase_analysis = phase_analysis
        
        if fig is None:
            self.fig, self.ax = subplots(num=fig_name)
        else:
            self.fig, self.ax = fig, fig.gca()
        self.cax = None
        self.uax = None
        self.message_out = message_out
        self.cmap = cmap

        self.fig.subplots_adjust(left=0.11, top=0.95, right=0.98, bottom=0.1)
        if show_raw:
            self.fig.subplots_adjust(top=0.85)
            self.uax = self.fig.add_axes([0.11, 0.85, 0.87, 0.13])
            self.uax.xaxis.set_major_formatter(NullFormatter())
            self.uax.yaxis.set_major_formatter(NullFormatter())
        
        elif show_colorbar:
            self.fig.subplots_adjust(right=0.9)
            self.cax = self.fig.add_axes([0.92, 0.1, 0.02, 0.85])
            self.cax.xaxis.set_major_formatter(NullFormatter())
            self.cax.yaxis.set_major_formatter(NullFormatter())
            self.cax.set_ylim(0, 1)

        self.sdpi = 1000.

        self.ax.set_xlabel('Time [s]', fontsize=font_size)
        self.ax.set_ylabel('Frequency [Hz]', fontsize=font_size)
        
 
        for label in (self.ax.get_xticklabels() + self.ax.get_yticklabels()):
            label.set_fontsize(font_size) # Size here overrides font_prop

        self.methods = 'time', 'freq.', 'sparse'

        self.sgamma = sgamma
        self.stau   = stau
        self.method = method
        
        #set parameters of the gamma slider
        sgamma.setRange(0.001*self.sdpi, 1*self.sdpi)
        sgamma.setValue(.5*self.sdpi)
        sgamma.setTracking(True)
        sgamma.setTickPosition(QSlider.NoTicks)
        sgamma.setSingleStep(.1)

        #set parameters of the NFFT slider
        self.tau0 = .6 #dimensionaless parameter for time resolution between 0 and 1
        stau.setRange(.3*self.sdpi, .8*self.sdpi)
        stau.setValue(self.tau0*self.sdpi)
        stau.setTracking(True)
        stau.setTickPosition(QSlider.NoTicks)
        stau.setSingleStep(1)

        #set parameters of the method menu
        
        self.method.addItem("STFT")
        self.method.addItem("SFFT")
        self.method.addItem("Sparse")
        self.DFT_backend = 'time'
        ind_method = [i for i, n in enumerate(self.methods) if n == self.DFT_backend][0]
        self.method.setCurrentIndex( ind_method)
        
        #default time and frequency range
        self.t_range = np.nan, np.nan
        self.f_range = 0, 2e5

        matplotlib.rcParams['xtick.direction'] = 'in'
        matplotlib.rcParams['ytick.direction'] = 'in'
        
        #time trace of the mhd mode 
        self.plt_trace, = self.ax.plot([], [], 'gray', zorder=99)

        c  = 1-np.array(plt.cm.get_cmap( cmap)(0))[:3] #inverse of the lowest color in the colormap
        self.plt_plasma_freq_n1, = self.ax.plot([], [], c=c, zorder=99, lw=.5)
        #self.plt_plasma_freq_n2, = self.ax.plot([], [], c=c, zorder=99, lw=.5,ls='--')

        if allow_selector:
            rectprops = dict(facecolor='gray', edgecolor='black', alpha=0.5, fill=True, zorder=99)
            self.RS1 = RectangleSelector(self.ax, self.line_select_callback, 
                                        drawtype='box', useblit=True, 
                                        button=[1, ], # don't use middle button
                                       minspanx=5, minspany=5, rectprops=rectprops, 
                                       spancoords='pixels')
        
    def __del__(self):
       pass

    def reset(self):
        #clear axis, remove images from the previous discharge
        if not self.initialized:
           return 

        self.initialized=False
    
        self.dR_corr = 0
        self.dZ_corr = 0  


        for artist in self.ax.get_images()+self.ax.lines+self.ax.artists:
            artist.remove()
            del artist
        self.ax.figure.canvas.draw_idle()
  
        if self.show_raw:
            for artist in self.uax.lines+self.uax.artists:
                artist.remove()
            self.uax.figure.canvas.draw_idle()
            
        self.plot_description.remove()
       
    def init_plot(self, data, window=None, tmin=None, tmax=None, 
                  fmin0=-np.infty, fmax0=np.infty, description='', mode_range=None, mode_num_lbl='' ):
        
        #object of the spectrogram
        
        if window is None and not hasattr(self, 'window'):
            window = 'gauss'
        elif window is None:
            window = self.window
            
        self.window = window
        self.description = description
        self.nt = data['tvec'].size
        self.dt = (data['tvec'][-1]-data['tvec'][0])/(self.nt-1)  #BUG not valid for not equally spaced time vectors!

        self.nfft0 = self.get_nfft(self.tau0)
        gamma = self.sgamma.value()/self.sdpi
        self.stft_img =  STFTImage(self.ax, gamma=gamma, cmap=self.cmap, 
                                   phase_analysis=self.phase_analysis, 
                                   colorbar_ax=self.cax, mrange=mode_range, 
                                   mode_num_lbl=mode_num_lbl)
        
        if tmin is None: tmin =  data['tvec'][0]
        if tmax is None: tmax =  data['tvec'][-1]

        #clear the image
        for im in self.ax.get_images():
            im.remove()
            del im

        self.plt_trace.set_data([], [])

        self.plt_plasma_freq_n1.set_visible(False)
        #self.plt_plasma_freq_n2.set_visible(False)
     
        if 'freq' in data:
            self.plt_plasma_freq_n1.set_data(data['freq_tvec'], data['freq'])
            #self.plt_plasma_freq_n2.set_data(data['freq_tvec'], data['freq']*2)
            self.plt_plasma_freq_n1.set_visible(True)
            #self.plt_plasma_freq_n2.set_visible(True)

   
        if self.show_raw:

            #clear it before use
            for l in self.uax.get_lines():
                l.remove()
                del l
            
            #object of the plot with raw data
            self.data_plot = DataPlot(self.uax, lw=.5, c='k')
            self.data_plot.prepare(data, tmin, tmax)

            #object of the both + STFT algorithms 
            self.stft_disp = STFTDisplay(data, self.stft_img, self.ax, self.data_plot, 
                        nfft=self.nfft0, win=window, phase_analysis=self.phase_analysis, 
                        colorbar_ax=self.cax, method=self.DFT_backend)
            
            self.cid3 = self.uax.callbacks.connect('xlim_changed', 
                                            self.stft_disp.plot_ax_update)

        else:
            self.stft_disp = STFTDisplay(data, self.stft_img, self.ax, nfft=self.nfft0, 
                            win=window, phase_analysis=self.phase_analysis, 
                            colorbar_ax=self.cax, method=self.DFT_backend)
    
        try:
            x, y, Z = self.stft_disp.__call__(tmin, tmax, fmin0, fmax0)
            if np.size(y):
            	fmax0 = min(fmax0, y[-1])
            	fmin0 = max(fmin0, y[0])
                
        except Exception as e:
            logger.error('Error: ', exc_info=True)
            logger.error( traceback.format_exc())
            return 

        if x.size < 2 or y.size < 2:
            self.ax.set_ylim( fmin0, fmax0)
            self.ax.set_xlim(tmin, tmax)
            logger.error('time range or frequency range is too small')
            raise

        self.ax.set_ylim(y[0], y[-1])
        self.ax.set_xlim(tmin, tmax)
        self.stft_img.prepare(x[0], x[-1], y[0], y[-1], Z)
        
        if self.initialized:
            self.plot_description.set_text(description)
        else:
            self.plot_description = self.ax.text(1.006, .05, description, rotation='vertical', 
                        transform=self.ax.transAxes, verticalalignment='bottom', 
                        size='xx-small', backgroundcolor='none', zorder=100)

        self.cid1 = self.ax.callbacks.connect('xlim_changed', self.stft_disp.ax_update)
        self.cid2 = self.ax.callbacks.connect('ylim_changed', self.stft_disp.ax_update)
        self.cid4 = self.fig.canvas.mpl_connect('resize_event', self.stft_disp.ax_update)

        self.ax.xaxis.set_minor_locator(AutoMinorLocator())
        self.ax.yaxis.set_minor_locator(AutoMinorLocator())

        self.cid5 = self.fig.canvas.mpl_connect('scroll_event', self.MouseWheelInteraction)
        self.cid6 = self.fig.canvas.mpl_connect('button_press_event', self.MouseClickInteraction)
        self.cid7 = self.fig.canvas.mpl_connect('button_release_event', self.MouseClickInteraction) #BUG if the mouse released out of the axis/figure, it will not return to the original state
        self.MWI_id = 0
        self.MCI_id = 0

        self.stft_disp.set_yticks()
        self.initialized=True
    
 
    def line_select_callback(self, eclick, erelease):
        'eclick and erelease are the press and release events'

        x1, y1 = eclick.xdata, eclick.ydata
        x2, y2 = erelease.xdata, erelease.ydata

        self.t_range = np.sort((x1, x2))
        self.f_range = np.sort((y1, y2))
        
        self.RS1.to_draw.set_visible(True)
        self.RS1.canvas.draw()

    def select_DFT_backend(self, label):
        self.DFT_backend = self.methods[label]
        if not self.initialized:  return 
    
        self.stft_disp.method = self.methods[label]
        
        xstart, xend, ystart, yend = self.stft_img.ax.axis()
        
        self.stft_disp.ax_update(self.stft_img.ax, force_update=True)
        
        self.apply_slider(0)
        
    def set_mode_range(self, mode_min, mode_max):
        self.stft_img.mode_range = mode_min, mode_max
        self.stft_img.init_mode_shuffle()
        self.stft_img.prepare_colorbar()
        self.stft_disp.ax_update(force_update=True)  #recalculate
        
    def updateGammaSlider(self, val):
        self.sgamma.setValue(val*self.sdpi)

    def updateTauSlider(self, val):
        self.stau.setValue(val*self.sdpi)
 
    def get_nfft(self, x):
        if self.DFT_backend in ('freq.', 'sparse'):  x = 1.2-x #Guess to hava roughly similar size of window
        
        tau = self.nt**x*self.dt
        nfft = tau/self.dt if self.window == 'gauss' else 2**ceil(np.log2(tau/self.dt))
        return nfft
    
    def get_tau0(self):
        return self.stau.value()/self.sdpi
    
    def apply_slider(self, val):
        if not self.initialized:  return 

        gamma = self.sgamma.value()/self.sdpi
 
        self.tau0 =  self.get_tau0()
        nfft = self.get_nfft( self.tau0)
        
        xstart, xend, ystart, yend = self.stft_img.im.get_extent()
        if nfft!= self.stft_disp.nfft:
            self.stft_disp.nfft = nfft
            #NOTE tau is not exatly equal to the width (sigma) of the gaussin window
            self.message_out('nfft: %d  window: %.2gms'%(nfft, nfft*self.dt*1000), 1000)
            x, y, Z = self.stft_disp.__call__(xstart, xend, ystart, yend)
            self.stft_img.z = Z
        if gamma!= self.stft_img.gamma:
            self.stft_img.gamma = gamma
                
        self.stft_img.update_image(xstart, xend, ystart, yend, self.stft_img.z)     

        self.stft_disp.ax_update(self.stft_img.ax)


    def calc_sig_power(self, tvec_f, fvec, df):
        #cauclated a power of signal at given frequency fvec, evolving in time 
        #BUG make it better?
        T = time.time()
                
        x = self.stft_disp.tvec
        y = self.stft_disp.signal
        ind = slice(*x.searchsorted(tvec_f[[0, -1]]))
        x = x[ind]
        if y.ndim > 1: 
            y = y[:, 0]
            
        y = np.single(np.copy(y[ind]))

        dt = (x[-1]-x[0])/x.size
        f = np.interp(x, tvec_f, fvec)
        corr = np.exp((dt*1j*2*np.pi)*np.cumsum(np.single(f)))

        #detrend 
        y -= (y[-1]-y[0])/(x[-1]-x[0])*(x-x[0])+y[0]

        f_nq = .5/dt 
        b, a = butter(5, df/f_nq)

        ylow = filtfilt(b, a, y/corr)
        
        return x, np.abs(ylow)


    def trace_spectrogram(self, f0, t0, delta_f ):
        
        tvec = self.stft_disp.x
        fvec = self.stft_disp.y 
        if self.phase_analysis:
            spect = np.copy(self.stft_img.z[0])
            it0 = tvec.searchsorted(t0)
            if0 = fvec.searchsorted(f0)
            spect[self.stft_img.z[1] != self.stft_img.z[1][if0, it0]] = 0 #select only one mode number
        else:
            spect = self.stft_img.z

        F = [f0]
        it0 = tvec.searchsorted(t0)

        noise_level = np.median(spect)*2
        #forward pass
        for k in range(it0, tvec.size):
            try:
                w = np.exp(-(fvec-F[-1])**2/(2*delta_f**2))    
            except:
                print('delta_f, f0, t0, F' , delta_f, f0, t0, F)
                raise

            ifmax = np.argmax(spect[:, k]*w)
            f = fvec[ifmax]
            A = -(spect[:, k]*w)[ifmax]
            
            #fine tuning 
            w = np.exp(-(fvec-f)**2/(2*delta_f**2))
            A, f = min_fine(fvec, -spect[:, k]*w)

            #if it will lost the mode
            if -A < noise_level:  f = F[-1]

            F.append(f)

        F = F[1:]
        F.reverse()

        #backward pass
        for k in range(it0-1, -1, -1):
            w = np.exp(-(fvec-F[-1])**2/(2*delta_f**2))    
            ifmax = np.argmax(spect[:, k]*w)
            f = fvec[ifmax]
            A = -(spect[:, k]*w)[ifmax]
            
            #fine tuning 
            w = np.exp(-(fvec-f)**2/(2*delta_f**2))    
            A, f = min_fine(fvec, -spect[:, k]*w)
                    
            #if it will lost the mode
            if -A < noise_level:  f = F[-1]

            F.append(f)
            
        F.reverse()
        
        tvec = np.asarray(tvec)
        self.plt_trace.set_data(tvec, F )
        self.fig.canvas.draw_idle()
        
        return tvec, np.array(F)


    def MouseClickInteraction(self, event):

        if self.MCI_id == id(event):
            return #sometimes the same event arrives several times :(
        self.MCI_id = id(event)
        
        QApplication.setOverrideCursor(Qt.WaitCursor)
        
        # trace the mode
        if  event.button == 1 and event.dblclick: #left button, double click
            
            self.mode_traced = True

            t0, f0 = event.xdata, event.ydata

            tau = self.get_nfft(self.tau0)*self.dt
            
            delta_f =  2*np.pi/tau
            
            if not f0 is None:
                T, F = self.trace_spectrogram(f0, t0, delta_f )
            
                tvec_power, Power = self.calc_sig_power( T, F , delta_f)
                self.mode_trace = tvec_power, Power
                xmin, xmax = self.stft_disp.im_ax.get_xlim()
                if self.show_raw:
                    self.data_plot.update_plot(xmin, xmax, tvec_power, Power)
                    
                self.fig.canvas.draw()   
            
        # remove the trace of the the mode
        elif  event.button == 3 and event.dblclick: #   right button, double click

            self.mode_traced = False
            
            self.plt_trace.set_data([], [] )
            xmin, xmax = self.stft_disp.im_ax.get_xlim()
            if self.show_raw:
                self.data_plot.update_plot(xmin, xmax, self.stft_disp.tvec, self.stft_disp.signal )
            self.fig.canvas.draw_idle()
               
        elif  event.button == 2: #middle botton - panning 
            
            if event.name  == 'button_press_event':
                self.panning_start_pos = event.xdata, event.ydata
 
            if event.name  == 'button_release_event' and hasattr(self, 'panning_start_pos'):

                curr_xlim = self.ax.get_xlim()
                curr_ylim = self.ax.get_ylim()

                #mouse was released out of the window
                if event.xdata is None:
                    event.xdata = float(event.x)/self.stft_disp.width *np.diff(curr_xlim)+curr_xlim[0]
                    event.ydata = float(event.y)/self.stft_disp.height*np.diff(curr_ylim)+curr_ylim[0]

                panning_shift_x = event.xdata-self.panning_start_pos[0]
                panning_shift_y = event.ydata-self.panning_start_pos[1]

                self.stft_disp.prevent_update = True  #trick to save one update :( 
                self.ax.set_xlim(curr_xlim-panning_shift_x )
                self.stft_disp.prevent_update = False #trick to save one update :( 
                self.ax.set_ylim(curr_ylim-panning_shift_y )
        
        elif  event.button == 1 and hasattr(self.stft_img, 'im'): #left button. write frequency and time at this point 
            xstart, xend, ystart, yend = self.stft_img.im.get_extent()

            if self.ax == event.inaxes:
                x, y = event.xdata, event.ydata

                if self.phase_analysis:
                    ix = int(round((x-xstart)/(xend-xstart)*(self.stft_img.z[1].shape[1]-1)))
                    iy = int(round((y-ystart)/(yend-ystart)*(self.stft_img.z[1].shape[0]-1))) 
                    try:
                        f = self.stft_img.z[1].T[ix, iy]
                        n = self.stft_img.N_mode[0][self.stft_img.N_mode[1]==f]
                        self.message_out('t: %.5fs f: %.3gkHz, mode num.:%d'%(x, y/1e3, n), 1000)
                    except:
                        pass
                else:
                    ix = min(max(0, (x-xstart)/(xend-xstart)), 1)
                    ix = int(round(ix*(self.stft_img.z.shape[1]-1)))
                    
                    iy = min(max(0, (y-ystart)/(yend-ystart)), 1)
                    iy = int(round(iy*(self.stft_img.z.shape[0]-1)))
                    try:
                        f = self.stft_img.z[iy, ix]
                        self.message_out('t: %.5fs f: %.3gkHz val:%.2e'%(x, y/1e3, f), 1000)
                    except:
                        pass
                    
            #select only a single mode number
            if self.phase_analysis:
                
                if hasattr(self, 'orig_z') and not self.orig_z is None: #draw original image
                    self.stft_img.update_image(xstart, xend, ystart, yend, self.orig_z)   
                    self.orig_z = None

                if self.cax == event.inaxes and event.name  == 'button_press_event':
                    #for older matplotlib
                    if all(np.array(self.cax.get_ylim()) == np.array((0.,1.))):
                        n = self.stft_img.N_mode[1][int(event.ydata*len(self.stft_img.N_mode[0]))]
                    else: #new matploltib
                        n_ = int(round(event.ydata))
                        n = self.stft_img.N_mode[1][self.stft_img.N_mode[0]==n_]
               
                    self.orig_z = np.copy(self.stft_img.z[0]), np.copy(self.stft_img.z[1])
                    self.stft_img.z[0][self.stft_img.z[1] != n] = 0
                    self.stft_img.update_image(xstart, xend, ystart, yend, self.stft_img.z)  
        
        QApplication.restoreOverrideCursor()            


    def MouseWheelInteraction(self, event):
        
        factor = 0.8  #zooming/unzooming factor 

        if self.MWI_id == id(event):
            return #sometimes the same event arrives several times :(
        
        self.MWI_id = id(event)

        QApplication.setOverrideCursor(Qt.WaitCursor)
    
        if event.inaxes == self.ax:
            #  zoom by mouse wheel
            curr_xlim = self.ax.get_xlim()
            curr_ylim = self.ax.get_ylim()
            new_width = (curr_xlim[1]-curr_xlim[0])*factor**event.step
            new_height= (curr_ylim[1]-curr_ylim[0])*factor**event.step

            relx = (curr_xlim[1]-event.xdata)/(curr_xlim[1]-curr_xlim[0])
            rely = (curr_ylim[1]-event.ydata)/(curr_ylim[1]-curr_ylim[0])

            self.stft_disp.prevent_update = True #trick to save one update :( 
            self.ax.set_ylim([event.ydata-new_height*(1-rely), event.ydata+new_height*rely])
            self.stft_disp.prevent_update = False
            self.ax.set_xlim([event.xdata- new_width*(1-relx), event.xdata+ new_width*relx])

        elif self.show_raw and event.inaxes == self.uax:
            #  zoom by mouse wheel in upper plot 
            curr_xlim = self.uax.get_xlim()
            new_width = (curr_xlim[1]-curr_xlim[0])*factor**event.step

            relx = (curr_xlim[1]-event.xdata)/(curr_xlim[1]-curr_xlim[0])
            xstart = event.xdata- new_width*(1-relx)
            xend = event.xdata+ new_width*(relx)
            self.ax.set_xlim([xstart, xend])

        self.apply_slider(0)
        
        QApplication.restoreOverrideCursor()


class STFTImage():

    def __init__(self, ax, gamma=.5, cmap='gnuplot2', phase_analysis=False, 
                 colorbar_ax=None, mrange=None, mode_num_lbl=''):

        self.ax = ax
        self.colorbar_ax = colorbar_ax
        self.gamma = gamma
        self.cmap = cmap
        self.phase_analysis = phase_analysis
        self.mode_num_lbl = mode_num_lbl

        self.mode_range = mrange
        
        if self.phase_analysis: 
            self.init_mode_shuffle()

    def init_mode_shuffle(self):
        N_min, N_max = self.mode_range
    
        N = np.arange(N_min, N_max+1, dtype='int8')
        np.random.seed(random_seed)
        np.random.shuffle(N)
        self.N_mode = np.arange(N_min, N_max+1, dtype='int8'), N[::-1]

    def prepare(self, xstart, xend, ystart, yend, z):
        #show spectrogram
  
        if self.phase_analysis:  #phase core
            self.prepare_colorbar()

        self.im = self.ax.imshow(np.zeros((2, 2)), origin='lower', extent=(xstart, xend, ystart, 
                   yend), aspect='auto', cmap=self.cmap, interpolation='gaussian', rasterized=True)#
        
        self.update_image(xstart, xend, ystart, yend, z)

    
    def prepare_colorbar(self):
    
        N_min = self.N_mode[0].min()
        N_max = self.N_mode[0].max()
        
        self.colors = colorize(np.linspace(0, 1, N_max-N_min+1, endpoint=False) )
        self.cmap = matplotlib.colors.ListedColormap(self.colors(self.N_mode[1], 1))

        bounds = np.linspace(N_min-.5, N_max+.5, N_max-N_min+2  )
        if not self.colorbar_ax is None:
            self.colorbar_ax.cla()
            self.cb = matplotlib.colorbar.ColorbarBase(self.colorbar_ax, 
                        cmap=self.cmap , boundaries=bounds, extendfrac='auto', 
                        ticks=np.arange(N_min, N_max+1), spacing='uniform')
            self.colorbar_ax.set_xlabel(self.mode_num_lbl, fontsize='small')
            self.colorbar_ax.xaxis.set_label_position('top') 
            self.colorbar_ax.tick_params(labelsize=8) 

       
    def GammaTransform(self, z, typ='log'):
        if type(z) is tuple:
            return self.GammaTransform(z[0], typ=typ), z[1] 
        if np.iscomplexobj(z):
            z = np.abs(z)
   
        if typ == 'log' and z.size > 1:
            #transformation by log
            ind = np.random.randint(0, z.size-1, size=10000)
            vmin = mquantiles(z.flat[ind], .5)[0]
            z_norm = np.tan(np.pi/2.000001*self.gamma)*np.nanmean(z)
            out = np.log1p(np.maximum(0, z-vmin*np.sqrt(self.gamma*2))/z_norm)

            return out
      
        if typ == 'gamma':#standard definition of gamma
            return z**self.gamma

        
    def update_image(self, xstart, xend, ystart, yend, z):
        #update spectrogram image
        self.z = z
        z = self.GammaTransform(self.z)
        t2 = time.time()

        if self.phase_analysis:  #phase core
            N_min, N_max = -6, 5
            if not self.mode_range is None:
                N_min, N_max = self.mode_range

            z = self.colors( z[1], z[0])

        self.im.axes = self.ax
        try:
            self.im.set_data(z)
        except:
            logger.error(  'error self.im.set_data(z)')
  
        self.im.set_extent((xstart, xend, ystart, yend))
 
        t3 = time.time()

        self.update_image_clim(xstart, xend, ystart, yend)

   
    def update_image_clim(self, xstart, xend, ystart, yend):
        E = self.im.get_extent()

        self.ax.autoscale_view(tight=True)
        
        if self.phase_analysis:
            self.im.autoscale()
        else:
            z = self.im.get_array()

            iy1, iy2 = np.int_((np.array((xstart, xend))-E[0])/(E[1]-E[0])*z.shape[1])
            ix1, ix2 = np.int_((np.array((ystart, yend))-E[2])/(E[3]-E[2])*z.shape[0])
 
            #choose optimal value to use the whole colorbar range efficiently 
            vmax1 = mquantiles(z[ix1:ix2, iy1:iy2].max(0), 0.95)[0]
            vmax2 = mquantiles(z[ix1:ix2, iy1:iy2].max(1), 0.95)[0]
            vmax = min(vmax1, vmax2)
   
            self.im.set_clim(0, vmax)
        self.ax.figure.canvas.draw_idle()
        
        
class DataPlot():

    def __init__(self, ax=None, **kwarg):
        self.kwarg = kwarg
        self.ax_raw = ax
        self.tmin = -np.infty
        self.tmax = np.infty

    def prepare(self, data, tmin, tmax):
        self.x = data['tvec']

        if data['signal'].ndim == 1:
            self.y = data['signal']
        elif data['signal'].ndim == 2:
            self.y = data['signal'][0]

            
        self.data_timetrace, = self.ax_raw.plot([], [], antialiased=True, **self.kwarg)
        self.update_plot(tmin, tmax)
        
    def update_plot(self, xstart, xend, xnew=None, ynew=None):
        #dims = self.ax.axesPatch.get_window_extent().bounds
         
        xlim =  self.ax_raw.get_xlim()

        if xstart == self.tmin and xend == self.tmax and xnew is None and ynew is None:
            return

        self.tmin = xstart
        self.tmax = xend
        
        dims = self.ax_raw.get_window_extent().bounds

        width = int(dims[2] + 0.5)
        
        if xnew is not None:
            self.x = xnew
        
        if ynew is not None:
            self.y = ynew
            
        istart, iend = self.x.searchsorted((xstart, xend))
        
        #random sampling
        random.seed(0)
        ind = np.linspace(istart+1, iend-2, width*20)+np.random.rand(width*20)*2-1
        ind = np.int_(np.unique(np.round(ind)))

        try:
            x = self.x[ind]
            y = self.y[ind].real
        except Exception as e:
            logger.error('update_plot error %s %d %d %d', str(e), len(self.x), len(ind), len(self.y))
            x , y = [0, ], [0, ]
            
        self.data_timetrace.set_data(x, y)
        self.ax_raw.set_xlim(xstart, xend)
        self.ax_raw.set_ylim(np.amin(y), np.amax(y))


class STFTDisplay():


    def __init__(self, data, image, im_ax, data_plot=None, method='time', 
                 nfft=2**13, win='gauss', phase_analysis=False, colorbar_ax=None):
 
        dims = im_ax.get_window_extent().bounds

        self.width  = int(dims[2] + 0.5) 
        self.height = int(dims[3] + 0.5)
        
        self.nfft = nfft
        self.phase_analysis = phase_analysis
        self.signal = data['signal']
        self.tvec = data['tvec']
        if 'phase_correction' in data and not data['phase_correction'] is None:
            self.phase_correction =  data['phase_correction']
        if 'Phi' in data:  #toroidal or polodial angle 
            self.Phi =  data['Phi']
        self.x = None
        self.y = None
        self.image = image
        self.data_plot = data_plot
        self.im_ax = im_ax
        self.colorbar_ax = colorbar_ax
        self.method = method
        self.dt = (self.tvec[-1]-self.tvec[0])/(self.tvec.size-1)

        #fast guess of the lowest dt (tvec can have variable sampling frequency!), 
        #but 5% error is there due to a single precision!!!

        self.window = win
        self.output = {}
        self.prevent_update = False
     

    def __call__(self, xstart, xend, ystart, yend):
        #compute the signal transformation
        QApplication.setOverrideCursor(Qt.WaitCursor)

        t = time.time()
        
        if xend   < self.tvec[0]:  xend   = self.tvec[0]+0.1
        if xstart > self.tvec[-1]: xstart = self.tvec[-1]-0.1
        if xstart > xend: xstart, xend = xend, xstart

        if self.method == 'time':
            A, fv, tv = stft(self.tvec, self.signal, self.nfft, resolution=self.width, 
                    window=self.window, tmin=xstart, tmax=xend, fmin=ystart, 
                    fmax=yend, pass_DC=False, complex_spectrum=self.phase_analysis)
            if fv.size > 10*self.height: #downsample is if it too large 
                A = A[:, ::fv.size//(10*self.height)]
                fv = fv[::fv.size//(10*self.height)]

        elif self.method == 'freq.':
            A, fv, tv = sfft(self.tvec, self.signal, self.nfft, resolution=self.height, 
                        window=self.window, fmin=ystart, fmax=yend, tmin=xstart, 
                        tmax=xend, pass_DC=False, complex_spectrum=self.phase_analysis)
            
        elif self.method == 'sparse':
            if self.phase_analysis:
                raise Exception('Not yet implemented')
            A, fv, tv = sstft(self.tvec, self.signal, self.nfft, xstart, xend, ystart, 
                           yend, self.width, self.height, zoom=4)

        #correction of the phase for the balooning coils
        if hasattr(self, 'phase_correction') and self.phase_analysis:
 
            correction = [np.interp(fv, self.phase_correction[:, i*2], self.phase_correction[:, i*2+1])\
                        for i in range(self.phase_correction.shape[1]//2)]
            correction = np.vstack(correction).T
            A*= np.exp(1j*correction)
       
        logger.info( 'FFT done in:  %5.2fs', (time.time() - t))
         
        #mode number estimation
        if self.phase_analysis:
            absA = np.abs(A)
            A/= absA; angA = A

            Phi = np.exp(-1j*self.Phi)
            n = self.image.N_mode[0]

            cmplxPhi = np.complex64(Phi[None, :]**n[:, None])
            
            ##go throught all mode numbers and find the one which gives the closest values of angles in all coils
            ang_norms = np.tensordot(angA, np.conj(cmplxPhi)/len(Phi), [2, 1])
            ang_norms = np.abs(ang_norms) #best value of 'n' will give ang_norms closest to circle abs(x) = 1
  
            N = self.image.N_mode[1]
            n_mode = N[np.argmax(ang_norms, -1)]
            absA = np.mean(absA, -1)
            A = absA.T, n_mode.T
        else:
            A = A.T
     
        self.x = tv
        self.y = fv

        self.output['fvec'] = fv
        self.output['tvec'] = tv
        self.output['spect'] = A

        QApplication.restoreOverrideCursor()

        return tv, fv, A
    
            
    def set_yticks(self):
        #plot yticks in Hz, kHz, MHz
             
        self.im_ax.yaxis.set_major_locator(MaxNLocator(6))
        yticks = self.im_ax.get_yticks()
        
        if yticks.max() < 1e3:
            unit, fact =  'Hz', 1
        elif yticks.max()  < 1e6:
            unit, fact = 'kHz', 1e3
        else:
            unit, fact = 'MHz', 1e6
            
        ndig = max(0, int(np.ceil(-np.log10(np.mean(np.diff(yticks)))) + np.log10(fact)))
        
        # fixing yticks with matplotlib.ticker "FixedLocator"
        self.im_ax.yaxis.set_major_locator(FixedLocator(yticks))

        self.im_ax.set_yticklabels([('%.'+str(ndig)+'f')%f for f in yticks/fact]) 
        self.im_ax.set_ylabel('Frequency [%s]'%unit, fontsize=font_size)


    def ax_update(self, ax=None, fine_resolution=False, force_update=False):
        #called when xlim or ylim was changed or window was resized 
        if self.prevent_update: return 
        ax = self.im_ax

        ax.set_autoscale_on(False) # Otherwise, infinite loop

        #Get the number of points from the number of pixels in the window
        dims = ax.get_window_extent().bounds

        width  = int(dims[2] + 0.5) 
        height = int(dims[3] + 0.5)
        
        if width!= self.width or height!= self.height:
            force_update = True
        
        self.width  = width
        self.height = height

        #increase resolution
        if fine_resolution and self.method != 'sparse':
            self.width  *= 2
            self.height *= 2

        #Get the range for the new area
        xstart, ystart, xdelta, ydelta = ax.viewLim.bounds
        xend = xstart + xdelta
        yend = ystart + ydelta
     
        #set margins given by time and frequency
        if xstart < self.tvec[0] or xend > self.tvec[-1]:
            xstart = max(xstart, self.tvec[0] )
            xend   = min(xend  , self.tvec[-1])
            _xstart,_xend = ax.get_xlim()
            if not (_xstart == xstart and _xend == xend):
                ax.set_xlim(xstart, xend)
         
        ymin = -.5/self.dt if self.signal.dtype == 'csingle' else 0
        if ystart < ymin or yend > .5/self.dt:
            ystart = max(ystart, ymin )
            yend   = min(yend  , .5/self.dt)
            ax.set_ylim(ystart, yend)

        if  xstart > xend: xstart, xend = xend, xstart

        self.set_yticks()

        # Update the image object with our new data and extent
        extent = self.image.im.get_extent()
        if self.data_plot is not None:
            self.data_plot.update_plot(xstart, xend)
        
        do_not_update = False
        if not force_update:
            #skip when nothing important was changed
            if self.method == 'time' and ystart >= extent[2] and \
                yend <= extent[3] and xstart==extent[0] and xend==extent[1]:
                do_not_update = True 
    
            if self.method == 'freq.' and ystart == extent[2] and \
                yend == extent[3] and xstart>=extent[0] and xend<=extent[1]:
                do_not_update = True  
            
            if self.method == 'sparse' and ystart == extent[2] and \
                yend == extent[3] and xstart==extent[0] and xend==extent[1]:
                do_not_update = True  

        if do_not_update:
            self.image.update_image_clim(xstart, xend, ystart, yend)
        else:
            x, y, z = self.__call__(xstart, xend, ystart, yend)
            self.image.update_image(xstart, xend, ystart, yend, z)

    def plot_ax_update(self, ax):
        xstart, xend = ax.get_xlim()
        if self.im_ax.get_xlim() != ax.get_xlim():
            self.im_ax.set_xlim(xstart, xend)
        

def extract_harmonics(tvec_data,t_range, f_range,n_harm_max,  cross_tvec_signal=None):
    
    from scipy.signal import get_window
    
    T = time.time()
    
    #not all signal must have teh same length and time vector
    fastest = np.argmax([(len(t)-1)/(t[-1]-t[0]) for t, s in tvec_data])
 
    nfft = next_fast_len(tvec_data[fastest][0].size)
    tvec = tvec_data[fastest][0]
    
    tvec = np.linspace(tvec[0], tvec[-1], tvec.size) #equally spaced vector!!

    
    cross_signal = None
    
    if cross_tvec_signal is not None:     
        cross_signal = np.interp(tvec, *cross_tvec_signal)
        win = get_window('hann', len(cross_signal), fftbins=False)   
        fdata0 = np.fft.rfft(cross_signal*win, nfft)[1:]
    
    #=================   use a stupid and robust method to get the amplitudes
    fvec = np.linspace(0, .5*nfft/np.diff(t_range)[0], nfft//2+1, endpoint=False)
    
    f0 = np.mean(f_range)
    df = np.diff(f_range)[0]
    
    if0 = fvec[1:].searchsorted(f0)
    idf  = fvec[1:].searchsorted(df)
    
    ind_f = [slice(*(if0*(n+1)+np.r_[-idf//2, idf//2+1])) for n in range(n_harm_max)]
    ind_f = [ind for ind in ind_f if np.size(fvec[ind])>0]
    n_harm = len(ind_f)
    
    fft_sig_n = [[] for n in range(n_harm)]
    fft_signoise = []
    
    

    data = []
    win = get_window('hann', len(tvec), fftbins=False)   
    for i, (t, sig) in enumerate(tvec_data):
        if t.size!= len(tvec):
            sig = np.interp(tvec, t, sig)
        
        data.append(sig)
        try:
            fdata = np.fft.rfft(sig*win, nfft)[1:]
        except Exception as e:
            print('fft error',e, sig.shape, win.shape, len(t), len(tvec))
            
        for n in range(n_harm):
            fft_sig_n[n].append(fdata[ind_f[n]])

        try:
            normalization = 1/np.mean(win)/len(win)*np.sqrt(ind_f[0].stop-ind_f[0].start)*4/3.
        except:
            print( 'normalization', ind_f, n_harm, if0, df, [slice(*(if0*(n+1)+np.r_[-df//2, (df+1)//2])) for n in range(n_harm)])
            raise

        fft_signoise.append(np.median(np.abs(fdata))*normalization)
 
    #absolute value of the perturbation 
    fft_signoise = np.array(fft_signoise)

    amplitude = []
    for n in range(n_harm):
        if cross_signal is None:
            A = np.array([np.linalg.norm(fs)/np.sum(win) for fs in fft_sig_n[n]])
            #substract contribution of the background random noise 
            A = np.sqrt(np.maximum(0, A**2-fft_signoise**2))
        else:
            X = fdata0[ind_f[n]]/np.linalg.norm(fdata0[ind_f[n]])
            A = np.array([np.abs(np.vdot(X, fs))/np.sum(win) for fs in fft_sig_n[n]])
            
        amplitude.append(A*3./2.)


    offset = np.array([np.asarray(sig).mean() for sig in data])
    corrupted = (fft_signoise == 0)#|(amplitude[0] < np.median(amplitude[0])/1e3)


    
    # use a advaced and more sensitive method. But only well coherent modes can be analyzed 

    cross_sig_num = None
    if cross_signal is None:
        #TODO use SVD to dinf the dominant component in the requested frequency rage? 
        #measure a cross phase with respect to the strongest signal
        cross_sig_num = np.argmax(amplitude[0]/(fft_signoise+np.nanmean(fft_signoise)*1e-6+corrupted))      
        cross_signal = data[cross_sig_num]
        

        
    dt = (tvec[-1]-tvec[0])/(len(tvec)-1)
    f_nq = .5/dt 
    cross_signal = detrend(cross_signal)
 
    # filter crosscorrelation signal in the selected frequency window
    fcross_signal = np.fft.rfft(cross_signal, nfft)
    
    df = f_nq/(nfft//2+1)
    ifmin = int(max(1, round(f_range[0]/df)))
    ifmax = int(min(len(cross_signal)//2+1, round(f_range[1]/df)))
    #sharp rectangular window in fourier space

    fcross_signal[:ifmin] = 0
    fcross_signal[ifmax:] = 0
    
    cross_signal_fband = np.fft.irfft(fcross_signal)[:len(cross_signal)]
    
    #get an envelope
    from scipy.signal import argrelmax
    i = argrelmax(np.abs(cross_signal_fband), mode='wrap')[0][1:-1]
    try:
        envelope = np.interp(tvec, tvec[i], np.abs(cross_signal_fband)[i])
    except:
        envelope = np.abs(cross_signal_fband)
    

    #normalize the aplitude
    b, a = butter(3, f0/f_nq*.5, 'lowpass')
    padlen = min(len(cross_signal)//2, 100)
    if not np.any(np.abs(np.roots(a))>=1):   envelope = filtfilt(b, a, envelope, padtype='even', padlen=padlen)

    cross_signal = np.copy(cross_signal_fband)
    cross_signal/= np.maximum(envelope, envelope.max()/20)[:len(cross_signal)]
    #signal whitening 
    cross_signal-= np.mean(cross_signal)
    cross_signal/= np.std(cross_signal)*np.sqrt(2)
    

    #get a 90 degrees shifted complex part 
    cmplx_sig = hilbert(cross_signal, nfft)[:tvec.size]

    amplitude2 = np.zeros((n_harm, len(data)))
    phi2       = np.zeros((n_harm, len(data)))
    #demodulate signal and find the first component

    
    complex_harm = np.zeros((n_harm, tvec.size), dtype='complex')
    complex_harm[0] = cmplx_sig
    for i in range(1, n_harm):
        complex_harm[i] = complex_harm[i-1]*cmplx_sig
    #make it more othogonal
    complex_harm, _ = np.linalg.qr(complex_harm.T)
            
    retro = np.zeros((tvec.size, len(data)))
    error = np.zeros(len(data))
    
    #remove issues with boundary effects
    win = get_window('hann', tvec.size, fftbins=False).astype('single')   
    win /= np.sum(win)/np.sqrt(len(win))

    for j, sig in enumerate(data):
        mask = np.ma.getmask(sig)
        #skip if more than 20% is corrupted
        if np.any(mask) and np.sum(mask)/len(mask) > 0.2:
           continue
  
        #perform a least squares fit of the orthonormal complex prototype to the measured signal
        sig_ = win*(sig-offset[j])
        #harm = 2*np.dot(win*(sig-offset[j]), np.conj(complex_harm))   
        harm = np.array([2*(sig_*np.conj(ch)).sum() for ch in complex_harm.T])
        amplitude2[:, j] = np.abs(harm)
        phi2[:, j] = np.angle(harm)
        #retrofit of the measured signal 
        retro[:, j] = np.dot(complex_harm, harm*np.sqrt(len(win))).real
        retro[:, j] += offset[j]
        #this might be wrong if there is another non-random signal
        error[j] = (sig-retro[:, j]).std()
    
    error *= np.sum(win**2)/np.sqrt(len(tvec))
        
    corrupted |= amplitude2[0] == 0   
    
    return n_harm, amplitude, tvec, cross_signal, amplitude2, phi2, retro, offset, error, corrupted, cross_sig_num


class RadialViewer:

    initialized = False
    #plot radial profile of the aplitude and phase of the diagnostic

    def __init__(self, parent, fig, n_harm, rho_lbl):

        self.parent = parent
        self.fig = fig
        self.fig.clf()
        self.n_harm_max = n_harm
        self.group = None
        self.rho_lbl = rho_lbl
        self.diag = None
        
        self.ax1 = self.fig.add_subplot(211)
        self.ax2 = self.fig.add_subplot(212, sharex= self.ax1)

        self.fig.subplots_adjust(hspace=0.05, wspace = 0.10)
        self.plot1 = []
        self.plot2 = []
        self.plot1_ = []
        self.plot2_ = []
        self.plot2_s = []

        c = 'b', 'g', 'r', 'y', 'm'

        for i in range(self.n_harm_max):
            self.plot1.append(self.ax1.plot([], [], 'o--', mfc='none', mec=c[i%len(c)]
                            , zorder=self.n_harm_max-i)[0])
            self.plot1_.append(self.ax1.errorbar(0, 0, 0, fmt='s-'+c[i%len(c)], 
                            capsize=0, zorder=self.n_harm_max-i, 
                            label='%d. harmonic'%(i+1)))
            self.plot2_.append(self.ax2.plot([], [], c=c[i%len(c)], zorder=self.n_harm_max-i)[0])
            self.plot2_s.append(self.ax2.scatter([], [], edgecolors=c[i%len(c)], 
                            s=50, zorder=99-i, linewidths=.5))
            
        #steady state profile
        self.plt_harm0, = self.ax1.plot([], [], '-k')
            
        self.m = 1  #mode m-number 

        self.plot2_theta, = self.ax2.plot([], [], 'kx-', zorder=99, label=r'Ideal m=%d mode'%self.m)
        #self.plot1_map_ampl, = self.ax1.plot([], [], 'k:')  #mapping correction factor

        self.leg1 = self.ax1.legend(loc='best', fancybox=True, prop={'size':12})
        self.leg1.get_frame().set_alpha(0.7)
        
        self.leg2  = self.ax2.legend(loc='best', fancybox=True, prop={'size':12})
        self.leg2.get_frame().set_alpha(0.7)
        
        self.q_surf_names = '1/1', '4/3', '3/2', '2/1'
        self.q_surf_vals = 1, 4/3., 3/2., 2

        nsurf = len(self.q_surf_names)
        
        self.plot_qsurf1 = [self.ax1.axvline(x=np.nan, c='k', ls='--', lw=.5, dashes=(5, 5)) for i in range(nsurf)]
        self.plot_qsurf2 = [self.ax2.axvline(x=np.nan, c='k', ls='--', lw=.5, dashes=(5, 5)) for i in range(nsurf)]
        self.text_qsurf  = [self.ax2.text(10, 10, '', fontsize=8, rotation=90, zorder=99, backgroundcolor='none') for i in range(nsurf)]
        for txt in self.text_qsurf: txt.set_clip_on(True)
 
        self.amplitude = []
        self.phi = None
        
        self.rho = None
        self.t_range = -np.infty, np.infty
        self.f_range = None

        self.ax1.set_ylabel('Perturbation $\delta$', fontsize=font_size)
        self.ax2.set_ylabel('Cross phase $\phi$', fontsize=font_size, labelpad=-4)
        rho_lbl = self.rho_lbl
        if self.rho_lbl[:3] == 'rho': rho_lbl = '\\'+rho_lbl[:4]+'{'+rho_lbl[4:]+'}'
        rho_lbl = '$'+ rho_lbl  +'$'
        
        self.ax2.set_xlabel(rho_lbl, fontsize=font_size+4)
        
        y_tick = np.linspace(-np.pi*2, np.pi*2, 9)
        y_label = [r"$-2\pi$", r"$-3\pi/2$", r"$-\pi$", r"$-\pi/2$", "$0$", 
                   r"$\pi/2$", r"$\pi$", r"$3\pi/2$", r"$2\pi$"]
        
        self.ax2.set_yticks(y_tick)
        self.ax2.set_yticklabels(y_label, fontsize=18) 
        self.multi = MultiCursor(self.fig.canvas, (self.ax1, self.ax2), color='k', lw=1)

        self.ax2.grid(True)
        self.ax1.grid(True)

        self.ax1.ticklabel_format(scilimits=(-2, 2), axis='y')
        self.ax1.set_xlim(-1, 1)
        self.ax2.set_ylim(-4, 4)
        
        for label in ( self.ax1.get_yticklabels()):
            label.set_fontsize(font_size) # Size here overrides font_prop
        self.ax1.xaxis.offsetText.set_visible(False)
        for label in (self.ax1.get_xticklabels()):
            label.set_visible(False)
            
        for label in (self.ax2.get_xticklabels() + self.ax2.get_yticklabels()):
            label.set_fontsize(font_size) # Size here overrides font_prop
        
        self.plot_description = self.ax2.text(1.008, .05, '', rotation='vertical', 
                        transform=self.ax2.transAxes, verticalalignment='bottom', 
                        size='xx-small', backgroundcolor='none')

    def reset(self):
        #reset when the shot number is changed
        if self.initialized:
            self.fig.clf()
            self.__init__( self.parent, self.fig, self.n_harm_max, self.rho_lbl)
        
    def set_data(self, diag, rho, theta0, signal, R, Z, Phi, mag_coord_diag, mag_coord_data, qsurfs, units):

        rho[np.isnan(rho)|(R<1)] = 2
        self.qsurfs = qsurfs
        self.diag = diag
        
        rho_ind = np.argsort(rho)
        self.rho = rho[rho_ind]
        
        self.data = [signal[i] for i in rho_ind]
        self.cross_signal = None
        self.initialized = True
        
        self.theta_diag = mag_coord_diag['theta_geo'][rho_ind]
        self.theta_star_diag = mag_coord_diag['theta_star'][rho_ind]
 
        self.R = R[rho_ind]
        self.Z = Z[rho_ind]
        self.theta0 = theta0[rho_ind]
        self.Phi0 = Phi #toroidal angle [rad], assume a single number?

        self.units = units

    def set_cross_signal(self, cross_signal):
        self.cross_signal =  cross_signal
        
    def change_m_number(self, ind):

        m = self.parent.m_numbers[ind]
        if m != self.m:
            self.m = m
            self.update()
            self.parent.tomo_m_num.setCurrentIndex(ind)
            self.leg2.get_texts()[0].set_text(r'Ideal m=%d mode'%self.m)

    def update(self):

        description = self.description+'  from %.2f to %.2fs, from %.1f to %.1fkHz'%(self.t_range[0], 
                                    self.t_range[1], self.f_range[0]/1e3, self.f_range[1]/1e3)
        self.plot_description.set_text(description)
        
        
        out = extract_harmonics(self.data,self.t_range,self.f_range,
                                self.n_harm_max,  self.cross_signal)

        n_harm, amplitude, self.tvec, cross_signal, self.amplitude2, self.phi2, self.retro, offset, error, corrupted, cross_sig_num = out

        #map poloidal theta of the diagnostic to the theta_star (WARNING theta star from CLISTE is not very accurate!!)
        theta_star0 = np.zeros_like(self.theta0)
        theta0 = np.unwrap(self.theta0, discont=1.5*np.pi)
        for i, (t0, t, ts) in enumerate(zip(theta0, self.theta_diag, self.theta_star_diag)):
            theta_star0[i] = np.interp(t0, np.r_[t-2*np.pi, t, t+2*np.pi], np.r_[ts-2*np.pi, ts, ts+2*np.pi])*self.m #use periodicity 
        theta_star0 -= np.nanmean(theta_star0)

        #make a nice plot of phi, if possible  without 2pi jumps
        ind = ~corrupted&(np.abs(self.rho)<1)

        
        try:
            phi = np.copy(self.phi2[0][ind])
            #keep  pi jumps sign consistent with theta_star0
            phi = np.unwrap(phi-theta_star0[ind])+theta_star0[ind]

            #rought shift with respect to the strongest
            largest_values = amplitude[0][ind] > mquantiles(amplitude[0][ind], .9)/2
            if not any(largest_values):
                largest_values = np.argmax(amplitude[0][ind])
            phi -= np.median((phi-theta_star0[ind])[largest_values]) 
            #shift with respect to theta_star0, ignoring pi jumps 
            weights = amplitude[0][ind]*(1-np.abs(self.rho[ind]))
            phi-= np.average((phi-theta_star0[ind]+np.pi/2)%np.pi - np.pi/2, weights=weights)
            #shift it by k*pi with respect to the strongest channels
            phi -= np.pi*round(np.average((phi-theta_star0[ind]), weights=weights)/np.pi)
        except:
            
            traceback.format_exc()
            print( 'phi, weights', phi, self.phi2[0][ind], ind, amplitude[0][ind])
     

        self.phi2[0][ind] = phi 
        self.amplitude2[:, corrupted] = np.nan
        theta_star0[np.abs(self.rho)>=1.1] = np.nan
        offset[corrupted] = np.nan
        self.phi2[0][corrupted] = np.nan 

        #plot amplitude
        for n in range(n_harm):
            amplitude[n][corrupted] = np.nan
            plotline, caplines, barlinecols =  self.plot1_[n]
            # Replot the data first
            x, y = self.rho, self.amplitude2[n]
            plotline.set_data(x, y)
            #yerr = np.array([np.std(sig)/np.sqrt(len(tvec_)/2.) for tvec_, sig in self.data])
            yerr = error#[n]
            # Find the ending points of the errorbars
            error_positions = (x, y-yerr), (x, y+yerr)
            # Update the error bars
            barlinecols[0].set_segments(zip(zip(x, y-yerr), zip(x, y+yerr)))

        for n in range(n_harm, len(self.plot1_)):
            plotline, caplines, barlinecols =  self.plot1_[n]
            plotline.set_data(np.nan, np.nan)
            barlinecols[0].set_segments(zip(zip(x*np.nan, x*np.nan), zip(x*np.nan, x*np.nan)))
            
        self.plot1[0].set_data(self.rho, amplitude[0])
        for n in range(1, len(self.plot1)):
            self.plot1[n].set_data(np.nan, np.nan)

        self.plt_harm0.set_data(self.rho, offset)

        #plot a phase pof ECE diag
        self.plot2_theta.set_data(self.rho, theta_star0)
        
        #plot a phase of the first harmonic
        n = 0
        self.plot2_[n].set_data(self.rho, self.phi2[n])
        c = self.plot2_s[n].get_edgecolor()
        w = self.amplitude2[n]/np.nanmax(self.amplitude2[n])
        w[np.isnan(w)] = 0
        self.plot2_s[n].set_facecolor((c.T*w+1-w).T)
        self.plot2_s[n].set_offsets(np.c_[self.rho, self.phi2[n]])

        ymax = 0
        for n in range(n_harm):
            ii = ind if any(ind) else  ~corrupted
            ymax = max(ymax, max(amplitude[n][ii]))
            ymax = max(ymax, max(self.amplitude2[n][ii]))

        #plot expected positions of the resonance surfaces 
        for i, name in enumerate(self.q_surf_names):
            self.plot_qsurf1[i].set_xdata(self.qsurfs[i])
            self.plot_qsurf2[i].set_xdata(self.qsurfs[i]) 
            self.text_qsurf[i].set_text('$q=%s$'%name)
            if np.isnan(self.qsurfs[i]):
                self.text_qsurf[i].set_visible(False)
            else:
                self.text_qsurf[i].set_visible(True)
                self.text_qsurf[i].set_x(self.qsurfs[i])
            f = .2
            y1, y2 = self.ax2.get_ylim()        
            self.text_qsurf[i].set_y(f*y2+(1-f)*y1)


        self.ax1.set_ylim(0, 1.1*ymax )
        self.ax1.figure.canvas.draw_idle()
        self.ax2.figure.canvas.draw_idle()


class Diag2DMapping(object):

    keyCtrl = False
    keyShift = False
    initialized = False
    loaded = False
    n_theta = 10
    n_rho = 10
    n_harm = 3
    ECH_gyrotron_locations = {}
    initialise2D = False
    cross_signal = None

    def __init__(self, parent, fig, n_contour, remback_button, show_ecei):

        self.parent = parent
        self.fig = fig
        self.fig.subplots_adjust(right=0.82, top=.95)
        self.ax = self.fig.add_subplot(111, label='2Dmap')
        self.plot_ece_active, = self.ax.plot([], [], 'wo', zorder=100, markeredgecolor='k')
        self.plot_ece_ignored, = self.ax.plot([], [], 'o' , zorder=100, mfc='none', markeredgecolor='k')
        self.ech_location, = self.ax.plot([], [], 'xk', zorder=100)

        self.plot_ecei, = self.ax.plot([], [], 'k.' , zorder=100, markeredgecolor='k')
        self.plot_ecei_region, = self.ax.plot([], [], 'k-' , zorder=100, lw=.5)

        X = np.zeros((1, self.n_rho))
        self.plot_mag_rho   = self.ax.plot(X, X, 'k--', lw=.5, zorder=99, dashes=(5, 5))
        X = np.zeros((1, self.n_theta))
        self.plot_mag_theta = self.ax.plot(X, X, 'k--', lw=.5, zorder=99, dashes=(5, 5))
        
        self.n_contour = n_contour
        self.remback_button = remback_button
        self.show_ecei = show_ecei
        
        #load shape of the tokamak walls
        if self.parent.tokamak == "AUG":
            import aug_sfutils as sf
            gc_d = sf.getgc()
            for gcc in gc_d.values():
                self.ax.plot(gcc.r, gcc.z, 'k', lw=.5)
        elif self.parent.tokamak == "DIIID":
            from loaders_DIIID import map_equ
            gc_r, gc_z = map_equ.get_gc()
        else:
            gc_r, gc_z =  {},{}
        try:
            for key in gc_r:
                self.ax.plot(gc_r[key], gc_z[key], 'k', lw=.5)
        except:
            logger.error('Loading of the  vessel shape failed!!')
            logger.error( traceback.format_exc())


        for label in (self.ax.get_xticklabels() + self.ax.get_yticklabels()):
            label.set_fontsize(font_size) # Size here overrides font_prop

        self.cbar_ax = self.fig.add_axes([0.85, 0.1, 0.04, 0.85], label='cbar')
        self.cbar_ax.xaxis.set_major_formatter(NullFormatter())
        self.cbar_ax.yaxis.set_major_formatter(NullFormatter())
        self.cbar_ax.tick_params(labelsize=font_size) 

        self.plot_description = self.ax.text(1.008, .05, '', rotation='vertical', 
                transform=self.ax.transAxes, verticalalignment='bottom', 
                size='xx-small', backgroundcolor='none')

        self.ax.axis('equal')
        self.ax.set_xlabel('R [m]', fontsize=font_size)
        self.ax.set_ylabel('z [m]', fontsize=font_size)
        self.ax.axis([1.4, 2, -0.25, 0.35]) #BUG should not be fixed
        self.shift_phi = 0
        self.m = 1
        self.lim = 0.5
        self.substract = False
        self.dR = 0
        self.dZ = 0
        self.sxr_tvec = None
        self.sxr_emiss = None
        
        self.t_range = -np.infty, np.infty
        self.f_range = None


        self.cid = self.fig.canvas.mpl_connect('scroll_event', self.MouseWheelInteraction)
        self.cid2 = fig.canvas.mpl_connect('key_press_event',  self.onKeyPress)
        self.cid3 = fig.canvas.mpl_connect('key_release_event', self.onKeyRelease)
    
    def reset(self):
        if self.initialized:
            self.fig.clf()
            self.__init__( self.parent, self.fig, self.n_contour, self.remback_button, self.show_ecei)

            
    def set_data(self,rho,R, Z, theta0,  Phi, data, mag_coord_diag, mag_coord_data,
                 units, n_harm,  shot, use_LFS_data, phase_locked):
          
          
        self.phase_locked = phase_locked
        self.units = units
        self.initialized = False
        self.shot = shot
        self.phi0 = 0 #used for phase locking of the mode 
        self.n_harm = n_harm
        self.Phi0 = Phi

 
        out = extract_harmonics(data,self.t_range, self.f_range,n_harm, cross_tvec_signal=self.cross_signal)
        n_harm, amplitude, tvec, cross_signal, amplitude2, phi2, retro, offset, error, corrupted, cross_sig_num = out
     
        
        rho[np.isnan(rho)|(R<1)] = 2   
        #nlfs = np.sum((np.abs(rho)<=1)&(rho>=0)&np.isfinite(amplitude2[0]))
        #nhfs = np.sum((np.abs(rho)<=1)&(rho<=0)&np.isfinite(amplitude2[0]))
        
        #rho_sign = 1 if nhfs <= nlfs else -1
        rho_sign = 1 if use_LFS_data else -1
            
        ind = (np.abs(rho)<=1)&(rho_sign*rho>=0)&np.isfinite(amplitude2[0])
        
        self.plot_ece_active.set_data( R[ ind], Z[ ind])
        self.plot_ece_ignored.set_data(R[~ind], Z[~ind])

        sind = np.argsort(np.abs(rho))
        ind = sind[ind[sind]]
        
        assert len(ind) > 2, 'Not enough signals inside of separatrix to make the mapping!'  
    
        self.tvec = tvec
        self.rho = rho[ind]
        self.retro =  retro[:, ind]
        self.theta0 = theta0[ind]%(2*np.pi)
        self.theta = mag_coord_diag['theta_geo'][ind]
        self.theta_star = mag_coord_diag['theta_star'][ind]
        
    
   
        #isoflux and isotheta surfaces, just for plotting 
        self.isoflux_R = mag_coord_data['R'][self.n_rho-1::self.n_rho]
        self.isoflux_Z = mag_coord_data['Z'][self.n_rho-1::self.n_rho]
        
        t = np.linspace(0, 2*np.pi, self.n_theta, endpoint=False)
        self.isotheta_R = np.array([np.interp(t, ts, r) for ts, r in zip(mag_coord_data['theta_star'], mag_coord_data['R'])])
        self.isotheta_Z = np.array([np.interp(t, ts, z) for ts, z in zip(mag_coord_data['theta_star'], mag_coord_data['Z'])])

        self.ax.axis([self.isoflux_R[-1].min(), self.isoflux_R[-1].max(), 
                      self.isoflux_Z[-1].min(), self.isoflux_Z[-1].max()]) #BUG should not be fixed
  
 
        Rmag_ece = mag_coord_diag['R'][ind]
        Zmag_ece = mag_coord_diag['Z'][ind]
        
        n_theta = Rmag_ece.shape[1]
        
        #add point on axis
        self.Rmag_ece = np.c_[Rmag_ece[0].mean()*np.ones(n_theta), Rmag_ece.T].T
        self.Zmag_ece = np.c_[Zmag_ece[0].mean()*np.ones(n_theta), Zmag_ece.T].T
     
        #load sxr data if availible and and ECH resonance location
        self.load_sxr()
        self.load_ech_resonances()
        
        self.loaded = True

    def set_data2D(self,R, Z, Phi,  data):
        #plot data in 2D grid 
 
        out = extract_harmonics(data, self.t_range, self.f_range, self.n_harm,
                              cross_tvec_signal=self.cross_signal)
        n_harm, amplitude, tvec, cross_signal, amplitude2, phi2, \
                          retro, offset, error, corrupted, cross_sig_num = out
        
                
        nrad = 8
        npol = R.size//nrad
        self.retro2D = retro.reshape(len(tvec), npol,nrad)
        self.invalid2D = corrupted.reshape(npol, nrad)
        self.offset2D = offset.reshape(npol, nrad)
       

        self.plot_ecei.set_data(R[~corrupted],Z[~corrupted])


        R = R.reshape(npol, nrad)
        Z = Z.reshape(npol, nrad)

        self.plot_ecei_region.set_data(np.r_[R[:,0],R[::-1,-1], R[0,0]],
                                       np.r_[Z[:,0],Z[::-1,-1], Z[0,0]])
        
        self.tvec2D = tvec
        self.coord2D = R, Z, Phi
        self.initialise2D = True
        
        

        
    def set_cross_signal(self, cross_signal):
        self.cross_signal =  cross_signal
        
    def change_m_number(self, ind):

        m = self.parent.m_numbers[ind]
        if m != self.m:
            self.m = m
            self.update()
            self.parent.tomo_m_num.setCurrentIndex(ind)
      
        
    def load_sxr(self):
        try:
            paths = ['Emissivity_%d.npz'%self.shot, 
                     os.path.expanduser('~/tomography/tmp/Emissivity_%d.npz'%self.shot),  ]
            for path in paths:
                if os.path.isfile(path):
                    break

            emiss = load(path, allow_pickle=True)

            if self.tvec[-1] > emiss['tvec'][-1] or  self.tvec[0] < emiss['tvec'][0]:
                raise Exception('out of SXR tomography %.3f - %.3fs'%(emiss['tvec'][0], emiss['tvec'][-1]))
                        
            if  self.sxr_emiss is None or not all(emiss['tvec'] == self.sxr_tvec):
                self.sxr_tvec = emiss['tvec']
                print( 'loading SXR from %.4f to %.4f'%(self.sxr_tvec[0], self.sxr_tvec[-1]))
                self.sxr_emiss = np.single(emiss['gres'])*emiss['gres_norm'][None, None, :]
                self.sxr_r = emiss['rvec']
                self.sxr_z = emiss['zvec']

        except Exception as e:
            self.sxr_emiss = None
    
    def load_ech_resonances(self):
         
        #plot ECH heating locations, only for DIII-D
        if self.shot not in self.ECH_gyrotron_locations:
            self.ECH_gyrotron_locations[self.shot] = []
            try:
                logger.info('Fetching ECH location' )
                MDSconn = self.parent.MDSconn
                MDSconn.openTree('AOT', self.shot) 
                for i in range(1,7):
                   try:         
                       RECH = MDSconn.get(r'_x=\AOT::TOP.TORAY.TORAY%d.PEAK:R'%i).data()
                       ZECH = MDSconn.get(r'_x=\AOT::TOP.TORAY.TORAY%d.PEAK:Z'%i).data()
                       time = MDSconn.get('dim_of(_x)').data() / 1e3
                   except:
                      continue
    
                   self.ECH_gyrotron_locations[self.shot].append([time, RECH, ZECH])
            except Exception as e:
                pass
            
        if len(self.ECH_gyrotron_locations.get(self.shot, [])):
            #plot the ECH locations 
            R,Z = [], []
            for t,r,z in self.ECH_gyrotron_locations[self.shot]:
                it = np.argmin(np.abs(t-np.mean(self.t_range)))
                R.append(r[it])
                Z.append(z[it])
            logger.info( 'Replot ECH location')
            self.ech_location.set_data(R,Z)
 
        
    def prepare(self):
        if not self.loaded: return 

        if self.substract:
            self.cmap = 'seismic'
        else:
            from copy import copy
            self.cmap = copy(plt.cm.get_cmap('nipy_spectral'))
            self.cmap._init()
            self.cmap.set_under('w')

        dT = 1/np.array(self.f_range)[::-1]
        #select a time interval used for mapping on theta_star
        i_start = self.tvec.searchsorted(self.tvec[self.tvec.size//2])#-0.5*mean(dT))
        self.t_start = self.tvec[i_start]
        r0 = self.retro[i_start]
        ind_tend = slice(*self.tvec.searchsorted(self.t_start+dT))
        i_end = ind_tend.start+np.argmin(np.sum((self.retro[ind_tend]-r0[None, :])**2, 1))
        
        #second pass to get a m-th minimum more accurately
        dT = (self.tvec[i_end]-self.t_start)*np.abs(self.m)
        dT = dT*np.array((1-0.5/np.abs(self.m), 1+0.5/np.abs(self.m)))
        ind_tend = slice(*self.tvec.searchsorted(self.t_start+dT))
        i_end = ind_tend.start+np.argmin(np.sum((self.retro[ind_tend]-r0[None, :])**2, 1))
        ind_t = slice(i_start, i_end+1)

        self.mode_Te = self.retro[ind_t]
        self.f0 = 1/(self.tvec[ind_t][-1]-self.tvec[ind_t][0])*np.abs(self.m)

        #position of the measurements in the theta star coordinates
        self.thetaStar0 = [np.interp(t0, t, ts) for t0, t, ts in zip(self.theta0, self.theta, self.theta_star)]
        self.thetaStar0 = (np.array(self.thetaStar0) + np.pi)%(2*np.pi) - np.pi

        #angular shift due to the time 
        self.phi = np.linspace(0, 2*np.pi, len(self.mode_Te), endpoint=True)
        
        #calculate vmin,vmax, levels
        self.UpdateLim(None)

        if self.substract: 
            mode_Te = self.mode_Te-self.mode_Te.mean(1)[:, None]

        #compensate the phase shift of the Te profile in a different timepoints
        #for a better visualization of the rotating mode 
        if self.phase_locked:
            if hasattr(self, 'rmax'):
                imax = np.argmin(bp.abs(self.rho-self.rmax))  #preselected position
            else:
                imax = np.argmax(np.std(self.mode_Te, 0))

            m  = np.mean(self.mode_Te[:, imax])
            Ac = np.sum((self.mode_Te[:, imax]-m)*np.cos(self.phi*self.m))
            As = np.sum((self.mode_Te[:, imax]-m)*np.sin(self.phi*self.m))
            self.phi0 = np.arctan2(As, Ac)/np.abs(self.m)
            print(' ECE phase shifted by %ddeg'%np.rad2deg(self.phi0 ))

        self.initialized = True


    def shift_phase(self, shift_phi):

        self.shift_phi = shift_phi%(2*np.pi)#periodicity -         
        self.time = self.t_start + self.shift_phi/(2*np.pi*self.f0)*np.abs(self.m)
        try:
            description = '#%d  at %.6fs, $\phi_0$ = %d, f$_0$ = %.4fkHz and m=%d'%(self.shot, \
                self.time, np.rad2deg(np.median(self.Phi0)), self.f0/1e3, self.m)
        except:
            description = ''
            print('Error: ', (self.shot, self.time, np.rad2deg(np.median(self.Phi0)), self.f0/1e3, self.m, shift_phi, self.t_start))
        self.plot_description.set_text(description)


    def shiftTe(self):
        nr, nt = self.theta_star.shape
        mode_Te = np.zeros((nr+1, nt))
        mode_Te[0] = self.mode_Te[:, 0].mean()  #fill the core by zero aplitude mode
        shift_phi = self.shift_phi + self.phi0
        
        #iterate over radial positions of the measurements
        for i, (t0, te, t) in enumerate(zip(self.thetaStar0, self.mode_Te.T, self.theta_star)):
            mode_Te[i+1] = np.interp((np.sign(self.m)*(t-t0)+shift_phi)%(2*np.pi), self.phi, te)

        if self.substract: 
            mode_Te -= mode_Te.mean(1)[:, None]

        return mode_Te
        
    
    def clear_contours(self, ax=None):
        #use self.ax if no other ax is provided
        if ax is None: ax = self.ax
        
        #remove contours from a previous plotting 
        for element in ('Te_contourf', 'Te_contour', 'Te2D_contour', 'Te2D_contourf'):
            if hasattr(ax, element):
                for coll in getattr(ax, element).collections:
                    try:    ax.collections.remove(coll)
                    except ValueError: pass#Everything is not removed for some reason!    

    
    def update(self, update_cax=True, update_mag=True, animate=False, ax = None, filled_contours = True):

        if not self.initialized: return 
 
       
        #use self.ax if no other ax is provided
        if ax is None: ax = self.ax
        self.clear_contours(ax=ax)
        self.shift_phase(self.shift_phi)

 
        if update_mag:
            for i, p in enumerate(self.plot_mag_rho):
                p.set_data(self.isoflux_R[i], self.isoflux_Z[i])
            for i, p in enumerate(self.plot_mag_theta):
                p.set_data(self.isotheta_R[:, i], self.isotheta_Z[:, i])

        collections = (self.plot_description, )
        
        
        #plot SXR data if availible 
        if self.sxr_emiss is not None and filled_contours:
            R, Z = np.meshgrid(self.sxr_r+self.dR, self.sxr_z+self.dZ)  #BUG!!!
            n = (0, 1, 1, 2, 3, 3)[np.abs(self.m)]  #BUG just guess of the most common mode!!
            dT = 0
            if self.parent.tokamak == "AUG":
                dT = (-40.5/360.)/self.f0*n/self.m  #BUG hardcoded toroidal shift between SXR and ECE
            if self.parent.tokamak == "DIIID":
                dT = (69/360.)/self.f0*n/self.m  #BUG hardcoded toroidal shift between SXR and ECE

            it = np.argmin(np.abs(self.sxr_tvec+dT-self.time))
            E = np.maximum(self.sxr_emiss[:, :, it], 0)
  
            i1, i2 = self.sxr_tvec.searchsorted(self.t_start + np.r_[0, 1]/self.f0*np.abs(self.m))
            
            vmin_sxr = 0
            vmax_sxr = mquantiles(self.sxr_emiss[:, :, i1:i2].max((0, 1)), 0.95)[0]
      
            self.units = 'W/m$^3$'
            prefix = 'k'
            
            levels_sxr=np.linspace(vmin_sxr, vmax_sxr, self.n_contour)
            self.SXR_contour = ax.contourf(R, Z, E/1e3, levels_sxr/1e3, cmap='nipy_spectral',
                                        vmin=vmin_sxr/1e3, vmax=vmax_sxr/1e3, extend='max')
            collections += (self.SXR_contour.collections, )
            
            
            
        mode_Te = self.shiftTe()/self.fact

        if not filled_contours or self.sxr_emiss is not None:
            if self.substract:
                levels=np.linspace(self.vmin, self.vmax, 15)
            ax.Te_contour = ax.contour(self.Rmag_ece, self.Zmag_ece, mode_Te, extend='both', linewidths = .5, 
                            vmin=self.vmin, vmax=self.vmax, levels=levels, colors='k')
            collections += (ax.Te_contour.collections, )

        else:
            ax.Te_contourf = ax.contourf(self.Rmag_ece, self.Zmag_ece, mode_Te, extend='both', 
                            vmin=self.vmin, vmax=self.vmax, levels=self.levels, cmap=self.cmap)
            
            collections += (ax.Te_contourf.collections, )
            if not self.substract:
                ax.Te_contour = ax.contour(ax.Te_contourf, colors='k', linewidths=0.3 )
                collections += (ax.Te_contour.collections, )

            if update_cax:
                #BUG update colorbar by creating a new one :( 
                self.cbar_ax.cla()
                cb = self.fig.colorbar(ax.Te_contourf, cax=self.cbar_ax )
                tick_locator = MaxNLocator(nbins=7)
                cb.locator = tick_locator
                cb.update_ticks()
                cb.set_label('T$_e$ ['+self.prefix+self.units+']', labelpad=4, fontsize=font_size)


        #plot ECEI data
        if self.initialise2D and self.show_ecei.isChecked():
            self.plot_ecei.set_visible(True)
            self.plot_ecei_region.set_visible(True) 


            R,Z,Phi = self.coord2D
            #BUG just guess of the most common mode!!
            n = (0, 1, 1, 2, 3, 3)[np.abs(self.m)]  

            #toroidal shift betwen ECE and ECEI results in some time shift between the data            
            dPhi = -Phi + self.Phi0    
            dT = dPhi/(self.f0*2*np.pi)*n/self.m  
            t = time.time()

            Te2D = interp1d(self.tvec2D, self.retro2D, axis=0,
                              copy=False, assume_sorted=True)(self.time-dT)
            #create masked array 
            Te2D = Te2D/self.fact
            if self.substract: 
                Te2D -= self.offset2D/self.fact

            #interpolate corrupted points
            if np.any(self.invalid2D):
                from scipy.interpolate import griddata
                Te2D[self.invalid2D] = griddata((R[~self.invalid2D], Z[~self.invalid2D]),
                                           Te2D[~self.invalid2D], (R[self.invalid2D],Z[self.invalid2D]))
       
    
            ax.Te2D_contourf = ax.contourf(R, Z, Te2D, extend='both', zorder = 10,
                                    vmin=self.vmin, vmax=self.vmax, 
                                    levels=self.levels, cmap=self.cmap)
            
            collections += (ax.Te2D_contourf.collections, )
            
            if not self.substract:
                ax.Te2D_contour = ax.contour(ax.Te2D_contourf, colors='k', 
                                            linewidths=0.3,zorder = 10,)
                collections += (ax.Te2D_contour.collections,)
        
        self.plot_ecei.set_visible(self.show_ecei.isChecked())
        self.plot_ecei_region.set_visible(self.show_ecei.isChecked()) 

                
        
        if not animate:
            self.fig.canvas.draw_idle()

        return collections


    def MouseWheelInteraction(self, event):

        if event.inaxes == self.ax:
            steps = 360 if self.keyCtrl else 36
            self.shift_phase(self.shift_phi + event.step*2*np.pi/steps)
            self.update(update_cax=False, update_mag=False)

    def onKeyPress(self, event):
        steps = 36
        if 'control' == event.key:
            self.keyCtrl=True
            steps = 360

        if 'shift' == event.key:
            self.keyShift=True
            
        if 'left' == event.key:
            self.shift_phase(self.shift_phi - 2*np.pi/steps)
            self.update(update_cax=False, update_mag=False)
            
        if 'right' == event.key:
            self.shift_phase(self.shift_phi + 2*np.pi/steps)
            self.update(update_cax=False, update_mag=False)

    def onKeyRelease(self, event):
        if 'control' == event.key:
            self.keyCtrl=False
        if 'shift' == event.key:
            self.keyShift=False

    def UpdateModeM(self, ind, animate=False):
        m = self.parent.m_numbers[ind]
        if m != self.m:
            self.m = m
            self.prepare()
            self.update(update_cax=False, animate=animate)
            self.parent.radial_m_num.setCurrentIndex(ind)

    def UpdateLim(self, ind=None):
        if ind is not  None:
            self.lim = ind/100.
            
        mode_Te = self.shiftTe()

        
        if np.abs(mode_Te).max() < 1e3: 
            self.prefix, self.fact = '' , 1e0
        elif np.abs(mode_Te).max() < 1e6:
            self.prefix, self.fact = 'k', 1e3
        else:
            self.prefix, self.fact = 'M', 1e6

        vmax = mode_Te.max()/self.fact
        vmin = mode_Te.min()/self.fact

        if self.substract:
            self.vmax = max(vmax, -vmin)
            self.vmin = -self.vmax
            n_contour = self.n_contour*2
        else:
            self.vmin = (vmax-vmin)*self.lim+vmin
            self.vmax = vmax 
            n_contour = self.n_contour

        self.levels = np.linspace(self.vmin, self.vmax, n_contour)
        
        
        self.update()
    
    def UpdatePlotType(self, ind):
        self.substract = self.remback_button.isChecked()
        self.prepare()
        self.update()

        
class GetInfoThread(QThread):
    
    finished = pyqtSignal(str)

    def __init__(self, parent, fun, *args):
        QThread.__init__(self, parent)
        self.parent = parent
        self.fun = fun
        self.args = args

    def run(self):
        while not self.parent.eqm_ready:
            time.sleep(0.1)
        
        try:
            info = self.fun(*self.args)
        except Exception as e:
            print( 'info err', e, self.fun)
            print( traceback.format_exc())
            info = ''
        self.finished.emit(info)


class MainGUI(QMainWindow):

    def __init__(self, shot=None, diag=None, group=None, signal=None, diag_phase=None, signal_phase=None, tmin=np.nan, tmax=np.nan, fmin=0, fmax=np.inf, parent=None ):
        QMainWindow.__init__(self, parent)
        
        self.setWindowTitle('Data preprocessing tool')
        
        self.shot = shot        
        self.diag = diag
        self.diag_group = group
        self.signal = signal
        self.diag_phase = diag_phase
        self.signal_phase = signal_phase

        self.data_loader = None
        self.data_loader_radial = None
        self.data_loader_phase = None
        self.signal_radial = None
        
        self.load_config()
        self.parent = parent
        self.MDSconn = None
        if self.tokamak == "AUG":
            #import map_equ_20200306 as map_equ
            #self.eqm = map_equ.equ_map()
            pass
        elif self.tokamak == "DIIID":
            from loaders_DIIID import map_equ
            import MDSplus as mds
            self.MDSconn = mds.Connection(self.mds_server )
            self.eqm = map_equ.equ_map(self.MDSconn)
        elif self.tokamak == "NSTX":
            from loaders_NSTX import map_equ
            import MDSplus as mds
            self.MDSconn = mds.Connection(self.mds_server )
            self.eqm = map_equ.equ_map(self.MDSconn)
        else:
            raise Exception("tokamak %s is not supported yet"%self.tokamak)
        
        #correction of the error in the plasma or diagnostic  position
        self.dR_corr = 0
        self.dZ_corr = 0        
        self.m_numbers = np.r_[-4:0, 1:8]
        self.n_numbers = np.r_[-3:0, 1:4]
 
        path = os.path.dirname(os.path.realpath(__file__))
        current_path = os.getcwd()
        os.chdir(path)

        #read a modules for data loading 
        filelist = os.listdir('loaders_'+self.tokamak)  
        filelist = [file[:-3] for file in filelist if file.endswith('.py')] 
        data_rutines = __import__('loaders_'+self.tokamak, fromlist=filelist)

        diag_loaders = {}
        for k, v in data_rutines.__dict__.items():
            if k.startswith('_'): continue
            check, loader = None, None
        
            for k2, v2 in v.__dict__.items():
                if k2.startswith('loader_'): #find loader
                    loader = v2
                if k2.startswith('check'):  #find a shotfile check routine 
                    check = v2
            if loader is not None:
                diag_loaders[k] = loader, check
    
        os.chdir(current_path)  
        
        from collections import OrderedDict
        self.diag_loaders = OrderedDict(sorted(diag_loaders.items(), key=lambda t: t[0]))
   
        self.tables_names = []
        self.curr_tab = 0
        self.movie_saved = 0

        self.create_menu()
        self.create_status_bar()

        self.create_main_frame()
        self.create_spec_table() 
        self.create_phase_table()
        self.create_radialplot_table()

        if self.tokamak != 'NSTX':
            self.create_2Dmap_table()
            self.create_rototomo_table()
        
        if not self.shot is None:
            self.shot_changed(self.shot, init=True)
            
        self.SpecWin.f_range = (fmin, fmax)
        self.SpecWin.t_range = (tmin, tmax)

        self.main_tab.currentChanged.connect(self.updatePanels)
        #try to load data if specified over commadline
        self.change_diagnostics()
        self.change_signal()
        
        #try to load data if specified over commadline

        self.change_diagnostics_phase()
        self.change_signal_phase()

        self.closeEvent = self.__del__
        
            #######################################
        #from time import time
        #t = time()
        #from multiprocessing import  Pool

        #def mds_load(tmp):
            #mds_server,  TDI = tmp
            #MDSconn = mds.Connection(mds_server )
            #data = []
            #for tdi in TDI:
                #try:
                    #data.append(MDSconn.get(tdi).data())
                #except:
                    #data.append([])
            #return data
        #group = 'LFS'
        #names = ['0301', '0302', '0303', '0304', '0305', '0306', '0307', '0308', '0501', '0502', '0503', '0504', '0505', '0506', '0507', '0508', '0701', '0702', '0703', '0704', '0705', '0706', '0707', '0708', '1101', '1102', '1103', '1104', '1105', '1106', '1107', '1108', '1201', '1202', '1203', '1204', '1205', '1206', '1207', '1208', '1301', '1302', '1303', '1304', '1305', '1306', '1307', '1308', '1501', '1502', '1503', '1504', '1505', '1506', '1507', '1508', '1701', '1702', '1703', '1704', '1705', '1706', '1707', '1708', '1901', '1902', '1903', '1904', '1905', '1906', '1907', '1908', '2101', '2102', '2103', '2104', '2105', '2106', '2107', '2108']
        #shot = 180180
        #numTasks = 8
        #TDI = [f'PTDATA2("{group+n}", {shot})' for n in names]
        #TDI = np.array_split(TDI, numTasks)
        #args = [("localhost", tdi) for tdi in TDI]
        
        #out = []
        #pool = Pool()
        #for o in pool.map(mds_load,args):
            #out.extend(o)
        #pool.close()
        #pool.join()
        #print(time()-t)
        
        
    
    def __del__(self, event=None):
        if hasattr(self, 'roto_tomo'):
            del self.roto_tomo
      
 
    def updatePanels(self, panel_ind):
        new_panel = self.tables_names[panel_ind]
        prew_panel = self.tables_names[self.curr_tab]
        #TODO write it in more efficient way 
        self.curr_tab = panel_ind
        QApplication.setOverrideCursor(Qt.WaitCursor)

        method = self.method.currentIndex()
        tau = self.stau.value()/self.SpecWin.sdpi
        
        self.rho = None
        #default values
        selected_trange = self.SpecWin.t_range
        selected_frange = self.SpecWin.f_range

        if prew_panel == 'Spectrogram':
            tau = self.stau.value()/self.SpecWin.sdpi
            method = self.method.currentIndex()
            self.trange = self.SpecWin.ax.get_xlim()
            self.frange = self.SpecWin.ax.get_ylim()
            
        if prew_panel == 'Cross-phaseogram':
            tau = self.stau_phase.value()/self.SpecWin_phase.sdpi
            method = self.method_phase.currentIndex()
            self.trange = self.SpecWin_phase.ax.get_xlim()
            self.frange = self.SpecWin_phase.ax.get_ylim()
            selected_trange = self.SpecWin_phase.t_range
            selected_frange = self.SpecWin_phase.f_range

        if prew_panel == 'Radial Profile':
            selected_trange = self.radial_view.t_range 
            selected_frange = self.radial_view.f_range 
            
        if prew_panel == '2D Te':
            selected_trange = self.Te2Dmap.t_range 
            selected_frange = self.Te2Dmap.f_range 
         
        if prew_panel == '2D SXR' and self.roto_tomo.initialized:
            selected_trange = self.roto_tomo.tmin, self.roto_tomo.tmax
            selected_frange = self.roto_tomo.fmin, self.roto_tomo.fmax
            
        #check if setting are valid
                    
        if new_panel in ['Radial Profile', '2D Te', '2D SXR']\
            and (not all(np.isfinite(selected_trange)) or not all(np.isfinite(selected_frange))):

            QMessageBox.warning(self, "Window not selected", "You must select the area for analysis\n in the spectrogram by a left mouse button ", QMessageBox.Ok)
            QApplication.restoreOverrideCursor()
            return 
        
        

        if  self.data_loader is None and new_panel in ['Radial Profile']:
            QMessageBox.warning(self, "Select diagnostic", "Diagnostic in the spectrogram panel is not selected", QMessageBox.Ok)
            QApplication.restoreOverrideCursor()
            return 
     

  
        #setup the new panel based on the previous one
        
        if new_panel == 'Spectrogram':
            self.SpecWin.ax.set_ylim(self.frange)
            self.SpecWin.ax.set_xlim(self.trange)
            self.SpecWin.t_range = selected_trange
            self.SpecWin.f_range = selected_frange
            
        if new_panel == 'Cross-phaseogram':
            self.SpecWin_phase.ax.set_ylim(self.frange)
            self.SpecWin_phase.ax.set_xlim(self.trange)
            self.SpecWin_phase.t_range = selected_trange
            self.SpecWin_phase.f_range = selected_frange

        if new_panel == 'Radial Profile':
            self.radial_view.t_range = selected_trange
            self.radial_view.f_range = selected_frange
            
            if self.cross_signal.currentIndex() == -1:
                self.cross_diagnostics.setCurrentIndex( self.cb_diagnostics.currentIndex())
                self.cross_sig_group.setCurrentIndex( self.cb_sig_group.currentIndex() )
                self.cross_signal.setCurrentIndex( 0)
            
            #if prew_panel != '2D Te':
            self.change_signal_radial(self.cross_signal.currentIndex(), update=True)

        if new_panel == '2D Te':
    
                
            self.Te2Dmap.t_range = selected_trange
            self.Te2Dmap.f_range = selected_frange
            


                
            self.change_signal_2Dmapping(update=True)
            self.load_ECEI_data()     
            #BUG loaded it only if requested!!
            self.Te2Dmap.prepare()
            self.Te2Dmap.update()
  
            self.save_anim_action.setEnabled(True)
                
        elif new_panel == '2D SXR':
            
            tmin, tmax = selected_trange
            fmin, fmax = selected_frange
            
            if self.SpecWin.initialized:
                sig0 = self.SpecWin.stft_disp.signal
                tvec0 = self.SpecWin.stft_disp.tvec
                if self.roto_tomo.initialized and self.roto_tomo.tmin == tmin and self.roto_tomo.tmax == tmax\
                    and self.roto_tomo.fmin == fmin and self.roto_tomo.fmax == fmax and self.shot == self.roto_tomo.shot:
                        pass  #if nothing has changed, just continue 
                elif self.shot is not None:
  
                    #update also tomography
                    try:
                        self.roto_tomo.prepare_tomo(self.tokamak, self.shot, tmin, tmax, fmin, fmax, self.eqm, tvec0, sig0)
                    except:
                        print('rotation tomography failed')
                        print( traceback.format_exc())

                        QMessageBox.warning(self, "rotation tomography ", traceback.format_exc(), QMessageBox.Ok)
                        return 

                    if self.roto_tomo.showTe:
                        logger.info('prepare ECE')
                        #update also 2D Te 
                        curr_tab = self.curr_tab
                        self.updatePanels(self.tables_names.index('2D Te'))
                        self.curr_tab = curr_tab
            else:
                QMessageBox.warning(self, "Load data issue", "First select frequency and time\n range in the Spectrogram", QMessageBox.Ok)
                return
            self.save_anim_action.setEnabled(True)
        else:
            self.save_anim_action.setEnabled(False)
        QApplication.restoreOverrideCursor()

        
    def save_plot(self, no_gui=False):

        canvas = None
        if not (self.SpecWin.initialized or self.SpecWin_phase.initialized ):
            print('spectrum or cross-spectrum is not initialised')
            return
        
        fine_resolution = False
        if self.curr_tab == 0:
            #increase the resolution before saving
            self.SpecWin.stft_disp.ax_update(self.SpecWin.stft_disp.im_ax, \
                            fine_resolution=fine_resolution, force_update=True)
 
            canvas = self.canvas
            canvas.draw_idle()
            pre = 'spectrum'

        if self.curr_tab == 1:
            #increase the resolution before saving
            self.SpecWin_phase.stft_disp.ax_update(self.SpecWin_phase.stft_disp.im_ax, 
                                             fine_resolution=fine_resolution, force_update=True)
            
            canvas = self.canvas_phase
            canvas.draw_idle()
            pre = 'cross_spectrum'

        if self.curr_tab == 2:
            canvas = self.canvas_radial
            pre = 'radial_profile'

        if self.curr_tab == 3:
            canvas = self.canvas_2Dmap
            pre = '2D_profile'
# git
        if self.curr_tab == 4:
            canvas = self.rototomo_canvas
            pre = '2D_SXR'

        if no_gui:
            path = pre+'_'+str(self.shot)
        else:
            save_gui = QFileDialog(self)
            save_gui.setDefaultSuffix( '.pdf' )
            #save_gui.setConfirmOverwrite(True )

            file_choices =  "PDF (*.pdf);;EPS (*.eps);;SVG (*.svg)"
            selected_filter = u''
            path = save_gui.getSaveFileName(self, 
                            'Save file', pre+'_'+str(self.shot), 
                            file_choices, selected_filter)
            if isinstance(path, tuple): #correction for python 3
                selected_filter = path[1] 
                path = path[0]  

            if path == '': return 
            if path[-4]!= '.': path+= str(selected_filter)[-5:-1]
            
        print( 'path: ', path)

        if not canvas is None:
            canvas.print_figure(path, dpi=self.dpi*4)
            self.statusBar().showMessage('Saved to %s' % path, 5000)
    
    
    def save_anim(self):

        T = time.time()
        print('save_anim',T - self.movie_saved )
        if T - self.movie_saved < 10:
            print('Prevent multiple calls of save_anim')
            return
        
        if self.tables_names[self.curr_tab] == '2D Te':
            t = self.radial_view.t_range[0]
            obj = self.Te2Dmap
            suffix = 'Te'

        elif self.tables_names[self.curr_tab] == '2D SXR':
            t = self.roto_tomo.tmin
            obj = self.roto_tomo
            suffix = 'SXR'
        else:
            raise Exception('animation not supported')

        import matplotlib.animation as animation
        theta = np.linspace(0, 2*np.pi/obj.m, 100//np.abs(obj.m), endpoint=False)
        
        def animate(t):
            sys.stdout.write("\r Writing movie: %2.0f%%" %(t/(2*np.pi/obj.m)*100))
            sys.stdout.flush()
            self.statusBar().showMessage("Writing movie: %2.0f%%" %(t/(2*np.pi/obj.m)*100), 1000)

            obj.shift_phase(t)
            out = obj.update(update_mag=False, update_cax=False, animate=True)
            return out
        
        save_gui = QFileDialog(self)
        save_gui.setDefaultSuffix( '.mp4' )
                
        path = save_gui.getSaveFileName(self, 'Save movie', 
                        'movie_%d_%.4fs_%s.mp4'%(self.shot, t, suffix), "MP4 file (*.mp4)")

        if isinstance(path, tuple): #correction for python 3
            path = path[0]
 
        print('Movie path', path)

        ani = animation.FuncAnimation(obj.fig, animate, theta, interval=25, blit=False, repeat =False) #BUG blit true? 

        writer = animation.FFMpegWriter(fps=15, bitrate=2000)
        try:
            ani.save(path, writer=writer, dpi=self.dpi)
            print( '\n Writing finished in %.1fs'%(time.time()-T))
            del ani
            self.statusBar().showMessage('Saved to %s' % path, 5000)
            self.movie_saved = time.time()
            return
        except:
            print('Animation failed')
            print( traceback.format_exc())

            QMessageBox.warning(self, "Animation failed", traceback.format_exc(), QMessageBox.Ok)
        


    def save_data(self):

        file_choices =  "NPZ (*.npz)"
        
        prew_panel = self.tables_names[self.curr_tab]
        if self.tables_names[self.curr_tab] == 'Cross-phaseogram':
            SpecWin = self.SpecWin_phase
            data_loader = self.data_loader_phase
            diag_group = self.diag_phase
            signal = ''
        else:
            SpecWin = self.SpecWin
            data_loader = self.data_loader
            diag_group = self.diag_group
            signal =  self.signal
            
        try:
            description = data_loader.get_description(diag_group, signal)
        except Exception as e:
            print( 'save_data:', e)
            description = ''
            
        name = description.replace(' ', '_')            
        path = QFileDialog.getSaveFileName(self, 'Save file', name+'.npz', file_choices)
        
        if isinstance(path, tuple): #correction for python 3
            path = path[0]
              
        if path == '': return 

        out = SpecWin.stft_disp.output
        if hasattr(SpecWin.stft_img, 'N_mode'):
            out['mode colorbar'] = SpecWin.stft_img.N_mode
        
        try:
            if self.radial_view.initialized:
                out['radial_prof_rho'] = self.radial_view.rho
                out['radial_prof_Ampl'] = self.radial_view.amplitude
                out['radial_prof_phase'] = self.radial_view.phi
                out['radial_prof_time'] = self.radial_view.t_range
                out['radial_prof_freq'] = self.radial_view.f_range
        except Exception as e:
            print( 'save_data2:', e)
            pass
        
        try:
            if hasattr(self, Te2Dmap) and self.Te2Dmap.initialized:
                out['Te2Dmap_R'] = self.Te2Dmap.Rmag_ece
                out['Te2Dmap_Z'] = self.Te2Dmap.Zmag_ece
                out['Te2DmapTe'] = self.Te2Dmap.shiftTe()
                out['Te2DmapECE_R'] = self.Te2Dmap.R
                out['Te2DmapECE_z'] = self.Te2Dmap.Z
         
        except Exception as e:
            print( 'save_data3:', e)
            pass
        
        if hasattr(SpecWin, 'data_plot'):
            x = SpecWin.data_plot.x
            y = SpecWin.data_plot.y
            xmin, xmax = SpecWin.ax.get_xlim()

            istart, iend  = x.searchsorted((xmin, xmax))
            out['signal'] = y[istart:iend].astype('single', copy=False)
            out['signal_tvec'] = x[istart:iend].astype('single', copy=False)
            
        x, y = SpecWin.plt_trace.get_data()
        out['mode_time'] = x 
        out['mode_freq'] = y

        try:
            info = data_loader.signal_info(diag_group, signal, out['tvec'].mean())
        except Exception as e:
            print( 'save_data:', e)
            info = ''
      
        np.savez_compressed(path, description=description, info=info, **out)
        self.statusBar().showMessage('Saved to %s' % path, 2000)
     
         
    def shot_changed(self, shot, init=False):

        try:
            shot = int(shot)
        except:
            return
        
        if not init and self.shot == shot:
            return 
        
        self.shot = shot

        self.shot_line_phase.setText(str(self.shot))
        self.shot_line.setText(str(self.shot))

        self.eqm_ready = False
        if self.tokamak != 'AUG':
            self.eqm.Close()
            
        #clear data from the previous discharge
        if hasattr(self,  'data_loader_ECEI'):
            del self.data_loader_ECEI
        if hasattr(self,  'data_loader_radial'):
            del self.data_loader_radial
        if hasattr(self,  'data_loader_2Dmap'):
            del self.data_loader_2Dmap
        if hasattr(self,  'data_loader_phase'):
            del self.data_loader_phase
        if hasattr(self,  'data_loader'):
            del self.data_loader
     
        self.SpecWin.reset()
        self.SpecWin_phase.reset()

        #set mode number back to one
        self.radial_m_num.setCurrentIndex(np.where(self.m_numbers == 1)[0])

        self.cross_signal.setCurrentIndex( -1)
        self.cb_diagnostics.clear()
        self.cb_diagnostics.addItem(' ') 
        for n in self.diag_loaders.keys():
            if len(self.diag_loaders[n])==1 or self.diag_loaders[n][1](self.shot):
                self.cb_diagnostics.addItem(str(n)) 
        
        self.cross_diagnostics.clear()
        self.cross_diagnostics.addItem(' ') 
        for n in self.diag_loaders.keys():
            if len(self.diag_loaders[n])==1 or self.diag_loaders[n][1](self.shot):
                self.cross_diagnostics.addItem(str(n))

        self.cb_diagnostics_phase.clear()
        self.cb_diagnostics_phase.addItem(' ') 

        for n, l in self.diag_loaders.items():
            if len(l)==1 or l[1](self.shot):
                if l[0].tor_mode_num or l[0].pol_mode_num:  #load only diagnostics supporting phase measurements
                    self.cb_diagnostics_phase.addItem(str(n)) 

        #clear and prepare the plots
        self.radial_view.reset()
        self.Te2Dmap.reset()


        def init_equlibrium():
            self.eqm_ready = True
            print( 'Initialise equilibrium: ', self.eq_diag, self.shot)

            if self.tokamak == 'AUG':
                import aug_sfutils as sf
                self.eqm = sf.EQU(self.shot, diag=self.eq_diag, exp=self.eq_exp, ed=self.eq_ed)
                if not self.eqm.sf.status:
                    self.eqm = sf.EQU(self.shot, diag='EQI')
                    if not self.eqm.sf.status:
                        print( """standard eq. do not exist!!! use real time equilibrium""")
                        self.eqm_ready = False

            elif self.tokamak == 'DIIID':
                try:
                    self.MDSconn.openTree('D3D', shot)
                    self.BRIEF = self.MDSconn.get(r'\D3D::TOP.COMMENTS:BRIEF').data()
                    self.CONFIG = self.MDSconn.get(r'\D3D::TOP.COMMENTS:CONFIG').data()
                    self.MDSconn.closeTree('D3D', shot)
                    if not isinstance(self.BRIEF,str):
                        self.BRIEF = self.BRIEF.decode()
                    print('Experiment name: ', self.BRIEF)
                except:
                    pass
                
                if not self.eqm.Open(self.shot, diag=self.eq_diag, exp=self.eq_exp, ed=self.eq_ed):
                    print( """Equlibrium shotfile: diag=%s, exp=%s, ed=%d do not exist!!!
                    standard shotfile will be used"""%(self.eq_diag, self.eq_exp, self.eq_ed))
                    if not self.eqm.Open(self.shot, diag='EFIT01'):
                        print( """standard eq. do not exist!!! use real time equilibrium""")
                        if not self.eqm.Open(self.shot, diag='EFITRT01'):
                            self.eqm_ready = False
            elif self.tokamak == 'NSTX':
                   if not self.eqm.Open(self.shot, diag='EFIT02', exp=self.eq_exp):
                       print( """EFIT02 was not found, try EFIT01""")
                       if not self.eqm.Open(self.shot, diag='EFIT01'):
                           eqm_ready = False
                            
            if self.eqm_ready and self.tokamak in ['DIIID','NSTX']:
                self.eqm._read_scalars()
                self.eqm._read_profiles()
                self.eqm._read_pfm()
                self.eqm.read_ssq()
                            

        self.statusBar().showMessage('Loading equilibrium ...', 5000 )

        #initialise equilibrium, save some time later! 
        init_equlibrium()

    def shot_phase_changed(self):

        if not str(self.shot_line_phase.text()).strip().isdigit(): 
            print( 'It is not a valid shot number')
            return  
        self.shot_changed(self.shot_line_phase.text())

    def mode_range_changed(self):

        try:
            num_min = int(self.node_number_min.text())
            num_max = int(self.node_number_max.text())
        except:
            print( 'not a valid number! ', self.node_number_min.text(), self.node_number_max.text())

        self.SpecWin_phase.set_mode_range(num_min, num_max)
    
    def shot_spect_changed(self):

        if not str(self.shot_line.text()).strip().isdigit(): 
            print( 'It is not a valid shot number')
            return 
        self.shot_changed(self.shot_line.text())

        
    def dRdZ_changed(self):
        #change correction of horizontal and vertical position of the magnetic equilibrium

        try:
            dR = self.QLinedR.text()
            dR  = 0 if dR == '' else float(dR)
            dZ = self.QLinedZ.text()
            dZ  = 0 if dZ == '' else float(dZ)
        except Exception as e:
            print('No a valid number')
            return
    
        if dR!= self.dR_corr or dZ != self.dZ_corr:
            self.dR_corr = dR
            self.dZ_corr = dZ
            self.updatePanels(self.curr_tab)
            self.Te2Dmap.dR = dR
            self.Te2Dmap.dZ = dZ


    def change_diagnostics(self, id=None):

        if self.shot is None or (self.diag is None and id is None) or (id is not None and id  <= 0):
            return
    
        QApplication.setOverrideCursor(Qt.WaitCursor)

        if id is not None:
            self.diag = str(self.cb_diagnostics.itemText(id))
        
        if hasattr(self, 'data_loader'):
            #del self.data_loader  #free memory of the previous data_loader
            self.data_loader = None
        try:
            if self.diag not in self.diag_loaders:
                raise Exception('Loader for diagnostic %s was not found'%self.diag)
            self.data_loader = self.diag_loaders[self.diag][0](self.shot, 
                        eqm=self.eqm, rho_lbl=self.rho_lbl, MDSconn=self.MDSconn)
        except Exception as e:
            print('error in loading')
            print( traceback.format_exc())
            QMessageBox.warning(self, "Loading problem", str(e), QMessageBox.Ok)
            return       
        
        QApplication.restoreOverrideCursor()

        self.diag_groups = self.data_loader.get_signal_groups()
        self.cb_sig_group.clear()
        self.cb_signal.clear()

        if len(self.diag_groups) > 1 :
            self.cb_sig_group.addItem(' ') 
        for n in self.diag_groups:
            self.cb_sig_group.addItem(str(n)) 
            
            
    def change_diagnostics_radial(self, id):

        if id <= 0 or self.shot is None: return 
        
        self.diag_radial = str(self.cross_diagnostics.itemText(id))
        QApplication.setOverrideCursor(Qt.WaitCursor)
        if hasattr(self, 'data_loader_radial'):
            #del self.data_loader_radial   #free memory of the previos data_loader
            self.data_loader_radial = None
        try:
            self.data_loader_radial = self.diag_loaders[self.diag_radial][0](self.shot, 
                        ed=0, eqm=self.eqm, rho_lbl=self.rho_lbl, MDSconn=self.MDSconn)
        except Exception as e:
            print( traceback.format_exc())
            QMessageBox.warning(self, "Loading problem", str(e), QMessageBox.Ok)
            return 
            
        self.diag_groups_radial = self.data_loader_radial.get_signal_groups()
        self.cross_sig_group.clear()
        self.cross_signal.clear()
        QApplication.restoreOverrideCursor()

        if len(self.diag_groups_radial) > 1 :
            self.cross_sig_group.addItem(' ') 
        for n in self.diag_groups_radial:
            self.cross_sig_group.addItem(str(n)) 
   

    def change_diagnostics_phase(self, id=None):

        if self.shot is None or (self.diag_phase is None and id is None) or (id is not None and id  <= 0):
            return
             
        if id is not None:
            self.diag_phase =   str(self.cb_diagnostics_phase.itemText(id))
            
        loader = self.diag_loaders[self.diag_phase][0]
        if hasattr(self, 'data_loader_phase'):
            self.data_loader_phase = None
            #del self.data_loader_phase
            
        try:
            self.data_loader_phase = loader(self.shot, 
                        ed=0, eqm=self.eqm, rho_lbl=self.rho_lbl, MDSconn=self.MDSconn)
        except:
            print( traceback.format_exc())
            QMessageBox.warning(self, "Loading problem", str(e), QMessageBox.Ok)
            return 
        self.sig_names_phase = self.data_loader_phase.get_names_phase()
        #print('change diag phase', self.sig_names_phase)
        try:
            self.cb_signal_phase.clear()
        except:
            pass

        if np.size(self.sig_names_phase) > 1:
            self.cb_signal_phase.addItem(' ') 
        for n in self.sig_names_phase:
            self.cb_signal_phase.addItem(str(n)) 
            
            
    def change_sig_group(self, ind):
        
        if ind == 0 and len(self.diag_groups) > 1: return 
        elif len(self.diag_groups) >1: ind-= 1

        if  len(self.diag_groups) == 0:
            logger.error('Diagnostic is not availible')
            return
            
        self.diag_group = self.diag_groups[ind]
        self.sig_names = self.data_loader.get_names(self.diag_group)

        self.cb_signal.clear()
        self.cb_signal.addItem(' ') 
        for n in self.sig_names:
            self.cb_signal.addItem(str(n)) 
            
            
    def change_sig_group_radial(self, ind):
                
        if ind <= 0 and len(self.diag_groups_radial) > 1: return 
    
        elif len(self.diag_groups_radial) > 1: ind-= 1

        self.diag_group_radial = self.diag_groups_radial[ind]
        self.sig_names_radial = self.data_loader_radial.get_names(self.diag_group_radial)

        self.cross_signal.clear()
        self.cross_signal.addItem('  ') 
        for n in self.sig_names_radial:
            self.cross_signal.addItem(str(n)) 

            
    def change_signal(self, num=None):
        
        if num is not None and num <=0  or self.data_loader is None: 
            return 
        
        elif num is not None:
            self.signal = self.sig_names[num-1]

        xlims = self.SpecWin.ax.get_xlim()
        ylims = self.SpecWin.ax.get_ylim()
        
        self.statusBar().showMessage('Loading ...', 1000 )
        
        QApplication.setOverrideCursor(Qt.WaitCursor)
        #http://stackoverflow.com/questions/8218900/how-can-i-change-the-cursor-shape-with-pyqt
        
        T = time.time()
        try:
            tvec, sig = self.data_loader.get_signal(self.diag_group, self.signal)
            if np.size(tvec) < 2:
                raise Exception('Too short time vector len(time) = %d'%np.size(tvec))
            if sig is None: raise

        except Exception as e:
            logger.error( traceback.format_exc())

            QMessageBox.warning(self, "Loading problem", "Check if the signal exist", QMessageBox.Ok)
            return 
 
        try:
            description = self.data_loader.get_description(self.diag_group, self.signal)
        except Exception as e:
            logger.error( traceback.format_exc())
            description = ''

        data = {'tvec':tvec, 'signal':sig}
        
        if time.time()-T > 0.1:
            logger.info(  self.diag+' data loaded in %5.2fs', (time.time()-T))
     
        if self.show_plasma_freq and self.data_loader.radial_profile and self.eqm_ready:
            try:
                rho = self.data_loader.get_rho(self.diag_group, (self.signal, ), 
                                    np.mean(xlims), dR=self.dR_corr, dZ=self.dZ_corr)[0]
                data['freq_tvec'], data['freq'] = self.data_loader.get_plasma_freq(rho)
            except Exception as e:
                #print( traceback.format_exc())
                data['freq_tvec'] = data['freq'] = np.nan
        
        #for NSTX, plot frequency of the modes on q = 1 to q = 4
        if self.show_plasma_freq and hasattr(self.data_loader,'get_plasma_freq_q') and self.eqm_ready:
            try:
                Q = 1,2,3,4
                data['freq_tvec'], data['freq'] = self.data_loader.get_plasma_freq_q(Q)
            except Exception as e:
                print( traceback.format_exc())
                data['freq_tvec'] = data['freq'] = np.nan
             
            
        
        if  self.SpecWin.initialized:
            self.SpecWin.init_plot(data, tmin=xlims[0], tmax=xlims[1], fmin0=ylims[0], 
                            fmax0=ylims[1], description=description)
        else:            
            self.SpecWin.init_plot(data, window=self.win_fun, description=description)
        self.canvas.draw_idle()

        self.InfoThread = GetInfoThread(self, self.data_loader.signal_info, self.diag_group, self.signal, np.mean(xlims))
        self.InfoThread.finished.connect(self.statusBar().showMessage)
        self.InfoThread.start()
        
        QApplication.restoreOverrideCursor()
        

    def change_signal_radial(self, num, update=False):

        if num <= 0 and not update:  return 

        self.statusBar().showMessage('Loading ...' )
        QApplication.setOverrideCursor(Qt.WaitCursor)

        self.signal_radial = self.sig_names_radial[num-1] if num > 0 else None
        tmin, tmax = self.radial_view.t_range
        t0 = (tmin+tmax)/2
        
        if not self.data_loader.radial_profile:
            QMessageBox.warning(self, "Radial profile is not supported", "Radial profile from this diagnostic is not avalible", QMessageBox.Ok)
            return 
        
        #do not load data again, if is not necessary 
        if  not hasattr(self, 'data_loader_radial') or self.data_loader.__class__.__name__ == self.data_loader_radial.__class__.__name__:
            self.data_loader_radial  = self.data_loader
        
        if update:
            T = time.time()

            if self.diag_group is None: return 

            #load data and if diag_group is "all" it will load everything 
            group = self.diag_group
            if self.diag in ['ECE']:
                names = self.data_loader.names  #load all channels 
            else:
                names = self.data_loader.get_names(group)

            data = self.data_loader.get_signal(group, names, calib=True, tmin=tmin, tmax=tmax)
            

            rho, theta, R, Z = self.data_loader.get_rho(group, names, t0, dR=self.dR_corr, dZ=self.dZ_corr)

           
            try:
                Phi = self.data_loader.get_phi_tor()
            except:
                Phi = 0
            

            
            #get a smooth calibration for the stationary profile 
            if self.diag == 'ECE' and hasattr(self.data_loader, 'get_Te0'):
                out = self.data_loader.get_Te0(tmin, tmax, self.dR_corr, self.dZ_corr)
                if not out is None:  #tha case whan CEC and IDA are not avalible 
                    for t, d in zip(out[1], data):
                        m = np.mean(d[1])
                        if m == 0: continue
                        d[1] = d[1]*(t/m)
                        

            mag_coord_data, mag_coord_diag = self.get_mag_eq_coordinates(self.data_loader, t0, rho)


            qsurfs = self.radial_view.q_surf_vals
            rho_qsurfs = self.data_loader.get_q_surfs(t0, qsurfs)
            units = self.data_loader.units
          
            self.radial_view.set_data(self.diag, rho, theta, data, R, Z, Phi, mag_coord_diag, mag_coord_data, rho_qsurfs, units)
            self.radial_view.group = self.diag_group
            

            
            if time.time()-T  > 0.5:
                logger.info('Radial profile data loaded in %5.2fs', (time.time()-T))


        if not self.signal_radial is None:
            cross_sig = self.data_loader_radial.get_signal(self.diag_group_radial, 
                                    self.signal_radial, calib=False, tmin=tmin, tmax=tmax)

            self.radial_view.set_cross_signal(cross_sig)
            info = self.data_loader_radial.signal_info(self.diag_group_radial, 
                                    self.signal_radial, t0)
            InfoThread = GetInfoThread(self, self.data_loader_radial.signal_info, 
                            self.diag_group_radial, self.signal_radial, t0)

        else:
            InfoThread = GetInfoThread(self, self.data_loader.signal_info, self.diag_group, self.signal, t0)

        InfoThread.finished.connect(self.statusBar().showMessage)

        InfoThread.start()
  

        try:
            description = self.data_loader.get_description(self.diag_group, self.signal)
        except Exception as e:
            print( traceback.format_exc())
            description = ''
            
        self.radial_view.description = description
        self.radial_view.update()
        QApplication.restoreOverrideCursor()

    
    def get_mag_eq_coordinates(self, loader, time, rho ):
        rho_grid = np.linspace(0.01, 1, 100)
        magr, magz, theta_geo, theta_star = loader.mag_theta_star(time, rho_grid, dR=self.dR_corr, dZ=self.dZ_corr)
        
        #r = np.hypot(magr-magr[(0, )], magz-magz[(0, )])
        #dr = np.gradient(r, axis=0)/np.gradient(rho_grid)[:, None]
        #dt = np.gradient(theta_star,axis=1)/np.gradient(theta)[None]
        #J = dt*dr  #jacobian of the transformation
        
        
        mag_coord_data = {'R':magr, 'Z':magz, 'theta_geo': theta_geo, 'theta_star': theta_star}
        
        #evalute it at rho locations of the diagnostics measurements
        rho_interp = np.clip(np.abs(rho), *rho_grid[[0,-1]])

        theta_geo    = np.tile(theta_geo, (len(rho_interp), 1))            
        magr = interp1d(rho_grid, magr , axis=0, assume_sorted=True)(rho_interp)
        magz = interp1d(rho_grid, magz , axis=0, assume_sorted=True)(rho_interp)
        theta_star = interp1d(rho_grid, theta_star, axis=0, assume_sorted=True)(rho_interp)
        #diag_J        = interp1d(rho_grid, J    , axis=0, assume_sorted=True)(rho_interp)

        mag_coord_diag = {'R':magr, 'Z':magz, 'theta_geo': theta_geo, 'theta_star': theta_star}
        
        return  mag_coord_data, mag_coord_diag


    def change_signal_2Dmapping(self, update=False):
  
        
        #currently only ECE  and ECEI diags are supported
        diag = 'ECE'

        #if num <= 0 and not update:  return 
    
        #self.signal_radial = self.sig_names_radial[num-1] if num > 0 else None

    
        T = time.time()

        self.statusBar().showMessage('Loading ...' )
        QApplication.setOverrideCursor(Qt.WaitCursor)
        
        tmin, tmax = self.Te2Dmap.t_range
              
        if not np.isfinite(tmax-tmin):
            QMessageBox.warning(self, "Frequency and time window not selected", "You must select the area for analysis\n in the spectrogram by a left mouse button ", QMessageBox.Ok)
            return 
        

        #do not load data again, if is not necessary 
        if hasattr(self, 'data_loader_2Dmap'):
            pass
        elif self.data_loader.__class__.__name__ == 'loader_'+diag:
            self.data_loader_2Dmap = self.data_loader
        else:

            try:
                self.data_loader_2Dmap = self.diag_loaders[diag][0](self.shot, 
                            eqm=self.eqm, rho_lbl=self.rho_lbl, MDSconn=self.MDSconn)
            except Exception as e:
                print('error in loading')
                print( traceback.format_exc())
                QMessageBox.warning(self, "Loading problem", str(e), QMessageBox.Ok)
                return                
        
        group = 'all'
        t0 = (tmin+tmax)/2
        #load signals for cross-correlation, for now use the same as in radial_plot
        if not self.signal_radial is None:
            cross_sig = self.data_loader_radial.get_signal(self.diag_group_radial, 
                                    self.signal_radial, calib=False, tmin=tmin, tmax=tmax)
            
            self.Te2Dmap.set_cross_signal(cross_sig)
            info = self.data_loader_radial.signal_info(self.diag_group_radial, 
                                    self.signal_radial, t0)
            InfoThread = GetInfoThread(self, self.data_loader_radial.signal_info, 
                            self.diag_group_radial, self.signal_radial, t0)

        else:
            #TODO what was this doing?
            InfoThread = GetInfoThread(self, self.data_loader_2Dmap.signal_info, group, self.signal, t0)



        if update:
            #load all ECE channels 
            names = self.data_loader_2Dmap.names 

            units = self.data_loader_2Dmap.units


            data = self.data_loader_2Dmap.get_signal(group, names, calib=True, tmin=tmin, tmax=tmax)
      
            rho, theta, R, Z = self.data_loader_2Dmap.get_rho(group, names, t0, 
                        dR=self.dR_corr, dZ=self.dZ_corr)

            Phi = self.data_loader_2Dmap.get_phi_tor()
         
            
            #get a smooth calibration for the stationary profile 
            if diag == 'ECE' and hasattr(self.data_loader_2Dmap, 'get_Te0'):
                out = self.data_loader_2Dmap.get_Te0(tmin, tmax, self.dR_corr+1e-4, self.dZ_corr)
                #if steady state profiles are availible, crosscalibrate fast data
                if out is not None:
                    for t, d in zip(out[1], data):
                        m = np.mean(d[1])
                        if m == 0: continue
                        d[1] = d[1]*(t/m)

            #get magnetic coordinates for poloidal mapping of the data
            mag_coord_data, mag_coord_diag = self.get_mag_eq_coordinates(self.data_loader_2Dmap, t0, rho)

            
            self.Te2Dmap.set_data(rho,R, Z, theta, Phi,data, mag_coord_diag, mag_coord_data,
                        units, self.Te2D_n_harm,  self.shot, self.use_LFS_data,
                             self.phase_locked_tomo)
          
          
            
            if time.time()-T  > 0.1:
                logger.info('2D profile data loaded in %5.2fs', (time.time()-T))


        


        InfoThread.finished.connect(self.statusBar().showMessage)
        InfoThread.start()

        try:
            description = self.data_loader_2Dmap.get_description(group, 'all')
        except Exception as e:
            print( traceback.format_exc())
            description = ''
        
        #update plot description
        self.Te2Dmap.description = description
        #self.Te2Dmap.update()
        QApplication.restoreOverrideCursor()

    def load_ECEI_data(self):

        if not self.tomo_plot_ecei.isChecked():
            #clean the ECEI overlay
            self.Te2Dmap.update()
            return

        #special case, load data 2D image with ECE
        if self.tokamak != 'DIIID':
            QMessageBox.warning(self, "Not Implemented",
                        'If you need it, follow DIII-D example and implement it by yourself. ', QMessageBox.Ok)
            self.tomo_plot_ecei.setChecked(False)
            return   
        
        
        T = time.time()

        self.statusBar().showMessage('Loading ...' )
        QApplication.setOverrideCursor(Qt.WaitCursor)
        
        tmin, tmax = self.Te2Dmap.t_range

        #do not load data again, if is not necessary 
        if  hasattr(self, 'data_loader_ECEI'):
            pass
        elif  not hasattr(self, 'data_loader_ECEI') and self.data_loader.__class__.__name__ == 'loader_ECEI':
            self.data_loader_ECEI = self.data_loader
        else:
            self.data_loader_ECEI = self.diag_loaders['ECEI'][0](self.shot, 
                            eqm=self.eqm, rho_lbl=self.rho_lbl, MDSconn=self.MDSconn)            
            
        if not (self.data_loader_ECEI.LFSGOOD or self.data_loader_ECEI.HFSGOOD):
            QMessageBox.warning(self, "NO ECEI data",
                        'ECEI data are not availible for this shot', QMessageBox.Ok)
            return 

        #teh HFS system on DIII-D is broken, but can be added in the future
        group = 'LFS'
        names = self.data_loader_ECEI.get_names(group)
        data = self.data_loader_ECEI.get_signal(group, names, calib=True, tmin=tmin, tmax=tmax)
        
        #make sure that chnage in data_loader will not affect data cached in data_loader_2Dmap

        t0 = (tmin+tmax)/2            
        rho, theta, R, Z = self.data_loader_ECEI.get_rho(group, names, t0, 
                    dR=self.dR_corr, dZ=self.dZ_corr)
        Phi = self.data_loader_ECEI.get_phi_tor()
 
        self.Te2Dmap.set_data2D(R, Z, Phi,  data)
        self.Te2Dmap.update()
        
        
        logger.info('ECEI profile data loaded in %5.2fs', (time.time()-T))

        
        
        QApplication.restoreOverrideCursor()
        
        
        
    def change_signal_phase(self, num=None):

        if num is not None and (num < 0 or num == 0 and len(self.sig_names_phase) > 1) or self.data_loader_phase is None: 
            return 
        elif num is not None:
            self.signal_phase = self.sig_names_phase[num-1]

        self.statusBar().showMessage('Loading ...' )
        QApplication.setOverrideCursor(Qt.WaitCursor)

        xlims = self.SpecWin_phase.ax.get_xlim()
        ylims = self.SpecWin_phase.ax.get_ylim()
        T = time.time()
        try:
            tvec, sig = self.data_loader_phase.get_signal_phase(self.signal_phase)
            if sig is None: raise
        except Exception as e:
            print( traceback.format_exc())

            QMessageBox.warning(self, "Loading problem", "Check if the signal exist", QMessageBox.Ok)
            return 
        
        if time.time()-T > 0.1:
            logger.info(  self.diag_phase+' data loaded in %5.2fs', (time.time()-T))
            
        self.statusBar().showMessage('')
        #option to set mode range from config file 
        if self.mode_range is None:
            mode_range = self.balooning_coils_mode_range
        else:
            mode_range = self.data_loader_phase.mode_range

        if not 'Toroidal' in self.signal_phase and hasattr( self.data_loader_phase, 'get_theta_pol'):
            phi = self.data_loader_phase.get_theta_pol( self.signal_phase, np.mean(xlims), 
                                                       rhop=self.rho_pol_mode)
            lbl = 'M'
        else:
            phi = self.data_loader_phase.get_phi_tor( self.signal_phase)
            lbl = 'N'
       
        phase_corr = self.data_loader_phase.get_phase_corrections( self.signal_phase)
        data = {'tvec':tvec, 'signal':sig, 'Phi':phi, 'phase_correction':phase_corr}
  

        #for NSTX, plot frequency of the modes on q = 1 to q = 4
        if self.show_plasma_freq and hasattr(self.data_loader,'get_plasma_freq_q') and self.eqm_ready:
            try:
                Q = 1,2,3,4
                data['freq_tvec'], data['freq'] = self.data_loader.get_plasma_freq_q(Q)
            except Exception as e:
                print( traceback.format_exc())
                data['freq_tvec'] = data['freq'] = np.nan
            
        
        try:
            description = self.data_loader_phase.get_description(   self.signal_phase , '')
        except:
            print( traceback.format_exc())

        if  self.SpecWin_phase.initialized:
            self.SpecWin_phase.init_plot(data, tmin=xlims[0], tmax=xlims[1], 
                        fmin0=ylims[0], fmax0=ylims[1], mode_range=mode_range, mode_num_lbl=lbl)
        else:
            self.SpecWin_phase.init_plot(data, mode_range=mode_range, mode_num_lbl=lbl, description=description)
        self.canvas_phase.draw()
        QApplication.restoreOverrideCursor()

        InfoThread = GetInfoThread(self, self.data_loader_phase.signal_info, 
                        self.diag_group, self.signal_phase, (xlims[1]+xlims[0])/2)
        InfoThread.finished.connect(self.statusBar().showMessage)
        InfoThread.start()

        
    def on_about(self):
        #msg = """ This is a tool for the preparation of the input signals for the tomography. 
       
        #Basic data processing 
        #*** Selecting proper time interval
        #*** Averadging by box car method and downsampling
        #*** Removing of the corupted channels
        
        #Advaced data processing 
        #*** Tool for design of the FIR filter
        #*** Compressed coding by complex SVD
        #*** Tool for verification of the results
        #"""
        msg = 'Simple universal tool to visualize a spectral information, \n in case of troubles contact todstrci@ipp.mpg.de'
        
        QMessageBox.about(self, "About this tool", msg.strip())
        
        
    def on_manual(self):

        path = os.path.dirname(os.path.realpath(__file__))
        text = open(path+os.sep+'README.md').read()
        
        help_win = QDialog()

        self.textEdit = QTextEdit(help_win)
        self.textEdit.setGeometry(QRect(0, 0, 750, 500))
        self.textEdit.setObjectName("textEdit")
        self.textEdit.setReadOnly(True)
        
        self.textEdit.setText(text)
        help_win.exec_()

        return

    
    def on_pick(self, event):
        # The event received here is of the type
        # matplotlib.backend_bases.PickEvent
        #
        # It carries lots of information, of which we're using
        # only a small amount here.
        # 
        box_points = event.artist.get_bbox().get_points()
        msg = "You've clicked on a bar with coords:\n %s" % box_points
        
        QMessageBox.information(self, "Click!", msg)
    
    
    def create_main_frame(self):

        self.dpi = 100
        self.cWidget = QWidget(self)
        self.gridLayout = QGridLayout(self.cWidget)
        self.setWindowTitle('Spectra viewer')
        self.Expand = QSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        self.Fixed  = QSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed)
        self.main_tab = QTabWidget(self)
        self.main_tab.setTabPosition(QTabWidget.North)
        self.gridLayout.addWidget(self.main_tab, 0, 0, 1, 2)
        self.setCentralWidget(self.cWidget)

        self.resize(600, 900)
        self.setCenter()
        
        
    def create_spec_table(self):
        
        self.tab_widget = QWidget(self)
        self.tables_names.append('Spectrogram')
        self.main_tab.insertTab(len(self.tables_names)-1, self.tab_widget, self.tables_names[-1])
        self.tab_widget.setSizePolicy(self.Expand)
        self.verticalLayout = QVBoxLayout(self.tab_widget)
        self.horizontalLayout = QHBoxLayout()

        self.verticalLayout.addLayout(self.horizontalLayout)

        self.setCenter()

        self.fig_spec = Figure(dpi=self.dpi)
        self.canvas = FigureCanvas(self.fig_spec)
        self.canvas.setParent(self.cWidget)

        self.canvas.setSizePolicy(self.Expand)
        self.canvas.updateGeometry()
    
        c =  self.cWidget.palette().color(QPalette.Base)
        
        self.fig_spec.patch.set_facecolor((c.red()/255., c.green()/255., c.blue()/255.))
        
        self.stau = QSlider(Qt.Horizontal)
        self.sgamma = QSlider(Qt.Horizontal)
        self.method  = QComboBox(self.cWidget)
        
        self.SpecWin =  SpectraViewer( self.sgamma, self.stau, self.method, 
                    self.statusBar().showMessage, show_raw=self.show_raw, allow_selector=True, 
                        fig= self.fig_spec, cmap=self.spect_cmap)

        # Bind the 'pick' event for clicking on one of the bars

        self.canvas.mpl_connect('pick_event', self.on_pick)

        # Create the navigation toolbar, tied to the canvas

        self.mpl_toolbar = NavigationToolbar(self.canvas, self.cWidget)
        
        self.shot_line = QLineEdit(self.cWidget)
        if not self.shot is None: self.shot_line.setText(str(self.shot))
        self.shot_line.setToolTip('Write shot number')
        self.shot_line.setFixedWidth(65)

        self.cb_diagnostics  = QComboBox(self.cWidget)
        self.cb_diagnostics.setToolTip('Select diagnotics')
        self.cb_diagnostics.setFixedWidth(100)

        self.cb_sig_group  = QComboBox(self.cWidget)
        self.cb_sig_group.setToolTip('Select signal group')
        self.cb_sig_group.setFixedWidth(100)

        self.cb_signal  = QComboBox(self.cWidget)
        self.cb_signal.setToolTip('Select signal for the spectrogram')
        self.cb_signal.setFixedWidth(80)

        label1 = QLabel('Shot:')
        label2 = QLabel('Diag:')
        label3 = QLabel('Group:')
        label4 = QLabel('Sig:')

        self.shot_line.editingFinished.connect(self.shot_spect_changed)
        self.cb_diagnostics.currentIndexChanged.connect(self.change_diagnostics)
        self.cb_sig_group.currentIndexChanged.connect(self.change_sig_group)
        self.cb_signal.currentIndexChanged.connect(self.change_signal)

        self.method.setFixedWidth(70)

        self.method.currentIndexChanged.connect(self.SpecWin.select_DFT_backend)

        self.method.setToolTip('Choose algorithm')
        self.stau.setToolTip('Set width of window for FFT')
        self.sgamma.setToolTip('Set contrast')

        slider_label2 = QLabel('TAU:')
        slider_label1 = QLabel('GAMMA:')
        combo_label3 = QLabel('Method:')

        self.sgamma.valueChanged.connect(self.SpecWin.apply_slider)
        self.stau.valueChanged.connect(self.SpecWin.apply_slider)
        
        #
        # Layout with box sizers
        # 

        hboxUp = QHBoxLayout()
        hboxDown = QHBoxLayout()        

        hboxUpL = QHBoxLayout()
        hboxUpR = QHBoxLayout()

        hboxUpL.addWidget(self.mpl_toolbar)
        hboxUpL.setAlignment(self.mpl_toolbar, Qt.AlignLeft)

        for l, w in ((slider_label1, self.sgamma), (slider_label2, self.stau)):
            hbox = QHBoxLayout()
            hbox.addWidget(l)
            hbox.addWidget(w)
            hbox.setAlignment(l, Qt.AlignRight)
            hbox.setAlignment(w, Qt.AlignLeft)
            hboxUpR.addLayout(hbox)
            
        hboxUp.addLayout(hboxUpL)
        hboxUp.addLayout(hboxUpR)

        for l, w in [(label1, self.shot_line), 
                     (label2, self.cb_diagnostics), 
                     (label3, self.cb_sig_group), 
                     (label4, self.cb_signal), 
                     (combo_label3, self.method)]:

            hbox = QHBoxLayout()
            hbox.addWidget(l)
            hbox.addWidget(w)
            hbox.setAlignment(l, Qt.AlignRight)
            hbox.setAlignment(w, Qt.AlignLeft)
            hboxDown.addLayout(hbox)
        
        vbox = QVBoxLayout()
        vbox.addWidget(self.canvas)

        vbox.addLayout(hboxUp)
        vbox.addLayout(hboxDown)

        self.horizontalLayout.addLayout(vbox)


    def create_radialplot_table(self):
        
        self.tab_widget_radial = QWidget(self)
  
        self.tables_names.append('Radial Profile')
        self.main_tab.insertTab(len(self.tables_names)-1, self.tab_widget_radial, 
                                self.tables_names[-1])
      
        self.tab_widget_radial.setSizePolicy(self.Expand)
        self.verticalLayout_radial = QVBoxLayout(self.tab_widget_radial)
        self.horizontalLayout_radial = QHBoxLayout()

        self.verticalLayout_radial.addLayout(self.horizontalLayout_radial )

        self.dpi = 100
        self.fig_radial = Figure(dpi=self.dpi)
        self.canvas_radial = FigureCanvas(self.fig_radial)
        self.canvas_radial.setParent(self.cWidget)
        self.canvas_radial.setSizePolicy(self.Expand)
        self.canvas_radial.updateGeometry()

        c =  self.cWidget.palette().color(QPalette.Base)

        self.fig_radial.patch.set_facecolor((c.red()/255., c.green()/255., c.blue()/255.))

        self.radial_view = RadialViewer(self, self.fig_radial, self.Te2D_n_harm, self.rho_lbl)

        self.canvas_radial.mpl_connect('pick_event', self.on_pick)
        
        self.mpl_toolbar_radial = NavigationToolbar(self.canvas_radial, self.cWidget)
        
        self.cross_diagnostics  = QComboBox(self.cWidget)
        self.cross_diagnostics.setToolTip('Select diagnotics for the crosscorelation')
        self.cross_diagnostics.setFixedWidth(80)

        self.cross_sig_group  = QComboBox(self.cWidget)
        self.cross_sig_group.setToolTip('Select signal group for the crosscorelation')
        self.cross_sig_group.setFixedWidth(80)

        self.cross_signal  = QComboBox(self.cWidget)
        self.cross_signal.setToolTip('Select signal for the crosscorelation')
        self.cross_signal.setFixedWidth(80)

        label_mode = QLabel('M mode:')
        self.radial_m_num  = QComboBox(self.cWidget)
        self.radial_m_num.setToolTip('You can specify poloidal number of the mode')
        self.radial_m_num.setFixedWidth(50)

        for n in self.m_numbers:  self.radial_m_num.addItem(str(n)) 
        self.radial_m_num.setCurrentIndex(np.where(self.m_numbers == 1)[0])

        label2 = QLabel('Diag:')
        label3 = QLabel('Group:')
        label4 = QLabel('Sig:')

        self.cross_diagnostics.currentIndexChanged.connect(self.change_diagnostics_radial)
        self.cross_sig_group.currentIndexChanged.connect(self.change_sig_group_radial)
        self.cross_signal.currentIndexChanged.connect(self.change_signal_radial)
        self.radial_m_num.currentIndexChanged.connect(self.radial_view.change_m_number)
    
        labeldR = QLabel('dR:')
        labeldZ = QLabel('dZ:')

        self.QLinedR = QLineEdit(self.cWidget)
        self.QLinedZ = QLineEdit(self.cWidget)
        
        self.load_config()
        self.QLinedR.setText(str(self.dR_corr))
        self.QLinedZ.setText(str(self.dZ_corr))

        self.QLinedR.setToolTip('horizontal shift of the magnetic equilibrium')
        self.QLinedR.setFixedWidth(40)
        
        self.QLinedZ.setToolTip('vertical shift of the magnetic equilibrium')
        self.QLinedZ.setFixedWidth(40)
        
        self.QLinedR.editingFinished.connect(self.dRdZ_changed)
        self.QLinedZ.editingFinished.connect(self.dRdZ_changed)

        #
        # Layout with box sizers
        # 
        hboxUp = QHBoxLayout()
        hboxDown = QHBoxLayout()        

        hboxUp.addWidget(self.mpl_toolbar_radial)
        hboxUp.setAlignment(self.mpl_toolbar_radial, Qt.AlignLeft)    

        for l, w in [ (label2, self.cross_diagnostics), (label3, 
                  self.cross_sig_group), (label4, self.cross_signal), (label_mode, 
                  self.radial_m_num), (labeldR, self.QLinedR), (labeldZ, self.QLinedZ)]:
            hbox = QHBoxLayout()
            hbox.addWidget(l)
            hbox.addWidget(w)
            hbox.setAlignment(l, Qt.AlignRight)
            hbox.setAlignment(w, Qt.AlignLeft)
            hboxDown.addLayout(hbox)
        
 
        vbox = QVBoxLayout()
        vbox.addWidget(self.canvas_radial)

        vbox.addLayout(hboxUp)
        vbox.addLayout(hboxDown)
        self.horizontalLayout_radial.addLayout(vbox)
        
       
    def create_phase_table(self):

        self.tab_widget_phase = QWidget(self)
        
        self.tables_names.append('Cross-phaseogram')
                                
        self.main_tab.insertTab(len(self.tables_names)-1, self.tab_widget_phase, 
                                    self.tables_names[-1])      

        self.tab_widget_phase.setSizePolicy(self.Expand)
        self.verticalLayout_phase = QVBoxLayout(self.tab_widget_phase)
        self.horizontalLayout_phase = QHBoxLayout()

        self.verticalLayout_phase.addLayout(self.horizontalLayout_phase )

        self.dpi = 100
        self.fig_phase = Figure(dpi=self.dpi)
        self.canvas_phase = FigureCanvas(self.fig_phase)
        self.canvas_phase.setParent(self.cWidget)
        self.canvas_phase.setSizePolicy(self.Expand)
        self.canvas_phase.updateGeometry()

        c =  self.cWidget.palette().color(QPalette.Base)

        self.fig_phase.patch.set_facecolor((c.red()/255., c.green()/255., c.blue()/255.))
        
        # Create the mpl Figure and FigCanvas objects. 
        #        
        # Since we have only one plot, we can use add_axes 
        # instead of add_subplot, but then the subplot
        # configuration tool in the navigation toolbar wouldn't
        # work.
        
        self.stau_phase = QSlider(Qt.Horizontal)
        self.sgamma_phase = QSlider(Qt.Horizontal)
        self.method_phase  = QComboBox(self.cWidget)
        self.SpecWin_phase =  SpectraViewer( self.sgamma_phase, self.stau_phase 
                        , self.method_phase, self.statusBar().showMessage, show_raw=False, cmap=self.spect_cmap, 
                        allow_selector=True, fig= self.fig_phase, phase_analysis=True, show_colorbar=True)

        # Bind the 'pick' event for clicking on one of the bars

        self.canvas_phase.mpl_connect('pick_event', self.on_pick)

        # Create the navigation toolbar, tied to the canvas

        self.mpl_toolbar_phase = NavigationToolbar(self.canvas_phase, self.cWidget)
        
        self.shot_line_phase = QLineEdit(self.cWidget)
        if not self.shot is None: self.shot_line_phase.setText(str(self.shot))
        self.shot_line_phase.setToolTip('Write AUG shot number')
        self.shot_line_phase.setFixedWidth(65)

        self.cb_diagnostics_phase = QComboBox(self.cWidget)
        self.cb_diagnostics_phase.setToolTip('Select diagnotics')
        self.cb_diagnostics_phase.setFixedWidth(80)

        self.cb_signal_phase  = QComboBox(self.cWidget)
        self.cb_signal_phase.setToolTip('Select signals for the phaseogram')
        self.cb_signal_phase.setFixedWidth(80)
        
        self.node_number_min = QLineEdit(self.cWidget)
        self.node_number_min.setText(str(self.mode_range[0]))
        self.node_number_min.setToolTip('Lower limit for the investigated range of the mode numbers')
        self.node_number_min.setFixedWidth(35)

        self.node_number_max = QLineEdit(self.cWidget)
        self.node_number_max.setText(str(self.mode_range[1]))
        self.node_number_max.setToolTip('U limit for the investigated range of the mode numbers')
        self.node_number_max.setFixedWidth(35)
     
        label1 = QLabel('Shot:')
        label2 = QLabel('Diag:')
        label4 = QLabel('Sig:')
        label_min = QLabel('  Mode num. from')
        label_max = QLabel('to ')
        
        self.shot_line_phase.editingFinished.connect(self.shot_phase_changed)
        self.node_number_min.editingFinished.connect(self.mode_range_changed)
        self.node_number_max.editingFinished.connect(self.mode_range_changed)
        self.cb_diagnostics_phase.currentIndexChanged.connect(self.change_diagnostics_phase)
        self.cb_signal_phase.currentIndexChanged.connect(self.change_signal_phase)

        self.method_phase.setFixedWidth(70)
        self.method_phase.currentIndexChanged.connect(self.SpecWin_phase.select_DFT_backend)
        self.method_phase.setToolTip('Choose algorithm')
        self.stau_phase.setToolTip('Set width of window for FFT')
        self.sgamma_phase.setToolTip('Set contrast')

        slider_label2 = QLabel('TAU:')
        slider_label1 = QLabel('GAMMA:')
        combo_label3  = QLabel('Method:')

        self.sgamma_phase.valueChanged.connect(self.SpecWin_phase.apply_slider)
        self.stau_phase.valueChanged.connect(self.SpecWin_phase.apply_slider)

        # Layout with box sizers
 
        hboxUp = QHBoxLayout()
        hboxDown = QHBoxLayout()        

        hboxUpL = QHBoxLayout()
        hboxUpL.addWidget(self.mpl_toolbar_phase)
        hboxUpL.setAlignment(self.mpl_toolbar_phase, Qt.AlignLeft)
        
        hboxUpR = QHBoxLayout()

        for l, w in ((slider_label1, self.sgamma_phase), (slider_label2, self.stau_phase)):
            hbox = QHBoxLayout()
            hbox.addWidget(l)
            hbox.addWidget(w)
            hbox.setAlignment(l, Qt.AlignRight)
            hbox.setAlignment(w, Qt.AlignLeft)
            hboxUpR.addLayout(hbox)

        hboxUp.addLayout(hboxUpL)
        hboxUp.addLayout(hboxUpR)

        for L, W, k in [(label1, self.shot_line_phase, True), 
                     (label2, self.cb_diagnostics_phase, True), 
                     (label4, self.cb_signal_phase, True), 
                     ((label_min, label_max), (self.node_number_min, self.node_number_max), False), 
                     (combo_label3, self.method_phase, True)]:
            
            hbox = QHBoxLayout()
            if k: L, W = [L], [W]
            for l, w in zip(L, W):
                hbox_inner = QHBoxLayout()
                hbox_inner.addWidget(l)
                hbox_inner.addWidget(w)
                hbox_inner.setAlignment(l, Qt.AlignRight)
                hbox_inner.setAlignment(w, Qt.AlignLeft)
                hbox.addLayout(hbox_inner)

            hboxDown.addLayout(hbox)
        
        vbox = QVBoxLayout()
        vbox.addWidget(self.canvas_phase)

        vbox.addLayout(hboxUp)
        vbox.addLayout(hboxDown)

        self.horizontalLayout_phase.addLayout(vbox)

        
    def create_2Dmap_table(self):
        
        self.tab_widget_tomo = QWidget(self)
        self.tables_names.append('2D Te')
        self.main_tab.insertTab(len(self.tables_names)-1, self.tab_widget_tomo, self.tables_names[-1])
        self.tab_widget_tomo.setSizePolicy(self.Expand)
        self.verticalLayout_tomo = QVBoxLayout(self.tab_widget_tomo)
        self.horizontalLayout_tomo  = QHBoxLayout()

        self.verticalLayout_tomo.addLayout(self.horizontalLayout_tomo )
        self.setCenter()

        self.fig_tomo = Figure(dpi=self.dpi)
        self.canvas_2Dmap = FigureCanvas(self.fig_tomo)
        self.canvas_2Dmap.setParent(self.cWidget)
        self.canvas_2Dmap.setFocusPolicy( Qt.WheelFocus )
        self.canvas_2Dmap.setFocus()
        self.canvas_2Dmap.setSizePolicy(self.Expand)
        self.canvas_2Dmap.updateGeometry()

        c =  self.cWidget.palette().color(QPalette.Base)
        
        self.fig_tomo.patch.set_facecolor((c.red()/255., c.green()/255., c.blue()/255.))
        self.tab_widget_tomo.setToolTip('Use a mouse wheel to rotate the mode, \n Ctrl+wheel for a fine shift')
        
        self.tomo_limit = QSlider(Qt.Horizontal)
        self.tomo_limit.setRange(0, 99)
        self.tomo_limit.setValue(70)
        self.tomo_limit.setTracking(True)
        self.tomo_limit.setTickPosition(QSlider.NoTicks)
        self.tomo_limit.setSingleStep(.05)

        self.canvas_2Dmap.mpl_connect('pick_event', self.on_pick)

        # Create the navigation toolbar, tied to the canvas

        self.mpl_toolbar_2D = NavigationToolbar(self.canvas_2Dmap, self.cWidget)
        self.mpl_toolbar_2D.setMaximumWidth(150)

        label2 = QLabel('Show ECEI')
        self.tomo_plot_ecei= QCheckBox(self.cWidget)
        self.tomo_plot_ecei.setToolTip('Add overlay with ECEI data')
        

        label1 = QLabel('M mode:')
        self.tomo_m_num  = QComboBox(self.cWidget)
        self.tomo_m_num.setToolTip('You must specify poloidal number of the mode')
        self.tomo_m_num.setFixedWidth(50)

        label3 = QLabel('Remove background')
        self.tomo_plot_substract = QCheckBox(self.cWidget)
        self.tomo_plot_substract.setToolTip('Substract the time averaged signal')
        
        for n in self.m_numbers:  self.tomo_m_num.addItem(str(n)) 
        self.tomo_m_num.setCurrentIndex( np.where(self.m_numbers == 1)[0])

        self.Te2Dmap = Diag2DMapping(self, self.fig_tomo, self.n_contour_2Dprof,
                                     self.tomo_plot_substract, self.tomo_plot_ecei)

        self.tomo_m_num.currentIndexChanged.connect( self.Te2Dmap.UpdateModeM)
        self.tomo_plot_substract.stateChanged.connect(self.Te2Dmap.UpdatePlotType)
        self.tomo_plot_ecei.stateChanged.connect(self.load_ECEI_data)
   
        label4 = QLabel('Scale')
        self.tomo_limit.valueChanged.connect(self.Te2Dmap.UpdateLim)
        self.tomo_limit.setToolTip('Set a lower limit for the Te plot')
        self.tomo_limit.setFixedWidth(50)

        # Layout with box sizers

        hboxC = QHBoxLayout()
        hboxL = QHBoxLayout()
        hboxR = QHBoxLayout()
        
        hboxL.addWidget(self.mpl_toolbar_2D)
        hboxL.setAlignment( self.mpl_toolbar_2D, Qt.AlignLeft)

        for l, w in [(label2, self.tomo_plot_ecei), (label1, self.tomo_m_num), (label3, 
                  self.tomo_plot_substract), (label4, self.tomo_limit) ]:
            hbox = QHBoxLayout()
            hbox.addWidget(l)
            hbox.setAlignment(l, Qt.AlignRight)
            hbox.addWidget(w)
            hbox.setAlignment(w, Qt.AlignLeft)
            hboxR.addLayout(hbox)
            
        hboxC.addLayout(hboxL)
        hboxC.addLayout(hboxR)
        
        vbox = QVBoxLayout()
        vbox.addWidget(self.canvas_2Dmap)
        vbox.addLayout(hboxC)
        self.horizontalLayout_tomo.addLayout(vbox)

       
    def create_rototomo_table(self):
        
        try:
            from roto_tomo.roto_tomo import Roto_tomo
        except Exception as e:
            print('import of the rotation tomography failed: '+str(e))
            print( traceback.format_exc())

            return 
        
        self.tab_widget_rototomo = QWidget(self)
        self.tables_names.append('2D SXR')
        self.main_tab.insertTab(len(self.tables_names)-1, self.tab_widget_rototomo, self.tables_names[-1])
        self.tab_widget_rototomo.setSizePolicy(self.Expand)
        self.verticalLayout_rototomo = QVBoxLayout(self.tab_widget_rototomo)
        self.horizontalLayout_rototomo  = QHBoxLayout()

        self.verticalLayout_rototomo.addLayout(self.horizontalLayout_rototomo )
        
        self.setCenter()

        self.fig_rototomo = Figure(dpi=self.dpi)
        self.rototomo_canvas = FigureCanvas(self.fig_rototomo)
        
        self.rototomo_canvas.setParent(self.cWidget)
        self.rototomo_canvas.setFocusPolicy( Qt.WheelFocus )
        self.rototomo_canvas.setFocus()        

        self.rototomo_canvas.setSizePolicy(self.Expand)
        self.rototomo_canvas.updateGeometry()

        c =  self.cWidget.palette().color(QPalette.Base)
        
        self.fig_rototomo.patch.set_facecolor((c.red()/255., c.green()/255., c.blue()/255.))
        self.rototomo_canvas.setToolTip('Use a mouse wheel to rotate the mode, \n Ctrl+wheel for a fine shift')
        
        self.rototomo_limit = QSlider(Qt.Horizontal)
        self.rototomo_limit.setRange(1, 100)
        self.rototomo_limit.setValue(100)
        self.rototomo_limit.setTracking(True)
        self.rototomo_limit.setTickPosition(QSlider.NoTicks)
        self.rototomo_limit.setSingleStep(.05)
           
        self.rototomo_reg = QSlider(Qt.Horizontal)
        self.rototomo_reg.setRange(1, 100)
        self.rototomo_reg.setValue(70)
        self.rototomo_reg.setTracking(True)
        self.rototomo_reg.setTickPosition(QSlider.NoTicks)
        self.rototomo_reg.setSingleStep(.01)

        self.rototomo_canvas.mpl_connect('pick_event', self.on_pick)
        
        # Create the navigation toolbar, tied to the canvas

        self.mpl_toolbar_rototomo = NavigationToolbar(self.rototomo_canvas, self.cWidget)
        items = 'Home', 'Pan', 'Zoom'
        self.mpl_toolbar_rototomo.toolitems = [t for t in NavigationToolbar2QT.toolitems if t[0] in items]
        self.mpl_toolbar_rototomo.setMaximumWidth(150)

        grid = QGridLayout()

        OptionsFrame = QGroupBox("Options" )
        OptionsFrame.setContentsMargins(0, 10, 0, 0)
        
        sizePolicy = QSizePolicy( QSizePolicy.Fixed, QSizePolicy.Expanding)
        OptionsFrame.setSizePolicy(sizePolicy)
        ModeFrame = QGroupBox("MHD Mode:" )
        ModeFrame.setContentsMargins(0, 20, 0, 0)
        label_m = QLabel('M:')
        self.rototomo_m_num  = QComboBox(self.cWidget)
        self.rototomo_m_num.setToolTip('You must specify poloidal number of the mode')
        self.rototomo_m_num.setFixedWidth(50)
        
        label_n = QLabel('N:')
        self.rototomo_n_num  = QComboBox(self.cWidget)
        self.rototomo_n_num.setToolTip('You must specify toroidal number of the mode')
        self.rototomo_n_num.setFixedWidth(50)
        
        label_bcg = QLabel('Background')
        self.rototomo_plot_bcg = QCheckBox(self.cWidget)
        self.rototomo_plot_bcg.setToolTip('Add the time averaged background')
        self.rototomo_plot_bcg.setChecked(True)
        
        label_TeOver = QLabel('T<sub>e</sub> contours')
        self.rototomo_TeOver = QCheckBox(self.cWidget)
        self.rototomo_TeOver.setToolTip('Overlayed by electron temperature contours')
        
        label_mag_flx = QLabel('Flux surfaces')
        self.rototomo_mag_flx = QCheckBox(self.cWidget)
        self.rototomo_mag_flx.setToolTip('Overlayed by magnetic flux surfaces')
        self.rototomo_mag_flx.setChecked(True)


        for m in self.m_numbers:  self.rototomo_m_num.addItem(str(m)) 
        self.rototomo_m_num.setCurrentIndex( np.where(self.m_numbers == 1)[0])
        for n in self.n_numbers:  self.rototomo_n_num.addItem(str(n)) 
        self.rototomo_n_num.setCurrentIndex( np.where(self.n_numbers == -1)[0])
   
        label_lim = QLabel('Scale')
        self.rototomo_limit.setToolTip('Set a lower limit for the colorscale')
        self.rototomo_limit.setFixedWidth(120)

        label_reg = QLabel('Regularization')
        self.rototomo_reg.setToolTip('Set a regularization level')
        self.rototomo_reg.setFixedWidth(120)
        self.rototomo_button = QPushButton('Set data', self)
        self.rototomo_button.clicked.connect(self.handleButtonData)
        self.rototomo_button.setFixedWidth(130)
        self.retrofit_button = QPushButton('Retrofit', self)
        self.retrofit_button.setFixedWidth(130)

        # Layout with box sizers
 
        vbox1 = QVBoxLayout()
        vbox1.setAlignment( Qt.AlignTop)

        vbox1.addWidget(ModeFrame)
        vbox1_up =  QVBoxLayout()
        ModeFrame.setLayout(vbox1_up)

        for l, w, frame in [(label_m, self.rototomo_m_num, vbox1_up), 
                    (label_n,    self.rototomo_n_num, vbox1_up), 
                    (label_bcg, self.rototomo_plot_bcg, vbox1), 
                    (label_TeOver, self.rototomo_TeOver, vbox1), 
                    (label_mag_flx, self.rototomo_mag_flx, vbox1)]:
            subbbox = QHBoxLayout()
            subbbox.addWidget(l)
            subbbox.setAlignment(l, Qt.AlignLeft)
            subbbox.addWidget(w)
            subbbox.setAlignment(w, Qt.AlignRight)
            frame.addLayout(subbbox)

        for l, w in [(label_lim,  self.rototomo_limit), 
                    (label_reg,  self.rototomo_reg)]:
            subvbox = QVBoxLayout()
            subvbox.addWidget(l)
            subvbox.setAlignment(l, Qt.AlignBottom)
            subvbox.addWidget(w)
            subvbox.setAlignment(w, Qt.AlignTop)
            vbox1.addLayout(subvbox)

        OptionsFrame.setLayout(vbox1)
        outer_vbox = QVBoxLayout()
        outer_vbox.addWidget(OptionsFrame )
        outer_vbox.addWidget(self.rototomo_button)
        outer_vbox.addWidget(self.retrofit_button)
        outer_vbox.addWidget(self.mpl_toolbar_rototomo)
        outer_vbox.setAlignment( self.mpl_toolbar_rototomo, Qt.AlignLeft)

        self.horizontalLayout_rototomo.addLayout(outer_vbox )
        self.horizontalLayout_rototomo.addWidget(self.rototomo_canvas)

        self.roto_tomo =  Roto_tomo(self, self.fig_rototomo, self.rototomo_m_num, 
                        self.rototomo_n_num, self.rototomo_plot_bcg, 
                        self.rototomo_TeOver, self.rototomo_mag_flx, self.rototomo_limit, 
                        self.rototomo_reg, self.rtomo_n_harm, self.rtomo_n_svd, 
                        self.eqm, self.rho_lbl, self.rtomo_show_contours)

        self.rototomo_m_num.currentIndexChanged.connect(self.roto_tomo.update_node_m_number)
        self.rototomo_n_num.currentIndexChanged.connect(self.roto_tomo.update_node_n_number)                    
        self.rototomo_plot_bcg.stateChanged.connect(self.roto_tomo.set_substract)
        self.rototomo_TeOver.stateChanged.connect(self.roto_tomo.TeOverplot)
        self.rototomo_mag_flx.stateChanged.connect(self.roto_tomo.show_magflx)
        self.rototomo_limit.valueChanged.connect(self.roto_tomo.set_plot_lim)
        self.rototomo_reg.valueChanged.connect(self.roto_tomo.set_reg)
        self.retrofit_button.clicked.connect(self.roto_tomo.show_retrofit)
        self.updatePanels(0)
            
        
    def handleButtonData(self):

        from roto_tomo.roto_tomo import  DataSettingWindow
        win = DataSettingWindow(self, self.roto_tomo)
        win.show()

        
    def create_status_bar(self):

        self.status_text = QLabel("")
        self.statusBar().addWidget(self.status_text, 1)

        
    def create_menu(self):

        self.file_menu = self.menuBar().addMenu("&File")
        
        save_plot_action = self.create_action("&Save plot", 
            shortcut="Ctrl+S", slot=self.save_plot, 
            tip="Save the plot")
        
        self.save_anim_action = self.create_action("Save &movie", 
            shortcut="Ctrl+M", slot=self.save_anim, 
            tip="Save the animation", checkable=False)
        
        self.save_anim_action.setEnabled(False)
        
        save_data_action = self.create_action("Save &data", 
            shortcut="Ctrl+D", slot=self.save_data, 
            tip="Save the data")
        
        quit_action = self.create_action("&Quit", slot=self.close, 
            shortcut="Ctrl+Q", tip="Close the application")
        
        self.add_actions(self.file_menu, 
            (save_plot_action, save_data_action, self.save_anim_action, 
             None, quit_action))

        self.help_menu = self.menuBar().addMenu("&Help")
        
        manual_action = self.create_action("&Manual", 
            shortcut='F1', slot=self.on_manual, 
            tip='About this tool')
        
        about_action = self.create_action("&About", 
             slot=self.on_about, 
            tip='About this tool')
        
        self.add_actions(self.help_menu, (about_action, manual_action))


    def add_actions(self, target, actions):

        for action in actions:
            if action is None:
                target.addSeparator()
            else:
                target.addAction(action)
                
                
    def setCenter(self):

        frect = QDesktopWidget.frameGeometry(self)
        frect.moveCenter(QDesktopWidget().availableGeometry(self).center());
        self.move(frect.topLeft())


    def create_action(  self, text, slot=None, shortcut=None, 
                        icon=None, tip=None, checkable=False):

        action = QAction(text, self)
        if icon is not None:
            action.setIcon(QIcon(":/%s.png" % icon))
        if shortcut is not None:
            action.setShortcut(shortcut)
        if tip is not None:
            action.setToolTip(tip)
            action.setStatusTip(tip)
        if slot is not None:
            action.triggered.connect(slot)
        if checkable:
            action.setCheckable(True)
        return action
    
 
    def load_config(self):
        
        config = configparser.RawConfigParser()
        
        cfg_path = os.path.expanduser('~/pyspecviewer.cfg')
        
        #if the personal cfg file does not exust, copy one from the code directory
        if not os.path.isfile(cfg_path):
            #any better way how to distinguis in betwen different instalations? 
            if os.path.exists('/fusion/projects/codes/pyspecview'):
                tok_name = 'DIIID'
            elif os.path.exists('/afs/ipp-garching.mpg.de/'):
                tok_name = 'AUG'
            elif os.path.exists('/p/'):
                tok_name = 'NSTX'
            else:
                tok_name = 'DIIID'
                print('Local installation of pySpecview - DIII-D?')
                
            from shutil import copyfile
            path = os.path.dirname(os.path.realpath(__file__))
            copyfile( path+os.sep+'pyspecviewer_'+tok_name+'.cfg', cfg_path)
  
        config.read(cfg_path)
    
        self.tokamak = config.get('basic', 'tokamak', fallback = 'AUG')
        self.mds_server = config.get('basic', 'mds_server', fallback = None)

        self.spect_cmap = config.get('spectrogram', 'spect_cmap')
        self.win_fun = config.get('spectrogram', 'win_fun')
        self.show_plasma_freq = config.get('spectrogram', 'show_plasma_freq').strip() == 'True'
        self.show_raw = config.get('spectrogram', 'show_raw').strip() == 'True'

        self.Te2D_n_harm = int(config.get('Diag2DMapping', 'Te2D_n_harm'))
        self.n_contour_2Dprof = int(config.get('Diag2DMapping', 'n_contour'))

        self.eq_diag = config.get('equilibrium', 'eq_diag')
        self.eq_ed = int(config.get('equilibrium', 'eq_ed'))
        self.eq_exp = config.get('equilibrium', 'eq_exp')
        self.rho_lbl = config.get('equilibrium', 'radial_coord').strip()
        self.dR_corr = float(config.get('equilibrium', 'dR'))
        self.dZ_corr = float(config.get('equilibrium', 'dz'))

        self.use_LFS_data = config.get('Diag2DMapping', 'use_LFS').strip() == 'True'
        self.phase_locked_tomo = config.get('Diag2DMapping', 'phase_locked').strip() == 'True'
        self.rtomo_n_svd  =  int(config.get('roto_tomo', 'n_svd', fallback=4))
        self.rtomo_n_harm =  int(config.get('roto_tomo', 'n_harm', fallback = 5))
        self.rtomo_show_contours =  config.get('roto_tomo', 'show_contours', fallback = 'False') == 'True'

 
            
        try:
            self.balooning_coils_mode_range =  np.int_(config.get('spectrogram', 'balooning_coils_mode_range').split(','))
        except:
            self.balooning_coils_mode_range = None
 
        self.mode_range = np.int_(config.get('spectrogram', 'mode_range', fallback='-6,5').split(','))
        self.rho_pol_mode = float(config.get('spectrogram', 'rho_pol_mode', fallback = 0.2))
       


if __name__ == '__main__':
 
    path = os.path.dirname(os.path.realpath(__file__))
    text = open(path+os.sep+'README.md').read()

    parser = argparse.ArgumentParser( usage=text)
    
    parser.add_argument('--shot', metavar='S', type=int, help='optional shot number', default=None)
    parser.add_argument('--diag', metavar='D', type=str, help='optional diag name', default=None)
    parser.add_argument('--group', metavar='G', type=str, help='optional group name', default=None)
    parser.add_argument('--sig', metavar='I', type=str, help='optional signal name', default=None)
    
    parser.add_argument('--diag_phase', metavar='P', type=str, help='optional group crossphase name', default=None)
    parser.add_argument('--signal_phase', metavar='J', type=str, help='optional signal crossphase name', default=None)

    parser.add_argument('--tmin', type=np.double, help='optional start time', default=np.nan)
    parser.add_argument('--tmax', type=np.double, help='optional end time', default=np.nan)
    parser.add_argument('--fmin', type=np.double, help='optional min. selected freq. range', default=0)
    parser.add_argument('--fmax', type=np.double, help='optional max. selected freq. range', default=np.inf)

    parser.add_argument('--save_fig', help='Save figure without openning GUI', action='store_true')

    args = parser.parse_args()

    app = QApplication(sys.argv)
    app.setStyle(QStyleFactory.create("plastique"))
    #set window icon
    app.setWindowIcon(QIcon(path+os.sep+'icon.png'))
    
    if os.name == 'nt':
        new_font = app.font()
        new_font.setPointSize(  12 )
        app.setFont( new_font )
        
    assert not (args.tmin > args.tmax)
    assert not (args.fmin > args.fmax)
    
    #try:
    form = MainGUI(args.shot, args.diag, args.group, args.sig, args.diag_phase, args.signal_phase, args.tmin, args.tmax, args.fmin, args.fmax )

    if args.save_fig:
        form.curr_tab = 1
        form.save_plot(no_gui=True)
    else:
        form.show()
        app.exec_()
    #except Exception as e:
       #print(e)
       #exit()
