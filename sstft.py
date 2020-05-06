#!/usr/bin/env python 
# -*- coding: utf-8 -*-

import numpy as np
from time import time 
from scipy.ndimage.interpolation import zoom as img_zoom
from scipy.fftpack import rfft,ifft
import scipy.sparse
from multiprocessing import cpu_count

import pyfftw
pyfftw.interfaces.cache.enable()
pyfftw.interfaces.cache.set_keepalive_time(30)

try:
    from scipy.fftpack import next_fast_len 
except:
    try:
        from scipy.fftpack.helper import next_fast_len
    except:
        from scipy.signal.signaltools import _next_regular as next_fast_len
          
"""

sparse short time transformation 
Sparse time-frequency representations
Timothy J. Gardner and Marcelo O. Magnasco

"""

__all__ = ["sstft"]


def histogram2d_fast(x, y,weights,xmin,xmax,ymin,ymax, bins=10):
    #fast 2d histogram based on coo_matrix function
    #dataout of the boundaries are wrapper around
    try:
        nx, ny = bins
    except TypeError:
        nx = ny = bins

    dx = (xmax - xmin) / (nx - 1.)
    dy = (ymax - ymin) / (ny - 1.)

    ind = (x>=xmin)&(x<=xmax)&(y>=ymin)&(y<=ymax)
    x,y,weights = x[ind],y[ind],weights[ind]

    # Basically, this is just doing what digitize does with one less copy
    
    x -= xmin
    x /= dx
    y -= ymin
    y /= dy
    x = x.astype('int16')
    y = y.astype('int16')
    
    # Now, we'll exploit a sparse coo_matrix to build the 2D histogram...
    grid = scipy.sparse.coo_matrix((weights, (x,y)), shape=(nx, ny)).toarray()

    return grid, np.linspace(xmin, xmax, nx), np.linspace(ymin, ymax, ny)


def close2pow(n,round_fun=np.ceil):

    if n < 1: return 1        
    return 2**int(round_fun(np.log2(n)))


def sstft(tvec,sig,nfft,tmin,tmax,fmin,fmax,width=500,height=500,zoom=8):

    #print 'computing sparse STFT, be patient!'
    f_nfft = nfft
    if sig.dtype in [np.cdouble, np.csingle]:
        raise Exception('Complex signals are not supported')

    nfft = next_fast_len(int(nfft)//8)
    zoom = int(zoom)
    iimin, iimax = tvec.searchsorted((tmin,tmax))
    tmin = tvec[iimin]
    tmax = tvec[iimax-1]
    
    dt = (tmax-tmin)/(iimax-iimin)
    
    sigma = 1.5
    dtype = 'single'

    fmax = min(fmax, 1/dt/2)

    sig = sig[iimin:iimax].astype(dtype ,copy=False)
  
    nf = next_fast_len(len(sig))
    
    #fast trick how to get a full fft from a real rfft 
    fx = rfft(sig,nf)  #use scipy.fftpack
    
    if nf%2 == 0:
        fx = np.r_[fx[0],fx[1:-1:2]+1j*fx[2:-1:2],fx[-1],fx[-3: :-2]-1j*fx[-2:1:-2]]
    else:
        fx = np.r_[fx[0],fx[1:-1:2]+1j*fx[2:  :2],       fx[-2:0:-2]-1j*fx[  :1:-2]]

    nf = len(fx)
    fx[0] = 0  #remove offset

    #factor for interpoaltion by FFT
    fft_zoom = max(1,close2pow(height*zoom/((len(sig)*nfft//2)//nf)//2,np.round))
    
    nalign = pyfftw.simd_alignment
    # input signal
    out_step = pyfftw.n_byte_align_empty((2,nfft*fft_zoom), nalign,dtype=np.complex64)
    # output signal
    in_step  = pyfftw.n_byte_align_empty((2,nfft*fft_zoom), nalign,dtype=np.complex64)
    fft_backward = pyfftw.FFTW(in_step,out_step, direction='FFTW_BACKWARD',
                flags=['FFTW_ESTIMATE','FFTW_DESTROY_INPUT'],axes=(-1,), threads=cpu_count()/2)
        
    ifmin = int(max(nf*fmin*(2*dt),  0))//2
    ifmax = int(min(nf*fmax*(2*dt), nf))//2
    
    fft_step = max(1,(ifmax-ifmin)//(height*zoom))
    fvec = np.arange(ifmin,ifmax,fft_step)/dt/2./nf*2

    t = np.arange(nfft)
    w = 2*np.pi*np.arange(nfft)/float(nfft)

    #fourier transformation of the gaussian window
    fgwin1 = np.exp(-(w-w[nfft//2])**2*(sigma**2)/2.)/sigma*np.sqrt(nfft)
    fgwin1 = fgwin1.astype(dtype)

    #fourier transformation of the gaussian window and derivation
    fgwin2 = (t[nfft//2]-t)*np.exp(-(w-w[nfft//2])**2*(sigma**2)/2.)*sigma/np.sqrt(nfft)*6.3
    fgwin2 = fgwin2.astype(dtype)*-1j

    nt = (len(sig)*nfft//2)//nf*2*fft_zoom
    w_int   = np.zeros((len(fvec),nt),dtype=dtype)
    t_int   = np.zeros((len(fvec),nt),dtype=dtype)
    weights = np.zeros((len(fvec),nt),dtype=dtype)

    t = np.linspace(0,nfft,nfft*fft_zoom,endpoint=False)[:nt]
    w = fvec*dt*(2*np.pi)*nf/nfft

    #calculate STFT in fourier domain and also chi and eta 
    for i,ii in enumerate(range(ifmin,ifmax,fft_step)):

        ind = slice(max(0,ii-nfft//2)-ii+nfft//2,nfft//2-ii+min(nf,ii+nfft//2))
        n = ind.stop-ind.start

        in_step[:, :n]  = fx[max(0,ii-nfft//2):min(ii+nfft//2,nf)]
        in_step[0, :n] *= fgwin1[ind]
        in_step[1, :n] *= fgwin2[ind]
        in_step[:, n:]  = 0

        fft_backward()
        chi, eta = out_step[:,:nt]

        np.abs(chi,out=weights[i])
        chi[weights[i] == 0] = 1e-6  #remove 0/0 division

        ratio = eta/chi
        w_int[i] = w[i]+ratio.imag/sigma**2
        t_int[i] = t+ratio.real
 
    #interpolate in frequency and time to reduce noise in the histogram
    zoom_fact = close2pow(width*zoom*1./len(fvec)),close2pow(height*zoom*1./nt)

    order = 1
    if zoom_fact != (1,1):
        weights_ = abs(img_zoom(weights,zoom_fact,output=np.single,order=order))
        w_int = img_zoom(w_int*weights,zoom_fact,output=np.single,order=order)
        w_int/= weights_
        t_int = img_zoom(t_int*weights,zoom_fact,output=np.single,order=order)
        t_int/= weights_
    else:
        weights_ = weights

    hist,x,y = histogram2d_fast(t_int.ravel(),w_int.ravel(),weights_.ravel(),
                            t[0],t[-1],w[0],w[-1], bins=(width,height))

    return hist, np.linspace(fvec[0],fvec[-1],len(y)), np.linspace(tmin,tmax,len(x))


if __name__ == '__main__':

    import matplotlib.pylab as plt

    data+= np.sin(tvec*50000*(1+.02*np.sin(tvec*10)))*2

    T = time()

    spec,f,t  = sstft(tvec,data,2**10,tvec[0],tvec[-1],1e3,15e3,zoom=8)
