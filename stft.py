import numpy as np
from scipy.signal import get_window
from time import time
from scipy.fftpack import fftfreq, fftshift
import pyfftw
pyfftw.interfaces.cache.enable()
pyfftw.interfaces.cache.set_keepalive_time(30)

from multiprocessing import cpu_count

#loading foe manu different version of scipy 
try:
    from scipy.fftpack import next_fast_len 
except:
    try:
        from scipy.fftpack.helper import next_fast_len
    except:
        from scipy.signal.signaltools import _next_regular as next_fast_len
        
__all__ = ["stft"]

def stft(tvec, signal,nfft,resolution=1000,window='gauss',fmin=-np.infty,fmax=np.infty,
         tmin=-np.infty, tmax=np.infty,pass_DC=True,complex_spectrum=False):
    """Short time Fourier Tranform. in time domain
        x - real or complex, transformation done along first axis
        
        nfft - length of the STFT window, it can be any number, not just power of 2!
        resolution - horizontal resolution of the spectrogram
        window - windowing function
        dt  - sampling time
        fmin - minimal frequency
        fmax - maximal frequency
        tmin - minimal time
        tmax - maximal time
        pass_DC - remove constant background
        complex_spectrum - output is complex spectrogram        
    """
     
    if signal.dtype in [np.double, np.cdouble]:
        raise Exception('use dtype single or csingle')

    complex_sig =  signal.dtype == np.csingle

    iimin, iimax = tvec.searchsorted((tmin,tmax)) #BUG SLOW for large signals 
    iimax-= 1
    n = iimax-iimin
    
    dt = (tvec[iimax]-tvec[iimin])/n  #BUG will be wrong for not equally spaced time vectors 

    f_nfft = nfft

    nfft = next_fast_len( int(nfft) )
    fvec = fftfreq(nfft,dt)
    if not complex_sig:
        fvec = fvec[:nfft//2]
        sig_dtype = 'single'
    else:
        fvec = fftshift(fvec)
        sig_dtype = 'csingle'

    ifmin,ifmax = fvec.searchsorted([fmin,fmax])
    ifmax = np.minimum(ifmax+1, len(fvec))

    # input signal
    nalign = pyfftw.simd_alignment
    sig = pyfftw.n_byte_align_empty(signal.shape[1:]+(nfft,), nalign,dtype=sig_dtype)
    # output signal
    nf = nfft if complex_sig else nfft//2+1
    out = pyfftw.n_byte_align_empty(signal.shape[1:]+(nf  ,), nalign,dtype=np.complex64)
    
    fft_forward = pyfftw.FFTW(sig,out, direction='FFTW_FORWARD',
                flags=['FFTW_ESTIMATE','FFTW_DESTROY_INPUT'],axes=(-1,),threads=cpu_count()//2)
    
    fft_step = max(1,n//resolution)
    dtype =  np.complex64 if complex_spectrum else np.single
    spec = np.empty((int(n/fft_step),ifmax-ifmin)+signal.shape[1:],dtype=dtype)
    win = None

    for i in range(int(n//fft_step)):
        imin = int(max(0,iimin+np.floor(i*fft_step-nfft//2)))
        imax = int(min(len(signal),iimin+np.floor(i*fft_step+nfft//2)))        
        
        if np.size(win) != imax-imin:
            if window == 'gauss':
                win = get_window((window, f_nfft*(imax-imin)//8//nfft),imax-imin, fftbins=False)
            else:
                win = get_window(window, imax-imin, fftbins=False)
                
            win = win.astype('single')
            
        #implicit conversion from int (or anything else) to single!
        if pass_DC:
            sig[...,:imax-imin] = signal[imin:imax].T
        else:
            sig[...,:imax-imin] = signal[imin:imax].T
            sig[...,:imax-imin] -= np.mean(signal[imin:imax],0)[None].T
            
        sig[...,imax-imin:] = 0
        sig[...,:imax-imin]*= win 

        #the main step, FFT
        fft_forward()  
        
        #prepare output spectrum
        if complex_sig: 
            #to get monotonously increasing frequency
            ind_neg = slice(min(nfft,ifmin+(nfft+1)//2),min(ifmax+(nfft+1)//2,nfft))
            ind_pos = slice(max(0,ifmin-(nfft)//2),max(ifmax-(nfft)//2,0))
            nneg = ind_neg.stop-ind_neg.start
            npos = ind_pos.stop-ind_pos.start
            if complex_spectrum:
                spec[i,:nneg] = out.T[ind_neg]
                spec[i,nneg:] = out.T[ind_pos]
            else:#compute abs(fft)
                np.abs(out.T[ind_neg],out=spec[i,:nneg])
                np.abs(out.T[ind_pos],out=spec[i,nneg:])
        else:
            if complex_spectrum:
                spec[i] = out.T[ifmin:ifmax]
            else:
                np.abs(out.T[ifmin:ifmax],out=spec[i])  #compute abs(fft)

    return spec,fvec[ifmin:ifmax],tvec[iimin:iimax:fft_step][:int(n//fft_step)]


if __name__ == '__main__':


    import matplotlib.pylab as plt
 
    N = 62*1600
    x = np.zeros((N))
    np.random.seed(1000)

    x[N//4] = 1
    x = (np.fft.fft((x),axis=0))

    x = np.single(x)
    tvec = np.arange(len(x))
    T = time()

    spec,f,t = stft(tvec, x,2**10,resolution=10000,window='hann',fmin=-.1,fmax=.1,
            tmin=-np.infty, tmax=np.infty,pass_DC=True, complex_spectrum=True)
    
    plt.imshow(abs(spec.T),origin='lower', extent=(t[0],t[-1],f[0],f[-1]),aspect='auto', interpolation='nearest')
    plt.colorbar()
    plt.show()

    plt.figure()
    plt.imshow(angle(spec.T),origin='lower', extent=(t[0],t[-1],f[0],f[-1]),aspect='auto', interpolation='nearest')
    plt.colorbar()
    plt.show()
