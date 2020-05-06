from scipy.signal import get_window
import numpy as np

from scipy.fftpack import fftfreq, rfft, ifft, fftshift
from time import time 
import matplotlib.pylab as plt
from multiprocessing import cpu_count

try:
    from scipy.fftpack import next_fast_len 
except:
    try:
        from scipy.fftpack.helper import next_fast_len
    except:
        from scipy.signal.signaltools import _next_regular as next_fast_len

__all__ = ["sfft"]
import pyfftw
pyfftw.interfaces.cache.enable()
pyfftw.interfaces.cache.set_keepalive_time(30)


def sfft(tvec, x,nfft,resolution=1000,window='hann',fmin=0,fmax=np.infty,
         tmin=-np.infty, tmax=np.infty,pass_DC=True,complex_spectrum=False):
    """Short time Fourier Tranform. in the frequency domain
        done along 1. axis!
    """
    if x.dtype in [np.cdouble, np.csingle]:
        raise Exception('Complex signals are not supported yet')

    f_nfft = nfft
    nfft = next_fast_len( int(nfft) )

    iimin = tvec.searchsorted(tmin)
    iimax = tvec.searchsorted(tmax)

    tmin = tvec[iimin]
    tmax = tvec[iimax-1]
    dt = (tmax-tmin)/(iimax-iimin)

    fmax = min(fmax, 1./dt/2)

    nt = next_fast_len(len(x))

    x = x[iimin:iimax]

    complex_sig =  x.dtype == np.csingle

    if not complex_sig:
        sig_dtype = 'single'
        nf =  nt//2+1
    else:
        sig_dtype = 'csingle'
        nf =  nt

    nalign = pyfftw.simd_alignment
    # input signal
    sig = pyfftw.n_byte_align_empty(x.shape[1:]+(nt,), nalign,dtype=sig_dtype)
    # output signalsignal
    out = pyfftw.n_byte_align_empty(x.shape[1:]+(nf,), nalign,dtype=np.complex64)
    fft_forward = pyfftw.FFTW(sig,out, direction='FFTW_FORWARD',
                flags=['FFTW_ESTIMATE','FFTW_DESTROY_INPUT'],axes=(-1,),threads=cpu_count()//2)

    # input signal
    out_step = pyfftw.n_byte_align_empty(x.shape[1:]+(nfft,), nalign,dtype=np.complex64)
    # output signal
    in_step  = pyfftw.n_byte_align_empty(x.shape[1:]+(nfft,), nalign,dtype=np.complex64)
    fft_backward = pyfftw.FFTW(in_step,out_step, direction='FFTW_BACKWARD',
                flags=['FFTW_ESTIMATE','FFTW_DESTROY_INPUT'],axes=(-1,),threads=cpu_count()//2)
    
    sig[...,:len(x)] = x
    sig[...,:len(x)] -= x.mean(0).T
    sig[...,len(x):] = 0

    fft_forward()
    
    del sig
    
    ifmin = int(max(nf*fmin*(2*dt),-nf))
    ifmax = int(min(nf*fmax*(2*dt), nf))

    fft_step = max(1,(ifmax-ifmin)//resolution//2)
    fvec = np.arange(ifmin,ifmax,fft_step)/dt/2./nf
    ntvec = (len(x)*nfft//2)//nf
    tvec = np.linspace(tmin,len(x)*tmax/nf,ntvec,endpoint=False)

    dtype =np.complex64 if complex_spectrum else  np.single
    spec = np.empty((ntvec, len(fvec))+x.shape[1:],dtype=dtype)

    if window == 'gauss':
        window = window,f_nfft/8.

    win = get_window(window, nfft, fftbins=False).astype('single')

    for i,ii in enumerate(range(ifmin,ifmax,fft_step)):
        
        L = min(ii+nfft//2,nf) - max(0,ii-nfft//2)
        in_step[...,:L] = out[...,max(0,ii-nfft//2):min(ii+nfft//2,nf)]
        in_step[...,L:] = 0
        in_step[...,:L] *= win[max(0,ii-nfft//2)-ii+nfft//2:nfft//2-ii+min(nf,ii+nfft//2)]

        #main step FFT 
        fft_backward()
        
        if complex_spectrum:
            spec[:,i] =  out_step.T[:ntvec]
        else:
            spec[:,i] = np.abs(out_step.T[:ntvec])

    if not pass_DC:
        spec = spec[:,1:]
        fvec = fvec[1:]

    return spec,fvec,tvec


if __name__ == "__main__":

    import matplotlib.pylab as plt

    N = 2**14+1
    x = np.zeros((N,10))
    np.random.seed(1000)
    x[N//4] = 1
    x = np.real(np.fft.fft((x),axis=0))
    x += np.random.randn(*x.shape)

    tvec = np.arange(len(x))
    T = time()

    spec,f,t = sfft(tvec, x,2**15+1,resolution=1000,window='hann',dt=1.,fmin=0,fmax=np.infty,
            tmin=-np.infty, tmax=np.infty,pass_DC=True)   
