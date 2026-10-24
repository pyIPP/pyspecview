"""
TokEye: segmentation of coherent modes in spectrograms by a convolutional neural network
https://github.com/PlasmaControl/TokEye

The network tokeye.onnx is not a part of the repository, copy it in this directory
from https://huggingface.co/nc1/big_tf_unet/resolve/main/tokeye.onnx
It runs on CPU by onnxruntime (pip install onnxruntime).
"""

import os
import numpy as np
from multiprocessing import cpu_count

__all__ = ["tokeye", "tokeye_model"]

tokeye_model = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'tokeye.onnx')
session = None


def tokeye(spect, size=512):
    """Probability of a coherent mode in each pixel of the spectrogram
        spect - real spectrogram (frequency, time), the lowest frequency in the first row
        size - spectrogram is averaged down to between size and 2*size pixels
               in each axis, the network was trained on 512x512 images
        returns array of the same shape as spect with values between 0 and 1
    """
    global session
    if session is None:
        try:
            import onnxruntime
        except ImportError:
            raise Exception('TokEye needs onnxruntime: pip install onnxruntime')
        options = onnxruntime.SessionOptions()
        options.intra_op_num_threads = max(1, cpu_count()//2)
        options.enable_cpu_mem_arena = False
        session = onnxruntime.InferenceSession(tokeye_model, options, providers=['CPUExecutionProvider'])

    nf, nt = spect.shape
    kf, kt = max(1, nf//size), max(1, nt//size)
    x = spect[:nf//kf*kf, :nt//kt*kt].reshape(nf//kf, kf, nt//kt, kt).mean((1, 3))
    if min(x.shape) < 64:
        raise Exception('TokEye needs at least 64x64 pixels, the spectrogram has %dx%d'%x.shape)

    x = np.log1p(np.maximum(np.nan_to_num(x), 0))
    x = np.clip(x, *np.percentile(x, (1, 99)))
    x = (x-x.mean())/(x.std()+1e-6)

    #channel 0: coherent modes, channel 1: transient events
    p = session.run(None, {'spectrogram': np.single(x)[None, None]})[0][0, 0]

    #back to the original size
    p = np.repeat(np.repeat(p, kf, 0), kt, 1)
    return np.pad(p, ((0, nf-p.shape[0]), (0, nt-p.shape[1])), mode='edge')
