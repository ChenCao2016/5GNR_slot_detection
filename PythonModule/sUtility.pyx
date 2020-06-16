#signal processing utilities, use cython to build

import scipy.signal as signal
import numpy as np
cimport numpy as np
cimport cython
@cython.boundscheck(False)
@cython.wraparound(False)  

#self_correlate: correlation of two segments in data with an interval 
#data: input complex signal
#cLength: segment length
#cInterval: interval between two segment
#normalize: normalize with first segment self-correlation
def self_correlate(np.ndarray[np.complex_t, ndim = 1] data, long cLength, long cInterval, normalize = False):

    cdef long rlength = data.shape[0] - cInterval - cLength

    if rlength <= 0:
        raise Exception("data size is shorter than correlation length")

    cdef np.ndarray[np.complex_t, ndim = 1] sumProduct = np.zeros(rlength, dtype=np.complex)
    cdef np.ndarray[np.complex_t, ndim = 1] sumProductReal = np.zeros(rlength, dtype=np.complex)

    cdef long n = 0

    cdef np.ndarray[np.complex_t, ndim = 1] temp = signal.correlate(data[0: cLength],
                                                                    data[cInterval: cInterval+cLength],
                                                                    mode="valid")

    sumProduct[n] = temp[0]

    if normalize:
        temp = signal.correlate(data[0: cLength],
                                data[0: cLength],
                                mode="valid")

        sumProductReal[n] = temp[0]

    n = n + 1

    cdef np.float_t sumProductRealCross,sumProductImageCross,sumProductRealSelf

    while n < rlength:

        sumProductRealCross = sumProduct[n-1].get('real')
        sumProductImageCross = sumProduct[n-1].get('imag')
        

        sumProductRealCross +=  + data[n+cLength].get('real') * data[n+cInterval+cLength].get('real') \
                                + data[n+cLength].get('imag') * data[n+cInterval+cLength].get('imag')

        sumProductImageCross += - data[n+cLength].get('real') * data[n+cInterval+cLength].get('imag') \
                                + data[n+cLength].get('imag') * data[n+cInterval+cLength].get('real')

        sumProductRealCross -=  + data[n-1].get('real') * data[n-1+cInterval].get('real') \
                                + data[n-1].get('imag') * data[n-1+cInterval].get('imag')
        
        sumProductImageCross -= - data[n-1].get('real') * data[n-1+cInterval].get('imag') \
                                + data[n-1].get('imag') * data[n-1+cInterval].get('real')

        sumProduct[n] = {'real': sumProductRealCross, 'imag': sumProductImageCross}

        if normalize:

            sumProductRealSelf = sumProductReal[n-1].get('real')

            sumProductRealSelf +=   + data[n+cLength].get('real') * data[n+cLength].get('real') + data[n+cLength].get('imag')*data[n+cLength].get('imag')

            sumProductRealSelf -=   + data[n-1].get('real') * data[n-1].get('real') + data[n-1].get('imag') * data[n-1].get('imag')

            sumProductReal[n] = {'real': sumProductRealSelf, 'imag': 0.0}



        # sumProduct[n] = sumProduct[n-1] - data[n-1] * np.conj(data[n-1+cInterval]) + data[n+cLength] * np.conj(data[n+cInterval+cLength])
        # sumProductReal[n] = sumProductReal[n-1] - data[n-1] * np.conj(data[n-1]) + data[n+cLength] * np.conj(data[n+cLength])

        n = n + 1


    cdef np.ndarray[np.complex_t] corr
    
    if normalize:
        corr = sumProduct/sumProductReal
    else:
        corr = sumProduct

    return corr
