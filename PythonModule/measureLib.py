import scipy.signal as signal
import numpy as np
from sUtility import *

CONSTELLATION = {
    "QPSK" : [1+1j,1-1j,-1+1j,-1-1j]/np.sqrt(2),
    "QAM"  : [
                4/3 + 4/3j, 4/3 - 4/3j, -4/3 + 4/3j, -4/3 - 4/3j,
                4/9 + 4/9j, 4/9 - 4/9j, -4/9 + 4/9j, -4/9 - 4/9j,
                4/3 + 4/9j, 4/3 - 4/9j, -4/3 + 4/9j, -4/3 - 4/9j,
                4/9 + 4/3j, 4/9 - 4/3j, -4/9 + 4/3j, -4/9 - 4/3j,
                ]/np.sqrt(2)
}


def channel_estimate(data,RE):
    
    assert np.size(data) == np.size(RE)

    res = data*np.conj(RE)

    tap = int(np.size(RE)/np.size(RE[np.absolute(RE) == 0]))

    # A simple average filter
    tapWin = 6
    L = int(tap * tapWin)
    b = np.ones(L) / tapWin
    a = [1]
    h = signal.lfilter(b, a, res)
    for i in range(tapWin - 1):
        h[i*tap: (i+1)*tap] = h[i*tap: (i+1)*tap] * tapWin /(i+1)

    return h

def evm_estimate(data, cons):

    ex = least_find(data,CONSTELLATION[cons])
    error = data - ex
    res = np.real(error*np.conj(error))/np.real(ex*np.conj(ex))

    return res






    