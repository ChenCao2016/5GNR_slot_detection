import scipy.signal as signal
import numpy as np
from sUtility import *

from matplotlib import pyplot

CONSTELLATION = {
    "QPSK" : [1+1j,1-1j,-1+1j,-1-1j]/np.sqrt(2),
    "QAM"  : [
                4/3 + 4/3j, 4/3 - 4/3j, -4/3 + 4/3j, -4/3 - 4/3j,
                4/9 + 4/9j, 4/9 - 4/9j, -4/9 + 4/9j, -4/9 - 4/9j,
                4/3 + 4/9j, 4/3 - 4/9j, -4/3 + 4/9j, -4/3 - 4/9j,
                4/9 + 4/3j, 4/9 - 4/3j, -4/9 + 4/3j, -4/9 - 4/3j,
                ]/np.sqrt(2),
    "64QAM" : np.array([-15.        -15.j        , -15.        -10.71428571j,
       -15.         -6.42857143j, -15.         -2.14285714j,
       -15.         +2.14285714j, -15.         +6.42857143j,
       -15.        +10.71428571j, -15.        +15.j        ,
       -10.71428571-15.j        , -10.71428571-10.71428571j,
       -10.71428571 -6.42857143j, -10.71428571 -2.14285714j,
       -10.71428571 +2.14285714j, -10.71428571 +6.42857143j,
       -10.71428571+10.71428571j, -10.71428571+15.j        ,
        -6.42857143-15.j        ,  -6.42857143-10.71428571j,
        -6.42857143 -6.42857143j,  -6.42857143 -2.14285714j,
        -6.42857143 +2.14285714j,  -6.42857143 +6.42857143j,
        -6.42857143+10.71428571j,  -6.42857143+15.j        ,
        -2.14285714-15.j        ,  -2.14285714-10.71428571j,
        -2.14285714 -6.42857143j,  -2.14285714 -2.14285714j,
        -2.14285714 +2.14285714j,  -2.14285714 +6.42857143j,
        -2.14285714+10.71428571j,  -2.14285714+15.j        ,
         2.14285714-15.j        ,   2.14285714-10.71428571j,
         2.14285714 -6.42857143j,   2.14285714 -2.14285714j,
         2.14285714 +2.14285714j,   2.14285714 +6.42857143j,
         2.14285714+10.71428571j,   2.14285714+15.j        ,
         6.42857143-15.j        ,   6.42857143-10.71428571j,
         6.42857143 -6.42857143j,   6.42857143 -2.14285714j,
         6.42857143 +2.14285714j,   6.42857143 +6.42857143j,
         6.42857143+10.71428571j,   6.42857143+15.j        ,
        10.71428571-15.j        ,  10.71428571-10.71428571j,
        10.71428571 -6.42857143j,  10.71428571 -2.14285714j,
        10.71428571 +2.14285714j,  10.71428571 +6.42857143j,
        10.71428571+10.71428571j,  10.71428571+15.j        ,
        15.        -15.j        ,  15.        -10.71428571j,
        15.         -6.42857143j,  15.         -2.14285714j,
        15.         +2.14285714j,  15.         +6.42857143j,
        15.        +10.71428571j,  15.        +15.j        ])/13.887301496588274
}


def channel_estimate(data,validRE,rbNum,rbOffset):

    # validRE = RE[12*rbOffset:12*(rbOffset + rbNum)]
    # validRE = RE
    
    assert np.size(data) == np.size(validRE)

    res = data*np.conj(validRE)

    tap = int(np.size(validRE)/np.size(validRE[np.absolute(validRE) == 0]))

    # A simple average filter
    tapWin = 1
    L = int(tap * tapWin)
    b = np.ones(L) / tapWin
    a = [1]
    h = signal.lfilter(b, a, res)
    for i in range(tapWin - 1):
        h[i*tap: (i+1)*tap] = h[i*tap: (i+1)*tap] * tapWin /(i+1)


    return h

def evm_estimate(data, cons):

    icons = []
    if type(cons) is list:
        for i in cons:
            icons = np.concatenate((icons, CONSTELLATION[i]))
    else:
        icons = CONSTELLATION[cons]

    ex = least_find(data, icons)
    error = data - ex
    res = np.sqrt(np.real(error*np.conj(error)))/np.sqrt(np.real(ex*np.conj(ex)))

    return res






    