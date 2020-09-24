import time
from scipy.fft import fft, fftfreq, fftshift, ifft
import scipy.signal as signal

import numpy as np
import array

from matplotlib import pyplot

from sUtility import *
from DMRSgeneration import *
import NRparameters

from measureLib import *

import configure


if __name__ == "__main__":

    #-----------------------------------------------------------------
    # load data
    filename = configure.file_name

    f = open(filename, 'rb')
    packData = f.read()
    f.close()
    IQsample = np.array(array.array('f',packData),dtype=np.complex)
    data = np.add(IQsample[0::2], 1j*IQsample[1::2])
    #-----------------------------------------------------------------

    #-------------------------------------------------------------------
    # parameters
    # FR2 example
    class para:
        dmrsSymb = 3
        numerology = 3
        rbNum = 66
        rbOffset = 0
        maxRb = 66
        FFTsizeExp = 10
        waveformSamplingRate = 2457600000
        scs = (1<<numerology) * 15000
        expectedSamplingRate = scs << FFTsizeExp
        samplingRateRatio = waveformSamplingRate/expectedSamplingRate
    # FR1 example
    class para:
        dmrsSymb = 3
        numerology = 1
        rbNum = 12
        rbOffset = 0
        maxRb = 273
        FFTsizeExp = 12
        waveformSamplingRate = 240000000
        scs = (1<<numerology) * 15000
        expectedSamplingRate = scs << FFTsizeExp
        samplingRateRatio = waveformSamplingRate/expectedSamplingRate

    class para:
        dmrsSymb = configure.dmrs_symbol
        numerology = configure.numerology
        rbNum = configure.rb_number
        rbOffset = configure.rb_offset
        maxRb = configure.max_rb_number
        FFTsizeExp = 12
        waveformSamplingRate = configure.sampling_rate_hz
        scs = (1<<numerology) * 15000
        expectedSamplingRate = scs << FFTsizeExp
        samplingRateRatio = waveformSamplingRate/expectedSamplingRate

    #-------------------------------------------------------------------
    # IQdata info and re-sample

    numSamples = np.size(data)

    print("Capture Length: " + str(numSamples/240))
    print("Number of samples: " + str(numSamples))

    print("Ratio: " + str(para.samplingRateRatio))
    print("Sampling Rate: " + str(int(para.expectedSamplingRate)))
    data = signal.resample(data, int(numSamples/para.samplingRateRatio))
    print("Number of samples: " + str(np.size(data)))

    #-------------------------------------------------------------------
    # CP self correlation

    p = NRparameters.unit(para.numerology,para.expectedSamplingRate)

    start = time.time()
    c = np.absolute(self_correlate(data, p.normal_cp_sample, p.symbol_sample, True))
    print("self_correlation time: " + str(time.time()-start))

    handle1 = pyplot.figure()
    pyplot.plot(c)

    #---------------------------------------------------------------------
    # frequency offset estimation, resolution (-pi,pi)

    freqOffset = 0
    
    for i in range(p.symbol_per_slot):

        cp_start = np.argmax(c[i*(p.symbol_sample + p.normal_cp_sample): (i+1)*(p.symbol_sample + p.normal_cp_sample)])
        phase_shift = np.multiply( data[cp_start:cp_start+p.normal_cp_sample], np.conj(data[cp_start + p.symbol_sample:cp_start+p.symbol_sample+p.normal_cp_sample]))
        phase_shift_ave = np.average(phase_shift)
        freqOffset = freqOffset + np.angle(phase_shift_ave)/p.symbol_sample

    freqOffset = freqOffset/p.symbol_per_slot/2

    print("Estimated frequency offset: " + str(freqOffset/2/np.pi*para.expectedSamplingRate) + "Hz")

    #freqOffset = 2830 / para.expectedSamplingRate * 2 * np.pi
    #---------------------------------------------------------------------
    #frequency offset compensation

    freqComp = np.sin(freqOffset*np.array(range(np.size(data)))) - 1j* np.sin(np.pi/2 + freqOffset*np.array(range(np.size(data))))
    data = np.multiply(data,freqComp)

    phase_shift = np.multiply( data[cp_start:cp_start+p.normal_cp_sample], np.conj(data[cp_start + p.symbol_sample:cp_start+p.symbol_sample+p.normal_cp_sample]))
    phase_shift_ave = np.average(phase_shift)
    freqOffset = np.angle(phase_shift_ave)/p.symbol_sample

    print("Frequency offset adjust: " + str(freqOffset/np.pi*para.expectedSamplingRate) + "Hz")

    #---------------------------------------------------------------------
    # uplink DMRS correlation

    a = UL_DMRS()
    a.transformPrecoding = configure.DFTs_OFDM
    a.N_Id_n_SCID = configure.N_Id_n_SCID
    a.antenna_port = configure.antenna_port
    a.n_RS_ID = configure.n_RS_ID


    cFloor = []
    cRoof = []
    slot = []

    print("uplink DMRS correlation result:")

    for i in range(p.slot_per_frame):

        RE = a.REvalue(para.maxRb,para.rbNum,i,para.dmrsSymb) #RE value on constellation is unit power

        reArrange = a.fillShift(RE,1<<para.FFTsizeExp,para.maxRb,para.rbNum,para.rbOffset)

        data  = (data - np.mean(data))/ np.std(data)
        refSignal = ifft(reArrange)
        
        c = np.absolute(signal.correlate(data,refSignal,mode = 'valid'))/30*(para.rbNum/12/4)

        cRoof.append(np.max(c))
        slot.append(i)

        print(str(i) + " : " +  str(np.max(c)))
    

    #---------------------------------------------------------------------
    # DMRS channel estimation
    
    RE = a.REvalue(para.maxRb,para.rbNum,0,para.dmrsSymb) #RE value on constellation is unit power
    reArrange = a.fillShift(RE,1<<para.FFTsizeExp,para.maxRb,para.rbNum,para.rbOffset)
    data  = (data - np.mean(data))/ np.std(data)
    refSignal = ifft(reArrange)
    c = np.absolute(signal.correlate(data,refSignal,mode = 'valid'))/30

    pyplot.plot(c[::])
    pyplot.title("uplink DMRS correlation result")
    handle1.show()

    symbol_start = np.argmax(c)
    
    res = fft(data[symbol_start:symbol_start + p.symbol_sample])/p.symbol_sample

    DMRS_symbol = a.ifillShift(res,1<<para.FFTsizeExp,para.maxRb,para.rbNum,para.rbOffset)

    handle2 = pyplot.figure()
    pyplot.plot(np.absolute(res))
    pyplot.title("DMRS symbol spectrum")
    handle2.show()

    h = channel_estimate(DMRS_symbol, RE, para.rbNum, para.rbOffset)
    
    handle2 = pyplot.figure()
    pyplot.plot(np.angle(h))
    pyplot.title("estimated channel phase")
    handle2.show()

    handle3 = pyplot.figure()
    pyplot.plot(np.abs(h),'-x')
    pyplot.title("estimated channel magitude")
    handle3.show()


    #---------------------------------------------------------------------
    # DMRS symbol EVM

    csymbol = DMRS_symbol/h

    RE_cSymbol = csymbol[np.real(RE) != 0]
    res = evm_estimate(RE_cSymbol,'QPSK')
    handle3 = pyplot.figure()
    pyplot.plot(res,'-x')
    pyplot.title("RE subcarrier EVM, in DMRS symbol")
    handle3.show()

    print("RE EVM: " + str(np.average(res*100)) + "%")

    Data_cSymbol = csymbol[np.real(RE) == 0]
    res = evm_estimate(Data_cSymbol, configure.modulation)
    handle3 = pyplot.figure()
    pyplot.plot(res,'-x')
    pyplot.title("Data subcarrier EVM, in DMRS symbol")
    handle3.show()

    print("Data EVM: " + str(np.average(res*100)) + "%")

    handle4 = pyplot.figure()
    pyplot.scatter(np.real(RE_cSymbol),np.imag(RE_cSymbol))
    pyplot.scatter(np.real(Data_cSymbol),np.imag(Data_cSymbol))
    pyplot.scatter(np.real(CONSTELLATION['QPSK']),np.imag(CONSTELLATION['QPSK']))
    pyplot.scatter(np.real(CONSTELLATION[configure.modulation]),np.imag(CONSTELLATION[configure.modulation]))
    pyplot.ylim(-2,+2)
    pyplot.xlim(-2,+2)
    pyplot.title("DMRS symbol constellation")
    handle4.show()

    #-------------------------------------------------------------------------
    # Data symbol EVM

    evm_all = []

    for index in range(2,3):

        index_offset = index - para.dmrsSymb

        data_symbol_start = symbol_start + (p.symbol_sample + p.normal_cp_sample) * index_offset  # next symbol

        res = fft(data[data_symbol_start:data_symbol_start + p.symbol_sample])/p.symbol_sample

        data_symbol = a.ifillShift(res,1<<para.FFTsizeExp,para.maxRb,para.rbNum,para.rbOffset)

        csymbol = data_symbol/h

        if configure.DFTs_OFDM:
            csymbol = ifft(csymbol)*np.sqrt(p.sc_per_rb*para.rbNum) 

        res = evm_estimate(csymbol,[configure.modulation])
        evm_all.append(np.copy(res))

        print("symbol " + str(index) + " data EVM: " + str(np.average(res*100)) + "%")

        handle4 = pyplot.figure()
        pyplot.grid()
        pyplot.scatter(np.real(csymbol),np.imag(csymbol))
        pyplot.ylim(-2,+2)
        pyplot.xlim(-2,+2)
        pyplot.title("symbol " + str(index) + " constellation")
        handle4.show()

        handle4 = pyplot.figure()
        pyplot.grid()
        pyplot.plot(res)
        pyplot.title("symbol " + str(index) + " subcarrier EVM")
        handle4.show()

    #-------------------------------------------------------------------------
    

    input()