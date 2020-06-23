import time
from scipy.fft import fft, fftfreq, fftshift, ifft
import scipy.signal as signal

import numpy as np
import array

from matplotlib import pyplot

from scpi.lp_socket import *
from IQsamplesIO import *

from sUtility import *
from DMRSgeneration import *
import NRparameters

from measureLib import *


if __name__ == "__main__":


    if True:
        logger.logLevel = 4
        logger.enablePrint = True
        logger.start()
        lp = connect_tester("10.201.13.138")
        interface = IO(lp,logger)
        #lp.exec("CHAN1;VSA11;INIT;")
        lp.query("*WAI;ERR:ALL?")
        IQsample = interface.readIQsample(1)
        data = np.array(IQsample[1]) + 1j*np.array(IQsample[2])
        logger.stop()
    else:
        #-----------------------------------------------------------------
        # load data
        #filename = "qcom_mu3_bw100_66@0.dat"
        filename = "0.dat"
        #filename = "mu3_bw100_qpsk_tfpc0_rbo000_rbd066.dat"

        f = open(filename, 'rb')
        packData = f.read()
        f.close()
        IQsample = np.array(array.array('f',packData))
        data = IQsample[0::2] + 1j* IQsample[1::2]
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
    
    for i in range(p.symbol_per_slot*2):

        cp_start = np.argmax(c[i*(p.symbol_sample + p.normal_cp_sample): (i+1)*(p.symbol_sample + p.normal_cp_sample)])
        phase_shift = np.multiply( data[cp_start:cp_start+p.normal_cp_sample], np.conj(data[cp_start + p.symbol_sample:cp_start+p.symbol_sample+p.normal_cp_sample]))
        phase_shift_ave = np.average(phase_shift)
        freqOffset = freqOffset + np.angle(phase_shift_ave)/p.symbol_sample

    freqOffset = freqOffset/p.symbol_per_slot/2

    print("Estimated frequency offset: " + str(freqOffset/2/np.pi*para.expectedSamplingRate) + "Hz")

    #---------------------------------------------------------------------
    # frequency offset compensation

    freqComp = np.sin(freqOffset*np.array(range(np.size(data)))) - 1j* np.sin(np.pi/2 + freqOffset*np.array(range(np.size(data))))
    data = np.multiply(data,freqComp)

    phase_shift = np.multiply( data[cp_start:cp_start+p.normal_cp_sample], np.conj(data[cp_start + p.symbol_sample:cp_start+p.symbol_sample+p.normal_cp_sample]))
    phase_shift_ave = np.average(phase_shift)
    freqOffset = np.angle(phase_shift_ave)/p.symbol_sample

    print("Frequency offset adjust: " + str(freqOffset/np.pi*para.expectedSamplingRate) + "Hz")

    #---------------------------------------------------------------------
    # uplink DMRS correlation

    a = UL_DMRS()
    a.transformPrecoding = False
    a.N_Id_n_SCID = 0
    a.antenna_port = 0


    cFloor = []
    cRoof = []
    slot = []

    for i in range(p.slot_per_frame):

        #RE value should be always generated with max number RB 
        RE = a.REvalue(p.sc_per_rb*para.maxRb,i,para.dmrsSymb)

        reArrange = a.fillShift(RE,1<<para.FFTsizeExp,para.maxRb,para.rbNum,para.rbOffset)

        data  = (data - np.mean(data))/ np.std(data)
        refSignal = ifft(reArrange)
        
        c = np.absolute(signal.correlate(data,refSignal,mode = 'valid'))/30*(273/12/4)

        if np.max(c) > 0.7:
            cRoof.append(np.max(c))
            slot.append(i)
            pyplot.plot(c[::])
        else:
            cFloor.append(np.max(c))

        print(str(i) + " : " +  str(np.max(c)))
    
    handle1.show()
    

    #---------------------------------------------------------------------
    # DMRS channel estimation
    
    RE = a.REvalue(p.sc_per_rb*para.maxRb,slot[0],para.dmrsSymb) #RE value on constellation is unit power

    reArrange = a.fillShift(RE,1<<para.FFTsizeExp,para.maxRb,para.rbNum,para.rbOffset)

    data  = (data - np.mean(data))/ np.std(data)
    refSignal = ifft(reArrange)
    c = np.absolute(signal.correlate(data,refSignal,mode = 'valid'))/30
    symbol_start = np.argmax(c)
    
    res = fft(data[symbol_start:symbol_start + p.symbol_sample]) 

    DMRS_symbol = a.ifillShift(res,1<<para.FFTsizeExp,para.maxRb,para.rbNum,para.rbOffset)

    h = channel_estimate(DMRS_symbol, RE, para.rbNum, para.rbOffset)
    
    handle2 = pyplot.figure()
    pyplot.plot(np.angle(h))
    handle2.show()

    handle3 = pyplot.figure()
    pyplot.plot(np.absolute(h),'-x')
    handle3.show()


    #---------------------------------------------------------------------
    # DMRS symbol EVM

    # csymbol = DMRS_symbol/h

    # RE_cSymbol = csymbol[np.real(RE) != 0]
    # res = evm_estimate(RE_cSymbol,'QPSK')
    # handle3 = pyplot.figure()
    # pyplot.plot(res,'-x')
    # handle3.show()

    # print("RE EVM: " + str(np.average(res*100)) + "%")


    # Data_cSymbol = csymbol[np.real(RE) == 0]
    # res = evm_estimate(Data_cSymbol,'QAM')
    # handle3 = pyplot.figure()
    # pyplot.plot(res,'-x')
    # handle3.show()

    # print("data EVM: " + str(np.average(res*100)) + "%")

    # handle4 = pyplot.figure()
    # pyplot.scatter(np.real(csymbol),np.imag(csymbol))
    # pyplot.scatter(np.real(CONSTELLATION['QPSK']),np.imag(CONSTELLATION['QPSK']))
    # pyplot.scatter(np.real(CONSTELLATION['QAM']),np.imag(CONSTELLATION['QAM']))
    # handle4.show()

    #-------------------------------------------------------------------------
    # all symbol EVM

    evm_all = []

    for index in range(p.symbol_per_slot):

        index_offset = index - para.dmrsSymb

        data_symbol_start = symbol_start + (p.symbol_sample + p.normal_cp_sample) * index_offset   # next symbol

        print(data_symbol_start)

        res = fft(data[data_symbol_start:data_symbol_start + p.symbol_sample]) 

        data_symbol = a.ifillShift(res,1<<para.FFTsizeExp,para.maxRb,para.rbNum,para.rbOffset)

        csymbol = data_symbol/h

        res = evm_estimate(csymbol,['64QAM','QPSK'])
        evm_all.append(np.copy(res))

        print("symbol " + str(index) + " data EVM: " + str(np.average(res*100)) + "%")

        # handle3 = pyplot.figure()
        # pyplot.plot(res,'-x')
        # handle3.show()
        # handle4 = pyplot.figure()
        # pyplot.scatter(np.real(csymbol),np.imag(csymbol))
        # handle4.show()

    handle1 = pyplot.figure()
    pyplot.subplot(211)
    cs = pyplot.contourf(np.array(evm_all))
    pyplot.subplot(212)
    cs = pyplot.contourf(10*np.log10(np.abs(np.array(evm_all))))
    pyplot.subplots_adjust(bottom=0.1, right=0.8, top=0.9)
    cax = pyplot.axes([0.85, 0.1, 0.03, 0.8])       
    pyplot.colorbar(cax=cax)
    pyplot.grid(True)
    handle1.show()

    input()
