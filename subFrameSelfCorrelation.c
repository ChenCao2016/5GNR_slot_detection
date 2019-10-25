#include "subFrameSelfCorrelation.h"


int subframeCorrelate(double *correlationArray, 
                        float * data, 
                        unsigned long long numSamples, 
                        unsigned long samplingRateRatio, 
                        unsigned long numerology,
                        int mode)
{

    const unsigned long FFTsizeExp = 12;
    const unsigned long numSamplePerSymbol = 4096;
    const unsigned long numSamplePerCP = 144 * 2;
    unsigned long numSamplePerCP0 = numSamplePerCP + (16 << (1 + numerology));

    unsigned long numSamplePerSymbolRev = numSamplePerSymbol * samplingRateRatio << 1;
    unsigned long numSamplePerCPrev = numSamplePerCP * samplingRateRatio << 1;
    unsigned long numSamplePerCP0Rev = numSamplePerCP0 * samplingRateRatio << 1;
    unsigned long stepSize = samplingRateRatio << 1;

    unsigned long correlationLength = (mode == 0)? numSamplePerCP0Rev : numSamplePerCPrev;

    double sumProductRealCross = 0;
    double sumProductImageCross = 0;
    double sumProductRealSelf = 0;

    unsigned long n = 0;


    for (unsigned long i = 0; i < correlationLength; i += stepSize){

        sumProductRealCross += data[i] * data[i + numSamplePerSymbolRev] 
                                + data[(i + 1)]*data[(i + 1) + numSamplePerSymbolRev];

        sumProductImageCross += - data[i] * data[(i + 1) + numSamplePerSymbolRev] 
                                + data[(i + 1)] * data[i + numSamplePerSymbolRev];

        sumProductRealSelf +=  data[i] * data[i] + data[(i + 1)]*data[(i + 1)];

        correlationArray[n] = 0;

        n = n + 1;

    }

    for (unsigned long long i = correlationLength; i < (numSamples - numSamplePerSymbolRev - 1); i += stepSize){

        sumProductRealCross += data[i] * data[i + numSamplePerSymbolRev] 
                                + data[(i + 1)]*data[(i + 1) + numSamplePerSymbolRev];

        sumProductImageCross += - data[i] * data[(i + 1) + numSamplePerSymbolRev] 
                                + data[(i + 1)] * data[i + numSamplePerSymbolRev];

        sumProductRealCross -= data[i - correlationLength] * data[i - correlationLength + numSamplePerSymbolRev] 
                                + data[(i + 1) - correlationLength]*data[(i + 1) - correlationLength + numSamplePerSymbolRev];
        
        sumProductImageCross -= - data[i - correlationLength] * data[(i + 1) - correlationLength + numSamplePerSymbolRev] 
                                + data[(i + 1)  - correlationLength] * data[i - correlationLength + numSamplePerSymbolRev];        

        sumProductRealSelf +=  data[i] * data[i] + data[(i + 1)]*data[(i + 1)];

        sumProductRealSelf -=  data[i - correlationLength] * data[i - correlationLength] 
                                + data[(i + 1) - correlationLength]*data[(i + 1) - correlationLength];

        correlationArray[n] = (sumProductRealCross*sumProductRealCross + sumProductImageCross*sumProductImageCross)
                                /(sumProductRealSelf * sumProductRealSelf);
    
        n = n + 1;
    }        



    return 1;
  
}



