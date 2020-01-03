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
    unsigned long stepSize = (samplingRateRatio << 1);

    unsigned long correlationLength = (mode == 0)? numSamplePerCP0Rev : numSamplePerCPrev;

    double sumProductRealCross = 0;
    double sumProductImageCross = 0;
    double sumProductRealSelf = 0;

    unsigned long long n = 0;


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


int subframeCorrelate2(double *correlationArray, 
                        float * data, 
                        unsigned long long numSamples, 
                        unsigned long samplingRateRatio, 
                        unsigned long numerology,
                        int mode)
{

    const unsigned long FFTsizeExp = 12;                                                        //4096
    const unsigned long numSamplePerSymbol = 4096;                                              //FFT size
    const unsigned long numSamplePerCP = 144 * 2;                                               //CP length
    unsigned long numSamplePerCP0 = numSamplePerCP + (16 << (1 + numerology));                  //First symbol CP length

    unsigned long numSamplePerSymbolRev = numSamplePerSymbol * samplingRateRatio << 1;          //Real + image
    unsigned long numSamplePerCPrev = numSamplePerCP * samplingRateRatio << 1;                  //Real + image
    unsigned long numSamplePerCP0Rev = numSamplePerCP0 * samplingRateRatio << 1;                //Real + image
    unsigned long stepSize = (samplingRateRatio << 1);                                          //Real + image

    unsigned long correlationLength = (mode == 0)? numSamplePerCP0Rev : numSamplePerCPrev;      //Correlation on first symbol CP

    double sumProductRealCross = 0;
    double sumProductImageCross = 0;
    double sumProductRealSelf = 0; 

    unsigned long bufferLength = correlationLength/stepSize + 100;
    double * sumProductRealCrossBuffer = (double *) malloc(sizeof(double) * bufferLength);
    double * sumProductImageCrossBuffer = (double *) malloc(sizeof(double) * bufferLength);
    double * sumProductRealSelfBuffer = (double *) malloc(sizeof(double) * bufferLength);

    unsigned long long n = 0;
    unsigned long j = 0;

    for (unsigned long i = 0; i < correlationLength; i += stepSize){

        sumProductRealCrossBuffer[j] =  + data[i] * data[i + numSamplePerSymbolRev] 
                                        + data[(i + 1)] * data[(i + 1) + numSamplePerSymbolRev];
        sumProductRealCross +=  sumProductRealCrossBuffer[j];

        sumProductImageCrossBuffer[j] = - data[i] * data[(i + 1) + numSamplePerSymbolRev] 
                                        + data[(i + 1)] * data[i + numSamplePerSymbolRev];
        sumProductImageCross += sumProductImageCrossBuffer[j];


        sumProductRealSelfBuffer[j] =   + data[i] * data[i] 
                                        + data[(i + 1)] * data[(i + 1)];
        sumProductRealSelf += sumProductRealSelfBuffer[j];

        correlationArray[n] = 0;

        n = n + 1;
        if (j == bufferLength - 1){
            j = 0;
        }
        else{
            j = j + 1;
        }

    }

    j = 0;

    for (unsigned long long i = correlationLength; i < (numSamples - numSamplePerSymbolRev - 1); i += stepSize){

        sumProductRealCross -= sumProductRealCrossBuffer[j];
        sumProductRealCrossBuffer[j] =  + data[i] * data[i + numSamplePerSymbolRev] 
                                        + data[(i + 1)] * data[(i + 1) + numSamplePerSymbolRev];
        sumProductRealCross +=  sumProductRealCrossBuffer[j];

        sumProductImageCross -= sumProductImageCrossBuffer[j];
        sumProductImageCrossBuffer[j] = - data[i] * data[(i + 1) + numSamplePerSymbolRev] 
                                        + data[(i + 1)] * data[i + numSamplePerSymbolRev];
        sumProductImageCross += sumProductImageCrossBuffer[j];

        sumProductRealSelf -= sumProductRealSelfBuffer[j];
        sumProductRealSelfBuffer[j] =   + data[i] * data[i] 
                                        + data[(i + 1)] * data[(i + 1)];
        sumProductRealSelf += sumProductRealSelfBuffer[j];

        correlationArray[n] = (sumProductRealCross*sumProductRealCross + sumProductImageCross*sumProductImageCross)
                                /(sumProductRealSelf * sumProductRealSelf);
    
        n = n + 1;
        
        if (j == bufferLength - 1){
            j = 0;
        }
        else{
            j = j + 1;
        }

    }        

    free(sumProductRealCrossBuffer);
    free(sumProductImageCrossBuffer);
    free(sumProductRealSelfBuffer);

    return 1;
}


int subframeCorrelate3(double *correlationArray,
                        double *phaseShiftArray,
                        float * data, 
                        unsigned long long numSamples, 
                        unsigned long samplingRateRatio, 
                        unsigned long numerology,
                        int mode)
{

    const unsigned long FFTsizeExp = 12;                                                        //4096
    const unsigned long numSamplePerSymbol = 4096;                                              //FFT size
    const unsigned long numSamplePerCP = 144 * 2;                                               //CP length
    unsigned long numSamplePerCP0 = numSamplePerCP + (16 << (1 + numerology));                  //First symbol CP length

    unsigned long numSamplePerSymbolRev = numSamplePerSymbol * samplingRateRatio << 1;          //Real + image
    unsigned long numSamplePerCPrev = numSamplePerCP * samplingRateRatio << 1;                  //Real + image
    unsigned long numSamplePerCP0Rev = numSamplePerCP0 * samplingRateRatio << 1;                //Real + image
    unsigned long stepSize = (samplingRateRatio << 1);                                          //Real + image

    unsigned long correlationLength = (mode == 0)? numSamplePerCP0Rev : numSamplePerCPrev;      //Correlation on first symbol CP

    double sumProductRealCross = 0;
    double sumProductImageCross = 0;
    double sumProductRealSelf = 0; 

    unsigned long bufferLength = correlationLength/stepSize + 100;
    double * sumProductRealCrossBuffer = (double *) malloc(sizeof(double) * bufferLength);
    double * sumProductImageCrossBuffer = (double *) malloc(sizeof(double) * bufferLength);
    double * sumProductRealSelfBuffer = (double *) malloc(sizeof(double) * bufferLength);

    unsigned long long n = 0;
    unsigned long j = 0;

    for (unsigned long i = 0; i < correlationLength; i += stepSize){

        sumProductRealCrossBuffer[j] =  + data[i] * data[i + numSamplePerSymbolRev] 
                                        + data[(i + 1)] * data[(i + 1) + numSamplePerSymbolRev];

        sumProductRealCross +=          sumProductRealCrossBuffer[j];

        sumProductImageCrossBuffer[j] = - data[i] * data[(i + 1) + numSamplePerSymbolRev] 
                                        + data[(i + 1)] * data[i + numSamplePerSymbolRev];

        sumProductImageCross +=         sumProductImageCrossBuffer[j];


        sumProductRealSelfBuffer[j] =   + data[i] * data[i] 
                                        + data[(i + 1)] * data[(i + 1)];

        sumProductRealSelf +=           sumProductRealSelfBuffer[j];

        correlationArray[n] =           0;
        phaseShiftArray[n] =            0;

        n = n + 1;
        if (j == bufferLength - 1){
            j = 0;
        }
        else{
            j = j + 1;
        }

    }

    j = 0;

    for (unsigned long long i = correlationLength; i < (numSamples - numSamplePerSymbolRev - 1); i += stepSize){

        sumProductRealCross -=          sumProductRealCrossBuffer[j];
        sumProductRealCrossBuffer[j] =  + data[i] * data[i + numSamplePerSymbolRev] 
                                        + data[(i + 1)] * data[(i + 1) + numSamplePerSymbolRev];
        sumProductRealCross +=          sumProductRealCrossBuffer[j];

        sumProductImageCross -=         sumProductImageCrossBuffer[j];
        sumProductImageCrossBuffer[j] = - data[i] * data[(i + 1) + numSamplePerSymbolRev] 
                                        + data[(i + 1)] * data[i + numSamplePerSymbolRev];
        sumProductImageCross +=         sumProductImageCrossBuffer[j];

        sumProductRealSelf -=           sumProductRealSelfBuffer[j];
        sumProductRealSelfBuffer[j] =   + data[i] * data[i] 
                                        + data[(i + 1)] * data[(i + 1)];
        sumProductRealSelf +=           sumProductRealSelfBuffer[j];

        correlationArray[n] =           (sumProductRealCross*sumProductRealCross + sumProductImageCross*sumProductImageCross)
                                        /(sumProductRealSelf * sumProductRealSelf);

        phaseShiftArray[n] =            sumProductImageCross/sumProductRealCross;
    
        n = n + 1;
        
        if (j == bufferLength - 1){
            j = 0;
        }
        else{
            j = j + 1;
        }

    }        

    free(sumProductRealCrossBuffer);
    free(sumProductImageCrossBuffer);
    free(sumProductRealSelfBuffer);

    return 1;
}

int subframeStartRegionCandidate(unsigned long long *region,
                                    unsigned long long *IndexPeak,
                                    unsigned long *numCandidate, 
                                    double *correlationArray,
                                    unsigned long long correlationArrayLength, 
                                    double threshold,
                                    unsigned long numCandidateLimit)
{
    double maxCorrelation = 0;

    for (unsigned long i = 0; i < correlationArrayLength; i++){
        if (correlationArray[i] > maxCorrelation ) {
            maxCorrelation =  correlationArray[i];
        }
    }

    int status = 0;

    double absThreshold = threshold * maxCorrelation;

    unsigned long num = 0;

    unsigned long long indexMax = 0;
    double maxValue = 0;

    for (unsigned long long i = 0; i < correlationArrayLength; i++){

        if (correlationArray[i] > absThreshold ) {
            
            if (status == 0){
                status = 1;
                region[(num << 1)] = i;
            }

            if (correlationArray[i] > maxValue){
                maxValue = correlationArray[i];
                indexMax = i;
            }
            
        }
        else{

            if (status == 1){                
                status = 0;
                region[(num << 1) + 1] = i;
                num ++;

                (*IndexPeak) = indexMax;
                indexMax = 0;
                maxValue = 0;
            }
        }

        if (num > numCandidateLimit){
            break;
        }
    }

    if (status == 1){                
        status = 0;
        region[num + 1] = correlationArray[correlationArrayLength - 1];
        num ++;
    }   

    (*numCandidate) = num;

    return 1;  
}

int frequencyOffsetEstimate(double *frequenceOffset,
                            float *data,
                            unsigned long long IndexPeak,
                            unsigned long numCandidate,
                            unsigned long samplingRateRatio,
                            unsigned long numerology)
{

    int ret = 0;

    const unsigned long FFTsizeExp = 12;
    const unsigned long numSamplePerSymbol = 4096;
    const unsigned long numSamplePerCP = 144 * 2;
    unsigned long numSamplePerCP0 = numSamplePerCP + (16 << (1 + numerology));

    unsigned long numSamplePerSymbolRev = numSamplePerSymbol * samplingRateRatio << 1;
    unsigned long numSamplePerCPrev = numSamplePerCP * samplingRateRatio << 1;
    unsigned long numSamplePerCP0Rev = numSamplePerCP0 * samplingRateRatio << 1;

    unsigned long long startIndex = (IndexPeak * samplingRateRatio << 1) - numSamplePerCP0Rev;

    if (startIndex >= 0){

        double sumProductRealCross =  0;
        double sumProductImageCross = 0;

        for (unsigned long long i = startIndex; i <numSamplePerCP0Rev; i = i + 2){

            sumProductRealCross +=  + data[i] * data[i + numSamplePerSymbolRev] 
                                    + data[(i + 1)]*data[(i + 1) + numSamplePerSymbolRev];

            sumProductImageCross += - data[i] * data[(i + 1) + numSamplePerSymbolRev] 
                                    + data[(i + 1)] * data[i + numSamplePerSymbolRev];

        }

        (*frequenceOffset) =         sumProductImageCross/sumProductRealCross;

        ret = 1;
    }
    else{
        ret = -1;
    }

    return ret;
}



