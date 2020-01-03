#include "DMRScorrelation.h"
#include <stdio.h>


int DMRScorrelate(double *correlationArray,
                        double *maxCorrelation,
                        unsigned long long *maxIndex,
                        double * DMRS,
                        float * data, 
                        unsigned long long numSamples, 
                        unsigned long samplingRateRatio,
                        int mode)
{

    unsigned long stepSize = 1;

    if (mode == 1){
        stepSize = samplingRateRatio;
    }


    double maxValue = 0;
    unsigned long long maxValueIndex = 0;
    unsigned long long n = 0;


    for (unsigned long long j = 0; j < (numSamples -(4096*samplingRateRatio << 1)); j += (stepSize << 1)){
        
        double sumProductRealCross = 0;
        double sumProductImageCross = 0;
        double sumProductRealSelf = 0;

        for (unsigned long i = 0; i < (4096 << 1); i += 2){

            sumProductRealCross +=  + DMRS[i] * data[j + i*samplingRateRatio] 
                                    + DMRS[i + 1]*data[j + i*samplingRateRatio + 1];

            sumProductImageCross += - DMRS[i] * data[j + i*samplingRateRatio + 1] 
                                    + DMRS[i + 1] * data[j + i*samplingRateRatio];

            sumProductRealSelf +=   + data[j + i*samplingRateRatio]  * data[j + i*samplingRateRatio] 
                                    + data[j + i*samplingRateRatio + 1]*data[j + i*samplingRateRatio + 1];

        }

        correlationArray[n] =   (sumProductRealCross*sumProductRealCross + sumProductImageCross*sumProductImageCross)
                                /(sumProductRealSelf * sumProductRealSelf);

        if (maxValue < correlationArray[n]){
            maxValue = correlationArray[n];
            maxValueIndex = n;
        }        

        
        n = n + 1;
    }
       
    (*maxCorrelation) = maxValue;
    (*maxIndex) = maxValueIndex; 

    return 1;
  
}


int DMRScorrelate2(double *correlationArray,
                        double *maxCorrelation,
                        unsigned long long *maxIndex,
                        double * DMRS,
                        float * data, 
                        unsigned long long numSamples, 
                        unsigned long samplingRateRatio,
                        int mode)
{

    unsigned long stepSize = 1;

    if (mode == 1){
        stepSize = samplingRateRatio;
    }


    double maxValue = 0;
    unsigned long long maxValueIndex = 0;
    unsigned long long n = 0;

    
    double correlationSum = 0;

    for (unsigned long long j = 0; j < (numSamples - (4096*samplingRateRatio << 1)); j += (stepSize << 1)){
        
        double sumProductRealCross = 0;
        double sumProductImageCross = 0;
        double sumProductRealSelf = 0;
      

        for (unsigned long i = 0; i < (4096 << 1); i += 2){

            sumProductRealCross +=  + DMRS[i] * data[j + i*samplingRateRatio] 
                                    + DMRS[i + 1]*data[j + i*samplingRateRatio + 1];

            sumProductImageCross += - DMRS[i] * data[j + i*samplingRateRatio + 1] 
                                    + DMRS[i + 1] * data[j + i*samplingRateRatio];

            sumProductRealSelf +=   + data[j + i*samplingRateRatio]  * data[j + i*samplingRateRatio] 
                                    + data[j + i*samplingRateRatio + 1]*data[j + i*samplingRateRatio + 1];

        }

        correlationArray[n] =   (sumProductRealCross*sumProductRealCross + sumProductImageCross*sumProductImageCross)
                                /(sumProductRealSelf * sumProductRealSelf);

        correlationSum = correlationSum + correlationArray[n];

        if (maxValue < correlationArray[n]){
            maxValue = correlationArray[n];
            maxValueIndex = n;
        }        
      
        n = n + 1;
    }
       
    (*maxCorrelation) = maxValue*n/correlationSum;
    (*maxIndex) = maxValueIndex; 

    return 1;
  
}

int DMRScandidateRegionShift(unsigned long long *regionOut, 
                                    unsigned long long *regionIn,
                                    unsigned long numCandidate, 
                                    unsigned long DMRSsymbolIndex)
{

    unsigned long symbolLength = 4384; //4096 + 144*2 normal cyclic prfix, l!= 0 or l!=7<<numberology
    unsigned long shift = symbolLength * DMRSsymbolIndex;

    for (unsigned long i = 0; i < (numCandidate<<1); i++){
        regionOut[i] = regionIn[i] + shift; 
    }    

    return 1;  
}



