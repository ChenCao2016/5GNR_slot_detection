
#ifndef __SUBFRAMESELFCORRELATION_H
#define __SUBFRAMESELFCORRELATION_H

int subframeCorrelate(double *correlationArray, 
                        float * data, 
                        unsigned long long numSamples, 
                        unsigned long samplingRateRatio, 
                        unsigned long numerology,
                        int mode);

#endif



