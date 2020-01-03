
#ifndef __SUBFRAMESELFCORRELATION_H
#define __SUBFRAMESELFCORRELATION_H

int subframeCorrelate(double *correlationArray, 
                        float * data, 
                        unsigned long long numSamples, 
                        unsigned long samplingRateRatio, 
                        unsigned long numerology,
                        int mode);

int subframeCorrelate2(double *correlationArray, 
                        float * data, 
                        unsigned long long numSamples, 
                        unsigned long samplingRateRatio, 
                        unsigned long numerology,
                        int mode);

int subframeCorrelate3(double *correlationArray,
                        double *phaseShiftArray,
                        float * data, 
                        unsigned long long numSamples, 
                        unsigned long samplingRateRatio, 
                        unsigned long numerology,
                        int mode);

int subframeStartRegionCandidate(unsigned long long *region,
                                    unsigned long long *IndexPeak, 
                                    unsigned long *numCandidate, 
                                    double *correlationArray,
                                    unsigned long long correlationArrayLength, 
                                    double threshold,
                                    unsigned long numCandidateLimit);

#endif



