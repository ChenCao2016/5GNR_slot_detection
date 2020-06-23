#ifndef __DMRSCORRELATION_H
#define __DMRSCORRELATION_H


int DMRScorrelate(double *correlationArray,
                        double *maxCorrelation,
                        unsigned long long *maxIndex,
                        double * DMRS,
                        float * data, 
                        unsigned long long numSamples, 
                        unsigned long samplingRateRatio,
                        int mode);

int DMRScorrelate2(double *correlationArray,
                        double *maxCorrelation,
                        unsigned long long *maxIndex,
                        double * DMRS,
                        float * data, 
                        unsigned long long numSamples, 
                        unsigned long samplingRateRatio,
                        int mode);

int DMRScandidateRegionShift(unsigned long long *regionOut, 
                                    unsigned long long *regionIn,
                                    unsigned long numCandidate, 
                                    unsigned long DMRSsymbolIndex);

#endif

