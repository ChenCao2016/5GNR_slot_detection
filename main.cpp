#include <iostream>
#include <fstream>
#include <string>
#include <chrono>
#include <math.h>

extern "C" {
#include "subFrameSelfCorrelation.h"
#include "DMRSgeneration.h"
#include "DMRScorrelation.h"
}

#include "fftw3.h"

using namespace std;

int main(){

    /*----------------------------------------------------------------------------------------------------------------
    Load waveform file
    -----------------------------------------------------------------------------------------------------------------*/
    std::ofstream result;
    unsigned long long start, end;

    start = chrono::duration_cast<chrono::milliseconds>(chrono::system_clock::now().time_since_epoch()).count();
    //std::ifstream fin("qcom_mu3_bw100_66@0.dat", std::ios::binary);
    //std::ifstream fin("nr5g_ul_fs2457p6_mu2_bw100_mod_QPSK_tfpc0_rbo000_rbd132_cc01_sf01.dat", std::ios::binary);
    //std::ifstream fin("0.dat", std::ios::binary);
    std::ifstream fin("mu3_bw100_qpsk_tfpc0_rbo000_rbd066.dat", std::ios::binary);

    unsigned long slot = 0;
    unsigned long dmrsSymb = 3;
    unsigned long numerology = 3;
    unsigned long rbNum = 66;

    const long long FFTsizeExp = 12;
    unsigned long waveformSamplingRate = 2457600000;
    unsigned long subcarrierSpacing = (1 << numerology) * 15000;
    unsigned long expectedSamplingRate = subcarrierSpacing << FFTsizeExp;
    unsigned long samplingRateRatio = waveformSamplingRate/expectedSamplingRate;

    std::cout << "Sampling Rate Ratio: " << samplingRateRatio << "\n";

    fin.seekg(0, std::ios::end);
    unsigned long long numSamples = fin.tellg() / sizeof(float);
    std::cout << "Capture Length: " << (numSamples>>1)/2457.6 << "us\n";
    std::cout << "Number of samples: " << (numSamples>>1) << "\n";
    fin.seekg(0, std::ios::beg);

    float *data;
    data = new float[numSamples];

    float f;
    for (unsigned long long i = 0; i <numSamples; i++){
        fin.read(reinterpret_cast<char*>(&f), sizeof(float));
        data[i] = f;
    } 

    end = chrono::duration_cast<chrono::milliseconds>(chrono::system_clock::now().time_since_epoch()).count();
    std::cout << "Load data: " << (end - start) << "ms\n";  

    /*----------------------------------------------------------------------------------------------------------------
    CP self correlation for subframe detection
    -----------------------------------------------------------------------------------------------------------------*/
    double *correlation;
    double *phaseShiftArray;
    unsigned long long correlationArrayLength = (numSamples>>1)/samplingRateRatio;
    std::cout << "Correlation Array Length: " << correlationArrayLength << "\n";
    correlation = new double[correlationArrayLength]; //CP correlation array
    phaseShiftArray = new double[correlationArrayLength];


    start = chrono::duration_cast<chrono::milliseconds>(chrono::system_clock::now().time_since_epoch()).count();
    //subframeCorrelate(correlation,data,numSamples,samplingRateRatio,numerology, 0);
    //subframeCorrelate2(correlation,data,numSamples,samplingRateRatio,numerology, 0);
    subframeCorrelate3(correlation,phaseShiftArray,data,numSamples,samplingRateRatio,numerology, 0);
    end = chrono::duration_cast<chrono::milliseconds>(chrono::system_clock::now().time_since_epoch()).count();
    std::cout << "Self Correlation: " << (end - start) << "ms\n";  

    unsigned long long *region;
    unsigned long long *IndexPeak;
    unsigned long numCandidateLimit = 10;
    region = new unsigned long long [numCandidateLimit << 1];
    IndexPeak = new unsigned long long [numCandidateLimit];
    unsigned long numCandidate = 0;

    start = chrono::duration_cast<chrono::milliseconds>(chrono::system_clock::now().time_since_epoch()).count();
    subframeStartRegionCandidate(region, IndexPeak, &numCandidate, correlation, correlationArrayLength, 0.95, numCandidateLimit);    
    end = chrono::duration_cast<chrono::milliseconds>(chrono::system_clock::now().time_since_epoch()).count();
    std::cout << "Find subframe candidate: " << (end - start) << "ms\n";
    std::cout << "Number of subframe candidate region: " << numCandidate << "\n";
    for (unsigned long i = 0; i < (numCandidate<<1); i += 2){
        std::cout << "["<< region[i] <<"," << region[i + 1] <<"]\n"; 
    }

    double max = 0;   
    for (unsigned long long i = 0; i < correlationArrayLength; i++){
        if (correlation[i] > max ) {
            max =  correlation[i];
        }
    }
    std::cout << "Max: " << max << "\n";

    /*----------------------------------------------------------------------------------------------------------------
    DMRS cross correlation
    -----------------------------------------------------------------------------------------------------------------*/
    unsigned long long *regionOut;
    regionOut = new unsigned long long [numCandidate << 1];
    DMRScandidateRegionShift(regionOut, region, numCandidate, dmrsSymb);

    std::cout << "DMRS symbol candidate\n"; 
    for (unsigned long i = 0; i < (numCandidate<<1); i += 2){
        std::cout << "["<< regionOut[i] <<"," << regionOut[i + 1] <<"]\n"; 
    }


    unsigned long candidate = 0;
    for (slot = 0; slot < 40; slot += 2){
        start = chrono::duration_cast<chrono::milliseconds>(chrono::system_clock::now().time_since_epoch()).count();
        float *re;
        float *reRearrange;
        re = new float[(12 * rbNum) << 1];
        reRearrange = new float[4096 << 1];
        REvalue(re, 12 * rbNum, slot, dmrsSymb);    
        fftShift(reRearrange, re, rbNum);
        end = chrono::duration_cast<chrono::milliseconds>(chrono::system_clock::now().time_since_epoch()).count();
        std::cout << "Generate DMRS sequence: " << (end - start) << "ms\n";    
        
        result.open ("REvalue1.csv");

        for (unsigned long i = 0; i < (12 * rbNum) << 1; i++){
            result  << re[i] << "\n";
        }
        result.close();

        result.open ("reRearrange1.csv");

        for (unsigned long i = 0; i < (4096 << 1); i++){
            result  << reRearrange[i] << "\n";
        }
        result.close();


        delete [] re;


        start = chrono::duration_cast<chrono::milliseconds>(chrono::system_clock::now().time_since_epoch()).count();
        fftw_complex *in, *out;
        fftw_plan p;
        int N = 4096;
        in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N);
        out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N);
        for (unsigned long int i = 0; i < N; i++){
            in[i][0] = reRearrange[i << 1];
            in[i][1] = reRearrange[(i << 1) + 1];
        }
        delete [] reRearrange;
        p = fftw_plan_dft_1d(N, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
        fftw_execute(p);
        fftw_destroy_plan(p);
        end = chrono::duration_cast<chrono::milliseconds>(chrono::system_clock::now().time_since_epoch()).count();
        std::cout << "IDFT: " << (end - start) << "ms\n";

        start = chrono::duration_cast<chrono::milliseconds>(chrono::system_clock::now().time_since_epoch()).count();
        
        double *DMRS;
        DMRS = new double [4096 << 1];

        for (unsigned long i = 0; i < 4096; i++){
            DMRS[(i << 1)] = out[i][0];
            DMRS[(i << 1) + 1] = out[i][1];
        } 

        result.open ("DMRS1.csv");

        for (unsigned long i = 0; i < (4096 << 1); i++){
            result  << DMRS[i] << "\n";
        }
        result.close();


        unsigned long rangeLength = regionOut[candidate*2 + 1] - regionOut[candidate*2];

        double *correlationArray;
        correlationArray = new double [rangeLength];
        double maxCorrelation = 0;
        unsigned long long maxIndex = 0;

        DMRScorrelate(correlationArray,
                                &maxCorrelation,
                                &maxIndex,
                                DMRS,
                                &data[regionOut[candidate*2] * samplingRateRatio * 2],
                                unsigned long long ((rangeLength + 4096) * samplingRateRatio * 2), 
                                samplingRateRatio,
                                1);

        end = chrono::duration_cast<chrono::milliseconds>(chrono::system_clock::now().time_since_epoch()).count();
        std::cout << "DMRS cross correlation: " << (end - start) << "ms\n";
        std::cout << "SLOT: " << slot << "\n";
        std::cout << "DMRS cross correlation: " << maxCorrelation << "," << maxIndex << "\n";

        if (slot == 10){
            result.open ("DMRScorrelation.csv");

            for (unsigned long i = 0; i < rangeLength; i++){
                result  << correlationArray[i] << "\n";
            }
            result.close();
        }


        delete [] DMRS;
        delete [] correlationArray;

    }


    
    std::cout << "I am done\n"; 
    getchar();

}
