//TS 138 211 V15.3.0 P-58

#include "DMRSgeneration.h"

const unsigned long Nc = 1600;
const unsigned long N_Id_n_SCID = 0;
const unsigned long n_SCID = 0;

unsigned long c_init(unsigned long slotNumInFrame, unsigned long symbolNumInSlot){

    unsigned long res = (((14 * slotNumInFrame + symbolNumInSlot + 1)*((N_Id_n_SCID << 1) + 1) << 17) + (N_Id_n_SCID << 1) + n_SCID) % (1 << 31);

    return res;
}

void GoldSequence(unsigned long *sequence, unsigned long length, unsigned long init){

    unsigned long *x1, *x2;
    x1 = (unsigned long *) malloc(sizeof(unsigned long) * (31 + Nc + length));
    x2 = (unsigned long *) malloc(sizeof(unsigned long) * (31 + Nc + length));

    for (unsigned long i = 0; i<31; i++){

        x2[i] = (init >> i) & 0x1;

        if (i == 0){
            x1[i] = 1;
        }
        else{
            x1[i] = 0;
        }
    }

    for (unsigned long i = 0; i< Nc + length; i++){
        x1[i + 31] = (x1[i + 3] + x1[i]) & 0x1;
        x2[i + 31] = (x2[i + 3] + x2[i + 2] + +x2[i+1] + x2[i]) & 0x1;
    }

    for (unsigned long i = 0; i< length; i++){
        sequence[i] = (x1[i + Nc] + x2[i + Nc]) & 0x1;
    }    

    free(x1);
    free(x2);

} 

void Rsequence(float* r, unsigned long length, unsigned long slotNumInFrame, unsigned long symbolNumInSlot){

    unsigned long goldLength = (length << 1) + 1;

    unsigned long * goldSequence;
    goldSequence = (unsigned long*) malloc(sizeof(unsigned long) * goldLength);

    unsigned long init = c_init(slotNumInFrame, symbolNumInSlot);

    GoldSequence(goldSequence, goldLength, init);


    for (unsigned long i = 0; i <(length << 1); i+=2){

        r[i] = (1 - (signed int)(goldSequence[i] << 1)) / 1.414;
        r[i + 1] = (1 - (signed int)(goldSequence[i + 1] << 1)) / 1.414;

    }

    free(goldSequence);

}


//Configuration type 1, Antenna port 0
void REvalue(float* resourceElement, unsigned long length, unsigned long slotNumInFrame, unsigned long symbolNumInSlot){

    unsigned long int rSequenceLength = (length >> 1) + 1;

    float *r;

    r = (float *)malloc(sizeof(float) * (rSequenceLength << 1));

    Rsequence(r, rSequenceLength, slotNumInFrame, symbolNumInSlot);


    unsigned long k = 0;
    for (unsigned long i = 0; i < ((length - 2) >> 2) + 1; i++){
        

        resourceElement[k] = r[(i << 1) << 1];
        resourceElement[k + 1] = r[((i << 1) << 1) + 1];

        k = k + 2;

        resourceElement[k] = 0;
        resourceElement[k + 1] = 0;

        k = k + 2;

        resourceElement[k] = r[((i << 1) + 1) << 1];
        resourceElement[k + 1] = r[(((i << 1) + 1) << 1) + 1];

        k = k + 2;

        resourceElement[k] = 0;
        resourceElement[k + 1] = 0;     


        k = k + 2;  

    }

    free(r);

}


void fftShift(float* out, float* in, unsigned long rbNum){

    unsigned long  mcount = 0;
    for (unsigned long i = 0; i < ((6 * rbNum) << 1); i += 2 ){
        out[i] = in[(mcount + 6 * rbNum) << 1];
        out[i + 1] = in[((mcount + 6 * rbNum)<<1) + 1];
        mcount = mcount + 1;
    }
  
    for (unsigned long i = ((6 * rbNum) << 1); i < ((4096 - 6 * rbNum) <<1); i += 2){
        out[i] = 0;
        out[i + 1] = 0;
    } 

    mcount = 0;
    for (unsigned long i = ((4096 - 6 * rbNum) << 1); i < (4096 << 1); i += 2){
        out[i] = in[mcount << 1];
        out[i + 1] = in[(mcount << 1) + 1];
        mcount = mcount + 1;
    }           


}