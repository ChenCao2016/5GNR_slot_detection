
#ifndef __DMRSGENERATION_H
#define __DMRSGENERATION_H

unsigned long c_init(unsigned long slotNumInFrame, unsigned long symbolNumInSlot);

void GoldSequence(unsigned long *sequence, unsigned long length, unsigned long c_init);

void Rsequence(float* r, unsigned long length, unsigned long slotNumInFrame, unsigned long symbolNumInSlot);

void REvalue(float* resourceElement, unsigned long length, unsigned long slotNumInFrame, unsigned long symbolNumInSlot);

void fftShift(float* out, float* in, unsigned long rbNum);

#endif

