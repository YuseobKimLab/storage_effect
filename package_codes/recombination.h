#ifndef __RECOMBINATION_H_
#define __RECOMBINATION_H_

#include "random.h"
#include "indiv.h"
#include "printBit.h"

unsigned int* Get_Recombination_Template(int nRecom, int nSeq, int nSeq_block);
struct Individual Recombination(struct Individual* recombinant, struct Individual* ind1, struct Individual* ind2, double recomRate, int nSeq, int nSeq_block,  unsigned int* nonsyn, double* delta, double pheno_wild);




#endif
