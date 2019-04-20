#ifndef __INDIV_H_
#define __INDIV_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "option.h"

struct Individual{
	unsigned int* seq;
	double pheno;

};

struct Individual New_Individual(struct Individual ind, int nSeq_block, double wild_pheno);
struct Individual CopyIndividual(struct Individual dest, struct Individual sourc, int nSeq_block);
void UpdatePhenotype (struct Individual* ind, unsigned int* nonsyn, double* delta, double pheno_wild, int nSeq_block);



#endif
