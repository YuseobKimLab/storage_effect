#ifndef __POP_H_
#define __POP_H_

#include <stdlib.h>
#include <stdio.h>
#include "indiv.h"
#include "printBit.h"
#include "option.h"

struct Subpopulation {
	struct Individual* indiv;
};



typedef struct Subpopulation* Population;

Population Alloc_Pop(int nDeme, int popsize[], int nSeq_block);
void Init_Pop(Population pop, int nDeme, int popsize[], int nSeq_block, double pheno_wild);


void PrintPop2(Population pop, int nDeme, int popsize[], int nSeq_block);
void PrintPop3(Population pop, int nDeme, int popsize[], int nSeq_block);



///////////////////////ERASE HERE//////////////////
void PrintBinary(unsigned int input);
///////////////////////ERASE HERE//////////////////


#endif
