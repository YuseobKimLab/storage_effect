#ifndef __REPRODUCTION_H_
#define __REPRODUCTION_H_

#include "subpopReproduction.h"
#include "pop.h"



Population Reproduction_Control(Population currentPop, Population pop1, Population pop2, unsigned int* nonsyn, double* delta, double pheno_wild, double** selOpt, int gen,listP* phenoTable,para_t p);

Population Reproduction_Refuge(Population currentPop, Population pop1, Population pop2,unsigned int* nonsyn, double* delta, double pheno_wild, double** selOpt, int gen,listP* phenoTable,para_t p);

Population Reproduction_Seedbank(Population currentPop, Population pop1, Population pop2,unsigned int* nonsyn, double* delta, double pheno_wild, double** selOpt, int gen,listP* phenoTable,para_t p);

Population Reproduction_Neutral(Population currentPop, Population pop1, Population pop2, unsigned int* nonsyn, double* delta, double pheno_wild, double** selOpt, int gen,listP* phenoTable,para_t p);

#endif
