#ifndef __SUBPOPREPRODUCTION_H
#define __SUBPOPREPRODUCTION_H


#include "pop.h"
#include "phenoList.h"
#include "recombination.h"
#include "parameter.h"
#include "random.h"



double FitOfPheno(double pheno, double sigma_s, double selOpt);

void subpop_Selection_Drift(struct Subpopulation parentSubpop, struct Subpopulation offspringSubpop, int popNum, int popSize, unsigned int* nonsyn, double* delta, double pheno_wild, double selOpt_gen,  listP* phenoTable,para_t p);

void subpop_noSelection_Drift(struct Subpopulation parentSubpop, struct Subpopulation offspringSubpop, int popNum, int popSize, unsigned int* nonsyn, double* delta, double pheno_wild, double selOpt_gen,  listP* phenoTable,para_t p);


void subpop_noSelection_noDrift(struct Subpopulation parentSubpop, struct Subpopulation offspringSubpop, int popNum, int popSize, unsigned int* nonsyn, double* delta, double pheno_wild, double selOpt_gen,  listP* phenoTable,para_t p);

#endif
