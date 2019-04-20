#ifndef __MUTATION_H_
#define __MUTATION_H_


#include <stdio.h>
#include "pop.h"
#include "random.h"
#include "option.h"
#include "printBit.h"
#include "option.h"
#include "parameter.h"

//void Mutation(Population pop, para_t p);

void Mutation(Population pop, para_t p, unsigned int* nonsyn, double* delta, double pheno_wild);


#endif
