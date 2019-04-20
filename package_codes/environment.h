#ifndef __ENVIRONMENT_H
#define __ENVIRONMENT_H

#include <math.h>
#include <string.h>
#include "random.h"
#include "printBit.h"
#include "option.h"

unsigned int* Alloc_NonSynSite(int nSeq_block);
double* Alloc_Delta(int nNon);
double** Alloc_SelOpt(int nDeme, int period);

void Assign_NonsynSite_Random(unsigned int* nonsynSite, int nSeq_block, int nNon, int nSyn);
void Assign_NonsynSite_nonRandom(unsigned int* nonsynSite, int nSeq_block, int nNon, int nSyn);


void Assign_Delta(double* delta, int nNon, double sigma_m);
void Print_Delta(double* delta, int nNon, int try, FILE* g_outfile);

void Assign_SelOpt(double** selOpt, int mode, int period, int nDeme, double selMag, double optPheno);
void Disturb_SelOpt(double** selOpt,double** selOpt_origianl, double epsilon_half, int period, int nDeme);

#endif
