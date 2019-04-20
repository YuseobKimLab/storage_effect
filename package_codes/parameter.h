#ifndef __PARAMETER_H_
#define __PARAMETER_H_

#include "option.h"
#include <stdio.h>
#include <stdlib.h>

extern FILE *g_outfile;

typedef struct {
	int nDeme;
	int nTry;
	int nGen;
	int burnIn;

	int nInd1;
	int nInd2;
	int nIndT;
	int* popsize;

	int nSeq_block;
	int nSeq;
	int nNon;
	int nSyn;

	int mig;
	int period;

	double Nu;
	double recomRate;

	double optPheno;
	double epsilon_half;
	double sigma_s;
	double sigma_m;

	long seed;

}para_t;




para_t ReadParameter(const char* inpname);
para_t PrintParameter_screen(para_t p);
para_t PrintParameter_outfile(para_t p, FILE* outfile);


#endif
