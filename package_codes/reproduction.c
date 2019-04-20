#include "reproduction.h"


Population Reproduction_Control(Population currentPop, Population pop1, Population pop2, unsigned int* nonsyn, double* delta, double pheno_wild, double** selOpt, int gen,listP* phenoTable,para_t p)
{
    Population parentPop;
    Population offspringPop;
    int i;
    double fit1, pheno1;
    int gen_period = gen % (p.period);


    if (currentPop == pop1)								//pop1 = parentPop, pop2 = offspringPop
	{
		parentPop = pop1;
		offspringPop = pop2;
	}
	else if (currentPop == pop2)						//pop2 = parentPop, pop1 = offspringPop
	{
		parentPop = pop2;
		offspringPop = pop1;
	}
                                                                                #if PRINT_SELECTION
                                                                                printf("-------parent pop-------\n");
                                                                                PrintPop3(parentPop, p.nDeme, p.popsize, p.nSeq_block);
                                                                                printf("-------off pop-------\n");
                                                                                PrintPop3(offspringPop, p.nDeme, p.popsize, p.nSeq_block);
                                                                                #endif


    i = 0;
    subpop_Selection_Drift(parentPop[i], offspringPop[i],i,  p.popsize[i], nonsyn, delta,  pheno_wild,  selOpt[i][gen_period],   phenoTable, p);
    i = 1;
    subpop_Selection_Drift(parentPop[i], offspringPop[i],i,  p.popsize[i], nonsyn, delta,  pheno_wild,  selOpt[i][gen_period],   phenoTable, p);

                                                                                #if PRINT_SELECTION
                                                                                printf("-------after reprod (parent )-------\n");
                                                                                PrintPop3(parentPop, p.nDeme, p.popsize, p.nSeq_block);
                                                                                printf("-------after reprod (offspring )-------\n");
                                                                                PrintPop3(offspringPop, p.nDeme, p.popsize, p.nSeq_block);
                                                                                #endif
    return offspringPop;
}



Population Reproduction_Refuge(Population currentPop, Population pop1, Population pop2,unsigned int* nonsyn, double* delta, double pheno_wild, double** selOpt, int gen,listP* phenoTable,para_t p)
{
    Population parentPop;
    Population offspringPop;
    int i;
    double fit1, pheno1;
    int gen_period = gen % (p.period);


    if (currentPop == pop1)								//pop1 = parentPop, pop2 = offspringPop
    {
        parentPop = pop1;
        offspringPop = pop2;
    }
    else if (currentPop == pop2)						//pop2 = parentPop, pop1 = offspringPop
    {
        parentPop = pop2;
        offspringPop = pop1;
    }
                                                                                #if PRINT_SELECTION
                                                                                printf("-------parent pop-------\n");
                                                                                PrintPop3(parentPop, p.nDeme, p.popsize, p.nSeq_block);
                                                                                printf("-------off pop-------\n");
                                                                                PrintPop3(offspringPop, p.nDeme, p.popsize, p.nSeq_block);
                                                                                #endif


    i = 0;
    subpop_Selection_Drift(parentPop[i], offspringPop[i],i,  p.popsize[i], nonsyn, delta,  pheno_wild,  selOpt[i][gen_period],  phenoTable, p);
    i = 1;
    subpop_noSelection_Drift(parentPop[i], offspringPop[i],i,  p.popsize[i], nonsyn, delta,  pheno_wild,  selOpt[i][gen_period],  phenoTable, p);

                                                                                #if PRINT_SELECTION
                                                                                printf("-------after reprod (parent )-------\n");
                                                                                PrintPop3(parentPop, p.nDeme, p.popsize, p.nSeq_block);
                                                                                printf("-------after reprod (offspring )-------\n");
                                                                                PrintPop3(offspringPop, p.nDeme, p.popsize, p.nSeq_block);
                                                                                #endif
    return offspringPop;
}


Population Reproduction_Seedbank(Population currentPop, Population pop1, Population pop2,unsigned int* nonsyn, double* delta, double pheno_wild, double** selOpt, int gen,listP* phenoTable,para_t p)
{
    Population parentPop;
    Population offspringPop;
    int i;
    double fit1, pheno1;
    int gen_period = gen % (p.period);


    if (currentPop == pop1)								//pop1 = parentPop, pop2 = offspringPop
    {
        parentPop = pop1;
        offspringPop = pop2;
    }
    else if (currentPop == pop2)						//pop2 = parentPop, pop1 = offspringPop
    {
        parentPop = pop2;
        offspringPop = pop1;
    }
                                                                                #if PRINT_SELECTION
                                                                                printf("-------parent pop-------\n");
                                                                                PrintPop3(parentPop, p.nDeme, p.popsize, p.nSeq_block);
                                                                                printf("-------off pop-------\n");
                                                                                PrintPop3(offspringPop, p.nDeme, p.popsize, p.nSeq_block);
                                                                                #endif


    i = 0;
    subpop_Selection_Drift(parentPop[i], offspringPop[i],i,  p.popsize[i], nonsyn, delta,  pheno_wild,  selOpt[i][gen_period],   phenoTable, p);
    i = 1;
    subpop_noSelection_noDrift(parentPop[i], offspringPop[i],i,  p.popsize[i], nonsyn, delta,  pheno_wild,  selOpt[i][gen_period],   phenoTable, p);

                                                                                #if PRINT_SELECTION
                                                                                printf("-------after reprod (parent )-------\n");
                                                                                PrintPop3(parentPop, p.nDeme, p.popsize, p.nSeq_block);
                                                                                printf("-------after reprod (offspring )-------\n");
                                                                                PrintPop3(offspringPop, p.nDeme, p.popsize, p.nSeq_block);
                                                                                #endif
    return offspringPop;
}



Population Reproduction_Neutral(Population currentPop, Population pop1, Population pop2, unsigned int* nonsyn, double* delta, double pheno_wild, double** selOpt, int gen,listP* phenoTable,para_t p)
{
    Population parentPop;
    Population offspringPop;
    int i;
    double fit1, pheno1;
    int gen_period = gen % (p.period);


    if (currentPop == pop1)								//pop1 = parentPop, pop2 = offspringPop
	{
		parentPop = pop1;
		offspringPop = pop2;
	}
	else if (currentPop == pop2)						//pop2 = parentPop, pop1 = offspringPop
	{
		parentPop = pop2;
		offspringPop = pop1;
	}
                                                                                #if PRINT_SELECTION
                                                                                printf("-------parent pop-------\n");
                                                                                PrintPop3(parentPop, p.nDeme, p.popsize, p.nSeq_block);
                                                                                printf("-------off pop-------\n");
                                                                                PrintPop3(offspringPop, p.nDeme, p.popsize, p.nSeq_block);
                                                                                #endif


    i = 0;
    subpop_noSelection_Drift(parentPop[i], offspringPop[i],i,  p.popsize[i], nonsyn, delta,  pheno_wild,  selOpt[i][gen_period],   phenoTable, p);
    i = 1;
    subpop_noSelection_Drift(parentPop[i], offspringPop[i],i,  p.popsize[i], nonsyn, delta,  pheno_wild,  selOpt[i][gen_period],   phenoTable, p);

                                                                                #if PRINT_SELECTION
                                                                                printf("-------after reprod (parent )-------\n");
                                                                                PrintPop3(parentPop, p.nDeme, p.popsize, p.nSeq_block);
                                                                                printf("-------after reprod (offspring )-------\n");
                                                                                PrintPop3(offspringPop, p.nDeme, p.popsize, p.nSeq_block);
                                                                                #endif
    return offspringPop;
}
