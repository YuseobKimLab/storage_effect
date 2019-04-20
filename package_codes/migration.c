#include "migration.h"

void Migration(Population pop, para_t p)
{
	#if PRINT_WHERE
	printf("\n\n");
	printf("+-------------------------+\n");
	printf("|         Migration       |\n");
	printf("+-------------------------+ \n\n");
	#endif

	int nSeq_block = p.nSeq_block;
	int nDeme = p.nDeme;
	int* popsize = p.popsize;
	int mig = p.mig;

	struct Individual tmp;
	tmp.seq = (unsigned int*) malloc(sizeof(unsigned int) * nSeq_block);

	int i;

	if (nDeme == 2)
	{
		for (i=0; i < mig; i++)
		{
			tmp = CopyIndividual(tmp, pop[0].indiv[i], nSeq_block);
			pop[0].indiv[i] = CopyIndividual(pop[0].indiv[i], pop[1].indiv[i], nSeq_block);
			pop[1].indiv[i] = CopyIndividual(pop[1].indiv[i] ,tmp ,nSeq_block);
		}
	}
	free(tmp.seq);
}
