#include "pop.h"





void PrintPop2(Population pop, int nDeme, int popsize[], int nSeq_block)
{
	int dm, i, s;
	//int count = 0;
	int seg_print;


	for (dm = 0; dm < nDeme; dm++)
	{

		printf(" > SubPopulation %d \n", dm);
		for (i = 0 ; i < popsize[dm] ; i++)
		{
			//sIdx_check = 0;
			printf("indiv %5d : ",  i);

			for (s = 0; s < nSeq_block; s++)
			{
				if(s==0)
					printf("%p ", &(pop[dm].indiv[i].seq[0]));

				seg_print = pop[dm].indiv[i].seq[s];
				PrintBinary(seg_print);
				printf("(%-15d)", seg_print);
				printf(" ");
				//count++;

			}
			printf("\n");
		}
	}
	printf("\n\n");

}


void PrintPop3(Population pop, int nDeme, int popsize[], int nSeq_block)
{
	int dm, i, s;
	//int count = 0;
	int seg_print;


	for (dm = 0; dm < nDeme; dm++)
	{

		printf(" > SubPopulation %d \n", dm);
		for (i = 0 ; i < popsize[dm] ; i++)
		{
			//sIdx_check = 0;
			printf("indiv %5d : ",  i);

			for (s = 0; s < nSeq_block; s++)
			{
				if(s==0)
					printf("%p(%6f)", &(pop[dm].indiv[i].seq[0]),  pop[dm].indiv[i].pheno);

				seg_print = pop[dm].indiv[i].seq[s];
				PrintBinary(seg_print);
				printf("(%-15d)", seg_print);
				printf(" ");
				//count++;

			}
			printf("\n");
		}
	}
	printf("\n\n");

}

Population Alloc_Pop(int nDeme, int popsize[], int nSeq_block)
{
	#if PRINT_WHERE
	printf("\n\n");
	printf("+-------------------------+\n");
	printf("|       Alloc_Pop         |\n");
	printf("+-------------------------+ \n\n");
	#endif


	int dm, i;

	Population pop = (struct Subpopulation*) malloc (sizeof(struct Subpopulation) * nDeme);
    for (dm=0; dm < nDeme; dm++)
	{
        pop[dm].indiv = (struct Individual*) malloc (sizeof(struct Individual) * popsize[dm]);

        for(i = 0; i < popsize[dm]; i++)
        {
            pop[dm].indiv[i].seq = (unsigned int*) malloc(nSeq_block * sizeof(unsigned int));
			memset(pop[dm].indiv[i].seq,  0, sizeof(unsigned int) * nSeq_block);
			pop[dm].indiv[i].pheno = 0;

        }
    }


    return pop;


}


void Init_Pop(Population pop, int nDeme, int popsize[], int nSeq_block, double pheno_wild)
{
	#if PRINT_WHERE
	printf("\n\n");
	printf("+-------------------------+\n");
	printf("|        InitPop          |\n");
	printf("+-------------------------+ \n\n");
	#endif

	int dm, i;

    for (dm=0; dm < nDeme; dm++)
	{
        pop[dm].indiv = (struct Individual*) malloc (sizeof(struct Individual) * popsize[dm]);

        for(i = 0; i < popsize[dm]; i++)
        {
            pop[dm].indiv[i].seq = (unsigned int*) malloc(nSeq_block * sizeof(unsigned int));
			memset(pop[dm].indiv[i].seq,  0, sizeof(unsigned int) * nSeq_block);
			pop[dm].indiv[i].pheno = pheno_wild;

        }
    }


   // return pop;
}
