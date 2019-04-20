#include "mutation.h"


void Mutation(Population pop, para_t p, unsigned int* nonsyn, double* delta, double pheno_wild)
//1114 version
{
	//mutation occurs regardless of the syn / nonsyn
	#if PRINT_WHERE
	printf("\n\n");
	printf("+-------------------------+\n");
	printf("|         Mutation        |\n");
	printf("+-------------------------+ \n\n");
	#endif

	double Nu = p.Nu;
	int nSeq_block = p.nSeq_block;
	int nIndT = p.nIndT;
	int* popsize = p.popsize;

	int nMut;
	int mut_indiv, mut_position;
	unsigned int seg, binary=0;
	int s, m;
	struct Individual* ranIdv;



	nMut = poidev(Nu * _block_ * nSeq_block, &gseed);
		//printf("nMut : %d\n", nMut);


	for ( m = 0; m < nMut ; m++)
	{
		#ifdef _MUT12_
			//mutation occurs in both subpopulation
			mut_indiv = (int) (ran1(&gseed) * nIndT);
			if (mut_indiv < popsize[0])
				ranIdv = &(pop[0].indiv[mut_indiv]);
			else
				ranIdv = &(pop[1].indiv[mut_indiv-popsize[0]]);
		#elif defined(_MUT1_)
			//mutation occurs in only sub0
			mut_indiv = (int) (ran1(&gseed) * popsize[0]);
			ranIdv = &(pop[0].indiv[mut_indiv]);
		#else
			printf("++++++++++++++++++++++++++++++\n+[Warning] Please compile with option+\n++++++++++++++++++++++++++++++\n");
		#endif


		mut_position = (int) (ran1(&gseed) * nSeq_block * _block_ );

											#if PRINT_MUTATION
											if (mut_indiv < popsize[0]) {printf("sub0, ind %d" ,mut_indiv);}
											else {printf("sub1, ind %d" ,mut_indiv-popsize[0]);}
											printf(" 	block %d site %d =>  ", mut_position/_block_, mut_position%_block_);
											#endif


		seg = ranIdv->seq[mut_position/_block_];
		seg = seg ^ (1 << ((mut_position%_block_)));
		ranIdv->seq[mut_position/_block_] = seg;

		UpdatePhenotype(ranIdv, nonsyn, delta, pheno_wild, p.nSeq_block);
	}
}
