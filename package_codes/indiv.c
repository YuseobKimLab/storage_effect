#include "indiv.h"




struct Individual New_Individual(struct Individual ind, int nSeq_block, double wild_pheno)
{

    ind.seq = (unsigned int*) malloc(nSeq_block * sizeof(unsigned int));
    memset(ind.seq,  0, sizeof(unsigned int) * nSeq_block);

    ind.pheno = wild_pheno;

    return ind;
}



struct Individual CopyIndividual(struct Individual dest, struct Individual sourc, int nSeq_block)
{
	memcpy(dest.seq, sourc.seq, (sizeof(int)) * nSeq_block);
	//indiv0 : destination / indiv1 : source

    dest.pheno = sourc.pheno;
	return dest;

}


void UpdatePhenotype (struct Individual* ind, unsigned int* nonsyn, double* delta, double pheno_wild, int nSeq_block)
{
    int i, count, delta_idx=0;
	double z = pheno_wild;
	unsigned int seg, template_seg;

    unsigned int* seq = ind->seq;
                                                                    #if PRINT_MUTATION
                                                                    printf("%lf ", z);
                                                                    #endif
	for (i = 0; i < nSeq_block; i++)
	{
		count = 0;
		seg = seq[i];
		template_seg = nonsyn[i];


		while ( count < _block_)
		{
			if ((template_seg & (1 << count))?1:0)			//if this site is nonsyn site
			{
				if ((seg & (1 << count))? 1: 0)				//if the allele is mutant type
				{
					z+= delta[delta_idx];
                                                                    #if PRINT_MUTATION
                                                                    printf("%lf ", delta[delta_idx]);
                                                                    #endif

				}
																	#if PRINT_MUTATION
																	else{printf("0(%lf) ", delta[delta_idx]);}
																	#endif
				delta_idx++;

			}

			count++;

		}


	}

	ind->pheno = z;
                                                                        #if PRINT_MUTATION
                                                                        printf("  => %lf \n", ind->pheno);
                                                                        #endif
    return;
}
