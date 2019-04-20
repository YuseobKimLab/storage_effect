#include "environment.h"

unsigned int* Alloc_NonSynSite(int nSeq_block)
{
    // To store the information about the site type (whether it is synonymous or nonsynonymous)
    unsigned int* nonsynSite;
    nonsynSite = (unsigned int*) malloc ( nSeq_block * sizeof(unsigned int));

    return nonsynSite;
}

void Assign_NonsynSite_Random(unsigned int* nonsynSite, int nSeq_block, int nNon, int nSyn)
{
	//make n sites to be nonsynonymous sites
	//mark nonsynonymous site as 1
	//return template

	int i, j, bit, count = 0;
	int mut_position;
	int mut_site, mut_block;
	unsigned int seg;
	int nSeq = nSeq_block * _block_;

	memset(nonsynSite, 0, sizeof(unsigned int) * nSeq_block);

	while (count < nNon)
	{
		mut_position = (int) (ran1(&gseed)*nSeq);

		mut_block = mut_position / _block_;
		mut_site = mut_position % _block_;

		seg = nonsynSite[mut_block];
		bit = (seg &  1 << (mut_site))?1:0;
																		#if PRINT_NONSYNSITE
																		printf("\nnonNum = %d (block %2d site %2d  (%d) ) >>%d<<\n" , count, mut_block, mut_site,mut_position, bit);
																		#endif
		if ( bit == 0 )
		{
			seg += ( 1 << (mut_site) );
			count++;
		}

		nonsynSite[mut_block] = seg;




	}

                                                                        #if PRINT_NONSYNSITE
                                                                        for (j = 0 ; j < nSeq_block; j++)
                                                                        {
                                                                            PrintBinary_asterik(nonsynSite[j]);
                                                                            printf(" ");
                                                                        }
                                                                        #endif

}

void Assign_NonsynSite_nonRandom(unsigned int* nonsynSite, int nSeq_block, int nNon, int nSyn)
{
  int i, j;
  int non_count = 0;
  memset(nonsynSite, 0, sizeof(unsigned int) * nSeq_block);


    if (nNon != nSyn)
    {
        //printf("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");         printf("+  [Warning] nNon != nSyn -> Please use Assign_NonsynSite_nonRandom instead  +\n"); printf("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
        //return ;
        for (i = 0; i < nSeq_block; i++)
        {
            for (j = 0; j < _block_/2; j++)
            {
              nonsynSite[i] = nonsynSite[i] << 2;
                if (non_count < nNon)
                {
                  nonsynSite[i] += 2;
                  non_count++;}
                else
                  {nonsynSite[i] += 0;}




            }


        }

    }

	//make n sites to be nonsynonymous sites
	//mark nonsynonymous site as 1
	//return template

    else
    {
    for (i = 0; i < nSeq_block; i++)
    {
        for (j = 0; j < _block_; j++)
        {
            nonsynSite[i] += 2;
            nonsynSite[i] = nonsynSite[i] << 2;
        }
        nonsynSite[i] += 2;

    }
                                                                                #if PRINT_NONSYNSITE
                                                                                for (j = 0 ; j < nSeq_block; j++)
                                                                                {
                                                                                    PrintBinary_asterik(nonsynSite[j]);
                                                                                    printf(" ");
                                                                                }
                                                                                #endif
    }

}


double* Alloc_Delta(int nNon)
{
    double* delta;
    delta = (double*) malloc (sizeof(double) * nNon);

    return delta;
}

void Assign_Delta(double* delta, int nNon, double sigma_m)
{
    int i;

	memset(delta, 0, sizeof(double)* nNon);
	for (i = 0; i< nNon; i++)
	{
        #ifndef _EXP_
        delta[i] = sigma_m * gasdev(&gseed);
        #else
        delta[i] = sigma_m * expdev(&gseed) * pow(-1, (int) (ran1(&gseed) * 100 ));
        #endif

	}
												#if PRINT_DELTA
												printf("delta : " );
												for (i = 0; i<nNon; i++){ printf("%lf ", delta[i]);}
												printf("\n");
												#endif
}

void Print_Delta(double* delta, int nNon, int try, FILE* outfile)
{
    int i;
    fprintf(outfile, "Try%d\t", try);
    for (i = 0; i < nNon; i++)
    {
        fprintf(outfile, "%lf ", delta[i]);
    }
    fprintf(outfile, "\n");

}

double** Alloc_SelOpt(int nDeme, int period)
{
    double** selOpt;
    int i;

    selOpt = (double**) malloc (sizeof(double*) * nDeme);
	for (i = 0; i < nDeme ;i++)
	{
		selOpt[i] = (double*) malloc (sizeof(double) * period) ;
	}

    return selOpt;

}

void Assign_SelOpt(double** selOpt, int mode, int period, int nDeme, double selMag, double optPheno)
{
    int i, j;
	int r = gseed % period;
	double s_t;
	int half_period = period / 2;


	for (i = 0; i < nDeme; i++)
	{
		memset(selOpt[i], 0, sizeof(double) * period);
	}

    if (mode == 1 )	// binary
	{

		for (i = 0; i < nDeme; i++)
		{
			for (j = 0; j < period; j++)
			{
				if (i ==0)
				{

					if ( j < half_period )
						selOpt[i][j] =  optPheno;
					else
						selOpt[i][j] = - optPheno;
				}
				else
				{
					if ( j < half_period )
						selOpt[i][j] = selMag * optPheno;
					else
						selOpt[i][j] = - (selMag * optPheno);
				}

			}
		}
	}
	else if (mode == 2) // sine fuction
	{
		for (i = 0; i < nDeme; i++)
		{
			for (j = 0; j < period; j++)
			{
				s_t =  optPheno * sin(2*PI*((double)(i+r)/period));
				if (i ==0)
				{
					selOpt[i][j] =  s_t;
				}
				else
				{
					selOpt[i][j] =  selMag  * s_t;
				}

			}

		}
	}
											#if PRINT_SELOPT
											for (j = 0; j < period; j++){    printf("%dth  : %lf\t%lf\n", j, selOpt[0][j], selOpt[1][j]);	}
											#endif

}

void Disturb_SelOpt(double** selOpt,double** selOpt_origianl, double epsilon_half, int period, int nDeme)
{
    int i,j;
	int half_period = period/2;
	double e1, e2;         //flucaution disturbance during season 1 and 2

	e1 = (ran1(&gseed) * 2 * epsilon_half) - epsilon_half;	//random from (-epsilon_half, epsilon_half
	e2 = (ran1(&gseed) * 2 * epsilon_half) - epsilon_half;

    for (i = 0; i < nDeme ; i++)
	{
		for (j = 0; j < period; j++)
		{

			if ( j < half_period )
				selOpt[i][j] = selOpt_origianl[i][j] + e1;
			else
				selOpt[i][j] = selOpt_origianl[i][j] + e2;


		}
	}

    #if PRINT_DIST_SELOPT
	printf("=====during this period, ===== \n" );
	printf(" ranE1 : %lf, ranE2 : %lf  \n", e1, e2);

	for (j = 0; j < period; j++)
	{ printf("%dth  : %lf\t%lf\n", j, selOpt[0][j], selOpt[1][j]); }
	printf("---------------------\n\n");
	#endif


}
