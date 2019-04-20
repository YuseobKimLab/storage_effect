#include "parameter.h"


para_t ReadParameter(const char* inpname)
{

	int i;
    FILE *infile ;

    para_t p;

	if ( (infile = fopen(inpname, "r")) == NULL)
    {
        printf("Cannot open input file\n");
        exit(-1);           //* <stlib.h> should be included. close every file when there is an error
    }
	//read parameters from outFile

	fscanf(infile, "%*s%d", &(p.nDeme));
	fscanf(infile, "%*s%d", &(p.nTry));
	fscanf(infile, "%*s%d", &(p.nGen));

	/*fscanf(infile, "%*s%d", &(p.nInd1));
	fscanf(infile, "%*s%d", &(p.nInd2));
	p.nIndT = p.nInd1 + p.nInd2;
	*/

	p.nIndT = 0;
	p.popsize = (int*) malloc (sizeof(int) * p.nDeme);
	for (i =0 ; i < p.nDeme; i++)
	{
		fscanf(infile, "%*s%d", &(p.popsize[i]));
		//vi//printf("> %d : %d ", i, p.popsize[i]);
		p.nIndT += p.popsize[i];
	}



	fscanf(infile, "%*s%d", &(p.nSeq_block));
	fscanf(infile, "%*s%d", &(p.nNon));
	p.nSeq = p.nSeq_block * _block_;
	p.nSyn = p.nSeq - p.nNon;



	fscanf(infile, "%*s%d", &(p.mig));
	fscanf(infile, "%*s%d", &(p.period));

	fscanf(infile, "%*s%lf", &(p.Nu));
	fscanf(infile, "%*s%lf", &(p.recomRate));



	fscanf(infile, "%*s%lf", &(p.optPheno));
	fscanf(infile, "%*s%lf", &(p.epsilon_half));
	fscanf(infile, "%*s%lf", &(p.sigma_s));
	fscanf(infile, "%*s%lf", &(p.sigma_m));

	fscanf(infile, "%*s%ld", &(p.seed));
	fscanf(infile, "%*s%d", &(p.burnIn));

    return p;
}



para_t PrintParameter_screen(para_t p)
{
    int i;

    #if defined(_SBC_)
	printf( " << SEEDBANK - CONTROL    MODE >> \n");
	#elif defined(_SBE_)
	printf( " << SEEDBANK - EXPERIMNET MODE >> \n");
	#elif defined(_SBN_)
	printf( " << SEEDBANK - NEUTRAL    MODE >> \n");
	#elif defined(_RSC_)
	printf( " << REFUGE   - CONTROL    MODE >> \n");
	#elif defined(_RSE_)
	printf( " << REFUGE   - EXPERIMNET MODE >> \n");
	#elif defined(_RSN_)
	printf( " << REFUGE   - NEUTRAL    MODE >> \n");
	#else
	printf("++++++++++++++++++++++++++++++\n+[Warning] Please compile with option+\n++++++++++++++++++++++++++++++\n");
	#endif

	#if defined(_EXP_)
	printf( " >>  Mutation Exponential << \n");
	#endif

	printf( "nDeme..... %d\n", p.nDeme);
	printf( "nTry...... %d\n", p.nTry);
	printf( "nGen...... %d\n", p.nGen);

	printf( "nInd1_2_T. ");
	for (i =0 ; i < p.nDeme; i++)
	{
		printf("%d ", p.popsize[i]);
	}
	printf("%d \n", p.nIndT);


	printf( "N_S_Total. %d %d %d\n", p.nNon, p.nSyn, p.nSeq);

	printf( "mig....... %d\n", p.mig);
	printf( "period.... %d\n", p.period);


	printf( "Nu........ %lf\n", p.Nu);
	printf( "recomRate. %lf \n", p.recomRate);

	printf( "optPheno.. %lf\n", p.optPheno);
	printf( "epsilon_half... %lf\n", p.epsilon_half);
	printf( "sigma_s... %lf\n", p.sigma_s);
	printf( "sigma_m... %lf\n", p.sigma_m);

	printf( "seed...... %ld\n", p.seed);
	printf( "burnIn.... %d\n", p.burnIn);

	printf("==================================\n\n");
}


para_t PrintParameter_outfile(para_t p, FILE* outfile)
{
    int i;

	#if defined(_SBC_)
	fprintf(outfile, " << SEEDBANK - CONTROL    MODE >> \n");
	#elif defined(_SBE_)
	fprintf(outfile, " << SEEDBANK - EXPERIMNET MODE >> \n");
	#elif defined(_SBN_)
	fprintf(outfile, " << SEEDBANK - NEUTRAL    MODE >> \n");
	#elif defined(_RSC_)
	fprintf(outfile, " << REFUGE   - CONTROL    MODE >> \n");
	#elif defined(_RSE_)
	fprintf(outfile, " << REFUGE   - EXPERIMNET MODE >> \n");
	#elif defined(_RSN_)
	fprintf(outfile, " << REFUGE   - NEUTRAL    MODE >> \n");
	#else
	fprintf(outfile, "++++++++++++++++++++++++++++++\n+[Warning] Please compile with option+\n++++++++++++++++++++++++++++++\n");
	#endif


	#if defined(_EXP_)
	fprintf(outfile,  " >>  Mutation Exponential << \n");
	#endif

	fprintf(outfile, "nDeme %d\n", p.nDeme);
	fprintf(outfile, "nTry %d\n", p.nTry);
	fprintf(outfile, "nGen %d\n", p.nGen);

	//fprintf(outfile, "nInd1_2_Total %d %d %d\n", p.nInd1,  p.nInd1, p.nIndT);
	fprintf(outfile, "nInd1_2_T. ");
	for (i =0 ; i < p.nDeme; i++)
	{
		fprintf(outfile,"%d ", p.popsize[i]);
	}
	fprintf(outfile,"%d \n", p.nIndT);

	fprintf(outfile, "N_S_Total %d %d %d\n", p.nNon, p.nSyn, p.nSeq);

	fprintf(outfile, "mig %d\n", p.mig);
	fprintf(outfile, "period %d\n", p.period);


	fprintf(outfile, "Nu %lf\n", p.Nu);
	fprintf(outfile, "recomRate %lf\n", p.recomRate);


	fprintf(outfile, "optPheno %lf\n", p.optPheno);
	fprintf(outfile, "epsilon_half %lf\n", p.epsilon_half);
	fprintf(outfile, "sigma_s %lf\n", p.sigma_s);
	fprintf(outfile, "sigma_m. %lf\n", p.sigma_m);

	fprintf(outfile, "seed %ld\n", p.seed);
	fprintf(outfile,"burnIn.... %d\n", p.burnIn);

	fprintf(outfile,"==================================\n\n");


	return p;
}
