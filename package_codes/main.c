#include <stdio.h>
#include <time.h>
#include "mutation.h"
#include "migration.h"
#include "parameter.h"
#include "random.h"
#include "pop.h"
#include "environment.h"
#include "phenoList.h"
#include "recombination.h"
#include "subpopReproduction.h"
#include "reproduction.h"
#include "writeFile.h"




long gseed;

void main(int argc, char* argv[])
{

    const char* inpname;
    const char* out_head;
    FILE *outfile, *freqfile, *freqfile1, *freqfile2, *phenofile, *samplefile;

    int try , gen;
	int  sampling_if_negative=0;
    int nGen, burnIn;
    double pheno_wild;

    Population pop, pop1, pop2;

    unsigned int* nonsynSite;
	double* delta;
	double** selOpt;
	double** selOpt_original;



    /*time.h*/
    clock_t before;
    double clockResult;
    double estimated_remaining_time;

    before = clock();

    //Get Input File Name and Output Files name
    if(argc==3){
        inpname = argv[1];
		out_head = argv[2];
    }
	else{
		printf("./PROGRAM_NAME inpname outname");
	}

    outfile = Open_File_Write(outfile, out_head, -1, "out");

    //------------------------------

    //Init_Simulation
	para_t parameter = ReadParameter(inpname);
    PrintParameter_screen(parameter);
    PrintParameter_outfile(parameter, outfile);

	gseed = (-1) * parameter.seed;
	nGen = parameter.nGen;
    burnIn = parameter.burnIn;

    nonsynSite = Alloc_NonSynSite(parameter.nSeq_block);
	delta = Alloc_Delta(parameter.nNon);
	selOpt = Alloc_SelOpt(parameter.nDeme, parameter.period);
	selOpt_original = Alloc_SelOpt(parameter.nDeme, parameter.period);

    pop1 = Alloc_Pop(parameter.nDeme, parameter.popsize, parameter.nSeq_block);
    pop2 = Alloc_Pop(parameter.nDeme, parameter.popsize, parameter.nSeq_block);
    pop = pop1;

    listP* phenoList = Create_PhenoList();




    //Start New Trial
    try = 0 ;
    while ( try < parameter.nTry )
    {
        if (try > 0)
        {
            clockResult = (double)(clock()-before)/CLOCKS_PER_SEC;
            estimated_remaining_time = clockResult * ((parameter.nTry-try) / try);
            estimated_remaining_time /= 60;

            printf("> Try %3d | %.2lf min(s) | Estimated remaining time : %.2lf min(s) \n", try, clockResult/(60), estimated_remaining_time);
        }
        else
        {
            printf("> Try %3d | ------------------------------------------- \n", try);
        }


        Assign_SelOpt(selOpt_original, 1, parameter.period, parameter.nDeme,  1, parameter.optPheno);
        Assign_NonsynSite_nonRandom(nonsynSite , parameter.nSeq_block, parameter.nNon, parameter.nSyn);
        Assign_Delta( delta, parameter.nNon, parameter.sigma_m);

        pheno_wild = ran1(&gseed) * (parameter.optPheno * 2) - parameter.optPheno;
                //printf("phenowild = %lf\n", pheno_wild);
        Init_Pop(pop, parameter.nDeme, parameter.popsize, parameter.nSeq_block, pheno_wild);
//PrintPop3(pop1, parameter.nDeme, parameter.popsize, parameter.nSeq_block);
                //printf("=========================\n");



        freqfile = Open_File_Write(freqfile, out_head, try, "freq");
        freqfile1 = Open_File_Write(freqfile1, out_head, try, "freq1");
        freqfile2 = Open_File_Write(freqfile2, out_head, try, "freq2");
        #ifndef _NOPHENO_
        phenofile = Open_File_Write(outfile, out_head, try, "pheno");
        #endif
        samplefile = Open_File_Write(outfile, out_head, try, "sample");

        Print_NSinfo(nonsynSite, parameter.nSeq_block, freqfile, freqfile1, freqfile2);

        //Start New Generation
        gen = 0;

        while ( gen <  nGen + burnIn)
        {
            Init_PhenoList(phenoList);

            if ((gen + burnIn) % (_PRINTOUT_EVERY_) == 0)
            {
                clockResult = (double)(clock()-before)/CLOCKS_PER_SEC;
                printf("     (%3d)| Generation %15d / time : %lf min\n", try, gen-burnIn , clockResult/60);

            }


            if ( gen % parameter.period == 0)
            {
				Disturb_SelOpt(selOpt, selOpt_original, parameter.epsilon_half, parameter.period, parameter.nDeme);
            }

            Mutation(pop, parameter, nonsynSite, delta, pheno_wild);
                    //PrintPop3(pop1, parameter.nDeme, parameter.popsize, parameter.nSeq_block);
                    //printf("=========================\n");


            Migration(pop, parameter);
                    //PrintPop3(pop1, parameter.nDeme, parameter.popsize, parameter.nSeq_block);
                    //printf("=========================\n");

            #if defined(_SELC_)
                pop = Reproduction_Control(pop, pop1, pop2, nonsynSite, delta, pheno_wild, selOpt, gen, phenoList, parameter);
            #elif defined(_SELSB_)
                pop = Reproduction_Seedbank(pop, pop1, pop2, nonsynSite, delta, pheno_wild, selOpt, gen, phenoList, parameter);
            #elif defined(_SELRS_)
                pop = Reproduction_Refuge(pop, pop1, pop2, nonsynSite, delta, pheno_wild, selOpt, gen, phenoList, parameter);
            #elif defined(_SELN_)
                pop = Reproduction_Neutral(pop, pop1, pop2, nonsynSite, delta, pheno_wild, selOpt, gen, phenoList, parameter);
            #else
                printf("++++++++++++++++++++++++++++++++++++++++++++++++\n"); printf("+     [Warning] Please compile with option     +\n"); printf("++++++++++++++++++++++++++++++++++++++++++++++++\n");
            #endif

               // printf("---------pop  -----------\n");
               // PrintPop3(pop, parameter.nDeme, parameter.popsize, parameter.nSeq_block);
                //printf("=========================\n");

            if ( gen > burnIn )  //after burnIn
            {
                Write_Freq_File(pop, gen - burnIn, parameter, freqfile1, freqfile2, freqfile);
                #ifndef _NOPHENO_
                Print_Pheno_Outfile(phenoList,  gen - burnIn,selOpt[0][gen % parameter.period], phenofile);		//#phenoHere
                #endif
                //if (gen % _sampling_generation_ == 0)
                if (sampling_if_negative < 0)
                {
                        Write_Sampled_Sequence(pop, gen - burnIn,selOpt[0][gen % parameter.period], parameter,samplefile );
                }

            }

            gen++;
		        sampling_if_negative++;
            if ((gen + parameter.period) % _sampling_generation_ == 0)
            {
	            sampling_if_negative = - parameter.period;
            }


        }




        Print_Delta(delta, parameter.nNon, try,  outfile);
		try++;

        fflush(outfile);
        fclose(freqfile);
        fclose(freqfile1);
        fclose(freqfile2);
        #ifndef _NOPHENO_
        fclose(phenofile);
        #endif

    }

    /*time.h*/
	clockResult = (double)(clock()-before)/CLOCKS_PER_SEC;
	printf("\n totalTime = %lf hour(s)\n", clockResult / (60 * 60));
	fclose(outfile);
	free(parameter.popsize);



}
