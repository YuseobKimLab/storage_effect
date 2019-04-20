#include "subpopReproduction.h"

double FitOfPheno(double pheno, double sigma_s, double selOpt)
{
    double fit;

    fit = exp(-(pheno - selOpt)*(pheno - selOpt) / (2.0 * sigma_s * sigma_s));
                                                                            	#if PRINT_FIT
                                                                            	printf("z = %lf / selOpt = %lf -> exp(%lf) = %lf ", pheno, selOpt, (-(pheno - selOpt)*(pheno - selOpt) / (2.0 * sigma_s * sigma_s)), fit);
                                                                            	#endif

	return fit;


}



void subpop_Selection_Drift(struct Subpopulation parentSubpop, struct Subpopulation offspringSubpop, int popNum, int popSize, unsigned int* nonsyn, double* delta, double pheno_wild, double selOpt_gen,  listP* phenoTable,para_t p)
{
    #if PRINT_WHERE
    printf("\n\n");
    printf("+-------------------------+\n");
    printf("|   Selection - S O D O   |\n");
    printf("+-------------------------+ \n\n");
    #endif

    int count;
    int rInd1;
    int rInd2;

    double fit1;
    double pheno1;

    int count_ran = 0;

    double recomRate = p.recomRate;

    struct Individual gamate = New_Individual(gamate, p.nSeq_block, pheno_wild);

    while (count < popSize)
    {
        rInd1 = (int) (popSize * ran1(&gseed));
        rInd2 = (int) (popSize * ran1(&gseed));
        while (rInd2 == rInd1) { rInd2 = (int) (popSize * ran1(&gseed)); }
        count_ran++;

        if  ( recomRate == 0 )
           gamate = CopyIndividual(gamate, parentSubpop.indiv[rInd1],  p.nSeq_block);
        else
           gamate = Recombination(&gamate, &(parentSubpop.indiv[rInd1]), &(parentSubpop.indiv[rInd2]), recomRate, p.nSeq, p.nSeq_block, nonsyn, delta, pheno_wild);

        fit1 = FitOfPheno(gamate.pheno, p.sigma_s, selOpt_gen);

        if ( ran1(&gseed) < fit1 )
        {
                                                                                #if PRINT_SELECTION
                                                                                int b;  printf("   rInd1 : %d  / rInf2 : %d / ga pheno : %f / ga fit : %f \n", rInd1, rInd2, gamate.pheno, fit1);       printf("   ga ");for (b = 0; b < p.nSeq_block; b++){PrintBinary(gamate.seq[b]);  printf(" ");}   printf("\n");
                                                                                #endif

           offspringSubpop.indiv[count] = CopyIndividual(offspringSubpop.indiv[count], gamate, p.nSeq_block);
           #ifndef _NOPHENO_
           Count_Pheno(phenoTable, popNum, gamate.pheno, p.nSeq_block);
           #endif

           count += 1;
           if (count == popSize) { break; }

        }


    }

    free(gamate.seq);
    #if PRINT_SELECTION
    printf(" - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n");
    #endif

 //   fprintf(yofile, "pop%d - %d\n", popNum, count_ran);

}


void subpop_noSelection_Drift(struct Subpopulation parentSubpop, struct Subpopulation offspringSubpop, int popNum, int popSize, unsigned int* nonsyn, double* delta, double pheno_wild, double selOpt_gen,  listP* phenoTable,para_t p)
{
    #if PRINT_WHERE
    printf("\n\n");
    printf("+-------------------------+\n");
    printf("|   Selection - S X D O   |\n");
    printf("+-------------------------+ \n\n");
    #endif

    int count;
    int rInd1;
    int rInd2;

    double fit1;
    double pheno1;

    int count_ran = 0;

    double recomRate = p.recomRate;

    struct Individual gamate = New_Individual(gamate, p.nSeq_block, pheno_wild);

    while (count < popSize)
    {
        rInd1 = (int) (popSize * ran1(&gseed));
        rInd2 = (int) (popSize * ran1(&gseed));
        while (rInd2 == rInd1) { rInd2 = (int) (popSize * ran1(&gseed)); }
        count_ran++;

        if  ( recomRate == 0 )
           gamate = CopyIndividual(gamate, parentSubpop.indiv[rInd1],  p.nSeq_block);
        else
           gamate = Recombination(&gamate, &(parentSubpop.indiv[rInd1]), &(parentSubpop.indiv[rInd2]), recomRate, p.nSeq, p.nSeq_block, nonsyn, delta, pheno_wild);
                                                                                                                                           #if PRINT_SELECTION
                                                                                                                                           int b;  printf("   rInd1 : %d  / rInf2 : %d / ga pheno : %f / ga fit : %f \n", rInd1, rInd2, gamate.pheno, fit1);       printf("   ga ");for (b = 0; b < p.nSeq_block; b++){PrintBinary(gamate.seq[b]);  printf(" ");}   printf("\n");
                                                                                                                                           #endif
            offspringSubpop.indiv[count] = CopyIndividual(offspringSubpop.indiv[count], gamate, p.nSeq_block);
            #ifndef _NOPHENO_
            Count_Pheno(phenoTable, popNum, gamate.pheno, p.nSeq_block);
            #endif


            count += 1;
            if (count == popSize) { break; }
        }

    free(gamate.seq);
    #if PRINT_SELECTION
    printf(" - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n");
    #endif

}



void subpop_noSelection_noDrift(struct Subpopulation parentSubpop, struct Subpopulation offspringSubpop, int popNum, int popSize, unsigned int* nonsyn, double* delta, double pheno_wild, double selOpt_gen,  listP* phenoTable,para_t p)
{
    #if PRINT_WHERE
    printf("\n\n");
    printf("+-------------------------+\n");
    printf("|   Selection - S X D X   |\n");
    printf("+-------------------------+ \n\n");
    #endif

    int count, j;
    int randInd;
    int tmp;
    double pheno;

    int* ranArray;

    ranArray = (int*) malloc (popSize * sizeof(int) );		//for random sorting of second subpopulation
    for ( j = 0; j < popSize; j++)
    {
        ranArray[j] = j;
    }

    for (count = 0 ; count < popSize; count++)
    {
        //randomly arrange index number
        //randInd = ran1(&gseed) * (popSize - count) + count ;
        randInd = (int) (ran1(&gseed) * (popSize - count)) + count ;


        tmp = ranArray[count];
        ranArray[count] = ranArray[randInd];
        ranArray[randInd] = tmp;

        offspringSubpop.indiv[count] = CopyIndividual(offspringSubpop.indiv[count], parentSubpop.indiv[ranArray[count]], p.nSeq_block);
        #ifndef _NOPHENO_
        Count_Pheno(phenoTable, popNum, offspringSubpop.indiv[count].pheno, p.nSeq_block);
        #endif

    /*    if (count == popSize)
        {
            break;
        }
        */
    }

                                                                                #if PRINT_SELECTION
                                                                                int b;
                                                                                /*printf("   rInd : %d  / fit : %f / pointer %p / seq 0 pointer %p ", rInd1, fit1, &(parentSubpop.indiv[rInd1]), &(parentSubpop.indiv[rInd1].seq[0]));*/
                                                                                printf("random order :");
                                                                                for (b = 0; b < popSize; b++)
                                                                                {
                                                                                    printf(" %d", ranArray[b]);
                                                                                }
                                                                                printf("\n");

                                                                                for (count = 0 ; count < popSize; count++)
                                                                                {

                                                                                printf("ga%2d : ",  count);
                                                                                for (b = 0; b < p.nSeq_block; b++){PrintBinary(offspringSubpop.indiv[count].seq[b]);  printf(" ");printf("(%-15d)", offspringSubpop.indiv[count].seq[b]);}   printf("\n");
                                                                                }
                                                                                #endif

    free(ranArray);

}
