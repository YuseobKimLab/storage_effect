#include "recombination.h"

unsigned int* Get_Recombination_Template(int nRecom, int nSeq, int nSeq_block)
{
    // making a recombination template (where is from parent 1 or 2; 0 from parent 1, 1 from parent 2)

    int i, n, b;
    int recomPosition, recomSite, recomBlock;

    unsigned int* template_recom;
    unsigned int template;

    template_recom = (unsigned int*) malloc (nSeq_block * sizeof(unsigned int));
	for (i = 0; i < nSeq_block; i++)
	{
		template_recom[i] = 0;
	}
                                                        #if PRINT_RECOM
                                                        nRecom = 2;
                                                        printf(" +----------------------------------------------------------------------------+\n");
                                                        #endif
   for (n = 0; n < nRecom; n++)
   {
       recomPosition = (int) (ran1(&gseed) * (nSeq-1)) + 1;    // if there is site 1,2,3
                                                               // recombination can happen between 1-2 and 2-3

       recomSite  =  recomPosition % _block_;
       recomBlock =  nSeq_block -1 - (recomPosition / _block_);
                                                                               #if PRINT_RECOM
                                                                               printf(" | recomPosition : %d (%d - %d)\n" , recomPosition, recomBlock, recomSite);
                                                                               #endif

       for (b = 0; b < nSeq_block; b++)
       {
           if (b > recomBlock)             // if recombination occurs at block 2 site 1,
                                           // all sites in block 3, 4, 5 ... will have parent 2's allele
           {
               template = ((1 << (_block_)) - 1);
               template_recom[b] = ( template_recom[b] ^ template );
           }
           else if (b == recomBlock)       // in block 2, where recombination occurs,
                                           // sites 2, 3, 4, ... will have parent 2's allele
           {
               template = ( 1 << recomSite )  -1;
               template_recom[b] = ( template_recom[b] ^ template );

           }


       }
                                                                               #if PRINT_RECOM
                                                                               printf(" |    "); for (b = 0; b < nSeq_block; b++){PrintBinary(template_recom[b]);  printf(" ");}   printf("\n");
                                                                               #endif


   }
                                                                               #if PRINT_RECOM
                                                                               printf(" +----------------------------------------------------------------------------+\n");
                                                                               #endif

    return template_recom;
}

struct Individual Recombination(struct Individual* recombinant, struct Individual* ind1, struct Individual* ind2, double recomRate, int nSeq, int nSeq_block,  unsigned int* nonsyn, double* delta, double pheno_wild)
{
    unsigned int* template_recom;
    int nRecom;
    int b;

    nRecom = poidev(recomRate * (nSeq-1) , &gseed);
                                                                                #if PRINT_RECOM
                                                                                nRecom = 2;
                                                                                printf("nRecomb %d ", nRecom);
                                                                                #endif
    if ( nRecom == 0 )
    {
        *recombinant = CopyIndividual(*recombinant, *ind1,  nSeq_block);
        return *recombinant;

    }

    template_recom = Get_Recombination_Template(nRecom, nSeq, nSeq_block);
    // make recombinant using template

                                                                                #if PRINT_RECOM
                                                                                printf(" | p1 "); for (b = 0; b < nSeq_block; b++){PrintBinary(ind1->seq[b]);  printf(" ");}   printf("\n");
                                                                                printf(" | tm "); for (b = 0; b < nSeq_block; b++){PrintBinary_asterik(template_recom[b]);  printf(" ");}   printf("\n");
                                                                                printf(" | p2 "); for (b = 0; b < nSeq_block; b++){PrintBinary(ind2->seq[b]);  printf(" ");}   printf("\n");
                                                                                #endif

    for (b = 0; b < nSeq_block; b++)
    {
        recombinant->seq[b] = ((~template_recom[b]) & ind1->seq[b]) + ((template_recom[b]) & ind2->seq[b]);
    }

                                                                                #if PRINT_RECOM
                                                                                printf(" | \n | re ");for (b = 0; b < nSeq_block; b++){PrintBinary(recombinant->seq[b]);  printf(" ");}   printf("\n");
                                                                                printf(" +----------------------------------------------------------------------------+\n");
                                                                                #endif


    UpdatePhenotype(recombinant, nonsyn, delta, pheno_wild, nSeq_block);
    free(template_recom);

    return *recombinant;
}
