#include "writeFile.h"

FILE* Open_File_Write(FILE* file, const char* file_head, int try, const char* file_end)
{

    char file_name[100];

    if (try == -1)
    {
        sprintf(file_name,  "%s.%s", file_head, file_end);
    }
    else
    {
        sprintf(file_name,  "%s_try%d.%s", file_head, try, file_end);
    }

    file = fopen(file_name, "w");

    return file;


}


void Write_Freq_File(Population pop, int gen, para_t p, FILE* freq1, FILE* freq2, FILE* freq)
{
	//was GetHet originally

	#if PRINT_WHERE
	printf("\n\n==========PrintFreq==========\n\n");
	#endif
	#if PRINT_WHERE_SHORT
	printf("==========PrintFreq==========\n");
	#endif


	int i, j, k, check = 1;
	char bit;
	unsigned int input;
	int countAllele[_block_];
	int countAllele2[_block_];
	fprintf(freq,"%d ", gen);
	fprintf(freq1,"%d ", gen);
	fprintf(freq2,"%d ", gen);




	for (j = 0; j < p.nSeq_block; j++)
	{
		memset(countAllele, 0, sizeof(int)*_block_);
		memset(countAllele2, 0, sizeof(int)*_block_);


		for (k = 0; k < p.popsize[0]; k++)				//subpop0
		{
			input = pop[0].indiv[k].seq[j];

			for (i = _block_-1; i>=0 ; i--)
			{
				bit = ( input & ( 1 << i ))?1:0;
				if (bit ==1)
					countAllele[i]++;
			}
		}
		for (k = 0; k < p.popsize[1]; k++)				//subpop1
		{
			input = pop[1].indiv[k].seq[j];

			for (i = _block_-1; i>=0 ; i--)
			{
				bit = ( input & ( 1 << i ))?1:0;
				if (bit ==1)
					countAllele2[i]++;
			}
		}



		for (i = _block_-1; i>=0 ; i--)
		{
			fprintf(freq1, "%d ", countAllele[i]);
		}

		for (i = _block_-1; i>=0 ; i--)
		{
			fprintf(freq2, "%d ", countAllele2[i]);
		}


		for (i = _block_-1; i>=0 ; i--)
		{
			fprintf(freq, "%d ", countAllele[i] + countAllele2[i]);
		}




	}

	fprintf(freq,"\n");
	fprintf(freq1,"\n");
	fprintf(freq2,"\n");


}


void Write_Sampled_Sequence(Population pop, int gen, double selOpt, para_t p, FILE* file)
{
    int i,j, rInd, popNum, indNum;
    int nIndT = p.nIndT;
    int nSeq_block = p.nSeq_block;
    int nNon = p.nNon;
    int* popsize = p.popsize;

    for (i = 0; i < _nSample_ ; i++)
    {
        rInd = (int) (ran1(&gseed) * nIndT);
        if (rInd < popsize[0])
        {
            popNum = 0;
            indNum = rInd;
        }
        else
        {   popNum = 1;
            indNum = rInd-popsize[0];

        }


        fprintf(file, "%d\t%lf\t%d\t%lf\t", gen, selOpt, i, pop[popNum].indiv[indNum].pheno);
        for (j = 0; j < nSeq_block ; j++)
        {
             PrintBinary_out(pop[popNum].indiv[indNum].seq[j], file );

        }
        fprintf(file, "\n");

    }
}
