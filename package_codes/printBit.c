//print seq
#include "printBit.h"

void PrintBinary(unsigned int input)
{
	int i, check = 1;
	char bit;

	for (i=_block_-1; i>=0; i--)
	{
		bit = (input & (1 << i))?1:0;
		//if (check)
		//{
			if (bit ==1)
			{
				check = 0;
				printf("%d", bit);
			}
		//}
		else printf("%d", bit);
	}
	//printf("\n");
}

void PrintBinary_out(unsigned int input, FILE* file )
{
	int i, check = 1;
	char bit;

	for (i=_block_-1; i>=0; i--)
	{
		bit = (input & (1 << i))?1:0;
		//if (check)
		//{
			if (bit ==1)
			{
				check = 0;
				fprintf(file,"%d", bit);
			}
		//}
		else fprintf(file,"%d", bit);
	}
	//printf("\n");
}


void PrintSeq_non_out(unsigned int* seq,unsigned int* nonsyn, FILE* file , int p_nSeq_block)	//has changed(180828) -> int p_nSeq_block instead of para_t p
{
	int i, j,  check = 1;
	char bit;
	unsigned int seg, template_seg;
	int count = 0;


	for (j = 0; j < p_nSeq_block; j++)
	{
		seg = seq[j];
		template_seg = nonsyn[j];
		count = 0;

		for (i=_block_-1; i>=0; i--)
		{


			if ((template_seg & (1 << i))?1:0)			//if this site is nonsyn site
			{
				//printf("block %d %d %d %lf", i,  count, delta_idx, z);
				bit =(seg & (1 << i))? 1: 0;

				if (bit ==1)
				{
					check = 0;
					fprintf(file,"%d", bit);
				}
				else
					fprintf(file,"%d", bit);



			}

		}

		fprintf(file,"/");
	}



	//printf("\n");
}

void PrintSeq_all_out(unsigned int* seq,unsigned int* nonsyn, FILE* file , int p_nSeq_block)	//has changed(180828) -> int p_nSeq_block instead of para_t p
{
	int i, j,  check = 1;
	char bit;
	unsigned int seg, template_seg;
	int count = 0;


	for (j = 0; j < p_nSeq_block; j++)
	{
		seg = seq[j];
		template_seg = nonsyn[j];
		count = 0;

		for (i=_block_-1; i>=0; i--)
		{



				//printf("block %d %d %d %lf", i,  count, delta_idx, z);
				bit =(seg & (1 << i))? 1: 0;

				if (bit ==1)
				{
					check = 0;
					fprintf(file,"%d", bit);
				}
				else
					fprintf(file,"%d", bit);





		}

		fprintf(file," ");
	}



	//printf("\n");
}


void PrintBinary_asterik(unsigned int input)
{
	int i, check = 1;
	char bit;

	for (i=_block_-1; i>=0; i--)
	{
		bit = (input & (1 << i))?1:0;
		//if (check)
		//{
			if (bit ==1)
			{
				check = 0;
				printf("*");
			}
		//}
		else printf("-");
	}
	//printf("\n");
}


void Print_NSinfo(unsigned int* seq, int p_nSeq_block, FILE* freq1, FILE* freq2, FILE* freq)
{

	int i, j;
	char bit;
	int count_n = 0;
	int count_s = 0;
	unsigned int seg;


	fprintf(freq,"gen ");
	fprintf(freq1,"gen ");
	fprintf(freq2,"gen ");


	for (j = 0; j < p_nSeq_block; j++)
	{
		seg = seq[j];

		for (i=_block_-1; i>=0; i--)
		{
			//printf("block %d %d %d %lf", i,  count, delta_idx, z);
				bit =(seg & (1 << i))? 1: 0;

				if (bit ==1)
				{
					fprintf(freq,"non%d ", count_n);
					fprintf(freq1,"non%d ", count_n);
					fprintf(freq2,"non%d ", count_n);
					count_n++;
				}
				else
				{
					fprintf(freq,"syn%d ", count_s);
					fprintf(freq1,"syn%d ", count_s);
					fprintf(freq2,"syn%d ", count_s);
					count_s++;
				}
		}
	}

	fprintf(freq,"\n");
	fprintf(freq1,"\n");
	fprintf(freq2,"\n");

}
