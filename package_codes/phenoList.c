#include "phenoList.h"


listP* Create_PhenoList()
{
	listP* phenoList = (listP*)malloc(sizeof(listP));
	phenoList->count = 0;
	phenoList->head = NULL;

	return phenoList;
}

void Init_PhenoList(listP* lptr)
{
	//initialize the list
	if (lptr->head != NULL)
	{
		nptrP tmp1 = lptr->head;
		nptrP tmp2;

		while (tmp1 != NULL)
		{
			tmp2 = tmp1;
			tmp1 = tmp2->next;

			free(tmp2);

		}

	}

	lptr->count = 0;
	lptr->head = NULL;

}

void Insert_Pheno(listP* lptr, int popnum, int numPheno, double pheno, int position, int p_nSeq_block)
{

	//insert value to the proper postion
	if ( position < 1 || position > (lptr->count)+1 )
	{
		printf("Position Out of Bound\n");
		return;
	}


	nptrP new_nptr = (nodeP*)malloc(sizeof(nodeP));
	if (popnum == 0)
	{
		new_nptr->numPheno1 = numPheno;
		new_nptr->numPheno2 = 0;
	}

	else
	{
		new_nptr->numPheno2 = numPheno;
		new_nptr->numPheno1 = 0;
	}

	new_nptr->pheno = pheno;


	if (position == 1)
	{
		new_nptr->next = lptr->head;
		lptr->head = new_nptr;
	}
	else
	{
		nptrP tmp = lptr->head;
		int i;
		for (i = 1; i <position-1; i++)
		{
			tmp = tmp->next;
		}
		new_nptr->next = tmp->next;
		tmp->next = new_nptr;
	}
	lptr->count++;
}

void Print_Pheno_Outfile(listP* lptr, int gen, double selOpt, FILE *phenofile)
{
	double seq;
	double pheno;
	int i, rank=0;


	nptrP tmp = lptr->head;
	while(tmp != NULL)
	{
		fprintf(phenofile, "%d %lf %lf %d %d %d \n", gen, selOpt, tmp->pheno, tmp->numPheno1 + tmp->numPheno2, tmp->numPheno1, tmp->numPheno2);
		tmp = tmp->next;

	}

	// sel opt 추가해서 인쇄하기
}

void Count_Pheno(listP* lptr, int popnum, double pheno, int p_nSeq_block)
{
	//transverse the list and
	//find the first position of the value (first form head)
	//if not exist, return 0

	nptrP tmp = lptr->head;
	int i = 1;
	int j;
	int same = 1;
	while ( tmp != NULL )
	{
		same = 1;
		for ( j = 0; j < p_nSeq_block; j++)
		{
			if (pheno != tmp->pheno)
			{
				same = 0;
			}
		}

		if (same == 1 )
		{
			if (popnum == 0)
				tmp->numPheno1++;
			else
				tmp->numPheno2++;
			break;
		}

		i++;
		tmp = tmp->next;
	}

	if (i > lptr->count )
	{
		Insert_Pheno(lptr, popnum, 1, pheno, i,  p_nSeq_block);		//inside this funtion, we reset the other population as 0
	}

}
