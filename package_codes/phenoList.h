#ifndef __PHENOLIST_H_
#define __PHENOLIST_H_
//linked list for ranking phenotype
// -> modified linked list including int array for sequece info

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

typedef struct _nodeP
{
	struct _nodeP *next;
	double pheno;
	int numPheno1;
	int numPheno2;
}nodeP;


typedef nodeP* nptrP;

typedef struct _listP
{
	int count;
	nptrP head;
}listP;

listP* Create_PhenoList();
void Init_PhenoList(listP* lptr);
void Insert_Pheno(listP* lptr, int popnum, int numPheno, double pheno, int position, int p_nSeq_block);
void Print_Pheno_Outfile(listP* lptr, int gen, double selOpt, FILE *phenofile);
void Count_Pheno(listP* lptr, int popnum, double pheno, int p_nSeq_block);


#endif
