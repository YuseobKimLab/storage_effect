#ifndef __PRINTBIT_H_
#define __PRINTBIT_H_

//print seq
#include "option.h"
#include <stdio.h>


void PrintBinary(unsigned int input);
void PrintBinary_out(unsigned int input, FILE* file );
void PrintSeq_non_out(unsigned int* seq,unsigned int* nonsyn, FILE* file , int p_nSeq_block);
void PrintSeq_all_out(unsigned int* seq,unsigned int* nonsyn, FILE* file , int p_nSeq_block);
void PrintBinary_asterik(unsigned int input);
void Print_NSinfo(unsigned int* seq, int p_nSeq_block, FILE* freq1, FILE* freq2, FILE* freq);


#endif
