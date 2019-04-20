#ifndef __WRITEFILE_H_
#define __WRITEFILE_H_

#include <stdio.h>
#include <string.h>
#include "parameter.h"
#include "pop.h"
#include "random.h"
#include "printBit.h"


FILE* Open_File_Write(FILE* file, const char* file_head, int try, const char* file_end);

void Write_Freq_File(Population pop, int gen,para_t p, FILE* freq1, FILE* freq2, FILE* freq);
void Write_Sampled_Sequence(Population pop, int gen, double selOpt, para_t p, FILE* file);
#endif
