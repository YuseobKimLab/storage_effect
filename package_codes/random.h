#ifndef __RANDOM_H_
#define __RANDOM_H_


#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)
#define PI 3.141592654

#include <math.h>


float gammln(float xx);
float poidev(float xm, long *idum);
float bnldev(float pp, int n, long *idum);
float gasdev(long *idum);
float expdev(long *idum);
float ran1(long *idum);

extern long gseed;

#endif
