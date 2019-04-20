#include "random.h"



float ran1(long *idum)
{
	int j;
	long k;
	static long iy=0;
	static long iv[NTAB];
	float temp;

	if (*idum <= 0 || !iy) {
		if (-(*idum) < 1) *idum=1;
		else *idum = -(*idum);
		for (j=NTAB+7;j>=0;j--) {
			k=(*idum)/IQ;
			*idum=IA*(*idum-k*IQ)-IR*k;
			if (*idum < 0) *idum += IM;
			if (j < NTAB) iv[j] = *idum;
		}
		iy=iv[0];
	}
	k=(*idum)/IQ;
	*idum=IA*(*idum-k*IQ)-IR*k;
	if (*idum < 0) *idum += IM;
	j=iy/NDIV;
	iy=iv[j];
	iv[j] = *idum;
	if ((temp=AM*iy) > RNMX) return RNMX;
	else return temp;
}


//gammlm
float gammln(float xx)
{
	double x, y, tmp, ser;
	static double cof[6] = {76.18009173,-86.50532033,24.01409822,-1.231739516,0.120858003e-2,-0.536382e-5};
	int j;
	
	y = x = xx;
	tmp = x +5.5;
	tmp -= (x+ 0.5) * log(tmp)	;
	ser = 1.000000000190015;
	for (j=0;j<=5;j++)	ser += cof[j]/++y;
	return -tmp+log(2.5066282746310005*ser/x);
}
	
//poidev 
#define PI 3.141592654
float poidev(float xm, long *idum)
{
	float gammln(float xx);
	float ran1(long *idum);
	static float sq, alxm, g, oldm=(-1.0);
	float em, t, y;
	
	if (xm < 12.0)
	{
		if (xm != oldm)
		{
			oldm = xm;
			g = exp(-xm);
		}
		em = -1;
		t = 1.0;
		do{
			++em;
			t *= ran1(idum);
		}
		while (t>g);
	}
	else
	{
		if (xm != oldm)
		{
			oldm = xm;
			sq = sqrt(2.0*xm);
			alxm = log(xm);
			g = xm*alxm - gammln(xm+1.0);
		}
		do {
			do {
					y = tan(PI*ran1(idum));
					em = sq*y+xm;
			}while (em < 0.0);
			em = floor(em);
			t = 0.9*(0.1+y*y)*exp(em*alxm-gammln(em+1.0)-g);
		}while (ran1(idum)>t);
	}
	return em;
}


float bnldev(float pp, int n, long *idum)
{
	int j;
	static int nold=(-1);
	float am,em,g,angle,p,bnl,sq,t,y;
	static float pold=(-1.0),pc,plog,pclog,en,oldg;
	
	p=(pp <= 0.5 ? pp : 1.0-pp);
	am=n*p;
	if (n < 25) {
		bnl=0.0;
		for (j=1;j<=n;j++)
			if (ran1(idum) < p) bnl += 1.0;
	} else if (am < 1.0) {
		g=exp(-am);
		t=1.0;
		for (j=0;j<=n;j++) {
			t *= ran1(idum);
			if (t < g) break;
		}
		bnl=(j <= n ? j : n);
	} else {
		if (n != nold) {
			en=n;
			oldg=gammln(en+1.0);
			nold=n;
		} if (p != pold) {
			pc=1.0-p;
			plog=log(p);
			pclog=log(pc);
			pold=p;
		}
		sq=sqrt(2.0*am*pc);
		do {
			do {
				angle=PI*ran1(idum);
				y=tan(angle);
				em=sq*y+am;
			} while (em < 0.0 || em >= (en+1.0));
			em=floor(em);
			t=1.2*sq*(1.0+y*y)*exp(oldg-gammln(em+1.0)
				-gammln(en-em+1.0)+em*plog+(en-em)*pclog);
		} while (ran1(idum) > t);
		bnl=em;
	}
	if (p != pp) bnl=n-bnl;
	return bnl;
	
}

float gasdev(long *idum)
{
	static int iset = 0;
	static float gset;
	float fac, rsq, v1, v2;
	
	if (*idum < 0) iset=0; 
	if (iset == 0) 
	{
		 do 
		 {
			v1=2.0*ran1(idum)-1.0; // pick two uniform numbers in the square extending from -1 to +1 in each direction, 
			v2=2.0*ran1(idum)-1.0; 
			rsq=v1*v1+v2*v2; // see if they are in the unit circle,
		}
		while (rsq >= 1.0 || rsq == 0.0); // and if they are not, try again. 
		fac=sqrt(-2.0*log(rsq)/rsq); /* Now make the Box-Muller transformation to get two normal deviates. Return one and save the other for next time. */
		gset=v1*fac;
		iset=1; // Set ag. 
		return (double) v2*fac; 
	} 
	else 
	{ // We have an extra deviate handy, 
		iset=0; // so unset the ag, 
		return (double) gset; // and return it. 
	} 
}
		
		

float expdev(long *idum)
{
	float ran1(long *idum);
	float dum;

	do
		dum=ran1(idum);
	while (dum == 0.0);
	return -log(dum);
}