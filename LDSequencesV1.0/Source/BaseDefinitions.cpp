/***************************************************************************
					   BaseDefinitions.cpp  -  description
						 -------------------
		begin				: Wed Apr 25 2001
		copyright			: (C) 2001 by Reinhold Kainhofer
		email				: reinhold@kainhofer.com
 ***************************************************************************/

/***************************************************************************
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 ***************************************************************************/

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "BaseDefinitions.h"

long ERFCMAXREK=7;

Zahl AsymErfc(Zahl zsq) {
	Zahl result=1;
	Zahl am=1;
	for (long m=1; m<ERFCMAXREK; m++) {
		am*=((-1)*(2.*m-1.)/(2.*zsq));
		result+=am;
	}
	return result;
}

/* ***************************************************************
 *  writing out the results (arrays) in Mathematica format
 ****************************************************************/

void WriteDoubleArray(ofstream &str, long len, double* vals, char* endchar, long dim) {
	long k;
	if (dim<=1) {
		str<<"{ "<<vals[0];
		for (k=1;k<len;k++) str<<", "<<vals[k];
		str<<"}"<<endchar<<"\n";
	} else {
		str<<"{ ";
		for (k=0;k<len-1;k++) {
			WriteDoubleArray(str, len, &vals[k*power(len, dim-1)],",", dim-1);
		}
		char tmp[30];
		strcpy(tmp, "}");
		strcat(tmp, endchar);
		strcat(tmp,"\n");
		WriteDoubleArray(str, len, &vals[(len-1)*power(len, dim-1)], tmp, dim-1);
	}
}


/* ============================================================ */
// General mathematical functions.

long binom(long i, long k) {
	if (k>i||k<0) return 0;
	long j, mx=max(k, i-k);
	Zahl res=1.;
	
	for (j=mx+1; j<=i; j++) res=res* (Zahl)j / (Zahl)(i-j+1) ;
//	for (j=1; j<=i-mx; j++) res/=j;
	return (long)(res+0.005);
};

long power(long m, long n) {
	long res=1;
	for (long i=1;i<=n;i++) res*=m;
	return res;
}
double power(double m, long n) {
	double res=1;
	for (long i=1;i<=n;i++) res*=m;
	return res;
}

double vectabs(double a1, double a2, double a3) {
	return sqrt(a1*a1+a2*a2+a3*a3);
}


bool IsPrimitiveRoot(long proot, long base) {
	bool res=true;
	long j=0;
	long tmp=1;
	while (res && (j<base-2)) {
		j++;
		tmp=(tmp*proot % base);
		if (tmp==1) res=false;
	}
	if (res) {
		tmp=(tmp*proot % base);
		if (tmp!=1) res=false;
	}
	return res;
}


double Phi(long base, long x){
	long tmp=x;
	double res=0;
	long j=1, rem=0;

	while (tmp>0) {
			rem=tmp % base;
			tmp/=base;
			res+=(rem*pow(base, -j));
			j++;
	}
	return res;
}


double CoeffToDouble(long base, long* cfs, long cflen) {
	double res=0.;
	for (long i=cflen-1; i>=0; i--) {
		res+=(double)(cfs[i]);
		res/=(double)base;
	}
	return res;
}

double CoeffToDoubleModified(long base, long* cfs, long cflen, long modfact) { //TODO: use modfact somewhere!!!!!!
	double res=0.;
	for (long i=cflen-1; i>=0; i--) {
		res+=(double)(cfs[i]);
		res/=(double)base;
	}
	return res;
}


long LongToCoeff(long nr, long base, long*cfs, long cflen) {
	long nn=nr;
	long maxcoeff=0;

	for (long i=0; i<cflen; i++) cfs[i]=0;
	while(nn>0 && maxcoeff<cflen) {
		cfs[maxcoeff]=nn%base;
		nn=(long)nn/base;
		maxcoeff++;
	}
	return maxcoeff;
}
