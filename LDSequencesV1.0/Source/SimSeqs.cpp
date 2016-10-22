/***************************************************************************
						SimSeqs.cpp  -  description
						-------------------
		begin				: Sun May 28 2000
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
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <iostream.h>
#include <fstream>
#include "SimSeqs.h"
#include "matrix.h"
#include "BaseDefinitions.h"


unsigned int W[53][32]; // This field is filled with values as given by Siegl



long kill_last_zeros(long *n) {		//	F"ur LDSqSobol, LDSqSobol+
	long rest,hilf;

	hilf=0;
	do {
		rest=(*n) % 2;
		(*n) /= 2;
		hilf++;
	}	while (rest!=1);
	return(hilf-1);
}


	/**	this function is needed for Atanassov's modified Halton sequence.
		It calculates the matrix A of values a_ij with
		prts[i]^{a_ij}=bases[j] (mod bases[i])
		@author: Reinhold Kainhofer */
bool CreateModExponents(long*A, long dim, long*bases, long*prts){
	long tmp=0;
	long val=0;
	long bs=0;
	bool noerror=true;
	for (long i=0; i<dim; i++) {
		for (long j=0; j<dim; j++) {
			if (i==j) A[i*dim+i]=0;
			else {
				bs=bases[j] % bases[i]; // compare with this var
				val=prts[i] % bases[i]; // initial value for a_ij=1
				tmp=1;  // holds a_ij
				while( (val != bs) && (tmp<=bases[i]) ) {
					tmp++;
					val=(val*prts[i]) % bases[i];
				}
				if (val==bs) A[i*dim+j]=tmp;
				else {
					noerror=false;
					A[i*dim+j]=0;
				}
			}
		}
	}
	return noerror;
}

long expand(long*A, bool*rowsUsed, long dim, long dimleft, long dimtocalc) {
	long fact=-1;
	long res=0;
	for (long i=0; i<dimtocalc; i++) {
		if (!rowsUsed[i]) {
			fact*=(-1);
			rowsUsed[i]=true;

			res+=fact*A[dim*i + dimleft]*expand(A, rowsUsed, dim, dimleft-1, dimtocalc);
			rowsUsed[i]=false;
		}
	}
	return res;
}

	/**	calculate the determinant via the expansion according to the last column.
		This is really slow, but never mind...
		flag the rows used so that we don't need another array for the cofactor.
		This saves a lot of mem... */
long CalcDeterminant(long*A, long dim, long dimtocalc) {
	// Calculate the determinant of A_i, which is the submatrix
	// of the first i rows and columns of A
	Matrix mx((int)dimtocalc);
	for (int i=0; i<dimtocalc; i++) {
		for (int j=0; j<dimtocalc; j++) {
			mx.Set(i+1,j+1,(Zahl)A[i*dim+j]);
		 }
	}
	long res=(long)rint(mx.Det());
	return res;
}

	/**	this function calculates the entries k_i for Atanassov's
		modified Halton sequence. The matrix A needs to have a determinant
		of 1, so this is used to determine k_i.
		@author: Reinhold Kainhofer  */
void MakeDet1(long*A, long dim, long i) {
	// set a_ii =0, so the new a_ii then has to be equal 1-det(A_i),
	// since det(A_{i-1}) == 1 by induction. A_i ist the matrix of the
	// first i rows and columns
	long determ=CalcDeterminant(A, dim,i);
	A[i*dim+i]=1-determ;
	return;
}

/* ***************************************************************
 *           reading in the primes
 ****************************************************************/

long ReadPrimes(long* primes, long m) {
	ifstream stream(primefile);
	if (!stream) return 0;
	for (long i=0;i<m;i++) {
		stream>>primes[i];
	}	return m;
}

long NextPrime(long number) {
	long pr=0;
	if (number<50) {
	// for small numbers don't bother to look into the file, explicitely check
		pr=number;
		if (pr%2==0 && pr>2) pr++;
		switch (pr) {
			case 1: return 2;
			case 3: return 3;
			case 5: return 5;
			case 7: return 7;
			case 9:
			case 11: return 11;
			case 13: return 13;
			case 15:
			case 17: return 17;
			case 19: return 19;
			case 21:
			case 23: return 23;
			case 25:
			case 27:
			case 29: return 29;
			case 31: return 31;
			case 33:
			case 35:
			case 37: return 37;
			case 39:
			case 41: return 41;
			case 43: return 43;
			case 45:
			case 47: return 47;
			case 49:
			case 51: return 51;
			default: return 2;
		}
	} else {
		ifstream stream(primefile);
		if (!stream) {
			cout<< "Error reading the primes file \"" << primefile<<"\n";
			return -1;
		}
		while (pr<number) stream>>pr;
	}
	return (long)pr;
}

long*InitGenericPrimes(long dm) {
	long*res=new long[dm];
	ReadPrimes(res, dm);
	return res;
}
long*InitGenericInts(long dm) {
	long*res=new long[dm];
	for (long i=0; i<dm; i++) res[i]=i/*+1*/;
	return res;
}

/* ================================================================= */


longpoly::longpoly(long cflen, long bs) {
	coefflen=cflen;
	coeff=new long[coefflen];
	memset(&coeff[0], 0,	coefflen*sizeof(long));
	deg=0;
	base=bs;
}

longpoly::~longpoly() {
	Free(coeff);
}

long longpoly::GetDegree(){
	return deg;
}
long longpoly::SetDegree(long degree) {
	long i;
	long*oldcoeff=coeff;
	if (degree>coefflen) {
		coeff=new long[degree+5];
		memset(&coeff[0], 0,	(degree+5)*sizeof(long));
		for (i=0;i<deg; i++) coeff[i]=oldcoeff[i];
		Free(oldcoeff);
		coefflen=degree+5;
	}
	for (i=deg; i<degree; i++) coeff[i]=0;
	deg=degree;
	return deg;
}
long longpoly::SetCoeff(long k, long val) {
	if(k>=deg)
		SetDegree(k);
	if (base>0) coeff[k]=val%base;
	else coeff[k]=val;
	if (val==0 && deg==k) {
		while (coeff[deg]==0) deg--;
	}
	return k;
}
long longpoly::GetCoeff(long what) {
	if (what<=deg) return coeff[deg];
	else return 0;
}

//longpoly* longpoly::operator *(longpoly& poly1) {
longpoly* longpoly::times(longpoly& poly1) {
	longpoly*result=(longpoly*)new longpoly( GetDegree() + poly1.GetDegree() );
	long coeffk;
	long i,k;

	for (k=0; k<result->GetDegree(); k++) {
		coeffk=0;
		for (i=0; i<k; i++) {
			coeffk+= (GetCoeff(i) * poly1.GetCoeff(k-i) );
			if (base>0) coeffk%=base;
		}
		result->SetCoeff(k, coeffk);
	}
	return result;
}



// ****************************************************************************
// LDSqBase is the base class for all the sequences.
// For every sequence, only the InitData(int, double), SetNames(), ExitData(),
// CalculateNextElement(double*, long) methods need to be overridden.
// ****************************************************************************

LDSqBase::LDSqBase(long*b, long dm, long iterations,
		long genau, double genau1, char* ex,char* nm) {
	qmc=true;
	numbers=NULL;
	bases=NULL;
	created=false;
	state=-44891368;
//	InitMethod(b, dm, iterations, genau, genau1);
}
void LDSqBase::InitMethod(long*b, long dm, long iterations,
		long genau, double genau1, char* ex,char* nm) {
	numbers=NULL;
	SetNames(ex,nm);
	dim=dm;
	if (b!=NULL) {
		bases=new long[dm];
		for (long i=0;i<dm;i++) {bases[i]=b[i];}
	} else {bases=InitGenericBases(dm);}
	len=iterations;
	lastnr=0;
	InitData(genau, genau1);
}

LDSqBase::~LDSqBase(){
	Free(numbers);
	ExitData();
	Free(bases);
}

long LDSqBase::CalculateNextElement(long nr, double*buffer,
		long bufflen) {return -1;}
long LDSqBase::NextElement(double*buffer, long bufflen, long nr) {
	long Nr;
	if (buffer==NULL || bufflen==0) return 0;
	if (nr<0) Nr=lastnr+1; else Nr=nr;
	lastnr=Nr;
	if (created) {
		long to=min(dim, bufflen);
		for (long i=0; i<to; i++) {
			buffer[i]=numbers[Nr*dim+i];
		}
		return min(dim, bufflen);
	} else return CalculateNextElement(Nr, buffer, bufflen);
}

	/**	Set the bases for the sequence (or change them if they were already given
		in the constructor */
long LDSqBase::SetBases(long*b, long dm, long genau, double genau1){
	bases=b;
	dim=dm;
	return InitData(genau, genau1);
}

	/** Set the name of the funciton and the default extension */
void LDSqBase::SetNames(char* ex, char* nm){
	name=nm;
	ext=ex;
}

bool LDSqBase::CreateNumbers() {
	Free(numbers);
	numbers=new double[dim*len];
	memset(&numbers[0], 0,	dim*len*sizeof(Zahl));
	if (numbers==NULL) return false;
	for (long i=0; i<len; i++) {
		CalculateNextElement(i, &numbers[i*dim], dim);
	}
	ExitData();
	created=true;
	return true;
}

	/** initialize the method specific data */
long LDSqBase::InitData(long genau, double genau1){
	return dim;
}


// ****************************************************************************
// MonteCarlo implements Pseudo Random Numbers
// ****************************************************************************
LDSqMonteCarlo::LDSqMonteCarlo(long tp, long status, char* ex,char* nm): LDSqBase(NULL,0,0,0,0){
	qmc=false;
	iy=0;
	type=tp;
	state=status;
	InitMethod(NULL, 0,0,0,0, ex, nm);
}

#define IM1 2147483399
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)
double LDSqBase::ran2(long *idum) {
/* Long perios (>2x10^18) random number gnerator of L'Ecuyer with
	 Bays-Durham shuffle and added safeguards. Returns a uniform
	 random deviate between 0.0and 1.0 (exclusive of the endpoint
	 values). Call with idum a negative integer to initialize;
	 thereafter, do not alter idum between successive deviates in a
	 sequence. RNMX should approximate the latgest floating value
	 that is less than 1 */
	int j;
	long k;
	static long idum2=123456789;
	static long iy=0;
	static long iv[NTAB];
	double temp;

	if (*idum<=0) {	// initialize
		if (-(*idum)<1) *idum=1;	// be sure to prevent idum=0
		else *idum=-(*idum);
		idum2=(*idum);
		for (j=NTAB+7; j>=0; j--) {	// load the shuffle table (after 8 warm-ups).
			k=(*idum)/IQ1;
			*idum=IA1*(*idum-k*IQ1)-k*IR1;
			if (*idum<0) *idum+=IM1;
			if (j<NTAB) iv[j]=*idum;
		}
		iy=iv[0];
	}
	k=(*idum)/IQ1; // start here when not initializeng
	*idum=IA1*(*idum-k*IQ1)-k*IR1;	//compute idum=(IA1*idum) % IM1 without overflows by Schrage's method.
	if (*idum<0) *idum+=IM1;
	k=idum2/IQ2;
	idum2=IA2*(idum2-k*IQ2)-k*IR2;	// compute idum2=(IA12*idum)%IM2 likewise
	if (idum2<0) idum2+=IM2;
	j=iy/NDIV; // will be in the range 0.. NTAB-1
	iy=iv[j]-idum2;	// here idum is shufled, idum and idum2 are combined to generate output
	iv[j]=*idum;
	if (iy<1) iy+=IMM1;
	if ((temp=AM*iy)>RNMX) return RNMX;	// because users don't expect endpoint values.
	else return temp;
}


#define IA 16807
#define IM 2147483647
#define IQ 127773
#define IR 2836
double LDSqBase::ran1(long *idum) {	
	long k, j=0;
	double temp;

	if (*idum<=0 || !iy) {
	 	if (-(*idum)<1) *idum=1;
		else *idum=-(*idum);
		for (j=NTAB+7;j>=0;j--) {
		 	k=(*idum)/IQ;
			*idum=IA*(*idum-k*IQ)-IR*k;
			if (*idum<0) *idum += IM;
			if (j<NTAB) iv[j]=*idum;
		}
		iy=iv[0];
	}
	k=(*idum)/IQ;
	*idum=IA*(*idum-k*IQ)-IR*k;
	if (*idum<0) *idum += IM;
	j=iy/NDIV;
	iy=iv[j];
	iv[j]= *idum;
	if ((temp=AM*iy)>RNMX) return RNMX;
	else return temp;
} //ran1

#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC (1.0/MBIG)
// according to Knuth, any large MBIG and any smaller (but still
// large) MSEED can be substituted for the above values.
double LDSqBase::ran3(long*idum) {
/* returns a uniform random deviate between 0.0 and 1.0. Set
	 idum to any negative value to initialize or reinitialize
	 the sequence. */
	static int inext, inextp;
	static long ma[56]; // the value 56 (range ma[1..55]) is special and should not be modified; see Knuth
	static int iff=0;
	long mj,mk;
	int i,ii,k;

	if (*idum<0||iff==0) { // initialization
		iff=1;
		mj=MSEED-(*idum<0 ? -*idum : *idum); // initialize ma[55] using the seed idum and the large number MSEED.
		mj %= MBIG;
		ma[55]=mj;
		mk=1;
		for (i=1; i<=54; i++) { // not initialize the rest of the table, in a slightly random order, with numbers that are not especially random.
			ii=(21*i)%55;
			ma[ii]=mk;
			mk=mj-mk;
			if (mk<MZ) mk+= MBIG;
			mj=ma[ii];
		}
		for (k=1; k<=4;k++) // we randomize them by "warming up the generator"
			for (k=1; i<=55; i++) {
				ma[i]-=ma[1+(i+30)%55];
				if(ma[i]<MZ) ma[i]+=MBIG;
			}
		inext=0;// prepare inidices for our first generatred number.
		inextp=31;	// the constant 31 is special; see Knuth
		*idum=1;
	}
	// here is where we start, except on initialization.
	if (++inext==56) inext=1; // increment inext and inextp, wrapping around 56 to 1
	if (++inextp==56) inextp=1;
	mj=ma[inext]-ma[inextp];	// generate a new random nubmer subtractively
	if (mj<MZ) mj+=MBIG;	// be sure that it is in range
	ma[inext]=mj;	// store it
	return mj*FAC; // and output the derived uniform deviate
}




// ****************************************************************************
// LDSqHalton implements the Halton sequence, and thus also the VanDerCorput
// sequence, which is just the Halton sequence in one dimension
// ****************************************************************************


long LDSqHalton::InitData(long genau, double genau1) {
	ExitData();
	prec=genau;

	p=new double[dim*(prec+2)];
	q=new double[dim*(prec+2)];
	memset(&p[0], 0,	dim*(prec+2)*sizeof(double));
	memset(&q[0], 0,	dim*(prec+2)*sizeof(double));

	for (long i=0;i<dim;i++) p[i*(prec+2)]=1;

	for (long j=0;j<=prec;j++) {
		for (long i=0;i<dim;i++) {
			p[i*(prec+2)+j+1]=p[i*(prec+2)+j]/(double)bases[i];
			q[i*(prec+2)+j]=(1.0-p[i*(prec+2)+j])*(1.0-genau1);
		}
	}
	return dim;
}

/** Assuming buffer already holds the old sequence */
long LDSqHalton::CalculateNextElement(long nr, double*buffer, long bufflen) {
	long i,d=min(bufflen,dim);
//	if (bufflen>dim) cout<<"Error: bufflen("<<bufflen<<")>dim("<<dim<<")\n";
	if ((nr!=lastnr) || ((nr%10000)==0)) {
		for (i=0; i<d;i++) {
			buffer[i]=Phi(bases[i], nr);
			if(buffer[i]>1|| buffer[i]<0) {
				cout<<"error: "<<buffer[i]<<"(b="<<bases[i]<<", i="<<i<<")\n";
			}
		}
	} else {
		for (i=0;i<d;i++) {
			long k=1;
			while (buffer[i]>=q[i*(prec+2)+k] && (k<prec)) {k++;}
			if (k==1) {buffer[i]+=(double)(1./bases[i]);}
			else {buffer[i]+=(p[i*(prec+2)+k-1]+p[i*(prec+2)+k]-1);}
			if (buffer[i]<0 || buffer[i]>=1) {
				buffer[i]=0;
			}
		}
	}
	for (i=d; i<bufflen; i++) {
		buffer[i]=ran2(&state);
	}
	return d;
}


long LDSqHalton::ExitData() {
	LDSqBase::ExitData();
	Free(p);
	Free(q);
	return 0;

}

// ****************************************************************************
// LDSqAtanassov implements the modified Halton sequence introduced by Atanassov.
// If no primitive roots are given, the lowest possible values are taken.
// ****************************************************************************


long LDSqAtanassov::InitData(long genau, double genau1) {
	long*tmp;
	long i;

	tmp=proots;
		// allocate memory for the primitive roots
	proots=new long[dim];
	memset(&proots[0], 0, dim*sizeof(long));
	if (tmp==NULL) GeneratePrimitiveRoots(bases, proots, dim);
	else { // just copy them over
		for (i=0; i<dim; i++) proots[i]=tmp[i];
	}
		// and the same with the modifiers.
	tmp=modifiers;
		// allocate memory for the primitive roots
	modifiers=new long[dim];
	memset(&modifiers[0], 0, dim*sizeof(long));
	if (tmp==NULL) GenerateModifiers(bases, proots, modifiers, dim);
	else { // just copy them over
		for (i=0; i<dim; i++) modifiers[i]=tmp[i];
	}
	factors=new long[dim];
	memset(&factors[0], 0, dim*sizeof(long));
	for (i=0; i<dim; i++) factors[i]=power(proots[i], modifiers[i])%bases[i];
	return dim;
}

long LDSqAtanassov::CalculateNextElement(long nr, double*buffer, long bufflen) {
	long* cfs=new long[MAXPREC];
	memset(&cfs[0], 0, MAXPREC*sizeof(long));
	long i,d=min(bufflen, dim);
//	if (bufflen>dim) cout<<"Error: bufflen("<<bufflen<<")>dim("<<dim<<")\n";
	for (i=0; i<d; i++) {
		LongToCoeff(nr, bases[i], cfs, MAXPREC);
		buffer[i]=CoeffToDoubleModified(bases[i], cfs, MAXPREC, factors[i]);
	}
	for (i=d; i<bufflen; i++) {
		buffer[i]=ran2(&state);
	}
	Free(cfs);
	return dim;
}


long LDSqAtanassov::ExitData() {
	LDSqBase::ExitData();
	Free(factors);
	Free(modifiers);
	Free(proots);
	return 0;

}

bool LDSqAtanassov::GeneratePrimitiveRoots(long*bases, long*prts, long dim){
	bool res=true;
	bool found;
	long i,j;
	for (i=0; i<dim; i++) {
		found=false;
		j=0;
		while((!found)&&(j<bases[i])) {
			j++;
			found=IsPrimitiveRoot(j, bases[i]);
		}
		prts[i]=j;
		res=(res && (j<bases[i]));
	}
	return res;
}

bool LDSqAtanassov::GenerateModifiers(long*bases, long*prts, long*modfs, long dim){
	long*A = new long[dim*dim];
	memset(&A[0], 0, dim*dim*sizeof(long));
	long i;
	if (CreateModExponents(A, dim, bases, prts)) {
		for (i=0; i<dim; i++) {
			 /* Calculate a_ii such that det(A)=1	=> the system of
					diophantine equations then has integer solutions
					(see proof in the thesis or in Atanassov's paper)*/
			 MakeDet1(A, dim, i);
			 modfs[i]=A[i*dim+i];
		}
		Free(A);
		return true;
	} else {
		Free(A);
		return false;
	}
}


// ****************************************************************************
// LDSqHammersley implements the Hammersley sequence, which is the Halton Sequence
// with one more dimension of the form n/N */
// ****************************************************************************

/** Assuming buffer already holds the old element */
long LDSqHammersley::CalculateNextElement(long nr, double*buffer, long bufflen) {
	dim--;
	long lngth=min(dim+1, bufflen);
	long i;

		// Reset to 0 so that every timestep uses the same set
	if ((nr%len) ==0) {
		for (i=0;i<lngth; i++) {
			buffer[i]=0;
		}
	}

		// Shift the elements to get rid of the additional entry
	for (i=Npos; i<dim; i++) buffer[i]=buffer[i+1];
		// then call the Halton sequence
	long res=LDSqHalton::CalculateNextElement(nr, buffer, bufflen-1);
		// and shift the elements back
	for (i=dim-1; i>=Npos; i--) buffer[i+1]=buffer[i];
	buffer[Npos]=(double)((double)(nr%len))/(double)len;
	dim++;
	return res;
}



// ****************************************************************************
// LDSqSobol implements the Sobol sequence, which is more or less a permuation
// of the VanDerCorpus sequence in base two
// ****************************************************************************
Zahl mmin(Zahl x, Zahl y) {
	if (x<=y) return x;
	else return y;
}
long lmin(long x, long y) {
	if (x<=y) return x;
	else return y;
}

long LDSqSobol::InitData(long genau, double genau1) {
	ExitData();
	dirnum_erz();
//	long r=lfloor(log(MAXLONG)/log(3));
	V=new long[lmin(dim, 51)];
	memset(&V[0], 0, lmin(dim, 51)*sizeof(long));
	return 0;
}

long LDSqSobol::CalculateNextElement(long nr, double*buffer, long bufflen) {
	long i, n=nr+1, ki, d=lmin(bufflen, lmin(dim, 51));
//	if (bufflen>dim) cout<<"Error: bufflen("<<bufflen<<")>dim("<<dim<<")\n";
	ki=kill_last_zeros(&n);
	for (i=0;i<d; i++) {
		V[i]=V[i]^W[i+1][ki+1];
		buffer[i]=(double)V[i]/konst;
	}
	for (i=d;i<bufflen; i++) {
		buffer[i]=ran2(&state);
	}
	return bufflen;
}

// the direction numbers for Sobol's sequence up to dimension 51
// are defined in the file SobolDirNum.i
void LDSqSobol::dirnum_erz() {
#include "SobolDirNum.i"
}

long LDSqSobol::ExitData() {
	LDSqBase::ExitData();
	Free(V);
	return 0;
}



// ****************************************************************************
// LDSqFaure implements the Faure sequence
// ****************************************************************************

LDSqFaure::LDSqFaure(long*b, long dm,long iterations, long genau, double genau1,
		long fstp, char* ex,char* nm):
		LDSqBase(b, dm, iterations, genau, genau1, ex, nm){
	r=0;
	base=0;
	c=NULL;
	fstep=power(base, fstp);
	InitMethod(b, dm, iterations, genau, genau1, ex, nm);
}


long LDSqFaure::InitData(long genau, double genau1) {
	// as base use the smallest prime larger than the dimension
	base=NextPrime(dim);
	if (base<0) return -1;
	r=(long)(log(MAXLONG)/log(base));
	c=new long[r*r];
	memset(&c[0], 0, r*r*sizeof(long));
	CreateDirMx(&c[0], r, base);
	return dim;
}


/** Create the Pascal matrix. If you want a different uppre triangular matrix,
		just overload this function in your subclass and create your own direction
		numbers
*/
void LDSqFaure::CreateDirMx(long*cc, long len, long bs) {
	if (base>0) {
		// Create the pascal matrix mod base into c
		long i, j;
		for (j=0; j<r; j++) c[j]=1;
		for (i=1; i<r; i++) {
			for (j=i*r+i; j<(i+1)*r; j++) {
				c[j]=(c[j-1]+c[j-1-r]) % base;
			}
		}
	}
}

void LDSqFaure::cyMultiply(long*cc, long*yold, long*ynew, long dm, long coefflen, long bs){
/* y=c.yold */
	for (long i=0; i<coefflen; i++) {
		ynew[i]=0;
		long first=i*dm;
		for (long j=i; j<coefflen; j++) {
			ynew[i]+=(yold[j]*cc[first+j]);
		}
		ynew[i]%=bs;
	}
}

long LDSqFaure::CalculateNextElement(long nr, double*buffer, long bufflen) {
	if (base<0) return 0;
	long d=min(dim, bufflen), k;
//	if (bufflen>dim) cout<<"Error: bufflen("<<bufflen<<")>dim("<<dim<<")\n";
	long* a=new long[r];
	memset(&a[0], 0, r*sizeof(long));
	long* aold=new long[r];
	memset(&aold[0], 0, r*sizeof(long));

		// Darstellung von n-1 in Basis base

	long tmp=nr+fstep;
	long j=LongToCoeff(tmp, base, &a[0], r);

		// daraus x1 erstellen (wie bei van der corput)
	buffer[0]=CoeffToDouble(base, a, j);
	for (k=1; k<d;k++){
		memcpy(&aold[0], &a[0], j*sizeof(long));
			// mit c multiplizieren ergibt neues y
		cyMultiply(c, aold, a, r, j, base);
			// sind Koeffizienten von xk in Basis q
		buffer[k]=CoeffToDouble(base, a, j);
	}
	Free(a);
	Free(aold);
	for (k=d;k<bufflen; k++) {
		buffer[k]=ran2(&state);
	}

	return d;
}

long LDSqFaure::ExitData() {
	LDSqBase::ExitData();
	Free(c);
	return 0;
}


// ****************************************************************************
// LDSqNetz implements (0,s)-Sequences using hyperderivatives and Lecot's
// recursive construction
// ****************************************************************************

LDSqNetz::LDSqNetz(long*b, long bas, long dm,long iterations, long lamb,
		char* ex,char* nm):LDSqBase(b, dm, iterations,0,0, ex, nm){
	base=bas;
	lambda=lamb;
	ncurr=0;
	nprime=0;
	InitMethod(b, dm, iterations, 0,0, ex, nm);
}
long LDSqNetz::InitData(long genau, double genau1){
	long p=0, bl=0, ind=0;
	long indexoffset=0;
	long i=0, j=0, l=0, m=0, tmp=0;
	
	nu=power((long)base,(long)lambda);
	ncmax=(long)(log(MAXLONG)/log(base)) + 1; // Worst case, but memory does no matter that much...
	ncoeff=new long[ncmax];
	memset(&ncoeff[0], 0, ncmax*sizeof(long));
	
	// just store the coefficients in an nu × lambda rectangle.
	// Many of them are 0, but not much less than half of them, so
	// there is not much to gain by trying to save only the non-zero
	// ones. This would just require an enormous calculation effort...
	ypjilen=(long)(nu*lambda*dim);
	ypji=new long[ypjilen];
	memset(&ypji[0], 0, ypjilen*sizeof(long));
	for (p=0;p<base;p++) {
		for (i=0; i<dim; i++) {
			setypji(p,0,i,p%base);
		} //i=0..dim
	} //p=0..b

	for (l=1; l<lambda; l++) {
		for (m=1; m<base; m++) {
			bl=power(base,l);
			indexoffset=m*bl*dim*lambda;
			for (p=m*bl; p<(m+1)*bl; p++) {
				ind=getindex(p, 0, 0);
				for (j=0; j<l; j++) {
					for (i=0; i<dim; i++) {
						tmp=ypji[ind-indexoffset]; //getypji(p-m*bl, j,i)
						tmp+=(long)(binom(l,j))*(long)(power(bases[i], l-j))*(long)(m);

						setypji(ind, tmp);
						ind++;
					} //i=0...s (Dimension)
				} //j=0...l
				for (i=0;i<dim; i++) {
					setypji(ind, m);
					ind++;
				}
			} //p=m*b^l...(m+1)*b^l
		} // m=1...b-1
	} // l=1...lambda
	return 0;
};

long LDSqNetz::getindex(long p, long j, long i) {
	return i+dim*(p*lambda+j);
}
long LDSqNetz::getypji(long index) {
	if (0<=index<ypjilen) return ypji[index];
	else return 0;
}
long LDSqNetz::getypji(long p, long j, long i) {
	long index=getindex(p,j,i);
	if (0<=index<ypjilen) return ypji[index];
	else return 0;
}
long LDSqNetz::setypji(long index, long value) {
	long oldval;
	if (0<=index<ypjilen) {
		oldval=ypji[index];
		ypji[index]=value % base;
		return oldval;
	} else return -1;
}

long LDSqNetz::setypji(long p, long j, long i, long value) {
	long index=getindex(p,j,i);
	return setypji(index, value);
}
long LDSqNetz::CalculateNextElement(long nr, double*buffer, long bufflen){
	long k, ypjinew=0, i=0, j=0;
	long d=min(dim, bufflen);
//	if (bufflen>dim) cout<<"Error: bufflen("<<bufflen<<")>dim("<<dim<<")\n";
	long index=0;
long ypjinew1;	
long ss;
long ss1;
long ss2;
long ss3;
	for (i=0; i<d; i++) buffer[i]=0;
	if (nr>=nu) {
		if ((ncurr*nu>nr)||(nr>=(ncurr+1)*nu)) {
			ncurr=(long)nr/nu;
			nprime=LongToCoeff(ncurr, base, ncoeff, ncmax);
		}
		for (j=0; j<lambda+nprime; j++){
			for (i=0; i<d; i++) {
				if (j>=lambda) ypjinew=0;
				else ypjinew=getypji(nr-ncurr*nu, j, i);

				for (k=max(j, lambda); k<lambda+nprime; k++) {
/*			long ss=(long)ncoeff[k-lambda];
					long ss1=(long)power(bases[i], k-j);
					long ss2=(long)binom(k,j);
					long ss3=ss1*ss2*ss;*/
					ss=(long)ncoeff[k-lambda];
					ss1=(long)power(bases[i], k-j);
					ss2=(long)binom(k,j);
					ss3=ss1*ss2*ss;
					ypjinew+=(long)((long)binom(k,j)*(long)power(bases[i],k-j)*
							(long)ncoeff[k-lambda]);
				} //k=max()...lambda+nprime
				ypjinew1=ypjinew % base;
				buffer[i]+=(double)((double)ypjinew1/(double)power(base,j+1));
			} //i=1..s
		} //j=1...lambda+nprime+1
	} else {
		index=getindex(nr, 0,0);
		for (j=0; j<lambda; j++) {
			for (i=0; i<d; i++) {
				buffer[i]+=(double)((double)ypji[index]/(double)power(base, j+1));
				index++;
			}
		}
	}
	for (i=d; i<bufflen; i++) {
		buffer[i]=ran2(&state);
	}
	return 0;
};

long LDSqNetz::ExitData(){
	LDSqBase::ExitData();
	Free(ncoeff);
	Free(ypji);
	return 0;
};


// ****************************************************************************
// LDSqNiederreiter implements (0,s)-Sequences using Niederreiter's algorithm
// using monic polynomials
// ****************************************************************************
LDSqNiederreiter::LDSqNiederreiter(long*b, long bas, long dm,long iterations,
		char* ex,char* nm):LDSqBase(b, dm, iterations,0,0, ex, nm){
	basis=bas;
	InitMethod(b, dm, iterations, 0,0, ex, nm);
}

long LDSqNiederreiter::ExitData(){
	LDSqBase::ExitData();
	return 0;
};


void LDSqNiederreiter::pol_mult(polynom ap, polynom bp, polynom* cp) {
	int i,j,k,hilf;
	// multipliziert zwei Polynome cp=ap*bp

	cp->deg= ap.deg+bp.deg;

	// obere Koeffizienten mit Nullen fuellen
	for (i= ap.deg+1; i <= cp->deg; i++) ap.coeff[i]= 0;
	for (j= bp.deg+1; j <= cp->deg; j++) bp.coeff[j]= 0;

	for (k=0; k <= cp->deg; k++) {
		cp->coeff[k]=0;
		for (i=0; i <= k; i++) {
			hilf= ap.coeff[i]*bp.coeff[k-i];
			cp->coeff[k]= (cp->coeff[k]+hilf) % basis;
		}
	}
}


int LDSqNiederreiter::entw_in_basis (int n, int *a)	{
// n in basis entwickeln
	int i=-1,r,mm;
	int q;

	do {
		q=n / basis;
		r=n % basis;
		a[++i]=r;
		n=q;
	} while (n != 0);
	mm = i;
	for (i=mm+1;i<=log_N;i++) a[i]=0;
	return(mm);
}


void LDSqNiederreiter::gen_poly(polynom* pol ) {
// generate monic polynomials: x, x+1, x+2, ... , x+BASE-1,... *
	int i,j,n=0,anz,d,is_irreducible;
	polynom a,b,c;
bool stop;

	anz=0;
	do {
//	for(n=0;;n++) {
		a.deg=entw_in_basis(n++,a.coeff);
		a.coeff[++a.deg]=1;
		is_irreducible=0;
		for(i=0;;i++) {
			b.deg=entw_in_basis(i,b.coeff);
			b.coeff[++b.deg]=1;
			if (b.deg>a.deg/2) { is_irreducible=1; break; }
			for(j=0;j<=a.deg;j++) c.coeff[j]=a.coeff[j];
			c.deg=a.deg;
			do {
				d=c.coeff[c.deg];
				for (j=0;j<=b.deg;j++)
					c.coeff[c.deg-j]=(c.coeff[c.deg-j]-d*b.coeff[b.deg-j])%basis;
				while((c.deg>=0)&&(c.coeff[c.deg]==0))
					c.deg--;
			}	while(c.deg>=b.deg);
			if (c.deg==(-1)) { is_irreducible=0; break; }
		}
		if (is_irreducible==1) {
			for(j=0;j<=a.deg;j++) pol[anz].coeff[j]=a.coeff[j];
			pol[anz].deg=a.deg;
			anz++;
//			if (anz>min(dim, maxdim)) break;
		}
		stop=anz<dim && anz<maxdim;
	} while (stop);//(anz<min(dim, maxdim));
}


// ************************************************************************
// *																																			*
// *	Erzeugung des nummer-ten NIEDER Vektors x der Dimension dimension	 *
// *																																			*
// ************************************************************************
long LDSqNiederreiter::CalculateNextElement (long nr, double *buffer, long bufflen) {
// x->buffer
// nummer->nr
	long d=min(min(dim, bufflen), maxdim);
//	if (bufflen>dim) cout<<"Error: bufflen("<<bufflen<<")>dim("<<dim<<")\n";

	int_vektor a;
	int rt,it,ht,jt,D=0;
	int pot,qn;

	entw_in_basis(nr,a);
	for (it=0;it<d;it++) { // it=0..d-1, not it=1..d, as in Siegl
		qn=0;
		pot=max_num / basis;
		for (jt=1;jt<=log_N;jt++) {
			D=0;
			for (rt=0;rt<log_N;rt++) {
				ht=(c[jt][rt][it]*a[rt]) % basis;
				D+=ht;
				D%=basis;
			};
			qn+=D*pot;
			pot/=basis;
		}
		buffer[it]=(double)qn/max_num;
	};
	for (it=d; it<bufflen; it++) buffer[it]=ran2(&state);
	return d;
}


long LDSqNiederreiter::InitData(long genau, double genau1){
// genarates the elements C_ij acc. B-F_N *
// variables used coincide with those in the paper *
	int i,j,k,l,m,e,q,r,u,hilfi;
	int_vektor v;
	polynom hilfp, b;

	log_N=0;
	do {
		log_N++;
		max_num=1;
		for(i=0;i<log_N;i++)
			max_num*=basis;
	} while((double)max_nieder/max_num>basis);

long rr=min(dim, maxdim);
	gen_poly(monic_pols);
	for (i=0; i<rr; i++) {
		e=monic_pols[i].deg;
		j=0;
		q=(-1);
		u=e;
		b.deg=0;
		b.coeff[0]=1;
		do {
			j++;
			if (u==e) {
				q++;
				u=0;
				pol_mult(b,monic_pols[i],&hilfp);
				m=hilfp.deg;
				for (k=0;k<m;k++)
					b.coeff[k]=(basis-hilfp.coeff[k]) % basis;

				b.coeff[m]=1;
				b.deg=m;

				for (k=0; k < m-1; k++) v[k]= 0;
				v[m-1]= 1;
				for (k=m; k <= (log_N+e-2);k++) {
					v[k]=0;
					for (l=1; l <= m; l++) {
						hilfi=(b.coeff[m-l]*v[k-l]) % basis;
						v[k]= (v[k]+ hilfi) % basis;
					}
				}
			}
			for (r=0;r<log_N;r++)
				c[j][r][i]= v[r+u];
			u++;

		} while (j<log_N);
	}

	// don't start at 0, since too many points are near 0!!!
	lastnr=1;
	for(i=0;i<log_N-1;i++)
		lastnr*=basis;
	return dim;
}
