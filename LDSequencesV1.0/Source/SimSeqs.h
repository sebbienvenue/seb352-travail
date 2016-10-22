/***************************************************************************
						SimSeqs.h  -  description
						-----------------
		begin				: Sun May 29 2000
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


#ifndef SIMSEQS_H
#define SIMSEQS_H
//#include <stdlib.h>
//#include <stdio.h>
//#include <math.h>
#include "BaseDefinitions.h"

/** These constants define bitwise masks for variables that describe which sequence to use for your calculation.
e.g. 107 means MC, Halton, Niederreiter, Faure and Atanassov */
#define USEMC 1
#define USEHALTON 2
#define USENET 4
#define USENIEDER 8
#define USESOBOL 16
#define USEFAURE 32
#define USEATANASSOV 64
#define USEHAMMERSLEY 128
#define USENALPHA 256
#define USEMAX USENALPHA


	// constants for the recursive construction of the Halton sequence
#define GENAU		1E-9
#define BGENAU	 29
#define MAXLONG 2147483647
#define MAXPREC 30
#define primefile "primes.dat"

/** some numbes for the nAlpha sequence:
	NAlphaOnes	 is [0;1,1,1,1,1,1,1,1....]
	NAlphaTwos	 is [0;2,2,2,2,2,2,2,2....]
	NAlphaThrees is [0;3,3,3,3,3,3,3,3....]
	NAlphaFours	is [0;4,4,4,4,4,4,4,4....]
*/
#define NAlphaOnes .6180339887498948482045969343656391177203091798058
#define NAlphaTwos .4142135623730950488016887242096980785696718753769
#define NAlphaThrees .3027756377319946465596106337352479731256482869226
#define NAlphaFours .2360679774997896964091736687312762354406183596115

/* *********************************************************** */
/* Konstanten, Typen																				 */
/* *********************************************************** */
//TODO
#define maxanz 310
#define max_nieder 900000000 // Max. Anz.der Niederr. Vektoren
#define maxdim 310	// max. dimension for the sequence. set this to a high value to allow 400-dim. sequences as we need them.
#define MAX_I_NORMAL 1025

typedef double double_vektor[maxanz];
typedef int		int_vektor[maxanz];
typedef struct
				{
					int				deg;						// degree of the polynomial
					int_vektor coeff;					// coefficients
																		 // p(x)= coeff[0]+coeff[1]*x+ ...
				} polynom;
typedef polynom polarr[maxdim+4];
//endTODO */

/** This helper class describes polynomials with long coefficients.
	I know, this should probably be implemented as a Template class
	so that arbitrary data types are allowed. I'm just too lazy for now ;-)
*/
class longpoly {
protected:
			/** stores the degree of the polynomial
			*/
	long deg;
			/** if we are in the Ring Z_n, this is the base for all the modulo operations.
			*/
	long base;
public:
		/** initializes the polynomial with the given number of coefficients (order, default is 15)
				and in base bs (default is 0). If bs is set <=0, this means: "Don't do any modulo
				operations on the coefficients",i.e. we are in Z
		*/
	longpoly(long cflen=15, long bs=0);
	~longpoly();
			/** Returns the degree of the polynomial (memory allocated for the coefficients)
			*/
	long GetDegree();
			/** Sets a new possible degree for the polynomial by allocating sufficient memory (old coefficients are preserved)
			*/
	long SetDegree(long deg);
			/** Sets the k-th coefficient to val
			*/
	long SetCoeff(long k, long val);
			/** Gets the k-th coefficient
			*/
	long GetCoeff(long k);

			/** multiply this polynomial with another one (poly1). A new longpoly object is returned, which has to be free manually / by the garbage collection
			*/
	longpoly* times(longpoly& poly1);
//	longpoly* operator *(longpoly& poly1);

protected:
			/** this field stores the coefficients of the polynomial
			*/
	long*coeff;
			/** Size of memory allocated for the coefficients
			*/
	long coefflen;
};


//c lass LDSqBase;
		/** This function reads in the first m primes from the file "primes.dat" in the working directory
				and stores them to primes, which must be large enough to hold all m primes
		*/
long ReadPrimes(long* primes, long m);
		/** Returns the smalles prime larger than or equal to "current"
		*/
long NextPrime(long current);
		/** returns a field of length dm filled with the smalles prime numbers. Use this e.g. as bases for the Halton sequence.
				The long* field has to be delete[]'ed manually.
		*/
long*InitGenericPrimes(long dm);
		/** returns a field of length dm filled with the smalles Integer numbers. Use this e.g. as bases for the Niederreiter sequence in a base p.
				The long* field has to be delete[]'ed manually.
		*/
long*InitGenericInts(long dm);




#define NTAB 32

/**	This is the base class for all the simulation sequences including
		 Monte Carlo (just random numbers).
		 If you intend to use MC / QMC sequences in your program, use an object of LDSqBase* type,
		 so that you can later determine which sequence you really want.
		 To initialize the sequence, just initialize a subclass and type-cast it to LDSqBase*.
		 To get the next element, use the method NextElement(double*buffer, long bufflen, long nr=-1);
		 To get the name of the sequence use GetName() or GetExtension()
		 You probably won't need any other methods.

		 If you decide you want any other sequence, just use a diffent constructor, you don't have to change anything else.
		 e.g.
		 <pre>
		 void docalc() {
			 // definition of vars, memory for buffer, get bases etc.
			 // ...
			 LDSqBase*numbers;
			 numbers=(LDSqBase*)new MonteCarlo();
			 //numbers=(LDSqBase*)new Halton(bases, 5); // uncomment this if you want to use the Halton -sequence. This is the only place where you need to change something if you want a different sequence.

			 numbers->NextElement(buffer, bufflen); // This is the same, whatever sequence you use!!!
			 // .. do whatever you want
			 numbers->NextElement(buffer, bufflen, 271); // returns the 271. element of the sequence

			 //.. some more code,
			 delete numbers; // This is the same for all sequences, too
			 //free mem etc.
		 }
		 </pre>

		 So, you see, to use different sequences, just initialize a differen object, but since
		 all use LDSqBase as base class (which supplies the framework and the virtual functions
		 that are implemented by the other sequence classes), the rest of the program, where you
		 really work with the sequence just stays the same.
 @author Reinhold Kainhofer
*/
class LDSqBase {
protected:
			/** Strings, which contain the name and the 3-letter name of the sequence, e.g. name="Halton", ext="hal"
			*/
	char* name, *ext;
			/** number of dimensions (also the length of long* bases)
			*/
	long dim;
			/** Number of points of this point set (for sequences this has no important meaning, but for the Hammersley set this is vital)
			*/
	long len;
			/** long-field that stores the bases (if needed), e.g. for the Halton/Hammersley sequence
			*/
	long* bases;
			/** If all the numbers should be created at the beginning, this field stores all the elements. This
					can be quite effective if you need the same elements quite often in your code, since they are created
					just once at the beginning and not on the fly when needed (as it is the default
			*/
	double* numbers;
			/** flag to show whether the elements are pre-created into the field numbers already and just need
					to be copied from there. If false, the element needs to be calculated.
			*/
	bool created;
			/* needed for the Pseudo random generators. Stores state and some other coefficients for shuffling
			*/
	long state, iy, iv[NTAB], type;
public:
			/** store the index of the last created element, to use lastnr+1 as default... Some sequences
					can be calculated recursively (like the Halton sequences), so compare if the number needed
					is lastnr+1 (if yes, use the previous element to create the next one much faster
			*/
	long lastnr;

public:
			/** Determines whether the sequence is a Low Discrepancy Sequence or not. This is a flag
					that is set once in the constructor. Your Program can check
					if (numberObject->qmc) {Do for QMC-sequence} else {Do for pseudo-random sequence}
			*/
	/*const*/ bool qmc;
			/** Constructor, initializes the sequence.
					@param b Pointer to an array of longs containing the bases for the sequence (if needed, NULL or 0 otherwise). If the sequence needs bases, but NULL is given, it calls InitGenericBases(..) to create default bases (usually the lowest dim prime numbers)
					@param dm dimension of the sequence.
					@param iterations Number of points needed (e.g. used for the Hammesley-point set)
					@param genau Needed for the Halton-sequence (and maybe for some other sequences as well)
					@param genau1 like genau
					@param ex extension returned by GetExtension(), usually you leave this to its default value
					@param nm Name of the sequence, returned by GetName(), usually you leave this to its default value
			*/
	LDSqBase(long*b=NULL, long dm=0,long iterations=0,
			long genau=0, double genau1=0., char* ex="",char* nm="");
	virtual ~LDSqBase();

			/** Returns the dimension of the LD sequence (stored in the member variable dim)
			*/
	virtual long GetDimension() {return dim;};
			/** Returns len (the number of points of this point set or sequence)
			*/
	virtual long GetLength() {return len;};
	virtual char* GetName() {return name;};
	virtual char* GetExtension() {return ext;};
			/** Set the bases for the sequence (or change them if they
					were already given in the constructor)
			*/
	virtual long SetBases(long*b, long dm, long genau, double genau1);
			/** Returns the bases used for this sequence into the buffer b, at most dm of them and returns
					how many bases have been copied to b
			*/
	virtual long GetBases(long*b, long dm) {
			for (long i=0;i<min(dm, dim);i++) b[i]=bases[i]; return min(dm, dim);}
			/** Writes the nr-th element (if nr=-1 or not given, it writes the next element) of the sequence in a buffer.
					If bufflen>dim, the remaining dimensions are pseudo-random numbers
					(i.e. the sequence is a hybrid sequence). This is the method you should
					call from within your program.

					If the number have been pre-created (because CreateNumbers() has been called to create
					them), it just copies the appropriate element to the buffer, otherwise it will really
					calculate it.
			*/
	virtual long NextElement(double*buffer, long bufflen, long nr=-1);
			/** Creates all the numbers (len defines how many are needed) and stores them in the array
					numbers, so that NextElement(..) can just copy them from the array and does not have
					to calcualate them.
			*/
	virtual bool CreateNumbers();
			/** Changes the name of the sequence and the default extension (return value of GetName() and GetExtension()).
					They are initialized by the constructor
			*/
	virtual void SetNames(char* ex="",char* nm="");
			/** Sets the seed for the ran1, ran2 and ran3 pseudo-number generators
			*/
	virtual void SetSeed(long int status) {state=status;}
			/** ran1 pseudo-Random number generator from the "Numerical Recipies in C":
					"minimal" random number generator of Park and Miller with Bays-Durham shuffle and added safeguards
			*/
	double ran1(long*idum);
			/** ran2 improved pseudo-Random number generator from the "Numerical Recipies in C"
					long period rng of L'Ecuyer with Bays-Durham shuffle and added safeguards
			*/
	double ran2(long*idum);
			/** ran3 improved pseudo-Random number generator from the "Numerical Recipies in C"
					Knuth's subtractive method rng
			*/
	double ran3(long*idum);

protected:
			/** Initializes the method, allocates memory needed,sets default values for the
					variables needed and calls InitData to do class-specific initialization...
					@param b Pointer to an array of longs containing the bases for the sequence (if needed, NULL or 0 otherwise)
					@param dm dimension of the sequence.
					@param iterations Number of points needed (e.g. used for the Hammesley-point set)
					@param genau Needed for the Halton-sequence (and maybe for some other sequences as well)
					@param genau1 like genau
					@param ex extension returned by GetExtension(), usually you leave this to its default value
					@param nm Name of the sequence, returned by GetName(), usually you leave this to its default value
			*/
	virtual void InitMethod(long*b=NULL, long dm=0,long iterations=0,
			long genau=0, double genau1=0., char* ex="",char* nm="");
			/** Calculates the next element and stores it to the buffer. If bufflen>dim, the remaining dimensions are pseudo-random numbers (i.e. the sequence is a hybrid sequence)
					This function should never be called directly! Use NextElement(..) instead!!!
					@see NextElement()
			*/
	virtual long CalculateNextElement(long nr, double*buffer, long bufflen);
			/** initialize the method specific data
			*/
	virtual long InitData(long genau, double genau1);
			/** Do class-specific freeing of memory etc.
					Called by the destructor
			*/
	virtual long ExitData() { return 0;}
			/** If no array of bases was given in the constructor, but the sequence needs them, this
					creates generic bases, e.g. most of the time these are the dm lowest prime numbers,
					but you can override this by reimplementing InitGenericBases in your subclass
			*/
	virtual long*InitGenericBases(long dm) {return NULL;};
};



/** These gives the Monte Carlo Method (just pseudo-random numbers), defines constant for the LDSqMonteCarlo class
*/
#define MC1 1	// "minimal" random number generator of Park and Miller with Bays-Durham shuffle and added safeguards
#define MC2 2	// long period rng of L'Ecuyer with Bays-Durham shuffle and added safeguards
#define MC3 3	// Knuth's subtractive method rng

/** This class implements pseudo-random numbers. To create them it uses either
		ran1, ran2 or ran3 depending on the value for tp that you pass to the constructor.
		@author Reinhold Kainhofer
*/
class LDSqMonteCarlo : public LDSqBase	{
public:
			/** Constructor, initializes the Monte Carlo sequence.
					@param tp defines whether to use ran1, ran2 or ran3 (use the constants MC1, MC2 and MC3 for this)
					@param status initial seed for the RNGs
			*/
	LDSqMonteCarlo(long tp=MC1, long status=-44891368, char* ex="mc",char* nm="Monte Carlo");
protected:
			/** Generates an bufflen-dimensional pseudo-random vector and stores it to buffer
			*/
	virtual long CalculateNextElement(long nr, double*buffer, long bufflen) {
		switch (type) {
			case MC2:	for (long i=0;i<bufflen;i++) { buffer[i]=ran2(&state);}; break;
			case MC3:	for (long i=0;i<bufflen;i++) { buffer[i]=ran3(&state);}; break;
			case MC1:default: for (long i=0;i<bufflen;i++) { buffer[i]=ran1(&state);}; break;
		}
		return dim;
	}
	virtual long InitData(long genau, double genau1){return dim;};
	virtual long ExitData(){return 0;};
};


/** Afflerbach is a special Pseudorandomnumber generator that uses a linear congruence operator
		(quite fast).
		@author Reinhold Kainhofer (code taken from Thomas Siegl's programs)
*/
class LDSqAfflerbach : public LDSqMonteCarlo {
			/** Store the old state (don't use status, since then we could not use ran1, ran2 or ran3 any
					more in this class
			*/
	long altx;
public:
			/** Sets the state "altx" to a generic value of 12785
			*/
	LDSqAfflerbach(char* ex="aff", char*nm="Afflerbach"):LDSqMonteCarlo(MC1,0,ex,nm) {altx=12785;}
protected:
	virtual long CalculateNextElement(long nr, double*buffer, long bufflen) {
		for (long i=0;i<bufflen;i++) buffer[i]=createNumber();
		return dim;
	}
private:
			/** This function implements the linear congruence operator
			*/
	virtual double createNumber(){
		long xx=(532393*altx+1);
		while (xx>=268435456) xx-=268435456;
		altx=xx;
		return ((double)xx/(double)268435456);
	}
};


/** This class implements the N *Alpha Sequence. No check on the independence of the alphas is done!!!
*/
class LDSqNAlpha : public LDSqBase {
public:
	LDSqNAlpha(Zahl*alph=NULL, long dm=0, long iterations=0, char*ex="nalpha", char*nm="[n Alpha]"):
			LDSqBase(NULL, dm, iterations, 0,0,ex,nm) {
			alpha=alph;
			InitMethod(NULL, dm, iterations, 0,0,ex,nm); }
protected:
			/** Here we should create dm independent bases into alpha. Not done yet!!!
			*/
	virtual long*InitGenericBases(long dm) {return NULL;/*TODO!!!*/};
			/** Create the nr-th element by a simple multiplication nr*alpha (mod 1)
			*/
	virtual long CalculateNextElement(long nr, double*buffer, long bufflen) {
		long d=min(bufflen,dim);
		for (long i=0; i<d;i++) buffer[i]=(nr*alpha[i])-floor((nr*alpha[i]));
		return d; }
			/** This double field stores the alpha-values
			*/
	Zahl*alpha;
};

/** This class implements the Halton Sequence in the bases stored in b passed to
		the constructor. The VanDerCorput Sequence is just a special case with dim=1.
*/
class LDSqHalton : public LDSqBase	{
public:
	LDSqHalton(long*b=NULL, long dm=0,long iterations=0, long genau=BGENAU,
			double genau1=GENAU, char* ex="hal",char* nm="Halton"):
			LDSqBase(b, dm, iterations, genau, genau1, ex, nm){
			p=NULL; q=NULL; lastnr=0; InitMethod(b, dm, iterations, genau, genau1,	ex, nm);}
protected:
			/** Calculates the next element from the current element (or if we need a different element
					it calculates it using the function Phi. If bufflen>dim=length(bases), the remaining
					dimensions are filled with pseudo-random numbers.
					@see Phi()
			*/
	virtual long CalculateNextElement(long nr, double*buffer, long bufflen);
			/** Precalculate some values needed to create the numbers much faster (recursively)
			*/
	virtual long InitData(long genau, double genau1);
	virtual long ExitData();
	virtual long*InitGenericBases(long dm) {return InitGenericPrimes(dm);};
			/** Precision (in powers of 1/2). We need this many precreated numbers 1/2^i and the other corresponding value 1-1/2^1 in p and q
			*/
	long prec;
			/** Stores 1/2^i and 1-1/2^i to speed up the calculation (recursive!!!)
			*/
	double *q,*p;
};



/** This class implements Atanassov's modified Halton sequence
*/
class LDSqAtanassov : public LDSqBase	{
protected: // Protected attributes
			/** stores the factors of Atanassov's definition
			*/
	long* factors;
			/** stores the primitive roots needed to create the sequence
			*/
	long* proots;
			/** stores the modifiers	needed in Atanassov's definition
			*/
	long* modifiers;
public:
	LDSqAtanassov(long*b=NULL, long dm=NULL,long iterations=0, long genau=BGENAU,
			double genau1=GENAU, char* ex="ata",char* nm="Atanassov", long*prts=NULL,long*mdfs=NULL):
			LDSqBase(b, dm, iterations, genau, genau1, ex, nm){
			proots=prts; modifiers=mdfs; InitMethod(b, dm, iterations, genau, genau1,	ex, nm);}
				/** writes the next element(s) of the sequence in a buffer. */
protected:
	virtual long CalculateNextElement(long nr, double*buffer, long bufflen);
	virtual long InitData(long genau=0, double genau1=0);
	virtual long ExitData();
	virtual long*InitGenericBases(long dm) {return InitGenericPrimes(dm);};
private:
			/** For Atanassov's modified Halton-Sequence we need a bunch of primitive roots. They are calculated here.
			*/
	bool GeneratePrimitiveRoots(long*bases, long*prts, long dim);
			/** The modifiers are calculated from the primite roots using this function
			*/
	bool GenerateModifiers(long*bases, long*prts, long*modfs, long dim);
};

/** This class implements the Hammersley Set by creathing the Halton-sequence in
		dim-1 dimensions and adding the additional dimension of the form k/N. Here the
		parameter len is really important!!!
*/
class LDSqHammersley : public LDSqHalton	{
protected:
		 /** npos gibt an, an welcher Stelle das n/iterations eingefuegt werden
				soll. Der Index is 0-basierend, negative Werte sind von rechts
				gezaehlt bzw. der Wert ist modulo dim
		 */
	int Npos;
public:
		 /** npos gibt an, an welcher Stelle das n/iterations eingefuegt werden
				soll. Der Index is 0-basierend, negative Werte sind von rechts
				gezaehlt bzw. der Wert ist modulo dim
		 */
	LDSqHammersley(long*b=NULL, long dm=0,long iterations=0, long genau=0,
			double genau1=0., int npos=1, char* ex="ham",char* nm="Hammersley"):
			LDSqHalton(b,dm-1, iterations, genau, genau1,ex,nm) {
			Npos=((npos%dm)+dm)%dm;dim++;};
				/** write the next elements of the sequence in a buffer. */
protected:
	virtual long CalculateNextElement(long nr, double*buffer, long bufflen);
};


#define konst 1073741824
/** This class implements the Sobol Sequence, direction numbers are only available
		for 51 dimensions, so we can only provide that many dimensions
*/
class LDSqSobol : public LDSqBase	{
protected:
	long*V;
public:
	LDSqSobol(long*b=NULL, long dm=0,long iterations=0, long genau=0,
			double genau1=0., char* ex="sob",char* nm="Sobol"):
			LDSqBase(b, dm, iterations, genau, genau1, ex,nm){
			V=NULL;InitMethod(b, dm, iterations, genau, genau1,ex,nm);}
				// write the next element(s) of the sequence in a buffer.
protected:
	virtual long CalculateNextElement(long nr, double*buffer, long bufflen);
	virtual long InitData(long genau, double genau1);
	virtual long ExitData();
private:
	void dirnum_erz();
};



/** Implementation of the Faure low discrepancy sequence, use a linear field c
		to store the pascal matrix and multiply it on the previous coefficients
		by matrix-multiplication to get the new coefficients
*/
class LDSqFaure : public LDSqBase	{
protected:
	long r, base;
	long*c;
public:
	// fstep contains base^fstp
	long fstep;
public:
	LDSqFaure(long*b=NULL, long dm=0,long iterations=0, long genau=0,
			double genau1=0., long fstp=0, char* ex="fau",char* nm="Faure");
protected:
	virtual long CalculateNextElement(long nr, double*buffer, long bufflen);
	virtual long InitData(long genau, double genau1);
	virtual long ExitData();
private:
		/** multiply the current vector of digits of the expansion with the Pascal Matrix to get
				the vector containing the digits for the next dimension
		*/
	void cyMultiply(long*cc, long*yold, long*ynew, long dm, long coefflen, long bs);
protected:
	virtual void CreateDirMx(long*cc, long len, long bs);
};

/** Implementation of the (0,s) nets, according to an algorithm given by Lecot
		for Niederreiter's construction using hyperderivatives of polynomials
*/
class LDSqNetz : public LDSqBase	{
protected:
	long base, nu, lambda;
	long nprime, ncurr, ncmax;
	long*ncoeff, *ypji;
	long ypjilen;
public:
	LDSqNetz(long*b=NULL, long bas=5, long dm=0,long iterations=0,
			long lamb=5, char* ex="ts" ,char* nm="(t,s)-Sequence");
protected:
	virtual long CalculateNextElement(long nr, double*buffer, long bufflen);
	virtual long InitData(long genau, double genau1);
	virtual long ExitData();
			/** This net sequences just need dm different integers, no need for them to be prime
					numbers, so use the dm lowest integers instead to minimize the discrepancy
			*/
	virtual long*InitGenericBases(long dm) {return InitGenericInts(dm);};
private:
	long getypji(long p, long j, long i);
	long getypji(long index);
	long setypji(long p, long j, long i, long value);
	long setypji(long index, long value);
	long getindex(long p, long j, long i);
};

/** Implementation of (0,s) nets, according to an algorithm given by Niederreiter
		using monic polynomials
*/
class LDSqNiederreiter : public LDSqBase	{
protected:
	long basis;
	long maxnum,log_N;
	int_vektor q;
	long max_num;
	long c[32][32][maxanz]; //	c_ij acc. B-F-N Niederr
	polarr monic_pols;
//	longpoly**MonicPolys;
//	long nu, lambda;
//	long nprime, ncurr, ncmax;
//	long*ncoeff, *ypji;
//	long ypjilen;
public:
	LDSqNiederreiter(long*b=NULL, long bas=5, long dm=0,long iterations=0,
			char* ex="nie" ,char* nm="Niederreiter (t,s)-Sequence");
protected:
	virtual long CalculateNextElement(long nr, double*buffer, long bufflen);
	virtual long InitData(long genau, double genau1);
	virtual long ExitData();
private:
	void pol_mult(polynom ap, polynom bp, polynom* cp);
	int entw_in_basis (int n, int *a);
	void gen_poly(polynom* pol );
};



#endif
