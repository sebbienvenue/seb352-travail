/***************************************************************************
						BaseDefinitions.h  -  description
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

#ifndef BASEDEFINITIONS_H
#define BASEDEFINITIONS_H

//#include <stdlib.h>
#include <stdio.h>
#include <fstream.h>

	/**	We want to use double, but who knows, maybe sometimes we need long double or just float. By using our own type we only have to change this here and then to recompile
	*/
#define Zahl double
#define min(a,b) (a<b)?a:b
#define epsilon 0.0000000001
	/**	approximate value for the square root of pi
	*/
#define M_SQRTPI 1.772453850905516027298167483341145182798
#define M_1SQRT2PI 0.398942280401432677939946059934	/* 1/sqrt(1*pi)*/


	/**	Defines how many term in the asymptotic expansion of erfc=1-erf should be used for the calculation
	*/
extern long ERFCMAXREK;

	/**	Calculates the Erfc by means of its asymptotic expansion.
		@param zsq The variable z^2 in the expansion.
	*/
Zahl AsymErfc(Zahl zsq);
	/**	Frees the array-pointer ptr and sets it to NULL afterwards
	*/
#define Free(ptr) delete[] ptr;ptr=NULL;

	/**	Writes the "dim"-dimensional array vals in Mathematica-List format to the file str
		@param endchar defines a string which should be written out after the list, e.g. ";" or "//MatrixForm"
		@param len
		@param vals The dim-dimensional array of double-Values
		@param endchar String that is appended to the file after the List is written out.
		@param dim
	*/
void WriteDoubleArray(ofstream &str, long len, double* vals, char* endchar, long dim);

#define lfloor(f) (long)floor(f)



	// general math functions
	/**	Returns m^n, where n and m are long variables. The result, of course, is also of type long.
	*/
long power(long m, long n);
	/**	Returns m^n where n and m are double variable. The result is of type double, too.
	*/
double power(double m, long n);
	/**	Returns the binomial coefficient (i k)
	*/
long binom(long i, long k);
	/**	Returns the length of the 3-dimensional vector (a1, a2, a3) (just sqrt(a1*a1 + a2*a2 + a3*a3) )
	*/
double vectabs(double a1, double a2, double a3);
	/**	Macro, returns the minimum of two given values
	*/
#define min(x,y) (x<y)?x:y
	/**	Macro, returns the maximum of two given values
	*/
#define max(x,y) (x>y)?x:y


	/**	The digit inversion function. This inverts x (written down in base "base") at the comma. E.g.
		Phi(3, 25): 25=2*9+2*3+1=221_3 => 0.112_3=1/3 + 1/9 + 2/27
	*/
double Phi(long base, long x);

	/**	Expands n in base "base".
		@param nr the number to develop in base
		@param base Which base to use for the expansion
		@param cfs long-array where the coefficients should be stored
		@param cflen length of the buffer cfs
		@return Returns the number of coefficients (lowest n s.t. nr<base^n) or cflen, depending in which one is smaller
	*/
long LongToCoeff(long nr, long base, long* cfs, long cflen);
	/**	returns the double-representation of the cofficients in base "base" given in cfs. E.g.
		cfs={1,2,1,2,1}, base=3 => returns 1/3 + 2/9 + 1/27 + 2/81 + 1/243
		This is the modified version, needed by Atanassov's sequence. (Not quite implemented yet ...)
	*/
double CoeffToDoubleModified(long base, long* cfs, long cflen, long modfact);
	/**	returns the double-representation of the cofficients in base "base" given in cfs. E.g.
		cfs={1,2,1,2,1}, base=3 => returns 1/3 + 2/9 + 1/27 + 2/81 + 1/243
	*/
double CoeffToDouble(long base, long* cfs, long cflen);
	/**	Checks if proot is a primitive root modulo base
	*/
bool IsPrimitiveRoot(long proot, long base);

	/**	calculates the determinant of the (upper left) i*i submatrix of the dim*dim Matrix A
	*/
long CalcDeterminant(long*A, long dim, long i);



#endif
