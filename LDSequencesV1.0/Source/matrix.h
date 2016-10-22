/***************************************************************************
						matrix.h  -  description
						---------------
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

// Matrix.h
#ifndef matrix_h							// Sentry, use file only if it's not already included.
#define matrix_h
#include <math.h>
#include "BaseDefinitions.h"

#define maxElem 20
//typedef int bool;
//typedef long double Zahl;
#define true 1
#define false 0
#define Pi 3.14159265358979323846264338328


//////////////////////////////////////////////////////////////////////////////////////////////
// Matrix-Objekt
//

class Matrix
{
	public:				 //Konstruktoren
		Matrix(int zSz=3,int sSz=-1, Zahl diag=1, Zahl other=0);
		Matrix(Matrix& x);
		~Matrix();

	protected:			//Speicherverwaltung
		bool AllocMxMem(int zSz,int sSz);
		bool FreeMxMem();
	protected:			//Daten
		Zahl *Elements;
		int zSize,sSize;

	public:		//Selektoren und Eigenschaften
		Zahl Get(const int i,const int j) {return (((i>zSize)||(j>sSize)||(zSize<=0)||(sSize<=0))?0:Elements[(i-1)*sSize+j-1]);}
		Zahl Set(int i,int j, const Zahl val);
		int GetSize(void) {return zSize;}
		int GetZSize(void) {return zSize;}
		int GetSSize(void) {return sSize;}
		bool Quadratic(void) {return sSize==zSize;}

/*public:				 //Operatoren
		Matrix&	operator =	(Matrix& x);
		Matrix&	operator *	(Matrix& x);
		Matrix&	operator +	(Matrix& x);
		Matrix&	operator - (Matrix& x);
		Matrix&	operator * (const Zahl x);*/

	public:				 //Rearrangements
		bool TauscheZeilen(int n, int m);
		bool TauscheSpalten(int n, int m);

/*	public:				 //Normen
		Zahl Schurnorm();
		Zahl ZSNorm();
		Zahl SSNorm();
*/
	public:				 //Determinante
		Zahl Det();
		protected:
			int GetBestPivot(int Col);

/*	public:				 //Datei Eingabe/Ausgabe
		bool Open(char* filename);
		bool Write(ofstream *str);

	public:				 //spezielle Matrizentypen
		bool Kronecker(float angle);
		bool WalshHardamard(void);
		bool Zielke(Zahl a);
		bool Hilbert(void);

	public:				 //Verfahren zur Inversion
		bool GaussJordan(void);
		bool Strassen(void);
		bool Schulz(Matrix& inv, char *fn="");
		bool LRInversion(void);

	protected:			//Teilmatrizen (+Multiplikation)
		Matrix GetTmx(int flag);
		bool SetMxFromTmx(Matrix& B11, Matrix& B12, Matrix& B21, Matrix& B22);
*/
};

#endif // Vektor_h
