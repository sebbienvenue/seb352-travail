/***************************************************************************
						matrix.cpp  -  description
						-----------------
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

// Matrix.cpp

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "matrix.h"

#define Abbruch 1e-25
#define maxRek 100
#define LeftTop 1
#define LT Lefttop
#define RightTop 2
#define RT RightTop
#define LeftBottom 3
#define LB LeftBottom
#define RightBottom 4
#define RB RightBottom



//////////////////////////////////////////////////////////
//	Hilfsfunktionen
/*Zahl max(Zahl value1, Zahl value2)
{	return ((value1>value2)?value1:value2);}*/




/******************************************************************
 *					MATRIX-OBJEKT
 ******************************************************************/



//////////////////////////////////////////////////////////
// Konstruktoren / Destruktor

Matrix::Matrix(int zSz,int sSz, Zahl diag, Zahl other) {
	int n, m;
	if (sSz<=0) sSz=zSz;

	AllocMxMem(zSz,sSz);
	for (n=1; n<=zSize; n++)
		for (m=1; m<=sSize;m++)
			if (m==n) Set(n,m,diag); else Set(n,m,other);
}

Matrix::Matrix(Matrix& x) {
	int i,j;
	AllocMxMem(x.GetZSize(),x.GetSSize());
	for (i=1;i<=zSize;i++)
		for(j=1;j<=sSize;j++)
			Set(i,j,x.Get(i,j));
}

Matrix::~Matrix() {
	FreeMxMem();
}



//////////////////////////////////////////////////////////
// Speicherverwaltung

bool Matrix::AllocMxMem(int zSz,int sSz) {
	zSize=zSz;
	sSize=sSz;

	if (!(Elements = new Zahl[zSize*sSize]) ) {
		zSize=0;
		sSize=0;
		Elements=NULL;
	}
	memset(&Elements[0], 0, zSize*zSize*sizeof(Zahl));
	return Elements==NULL;
}

bool Matrix::FreeMxMem() {
	Free(Elements);
	zSize=0;
	sSize=0;
	return true;
}



//////////////////////////////////////////////////////////
// Selektoren

Zahl Matrix::Set(int i,int j,const Zahl val) {
	Zahl rv=Get(i,j);
	if (i>zSize||j>sSize) return 0;
	Elements[(i-1)*sSize+j-1]=val;
	return rv;
}



//////////////////////////////////////////////////////////
// Operatoren (=, *, +, -)
/*
Matrix& Matrix::operator =(Matrix& x) {
	int i,j;
	FreeMxMem();
	AllocMxMem(x.zSize,x.sSize);
	for (i=1;i<=zSize;i++)
		for (j=1;j<=sSize;j++)
			Set(i,j,x.Get(i,j));
	return *this;
}

Matrix& Matrix::operator *(Matrix& x) {
	Matrix m=Matrix(zSize,x.sSize);
	int i,j,k;
	Zahl el;
	if(sSize!=x.zSize) return Matrix(0);

	for (i=1;i<=zSize;i++)
		for (j=1;j<=x.sSize;j++) {
			el=0;
			for (k=1;k<=sSize;k++) {
				el+=Get(i,k)*x.Get(k,j);
			}
			m.Set(i,j,el);
		}
	return m;
}

Matrix Matrix::operator +(Matrix& x) {
	if ((zSize!=x.zSize)||(sSize!=x.sSize)) return Matrix(0);
	Matrix m=Matrix(zSize, sSize);
	int i,j;
	for (i=1;i<=zSize;i++)
		for (j=1;j<=sSize;j++)
			m.Set(i,j,Get(i,j)+x.Get(i,j));
	return m;
}*/
/*
		class Matrix {
		public:
			Matrix(unsigned rows, unsigned cols);
			double& operator() (unsigned row, unsigned col);
			double	operator() (unsigned row, unsigned col) const;
			// ...
		 ~Matrix();															// Destructor
			Matrix(const Matrix& m);							 // Copy constructor
			Matrix& operator= (const Matrix& m);	 // Assignment operator
			// ...
		private:
			unsigned rows_, cols_;
			double* data_;
		};

		inline
		Matrix::Matrix(unsigned rows, unsigned cols)
			: rows_ (rows),
				cols_ (cols),
				data_ (new double[rows * cols])
		{
			if (rows == 0 || cols == 0)
				throw BadIndex("Matrix constructor has 0 size");
		}

		inline
		Matrix::~Matrix()
		{
			Free(data_);
		}

		inline
		double& Matrix::operator() (unsigned row, unsigned col)
		{
			if (row >= rows_ || col >= cols_)
				throw BadIndex("Matrix subscript out of bounds");
			return data_[cols_*row + col];
		}

		inline
		double Matrix::operator() (unsigned row, unsigned col) const
		{
			if (row >= rows_ || col >= cols_)
				throw BadIndex("const Matrix subscript out of bounds");
			return data_[cols_*row + col];
		}

*/
/*Matrix Matrix::operator -(Matrix& x)
{
	if ((zSize!=x.zSize)||(sSize!=x.sSize)) return Matrix(0);
	Matrix m=Matrix(zSize, sSize);
	int i,j;
	for (i=1;i<=zSize;i++)
		for (j=1;j<=sSize;j++)
			m.Set(i,j,Get(i,j)-x.Get(i,j));
	return m;
}

Matrix Matrix::operator *(Zahl x)
{
	Matrix m=Matrix(zSize, sSize);
	int i,j;
	for (i=1;i<=zSize;i++)
		for (j=1;j<=sSize;j++)
			m.Set(i,j,Get(i,j)*x);
	return m;
}*/



//////////////////////////////////////////////////////////
// Rearrangements

bool Matrix::TauscheZeilen(int m, int n) {
	int k;
	Zahl x;
	for (k=1; k<=sSize; k++) {
		x=Get(n,k);
		Set(n,k,Get(m,k));
		Set(m,k,x);
	}
	return true;
}

bool Matrix::TauscheSpalten(int n, int m) {
	int k;
	Zahl x;
	for (k=1; k<=zSize; k++) {
		x=Get(k,n);
		Set(k,n,Get(k,m));
		Set(k,m,x);
	}
	return true;
}



//////////////////////////////////////////////////////////
// Normen
/*
Zahl Matrix::Schurnorm() {
	int n, m;
	Zahl N=0;
	for (n=1; n<=zSize; n++)
		for (m=1; m<=sSize; m++)
			N+=(Get(n,m)*Get(n,m));
	return sqrt(N);
}

Zahl Matrix::ZSNorm() {
	int n, m;
	Zahl N=0;
	Zahl N1=0;

	for (n=1; n<=zSize; n++) {
		for (m=1; m<=sSize; m++)
			N+=fabs(Get(n,m));
		N1=max(N, N1);
		N=0;
	}
	return N1;
}

Zahl Matrix::SSNorm() {
	int n, m;
	Zahl N=0;
	Zahl N1=0;

	for (n=1; n<=sSize; n++) {
		for (m=1; m<=zSize; m++)
			N+=fabsl(Get(m,n));
		N1=max(N, N1);
		N=0;
	}
	return N1;
}
*/


//////////////////////////////////////////////////////////
// Determinante
Zahl Matrix::Det(void) {
	if (!Quadratic()) return 0;
	Matrix m(*this);
	Zahl d=1;
	int i,j,k,p=0;

	// LR-Zerlegung

	for (k=1;k<=zSize-1;k++) {
		p=m.GetBestPivot(k);
		if (p==-1) return 0;
		if (p!=k) {
			d*=(-1);
			m.TauscheZeilen(k,p);
		}
		d*=m.Get(k, k);
		for(i=k+1;i<=zSize;i++) {
			m.Set(i,k,m.Get(i,k)/m.Get(k,k));
			for(j=k+1;j<=zSize;j++)
				m.Set(i,j, m.Get(i,j)-m.Get(i,k)*m.Get(k,j));
		}
	}
	d*=m.Get(zSize,zSize);

	return d;
}

int Matrix::GetBestPivot(int Col) {
	Zahl s,max=0,q=0;
	int i,j,p=0;

	for(i=Col;i<=zSize;i++) 	{
		s=0;
		for(j=Col;j<=sSize;j++)
			s+=fabsl(Get(i,j));
		if(s==0) return -1;
		q=fabsl(Get(i,Col))/s;
		if (q>max) {
			max=q;
			p=i;
		}
	}
	if (max==0) return -1;
	return p;
}



//////////////////////////////////////////////////////////
// Datei Eingabe/Ausgabe

/*bool Matrix::Open(char* filename) {
	ifstream stream(filename);
	int zsize=0,ssize=0, i,j;
	char Buffer[255];
	Zahl el;

	if (!stream) return false;
	i=fscanf(stream,"Matrix %i %i ", &zsize,&ssize);
	if (i==1) ssize=zsize;
	if (i<1) return false;
	if ((zsize<=0)||(ssize<=0)) return false;

	FreeMxMem();
	AllocMxMem(zsize,ssize);

	fscanf(stream," { ");
	fscanf(stream," {{ ");
	for (i=1;i<=zsize;i++) {
		fscanf(stream, " { ");
		fscanf(stream, " ");

		for (j=1;j<=ssize;j++) {
			if (fscanf(stream," %Lg ", &el)==1) Set(i,j, el); else Set(i,j,0);
			fscanf(stream," , ");
		}
		fscanf(stream," } ");
		fscanf(stream," , ");
	}
	fclose(stream);
	return true;
}

bool Matrix::Write(ofstream *str) {
	int sz=GetSize(),i,j;
	str<<"{";
	for(i=1;i<=sz;i++)	 {
		if (i!=1) str<<" ";
		str<<"{";
		for (j=1;j<=sz;j++) {
			str<<Get(i,j);
			if (j!=sz) str<<",";
		}
		str<<"}";
		if (i==sz) str<<"}"; else str<<",";
		str<<"\n"<<flush;
	}
	return true;
}



//////////////////////////////////////////////////////////
// spezielle Matrizentypen

bool Matrix::Kronecker(float angle) {
	if (!Quadratic()) return false;

	switch(zSize) {
		case 1:
			Set(1,1,1);
			return true;
		case 2:
			Set(1,1,cos(angle*Pi/180));
			Set(1,2,sin(angle*Pi/180));
			Set(2,1,Get(1,2));
			Set(2,2,-Get(1,1));
			return true;
		default:
			float lg2=pow(2, (int)(log(zSize)/log(2)));
			if (lg2!=zSize) return false;
			Matrix m(zSize/2);
			if (!m.Kronecker(angle)) return false;
			return SetMxFromTmx(m*cos(angle*Pi/180),m*sin(angle*Pi/180),m*sin(angle*Pi/180),m*(-cos(angle*Pi/180)));
	}
}

bool Matrix::WalshHardamard(void) {
	if (!Quadratic()) return false;

	if (zSize==1) {
		Set(1,1,1);
		return true;
	}
	float lg2=pow(2, (int)(log(zSize)/log(2)));
	if (lg2!=zSize) return false;
	Matrix m(zSize/2);
	if (!m.WalshHardamard()) return false;
	return SetMxFromTmx(m,m,m,m*(-1));
}

bool Matrix::Zielke(Zahl a) {
	if (!Quadratic()) return false;
	int n,m;
	for (m=1;m<=zSize;m++) {
		for (n=1;n<=zSize-m;n++)
			Set(m,n,a+1);
		for (n=zSize-m+1;n<=zSize;n++)
			Set(m,n,a);
	}
	Set(zSize, zSize, a-1);
	return true;
}

bool Matrix::Hilbert(void) {
	int i,j;
	for (i=1;i<=zSize;i++)
		for (j=1;j<=zSize;j++)
			Set(i,j,1/(Zahl)(i+j-1));
	return (zSize>0);
}



//////////////////////////////////////////////////////////
// Verfahren zur Inversion

bool Matrix::GaussJordan(void) {
	if (!Quadratic()) return false;

	int p[maxElem];
	int i,j,k;
	Zahl piv=1;

	for (k=1; k<=zSize; k++) {
		p[k]=GetBestPivot(k);
		if (p[k]==-1) return false;
		if (p[k]!=k) TauscheZeilen(k,p[k]);

		piv=Get(k,k);
		if (piv==0) return false;

		for (j=1; j<=zSize; j++) {					//Spalten durchlaufen
			if(j!=k) {
				Set(k,j,Get(k,j)/(-piv));								 //Elemente der Pivotzeile
			 	for (i=1; i<=zSize; i++)									//Spalten durchlaufen
					if (i!=k) Set(i,j,Get(i,j)+Get(i,k)*Get(k,j));	//andere Elemente
			}
		}

		for (i=1; i<=zSize; i++)
			Set(i,k,Get(i,k)/piv);											//Elemente der Pivotspalte
		Set(k,k,1/piv);															 //Pivotelement
	}

	// Spaltenvertauschungen
	for (k=zSize-1; k>=1;k--) {
		if (p[k]!=k) TauscheSpalten(k,p[k]);
	}

	return true;
}

bool Matrix::Strassen(void) {
	if (!Quadratic()) return false;

	Matrix M[3], B[2][2];
	Zahl x=0;
	Zahl y=0;

//	if(Det()==0) return false;

	switch (zSize) {
		case 1:
				if (Get(1,1)==0) return false;
				Set(1,1,1/Get(1,1));
				return true;
		case 2:
				x=Get(1,1)*Get(2,2)-Get(1,2)*Get(2,1);
				if (x==0) return false;
				y=Get(2,2)/x;
				Set(2,2,Get(1,1)/x);
				Set(1,1,y);
				Set(1,2,Get(1,2)/(-x));
				Set(2,1,Get(2,1)/(-x));
				return true;
		default:
				M[0]=GetTmx(LeftTop);
				if (!M[0].Strassen()) return false;
				M[1]=M[0]*GetTmx(RT);
				M[2]=GetTmx(LB)*M[1]-GetTmx(RightBottom);
				if (!M[2].Strassen()) return false;
				B[0][1]=M[1]*M[2];
				B[1][0]=M[2]*GetTmx(LB)*M[0];
				B[0][0]=M[0]-M[1]*B[1][0];
				B[1][1]=(M[2]*(-1));
				return SetMxFromTmx(B[0][0], B[0][1], B[1][0], B[1][1]);
	}
}

bool Matrix::Schulz(Matrix& inv,char *fn) {
	if (!Quadratic()) return false;

//	Matrix inv(*this);
	Matrix res;
	int n=0;
	int i,j;
	Zahl oldres,newres;
	int higher=0;


	res=Matrix(zSize,sSize)-(*this)*inv;
	if (res.zSize==0) return false;
	newres=res.Schurnorm();

ofstream stream(fn);
stream<<"Folgende Matrix wurde mit dem Schulz-Algorithmus invertiert:\n";
Write(stream);
stream<<"\nDie Schurnorm des Residuums nach den einzelnen Iterationsschritten war folgendes:\nRes={";

	if ( (res.SSNorm()>=1) && (res.ZSNorm()>=1) && (res.Schurnorm()>=1)) return false;


	while ((newres>Abbruch) && (n<maxRek) && (higher<10) ) {
stream<<res.Schurnorm()<<",";
		n++;
		inv=inv*(Matrix(zSize)+res);
//		inv=inv*(Matrix(zSize)*2-(*this)*inv);
//		res=res*res;							// sollte (theoretisch!!!) funktionieren statt nächster Zeile
		res=Matrix(zSize)-(*this)*inv;
		oldres=newres;
		newres=res.Schurnorm();
		if (newres>oldres) higher++;
	}
stream<<res.Schurnorm()<<"}\n\n\nDie Inverse:\n";
inv.Write(stream);
stream<<"\n\nUnd die Residuenmatrix:\n";
res.Write(stream);
stream.close();

	for (i=1;i<=zSize;i++)
		for (j=1;j<=sSize;j++)
			Set(i,j,inv.Get(i,j));

	return true;
//	return ((n<maxRek)||(res.Schurnorm()<=0.0001));
}

bool Matrix::LRInversion(void) {
	if (!Quadratic()) return false;

	int p[maxElem];
	int k, i, j, z;
	Zahl h;
	Matrix inv(zSize, zSize);

	// LR-Zerlegung
	for (k=1;k<=zSize-1;k++) {
		p[k-1]=GetBestPivot(k);

		if (p[k-1]==-1) return false;

		if(p[k-1]!=k) TauscheZeilen(k,p[k-1]);

		for(i=k+1;i<=zSize;i++) {
			if (Get(k,k)==0) return false;
			Set(i,k, Get(i,k)/Get(k,k));
			for(j=k+1;j<=zSize;j++)
				Set(i,j,Get(i,j)-Get(i,k)*Get(k,j));
		}
	}

		// Vorwärtseinsetzen
	for(k=1;k<zSize;k++) {		 //Zeilenvertauschungen, lieáen sich noch vereinfachen
		inv.TauscheZeilen(p[k-1],k);
	}


	for(z=1;z<=zSize;z++) {
		for (i=1;i<=zSize;i++)
			for(j=1;j<i;j++)

				inv.Set(i,z,inv.Get(i,z)-Get(i,j)*inv.Get(j,z));

		// Rückwärtseinsetzen

		for(i=zSize;i>=1;i--) {
			for(k=i+1;k<=zSize;k++)
				inv.Set(i,z, inv.Get(i,z)+Get(i,k)*inv.Get(k,z));
			if (Get(i,i)==0) return false;
			inv.Set(i,z, inv.Get(i,z)/Get(i,i)*(-1));
		}

	}
	for (i=1;i<=zSize;i++)
		for (j=1;j<=zSize;j++)
			Set(i,j,-inv.Get(i,j));

	return true;
}


//////////////////////////////////////////////////////////
// Teilmatrizen

Matrix Matrix::GetTmx(int flag) {
	int z=1,s=1,ze=(zSize+1)/2,se=(sSize+1)/2,i,j;
	switch(flag) {
		case RT:
			s=(sSize+1)/2+1;
			se=sSize;
			break;
		case RB:
			s=(sSize+1)/2+1;
			se=sSize;
		case LB:
			z=(zSize+1)/2+1;
			ze=zSize;
	}
	Matrix Tmx(ze-z+1,se-s+1);
	for (i=z;i<=ze;i++)
		for(j=s;j<=se;j++)
			Tmx.Set(i-z+1,j-s+1,Get(i,j));

	return Tmx;
}

bool Matrix::SetMxFromTmx(Matrix& B11, Matrix& B12, Matrix& B21, Matrix& B22) {
	int i,j;

	if ((B11.zSize!=B12.zSize)||(B11.sSize!=B21.sSize)||(B21.zSize!=B22.zSize)||(B12.sSize!=B22.sSize)) return false;

	FreeMxMem();
	AllocMxMem(B11.zSize+B22.zSize,B11.sSize+B22.sSize);

	for (i=1;i<=B11.zSize;i++)
		for (j=1;j<=B11.sSize;j++)
			Set(i,j,B11.Get(i,j));

	for (i=1;i<=B12.zSize;i++)
		for (j=1;j<=B12.sSize;j++)
			Set(i,j+B11.sSize,B12.Get(i,j));

	for (i=1;i<=B21.zSize;i++)
		for (j=1;j<=B21.sSize;j++)
			Set(i+B11.zSize,j,B21.Get(i,j));

	for (i=1;i<=B22.zSize;i++)
		for (j=1;j<=B22.sSize;j++)
			Set(i+B11.zSize,j+B11.sSize,B22.Get(i,j));

	return true;
}

*/
