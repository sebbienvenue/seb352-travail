/***************************************************************************
						Distributions.h  -  description
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


#ifndef DISTRIBUTIONS_H
#define DISTRIBUTIONS_H
#include <math.h>
#include <string>
#include "BaseDefinitions.h"

#define RANGE_NONE				-1
#define RANGE_MININF_INF		0
#define RANGE_0_INF				1
#define RANGE_0_1				2
#define RANGE_PARA1_PARA2	3
#define RANGE_0_PARA1			4
#define RANGE_PARA1_INF		5
#define RANGE_PARA2_INF		6
#define RANGE_PARA3_INF		7
#define RANGE_PARA1P1_INF	8

#define INFINITY 10e+300
#define inf INFINITY
#define MOMENT_NONE INFINITY


#define USE_EXP_DIST				1
#define USE_GAMMA_DIST			2
#define USE_NORMAL_DIST			3
#define USE_LOGNORMAL_DIST	4
#define USE_PARETO0_DIST			5
#define USE_WEIBULL_DIST			6


// Mathematical function to be defined somewhere else...
//Zahl IncompleteBeta( Zahl p1, Zahl p2, Zahl p3) {return 0;}//TODO
//Zahl IncompleteGamma(Zahl p1,Zahl p2) {return 0;} //TODO
#define IncompleteBeta( p1, p2, p3)	0
#define IncompleteGamma(p1,p2) 0
#define tgamma(x) exp(lgamma(x))




/* =================================================================================== */




	/** This is the base class for all distributions
		*@author Reinhold Kainhofer
	*/

class Distribution {
protected:
	Zahl param[4];
	int paramnr;
public:
	Distribution(int paramnumber=0, Zahl para=0, Zahl para1=0){param[0]=para; param[1]=para1; paramnr=paramnumber;};
	virtual ~Distribution(){};
		/**	Returns a string specifying the name of this distribution
		*/
	virtual string DistributionName(){return "Generic";};
	virtual string DistributionShortName(){ return DistributionName().substr(0,3);}
		/**	returns the value f(y) of the distribution density at y
		*/
	virtual Zahl pdf(Zahl y) {return pdf(y, param[0]);};
	virtual Zahl pdf(Zahl y, Zahl para) {return pdf(y, para, param[1]);};
	virtual Zahl pdf(Zahl y, Zahl para, Zahl para1)=0;
		/**	transforms the uniform variate y to a variate distributed with this distribution
		*/
	virtual Zahl Transform(Zahl y) {return Transform(y,param[0]);};	// transforms uniformly distributed variate into a variable distributed according to this distribution
	virtual Zahl Transform(Zahl y, Zahl para) {return Transform(y,para, param[1]);};	// transforms uniformly distributed variate into a variable distributed according to this distribution
	virtual Zahl Transform(Zahl y, Zahl para, Zahl para1) =0;
		/**	transforms the uniform variate y to a variate distributed with this distribution, but only between 0 and cut, so the density needs to be rescaled
		*/
	virtual inline Zahl TransformCut(Zahl y, Zahl cut) {return TransformCut(y, cut, param[0]);};
	virtual inline Zahl TransformCut(Zahl y, Zahl cut, Zahl para) {return TransformCut(y, cut, para, param[1]);};
	virtual inline Zahl TransformCut(Zahl y, Zahl cut, Zahl para, Zahl para1) {Zahl res=Transform(y*cdf(cut, para, para1), para, para1); if ((res>cut) && (res-cut<0.005)) res=cut; return res;};
		/**	Sets the parameter for the distribution
		*/
	virtual inline void SetParameter(Zahl para, int nr=0) {if (0<=nr && nr<paramnr) param[nr]=para;};
		/**	Returns the parameter of this distribution
		*/
	virtual inline Zahl GetParameter(int nr=0){if (0<=nr && nr<paramnr) return param[nr]; else return 0;};
		/**	returns the distribution function F(y)
		*/
	virtual Zahl cdf(Zahl y) {return cdf(y,param[0]);};
	virtual Zahl cdf(Zahl y, Zahl para) {return cdf(y,para, param[1]);};
	virtual Zahl cdf(Zahl y, Zahl para, Zahl para1) =0;
		/**	returns 1-F(y)
		*/
	virtual Zahl cdfc(Zahl y) {return cdfc(y, param[0]);}
	virtual Zahl cdfc(Zahl y, Zahl para) {return cdfc(y, para, param[1]);}
	virtual Zahl cdfc(Zahl y, Zahl para, Zahl para1) {return 1-cdf(y, para, para1);}

	virtual int Range()=0;
	virtual int ParameterRange(int paranr) =0;
	virtual Zahl mgf(Zahl y, Zahl para, Zahl para1)=0;
	virtual Zahl Moment(int k, Zahl para, Zahl para1)=0;
	virtual Zahl Mean() {return Mean(param[0], param[1]);};
	virtual Zahl Mean(Zahl para, Zahl para1) {return Moment(1,para,para1);};
	virtual Zahl Variance(Zahl para, Zahl para1) {return Moment(2,para,para1);};
	virtual Zahl StdDeviation(Zahl para, Zahl para1) {Zahl vari=Variance(para,para1); if (vari==MOMENT_NONE) return MOMENT_NONE; else return sqrt(vari);}
	virtual Zahl Skewness(Zahl para, Zahl para1) {return Moment(3,para,para1);};
	virtual Zahl Kurtosis(Zahl para, Zahl para1) {return Moment(4,para,para1);};
	virtual Zahl LogLikelihood(Zahl y, Zahl para, Zahl para1)=0;
	virtual Zahl CharacteristicFunction(Zahl y, Zahl para, Zahl para1)=0;
	virtual Zahl Median(Zahl para, Zahl para1)=0;
	virtual Zahl Q1(Zahl para, Zahl para1)=0;
	virtual Zahl Q3(Zahl para, Zahl para1)=0;
	virtual Zahl Mode(Zahl para, Zahl para1)=0;
	virtual Zahl qMean(Zahl para, Zahl para1)=0;
	virtual Zahl qMode(Zahl para, Zahl para1)=0;
};
//TODO: Discrete distributions
class DistDiscrete : public Distribution {
};


class DistTemplate : public Distribution {
public:
	DistTemplate(Zahl para=0, Zahl para1=0) : Distribution(0,para, para1){}
	virtual ~DistTemplate(){};
	virtual string DistributionName(){return "Template";};
	virtual Zahl pdf(Zahl y, Zahl para, Zahl para1) {return 0;} //TODO
	virtual Zahl Transform(Zahl y, Zahl para, Zahl para1)	{return 0;} //TODO
	virtual inline Zahl TransformCut(Zahl y, Zahl cut, Zahl para, Zahl para1) {return Transform(y*cdf(cut, para, para1), para, para1);};
	virtual Zahl cdf(Zahl y, Zahl para, Zahl para1)	{return 0;} //TODO
	virtual Zahl cdfc(Zahl y, Zahl para, Zahl para1) {return 1-cdf(y, para, para1);}

	virtual int Range() {return 0;} //TODO
	virtual int ParameterRange(int paranr) {return 0;} //TODO;
	virtual Zahl mgf(Zahl y, Zahl para, Zahl para1) {return 0;} //TODO
	virtual Zahl Moment(int k, Zahl para, Zahl para1) {return 0;} //TODO
	virtual Zahl Mean(Zahl para, Zahl para1) {return Moment(1,para,para1);};
	virtual Zahl Variance(Zahl para, Zahl para1) {return Moment(2,para,para1);};
	virtual Zahl Skewness(Zahl para, Zahl para1) {return Moment(3,para,para1);};
	virtual Zahl Kurtosis(Zahl para, Zahl para1) {return Moment(4,para,para1);};
	virtual Zahl LogLikelihood(Zahl y, Zahl para, Zahl para1) {return 0;} //TODO
	virtual Zahl CharacteristicFunction(Zahl y, Zahl para, Zahl para1) {return 0;} //TODO
	virtual Zahl Median(Zahl para, Zahl para1) {return 0;} //TODO
	virtual Zahl Q1(Zahl para, Zahl para1) {return 0;} //TODO
	virtual Zahl Q3(Zahl para, Zahl para1) {return 0;} //TODO
	virtual Zahl Mode(Zahl para, Zahl para1) {return 0;} //TODO
	virtual Zahl qMean(Zahl para, Zahl para1) {return 0;} //TODO
	virtual Zahl qMode(Zahl para, Zahl para1) {return 0;} //TODO
};


class DistBeta : public Distribution {
public:
	DistBeta(Zahl para=0, Zahl para1=0, Zahl para2=0, Zahl para3=0) : Distribution(4,para, para1){ param[2]=para2; param[3]=para3;}
	virtual ~DistBeta(){};
	virtual string DistributionName(){return "Beta";}
	virtual Zahl pdf(Zahl y, Zahl para, Zahl para1) {return tgamma(param[2]+param[3])/(tgamma(param[2])*tgamma(param[3])*pow(param[1]-param[0], param[2]+param[3]-1)) *pow(y-param[0], param[2]-1)*pow(param[1]-y, param[3]-1);}
	virtual Zahl Transform(Zahl y, Zahl para, Zahl para1)	{return 0;} //TODO numerically???
	virtual Zahl cdf(Zahl y, Zahl para, Zahl para1)	{return IncompleteBeta( (y-para)/(para1-para), param[2], param[3]);}

	virtual int Range() {return RANGE_PARA1_PARA2;}
	virtual int ParameterRange(int paranr)	{ switch (paranr) {
		case 1: case 2: return RANGE_MININF_INF;
		case 3: case 4: return RANGE_0_INF;
		default: return RANGE_NONE;} }
	virtual Zahl mgf(Zahl y, Zahl para, Zahl para1) {return 0;} //TODO
	virtual Zahl Moment(int k, Zahl para, Zahl para1) {return 0;} //TODO numerically
	virtual Zahl Mean(Zahl para, Zahl para1) {return (para*param[3] + para1*param[2])/(param[2]+param[3]);}
	virtual Zahl Variance(Zahl para, Zahl para1) {return (param[2]*param[3]*(para1-para)*(para1-para) ) / ( (param[2]+param[3]+1)*(param[2]+param[3])*(param[2]+param[3]) );}
	virtual Zahl StdDeviation(Zahl para, Zahl para1) {return Moment(2, para,para1);}//TODO
	virtual Zahl Skewness(Zahl para, Zahl para1) {
				 Zahl C=param[2],D=param[3];
				 return (2*C*D*(D-C))	/	( pow(C+D, 3) * (C+D+1) * (C+D+2) * sqrt(pow((C*D)/( (C+D)*(C+D)*(C+D+1)), 3)) ); }
	virtual Zahl Kurtosis(Zahl para, Zahl para1) {
				 Zahl C=param[2],D=param[3];
				 return 3*( ( (C*C*(D+2) + 3*D*D + C*D*(D-2))*(C+D+1)) / (C*D*(C+D+2)*(C+D+3) ) -1);}
	virtual Zahl LogLikelihood(Zahl y, Zahl para, Zahl para1) {return 0;} //TODO
	virtual Zahl CharacteristicFunction(Zahl y, Zahl para, Zahl para1) {return 0;} //TODO
	virtual Zahl Median(Zahl para, Zahl para1) {return 0;} //TODO no simple closed form
	virtual Zahl Q1(Zahl para, Zahl para1) {return 0;} //TODO no simple closed form
	virtual Zahl Q3(Zahl para, Zahl para1) {return 0;} //TODO no simple closed form
	virtual Zahl Mode(Zahl para, Zahl para1) {Zahl A=para, B=para1, C=param[2],D=param[3];
				 if (C==1 && D==1) return -inf;
				 else return (A*(D-1) + B*(C-1))/(C+D-2); }
	virtual Zahl qMean(Zahl para, Zahl para1) {return 0;} //TODO no simple closed form
	virtual Zahl qMode(Zahl para, Zahl para1) {return 0;} //TODO no simple closed form
};

// Beta(0,1,C,1) is often called the Power-function distribution.
//class DistPowerFunction : public DistBeta {

//TODO: Discrete distributions
//class DistBinomial : public DistDiscrete {
//	DistDiscrete(Zahl para=0.5,

// The Bradford distribution has been used to model the distribution of references among several sources
class DistBradford : public Distribution {
public:
	DistBradford(Zahl para=0, Zahl para1=1, Zahl para2=5) : Distribution(3,para, para1){param[2]=para2;}
	virtual ~DistBradford(){}
	virtual string DistributionName(){return "Bradford";}
	virtual Zahl pdf(Zahl y, Zahl para, Zahl para1) { return param[2]/( (param[2]*(y-para) + para1 - para) * log(param[2]+1) );}
	virtual Zahl Transform(Zahl y, Zahl para, Zahl para1) { return (para*(param[2]+1) - para1 + (para1-para)*pow(param[2]+1, y) ) / param[2];}
//	virtual inline Zahl TransformCut(Zahl y, Zahl cut, Zahl para, Zahl para1) {return Transform(y*cdf(cut, para, para1), para, para1);};
	virtual Zahl cdf(Zahl y, Zahl para, Zahl para1) { return log(1+ (param[2]*(y-para))/(para1-para) )	/ log(param[2]+1);}

	virtual int Range() {return RANGE_PARA1_PARA2;};
	virtual int ParameterRange(int paranr)	{switch (paranr) {
				 case 1: case 2: return RANGE_MININF_INF;
				 case 3: return RANGE_0_INF;
				 default: return RANGE_NONE;} }
	virtual Zahl mgf(Zahl y, Zahl para, Zahl para1) {return 0;} //TODO, nicht angegegben
	virtual Zahl Moment(int k, Zahl para, Zahl para1) {return 0;} //TODO, nicht angegeben
	virtual Zahl Mean(Zahl para, Zahl para1) {
				Zahl A=para, B=para1, C=param[2], k;
				 k=log(C+1);
				 return (C*(B-A) + k*( A*C+A-B) ) / (C*k);}
	virtual Zahl Variance(Zahl para, Zahl para1) {
				 Zahl A=para, B=para1, C=param[2], k;
				 k=log(C+1);
				 return (B-A)*(B-A) *( C*(k-2) + 2*k)	/ (2*C*k*k);}
	virtual Zahl StdDeviation(Zahl para, Zahl para1) {return Moment(2, para,para1);};//TODO
	virtual Zahl Skewness(Zahl para, Zahl para1) {
				Zahl C=param[2], k;
				 k=log(C+1);
				 return M_SQRT2* (12*C*C - 9*k*C*(C+2) + 2*k*k*(C*(C+3)+3) ) / (sqrt(C*(C*(k-2)+2*k))*(3*C*(k-2)+6*k) );}
	virtual Zahl Kurtosis(Zahl para, Zahl para1) {
				 Zahl C=param[2], k;
				 k=log(C+1);
				 return (C*C*C*(k-3)*(k*(3*k-16)+24) + 12*k*C*C*(k-4)*(k-3) + 6*C*k*k*(3*k-14) + 12*k*k*k) / (3*C*pow(C*(k-2) + 2*k, 2) );}
	virtual Zahl LogLikelihood(Zahl y, Zahl para, Zahl para1) {return 0;}// TODO, nicht angegeben
	virtual Zahl CharacteristicFunction(Zahl y, Zahl para, Zahl para1) {return 0;} //TODO, nicht angegeben
	virtual Zahl Median(Zahl para, Zahl para1) {return (para*(param[2]+1) - para1 + (para1-para)*sqrt(param[2]+1)) / param[2];}
	virtual Zahl Q1(Zahl para, Zahl para1) { return (para*(param[2]+1)-para1+(para1-para)*sqrt(sqrt(param[2]+1)) ) / param[2];}
	virtual Zahl Q3(Zahl para, Zahl para1) { return (para*(param[2]+1)-para1+(para1-para)*sqrt(sqrt(pow(param[2]+1, 3) )) ) / param[2];}
	virtual Zahl Mode(Zahl para, Zahl para1) {return para;}
	virtual Zahl qMean(Zahl para, Zahl para1) {Zahl C=param[2]; return log(C/log(C+1)) / log(C+1);}
	virtual Zahl qMode(Zahl para, Zahl para1) {return 0;}
};


// class DistBurr //TODO

// U,V ~ Normal(0,1)	-> U/V~Cauchy(0,1)
// Z~Cauchy -> W=1/(a+b*Z) ~Cauchy
// The Cauchy distribution is sometimes called the Lorentz distribution

class DistCauchy : public Distribution {
public:
	DistCauchy(Zahl para=0, Zahl para1=1) : Distribution(2,para, para1){}
	virtual ~DistCauchy(){}
	virtual string DistributionName(){return "Cauchy";};
	virtual Zahl pdf(Zahl y, Zahl para, Zahl para1) {return 1/(M_PI*para1*(1+pow((y-para)/para1, 2)) ); };
	virtual Zahl Transform(Zahl y, Zahl para, Zahl para1)	{return para + para*tan(M_PI*(y-0.5));};
//	virtual inline Zahl TransformCut(Zahl y, Zahl cut, Zahl para, Zahl para1) {return Transform(y*cdf(cut, para, para1), para, para1);};
	virtual Zahl cdf(Zahl y, Zahl para, Zahl para1)	{return 0.5+1/M_PI*atan( (y-para)/para1 );};

	virtual int Range() { return RANGE_MININF_INF;};
	virtual int ParameterRange(int paranr)	{ if (paranr==0) return RANGE_MININF_INF; if (paranr==1) return RANGE_0_INF; return RANGE_NONE;};
	virtual Zahl mgf(Zahl y, Zahl para, Zahl para1) {return MOMENT_NONE;}
	virtual Zahl Moment(int k, Zahl para, Zahl para1) {return MOMENT_NONE;};
	virtual Zahl LogLikelihood(Zahl y, Zahl para, Zahl para1) {return 0;} //TODO, nicht angegeben
	virtual Zahl CharacteristicFunction(Zahl y, Zahl para, Zahl para1) {return 0;} // TODO: nicht angegeben
	virtual Zahl Median(Zahl para, Zahl para1) {return para;};
	virtual Zahl Q1(Zahl para, Zahl para1) {return para-para1;};
	virtual Zahl Q3(Zahl para, Zahl para1) {return para+para1;};
	virtual Zahl Mode(Zahl para, Zahl para1) {return para;};
	virtual Zahl qMean(Zahl para, Zahl para1) {return MOMENT_NONE;};
	virtual Zahl qMode(Zahl para, Zahl para1) {return 0.5;};
};
#define DistLorentz DistChauchy

/* =================================================================================== */






/** This class describes the exponential distribution
	*@author Reinhold Kainhofer
	*/
class DistExp : public Distribution {
public:
	/** Initialize the distribution, if no parameter is given, use the natural assumption Lambda=1,
			para is lambda, para1 is the shift, so this class also describes the Exponential Tail distribution.
			Note that in [McL] B=1/lambda=1/para, also A is para1, while 1/B is para here, so their order is switched, too
	*/
	DistExp(Zahl lamb=1, Zahl para=0):Distribution(2, lamb, 0) {}
	virtual ~DistExp() {}
	virtual string DistributionName(){return "Exponential";};
	/** Returns the density of the exponential distribution evaluated at y:
			Lambda*exp(-Lambda*y)
	*/
	virtual Zahl pdf(Zahl y, Zahl para, Zahl para1) { return para*exp(para*(para1-y) );}
	/** Transforms the uniform variate y into an exponentially transformed one with parameter para, which can be different from the parameter stored inside the class
	*/
	virtual Zahl Transform(Zahl y, Zahl para, Zahl para1) {return para1-log(1-y)/para;};
	virtual Zahl cdf(Zahl y, Zahl para, Zahl para1) {return 1-exp(para*(para1-y));}
	virtual Zahl cdfc(Zahl y, Zahl para, Zahl para1) {return exp(para*(para1-y) );}
	virtual int Range() {return RANGE_PARA2_INF;}
	virtual int ParameterRange(int paranr) {if (paranr==1) return RANGE_0_INF; if (paranr==2) return RANGE_MININF_INF; return RANGE_NONE;}
	virtual Zahl mgf(Zahl y, Zahl para, Zahl para1) {return 0;} //TODO
	virtual Zahl Moment(int k, Zahl para, Zahl para1) {return 0;} //TODO
	virtual Zahl Mean(Zahl para, Zahl para1) {return 1/para+para1;}
	virtual Zahl Variance(Zahl para, Zahl para1) {return 1/(para*para);}
	virtual Zahl Skewness(Zahl para, Zahl para1) {return 2;}
	virtual Zahl Kurtosis(Zahl para, Zahl para1) {return 6;}
	virtual Zahl LogLikelihood(Zahl y, Zahl para, Zahl para1) {return 0;} //TODO
	virtual Zahl CharacteristicFunction(Zahl y, Zahl para, Zahl para1) {return 0;} //TODO
	virtual Zahl Median(Zahl para, Zahl para1) {return para1+ log(2)/para;}
	virtual Zahl Q1(Zahl para, Zahl para1) {return para1+log(4/3)/para;}
	virtual Zahl Q3(Zahl para, Zahl para1) {return para1+log(4)/para;}
	virtual Zahl Mode(Zahl para, Zahl para1) {return para1;}
	virtual Zahl qMean(Zahl para, Zahl para1) {return (M_E-1)/M_E;}
	virtual Zahl qMode(Zahl para, Zahl para1) {return 0;}
};




// Gamma(A,B,C) where C is an integer, is the Erlang distributions
// Gamma(A,B,1) is the exponential distribution
// Gamma(0,2,nu/2) is the Chi-square distribution with nu degrees of freedom

// if Z_1 ~Gamma(A,B,C_1) and Z_2~Gamma(A,B,C_2)	=> (Z_1+Z_2)~Gamma(A,B,C_1+C_2)
// Z_1, ... Z_nu ~Normal(0,1)	=> W=sum_{k=1}^{nu} Z_k^2 ~Gamma(0,2,nu/2)
// Z_1, ... Z_n ~Exp(lambda, A)	=> W=sum_{k=1}^{nu} Z_k ~Erlang(A, 1/lambda, n)

/** This class describes the Gamma distribution with parameters alpha and beta
	*@author Reinhold Kainhofer
	*/
class DistGamma : public Distribution {
public:
	/** Initializes the Gamma distribution with the given parameters alpha and beta
	*/
	DistGamma(Zahl al=1, Zahl be=1, Zahl A=0):Distribution(3, al, be) {param[2]=A;}
	virtual ~DistGamma() {}
	virtual string DistributionName(){return "Gamma";}
	/** returns the density of the gamma distribution
			para1=beta=1/B in [McL]
			para=k=C in [McL]
			param[2] is the shift
	*/
	virtual Zahl pdf(Zahl y, Zahl para, Zahl para1) { return 1/(para1*tgamma(para))*pow( (y-param[2])*para1,para-1)*exp( (param[2]-y)*para1); }
	/** transforms the uniform variate y into a gamma distributed one
	*/
	virtual Zahl Transform(Zahl y, Zahl para, Zahl para1);
	/** transforms the uniform variate y into a gamma distributed one, which is only in the	interval [0,cut), so the density needs to be rescaled
	*/
	virtual Zahl TransformCut(Zahl y, Zahl cut, Zahl para, Zahl para1);
	virtual Zahl cdf(Zahl y, Zahl para, Zahl para1) {return IncompleteGamma(para, (y-param[2])*para1);}

	virtual int Range() {return RANGE_PARA3_INF;}
	virtual int ParameterRange(int paranr) {switch (paranr) {
				case 1: return RANGE_MININF_INF;
				case 2: case 3: return RANGE_0_INF;
				default: return RANGE_NONE;} }
	virtual Zahl mgf(Zahl y, Zahl para, Zahl para1) {return 0;} //TODO
	virtual Zahl Moment(int k, Zahl para, Zahl para1) {return 0;} //TODO
	virtual Zahl Mean(Zahl para, Zahl para1) {return param[2] + para/para1;}
	virtual Zahl Variance(Zahl para, Zahl para1) {return para/(para1*para1);}
	virtual Zahl Skewness(Zahl para, Zahl para1) {return 2/sqrt(para);}
	virtual Zahl Kurtosis(Zahl para, Zahl para1) {return 6/para;}
	virtual Zahl LogLikelihood(Zahl y, Zahl para, Zahl para1) {return 0;} //TODO
	virtual Zahl CharacteristicFunction(Zahl y, Zahl para, Zahl para1) {return 0;} //TODO
	virtual Zahl Median(Zahl para, Zahl para1) {return 0;} //TODO, no simple closed form
	virtual Zahl Q1(Zahl para, Zahl para1) {return 0;} //TODO, no simple closed form
	virtual Zahl Q3(Zahl para, Zahl para1) {return 0;} //TODO, no simple closed form
	virtual Zahl Mode(Zahl para, Zahl para1) {return param[2]+ (para-1)/para1;}
	virtual Zahl qMean(Zahl para, Zahl para1) {return IncompleteGamma(param[0], param[0]);}
	virtual Zahl qMode(Zahl para, Zahl para1) {return IncompleteGamma(param[0], param[0]-1);}
};




class DistLaplace : public DistExp {
public:
	DistLaplace(Zahl lamb=1, Zahl para=0):DistExp(lamb, para) {}
	virtual ~DistLaplace(){}
	virtual string DistributionName(){return "Laplace";};
	virtual Zahl pdf(Zahl y, Zahl para, Zahl para1) {return 1/2.*DistExp::pdf(fabs(y), para, para1);}
	virtual Zahl Transform(Zahl y, Zahl para, Zahl para1) { Zahl u=2*y-1; if (u<0) return -DistExp::Transform(-u,para,para1); else return DistExp::Transform(u,para,para1);}
	virtual Zahl cdf(Zahl y, Zahl para, Zahl para1) {if (y<para1) return 0.5*exp(para*(y-para1)); else return 1-0.5*exp(para*(para1-y));}

	virtual int Range() {return RANGE_MININF_INF;}
	virtual Zahl mgf(Zahl y, Zahl para, Zahl para1) {return 0;} //TODO
	virtual Zahl Moment(int k, Zahl para, Zahl para1) {return 0;} //TODO
	virtual Zahl Mean(Zahl para, Zahl para1) {return para1;}
	virtual Zahl Variance(Zahl para, Zahl para1) {return 2/(para*para);}
	virtual Zahl Skewness(Zahl para, Zahl para1) {return 0;}
	virtual Zahl Kurtosis(Zahl para, Zahl para1) {return 3;}
	virtual Zahl LogLikelihood(Zahl y, Zahl para, Zahl para1) {return 0;} //TODO
	virtual Zahl CharacteristicFunction(Zahl y, Zahl para, Zahl para1) {return 0;} //TODO
	virtual Zahl Median(Zahl para, Zahl para1) {return para1;}
	virtual Zahl Q1(Zahl para, Zahl para1) {return para1-log(2)/para;}
	virtual Zahl Q3(Zahl para, Zahl para1) {return para1+log(2)/para;}
	virtual Zahl Mode(Zahl para, Zahl para1) {return para1;}
	virtual Zahl qMean(Zahl para, Zahl para1) {return 0.5;}
	virtual Zahl qMode(Zahl para, Zahl para1) {return 0.5;}
};
#define DistDoubleExponential DistLaplace
#define DistBilateralExponential DistLaplace



/* ================================================================================== */




/** This class describes the Normal (Gauss) distribution
	*@author Reinhold Kainhofer
	*/
class DistNormal : public Distribution {
public:
	DistNormal(Zahl mm=0, Zahl ss=1) : Distribution(2, mm, ss) {}
	virtual ~DistNormal() {}
	virtual string DistributionName(){return "Normal (Gauss)";};
	virtual Zahl pdf(Zahl y, Zahl para, Zahl para1) { return M_1SQRT2PI/para1*exp(-(y-para)*(y-para)/(2*para1*para1));};
	virtual Zahl Transform(Zahl y, Zahl mm, Zahl ss);
				/** The distribution function needs to be calculated numerically. Not implemented...
				*/
	virtual Zahl cdf(Zahl y, Zahl para, Zahl para1);
	virtual Zahl cdfc(Zahl y, Zahl para, Zahl para1) {return cdf(-y, para, para1);}

	virtual int Range() {return RANGE_MININF_INF;}
	virtual int ParameterRange(int paranr) {if (paranr==1) return RANGE_MININF_INF; if (paranr==2) return RANGE_0_INF; return RANGE_NONE;}
	virtual Zahl mgf(Zahl y, Zahl para, Zahl para1) {return 0;} //TODO
	virtual Zahl Moment(int k, Zahl para, Zahl para1) {return 0;} //TODO
	virtual Zahl Mean(Zahl para, Zahl para1) {return para;}
	virtual Zahl Variance(Zahl para, Zahl para1) {return para1*para1;}
	virtual Zahl Skewness(Zahl para, Zahl para1) {return 0;}
	virtual Zahl Kurtosis(Zahl para, Zahl para1) {return 0;}
	virtual Zahl LogLikelihood(Zahl y, Zahl para, Zahl para1) {return 0;} //TODO
	virtual Zahl CharacteristicFunction(Zahl y, Zahl para, Zahl para1) {return 0;} //TODO
	virtual Zahl Median(Zahl para, Zahl para1) {return para;}
	virtual Zahl Q1(Zahl para, Zahl para1) {return para-0.6745*para1;}
	virtual Zahl Q3(Zahl para, Zahl para1) {return para+0.6745*para1;}
	virtual Zahl Mode(Zahl para, Zahl para1) {return para;}
	virtual Zahl qMean(Zahl para, Zahl para1) {return 0.5;}
	virtual Zahl qMode(Zahl para, Zahl para1) {return 0.5;}
};
#define DistGauss DistNormal

// Z~Normal(A,B)	=>	|Z| ~FoldedNormal(A,B)
class DistFoldedNormal : public DistNormal {
public:
	DistFoldedNormal(Zahl mm=0, Zahl ss=1) : DistNormal(mm, ss) {}
	virtual ~DistFoldedNormal() {}
	virtual string DistributionName(){return "Folded Normal";};
	virtual Zahl pdf(Zahl y, Zahl para, Zahl para1) { return M_SQRT2/(para1*M_SQRTPI)*cosh(para*y/(para1*para1) )*exp(-0.5*(y*y+para*para)/(para1*para1));};
	virtual Zahl Transform(Zahl y, Zahl mm, Zahl ss) {return 0;} //TODO
				/** The distribution function needs to be calculated numerically. Not implemented...
				*/
	virtual Zahl cdf(Zahl y, Zahl para, Zahl para1) { return DistNormal::cdf(y, para, para1)-DistNormal::cdf(-y, para, para1);}
	virtual Zahl cdfc(Zahl y, Zahl para, Zahl para1) {return 1-cdf(y, para, para1);}

	virtual int Range() {return RANGE_0_INF;}
	virtual int ParameterRange(int paranr) {if (paranr==1) return RANGE_MININF_INF; if (paranr==2) return RANGE_0_INF; return RANGE_NONE;}
	virtual Zahl mgf(Zahl y, Zahl para, Zahl para1) {return 0;} //TODO
	virtual Zahl Moment(int k, Zahl para, Zahl para1) {return 0;} //TODO
	virtual Zahl Mean(Zahl para, Zahl para1) {return 0;} //TODO
	virtual Zahl Variance(Zahl para, Zahl para1) {return 0;} //TODO
	virtual Zahl Skewness(Zahl para, Zahl para1) {return 0;} //TODO
	virtual Zahl Kurtosis(Zahl para, Zahl para1) {return 0;} //TODO
	virtual Zahl LogLikelihood(Zahl y, Zahl para, Zahl para1) {return 0;} //TODO
	virtual Zahl CharacteristicFunction(Zahl y, Zahl para, Zahl para1) {return 0;} //TODO
	virtual Zahl Median(Zahl para, Zahl para1) {return 0;} //TODO, no simple closed form
	virtual Zahl Q1(Zahl para, Zahl para1) {return 0;} //TODO, no simple closed form
	virtual Zahl Q3(Zahl para, Zahl para1) {return 0;} //TODO, no simple closed form
	virtual Zahl Mode(Zahl para, Zahl para1) {return 0;} //TODO, no simple closed form
	virtual Zahl qMean(Zahl para, Zahl para1) {return 0;} //TODO, no simple closed form
	virtual Zahl qMode(Zahl para, Zahl para1) {return 0;} //TODO, no simple closed form
};

// Special case of Chi and Folded Normal distribution
class DistHalfNormal : public DistNormal {
public:
	DistHalfNormal(Zahl mm=0, Zahl ss=1) : DistNormal(mm, ss) {}
	virtual ~DistHalfNormal() {}
	virtual string DistributionName(){return "Half Normal";};
	virtual Zahl pdf(Zahl y, Zahl para, Zahl para1) {return 2*DistNormal::pdf(y,para,para1);}
	virtual Zahl Transform(Zahl y, Zahl mm, Zahl ss) {return DistNormal::Transform(y/2.+0.5, mm, ss);}
				/** The distribution function needs to be calculated numerically. Not implemented...
				*/
	virtual Zahl cdf(Zahl y, Zahl para, Zahl para1) { return 2*DistNormal::cdf(y,para,para1) - 1;}
	virtual Zahl cdfc(Zahl y, Zahl para, Zahl para1) {return 2*DistNormal::cdf(y,para,para1);}

	virtual int Range() {return RANGE_PARA1_INF;}
	virtual Zahl mgf(Zahl y, Zahl para, Zahl para1) {return 0;} //TODO
	virtual Zahl Moment(int k, Zahl para, Zahl para1) {return 0;} //TODO
	virtual Zahl Mean(Zahl para, Zahl para1) {return para+para1*M_SQRT2/M_SQRTPI;}
	virtual Zahl Variance(Zahl para, Zahl para1) {return para1*para1*(1-2/M_PI);}
	virtual Zahl Skewness(Zahl para, Zahl para1) {return M_SQRT2*(4-M_PI) / sqrt(pow(M_PI-2, 2));}// TODO: return precalculated value
	virtual Zahl Kurtosis(Zahl para, Zahl para1) {return 8*(M_PI-3)/ ( (M_PI-2)*(M_PI-2) );} // TODO: return precalculated value
	virtual Zahl LogLikelihood(Zahl y, Zahl para, Zahl para1) {return 0;} //TODO
	virtual Zahl CharacteristicFunction(Zahl y, Zahl para, Zahl para1) {return 0;} //TODO
	virtual Zahl Median(Zahl para, Zahl para1) {return para+0.6745*para1;}
	virtual Zahl Q1(Zahl para, Zahl para1) {return para+0.3286*para1;}
	virtual Zahl Q3(Zahl para, Zahl para1) {return para+1.150*para1;}
	virtual Zahl Mode(Zahl para, Zahl para1) {return para;}
	virtual Zahl qMean(Zahl para, Zahl para1) {return 0.5751;}
	virtual Zahl qMode(Zahl para, Zahl para1) {return 0;}
};


/** This class describes the LOG-Normal (Gauss) distribution
	*@author Reinhold Kainhofer
	*/
class DistLogNormal : public DistNormal {
public:
	DistLogNormal(Zahl mm=0, Zahl ss=1):DistNormal(mm,ss) {};
	virtual ~DistLogNormal() {};
	virtual string DistributionName(){return "LogNormal";};
	virtual Zahl pdf(Zahl y, Zahl para, Zahl para1) { return DistNormal::pdf(log(y), para,para1)/y;};
	virtual Zahl Transform(Zahl y, Zahl mm, Zahl ss) {return exp(DistNormal::Transform(y,mm,ss));};
	virtual Zahl cdf(Zahl y, Zahl para, Zahl para1) { return DistNormal::cdf(log(y), para, para1);}

	virtual int Range() {return RANGE_0_INF;}
	virtual Zahl mgf(Zahl y, Zahl para, Zahl para1) {return 0;} //TODO
	virtual Zahl Moment(int k, Zahl para, Zahl para1) {return 0;} //TODO
	virtual Zahl Mean(Zahl para, Zahl para1) {return exp(para+para1*para1/2);}
	virtual Zahl Variance(Zahl para, Zahl para1) {return exp(2*para+para1*para1) * (exp(para1*para1) - 1);}
	virtual Zahl Skewness(Zahl para, Zahl para1) {return (M_E+2) *sqrt(M_E-1);} // TODO: return precalculated value
	virtual Zahl Kurtosis(Zahl para, Zahl para1) {return M_E*M_E*M_E*M_E + 3*M_E*M_E*M_E + 3*M_E*M_E - 6;} //TODO: return precalculated value
	virtual Zahl LogLikelihood(Zahl y, Zahl para, Zahl para1) {return 0;} //TODO
	virtual Zahl CharacteristicFunction(Zahl y, Zahl para, Zahl para1) {return 0;} //TODO
	virtual Zahl Median(Zahl para, Zahl para1) {return exp(para);}
	virtual Zahl Q1(Zahl para, Zahl para1) {return exp(para-0.6745*para1);}
	virtual Zahl Q3(Zahl para, Zahl para1) {return exp(para+0.6745*para1);}
	virtual Zahl Mode(Zahl para, Zahl para1) {return exp(para-para1*para1);}
	virtual Zahl qMean(Zahl para, Zahl para1) {return 0;} //TODO, no simple closed form
	virtual Zahl qMode(Zahl para, Zahl para1) {return 0;} //TODO, no simple closed form
};
#define DIstCobbDouglas DistLogNormal
#define DistAntiLogNormal DistLogNormal


/** This class describes the Pareto distribution
	@author Reinhold Kainhofer
	*/
class DistPareto : public Distribution {
public:
	DistPareto(Zahl aa=1, Zahl bb=1) : Distribution(2, aa, bb) {}
	virtual ~DistPareto() {}
	virtual string DistributionName(){return "Pareto";};
	virtual Zahl pdf(Zahl y, Zahl para, Zahl para1) { return y<para1 ? 0 : (para/para1)*pow(para1/y, para+1); };
	virtual Zahl Transform(Zahl y, Zahl aa, Zahl bb) {return bb*(pow(y, -1/aa));};
	virtual Zahl cdf(Zahl y, Zahl para, Zahl para1) {return 1-pow(para1/y, para);}
	virtual Zahl cdfc(Zahl y, Zahl para, Zahl para1) {return pow(para1/y, para);}

	virtual int Range() {return RANGE_PARA2_INF;} //TODO: is this right???
	virtual int ParameterRange(int paranr) { if (paranr>0 && paranr<=2) return RANGE_0_INF; else return RANGE_NONE;}
	virtual Zahl mgf(Zahl y, Zahl para, Zahl para1) {return 0;} //TODO
	virtual Zahl Moment(int k, Zahl para, Zahl para1) {if (k>=para) return MOMENT_NONE; else {return 0; }} //TODO
	virtual Zahl Mean(Zahl para, Zahl para1) {if (1<para) return (para*para1/(para-1)); else return MOMENT_NONE;}
	virtual Zahl Variance(Zahl para, Zahl para1) {if (2<para) return (para*para1*para1)/( (para-2)*(para-1)*(para-1)); else return MOMENT_NONE;}
	virtual Zahl Skewness(Zahl para, Zahl para1) {if (3<para) return ( 2*(para+1)*sqrt(para-2) ) / ((para-3)*sqrt(para)); else return MOMENT_NONE;}
	virtual Zahl Kurtosis(Zahl para, Zahl para1) {if (4<para) return 6/para*(para*para*para+ para*para-6*para-2) / (para*para-7*para+12); else return MOMENT_NONE;}
	virtual Zahl LogLikelihood(Zahl y, Zahl para, Zahl para1) {return 0;} //TODO
	virtual Zahl CharacteristicFunction(Zahl y, Zahl para, Zahl para1) {return 0;} //TODO
	virtual Zahl Median(Zahl para, Zahl para1) {return para1*pow(2., 1/para);}
	virtual Zahl Q1(Zahl para, Zahl para1) {return para1*pow(4./3., 1/para);}
	virtual Zahl Q3(Zahl para, Zahl para1) {return para1*pow(4., 1/para);}
	virtual Zahl Mode(Zahl para, Zahl para1) {return para1;}
	virtual Zahl qMean(Zahl para, Zahl para1) {return 1-pow((para-1)/para, 1/para);}
	virtual Zahl qMode(Zahl para, Zahl para1) {return 0;}
};


/** This class describes the 0-Point Pareto distribution
	@author Reinhold Kainhofer
	*/
class DistPareto0 : public DistPareto {
public:
	DistPareto0(Zahl aa=1, Zahl bb=1) : DistPareto(aa, bb) {}
	virtual ~DistPareto0() {}
	virtual string DistributionName(){return "Pareto (0-Point)";}
	virtual Zahl pdf(Zahl y, Zahl para, Zahl para1) { return (para/para1)*pow((1+y/para1), -para-1); };
	virtual Zahl Transform(Zahl y, Zahl aa, Zahl bb) {return bb*(pow(1-y, -1/aa)-1);};
//	virtual Zahl TransformCut(Zahl y, Zahl cut, Zahl para); // TODO
	virtual Zahl cdf(Zahl y, Zahl para, Zahl para1) {return 1-pow(1+y/para1, -para);}
	virtual Zahl cdfc(Zahl y, Zahl para, Zahl para1) {return pow(1+y/para1, -para);}

	virtual int Range() {return RANGE_0_INF;}
	virtual Zahl mgf(Zahl y, Zahl para, Zahl para1) {return 0;} //TODO
	virtual Zahl Moment(int k, Zahl para, Zahl para1) {return 0;} //TODO
	virtual Zahl Mean(Zahl para, Zahl para1) {return Moment(1,para,para1);};
	virtual Zahl Variance(Zahl para, Zahl para1) {return Moment(2,para,para1);};
	virtual Zahl Skewness(Zahl para, Zahl para1) {return Moment(3,para,para1);};
	virtual Zahl Kurtosis(Zahl para, Zahl para1) {return Moment(4,para,para1);};
	virtual Zahl LogLikelihood(Zahl y, Zahl para, Zahl para1) {return 0;} //TODO
	virtual Zahl CharacteristicFunction(Zahl y, Zahl para, Zahl para1) {return 0;} //TODO
	virtual Zahl Median(Zahl para, Zahl para1) {return 0;} //TODO
	virtual Zahl Q1(Zahl para, Zahl para1) {return 0;} //TODO
	virtual Zahl Q3(Zahl para, Zahl para1) {return 0;} //TODO
	virtual Zahl Mode(Zahl para, Zahl para1) {return 0;} //TODO
	virtual Zahl qMean(Zahl para, Zahl para1) {return 0;} //TODO
	virtual Zahl qMode(Zahl para, Zahl para1) {return 0;} //TODO
};


/** This class describes the uniform distribution
	*@author Reinhold Kainhofer
	*/
class DistUniform : public Distribution {
public:
	DistUniform(Zahl a=0, Zahl b=1):Distribution(2,a,b) {}
	virtual ~DistUniform(){}
	virtual string DistributionName(){return "Uniform";};
	virtual Zahl pdf(Zahl y, Zahl para, Zahl para1) {return 1./(para1-para);}
	virtual Zahl Transform(Zahl y, Zahl para, Zahl para1) { return para+y*(para1-para);}
	virtual Zahl cdf(Zahl y, Zahl para, Zahl para1) {return max(0, min(1, (y-para)/(para1-para) ));}
	virtual int Range() {return RANGE_PARA1_PARA2;}
	virtual int ParameterRange(int paranr)	{if (0<paranr && paranr<=2 ) return RANGE_MININF_INF; else return RANGE_NONE;}
	virtual Zahl mgf(Zahl y, Zahl para, Zahl para1) {return 0;} //TODO, nicht angegeben
	virtual Zahl Moment(int k, Zahl para, Zahl para1) { return 0;} //TODO: nicht angegeben
	virtual Zahl Mean(Zahl para, Zahl para1) {return (para+para1)/2;}
	virtual Zahl Variance(Zahl para, Zahl para1) {return (para1-para)*(para1-para)/12;}
	virtual Zahl Skewness(Zahl para, Zahl para1) {return 0;}
	virtual Zahl Kurtosis(Zahl para, Zahl para1) {return -1.2;}
	virtual Zahl LogLikelihood(Zahl y, Zahl para, Zahl para1) {return 0;} //TODO
	virtual Zahl CharacteristicFunction(Zahl y, Zahl para, Zahl para1) {return 0;} //TODO
	virtual Zahl Median(Zahl para, Zahl para1) {return (para+para1)/2;}
	virtual Zahl Q1(Zahl para, Zahl para1) {return (3*para+para1)/4;}
	virtual Zahl Q3(Zahl para, Zahl para1) {return (para+3*para1)/4;}
	virtual Zahl Mode(Zahl para, Zahl para1) {return MOMENT_NONE;}
	virtual Zahl qMean(Zahl para, Zahl para1) {return 0.5;}
	virtual Zahl qMode(Zahl para, Zahl para1) {return MOMENT_NONE;}
};
#define DistRectangular DistUniform

// y~Weibull(C,B,A) => X=((y-A)/B)^C ~Exp(1,0)
/** This class describes the Weibull distribution, which is a generalization of the Exp dist, so maybe I should make DistExp make its base classe??? TODO
		@author Reinhold Kainhofer
	*/
class DistWeibull : public Distribution {
public:
	DistWeibull(Zahl aa=1, Zahl bb=1, Zahl A=0) : Distribution(3, aa, bb) {param[2]=A;}
	virtual ~DistWeibull() {}
	virtual string DistributionName(){return "Weibull";}
	virtual Zahl pdf(Zahl y, Zahl para, Zahl para1) { return (para/para1)*pow((y-param[2])/para1, para-1)*exp(-pow((y-param[2])/para1, para) ); };
	virtual Zahl Transform(Zahl y, Zahl mm, Zahl ss) {return param[2]+pow(-log(1-y), 1/mm)*ss;}; //Calculated by use of the	quantile function
	virtual Zahl cdf(Zahl y, Zahl para, Zahl para1) {return 1-exp(-pow((y-param[2])/para1, para) );}
	virtual Zahl cdfc(Zahl y, Zahl para, Zahl para1) {return exp(-pow((y-param[2])/para1, para) );}

	virtual int Range() {return RANGE_PARA1_INF;}
	virtual int ParameterRange(int paranr) {switch (paranr) {
				case 1: case 2: return RANGE_0_INF;
				case 3: return RANGE_MININF_INF;
				default: return RANGE_NONE;}}
	virtual Zahl mgf(Zahl y, Zahl para, Zahl para1) {return 0;} //TODO
	virtual Zahl Moment(int k, Zahl para, Zahl para1) {return 0;} //TODO
	virtual Zahl Mean(Zahl para, Zahl para1) {return param[2]+ para1*tgamma((para+1)/para); }
	virtual Zahl Variance(Zahl para, Zahl para1) {return para1*para1* ( tgamma((para+2)/para) - pow(tgamma((para+1)/para), 2));}
	virtual Zahl Skewness(Zahl para, Zahl para1) {
				Zahl g1=tgamma((para+1)/para), g2=tgamma((para+2)/para), g3=tgamma((para+1)/para);
				return (2*pow(g1, 3) - 3*g1*g2 + g3 ) / pow( g2-g1*g1, 1.5);}
	virtual Zahl Kurtosis(Zahl para, Zahl para1) {
				Zahl g1=tgamma((para+1)/para), g2=tgamma((para+2)/para), g3=tgamma((para+1)/para), g4=tgamma((para+4)/para);
				return (-3*g1*g1*g1*g1 + 6*g1*g1*g2 - 4*g1*g3 + g4 ) / pow(g1*g1-g2, 2) -3;}
	virtual Zahl LogLikelihood(Zahl y, Zahl para, Zahl para1) {return 0;} //TODO
	virtual Zahl CharacteristicFunction(Zahl y, Zahl para, Zahl para1) {return 0;} //TODO
	virtual Zahl Median(Zahl para, Zahl para1) {return param[2] + para1*pow(M_LOG2E, 1/para);}
	virtual Zahl Q1(Zahl para, Zahl para1) {return param[2]+para1*pow(log(4/3.), 1/para);}
	virtual Zahl Q3(Zahl para, Zahl para1) {return param[2]+para1*pow(log(4.), 1/para);}
	virtual Zahl Mode(Zahl para, Zahl para1) {if (para<=1) return param[2]; else return param[2]+para1*pow((para-1)/para, 1/para);}
	virtual Zahl qMean(Zahl para, Zahl para1) {return 1-exp(-pow(tgamma(1+1/para), para));}
	virtual Zahl qMode(Zahl para, Zahl para1) {return 1-exp(1/para-1);}
};
#define DistFrechet DistWeibull
//TODO
// Weibull(1AB) is the exponential distribution
//Weibull (2,1,0) is the standard Rayleigh distribution


/** This initializes a Distribution and returns a pointer to it
	 @param type defines which distribution is created, constants defined in Distributions.h
	 @param para First parameter
	 @param para1 Second parameter (ignored if distribution has just one param)
*/
Distribution*CreateDistribution(long type, Zahl para, Zahl para1);

#endif
