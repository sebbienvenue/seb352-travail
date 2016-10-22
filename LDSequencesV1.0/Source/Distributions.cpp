/***************************************************************************
						Distributions.cpp  -  description
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


#include "Distributions.h"

Distribution*CreateDistribution(long type, Zahl para, Zahl para1) {
	Distribution *dist;
	switch (type) {
		 case USE_EXP_DIST: dist=new DistExp(para, para1); break;
		 case USE_GAMMA_DIST: dist=new DistGamma(para, para1); break;
		 case USE_NORMAL_DIST: dist=new DistNormal(para, para1); break;
		 case USE_LOGNORMAL_DIST: dist=new DistLogNormal(para, para1); break;
		 case USE_PARETO0_DIST: dist=new DistPareto0(para, para1); break;
		 case USE_WEIBULL_DIST: dist=new DistWeibull(para, para1); break;
		 default: cerr<<"Unrecognized distribution "<<type<<" for t, using exponential distribution."<<endl;
				dist=new DistExp(para, para1);
				break;
	}
	return dist;
}

	/** Transform the number y \in [0,1) to an exponentially distributed variable with parameter para*/
Zahl DistGamma::Transform(Zahl y, Zahl para, Zahl para1) {	return 0; }//TODO

	/** Transform the number y \in [0,1) to a random variable which has the (rescaled) density of
	the exponential distribution, but lies only in the interval [0,cut). */
Zahl DistGamma::TransformCut(Zahl y, Zahl cut, Zahl para, Zahl para1){	return 0; }//TODO

	/** coefficients for Hasting's approximation */
#define QUANT_NORM_A0 2.515517
#define QUANT_NORM_A1 0.802853
#define QUANT_NORM_A2 0.010328
#define QUANT_NORM_B1 1.432788
#define QUANT_NORM_B2 0.189269
#define QUANT_NORM_B3 0.001308
	/**	Transform the number y \in [0,1) to a normal distributed variable with N(mm, ss) by first transforming it to N(0,1)
		The quantile function is approximated by a broken rational function. folr 0.5<=y<1 use Hasting's approximation  */
Zahl DistNormal::Transform(Zahl y, Zahl mm, Zahl ss) {
	Zahl t;
	Zahl p=y;
	if (y>=0.5) p=1-y;	// For 0.5<y<1, use 1-y, which is <0.5 so that we can apply the approximation for Q^{-1}, in the end we have to distinguish the two cases again...
	t=sqrt(-2*log(p));
	Zahl enumerator=QUANT_NORM_A0 + t*(QUANT_NORM_A1 + t*QUANT_NORM_A2);
	Zahl denominator=1 + t*(QUANT_NORM_B1+ t*(QUANT_NORM_B2+t*QUANT_NORM_B3));
	Zahl n01variate=t-enumerator/denominator;
	if (y<0.5)	n01variate=-n01variate;
	// now that we have a N(0,1) variate, just stretch and translate it to get a N(mu, sigma) variate
	return n01variate*ss+mm;
}

	/** coefficients for the polynomial approximation of the cdf */
#define CDF_NORM_P 0.2316419
#define CDF_NORM_B1 0.319381530
#define CDF_NORM_B2 -0.356563782
#define CDF_NORM_B3 1.781477937
#define CDF_NORM_B4 -1.821255978
#define CDF_NORM_B5 1.330274429
	/** Approximates the cdf of the normal distribution by a 5th order polynomial, epsilon<7.5x10^{-8} */
Zahl DistNormal::cdf(Zahl y, Zahl para, Zahl para1) {
	Zahl yy=(y-para)/para1; // rescale it to an N(0,1) distribution
	Zahl y1=yy;
	if (y1<0) y1=-yy;
	Zahl t=1/(1+CDF_NORM_P*y1);
	Zahl result=DistNormal::pdf(y1, 0, 1)*t*(CDF_NORM_B1 + t*(CDF_NORM_B2 + t*(CDF_NORM_B3 + t*(CDF_NORM_B4 + t*CDF_NORM_B5))));
	if (yy>0) result=1-result;
	return result;
}
