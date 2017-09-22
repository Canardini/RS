#include<Math/AnalyticalFormulas/Black.h>
#include<cmath>

double Black::black(int sign, double E, double strike, double df,
	double time, double vol)
{
	double d1, d2, pv;
	vol = (vol<1e-6) ? 1e-6 : vol;
	if (time*vol <= 0) {
		pv = (sign*(E - strike) > 0) ? (sign*(E - strike)) : 0;
		pv *= df;
		return pv;
	}

	if (strike < 0) {
		pv = (sign>0 ? df*(E - strike) : 0);
		return pv;
	}

	d1 = (log(E / strike) + (vol*vol / 2)*time) / (vol*sqrt(time));
	d2 = d1 - vol*sqrt(time);
	pv = df*(sign*E*gauss_prob(sign*d1) - sign*strike*gauss_prob(sign*d2));

	return pv;
}

double Black::black_delta(int sign, double E, double strike, double df,
	double time, double vol)
{
	double d1, delta;

	vol = (vol<1e-6) ? 1e-6 : vol;
	d1 = (log(E / strike) + (vol*vol / 2)*time) / (vol*sqrt(time));
	if (sign==1)
		delta = df*gauss_prob(d1);
	else
		delta = df*(gauss_prob(d1) - 1);
	return delta;
}

double Black::black_gamma(int sign, double E, double strike, double df,
	double time, double vol)
{
	double d1, dd1, gamma;

	vol = (vol<1e-6) ? 1e-6 : vol;
	d1 = (log(E / strike) + (vol*vol / 2)*time) / (vol*sqrt(time));
	dd1 = 1 / (E*vol*sqrt(time));
	gamma = df*exp(-d1*d1 / 2)*dd1 / sqrt(2 * PI);
	return gamma;
}

double Black::black_vega(double E, double strike, double df,
	double time, double vol)
{
	double d1, vega;
	vol = (vol<1e-6) ? 1e-6 : vol;
	d1 = (log(E / strike) + (vol*vol / 2)*time) / (vol*sqrt(time));
	vega = df*E*exp(-d1*d1 / 2)*sqrt(time);
	return vega;
}

double Black::black_normal(int sign, double E, double strike, double df,
	double time, double vol)
{
	double d1, pv;

	if (time*vol <= 0) {
		pv = (sign*(E - strike) > 0) ? (sign*(E - strike)) : 0;
		pv *= df;
		return pv;
	}


	d1 = (strike - E) / (vol*sqrt(time));
	if (sign==1)
		pv = df*((E - strike)*(1 - gauss_prob(d1)) + (vol*sqrt(time) / sqrt(2 * PI))*exp(-d1*d1 / 2));
	else
		pv = df*((strike - E)*gauss_prob(d1) + (vol*sqrt(time) / sqrt(2 * PI))*exp(-d1*d1 / 2));
	return pv;
}



 

double Black::gauss_prob(double x)
/* gaussian density function */
{

	double        GtimesP, Prob, t, absx, G;
	double        p = .2316419;
	double        b[] = { .319381530, -.356563782, 1.781477937, -1.821255978, 1.330274429 };

	//Suzuki put in a new erf to calculate normal distribution 11/06/12
	//Old function is about 10% faster but the result is poor, so we suppress (flase)

	if (false){
		absx = fabs(x);

		if (absx >  5.5)
			GtimesP = 0;
		else {
			t = 1 / (1 + p*absx);
			G = .3989422804014 * exp(-absx*absx / 2);
			GtimesP = G*((((b[4] * t + b[3])*t + b[2])*t + b[1])*t + b[0])*t;
		}

		Prob = (x >= 0) ? (1 - GtimesP) : GtimesP;

	}
	else{
		//new function
		//Prob =  0.5 + 0.5*erf(x/sqrt(2.0));
		Prob = 0.5 + 0.5*erf(x*0.7071067811865474617150084668538);
	}
	return Prob;
}





