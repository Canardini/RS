#ifndef BLACK_H
#define BLACK_H


namespace Black{
	const double PI = 3.14159265358979;

double black(int sign, double E, double strike, double df,
		double time, double vol);

double black_delta(int sign, double E, double strike, double df,
	double time, double vol);

double black_gamma(int sign, double E, double strike, double df,
	double time, double vol);

double black_vega(double E, double strike, double df,
	double time, double vol);

double black_normal(int sign, double E, double strike, double df,
	double time, double vol);

double gauss_prob(double x);
}










#endif // !BLACK_H
