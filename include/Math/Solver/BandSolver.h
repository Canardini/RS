#ifndef BAND_SOLVER_H
#define BAND_SOLVER_H

#pragma once

#include <vector>
#include <Math/Matrix/Matrix.h>
#include <boost/function.hpp>

// Numerical Recipes

#define SWAP(a,b) {dum=(a);(a)=(b);(b)=dum;}

template<class T>
inline T SIGN(const T &a, const T &b)
{
	return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);
}

inline float SIGN(const float &a, const double &b)
{
	return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);
}

inline float SIGN(const double &a, const float &b)
{
	return (float)(b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a));
}

namespace BandSolver{

	void bandec(dMatrix &a, int n, const int m1, const int m2, dMatrix & al, std::vector<int> & indx, double & d);
	void dbanbks(const dMatrix & a, const int n, const int m1, const int m2, const dMatrix & al, std::vector<int> & indx, dVector & b);
	void dtridag(const dVector & a, const dVector & b, const dVector & c, const dVector & res, dVector & u);

}

namespace Optimizer{

	void lnsrch(std::vector<double> &xold, const double fold, std::vector<double> & g, std::vector<double> &p,
		std::vector<double> &x, double &f, const double stpmax, bool &check, boost::function<double(const std::vector<double> &)> inFunc);
	void dfpmin(std::vector<double> &p, const double gtol, int &iter, double &fret, boost::function<double(const std::vector<double> &)> inFunc);
	void gradient(std::vector<double> & inX, double fold, boost::function<double(const std::vector<double> &)> inFunc, std::vector<double> & outGradX);
}

namespace RootFinder{
	double zbrent(boost::function<double(double &)> inFunc, const double x1, const double x2, const double tol, double &fvalue);

}

	
#endif