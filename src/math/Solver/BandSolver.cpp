#include <Math/Solver/BandSolver.h>
#include <cmath>
#include <cassert>
#include <limits>
#include<iostream>
#include <iomanip>
void BandSolver::dbanbks(const dMatrix & a,const int n,const int m1, const int m2, const dMatrix & al, std::vector<int> &indx, dVector & b){
	int  i, k, l;
	int mm;
	double dum;

	mm = m1 + m2 + 1;
	l = m1;
	for (k = 0; k < n; k++) {
		i = indx[k]-1;
		if (i != k) SWAP(b[k], b[i])
		if (l < n) l++;
		for (i = k + 1; i < l; i++) b[i] -= al[k][i-k-1] * b[k];
	}
	l = 1;
	for (i = n-1; i >= 0; i--) {
		dum = b[i];
		for (k = 1; k < l; k++) dum -= a[i][k] * b[k + i];
		b[i] = dum / a[i][0];
		if (l < mm) l++;
	}

}

void  BandSolver::dtridag(const dVector & a, const dVector & b, const dVector & c, const dVector & r, dVector & u)
{
	int j;
	int n=a.size();
	double bet;
	dVector gam(n);

	assert(b[0] != 0.0);
	u[0] = r[0] / (bet = b[0]);
	for (j = 1; j < n; j++) {
		gam[j] = c[j - 1] / bet;
		bet = b[j] - a[j] * gam[j];
		assert(bet != 0.0);
		u[j] = (r[j] - a[j] * u[j - 1]) / bet;
	}
	for (j = (n - 2); j >= 0; j--)
		u[j] -= gam[j + 1] * u[j + 1];
}

void BandSolver::bandec(dMatrix &a, int n, const int m1, const int m2, dMatrix & al, std::vector<int> & indx, double & d){

	const double TINY = 1.0e-40;
	int i, j, k, l, mm;
	double dum;

	mm = m1 + m2 + 1;
	l = m1;
	for (i = 0; i<m1; i++) {
		for (j = m1 - i; j<mm; j++) a[i][j - l] = a[i][j];
		l--;
		for (j = mm - l - 1; j<mm; j++) a[i][j] = 0.0;
	}
	d = 1.0;
	l = m1;
	for (k = 0; k<n; k++) {
			dum = a[k][0];
		i = k;
		if (l<n) l++;
		for (j = k + 1; j<l; j++) {
			if (abs(a[j][0]) > abs(dum)) {
				dum = a[j][0];
				i = j;
			}
		}
		indx[k] = i + 1;
		if (dum == 0.0) a[k][0] = TINY;

		if (i != k) {
			d = -d;
			for (j = 0; j<mm; j++) SWAP(a[k][j], a[i][j]);
		}
		for (i = k + 1; i<l; i++) {
			
			dum = a[i][0] / a[k][0];
			al[k][i - k - 1] = dum;
			for (j = 1; j<mm; j++) a[i][j - 1] = a[i][j] - dum*a[k][j];
			a[i][mm - 1] = 0.0;
		}
	}
}



void Optimizer::lnsrch(std::vector<double> &xold, const double fold, std::vector<double> & g, std::vector<double> &p,
	std::vector<double> &x, double &f, const double stpmax, bool &check, boost::function<double(const std::vector<double> &)> inFunc){
	const double ALF = 1.0e-4, TOLX = 1.0e-7;
	double a, alam, alam2 = 0.0, alamin, b, disc, f2 = 0.0;
	double rhs1, rhs2, slope = 0.0, sum = 0.0, temp, test, tmplam;
	int i, n = xold.size();
	check = false;
	for (i = 0; i<n; i++) sum += p[i] * p[i];
	sum = sqrt(sum);
	if (sum > stpmax)
	for (i = 0; i < n; i++)
		p[i] *= stpmax / sum;
	for (i = 0; i < n; i++)
		slope += g[i] * p[i];
	if (slope >= 0.0){
		check = false;
		return;
	}
	test = 0.0;
	for (i = 0; i<n; i++) {
		temp = abs(p[i]) / std::max(abs(xold[i]), 1.0);
		if (temp > test) test = temp;
	}
	alamin = TOLX / test;
	alam = 1.0;
	for (;;) {
		for (i = 0; i < n; i++) x[i] = xold[i] + alam*p[i];
		f = inFunc(x);
		if (alam < alamin) {
			for (i = 0; i < n; i++) x[i] = xold[i];
			check = true;
			return;
		}
		else if (f <= fold + ALF*alam*slope) return;
		else {
			if (alam == 1.0)
				tmplam = -slope / (2.0*(f - fold - slope));
			else {
				rhs1 = f - fold - alam*slope;
				rhs2 = f2 - fold - alam2*slope;
				a = (rhs1 / (alam*alam) - rhs2 / (alam2*alam2)) / (alam - alam2);
				b = (-alam2*rhs1 / (alam*alam) + alam*rhs2 / (alam2*alam2)) / (alam - alam2);
				if (a == 0.0) tmplam = -slope / (2.0*b);
				else {
					disc = b*b - 3.0*a*slope;
					if (disc < 0.0) tmplam = 0.5*alam;
					else if (b <= 0.0) tmplam = (-b + sqrt(disc)) / (3.0*a);
					else tmplam = -slope / (b + sqrt(disc));
				}
				if (tmplam>0.5*alam)
					tmplam = 0.5*alam;
			}
		}
		alam2 = alam;
		f2 = f;
		alam = std::max(tmplam, 0.1*alam);
	}
}

void Optimizer::dfpmin(std::vector<double> &p, const double gtol, int &iter, double &fret, boost::function<double(const std::vector<double> &)> inFunc){
	const int ITMAX = 200;
	const double EPS = 3.0e-8;
	const double TOLX = 4 * EPS, STPMX = 100.0;
	bool check;
	double den, fac, fad, fae, fp, stpmax, sum = 0.0, sumdg, sumxi, temp, test;
	int n = p.size();
	std::vector<double> dg(n), g(n), hdg(n), pnew(n), xi(n);
	std::vector<std::vector<double>> hessin(n, std::vector<double>(n));
	fp = inFunc(p);
	Optimizer::gradient(p, fp, inFunc, g);
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) hessin[i][j] = 0.0;
		hessin[i][i] = 1.0;
		xi[i] = -g[i];
		sum += p[i] * p[i];
	}
	stpmax = STPMX*std::max(sqrt(sum), double(n));
	for (int its = 0; its < ITMAX; its++) {
		iter = its;
		Optimizer::lnsrch(p, fp, g, xi, pnew, fret, stpmax, check, inFunc);
		if (!check){
			return;
		}
		fp = fret;
		for (int i = 0; i < n; i++) {
			xi[i] = pnew[i] - p[i];
			p[i] = pnew[i];
		}
		test = 0.0;
		for (int i = 0; i<n; i++) {
			temp = abs(xi[i]) / std::max(abs(p[i]), 1.0);
			if (temp > test) test = temp;
		}
		if (test < TOLX)
			return;
		for (int i = 0; i < n; i++) dg[i] = g[i];
		Optimizer::gradient(p, fp, inFunc, g);
		test = 0.0;
		den = std::max(fret, 1.0);
		for (int i = 0; i<n; i++) {
			temp = abs(g[i])*std::max(abs(p[i]), 1.0) / den;
			if (temp > test) test = temp;
		}
		if (test < gtol)
			return;
		for (int i = 0; i < n; i++)
			dg[i] = g[i] - dg[i];
		for (int i = 0; i < n; i++) {
			hdg[i] = 0.0;
			for (int j = 0; j < n; j++) hdg[i] += hessin[i][j] * dg[j];
		}
		fac = fae = sumdg = sumxi = 0.0;
		for(int i = 0; i<n; i++) {
			fac += dg[i] * xi[i];
			fae += dg[i] * hdg[i];
			sumdg += (dg[i])*(dg[i]);
			sumxi += xi[i] * xi[i];
		}
		if (fac > sqrt(EPS*sumdg*sumxi)) {
			fac = 1.0 / fac;
			fad = 1.0 / fae;
			for (int i = 0; i < n; i++) dg[i] = fac*xi[i] - fad*hdg[i];
			for (int i = 0; i < n; i++) {
				for (int j = i; j < n; j++) {
					hessin[i][j] += fac*xi[i] * xi[j]
						- fad*hdg[i] * hdg[j] + fae*dg[i] * dg[j];
					hessin[j][i] = hessin[i][j];
				}
			}
		}
		for (int i = 0; i < n; i++) {
			xi[i] = 0.0;
			for (int j = 0; j < n; j++) xi[i] -= hessin[i][j] * g[j];
		}
	}
	//throw("too many iterations in dfpmin");
}



void Optimizer::gradient(std::vector<double> & inX, double inFold, boost::function<double(const std::vector<double> &)> inFunc, std::vector<double> & outGradX){
	const double EPS = 1e-8;
	int n = inX.size();
	std::vector<double> xh = inX;
	double fold = inFold;
	//double fold = inFunc(inX);
	for (int j = 0; j<n; j++) {
		double temp = inX[j];
		double h = EPS*abs(temp);
		if (h == 0.0) h = EPS;
		xh[j] = temp + h; 
			h = xh[j] - temp;
		double fh = inFunc(xh);
		xh[j] = temp;
		outGradX[j] = (fh - fold) / h;
	}

}

double RootFinder::zbrent(boost::function<double(double &)> inFunc, const double x1, const double x2, const double tol,double &fvalue){
	const int ITMAX = 100;
	const double EPS = 1e-8;
	int iter;
	double a = x1, b = x2, c = x2, d, e, min1, min2;
	double fa = inFunc(a), fb = inFunc(b), fc, p, q, r, s, tol1, xm;

	if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0)){
		//error
	}
	fc = fb;
	for (iter = 0; iter<ITMAX; iter++) {
		if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
			c = a;
			fc = fa;
			e = d = b - a;
		}
		if (fabs(fc) < fabs(fb)) {
			a = b;
			b = c;
			c = a;
			fa = fb;
			fb = fc;
			fc = fa;
		}
		tol1 = 2.0*EPS*fabs(b) + 0.5*tol;
		xm = 0.5*(c - b);
		if (fabs(xm) <= tol1 || fb == 0.0){
			fvalue = fb;
			return b;
		}
		if (fabs(e) >= tol1 && fabs(fa) > fabs(fb)) {
			s = fb / fa;
			if (a == c) {
				p = 2.0*xm*s;
				q = 1.0 - s;
			}
			else {
				q = fa / fc;
				r = fb / fc;
				p = s*(2.0*xm*q*(q - r) - (b - a)*(r - 1.0));
				q = (q - 1.0)*(r - 1.0)*(s - 1.0);
			}
			if (p > 0.0) q = -q;
			p = fabs(p);
			min1 = 3.0*xm*q - fabs(tol1*q);
			min2 = fabs(e*q);
			if (2.0*p < (min1 < min2 ? min1 : min2)) {
				e = d;
				d = p / q;
			}
			else {
				d = xm;
				e = d;
			}
		}
		else {
			d = xm;
			e = d;
		}
		a = b;
		fa = fb;
		if (fabs(d) > tol1)
			b += d;
		else
			b += SIGN(tol1, xm);
		fb = inFunc(b);
	}
	//"Maximum number of iterations exceeded in zbrent");
	return 0.0;
}




