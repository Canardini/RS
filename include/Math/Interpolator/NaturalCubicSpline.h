#ifndef NATURAL_CUBIC_SPLINE_H
#define NATURAL_CUBIC_SPLINE_H


#pragma once 
#include <Math\Interpolator\Interpolator.h>

class NaturalCubicSpline : public Interpolator {

public:
	NaturalCubicSpline(){};
	NaturalCubicSpline(const std::vector<double> & inX, const std::vector<double> &inY);
	virtual double calcY(double inX);
	void setSplineData();
	virtual void updateInterpolator(){ setSplineData(); };
	virtual boost::shared_ptr<Interpolator> clone() const;
private:
	std::vector<double> mCoeff;
};

#endif
