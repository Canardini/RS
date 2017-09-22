#ifndef TWO_D_NATURAL_CUBIC_SPLINE_H
#define TWO_D_NATURAL_CUBIC_SPLINE_H


#pragma once 
#include <Math\Interpolator\2DInterpolator.h>
#include<Math\Interpolator\NaturalCubicSpline.h>

class TwoDNaturalCubicSpline : public TwoDInterpolator
{
public:
	TwoDNaturalCubicSpline(){};
	TwoDNaturalCubicSpline(const std::vector<double> & inX, const std::vector<double> &inY, const std::vector<std::vector<double>> &inZ);
	virtual double calcZ(double inX,double inY);
	virtual boost::shared_ptr<TwoDInterpolator> clone() const;
	void virtual setSplineData();

private:
};

#endif