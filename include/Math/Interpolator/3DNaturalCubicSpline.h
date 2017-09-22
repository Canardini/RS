#ifndef THREE_D_NATURAL_CUBIC_SPLINE_H
#define THREE_D_NATURAL_CUBIC_SPLINE_H


#pragma once 
#include<Math\Interpolator\2DNaturalCubicSpline.h>
#include <Math\Interpolator\3DInterpolator.h>

class ThreeDNaturalCubicSpline : public ThreeDInterpolator  {

public:
	ThreeDNaturalCubicSpline(){};
	ThreeDNaturalCubicSpline(const std::vector<double> & inX, const std::vector<double> &inY,
	const std::vector<double> &inZ, const std::vector<std::vector<std::vector<double>>> &inXYZ);
	double calcXYZ(double inX, double inY, double inZ);
	virtual boost::shared_ptr<ThreeDInterpolator> clone() const;

private:

};

#endif