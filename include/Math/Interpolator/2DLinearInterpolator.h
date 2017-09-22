#ifndef TWO_D_LINEAR_INTERPOLATOR_H
#define TWO_D_LINEAR_INTERPOLATOR_H


#pragma once 
#include <Math\Interpolator\2DInterpolator.h>
#include<Math\Interpolator\LinearInterpolation.h>

class TwoDLinearInterpolator : public TwoDInterpolator
{

public:
	TwoDLinearInterpolator(){};
	TwoDLinearInterpolator(const std::vector<double> & inX, const std::vector<double> &inY, const std::vector<std::vector<double>> &inZ);
	virtual double calcZ(double inX, double inY);
	void virtual setSplineData();
	virtual boost::shared_ptr<TwoDInterpolator> clone() const;

private:
};

#endif