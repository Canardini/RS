#ifndef LINEAR_INTERPOLATION_H
#define LINEAR_INTERPOLATION_H


#pragma once 
#include <Math\Interpolator\Interpolator.h>

class LinearInterpolator : public Interpolator {

public:
	LinearInterpolator(){};
	LinearInterpolator(const std::vector<double> & inX, const std::vector<double> &inY);
	virtual double calcY(double inX);
	virtual void updateInterpolator(){};
	virtual boost::shared_ptr<Interpolator> clone() const;
private:
	
};

#endif
