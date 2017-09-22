#ifndef THREE_D_LINEAR_INTERPOLATOR_H
#define THREE_D_LINEAR_INTERPOLATOR_H


#pragma once 
#include<Math\Interpolator\2DLinearInterpolator.h>
#include <Math\Interpolator\3DInterpolator.h>

class ThreeDLinearInterpolator : public ThreeDInterpolator  {

public:
	ThreeDLinearInterpolator(){};
	ThreeDLinearInterpolator(const std::vector<double> & inX, const std::vector<double> &inY,
		const std::vector<double> &inZ, const std::vector<std::vector<std::vector<double>>> &inXYZ);
	double calcXYZ(double inX, double inY, double inZ);
	virtual boost::shared_ptr<ThreeDInterpolator> clone() const ;

private:
};

#endif