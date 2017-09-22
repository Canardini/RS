#ifndef THREE_D_INTERPOLATOR_H
#define THREE_D_INTERPOLATOR_H

#pragma once

#include<vector>
#include <boost/make_shared.hpp>
#include <Math\Interpolator\2DInterpolatorBridge.h>

class ThreeDInterpolator{

public:
	ThreeDInterpolator(const std::vector<double> & inX, const std::vector<double> &inY,
		const std::vector<double> &inZ, const std::vector<std::vector<std::vector<double>>> &inXYZ)
		:mX(inX), mY(inY), mZ(inZ), mXYZ(inXYZ)
	{
	};
	ThreeDInterpolator(){};
	virtual double calcXYZ(double inX, double inY, double inZ) = 0;
	virtual boost::shared_ptr<ThreeDInterpolator> clone() const = 0;
	virtual ~ThreeDInterpolator(){};

protected:
	std::vector<double> mX;
	std::vector<double> mY;
	std::vector<double> mZ;
	std::vector<std::vector<std::vector<double>>>  mXYZ;
	std::vector<std::vector<InterpolatorBridge>> mInterpolator;
 
};

#endif  