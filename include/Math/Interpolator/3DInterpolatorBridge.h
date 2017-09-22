#ifndef THREE_D_INTERPOLATOR_BRIDGE_H
#define THREE_D_INTERPOLATOR_BRIDGE_H


#pragma once
#include <Math\Interpolator\3DInterpolator.h>
#include <boost/make_shared.hpp>

class ThreeDInterpolatorBridge{

public:

	ThreeDInterpolatorBridge(const ThreeDInterpolatorBridge & original);
	ThreeDInterpolatorBridge(const ThreeDInterpolator & innerInterpolator);
	ThreeDInterpolatorBridge& operator=(const ThreeDInterpolatorBridge& original);
	double calcXYZ(double inX, double inY, double inZ);
	ThreeDInterpolatorBridge(){};

private:
	boost::shared_ptr<ThreeDInterpolator> mInterpolatorPtr;
};

#endif