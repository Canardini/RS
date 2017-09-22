#ifndef TWO_D_INTERPOLATOR_BRIDGE_H
#define TWO_D_INTERPOLATOR_BRIDGE_H


#pragma once
#include <Math\Interpolator\2DInterpolator.h>
#include <boost/make_shared.hpp>

class TwoDInterpolatorBridge{

public:

	TwoDInterpolatorBridge(const TwoDInterpolatorBridge & original);
	TwoDInterpolatorBridge(const TwoDInterpolator & innerInterpolator);
	TwoDInterpolatorBridge& operator=(const TwoDInterpolatorBridge& original);
	double calcZ(double inX, double inY);
	TwoDInterpolatorBridge(){};

	void findNeighboursPoints(std::vector<int> & outNeighbours, int &Num, double inX, double inY);
	void shiftXYZ(int i, int j,double shift){ mInterpolatorPtr->shiftXYZ(i, j,shift); };
	void setSplineData(){ mInterpolatorPtr->setSplineData(); };
	int getsizeX(){ return mInterpolatorPtr->getsizeX(); };
	int getsizeY(){ return mInterpolatorPtr->getsizeY(); };
private:
	boost::shared_ptr<TwoDInterpolator> mInterpolatorPtr;
};

#endif