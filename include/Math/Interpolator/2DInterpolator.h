#ifndef TWO_D_INTERPOLATOR_H
#define TWO_D_INTERPOLATOR_H

#pragma once

#include<vector>
#include <boost/make_shared.hpp>
#include <Math\Interpolator\InterpolatorBridge.h>

class TwoDInterpolator{

public:
	TwoDInterpolator(const std::vector<double> & inX, const std::vector<double> &inY,
					 const std::vector<std::vector<double>> &inZ)
					 :mX(inX), mY(inY), mZ(inZ)
	{
	};
	TwoDInterpolator(){};
	void virtual setSplineData()=0;
	virtual double calcZ(double inX, double inY) = 0;
	virtual boost::shared_ptr<TwoDInterpolator> clone() const = 0;
	virtual ~TwoDInterpolator(){};

	void findNeighboursPoints(std::vector<int> & outNeighbours, int &Num, double inX, double inY);
	int  locate(const std::vector<double> & xx, double inX0);
	void shiftXYZ(int i, int j, double shift);
	int getsizeX(){ return mX.size(); };
	int getsizeY(){ return mY.size(); };


protected:
	std::vector<double> mX;
	std::vector<double> mY;
	std::vector<std::vector<double>>  mZ;
	std::vector<InterpolatorBridge> mInterpolator;
};

#endif  