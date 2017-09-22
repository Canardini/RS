#ifndef INTERPOLATOR_BRIDGE_H
#define INTERPOLATOR_BRIDGE_H


#pragma once

#include <Math\Interpolator\Interpolator.h>
#include <boost/make_shared.hpp>

class InterpolatorBridge{

public :
	
	InterpolatorBridge(const InterpolatorBridge & original);
	InterpolatorBridge(const Interpolator & innerInterpolator);
	InterpolatorBridge(){};
	InterpolatorBridge& operator=(const InterpolatorBridge& original);
	double calcY(double inX);
	virtual void updateInterpolator(){ mInterpolatorPtr->updateInterpolator(); };
	void setXY(const std::vector<double> & inX, const std::vector <double> & inY){ mInterpolatorPtr->setXY(inX, inY);};
private :
	boost::shared_ptr<Interpolator> mInterpolatorPtr;
};

#endif