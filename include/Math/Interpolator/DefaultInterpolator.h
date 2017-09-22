#ifndef DEFAULT_INTERP_H
#define DEFAULT_INTERP_H


#pragma once 
#include <Math\Interpolator\Interpolator.h>

class DefaultInterpolator : public Interpolator {

public:

	DefaultInterpolator() {};
	
	virtual double calcY(double inX){ return 0; };
	
	virtual void updateInterpolator(){};
	virtual boost::shared_ptr<Interpolator> clone() const 	{
		boost::shared_ptr<DefaultInterpolator> lPtr = boost::make_shared<DefaultInterpolator>(*this);
		return lPtr;
	};
	virtual ~DefaultInterpolator(){};

};

#endif
