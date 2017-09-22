#ifndef CONSTANT_VOLATILITY_MODEL_H
#define CONSTANT_VOLATILITY_MODEL_H


#pragma once

#include<Model/VolatilityModel/VolatilityModel.h>

class ConstantVolatility : public VolatilityModel{

public:
	ConstantVolatility(double inVol) :VolatilityModel(0.0),mVolatility(inVol){};

	virtual boost::shared_ptr<VolatilityModel> clone() const;
	virtual double getVolatility(){return mVolatility;};
	
	virtual ~ConstantVolatility(){};
private:
	double mVolatility;
};

#endif