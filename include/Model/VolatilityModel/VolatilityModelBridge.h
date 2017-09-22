#ifndef VOLATILITY_BRIDGE_H
#define VOLATILITY_BRIDGE_H


#pragma once

#include <Model\VolatilityModel\VolatilityModel.h>
#include <boost/make_shared.hpp>

class VolatilityModelBridge{

public:

	VolatilityModelBridge(const VolatilityModelBridge & original);
	VolatilityModelBridge(const VolatilityModel & innerIRModel);
	VolatilityModelBridge& operator=(const VolatilityModelBridge& original);

	double getVolatility(){ return mVolatilityModelPtr->getVolatility(); };
	void setEntity(double inEntity){ mVolatilityModelPtr->setEntity(inEntity); };
	void setIndex(int inIndex){ mVolatilityModelPtr->setIndex(inIndex); };
	VolatilityModelBridge(){};

private:
	boost::shared_ptr<VolatilityModel> mVolatilityModelPtr;
};

#endif