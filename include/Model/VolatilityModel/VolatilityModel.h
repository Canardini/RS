#ifndef VOLATILITY_MODEL_H
#define VOLATILITY_MODEL_H

#pragma once

#include<vector>
#include <boost/make_shared.hpp>

class VolatilityModel{

public:
	VolatilityModel(double inEntity):mEntity(inEntity){};

	virtual double getVolatility() = 0;
	virtual boost::shared_ptr<VolatilityModel> clone() const = 0;
	void setEntity(double inEntity){ mEntity = inEntity; };
	void setIndex(int inIndex){ mIndex = inIndex; };

protected:
	int mIndex;
	double mEntity;
};



#endif