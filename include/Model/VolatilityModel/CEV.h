#ifndef CEV_MODEL_H
#define CEV_MODEL_H


#pragma once

#include<Model/VolatilityModel/VolatilityModel.h>

class CEV : public VolatilityModel{

public:
	CEV(double inAlpha, const std::vector<double> &inLambdas, double inEntity);

	virtual boost::shared_ptr<VolatilityModel> clone() const;
	virtual double getVolatility();
	virtual ~CEV(){};

private:
	std::vector<double> mLambdas;
	double mAlpha;
};

#endif