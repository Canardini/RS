#ifndef DD_MODEL_H
#define DD_MODEL_H


#pragma once
#include<Model/VolatilityModel/VolatilityModel.h>

class DisplacedDiffusion : public VolatilityModel{

public:
	DisplacedDiffusion(const std::vector<double> &inLambdas, 
		const std::vector<double> &inBetas, double inEntity, const std::vector<double>& inShift);

	
	virtual boost::shared_ptr<VolatilityModel> clone() const;
	virtual double getVolatility();
	virtual ~DisplacedDiffusion(){};

private:
	std::vector<double> mLambdas;
	std::vector<double> mBetas;
	std::vector<double> mShift;
};





#endif