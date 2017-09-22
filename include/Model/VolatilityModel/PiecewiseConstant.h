#ifndef PIECEWISE_CONSTANT_H
#define PIECEWISE_CONSTANT_H

#pragma once

#include<Model/VolatilityModel/VolatilityModel.h>

class PiecewiseConstant : public VolatilityModel{

public:
	PiecewiseConstant(const std::vector<double> &inLambdas );

	virtual boost::shared_ptr<VolatilityModel> clone() const;
	virtual double getVolatility();
	virtual ~PiecewiseConstant(){};

private:
	std::vector<double> mLambdas;
 
};

#endif