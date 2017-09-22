#ifndef SWAPTION_H
#define SWAPTION_H

#pragma once
#include<Product\InterestRateSwap.h>

class Swaption :public InterestRateSwap
{
public:

	Swaption(const MarkovianIRModelBridge & mMarkovianIRModel,
		int inNumberOfStates,
		std::vector<double> &inMarkovianStates,
		double inValueDate,
		const std::vector<IRSLegInfo> & inIRSLegInfo,
		std::vector<bgreg::date> &inExerciseDates
		);

	virtual boost::shared_ptr<MarkovianProduct> clone() const;
	virtual double getPrice();

private:
	std::vector<bgreg::date> mExerciseDates;
	int mNumberOfExerciseDates;
	std::vector<double> mDoubleExerciseDates;
};
#endif 