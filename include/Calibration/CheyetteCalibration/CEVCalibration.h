#ifndef CEV_CALIBRATION_H
#define CEV_CALIBRATION_H

#pragma once
#include<vector>
#include <boost/date_time/gregorian/gregorian_types.hpp>
#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>
#include<Product/InterestRateSwap.h>
#include <Math/PDE/Grid.h>
#include<Math/AnalyticalFormulas/Black.h>

namespace bgreg = boost::gregorian;

class CEVCalibration{

public:
	CEVCalibration(const std::vector<double> &inStripVolatilities,
				   const std::vector<bgreg::date> & inExercisesDates,
				   const InterestRateSwap &inIRS,
				   const bgreg::date & inToday,
				   double inAlpha,
				   TimeGrid & inTimeGrid,
				   int inNumberOfSimulations);

	void runCalibration();
	void getLambdas(std::vector<double> &outLambdas);
	void monteCarloCorrection();

	void getBounds(double &inEX, double &inEY, double &inVarX, double &inVarY);

	std::vector<double> getMarketPrices(){ return mMarketPrices; };
	std::vector<double> getFWPrices(){ return mFWPrices; };
	std::vector<double> getSWR(){ return mSwapRates; };

private:

	std::vector<double> mStripVolatilities;
	std::vector<bgreg::date>  mExercisesDates;
	InterestRateSwap mIRS;
	int mNumberOfVols;
	bgreg::date mToday;
	std::vector<double> mLambdas;
	std::vector<double> mSwapRates;
	std::vector<double> mMarketPrices;
	double mAlpha;
	TimeGrid mTimeGrid;
	int mNumberOfSimulations;
	std::vector<double> mFWPrices;
};

#endif