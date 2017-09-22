#ifndef DD_CALIBRATION_H
#define DD_CALIBRATION_H

#pragma once
#include<vector>
#include<Product/InterestRateSwap.h>
#include <Math/Solver/DifferentialEvolution.h>
#include <boost/date_time/gregorian/gregorian_types.hpp>
#include<Math/AnalyticalFormulas/Black.h>
#include <Math/PDE/Grid.h>
#include <boost/bind.hpp>
#include <limits>
namespace bgreg = boost::gregorian;

class DDCalibration{

public:

	
	struct shiftedLogNormalfunction{
		double mSwapRate;
		double mFxdRate1;
		double mFxdRate2;
		double mFxdRate3;
		double mBSPrice1;
		double mBSPrice2;
		double mBSPrice3;
		double mW1;
		double mW2;
		double mW3;
		double mAnnuity;
		double mMaturity;
		double mVol1;
		double mVol2;
		double mVol3;
		int mCallOrPut;
		double operator()(const std::vector<double> & inX){
			return SLNObjectiveFunction(inX);
		}
		double SLNObjectiveFunction(const std::vector<double> & inX){
			double v1, v2, v3;
			double eta, beta, ratio, vol, lambda;
			double res1,res2,res3;
			lambda = inX[0];
			eta = inX[1];
			beta = (1.0 - eta)*mSwapRate;
			vol = fabs(eta*lambda);
			if (eta == 0.0){
				ratio = 0.0;
			}
			else{
				ratio = beta / eta;
			}
			v1 = Black::black(mCallOrPut, mSwapRate + ratio, mFxdRate1 + ratio, mAnnuity, mMaturity, vol);
			v2 = Black::black(mCallOrPut, mSwapRate + ratio, mFxdRate2 + ratio, mAnnuity, mMaturity, vol);
			v3 = Black::black(mCallOrPut, mSwapRate + ratio, mFxdRate3 + ratio, mAnnuity, mMaturity, vol);

			if ((mBSPrice1 == 0) | (mW1 == 0)){
				res1 = 0.0;
			}
			else{
				res1 = mW1*pow((v1 - mBSPrice1) / mBSPrice1, 2);
			}
			if ((mBSPrice2 == 0) | (mW2 == 0)){
				res2 = 0.0;
			}
			else{
				res2 = mW2*pow((v2 - mBSPrice2) / mBSPrice2, 2);
			}
			if ((mBSPrice3 == 0) | (mW3 == 0)){
				res3 = 0.0;
			}
			else{
				res3 = mW3*pow((v3 - mBSPrice3) / mBSPrice3, 2);
			}
			return res1+res2+res3;
		}
	};

	void fillSLNStructure(shiftedLogNormalfunction & inObj, double inSwapRate,
		double inFxdRate1,
		double inFxdRate2,
		double inFxdRate3,
		double inBSPrice1,
		double inBSPrice2,
		double inBSPrice3,
		double inW1,
		double inW2,
		double inW3,
		double inAnnuity,
		double inMaturity,
		double inVol1,
		double inVol2,
		double inVol3,
		int inCallOrPut){
		inObj.mSwapRate = inSwapRate;
		inObj.mFxdRate1 = inFxdRate1;
		inObj.mFxdRate2 = inFxdRate2;
		inObj.mFxdRate3 = inFxdRate3;
		inObj.mBSPrice1 = inBSPrice1;
		inObj.mBSPrice2 = inBSPrice2;
		inObj.mBSPrice3 = inBSPrice3;
		inObj.mW1 = inW1;
		inObj.mW2 = inW2;
		inObj.mW3 = inW3;
		inObj.mAnnuity = inAnnuity;
		inObj.mMaturity = inMaturity;
		inObj.mVol1 = inVol1;
		inObj.mVol2 = inVol2;
		inObj.mVol3 = inVol3;
		inObj.mCallOrPut = inCallOrPut;
	};

	DDCalibration(const std::vector<double> &inMinVols, const std::vector<double> &inCenteredVols,
		const std::vector<double> &inPlusVols,
		const std::vector<double> &inStrikes,
		const std::vector<bgreg::date> & inExercisesDates,
		const InterestRateSwap &inIRS,
		const TimeGrid &inTimeGrid,
		const bgreg::date & inToday);

	void runCalibration();
	void runLambdaCalibration();
	void runBetaCalibration();
	void runMarketLambdaBetaCalibration();
	double RectangularGammaIntegral(bool inUsePivot, double SWRate, double lambda, double gamma, int flag, const bgreg::date &lT1, const bgreg::date &lT2, int n);
	void getLambdas(std::vector<double> &outLambdas);
	void getBetas(std::vector<double> &outBetas);
	void getSWRates(std::vector<double> &outSWRates);
	double getAverageY(double inT);
	void getAverageX0(double inT);
	void getAverageX(double inT);
	double getApproximatedVariance(double inT);
	double swapRateFunction(double inX);
	double integrator(const std::vector<double> & inX);
	void advancedCalibration(int inIndex);

	std::vector<double> getCenteredMarketPrice(){ return mCenteredMarketPrice; };
	std::vector<double> getFWPrices(){ return mFWPrices; };
	std::vector<double> getMarketErrors(){ return mMarketErrors; };



private:

	std::vector<double> mCenteredVols;
	std::vector<double> mMinVols;
	std::vector<double> mPlusVols;
	std::vector<bgreg::date>  mExercisesDates;
	std::vector<double>  mDoubleExercisesDates;
	InterestRateSwap mIRSpivot;
	InterestRateSwap mIRS;
	int mNumberOfVols;
	bgreg::date mToday;
	std::vector<double> mLambdas;
	std::vector<double> mBetas;
	std::vector<double> mMarketLambdas;
	std::vector<double> mMarketErrors;
	std::vector<double> mMarketBetas;
	std::vector<double> mSwapRate;
	std::vector<double> mCenteredMarketPrice;
	std::vector<double> mMinusMarketPrice;
	std::vector<double> mPlusMarketPrice;
	std::vector<double> mCenteredCalibratedPrice;
	std::vector<double> mMinusCalibratedPrice;
	std::vector<double> mPlusCalibratedPrice;
	std::vector<double> mStrikes;
	int mNumSteps;
	double mKappa;
	std::vector<double> mTempMarkov;
	double mtempSwapRate;
	TimeGrid mTimeGrid;
	double mAdvancedLambda;
	double mAdvancedBeta;
	std::vector<double> mAdvancedMarketLambdas;
	std::vector<double> mAdvancedMarketBetas;
	std::vector<double> mFWPrices;

	int mCompteur;

};

#endif