#ifndef MARKOVIAN_IR_MODEL_H
#define MARKOVIAN_IR_MODEL_H

#pragma once
#include <boost/make_shared.hpp>
#include <Model/VolatilityModel/VolatilityModelBridge.h>
#include<YieldCurve\DiscountCurve.h>

class MarkovianIRModel {

public:
	typedef boost::shared_ptr<DiscountCurve> tDiscountCurvePtr;
	typedef boost::shared_ptr<DiscountCurve> tForwardCurvePtr;

	MarkovianIRModel(const DiscountCurve & inDiscountCurve,
					 const DiscountCurve & inForwardCurve,
					 std::vector<double> inMarkovianFactors);

	double getMarkovianFactors(int index){ return mMarkovianFactors[index]; };
	int getnumberOfFactors(){ return mNumberOfFactors; };
	double  getInitDiscountValue(double inMaturity);
	double getInitDiscountValueByDate(bgreg::date & inMaturity);
	double  getInitFWValue(double inT1, double inT2);
	double getInitFWValueByDate(bgreg::date & inMaturity);
	int locate(const std::vector<double> & xx, double inX0);
	virtual boost::shared_ptr<MarkovianIRModel> clone() const = 0;

	virtual double getDiscountValue(double inValueDate, double inMaturity) = 0;
	virtual double getDiscountValueByDate(double inValueDate, bgreg::date & inMaturity) = 0;

	virtual double getFowardValue(double inValueDate, double inT1, double inT2) = 0;
	virtual double getFowardValueByDate(double inValueDate, bgreg::date & inMaturity) = 0;

	virtual double getGtT(const bgreg::date & inValueDate, const bgreg::date & inMaturity)=0;
	virtual double getGtT(double inValueDate, double inMaturity) = 0;

	
	void updateMarkovianFactors(std::vector<double> & inFactors);

	double getGDiscount(){ return mGDiscount; };
	double getGFWStart(){ return mGFWStart; };
	double getGFWEnd(){ return mGFWEnd; };
	double getFWDFSpread(){ return mFWDFSpread; };
	double getIndexAcc(){ return mIndexAcc; };
	virtual double updateForwardValue(double delta, double FW, double OISspread, std::vector<double> &dx, double G1, double G2) = 0;
	virtual double updateDiscountValue(double DF, std::vector<double> &dx, double G) = 0;
	virtual double integralGG(double inTa, double inTb, double inTk, double inTj) = 0;
	virtual double getDeterministicDrift(double inValueDate) = 0;
	virtual double getInitSource(double inDate)=0;
	bgreg::date getIndexEnd(){ return  mIndexEndDate; };
	virtual double getFWDFdiff(double OIS, double delta, double G1, double G2) = 0;
	virtual void updateMeanReversion(double inKappa) = 0;
	virtual double getMeanReversion()=0;
	void setIndexEnd(bgreg::date & inDate){ mIndexEndDate= mForwardCurve->getEndDate(inDate); };
	virtual ~MarkovianIRModel(){};
    void setVolAndGrid(const std::vector<double> & inVolGrid, const std::vector<double> & inVols){
		mOptionalVols = inVols;
		mOptionalVolTime = inVolGrid;
	};
protected:
	boost::shared_ptr<DiscountCurve> mDiscountCurve;
	boost::shared_ptr<DiscountCurve> mForwardCurve;

	std::vector<double> mMarkovianFactors;
	int mNumberOfFactors;
	double mGDiscount;
	double mGFWStart;
	double mGFWEnd;
	double mFWDFSpread;
	double mIndexAcc;
	bgreg::date mIndexEndDate;

	std::vector<double> mOptionalVolTime;
	std::vector<double> mOptionalVols;
};

#endif