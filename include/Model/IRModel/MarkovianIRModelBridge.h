#ifndef IRMODEL_BRIDGE_H
#define IRMODEL_BRIDGE_H


#pragma once

#include <Model\IRModel\MarkovianIRModel.h>
#include <boost/make_shared.hpp>

class MarkovianIRModelBridge{

public:

	MarkovianIRModelBridge(const MarkovianIRModelBridge & original);
	MarkovianIRModelBridge(const MarkovianIRModel & innerIRModel);
	MarkovianIRModelBridge& operator=(const MarkovianIRModelBridge& original);

	double getDiscountValue(double inValueDate, double inMaturity);
	double getDiscountValueByDate(double inValueDate, bgreg::date & inMaturity);
	double getFowardValue(double inValueDate, double inT1, double inT2);
	double getFowardValueByDate(double inValueDate, bgreg::date & inMaturity) ;
	double getGtT(double inValueDate, double inMaturity){
		return mMarkovianIRModelPtr->getGtT(inValueDate, inMaturity);
	};

	double getGtT(const bgreg::date & inValueDate, const bgreg::date & inMaturity);
	double getInitDiscountValue(double inMaturity);
	double getInitDiscountValueByDate(bgreg::date & inMaturity);
	double getInitFWValueByDate(bgreg::date & inStartDate);
	void   updateMarkovianFactors(std::vector<double> & inFactors);
	double getGDiscount(){ return mMarkovianIRModelPtr->getGDiscount(); };
	double getGFWStart(){ return mMarkovianIRModelPtr->getGFWStart(); };
	double getGFWEnd(){ return mMarkovianIRModelPtr->getGFWEnd(); };
	double getFWDFSpread(){ return mMarkovianIRModelPtr->getFWDFSpread(); };
	double getIndexAcc(){ return mMarkovianIRModelPtr->getIndexAcc(); };
	bgreg::date getIndexEnd(){ return  mMarkovianIRModelPtr->getIndexEnd(); };
	double updateForwardValue(double delta, double FW, double OISspread, std::vector<double> &dx, double G1, double G2);
	double updateDiscountValue(double DF, std::vector<double> &dx, double G){ return mMarkovianIRModelPtr->updateDiscountValue(DF, dx, G); };
	double getInitSource(double inDate){ return  mMarkovianIRModelPtr->getInitSource(inDate); };
	virtual double getFWDFdiff(double inOIS, double inDelta, double inG1, double inG2){
		return mMarkovianIRModelPtr->getFWDFdiff(inOIS,inDelta, inG1,inG2);
	};
	void updateMeanReversion(double inKappa){ return mMarkovianIRModelPtr->updateMeanReversion(inKappa); };
	double getMeanReversion(){ return mMarkovianIRModelPtr->getMeanReversion(); };

	virtual double integralGG(double inTa, double inTb, double inTk, double inTj){ return mMarkovianIRModelPtr->integralGG(inTa,inTb,inTk,inTj); };
	
	virtual double getDeterministicDrift(double inValueDate){ return mMarkovianIRModelPtr->getDeterministicDrift(inValueDate); };
	void setVolAndGrid(const std::vector<double> & inVolGrid, const std::vector<double> & inVols){
		mMarkovianIRModelPtr->setVolAndGrid(inVolGrid, inVols);
	}
	MarkovianIRModelBridge(){};

private:
	boost::shared_ptr<MarkovianIRModel> mMarkovianIRModelPtr;
};

#endif