#ifndef ONE_FACTOR_CHEYETTE_H
#define ONE_FACTOR_CHEYETTE_H

#pragma once
#include <Model/IRModel/MarkovianIRModel.h>

class OneFactorCheyette : public MarkovianIRModel
{

public:

	OneFactorCheyette(const DiscountCurve & inDiscountCurve,
					  const DiscountCurve & inForwardCurve,
		              std::vector<double> inMarkovianFactors,
					  double inKappa);


	virtual boost::shared_ptr<MarkovianIRModel> clone() const;

	virtual double getDiscountValue(double inValueDate, double inMaturity) ;
	double getOldDiscountValueByDate(double inValueDate, bgreg::date & inMaturity);
	virtual double getInitSource(double inDate){ return mDiscountCurve->getInstFWrate(inDate); };
	virtual double getFowardValue(double inValueDate, double inT1, double inT2);
	double getOldFowardValueByDate(double inValueDate, bgreg::date & inMaturity);
	virtual double getGtT(double inValueDate, double inMaturity);
	double getGtT(double inValueDate, const bgreg::date & inMaturity);
	virtual double getGtT(const bgreg::date & inValueDate, const bgreg::date & inMaturity);
 
	

	double updateForwardValue(double delta, double FW, double OISspread, std::vector<double> &dx, double G1, double G2);
	virtual double getFWDFdiff(double inOIS, double inDelta, double inG1, double inG2);
	virtual double updateDiscountValue(double DF, std::vector<double> &dx, double G);
	virtual double integralGG(double inTa, double inTb, double inTk, double inTj);
	virtual void updateMeanReversion(double inKappa){ setKappa(inKappa); };
	virtual double getMeanReversion(){ return mKappa; };
	void setKappa(double inKappa){ mKappa = inKappa; };
	double getDiscountValueByDate(double inValueDate, bgreg::date & inMaturity);
	double getFowardValueByDate(double inValueDate, bgreg::date & inT1);
	virtual double getDeterministicDrift(double inValueDate);
	double getY(double inValueDate);
	double integralForY(double v, double w, double inValueDate);
	virtual ~OneFactorCheyette(){};

private:
	double mKappa;


};

#endif
