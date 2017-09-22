#ifndef INTEREST_RATE_SWAP_H
#define INTEREST_RATE_SWAP_H

#pragma once
#include<Product\MarkovianProduct.h>

class InterestRateSwap :public MarkovianProduct
{
public:

	enum ProductType{Swap,Swaption};

	struct IRSLegInfo{
		bool isFixedOrFloat; 
		bool isPayer;
		std::string currency;
		double NotionalByDef;
		std::string ReferenceRate;
		double Spread;
		double FixedRate;
		std::string Interval;
	    std::string Daycount;
		std::string SettlementAdjust;
		std::string BusinessCalendar;
		bgreg::date SpotDate;
		bgreg::date TerminationDate;
		bgreg::date TodayDate;
		bgreg::date ForwardStartDate;
		std::vector<bgreg::date> ExercisesDates;
		std::string AmmortizationType;
		
	};

	InterestRateSwap(const MarkovianIRModelBridge & mMarkovianIRModel,
					 int inNumberOfStates,
					 std::vector<double> &inMarkovianStates,
		             double inValueDate,
					 const std::vector<IRSLegInfo> & inIRSLegInfo
					 );

	virtual boost::shared_ptr<MarkovianProduct> clone() const;
	double getLegPV(int index);
	void setConstantVariables();
	void  setInitStartDates(int inIndex);
	void  setInitEndDates(int inIndex);
	void  setInitAccPeriods(int inIndex);
	void  setDiscountFactors(int inIndex);
	void  setRates(int inIndex);
	void  setInitNotional(int inIndex);
	void setNextStartDates(int inIndex);
	void setNumberOfFloatingLegs();
	void setAccAdjust(int inIndex);
	double updateIndexPrice(std::vector<double> &inMarkovianStates, int inIndex);
	double updatePrice(std::vector<double> &inMarkovianStates);
	void setAccretingNotional(double inRate=0.0);
	void setForwardDate(int inIndex);
	void setForwardDate(const bgreg::date & inDate);
	void changeStructure(const bgreg::date & inDate, double inValueDate);
	double getSwapRate(){ return mSwapRate; };
	virtual double getKeyEntity(){ return getSwapRate(); };
	double getAnnuity(){ return mAnnuity; };
	bool isPayer(){
		if ((mIRSLegInfo[0].isFixedOrFloat == false) & (mIRSLegInfo[0].isPayer == false)){ return true; }
		else{ return false; }
	}
	void updateMeanReversion(double inKappa){ mMarkovianIRModel.updateMeanReversion(inKappa); };
	double getMeanReversion(){ return mMarkovianIRModel.getMeanReversion(); };
	double ClosedFormIntegral(const bgreg::date &inValueDate, const bgreg::date &inTa, const bgreg::date &inTb);
	void swapRateSensitivity(double inX,double inY,const bgreg::date &inValueDate, bool inBothSensitivities, double & outDx, double & outDxx);
	void swapRateSensitivity(const bgreg::date &inValueDate, bool inBothSensitivities, double & outDx, double & outDxx);

	virtual double getPrice();

	double getFixedRate(){ return mDefaultFixedRate; };

private:

	std::vector<IRSLegInfo> mIRSLegInfo;
	double mSwapRate;
	double mAnnuity;
	double mDefaultFixedRate;
	int mNumberOfLegs;
	int mNumberOfFloatingLegs;
	std::vector<std::vector<double>> mLegAccPeriods;
	std::vector<std::vector<double>> mDiscountFactors;
	std::vector<std::vector<double>> mGDiscount;
	std::vector<std::vector<double>> mGFWStart;
	std::vector<std::vector<double>> mGFWEnd;
	std::vector<std::vector<double>> mRates;
	std::vector<std::vector<double>> mFWDFRates;
	std::vector<std::vector<double>> mIndexAcc;
	std::vector<std::vector<double>> mFWDFspreads;
	std::vector<std::vector<double>> mNotionals;
	std::vector<std::vector<bgreg::date>> mStartDates;
	std::vector<std::vector<bgreg::date>> mIndexStartDates;
	std::vector<std::vector<bgreg::date>> mIndexEndDates;
	std::vector<std::vector<bgreg::date>> mEndDates;
	std::vector<std::vector<std::string>> mEndStrDates;
	std::vector<std::vector<double>> mDoubleStartDates;
	std::vector<std::vector<double>> mDoubleEndDates;
	std::vector<double> mAccAdjust;
	std::vector<int> mNextStartDates;
	std::vector<double> mLegPVs;
};
#endif 