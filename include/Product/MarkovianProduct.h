#pragma once


#include<vector>
#include <Model/IRModel/MarkovianIRModelBridge.h>
#include <Calendar\Calendar.h>
class MarkovianProduct 
{

public:
	MarkovianProduct(const MarkovianIRModelBridge & mMarkovianIRModel,
		             int inNumberOfStates,
					 std::vector<double> &inMarkovianStates,
					 double inValueDate);
	int getNumberOfFactors(){ return mNumberOfFactors; }
	virtual double getPrice()=0;
	virtual double updatePrice(std::vector<double> & inFactors) = 0;
	std::vector<double> getMarkovianFactors(){ return mMarkovianFactors; }
	virtual boost::shared_ptr<MarkovianProduct> clone() const = 0;
	void updateMarkovianFactors(std::vector<double> & inFactors);
	double getInitSource(double inDate){ return mMarkovianIRModel.getInitSource(inDate); };
	double getDrifAdjust(double inDate){ return mMarkovianIRModel.getDeterministicDrift(inDate); };
	double getGtT(double inDate, double inMaturity){ return mMarkovianIRModel.getGtT(inDate, inMaturity); };

	bool getIsExercisable(){ return mIsExercisable; };
	void setValueDate(double inValueDate){ mValueDate = inValueDate; };
	virtual double getKeyEntity()=0;
	virtual void changeStructure(const bgreg::date & inDate, double inValueDate) = 0;
	void setVolAndGrid(const std::vector<double> & inVolGrid, const std::vector<double> & inVols){
		mMarkovianIRModel.setVolAndGrid(inVolGrid, inVols);
	}
	virtual ~MarkovianProduct();

protected:
	int mNumberOfFactors;
	double mValueDate;
	std::vector<double> mMarkovianFactors;
	MarkovianIRModelBridge mMarkovianIRModel;
	bool mIsExercisable;

};

