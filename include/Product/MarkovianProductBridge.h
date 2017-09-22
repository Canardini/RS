#ifndef MARKOVIAN_PRODUCT_BRIDGE_H
#define MARKOVIAN_PRODUCT_BRIDGE_H


#pragma once

#include <Product\MarkovianProduct.h>
#include <boost/make_shared.hpp>

class MarkovianProductBridge{

public:

	MarkovianProductBridge(const MarkovianProductBridge & original);
	MarkovianProductBridge(const MarkovianProduct & innerIRModel);
	MarkovianProductBridge& operator=(const MarkovianProductBridge& original);

	void updateMarkovianFactors(std::vector<double> & inFactors);
	double getPrice();
	double updatePrice(std::vector<double> & inFactors){ return mMarkovianProductPtr->updatePrice(inFactors); };
	bool getIsExercisable(){ return mMarkovianProductPtr->getIsExercisable(); };
	double getInitSource(double inDate){ return mMarkovianProductPtr->getInitSource(inDate); };
	double getDrifAdjust(double inDate){ return mMarkovianProductPtr->getDrifAdjust(inDate); };
	double  getGtT(double inDate, double inMaturity){ return mMarkovianProductPtr->getGtT(inDate,inMaturity); };
	virtual void changeStructure(const bgreg::date & inDate, double inValueDate) { mMarkovianProductPtr->changeStructure(inDate, inValueDate); };
	double getKeyEntity(){ return mMarkovianProductPtr->getKeyEntity(); };
	void setValueDate(double inValueDate){ mMarkovianProductPtr->setValueDate(inValueDate); };
	void setVolAndGrid(const std::vector<double> & inVolGrid, const std::vector<double> & inVols){
		mMarkovianProductPtr->setVolAndGrid(inVolGrid, inVols);
	}
	MarkovianProductBridge(){};

private:
	boost::shared_ptr<MarkovianProduct> mMarkovianProductPtr;
};

#endif