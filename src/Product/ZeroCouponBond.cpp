#include<Product\ZeroCouponBond.h>

ZeroCouponBond::ZeroCouponBond(const MarkovianIRModelBridge & inMarkovianIRModel,
	int inNumberOfStates,
	std::vector<double> &inMarkovianStates,
	double inValueDate,
	double inMaturity)
	: MarkovianProduct(inMarkovianIRModel, inNumberOfStates, inMarkovianStates, inValueDate), mMaturity(inMaturity)
{
}

boost::shared_ptr<MarkovianProduct> ZeroCouponBond::clone() const{

	boost::shared_ptr<ZeroCouponBond> lPtr = boost::make_shared<ZeroCouponBond>(*this);
	return lPtr;
}

double ZeroCouponBond::getPrice(){
	return mMarkovianIRModel.getDiscountValue(mValueDate, mMaturity);
}

double ZeroCouponBond::updatePrice(std::vector<double> & inFactors){
	updateMarkovianFactors(inFactors);
	return getPrice();
}


