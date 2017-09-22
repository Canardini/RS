#ifndef ZERO_COUPON_BOND_H
#define ZERO_COUPON_BOND_H

#pragma once
#include<Product\MarkovianProduct.h>

class ZeroCouponBond :public MarkovianProduct
{
public :
	ZeroCouponBond(const MarkovianIRModelBridge & mMarkovianIRModel,
				   int inNumberOfStates,
		           std::vector<double> &inMarkovianStates,
		           double inValueDate,
		           double inMaturity);
	virtual boost::shared_ptr<MarkovianProduct> clone() const;
	virtual double getKeyEntity(){ return 0.0;};
	virtual double getPrice();
	virtual double updatePrice(std::vector<double> & inFactors);
	virtual void changeStructure(const bgreg::date & inDate, double inValueDate){};
private :
	double mMaturity;
};
#endif // !ZERO_COUPON_BOND_H
