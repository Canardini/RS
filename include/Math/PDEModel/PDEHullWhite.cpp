#include <Math/PDEModel/PDEHullWhite.h>

PDEHullWhite::PDEHullWhite(double inKappa, double inf0t, const VolatilityModelBridge & inVolatilityModel)
:PDEmodel(inVolatilityModel), mKappa(inKappa)
{
	mConvections.resize(1);
	mDiffusions.resize(1);
	mSourceTerms.resize(1);
	mFactors.resize(1);
	mInitSource = inf0t;

}

boost::shared_ptr<PDEmodel> PDEHullWhite::clone() const{

	return boost::make_shared<PDEHullWhite>(*this);

}


void PDEHullWhite::updateConvections(const std::vector<double> & inFactors)
{
	mConvections[0] = -mKappa*inFactors[0];

}
void PDEHullWhite::updateDiffusions(const std::vector<double> & inFactors)
{
	mDiffusions[0] = mVolatilityModel.getVolatility();
}
void PDEHullWhite::updateMixedDerivatives(const std::vector<double> & inFactors)
{

}
void PDEHullWhite::updateSourceTerms(const std::vector<double> & inFactors)
{
	if (mFactors[0] != inFactors[0]){
		mSourceTerms[0] = (mInitSource + inFactors[0]);
		mFactors[0] = inFactors[0];
	}
}

void PDEHullWhite::updateConvections(const std::vector<double> & inFactors, int index){
	updateConvections(inFactors);
}

void PDEHullWhite::updateDiffusions(const std::vector<double> & inFactors, int index){
	updateDiffusions(inFactors);
}
void PDEHullWhite::updateSourceTerms(const std::vector<double> & inFactors, int index){
}

 
