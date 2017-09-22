#include <Math/PDEModel/PDE1FCheyette.h>

PDE1FCheyette::PDE1FCheyette(double inKappa, double inf0t, const VolatilityModelBridge & inVolatilityModel)
:PDEmodel(inVolatilityModel), mKappa(inKappa)
{
	mConvections.resize(2);
	mDiffusions.resize(2);
	mSourceTerms.resize(2);
	mFactors.resize(2);
	mInitSource = inf0t;

}

boost::shared_ptr<PDEmodel> PDE1FCheyette::clone() const{

	return boost::make_shared<PDE1FCheyette>(*this);

}


void PDE1FCheyette::updateConvections(const std::vector<double> & inFactors)
{
		double lVol = mVolatilityModel.getVolatility();
		mConvections[0] = inFactors[1] - mKappa*inFactors[0];
		mConvections[1] = lVol*lVol - 2*mKappa*inFactors[1];
}
void PDE1FCheyette::updateDiffusions(const std::vector<double> & inFactors)
{
		mDiffusions[0] = mVolatilityModel.getVolatility();
		mDiffusions[1] = 0.0;

}
void PDE1FCheyette::updateMixedDerivatives(const std::vector<double> & inFactors)
{

}
void PDE1FCheyette::updateSourceTerms(const std::vector<double> & inFactors)
{
	if (mFactors[0] != inFactors[0]){
		mSourceTerms[0] = (mInitSource + inFactors[0]);
		mSourceTerms[1] = mSourceTerms[0];
		mFactors[0] = inFactors[0];
	}
}

void PDE1FCheyette::updateConvections(const std::vector<double> & inFactors, int index){
	if (index == 0){
		mConvections[0] = inFactors[1] - mKappa*inFactors[0];
	}
	else{
		double lVol = mVolatilityModel.getVolatility();
		mConvections[1] = lVol*lVol - 2 * mKappa*inFactors[1];
	}
}

void PDE1FCheyette::updateDiffusions(const std::vector<double> & inFactors, int index){
	if (index == 0){
		mConvections[0] = inFactors[1] - mKappa*inFactors[0];
	}
	else{
		double lVol = mVolatilityModel.getVolatility();
		mConvections[1] = lVol*lVol - 2 * mKappa*inFactors[1];
	}
}
void PDE1FCheyette::updateSourceTerms(const std::vector<double> & inFactors, int index){
}

PDE1FCheyette::~PDE1FCheyette()
{
}
