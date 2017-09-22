#include <Math/PDEModel/PDEPC1FCheyette.h>

PDEPC1FCheyette::PDEPC1FCheyette(double inKappa, double inf0t, const VolatilityModelBridge & inVolatilityModel)
:PDEmodel(inVolatilityModel), mKappa(inKappa)
{
	mConvections.resize(1);
	mDiffusions.resize(1);
	mSourceTerms.resize(1);
	mFactors.resize(1);
	mInitSource = inf0t;

}

boost::shared_ptr<PDEmodel> PDEPC1FCheyette::clone() const{

	return boost::make_shared<PDEPC1FCheyette>(*this);

}


void PDEPC1FCheyette::updateConvections(const std::vector<double> & inFactors)
{
	mConvections[0] = mDriftAdjust-mKappa*inFactors[0];
	if (mFactors[0] != inFactors[0]){
		mFactors[0] = inFactors[0];
	}

}
void PDEPC1FCheyette::updateDiffusions(const std::vector<double> & inFactors)
{
	mDiffusions[0] = mVolatilityModel.getVolatility();
	if (mFactors[0] != inFactors[0]){
		mFactors[0] = inFactors[0];
	}
}
void PDEPC1FCheyette::updateMixedDerivatives(const std::vector<double> & inFactors)
{

}
void PDEPC1FCheyette::updateSourceTerms(const std::vector<double> & inFactors)
{
		mSourceTerms[0] = (mInitSource + inFactors[0]);
		
		if (mFactors[0] != inFactors[0]){
			mFactors[0] = inFactors[0]; 
		}
}

void PDEPC1FCheyette::updateConvections(const std::vector<double> & inFactors, int index){
	updateConvections(inFactors);
}

void PDEPC1FCheyette::updateDiffusions(const std::vector<double> & inFactors, int index){
	updateDiffusions(inFactors);
}
void PDEPC1FCheyette::updateSourceTerms(const std::vector<double> & inFactors, int index){
}


