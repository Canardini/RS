#include<Model/VolatilityModel/DisplacedDiffusion.h>

DisplacedDiffusion::DisplacedDiffusion(const std::vector<double> &inLambdas,
	const std::vector<double> &inBetas, double inEntity, const std::vector<double>& inShift) :
	VolatilityModel(inEntity), mLambdas(inLambdas), mBetas(inBetas), mShift(inShift)
{
}

boost::shared_ptr<VolatilityModel> DisplacedDiffusion::clone() const{

	boost::shared_ptr<DisplacedDiffusion> lPtr = boost::make_shared<DisplacedDiffusion>(*this);
	return lPtr;
}

double DisplacedDiffusion::getVolatility(){
	return mLambdas[mIndex] *(mBetas[mIndex]*mEntity+(1-mBetas[mIndex])*mShift[mIndex]);
}