#include<Model/VolatilityModel/CEV.h>


CEV::CEV(double inAlpha, const std::vector<double> &inLambdas, double inEntity)
:VolatilityModel(inEntity),mLambdas(inLambdas), mAlpha(inAlpha)
{
}


boost::shared_ptr<VolatilityModel> CEV::clone() const{

	boost::shared_ptr<CEV> lPtr = boost::make_shared<CEV>(*this);
	return lPtr;
}


double CEV::getVolatility(){
	return mLambdas[mIndex] * pow(fabs(mEntity), mAlpha);
}