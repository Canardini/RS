#include<Model/VolatilityModel/PiecewiseConstant.h>


PiecewiseConstant::PiecewiseConstant(const std::vector<double> &inLambdas)
:VolatilityModel(0.0), mLambdas(inLambdas) 
{
}


boost::shared_ptr<VolatilityModel> PiecewiseConstant::clone() const{

	boost::shared_ptr<PiecewiseConstant> lPtr = boost::make_shared<PiecewiseConstant>(*this);
	return lPtr;
}


double PiecewiseConstant::getVolatility(){
	return mLambdas[mIndex] ;
}