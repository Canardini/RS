#include<Model/VolatilityModel/ConstantVolatility.h>


boost::shared_ptr<VolatilityModel> ConstantVolatility::clone() const{

	boost::shared_ptr<ConstantVolatility> lPtr = boost::make_shared<ConstantVolatility>(*this);
	return lPtr;
}