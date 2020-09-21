#include<Model\VolatilityModel\VolatilityModelBridge.h>




VolatilityModelBridge::VolatilityModelBridge(const VolatilityModelBridge& inOriginal){
	mVolatilityModelPtr = inOriginal.mVolatilityModelPtr->clone();
}

VolatilityModelBridge::VolatilityModelBridge(const VolatilityModel& inInnerIRModel){
	mVolatilityModelPtr = inInnerIRModel.clone();
}

VolatilityModelBridge& VolatilityModelBridge::operator=(const VolatilityModelBridge& inOriginal){

	if (this != &inOriginal){
		mVolatilityModelPtr = inOriginal.mVolatilityModelPtr->clone();
	}
	return *this;
}

