#include <Math\PDEModel\PDEModelBridge.h>




PDEModelBridge::PDEModelBridge(const PDEModelBridge& inOriginal){
	mPDEModelPtr = inOriginal.mPDEModelPtr->clone();
}

PDEModelBridge::PDEModelBridge(const PDEmodel& inInnerPDEModel){
	mPDEModelPtr = inInnerPDEModel.clone();
}

PDEModelBridge& PDEModelBridge::operator=(const PDEModelBridge& inOriginal){

	if (this != &inOriginal){
		mPDEModelPtr = inOriginal.mPDEModelPtr->clone();
	}
	return *this;
}

