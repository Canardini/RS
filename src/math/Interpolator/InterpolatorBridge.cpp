#include<Math\Interpolator\InterpolatorBridge.h>



InterpolatorBridge::InterpolatorBridge(const InterpolatorBridge& inOriginal){
	mInterpolatorPtr = inOriginal.mInterpolatorPtr->clone();
}

InterpolatorBridge::InterpolatorBridge(const Interpolator& inInnerInterpolator){
	mInterpolatorPtr = inInnerInterpolator.clone();
}

InterpolatorBridge& InterpolatorBridge::operator=(const InterpolatorBridge& inOriginal){

	if (this != &inOriginal){
		mInterpolatorPtr = inOriginal.mInterpolatorPtr->clone();
	}
	return *this;
}

double InterpolatorBridge::calcY(double inX){
	return mInterpolatorPtr->calcY(inX);
}