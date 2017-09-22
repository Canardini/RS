#include<Math\Interpolator\3DInterpolatorBridge.h>

ThreeDInterpolatorBridge::ThreeDInterpolatorBridge(const ThreeDInterpolatorBridge& inOriginal){
	mInterpolatorPtr = inOriginal.mInterpolatorPtr->clone();
}

ThreeDInterpolatorBridge::ThreeDInterpolatorBridge(const ThreeDInterpolator& inInnerInterpolator){
	mInterpolatorPtr = inInnerInterpolator.clone();
}

ThreeDInterpolatorBridge& ThreeDInterpolatorBridge::operator=(const ThreeDInterpolatorBridge& inOriginal){

	if (this != &inOriginal){
		mInterpolatorPtr = inOriginal.mInterpolatorPtr->clone();
	}
	return *this;
}

double ThreeDInterpolatorBridge::calcXYZ(double inX, double inY,double inZ){
	return mInterpolatorPtr->calcXYZ(inX, inY,inZ);
}