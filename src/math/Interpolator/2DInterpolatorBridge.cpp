#include<Math\Interpolator\2DInterpolatorBridge.h>

TwoDInterpolatorBridge::TwoDInterpolatorBridge(const TwoDInterpolatorBridge& inOriginal){
	mInterpolatorPtr = inOriginal.mInterpolatorPtr->clone();
}

TwoDInterpolatorBridge::TwoDInterpolatorBridge(const TwoDInterpolator& inInnerInterpolator){
	mInterpolatorPtr = inInnerInterpolator.clone();
}

TwoDInterpolatorBridge& TwoDInterpolatorBridge::operator=(const TwoDInterpolatorBridge& inOriginal){

	if (this != &inOriginal){
		mInterpolatorPtr = inOriginal.mInterpolatorPtr->clone();
	}
	return *this;
}

double TwoDInterpolatorBridge::calcZ(double inX, double inY){
	return mInterpolatorPtr->calcZ(inX,inY);
}

void TwoDInterpolatorBridge::findNeighboursPoints(std::vector<int> & outNeighbours, int &Num, double inX, double inY){
	mInterpolatorPtr->findNeighboursPoints(outNeighbours, Num, inX, inY);
}