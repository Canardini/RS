#include<Math\Interpolator\2DLinearInterpolator.h>


TwoDLinearInterpolator::TwoDLinearInterpolator(const std::vector<double> & inX, const std::vector<double> &inY, const std::vector<std::vector<double>> & inZ)
: TwoDInterpolator(inX, inY, inZ)
{
	setSplineData();
}

void TwoDLinearInterpolator::setSplineData(){
	mInterpolator.resize(0);
	for (size_t i = 0; i < mX.size(); i++){
		LinearInterpolator lCubicSpline(mY, mZ[i]);
		mInterpolator.push_back(lCubicSpline);
	}
}

double TwoDLinearInterpolator::calcZ(double inX, double inY){
	std::vector<double> lValY(mX.size());

	for (size_t i = 0; i < mX.size(); i++){
		lValY[i] = mInterpolator[i].calcY(inY);
	}

	LinearInterpolator lSpline(mX, lValY);
	return lSpline.calcY(inX);
}

boost::shared_ptr<TwoDInterpolator> TwoDLinearInterpolator::clone() const{

	boost::shared_ptr<TwoDLinearInterpolator> lPtr = boost::make_shared<TwoDLinearInterpolator>(*this);
	return lPtr;
}