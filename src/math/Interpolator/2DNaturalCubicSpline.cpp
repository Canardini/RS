#include<Math\Interpolator\2DNaturalCubicSpline.h>


TwoDNaturalCubicSpline::TwoDNaturalCubicSpline(const std::vector<double> & inX, const std::vector<double> &inY, const std::vector<std::vector<double>> & inZ)
:TwoDInterpolator(inX, inY, inZ)
{
	setSplineData();
}


void TwoDNaturalCubicSpline::setSplineData(){
	mInterpolator.resize(0);
	for (size_t i = 0; i < mX.size(); i++){
		NaturalCubicSpline lCubicSpline(mY, mZ[i]);
		mInterpolator.push_back(lCubicSpline);
	}
}

double TwoDNaturalCubicSpline::calcZ(double inX, double inY){
	std::vector<double> lValY(mX.size());

	for (size_t i = 0; i < mX.size(); i++){
		lValY[i] = mInterpolator[i].calcY(inY);
	}

	NaturalCubicSpline lSpline(mX, lValY);

	return lSpline.calcY(inX);
}

boost::shared_ptr<TwoDInterpolator> TwoDNaturalCubicSpline::clone() const{

	boost::shared_ptr<TwoDNaturalCubicSpline> lPtr = boost::make_shared<TwoDNaturalCubicSpline>(*this);
	return lPtr;
}
