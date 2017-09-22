#include<Math\Interpolator\3DNaturalCubicSpline.h>


ThreeDNaturalCubicSpline::ThreeDNaturalCubicSpline(const std::vector<double> & inX, const std::vector<double> &inY,
	const std::vector<double> &inZ, const std::vector<std::vector<std::vector<double>>> &inXYZ)
	:ThreeDInterpolator(inX,inY,inZ,inXYZ)
{
	for (size_t i = 0; i < mX.size(); i++){
		mInterpolator.push_back(std::vector<InterpolatorBridge>(mY.size()));
		for (size_t j = 0; j < mY.size(); j++){
			NaturalCubicSpline lCubicSpline(inZ, inXYZ[i][j]);
			mInterpolator[i][j] = lCubicSpline;
		}
	}
}

double ThreeDNaturalCubicSpline::calcXYZ(double inX, double inY, double inZ){
	std::vector<std::vector<double>> lValZ(mX.size(),std::vector<double>(mY.size()));

	for (size_t i = 0; i < mX.size(); i++){
		for (size_t j = 0; j < mY.size(); j++){
			lValZ[i][j] = mInterpolator[i][j].calcY(inZ);
		}
	}

	TwoDNaturalCubicSpline lSpline(mX, mY, lValZ);
	return lSpline.calcZ(inX,inY);
}

boost::shared_ptr<ThreeDInterpolator> ThreeDNaturalCubicSpline::clone() const{

	boost::shared_ptr<ThreeDNaturalCubicSpline> lPtr = boost::make_shared<ThreeDNaturalCubicSpline>(*this);
	return lPtr;
}