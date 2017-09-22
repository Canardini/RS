#include<Math\Interpolator\LinearInterpolation.h>




LinearInterpolator::LinearInterpolator(const std::vector<double> & inX, const std::vector <double> & inY)
:Interpolator(inX, inY)
{
}

boost::shared_ptr<Interpolator> LinearInterpolator::clone() const{

	boost::shared_ptr<LinearInterpolator> lPtr = boost::make_shared<LinearInterpolator>(*this);
	return lPtr;
}


double LinearInterpolator::calcY(double inX){

	double res;
	if ((mNum == 1) || (inX <= mX[0])){
		return mY[0];
	}
	else if (mX[mNum-1] <= inX){
		return mY[mNum - 1];
	}
	int left, right, mID;
	left = 0;
	right = mNum - 1;
	mID = 0;

	do{
		mID = (left + right) / 2;
		if ((mX[mID] < inX) &(inX <= mX[mID + 1])) {
			break;
		}
		else if ((mX[mID - 1] < inX) & (inX <= mX[mID])){
			mID = mID - 1;
			break;
		}
	 
		if (inX < mX[mID]){
			right = mID - 1;
		}
		else{
			left = mID + 1;
		}
	}while (left <= right);

	res = (inX - mX[mID]) * (mY[mID + 1] - mY[mID]) / (mX[mID + 1] - mX[mID]) + mY[mID];
	return res;
}
