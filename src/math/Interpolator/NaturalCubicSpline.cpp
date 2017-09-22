#include<Math\Interpolator\NaturalCubicSpline.h>




NaturalCubicSpline::NaturalCubicSpline(const std::vector<double> & inX, const std::vector <double> & inY)
	:Interpolator(inX, inY)
{
	setSplineData();
}

boost::shared_ptr<Interpolator> NaturalCubicSpline::clone() const{

	boost::shared_ptr<NaturalCubicSpline> lPtr = boost::make_shared<NaturalCubicSpline>(*this);
	return lPtr;
}



void NaturalCubicSpline::setSplineData(){
	mCoeff.resize(mNum);

	mCoeff[0] = 0;
	mCoeff[mNum - 1] = 0;

	std::vector<double> u(mNum - 1);
	u[0] = 0;

	double sig, p;

	for (int i = 1; i < mNum - 1; i++){
		sig = (mX[i] - mX[i - 1]) / (mX[i + 1] - mX[i - 1]);
		p = sig * mCoeff[i - 1] + 2;
		mCoeff[i] = (sig - 1) / p;
		u[i] = (mY[i + 1] - mY[i]) / (mX[i + 1] - mX[i]) - (mY[i] - mY[i - 1]) / (mX[i] - mX[i - 1]);
		u[i] = (6.0 * u[i] / (mX[i + 1] - mX[i - 1]) - sig * u[i - 1]) / p;
	}

	for (int i = 0; i < mNum - 1; i++){
		mCoeff[mNum - 2 - i] = mCoeff[mNum - 2 - i] * mCoeff[mNum - 2 - i + 1] + u[mNum - 2 - i];
	}
}


double NaturalCubicSpline::calcY(double inX){

	int k;
	double h, a, b;
	int klo = 0;
	int khi = mNum - 1;
	double res;
	do{
		k = (khi + klo) / 2;
		if (inX<mX[k]){
			khi = k;
		}
		else{
			klo = k;
		}
	} while ((khi - klo)>1);

	h = mX[khi] - mX[klo];
	a = (mX[khi] - inX) / h;
	b = (inX - mX[klo]) / h;
	res = a * mY[klo] + b * mY[khi] + ((a * a * a - a) * mCoeff[klo] + (b * b * b - b) * mCoeff[khi]) * (h * h) / 6.0;

	return res;
}
