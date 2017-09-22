#include<Model\IRModel\OneFactorCheyette.h>

OneFactorCheyette::OneFactorCheyette(const DiscountCurve & inDiscountCurve,
	const DiscountCurve & inForwardCurve,
	std::vector<double> inMarkovianFactors,
	double inKappa):
	MarkovianIRModel(inDiscountCurve, inForwardCurve, inMarkovianFactors), mKappa(inKappa)
{
}



boost::shared_ptr<MarkovianIRModel> OneFactorCheyette::clone() const{

	boost::shared_ptr<OneFactorCheyette> lPtr = boost::make_shared<OneFactorCheyette>(*this);
	return lPtr;
}

double OneFactorCheyette::getDiscountValue(double inValueDate, double inMaturity){

	double lPT = mDiscountCurve->getDFbyDouble(inMaturity);
	double lPt = mDiscountCurve->getDFbyDouble(inValueDate);
	mGDiscount = getGtT(inValueDate, inMaturity);

	return (lPT / lPt)*exp(-mMarkovianFactors[0] * mGDiscount - 0.5*mGDiscount*mGDiscount*mMarkovianFactors[1]);
}
double OneFactorCheyette::getFowardValue(double inValueDate, double inT1,double inT2){

	// If FWCurve!=DFCurve, spread are constant.
	mFWDFSpread = mForwardCurve->getFWbyDouble(inT1, inT2) - mDiscountCurve->getFWbyDouble(inT1, inT2);
	mIndexAcc = inT2 - inT1;
	double lPT1 = mDiscountCurve->getDFbyDouble(inT1);
	double lPT2 = mDiscountCurve->getDFbyDouble(inT2);
	mGFWStart= getGtT(inValueDate, inT1);
	mGFWEnd = getGtT(inValueDate, inT2);

	// redundant : lpT1/lPT2=1+deltaOIS
	double lRes = (1 / mIndexAcc)*((lPT1 / lPT2)*exp(-(mGFWStart - mGFWEnd)*mMarkovianFactors[0] - 0.5*(mGFWStart*mGFWStart - mGFWEnd *mGFWEnd)*mMarkovianFactors[1]) - 1.0);
	return lRes + mFWDFSpread;
}

double OneFactorCheyette::getOldFowardValueByDate(double inValueDate, bgreg::date & inT1){
	setIndexEnd(inT1);
	double lT1 = mForwardCurve->getDoubleStartDate(inT1);
	double lT2 = mForwardCurve->getDoubleStartDate(mIndexEndDate);// je fais le crevard, c est bien StartDate;
	return getFowardValue(inValueDate, lT1, lT2);
}



double OneFactorCheyette::getOldDiscountValueByDate(double inValueDate, bgreg::date & inMaturity){
	if (inValueDate != 0.0){
		double lMaturity = mDiscountCurve->getDoubleDate(inMaturity);
		return getDiscountValue(inValueDate, lMaturity);
	}
	return mDiscountCurve->getDoubleDate(inMaturity);
}

double OneFactorCheyette::getDiscountValueByDate(double inValueDate, bgreg::date & inMaturity){
	double lPT = mDiscountCurve->getDFbyDate(inMaturity);
	double lPt = 1.0;
	if (inValueDate > 0.0){
		lPt = mDiscountCurve->getDFbyDouble(inValueDate);
	}

	mGDiscount = getGtT(inValueDate, inMaturity);

	return (lPT / lPt)*exp(-mMarkovianFactors[0] * mGDiscount - 0.5*mGDiscount*mGDiscount*mMarkovianFactors[1]);

}

double OneFactorCheyette::getFowardValueByDate(double inValueDate, bgreg::date & inT1){
	setIndexEnd(inT1);
	// If FWCurve!=DFCurve, spread are constant.
	double lFWOIS = mDiscountCurve->getFWbyDate(inT1);
	double lFWLIBOR = mForwardCurve->getFWbyDate(inT1);
	mFWDFSpread = lFWLIBOR - lFWOIS;
	mIndexAcc = mDiscountCurve->getAccDate(inT1, mIndexEndDate);
	mGFWStart = getGtT(inValueDate, inT1);
	mGFWEnd = getGtT(inValueDate, mIndexEndDate);

	// redundant : lpT1/lPT2=1+deltaOIS
	double lRes = (1 / mIndexAcc)*((1 + mIndexAcc*lFWOIS)*exp(-(mGFWStart - mGFWEnd)*mMarkovianFactors[0] - 0.5*(mGFWStart*mGFWStart - mGFWEnd *mGFWEnd)*mMarkovianFactors[1]) - 1.0);
	return lRes + mFWDFSpread;
}

double  OneFactorCheyette::getFWDFdiff(double inOIS, double inDelta, double inG1, double inG2){
	return (1 / inDelta)*((1 + inDelta*inOIS)*(inG2 - inG1));
}

double OneFactorCheyette::getGtT(const bgreg::date & inValueDate, const bgreg::date & inMaturity){
	if ((mKappa == 0 )| (inValueDate > inMaturity)){
		return 0;
	}
	double lPeriod=DCC::Actual365Fixed::dayCountFraction(inValueDate, inMaturity);
	return (1 / mKappa)*(1 - exp(-mKappa*(lPeriod)));
}

double OneFactorCheyette::getGtT(double inValueDate, const bgreg::date & inMaturity){
	double lT = mDiscountCurve->getDoubleDate(inMaturity);
	if ((mKappa == 0) | (inValueDate > lT)){
		return 0;
	}
	double lPeriod = lT-inValueDate;
	return (1 / mKappa)*(1 - exp(-mKappa*(lPeriod)));
}

double OneFactorCheyette::getGtT(double inValueDate, double inMaturity){
	if ((mKappa == 0) | (inValueDate > inMaturity)){
		return 0;
	}
	return (1 / mKappa)*(1 - exp(-mKappa*(inMaturity - inValueDate)));
}

double  OneFactorCheyette::updateForwardValue(double delta, double FW, double OISspread, std::vector<double> &dx, double G1, double G2){


	double res;
	if (delta == 0){
		res = FW;
	}
	else{
		res = (1 / delta)*((1 + delta*FW)*exp(-0.5*dx[1]*(G1*G1 - G2*G2) - dx[0]*(G1 - G2)) - 1.0) + OISspread;
	}
	return res;
}

double OneFactorCheyette::updateDiscountValue(double DF, std::vector<double> &dx, double G){
	return DF*exp(-G*(dx[0]) - 0.5*(dx[1])*G*G);
}

double OneFactorCheyette::integralGG(double inTa, double inTb, double inTk, double inTj){ // Integral GtTGtTdt
	double res;

	res = inTb - inTa;
	res += (exp(-mKappa*(inTj - inTa)) - exp(-mKappa*(inTj - inTb))) / mKappa;
	res += (exp(-mKappa*(inTk - inTa)) - exp(-mKappa*(inTk - inTb))) / mKappa;
	res += (exp(-mKappa*(inTk + inTj - 2 * inTb)) - exp(-mKappa*(inTk + inTj - 2 * inTa))) / (2 * mKappa);
	res /= (mKappa*mKappa);

	return res;

}

double  OneFactorCheyette::getDeterministicDrift(double inValueDate){
	if (mOptionalVols.size() == 0){
		return 0.0;
	}
	return getY(inValueDate);
}

double OneFactorCheyette::integralForY(double v, double w, double inValueDate){
	double lRes = exp(-2 * mKappa*(inValueDate - w)) - exp(-2 * mKappa*(inValueDate - v));
	return lRes / (2 * mKappa);
}


double OneFactorCheyette::getY(double inValueDate){

	if (inValueDate <= 0.0){
		return 0.0;
	}

	int lSize = mOptionalVolTime.size() - 1;

	int lt = locate(mOptionalVolTime, inValueDate);
	if (lt == -1){
		return  mOptionalVols[0] * mOptionalVols[0] * integralForY(0, inValueDate, inValueDate);
	}
	double lVol(0);
	if (lt >= lSize){
		lVol += mOptionalVols[lSize] * mOptionalVols[lSize] * integralForY(mOptionalVolTime[lSize], inValueDate, inValueDate);

	}
	else{
		lVol += mOptionalVols[lt + 1] * mOptionalVols[lt + 1] * integralForY(mOptionalVolTime[lt], inValueDate, inValueDate);
	}

	for (int i = 0; i <= lt; i++){
		if (i == 0){
			lVol += mOptionalVols[i] * mOptionalVols[i] * integralForY(0.0, mOptionalVolTime[i], inValueDate);

		}
		else{
			lVol += mOptionalVols[i] * mOptionalVols[i] * integralForY(mOptionalVolTime[i - 1], mOptionalVolTime[i], inValueDate);
		}
	}

	return lVol;
}