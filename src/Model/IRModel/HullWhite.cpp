#include<Model\IRModel\HullWhite.h>

HullWhite::HullWhite(const DiscountCurve & inDiscountCurve,
	const DiscountCurve & inForwardCurve,
	std::vector<double> inMarkovianFactors,
	double inKappa, const std::vector<double> & inVols,
	const std::vector<double> & inVolTimeGrid) :
	MarkovianIRModel(inDiscountCurve, inForwardCurve, inMarkovianFactors), mKappa(inKappa),
	mVols(inVols), mVolTimeGrid(inVolTimeGrid)
{
}

boost::shared_ptr<MarkovianIRModel> HullWhite::clone() const{

	boost::shared_ptr<HullWhite> lPtr = boost::make_shared<HullWhite>(*this);
	return lPtr;
}


double HullWhite::getVolPart(double inValueDate, double inMaturity){

	if (inValueDate >= inMaturity){
		return 0.0;
	}
	double lVol(0);
	int lSize = mVolTimeGrid.size() - 1;

	int lt = locate(mVolTimeGrid, inValueDate);
	int lT = locate(mVolTimeGrid, inMaturity);

	/// inValueDate>inMaturity
	if ((lt == lSize)&(lT == lSize)){
		lVol = mVols[lT] * mVols[lT] * getfundamentalVols(inValueDate, inMaturity, inMaturity);
		return lVol;
	}

	if ((lt == -1)&(lT == -1)){
		lVol = mVols[0] * mVols[0] * getfundamentalVols(inValueDate, inMaturity, inMaturity);
		return lVol;
	}

	lVol += mVols[lt + 1] * mVols[lt + 1] * getfundamentalVols(inValueDate, mVolTimeGrid[lt + 1], inMaturity);

	if (lT >= lSize){
		lVol += mVols[lT] * mVols[lT] * getfundamentalVols(mVolTimeGrid[lT], inMaturity, inMaturity);
	}
	else{
		lVol += mVols[lT+1] * mVols[lT+1] * getfundamentalVols(mVolTimeGrid[lT], inMaturity, inMaturity);
	}


	for (int i = lt + 2; i <= lT;i++){
		lVol += mVols[i] * mVols[i] * getfundamentalVols(mVolTimeGrid[i-1], mVolTimeGrid[i], inMaturity);
	}
	return lVol;
}


double HullWhite::getdiffVol(double inMaturity){
	double eps(0.01);
	double lVtm = getVolPart(0, inMaturity - eps);
	double lVtp = getVolPart(0, inMaturity + eps);

	return (lVtp - lVtm) / (2 * eps);
}


double HullWhite::getfundamentalVols(double l, double u,double lMaturity){
	double lRes = exp(-2 * mKappa*lMaturity)*(exp(mKappa*u) - exp(mKappa*l))
		*(exp(mKappa*u) + exp(mKappa*l) - 4 * exp(mKappa*lMaturity)) + 2 * mKappa*(u - l);
	return lRes / (2 * mKappa*mKappa*mKappa);
}


double HullWhite::getInitSource(double inDate){
	double lf0t=mDiscountCurve->getInstFWrate(inDate);

	double ldiffV0t = getdiffVol(inDate);

	return 0.5*ldiffV0t + lf0t;
}


double HullWhite::getDiscountValue(double inValueDate, double inMaturity){
	/*if (inValueDate >= inMaturity){
		mGDiscount = 0;
		return 1.0;
	}*/
	double lPT = mDiscountCurve->getDFbyDouble(inMaturity);
	double lPt = mDiscountCurve->getDFbyDouble(inValueDate);
	mGDiscount = getGtT(inValueDate, inMaturity);
	
	double lVtT = getVolPart(inValueDate, inMaturity);
	double lV0T = getVolPart(0, inMaturity);
	double lV0t = getVolPart(0, inValueDate);

	return (lPT / lPt)*exp(0.5*(lVtT - lV0T + lV0t) - mGDiscount*mMarkovianFactors[0]);
}

double HullWhite::getFowardValue(double inValueDate, double inT1, double inT2){
	if (inValueDate >= inT1){
		mGFWStart = 0;
		mGFWEnd = 0;

		return 0.0;
	}
	// If FWCurve!=DFCurve, spread are constant.
	double lFWOIS = mDiscountCurve->getFWbyDouble(inT1, inT2);
	mFWDFSpread = mForwardCurve->getFWbyDouble(inT1, inT2) - lFWOIS;
	mIndexAcc = inT2 - inT1;
	double lPT1 = mDiscountCurve->getDFbyDouble(inT1);
	double lPT2 = mDiscountCurve->getDFbyDouble(inT2);
	mGFWStart = getGtT(inValueDate, inT1);
	mGFWEnd = getGtT(inValueDate, inT2);

	double lVtT1 = getVolPart(inValueDate, inT1);
	double lVtT2 = getVolPart(inValueDate, inT2);

	double lV0T1 = getVolPart(0, inT1);
	double lV0T2 = getVolPart(0, inT2);

	// redundant : lpT1/lPT2=1+deltaOIS
	double lRes = (1 / mIndexAcc)*((1 + mIndexAcc*lFWOIS)*exp(-(mGFWStart - mGFWEnd)*mMarkovianFactors[0] + 0.5*(lVtT1 - lVtT2 - lV0T1 + lV0T2))-1.0);
	return lRes + mFWDFSpread;
}



double HullWhite::getDiscountValueByDate(double inValueDate, bgreg::date & inMaturity){
	double lPT = mDiscountCurve->getDFbyDate(inMaturity);
	double lPt = 1.0;

	if (inValueDate > 0.0){
		lPt = mDiscountCurve->getDFbyDouble(inValueDate);
	}
	double lMaturity = mDiscountCurve->getDoubleDate(inMaturity);
	
	if (inValueDate >= lMaturity){
		return (lPT / lPt);
	}
	
	
	mGDiscount = getGtT(inValueDate, lMaturity);

	double lVtT = getVolPart(inValueDate, lMaturity);
	double lV0T = getVolPart(0, lMaturity);
	double lV0t = getVolPart(0, inValueDate);
	double lAdjuster = (lVtT - lV0T + lV0t);
	 
	return (lPT / lPt)*exp(0.5*lAdjuster - mGDiscount*mMarkovianFactors[0]);
}

double HullWhite::getFowardValueByDate(double inValueDate, bgreg::date & inT1){
	setIndexEnd(inT1);
	// If FWCurve!=DFCurve, spread are constant.
	double lFWOIS = mDiscountCurve->getFWbyDate(inT1);
	double lFWLIBOR = mForwardCurve->getFWbyDate(inT1);
	mFWDFSpread = lFWLIBOR - lFWOIS;
	mIndexAcc = mDiscountCurve->getAccDate(inT1, mIndexEndDate);
	mGFWStart = getGtT(inValueDate, inT1);
	mGFWEnd = getGtT(inValueDate, mIndexEndDate);

	double lT1 = mDiscountCurve->getDoubleDate(inT1);
	if (inValueDate >= lT1){
		return lFWLIBOR;
	}
	/*if (inValueDate >= lT1){
		mGFWStart = 0;
		mGFWEnd = 0;
		return 0.0;
	}*/
	double lT2 = lT1+mIndexAcc;

	double lVtT1 = getVolPart(inValueDate, lT1);
	double lVtT2 = getVolPart(inValueDate, lT2);

	double lV0T1 = getVolPart(0, lT1);
	double lV0T2 = getVolPart(0, lT2);

	// redundant : lpT1/lPT2=1+deltaOIS
	double lAdjuster = (lVtT1 - lVtT2 - lV0T1 + lV0T2);
	if (abs(lAdjuster) > 0.1){
		double temp = 12;
	}
	
	double lRes = (1 / mIndexAcc)*((1 + mIndexAcc*lFWOIS)*exp(-(mGFWStart - mGFWEnd)*mMarkovianFactors[0] + 0.5*lAdjuster) - 1.0);
	return lRes + mFWDFSpread;
}

double  HullWhite::getFWDFdiff(double inOIS, double inDelta, double inG1, double inG2){
	return (1 / inDelta)*((1 + inDelta*inOIS)*(inG2 - inG1));
}

double HullWhite::getGtT(const bgreg::date & inValueDate, const bgreg::date & inMaturity){
	if ((mKappa == 0) | (inValueDate > inMaturity)){
		return 0;
	}
	double lPeriod = DCC::Actual365Fixed::dayCountFraction(inValueDate, inMaturity); // a changer et checker l impact

	return (1 / mKappa)*(1 - exp(-mKappa*(lPeriod)));
}

double HullWhite::getGtT(double inValueDate, const bgreg::date & inMaturity){
	double lT = mDiscountCurve->getDoubleDate(inMaturity);
	if ((mKappa == 0) | (inValueDate > lT)){
		return 0;
	}
	double lPeriod = lT - inValueDate;
	return (1 / mKappa)*(1 - exp(-mKappa*(lPeriod)));
}

double HullWhite::getGtT(double inValueDate, double inMaturity){
	if ((mKappa == 0) | (inValueDate > inMaturity)){
		return 0;
	}
	return (1 / mKappa)*(1 - exp(-mKappa*(inMaturity - inValueDate)));
}

double  HullWhite::updateForwardValue(double delta, double FW, double OISspread, std::vector<double> &dx, double G1, double G2){
	 
	double res;
	if (delta == 0){
		res = FW;
	}
	else{
		res = (1 / delta)*((1 + delta*FW)*exp(- dx[0]*(G1 - G2)) - 1.0) + OISspread;
	}
	return res;
}

double HullWhite::updateDiscountValue(double DF, std::vector<double> &dx, double G){
	return DF*exp(-G*(dx[0]));
}

double HullWhite::integralGG(double inTa, double inTb, double inTk, double inTj){ // Integral GtTGtTdt
	double res(0);

	/*res = inTb - inTa;
	res += (exp(-mKappa*(inTj - inTa)) - exp(-mKappa*(inTj - inTb))) / mKappa;
	res += (exp(-mKappa*(inTk - inTa)) - exp(-mKappa*(inTk - inTb))) / mKappa;
	res += (exp(-mKappa*(inTk + inTj - 2 * inTb)) - exp(-mKappa*(inTk + inTj - 2 * inTa))) / (2 * mKappa);
	res /= (mKappa*mKappa);
*/
	return res;

}