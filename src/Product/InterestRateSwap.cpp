#include<Product\InterestRateSwap.h>

InterestRateSwap::InterestRateSwap(const MarkovianIRModelBridge & inMarkovianIRModel,
	int inNumberOfStates,
	std::vector<double> &inMarkovianStates,
	double inValueDate,
	const std::vector<IRSLegInfo> & inIRSLegInfo
	)
	: MarkovianProduct(inMarkovianIRModel, inNumberOfStates, inMarkovianStates, inValueDate), mIRSLegInfo(inIRSLegInfo), mNumberOfLegs(mIRSLegInfo.size())
{
	mStartDates.resize(mNumberOfLegs);
	mDoubleStartDates.resize(mNumberOfLegs);
	mEndStrDates.resize(mNumberOfLegs);
	mEndDates.resize(mNumberOfLegs);
	mDoubleEndDates.resize(mNumberOfLegs);
	mLegAccPeriods.resize(mNumberOfLegs);
	mNotionals.resize(mNumberOfLegs);
	mDiscountFactors.resize(mNumberOfLegs);
	mRates.resize(mNumberOfLegs);
	mNextStartDates.resize(mNumberOfLegs);
	mGDiscount.resize(mNumberOfLegs);
	mAccAdjust.resize(mNumberOfLegs);
	setNumberOfFloatingLegs();
	mGFWStart.resize(mNumberOfFloatingLegs);
	mGFWEnd.resize(mNumberOfFloatingLegs);
	mFWDFspreads.resize(mNumberOfFloatingLegs);
	mIndexAcc.resize(mNumberOfFloatingLegs);
	mFWDFRates.resize(mNumberOfFloatingLegs);
	mIndexStartDates.resize(mNumberOfFloatingLegs);
	mIndexEndDates.resize(mNumberOfFloatingLegs);
	mLegPVs.resize(mNumberOfLegs);
	setConstantVariables();
	if (mIRSLegInfo[0].AmmortizationType != "None"){
		setAccretingNotional();
	}
}

boost::shared_ptr<MarkovianProduct> InterestRateSwap::clone() const{

	boost::shared_ptr<InterestRateSwap> lPtr = boost::make_shared<InterestRateSwap>(*this);
	return lPtr;
}

void InterestRateSwap::changeStructure(const bgreg::date & inDate, double inValueDate){
	setValueDate(inValueDate);
	setForwardDate(inDate);
	for (int i = 0; i < mNumberOfLegs; i++){
		setNextStartDates(i);
		setAccAdjust(i);
	}
}

double InterestRateSwap::getLegPV(int inIndex){
	double lRes = 0.0;
	int lIndex = mNextStartDates[inIndex];
	for (size_t j = lIndex; j < mStartDates[inIndex].size(); j++){
		lRes += mNotionals[inIndex][j] * mDiscountFactors[inIndex][j] * mLegAccPeriods[inIndex][j] * mRates[inIndex][j];
	}

	if (mIRSLegInfo[inIndex].ForwardStartDate>mStartDates[inIndex][0]){
		lRes += mNotionals[inIndex][lIndex - 1] * mDiscountFactors[inIndex][lIndex - 1] * mAccAdjust[inIndex] * mRates[inIndex][lIndex-1];
	}
	return lRes;
}

void InterestRateSwap::setInitNotional(int inIndex){
	int n = mStartDates[inIndex].size();
	mNotionals[inIndex].resize(n);
	int lCoeff;
	(mIRSLegInfo[inIndex].isPayer ? lCoeff = -1 : lCoeff = 1);
	for (int i = 0; i < n; i++){
		mNotionals[inIndex][i] = lCoeff*mIRSLegInfo[inIndex].NotionalByDef;
	}
}

void InterestRateSwap::setAccretingNotional(double inRate){
	int lFixedIndex, lFloatIndex;//Asssuming two legs;
	if (mIRSLegInfo[0].isFixedOrFloat){
		lFixedIndex = 0;
		lFloatIndex = 1;
	}
	else{
		lFixedIndex = 1;
		lFloatIndex = 0;
	}

	double lRate = inRate;
	if (lRate == 0){
		lRate = mIRSLegInfo[lFixedIndex].FixedRate;
	}
	int n = mStartDates[lFixedIndex].size();
	for (int j = 0; j < n; j++){
		if (j != 0){
			mNotionals[lFixedIndex][j] = mNotionals[lFixedIndex][j-1] * (1 + mLegAccPeriods[lFixedIndex][j] * lRate);
		}
	}
	n = mStartDates[lFloatIndex].size();
	int lCompteur(0);
	int j(0);
	while (j < n){
		if (mEndDates[lFloatIndex][j] <= mEndDates[lFixedIndex][lCompteur]){
			mNotionals[lFloatIndex][j] = -mNotionals[lFixedIndex][lCompteur];
			j++;
		}
		else{
			lCompteur++;
		}
	}
}

void InterestRateSwap::setForwardDate(int inIndex){
	for (int i = 0; i < mNumberOfLegs; i++){
		mIRSLegInfo[i].ForwardStartDate = mIRSLegInfo[i].ExercisesDates[inIndex];
	}
}

void InterestRateSwap::setForwardDate(const bgreg::date & inDate){
	for (int i = 0; i < mNumberOfLegs; i++){
		mIRSLegInfo[i].ForwardStartDate = inDate;
	}
}

void InterestRateSwap::setConstantVariables(){

	for (int i = 0; i < mNumberOfLegs; i++){
		setInitStartDates(i);
		setInitEndDates(i);
		setInitAccPeriods(i);
		setNextStartDates(i);
		setAccAdjust(i);
		setInitNotional(i);
	}
}

double InterestRateSwap::getPrice(){
	double lPrice = 0;
	double lNumerator=0;
	int sign;
	mAnnuity = 0;
	for (int i = 0; i < mNumberOfLegs; i++){
		setDiscountFactors(i);
		setRates(i);
		mLegPVs[i] = getLegPV(i);
		sign = (mNotionals[i][0] > 0) - (mNotionals[i][0] < 0);
		if (mIRSLegInfo[i].isFixedOrFloat){
			mAnnuity = mLegPVs[i] / (sign*mRates[i][0]);
		}
		else{
			lNumerator += mLegPVs[i] * sign;
		}
		lPrice += mLegPVs[i];
	}
	mSwapRate = lNumerator / mAnnuity;
	return lPrice;
	//return lPrice-mLegPVs[1];

}

double InterestRateSwap::updatePrice(std::vector<double> &inMarkovianStates){
	double lPrice = 0;
	double lLegPrice = 0;
	double lNumerator = 0;
	int sign;
	for (int i = 0; i < mNumberOfLegs; i++){
		lLegPrice = updateIndexPrice(inMarkovianStates, i);
		sign = (mNotionals[i][0] > 0) - (mNotionals[i][0] < 0);

		if (mIRSLegInfo[i].isFixedOrFloat){
			mAnnuity = lLegPrice / (sign*mRates[i][0]);
		}
		else{
			lNumerator += lLegPrice*sign;
		}
		lPrice += lLegPrice;
	}
	mSwapRate = lNumerator / mAnnuity;

	return lPrice;
	//return lPrice-lLegPrice;

}

void InterestRateSwap::setNextStartDates(int inIndex){
	int n = mStartDates[inIndex].size();
	mNextStartDates[inIndex] = 0;
	for (int i = 0; i < n; i++){
		if (mStartDates[inIndex][i]>=mIRSLegInfo[inIndex].ForwardStartDate){
			mNextStartDates[inIndex] = i;
			break;
		}
		if (i == n - 1){
			mNextStartDates[inIndex] = n;
		}
	}
	
}

void InterestRateSwap::setAccAdjust(int inIndex){
	if (mNextStartDates[inIndex]>0){
		mAccAdjust[inIndex] = DCC::dayCountFraction(mIRSLegInfo[inIndex].ForwardStartDate, mEndDates[inIndex][mNextStartDates[inIndex]-1], mIRSLegInfo[inIndex].Daycount);
	}
	else{
		mAccAdjust[inIndex] = 0.0;
	}
}

void InterestRateSwap::setNumberOfFloatingLegs(){
	int lRes = 0;
	for (int i = 0; i < mNumberOfLegs; i++){
		if (mIRSLegInfo[i].isFixedOrFloat == false){
			lRes++;
		}
	}
	mNumberOfFloatingLegs = lRes;
}

double InterestRateSwap::updateIndexPrice(std::vector<double> &inMarkovianStates, int inIndex){

	std::vector<double> dx(inMarkovianStates.size());
	for (size_t i = 0; i < dx.size(); i++){
		dx[i] = inMarkovianStates[i] - mMarkovianFactors[i];
	}
	double lRes = 0.0;
	double lDF;
	int lIndex = mNextStartDates[inIndex];
	if (mIRSLegInfo[inIndex].isFixedOrFloat){

		for (size_t j = lIndex; j < mStartDates[inIndex].size(); j++){
			lDF = mMarkovianIRModel.updateDiscountValue(mDiscountFactors[inIndex][j], dx, mGDiscount[inIndex][j]);
			lRes += mNotionals[inIndex][j] * lDF* mLegAccPeriods[inIndex][j] * mRates[inIndex][j];
		}

		if (lIndex>0){
			lDF = mMarkovianIRModel.updateDiscountValue(mDiscountFactors[inIndex][lIndex - 1], dx, mGDiscount[inIndex][lIndex - 1]);
			lRes += mNotionals[inIndex][lIndex - 1] * lDF* mAccAdjust[inIndex] * mRates[inIndex][lIndex - 1];
		}
	}
	else{
		double lFW;
		for (size_t j = lIndex; j < mStartDates[inIndex].size(); j++){
			lDF = mMarkovianIRModel.updateDiscountValue(mDiscountFactors[inIndex][j], dx, mGDiscount[inIndex][j]);
			lFW = mMarkovianIRModel.updateForwardValue(mIndexAcc[0][j], mFWDFRates[0][j], mFWDFspreads[0][j], dx, mGFWStart[0][j], mGFWEnd[0][j]);
			lRes += mNotionals[inIndex][j] * lDF * mLegAccPeriods[inIndex][j] * lFW;
		}
		if (lIndex>0){
			lDF = mMarkovianIRModel.updateDiscountValue(mDiscountFactors[inIndex][lIndex - 1], dx, mGDiscount[inIndex][lIndex - 1]);
			lFW = mMarkovianIRModel.updateForwardValue(mIndexAcc[0][lIndex - 1], mFWDFRates[0][lIndex - 1], mFWDFspreads[0][lIndex - 1], dx, mGFWStart[0][lIndex - 1], mGFWEnd[0][lIndex - 1]);
			lRes += mNotionals[inIndex][lIndex - 1] * lDF* mAccAdjust[inIndex] * lFW;
		}
	}

	return lRes;
}

void InterestRateSwap::setInitStartDates(int inIndex){
	mStartDates[inIndex].resize(0);
	mDoubleStartDates[inIndex].resize(0);

	mStartDates[inIndex].push_back(mIRSLegInfo[inIndex].SpotDate);
	
	mDoubleStartDates[inIndex].push_back(DCC::Actual365Fixed::dayCountFraction(mIRSLegInfo[inIndex].TodayDate, mIRSLegInfo[inIndex].SpotDate));
	bgreg::date lNextDate = Calendar::getRespondDate(mStartDates[inIndex][0], mIRSLegInfo[inIndex].Interval, mIRSLegInfo[inIndex].BusinessCalendar, mIRSLegInfo[inIndex].SettlementAdjust);
	
	int lFlag = ((mIRSLegInfo[inIndex].TerminationDate) - lNextDate).days();
	int i=1;
	std::string lTenorType;
	lTenorType = Calendar::getTenorType(mIRSLegInfo[inIndex].Interval);
	std::string lUpdatedTenor;
	int lTenor = Calendar::getTenor(mIRSLegInfo[inIndex].Interval);
	while (true){
		lUpdatedTenor = boost::lexical_cast<std::string>(lTenor*i) + lTenorType;
		lNextDate = Calendar::getRespondDate(mStartDates[inIndex][0], lUpdatedTenor, mIRSLegInfo[inIndex].BusinessCalendar, mIRSLegInfo[inIndex].SettlementAdjust);
		lFlag = ((mIRSLegInfo[inIndex].TerminationDate) - lNextDate).days();
		if (lFlag <= 0){
			break;
		}
		mStartDates[inIndex].push_back(lNextDate);
		mDoubleStartDates[inIndex].push_back(DCC::Actual365Fixed::dayCountFraction(mIRSLegInfo[inIndex].TodayDate,lNextDate));
		i++;
	}
}

void InterestRateSwap::setInitEndDates(int inIndex){
	mEndDates[inIndex].resize(0);
	mEndStrDates[inIndex].resize(0);
	mDoubleEndDates[inIndex].resize(0);
	bgreg::date lNextDate;
	int lFlag ;
	int i = 1;
	std::string lTenorType;
	lTenorType = Calendar::getTenorType(mIRSLegInfo[inIndex].Interval);
	int lTenor = Calendar::getTenor(mIRSLegInfo[inIndex].Interval);
	std::string lUpdatedTenor;
	while (true){
		lUpdatedTenor = boost::lexical_cast<std::string>(lTenor*i)+lTenorType;
		lNextDate = Calendar::getRespondDate(mStartDates[inIndex][0], lUpdatedTenor, mIRSLegInfo[inIndex].BusinessCalendar, mIRSLegInfo[inIndex].SettlementAdjust);
		lFlag = ((mIRSLegInfo[inIndex].TerminationDate) - lNextDate).days();
		if (lFlag < 0){
			break;
		}
		mEndDates[inIndex].push_back(lNextDate);
		mEndStrDates[inIndex].push_back(to_iso_extended_string(lNextDate));
		mDoubleEndDates[inIndex].push_back(DCC::Actual365Fixed::dayCountFraction(mIRSLegInfo[inIndex].TodayDate, lNextDate));
		i++;
	}
	if (mEndDates[inIndex][mEndDates[inIndex].size() - 1] != mIRSLegInfo[inIndex].TerminationDate){
		mEndDates[inIndex].push_back(mIRSLegInfo[inIndex].TerminationDate);
		mEndStrDates[inIndex].push_back(to_iso_extended_string(mIRSLegInfo[inIndex].TerminationDate));
		mDoubleEndDates[inIndex].push_back(DCC::Actual365Fixed::dayCountFraction(mIRSLegInfo[inIndex].TodayDate, mIRSLegInfo[inIndex].TerminationDate));
	}
}

void InterestRateSwap::setInitAccPeriods(int inIndex){
	int n = mStartDates[inIndex].size();
	mLegAccPeriods[inIndex].resize(n);
	for (int i = 0; i < n; i++){
		mLegAccPeriods[inIndex][i] = DCC::dayCountFraction(mStartDates[inIndex][i], mEndDates[inIndex][i], mIRSLegInfo[inIndex].Daycount);
	}
}

void InterestRateSwap::setDiscountFactors(int inIndex){
	int n = mEndDates[inIndex].size();
	mDiscountFactors[inIndex].resize(n);
	mGDiscount[inIndex].resize(n);
	for (int i = 0; i < n; i++){
		mDiscountFactors[inIndex][i] = mMarkovianIRModel.getDiscountValueByDate(mValueDate, mEndDates[inIndex][i]);
		mGDiscount[inIndex][i] = mMarkovianIRModel.getGDiscount();
	}
}

void InterestRateSwap::setRates(int inIndex){
	int n = mStartDates[inIndex].size();
	mRates[inIndex].resize(n);
	
	if (mIRSLegInfo[inIndex].isFixedOrFloat == true){
		for (int i = 0; i < n; i++){
			mRates[inIndex][i] = mIRSLegInfo[inIndex].FixedRate;
		}
		mDefaultFixedRate = mRates[inIndex][0];
	}

	else{
		mGFWStart[0].resize(n);// inIndex n est pas forcement le bon index
		mGFWEnd[0].resize(n);
		mFWDFspreads[0].resize(n);
		mFWDFRates[0].resize(n);
		mIndexAcc[0].resize(n);
		mIndexStartDates[0].resize(n);
		mIndexEndDates[0].resize(n); 
		for (int i = 0; i < n; i++){
			mIndexStartDates[0][i] = mStartDates[inIndex][i];
			mRates[inIndex][i] = mMarkovianIRModel.getFowardValueByDate(mValueDate, mStartDates[inIndex][i]);
			mIndexEndDates[0][i] = mMarkovianIRModel.getIndexEnd();
			mGFWStart[0][i] = mMarkovianIRModel.getGFWStart();
			mGFWEnd[0][i] = mMarkovianIRModel.getGFWEnd();
			mFWDFspreads[0][i] = mMarkovianIRModel.getFWDFSpread();
			mFWDFRates[0][i] = mRates[inIndex][i] - mFWDFspreads[0][i];
			mIndexAcc[0][i] = mMarkovianIRModel.getIndexAcc();
		}
	}
}

// integral SWR(t,0,0)^2
double InterestRateSwap::ClosedFormIntegral(const bgreg::date &inValueDate, const bgreg::date &inTa, const bgreg::date &inTb){
	
	int startL, endL;
	int startF, endF;
	double P;
	double M = 0.0;
	double N = 0.0;
	double FWOIS;
	double SWR;
	double DiffFWOIS;

	int i, j, k, l;
	
	double ta = DCC::Actual365Fixed::dayCountFraction(inValueDate, inTa);
	double tb = DCC::Actual365Fixed::dayCountFraction(inValueDate, inTb);
	int lIndexL;
	int lIndexF;

	if (mIRSLegInfo[0].isFixedOrFloat==false){
		lIndexL = 0;
		lIndexF = 1;
		startL = mNextStartDates[0] ;
		endL = mStartDates[0].size()-1;
		startF = mNextStartDates[1];
		endF = mStartDates[1].size() - 1;
	}
	else{
		lIndexL = 1;
		lIndexF = 0;
		startL = mNextStartDates[1];
		endL = mStartDates[1].size() - 1;
		startF = mNextStartDates[0];
		endF = mStartDates[0].size() - 1;
	}
	double integraleNp = 0.0;
	double integraleNpMp = 0.0;
	double integraleMp = 0.0;
	double Ki[6];
	double Ti[6];
	double Kj[6];
	double Tj[6];
	SWR = mSwapRate;

	i = startL;
	if ((mAccAdjust[lIndexL] != 0)){
		Ti[3] = DCC::Actual365Fixed::dayCountFraction(inValueDate, mEndDates[lIndexL][i - 1]);
		Ti[4] = DCC::Actual365Fixed::dayCountFraction(inValueDate, mIndexEndDates[0][i - 1]);
		Ti[5] = DCC::Actual365Fixed::dayCountFraction(inValueDate, mIndexStartDates[0][i - 1]);
		FWOIS = mFWDFRates[0][i - 1];
		P = mDiscountFactors[lIndexL][i - 1];
		Ki[3] = (mAccAdjust[lIndexL])*(P)*(FWOIS + mFWDFspreads[0][i - 1]);
		N -= Ki[3];
		DiffFWOIS = mMarkovianIRModel.getFWDFdiff(FWOIS, mIndexAcc[0][i - 1], 0.0, 1.0);
		Ki[4] = -P*(mAccAdjust[lIndexL])*DiffFWOIS;
		Ki[5] = -Ki[4];
	}

	for (i = startL; i <= endL; i++) {
		FWOIS = mFWDFRates[0][i];
		P = mDiscountFactors[lIndexL][i];
		Ti[0] = DCC::Actual365Fixed::dayCountFraction(inValueDate, mEndDates[lIndexL][i]);
		Ti[1] = DCC::Actual365Fixed::dayCountFraction(inValueDate, mIndexEndDates[0][i]);
		Ti[2] = DCC::Actual365Fixed::dayCountFraction(inValueDate, mIndexStartDates[0][i]);
	
		Ki[0] = -((mLegAccPeriods[lIndexL][i])*(P)*(FWOIS + mFWDFspreads[0][i]));
		N -= Ki[0];

		DiffFWOIS = mMarkovianIRModel.getFWDFdiff(FWOIS, mIndexAcc[0][i], 0.0, 1.0);
		Ki[1] = P*(mLegAccPeriods[lIndexL][i])*DiffFWOIS;
		Ki[2] = -Ki[1];
		if ((i > startL) | (mAccAdjust[lIndexL] == 0)){
			Ki[3] = 0.0;
			Ki[4] = 0.0;
			Ki[5] = 0.0;
		}
		j = startL;
		if ((mAccAdjust[lIndexL] != 0)){
			Tj[3] = DCC::Actual365Fixed::dayCountFraction(inValueDate, mEndDates[lIndexL][j - 1]);
			Tj[4] = DCC::Actual365Fixed::dayCountFraction(inValueDate, mIndexEndDates[0][j - 1]);
			Tj[5] = DCC::Actual365Fixed::dayCountFraction(inValueDate, mIndexStartDates[0][j - 1]);
			FWOIS = mFWDFRates[0][j - 1];
			P = mDiscountFactors[lIndexL][j - 1];
			Kj[3] = (mAccAdjust[lIndexL])*(P)*(FWOIS + mFWDFspreads[0][j - 1]);
			DiffFWOIS = mMarkovianIRModel.getFWDFdiff(FWOIS, mIndexAcc[0][j - 1], 0.0, 1.0);
			Kj[4] = -P*(mAccAdjust[lIndexL])*DiffFWOIS;
			Kj[5] = -Kj[4];
		}

		for (j = startL; j <= endL; j++) {
			Tj[0] = DCC::Actual365Fixed::dayCountFraction(inValueDate, mEndDates[lIndexL][j]);
			Tj[1] = DCC::Actual365Fixed::dayCountFraction(inValueDate, mIndexEndDates[0][j]);
			Tj[2] = DCC::Actual365Fixed::dayCountFraction(inValueDate, mIndexStartDates[0][j]);
			FWOIS = mFWDFRates[0][j];
			P = mDiscountFactors[lIndexL][j];

			Kj[0] = -((mLegAccPeriods[lIndexL][j])*(P)*(FWOIS + mFWDFspreads[0][j]));
			DiffFWOIS = mMarkovianIRModel.getFWDFdiff(FWOIS, mIndexAcc[0][j], 0.0, 1.0);
			Kj[1] = P*(mLegAccPeriods[lIndexL][j])*DiffFWOIS;
			Kj[2] = -Kj[1];

			if ((j > startL) | (mAccAdjust[lIndexL] == 0)){
				Kj[3] = 0.0;
				Kj[4] = 0.0;
				Kj[5] = 0.0;
			}

			for (k = 0; k<6; k++){
				for (l = 0; l<6; l++){
					if ((Ki[k] != 0.0) && (Kj[l] != 0.0)){
						integraleNp += Ki[k] * Kj[l] * mMarkovianIRModel.integralGG(ta, tb, Ti[k], Tj[l]);
					}
				}
			}
		}

		j = startF;
		if ((mAccAdjust[lIndexF] != 0)){
			Tj[1] = DCC::Actual365Fixed::dayCountFraction(inValueDate, mEndDates[lIndexF][j - 1]);
			P = mDiscountFactors[lIndexL][j - 1];;
			Kj[1] = P*(mAccAdjust[lIndexF]);
		}

		for (j = startF; j <= endF; j++) {
			Tj[0] = DCC::Actual365Fixed::dayCountFraction(inValueDate, mEndDates[lIndexF][j]);
			P = mDiscountFactors[lIndexF][j];
			Kj[0] = -P*(mLegAccPeriods[lIndexF][j]);

			if ((j > startF) | (mAccAdjust[lIndexF] == 0)){
				Kj[1] = 0.0;
			}

			for (k = 0; k<6; k++){
				for (l = 0; l<2; l++){
					if ((Ki[k] != 0.0) && (Kj[l] != 0.0)){
						integraleNpMp += Ki[k] * Kj[l] * mMarkovianIRModel.integralGG(ta, tb, Ti[k], Tj[l]);
					}
				}
			}
		}
	}

	M = 0.0;
	i = startF;
	if ((mAccAdjust[lIndexF] != 0)){
		Ti[1] = DCC::Actual365Fixed::dayCountFraction(inValueDate, mEndDates[lIndexF][i - 1]);
		P = mDiscountFactors[lIndexF][i - 1];
		Ki[1] = P*mAccAdjust[lIndexF];
		M -= Ki[1];
	}
	for (i = startF; i <= endF; i++) {
		Ti[0] = DCC::Actual365Fixed::dayCountFraction(inValueDate, mEndDates[lIndexF][i]);
		P = mDiscountFactors[lIndexF][i];
		Ki[0] = -P*(mLegAccPeriods[lIndexF][i]);
		M -= Ki[0];
		if ((i > startF) | (mAccAdjust[lIndexF] == 0)){
			Ki[1] = 0.0;
		}

		j = startF;
		if ((mAccAdjust[lIndexF] != 0)){
			Tj[1] = DCC::Actual365Fixed::dayCountFraction(inValueDate, mEndDates[lIndexF][j - 1]);
			P = mDiscountFactors[lIndexF][j - 1];
			Kj[1] = P*mAccAdjust[lIndexF];
		}
		for (j = startF; j <= endF; j++) {
			Tj[0] = DCC::Actual365Fixed::dayCountFraction(inValueDate, mEndDates[lIndexF][j]);
			P = mDiscountFactors[lIndexF][j];
			Kj[0] = -P*(mLegAccPeriods[lIndexF][j]);

			if ((j > startF) | (mAccAdjust[lIndexF] == 0)){
				Kj[1] = 0.0;
			}
			for (k = 0; k<2; k++){
				for (l = 0; l<2; l++){
					if ((Ki[k] != 0.0) && (Kj[l] != 0.0)){
						integraleMp += Ki[k] * Kj[l] * mMarkovianIRModel.integralGG(ta, tb, Ti[k], Tj[l]);
					}
				}
			}
		}
	}
	return (integraleNp - 2 * SWR*integraleNpMp + SWR*SWR*integraleMp) / (M*M);
}

void InterestRateSwap::swapRateSensitivity(const bgreg::date &inValueDate, bool inBothSensitivities, double & outDx, double & outDxx){
	return swapRateSensitivity(0.0, 0.0, inValueDate, inBothSensitivities,outDx, outDxx);
}
void InterestRateSwap::swapRateSensitivity(double inX, double inY, const bgreg::date &inValueDate, bool inBothSensitivities, double & outDx, double & outDxx){
	int startL, endL;
	int startF, endF;
	int lIndexL;
	int lIndexF;
	std::vector<double> dx(2);
	dx[0] = inX;
	dx[1] = inY;

	double x = inX; //violation of oop, on suppose qu'il ya seulement deux facteurs.
	double y = inY;
	double G;
	double G1, G2;
	double P;
	double N = 0.0;
	double Np = 0.0;
	double Npp = 0.0;
	double M = 0.0;
	double Mp = 0.0;
	double Mpp = 0.0;
	double FWOIS;
	double DiffFWOIS;
	double DiffSQRFWOIS;
	int i;

	if (mIRSLegInfo[0].isFixedOrFloat == false){
		lIndexL = 0;
		lIndexF = 1;
		startL = mNextStartDates[0];
		endL = mStartDates[0].size() - 1;
		startF = mNextStartDates[1];
		endF = mStartDates[1].size() - 1;
	}
	else{
		lIndexL = 1;
		lIndexF = 0;
		startL = mNextStartDates[1];
		endL = mStartDates[1].size() - 1;
		startF = mNextStartDates[0];
		endF = mStartDates[0].size() - 1;
	}
	i = startL;
	if ((mAccAdjust[lIndexL] != 0)){
		G = mMarkovianIRModel.getGtT(inValueDate, mEndDates[lIndexL][i - 1]);
		G1 = mMarkovianIRModel.getGtT(inValueDate, mIndexStartDates[0][i - 1]);
		G2 = mMarkovianIRModel.getGtT(inValueDate, mIndexEndDates[0][i - 1]);
		FWOIS = mMarkovianIRModel.updateForwardValue(mIndexAcc[0][i - 1], mFWDFRates[0][i - 1], mFWDFspreads[0][i - 1], dx, G1, G2);
		P = mMarkovianIRModel.updateDiscountValue(mDiscountFactors[lIndexL][i - 1], dx, G);
		N += (mAccAdjust[lIndexL])*(P)*(FWOIS);
		Np -= G*P*(mAccAdjust[lIndexL])*(FWOIS);
		DiffFWOIS = mMarkovianIRModel.getFWDFdiff(FWOIS - mFWDFspreads[0][i - 1], mIndexAcc[0][i - 1], G1, G2);
		Np += P*(mAccAdjust[lIndexL])*DiffFWOIS;
		if (inBothSensitivities){
			Npp += G*G*P*(mAccAdjust[lIndexL])*(FWOIS);
			Npp -= 2 * G*P*(mAccAdjust[lIndexL])*DiffFWOIS;
			DiffSQRFWOIS = mMarkovianIRModel.getFWDFdiff(DiffFWOIS - 1 / ((mIndexAcc[0][i - 1])), (mIndexAcc[0][i - 1]), G1, G2);
			Npp += P*(mAccAdjust[lIndexL])*DiffSQRFWOIS;
		}
	}

	for (i = startL; i <= endL; i++) {
		G = mMarkovianIRModel.getGtT(inValueDate, mEndDates[lIndexL][i]);
		G1 = mMarkovianIRModel.getGtT(inValueDate, mIndexStartDates[0][i]);
		G2 = mMarkovianIRModel.getGtT(inValueDate, mIndexEndDates[0][i]);
		FWOIS = mMarkovianIRModel.updateForwardValue(mIndexAcc[0][i], mFWDFRates[0][i], mFWDFspreads[0][i], dx, G1, G2);
		P = mMarkovianIRModel.updateDiscountValue(mDiscountFactors[lIndexL][i],dx, G);
		N += (mLegAccPeriods[lIndexL][i])*(P)*(FWOIS);
		Np -= G*P*(mLegAccPeriods[lIndexL][i])*(FWOIS);
		DiffFWOIS = mMarkovianIRModel.getFWDFdiff(FWOIS - mFWDFspreads[0][i], mIndexAcc[0][i], G1, G2);
		Np += P*(mLegAccPeriods[lIndexL][i])*DiffFWOIS;
		if (inBothSensitivities){
			Npp += G*G*P*(mLegAccPeriods[lIndexL][i])*(FWOIS);
			Npp -= 2 * G*P*(mLegAccPeriods[lIndexL][i])*DiffFWOIS;
			DiffSQRFWOIS = mMarkovianIRModel.getFWDFdiff(DiffFWOIS - 1 / ((mIndexAcc[0][i])), (mIndexAcc[0][i]), G1, G2);
			Npp += P*(mLegAccPeriods[lIndexL][i])*DiffSQRFWOIS;
		}
	}
	i = startF;
	if ((mAccAdjust[lIndexF] != 0)){
		G = mMarkovianIRModel.getGtT(inValueDate, mEndDates[lIndexF][i - 1]);
		P = mMarkovianIRModel.updateDiscountValue(mDiscountFactors[lIndexF][i - 1], dx, G);
		M += (mAccAdjust[lIndexF])*P;
		Mp -= (mAccAdjust[lIndexF])*P*G;
		if (inBothSensitivities){
			Mpp += (mAccAdjust[lIndexF])*P*G*G;
		}
	}
	for (i = startF; i <= endF; i++) {
		G = mMarkovianIRModel.getGtT(inValueDate, mEndDates[lIndexF][i]);
		P = mMarkovianIRModel.updateDiscountValue(mDiscountFactors[lIndexF][i], dx, G);

		M += (mLegAccPeriods[lIndexF][i])*P;
		Mp -= (mLegAccPeriods[lIndexF][i])*P*G;
		if (inBothSensitivities){
			Mpp += (mLegAccPeriods[lIndexF][i])*P*G*G;
		}
	}
	if (M != 0.0){
		outDx = (M*Np - N*Mp) / (M*M);
	}
	else{
		outDx = 0.0;
	}

	if (inBothSensitivities){
		if (M != 0.0){
			outDxx = (Npp*M - 2 * Mp*Np - Mpp*N) / (M*M) + (2 * Mp*Mp*N) / (M*M*M);
		}
		else{
			outDxx = 0.0;
		}

	}
}






