#include<Calibration/CheyetteCalibration/DDCalibration.h>
#include <Math/Solver/BandSolver.h>
#include <boost/thread.hpp>

DDCalibration::DDCalibration(const std::vector<double> &inMinVols, 
	const std::vector<double> &inCenteredVols,
	const std::vector<double> &inPlusVols,
	const std::vector<double> &inStrikes,
	const std::vector<bgreg::date> & inExercisesDates,
	const InterestRateSwap &inIRS,
	const TimeGrid &inTimeGrid,
	const bgreg::date & inToday) :
	mCenteredVols(inCenteredVols), mMinVols(inMinVols), mPlusVols(inPlusVols), mStrikes(inStrikes), mExercisesDates(inExercisesDates), mIRSpivot(inIRS),
	mIRS(inIRS), mToday(inToday), mTimeGrid(inTimeGrid)
{
	mNumberOfVols = mCenteredVols.size();
	mLambdas.resize(mNumberOfVols);
	mBetas.resize(mNumberOfVols);
	mSwapRate.resize(mNumberOfVols);
	mCenteredMarketPrice.resize(mNumberOfVols);
	mMinusMarketPrice.resize(mNumberOfVols);
	mPlusMarketPrice.resize(mNumberOfVols);
	mCenteredCalibratedPrice.resize(mNumberOfVols);
	mMinusCalibratedPrice.resize(mNumberOfVols);
	mPlusCalibratedPrice.resize(mNumberOfVols);
	mMarketLambdas.resize(mNumberOfVols);
	mMarketErrors.resize(mNumberOfVols);
	mMarketBetas.resize(mNumberOfVols);
	mNumSteps =(int)(std::ceil(mTimeGrid[mTimeGrid.getGridSize()-1]*10));
	mKappa = 0.001; // to read from somewhere ( just for Improved DD)
	mDoubleExercisesDates.resize(mExercisesDates.size());
	mTempMarkov.resize(2); 
	mAdvancedMarketLambdas.resize(mNumberOfVols);
	mAdvancedMarketBetas.resize(mNumberOfVols);
	mFWPrices.resize(mNumberOfVols);

}

void DDCalibration::runCalibration(){
	runMarketLambdaBetaCalibration();
	runLambdaCalibration();
	runBetaCalibration();
}

void DDCalibration::getLambdas(std::vector<double> &outLambdas){
	int lNout = outLambdas.size();
	if (lNout == mNumberOfVols){
		for (int i = 0; i < lNout; i++){
			outLambdas[i] = mLambdas[i];
		}
	}
}

void DDCalibration::getBetas(std::vector<double> &outBetas){
	int lNout = outBetas.size();
	if (lNout == mNumberOfVols){
		for (int i = 0; i < lNout; i++){
			outBetas[i] = mBetas[i];
		}
	}
}

void DDCalibration::getSWRates(std::vector<double> &outSWRates){
	int lNout = outSWRates.size();
	if (lNout == mNumberOfVols){
		for (int i = 0; i < lNout; i++){
			outSWRates[i] = mSwapRate[i];
		}
	}
}

void DDCalibration::runMarketLambdaBetaCalibration(){
	double lMaturity;
	double lAnnuity;
	std::vector<double> lP(2);
	std::vector<double> lW(3);
	std::vector<double> oLower(2);
	std::vector<double> oUpper(2);
	 
	int inCallOrPut;
	int lIter=0;
	double lfret=0;
	lW[1] = 100.0;
	shiftedLogNormalfunction oLogNormalFunction;
	int N = std::max<std::size_t>(boost::thread::hardware_concurrency(), 1);
	for (int j = 0; j < mNumberOfVols; j++){
		mIRS.changeStructure(mExercisesDates[j], 0);
		mFWPrices[j]=mIRS.getPrice(); // doit etre optimize car on reconstruit la courbe de zero coupon
		mSwapRate[j] = mIRS.getSwapRate();
		lAnnuity = mIRS.getAnnuity();
		lMaturity = DCC::Actual365Fixed::dayCountFraction(mToday, mExercisesDates[j]);
		mDoubleExercisesDates[j] = lMaturity;
		if (mMinVols[j] == 0){
			lW[0] = 0.0;
		}
		else{
			lW[0] = 1.0;
		}
		if (mPlusVols[j] == 0){
			lW[2] = 0.0;
		}
		else{
			lW[2] = 1.0;
		}
		
		
		(mIRS.isPayer() ? inCallOrPut = 1 : inCallOrPut = -1);
		mMinusMarketPrice[j] = Black::black(inCallOrPut, mSwapRate[j], mStrikes[0], lAnnuity, lMaturity, mMinVols[j]);
		mCenteredMarketPrice[j] = Black::black(inCallOrPut, mSwapRate[j], mStrikes[1], lAnnuity, lMaturity, mCenteredVols[j]);
		mPlusMarketPrice[j] = Black::black(inCallOrPut, mSwapRate[j], mStrikes[2], lAnnuity, lMaturity, mPlusVols[j]);
		
		fillSLNStructure(oLogNormalFunction, mSwapRate[j], mStrikes[0], mStrikes[1], mStrikes[2], mMinusMarketPrice[j],
						 mCenteredMarketPrice[j], mPlusMarketPrice[j], lW[0], lW[1], lW[2], lAnnuity, lMaturity, 
						 mMinVols[j], mCenteredVols[j], mPlusVols[j], inCallOrPut);
	
		

		oLower[0] = 1e-12;
		oLower[1] = -0.1;
		oUpper[0] = 1.5;
		oUpper[1] = 1.2;
		DifferentialEvolution2 oDE;
		oDE.setNumberOfThreads(N);
		int res = oDE.Optimize(oLogNormalFunction, lP, lfret, oLower, oUpper);
		double tempVal = oLogNormalFunction(lP);
		//Optimizer::dfpmin(lP, 1.0e-15, lIter, lfret, oLogNormalFunction);
		mMarketErrors[j] = lfret;
		mMarketLambdas[j] = lP[0];
		mMarketBetas[j] = lP[1];

	}
}

void DDCalibration::runLambdaCalibration(){
	double lMaturity;
	double lIntegral;
	double sum = 0.0;
	bgreg::date lTi;
	bgreg::date lTim;
	for (int j = 0; j < mNumberOfVols; j++){
		lTim = mToday;
		mIRS.changeStructure(mExercisesDates[j], 0);
		mIRS.getPrice();
		lMaturity = DCC::Actual365Fixed::dayCountFraction(mToday, mExercisesDates[j]);
		sum = 0.0;
		for (int i = 0; i < j; i++){
			lTi = mExercisesDates[i];
			lIntegral = mIRS.ClosedFormIntegral(mToday, lTim, lTi);
			sum += mLambdas[i] * mLambdas[i] * pow(mSwapRate[i], 2)*lIntegral;
			lTim = lTi;
		}
		lTi = mExercisesDates[j];
		lIntegral = mIRS.ClosedFormIntegral(mToday, lTim, lTi);
		mLambdas[j] = mMarketLambdas[j] * mMarketLambdas[j] * mSwapRate[j] * mSwapRate[j] * lMaturity - sum;

		if ((lIntegral <= 0) | (mLambdas[j] < 0)){
			if (j>0){
				mLambdas[j] = mLambdas[j - 1];
			}
			else{
				mLambdas[j] = 0.0;
			}
		}
		else{
			mLambdas[j] /= pow(mSwapRate[j], 2)*lIntegral;
			mLambdas[j] = sqrt(mLambdas[j]); /// A regle le probleme de 1#indf
		}
		
	}
}

void DDCalibration::runBetaCalibration(){
		double lMaturity;
		int i, j;
		double sum = 0.0;
		double Integral;
		bgreg::date lTi;
		bgreg::date lTim;
		int lNumberOfSteps;
		for (j = 0; j < mNumberOfVols; j++){
			lTim = mToday;
			mIRS.changeStructure(mExercisesDates[j], 0);
			mIRS.getPrice();
			lMaturity = DCC::Actual365Fixed::dayCountFraction(mToday, mExercisesDates[j]);
			sum = 0.0;
			mNumSteps = 0;
			for (i = 0; i < j; i++){
				lTi = mExercisesDates[i];
				mIRSpivot.changeStructure(lTi, 0);
				mIRSpivot.getPrice();
				lNumberOfSteps = (int)(std::ceil(DCC::Actual365Fixed::dayCountFraction(lTim, lTi)*20));
				//std::cout << lNumberOfSteps << "\n";
				
				Integral = RectangularGammaIntegral(true, mSwapRate[i], mLambdas[i], mBetas[i], 1, lTim, lTi, lNumberOfSteps);
				sum = sum + Integral;
				lTim = lTi;
				mNumSteps += lNumberOfSteps;
			}
			
			
			lTi = mExercisesDates[j];
			lNumberOfSteps = (int)(std::ceil(DCC::Actual365Fixed::dayCountFraction(lTim, lTi) * 20));
			//std::cout << lNumberOfSteps << "\n";
			mNumSteps += lNumberOfSteps;
			//std::cout << mNumSteps << "\n";
		

			Integral = RectangularGammaIntegral(false, mSwapRate[j], mLambdas[j], 0.0, 1, lTim, lTi, lNumberOfSteps);
			
			sum = sum + Integral;
			Integral = RectangularGammaIntegral(false, 0.0, 0.0, 0.0, 2, mToday, lTi, mNumSteps);
			mBetas[j] = 2 * mMarketLambdas[j] * mMarketLambdas[j] * mSwapRate[j] * mMarketBetas[j] * Integral - sum;
			Integral = RectangularGammaIntegral(false, mSwapRate[j], mLambdas[j], 0.0, 0, lTim, lTi, lNumberOfSteps);
			
			if (Integral == 0){
				if (j > 0){
					mBetas[j] = mBetas[j-1];
				}
				else{
					mBetas[j] = 0.0;
				}
			}
			else{
				mBetas[j] = mBetas[j] / (Integral);
				if (mBetas[j] < -2){
					mBetas[j] =-2;
				}
				else if (mBetas[j]>2){
					mBetas[j] = 2.0;
				}
			}
				 
		}
}


double DDCalibration::RectangularGammaIntegral(bool inUsePivot,double SWRate, double lambda, double gamma, int flag, const bgreg::date &inT1, const bgreg::date &inT2, int n){
	double sum = 0.0;
	double diffk = 0.0;
	double diffi = 0.0;
	double diffkxx = 0.0;
	double diffsquared = 0.0;

	double T = DCC::Actual365Fixed::dayCountFraction(inT1, inT2);
	int days = (int)(365.0*T / (1.0*n));
	if (days == 0){
		days = 1;
	}

	bgreg::date  temp;
	bgreg::date  temp2;
	temp = inT1;
	double dcf = DCC::Actual365Fixed::dayCountFraction(temp, inT2);
	bgreg::days d(days);
	temp2 = temp + d;
	double dt = DCC::Actual365Fixed::dayCountFraction(temp, temp2);

	while (dcf>0){
		if (inUsePivot){
			mIRSpivot.swapRateSensitivity(temp, false, diffi, diffi);  //must be at x=0 y=0
		}
		else{
			mIRS.swapRateSensitivity(temp, false, diffi, diffi);
		}

		mIRS.swapRateSensitivity(temp, true, diffk, diffkxx);//must be at x=0 y=0

		if (flag == 0){
			diffsquared = 2.0*diffk*lambda*SWRate*(lambda*diffk*diffi);
		}
		else if (flag == 1){
			diffsquared = 2.0*diffk*lambda*lambda*SWRate*(diffkxx*SWRate + gamma*diffk*diffi);
		}
		else{
			diffsquared = diffk;
		}
		sum = sum + diffsquared*dt;
		temp = temp2;
		temp2 = temp + d;
		dcf = DCC::Actual365Fixed::dayCountFraction(temp2, inT2);
		dt = DCC::Actual365Fixed::dayCountFraction(temp, temp2);
	}
	dt = DCC::Actual365Fixed::dayCountFraction(temp, inT2);
	if (inUsePivot){
		mIRSpivot.swapRateSensitivity(temp, false, diffi, diffi);  //must be at x=0 y=0
	}
	else{
		mIRS.swapRateSensitivity(temp, false, diffi, diffi);
	}

	mIRS.swapRateSensitivity(temp, true, diffk, diffkxx);//must be at x=0 y=0
	if (flag == 0){
		diffsquared = 2.0*diffk*lambda*SWRate*(lambda*diffk*diffi);
	}
	else if (flag == 1){
		diffsquared = 2.0*diffk*lambda*lambda*SWRate*(diffkxx*SWRate + gamma*diffk*diffi);
	}
	else{
		diffsquared = diffk;
	}
	sum = sum + diffsquared*dt;
	return sum;
}


double DDCalibration::getAverageY(double inT){
	double sum = 0.0;

	for (int i = 0; i < mNumberOfVols; i++){
		if (mDoubleExercisesDates[i] <= inT){
			if (i == 0){
				sum += (mLambdas[i] * mLambdas[i] * mSwapRate[i] * mSwapRate[i])
					*(exp(2 * mKappa*(mDoubleExercisesDates[i])) - exp(2 * mKappa*(0))) / (2 * mKappa);
			}
			else{
				sum += (mLambdas[i] * mLambdas[i] * mSwapRate[i] * mSwapRate[i])
					*(exp(2 * mKappa*(mDoubleExercisesDates[i])) - exp(2 * mKappa*(mDoubleExercisesDates[i - 1]))) / (2 * mKappa);
			}
		}
		else{
			if (i == 0){
				sum += (mLambdas[i] * mLambdas[i] * mSwapRate[i] * mSwapRate[i])
					*(exp(2 * mKappa*(inT)) - exp(2 * mKappa*(0))) / (2 * mKappa);
			}
			else{
				sum += (mLambdas[i] * mLambdas[i] * mSwapRate[i] * mSwapRate[i])
					*(exp(2 * mKappa*(inT)) - exp(2 * mKappa*(mDoubleExercisesDates[i - 1]))) / (2 * mKappa);
			}
			break;
		}
	}
	return exp(-2 * mKappa*inT)*sum;
}

void DDCalibration::getAverageX0(double inT){
	mTempMarkov[1] = getAverageY(inT);
	/*boost::function<double(double &)> f = std::bind1st(std::mem_fun(&DDCalibration::swapRateFunction), this);
	mTempMarkov[0] = RootFinder::zbrent(f, -0.01, 0.01, 1e-8);*/
} 
double DDCalibration::swapRateFunction(double inX){
		mTempMarkov[0] = inX;
		mIRS.updatePrice(mTempMarkov);
		return mIRS.getSwapRate() - mtempSwapRate;
}

void DDCalibration::getAverageX(double inT){
	//double lSWRdiff;
	//double lSWRdiffSQR;
	//bgreg::days d((int)((ceil)(365.0*inT)));
	////mIRS.changeStructure(mExercisesDates[0], 0.0);
	////mIRS.getPrice();

	//double lVariance = getApproximatedVariance(inT);// A revoir car quand one change la structure de today at ti , tout change.
	//mIRS.swapRateSensitivity(mTempMarkov[0],mTempMarkov[1],mToday + d, true, lSWRdiff, lSWRdiffSQR);
	//double lDiffInverseSWR = -(lSWRdiffSQR) / (lSWRdiff*lSWRdiff*lSWRdiff);// a revoir
	//mTempMarkov[0] += 0.5*lDiffInverseSWR*lVariance;
}

double DDCalibration::getApproximatedVariance(double inT){

	
	double lIntegral;
	double sum = 0.0;
	bgreg::date lTi;
	bgreg::date lTim;
	for (int j = 0; j < mNumberOfVols; j++){
		lTim = mToday;	
		sum = 0.0;
		if (mDoubleExercisesDates[j] <= inT){
			lTi = mExercisesDates[j];
			lIntegral = mIRS.ClosedFormIntegral(mToday, lTim, lTi);
			sum += mLambdas[j] * mLambdas[j] * pow(mSwapRate[j], 2)*lIntegral;
			lTim = lTi;
		}
		else{
			bgreg::days d((int)((ceil)(365.0*inT)));
			lIntegral = mIRS.ClosedFormIntegral(mToday, lTim, mToday + d);
			sum += mLambdas[j] * mLambdas[j] * pow(mSwapRate[j], 2)*lIntegral;
			break;
		}
	}
	return sum;
}


double DDCalibration::integrator(const std::vector<double> & inX){
	int i =0;
	double lambda;
	double beta;
	double lambdaB = 0;
	double omegatilda = 0;
	double omegaB = 0;
	double bB = 0;
	double NewbB = 0;
	double NewDenominatorB = 0;

	double lSWRdiff, lSWRdiffSQR;
	double localSWRate;
	double localdiffSWRate;

	double dt;
	double lVol, lVoldiff;
	int compteur = mCompteur;
	

	int N = mTimeGrid.getGridSize();
	/*double testB = 0;
	double testBBarre = 0;
	double testLambda = 0;
	double testSum = 0.0;
	double lTestFormerMat =0.0;
	double lTestNewMat;*/
	std::vector<double> lPhij;

	for (int j = 0; j <=compteur; j++){
		if (j != compteur){
			mAdvancedLambda = mLambdas[j];
			mAdvancedBeta = mBetas[j];
		}
		else{
			mAdvancedLambda = inX[0];
			mAdvancedBeta = inX[1];
		}
		/*testSum = 0.0;
		for (int l = 0; l < j; l++){
			lTestNewMat=DCC::Actual365Fixed::dayCountFraction(mToday, mExercisesDates[l]);
			testSum += mLambdas[l] * mLambdas[l] * (lTestNewMat - lTestFormerMat);
			lTestFormerMat = lTestNewMat;
		}
		lTestNewMat = DCC::Actual365Fixed::dayCountFraction(mToday, mExercisesDates[j]);
		testSum = testSum*(lTestNewMat - lTestFormerMat);
		testSum += 0.5*mLambdas[j] * mLambdas[j] * (lTestNewMat - lTestFormerMat)*(lTestNewMat - lTestFormerMat);
		testSum *= mLambdas[j] * mLambdas[j];
		testBBarre += testSum;
		testB += mBetas[j] * testSum;
		testLambda += mLambdas[j] * mLambdas[j] * (lTestNewMat - lTestFormerMat);*/
	/*	mAdvancedLambda = inX[0];
		mAdvancedBeta = inX[1];*/
		do {
			i++;
			mIRS.changeStructure(mExercisesDates[compteur], mTimeGrid[i]);
			mIRS.getPrice();
			mtempSwapRate = mSwapRate[compteur];
			//getAverageX0(mTimeGrid[i]);
			//getAverageX(mTimeGrid[i]);
			
			bgreg::days d((int)(ceil(365.0*mTimeGrid[i])));
			mIRS.swapRateSensitivity(mTempMarkov[0], mTempMarkov[1], mToday + d, true, lSWRdiff, lSWRdiffSQR);
			mIRS.changeStructure(mExercisesDates[j], mTimeGrid[i]);
			mIRS.getPrice();
			mIRS.swapRateSensitivity(mTempMarkov[0], mTempMarkov[1], mToday + d, false, localdiffSWRate, localdiffSWRate);
			mIRS.updatePrice(mTempMarkov);
			localSWRate = mIRS.getSwapRate();
			lVol = mAdvancedLambda*(mAdvancedBeta*localSWRate + (1 - mAdvancedBeta)*mSwapRate[j]);
			lambda = lSWRdiff*lVol / mSwapRate[compteur];
			

			lVoldiff = mAdvancedBeta*mAdvancedLambda*localdiffSWRate;
			//beta = lVoldiff / (lVol*lSWRdiff);
			beta = lVoldiff / (lVol*lSWRdiff) + lSWRdiffSQR / (lSWRdiff*lSWRdiff);
			beta = beta*mSwapRate[compteur];
	

			dt = mTimeGrid[i] - mTimeGrid[i - 1];
			lambdaB += lambda*lambda*dt;

			lPhij.push_back(lambdaB);
			NewDenominatorB += lambda*lambda*dt*lPhij[i-1];
			NewbB += beta*lambda*lambda*dt*lPhij[i-1];
		} while ((mTimeGrid.getIsKeyDate(i) == false)&(i!=N-1));

		mAdvancedMarketLambdas[j] = sqrt(lambdaB / mTimeGrid[i]);
		mAdvancedMarketBetas[j] = NewbB / NewDenominatorB;
	}
	/*testB /= testBBarre;
	testLambda /= lTestNewMat;
	testLambda = sqrt(testLambda);*/
	return pow(((mAdvancedMarketLambdas[compteur] - mMarketLambdas[compteur]) / mMarketLambdas[compteur]), 2)
		+ pow(((mAdvancedMarketBetas[compteur] - mMarketBetas[compteur]) / mMarketBetas[compteur]), 2);

}

void DDCalibration::advancedCalibration(int inCompteur){
	boost::function<double(const std::vector<double> &)> inFunc  = boost::bind(&DDCalibration::integrator, this,_1);
	std::vector<double> lP(2);
	int lIter=0;
	
	mCompteur = inCompteur;
	lP[0] = mLambdas[mCompteur];
	lP[1] = mBetas[mCompteur];
	std::vector<double> oLower(2);
	std::vector<double> oUpper(2);
	double lLambdaEps = 0.01;
	double lBetaEps = 0.1;

	oLower[0] = lP[0] - lLambdaEps;
	oLower[1] = lP[1] - lBetaEps;
	oUpper[0] = lP[0] +lLambdaEps;
	oUpper[1] = lP[1] + lBetaEps;
	double lfret = 0.0;
	//double lfret = integrator(lP);
	
	Optimizer::dfpmin(lP, 1.0e-15, lIter, lfret, inFunc);
	//int N = std::max<std::size_t>(boost::thread::hardware_concurrency(), 1);
	//DifferentialEvolution2 oDE;
	//oDE.setNumberOfThreads(N);
	//int res = oDE.Optimize(inFunc, lP, lfret, oLower, oUpper);
	mLambdas[mCompteur] = lP[0];
	mBetas[mCompteur] = lP[1];
}