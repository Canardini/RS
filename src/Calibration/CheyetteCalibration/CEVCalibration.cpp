#include<Calibration/CheyetteCalibration/CEVCalibration.h>

CEVCalibration::CEVCalibration(const std::vector<double> &inStripVolatilities,
							   const std::vector<bgreg::date> & inExercisesDates,
							   const InterestRateSwap &inIRS,
							   const bgreg::date & inToday,
							   double inAlpha,
							   TimeGrid &inTimeGrid,
							   int inNumberOfSimulations) :
							   mStripVolatilities(inStripVolatilities), mExercisesDates(inExercisesDates), 
							   mIRS(inIRS), mToday(inToday), mAlpha(inAlpha),mTimeGrid(inTimeGrid), mNumberOfSimulations(inNumberOfSimulations)
{
	mNumberOfVols=mStripVolatilities.size();
	mLambdas.resize(mNumberOfVols);
	mSwapRates.resize(mNumberOfVols);
	mMarketPrices.resize(mNumberOfVols);
	mFWPrices.resize(mNumberOfVols);
}

void CEVCalibration::getLambdas(std::vector<double> &outLambdas){
	int lNout = outLambdas.size();
	if (lNout == mNumberOfVols){
		for (int i = 0; i < lNout; i++){
			outLambdas[i] = mLambdas[i];
		}
	}
}

void CEVCalibration::runCalibration(){
	std::vector<double> lSwapRate(mNumberOfVols);
	double lMaturity;
	double lIntegral;
	double sum = 0.0;
	bgreg::date lTi;
	bgreg::date lTim;
	int inCallOrPut;
	double lAnnuity;
	for (int j = 0; j < mNumberOfVols; j++){
		lTim = mToday;
		mIRS.changeStructure(mExercisesDates[j], 0);
		mFWPrices[j] = mIRS.getPrice();
		lSwapRate[j] = mIRS.getSwapRate();
		mSwapRates[j] = lSwapRate[j];
		lMaturity = DCC::Actual365Fixed::dayCountFraction(mToday, mExercisesDates[j]);
		(mIRS.isPayer() ? inCallOrPut = 1 : inCallOrPut = -1);
		lAnnuity = mIRS.getAnnuity();
		mMarketPrices[j] = Black::black(inCallOrPut, mSwapRates[j], mIRS.getFixedRate(), lAnnuity, lMaturity, mStripVolatilities[j]);
		sum = 0.0;
		for (int i = 0; i < j; i++){
			lTi = mExercisesDates[i];
			lIntegral = mIRS.ClosedFormIntegral(mToday, lTim, lTi);
			sum += mLambdas[i] * mLambdas[i] * pow(lSwapRate[i], 2 * mAlpha)*lIntegral;
			lTim = lTi;
		}
		lTi = mExercisesDates[j];
		lIntegral = mIRS.ClosedFormIntegral(mToday, lTim, lTi);
		mLambdas[j] = mStripVolatilities[j] * mStripVolatilities[j] * lSwapRate[j] * lSwapRate[j] * lMaturity - sum;
		if ((lIntegral <= 0) | (mLambdas[j] < 0)){
			if (j>0){
				mLambdas[j] = mLambdas[j - 1];
			}
			else{
				mLambdas[j] = 0.0;
			}
		}
		else{
			mLambdas[j] /= pow(lSwapRate[j], 2 * mAlpha)*lIntegral;
			mLambdas[j] = sqrt(mLambdas[j]);
		}
	}
}

void CEVCalibration::getBounds(double &inEX, double &inEY, double &inVarX, double &inVarY){
	double lXB;
	double lYB = 0.0;
	double lVarX = 0.0;
	double lKappaFactor;
	double lKappa = mIRS.getMeanReversion();
	double lMaturity = DCC::Actual365Fixed::dayCountFraction(mToday, mExercisesDates[mNumberOfVols-1]);
	double lTim = 0;
	double lTi;
	//Expectation of Y
	for (int j = 0; j < mNumberOfVols; j++){
		lTi = DCC::Actual365Fixed::dayCountFraction(mToday, mExercisesDates[j]);
		lKappaFactor = (exp(-2 * lKappa*(lMaturity - lTi)) - exp(-2 * lKappa*(lMaturity - lTim)));
		lKappaFactor /= 2 * lKappa;
		lYB += mLambdas[j] * pow(mSwapRates[j], mAlpha)*lKappaFactor;
		lVarX += mLambdas[j] * mLambdas[j] * pow(mSwapRates[j], 2 * mAlpha)*lKappaFactor;
		lTim = lTi;
	}
	lXB = (1 - exp(-lKappa*lMaturity))*lYB / lMaturity;

}

void CEVCalibration::monteCarloCorrection(){
	/*mRandomVariables.resize(mNumberOfSimulations);*/
	int mAntitheticNumber =(int)( mNumberOfSimulations / 2);
	boost::mt19937 rng; // I don't seed it on purpouse (it's not relevant)
	boost::normal_distribution<> nd(0.0, 1.0);
	boost::variate_generator<boost::mt19937&,	boost::normal_distribution<> > var_nor(rng, nd);
	int lCompteur = 0;
	int lTimeCompteur = 1;
	double xp, x(0);
	double yp, y(0);
	double kappa = 0.1;
	double vol;
	//double voldiff;
	double z1;
	double dt;

	double MCprice = 0;

	double integralR = 0;
	
	double lSwapRatediff = 0;
	double T = 0;
	double f0t;
	std::vector<double> vecf0t(0);
	std::vector<double> vecdt(0);
	std::vector<double> lMarkovStates(2);
	bgreg::days d;
	for (int i = 0; i < mNumberOfSimulations; i++){
		lCompteur = 0;
		lTimeCompteur = 0;
		x = 0.0;
		y = 0.0;
		integralR = 0;
		//mIRS.changeStructure(mTimeGrid.getKeyDate(0), 0.0);
		//double tempPrice=mIRS.getPrice();
		//lSWRate = mIRS.getSwapRate();

		//integralR += (x + 0.03)*dt;

		////integralR += (x + 0.0)*dt;
		//vol = mLambdas[0] * pow(abs(lSWRate), mAlpha);
		//if (abs(lSWRate) < 1e-8){
		//	voldiff = 0.0;
		//}
		//else{
		//	mIRS.swapRateSensitivity(mToday, false, lSwapRatediff, lSwapRatediff);
		//	voldiff = mLambdas[0] * mAlpha*(pow(abs(lSWRate), mAlpha - 1.0))*lSwapRatediff;
		//}
		T = 0;
		while (!mTimeGrid.getIsKeyDate(lTimeCompteur)){
			z1 = var_nor();
			vol = 0.2; 
			if (i == 0){
				dt = mTimeGrid[lTimeCompteur + 1] - mTimeGrid[lTimeCompteur];
				vecdt.push_back(dt);
				f0t = mIRS.getInitSource(mTimeGrid[lTimeCompteur + 1]);
				vecf0t.push_back(f0t);
			}
			
			xp = x + (y - kappa*x)*vecdt[lTimeCompteur] + vol*sqrt(vecdt[lTimeCompteur])*z1;
			//xp = x + (y - kappa*x)*dt + vol*sqrt(dt)*z1 + 0.5*dt*voldiff*vol*(z1*z1 - 1);
			yp = y + (vol*vol - 2 * kappa*y)*vecdt[lTimeCompteur];
			
			integralR += (xp + vecf0t[lTimeCompteur])*vecdt[lTimeCompteur];
			lTimeCompteur++;
			T += vecdt[lTimeCompteur];
			
			
			//integralR += (xp + 0.03)*dt;


		/*	lMarkovStates[0] = xp;
			lMarkovStates[1] = yp;
			lTimeCompteur++;
			dt = mTimeGrid[lTimeCompteur+1] - mTimeGrid[lTimeCompteur];
			mIRS.changeStructure(mTimeGrid.getKeyDate(lCompteur), mTimeGrid[lTimeCompteur]);
			SwapSample = mIRS.updatePrice(lMarkovStates);
			lSWRate = mIRS.getSwapRate();
			bgreg::days d((int)ceil(mTimeGrid[lTimeCompteur] * 365.0));
			vol = mLambdas[0] * pow(abs(lSWRate), mAlpha);
			integralR += (xp + mIRS.getInitSource(mTimeGrid[lTimeCompteur]))*dt;
			if (abs(lSWRate) < 1e-8){
				voldiff = 0.0;
			}
			else{
				mIRS.swapRateSensitivity(mToday+d, false, lSwapRatediff, lSwapRatediff);
				voldiff = mLambdas[0] * mAlpha*(pow(abs(lSWRate), mAlpha - 1.0))*lSwapRatediff;
			}*/
			x = xp;
			y = yp;
		}
		//integralR -= (xp + mIRS.getInitSource(mTimeGrid[lTimeCompteur]))*dt;
		//integralR -= (xp +0.03)*dt;

		//MCprice += exp(-integralR)*std::max(SwapSample, 0.0);
		//integralR += (xp + mIRS.getInitSource(mTimeGrid[lTimeCompteur]))*dt;

		MCprice += exp(-integralR) ;

	}
	MCprice /= mNumberOfSimulations;
	
}

