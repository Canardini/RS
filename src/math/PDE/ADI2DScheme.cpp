#include <Math/PDE/ADI2DScheme.h>




double f(double x, double y){
	return exp(x)*cos(x*y);
}

double fx(double x, double y){
	return exp(x)*(cos(x*y)-y*sin(x*y));
}
double fxx(double x, double y){
	return fx(x, y) + exp(x)*(-cos(x*y)*y*y - y*sin(x*y));
}


class FillingPayoffEngine{
public:
	FillingPayoffEngine(MarkovianProductBridge inMarkovianProduct);

private:
	MarkovianProductBridge mMarkovianProduct;
	int Threadnumber;
};


ADI2DScheme::ADI2DScheme(MarkovianProductBridge inMarkovianProduct,
	PDEModelBridge inPDEmodel,
	boost::shared_ptr<Grid> inGrid,
	TimeGrid &inTimeGrid,
	int inTimeSteps,
	double inSchemeType,
	PDVector inStencilConvections,
	PDVector inStencilDiffusions,
	PDVector inStencilMixedDerivatives,
	bool inUseCorrector,
	double inLambda,
	std::string & inProductType) :
	ADIScheme(inMarkovianProduct, inPDEmodel, inGrid, inTimeGrid, inTimeSteps, inSchemeType, inStencilConvections, inStencilDiffusions, inStencilMixedDerivatives, inUseCorrector, inLambda), mProductType(inProductType)
{
	mPDEDimension = 2;
	mObjectiveFunction = dMatrix(mGrid->getGridSize(0), dVector(mGrid->getGridSize(1)));
	//mObjectiveFunction = boost::numeric::ublas::matrix<double>(mGrid->getGridSize(0),mGrid->getGridSize(1));
	mEntities = dMatrix(mGrid->getGridSize(0), dVector(mGrid->getGridSize(1)));
	mExerciseFlags = std::vector<std::vector<int>>(mGrid->getGridSize(0), std::vector<int>(mGrid->getGridSize(1)));
	//mPrices = boost::numeric::ublas::matrix<double>(mGrid->getGridSize(0), mGrid->getGridSize(1));
	mPrices = dMatrix(mGrid->getGridSize(0), dVector(mGrid->getGridSize(1)));

	//mVol = dMatrix(mGrid->getGridSize(0), dVector(mGrid->getGridSize(1)));
	mPredictors.resize(1);
	if (false){

	//if ((mGrid->getIsUniformGrid(0)) & (mGrid->getIsUniformGrid(1))){
		mUniformdX = mGrid->getGridElement(0, 1) - mGrid->getGridElement(0, 0);
		mUniformdY = mGrid->getGridElement(1, 1) - mGrid->getGridElement(1, 0);
		if (mStencilConvections[0] == ThreePoints){
			std::vector<double> lTempDX(2, mUniformdX);
			mStencil1STX.resize(3);
			stencilFirstOrderCoefficients(lTempDX, true, ThreePoints, mStencil1STX);
		}
		else{
			std::vector<double> lTempDX(4, mUniformdX);
			mStencil1STX.resize(5);
			stencilFirstOrderCoefficients(lTempDX, true, FivePoints, mStencil1STX);
		}
		if (mStencilConvections[1] == ThreePoints){
			std::vector<double> lTempDY(2, mUniformdX);
			mStencil1STY.resize(3);
			stencilFirstOrderCoefficients(lTempDY, true, ThreePoints, mStencil1STY);
		}
		else{
			std::vector<double> lTempDY(4, mUniformdX);
			mStencil1STY.resize(5);
			stencilFirstOrderCoefficients(lTempDY, true, FivePoints, mStencil1STY);
		}
		if (mStencilDiffusions[0] == ThreePoints){
			std::vector<double> lTempDX(2, mUniformdX);
			mStencil2ndX.resize(3);
			stencilSecondOrderCoefficients(lTempDX, true, ThreePoints, mStencil2ndX);
		}
	}
	else{
		generateDX();
		generateFirstStencils();
		generateSecondStencils();
	}
	/*mNumberOfCPUs = std::max<std::size_t>(boost::thread::hardware_concurrency(),1);

	for (int i = 0; i < mNumberOfCPUs; i++){
		mMarkovianProductVector.push_back(mMarkovianProduct);
	}*/

}



void ADI2DScheme::fillOFwithPayoff(){
	std::vector<double> lStates(2);
	//dMatrix mFx = dMatrix(mGrid->getGridSize(0), dVector(mGrid->getGridSize(1)));
	//dMatrix mFxx = dMatrix(mGrid->getGridSize(0), dVector(mGrid->getGridSize(1)));
	//dMatrix mFDx = dMatrix(mGrid->getGridSize(0), dVector(mGrid->getGridSize(1)));
	//dMatrix mFDxx = dMatrix(mGrid->getGridSize(0), dVector(mGrid->getGridSize(1)));

	for (int i = 0; i < mGrid->getGridSize(0); i++){
		lStates[0] = mGrid->getGridElement(0, i);
		for (int j = 0; j < mGrid->getGridSize(1); j++){
			lStates[1] = mGrid->getGridElement(1, j);
			if ((i == 0) & (j == 0)){
				mMarkovianProduct.updateMarkovianFactors(lStates);
				// mObjectiveFunction[i][j] = mMarkovianProduct.getPrice();
				//mObjectiveFunction[i][j] = std::max(mMarkovianProduct.getPrice(), 0.0);
				//mObjectiveFunction(i,j) = mMarkovianProduct.getPrice();
				//mObjectiveFunction(i, j) = 1.0;

				//mObjectiveFunction(i,j) = std::max(mMarkovianProduct.getPrice(), 0.0);
				mObjectiveFunction[i][j] = std::max(mMarkovianProduct.getPrice(), 0.0);

				/*mObjectiveFunction[i][j] = f(lStates[0], lStates[1]);
				mFx[i][j] = fx(lStates[0], lStates[1]);
				mFxx[i][j] = fxx(lStates[0], lStates[1]);*/
				/*mFDx[i][j] = newDifferentiator2D(i,j,ThreePoints,true,false);
				mFDxx[i][j] = newDifferentiator2D(i, j,ThreePoints, true, false);*/

				//mObjectiveFunction[i][j] =1.0;

			}
			else{
				//mObjectiveFunction[i][j] = mMarkovianProduct.updatePrice(lStates);
				mObjectiveFunction[i][j] = std::max(mMarkovianProduct.updatePrice(lStates),0.0);
				//mObjectiveFunction(i, j) = std::max(mMarkovianProduct.updatePrice(lStates), 0.0);
				//mObjectiveFunction(i, j) = mMarkovianProduct.updatePrice(lStates);
				//mObjectiveFunction(i, j) = 1.0;

				/*mObjectiveFunction[i][j] = f(lStates[0], lStates[1]);
				mFx[i][j] = fx(lStates[0], lStates[1]);
				mFxx[i][j] = fxx(lStates[0], lStates[1]);*/

				//mObjectiveFunction[i][j] = 1.0;
			}
			//if (mObjectiveFunction(i, j) == 0.0){
			if (mObjectiveFunction[i][j] == 0.0){

				mExerciseFlags[i][j] = 0;
			}
			else{
				mExerciseFlags[i][j] = 1;
			}
		}
	}



	int lY = mGrid->locate(1, mAverageDistribution);
	double lIndexY;
	double lProbaTotal(0.0);
	double lProbaNm, lProbaN;
	double dm, d;
	int lSizeX = mGrid->getGridSize(0);
	(lY == mGrid->getGridSize(1) ? lIndexY = lY : lIndexY = lY + 1);
	if (lY == -1){ lY++; };
	for (int i = 0; i < lSizeX; i++){
		if ((mExerciseFlags[i][lIndexY] == 1.0) | (mExerciseFlags[i][lIndexY - 1] == 1.0)){
			if (i == 0){
				d = mGrid->getGridElement(0, i) - mAverageDistribution;
				d /= sqrt(4*mVarianceDistribution);
				lProbaN = Black::gauss_prob(d);
				lProbaTotal += lProbaN;
			}
			else if (i == lSizeX - 1){
				d = mGrid->getGridElement(0, i) - mAverageDistribution;
				d /= sqrt(4*mVarianceDistribution);
				lProbaN = 1 - Black::gauss_prob(d);
				lProbaTotal += lProbaN;
			}

			else{
				dm = mGrid->getGridElement(0, i - 1) - mAverageDistribution;
				dm /= sqrt(4*mVarianceDistribution);
				d = mGrid->getGridElement(0, i) - mAverageDistribution;
				d /= sqrt(4*mVarianceDistribution);
				lProbaNm = Black::gauss_prob(dm);
				lProbaN = Black::gauss_prob(d);
				lProbaTotal += lProbaN - lProbaNm;
			}
		}
	}
	mProbaComponents.push_back(lProbaTotal);
}

class LocalFunction{
	ADI2DScheme* mObj;
	int mIndex;

public:
	LocalFunction(ADI2DScheme* inObj, int inIndex) : mObj(inObj), mIndex(inIndex) {}
	void operator()(){
		mObj->fillOFwithPayoff(mIndex);
	}
};

class LocalFunction2{
	ADI2DScheme* mObj;
	int mIndex;

public:
	LocalFunction2(ADI2DScheme* inObj, int inIndex) : mObj(inObj), mIndex(inIndex) {}
	void operator()(){
		mObj->fillEntities(mIndex);
	}
};

void ADI2DScheme::multiThreadPayoff(){
	boost::thread_group lThreads;
	//boost::barrier lBarrier(mNumberOfCPUs);
	for (int i = 0; i < mNumberOfCPUs; ++i)
	{
		lThreads.create_thread(LocalFunction(this, i));
	}
	lThreads.join_all();
}

void ADI2DScheme::multiThreadEntity(){
	boost::thread_group lThreads;
	//boost::barrier lBarrier(mNumberOfCPUs);
	for (int i = 0; i < mNumberOfCPUs; ++i)
	{
		lThreads.create_thread(LocalFunction2(this, i));
	}
	lThreads.join_all();
}

void ADI2DScheme::fillOFwithPayoff(int inThread){
	std::vector<double> lStates(2);
	for (int i = 0; i < mGrid->getGridSize(0); i++){
		if (i%mNumberOfCPUs == inThread){
			lStates[0] = mGrid->getGridElement(0, i);
			for (int j = 0; j < mGrid->getGridSize(1); j++){
				lStates[1] = mGrid->getGridElement(1, j);
				if ((i == inThread) & (j == 0)){
					mMarkovianProductVector[inThread].updateMarkovianFactors(lStates);
					//mObjectiveFunction(i, j) = std::max(mMarkovianProductVector[inThread].getPrice(), 0.0);
					mObjectiveFunction[i][j] = std::max(mMarkovianProductVector[inThread].getPrice(), 0.0);

				}
				else{
					//mObjectiveFunction(i, j) = std::max(mMarkovianProductVector[inThread].updatePrice(lStates), 0.0);
					mObjectiveFunction[i][j] = std::max(mMarkovianProductVector[inThread].updatePrice(lStates), 0.0);

				}
			}
		}
	}


}

void ADI2DScheme::fillEntities(){
	std::vector<double> lStates(2);
	for (int i = 0; i < mGrid->getGridSize(0); i++){
		lStates[0] = mGrid->getGridElement(0, i);
		for (int j = 0; j < mGrid->getGridSize(1); j++){
			lStates[1] = mGrid->getGridElement(1, j);
			if ((i == 0) & (j == 0)){
				mMarkovianProduct.updateMarkovianFactors(lStates);
				mPrices[i][j] = mMarkovianProduct.getPrice();
				//mPrices(i, j) = mMarkovianProduct.getPrice();

				mEntities[i][j] = mMarkovianProduct.getKeyEntity();
			}
			else{
				mPrices[i][j] = mMarkovianProduct.updatePrice(lStates);
				//mPrices(i, j) = mMarkovianProduct.updatePrice(lStates);

				mEntities[i][j] = mMarkovianProduct.getKeyEntity();
			}
		}
	}
}



void ADI2DScheme::fillExercisePayOff(){
	std::vector<double> lStates(2);
	for (int i = 0; i < mGrid->getGridSize(0); i++){
		for (int j = 0; j < mGrid->getGridSize(1); j++){
			//mObjectiveFunction(i, j) = std::max(mObjectiveFunction(i, j), mPrices(i, j));
			//if (mObjectiveFunction(i, j) == mPrices(i, j)){
			mObjectiveFunction[i][j] = std::max(mObjectiveFunction[i][j], mPrices[i][j]);
			if (mObjectiveFunction[i][j] == mPrices[i][j]){
				mExerciseFlags[i][j] = 1;
			}
			else{
				mExerciseFlags[i][j] = 0;
			}
		}
	}

	int lY = mGrid->locate(1, mAverageDistribution);
	double lIndexY;
	double lProbaTotal(0.0);
	double lProbaNm, lProbaN;
	double dm, d;
	int lSizeX = mGrid->getGridSize(0);
	(lY == mGrid->getGridSize(1) ? lIndexY = lY : lIndexY = lY + 1);
	if (lY == -1){ lY++; };
	for (int i = 0; i < lSizeX; i++){
		if ((mExerciseFlags[i][lIndexY] == 1.0) | (mExerciseFlags[i][lIndexY - 1] == 1.0)){
			if (i == 0){
				d = mGrid->getGridElement(0, i) - mAverageDistribution;
				d /= sqrt(4*mVarianceDistribution);
				lProbaN = Black::gauss_prob(d);
				lProbaTotal += lProbaN;
			}
			else if (i == lSizeX - 1){
				d = mGrid->getGridElement(0, i) - mAverageDistribution;
				d /= sqrt(4*mVarianceDistribution);
				lProbaN = 1-Black::gauss_prob(d);
				lProbaTotal += lProbaN;
			}
			
			else{
				dm = mGrid->getGridElement(0, i - 1) - mAverageDistribution;
				dm /= sqrt(4*mVarianceDistribution);
				d = mGrid->getGridElement(0, i) - mAverageDistribution;
				d /= sqrt(4*mVarianceDistribution);
				lProbaNm = Black::gauss_prob(dm);
				lProbaN = Black::gauss_prob(d);
				lProbaTotal += lProbaN - lProbaNm;
			}
		}
	}
	mProbaComponents.push_back(lProbaTotal);
}


void ADI2DScheme::fillEntities(int inThread){
	std::vector<double> lStates(2);
	for (int i = 0; i < mGrid->getGridSize(0); i++){
		if (i%mNumberOfCPUs == inThread){
			lStates[0] = mGrid->getGridElement(0, i);
			for (int j = 0; j < mGrid->getGridSize(1); j++){
				lStates[1] = mGrid->getGridElement(1, j);
				if ((i == inThread) & (j == 0)){
					mMarkovianProductVector[inThread].updateMarkovianFactors(lStates);
					mEntities[i][j] = mMarkovianProductVector[inThread].getPrice();
					mEntities[i][j] = mMarkovianProductVector[inThread].getKeyEntity();
				}
				else{
					mEntities[i][j] = mMarkovianProductVector[inThread].updatePrice(lStates);
					mEntities[i][j] = mMarkovianProductVector[inThread].getKeyEntity();
				}
			}
		}
	}
}

void ADI2DScheme::generateVol(){
	std::vector<double> lStates(2);
	for (int i = 0; i < mGrid->getGridSize(0); i++){
		lStates[0] = mGrid->getGridElement(0, i);
		for (int j = 0; j < mGrid->getGridSize(1); j++){
			lStates[1] = mGrid->getGridElement(1, j);
			if ((i == 0) & (j == 0)){
				mMarkovianProduct.updateMarkovianFactors(lStates);
				mEntities[i][j] = mMarkovianProduct.getPrice();
				mEntities[i][j] = mMarkovianProduct.getKeyEntity();
			}
			else{
				mEntities[i][j] = mMarkovianProduct.updatePrice(lStates);
				mEntities[i][j] = mMarkovianProduct.getKeyEntity();
			}
		}
	}
}


//void ADI2DScheme::solve(int timeIndex){
//
//	int compteur = mTimeGrid.getNumOfKeyDates() - 1;
//	mMarkovianProduct.changeStructure(mTimeGrid.getKeyDate(compteur), mTimeGrid[timeIndex]);
//	fillOFwithPayoff();
//	double dt;
//	for (int i = timeIndex - 1; i >= 0; i--)
//	{
//		dt = mTimeGrid[i + 1] - mTimeGrid[i];
//		if (mTimeGrid.getIsKeyDate(i)){
//			compteur--;
//			mMarkovianProduct.changeStructure(mTimeGrid.getKeyDate(compteur), mTimeGrid[i]);
//		}
//		mMarkovianProduct.setValueDate(mTimeGrid[i]);
//		fillEntities();
//		mPDEmodel.setIndex(compteur);
//		mPDEmodel.setInitSource(mMarkovianProduct.getInitSource(mTimeGrid[i]));
//		solveState(0, 1, dt);
//		solveState(1, 0, dt);
//		applyBoundaryConditions(Numerical);
//
//		if (mTimeGrid.getIsKeyDate(i)){
//			//fillOFwithPayoff();
//		}
//	}
//}

void ADI2DScheme::solve(int compteur){

	int lcompteur = compteur;
	int ltimeIndex = mTimeGrid.getIndexKeyDates(compteur);
	mMarkovianProduct.changeStructure(mTimeGrid.getKeyDate(compteur), mTimeGrid[ltimeIndex]);
	mVarianceDistribution = mMarkovianProduct.getDrifAdjust(mTimeGrid[ltimeIndex]);
	mAverageDistribution = mMarkovianProduct.getGtT(0.0, mTimeGrid[ltimeIndex])*mVarianceDistribution;
	fillOFwithPayoff();
	double dt;
	
	for (int i = ltimeIndex - 1; i >= 0; i--)
	{
		dt = mTimeGrid[i + 1] - mTimeGrid[i];
		if (mTimeGrid.getIsKeyDate(i)){
			compteur--;
			mMarkovianProduct.changeStructure(mTimeGrid.getKeyDate(compteur), mTimeGrid[i]);
		}
		mMarkovianProduct.setValueDate(mTimeGrid[i]);
		fillEntities();
		mPDEmodel.setIndex(compteur);
		mPDEmodel.setInitSource(mMarkovianProduct.getInitSource(mTimeGrid[i]));
		solveState(0, 1, dt);
		solveState(1, 0, dt);
		applyBoundaryConditions(Numerical);

		if (mProductType != "European"){
			if (mTimeGrid.getIsKeyDate(i)){
				mVarianceDistribution = mMarkovianProduct.getDrifAdjust(mTimeGrid[i]);
				mAverageDistribution = mMarkovianProduct.getGtT(0.0, mTimeGrid[i])*mVarianceDistribution;
				fillExercisePayOff();
			}
		}
	}
}

void ADI2DScheme::solveThr(int timeIndex){

	int compteur = mTimeGrid.getNumOfKeyDates() - 1;
	for (int i = 0; i < mNumberOfCPUs; i++){
		mMarkovianProductVector[i].changeStructure(mTimeGrid.getKeyDate(compteur), mTimeGrid[timeIndex]);
	}
	mMarkovianProduct.changeStructure(mTimeGrid.getKeyDate(compteur), mTimeGrid[timeIndex]);
	std::clock_t start;
	start = std::clock();
	multiThreadPayoff();
	mDurations.push_back((std::clock() - start) / (double)CLOCKS_PER_SEC);
	double dt;
	for (int i = timeIndex - 1; i >= 0; i--)
	{
		dt = mTimeGrid[i + 1] - mTimeGrid[i];
		if (mTimeGrid.getIsKeyDate(i)){
			compteur--;
			mMarkovianProduct.changeStructure(mTimeGrid.getKeyDate(compteur), mTimeGrid[i]);
			for (int k = 0; k < mNumberOfCPUs; k++){
				mMarkovianProductVector[k].changeStructure(mTimeGrid.getKeyDate(compteur), mTimeGrid[k]);
			}
		}
		mMarkovianProduct.setValueDate(mTimeGrid[i]);
		for (int k = 0; k < mNumberOfCPUs; k++){
			mMarkovianProductVector[k].setValueDate(mTimeGrid[i]);
		}
		multiThreadEntity();
		mPDEmodel.setIndex(compteur);
		mPDEmodel.setInitSource(mMarkovianProduct.getInitSource(mTimeGrid[i]));
		solveState(0, 1, dt);
		solveState(1, 0, dt);
		applyBoundaryConditions(Numerical);

		if (mTimeGrid.getIsKeyDate(i)){
			//fillOFwithPayoff();
		}
	}
}

//void ADI2DScheme::solveThr(int timeIndex){
//
//	
//
//	/*typedef boost::shared_ptr<ThreadOperator> CalcP;
//	typedef boost::shared_ptr<boost::thread> ThreadP;
//	std::vector<ThreadP> lThreads;
//	std::vector<CalcP> lCalculatorVec;
//	for (int i = timeIndex - 1; i >= 0; i--)
//	{
//		for (int j = 0; j < mNumberOfCPUs; j++){
//			CalcP lCalp(new ThreadOperator(this, 0, 1, mTimeGrid[i + 1] - mTimeGrid[i],j));
//			lThreads.push_back(ThreadP(new boost::thread(boost::ref(*lCalp))));
//		}
//		for (int j = 0; j < mNumberOfCPUs; j++){
//			lThreads[i]->join();
//		}
//		applyBoundaryConditions(Numerical);
//	}*/
//
//
//}


void ADI2DScheme::generateCorrectorBandMatrix(dMatrix & outCorrector, int inStep){

}
void ADI2DScheme::applyBoundaryConditions(BoundaryConditions type){
	int lStarterX;
	int lStarterY;
	bool bFlagX = (mStencilConvections[0] == FivePoints) | (mStencilDiffusions[0] == FivePoints);
	bool bFlagY = (mStencilConvections[1] == FivePoints) | (mStencilDiffusions[1] == FivePoints);
	(bFlagX ? lStarterX = 2 : lStarterX = 1);
	(bFlagY ? lStarterY = 2 : lStarterY = 1);
	int nX = mObjectiveFunction.size();
	//int nX = mObjectiveFunction.size1();
	int nY = mObjectiveFunction[0].size();
	//int nY = mObjectiveFunction.size2();
	for (int j = lStarterY; j < nY-lStarterY; j++){
		if (bFlagX){
			mObjectiveFunction[1][j] = 2 * mObjectiveFunction[2][j] - mObjectiveFunction[3][j];
			//mObjectiveFunction(1,j) = 2 * mObjectiveFunction(2,j) - mObjectiveFunction(3,j);

			mObjectiveFunction[0][j] = 2 * mObjectiveFunction[1][j] - mObjectiveFunction[2][j];	
			//mObjectiveFunction(0,j) = 2 * mObjectiveFunction(1,j) - mObjectiveFunction(2,j);

			mObjectiveFunction[nX-2][j] = 2 * mObjectiveFunction[nX-3][j] - mObjectiveFunction[nX-4][j];
			//mObjectiveFunction(nX - 2,j) = 2 * mObjectiveFunction(nX - 3,j) - mObjectiveFunction(nX - 4,j);

			mObjectiveFunction[nX-1][j] = 2 * mObjectiveFunction[nX-2][j] - mObjectiveFunction[nX-3][j];
			//mObjectiveFunction(nX - 1,j) = 2 * mObjectiveFunction(nX - 2,j) - mObjectiveFunction(nX - 3,j);
		}
		else{
			mObjectiveFunction[0][j] = 2 * mObjectiveFunction[1][j] - mObjectiveFunction[2][j];
			//mObjectiveFunction(0,j) = 2 * mObjectiveFunction(1,j) - mObjectiveFunction(2,j);

			mObjectiveFunction[nX-1][j] = 2 * mObjectiveFunction[nX-2][j] - mObjectiveFunction[nX-3][j];
			//mObjectiveFunction(nX - 1,j) = 2 * mObjectiveFunction(nX - 2,j) - mObjectiveFunction(nX - 3,j);

		}
	}

	for (int i = 0; i < nX; i++){
		if (bFlagY){
			mObjectiveFunction[i][1] = 2 * mObjectiveFunction[i][2] - mObjectiveFunction[i][3];
			//mObjectiveFunction(i,1) = 2 * mObjectiveFunction(i,2) - mObjectiveFunction(i,3);

			mObjectiveFunction[i][0] = 2 * mObjectiveFunction[i][1] - mObjectiveFunction[i][2];
			//mObjectiveFunction(i,0) = 2 * mObjectiveFunction(i,1) - mObjectiveFunction(i,2);

			mObjectiveFunction[i][nY - 2] = 2 * mObjectiveFunction[i][nY - 3] - mObjectiveFunction[i][nY - 4];
			//mObjectiveFunction(i,nY - 2) = 2 * mObjectiveFunction(i,nY - 3) - mObjectiveFunction(i,nY - 4);

			mObjectiveFunction[i][nY - 1] = 2 * mObjectiveFunction[i][nY - 2] - mObjectiveFunction[i][nY - 3];
			//mObjectiveFunction(i,nY - 1) = 2 * mObjectiveFunction(i,nY - 2) - mObjectiveFunction(i,nY - 3);

		}
		else{
			mObjectiveFunction[i][0] = 2 * mObjectiveFunction[i][1] - mObjectiveFunction[i][2];
			//mObjectiveFunction(i,0) = 2 * mObjectiveFunction(i,1) - mObjectiveFunction(i,2);

			mObjectiveFunction[i][nY - 1] = 2 * mObjectiveFunction[i][nY - 2] - mObjectiveFunction[i][nY - 3];
			//mObjectiveFunction(i,nY - 1) = 2 * mObjectiveFunction(i,nY - 2) - mObjectiveFunction(i,nY - 3);

		}
	}

}

void ADI2DScheme::generateLowerBoundCondition(int index, PointDiscretization inPD, BoundaryConditions inBC, std::vector<double> outVec){
}
void ADI2DScheme::generateUpperBoundCondition(int index, PointDiscretization inPD, BoundaryConditions inBC, std::vector<double> outVec){
}


void ADI2DScheme::solveState(int index,int pivot, double dt){
	dMatrix A;
	dMatrix al;
	dVector dPivot(2);
	dVector b;
	dVector u;
	int Nindex, Npivot;
	double d;
	std::vector<int> indx;
	Nindex = mGrid->getGridSize(index);
	Npivot = mGrid->getGridSize(pivot);
	int lIndexPivot;
	int Npoints;
	int lStarter;
	bool bFlag = (mStencilConvections[pivot] == FivePoints) | (mStencilDiffusions[pivot] == FivePoints);
	(bFlag ? Npoints = 5 : Npoints = 3);
	(bFlag ? lStarter = 2 : lStarter = 1);

	if ((mStencilConvections[index] == FivePoints) | (mStencilDiffusions[index] == FivePoints)){
			A.resize(Nindex - 4, dVector(5));
			al.resize(Nindex - 4, dVector(2));
			b.resize(Nindex - 4);
			indx.resize(Nindex - 4);
	}
	else{
			A.resize(3, dVector(Nindex - 2));
			mPredictors[index] = dMatrix(Nindex-2, dVector(Npivot - (Npoints - 1))); // can be done before in the code
			b.resize(Nindex - 2);
			u.resize(Nindex - 2);
	}
		for (int i = 0; i < Npivot-(Npoints-1); i++){
			lIndexPivot = i + lStarter;
			mIndexPivot = lIndexPivot;
			dPivot[pivot] = mGrid->getGridElement(pivot, lIndexPivot);
			generatePredictorBandMatrix(A, dt, index, dPivot);;//generate the band and the rhs in the same time.
			generateRHSPredictor(b, dt, index, lIndexPivot);
			if ((mStencilConvections[index] == FivePoints)){

				BandSolver::bandec(A, Nindex - 4, 2, 2, al, indx, d);
				BandSolver::dbanbks(A, Nindex - 4, 2, 2, al, indx, b);
				for (int j = 2; j < Nindex - 2; j++){
					mObjectiveFunction[lIndexPivot][j] = b[j-2]; 
					//mObjectiveFunction(lIndexPivot,j) = b[j - 2];

					// obviously we assume the 5 points is for the second factor
				}
			}
			else{
				BandSolver::dtridag(A[0], A[1], A[2], b, u);
				Matrix::fillMatrixWithVector(mPredictors[index], u, i, u.size()); // to optimize, shitty way
				// obviously we assume the 5 points is for the second factor
			}

		}	
}

void ADI2DScheme::solveStateThr(int index, int pivot, double dt,int inNumberofThreads){
 
	
	
	int Nindex, Npivot;
	double d;
	std::vector<int> indx;
	Nindex = mGrid->getGridSize(index);
	Npivot = mGrid->getGridSize(pivot);
	int lIndexPivot;
	int Npoints;
	int lStarter;
	bool bFlag = (mStencilConvections[pivot] == FivePoints) | (mStencilDiffusions[pivot] == FivePoints);
	(bFlag ? Npoints = 5 : Npoints = 3);
	(bFlag ? lStarter = 2 : lStarter = 1);

	
	for (int i = 0; i < Npivot - (Npoints - 1); i++){
		if (i%mNumberOfCPUs == inNumberofThreads){
			dMatrix A;
			dVector b;
			dVector u;
			dMatrix al;
			dVector dPivot(2);
			if ((mStencilConvections[index] == FivePoints) | (mStencilDiffusions[index] == FivePoints)){
				A.resize(Nindex - 4, dVector(5));
				al.resize(Nindex - 4, dVector(2));
				b.resize(Nindex - 4);
				indx.resize(Nindex - 4);
			}
			else{
				A.resize(3, dVector(Nindex - 2));
				mPredictors[index] = dMatrix(Nindex - 2, dVector(Npivot - (Npoints - 1))); // can be done before in the code
				b.resize(Nindex - 2);
				u.resize(Nindex - 2);
			}
			lIndexPivot = i + lStarter;
			mIndexPivot = lIndexPivot;
			dPivot[pivot] = mGrid->getGridElement(pivot, lIndexPivot);
			generatePredictorBandMatrix(A, dt, index, dPivot);;//generate the band and the rhs in the same time.
			generateRHSPredictor(b, dt, index, lIndexPivot);
			if ((mStencilConvections[index] == FivePoints)){

				BandSolver::bandec(A, Nindex - 4, 2, 2, al, indx, d);
				BandSolver::dbanbks(A, Nindex - 4, 2, 2, al, indx, b);
				for (int j = 2; j < Nindex - 2; j++){
					mObjectiveFunction[lIndexPivot][j] = b[j-2]; 
					//mObjectiveFunction(lIndexPivot, j) = b[j - 2];

					// obviously we assume the 5 points is for the second factor
				}
			}
			else{
				BandSolver::dtridag(A[0], A[1], A[2], b, u);
				Matrix::fillMatrixWithVector(mPredictors[index], u, i, u.size()); // to optimize, shitty way
				// obviously we assume the 5 points is for the second factor
			}
		}
	}
}

//void ADI2DScheme::solveStatethr(int index, int pivot, double dt, int inThread){
//	dMatrix A;
//	dMatrix al;
//	dVector dPivot(2);
//	dVector b;
//	dVector u;
//	int Nindex, Npivot;
//	double d;
//	std::vector<int> indx;
//	Nindex = mGrid->getGridSize(index);
//	Npivot = mGrid->getGridSize(pivot);
//	int lIndexPivot;
//	int Npoints;
//	int lStarter;
//	bool bFlag = (mStencilConvections[pivot] == FivePoints) | (mStencilDiffusions[pivot] == FivePoints);
//	(bFlag ? Npoints = 5 : Npoints = 3);
//	(bFlag ? lStarter = 2 : lStarter = 1);
//
//	if ((mStencilConvections[index] == FivePoints) | (mStencilDiffusions[index] == FivePoints)){
//		A.resize(Nindex - 4, dVector(5));
//		al.resize(Nindex - 4, dVector(2));
//		b.resize(Nindex - 4);
//		indx.resize(Nindex - 4);
//	}
//	else{
//		A.resize(3, dVector(Nindex - 2));
//		mPredictors[index] = dMatrix(Nindex - 2, dVector(Npivot - (Npoints - 1))); // can be done before in the code
//		b.resize(Nindex - 2);
//		u.resize(Nindex - 2);
//	}
//	for (int i = 0; i < Npivot - (Npoints - 1); i++){
//		if (i%mNumberOfCPUs == inThread){
//			lIndexPivot = i + lStarter;
//			dPivot[pivot] = mGrid->getGridElement(pivot, lIndexPivot);
//			generatePredictorBandMatrix(A, dt, index, dPivot);
//			generateRHSPredictor(b, dt, index, lIndexPivot);
//			if ((mStencilConvections[index] == FivePoints)){
//
//				BandSolver::bandec(A, Nindex - 4, 2, 2, al, indx, d);
//				BandSolver::dbanbks(A, Nindex - 4, 2, 2, al, indx, b);
//				for (int j = 2; j < Nindex - 2; j++){
//					mObjectiveFunction[lIndexPivot][j] = b[j - 2];
//					// obviously we assume the 5 points is for the second factor
//				}
//			}
//			else{
//				BandSolver::dtridag(A[0], A[1], A[2], b, u);
//				Matrix::fillMatrixWithVector(mPredictors[index], u, i, u.size()); // to optimize, shitty way
//				// obviously we assume the 5 points is for the second factor
//			}
//		}
//	}
//}




//double ADI2DScheme::differentiator2D(int x, int y,bool isUniform,
//	PointDiscretization type, bool IsFirstOrSecondState, bool IsFirstOrSecondDerivative)
//{
//	double dDiff=0.0;
//	dVector dx;
//	int Npoints;
//	int lStarter;
//	(type == ThreePoints ? Npoints = 3 : Npoints = 5);
//	(type == ThreePoints ? lStarter = 1 : lStarter = 2);
//
//	dx.resize(Npoints - 1);
//	for (int i = 0; i < Npoints-1; i++){
//		if (IsFirstOrSecondState){
//			dx[i] = mGrid->getGridElement(0, x - lStarter+ i + 1) - mGrid->getGridElement(0, x - lStarter + i);
//		}
//		else{
//			dx[i] = mGrid->getGridElement(1, y -lStarter + i + 1) - mGrid->getGridElement(1,y -lStarter + i);
//		}
//	}
//	if (type == FivePoints){
//		dx[0] += dx[1];
//		dx[3] += dx[2];
//	}
//	dVector dStencil(Npoints);
//	if (IsFirstOrSecondDerivative){
//		stencilFirstOrderCoefficients(dx, isUniform, type, dStencil);
//	}
//	else{
//		stencilSecondOrderCoefficients(dx, isUniform, type, dStencil);
//	}
//
//	for (int i = 0; i < Npoints; i++){
//		if (IsFirstOrSecondState){
//			dDiff += dStencil[i] * mObjectiveFunction[x - lStarter + i][y];
//		}
//		else{
//			dDiff += dStencil[i] * mObjectiveFunction[x][y - lStarter + i];
//		}
//	}
//	return dDiff;
//}


double ADI2DScheme::newDifferentiator2D(int x, int y,
	PointDiscretization type, bool IsFirstOrSecondState, bool IsFirstOrSecondDerivative)
{
	double dDiff = 0.0;
	int Npoints;
	int lStarter;
	(type == ThreePoints ? Npoints = 3 : Npoints = 5);
	(type == ThreePoints ? lStarter = 1 : lStarter = 2);

	if (IsFirstOrSecondDerivative){

		for (int i = 0; i < Npoints; i++){
			if (IsFirstOrSecondState){
				dDiff += mFirstOrderStencils[0][x - lStarter][i] * mObjectiveFunction[x - lStarter + i][y];
				//dDiff += mFirstOrderStencils[0][x - lStarter][i] * mObjectiveFunction(x - lStarter + i,y);

			}
			else{
				dDiff += mFirstOrderStencils[1][y - lStarter][i] * mObjectiveFunction[x][y - lStarter + i];
				//dDiff += mFirstOrderStencils[1][y - lStarter][i] * mObjectiveFunction(x,y - lStarter + i);

			}
		}
	}
	else{

		for (int i = 0; i < Npoints; i++){
			if (IsFirstOrSecondState){
				dDiff += mSecondOrderStencils[x - lStarter][i] * mObjectiveFunction[x - lStarter + i][y];
				//dDiff += mSecondOrderStencils[x - lStarter][i] * mObjectiveFunction(x - lStarter + i,y);

			}
			else{
				dDiff += mSecondOrderStencils[x - lStarter][i] * mObjectiveFunction[x][y - lStarter + i];
				//dDiff += mSecondOrderStencils[x - lStarter][i] * mObjectiveFunction(x,y - lStarter + i);

			}
		}
	}
	return dDiff;
}

double ADI2DScheme::uniformDifferentiator2D(int x, int y,
	PointDiscretization type, bool IsFirstOrSecondState, bool IsFirstOrSecondDerivative)
{
	double dDiff = 0.0;
	int Npoints;
	int lStarter;
	(type == ThreePoints ? Npoints = 3 : Npoints = 5);
	(type == ThreePoints ? lStarter = 1 : lStarter = 2);

	if (IsFirstOrSecondDerivative){

		for (int i = 0; i < Npoints; i++){
			if (IsFirstOrSecondState){
				dDiff += mFirstOrderStencils[0][x - lStarter][i] * mObjectiveFunction[x - lStarter + i][y];
				//dDiff += mStencil1STX[i] * mObjectiveFunction(x - lStarter + i, y);

			}
			else{
				dDiff += mFirstOrderStencils[1][y - lStarter][i] * mObjectiveFunction[x][y - lStarter + i];
				//dDiff += mStencil1STY[i] * mObjectiveFunction(x, y - lStarter + i);

			}
		}
	}
	else{

		for (int i = 0; i < Npoints; i++){
			if (IsFirstOrSecondState){
				dDiff += mSecondOrderStencils[x - lStarter][i] * mObjectiveFunction[x - lStarter + i][y];
				//dDiff += mStencil2ndX[i] * mObjectiveFunction(x - lStarter + i, y);

			}
			else{
				dDiff += mSecondOrderStencils[x - lStarter][i] * mObjectiveFunction[x][y - lStarter + i];
				//dDiff += mStencil2ndX[i] * mObjectiveFunction(x, y - lStarter + i);

			}
		}
	}
	return dDiff;
}



void ADI2DScheme::generateRHSPredictor(dVector &outVector, double dt, int index, int pivot)
{
	int Nindex = mGrid->getGridSize(index);
	dVector stateFactors(2);
	dVector dX, dY;
	double dDiffConvectionX, dDiffConvectionY;
	double dDiffDiffusionX, dDiffDiffusionY;
	double dConvectionX, dDiffusionX;
	double dConvectionY, dDiffusionY;
	double dSource;
	double uij,vij;
	bool lIsXUniform=mGrid->getIsUniformGrid(0);
	bool lIsYUniform = mGrid->getIsUniformGrid(1);
	int lVariableIndex;
	int Npoints;
	int lStarter;
	int lStarterPivot;
	bool bFlag = (mStencilConvections[index] == FivePoints) | (mStencilDiffusions[index] == FivePoints);
	
	(bFlag ? Npoints = 5 : Npoints = 3);
	(bFlag ? lStarter = 2 : lStarter = 1);
	if (index == 0)
	{
		stateFactors[1] = mGrid->getGridElement(1, pivot);
		dDiffDiffusionY = 0;//temp
		dDiffusionY = 0;//temp
		for (int j = 0; j < Nindex - Npoints+1; j++){
			lVariableIndex = lStarter+ j;
			stateFactors[index] = mGrid->getGridElement(index, lVariableIndex);
			if (index != 0){
				mPDEmodel.setEntity(mEntities[mIndexPivot][lVariableIndex]);
			}
			else{
				mPDEmodel.setEntity(mEntities[lVariableIndex][mIndexPivot]);
			}
			mPDEmodel.updateConvections(stateFactors);
			mPDEmodel.updateDiffusions(stateFactors);
			mPDEmodel.updateSourceTerms(stateFactors); // redundant : update juste un facteur

			dConvectionX = mPDEmodel.getConvection(0);
			dDiffusionX = mPDEmodel.getDiffusion(0);
			dConvectionY = mPDEmodel.getConvection(1);
			//dDiffusionY = mPDEmodel->getDiffusion(1);
			dSource = mPDEmodel.getSourceTerms(0);
			
			//if (lIsXUniform & lIsYUniform){
			if (false){

				dDiffConvectionX = uniformDifferentiator2D(lVariableIndex, pivot, mStencilConvections[0], true, true);
				dDiffDiffusionX = uniformDifferentiator2D(lVariableIndex, pivot, mStencilDiffusions[0], true, false); // redundant calculation
				dDiffConvectionY = uniformDifferentiator2D(lVariableIndex, pivot, mStencilConvections[1], false, true);
				//dDiffDiffusionY = differentiator2D(lVariableIndex, pivot, lIsYUniform, mStencilDiffusions[1], false, false);
			}
			else{
				dDiffConvectionX = newDifferentiator2D(lVariableIndex, pivot, mStencilConvections[0], true, true);
				dDiffDiffusionX = newDifferentiator2D(lVariableIndex, pivot, mStencilDiffusions[0], true, false); // redundant calculation
				dDiffConvectionY = newDifferentiator2D(lVariableIndex, pivot, mStencilConvections[1], false, true);
				//dDiffDiffusionY = differentiator2D(lVariableIndex, pivot, lIsYUniform, mStencilDiffusions[1], false, false);
			}
			

			uij = mObjectiveFunction[lVariableIndex][pivot];
			//uij = mObjectiveFunction(lVariableIndex,pivot);

			outVector[j] = uij
				+ dt*(1 - mThetaScheme)*
				(-0.5*dSource*uij
				+ dConvectionX*dDiffConvectionX
				+ 0.5*dDiffusionX*dDiffusionX*dDiffDiffusionX)
				+ dt*(-0.5*dSource*uij
				+ dConvectionY*dDiffConvectionY
				+ 0.5*dDiffusionY*dDiffusionY*dDiffDiffusionY);
		}
	}
	else{
		bool bFlagPivot = (mStencilConvections[0] == FivePoints) | (mStencilDiffusions[0] == FivePoints);
		(bFlagPivot ? lStarterPivot = 2 : lStarterPivot = 1);
		stateFactors[0] = mGrid->getGridElement(0, pivot);
		dDiffDiffusionY = 0;//temp
		for (int j = 0; j < Nindex - Npoints + 1; j++){
			lVariableIndex = lStarter + j;
			stateFactors[index] = mGrid->getGridElement(index, lVariableIndex);
			if (index != 0){
				mPDEmodel.setEntity(mEntities[mIndexPivot][lVariableIndex]);
			}
			else{
				mPDEmodel.setEntity(mEntities[lVariableIndex][mIndexPivot]);
			}
			/*mPDEmodel.updateConvections(stateFactors);
			mPDEmodel.updateDiffusions(stateFactors);
			mPDEmodel.updateSourceTerms(stateFactors);*/

			mPDEmodel.updateConvections(stateFactors,index);
			mPDEmodel.updateDiffusions(stateFactors);
			

			dConvectionY = mPDEmodel.getConvection(index);
			dDiffusionY = mPDEmodel.getDiffusion(index);
			dSource = mPDEmodel.getSourceTerms(index);
			if (false){
			//	if (lIsXUniform & lIsYUniform){

				dDiffConvectionY = uniformDifferentiator2D(pivot, lVariableIndex, mStencilConvections[index], false, true);
			}
			else{
				dDiffConvectionY = newDifferentiator2D(pivot, lVariableIndex, mStencilConvections[index], false, true);

			}
			//dDiffDiffusionY = differentiator2D(pivot, lVariableIndex, lIsYUniform, mStencilDiffusions[index], false, false);
			vij = mPredictors[0][pivot-lStarterPivot][j];
			uij = mObjectiveFunction[pivot][lVariableIndex];
			//uij = mObjectiveFunction(pivot,lVariableIndex);

			outVector[j] = vij
				- dt*(mThetaScheme)*
				(-(1/(1.0*mPDEDimension))*dSource*uij
				+ dConvectionY*dDiffConvectionY
				+ 0.5*dDiffusionY*dDiffusionY*dDiffDiffusionY);	
		}

	}
}
