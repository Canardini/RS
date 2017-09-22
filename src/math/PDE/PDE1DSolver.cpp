#include <Math/PDE/PDE1DSolver.h>

PDE1DSolver::PDE1DSolver(MarkovianProductBridge inMarkovianProduct,
	PDEModelBridge  inPDEmodel,
	boost::shared_ptr<Grid> inGrid,
	TimeGrid & inTimeGrid,
	int inTimeSteps,
	double inSchemeType,
	PDVector inStencilConvections,
	PDVector inStencilDiffusions,
	PDVector inStencilMixedDerivatives,
	std::string & inProductType
	) :
	PDEScheme(inMarkovianProduct, inPDEmodel, inGrid, inTimeGrid, inTimeSteps, inSchemeType, inStencilConvections, inStencilDiffusions, inStencilMixedDerivatives), mProductType(inProductType)
{

	generateDX();
	generateFirstStencils();
	generateSecondStencils();
	mObjectiveFunction.resize(mGrid->getGridSize(0));
	mPrices.resize(mGrid->getGridSize(0));
}

void PDE1DSolver::generateDX(){ 
	mdX.resize(mGrid->getGridSize(0) - 2);
	for (int i = 0; i < mGrid->getGridSize(0) - 2; i++){
		mdX[i].resize(2);
		mdX[i][0] = mGrid->getGridElement(0, i + 1) - mGrid->getGridElement(0, i);
		mdX[i][1] = mGrid->getGridElement(0, i + 2) - mGrid->getGridElement(0, i + 1);
	}
}



void PDE1DSolver::generateFirstStencils(){
	 
	mFirstOrderStencils.resize(mGrid->getGridSize(0) - 2);
	bool lX = mGrid->getIsUniformGrid(0);
	for (int i = 0; i < mGrid->getGridSize(0) - 2; i++){
		mFirstOrderStencils[i].resize(3);
		stencilFirstOrderCoefficients(mdX[i], lX, ThreePoints, mFirstOrderStencils[i]);
	}
}

void PDE1DSolver::generateSecondStencils(){
	bool lX = mGrid->getIsUniformGrid(0);
	mSecondOrderStencils = std::vector<std::vector<double>>(std::vector<std::vector<double>>(mGrid->getGridSize(0) - 2, std::vector<double>(3)));
	for (int i = 0; i < mGrid->getGridSize(0) - 2; i++){
		stencilSecondOrderCoefficients(mdX[i], lX, ThreePoints, mSecondOrderStencils[i]);
	}
}

void PDE1DSolver::generatePredictorBandMatrix(dMatrix &outMatrix, double dt){

	double dConvection;
	double dDiffusion;
	double dSource;
	int Nindex = mGrid->getGridSize(0);
	bool bFlag = false;
	bool lIsUniformgrid = mGrid->getIsUniformGrid(0); //Assume 2d scheme

	dVector a(Nindex - 2);
	dVector b(Nindex - 2);
	dVector c(Nindex - 2);
	dVector dCoeffConvection(3);
	dVector dCoeffDiffusion(3);
	dVector dX(2);
	dVector stateFactors(1);

	for (int i = 0; i < Nindex - 2; i++){
		stateFactors[0] = mGrid->getGridElement(0, i + 1);
		mPDEmodel.updateConvections(stateFactors);
		mPDEmodel.updateDiffusions(stateFactors);
		mPDEmodel.updateSourceTerms(stateFactors);
		dConvection = mPDEmodel.getConvection(0); // to optimize
		dDiffusion = mPDEmodel.getDiffusion(0);
		dSource = mPDEmodel.getSourceTerms(0);

		a[i] = dConvection*mFirstOrderStencils[i][0] + 0.5*dDiffusion*dDiffusion*mSecondOrderStencils[i][0];
		a[i] *= -mThetaScheme*dt;

		b[i] = -dSource + dConvection*mFirstOrderStencils[i][1]
			+ 0.5*dDiffusion*dDiffusion*mSecondOrderStencils[i][1];
		b[i] *= -mThetaScheme*dt;
		b[i] += 1;

		c[i] = dConvection*mFirstOrderStencils[i][2] + 0.5*dDiffusion*dDiffusion*mSecondOrderStencils[i][2];
		c[i] *= -mThetaScheme*dt;

	}
	//Numerical Boundary conditions
	c[0] -= a[0];
	b[0] += 2 * a[0];
	a[Nindex - 3] -= c[Nindex - 3];
	b[Nindex - 3] += 2 * c[Nindex - 3];

	outMatrix[0] = a;
	outMatrix[1] = b;
	outMatrix[2] = c;
}

void PDE1DSolver::fillOFwithPayoff(){
	std::vector<double> lStates(1);
	for (int i = 0; i < mGrid->getGridSize(0); i++){
		lStates[0] = mGrid->getGridElement(0, i);
		if (i == 0){
			mMarkovianProduct.updateMarkovianFactors(lStates);
			//mObjectiveFunction[i] = std::max(mMarkovianProduct.getPrice(), 0.0);
			mObjectiveFunction[i] = mMarkovianProduct.getPrice();
			//mObjectiveFunction[i] = 1.0;


		}
		else{
			//mObjectiveFunction[i] = std::max(mMarkovianProduct.updatePrice(lStates), 0.0);
			//mObjectiveFunction[i] = mMarkovianProduct.getPrice();
			/*mMarkovianProduct.updateMarkovianFactors(lStates);
			mObjectiveFunction[i] = mMarkovianProduct.getPrice();*/

			mObjectiveFunction[i] = mMarkovianProduct.updatePrice(lStates);
			//mObjectiveFunction[i] = 1.0;

		}
	}
}

void PDE1DSolver::fillPrices(){
	std::vector<double> lStates(1);
	for (int i = 0; i < mGrid->getGridSize(0); i++){
		lStates[0] = mGrid->getGridElement(0, i);
		if (i == 0){
			mMarkovianProduct.updateMarkovianFactors(lStates);
			mPrices[i] = std::max(mMarkovianProduct.getPrice(), 0.0);
		}
		else{
			mPrices[i] = std::max(mMarkovianProduct.updatePrice(lStates), 0.0);
		}
	}
}

void PDE1DSolver::applyBoundaryConditions(){
	int nX = mObjectiveFunction.size();
	mObjectiveFunction[0] = 2 * mObjectiveFunction[1] - mObjectiveFunction[2];
	mObjectiveFunction[nX - 1] = 2 * mObjectiveFunction[nX - 2] - mObjectiveFunction[nX - 3];
}
void PDE1DSolver::solve(int compteur){

	int lcompteur = compteur;
	int ltimeIndex = mTimeGrid.getIndexKeyDates(compteur);
	mMarkovianProduct.changeStructure(mTimeGrid.getKeyDate(compteur), mTimeGrid[ltimeIndex]);
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
		mPDEmodel.setIndex(compteur);
		mPDEmodel.setInitSource(mMarkovianProduct.getInitSource(mTimeGrid[i]));
		mPDEmodel.setDriftAdjust(mMarkovianProduct.getDrifAdjust(mTimeGrid[i]));

		solveState(dt);
		applyBoundaryConditions();

		if (mProductType != "European"){
			if (mTimeGrid.getIsKeyDate(i)){
				fillPrices();
				fillExercisePayOff();
			}
		}
	}
}

void PDE1DSolver::fillExercisePayOff(){
	for (int i = 0; i < mGrid->getGridSize(0); i++){
		mObjectiveFunction[i] = std::max(mObjectiveFunction[i], mPrices[i]);
	}
}


void PDE1DSolver::solveState(double dt){
	dMatrix A;
	dMatrix al;
	dVector dPivot(2);
	dVector b;
	dVector u;
	int Nindex;

	int Npoints;
	Nindex = mGrid->getGridSize(0);
	int lStarter;
	bool bFlag = false;
	(bFlag ? Npoints = 5 : Npoints = 3);
	(bFlag ? lStarter = 2 : lStarter = 1);


	A.resize(3, dVector(Nindex - 2));
	b.resize(Nindex - 2);
	u.resize(Nindex - 2);

	generatePredictorBandMatrix(A, dt);;//generate the band and the rhs in the same time.
	generateRHSPredictor(b, dt);
	BandSolver::dtridag(A[0], A[1], A[2], b, u);
	for (int i = 0; i < Nindex - 2; i++){
		mObjectiveFunction[i + 1] = u[i];
	}
}

void PDE1DSolver::generateRHSPredictor(dVector &outVector, double dt)
{
	int Nindex = mGrid->getGridSize(0);
	dVector stateFactors(1);
	dVector dX, dY;
	double dDiffConvectionX;
	double dDiffDiffusionX;
	double dConvectionX, dDiffusionX;
	double dSource;
	double uij;
	bool lIsXUniform = mGrid->getIsUniformGrid(0);
	int lVariableIndex;
	int Npoints;
	int lStarter;

	bool bFlag = false;

	(bFlag ? Npoints = 5 : Npoints = 3);
	(bFlag ? lStarter = 2 : lStarter = 1);

	for (int j = 0; j < Nindex - Npoints + 1; j++){
		lVariableIndex = lStarter + j;
		stateFactors[0] = mGrid->getGridElement(0, lVariableIndex);

		mPDEmodel.updateConvections(stateFactors);
		mPDEmodel.updateDiffusions(stateFactors);
		mPDEmodel.updateSourceTerms(stateFactors); // redundant : update juste un facteur

		dConvectionX = mPDEmodel.getConvection(0);
		dDiffusionX = mPDEmodel.getDiffusion(0);
		dSource = mPDEmodel.getSourceTerms(0);


		dDiffConvectionX = newDifferentiator2D(lVariableIndex, mStencilConvections[0], true);
		dDiffDiffusionX = newDifferentiator2D(lVariableIndex, mStencilDiffusions[0], false);

		uij = mObjectiveFunction[lVariableIndex];

		outVector[j] = uij
			+ dt*(1 - mThetaScheme)*
			(-dSource*uij
			+ dConvectionX*dDiffConvectionX
			+ 0.5*dDiffusionX*dDiffusionX*dDiffDiffusionX);

	}
}

double PDE1DSolver::newDifferentiator2D(int x,
	PointDiscretization type, bool IsFirstOrSecondDerivative)
{
	double dDiff = 0.0;
	int Npoints;
	int lStarter;
	(type == ThreePoints ? Npoints = 3 : Npoints = 5);
	(type == ThreePoints ? lStarter = 1 : lStarter = 2);

	if (IsFirstOrSecondDerivative){

		for (int i = 0; i < Npoints; i++){
				dDiff += mFirstOrderStencils[x - lStarter][i] * mObjectiveFunction[x - lStarter + i];
		}
	}
	else{

		for (int i = 0; i < Npoints; i++){

			dDiff += mSecondOrderStencils[x - lStarter][i] * mObjectiveFunction[x - lStarter + i];
		}
	}
	return dDiff;
}