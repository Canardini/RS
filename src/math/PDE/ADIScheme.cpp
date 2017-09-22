#include <Math/PDE/ADIScheme.h>

ADIScheme::ADIScheme(MarkovianProductBridge inMarkovianProduct,
					 PDEModelBridge inPDEmodel,
	                 boost::shared_ptr<Grid> inGrid,
					 TimeGrid &inTimeGrid,
	                 int inTimeSteps,
	                 double inSchemeType,
					 PDVector inStencilConvections,
					 PDVector inStencilDiffusions,
					 PDVector inStencilMixedDerivatives,
					 bool inUseCorrector,
					 double inLambda):
					 PDEScheme(inMarkovianProduct, inPDEmodel, inGrid, inTimeGrid, inTimeSteps, inSchemeType, inStencilConvections, inStencilDiffusions, inStencilMixedDerivatives), mUseCorrector(inUseCorrector),
mLambda(inLambda)
{
}

void ADIScheme::generateDX(){
	mdX.resize(2);
	mdX[0].resize(mGrid->getGridSize(0) - 2);
	for (int i = 0; i < mGrid->getGridSize(0) - 2; i++){
		mdX[0][i].resize(2);
		mdX[0][i][0] = mGrid->getGridElement(0, i + 1) - mGrid->getGridElement(0, i);
		mdX[0][i][1] = mGrid->getGridElement(0, i + 2) - mGrid->getGridElement(0, i + 1);
	}

	mdX[1].resize(mGrid->getGridSize(1) - 4);
	for (int i = 0; i < mGrid->getGridSize(1) - 4; i++){
		mdX[1][i].resize(4);
		mdX[1][i][0] = mGrid->getGridElement(1, i + 2) - mGrid->getGridElement(1, i);
		mdX[1][i][1] = mGrid->getGridElement(1, i + 2) - mGrid->getGridElement(1, i + 1);
		mdX[1][i][2] = mGrid->getGridElement(1, i + 3) - mGrid->getGridElement(1, i + 2);
		mdX[1][i][3] = mGrid->getGridElement(1, i + 4) - mGrid->getGridElement(1, i + 2);
	}

}



void ADIScheme::generateFirstStencils(){

	
	mFirstOrderStencils.resize(2);
	mFirstOrderStencils[0].resize(mGrid->getGridSize(0) - 2);
	bool lX = mGrid->getIsUniformGrid(0);
	bool lY = mGrid->getIsUniformGrid(1);

	for (int i = 0; i < mGrid->getGridSize(0) - 2; i++){
		mFirstOrderStencils[0][i].resize(3);
		stencilFirstOrderCoefficients(mdX[0][i], lX, ThreePoints, mFirstOrderStencils[0][i]);
	}

	mFirstOrderStencils[1].resize(mGrid->getGridSize(1) - 4);
	for (int i = 0; i < mGrid->getGridSize(1) - 4; i++){
		mFirstOrderStencils[1][i].resize(5);
		stencilFirstOrderCoefficients(mdX[1][i], lY, FivePoints, mFirstOrderStencils[1][i]);
	}
}

void ADIScheme::generateSecondStencils(){

	bool lX = mGrid->getIsUniformGrid(0);
	bool lY = mGrid->getIsUniformGrid(1);
	mSecondOrderStencils = std::vector<std::vector<double>>( std::vector<std::vector<double>>(mGrid->getGridSize(0) - 2, std::vector<double>(3)));
	for (int i = 0; i < mGrid->getGridSize(0) - 2; i++){
		stencilSecondOrderCoefficients(mdX[0][i],lX, ThreePoints, mSecondOrderStencils[i]); // Erreur, pas LY mais LX
	}
}

 void ADIScheme::generatePredictorBandMatrix(dMatrix &outMatrix, double dt,int index,const std::vector<double> & pivot){

	 double dConvection;
	 double dDiffusion;
	 double dSource;
	 int Nindex = mGrid->getGridSize(index);
	 bool bFlag = (mStencilConvections[index] == FivePoints) | (mStencilDiffusions[index] == FivePoints); 
	 bool lIsUniformgrid = mGrid->getIsUniformGrid(0) & mGrid->getIsUniformGrid(1); //Assume 2d scheme
	 // 5-points
	 if (bFlag){
		 dVector stateFactors(pivot);
		 for (int j = 0; j < 5; j++){
			 for (int i = 0; i < Nindex - 4; i++)
			 {
				 stateFactors[index] = mGrid->getGridElement(index, i + 2);
				 if (index != 0){
					 mPDEmodel.setEntity(mEntities[mIndexPivot][i + 2]);
				 }
				 else{
					 mPDEmodel.setEntity(mEntities[i + 2][mIndexPivot]);
				 }
				 //mPDEmodel.updateConvections(stateFactors);
				 mPDEmodel.updateConvections(stateFactors,index);

				 mPDEmodel.updateDiffusions(stateFactors);
				 
				 dConvection = mPDEmodel.getConvection(index);
				 dDiffusion = mPDEmodel.getDiffusion(index);

				 if (true){
				//	 if (!lIsUniformgrid){

					 if (dConvection == 0){
						 outMatrix[i][j] = 0.0;
					 }
					 else{
						 outMatrix[i][j] = dConvection*mFirstOrderStencils[index][i][j];

					 }
					 if (dDiffusion != 0){
						 outMatrix[i][j] += 0.5*dDiffusion*dDiffusion*mSecondOrderStencils[i][j];
					 }
				 }
				 else{
					 if (dConvection == 0){
						 outMatrix[i][j] = 0.0;
					 }
					 else{
						 if (index == 0){
							 outMatrix[i][j] = dConvection*mStencil1STX[j];
						 }
						 else{
							 outMatrix[i][j] = dConvection*mStencil1STY[j];
						 }
					 }
					 if (dDiffusion != 0){
						 if (index == 0){
							 outMatrix[i][j] += 0.5*dDiffusion*dDiffusion*mStencil2ndX[j];
						 }
						 else{
							 outMatrix[i][j] += 0.5*dDiffusion*dDiffusion*mStencil2ndX[j]; // no diffusion in Y

						 }
					 }

				 }

				 outMatrix[i][j] *= -mThetaScheme*dt;
				 if (j==2){
					 mPDEmodel.updateSourceTerms(stateFactors);
					 dSource = mPDEmodel.getSourceTerms(index);
					 outMatrix[i][j] += (1 / (1.0*mPDEDimension))*dSource*mThetaScheme*dt;
					 outMatrix[i][j] += 1;
				 }
			 }
		 }

		 //Numerical boundary condition : LowerBound
		 outMatrix[0][2] += 2 * outMatrix[0][0] + 2 * outMatrix[0][1];
		 outMatrix[0][3] -= outMatrix[0][1] ;
		 outMatrix[0][4] -= outMatrix[0][0];

		 outMatrix[1][1] += 2 * outMatrix[1][0];
		 outMatrix[1][2] -= outMatrix[1][0];

		 //Numerical boundary condition : UpperBound
		 outMatrix[Nindex - 5][0] -= outMatrix[Nindex - 5][4];
		 outMatrix[Nindex - 5][1] -= outMatrix[Nindex - 5][3];
		 outMatrix[Nindex - 5][2] += 2 * outMatrix[Nindex - 5][3] + 2  * outMatrix[Nindex - 5][4];

		 outMatrix[Nindex - 6][2] -= outMatrix[Nindex - 6][4];
		 outMatrix[Nindex - 6][3] += 2*outMatrix[Nindex - 6][4];
   
	 }
	 else{ // 3-points
		 
		 dVector a(Nindex-2);
		 dVector b(Nindex-2);
		 dVector c(Nindex-2);
		 dVector dCoeffConvection(3);
		 dVector dCoeffDiffusion(3);
		 dVector dX(2);
		 dVector stateFactors(pivot);

		 for (int i = 0; i <Nindex-2; i++){
			 stateFactors[index] = mGrid->getGridElement(index, i+1);
			 if (index != 0){
				 mPDEmodel.setEntity(mEntities[mIndexPivot][i + 2]);
			 }
			 else{
				 mPDEmodel.setEntity(mEntities[i + 2][mIndexPivot]);
			 }
			 mPDEmodel.updateConvections(stateFactors);
			 mPDEmodel.updateDiffusions(stateFactors);
			 mPDEmodel.updateSourceTerms(stateFactors);
			 dConvection = mPDEmodel.getConvection(index); // to optimize
			 dDiffusion = mPDEmodel.getDiffusion(index);
			 dSource = mPDEmodel.getSourceTerms(index);
			 //if (!lIsUniformgrid){

			 if (true){
				 a[i] = dConvection*mFirstOrderStencils[index][i][0] + 0.5*dDiffusion*dDiffusion*mSecondOrderStencils[i][0];
				 a[i] *= -mThetaScheme*dt;

				 b[i] = -(1 / (1.0*mPDEDimension))*dSource + dConvection*mFirstOrderStencils[index][i][1]
					 + 0.5*dDiffusion*dDiffusion*mSecondOrderStencils[i][1];
				 b[i] *= -mThetaScheme*dt;
				 b[i] += 1;

				 c[i] = dConvection*mFirstOrderStencils[index][i][2] + 0.5*dDiffusion*dDiffusion*mSecondOrderStencils[i][2];
				 c[i] *= -mThetaScheme*dt;
			 }
			 else{
				 if (index == 0){
					 a[i] = dConvection*mStencil1STX[0] + 0.5*dDiffusion*dDiffusion*mStencil2ndX[0];
					 a[i] *= -mThetaScheme*dt;


					 b[i] = -(1 / (1.0*mPDEDimension))*dSource + dConvection*mStencil1STX[1]
						 + 0.5*dDiffusion*dDiffusion*mStencil2ndX[1];
					 b[i] *= -mThetaScheme*dt;
					 b[i] += 1;

					 c[i] = dConvection*mStencil1STX[2] + 0.5*dDiffusion*dDiffusion*mStencil2ndX[2];
					 c[i] *= -mThetaScheme*dt;
				 }
			 }
		 }
		 //Numerical Boundary conditions
		 c[0] -= a[0];
		 b[0] += 2*a[0];
		 a[Nindex - 3] -= c[Nindex - 3];
		 b[Nindex - 3] += 2 * c[Nindex - 3];
		 
		 outMatrix[0] = a;
		 outMatrix[1] = b;
		 outMatrix[2] = c;
	 }
 }

 
