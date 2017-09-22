#ifndef ADISCHEME_H
#define ADISCHEME_H
#include "PDEScheme.h"
#pragma once
class ADIScheme : public PDEScheme
{
public:
	
	
	ADIScheme(MarkovianProductBridge inMarkovianProduct,
		PDEModelBridge  inPDEmodel,
		boost::shared_ptr<Grid> inGrid,
		TimeGrid & inTimeGrid,
		int inTimeSteps,
		double inSchemeType,
		PDVector inStencilConvections,
		PDVector inStencilDiffusions,
		PDVector inStencilMixedDerivatives,
		bool inUseCorrector,
		double inLambda);

	void generatePredictorBandMatrix(dMatrix & outPredictor, double dt,int inStep,const std::vector<double> & pivot);
	virtual void generateCorrectorBandMatrix(dMatrix & outCorrector, int inStep)=0; 
	virtual void generateRHSPredictor(dVector &outVector, double dt, int index, int pivot)=0;
	virtual void generateRHSCorrector()=0;
	virtual void solveState(int index, int pivot, double dt) = 0;
	virtual void solve(int inIndextime)=0;
	virtual void fillOFwithPayoff() = 0;
	void generateDX();
	
	void generateFirstStencils();
	void generateSecondStencils();
	virtual ~ADIScheme(){};
protected :
	std::vector<dMatrix> mPredictors;
    std::vector<dMatrix> mCorrectors;
	std::vector<double> mDurations;

	std::vector<std::vector<std::vector<double>>> mdX;
	std::vector<std::vector<double>> mVol;
	std::vector<std::vector<std::vector<double>>> mFirstOrderStencils;
	std::vector<std::vector<double>> mSecondOrderStencils;
	dMatrix mEntities;
	bool mUseCorrector;
	double mLambda;
	int mIndexPivot;
	double mUniformdX;
	double mUniformdY;
	std::vector<double> mStencil1STX;
	std::vector<double> mStencil1STY;
	std::vector<double> mStencil2ndX;

};

#endif 

