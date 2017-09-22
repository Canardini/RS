#ifndef PDE_ONED_SOLVER_H
#define PDE_ONED_SOLVER_H
#include "PDEScheme.h"
#pragma once
class PDE1DSolver : public PDEScheme
{
public:


	PDE1DSolver(MarkovianProductBridge inMarkovianProduct,
		PDEModelBridge  inPDEmodel,
		boost::shared_ptr<Grid> inGrid,
		TimeGrid & inTimeGrid,
		int inTimeSteps,
		double inSchemeType,
		PDVector inStencilConvections,
		PDVector inStencilDiffusions,
		PDVector inStencilMixedDerivatives,
		std::string & inProductType
	);

	void generatePredictorBandMatrix(dMatrix & outPredictor, double dt);
	virtual void generateRHSPredictor(dVector &outVector, double dt);
	virtual void generateLowerBoundCondition(int index, PointDiscretization inPD, BoundaryConditions inBC, std::vector<double> outVec){};
	virtual void generateUpperBoundCondition(int index, PointDiscretization inPD, BoundaryConditions inBC, std::vector<double> outVec){};

    void solve(int inIndextime);
	void fillOFwithPayoff();

	void generateDX();
	void fillExercisePayOff();
	void applyBoundaryConditions();
	void solveState(double dt);
	void generateFirstStencils();
	void generateSecondStencils();
	void fillPrices();
	double newDifferentiator2D(int x,PointDiscretization type, bool IsFirstOrSecondDerivative);
	double getObjectiveFunctionValue(int index){ return mObjectiveFunction[index]; };
	virtual ~PDE1DSolver(){};
protected:
	std::vector<dMatrix> mPredictors;
	std::vector<dMatrix> mCorrectors;
	std::vector<double> mDurations;

	std::vector<std::vector<double>> mdX;



	std::vector<std::vector<double>> mVol;
	std::vector<std::vector<double>> mFirstOrderStencils;
	std::vector<std::vector<double>> mSecondOrderStencils;
	std::vector<double> mObjectiveFunction;
	std::vector<double> mPrices;
	std::string mProductType;
	dMatrix mEntities;
	bool mUseCorrector;
	double mLambda;
	int mIndexPivot;
	double mUniformdX;
 
};

#endif 

