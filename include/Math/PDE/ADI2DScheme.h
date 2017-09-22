#ifndef ADI2DSCHEME_H
#define ADI2DSCHEME_H


#pragma once
#include <Math/PDE/ADIScheme.h>
#include<Model\IRModel\MarkovianIRModelBridge.h>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/thread.hpp>
#include <boost/phoenix.hpp>
#include <boost/ref.hpp>
#include<Math/AnalyticalFormulas/Black.h>

class ADI2DScheme :
	public ADIScheme
{
public:


	

	ADI2DScheme(MarkovianProductBridge inMarkovianProduct,
				PDEModelBridge  inPDEmodel,
				boost::shared_ptr<Grid> inGrid,
				TimeGrid &inTimeGrid,
				int inTimeSteps,
				double inSchemeType,
				PDVector inStencilConvections,
				PDVector inStencilDiffusions,
				PDVector inStencilMixedDerivatives,
				bool inUseCorrector,
				double inLambda,
				std::string & inProductType);

	
	virtual void solveState(int index,int pivot, double dt) ;
	void  solveStateThr(int index, int pivot, double dt,int inNumberofThreads);
	virtual void solve(int index);
	void solveThr(int timeIndex);

	double differentiator2D(int x, int y,bool isUniform, PointDiscretization type,
		                    bool IsFirstOrSecondState, bool IsFirstOrSecondDerivative);
	double newDifferentiator2D(int x, int y, PointDiscretization type,
		bool IsFirstOrSecondState, bool IsFirstOrSecondDerivative);
	double uniformDifferentiator2D(int x, int y, PointDiscretization type,
		bool IsFirstOrSecondState, bool IsFirstOrSecondDerivative);
	virtual void generateRHSPredictor(dVector &outVector, double dt, int index, int pivot);
	virtual void generateRHSCorrector(){};
	//double getObjectiveFunctionValue(int index1, int index2){ return mObjectiveFunction(index1,index2); };
	double getObjectiveFunctionValue(int index1, int index2){ return mObjectiveFunction[index1][index2]; };

	virtual void generateCorrectorBandMatrix(dMatrix & outCorrector, int inStep);
	virtual void generateLowerBoundCondition(int index, PointDiscretization inPD, BoundaryConditions inBC, std::vector<double> outVec);
	virtual void generateUpperBoundCondition(int index, PointDiscretization inPD, BoundaryConditions inBC, std::vector<double> outVec);
	void applyBoundaryConditions(BoundaryConditions type);
	virtual void fillOFwithPayoff();
	void fillOFwithPayoff(int inThread);
	void multiThreadPayoff();
	void generateVol();
	void fillEntities();
	void fillExercisePayOff();
	void fillEntities(int inThread);
	void multiThreadEntity();
	std::vector<double> getProbaComponents(){ return mProbaComponents; };

	virtual ~ADI2DScheme(){};
private: 
	dMatrix mObjectiveFunction;
	dMatrix mPrices;

	/*boost::numeric::ublas::matrix<double> mObjectiveFunction;
	boost::numeric::ublas::matrix<double> mPrices;*/
	int mNumberOfCPUs;
	std::vector<MarkovianProductBridge> mMarkovianProductVector;
	std::string mProductType;
	std::vector<std::vector<int>> mExerciseFlags;

	double mAverageDistribution;
	double mVarianceDistribution;
	std::vector<double> mProbaComponents;

};

#endif

