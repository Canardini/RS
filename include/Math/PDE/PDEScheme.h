#ifndef PDE_SCHEME_H
#define PDE_SCHEME_H

#pragma once
#include <boost/make_shared.hpp>
#include <numeric>

#include <Math/Solver/BandSolver.h>
#include <Math/Matrix/Matrix.h>
#include <Product\MarkovianProductBridge.h>
#include <Math/PDE/Grid.h>
#include <Math\PDEModel\PDEModelBridge.h>

class PDEScheme
{
public:
	
	enum PointDiscretization { ThreePoints, FivePoints };
	enum BoundaryConditions { Numerical,Neumann,Dirichlet,Robbins };

	typedef boost::shared_ptr<Grid> tGridPtr;
	typedef boost::shared_ptr<TimeGrid> tTimeGridPtr;
	typedef std::vector<PointDiscretization> PDVector;


	virtual void generateLowerBoundCondition(int index, PointDiscretization inPD, BoundaryConditions inBC, dVector outVec)=0;
	virtual void generateUpperBoundCondition(int index, PointDiscretization inPD, BoundaryConditions inBC, dVector outVec)=0;

	void stencilFirstOrderCoefficients(const dVector & dx, bool IsUniform,
		                               PointDiscretization type, 
									   std::vector<double> &outVec);
	void stencilFirstOrderCoefficients(const std::vector<double> & dx, bool IsUniform,
		PointDiscretization type, int index, double & outVec);
	void stencilSecondOrderCoefficients(const dVector & dx, bool IsUniform,
		                                PointDiscretization type, 
										std::vector<double> &outVec);
	double differentiator(const dVector & u,
		                  const dVector & dx, 
						  int index,
		                  PointDiscretization type, 
		                  bool isFirstOrSecond);

	virtual ~PDEScheme(){};

protected:
	explicit PDEScheme(MarkovianProductBridge inMarkovianProduct,
		               PDEModelBridge inPDEmodel,
		               tGridPtr inGrid,
					   TimeGrid &inTimeGrid,
		               int inTimeSteps,
		               double inSchemeType,
					   PDVector inStencilConvections,
					   PDVector inStencilDiffusions,
					   PDVector inStencilMixedDerivatives);

	
	MarkovianProductBridge mMarkovianProduct;
	PDEModelBridge mPDEmodel;
	tGridPtr mGrid;
	TimeGrid mTimeGrid;
	
	PDVector mStencilConvections;
	PDVector mStencilDiffusions;
	PDVector mStencilMixedDerivatives;
	
	int mPDEDimension;
	int mTimeSteps;
	double mThetaScheme; 

};
#endif // !PDE_SCHEME_H