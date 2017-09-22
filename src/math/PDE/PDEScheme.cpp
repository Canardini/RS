#include <Math/PDE/PDEScheme.h>



PDEScheme::PDEScheme(MarkovianProductBridge inMarkovianProduct,
					 PDEModelBridge inPDEmodel,
	                 boost::shared_ptr<Grid> inGrid,
					 TimeGrid & inTimeGrid, 
	                 int inTimeSteps,
	                 double inThetaScheme,
					 PDVector inStencilConvections,
					 PDVector inStencilDiffusions,
					 PDVector inStencilMixedDerivatives) :
					 mMarkovianProduct(inMarkovianProduct), mPDEmodel(inPDEmodel),
					 mGrid(boost::make_shared<Grid>(*inGrid)), mTimeGrid(inTimeGrid), mTimeSteps(inTimeSteps), mThetaScheme(inThetaScheme),
					 mStencilConvections(inStencilConvections), mStencilDiffusions(inStencilDiffusions), mStencilMixedDerivatives(inStencilMixedDerivatives)

 {}

void PDEScheme::stencilFirstOrderCoefficients(const std::vector<double> & dx, bool IsUniform,
	PointDiscretization type, std::vector<double> &outVec)
{
	
	if (type == ThreePoints){
		outVec[0] = -(dx[1]/ (dx[0]*(dx[0]+dx[1])));
		outVec[1] = (dx[1]-dx[0])/(dx[0]*dx[1]);
		outVec[2] = (dx[0] / (dx[1] * (dx[0] + dx[1])));
	}
	else
	{
		if (IsUniform){
			outVec[0] = (1 / (12.0*dx[1]));
			outVec[1] = -8.0*outVec[0];
			outVec[2] = 0.0;
			outVec[3] = -outVec[1];
			outVec[4] = -outVec[0];
		}
		else{
			double dDpp, dDp, dDm, dDmm;
			double dDpp2, dDp2, dDm2, dDmm2;
			double dDpp4, dDp4, dDm4, dDmm4;
			double dGamma01, dGamma02;
			double dGamma11, dGamma12, dGamma32;
			double dGamma21, dGamma22, dGamma31;
			 
			dDpp = dx[3];
			dDp = dx[2];
			dDm = dx[1];
			dDmm = dx[0];

			dDpp2 = dDpp*dDpp;
			dDp2 = dDp*dDp;
			dDm2 = dDm*dDm;
			dDmm2 = dDmm*dDmm;

			dDpp4 = dDpp2*dDpp2;
			dDp4 = dDp2*dDp2;
			dDm4 = dDm2*dDm2;
			dDmm4 = dDmm2*dDmm2;

			dGamma01 = dDm4 - dDp4;
			dGamma02 = dDmm4 - dDpp4;
			dGamma11 = dDm4*dDp + dDp4*dDm;
			dGamma21 = 0.5*(dDp2*dDm4-dDm2*dDp4);
			dGamma12 = dDmm4*dDpp+dDmm*dDpp4;
			dGamma22 = 0.5*(dDpp2*dDmm4-dDmm2*dDpp4);
			dGamma32 = (1 / 6.0)*(dDpp4*dDmm4/dDpp+dDmm4*dDpp4/dDmm);
			dGamma31 = (1 / 6.0)*(dDp4*dDm4/dDp+dDm4*dDp4/dDm);

			double M0 = dGamma31*dGamma02 - dGamma32*dGamma01;
			double M1 = dGamma31*dGamma12 - dGamma32*dGamma11;
			double M2 = dGamma31*dGamma22 - dGamma32*dGamma21;
 
			double alpham = -dDp / (dDm*(dDm + dDp));
			double alpha0 = (dDp-dDm) / (dDm*dDp);
			double alphap = dDm / (dDp*(dDm + dDp));


			outVec[0] = -dGamma31*dDpp4 / M1;
			outVec[1] = (dDp4*dGamma32 - M2*alpham) / M1;
			outVec[2] = (-M0-alpha0*M2)/M1;
			outVec[3] = (-dGamma32*dDm4 - alphap*M2) / M1;
			outVec[4] = dGamma31*dDmm4/M1;
		}
	}
}

void PDEScheme::stencilFirstOrderCoefficients(const std::vector<double> & dx, bool IsUniform,
	PointDiscretization type, int index, double &outVec)
{

	if (type == ThreePoints){
		switch (index){
		case 0: 
					outVec = -(dx[1] / (dx[0] * (dx[0] + dx[1])));
					break;
		
		case 1: 
					outVec = (dx[1] - dx[0]) / (dx[0] * dx[1]);
					break;
		
		case 2: 
					outVec = (dx[0] / (dx[1] * (dx[0] + dx[1])));
					break;
		}
	}
	else
	{
		if (IsUniform){
			switch (index){
			case 0: 
				outVec = (1 / (12.0*dx[1]));
				break;
			case 1: 
						outVec = -2 / (dx[1] * 3.0);
						break;
			case 2: 
						outVec = 0.0;
						break;
			case 3: 
						outVec = 2 / (dx[1] * 3.0);
						break;
			case 4: 
						outVec = -(1 / (12.0*dx[1]));
						break;
			}
	
		}
		else{
			/*double dDpp, dDp, dDm, dDmm;
			double dDpp2, dDp2, dDm2, dDmm2;
			double dDpp4, dDp4, dDm4, dDmm4;
			double dGamma01, dGamma02;
			double dGamma11, dGamma12, dGamma32;
			double dGamma21, dGamma22, dGamma31;

			dDpp = dx[3];
			dDp = dx[2];
			dDm = dx[1];
			dDmm = dx[0];

			dDpp2 = dDpp*dDpp;
			dDp2 = dDp*dDp;
			dDm2 = dDm*dDm;
			dDmm2 = dDmm*dDmm;

			dDpp4 = dDpp2*dDpp2;
			dDp4 = dDp2*dDp2;
			dDm4 = dDm2*dDm2;
			dDmm4 = dDmm2*dDmm2;

			dGamma01 = dDm4 - dDp4;
			dGamma02 = dDmm4 - dDpp4;
			dGamma11 = dDm4*dDp + dDp4*dDm;
			dGamma21 = 0.5*(dDp2*dDm4 - dDm2*dDp4);
			dGamma12 = dDmm4*dDpp + dDmm*dDpp4;
			dGamma22 = 0.5*(dDpp2*dDmm4 - dDmm2*dDpp4);
			dGamma32 = (1 / 6.0)*(dDpp4*dDmm4 / dDpp + dDmm4*dDpp4 / dDmm);
			dGamma31 = (1 / 6.0)*(dDp4*dDm4 / dDp + dDm4*dDp4 / dDm);

			double M0 = dGamma31*dGamma02 - dGamma32*dGamma01;
			double M1 = dGamma31*dGamma12 - dGamma32*dGamma11;
			double M2 = dGamma31*dGamma22 - dGamma32*dGamma21;

			double alpham = -dDp / (dDm*(dDm + dDp));
			double alpha0 = (dDp - dDm) / (dDm*dDp);
			double alphap = dDm / (dDp*(dDm + dDp));


			outVec[0] = -dGamma31*dDpp4 / M1;
			outVec[1] = (dDp4*dGamma32 - M2*alpham) / M1;
			outVec[2] = (-M0 - alpha0*M2) / M1;
			outVec[3] = (-dGamma32*dDm4 - alphap*M2) / M1;
			outVec[4] = dGamma31*dDmm4 / M1;*/
		}
	}
}

void PDEScheme::stencilSecondOrderCoefficients(const std::vector<double> & dx, bool IsUniform,
	PointDiscretization type, std::vector<double> &outVec){

	if (type == ThreePoints){
		outVec[0] = (2 / (dx[0] * (dx[0] + dx[1])));
		outVec[1] = outVec[0] * dx[0] * (-1 / dx[0] - 1 / dx[1]);
		outVec[2] = outVec[0] * dx[0] / dx[1];
	}
	else{
		if (IsUniform){
			outVec[0] = -1 / (12 * dx[0] * dx[0]);
			outVec[1] = -16.0*outVec[0];
			outVec[2] = 30 * outVec[0];
			outVec[3] = -outVec[1];
			outVec[4] = outVec[0];
		}
	}
}

double PDEScheme::differentiator(const std::vector<double> & u,
	const std::vector<double> & dx, int index,
	PointDiscretization type, bool IsFirstOrSecond)
{
	double dDiff;
	int Npoints;
	(type == ThreePoints ? Npoints = 3 : Npoints = 5);

	dVector dStencil(Npoints);
	if (IsFirstOrSecond){
		stencilFirstOrderCoefficients(dx,mGrid->getIsUniformGrid(index), type,dStencil);
	}
	else{
		stencilSecondOrderCoefficients(dx, mGrid->getIsUniformGrid(index), type, dStencil);
	}
		
	dDiff = Matrix::dotProduct(dStencil, u);
	return dDiff;
}
