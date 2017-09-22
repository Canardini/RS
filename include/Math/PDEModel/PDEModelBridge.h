#ifndef PDE_MODEL_BRIDGE_H
#define PDE_MODEL_BRIDGE_H


#pragma once

#include <Math\PDEModel\PDEmodel.h>
#include <boost/make_shared.hpp>

class PDEModelBridge{

public:

	PDEModelBridge(const PDEModelBridge & original);
	PDEModelBridge(const PDEmodel & innerPDEModel);
	PDEModelBridge& operator=(const PDEModelBridge& original);

	void updateConvections(const std::vector<double> & inFactors){ mPDEModelPtr->updateConvections(inFactors); };
	void updateDiffusions(const std::vector<double> & inFactors){ mPDEModelPtr->updateDiffusions(inFactors); };
	void updateSourceTerms(const std::vector<double> & inFactors){ mPDEModelPtr->updateSourceTerms(inFactors); };

	void updateConvections(const std::vector<double> & inFactors, int index){ mPDEModelPtr->updateConvections(inFactors, index); };
	void updateDiffusions(const std::vector<double> & inFactors, int index){ mPDEModelPtr->updateDiffusions(inFactors, index); };
	void updateSourceTerms(const std::vector<double> & inFactors, int index){ mPDEModelPtr->updateSourceTerms(inFactors, index); };

	double getConvection(int index){ return  mPDEModelPtr->getConvection(index); };
	double getDiffusion(int index){ return mPDEModelPtr->getDiffusion(index); };
	double getSourceTerms(int index){ return mPDEModelPtr->getSourceTerms(index); };
	double getMixedDerivatives(int index1, int index2){ return mPDEModelPtr->getMixedDerivatives(index1, index2);};
	void setFactors(const std::vector<double> & inFactors){ mPDEModelPtr->setFactors(inFactors); };
	void setFactor(double inFactor, int index){ mPDEModelPtr->setFactor(inFactor,index); };
	void setInitSource(double inSource){ mPDEModelPtr->setInitSource(inSource); };
	void setEntity(double inEntity){ mPDEModelPtr->setEntity(inEntity); };
	void setIndex(int inIndex){ mPDEModelPtr->setIndex(inIndex); };
	void setDriftAdjust(double inDrift){ mPDEModelPtr->setDriftAdjust(inDrift); };
	PDEModelBridge(){};

private:
	boost::shared_ptr<PDEmodel> mPDEModelPtr;
};

#endif