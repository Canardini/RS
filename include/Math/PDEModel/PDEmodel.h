#ifndef PDE_MODEL_H
#define PDE_MODEL_H
#pragma once

#include <vector>
#include <Model/VolatilityModel/VolatilityModelBridge.h>
#include <boost/make_shared.hpp>

class PDEmodel
{
public:
	
	PDEmodel(const VolatilityModelBridge & inVolatilityModel) :mVolatilityModel(inVolatilityModel){};
	virtual void updateConvections(const std::vector<double> & inFactors)=0;
	virtual void updateDiffusions(const std::vector<double> & inFactors)=0 ;
	virtual void updateSourceTerms(const std::vector<double> & inFactors)=0 ;

	virtual void updateConvections(const std::vector<double> & inFactors, int index) = 0;
	virtual void updateDiffusions(const std::vector<double> & inFactors, int index) = 0;
	virtual void updateSourceTerms(const std::vector<double> & inFactors, int index) = 0;

	virtual boost::shared_ptr<PDEmodel> clone() const = 0;

	double getConvection(int index){ return mConvections[index]; }
	double getDiffusion(int index){ return mDiffusions[index];}
	double getSourceTerms(int index){ return mSourceTerms[index]; }
	double getMixedDerivatives(int index1,int index2 ){ return mMixedDerivatives[index1][index2]; }
	void setFactors(const std::vector<double> & inFactors){ mFactors = inFactors; };
	void setFactor(double inFactor,int index){ mFactors[index] = inFactor; };
	void setInitSource(double inSource){ mInitSource = inSource; };
	void setDriftAdjust(double inDrift){ mDriftAdjust = inDrift; };

	void setEntity(double inEntity){ mVolatilityModel.setEntity(inEntity);};
	void setIndex(int inIndex){ mVolatilityModel.setIndex(inIndex); };
	virtual ~PDEmodel(){};

protected:
	int mNumberOfFactors;
	std::vector<double> mFactors;
    std::vector<double> mConvections;
	std::vector<double> mDiffusions;
    std::vector<std::vector<double>> mMixedDerivatives;
	std::vector<double> mSourceTerms;
	double mInitSource;
	double mVol;
	double mDriftAdjust;
	VolatilityModelBridge mVolatilityModel;

};

#endif // !PDE_MODEL_H
