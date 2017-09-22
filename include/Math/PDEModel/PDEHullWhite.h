#ifndef PDE_HULL_WHITE_H
#define PDE_HULL_WHITE_H

#pragma once
#include <Math/PDEModel/PDEmodel.h>

class PDEHullWhite :public PDEmodel
{
public:
	PDEHullWhite(double inKappa, double inf0t, const VolatilityModelBridge & inVolatilityModel);

	boost::shared_ptr<PDEmodel> clone() const;
	virtual void updateConvections(const std::vector<double> & inFactors);
	virtual void updateConvections(const std::vector<double> & inFactors, int index);
	virtual void updateDiffusions(const std::vector<double> & inFactors);
	virtual void updateDiffusions(const std::vector<double> & inFactors, int index);
	virtual void updateMixedDerivatives(const std::vector<double> & inFactors);
	virtual void updateSourceTerms(const std::vector<double> & inFactors);
	virtual void updateSourceTerms(const std::vector<double> & inFactors, int index);
	void updateCheyetteVol(const std::vector<double> & inFactors);
	virtual ~PDEHullWhite(){};

private:
	double mKappa;

};

#endif // !PDE_CHEYETTE_H
