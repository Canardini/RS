#ifndef PDE_PC_CHEYETTE_H
#define PDE_PC_CHEYETTE_H


#pragma once
#include <Math/PDEModel/PDEmodel.h>

class PDEPC1FCheyette :public PDEmodel
{
public:
	PDEPC1FCheyette(double inKappa, double inf0t, const VolatilityModelBridge & inVolatilityModel);

	boost::shared_ptr<PDEmodel> clone() const;
	virtual void updateConvections(const std::vector<double> & inFactors);
	virtual void updateConvections(const std::vector<double> & inFactors, int index);
	virtual void updateDiffusions(const std::vector<double> & inFactors);
	virtual void updateDiffusions(const std::vector<double> & inFactors, int index);
	virtual void updateMixedDerivatives(const std::vector<double> & inFactors);
	virtual void updateSourceTerms(const std::vector<double> & inFactors);
	virtual void updateSourceTerms(const std::vector<double> & inFactors, int index);
	void updateCheyetteVol(const std::vector<double> & inFactors);
	virtual ~PDEPC1FCheyette(){};

private:
	double mKappa;


};

#endif // !PDE_CHEYETTE_H
