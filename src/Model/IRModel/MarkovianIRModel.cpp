#include<Model\IRModel\MarkovianIRModel.h>


MarkovianIRModel::MarkovianIRModel(const DiscountCurve & inDiscountCurve,
	const DiscountCurve & inForwardCurve,
	std::vector<double> inMarkovianFactors)
{
	mDiscountCurve = (inDiscountCurve.clone());
	mForwardCurve = (inForwardCurve.clone());
	mMarkovianFactors = inMarkovianFactors;
	mNumberOfFactors=mMarkovianFactors.size();

}

void MarkovianIRModel::updateMarkovianFactors(std::vector<double> & inFactors){
	for (int i = 0; i < mNumberOfFactors; i++){
		mMarkovianFactors[i] = inFactors[i];
	}
}

double MarkovianIRModel::getInitDiscountValue(double inMaturity){

	return mDiscountCurve->getDFbyDouble(inMaturity);
}


double MarkovianIRModel::getInitDiscountValueByDate(bgreg::date & inMaturity){

	return mDiscountCurve->getDFbyDate(inMaturity);
}


double MarkovianIRModel::getInitFWValue(double inT1, double inT2){
	return mForwardCurve->getFWbyDouble(inT1, inT2);
}


double MarkovianIRModel::getInitFWValueByDate(bgreg::date & inMaturity){
	return mForwardCurve->getFWbyDate(inMaturity);
}

int MarkovianIRModel::locate(const std::vector<double> & xx, double inX0){
	int ju, jm, jl;
	int mm = 1;
	int n = xx.size();
	if (inX0 > xx[n - 1]){
		return n - 1;
	}
	if (inX0 < xx[0]){
		return -1;
	}

	bool ascnd = (xx[n - 1] >= xx[0]);
	jl = 0;
	ju = n - 1;
	while (ju - jl > 1) {
		jm = (ju + jl) >> 1;
		if (inX0 >= xx[jm] == ascnd)
			jl = jm;
		else
			ju = jm;
	}

	return jl;
}


