#include<Model\IRModel\MarkovianIRModelBridge.h>



MarkovianIRModelBridge::MarkovianIRModelBridge(const MarkovianIRModelBridge& inOriginal){
	mMarkovianIRModelPtr = inOriginal.mMarkovianIRModelPtr->clone();
}

MarkovianIRModelBridge::MarkovianIRModelBridge(const MarkovianIRModel& inInnerIRModel){
	mMarkovianIRModelPtr = inInnerIRModel.clone();
}

MarkovianIRModelBridge& MarkovianIRModelBridge::operator=(const MarkovianIRModelBridge& inOriginal){
	if (this != &inOriginal){
		mMarkovianIRModelPtr = inOriginal.mMarkovianIRModelPtr->clone();
	}
	return *this;
}



double MarkovianIRModelBridge::getGtT(const bgreg::date & inValueDate, const bgreg::date & inMaturity){
	return mMarkovianIRModelPtr->getGtT(inValueDate, inMaturity);
}

double MarkovianIRModelBridge::getInitDiscountValue(double inMaturity){
	return mMarkovianIRModelPtr->getInitDiscountValue(inMaturity);
}

double MarkovianIRModelBridge::getInitDiscountValueByDate(bgreg::date & inMaturity){
	return mMarkovianIRModelPtr->getInitDiscountValueByDate(inMaturity);
}


double MarkovianIRModelBridge::getDiscountValue(double inValueDate, double inMaturity){
	return mMarkovianIRModelPtr->getDiscountValue(inValueDate,inMaturity);
}

double MarkovianIRModelBridge::getDiscountValueByDate(double inValueDate, bgreg::date & inMaturity){
	return mMarkovianIRModelPtr->getDiscountValueByDate(inValueDate, inMaturity);
}


double MarkovianIRModelBridge::getFowardValue(double inValueDate, double inT1, double inT2){
	return mMarkovianIRModelPtr->getFowardValue(inValueDate,inT1,inT2);
}

double MarkovianIRModelBridge::getFowardValueByDate(double inValueDate, bgreg::date & inT1){
	return mMarkovianIRModelPtr->getFowardValueByDate(inValueDate, inT1);
}

double MarkovianIRModelBridge::getInitFWValueByDate(bgreg::date & inStartDate){
	return mMarkovianIRModelPtr->getInitFWValueByDate(inStartDate);
}

void MarkovianIRModelBridge::updateMarkovianFactors(std::vector<double> & inFactors){
	return mMarkovianIRModelPtr->updateMarkovianFactors(inFactors);
}

double MarkovianIRModelBridge::updateForwardValue(double delta, double FW, double OISspread, std::vector<double> &dx, double G1, double G2)
{
	return mMarkovianIRModelPtr->updateForwardValue(delta, FW, OISspread, dx, G1, G2);
}