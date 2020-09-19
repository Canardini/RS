#include<Product\MarkovianProductBridge.h>




MarkovianProductBridge::MarkovianProductBridge(const MarkovianProductBridge& inOriginal){
	mMarkovianProductPtr = inOriginal.mMarkovianProductPtr->clone();
}

MarkovianProductBridge::MarkovianProductBridge(const MarkovianProduct& inInnerIRModel){
	mMarkovianProductPtr = inInnerIRModel.clone();
}

MarkovianProductBridge& MarkovianProductBridge::operator=(const MarkovianProductBridge& inOriginal){

	if (this != &inOriginal){
		mMarkovianProductPtr = inOriginal.mMarkovianProductPtr->clone();
	}
	return *this;
}

void MarkovianProductBridge::updateMarkovianFactors(std::vector<double> & inFactors){
	return mMarkovianProductPtr->updateMarkovianFactors(inFactors);
}

double MarkovianProductBridge::getPrice(){
	return mMarkovianProductPtr->getPrice();
}


