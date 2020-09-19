#include <Product\MarkovianProduct.h>


MarkovianProduct::MarkovianProduct(const MarkovianIRModelBridge & inMarkovianIRModel, int inNumberOfFactors, std::vector<double> &inMarkovianFactors, double inValueDate)
: mMarkovianIRModel(inMarkovianIRModel), mNumberOfFactors(inNumberOfFactors), mMarkovianFactors(inMarkovianFactors), mValueDate(inValueDate)
{

}

void MarkovianProduct ::updateMarkovianFactors(std::vector<double> & inFactors){
	for (int i = 0; i < mNumberOfFactors; i++){
		mMarkovianFactors[i] = inFactors[i];
	}
	mMarkovianIRModel.updateMarkovianFactors(inFactors);
}
MarkovianProduct::~MarkovianProduct()
{
}
