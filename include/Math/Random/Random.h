#ifndef RANDOM_RANDOM_H
#define RANDOM_RANDOM_H
/**
*  @file Random.h boostÇÃóêêîê∂ê¨ÇÃWrapper
*/
#pragma once

#include <random>

template<class Dist_, class Gen_ = std::mt19937>
class Random {
	Gen_ mGen;
	Dist_ mDst;
	int mIndex;
	int mMax;
public:
	Random(int inSeed = 1) : mGen(inSeed), mIndex(0), mMax(1) {}
	template<typename T1_>
	Random(int inSeed, T1_ inDstArg1) : mGen(inSeed), mDst(inDstArg1), mIndex(0), mMax(1) {}
	template<typename T1_, typename T2_>
	Random(int inSeed, T1_ inDstArg1, T2_ inDstArg2) : mGen(inSeed), mDst(inDstArg1, inDstArg2), mIndex(0), mMax(1) {}

	void reset(int inSeed = 1){
		mGen = Gen_(inSeed);
		mDst(mGen);
	}

	void split(int inMax, int inIndex){
		mIndex = std::max(inIndex, 0);
		mMax = std::max(inMax, 1);
	}

	typename Dist_::result_type operator()() {
		double lVal = 0.0;
		for (int i = 0; i<mMax; i++){
			if (i == mIndex){
				lVal = mDst(mGen);
			}
		}
		return lVal;
	}
};


template<class Gen_ = std::mt19937>
class NormalRandom {
	typedef std::normal_distribution<> Dist_;
	Gen_ mGen;
	Dist_ mDst;

	int mIndex;
	int mMax;

	double getRand(){
		double lVal = 0.0;
		for (int i = 0; i<mMax; i++){
			if (i == mIndex){
				lVal = mDst(mGen);
			}
		}
		return lVal;
	}

public:
	NormalRandom(int inSeed = 1) : mGen(inSeed), mDst(0, 1), mIndex(0), mMax(1) {}
	NormalRandom(double inMean, double inStdV, int inSeed = 1) : mGen(inSeed), mDst(inMean, inStdV), mIndex(0), mMax(1) {}

	void reset(int inSeed = 1){
		mGen = Gen_(inSeed);
		mDst(mGen);
	}

	void split(int inMax, int inIndex){
		mIndex = std::max(inIndex, 0);
		mMax = std::max(inMax, 1);
	}

	typename Dist_::result_type operator()() { return getRand(); }

	typename Dist_::result_type getRandomVariate() { return getRand(); }

	typename Dist_::result_type getRandomVariateInt() { return getRand(); }

	template<typename T_>
	void getNRandomVariates(T_ inIterStart, T_ inIterEnd){
		for (; inIterStart != inIterEnd; ++inIterStart){
			*inIterStart = getRand();
		}
	}

	template<typename T_>
	void getNRandomVariates(std::vector<T_>& outRandVector)
	{
		for (size_t i = 0; i<outRandVector.size(); i++){
			outRandVector[i] = getRand();
		}
	}

	template<typename T_>
	void operator()(int inNum, T_* outRandVector) {
		for (int i = 0; i<inNum; i++){
			outRandVector[i] = getRand();
		}
	}

	template<typename T_>
	void operator()(std::vector<T_>& outRandVector) {
		for (size_t i = 0; i<outRandVector.size(); i++){
			outRandVector[i] = getRand();
		}
	}


	template<typename T_>
	void operator()(std::vector<std::vector<T_>>& outRandMatrix) {
		for (std::size_t i = 0; i<outRandMatrix.size(); i++){
			for (size_t j = 0; j<outRandMatrix[i].size(); j++){
				outRandMatrix[i][j] = getRand();
			}
		}
	}

	void jump(int inStep){
		// this is a temporary implementation
		int l = 1;
		l << inStep;
		for (int i = 0; i < l; i++){
			getRand();
		}
	}
};


template<class Gen_ = std::mt19937>
class UniformRandom {

#if defined(__GNUC__) && (__GNUC__ < 4 || (__GNUC__ ==4 && __GNUC_MINOR__ < 5))
	typedef std::uniform_real<> Dist_;
#else
	typedef std::uniform_real_distribution<> Dist_;
#endif
	Gen_ mGen;
	Dist_ mDst;

	int mIndex;
	int mMax;

	double getRand(){
		double lVal = 0.0;
		for (int i = 0; i<mMax; i++){
			if (i == mIndex){
				lVal = mDst(mGen);
			}
		}
		return lVal;
	}

public:
	UniformRandom() : mGen(), mDst(), mIndex(0), mMax(1)  {}
	UniformRandom(int inSeed) : mGen(inSeed), mDst(), mIndex(0), mMax(1)  {}
	UniformRandom(double inStart, double inEnd) : mGen(), mDst(inStart, inEnd), mIndex(0), mMax(1)  {}
	UniformRandom(double inStart, double inEnd, int inSeed) : mGen(inSeed), mDst(inStart, inEnd), mIndex(0), mMax(1)  {}

	void reset(int inSeed = 1){
		mGen = Gen_(inSeed);
		mDst(mGen);
	}

	void split(int inMax, int inIndex){
		mIndex = std::max(inIndex, 0);
		mMax = std::max(inMax, 1);
	}

	typename Dist_::result_type operator()() { return getRand(); }

	typename Dist_::result_type getRandomVariate() { return getRand(); }

	template<typename T_>
	void getNRandomVariates(std::vector<T_>& outRandVector)
	{
		for (size_t i = 0; i<outRandVector.size(); i++){
			outRandVector[i] = getRand();
		}
	}

	template<typename T_>
	void getNRandomVariates(T_ inIterStart, T_ inIterEnd){
		for (; inIterStart != inIterEnd; ++inIterStart){
			*inIterStart = getRand();
		}
	}

	template<typename T_>
	void operator()(std::vector<T_>& outRandVector) {
		for (size_t i = 0; i<outRandVector.size(); i++){
			outRandVector[i] = getRand();
		}
	}

	template<typename T_>
	void operator()(std::vector<std::vector<T_>>& outRandMatrix) {
		for (size_t i = 0; i<outRandMatrix.size(); i++){
			for (size_t j = 0; j<outRandMatrix[i].size(); j++){
				outRandMatrix[i][j] = getRand();
			}
		}
	}

	void jump(int inStep){
		// this is a temporary implementation
		int l = 1;
		l = l << inStep;
		for (int i = 0; i < l; i++){
			getRand();
		}
	}
};

template<class Gen_ = std::mt19937>
class UniformIntRandom {
#if defined(__GNUC__) && (__GNUC__ < 4 || (__GNUC__ ==4 && __GNUC_MINOR__ < 5))
	typedef std::uniform_int<> Dist_;
#else
	typedef std::uniform_int_distribution<> Dist_;
#endif
	Gen_ mGen;
	Dist_ mDst;

	int mIndex;
	int mMax;

	Dist_::result_type getRand(){
		Dist_::result_type lVal = static_cast<Dist_::result_type>(0);
		for (int i = 0; i<mMax; i++){
			if (i == mIndex){
				lVal = mDst(mGen);
			}
		}
		return lVal;
	}

public:
	UniformIntRandom() : mGen(), mDst(0, 1), mIndex(0), mMax(1)  {}
	UniformIntRandom(int inStart, int inEnd) : mGen(), mDst(inStart, inEnd), mIndex(0), mMax(1)  {}
	UniformIntRandom(int inStart, int inEnd, int inSeed) : mGen(inSeed), mDst(inStart, inEnd), mIndex(0), mMax(1)  {}

	void reset(int inSeed = 1){
		mGen = Gen_(inSeed);
		mDst(mGen);
	}

	void split(int inMax, int inIndex){
		mIndex = std::max(inIndex, 0);
		mMax = std::max(inMax, 1);
	}

	typename Dist_::result_type operator()() { return getRand(); }

	typename Dist_::result_type getRandomVariateInt() { return getRand(); }

	template<typename T_>
	void getNRandomVariates(T_ inIterStart, T_ inIterEnd){
		for (; inIterStart != inIterEnd; ++inIterStart){
			*inIterStart = getRand();
		}
	}

	template<typename T_>
	void getNRandomVariates(std::vector<T_>& outRandVector)
	{
		for (size_t i = 0; i<outRandVector.size(); i++){
			outRandVector[i] = getRand();
		}
	}

	template<typename T_>
	void operator()(std::vector<T_>& outRandVector) {
		for (size_t i = 0; i<outRandVector.size(); i++){
			outRandVector[i] = getRand();
		}
	}

	template<typename T_>
	void operator()(std::vector<std::vector<T_>>& outRandMatrix) {
		for (size_t i = 0; i<outRandMatrix.size(); i++){
			for (size_t j = 0; j<outRandMatrix[i].size(); j++){
				outRandMatrix[i][j] = getRand();
			}
		}
	}
	void jump(int inStep){
		// this is a temporary implementation
		int l = 1;
		l << inStep;
		for (int i = 0; i < l; i++){
			getRand();
		}
	}
};

#endif /* RANDOM_RANDOM_H */
