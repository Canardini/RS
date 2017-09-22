#ifndef MATH_OPTIMIZATION_LEGACY_DIFFERENTIALEVOLUTION2_H
#define MATH_OPTIMIZATION_LEGACY_DIFFERENTIALEVOLUTION2_H

#pragma once

#include <vector>
#include <boost/function.hpp>
#include<Math/Random/Random.h>

class DifferentialEvolution2
{
	int mMaxIter;
	int mNP;
	bool mUseSymmetryShortening;
	boost::function<void(std::vector<double> &)> mSymmetryShortening;
	double mCrossOverRatio;
	double mScaleFactor;
	int mSeed;
	int mNumOfThreads;
	bool mOutputLog;

public:
	DifferentialEvolution2(int inMaxIter = 100l, int inNp = 200, bool inUseSymmetryShortening = false, const boost::function<void(std::vector<double> &)>& inSymmetryShortening = &(DifferentialEvolution2::nullFunc), double inCrossOverRatio = 0.5, double inScaleFactor = 0.8, int inSeed = 1011, int inNumOfThreads =1, bool inOutputLog = false) :
		mMaxIter(inMaxIter), mNP(inNp), mUseSymmetryShortening(inUseSymmetryShortening), mSymmetryShortening(inSymmetryShortening), mCrossOverRatio(inCrossOverRatio), mScaleFactor(inScaleFactor), mSeed(inSeed), mNumOfThreads(inNumOfThreads), mOutputLog(inOutputLog){}
	void setNumberOfThreads(int inNum){ mNumOfThreads = inNum; };
	DifferentialEvolution2(const DifferentialEvolution2 & inObj) :
		mMaxIter(inObj.mMaxIter), mNP(inObj.mNP), mUseSymmetryShortening(inObj.mUseSymmetryShortening), mSymmetryShortening(inObj.mSymmetryShortening), mCrossOverRatio(inObj.mCrossOverRatio), mScaleFactor(inObj.mScaleFactor), mSeed(inObj.mSeed), mNumOfThreads(inObj.mNumOfThreads), mOutputLog(inObj.mOutputLog){}
	int Optimize
		(boost::function<double(const std::vector<double> &)> inFunc,
		std::vector<double> &inArgument,
		double &outRet,
		const std::vector<double> &inLowerBound,
		const std::vector<double> &inUpperBound);

	~DifferentialEvolution2(){}
	static void nullFunc(std::vector<double>&);
};

#endif 