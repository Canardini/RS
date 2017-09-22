#include <vector>
#include <algorithm>

#include <boost/function.hpp>
#include <boost/thread.hpp>
#include <Math/Solver/DifferentialEvolution.h>


typedef std::vector<double> DoubleRow;
typedef std::vector<DoubleRow> DoubleMatrix;

class DifferentialEvolutionEngine {
	DoubleRow  mInitialObjectiveFunction;
	DoubleRow  mCrossOverObjectiveFunction;
	boost::function<double(const std::vector<double> &)> mFunc;
	DoubleMatrix  mInitialPopulation;
	DoubleMatrix  mCrossOverPopulation;
	int mNumberOfThreads;
	int mNP;
	int  mMaxIter;

	std::vector<UniformRandom<>> 	mInitialUniformRandom;
	std::vector<UniformIntRandom<>>  mTempRandom;
	std::vector<UniformIntRandom<>>  mRamdomArgsize;
	int mSeed;
	std::vector<double> mArgument;
	double mRet;
	const std::vector<double> mLowerBound;
	const std::vector<double> mUpperBound;
	double mCrossOverRatio;
	double mScaleFactor;
	bool mUseSymmetryShortening;
	boost::function<void(std::vector<double> &)> mSymmetryShortening;
	bool mOutputLog;

public:
	class LocalFunction{
		DifferentialEvolutionEngine* mObj;
		int mIndex;
		boost::barrier& mBarrier;
	public:
		LocalFunction(DifferentialEvolutionEngine* inObj, int inIndex, boost::barrier& inBarrier) : mObj(inObj), mIndex(inIndex), mBarrier(inBarrier) {}
		void operator()(){
			mObj->DEProcess(mIndex, mBarrier);
		}
	};

	DifferentialEvolutionEngine(DoubleRow inInitialObjectiveFunction,
		boost::function<double(const std::vector<double> &)> inFunc,
		DoubleMatrix inInitialPopulation,
		int inNumberOfThreads,
		int inNP,
		int inMaxiter,
		int inSeed,
		std::vector<double> &inArgument,
		double &outRet,
		const std::vector<double> &inLowerBound,
		const std::vector<double> &inUpperBound,
		double inCrossOverRatio,
		double inScaleFactor,
		bool inUseSymmetryShortening,
		boost::function<void(std::vector<double> &)> inSymmetryShortening,
		bool inOutputLog
		) : mInitialObjectiveFunction(inInitialObjectiveFunction), mCrossOverObjectiveFunction(inInitialObjectiveFunction), mFunc(inFunc), mInitialPopulation(inInitialPopulation), mCrossOverPopulation(inInitialPopulation), mNumberOfThreads(inNumberOfThreads), mNP(inNP), mMaxIter(inMaxiter), mSeed(inSeed),
		mArgument(inArgument), mRet(outRet), mLowerBound(inLowerBound), mUpperBound(inUpperBound), mCrossOverRatio(inCrossOverRatio), mScaleFactor(inScaleFactor), mUseSymmetryShortening(inUseSymmetryShortening), mSymmetryShortening(inSymmetryShortening), mOutputLog(inOutputLog)
	{
		for (int i = 0; i < mNumberOfThreads + 1; ++i)
		{
			if (i < mNumberOfThreads)
			{
				if (i == 0)
				{
					mTempRandom.push_back(UniformIntRandom<>(0, mNP - 1, mSeed));
				}
				else
				{
					mTempRandom.push_back(UniformIntRandom<>(0, mNP - 1, mTempRandom[0]()));
				}
				mRamdomArgsize.push_back(UniformIntRandom<>(0, static_cast<int>(inArgument.size() - 1), mTempRandom[0]()));
			}

			mInitialUniformRandom.push_back(UniformRandom<>(mTempRandom[0]()));
		}
		resetPopulation();
	}

	~DifferentialEvolutionEngine(){}
	void resetPopulation()
	{
		for (int i = 0; i < mNP; i++)
		{
			for (size_t j = 0; j < mArgument.size(); j++)
			{
				mInitialPopulation[i][j] = mLowerBound[j] + (mUpperBound[j] - mLowerBound[j])*mInitialUniformRandom[0]();
			}
			if (mUseSymmetryShortening)
			{
				mSymmetryShortening(mInitialPopulation[i]);
			}
		}

	}

	double getRet()
	{
		return mRet;
	};
	std::vector<double> getArg()
	{
		return mArgument;
	};

	void partialPermutation2(UniformIntRandom<> & UniformIntRandom_,
		int & lBasenumber,
		int & lDiffnumber1,
		int & lDiffnumber2,
		int lMutatednumber)
	{

		lBasenumber = UniformIntRandom_();
		while (lBasenumber == lMutatednumber){
			lBasenumber = UniformIntRandom_();
		}

		lDiffnumber1 = UniformIntRandom_();
		while ((lDiffnumber1 == lMutatednumber) || (lDiffnumber1 == lBasenumber)){
			lDiffnumber1 = UniformIntRandom_();
		}

		lDiffnumber2 = UniformIntRandom_();
		while ((lDiffnumber2 == lMutatednumber) || (lDiffnumber2 == lBasenumber) || (lDiffnumber2 == lDiffnumber1)){
			lDiffnumber2 = UniformIntRandom_();
		}

	}

	void initializationObj(int inThreadID){
		for (int i = 0; i < mNP; i++)
		{
			if (i % mNumberOfThreads == inThreadID)
			{
				mCrossOverObjectiveFunction[i] = mInitialObjectiveFunction[i] = mFunc(mInitialPopulation[i]);
			}
		}
	}

	void mutationCrossOverProcess(int inThreadID, boost::barrier& inBarrier)
	{
		for (int k = 0; k < mMaxIter; k++)
		{
			if ((k + 1) % 100 == 0 && mOutputLog)
			{
				std::vector<double>::iterator low = std::min_element(mInitialObjectiveFunction.begin(), mInitialObjectiveFunction.end());
				int lMinElement = static_cast<int>(low - mInitialObjectiveFunction.begin());
				std::stringstream lSS;
				lSS << std::setprecision(20);
				lSS << "thread number : " << inThreadID << std::endl;
				lSS << "reached  " << k + 1 << " simulations.  " << std::endl;
				lSS << "minimum value :  " << mInitialObjectiveFunction[lMinElement] << std::endl;
				for (size_t i = 0; i < mArgument.size(); i++)
				{
					lSS << "X" << "[" << i << "]= " << mInitialPopulation[lMinElement][i] << std::endl;
				}
				lSS << std::endl;
				//FELUtilityNS::logDebug("FELMathOptimizationNS::legacy::DifferentialEvolution2", lSS.str());
			}

			//Mutation & CrossOver
			for (int i = 0; i < mNP; i++)
			{
				if (i % mNumberOfThreads == inThreadID)
				{
					int lMutatednumber = i;
					int lBasenumber = 0;
					int lDiffnumber1 = 0;
					int lDiffnumber2 = 0;
					partialPermutation2(mTempRandom[inThreadID], lBasenumber, lDiffnumber1, lDiffnumber2, lMutatednumber);
					int r = mRamdomArgsize[inThreadID]();
					for (size_t j = 0; j < mArgument.size(); j++)
					{
						double u = mInitialUniformRandom[inThreadID + 1]();
						if ((j == r) || u <= mCrossOverRatio)
						{
							//lCrossOverPopulation[j] = std::max(mLowerBound[j], std::min(mUpperBound[j], mInitialPopulation[lBasenumber][j] + mScaleFactor*(mInitialPopulation[lDiffnumber1][j] - mInitialPopulation[lDiffnumber2][j])));
							mCrossOverPopulation[i][j] = mInitialPopulation[lBasenumber][j] + mScaleFactor*(mInitialPopulation[lDiffnumber1][j] - mInitialPopulation[lDiffnumber2][j]);

							if ((mCrossOverPopulation[i][j] > mUpperBound[j]) || (mCrossOverPopulation[i][j] < mLowerBound[j]))
							{
								mCrossOverPopulation[i][j] = mLowerBound[j] + (mUpperBound[j] - mLowerBound[j])*u;
							}
						}
						else
						{
							mCrossOverPopulation[i][j] = mInitialPopulation[i][j];
						}
					}
					if (mUseSymmetryShortening)
					{
						mSymmetryShortening(mCrossOverPopulation[i]);
					}
					double temp = mFunc(mCrossOverPopulation[i]);
					//selection
					if (mInitialObjectiveFunction[i] > temp)
					{
						mCrossOverObjectiveFunction[i] = temp;
					}
					else{
						mCrossOverPopulation[i] = mInitialPopulation[i];
					}
				}
			}

			inBarrier.wait();
			for (int i = 0; i < mNP; i++){
				if (i % mNumberOfThreads == inThreadID)
				{
					mInitialPopulation[i] = mCrossOverPopulation[i];
					mInitialObjectiveFunction[i] = mCrossOverObjectiveFunction[i];
				}
			}
			inBarrier.wait();
		}
	}

	void DEProcess(int inThreadID, boost::barrier& inBarrier){
		initializationObj(inThreadID);
		mutationCrossOverProcess(inThreadID, inBarrier);
	}

	void multiThread()
	{
		boost::thread_group lThreads;
		boost::barrier lBarrier(mNumberOfThreads);
		for (int i = 0; i < mNumberOfThreads; ++i)
		{
			lThreads.create_thread(LocalFunction(this, i, lBarrier));
		}

		lThreads.join_all();

		std::vector<double>::iterator low = std::min_element(mInitialObjectiveFunction.begin(), mInitialObjectiveFunction.end());
		int lMinElement = static_cast<int>(low - mInitialObjectiveFunction.begin());
		mArgument = mInitialPopulation[lMinElement];
		mRet = mInitialObjectiveFunction[lMinElement];
	}

};
	void DifferentialEvolution2::nullFunc(std::vector<double>&){}

	int DifferentialEvolution2::Optimize
		(boost::function<double(const std::vector<double> &)> inFunc,
		std::vector<double> &inArgument,
		double &outRet,
		const std::vector<double> &inLowerBound,
		const std::vector<double> &inUpperBound)
	{
		int lNumberOfThreads = mNumOfThreads;

		DoubleMatrix lInitialPopulation(mNP, DoubleRow(inArgument.size()));
		DoubleRow lInitialObjectiveFunction(mNP);

		DifferentialEvolutionEngine DEEngine(lInitialObjectiveFunction, inFunc, lInitialPopulation, lNumberOfThreads, mNP, mMaxIter, mSeed,
			inArgument, outRet, inLowerBound, inUpperBound,
			mCrossOverRatio, mScaleFactor, mUseSymmetryShortening, mSymmetryShortening, mOutputLog);

		DEEngine.multiThread();
		inArgument = DEEngine.getArg();
		outRet = DEEngine.getRet();
		return mMaxIter;
	}


