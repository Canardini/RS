#include <Math/PDE/Grid.h>
 
Grid::Grid(int inGridDimension,
	std::vector<double> &inLowerBounds,
	std::vector<double> &inUpperBounds,
	dVector &inCenterValues,
	std::vector<int> &inGridSizes):
	mGridDimension(inGridDimension),
	mLowerBounds(inLowerBounds),
	mUpperBounds(inUpperBounds),
	mCenterValues(inCenterValues),
	mGridSizes(inGridSizes)
{
	mIsUniformGrid.resize(mGridDimension);
}

Grid::Grid(int inGridDimension,
	std::vector<double> &inLowerBounds,
	std::vector<double> &inUpperBounds,
	dVector &inCenterValues,
	std::vector<int> &inGridSizes,
	double inAlpha) :
	mGridDimension(inGridDimension),
	mLowerBounds(inLowerBounds),
	mUpperBounds(inUpperBounds),
	mCenterValues(inCenterValues),
	mGridSizes(inGridSizes),
	mAlpha(inAlpha)
{
	mIsUniformGrid.resize(mGridDimension);
}

int Grid::getGridSize(int index){
	return mGridSizes[index];
}



int Grid::locate(int i, double inX0){
	int ju, jm, jl;
	int mm = 1;
	int n = mGridSizes[i];
	if (inX0 > mGrid[i][n - 1]){
		return n - 1;
	}
	if (inX0 < mGrid[i][0]){
		return -1;
	}

	bool ascnd = (mGrid[i][n - 1] >= mGrid[i][0]);
	jl = 0;
	ju = n - 1;
	while (ju - jl > 1) {
		jm = (ju + jl) >> 1;
		if (inX0 >= mGrid[i][jm] == ascnd)
			jl = jm;
		else
			ju = jm;
	}

	return jl;
}


void Grid::buildArcSinus(int inIndex){
	double lC1 = asinh((mLowerBounds[inIndex] - mCenterValues[inIndex]) / mAlpha);
	double lC2 = asinh((mUpperBounds[inIndex] - mCenterValues[inIndex]) / mAlpha);
	int lN = mGridSizes[inIndex];
	std::vector<double> oUniformGrid(lN);
	double dt = 1.0 / (mGridSizes[inIndex] - 1);
	double lVal;
	for (int j = 0; j < lN; j++){
		oUniformGrid[j]= dt*j;
	}
	mGrid.push_back(dVector(0));
	for (int j = 0; j < lN; j++){
		lVal = mCenterValues[inIndex] + mAlpha*sinh(lC2*oUniformGrid[j] + lC1*(1 - oUniformGrid[j]));
		if (lVal == 0.0){
			mZeroIndexes.push_back(j);
		}
		mGrid[0].push_back(lVal);
	}
	mIsUniformGrid[inIndex] = false;
}
   
void Grid::buildUniformGrid(void){
	double dt;
	for (int i = 0; i < mGridDimension; i++){
		mGrid.push_back(dVector(0));
		mIsUniformGrid[i] = true;
		dt = (mUpperBounds[i] - mLowerBounds[i]) / (mGridSizes[i]-1);
		for (int j = 0; j < mGridSizes[i]; j++){
				mGrid[i].push_back(mLowerBounds[i]+dt*j);
		}
	}
}

//void Grid::buildCenteredUniformGrid(void){
//	double dt;
//	double dCentreVal;
//	int halfIndex;
//	for (int i = 0; i < mGridDimension; i++){
//		if (mGridSizes[i] % 2 == 1){
//			mGrid.push_back(dVector(0));
//			dCentreVal = mCenterValues[i];
//			mIsUniformGrid[i] = true;
//			dt = 2*(mUpperBounds[i]) / (mGridSizes[i] - 1);
//			halfIndex = (mGridSizes[i] - 1) / 2;
//			mCenterIndexes.push_back(halfIndex);
//			for (int j = 0; j < halfIndex; j++){
//				mGrid[i].push_back(mLowerBounds[i] + dt*j);
//			}
//			mGrid[i].push_back(dCentreVal);
//			for (int j = halfIndex+1; j < mGridSizes[i]; j++){
//				mGrid[i].push_back(dCentreVal + dt*(j - halfIndex));
//			}
//
//		}
//	}
//}

void Grid::buildCenteredUniformGrid(void){

	for (int i = 0; i < mGridDimension; i++){
		buildCenteredUniformGrid(i);
	}

}

void Grid::buildCenteredUniformGrid(int inIndex){
	double dt;
	double dCentreVal;
	int halfIndex;
	double lVal;
	if (mGridSizes[inIndex] % 2 == 1){
		mGrid.push_back(dVector(0));
		dCentreVal = mCenterValues[inIndex];
		mIsUniformGrid[inIndex] = true;
		dt = 2 * (mUpperBounds[inIndex]) / (mGridSizes[inIndex] - 1);
		halfIndex = (mGridSizes[inIndex] - 1) / 2;
		mCenterIndexes.push_back(halfIndex);
		if (dCentreVal == 0.0){
			mZeroIndexes.push_back(halfIndex);
		}
		for (int j = 0; j < halfIndex; j++){
			lVal = mLowerBounds[inIndex] + dt*j;
			if (lVal == 0.0){
				mZeroIndexes.push_back(j);
			}
			mGrid[inIndex].push_back(lVal);
		}
		mGrid[inIndex].push_back(dCentreVal);
		for (int j = halfIndex + 1; j < mGridSizes[inIndex]; j++){
			lVal = dCentreVal + dt*(j - halfIndex);
			if (lVal == 0.0){
				mZeroIndexes.push_back(j);
			}
			mGrid[inIndex].push_back(lVal);
		}
	}
}
Grid::~Grid(void)
{
}


TimeGrid::TimeGrid(int inGridSize,
	double inLowerBound,
	double inUpperBound,
	const bgreg::date & inToday
	) :mGridSize(inGridSize), mLowerBound(inLowerBound), mUpperBound(inUpperBound), mToday(inToday)
{
	mTimeGrid.resize(mGridSize);
	mIsKeyDates.resize(mGridSize);
}

TimeGrid::TimeGrid(dVector inTimeGrid)
: mGridSize(inTimeGrid.size()), mLowerBound(inTimeGrid[0]), mUpperBound(inTimeGrid[mGridSize-1])
{
}

void TimeGrid::buildUniform()
{
	double dt = (mUpperBound - mLowerBound) / (mGridSize-1);
	for (int j = 0; j < mGridSize; j++){
		mTimeGrid[j]=mLowerBound + dt*j;
		mIsKeyDates[j] = false;
	}
}


void TimeGrid::setKeyDates(const std::vector<bgreg::date> & inDates){
	mNumberOfKeyDates= inDates.size();
	mKeyDates.resize(mNumberOfKeyDates);

	for (int i = 0; i < mNumberOfKeyDates; i++){
		mKeyDates[i] = inDates[i];
	}
}

void TimeGrid::insertExerciseDates(const std::vector<bgreg::date> & inDates){
	
	setKeyDates(inDates);
	mIndexKeyDates.resize(mNumberOfKeyDates);
	double lValue;
	int lIndex;
	int compteur=0;
	for (int i = 0; i <mNumberOfKeyDates; i++){
		lValue = DCC::Actual365Fixed::dayCountFraction(mToday, inDates[i]);
		dVector::iterator it = std::lower_bound(mTimeGrid.begin(), mTimeGrid.end(), lValue);
		lIndex = it - mTimeGrid.begin();
		if (abs(mTimeGrid[lIndex] - lValue)>1/365.0){
			mTimeGrid.insert(it, lValue);
			mIsKeyDates.push_back(false);
			mIsKeyDates[lIndex] = true;
			mIndexKeyDates[i] = lIndex;
			mGridSize++;
		}
		else{
			mTimeGrid[lIndex] = lValue;
			mIndexKeyDates[i] = lIndex;
			mIsKeyDates[lIndex] = true;
		}
	}
	mLowerBound = mTimeGrid[0];
	mUpperBound = mTimeGrid[mGridSize - 1];
	mIsKeyDates[mGridSize - 1] = true;
	mIndexKeyDates[mNumberOfKeyDates - 1] = mGridSize - 1;
}


