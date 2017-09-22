#ifndef GRID_H
#define GRID_H
#pragma once

#include <Math/Matrix/Matrix.h>
#include <Calendar/DayCountConvention.h>
#include<vector>
#include <boost/date_time/gregorian/gregorian_types.hpp>
namespace bgreg = boost::gregorian;

class Grid
{
public:
	explicit Grid(int inGridDimension,
		 dVector &inLowerBounds,
		 dVector &inUpperBounds,
		 dVector &inCenterValues,
		 std::vector<int> &inGridSizes);

	explicit Grid(int inGridDimension,
		dVector &inLowerBounds,
		dVector &inUpperBounds,
		dVector &inCenterValues,
		std::vector<int> &inGridSizes,
		double inAlpha);

	double getGridElement(int index1, int index2){return mGrid[index1][index2];};
	void buildUniformGrid();
	void buildCenteredUniformGrid();
	void buildCenteredUniformGrid(int inIndex);
	void buildArcSinus(int inIndex);
	int locate(int i,double inX0);

	bool getIsUniformGrid(int index){return mIsUniformGrid[index];};
	int getCenterIndexes(int index1){ return mCenterIndexes[index1]; };
	int getZeroIndexes(int index1){ return mZeroIndexes[index1]; };

	int getGridSize(int index);
	~Grid(void);
private:

	int mGridDimension;
	std::vector<bool>   mIsUniformGrid;
	std::vector<double> mLowerBounds;
	std::vector<double> mCenterValues;
	std::vector<int> mCenterIndexes;
	std::vector<int> mZeroIndexes;
	std::vector<double> mUpperBounds;
	std::vector<int> mGridSizes;
	std::vector<std::vector<double>> mGrid;
	double mAlpha;

};

class TimeGrid {
public:
	TimeGrid(int inGridSize,
		     double inLowerBound,
			 double inUpperBound, const bgreg::date & inToday
		);

	TimeGrid(dVector inTimeVector);
	void buildUniform();
	
	double operator[](int i){ return mTimeGrid[i]; };
	int getGridSize(){ return mGridSize; };
	bool getIsKeyDate(int inIndex){ return mIsKeyDates[inIndex]; };
	void setKeyDates(const std::vector<bgreg::date> & inDates);
	bgreg::date getKeyDate(int inIndex){ return mKeyDates[inIndex]; };
	void insertExerciseDates(const std::vector<bgreg::date> & inDates);
	int getNumOfKeyDates(){ return mNumberOfKeyDates; };
	int getIndexKeyDates(int inKeyDateNumber){ return mIndexKeyDates[inKeyDateNumber]; };

private :
	dVector mTimeGrid;
	int mGridSize;
	double mLowerBound;
	double mUpperBound;
	bgreg::date mToday;
	std::vector<bool> mIsKeyDates;
	std::vector<int> mIndexKeyDates;
	int mNumberOfKeyDates;
	std::vector<bgreg::date> mKeyDates;


};


#endif
