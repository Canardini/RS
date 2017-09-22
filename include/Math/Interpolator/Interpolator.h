#ifndef INTERPOLATOR_H
#define INTERPOLATOR_H

#pragma once

#include<vector>
#include <boost/make_shared.hpp>

class Interpolator{

public:
	Interpolator(const std::vector<double> & inX, const std::vector<double> &inY)
		:mX(inX), mY(inY), mNum(inX.size())
	{
	};
	Interpolator(){};

	int getNum(){ return mNum; };
	virtual double calcY(double inX) = 0;
	void setXY(const std::vector<double> & inX, const std::vector <double> & inY);
	virtual void updateInterpolator()=0;
	virtual boost::shared_ptr<Interpolator> clone() const = 0;
	virtual ~Interpolator(){};

protected:
	std::vector<double> mX;
	std::vector<double> mY;
	int mNum;
};

#endif 