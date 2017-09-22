#ifndef DISCOUNT_CURVE_H
#define DISCOUNT_CURVE_H

#pragma once
#include <boost/date_time/gregorian/gregorian_types.hpp>
#include <Calendar\DayCountConvention.h>
#include <boost/shared_ptr.hpp>
#include<Math/Interpolator/DefaultInterpolator.h>

#include<Math/Interpolator/InterpolatorBridge.h>

#include<Calendar\DayCountConvention.h>
#include <Calendar\Calendar.h>
#include <boost/make_shared.hpp>
namespace bgreg = boost::gregorian;


class DiscountCurve{

public:
	DiscountCurve(){};
	DiscountCurve(const bgreg::date & inTodayDate) :mTodayDate(inTodayDate){
		DefaultInterpolator lInt;
		InterpolatorBridge lIB(lInt);
		mInterpolator = lIB;
		};
	DiscountCurve(const InterpolatorBridge & inInterpolator, bgreg::date inTodayDate,
		std::string inHolidayCenter, std::string inHolidayRule, std::string inYearDay, std::string inForwardTenor);
	DiscountCurve(const DiscountCurve&  inDiscountCurveObj);

	virtual boost::shared_ptr<DiscountCurve> clone() const;

	
	virtual double getDFbyDate(const bgreg::date &inDate);
	virtual double getDFbyDouble(double inT);
	double getFWbyDouble(double inT1, double inT2);
	double getFWbyDate(const bgreg::date &inDate1, const bgreg::date &inDate2);
	double getDoubleDate(const bgreg::date &inDate1);
	double getInstFWrate(double inDate);
	double getFWbyDate(const bgreg::date &inDate1);
	double getDoubleEndDate(const bgreg::date &inDate1);
	bgreg::date getEndDate(const bgreg::date &inDate1);
	double getDoubleStartDate(const bgreg::date &inDate1);
	double getAccDate(const bgreg::date &inDate1, const bgreg::date &inDate2);

	
	void importTempCurve();
	virtual ~DiscountCurve(){};

protected :

	InterpolatorBridge  mInterpolator;
	bgreg::date mTodayDate;
	std::string mHolidayCenter;
	std::string mHolidayRule;
	std::string mYearDay;
	std::string mForwardTenor;

};

#endif 