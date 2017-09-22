#ifndef OLD_FORWARD_CURVE_H
#define OLD_FORWARD_CURVE_H

#pragma once
#include <boost/date_time/gregorian/gregorian_types.hpp>
#include <boost/shared_ptr.hpp>
#include <Calendar\DayCountConvention.h>
#include<Math/Interpolator/InterpolatorBridge.h>
#include <Calendar\Calendar.h>
namespace bgreg = boost::gregorian;


class ForwardCurve{

public:
	ForwardCurve(const InterpolatorBridge & inInterpolator, bgreg::date inSpotDate,
		         std::string inHolidayCenter, std::string inHolidayRule, std::string inYearDay,
		         std::string inForwardTenor,std::string mForwardName);
	double getDFbyDate(bgreg::date &inDate);
	double getDFbyDouble(double inT);
	double getFWbyDouble(double inT1, double inT2);
	double getFWbyDate(bgreg::date &inDate1, bgreg::date &inDate2);
	double getFWbyDate(bgreg::date &inDate1);

	double getDoubleEndDate(bgreg::date &inDate1);
	bgreg::date getEndDate(bgreg::date &inDate1);
	double getDoubleStartDate(bgreg::date &inDate1);
	void  importTempCurve();

private:
	InterpolatorBridge mInterpolator;
	bgreg::date mSpotDate;
	std::string mHolidayCenter;
	std::string mHolidayRule;
	std::string mYearDay;
	std::string mForwardTenor;
	std::string mForwardName;

};

#endif 