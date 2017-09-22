#ifndef	FORWARD_CURVE_H
#define FORWARD_CURVE_H

#pragma once
#include <boost/date_time/gregorian/gregorian_types.hpp>
#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>

#include <Calendar\DayCountConvention.h>
#include<Calendar\DayCountConvention.h>
#include <Calendar\Calendar.h>
#include <YieldCurve/DiscountCurve.h>


namespace bgreg = boost::gregorian;


class ForwardCurve:public DiscountCurve{

public:
	ForwardCurve(){};
	ForwardCurve(const bgreg::date & inSpotDate);
	ForwardCurve(const InterpolatorBridge & inInterpolator, bgreg::date inSpotDate,
		std::string inHolidayCenter, std::string inHolidayRule, std::string inYearDay, std::string inForwardTenor);
	ForwardCurve(const ForwardCurve&  inForwardCurveObj);

	virtual boost::shared_ptr<DiscountCurve> clone() const;

	virtual double getDFbyDate(bgreg::date &inDate);
	virtual double getDFbyDouble(double inT);
	double getFWbyDouble(double inT1, double inT2);
	double getFWbyDate(bgreg::date &inDate1, bgreg::date &inDate2);
	double getFWbyDate(bgreg::date &inDate1);
	double getDoubleDate(bgreg::date &inDate1);
	double getInstFWrate(double inDate);
	double getDoubleEndDate(bgreg::date &inDate1);
	bgreg::date getEndDate(bgreg::date &inDate1);
	double getDoubleStartDate(bgreg::date &inDate1);

	void importTempCurve();
	virtual ~ForwardCurve(){};

private:

};

#endif 