#ifndef DAY_COUNT_CONVENTION_H
#define DAY_COUNT_CONVENTION_H

#pragma once

#include <boost/date_time/gregorian/gregorian_types.hpp>

namespace bgreg = boost::gregorian;

namespace DCC{

	int dayCount(const bgreg::date& inStartDate, const bgreg::date& inEndDate);
	double dayCountFraction(const bgreg::date& inStartDate, const bgreg::date& inEndDate,std::string &inDayConvention);
	namespace Actual365
	{
		double dayCountFraction(const bgreg::date& inStartDate, const bgreg::date& inEndDate);
	};

	namespace Actual365Fixed
	{
		double dayCountFraction(const bgreg::date& inStartDate, const bgreg::date& inEndDate);
	};

	namespace Actual360
	{
		double dayCountFraction(const bgreg::date& inStartDate, const bgreg::date& inEndDate);
	};

	namespace Thirty360
	{
		int dayCount(const bgreg::date& inStartDate, const bgreg::date& inEndDate);
		double dayCountFraction(const bgreg::date& inStartDate, const bgreg::date& inEndDate);
	};
}

#endif