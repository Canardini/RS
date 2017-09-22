#ifndef SERIAL_H
#define SERIAL_H

#pragma once

#include <boost/date_time/gregorian/gregorian_types.hpp>

namespace bgreg = boost::gregorian;

namespace Serial{
	void dateFromSerial(int inSerialDate, int &outYear, int &outMonth, int &outDay);
	int serialFromDate(const bgreg::date &inDate);
	int serialFromDate(int inYear, int inMonth, int inDay);
	bgreg::date dateFromSerial(int inSerialDate);
}

#endif