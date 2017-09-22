#ifndef CALENDAR_H
#define CALENDAR_H


#include<map>
#include<vector>
#include <fstream>
#include <iostream>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/split.hpp>

#include "boost/date_time/gregorian/gregorian.hpp" 
namespace bgreg = boost::gregorian;

class Calendar{

public:

	Calendar(){};

	static bgreg::date getRespondDate(const bgreg::date &target,std::string & inTenor,
							   std::string &inCityCenter,std::string &inRule);

	static bgreg::date getNextBusinessDay(const bgreg::date &target, std::string &inCityCenter, int inNdays);
	static bgreg::date addDays(const bgreg::date &target, int inNdays);
	static bool isHolidayCenter(const bgreg::date &target, std::string &inCityCenter);
	static bool isHoliday(const bgreg::date &target, std::string &inCityCenter);
	static bgreg::date Calendar::addMonth(const bgreg::date &target, int inNmonths, int inFixedDays, bool isFebAdj, std::string &inCityCenter);
	
	
	static std::map<std::string, std::vector<bgreg::date>>mHolidayMap;
	
	static long importCalendar();
	
	static bgreg::date Calendar::addYear(const bgreg::date &target, int inNYears, int inFixedDays,
		bool isFebAdj, std::string &inCityCenter);
	static int getTenor(std::string & inTenor);
	static std::string Calendar::getTenorType(std::string & inTenor);
	static bool Calendar::isMonthEnd(const bgreg::date &target, std::string &inCityCenter);
};


#endif