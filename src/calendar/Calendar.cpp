#include<Calendar\Calendar.h>


int Calendar::getTenor(std::string & inTenor){
	return std::atoi(inTenor.c_str());
}

std::string Calendar::getTenorType(std::string & inTenor){
	int n = inTenor.size();

	for (int i = 0; i < n; i++){
		if(std::isalpha(inTenor[i])){
			return inTenor.substr(i);
		}
	}
	return inTenor;
}

std::map<std::string, std::vector<bgreg::date> > Calendar::mHolidayMap;

long Calendar::importCalendar(){
	std::ifstream infile("F:/rnd/RSCheyette/Calendar.csv", std::ios::in);
	std::string lcontent;
	std::vector<std::string> lCityCenters;
	std::vector<std::string> lDates;

	if (infile){ // Pick the city centers
		std::getline(infile, lcontent);
		boost::split(lCityCenters, lcontent, boost::is_any_of(","));
	}
 

	while (!infile.eof()){
		std::getline(infile, lcontent);
		std::vector<std::string> lDates;
		boost::split(lDates, lcontent, boost::is_any_of(","));
		if (lcontent != ""){
			for (size_t i = 0; i < lCityCenters.size(); i++){
				if (lDates[i] != ""){
					double test = 0.9;
					bgreg::date d(bgreg::from_us_string(lDates[i]));
					mHolidayMap[lCityCenters[i]].push_back(d);
				}
			}
		}
	}
	return 1;
}



bgreg::date Calendar::getRespondDate(const bgreg::date &target, std::string & inTenor,
	std::string &inCityCenter, std::string &inRule)
{
	std::string sTenorType = getTenorType(inTenor);

	bool lFlag = (sTenorType == "SN") |( sTenorType == "ON") | (sTenorType == "O/N");
	
	if(lFlag){
		return getNextBusinessDay(target, inCityCenter, 1);
	}
	bgreg::date lRes;
	lFlag = (sTenorType == "T/N");
	if (lFlag){
		lRes = getNextBusinessDay(target, inCityCenter, 1);
		lRes = getNextBusinessDay(lRes, inCityCenter, 1);
		return lRes;
	}

	lFlag = (sTenorType == "SW");
	if (lFlag){
		lRes = addDays(target,7+1);
		lRes = getNextBusinessDay(lRes, inCityCenter, -1);

		return lRes;
	}


		
	int amount = getTenor(inTenor);
	lFlag = (sTenorType == "W");
	if (lFlag){
		lRes = addDays(target, 7 * amount);
	}
	lFlag = (sTenorType == "M");
	if (lFlag){
		lRes =addMonth(target, amount, 0, false, inCityCenter);
	}

	lFlag = (sTenorType == "Y");
	if (lFlag){
		lRes = addYear(target, amount, 0, false, inCityCenter);
	}

	if (isHoliday(lRes, inCityCenter) == true){
		int lDirect = 1;
		if ((inRule == "PR") | ((inRule == "MDFL") & isMonthEnd(lRes, inCityCenter))){
			lDirect = -1;
		}
		
		lRes = getNextBusinessDay(lRes, inCityCenter, lDirect);
		
	}
	return lRes;
}

bgreg::date Calendar::getNextBusinessDay(const bgreg::date &target, std::string &inCityCenter, int inNdays){
	bgreg::date lRes;
	lRes = addDays(target, inNdays);
	if (isHoliday(lRes, inCityCenter) == true){
		lRes = getNextBusinessDay(lRes,inCityCenter, 1);
		
	}
	return lRes;
}

bool Calendar::isHolidayCenter(const bgreg::date &target, std::string &inCityCenter){
	bgreg::date::day_of_week_type oDayOfWeek = target.day_of_week();

	if ((oDayOfWeek == 0) | (oDayOfWeek == 6)){
		return true;
	}

	for (size_t i = 0; i < mHolidayMap[inCityCenter].size(); i++){
		if (target == mHolidayMap[inCityCenter][i]){
			return true;
		}
		else if (target<mHolidayMap[inCityCenter][i]){
			break;
		}
	}
	return false; 
}

bool Calendar::isHoliday(const bgreg::date &target, std::string &inCityCenter){
	std::vector<std::string> lCityCenters;
	boost::split(lCityCenters, inCityCenter, boost::is_any_of("+"));
	bool lTemp;
	for (size_t i = 0; i < lCityCenters.size(); i++){
		lTemp=isHolidayCenter(target, lCityCenters[i]);
		if (lTemp == true){
			return true;
		}
	}
	return false;
}

bool Calendar::isMonthEnd(const bgreg::date &target, std::string &inCityCenter){
	bgreg::date lRes;
	lRes = getNextBusinessDay(target, inCityCenter, 1);
	if (lRes.month() == target.month()){
		return false;
	}
	else{
		return true;
	}
}

bgreg::date Calendar::addDays(const bgreg::date &target, int inNdays){
	bgreg::date lRes;
	if (inNdays != 0){
		if (inNdays > 0){
			bgreg::days lDays(inNdays);
			lRes = target + lDays;
			return lRes;
		}

		else{
			bgreg::days lDays(-inNdays);
			lRes = target - lDays;
			return lRes;
		}
	}
	return target;
}


bgreg::date Calendar::addMonth(const bgreg::date &target, int inNmonths, int inFixedDays, 
	                           bool isFebAdj, std::string &inCityCenter){
	
	int lTempday;
	bool lIsFeb=false;
	if (inFixedDays == 0){
		lTempday = target.day();
	}
	else{
		lTempday = inFixedDays;
	}
	if (((target.month() + inNmonths) % 12 == 2) & (lTempday > 28)){
		lIsFeb = true;
	}

	bgreg::months l(inNmonths);
	bgreg::date lRes = target + l;
	if (isFebAdj & lIsFeb){
		
	}
	else{
		int days = target.day() - lRes.day();
		bgreg::days d(days);
		lRes = lRes + d;
	}
	return lRes;
}

bgreg::date Calendar::addYear(const bgreg::date &target, int inNyears, int inFixedDays,
	bool isFebAdj, std::string &inCityCenter){

	bgreg::years l(inNyears);

	return target + l;

}
