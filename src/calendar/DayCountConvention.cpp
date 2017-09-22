#include <Calendar\DayCountConvention.h>


double DCC::dayCountFraction(const bgreg::date& inStartDate, 
								  const bgreg::date& inEndDate,
	                              std::string &inDayConvention){

	if(inDayConvention=="Act365"){
		return DCC::Actual365::dayCountFraction(inStartDate,inEndDate);
	}
	else if (inDayConvention == "Act360"){
		return DCC::Actual360::dayCountFraction(inStartDate, inEndDate);
	}
	else if (inDayConvention == "30/360"){
		return DCC::Thirty360::dayCountFraction(inStartDate, inEndDate);
	}
	else{
		return DCC::Actual365Fixed::dayCountFraction(inStartDate, inEndDate);
	}

}

int DCC::dayCount(const bgreg::date& inStartDate, const bgreg::date& inEndDate) 
{
	return (inEndDate - inStartDate).days();
}

double DCC::Actual365::dayCountFraction(const bgreg::date& inStartDate, const bgreg::date& inEndDate) 
{
	const int iStartYear = inStartDate.year();
	const int iEndYear = inEndDate.year();
	const int iNumDaysStartYear = bgreg::gregorian_calendar::is_leap_year(iStartYear) ? 366 : 365;
	const int iNumDaysEndYear = bgreg::gregorian_calendar::is_leap_year(iEndYear) ? 366 : 365;

	return static_cast<double>(iNumDaysStartYear - inStartDate.day_of_year() + 1) / iNumDaysStartYear
		+ static_cast<double>(inEndDate.day_of_year() - 1) / iNumDaysEndYear
		+ static_cast<double>(iEndYear - iStartYear - 1);
}

double DCC::Actual365Fixed::dayCountFraction(const bgreg::date& inStartDate, const bgreg::date& inEndDate)
{
	return dayCount(inStartDate, inEndDate) / 365.0;
}

double DCC::Actual360::dayCountFraction(const bgreg::date& inStartDate, const bgreg::date& inEndDate)
{
	return dayCount(inStartDate, inEndDate) / 360.0;
}

int DCC::Thirty360::dayCount(const bgreg::date& inStartDate, const bgreg::date& inEndDate) 
{
	const int				iYears(inEndDate.year() - inStartDate.year());
	const int				iMonths(inEndDate.month() - inStartDate.month());
	const bgreg::greg_day	iStartDay(inStartDate.day());
	const bgreg::greg_day	iEndDate(inEndDate.day());
	const bgreg::greg_day	iModifiedStartDay((iStartDay == 31) ? bgreg::greg_day(30) : iStartDay);
	const bgreg::greg_day	iModifiedEndDay(((iEndDate == 31) & (iModifiedStartDay == 30)) ? bgreg::greg_day(30) : iEndDate);
	const int				iDays(iModifiedEndDay - iModifiedStartDay);

	return 360 * iYears + 30 * iMonths + iDays;
}

double DCC::Thirty360::dayCountFraction(const bgreg::date& inStartDate, const bgreg::date& inEndDate)
{
	return DCC::Thirty360::dayCount(inStartDate, inEndDate) / 360.0;
}
