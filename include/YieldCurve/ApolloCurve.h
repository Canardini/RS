#ifndef APOLLO_CURVE_H
#define APOLLO_CURVE_H


#include <boost/date_time/gregorian/gregorian_types.hpp>
#include "boost/date_time/gregorian/gregorian.hpp" 
#import "C:\Program Files (x86)\MCMRiskDev\MCMAnalyticsComLib.tlb" raw_interfaces_only
using namespace MCMAnalyticsComLib;
#include<vector>
#include<YieldCurve\DiscountCurve.h>
#include <Calendar\Serial.h>

namespace bgreg = boost::gregorian;

class ApolloCurve : public DiscountCurve {

public:
	ApolloCurve(){};
	ApolloCurve(const bgreg::date & inSpotDate, const std::string & inUserID, const std::string & inPassword, const std::string & inTradeBB,
				const std::vector<std::string> & inCurveName, long inCurveNumber);

	boost::shared_ptr<DiscountCurve> clone() const;
	ApolloCurve(const ApolloCurve & inApolloCurve);
	virtual double getDFbyDate(const bgreg::date &inDate);
	virtual double getDFbyDouble(double inDate);
	virtual ~ApolloCurve(){};
	void setCurveNumber(long inCurveNumber){ mCurveNumber = inCurveNumber; };
	void setConventions(std::string inHolidayCenter, std::string inHolidayRule, std::string inYearDay, std::string inForwardTenor);

private :

	std::string mUserID;
	std::string mPassword;
	std::string mTradeDB;
	std::vector<std::string> mCurveName;
	DbEngineXfacePtr mDbengine;
	long mCurveNumber;
};

#endif // !APOLLO_CURVE_H