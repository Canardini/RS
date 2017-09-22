#ifndef SUSANO_YIELD_CURVE_H
#define SUSANO_YIELD_CURVE_H


#include <boost/date_time/gregorian/gregorian_types.hpp>
#include "boost/date_time/gregorian/gregorian.hpp" 
#include<YieldCurve\DiscountCurve.h>
#include "SusanooCLR.h"

namespace bgreg = boost::gregorian;

class SusanoYieldCurve : public DiscountCurve {

public:
	SusanoYieldCurve(){};
	SusanoYieldCurve(const bgreg::date & inSpotDate,  const std::string & inTradeBB, const std::string & inTodayDate,
		const std::vector<int> & inCurveNumber);

	boost::shared_ptr<DiscountCurve> clone() const;
	SusanoYieldCurve(const SusanoYieldCurve & inSusanoYieldCurve);
	virtual double getDFbyDate(const bgreg::date &inDate);
	virtual double getDFbyDouble(double inDate){ return 0; };
	virtual ~SusanoYieldCurve(){};
	//void setCurveNumber(long inCurveNumber){ mCurveNumber = inCurveNumber; };
	void setConventions(std::string inHolidayCenter, std::string inHolidayRule, std::string inYearDay, std::string inForwardTenor);
	void setActiveCurve(int inActiveCurve){ mActiveCurve = inActiveCurve; };
	void setActiveQuote(const std::string  & inActiveQuote){ mActiveQuote = inActiveQuote; };

private:
	std::vector<int> mCurveNumber;
	std::string mTradeDB;
	SusanooAPIWrapper mSusanoo;
	int mActiveCurve;
	std::string mActiveQuote;
};

#endif // !APOLLO_CURVE_H