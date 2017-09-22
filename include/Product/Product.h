#pragma once
#include <boost/date_time/gregorian/gregorian_types.hpp>
class Product
{
public:
	Product();
private:
	boost::gregorian::date mValuationDate;
	~Product();
};

