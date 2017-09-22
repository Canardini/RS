#include <Calendar\Serial.h>


int Serial::serialFromDate(int inYear, int inMonth, int inDay)
{
	// As only VBA supports negative dates. See also remark in function declaration.

	// Excel/Lotus 123 have a bug with 29-02-1900. 1900 is not a
	// leap year, but Excel/Lotus 123 think it is...
	if ((inDay == 29) & (inMonth == 02) & (inYear == 1900))
	{
		return 60;
	}

	// DMY to Modified Julian calculation with an extra subtraction of 2415019.
	int nSerialDate =
		int((1461 * (inYear + 4800 + int((inMonth - 14) / 12))) / 4) +
		int((367 * (inMonth - 2 - 12 * ((inMonth - 14) / 12))) / 12) -
		int((3 * (int((inYear + 4900 + int((inMonth - 14) / 12)) / 100))) / 4) +
		inDay - 2415019 - 32075;

	if (nSerialDate < 61)
	{
		// Because of the 29-02-1900 bug, any serial date
		// under 60 is one off... Compensate. We need
		// to check for 61 here as the serial date calculation
		// is correct and returns 60 for 28 Feb 1900.
		--nSerialDate;
	}

	return (int)nSerialDate;
}

int Serial::serialFromDate(const bgreg::date &inDate)
{
	return serialFromDate(inDate.year(), inDate.month(), inDate.day());
}

bgreg::date Serial::dateFromSerial(int inSerialDate)
{
	int iYear, iMonth, iDay;
	dateFromSerial(inSerialDate, iYear, iMonth, iDay);
	bgreg::date oDate(iYear, iMonth, iDay);
	return oDate;
}

void Serial::dateFromSerial(int inSerialDate, int &outYear, int &outMonth, int &outDay)
{
	// Only VBA supports negative dates. See also remark in function declaration.

	// Excel/Lotus 123 have a bug with 29-02-1900. 1900 is not a
	// leap year, but Excel/Lotus 123 think it is...
	if (inSerialDate == 60)
	{
		outDay = 29;
		outMonth = 2;
		outYear = 1900;

		return;
	}
	else if (inSerialDate < 60)
	{
		// Because of the 29-02-1900 bug, any serial date
		// under 60 is one off... Compensate.
		++inSerialDate;
	}

	// Modified Julian to DMY calculation with an addition of 2415019
	int lTerm1 = inSerialDate + 68569 + 2415019;
	int lTerm2 = int((4 * lTerm1) / 146097);
	lTerm1 = lTerm1 - int((146097 * lTerm2 + 3) / 4);
	int lTerm3 = int((4000 * (lTerm1 + 1)) / 1461001);
	lTerm1 = lTerm1 - int((1461 * lTerm3) / 4) + 31;
	int lTerm4 = int((80 * lTerm1) / 2447);
	outDay = lTerm1 - int((2447 * lTerm4) / 80);
	lTerm1 = int(lTerm4 / 11);
	outMonth = lTerm4 + 2 - (12 * lTerm1);
	outYear = 100 * (lTerm2 - 49) + lTerm3 + lTerm1;
}
