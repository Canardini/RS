#ifndef APOLLO_VOLATILITY_CUBE_H
#define APOLLO_VOLATILITY_CUBE_H

#pragma once

//#include "stdafx.h"
#include "windows.h"
#include <sql.h>
#include <sqltypes.h>
#include <sqlext.h>
#include <string>
#include <vector>

#include <boost/date_time/gregorian/gregorian_types.hpp>
#include "boost/date_time/gregorian/gregorian.hpp" 
#include<Math/Interpolator/2DNaturalCubicSpline.h>
#include<Math/Interpolator/3DNaturalCubicSpline.h>
#include<Math/Interpolator/2DLinearInterpolator.h>
#include<Math/Interpolator/3DLinearInterpolator.h>

#include<Math/Interpolator/2DInterpolatorBridge.h>
#include<Math/Interpolator/3DInterpolatorBridge.h>



namespace bgreg = boost::gregorian;

class ApolloVolatilityCube{

public:
	ApolloVolatilityCube(std::string & inToday, std::string & inDB, std::string &inCurrency,std::string &inInterpolationMethod);
	std::string  convertDate(std::string & inToday1);
	std::map<double, int> convertDoubleVecIntoInt(const std::vector<double> & inX);
	std::vector<double> getSubSortedVec(const std::vector<double> & inX);
	std::vector<std::vector<std::vector<double>>> getCubeValue(std::vector<double> & inOriginalX, std::vector<double> & inOriginalY,
													           std::vector<double> & inOriginalZ, std::vector<double> & inXYZValue,
												               std::map<double, int> & inX, std::map<double, int> & inY, std::map<double, int> & inZ);
	std::vector<std::vector<double>> getSurfaceValue(std::vector<double> & inOriginalX, std::vector<double> & inOriginalY, 
		                                             std::vector<double> & inXYValue,
			                                         std::map<double, int> & inX, std::map<double, int> & inY);

	
	double getVol(bgreg::date &inExerciseDate, bgreg::date & inTerminationDate, double StikeAtmspread);

	void addATMhedge(bgreg::date &inExerciseDate, bgreg::date & inTerminationDate);

	void updateShiftedSpline();
	void shifteSpline(int i, int j, double shift);
	int getMonth(const int unit, const int count)
	{
		int number = 0;
		switch (unit)
		{
		case 4:
			break;
		case 30:
			number = count;
			break;
		case 40:
			number = count * 12;
			break;
		default:
			break;
		}
		return number;
	}
	double getHedgeVol(bgreg::date &inExerciseDate, bgreg::date & inTerminationDate, double StikeAtmspread);
	int getNumberOfVegaHedges(){ return mNumberOfVegaHedges; };
	int getsizeATMX(){ return m2dSpline.getsizeX();};
	int getsizeATMY(){ return m2dSpline.getsizeY();};


	bool isHedge(int i, int j){ return mIsATMHedged[i][j]; }
	


private:
	TwoDInterpolatorBridge m2dSpline;
	ThreeDInterpolatorBridge m3dSpline;
	
	std::vector<std::vector<bool>> mIsATMHedged;

	TwoDInterpolatorBridge m2dSplineShifted;
	ThreeDInterpolatorBridge m3dSplineShifted;

	bgreg::date mToday;
	int mNumberOfVegaHedges;
};

	class DBconn
	{
		//private: 
	protected:
		SQLHENV henv;
		SQLHDBC hdbc;
		SQLHSTMT hstmt;
		SQLRETURN rc;
		std::vector<std::string> dbname;
		std::string userID;
		std::string password;

	public:
		DBconn()
		{
			henv = NULL;
			hdbc = NULL;
			hstmt = NULL;
			rc = 0;
		}

		DBconn(SQLHENV mhenv, SQLHDBC mhdbc, SQLHSTMT mhstmt, SQLRETURN mrc, std::vector<std::string> mdbname, std::string musername, std::string mpassword)
		{
			henv = mhenv;
			hdbc = mhdbc;
			hstmt = mhstmt;
			rc = mrc;
			dbname = mdbname;
			userID = musername;
			password = mpassword;
		}

		DBconn(const DBconn &inDBconn) : henv(inDBconn.henv), hdbc(inDBconn.hdbc), hstmt(inDBconn.hstmt), rc(inDBconn.rc), dbname(inDBconn.dbname), userID(inDBconn.userID), password(inDBconn.password){}

		DBconn & operator = (const DBconn &inDBconn) {
			henv = inDBconn.henv;
			hdbc = inDBconn.hdbc;
			hstmt = inDBconn.hstmt;
			rc = inDBconn.rc;
			dbname = inDBconn.dbname;
			userID = inDBconn.userID;
			password = inDBconn.password;
			return *this;
		}

		void init() {
			/*	SQLHENV henv = NULL;
			SQLHDBC hdbc = NULL;
			SQLHSTMT hstmt = NULL;
			SQLRETURN rc = 0; */

			std::string connstr = "Dsn=";
			connstr += dbname[0] + ";";
			//connstr += ";UID=dbdvlp;PWD=W!eb4mcmc;";
			connstr += "UID=" + userID + ";";
			connstr += "PWD=" + password + ";";
			connstr += "DATABASE=" + dbname[1] + ";";
			connstr += "SRVR=" + dbname[2] + ";";
			std::string ws(connstr.begin(), connstr.end());
			SQLCHAR * ConnStrIn = (SQLCHAR *)ws.c_str();

			rc = SQLAllocHandle(SQL_HANDLE_ENV, SQL_NULL_HANDLE, &henv);
			rc = SQLSetEnvAttr(henv, SQL_ATTR_ODBC_VERSION, (void*)SQL_OV_ODBC3, 0);
			rc = SQLAllocHandle(SQL_HANDLE_DBC, henv, &hdbc);

			rc = SQLDriverConnect(hdbc, NULL, ConnStrIn, SQL_NTS, NULL, 0, NULL, SQL_DRIVER_NOPROMPT);
			if (rc != SQL_SUCCESS) {
				std::string errmsg = "Cannot connect DataBase: " + dbname[0] + "\n";
				errmsg += "Please check User ID and Password.";
				throw errmsg;
				//exit(1);
			}

			rc = SQLAllocStmt(hdbc, &hstmt);
		}

		SQLHSTMT gethstmt(){ return hstmt; }

		~DBconn() {};
	};


#endif
