/*
 * CRF_LogMath.cpp
 *
 * Copyright (c) 2010
 * Author: Jeremy Morris
 *
 */
#include "CRF_LogMath.h"

/*
 *  The various logAdd functions implemented in this file perform addition in log space.  They
 *  all follow the same basic algorithm:
 *
 *  LogAdd algorithm
 * public static double logAdd(double logX, double logY) {
 *     // 1. make X the max
 *     if (logY > logX) {
 *         double temp = logX;
 *         logX = logY;
 *         logY = temp;
 *     }
 *     // 2. now X is bigger
 *     if (logX == Double.NEGATIVE_INFINITY) {
 *         return logX;
 *     }
 *     // 3. how far "down" (think decibels) is logY from logX?
 *     //    if it's really small (20 orders of magnitude smaller), then ignore
 *     double negDiff = logY - logX;
 *     if (negDiff < -20) {
 *         return logX;
 *     }
 *     // 4. otherwise use some nice algebra to stay in the log domain
 *     //    (except for negDiff)
 *     return logX + java.lang.Math.log(1.0 + java.lang.Math.exp(negDiff));
 * }
*/

/*
 *  Implements logAdd algorithm between two numbers
 */
double CRF_LogMath::logAdd(double loga, double logb)
{
	double logx = loga;
	double logy = logb;
	if (logy>logx)
	{
		// Ensure that logx is always the larger of the two
		logy=loga;
		logx=logb;
	}
	double negDiff=logy-logx;
	/*if (negDiff < -20) {
		return logx;
	}*/
	double retVal=0;
	try {
		retVal = logx + logE(1.0 + expE(negDiff));
	}
	catch (exception& e) {
		string errstr="logAdd caught exception "+string(e.what())+" while taking exp";
		throw runtime_error(errstr);
	}
	return retVal;
}

/*
 * Implements logAdd as a summation over an entire vector of values
 */
double CRF_LogMath::logAdd(double* R, int idx) {
	// 1. Find the maximum element
	double max=R[0];
	for (int i=1; i<idx; i++) {
		if (R[i]>max) { max=R[i]; }
	}
	// 2. Create the summation
	double sum=0.0;
	for (int i=0; i<idx; i++) {
		try {
			sum += expE(R[i]-max);
		}
		catch (exception& e) {
			string errstr="logAdd caught exception: "+string(e.what())+" while creating summation";
			throw runtime_error(errstr);
		}
	}
	// 3. Apply the log function to the summation
	double lsum=0.0;
	try {
		lsum=logE(sum);
	}
	catch (exception& e) {
		string errstr="logAdd caught exception: "+string(e.what())+" while taking exp of summation";
		throw runtime_error(errstr);
	}
	return max+lsum;
}

/*
 * As above, but the index of the maximum value in the vector R is known before calling logAdd
 * (saving a search through the vector for the maximum value).
 */
double CRF_LogMath::logAdd(double* R, double max, int idx) {
	// Perform the log add without having to find the maximum value in the array R
	// 2. Create the summation
	double sum=0.0;
	for (int i=0; i<idx; i++) {
		try {
			sum += expE(R[i]-max);
		}
		catch (exception& e) {
			string errstr="logAdd caught exception: "+string(e.what())+" while creating summation";
			throw runtime_error(errstr);
		}
	}
	// 3. Apply the log function to the summation
	double lsum=0.0;
	try {
		lsum=logE(sum);
	}
	catch (exception& e) {
		string errstr="logAdd caught exception: "+string(e.what())+" while taking exp of summation";
		throw runtime_error(errstr);
	}
	return max+lsum;
}

/*
 * As above, but implemented as a vector instead of as an array
 */
double CRF_LogMath::logAdd(vector<double>* R, int idx) {
	// Logic for performing logAdd as a summation over an entire 1D vector
	// 1. Find the maximum element
	double max=R->at(0);
	for (int i=1; i<idx; i++) {
		if (R->at(i)>max) { max=R->at(i); }
	}
	// 2. Create the summation
	double sum=0.0;
	for (int i=0; i<idx; i++) {
		try {
			sum += expE(R->at(i)-max);
		}
		catch (exception& e) {
			string errstr="logAdd caught exception: "+string(e.what())+" while creating summation";
			throw runtime_error(errstr);
		}
	}
	// 3. Apply the log function to the summation
	double lsum=0.0;
	try {
		lsum=logE(sum);
	}
	catch (exception& e) {
		string errstr="logAdd caught exception: "+string(e.what())+" while taking exp of summation";
		throw runtime_error(errstr);
	}
	return max+lsum;
}

/*
 * As above, but implemented as a vector instead of as an array and the maximum value is known
 * before calling the function.
 */
double CRF_LogMath::logAdd(vector<double>* R, double max, int idx) {
	// Perform the log add without having to find the maximum value in the array R
	// 2. Create the summation
	double sum=0.0;
	for (int i=0; i<idx; i++) {
		try {
			sum += expE(R->at(i)-max);
		}
		catch (exception& e) {
			string errstr="logAdd caught exception: "+string(e.what())+" while creating summation";
			throw runtime_error(errstr);
		}
	}
	// 3. Apply the log function to the summation
	double lsum=0.0;
	try {
		lsum=logE(sum);
	}
	catch (exception& e) {
		string errstr="logAdd caught exception: "+string(e.what())+" while taking exp of summation";
		throw runtime_error(errstr);
	}
	return max+lsum;
}

/*
 * protected log function - throws exceptions when given a zero or when the result is infinity or NAN.
 */
double CRF_LogMath::logE(double a)
{
//	if (a==1) { return 0; }
	if (a==0) {
		string errstr="logE caught attempt to take the log of zero when taking log of: "+stringify(a);
		throw overflow_error(errstr);
	}
	double b = log(a);
	if (isnan(b) || isinf(b)) {
		string errstr="logE caught overflow error when taking log of: "+stringify(a);
		throw overflow_error(errstr);
	}
	return b;
}

/*
 * Protected exponential function.  Throws exceptions when the result would overflow or if the result
 * is NAN or infinity.
 */
double CRF_LogMath::expE(double a)
{
//	if (a == 0) { return 1; }
	if (a >= CRF_DBL_LN_MAX) {
		string errstr = "expE caught overflow error when taking exp of: "+stringify(a)+", larger than max value "+stringify(CRF_DBL_LN_MAX);
		throw overflow_error(errstr);
	}
	double b = exp(a);
	if ((isnan(b) || isinf(b))) {
		string errstr = "expE caught overflow error when taking exp of: "+stringify(a);
		throw overflow_error(errstr);
	}
	return b;
}

