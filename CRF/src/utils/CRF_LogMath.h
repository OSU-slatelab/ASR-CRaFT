#ifndef CRF_LOGMATH_H_
#define CRF_LOGMATH_H_
/*
 * CRF_LogMath.h
 *
 * Copyright (c) 2010
 * Author: Jeremy Morris
 *
 * Implements functions pertaining to log math operations.
 *   Specifically, implements various forms of addition in the log space and a protected
 *   exponential function, as well as constants for log math.
 */
#include <stdexcept>
#include <iostream>
#include <string>
#include <vector>

#include "math.h"
#include "float.h"
#include "CRF_Utils.h"

using namespace std;

namespace CRF_LogMath {
	const double CRF_DBL_LN_MAX = log(DBL_MAX);
	const double LOG0 = -1*DBL_MAX;
	double logAdd(double a, double b);
	double logAdd(double* R, int idx);
	double logAdd(double* R, double max, int idx);
	double logAdd(vector<double>* R, int idx);
	double logAdd(vector<double>* R, double max, int idx);
	double logE(double a);
	double expE(double a);
}

#endif /*CRF_LOGMATH_H_*/
