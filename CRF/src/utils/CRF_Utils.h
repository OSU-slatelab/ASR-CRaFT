#ifndef CRF_UTILS_H_
#define CRF_UTILS_H_
/*
 * CRF_Utils.h
 *
 * Copyright (c) 2010
 * Author: Jeremy Morris
 *
 * Contains code for various utility functions that are not part of other classes.
 */

#include <iostream>
#include <sstream>
#include <string>
#include <stdexcept>

// This code taken from the C++ FAQ Lite
//		http://www.parashift.com/c++-faq-lite/misc-technical-issues.html#faq-39.1

inline std::string stringify(double x)
 {
   std::ostringstream o;
   if (!(o << x))
     throw std::runtime_error("stringify(double)");
   return o.str();
 }

#endif /*CRF_UTILS_H_*/
