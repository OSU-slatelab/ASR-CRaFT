#include <iostream>
#include <stdexcept>

#include "float.h"

#include "CRF_LogMath.h"

using namespace std;
using namespace CRF_LogMath;

int main(int argc, const char* argv[])
{
	double result=0.0;
	try {
		result=logAdd(DBL_MAX*2,DBL_MAX*2);
	}
	catch(exception& e) {
		cerr << "LogMathTest caught exception: " << e.what() << endl;
	}
	cout << "Result: " << result << endl;
	
	double* test_array = new double[10];
	for (int i=0; i<10; i++) {
		test_array[i]=DBL_MAX*2;
	}
	try {
		result=logAdd(test_array,10);
	}
	catch(exception& e) {
		cerr << "LogMathTest caught exception: " << e.what() << endl;
	}
	cout << "Result: " << result << endl;
	
	try {
		result=logAdd(test_array,test_array[0],10);
	}
	catch(exception& e) {
		cerr << "LogMathTest caught exception: " << e.what() << endl;
	}
	cout << "Result: " << result << endl;
	try {
		result=logE(0);
	}
	catch(exception& e) {
		cerr << "LogMathTest caught exception: " << e.what() << endl;
	}
	cout << "Result: " << result << endl;
	try {
		result=logE(1.0/0);
	}
	catch(exception& e) {
		cerr << "LogMathTest caught exception: " << e.what() << endl;
	}
	cout << "Result: " << result << endl;
	try {
		result=expE(DBL_MAX);
	}
	catch(exception& e) {
		cerr << "LogMathTest caught exception: " << e.what() << endl;
	}
	cout << "Result: " << result << endl;
	try {
		result=expE(LOG0);
	}
	catch(exception& e) {
		cerr << "LogMathTest caught exception: " << e.what() << endl;
	}
	cout << "exp(LOG0) Result: " << result << endl;
	try {
		result=expE(DBL_MIN);
	}
	catch(exception& e) {
		cerr << "LogMathTest caught exception: " << e.what() << endl;
	}
	cout << "Result: " << result << endl;
	try {
		result=logE(DBL_MIN);
	}
	catch(exception& e) {
		cerr << "LogMathTest caught exception: " << e.what() << endl;
	}
	cout << "Result: " << result << endl;
	try {
		result=logAdd(LOG0,LOG0);
	}
	catch(exception& e) {
		cerr << "LogMathTest caught exception: " << e.what() << endl;
	}
	cout << "Result LOG0+LOG0: " << result << endl;
	try {
		result=expE(logAdd(LOG0,LOG0));
	}
	catch(exception& e) {
		cerr << "LogMathTest caught exception: " << e.what() << endl;
	}
	cout << "Result exp(LOG0+LOG0): " << result << endl;
}
