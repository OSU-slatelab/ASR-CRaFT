/*
 * CRF_MLFManager.h
 *
 * Copyright (c) 2010
 * Author: Jeremy Morris
 *
 * Contains the class definitions for CRF_MLFManager
 * Uses to read the MLF file format used by HTK.
 */
#ifndef CRF_MLFMANAGER_H_
#define CRF_MLFMANAGER_H_
#include "fst/fstlib.h"
#include <stdexcept>
#include <vector>

using namespace fst;
using namespace std;

/*
 * class CRF_MLFManager
 *
 * Reads the MLF file format used by HTK.  Uses the MLF to generate lattices in OpenFST format.
 *
 */

class CRF_MLFManager {
private:
	SymbolTable* symTab;
	vector< vector<int>* > transcripts;
	map<string,int > fnameTable;
	string getKey(string fname);
public:
	CRF_MLFManager(char* mlffile, char* olist, SymbolTable* symTab);
	virtual ~CRF_MLFManager();
	void readMLF(char* mlffile);
	StdVectorFst* getFst(string fname);
};

#endif /* CRF_MLFMANAGER_H_ */
