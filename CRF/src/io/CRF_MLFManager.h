/*
 * CRF_MLFManager.h
 *
 *  Created on: Sep 15, 2009
 *      Author: morrijer
 */

#ifndef CRF_MLFMANAGER_H_
#define CRF_MLFMANAGER_H_
#include "fst/fstlib.h"
#include <stdexcept>
#include <vector>

using namespace fst;
using namespace std;

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
