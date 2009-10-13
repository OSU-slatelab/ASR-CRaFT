/*
 * CRF_MLFManager.cpp
 *
 *  Created on: Sep 15, 2009
 *      Author: morrijer
 */

#include "CRF_MLFManager.h"
typedef StdArc::StateId StateId;
using namespace fst;

CRF_MLFManager::CRF_MLFManager(char* mlffile, char* olist, SymbolTable* symTab)
	: symTab(symTab){
	// TODO Auto-generated constructor stub
	this->readMLF(mlffile);
}

CRF_MLFManager::~CRF_MLFManager() {
	// TODO Auto-generated destructor stub
	cout << "in CRF_MLFManager destructor" << endl;
}

string CRF_MLFManager::getKey(string fname) {
	// Takes in a filename in HTK MLF format (e.g. "*/filename.lab")
	// Returns just the key portion (e.g. "filename")

	int lastc=fname.rfind(".");
	int firstc=fname.rfind("/");
	return(fname.substr(firstc+1,lastc-firstc-1));
}

void CRF_MLFManager::readMLF(char* mlffile) {
	std::ifstream ifile;
	ifile.open(mlffile);
	bool isMLF = false;
	bool nextIsFname = false;
	if (this->symTab == NULL ) {
		string errstr="SymbolTable undefined in CRF_MLFManager";
		throw runtime_error(errstr);
	}
	int count=0;
	if (ifile.is_open()) {
		while (!ifile.eof()) {
			std::string s;
			getline(ifile,s);
			if (s == string("#!MLF!#")) {
				isMLF=true;
				nextIsFname=true;
			}
			else {
				if (!isMLF) {
					string errstr="File "+string(mlffile)+" is not a wellformed MLF";
					throw runtime_error(errstr);
				}
				if (s == "") {
					// Do nothing on blank lines
				}
				else if (s[0]=='"' && s[s.size()-1]=='"') {
					string key=this->getKey(s);
					if (key == "") {
						string errstr="CRF_MLFManager error finding key in string: "+s;
						throw runtime_error(errstr);
					}
					//cout << "LABEL: " << s << "\t" << key << endl;
					this->transcripts.push_back(new vector<int>);
					this->fnameTable[key]=count;
				}
				else if (s[0]=='.' && s.size()<=1) {
					//cout << "EOS: " << s << endl;
					/*vector<int>::iterator itTrans;
					for(itTrans = this->transcripts.at(count)->begin(); itTrans != this->transcripts.at(count)->end(); itTrans++)
					{
						cout << *itTrans << " ";
					}
					cout << endl;*/
					count++;
				}
				else {
					try {
						//cout << s << "\t" << symTab->Find(s.c_str()) << endl;
						this->transcripts.at(count)->push_back(symTab->Find(s.c_str()));
					}
		    		catch (exception &e) {
		    			cerr << "Exception: " << e.what() << endl;
		    		}
				}
			}
		}
		ifile.close();
	}
	else {
		string errstr="Unable to read MLF from file "+string(mlffile);
		throw runtime_error(errstr);
	}

}

StdVectorFst* CRF_MLFManager::getFst(string fname) {
	string key = this->getKey(fname);
	if (key == "") {
		string errstr="Unable to acquire key from filename "+fname;
		throw runtime_error(errstr);
	}
	int idx=this->fnameTable[key];
	StdVectorFst* fst = new StdVectorFst();
	StateId startState = fst->AddState();
	fst->SetStart(startState);
	StateId prevState=startState;
	StateId curState;

	//cout << key << "\t";
	vector<int>::iterator itNode;
	for(itNode = this->transcripts.at(idx)->begin(); itNode != this->transcripts.at(idx)->end(); itNode++)
	{
		curState=fst->AddState();
		fst->AddArc(prevState,StdArc(*itNode,*itNode,0,curState));
		//cout << *itNode << " ";
		prevState=curState;
	}
	fst->SetFinal(curState,0);
	//cout << endl;
	return fst;
}
