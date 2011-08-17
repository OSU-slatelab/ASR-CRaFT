/*
 * CRF_InLabStream_RandPresent.cpp
 *
 * Copyright (c) 2010
 * Author: Jeremy Morris
 *
 */

#include "CRF_MLFManager.h"
typedef StdArc::StateId StateId;
using namespace fst;

/*
 * CRF_MLFManager constructor
 *
 * Input: mlffile - text file containing MLF formatted utterances
 *        olist - text file containing output tokens in OpenFst format
 *        symtab - OpenFst symbol table of output tokens
 *
 */
CRF_MLFManager::CRF_MLFManager(char* mlffile, char* olist, SymbolTable* symTab)
	: symTab(symTab){
	this->readMLF(mlffile);
}

/*
 * CRF_MLFManager destructor
 */
CRF_MLFManager::~CRF_MLFManager() {
	//cout << "in CRF_MLFManager destructor" << endl;
}

/*
 * CRF_MLFManager::getKey
 *
 * Input: fname - string containing the file name of the utterance being searched for (e.g. "filename.lab")
 *
 * Returns: the key portion of the filename (e.g. "filename.lab" -> "filename")
 */
string CRF_MLFManager::getKey(string fname) {

	int lastc=fname.rfind(".");
	int firstc=fname.rfind("/");
	return(fname.substr(firstc+1,lastc-firstc-1));
}

/*
 * CRF_MLFManager::readMLF
 *
 * Input: mlffile - string containing the name of the input mlffile
 *
 * Reads the MLF file and stores its contents in the transcripts vector.
 */
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

/*
 * CRF_MLFManager::getFst
 *
 * Input: fname - the filename of the sequence being requested
 *
 * Returns: OpenFst formated VectorFst of the transcript from the MLF
 *
 */
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
