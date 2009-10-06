#include "CRF_Model.h"



CRF_Model::CRF_Model(QNUInt32 num_labs)
	: nlabs(num_labs)
{
	this->useLog=true;
	this->useMask=false;
	this->init_present=0;
	this->node_type=STD_STATE;
}

CRF_Model::~CRF_Model()
{
	if (this->lambda != NULL ) {delete [] this->lambda;}
	if (this->lambdaAcc != NULL ) {delete [] this->lambdaAcc;}
	if (this->featureMap != NULL ) {delete this->featureMap;}
}

QNUInt32 CRF_Model::getNLabs()
{
	return this->nlabs;
}

//void CRF_Model::setRangeMachine(CRF_Range* rng)
//{
//	this->rangeMachine = rng;
//	QNUInt32 len=rng->getNumFtrs();
//	this->setLambda(new double[len],len);
//	//this->lambda = new float[len];
//}

//CRF_Range* CRF_Model::getRangeMachine()
//{
//	return this->rangeMachine;
//}

void CRF_Model::setFeatureMap(CRF_FeatureMap* map)
{
	this->featureMap = map;
	QNUInt32 len=map->getNumFtrFuncs();
	this->setLambda(new double[len],len);
	this->setLambdaAcc(new double[len]);
	this->resetLambda();
	this->setNodeType();
}



CRF_FeatureMap* CRF_Model::getFeatureMap()
{
	return this->featureMap;
}

double* CRF_Model::getLambda()
{
	return this->lambda;
}

double* CRF_Model::getLambdaAcc()
{
	return this->lambdaAcc;
}

QNUInt32 CRF_Model::getLambdaLen()
{
	return this->lambda_len;
}

void CRF_Model::setLambda(double* lam, QNUInt32 lam_len)
{
	this->lambda=lam;
	this->lambda_len=lam_len;
}

void CRF_Model::setLambdaAcc(double* lam)
{
	this->lambdaAcc=lam;
}

void CRF_Model::resetLambda()
{
	for (QNUInt32 i=0; i<this->lambda_len; i++) {
		this->lambda[i]=0.0;
		this->lambdaAcc[i]=0.0;
	}
}


bool CRF_Model::writeToFile(const char* fname)
{
	std::ofstream ofile;
	ofile.open(fname);
	if (ofile.is_open()) {
		for (QNUInt32 i=0; i<this->lambda_len; i++) {
			ofile << lambda[i] << std::endl;
		}
		ofile.close();
		return true;
	}
	else {
		return false;
	}
}

bool CRF_Model::writeToFile(const char* fname, double* lam, QNUInt32 ll)
{
	std::ofstream ofile;
	ofile.open(fname);
	if (ofile.is_open()) {
		for (QNUInt32 i=0; i<ll; i++) {
			ofile << lam[i] << std::endl;
		}
		ofile.close();
		return true;
	}
	else {
		return false;
	}
}

bool CRF_Model::readFromFile(const char* fname)
{
	std::ifstream ifile;
	ifile.open(fname);
	if (ifile.is_open()) {
		for (QNUInt32 i=0; i<this->lambda_len; i++) {
			std::string s;
			getline(ifile,s);
			std::istringstream iss(s);
			iss >> std::dec >> lambda[i];
		//std::cout << lambda[i] << std::endl;
		}
		ifile.close();
		return true;
	}
	else {
		return false;
	}
}

bool CRF_Model::readAverageFromFile(const char* fname, int present)
{
	this->init_present=present;
	std::ifstream ifile;
	ifile.open(fname);
	if (ifile.is_open()) {
		for (QNUInt32 i=0; i<this->lambda_len; i++) {
			std::string s;
			getline(ifile,s);
			std::istringstream iss(s);
			iss >> std::dec >> lambdaAcc[i];
			if (present>0) {
				lambdaAcc[i]=lambdaAcc[i]*present;
			}
		}
		ifile.close();
		return true;
	}
	else {
		return false;
	}
}

QNUInt32 CRF_Model::getPresentations()
{
	return this->init_present;
}

void CRF_Model::setUseLog(bool isLog) {
	this->useLog=isLog;
	this->setNodeType();
}

void CRF_Model::setUseMask(bool isMasked) {
	this->useMask=isMasked;
	this->setNodeType();
}

void CRF_Model::setNodeType() {
	if (this->useLog) {
		if (this->getFeatureMap()->getNumStates()>1) {
			if (this->useMask) {
				this->node_type=STD_NSTATELOGMASKED;
			}
			else {
				this->node_type=STD_NSTATELOG;
			}
		}
		else {
			if (this->useMask) {
				this->node_type=STD_STATELOGMASKED;
			}
			else {
				this->node_type=STD_STATELOG;
			}
		}
	}
	else {
		if (this->getFeatureMap()->getNumStates()>1) {
			if (this->useMask) {
				this->node_type=STD_NSTATEMASKED;
			}
			else {
				this->node_type=STD_NSTATE;
			}
		}
		else {
			if (this->useMask) {
				this->node_type=STD_STATEMASKED;
			}
			else {
				this->node_type=STD_STATE;
			}
		}
	}

}

nodetype CRF_Model::getNodeType() {
	return this->node_type;
}
