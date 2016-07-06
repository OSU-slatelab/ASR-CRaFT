/*
 * CRF_GradBuilder.cpp
 *
 * Copyright (c) 2010
 * Author: Jeremy Morris
 *
 */

#include "CRF_Model.h"

/*
 * CRF_Model constructor
 *
 * Input: num_labs - number of possible labels in the CRF model
 */
CRF_Model::CRF_Model(QNUInt32 num_labs)
	: nlabs(num_labs)
{
	this->useLog=true;
	this->useMask=false;
	this->init_present=0;
	this->node_type=STD_STATE;

	//Added by Ryan, default values. It should be nlabs = nActualLabs * lab_max_dur.
	this->lab_max_dur = 1;
	this->nActualLabs = num_labs;
	this->model_type = STDFRAME;
	this->use_broken_class_label = false;
	this->init_iter = 0;
}

/*
 * CRF_Model destructor
 */
CRF_Model::~CRF_Model()
{
	if (this->lambda != NULL ) {delete [] this->lambda;}
	if (this->lambdaAcc != NULL ) {delete [] this->lambdaAcc;}
	if (this->featureMap != NULL ) {delete this->featureMap;}

	// added by Ryan
	if (this->gradSqrAcc != NULL ) {delete [] this->gradSqrAcc;}
}

/*
 * CRF_Model::getNLabs
 *
 * Returns: Number of labels defined for the CRF Model
 */
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

/*
 * CRF_Model::setFeatureMap
 *
 * Mutator function to set the feature map for the model.
 * As a side effect, also defines the size of the lambda vector, based on how many feature funcs
 * the feature map defines.
 */
void CRF_Model::setFeatureMap(CRF_FeatureMap* map)
{
	this->featureMap = map;
	QNUInt32 len=map->getNumFtrFuncs();
	this->setLambda(new double[len],len);
	this->setLambdaAcc(new double[len]);

	// Added by Ryan
	this->setGradSqrAcc(new double[len]);

	this->resetLambda();
	this->setNodeType();
}


/*
 * CRF_Mode::getFeatureMap
 *
 * Accessor function to retrieve the model's feature map
 */
CRF_FeatureMap* CRF_Model::getFeatureMap()
{
	return this->featureMap;
}

/*
 * CRF_Model::getLambda
 *
 * Accessor function to retrieve the model's lambda vector
 */
double* CRF_Model::getLambda()
{
	return this->lambda;
}

/*
 * CRF_Model::getLambdaAcc
 *
 * Accessor function to retrieve the model's lambda accumulator vector
 *
 * The accumulator vector is used to create average weights in stochastic gradient descent processing
 * (or other training methods that use weight averaging)
 */
double* CRF_Model::getLambdaAcc()
{
	return this->lambdaAcc;
}

/*
 * CRF_Model::getLambdaLen
 *
 * Accessor function to retrieve the length of the lambda vector
 */
QNUInt32 CRF_Model::getLambdaLen()
{
	return this->lambda_len;
}

// Added by Ryan
/*
 * CRF_Model::getGradSqrAcc
 *
 * Accessor function to retrieve the model's gradient square sum accumulator vector
 *
 * The accumulator vector is used for AdaGrad in stochastic gradient descent processing
 */
double* CRF_Model::getGradSqrAcc()
{
	return this->gradSqrAcc;
}

/*
 * CRF_Model::setLambda
 *
 * Mutator function to set the lambda vector and the length of the lambda vector
 */
void CRF_Model::setLambda(double* lam, QNUInt32 lam_len)
{
	this->lambda = lam;
	this->lambda_len = lam_len;
}

/*
 * CRF_Model::setLambdaAcc
 *
 * Mutator function to set the lambda accumulator vector
 * (Note that this vector will be treated as the same length as the lambda vector)
 */
void CRF_Model::setLambdaAcc(double* lam)
{
	this->lambdaAcc = lam;
}

// Added by Ryan
/*
 * CRF_Model::setGradSqrAcc
 *
 * Mutator function to set the accumulated gradient square sum vector
 * (Note that this vector will be treated as the same length as the lambda vector)
 */
void CRF_Model::setGradSqrAcc(double* grad)
{
	this->gradSqrAcc = grad;
}

/*
 * CRF_Model::resetLambda
 *
 * Mutator function to reset the lambda vector to zero.
 */
void CRF_Model::resetLambda()
{
	for (QNUInt32 i = 0; i < this->lambda_len; i++) {
		this->lambda[i] = 0.0;
		this->lambdaAcc[i] = 0.0;

		// Added by Ryan
		this->gradSqrAcc[i] = 0.0;
	}
}

/*
 * CRF_Model::writeToFile
 *
 * Inputs: fname - filename to write file to
 *
 * Returns: true if file can be written, false otherwise
 *
 * Writes the current lambda vector to the file fname
 */
bool CRF_Model::writeToFile(const char* fname)
{
	std::ofstream ofile;
	ofile.open(fname);
	if (ofile.is_open()) {
		for (QNUInt32 i=0; i<this->lambda_len; i++) {
			ofile << lambda[i] << std::endl;
		}

		// Added by Ryan
		if (ofile.bad() || ofile.fail()) {
			string errstr="CRF_Model::writeToFile() caught exception: errors when writing the weights to the file:\n";
			errstr += string(fname) + "\n";
			errstr += "Maybe because it is running out of space, or the file is deleted, or the permission is changed.";
			throw runtime_error(errstr);
		}

		ofile.close();
		return true;
	}
	else {

		// Added by Ryan
		string errstr="CRF_Model::writeToFile() caught exception: cannot open the file:\n";
		errstr += string(fname) + "\n";
		errstr += "Maybe because it is running out of space, or the file is deleted, or the permission is changed.";
		throw runtime_error(errstr);

		return false;
	}
}

/*
 * CRF_Model::writeToFile
 *
 * Inputs: fname - filename to write file to
 *         lam - lambda vector to write to file
 *         ll - length of lambda vector lam
 *
 * Returns: true if file can be written, false otherwise
 *
 * Writes the arbitrary lambda vector lam to the file fname
 */

bool CRF_Model::writeToFile(const char* fname, double* lam, QNUInt32 ll)
{
	std::ofstream ofile;
	ofile.open(fname);
	if (ofile.is_open()) {
		for (QNUInt32 i=0; i<ll; i++) {
			ofile << lam[i] << std::endl;
		}

		// Added by Ryan
		if (ofile.bad()) {
			string errstr="CRF_Model::writeToFile() caught exception: errors when writing the weights to the file:\n";
			errstr += string(fname) + "\n";
			errstr += "Maybe because it is running out of space, or the file is deleted, or the permission is changed.";
			throw runtime_error(errstr);
		}

		ofile.close();
		return true;
	}
	else {

		// Added by Ryan
		string errstr="CRF_Model::writeToFile() caught exception: cannot open the file:\n";
		errstr += string(fname) + "\n";
		errstr += "Maybe because it is running out of space, or the file is deleted, or the permission is changed.";
		throw runtime_error(errstr);

		return false;
	}
}

/*
 * CRF_Model::readFromFile
 *
 * Inputs: fname - file name to read lambda values from
 *
 * Returns: true if file can be read, false otherwise
 *
 * Reads values from file fname into lambda vector
 */
bool CRF_Model::readFromFile(const char* fname)
{
	std::ifstream ifile;
	ifile.open(fname);
	if (ifile.is_open()) {
		for (QNUInt32 i = 0; i < this->lambda_len; i++) {
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

/*
 * CRF_Model::readAverageFromFile
 *
 * Inputs: fname - file name to read weights from
 *         present - number of presentations used for averaging
 *
 * Returns: true if file can be read, false otherwise
 *
 * Reads the weights from file fname, multiplies them by the number of presentations present and
 * stores them in the lambda accumulator
 *
 * Used to start training from a previous set of accumulated weight values
 */
bool CRF_Model::readAverageFromFile(const char* fname, int present)
{
	this->init_present = present;
	std::ifstream ifile;
	ifile.open(fname);
	if (ifile.is_open()) {
		for (QNUInt32 i = 0; i < this->lambda_len; i++) {
			std::string s;
			getline(ifile,s);
			std::istringstream iss(s);
			iss >> std::dec >> lambdaAcc[i];
			if (present > 0) {
				lambdaAcc[i] = lambdaAcc[i] * present;

				// added by Ryan, just for debugging
				//std::cout << lambdaAcc[i] << std::endl;
			}
		}
		ifile.close();
		return true;
	}
	else {
		return false;
	}
}

// Added by Ryan
/*
 * CRF_Model::readGradSqrAccFromFile
 *
 * Inputs: fname - file name to read lambda values from
 *
 * Returns: true if file can be read, false otherwise
 *
 * Reads values from file fname into the accumulated square sum vector
 */
bool CRF_Model::readGradSqrAccFromFile(const char* fname)
{
	std::ifstream ifile;
	ifile.open(fname);
	if (ifile.is_open()) {
		for (QNUInt32 i = 0; i < this->lambda_len; i++) {
			std::string s;
			getline(ifile,s);
			std::istringstream iss(s);
			iss >> std::dec >> this->gradSqrAcc[i];
			//std::cout << this->gradSqrAcc[i] << std::endl;
		}
		ifile.close();
		return true;
	}
	else {
		return false;
	}
}

/*
 * CRF_Model::getPresentations
 *
 * Accessor function to the number of presentations used to train the lambda vector
 */
QNUInt32 CRF_Model::getPresentations()
{
	return this->init_present;
}

// Added by Ryan
/*
 * CRF_Model::setLabMaxDur
 *
 * Mutator function to set the maximum duration of labels
 *
 */
void CRF_Model::setLabMaxDur(QNUInt32 max_duration) {
	if (max_duration == 0)
	{
		string errstr="CRF_Model::setLabMaxDur() caught exception: the maximum duration of labels must be larger than 0.";
		throw runtime_error(errstr);
	}
	this->lab_max_dur=max_duration;
}

// Added by Ryan
/*
 * CRF_Model::getLabMaxDur
 *
 * Returns: The maximum duration of labels defined for the CRF Model
 */
QNUInt32 CRF_Model::getLabMaxDur()
{
	return this->lab_max_dur;
}

// Added by Ryan
/*
 * CRF_Model::setNActualLabs
 *
 * Mutator function to set the number of actual labels (without duration)
 *
 */
void CRF_Model::setNActualLabs(QNUInt32 num_actual_labs)
{
	this->nActualLabs = num_actual_labs;
}

// Added by Ryan
/*
 * CRF_Model::getNActualLabs
 *
 * Returns: the number of actual labels (without duration)
 *
 */
QNUInt32 CRF_Model::getNActualLabs()
{
	return this->nActualLabs;
}

// Added by Ryan
/*
 * CRF_Model::setModelType
 *
 * Mutator function to set the model type
 *
 */
void CRF_Model::setModelType(modeltype mtype)
{
	this->model_type = mtype;
}

// Added by Ryan
/*
 * CRF_Model::getModelType
 *
 * Accessor function to get the model type
 */
modeltype CRF_Model::getModelType()
{
	return this->model_type;
}

// Added by Ryan
/*
 * CRF_Model::setBrokenClassLabel
 *
 * Mutator function to set whether to use broken-class label
 *
 */
void CRF_Model::setBrokenClassLabel(bool useBrokenClass) {
	this->use_broken_class_label = useBrokenClass;
}

// Added by Ryan
/*
 * CRF_Model::ifUseBrokenClassLabel
 *
 * Accessor function to get whether to use broken-class label
 */
bool CRF_Model::ifUseBrokenClassLabel() {
	return this->use_broken_class_label;
}


// Added by Ryan
/*
 * CRF_Model::setInitIter
 *
 * Mutator function to set the starting training iteration.
 * It is non-zero only when the initial weight file is set.
 *
 */
void CRF_Model::setInitIter(QNUInt32 start_iter)
{
	this->init_iter = start_iter;
}

// Added by Ryan
/*
 * CRF_Model::getInitIter
 *
 * Returns: the starting training iteration.
 * It is non-zero only when the initial weight file is set.
 */
QNUInt32 CRF_Model::getInitIter()
{
	return this->init_iter;
}

/*
 * CRF_Model::setUseLog
 *
 * Mutator function to set the use logspace flag
 *
 */
void CRF_Model::setUseLog(bool isLog) {
	this->useLog=isLog;
	this->setNodeType();
}

/*
 * CRF_Model::setuseMask
 *
 * Mutator function to set the use mask flag
 */
void CRF_Model::setUseMask(bool isMasked) {
	this->useMask=isMasked;
	this->setNodeType();
}

/*
 * CRF_Model::setNodeType
 *
 * Mutator function to set the nodetype
 */
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

/*
 * CRF_Model::getNodeType
 *
 * Accessor function to get the node type
 */
nodetype CRF_Model::getNodeType() {
	return this->node_type;
}
