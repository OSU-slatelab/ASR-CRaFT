#ifndef CRFMODEL_H_
#define CRFMODEL_H_

#include "CRF.h"
#include "ftrmaps/CRF_FeatureMap.h"

class CRF_Model
{
protected:
	QNUInt32 nlabs;
	double* lambda;  // Holds the lambda weight values for the CRF
	double* lambdaAcc;  // Holds the lambda accumlation values for averaging
	QNUInt32 lambda_len; // Holds the size of the lambda vector
	//CRF_Range* rangeMachine; // Holds the range generating machine for the CRF
	CRF_FeatureMap* featureMap;
	QNUInt32 init_present;
	nodetype node_type;
public:
	CRF_Model(QNUInt32 num_labs);
	bool useLog;
	bool useMask;
	virtual ~CRF_Model();
	QNUInt32 getNLabs();
	//virtual void setRangeMachine(CRF_Range*);
	//virtual CRF_Range* getRangeMachine();
	virtual void setFeatureMap(CRF_FeatureMap*);
	virtual CRF_FeatureMap* getFeatureMap();
	virtual double* getLambda();
	virtual QNUInt32 getLambdaLen();
	virtual QNUInt32 getPresentations();
	virtual double* getLambdaAcc();
	virtual void setLambda(double*,QNUInt32);
	virtual void setLambdaAcc(double*);
	virtual void resetLambda();
	virtual bool writeToFile(const char* fname);
	virtual bool writeToFile(const char* fname,double* lam, QNUInt32 ll);
	virtual bool readFromFile(const char* fname);
	virtual bool readAverageFromFile(const char* fname, int present);
	virtual void setUseLog(bool isLog);
	virtual void setUseMask(bool isMasked);
	virtual void setNodeType();
	virtual nodetype getNodeType();
};

#endif /*CRF_H_*/
