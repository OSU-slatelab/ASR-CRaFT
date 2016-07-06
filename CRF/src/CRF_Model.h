#ifndef CRFMODEL_H_
#define CRFMODEL_H_
/*
 * CRF_Model.h
 *
 * Copyright (c) 2010
 * Author: Jeremy Morris
 *
 * Implements the CRF_Model class
 */

#include "CRF.h"
#include "ftrmaps/CRF_FeatureMap.h"

/*
 * class CRF_Model
 *
 * Holds the relevant information about the underlying CRF model, including the topology of the model
 * (as dictated by the node type), the interaction of the features with the model (as defined by the
 * feature map), and the current parameters of the model (as defined by the lambda vector).  Most
 * of the functions in this class are accessor and mutator functions for these elements.
 */
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

	//Added by Ryan
	QNUInt32 lab_max_dur;
	// TODO: set up different nActualLabs for different kinds of nodes.
	QNUInt32 nActualLabs;
	modeltype model_type;
	bool use_broken_class_label;
	QNUInt32 init_iter;
	double* gradSqrAcc;  // Holds the accumulated gradient square sums for AdaGrad

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
	virtual bool writeToFile(const char* fname, double* lam, QNUInt32 ll);
	virtual bool readFromFile(const char* fname);
	virtual bool readAverageFromFile(const char* fname, int present);
	virtual void setUseLog(bool isLog);
	virtual void setUseMask(bool isMasked);
	virtual void setNodeType();
	virtual nodetype getNodeType();

	//Added by Ryan
	virtual void setLabMaxDur(QNUInt32 max_duration);
	virtual QNUInt32 getLabMaxDur();
	virtual void setNActualLabs(QNUInt32 num_actual_labs);
	virtual QNUInt32 getNActualLabs();
	virtual void setModelType(modeltype mtype);
	virtual modeltype getModelType();
	virtual void setBrokenClassLabel(bool useBrokenClass);
	virtual bool ifUseBrokenClassLabel();
	virtual void setInitIter(QNUInt32 start_iter);
	virtual QNUInt32 getInitIter();
	virtual bool readGradSqrAccFromFile(const char* fname);
	virtual void setGradSqrAcc(double* grad);
	virtual double* getGradSqrAcc();
};

#endif /*CRF_H_*/
