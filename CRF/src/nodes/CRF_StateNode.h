#ifndef CRF_STATENODE_H_
#define CRF_STATENODE_H_

#include <string>
#include "fst/fstlib.h"
#include "../CRF.h"
#include "../CRF_Model.h"
#include <vector>

//typedef pair<uint,uint> StatePair;
//typedef pair<uint,uint> LabelPair;
//typedef pair< StatePair, LabelPair > StateNode;
/*
class CRF_ViterbiState {
public:
	uint ilabel;
	uint olabel;
	uint prevstate;
	uint nextstate;
	CRF_ViterbiState(uint ilab, uint olab, uint pstate, uint nstate) :
		ilabel(ilab), olabel(olab), prevstate(pstate), nextstate(nstate) {}
	bool operator == (const CRF_ViterbiState p) const {
		//return (ilabel==p.ilabel && olabel==p.olabel && prevstate==p.prevstate && nextstate==p.nextstate);
		//return (ilabel==p.ilabel && olabel==p.olabel && nextstate==p.nextstate);
		//return (ilabel==p.ilabel && nextstate==p.nextstate);
		//return (nextstate==p.nextstate && prevstate==p.prevstate);
		return (nextstate==p.nextstate);
	}

	bool operator < (const CRF_ViterbiState p) const {
		// definite order - order by prevstate, then by nextstate, then by ilabel, then by olabel
		//if (prevstate==p.prevstate) {
			//if (nextstate==p.nextstate) {
				//if (ilabel == p.ilabel) {
				//	return (olabel < p.olabel);
				//}
				//else {
				//	return (ilabel < p.ilabel);
				//}
			//}
			//else {
				return (nextstate < p.nextstate);
			//}
		//}
		//else {
		//	return (prevstate < p.prevstate);
		//}
	}
};
*/



using namespace fst;

class CRF_StateNode
{
protected:
	float* ftrBuf;
	QNUInt32 ftrBuf_size;
	QNUInt32 label;
	CRF_Model* crf_ptr;
	double* alphaArray;
	double* betaArray;
	double* alphaBetaArray;
	vector<double> alphaArrayAligned;
	vector<double> alphaArrayAlignedBase;
	vector<double> betaArrayAligned;
	vector<double> betaArrayAlignedBase;
	QNUInt32 alphaSize;
	double alphaScale;
	QNUInt32 nLabs;
	double* prevAlpha;
public:
	vector<double> gammaProbs;
	vector<double> gammaTransProbs;
	vector<double> alphaVector;
	vector<double> betaVector;
	vector<double> gammaProbsAligned;
	vector<double> gammaTransProbsAligned;
	vector<int> latticeNodes;
	vector<int> labLatticeNodes;
	vector<const StdArc*> viterbiStates;
	//vector<uint> viterbiStateIds;
	//vector<CRF_ViterbiState> viterbiStateIds;
	vector<uint> viterbiPhnIds;
	//vector<int> viterbiIlabs;
	//vector<int> viterbiOlabs;
	//vector<float> viterbiWeights;
	vector<int> viterbiPointers;
	//double viterbiPruningWeight;
	CRF_StateNode(float* fb, QNUInt32 sizeof_fb, QNUInt32 lab, CRF_Model* crf);
	virtual ~CRF_StateNode();
	virtual double computeTransMatrix();
	virtual double computeTransMatrixLog();
	virtual double computeAlpha(double* prev_alpha);
	virtual double computeFirstAlpha(double* prev_alpha);
	virtual double computeBeta(double* result_beta, double scale=1.0);
	virtual double* computeAlphaBeta(double Zx);
	virtual void setTailBeta();
	virtual double computeExpF(double* ExpF, double* grad, double Zx, double* prev_alpha, QNUInt32 prev_lab);
	virtual double computeSoftExpF(double* ExpF, double* grad, double Zx, double soft_Zx, double* prev_alpha, vector<double>* prevAlphaAligned, bool firstFrame);
	virtual double computeAlphaSum();
	virtual double computeAlphaAlignedSum();
	virtual void reset(float *fb, QNUInt32 sizeof_fb, QNUInt32 lab, CRF_Model* crf_in);
	virtual double* getAlpha();
	virtual double* getPrevAlpha();
	virtual double* getBeta();
	virtual double* getAlphaBeta();
	virtual vector<double>* getAlphaVector();
	virtual vector<double>* getBetaVector();
	virtual vector<double>* getAlphaAligned();
	virtual vector<double>* getBetaAligned();
	virtual vector<double>* getAlphaAlignedBase();
	virtual vector<double>* getBetaAlignedBase();
	virtual QNUInt32 getLabel();
	virtual double getAlphaScale();
	virtual double getTransValue(QNUInt32 prev_lab, QNUInt32 cur_lab);
	virtual double getStateValue(QNUInt32 cur_lab);
	virtual double getStateValue(QNUInt32 cur_lab, QNUInt32 cur_mix);
	virtual double getFullTransValue(QNUInt32 prev_lab, QNUInt32 cur_lab);
	static CRF_StateNode* createStateNode(float* fb, QNUInt32 sizeof_fb, QNUInt32 lab, CRF_Model* crf);
	virtual float *getFtrBuffer();
	virtual QNUInt32 getFtrBufferSize();
};

#endif /*CRF_STATENODE_H_*/
