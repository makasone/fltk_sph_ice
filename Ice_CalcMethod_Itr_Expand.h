//�v�Z���@��I�����鏃�����z�֐�

#ifndef _ICE_CALC_METHOD_ITR_EXPAND_
#define _ICE_CALC_METHOD_ITR_EXPAND_

#include "Ice_CalcMethod.h"

#include "Ice_ClusterMove.h"
#include "Ice_JudgeMove.h"
#include "Ice_Convolution.h"

#include "IceObject.h"
#include "Ice_SM.h"

using namespace std;

class Ice_CalcMethod_Itr_Expand: public Ice_CalcMethod
{
public:
	Ice_CalcMethod_Itr_Expand(const vector<Ice_SM*>& iceSM, Ice_ClusterMove* clusterMove, Ice_Convolution* convo);
	~Ice_CalcMethod_Itr_Expand();
	
	void SetObjMove(Ice_ClusterMove* clusterMove);
	void SetConvolution(Ice_Convolution* convo);

	void StepObjMove();
	void StepObjMoveDebug();

	void CopyOriginalObject(vector<Ice_SM>& copyObj);
	void DeleteCopyObject(vector<Ice_SM>& copyObj);
	void ReplaceCluster(Ice_SM copyObj[125]);

	void GetExpandeCluster();
	void ExpandeCluster();
	void ContractCluster();

private:
	vector<Ice_SM*> m_iceSM;

	//�^���v�Z���@�������N���X
	Ice_ClusterMove* m_iceMove;

	//�ŏI�������ʂ����߂�N���X
	Ice_Convolution* m_iceConvo;
};

#endif