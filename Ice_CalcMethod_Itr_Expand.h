//�v�Z���@��I�����鏃�����z�֐�

#ifndef _ICE_CALC_METHOD_ITR_EXPAND_
#define _ICE_CALC_METHOD_ITR_EXPAND_

#include "Ice_CalcMethod.h"

#include "Ice_ClusterMove.h"
#include "Ice_JudgeMove.h"
#include "Ice_Convolution.h"

#include "QueryCounter.h"
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

	void StepObjMoveTest();

	void CopyOriginalObject(vector<vector<unsigned>>& copyIndxes);
	void ReplaceCluster(const vector<vector<unsigned>>& copyIndxes);

	void GetExpandeCluster();
	void ExpandeCluster(vector<int>& searchFinishIndxes);
	void ExpandeCluster_Test(vector<int>& searchFinishIndxes);
	void ExpandCluster_Test2(vector<int>& searchFinishIndxes);

	void ContractCluster();

	void CalcVel();

private:
	vector<Ice_SM*> m_iceSM;

	//�^���v�Z���@�������N���X
	Ice_ClusterMove* m_iceMove;

	//�ŏI�������ʂ����߂�N���X
	Ice_Convolution* m_iceConvo;
};

#endif