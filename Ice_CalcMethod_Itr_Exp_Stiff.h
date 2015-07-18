//�N���X�^�g�唽���{臒l�w��

#ifndef _ICE_CALC_METHOD_ITR_EXP_STIFF_
#define _ICE_CALC_METHOD_ITR_EXP_STIFF_

#include "Ice_CalcMethod.h"

#include "Ice_ClusterMove.h"
#include "Ice_JudgeMove.h"
#include "Ice_Convolution.h"
#include "Ice_CalcStiffData.h"
#include "Ice_CalcStiffData_Summation.h"
#include "Ice_CalcStiffData_Average.h"
#include "Ice_CalcStiffData_StdDevision.h"
#include "Ice_CalcStiffData_CompareRigid.h"

#include "IceObject.h"
#include "Ice_SM.h"

using namespace std;

class Ice_CalcMethod_Itr_Exp_Stiff : public Ice_CalcMethod
{
public:
	//Ice_CalcMethod_Itr_Exp_Stiff(const vector<Ice_SM*>& iceSM, Ice_ClusterMove* clusterMove, Ice_Convolution* convo);
	Ice_CalcMethod_Itr_Exp_Stiff(const vector<Ice_SM*>& iceSM, Ice_SimuMethod* simuMethod, Ice_Convolution* convo);
	~Ice_CalcMethod_Itr_Exp_Stiff();
	
	//void SetObjMove(Ice_ClusterMove* clusterMove);
	void SetObjMove(Ice_SimuMethod* simuMethid);
	void SetConvolution(Ice_Convolution* convo);

	void StepObjMove();
	void StepObjMoveDebug();

	void CalcVel();

	void ExpandeCluster(vector<int>& searchFinishIndxes);

	void CopyOriginalObject(vector<vector<unsigned>>& copyIndxes);
	void ReplaceCluster(const vector<vector<unsigned>>& copyIndxes);

	void DebugStiffness();

private:
	vector<Ice_SM*> m_iceSM;

	//�^���v�Z���@�������N���X
	//Ice_ClusterMove* m_iceMove;
	Ice_SimuMethod* m_simuMethod;

	//�ŏI�������ʂ����߂�N���X
	Ice_Convolution* m_iceConvo;

	//臒l�Ɣ�r����l���v�����邽�߂̃N���X
	Ice_CalcStiffData* m_iceCalcStiff;
};

#endif