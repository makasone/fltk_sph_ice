//�v�Z���@��I�����鏃�����z�֐�

#ifndef _ICE_CALC_METHOD_ITERATION_
#define _ICE_CALC_METHOD_ITERATION_

#include "Ice_CalcMethod.h"

#include "Ice_ClusterMove.h"
#include "Ice_JudgeMove.h"
#include "Ice_Convolution.h"

#include "IceObject.h"
#include "Ice_SM.h"

using namespace std;

class Ice_CalcMethod_Iteration : public Ice_CalcMethod
{
public:
	Ice_CalcMethod_Iteration(const vector<Ice_SM*>& iceSM, Ice_ClusterMove* clusterMove, Ice_Convolution* convo);
	~Ice_CalcMethod_Iteration();
	
	void SetObjMove(Ice_ClusterMove* clusterMove);
	void SetConvolution(Ice_Convolution* convo);

	void StepObjMove();
	void StepObjMoveDebug();

private:
	vector<Ice_SM*> m_iceSM;

	//�^���v�Z���@�������N���X
	Ice_ClusterMove* m_iceMove;

	//�ŏI�������ʂ����߂�N���X
	Ice_Convolution* m_iceConvo;
};

#endif