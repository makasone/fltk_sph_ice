//�v�Z���@��I�����鏃�����z�֐�

#ifndef _ICE_CALC_METHOD_NORMAL_
#define _ICE_CALC_METHOD_NORMAL_

#include "Ice_CalcMethod.h"

#include "Ice_ClusterMove.h"

using namespace std;

class Ice_CalcMethod_Normal : public Ice_CalcMethod
{
public:
	Ice_CalcMethod_Normal(Ice_ClusterMove* clusterMove);
	~Ice_CalcMethod_Normal();

	void SetObjMove(Ice_ClusterMove* clusterMove);

	void StepObjMove();
	void StepObjMoveDebug();

private:
	//�^���v�Z���@�������N���X
	Ice_ClusterMove* m_iceMove;
};

#endif