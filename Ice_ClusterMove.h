//運動計算で用いる手法を選択する純粋仮想関数

#ifndef _ICE_CLUSTER_MOVE_
#define _ICE_CLUSTER_MOVE_

#include "Ice_JudgeMove.h"

class Ice_ClusterMove
{
public:
	virtual void SetJudgeMove(Ice_JudgeMove* judge) = 0;
	virtual Ice_JudgeMove* GetJudgeMove() = 0;

	virtual void StepObjMove() = 0;
	virtual void StepObjMoveItr() = 0;

	virtual void StepObjMoveDebug() = 0;
};

#endif