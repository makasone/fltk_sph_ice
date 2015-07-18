//�V�~�����[�V�����C�^���v�Z���@��I�����鏃�����z�֐�

#ifndef _ICE_SIMU_METHOD_
#define _ICE_SIMU_METHOD_

#include "Ice_JudgeMove.h"

class Ice_SimuMethod
{
public:
	virtual void SetJudgeMove(Ice_JudgeMove* judge) = 0;
	virtual Ice_JudgeMove* GetJudgeMove() = 0;

	virtual void StepObjMove() = 0;
	virtual void StepObjMoveItr() = 0;

	virtual void StepObjUpdate() = 0;
	virtual void StepObjUpdateItr() = 0;
	virtual void StepObjUpdateItrEnd() = 0;

	virtual void StepObjMoveDebug() = 0;
	virtual void StepObjMoveItrDebug() = 0;
};

#endif