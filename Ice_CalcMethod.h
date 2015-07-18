//�v�Z���@��I�����鏃�����z�֐��@�^���v�Z�̏��

#ifndef _ICE_CALC_METHOD_
#define _ICE_CALC_METHOD_

#include "Ice_SimuMethod.h"

class Ice_CalcMethod
{
public:
	virtual void SetObjMove(Ice_SimuMethod* simuMethod) = 0;

	virtual void StepObjMove() = 0;
	virtual void StepObjMoveDebug() = 0;
};

#endif