//計算方法を選択する純粋仮想関数　運動計算の上位

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