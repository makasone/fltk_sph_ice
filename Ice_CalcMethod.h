//計算方法を選択する純粋仮想関数

#ifndef _ICE_CALC_METHOD_
#define _ICE_CALC_METHOD_

#include "Ice_ClusterMove.h"

class Ice_CalcMethod
{
public:
	virtual void SetObjMove(Ice_ClusterMove* clusterMove) = 0;

	virtual void StepObjMove() = 0;
};

#endif