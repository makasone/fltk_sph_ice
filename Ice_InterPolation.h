//補間で用いる手法を選択する純粋仮想関数

#ifndef _ICE_INTER_POLATION_
#define _ICE_INTER_POLATION_

#include "Ice_InterPolationJudge.h"

class Ice_InterPolation
{
public:
	virtual void StepInterPolation() = 0;
	virtual void StepInterPolationItr() = 0;

	virtual void SetIntrpJudge(Ice_InterPolationJudge* judge) = 0;
	virtual Ice_InterPolationJudge* GetIntrpJudge() = 0;
};

#endif