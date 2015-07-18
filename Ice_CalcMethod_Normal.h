//計算方法を選択する純粋仮想関数

#ifndef _ICE_CALC_METHOD_NORMAL_
#define _ICE_CALC_METHOD_NORMAL_

#include "Ice_CalcMethod.h"
#include "Ice_SimuMethod.h"
#include "Ice_Convolution.h"

using namespace std;

class Ice_CalcMethod_Normal : public Ice_CalcMethod
{
public:
	Ice_CalcMethod_Normal(Ice_SimuMethod* simuMethod, Ice_Convolution* convo);
	~Ice_CalcMethod_Normal();

	void SetObjMove(Ice_SimuMethod* simuMethod){	m_simu = simuMethod;	}
	void SetConvolution(Ice_Convolution* convo){	m_convo = convo;	}

	void StepObjMove();

	void StepObjMoveDebug();

private:
	//運動計算方法を扱うクラス
	Ice_SimuMethod* m_simu;

	Ice_Convolution* m_convo;
};

#endif