#include "Ice_CalcMethod_Normal.h"

typedef Ice_CalcMethod_Normal CalcNormal;

CalcNormal::Ice_CalcMethod_Normal(Ice_SimuMethod* simuMethod, Ice_Convolution* convo)
{
	SetObjMove(simuMethod);
	SetConvolution(convo);
}

CalcNormal::~Ice_CalcMethod_Normal()
{
}

void CalcNormal::StepObjMove()
{
	m_simu->StepObjMove();
	m_convo->StepConvolution();
	m_simu->StepObjUpdate();
}

void CalcNormal::StepObjMoveDebug()
{
	m_simu->StepObjMoveDebug();
}