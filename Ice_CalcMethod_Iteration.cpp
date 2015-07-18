#include "Ice_CalcMethod_Iteration.h"

typedef Ice_CalcMethod_Iteration CalcIteration;

CalcIteration::Ice_CalcMethod_Iteration(const vector<Ice_SM*>& iceSM, Ice_SimuMethod* simuMethod, Ice_Convolution* convo)
{
	m_iceSM = iceSM;
	SetObjMove(simuMethod);
	SetConvolution(convo);
}

CalcIteration::Ice_CalcMethod_Iteration(Ice_SimuMethod* simuMethod, Ice_Convolution* convo)
{
	SetObjMove(simuMethod);
	SetConvolution(convo);
}

CalcIteration::~Ice_CalcMethod_Iteration()
{
}

void CalcIteration::SetObjMove(Ice_SimuMethod* simuMethod)
{
	m_iceSimu = simuMethod;
}

void CalcIteration::SetConvolution(Ice_Convolution* convo)
{
	m_iceConvo = convo;
}

void CalcIteration::StepObjMove()
{
	////初回
	//m_iceSimu->StepObjMove();

	////反復
	//for(int itr = 1; itr < Ice_SM::GetIteration(); itr++)
	//{
	//	m_iceConvo->StepConvolution();
	//	m_iceSimu->StepObjMoveItr();
	//}

	////速度算出
	//#pragma omp parallel for
	//for(int i = 0; i < IceObject::GetParticleNum(); ++i)
	//{	
	//	if(m_iceSimu->GetJudgeMove()->JudgeMove(i) == false){ continue;	}
	//	m_iceSM[i]->integrateIteration();
	//}

	//初回
	m_iceSimu->StepObjMove();
	m_iceConvo->StepConvolution();
	m_iceSimu->StepObjUpdateItr();

	//反復
	for(int itr = 1; itr < Ice_SM::GetIteration(); itr++){
		m_iceSimu->StepObjMoveItr();
		m_iceConvo->StepConvolution();
		m_iceSimu->StepObjUpdateItr();
	}

	//粒子位置，速度，角速度，姿勢を更新
	m_iceSimu->StepObjUpdateItrEnd();
}

void CalcIteration::StepObjMoveDebug()
{
	////初回
	//m_iceSimu->StepObjMoveDebug();

	////反復
	//for(int itr = 1; itr < Ice_SM::GetIteration(); itr++)
	//{
	//	m_iceConvo->StepConvolution();
	//	m_iceSimu->StepObjMoveItrDebug();
	//}

	////速度算出
	//#pragma omp parallel for
	//for(int i = 0; i < IceObject::GetParticleNum(); ++i)
	//{	
	//	if(m_iceSimu->GetJudgeMove()->JudgeMoveDebug(i) == false){ continue;	}
	//	m_iceSM[i]->integrateIteration();
	//}
}