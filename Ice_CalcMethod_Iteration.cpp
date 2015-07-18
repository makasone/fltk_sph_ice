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
	////����
	//m_iceSimu->StepObjMove();

	////����
	//for(int itr = 1; itr < Ice_SM::GetIteration(); itr++)
	//{
	//	m_iceConvo->StepConvolution();
	//	m_iceSimu->StepObjMoveItr();
	//}

	////���x�Z�o
	//#pragma omp parallel for
	//for(int i = 0; i < IceObject::GetParticleNum(); ++i)
	//{	
	//	if(m_iceSimu->GetJudgeMove()->JudgeMove(i) == false){ continue;	}
	//	m_iceSM[i]->integrateIteration();
	//}

	//����
	m_iceSimu->StepObjMove();
	m_iceConvo->StepConvolution();
	m_iceSimu->StepObjUpdateItr();

	//����
	for(int itr = 1; itr < Ice_SM::GetIteration(); itr++){
		m_iceSimu->StepObjMoveItr();
		m_iceConvo->StepConvolution();
		m_iceSimu->StepObjUpdateItr();
	}

	//���q�ʒu�C���x�C�p���x�C�p�����X�V
	m_iceSimu->StepObjUpdateItrEnd();
}

void CalcIteration::StepObjMoveDebug()
{
	////����
	//m_iceSimu->StepObjMoveDebug();

	////����
	//for(int itr = 1; itr < Ice_SM::GetIteration(); itr++)
	//{
	//	m_iceConvo->StepConvolution();
	//	m_iceSimu->StepObjMoveItrDebug();
	//}

	////���x�Z�o
	//#pragma omp parallel for
	//for(int i = 0; i < IceObject::GetParticleNum(); ++i)
	//{	
	//	if(m_iceSimu->GetJudgeMove()->JudgeMoveDebug(i) == false){ continue;	}
	//	m_iceSM[i]->integrateIteration();
	//}
}