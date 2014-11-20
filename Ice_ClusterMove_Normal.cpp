#include "Ice_ClusterMove_Normal.h"

typedef Ice_ClusterMove_Normal MoveNormal;

MoveNormal::Ice_ClusterMove_Normal(const vector<Ice_SM*>& iceSM)
{
	m_iceSM = iceSM;
}

MoveNormal::~Ice_ClusterMove_Normal()
{
}

void MoveNormal::SetJudgeMove(Ice_JudgeMove* judge)
{
	m_iceJudge = judge;	
}

Ice_JudgeMove* MoveNormal::GetJudgeMove()
{
	return m_iceJudge;
}

void MoveNormal::StepObjMove()
{
	#pragma omp parallel for
	for(int i = 0; i < IceObject::GetParticleNum(); ++i)
	{	
		if(m_iceJudge->JudgeMove(i) == false){	continue;	}

		m_iceSM[i]->UpdateCPU();				//運動計算
	}
}

void MoveNormal::StepObjMoveItr()
{
	#pragma omp parallel for
	for(int i = 0; i < IceObject::GetParticleNum(); ++i)
	{	
		if(m_iceJudge->JudgeMove(i) == false){	continue;	}

		m_iceSM[i]->ShapeMatchingIteration();		//現在の粒子位置を用いてSM法
	}	
}