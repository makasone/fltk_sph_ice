#include "Ice_SimuMethod_ShapeMatching.h"

typedef Ice_SimuMethod_SM MoveSM;


MoveSM::Ice_SimuMethod_SM(const vector<Ice_SM*>& iceSM)
{
	m_iceSM = iceSM;
}

MoveSM::Ice_SimuMethod_SM(const vector<ElasticObj*>& elasticObj, const vector<OrientedParticle*>& particles)
{
	m_elasticObj = elasticObj;
	m_vOrientedPrtes = particles;
}

MoveSM::~Ice_SimuMethod_SM()
{
}

void MoveSM::SetJudgeMove(Ice_JudgeMove* judge)
{
	m_iceJudge = judge;	
}

Ice_JudgeMove* MoveSM::GetJudgeMove()
{
	return m_iceJudge;
}

void MoveSM::StepObjMove()
{
	//マウスによるドラッグを反映させるために，無理やり値を更新
	OrientedParticleBaseElasticObject::CopyPrtToClstrPos(IceObject::GetParticleNum());

	//運動計算
	#pragma omp parallel for
	for(int i = 0; i < IceObject::GetParticleNum(); i++){	
		if(! m_iceJudge->JudgeMove(i)){	continue;	}

		m_elasticObj[i]->UpdateCluster_SM();
	}
}

void MoveSM::StepObjMoveItr()
{
	#pragma omp parallel for
	for(int i = 0; i < IceObject::GetParticleNum(); ++i){	
		if(! m_iceJudge->JudgeMove(i)){	continue;	}

		//m_iceSM[i]->ShapeMatchingIteration();
		m_elasticObj[i]->UpdateCluster_SM_Itr();
	}
}

//粒子情報の更新
void MoveSM::StepObjUpdate()
{
	//速度，角速度，姿勢を更新
	#pragma omp parallel for
	for(int i = 0; i < IceObject::GetParticleNum(); i++){
		m_vOrientedPrtes[i]->Update();
	}
}

void MoveSM::StepObjUpdateItr()
{

}

void MoveSM::StepObjUpdateItrEnd()
{
	////粒子位置，速度，角速度，姿勢を更新
	//#pragma omp parallel for
	//for(int i = 0; i < IceObject::GetParticleNum(); i++){
	//	//if(! m_iceJudge->JudgeMove(i)){	continue;	}

	//	m_elasticObj[i]->UpdateCluster_SM_Itr_End();
	//}

	//速度，角速度，姿勢を更新
	#pragma omp parallel for
	for(int i = 0; i < IceObject::GetParticleNum(); i++){
		m_vOrientedPrtes[i]->UpdatePosAndVel();
	}
}

//デバッグ
void MoveSM::StepObjMoveDebug()
{
	#pragma omp parallel for
	for(int i = 0; i < IceObject::GetParticleNum(); ++i)
	{	
		if(m_iceJudge->JudgeMoveDebug(i) == false){	continue;	}

		//m_iceSM[i]->UpdateCPU();
	}	
}

void MoveSM::StepObjMoveItrDebug()
{
	#pragma omp parallel for
	for(int i = 0; i < IceObject::GetParticleNum(); ++i)
	{	
		if(m_iceJudge->JudgeMoveDebug(i) == false){	continue;	}

		//m_iceSM[i]->ShapeMatchingIteration();
	}	
}