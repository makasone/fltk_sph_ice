#include "Ice_SimuMethod_DistanceConstraint.h"

typedef Ice_SimuMethod_DC MoveDC;


MoveDC::Ice_SimuMethod_DC(const vector<Ice_SM*>& iceSM)
{
	m_iceSM = iceSM;
}

MoveDC::Ice_SimuMethod_DC(const vector<ElasticObj*>& elasticObj, const vector<OrientedParticle*>& particles)
{
	m_elasticObj = elasticObj;
	m_vOrientedPrtes = particles;
}

MoveDC::~Ice_SimuMethod_DC()
{
}

void MoveDC::SetJudgeMove(Ice_JudgeMove* judge)
{
	m_iceJudge = judge;	
}

Ice_JudgeMove* MoveDC::GetJudgeMove()
{
	return m_iceJudge;
}

void MoveDC::StepObjMove()
{
	//マウスによるドラッグを反映させるために，無理やり値を更新
	OrientedParticleBaseElasticObject::CopyPrtToClstrPos(IceObject::GetParticleNum());

	//運動計算
	#pragma omp parallel for
	for(int i = 0; i < IceObject::GetParticleNum(); i++){	
		if(! m_iceJudge->JudgeMove(i)){	continue;	}

		m_elasticObj[i]->UpdateCluster_DC();
	}
}

void MoveDC::StepObjMoveItr()
{
	#pragma omp parallel for
	for(int i = 0; i < IceObject::GetParticleNum(); ++i){	
		if(! m_iceJudge->JudgeMove(i)){	continue;	}

		m_elasticObj[i]->UpdateCluster_DC_Itr();
	}
}

//粒子情報の更新
void MoveDC::StepObjUpdate()
{
	//速度，角速度，姿勢を更新
	#pragma omp parallel for
	for(int i = 0; i < IceObject::GetParticleNum(); i++){
		m_vOrientedPrtes[i]->UpdatePosAndVel();
	}
}

void MoveDC::StepObjUpdateItr()
{

}

void MoveDC::StepObjUpdateItrEnd()
{
	//速度，角速度，姿勢を更新
	#pragma omp parallel for
	for(int i = 0; i < IceObject::GetParticleNum(); i++){
		m_vOrientedPrtes[i]->UpdatePosAndVel();
	}
}

//デバッグ
void MoveDC::StepObjMoveDebug()
{
	#pragma omp parallel for
	for(int i = 0; i < IceObject::GetParticleNum(); ++i)
	{	
		if(m_iceJudge->JudgeMoveDebug(i) == false){	continue;	}

		//m_iceSM[i]->UpdateCPU();
	}	
}

void MoveDC::StepObjMoveItrDebug()
{
	#pragma omp parallel for
	for(int i = 0; i < IceObject::GetParticleNum(); ++i)
	{	
		if(m_iceJudge->JudgeMoveDebug(i) == false){	continue;	}

		//m_iceSM[i]->ShapeMatchingIteration();
	}	
}