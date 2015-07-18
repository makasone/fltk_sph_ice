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
	//�}�E�X�ɂ��h���b�O�𔽉f�����邽�߂ɁC�������l���X�V
	OrientedParticleBaseElasticObject::CopyPrtToClstrPos(IceObject::GetParticleNum());

	//�^���v�Z
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

//���q���̍X�V
void MoveSM::StepObjUpdate()
{
	//���x�C�p���x�C�p�����X�V
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
	////���q�ʒu�C���x�C�p���x�C�p�����X�V
	//#pragma omp parallel for
	//for(int i = 0; i < IceObject::GetParticleNum(); i++){
	//	//if(! m_iceJudge->JudgeMove(i)){	continue;	}

	//	m_elasticObj[i]->UpdateCluster_SM_Itr_End();
	//}

	//���x�C�p���x�C�p�����X�V
	#pragma omp parallel for
	for(int i = 0; i < IceObject::GetParticleNum(); i++){
		m_vOrientedPrtes[i]->UpdatePosAndVel();
	}
}

//�f�o�b�O
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