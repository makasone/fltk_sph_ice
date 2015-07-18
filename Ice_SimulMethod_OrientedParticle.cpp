#include "Ice_SimuMethod_ShapeMatching.h"

typedef Ice_SimuMethod_OP MoveSM;


MoveSM::Ice_SimuMethod_OP(const vector<Ice_SM*>& iceSM)
{
	m_iceSM = iceSM;
}

MoveSM::Ice_SimuMethod_OP(const vector<ElasticObj*>& elasticObj, const vector<OrientedParticle*>& particles, IceStructure* iceStrct)
{
	m_elasticObj = elasticObj;
	m_vOrientedPrtes = particles;
	m_iceStrct = iceStrct;
}

MoveSM::~Ice_SimuMethod_OP()
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
	////�^���v�Z
	//#pragma omp parallel for
	//for(int i = 0; i < IceObject::GetParticleNum(); ++i)
	//{	
	//	if(m_iceJudge->JudgeMove(i) == false){	continue;	}

	//	m_iceSM[i]->UpdateCPU();
	//}

	//�}�E�X�ɂ��h���b�O�𔽉f�����邽�߂ɁC�������l���X�V
	OrientedParticleBaseElasticObject::CopyPrtToClstrPos(IceObject::GetParticleNum());

	//�T���v�����O���q���X�V
	#pragma omp parallel for
	for(int i = 0; i < IceObject::GetParticleNum(); i++){
		if(! m_iceJudge->JudgeMove(i)){	continue;	}

		m_vOrientedPrtes[i]->Integrate();
	}

	//��T���v�����O���q���X�V
	#pragma omp parallel for
	for(int i = 0; i < IceObject::GetParticleNum(); i++){
		//��T���v�����O���q�̏����X�V
		if(m_iceJudge->JudgeMove(i)){	continue;	}

		//m_vOrientedPrtes[i]->Integrate_NotSampled();
		m_vOrientedPrtes[i]->Integrate_NotSampled2(m_iceStrct);
	}

	//�^���v�Z
	#pragma omp parallel for
	for(int i = 0; i < IceObject::GetParticleNum(); i++){	
		if(! m_iceJudge->JudgeMove(i)){	continue;	}

		m_elasticObj[i]->UpdateCluster_OP();
	}
}

void MoveSM::StepObjMoveItr()
{
	//#pragma omp parallel for
	//for(int i = 0; i < IceObject::GetParticleNum(); ++i){	
	//	if(! m_iceJudge->JudgeMove(i)){	continue;	}

	//	m_iceSM[i]->ShapeMatchingIteration();
	//}

	//�T���v�����O���q���X�V
	#pragma omp parallel for
	for(int i = 0; i < IceObject::GetParticleNum(); i++){
		if(! m_iceJudge->JudgeMove(i)){	continue;	}

		m_vOrientedPrtes[i]->Integrate_Itr();
	}

	//��T���v�����O���q���X�V
	#pragma omp parallel for
	for(int i = 0; i < IceObject::GetParticleNum(); i++){
		//��T���v�����O���q�̏����X�V
		if(m_iceJudge->JudgeMove(i)){	continue;	}

		//m_vOrientedPrtes[i]->Integrate_NotSampled_Itr();	//������g���ƁCUpdateItr�̒��g��ʏ�Ɣ������œ�ʂ�p�ӂ��Ȃ��Ƃ����Ȃ�
		m_vOrientedPrtes[i]->Integrate_NotSampled2_Itr(m_iceStrct);
	}

	//�^���v�Z
	#pragma omp parallel for
	for(int i = 0; i < IceObject::GetParticleNum(); i++){	
		if(! m_iceJudge->JudgeMove(i)){	continue;	}

		m_elasticObj[i]->UpdateCluster_OP();
	}
}

//���q���̍X�V
void MoveSM::StepObjUpdate()
{
	//�ʔ������ƂɁC�����p�̏����Q��g�ݍ��킹�����̂��ʏ�̉^���v�Z��\���ł����D
	//�܂�C���O��ς���ׂ�
	StepObjUpdateItr();
	StepObjUpdateItrEnd();
}

void MoveSM::StepObjUpdateItr()
{
	//���x�C�p���x�C�p�����X�V
	#pragma omp parallel for
	for(int i = 0; i < IceObject::GetParticleNum(); i++){
		if(! m_iceJudge->JudgeMove(i)){	continue;	}

		m_vOrientedPrtes[i]->UpdateOrientation();
	}
}

void MoveSM::StepObjUpdateItrEnd()
{
	//���q�ʒu�C���x�C�p���x�C�p�����X�V
	#pragma omp parallel for
	for(int i = 0; i < IceObject::GetParticleNum(); i++){
		m_vOrientedPrtes[i]->Update_ItrEnd();
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