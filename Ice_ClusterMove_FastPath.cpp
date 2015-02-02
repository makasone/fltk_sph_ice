#include "Ice_ClusterMove_Normal.h"

typedef Ice_ClusterMove_FastPath MoveFastPath;

MoveFastPath::Ice_ClusterMove_FastPath(const vector<Ice_SM*>& iceSM, Surf_SM* surfSM)
{
	m_iceSM = iceSM;
	m_surfSM = surfSM;
}

MoveFastPath::~Ice_ClusterMove_FastPath()
{
}

void MoveFastPath::SetJudgeMove(Ice_JudgeMove* judge)
{
	m_iceJudge = judge;	
}

Ice_JudgeMove* MoveFastPath::GetJudgeMove()
{
	return m_iceJudge;
}

void MoveFastPath::StepObjMove()
{
	m_surfSM->UpdatePrefixSum();							//prefixSum�̍X�V

	//�N���X�^�̃p�����[�^�X�V
	#pragma omp parallel for
	for(int i = 0; i < IceObject::GetParticleNum(); ++i)
	{	
		if(m_iceJudge->JudgeMove(i) == false){	continue;	}

		m_iceSM[i]->SetNowCm(m_surfSM->CalcCmSum(i));		//�d�S�̍X�V
		m_iceSM[i]->SetApq(m_surfSM->CalcApqSum(i));		//�ό`�s��̍X�V
		m_iceSM[i]->UpdateUsePathCPU();						//�^���v�Z
	}
}

void MoveFastPath::StepObjMoveItr()
{
	m_surfSM->UpdatePrefixSumItr();							//prefixSum�̍X�V

	//�e�N���X�^�̃f�[�^�X�V
	#pragma omp parallel for
	for(int i = 0; i < IceObject::GetParticleNum(); ++i)
	{	
		if(m_iceJudge->JudgeMove(i) == false){	continue;	}

		m_iceSM[i]->SetNowCm(m_surfSM->CalcCmSum(i));		//�d�S�̍X�V
		m_iceSM[i]->SetApq(m_surfSM->CalcApqSum(i));		//�ό`�s��̍X�V
		m_iceSM[i]->ShapeMatchingUsePath();					//���݂̈ʒu��SM�@
	}
}

void MoveFastPath::StepObjMoveDebug()
{
	m_surfSM->UpdatePrefixSum();							//prefixSum�̍X�V

	//�N���X�^�̃p�����[�^�X�V
	#pragma omp parallel for
	for(int i = 0; i < IceObject::GetParticleNum(); ++i)
	{	
		if(m_iceJudge->JudgeMoveDebug(i) == false){	continue;	}

		m_iceSM[i]->SetNowCm(m_surfSM->CalcCmSum(i));		//�d�S�̍X�V
		m_iceSM[i]->SetApq(m_surfSM->CalcApqSum(i));		//�ό`�s��̍X�V
		m_iceSM[i]->UpdateUsePathCPU();						//�^���v�Z
	}
}

void MoveFastPath::StepObjMoveItrDebug()
{
	m_surfSM->UpdatePrefixSumItr();							//prefixSum�̍X�V

	//�e�N���X�^�̃f�[�^�X�V
	#pragma omp parallel for
	for(int i = 0; i < IceObject::GetParticleNum(); ++i)
	{	
		if(m_iceJudge->JudgeMoveDebug(i) == false){	continue;	}

		m_iceSM[i]->SetNowCm(m_surfSM->CalcCmSum(i));		//�d�S�̍X�V
		m_iceSM[i]->SetApq(m_surfSM->CalcApqSum(i));		//�ό`�s��̍X�V
		m_iceSM[i]->ShapeMatchingUsePath();					//���݂̈ʒu��SM�@
	}
}