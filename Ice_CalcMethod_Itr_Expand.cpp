#include "Ice_CalcMethod_Itr_Expand.h"

typedef Ice_CalcMethod_Itr_Expand CalcIteration;

CalcIteration::Ice_CalcMethod_Itr_Expand(const vector<Ice_SM*>& iceSM, Ice_ClusterMove* clusterMove, Ice_Convolution* convo)
{
	m_iceSM = iceSM;
	SetObjMove(clusterMove);
	SetConvolution(convo);
}

CalcIteration::~Ice_CalcMethod_Itr_Expand()
{
}

void CalcIteration::SetObjMove(Ice_ClusterMove* clusterMove)
{
	m_iceMove = clusterMove;
}

void CalcIteration::SetConvolution(Ice_Convolution* convo)
{
	m_iceConvo = convo;
}

void CalcIteration::StepObjMove()
{
	//����
	m_iceMove->StepObjMove();

	vector<Ice_SM*> expandObj;

	//����
	for(int itr = 1; itr < Ice_SM::GetIteration(); itr++)
	{
		m_iceConvo->StepConvolution();
		m_iceMove->StepObjMoveItr();

		//�ό`�ʂ̑傫���N���X�^�����o
		GetExpandeCluster(expandObj);

		//�N���X�^�̉e���͈͂��g�債�čč\��
		ExpandeCluster(expandObj);
	}

	//�N���X�^�̍\�������ɖ߂�
	ContractCluster(expandObj);

	//���x�Z�o
	#pragma omp parallel for
	for(int i = 0; i < IceObject::GetParticleNum(); ++i)
	{	
		if(m_iceMove->GetJudgeMove()->JudgeMove(i) == false){	continue;	}
		m_iceSM[i]->integrateIteration();
	}
}

//�č\�z����N���X�^�����o
void CalcIteration::GetExpandeCluster(vector<Ice_SM*>& remakeObjes)
{
	for(int i = 0; i < IceObject::GetParticleNum(); ++i)
	{

	}
}

//�N���X�^�̉e���͈͂��g�債�čč\��
void CalcIteration::ExpandeCluster(vector<Ice_SM*>& remakeObjes)
{
	for(int i = 0; i < IceObject::GetParticleNum(); ++i)
	{

	}
}

//�N���X�^�̉e���͈͂��k�����čč\��
void CalcIteration::ContractCluster(vector<Ice_SM*>& remakeobjes)
{
	for(int i = 0; i < IceObject::GetParticleNum(); ++i)
	{

	}
}

//�f�o�b�O
void CalcIteration::StepObjMoveDebug()
{
	//����
	m_iceMove->StepObjMoveDebug();

	//����
	for(int itr = 1; itr < Ice_SM::GetIteration(); itr++)
	{
		m_iceConvo->StepConvolution();
		m_iceMove->StepObjMoveItrDebug();
	}

	//���x�Z�o
	#pragma omp parallel for
	for(int i = 0; i < IceObject::GetParticleNum(); ++i)
	{	
		if(m_iceMove->GetJudgeMove()->JudgeMoveDebug(i) == false){	continue;	}
		m_iceSM[i]->integrateIteration();
	}
}