#include "Ice_CalcMethod_Itr_Stiffness.h"

typedef Ice_CalcMethod_Itr_Stiffness CalcIteration;

CalcIteration::Ice_CalcMethod_Itr_Stiffness(const vector<Ice_SM*>& iceSM, Ice_ClusterMove* clusterMove, Ice_Convolution* convo)
{
	m_iceSM = iceSM;
	SetObjMove(clusterMove);
	SetConvolution(convo);

	//臒l�p�f�[�^�Z�o�N���X
	//m_iceCalcStiff = new Ice_CalcStiffData_Summation(iceSM);	//���ό`��
	m_iceCalcStiff = new Ice_CalcStiffData_Average(iceSM);		//���ϕω���
	//m_iceCalcStiff = new Ice_CalcStiffData_StdDevision(iceSM, m_iceMove->GetJudgeMove());		//���q�ʒu�̕��U
	//m_iceCalcStiff = new Ice_CalcStiffData_CompareRigid(iceSM, m_iceMove->GetJudgeMove());	//���̂Ƃ̍���
}

CalcIteration::~Ice_CalcMethod_Itr_Stiffness()
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

	//rigid�̍ŏI���ʁ@�d�͂𔽉f���Ă���o�[�W����
	m_iceCalcStiff->StepUpdate();

	//MSM�̍ŏI���ʎZ�o
	m_iceConvo->StepConvolution();

	//���ϕό`�ʑ���
	float threshold = m_iceCalcStiff->StepCalcData();

	//����	
	while(threshold > Ice_SM::GetItrStiffness())
	{
		//�e�N���X�^�ŉ^���v�Z
		m_iceMove->StepObjMoveItr();

		//rigid�̍ŏI����
		m_iceCalcStiff->StepUpdateItr();

		//MSM�̍ŏI���ʎZ�o
		m_iceConvo->StepConvolution();

		//臒l����
		threshold = m_iceCalcStiff->StepCalcData();
	}

	//���x�Z�o
	#pragma omp parallel for
	for(int i = 0; i < IceObject::GetParticleNum(); ++i)
	{	
		if(m_iceMove->GetJudgeMove()->JudgeMove(i) == false){	continue;	}
		m_iceSM[i]->integrateIteration();
	}
}

void CalcIteration::StepObjMoveDebug()
{
	//����
	m_iceMove->StepObjMoveDebug();

	//rigid�̍ŏI���ʁ@�d�͂𔽉f���Ă���o�[�W����
	m_iceCalcStiff->StepUpdate();

	//MSM�̍ŏI���ʎZ�o
	m_iceConvo->StepConvolution();

	//���ϕό`�ʑ���
	float threshold = m_iceCalcStiff->StepCalcData();
	int count = 0;

	//cout << __FUNCTION__ << " threshold = " << threshold << endl;

	//����	
	while(threshold > Ice_SM::GetItrStiffness())
	{
		//�e�N���X�^�ŉ^���v�Z
		m_iceMove->StepObjMoveItrDebug();

		//rigid�̍ŏI����
		m_iceCalcStiff->StepUpdateItr();

		//MSM�̍ŏI���ʎZ�o
		m_iceConvo->StepConvolution();

		//臒l����
		threshold = m_iceCalcStiff->StepCalcData();

		count++;

		cout << __FUNCTION__ << " thrshold:" << threshold << " count:" << count << endl;
	}

	//���x�Z�o
	#pragma omp parallel for
	for(int i = 0; i < IceObject::GetParticleNum(); ++i)
	{	
		if(m_iceMove->GetJudgeMove()->JudgeMoveDebug(i) == false){	continue;	}
		m_iceSM[i]->integrateIteration();
	}

//�f�o�b�O
	cout << __FUNCTION__ << " " << count << endl;

	string result = RESULT_DATA_PATH;
	result += "Itr_Stiffness_CountNum.txt";
	ofstream ofs(result, ios::app);
	ofs << count << endl;
}

//�d����臒l�𑪒肷�鏈���̃f�o�b�O
void CalcIteration::DebugStiffness()
{
	m_iceCalcStiff->StepCalcDataDebug();
}