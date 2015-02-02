#include "Ice_CalcMethod_Itr_Stiffness.h"

typedef Ice_CalcMethod_Itr_Stiffness CalcIteration;

CalcIteration::Ice_CalcMethod_Itr_Stiffness(const vector<Ice_SM*>& iceSM, Ice_ClusterMove* clusterMove, Ice_Convolution* convo)
{
	m_iceSM = iceSM;
	SetObjMove(clusterMove);
	SetConvolution(convo);

	//臒l�p�f�[�^�Z�o�N���X
	//m_iceCalcStiff = new Ice_CalcStiffData_Summation(iceSM);	//���ό`��
	m_iceCalcStiff = new Ice_CalcStiffData_StdDevision(iceSM, m_iceMove->GetJudgeMove());	//���q�ʒu�̕��U
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
	//����̉^���v�Z
	m_iceMove->StepObjMove();		//���̂܂܌ĂԂ���

	//�������邩���肷�邽�߂̒l
	float threshold = m_iceCalcStiff->StepCalcData();

	//����	
	while(threshold > Ice_SM::GetItrStiffness())
	{
		//�^���v�Z
		m_iceConvo->StepConvolution();
		m_iceMove->StepObjMoveItr();

		//�������邩���肷�邽�߂̒l
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

	//���ϕό`�ʑ���
	float threshold = m_iceCalcStiff->StepCalcData();
	int count = 0;

	//cout << __FUNCTION__ << " threshold = " << threshold << endl;

	//����	
	while(threshold > Ice_SM::GetItrStiffness())
	{
		//�^���v�Z
		m_iceConvo->StepConvolution();
		m_iceMove->StepObjMoveItrDebug();

		//���ϕό`�ʑ���
		threshold = m_iceCalcStiff->StepCalcData();

		//cout << __FUNCTION__ << " threshold = " << threshold << endl;

		count++;
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