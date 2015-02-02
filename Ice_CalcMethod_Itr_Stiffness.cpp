#include "Ice_CalcMethod_Itr_Stiffness.h"

typedef Ice_CalcMethod_Itr_Stiffness CalcIteration;

CalcIteration::Ice_CalcMethod_Itr_Stiffness(const vector<Ice_SM*>& iceSM, Ice_ClusterMove* clusterMove, Ice_Convolution* convo)
{
	m_iceSM = iceSM;
	SetObjMove(clusterMove);
	SetConvolution(convo);

	//閾値用データ算出クラス
	//m_iceCalcStiff = new Ice_CalcStiffData_Summation(iceSM);	//総変形量
	m_iceCalcStiff = new Ice_CalcStiffData_StdDevision(iceSM, m_iceMove->GetJudgeMove());	//粒子位置の分散
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
	//初回の運動計算
	m_iceMove->StepObjMove();		//そのまま呼ぶだけ

	//反復するか判定するための値
	float threshold = m_iceCalcStiff->StepCalcData();

	//反復	
	while(threshold > Ice_SM::GetItrStiffness())
	{
		//運動計算
		m_iceConvo->StepConvolution();
		m_iceMove->StepObjMoveItr();

		//反復するか判定するための値
		threshold = m_iceCalcStiff->StepCalcData();
	}

	//速度算出
	#pragma omp parallel for
	for(int i = 0; i < IceObject::GetParticleNum(); ++i)
	{	
		if(m_iceMove->GetJudgeMove()->JudgeMove(i) == false){	continue;	}
		m_iceSM[i]->integrateIteration();
	}
}

void CalcIteration::StepObjMoveDebug()
{
	//初回
	m_iceMove->StepObjMoveDebug();

	//平均変形量測定
	float threshold = m_iceCalcStiff->StepCalcData();
	int count = 0;

	//cout << __FUNCTION__ << " threshold = " << threshold << endl;

	//反復	
	while(threshold > Ice_SM::GetItrStiffness())
	{
		//運動計算
		m_iceConvo->StepConvolution();
		m_iceMove->StepObjMoveItrDebug();

		//平均変形量測定
		threshold = m_iceCalcStiff->StepCalcData();

		//cout << __FUNCTION__ << " threshold = " << threshold << endl;

		count++;
	}

	//速度算出
	#pragma omp parallel for
	for(int i = 0; i < IceObject::GetParticleNum(); ++i)
	{	
		if(m_iceMove->GetJudgeMove()->JudgeMoveDebug(i) == false){	continue;	}
		m_iceSM[i]->integrateIteration();
	}

//デバッグ
	cout << __FUNCTION__ << " " << count << endl;

	string result = RESULT_DATA_PATH;
	result += "Itr_Stiffness_CountNum.txt";
	ofstream ofs(result, ios::app);
	ofs << count << endl;
}