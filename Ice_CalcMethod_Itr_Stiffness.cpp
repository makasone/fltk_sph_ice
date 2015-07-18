#include "Ice_CalcMethod_Itr_Stiffness.h"

typedef Ice_CalcMethod_Itr_Stiffness CalcIteration;

//CalcIteration::Ice_CalcMethod_Itr_Stiffness(const vector<Ice_SM*>& iceSM, Ice_ClusterMove* clusterMove, Ice_Convolution* convo)
CalcIteration::Ice_CalcMethod_Itr_Stiffness(const vector<Ice_SM*>& iceSM, Ice_SimuMethod* simuMethod, Ice_Convolution* convo)
{
	m_iceSM = iceSM;
	//SetObjMove(clusterMove);
	SetObjMove(simuMethod);
	SetConvolution(convo);

	//閾値用データ算出クラス
	//m_iceCalcStiff = new Ice_CalcStiffData_Summation(iceSM);	//総変形量
	m_iceCalcStiff = new Ice_CalcStiffData_Average(iceSM);		//平均変化量
	//m_iceCalcStiff = new Ice_CalcStiffData_StdDevision(iceSM, m_iceMove->GetJudgeMove());		//粒子位置の分散
	//m_iceCalcStiff = new Ice_CalcStiffData_CompareRigid(iceSM, m_iceMove->GetJudgeMove());	//剛体との差分
}

CalcIteration::~Ice_CalcMethod_Itr_Stiffness()
{
}

//void CalcIteration::SetObjMove(Ice_ClusterMove* clusterMove)
void CalcIteration::SetObjMove(Ice_SimuMethod* simuMethod)
{
	//m_iceMove = clusterMove;
	m_iceSimu = simuMethod;
}

void CalcIteration::SetConvolution(Ice_Convolution* convo)
{
	m_iceConvo = convo;
}

void CalcIteration::StepObjMove()
{
	//初回
	//m_iceMove->StepObjMove();
	m_iceSimu->StepObjMove();

	//rigidの最終結果　重力を反映しているバージョン
	m_iceCalcStiff->StepUpdate();

	//MSMの最終結果算出
	m_iceConvo->StepConvolution();

	//平均変形量測定
	float threshold = m_iceCalcStiff->StepCalcData();

	//反復	
	while(threshold > Ice_SM::GetItrStiffness())
	{
		//各クラスタで運動計算
		//m_iceMove->StepObjMoveItr();
		m_iceSimu->StepObjMoveItr();

		//rigidの最終結果
		m_iceCalcStiff->StepUpdateItr();

		//MSMの最終結果算出
		m_iceConvo->StepConvolution();

		//閾値測定
		threshold = m_iceCalcStiff->StepCalcData();
	}

	//速度算出
	#pragma omp parallel for
	for(int i = 0; i < IceObject::GetParticleNum(); ++i)
	{	
		//if(m_iceMove->GetJudgeMove()->JudgeMove(i) == false){	continue;	}
		if(m_iceSimu->GetJudgeMove()->JudgeMove(i) == false){	continue;	}
		m_iceSM[i]->integrateIteration();
	}
}

void CalcIteration::StepObjMoveDebug()
{
	//初回
	//m_iceMove->StepObjMoveDebug();
	m_iceSimu->StepObjMoveDebug();

	//rigidの最終結果　重力を反映しているバージョン
	m_iceCalcStiff->StepUpdate();

	//MSMの最終結果算出
	m_iceConvo->StepConvolution();

	//平均変形量測定
	float threshold = m_iceCalcStiff->StepCalcData();
	int count = 0;

	//cout << __FUNCTION__ << " threshold = " << threshold << endl;

	//反復	
	while(threshold > Ice_SM::GetItrStiffness())
	{
		//各クラスタで運動計算
		m_iceSimu->StepObjMoveItrDebug();

		//rigidの最終結果
		m_iceCalcStiff->StepUpdateItr();

		//MSMの最終結果算出
		m_iceConvo->StepConvolution();

		//閾値測定
		threshold = m_iceCalcStiff->StepCalcData();

		count++;

		cout << __FUNCTION__ << " thrshold:" << threshold << " count:" << count << endl;
	}

	//速度算出
	#pragma omp parallel for
	for(int i = 0; i < IceObject::GetParticleNum(); ++i)
	{	
		if(m_iceSimu->GetJudgeMove()->JudgeMoveDebug(i) == false){	continue;	}
		m_iceSM[i]->integrateIteration();
	}

//デバッグ
	cout << __FUNCTION__ << " " << count << endl;

	string result = RESULT_DATA_PATH;
	result += "Itr_Stiffness_CountNum.txt";
	ofstream ofs(result, ios::app);
	ofs << count << endl;
}

//硬さの閾値を測定する処理のデバッグ
void CalcIteration::DebugStiffness()
{
	m_iceCalcStiff->StepCalcDataDebug();
}