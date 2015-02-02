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
	m_surfSM->UpdatePrefixSum();							//prefixSumの更新

	//クラスタのパラメータ更新
	#pragma omp parallel for
	for(int i = 0; i < IceObject::GetParticleNum(); ++i)
	{	
		if(m_iceJudge->JudgeMove(i) == false){	continue;	}

		m_iceSM[i]->SetNowCm(m_surfSM->CalcCmSum(i));		//重心の更新
		m_iceSM[i]->SetApq(m_surfSM->CalcApqSum(i));		//変形行列の更新
		m_iceSM[i]->UpdateUsePathCPU();						//運動計算
	}
}

void MoveFastPath::StepObjMoveItr()
{
	m_surfSM->UpdatePrefixSumItr();							//prefixSumの更新

	//各クラスタのデータ更新
	#pragma omp parallel for
	for(int i = 0; i < IceObject::GetParticleNum(); ++i)
	{	
		if(m_iceJudge->JudgeMove(i) == false){	continue;	}

		m_iceSM[i]->SetNowCm(m_surfSM->CalcCmSum(i));		//重心の更新
		m_iceSM[i]->SetApq(m_surfSM->CalcApqSum(i));		//変形行列の更新
		m_iceSM[i]->ShapeMatchingUsePath();					//現在の位置でSM法
	}
}

void MoveFastPath::StepObjMoveDebug()
{
	m_surfSM->UpdatePrefixSum();							//prefixSumの更新

	//クラスタのパラメータ更新
	#pragma omp parallel for
	for(int i = 0; i < IceObject::GetParticleNum(); ++i)
	{	
		if(m_iceJudge->JudgeMoveDebug(i) == false){	continue;	}

		m_iceSM[i]->SetNowCm(m_surfSM->CalcCmSum(i));		//重心の更新
		m_iceSM[i]->SetApq(m_surfSM->CalcApqSum(i));		//変形行列の更新
		m_iceSM[i]->UpdateUsePathCPU();						//運動計算
	}
}

void MoveFastPath::StepObjMoveItrDebug()
{
	m_surfSM->UpdatePrefixSumItr();							//prefixSumの更新

	//各クラスタのデータ更新
	#pragma omp parallel for
	for(int i = 0; i < IceObject::GetParticleNum(); ++i)
	{	
		if(m_iceJudge->JudgeMoveDebug(i) == false){	continue;	}

		m_iceSM[i]->SetNowCm(m_surfSM->CalcCmSum(i));		//重心の更新
		m_iceSM[i]->SetApq(m_surfSM->CalcApqSum(i));		//変形行列の更新
		m_iceSM[i]->ShapeMatchingUsePath();					//現在の位置でSM法
	}
}