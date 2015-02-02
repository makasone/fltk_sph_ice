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
	//初回
	m_iceMove->StepObjMove();

	vector<Ice_SM*> expandObj;

	//反復
	for(int itr = 1; itr < Ice_SM::GetIteration(); itr++)
	{
		m_iceConvo->StepConvolution();
		m_iceMove->StepObjMoveItr();

		//変形量の大きいクラスタを検出
		GetExpandeCluster(expandObj);

		//クラスタの影響範囲を拡大して再構成
		ExpandeCluster(expandObj);
	}

	//クラスタの構造を元に戻す
	ContractCluster(expandObj);

	//速度算出
	#pragma omp parallel for
	for(int i = 0; i < IceObject::GetParticleNum(); ++i)
	{	
		if(m_iceMove->GetJudgeMove()->JudgeMove(i) == false){	continue;	}
		m_iceSM[i]->integrateIteration();
	}
}

//再構築するクラスタを検出
void CalcIteration::GetExpandeCluster(vector<Ice_SM*>& remakeObjes)
{
	for(int i = 0; i < IceObject::GetParticleNum(); ++i)
	{

	}
}

//クラスタの影響範囲を拡大して再構成
void CalcIteration::ExpandeCluster(vector<Ice_SM*>& remakeObjes)
{
	for(int i = 0; i < IceObject::GetParticleNum(); ++i)
	{

	}
}

//クラスタの影響範囲を縮小して再構成
void CalcIteration::ContractCluster(vector<Ice_SM*>& remakeobjes)
{
	for(int i = 0; i < IceObject::GetParticleNum(); ++i)
	{

	}
}

//デバッグ
void CalcIteration::StepObjMoveDebug()
{
	//初回
	m_iceMove->StepObjMoveDebug();

	//反復
	for(int itr = 1; itr < Ice_SM::GetIteration(); itr++)
	{
		m_iceConvo->StepConvolution();
		m_iceMove->StepObjMoveItrDebug();
	}

	//速度算出
	#pragma omp parallel for
	for(int i = 0; i < IceObject::GetParticleNum(); ++i)
	{	
		if(m_iceMove->GetJudgeMove()->JudgeMoveDebug(i) == false){	continue;	}
		m_iceSM[i]->integrateIteration();
	}
}