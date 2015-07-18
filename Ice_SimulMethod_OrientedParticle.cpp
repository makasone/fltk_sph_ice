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
	////運動計算
	//#pragma omp parallel for
	//for(int i = 0; i < IceObject::GetParticleNum(); ++i)
	//{	
	//	if(m_iceJudge->JudgeMove(i) == false){	continue;	}

	//	m_iceSM[i]->UpdateCPU();
	//}

	//マウスによるドラッグを反映させるために，無理やり値を更新
	OrientedParticleBaseElasticObject::CopyPrtToClstrPos(IceObject::GetParticleNum());

	//サンプリング粒子を更新
	#pragma omp parallel for
	for(int i = 0; i < IceObject::GetParticleNum(); i++){
		if(! m_iceJudge->JudgeMove(i)){	continue;	}

		m_vOrientedPrtes[i]->Integrate();
	}

	//非サンプリング粒子を更新
	#pragma omp parallel for
	for(int i = 0; i < IceObject::GetParticleNum(); i++){
		//非サンプリング粒子の情報を更新
		if(m_iceJudge->JudgeMove(i)){	continue;	}

		//m_vOrientedPrtes[i]->Integrate_NotSampled();
		m_vOrientedPrtes[i]->Integrate_NotSampled2(m_iceStrct);
	}

	//運動計算
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

	//サンプリング粒子を更新
	#pragma omp parallel for
	for(int i = 0; i < IceObject::GetParticleNum(); i++){
		if(! m_iceJudge->JudgeMove(i)){	continue;	}

		m_vOrientedPrtes[i]->Integrate_Itr();
	}

	//非サンプリング粒子を更新
	#pragma omp parallel for
	for(int i = 0; i < IceObject::GetParticleNum(); i++){
		//非サンプリング粒子の情報を更新
		if(m_iceJudge->JudgeMove(i)){	continue;	}

		//m_vOrientedPrtes[i]->Integrate_NotSampled_Itr();	//これを使うと，UpdateItrの中身を通常と反復時で二通り用意しないといけない
		m_vOrientedPrtes[i]->Integrate_NotSampled2_Itr(m_iceStrct);
	}

	//運動計算
	#pragma omp parallel for
	for(int i = 0; i < IceObject::GetParticleNum(); i++){	
		if(! m_iceJudge->JudgeMove(i)){	continue;	}

		m_elasticObj[i]->UpdateCluster_OP();
	}
}

//粒子情報の更新
void MoveSM::StepObjUpdate()
{
	//面白いことに，反復用の処理２つを組み合わせたものが通常の運動計算を表現できた．
	//つまり，名前を変えるべき
	StepObjUpdateItr();
	StepObjUpdateItrEnd();
}

void MoveSM::StepObjUpdateItr()
{
	//速度，角速度，姿勢を更新
	#pragma omp parallel for
	for(int i = 0; i < IceObject::GetParticleNum(); i++){
		if(! m_iceJudge->JudgeMove(i)){	continue;	}

		m_vOrientedPrtes[i]->UpdateOrientation();
	}
}

void MoveSM::StepObjUpdateItrEnd()
{
	//粒子位置，速度，角速度，姿勢を更新
	#pragma omp parallel for
	for(int i = 0; i < IceObject::GetParticleNum(); i++){
		m_vOrientedPrtes[i]->Update_ItrEnd();
	}
}

//デバッグ
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