#include "Ice_CalcMethod_Iteration.h"

typedef Ice_CalcMethod_Iteration CalcIteration;

CalcIteration::Ice_CalcMethod_Iteration(const vector<Ice_SM*>& iceSM, Ice_ClusterMove* clusterMove, Ice_InterPolation* intrp)
{
	m_iceSM = iceSM;
	SetObjMove(clusterMove);
	SetIntrp(intrp);
}

CalcIteration::~Ice_CalcMethod_Iteration()
{
}

void CalcIteration::SetObjMove(Ice_ClusterMove* clusterMove)
{
	m_iceMove = clusterMove;
}

void CalcIteration::SetIntrp(Ice_InterPolation* intrp)
{
	m_iceInterPolation = intrp;
}

void CalcIteration::StepObjMove()
{
	//初回
	m_iceMove->StepObjMove();		//そのまま呼ぶだけ
	m_iceInterPolation->StepInterPolationItr();

	//反復
	for(int itr = 0; itr < Ice_SM::GetIteration(); itr++)
	{
		m_iceMove->StepObjMoveItr();
		m_iceInterPolation->StepInterPolationItr();
	}

	//速度算出
	#pragma omp parallel for
	for(int i = 0; i < IceObject::GetParticleNum(); ++i)
	{	
		if(m_iceMove->GetJudgeMove()->JudgeMove(i) == false){	continue;	}
		m_iceSM[i]->integrateIteration();
	}
}

//Weightな補間をしていないのに注意
void CalcIteration::InterPolationForCluster()
{
	unsigned sm_particleNum = IceObject::GetParticleNum();
	vector<unsigned> addParticleNum(sm_particleNum, 0);

	float* sldPos = Ice_SM::GetSldPosPointer();
	float* sldVel = Ice_SM::GetSldVelPointer();

	//sldの初期化
	for(int cIndx = 0; cIndx < sm_particleNum; cIndx++)
	{
		for(int dim = 0; dim < SM_DIM; dim++)
		{
			sldPos[cIndx*SM_DIM+dim] = 0.0f;
			sldVel[cIndx*SM_DIM+dim] = 0.0f;
		}
	}

	for(int cIndx = 0; cIndx < sm_particleNum; cIndx++)
	{
		//最終結果算出に用いるクラスタの判定
		if(m_iceInterPolation->GetIntrpJudge()->JudgeInterPolation(cIndx) == false){	continue;	}

		//クラスタが含んでいる各粒子の位置・速度を取得
		//各位置・速度をstaticな最終位置に足す
		for(int oIndx = 0; oIndx < m_iceSM[cIndx]->GetIndxNum(); oIndx++)
		{
			int pIndx = m_iceSM[cIndx]->GetParticleIndx(oIndx);

			if(pIndx == MAXINT){	continue;	}

			Vec3 pos = m_iceSM[cIndx]->GetVertexPos(oIndx);
			Vec3 vel = m_iceSM[cIndx]->GetVertexVel(oIndx);

			for(int dim = 0; dim < SM_DIM; dim++)
			{
				sldPos[pIndx*SM_DIM+dim] += pos[dim];
				sldVel[pIndx*SM_DIM+dim] += vel[dim];
			}

			//足した数をカウント
			addParticleNum[pIndx] += 1;
		}
	}

	//TODO::固体と液体で線形補完していないのに注意
	//足した数で平均し，粒子位置に反映
	for(int i = 0; i < sm_particleNum; i++)
	{
		int pIndx = i*SM_DIM;
		int smIndx = i*SM_DIM;

		float clusterNum = (float)addParticleNum[i];
		if(clusterNum == 0){	continue;	}

		//粒子位置に反映
		for(int dim = 0; dim < SM_DIM; dim++)
		{
			sldPos[pIndx+dim] = sldPos[smIndx+dim]/clusterNum;
			sldVel[pIndx+dim] = sldVel[smIndx+dim]/clusterNum;
		}
	}
}