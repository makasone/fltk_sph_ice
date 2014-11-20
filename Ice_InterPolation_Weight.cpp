#include "Ice_InterPolation_Weight.h"

typedef Ice_InterPolation_Weight IntrpWeight;

IntrpWeight::Ice_InterPolation_Weight(const vector<Ice_SM*>& iceSM, IceStructure* iceStrct)
{
	m_iceSM = iceSM;
	m_iceStrct = iceStrct;
}

IntrpWeight::~Ice_InterPolation_Weight()
{
}

void IntrpWeight::SetIntrpJudge(Ice_InterPolationJudge* judge)
{
	m_iceJudge = judge;	
}

Ice_InterPolationJudge* IntrpWeight::GetIntrpJudge()
{
	return m_iceJudge;
}

void IntrpWeight::StepInterPolation()
{
	unsigned sm_particleNum = IceObject::GetParticleNum();
	vector<float> deformationSum(sm_particleNum, 0.0f);

	float* sldPos = Ice_SM::GetSldPosPointer();
	float* sldVel = Ice_SM::GetSldVelPointer();

	float* s_sphPrtPos = IceObject::GetSPHHostPosPointer();
	float* s_sphPrtVel = IceObject::GetSPHHostVelPointer();

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
		if(m_iceJudge->JudgeInterPolation(cIndx) == false){	continue;	}

		//クラスタが含んでいる各粒子の位置・速度を取得
		//各位置・速度をstaticな最終位置に足す
		for(int oIndx = 0; oIndx < m_iceSM[cIndx]->GetIndxNum(); oIndx++)
		{
			int pIndx = m_iceSM[cIndx]->GetParticleIndx(oIndx);

			if(pIndx == MAXINT){	continue;	}

			float defAmount = m_iceSM[cIndx]->GetDefAmount();
			//defAmount = pow(defAmount, 2.0f);

			Vec3 pos = m_iceSM[cIndx]->GetVertexPos(oIndx) * defAmount;
			Vec3 vel = m_iceSM[cIndx]->GetVertexVel(oIndx) * defAmount;

			for(int dim = 0; dim < SM_DIM; dim++)
			{
				sldPos[pIndx*SM_DIM+dim] += pos[dim];
				sldVel[pIndx*SM_DIM+dim] += vel[dim];
			}

			//総変形量をカウント
			deformationSum[pIndx] += defAmount;
		}
	}

	//TODO::固体と液体で線形補完していないのに注意
	//足した数で平均し，粒子位置に反映
	for(int i = 0; i < sm_particleNum; i++)
	{
		int sphIndx = i*4;
		int smIndx = i*SM_DIM;

		float clusterNum = (float)deformationSum[i];
		if(clusterNum == 0){	continue;	}

		////デバッグ
		//if(clusterNum == 0)
		//{	
		//	for(int dim = 0; dim < SM_DIM; dim++)
		//	{
		//		s_sphPrtPos[sphIndx+dim] = 0.0f;
		//		s_sphPrtVel[sphIndx+dim] = 0.0f;
		//	}
		//	continue;	
		//}

		//粒子位置に反映
		for(int dim = 0; dim < SM_DIM; dim++)
		{
			s_sphPrtPos[sphIndx+dim] = sldPos[smIndx+dim]/clusterNum;
			s_sphPrtVel[sphIndx+dim] = sldVel[smIndx+dim]/clusterNum;
		}
	}

	//	//液体と固体の補間
	//	int sphIndx = pIndx*4;
	//	float intrpCff = m_iceObj->GetInterPolationCff(pIndx);
	//	double intrps = 1.0 - intrpCff;	//補間係数
	//
	//	for(int j = 0; j < 3; j++)
	//	{
	//		s_sphPrtVel[sphIndx+j] = vel[j] * intrpCff + s_sphPrtVel[sphIndx+j] * intrps;
	//		s_sphPrtPos[sphIndx+j] = pos[j] * intrpCff + s_sphPrtPos[sphIndx+j] * intrps;
	//	}
	//}
}

void IntrpWeight::StepInterPolationItr()
{
	unsigned sm_particleNum = IceObject::GetParticleNum();
	vector<float> deformationSum(sm_particleNum, 0.0f);

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
		if(m_iceJudge->JudgeInterPolation(cIndx) == false){	continue;	}

		//クラスタが含んでいる各粒子の位置・速度を取得
		//各位置・速度をstaticな最終位置に足す
		for(int oIndx = 0; oIndx < m_iceSM[cIndx]->GetIndxNum(); oIndx++)
		{
			int pIndx = m_iceSM[cIndx]->GetParticleIndx(oIndx);

			if(pIndx == MAXINT){	continue;	}

			float defAmount = m_iceSM[cIndx]->GetDefAmount();
			//defAmount = pow(defAmount, 2.0f);

			Vec3 pos = m_iceSM[cIndx]->GetVertexPos(oIndx) * defAmount;
			Vec3 vel = m_iceSM[cIndx]->GetVertexVel(oIndx) * defAmount;

			for(int dim = 0; dim < SM_DIM; dim++)
			{
				sldPos[pIndx*SM_DIM+dim] += pos[dim];
				sldVel[pIndx*SM_DIM+dim] += vel[dim];
			}

			//総変形量をカウント
			deformationSum[pIndx] += defAmount;
		}
	}

	//TODO::固体と液体で線形補完していないのに注意
	//足した数で平均し，粒子位置に反映
	for(int i = 0; i < sm_particleNum; i++)
	{
		int pIndx = i*SM_DIM;
		int smIndx = i*SM_DIM;

		float clusterNum = (float)deformationSum[i];
		if(clusterNum == 0){	continue;	}

		//粒子位置に反映
		for(int dim = 0; dim < SM_DIM; dim++)
		{
			sldPos[pIndx+dim] = sldPos[smIndx+dim]/clusterNum;
			sldVel[pIndx+dim] = sldVel[smIndx+dim]/clusterNum;
		}
	}
}