#include "Ice_Convolution_Weight.h"

typedef Ice_Convolution_Weight ConvoWeight;

ConvoWeight::Ice_Convolution_Weight(const vector<Ice_SM*>& iceSM, IceStructure* iceStrct)
{
	m_iceSM = iceSM;
	m_iceStrct = iceStrct;
	m_kernelDegree = 1.0f;
}

ConvoWeight::~Ice_Convolution_Weight()
{
}

void ConvoWeight::SetConvoJudge(Ice_ConvoJudge* judge)
{
	m_iceJudge = judge;	
}

Ice_ConvoJudge* ConvoWeight::GetConvoJudge()
{
	return m_iceJudge;
}

void ConvoWeight::StepConvolution()
{
	unsigned sm_particleNum = IceObject::GetParticleNum();
	vector<float> deformationSum(sm_particleNum, 0.0f);

	float* sldPos = Ice_SM::GetSldPosPointer();
	float* sldVel = Ice_SM::GetSldVelPointer();

	float* s_sphPrtPos = IceObject::GetSPHHostPosPointer();
	float* s_sphPrtVel = IceObject::GetSPHHostVelPointer();

	//sldの初期化
	Ice_SM::ResetFinalParamPointer(sm_particleNum);

	for(int cIndx = 0; cIndx < sm_particleNum; cIndx++)
	{
		//最終結果算出に用いるクラスタの判定
		if(m_iceJudge->JudgeConvolution(cIndx) == false){	continue;	}

		//クラスタが含んでいる各粒子の位置・速度を取得
		//各位置・速度をstaticな最終位置に足す
		for(int oIndx = 0; oIndx < m_iceSM[cIndx]->GetIndxNum(); oIndx++)
		{
			int pIndx = m_iceSM[cIndx]->GetParticleIndx(oIndx);
			if(pIndx == MAXINT){	continue;	}

			float defAmount = pow(m_iceSM[cIndx]->GetDefAmount(), m_kernelDegree);

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
		int smIndx = i*SM_DIM;

		float clusterNum = (float)deformationSum[i];
		if(clusterNum == 0){	continue;	}

		//固体の最終位置
		for(int dim = 0; dim < SM_DIM; dim++)
		{
			sldPos[smIndx+dim] /= clusterNum;
			sldVel[smIndx+dim] /= clusterNum;
		}
	}
}

void ConvoWeight::StepConvolutionDebug()
{
	unsigned sm_particleNum = IceObject::GetParticleNum();
	vector<float> deformationSum(sm_particleNum, 0.0f);

	float* sldPos = Ice_SM::GetSldPosPointer();
	float* sldVel = Ice_SM::GetSldVelPointer();

	float* s_sphPrtPos = IceObject::GetSPHHostPosPointer();
	float* s_sphPrtVel = IceObject::GetSPHHostVelPointer();

	//sldの初期化
	Ice_SM::ResetFinalParamPointer(sm_particleNum);

	for(int cIndx = 0; cIndx < sm_particleNum; cIndx++)
	{
		//最終結果算出に用いるクラスタの判定
		if(m_iceJudge->JudgeConvolution(cIndx) == false){	continue;	}

		//クラスタが含んでいる各粒子の位置・速度を取得
		//各位置・速度をstaticな最終位置に足す
		for(int oIndx = 0; oIndx < m_iceSM[cIndx]->GetIndxNum(); oIndx++)
		{
			int pIndx = m_iceSM[cIndx]->GetParticleIndx(oIndx);
			if(pIndx == MAXINT){	continue;	}

			float defAmount = pow(m_iceSM[cIndx]->GetDefAmount(), m_kernelDegree);

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
		int smIndx = i*SM_DIM;

		float clusterNum = (float)deformationSum[i];
		if(clusterNum == 0){	continue;	}

		//固体の最終位置
		for(int dim = 0; dim < SM_DIM; dim++)
		{
			sldPos[smIndx+dim] /= clusterNum;
			sldVel[smIndx+dim] /= clusterNum;
		}
	}
}