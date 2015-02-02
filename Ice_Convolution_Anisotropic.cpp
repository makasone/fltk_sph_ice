#include "Ice_Convolution_Anisotropic.h"

typedef Ice_Convolution_Anisotropic ConvoAnisotropict;

ConvoAnisotropict::Ice_Convolution_Anisotropic(const vector<Ice_SM*>& iceSM, IceStructure* iceStrct)
{
	m_iceSM = iceSM;
	m_iceStrct = iceStrct;
	m_dirVec = Vec3(0.0);
}

ConvoAnisotropict::~Ice_Convolution_Anisotropic()
{
}

void ConvoAnisotropict::SetConvoJudge(Ice_ConvoJudge* judge)
{
	m_iceJudge = judge;	
}

Ice_ConvoJudge* ConvoAnisotropict::GetConvoJudge()
{
	return m_iceJudge;
}

void ConvoAnisotropict::StepConvolution()
{
	float* sldPos = Ice_SM::GetSldPosPointer();
	float* sldVel = Ice_SM::GetSldVelPointer();

	float* s_sphPrtPos = IceObject::GetSPHHostPosPointer();
	float* s_sphPrtVel = IceObject::GetSPHHostVelPointer();

	unsigned sm_particleNum = IceObject::GetParticleNum();
	vector<float> similaritySum(sm_particleNum, 0.0f);
	
	int side = 17;
	m_dirVec = Vec3(0.0f, -1.0f, -1.0f);
	m_dirVec = Unit(m_dirVec);

	//前フレームの位置情報
	float* pre_sldPos = new float[sm_particleNum * SM_DIM];
	for(int pIndx = 0; pIndx < sm_particleNum; pIndx++)
	{
		for(int dim = 0; dim < SM_DIM; dim++)
		{
			pre_sldPos[pIndx*SM_DIM+dim] = sldPos[pIndx*SM_DIM+dim];
		}
	}

	//sldの初期化
	Ice_SM::ResetFinalParamPointer(sm_particleNum);

	for(int cIndx = 0; cIndx < sm_particleNum; cIndx++)
	{
		//最終結果算出に用いるクラスタの判定
		if(m_iceJudge->JudgeConvolution(cIndx) == false){	continue;	}

		Vec3 orgCm = m_iceSM[cIndx]->GetOrgCm();

		//クラスタが含んでいる各粒子の位置・速度を取得
		//各位置・速度をstaticな最終位置に足す
		for(int oIndx = 0; oIndx < m_iceSM[cIndx]->GetIndxNum(); oIndx++)
		{
			int pIndx = m_iceSM[cIndx]->GetParticleIndx(oIndx);
			if(pIndx == MAXINT){	continue;	}

			//ijiriら2012の手法
			float mass = 1.0f;
			Vec3 prtPos = m_iceSM[cIndx]->GetOrgPos(oIndx);
			Vec3 prtVec = Unit(prtPos-orgCm);
			float similality = mass * (dot(m_dirVec, prtVec) * dot(m_dirVec, prtVec) + 0.01f);

			Vec3 pos = m_iceSM[cIndx]->GetVertexPos(oIndx) * similality;
			Vec3 vel = m_iceSM[cIndx]->GetVertexVel(oIndx) * similality;

			for(int dim = 0; dim < SM_DIM; dim++)
			{
				sldPos[pIndx*SM_DIM+dim] += pos[dim];
				sldVel[pIndx*SM_DIM+dim] += vel[dim];
			}

			//総変形量をカウント
			similaritySum[pIndx] += similality;
		}
	}

	float timeStep = 0.01f;
	Vec3 gravity(0.0, -9.81, 0.0);

	//TODO::固体と液体で線形補完していないのに注意
	//足した数で平均し，粒子位置に反映
	for(int i = 0; i < sm_particleNum; i++)
	{
		int smIndx = i*SM_DIM;

		if(similaritySum[i] <= 0.0f){	continue;	}

		//固体の最終位置・速度
		for(int dim = 0; dim < SM_DIM; dim++)
		{
			sldPos[smIndx+dim] /= similaritySum[i];
			sldVel[smIndx+dim] /= similaritySum[i];

			//前フレームからの差分で速度更新
			//sldVel[smIndx+dim] += (sldPos[smIndx+dim] - pre_sldPos[smIndx+dim]) / timeStep + gravity[dim] * timeStep;

			//最終位置
			//sldPos[smIndx+dim] += timeStep * sldVel[smIndx+dim];
		}
	}

	delete[] pre_sldPos;
}

void ConvoAnisotropict::StepConvolutionDebug()
{

}