#include "Ice_Convolution_Normal.h"

typedef Ice_Convolution_Normal ConvoNormal;

ConvoNormal::Ice_Convolution_Normal(const vector<Ice_SM*>& iceSM, IceStructure* iceStrct)
{
	m_iceSM = iceSM;
	m_iceStrct = iceStrct;
}

ConvoNormal::~Ice_Convolution_Normal()
{
}

void ConvoNormal::SetConvoJudge(Ice_ConvoJudge* judge)
{
	m_iceJudge = judge;	
}

Ice_ConvoJudge* ConvoNormal::GetConvoJudge()
{
	return m_iceJudge;
}

void ConvoNormal::StepConvolution()
{
	//現在はm_iceStrctを用いずにsm法のデータm_iceSMだけで計算している
	unsigned pNum = IceObject::GetParticleNum();
	vector<unsigned> addParticleNum(pNum, 0);

	float* sldPos = Ice_SM::GetSldPosPointer();
	float* sldVel = Ice_SM::GetSldVelPointer();

	float* s_sphPrtPos = IceObject::GetSPHHostPosPointer();
	float* s_sphPrtVel = IceObject::GetSPHHostVelPointer();

	//sldの初期化
	Ice_SM::ResetFinalParamPointer(pNum);

	//単純な足しあわせ
	for(int cIndx = 0; cIndx < pNum; cIndx++)
	{
		//最終結果算出に用いるクラスタの判定
		if(m_iceJudge->JudgeConvolution(cIndx) == false){	continue;	}

		//クラスタが含んでいる各粒子の位置・速度をstaticな最終位置に足す
		Vec3 pos, vel;
		int pIndx, dim;

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

			//粒子数のカウント
			addParticleNum[pIndx] += 1;
		}
	}

	//平均値を固体位置に反映
	for(int pIndx = 0; pIndx < pNum; pIndx++)
	{
		int sphIndx = pIndx*4;
		int smIndx = pIndx*SM_DIM;

		//CtoPNum == PtoCNumより
		float clusterNum = (float)addParticleNum[pIndx];
		if(clusterNum <= 0.0f){	continue;	}

		//固体の最終位置
		for(int dim = 0; dim < SM_DIM; dim++)
		{
			sldPos[smIndx+dim] /= clusterNum;
			sldVel[smIndx+dim] /= clusterNum;
		}
	}
}

//テストコード　PtoCを毎フレーム作ってみたり
void ConvoNormal::StepConvolution2()
{
	unsigned pNum = IceObject::GetParticleNum();
	vector<unsigned> addParticleNum(pNum, 0);

	float* sldPos = Ice_SM::GetSldPosPointer();
	float* sldVel = Ice_SM::GetSldVelPointer();

	float* s_sphPrtPos = IceObject::GetSPHHostPosPointer();
	float* s_sphPrtVel = IceObject::GetSPHHostVelPointer();

	//並列処理に用いる変数
	Vec3 pos, vel;
	int oIndx, cIndx, coIndx, dim;

	//PtoCを作成
	vector<vector<pair<unsigned, unsigned>>> PtoC(pNum);

	//サイズ確保
	for(int pIndx = 0; pIndx < pNum; pIndx++)
	{
		PtoC[pIndx].resize(m_iceSM[pIndx]->GetNumVertices());
	}

	//PtoC=CtoPを利用し，クラスタ内で目当ての粒子を探索　並列可能
	#pragma omp parallel for private(oIndx, cIndx, coIndx)
	for(int pIndx = 0; pIndx < pNum; pIndx++)
	{
		if(m_iceJudge->JudgeConvolution(pIndx) == false){	continue;	}

		for(oIndx = 0; oIndx < m_iceSM[pIndx]->GetIndxNum(); oIndx++)
		{
			//pIndx番目のクラスタに含まれるcIndxという数字は，
			//クラスタに含まれる粒子であり近傍クラスタも示している．
			cIndx = m_iceSM[pIndx]->GetParticleIndx(oIndx);
			if(cIndx == MAXINT){	continue;	}

			//近傍クラスタの中からpIndx番目の粒子を探す
			coIndx = m_iceSM[cIndx]->SearchIndx(pIndx);
			if(coIndx == MAXINT){	continue;	}
			
			PtoC[pIndx][oIndx] = pair<unsigned, unsigned>(cIndx, coIndx);
		}
	}

	//sldの初期化
	Ice_SM::ResetFinalParamPointer(pNum);

	//PtoCの情報から各粒子の最終位置を決定
	#pragma omp parallel for private(oIndx, cIndx, coIndx, dim, pos, vel)
	for(int pIndx = 0; pIndx < pNum; pIndx++)
	{
		//最終結果算出に用いるクラスタの判定
		if(m_iceJudge->JudgeConvolution(pIndx) == false){	continue;	}

		//単純な平均
		for(oIndx = 0; oIndx < PtoC[pIndx].size(); oIndx++)
		{
			cIndx = PtoC[pIndx][oIndx].first;
			coIndx = PtoC[pIndx][oIndx].second;

			Vec3 pos = m_iceSM[cIndx]->GetVertexPos(coIndx);
			Vec3 vel = m_iceSM[cIndx]->GetVertexVel(coIndx);

			for(dim = 0; dim < SM_DIM; dim++)
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
	for(int pIndx = 0; pIndx < pNum; pIndx++)
	{
		int smIndx = pIndx*SM_DIM;

		float clusterNum = (float)addParticleNum[pIndx];
		if(clusterNum == 0){	continue;	}

		//粒子位置に反映
		for(int dim = 0; dim < SM_DIM; dim++)
		{
			sldPos[pIndx+dim] /= clusterNum;
			sldVel[pIndx+dim] /= clusterNum;
		}
	}
}

void ConvoNormal::StepConvolutionDebug()
{
	unsigned pNum = IceObject::GetParticleNum();
	vector<unsigned> addParticleNum(pNum, 0);

	float* sldPos = Ice_SM::GetSldPosPointer();
	float* sldVel = Ice_SM::GetSldVelPointer();

	float* s_sphPrtPos = IceObject::GetSPHHostPosPointer();
	float* s_sphPrtVel = IceObject::GetSPHHostVelPointer();

	//sldの初期化
	Ice_SM::ResetFinalParamPointer(pNum);

	//単純な足しあわせ
	for(int cIndx = 0; cIndx < pNum; cIndx++)
	{
		//最終結果算出に用いるクラスタの判定
		if(m_iceJudge->JudgeConvolution(cIndx) == false){	continue;	}

		//クラスタが含んでいる各粒子の位置・速度をstaticな最終位置に足す
		Vec3 pos, vel;
		int pIndx, dim;

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

			//粒子数のカウント
			addParticleNum[pIndx] += 1;
		}
	}

	//平均値を固体位置に反映
	for(int pIndx = 0; pIndx < pNum; pIndx++)
	{
		int sphIndx = pIndx*4;
		int smIndx = pIndx*SM_DIM;

		//CtoPNum == PtoCNumより
		float clusterNum = (float)addParticleNum[pIndx];
		if(clusterNum <= 0.0f){	continue;	}

		//固体の最終位置
		for(int dim = 0; dim < SM_DIM; dim++)
		{
			sldPos[smIndx+dim] /= clusterNum;
			sldVel[smIndx+dim] /= clusterNum;
		}
	}
}