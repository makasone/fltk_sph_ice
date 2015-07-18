#include "Ice_Convolution_Weight.h"

typedef Ice_Convolution_Weight ConvoWeight;

ConvoWeight::Ice_Convolution_Weight(const vector<Ice_SM*>& iceSM, IceStructure* iceStrct)
{
	m_iceSM = iceSM;
	m_iceStrct = iceStrct;
	m_kernelDegree = 1.0f;
}

ConvoWeight::Ice_Convolution_Weight(const vector<ElasticObj*>& elasticObj, const vector<OrientedParticle*>& particles, IceStructure* iceStrct)
{
	m_elasticObj = elasticObj;
	m_vOrientedPrtes = particles;
	m_iceStrct = iceStrct;
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
	////現在はm_iceStrctを用いずにsm法のデータm_iceSMだけで計算している
	//unsigned sm_particleNum = IceObject::GetParticleNum();
	//vector<float> deformationSum(sm_particleNum, 0.0f);

	//float* sldPos = Ice_SM::GetSldPosPointer();
	//float* sldVel = Ice_SM::GetSldVelPointer();

	//float* s_sphPrtPos = IceObject::GetSPHHostPosPointer();
	//float* s_sphPrtVel = IceObject::GetSPHHostVelPointer();

	////sldの初期化
	////Ice_SM::ResetFinalParamPointer(sm_particleNum);
	//float tempPos[7000 * 3] = {};
	//float tempVel[7000 * 3] = {};

	//for(int cIndx = 0; cIndx < sm_particleNum; cIndx++){
	//	//最終結果算出に用いるクラスタの判定
	//	if(m_iceJudge->JudgeConvolution(cIndx) == false){	continue;	}

	//	//クラスタが含んでいる各粒子の位置・速度を取得
	//	//各位置・速度をstaticな最終位置に足す
	//	for(int oIndx = 0; oIndx < m_iceSM[cIndx]->GetIndxNum(); oIndx++)
	//	{
	//		int pIndx = m_iceSM[cIndx]->GetParticleIndx(oIndx);
	//		if(pIndx == MAXINT){	continue;	}

	//		//重い処理みたい
	//		//float defAmount = pow(m_iceSM[cIndx]->GetDefAmount(), m_kernelDegree);
	//		//for(int i = 0; i < m_kernelDegree; i++){
	//		//	defAmount *= m_iceSM[cIndx]->GetDefAmount();
	//		//}
	//		float defAmount = m_iceSM[cIndx]->GetDefAmount();
	//		
	//		Vec3& pos = m_iceSM[cIndx]->GetVertexPos(oIndx) * defAmount;
	//		Vec3& vel = m_iceSM[cIndx]->GetVertexVel(oIndx) * defAmount;

	//		int sldIndx = pIndx * SM_DIM;
	//		for(int dim = 0; dim < SM_DIM; dim++){
	//			tempPos[sldIndx+dim] += pos[dim];
	//			tempVel[sldIndx+dim] += vel[dim];
	//		}

	//		//総変形量をカウント
	//		deformationSum[pIndx] += defAmount;
	//	}
	//}

	////TODO::固体と液体で線形補完していないのに注意
	////足した数で平均し，粒子位置に反映
	//for(int i = 0; i < sm_particleNum; i++)
	//{
	//	int smIndx = i*SM_DIM;

	//	float clusterNum = deformationSum[i];
	//	if(clusterNum <= 0.0f){	continue;	}

	//	//固体の最終位置
	//	for(int dim = 0; dim < SM_DIM; dim++){
	//		sldPos[smIndx+dim] = tempPos[smIndx+dim]/clusterNum;
	//		sldVel[smIndx+dim] = tempVel[smIndx+dim]/clusterNum;
	//	}
	//}

	unsigned pNum = IceObject::GetParticleNum();

	//最終位置決定
	//現在はm_iceStrctを用いずにsm法のデータm_iceSMだけで計算している
	vector<unsigned> addParticleNum(pNum, 0);
	vector<float> defSum(pNum, 0.0f);
	vector<Vec3> prtPoses(pNum, Vec3(0.0f));

	//単純な足しあわせ
	for(int cIndx = 0; cIndx < pNum; cIndx++){

		///最終結果算出に用いるクラスタの判定
		if(m_iceJudge->JudgeConvolution(cIndx) == false){	continue;	}

		for(int oIndx = 0; oIndx < m_elasticObj[cIndx]->GetIndxNum(); oIndx++){
			int pIndx = m_elasticObj[cIndx]->GetParticleIndx(oIndx);
			if(pIndx == MAXINT){	continue;	}

			Vec3 pos = m_elasticObj[cIndx]->GetVertexPos(oIndx);
			float defAmount = m_elasticObj[cIndx]->DefAmount();
			//float defAmount = m_elasticObj[cIndx]->DefAmount() * m_orientedObj[cIndx]->DefAmount();
			prtPoses[pIndx] += pos * defAmount;

			//総変形量をカウント
			defSum[pIndx] += defAmount;
		}
	}

	//平均値を固体位置に反映
	for(int pIndx = 0; pIndx < pNum; pIndx++){
		int smIndx = pIndx*SM_DIM;

		//CtoPNum == PtoCNumより
		float defAmount = defSum[pIndx];
		if(defAmount <= 0.0f){	continue;	}

		Vec3 pos = prtPoses[pIndx] / defAmount;
		m_vOrientedPrtes[pIndx]->PrdPos(pos);
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