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

	//sld�̏�����
	Ice_SM::ResetFinalParamPointer(sm_particleNum);

	for(int cIndx = 0; cIndx < sm_particleNum; cIndx++)
	{
		//�ŏI���ʎZ�o�ɗp����N���X�^�̔���
		if(m_iceJudge->JudgeConvolution(cIndx) == false){	continue;	}

		//�N���X�^���܂�ł���e���q�̈ʒu�E���x���擾
		//�e�ʒu�E���x��static�ȍŏI�ʒu�ɑ���
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

			//���ό`�ʂ��J�E���g
			deformationSum[pIndx] += defAmount;
		}
	}

	//TODO::�ő̂Ɖt�̂Ő��`�⊮���Ă��Ȃ��̂ɒ���
	//���������ŕ��ς��C���q�ʒu�ɔ��f
	for(int i = 0; i < sm_particleNum; i++)
	{
		int smIndx = i*SM_DIM;

		float clusterNum = (float)deformationSum[i];
		if(clusterNum == 0){	continue;	}

		//�ő̂̍ŏI�ʒu
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

	//sld�̏�����
	Ice_SM::ResetFinalParamPointer(sm_particleNum);

	for(int cIndx = 0; cIndx < sm_particleNum; cIndx++)
	{
		//�ŏI���ʎZ�o�ɗp����N���X�^�̔���
		if(m_iceJudge->JudgeConvolution(cIndx) == false){	continue;	}

		//�N���X�^���܂�ł���e���q�̈ʒu�E���x���擾
		//�e�ʒu�E���x��static�ȍŏI�ʒu�ɑ���
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

			//���ό`�ʂ��J�E���g
			deformationSum[pIndx] += defAmount;
		}
	}

	//TODO::�ő̂Ɖt�̂Ő��`�⊮���Ă��Ȃ��̂ɒ���
	//���������ŕ��ς��C���q�ʒu�ɔ��f
	for(int i = 0; i < sm_particleNum; i++)
	{
		int smIndx = i*SM_DIM;

		float clusterNum = (float)deformationSum[i];
		if(clusterNum == 0){	continue;	}

		//�ő̂̍ŏI�ʒu
		for(int dim = 0; dim < SM_DIM; dim++)
		{
			sldPos[smIndx+dim] /= clusterNum;
			sldVel[smIndx+dim] /= clusterNum;
		}
	}
}