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
	//���݂�m_iceStrct��p������sm�@�̃f�[�^m_iceSM�����Ōv�Z���Ă���
	unsigned pNum = IceObject::GetParticleNum();
	vector<unsigned> addParticleNum(pNum, 0);

	float* sldPos = Ice_SM::GetSldPosPointer();
	float* sldVel = Ice_SM::GetSldVelPointer();

	float* s_sphPrtPos = IceObject::GetSPHHostPosPointer();
	float* s_sphPrtVel = IceObject::GetSPHHostVelPointer();

	//sld�̏�����
	Ice_SM::ResetFinalParamPointer(pNum);

	//�P���ȑ������킹
	for(int cIndx = 0; cIndx < pNum; cIndx++)
	{
		//�ŏI���ʎZ�o�ɗp����N���X�^�̔���
		if(m_iceJudge->JudgeConvolution(cIndx) == false){	continue;	}

		//�N���X�^���܂�ł���e���q�̈ʒu�E���x��static�ȍŏI�ʒu�ɑ���
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

			//���q���̃J�E���g
			addParticleNum[pIndx] += 1;
		}
	}

	//���ϒl���ő̈ʒu�ɔ��f
	for(int pIndx = 0; pIndx < pNum; pIndx++)
	{
		int sphIndx = pIndx*4;
		int smIndx = pIndx*SM_DIM;

		//CtoPNum == PtoCNum���
		float clusterNum = (float)addParticleNum[pIndx];
		if(clusterNum <= 0.0f){	continue;	}

		//�ő̂̍ŏI�ʒu
		for(int dim = 0; dim < SM_DIM; dim++)
		{
			sldPos[smIndx+dim] /= clusterNum;
			sldVel[smIndx+dim] /= clusterNum;
		}
	}
}

//�e�X�g�R�[�h�@PtoC�𖈃t���[������Ă݂���
void ConvoNormal::StepConvolution2()
{
	unsigned pNum = IceObject::GetParticleNum();
	vector<unsigned> addParticleNum(pNum, 0);

	float* sldPos = Ice_SM::GetSldPosPointer();
	float* sldVel = Ice_SM::GetSldVelPointer();

	float* s_sphPrtPos = IceObject::GetSPHHostPosPointer();
	float* s_sphPrtVel = IceObject::GetSPHHostVelPointer();

	//���񏈗��ɗp����ϐ�
	Vec3 pos, vel;
	int oIndx, cIndx, coIndx, dim;

	//PtoC���쐬
	vector<vector<pair<unsigned, unsigned>>> PtoC(pNum);

	//�T�C�Y�m��
	for(int pIndx = 0; pIndx < pNum; pIndx++)
	{
		PtoC[pIndx].resize(m_iceSM[pIndx]->GetNumVertices());
	}

	//PtoC=CtoP�𗘗p���C�N���X�^���Ŗړ��Ă̗��q��T���@����\
	#pragma omp parallel for private(oIndx, cIndx, coIndx)
	for(int pIndx = 0; pIndx < pNum; pIndx++)
	{
		if(m_iceJudge->JudgeConvolution(pIndx) == false){	continue;	}

		for(oIndx = 0; oIndx < m_iceSM[pIndx]->GetIndxNum(); oIndx++)
		{
			//pIndx�Ԗڂ̃N���X�^�Ɋ܂܂��cIndx�Ƃ��������́C
			//�N���X�^�Ɋ܂܂�闱�q�ł���ߖT�N���X�^�������Ă���D
			cIndx = m_iceSM[pIndx]->GetParticleIndx(oIndx);
			if(cIndx == MAXINT){	continue;	}

			//�ߖT�N���X�^�̒�����pIndx�Ԗڂ̗��q��T��
			coIndx = m_iceSM[cIndx]->SearchIndx(pIndx);
			if(coIndx == MAXINT){	continue;	}
			
			PtoC[pIndx][oIndx] = pair<unsigned, unsigned>(cIndx, coIndx);
		}
	}

	//sld�̏�����
	Ice_SM::ResetFinalParamPointer(pNum);

	//PtoC�̏�񂩂�e���q�̍ŏI�ʒu������
	#pragma omp parallel for private(oIndx, cIndx, coIndx, dim, pos, vel)
	for(int pIndx = 0; pIndx < pNum; pIndx++)
	{
		//�ŏI���ʎZ�o�ɗp����N���X�^�̔���
		if(m_iceJudge->JudgeConvolution(pIndx) == false){	continue;	}

		//�P���ȕ���
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

			//�����������J�E���g
			addParticleNum[pIndx] += 1;
		}
	}

	//TODO::�ő̂Ɖt�̂Ő��`�⊮���Ă��Ȃ��̂ɒ���
	//���������ŕ��ς��C���q�ʒu�ɔ��f
	for(int pIndx = 0; pIndx < pNum; pIndx++)
	{
		int smIndx = pIndx*SM_DIM;

		float clusterNum = (float)addParticleNum[pIndx];
		if(clusterNum == 0){	continue;	}

		//���q�ʒu�ɔ��f
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

	//sld�̏�����
	Ice_SM::ResetFinalParamPointer(pNum);

	//�P���ȑ������킹
	for(int cIndx = 0; cIndx < pNum; cIndx++)
	{
		//�ŏI���ʎZ�o�ɗp����N���X�^�̔���
		if(m_iceJudge->JudgeConvolution(cIndx) == false){	continue;	}

		//�N���X�^���܂�ł���e���q�̈ʒu�E���x��static�ȍŏI�ʒu�ɑ���
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

			//���q���̃J�E���g
			addParticleNum[pIndx] += 1;
		}
	}

	//���ϒl���ő̈ʒu�ɔ��f
	for(int pIndx = 0; pIndx < pNum; pIndx++)
	{
		int sphIndx = pIndx*4;
		int smIndx = pIndx*SM_DIM;

		//CtoPNum == PtoCNum���
		float clusterNum = (float)addParticleNum[pIndx];
		if(clusterNum <= 0.0f){	continue;	}

		//�ő̂̍ŏI�ʒu
		for(int dim = 0; dim < SM_DIM; dim++)
		{
			sldPos[smIndx+dim] /= clusterNum;
			sldVel[smIndx+dim] /= clusterNum;
		}
	}
}