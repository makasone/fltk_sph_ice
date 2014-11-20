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
	//����
	m_iceMove->StepObjMove();		//���̂܂܌ĂԂ���
	m_iceInterPolation->StepInterPolationItr();

	//����
	for(int itr = 0; itr < Ice_SM::GetIteration(); itr++)
	{
		m_iceMove->StepObjMoveItr();
		m_iceInterPolation->StepInterPolationItr();
	}

	//���x�Z�o
	#pragma omp parallel for
	for(int i = 0; i < IceObject::GetParticleNum(); ++i)
	{	
		if(m_iceMove->GetJudgeMove()->JudgeMove(i) == false){	continue;	}
		m_iceSM[i]->integrateIteration();
	}
}

//Weight�ȕ�Ԃ����Ă��Ȃ��̂ɒ���
void CalcIteration::InterPolationForCluster()
{
	unsigned sm_particleNum = IceObject::GetParticleNum();
	vector<unsigned> addParticleNum(sm_particleNum, 0);

	float* sldPos = Ice_SM::GetSldPosPointer();
	float* sldVel = Ice_SM::GetSldVelPointer();

	//sld�̏�����
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
		//�ŏI���ʎZ�o�ɗp����N���X�^�̔���
		if(m_iceInterPolation->GetIntrpJudge()->JudgeInterPolation(cIndx) == false){	continue;	}

		//�N���X�^���܂�ł���e���q�̈ʒu�E���x���擾
		//�e�ʒu�E���x��static�ȍŏI�ʒu�ɑ���
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

			//�����������J�E���g
			addParticleNum[pIndx] += 1;
		}
	}

	//TODO::�ő̂Ɖt�̂Ő��`�⊮���Ă��Ȃ��̂ɒ���
	//���������ŕ��ς��C���q�ʒu�ɔ��f
	for(int i = 0; i < sm_particleNum; i++)
	{
		int pIndx = i*SM_DIM;
		int smIndx = i*SM_DIM;

		float clusterNum = (float)addParticleNum[i];
		if(clusterNum == 0){	continue;	}

		//���q�ʒu�ɔ��f
		for(int dim = 0; dim < SM_DIM; dim++)
		{
			sldPos[pIndx+dim] = sldPos[smIndx+dim]/clusterNum;
			sldVel[pIndx+dim] = sldVel[smIndx+dim]/clusterNum;
		}
	}
}