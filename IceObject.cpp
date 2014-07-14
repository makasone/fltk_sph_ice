#include "IceObject.h"

float* IceObject::s_sphPrtPos;
float* IceObject::s_sphPrtVel;

float* IceObject::m_fInterPolationCoefficience;

//�f�o�C�X�|�C���^
float* IceObject::sd_sldPrtPos;	
float* IceObject::sd_sldPrtVel;

float* IceObject::sd_ObjPrtPos;
float* IceObject::sd_ObjPrtVel;

int IceObject::sm_particleNum;
int IceObject::sm_tetraNum;
int IceObject::sm_clusterNum;


IceObject::IceObject(float* pos, float* vel, int pMaxNum, int cMaxNum, int tMaxNum)
{
	s_sphPrtPos = pos;
	s_sphPrtVel = vel;

	InitIceObj(pMaxNum, cMaxNum, tMaxNum);

}

IceObject::~IceObject()
{
}

//���ꂼ��̃N���X�E�ϐ��̏�����
void IceObject::InitIceObj(int pMaxNum, int cMaxNum, int tMaxNum)
{	cout << __FUNCTION__ << endl;
	//���̂̍\���̏�����
	m_iceStrct = new IceStructure(pMaxNum, cMaxNum, tMaxNum);

	//�^���v�Z���s��SM�N���X�^�̏�����

	//��ԏ����̂��߂̃p�����[�^�̏�����
	InitInterPolation();
}

//��ԏ����̂��߂̃p�����[�^�̏�����
void IceObject::InitInterPolation()
{
	m_fInterPolationCoefficience = new float[sm_particleNum];	//���`��ԌW��
	cout << __FUNCTION__ << " pnum = " << sm_particleNum << endl;

	for(int i = 0; i < sm_particleNum; ++i)
	{
		m_fInterPolationCoefficience[i] = 1.0f;
	}
}

//�ő̂̉^���v�Z
void IceObject::StepObjMove()
{
	//GPU��p�����N���X�^�̉^���v�Z
	//Ice_SM::UpdateGPU();

	//TODO::������GPU�ŏ�������D
	for(int i = 0; i < sm_particleNum; i++)
	{
		if(GetPtoCNum(i) == 0){	continue;	}
		m_iceMove[i]->CopyDeviceToInstance(i);
	}
}

//���̂ƌő̂̍ŏI�I�ȉ^���v�Z
void IceObject::StepInterPolation()
{
	Vec3 pos,vel;

	for(int i = 0; i < sm_particleNum; ++i)
	{
		if(GetParticleNum() <= i){	continue;	}	//�Z���݂̂̎����̂Ƃ��ɕK�v�ɂȂ�D
		if(GetPtoCNum(i) <= 0){		continue;	}

		//�ő̉^���̍ŏI�ʒu�v�Z
		CalcAverageCPU(i, pos, vel);

		//�t�̂ƌő̂̕��
		LinerInterPolationCPU(i, pos, vel);
	}
}

void IceObject::CalcAverageCPU(const int pIndx, Vec3& pos, Vec3& vel)
{
		//���ꂼ��̃x�N�g�������������ς��Ƃ�
		pos = Vec3(0.0, 0.0, 0.0);
		vel = Vec3(0.0, 0.0, 0.0);
		double shapeNum = 0.0;		//�N���X�^�̐�

		for(int j = 0; j < GetPtoCIndx(pIndx); ++j)
		{
			int jcIndx = GetPtoC(pIndx, j, 0);
			int joIndx = GetPtoC(pIndx, j, 1);

			if(jcIndx == -1 || joIndx == -1){	continue;	}

			pos += m_iceMove[jcIndx]->GetVertexPos(joIndx);
			vel += m_iceMove[jcIndx]->GetVertexVel(joIndx);

			shapeNum += 1.0;
		}

		int jpIndx = pIndx*4;
		

		//�N���X�^�̐��Ŋ���
		//�ǂ̃N���X�^�ɂ��܂܂�Ă��Ȃ��ꍇ�C�^����SPH�@�ɏ]��
		if(shapeNum != 0.0)
		{
			pos /= shapeNum;
			vel /= shapeNum;
		}		
		else
		{
			pos = Vec3(s_sphPrtPos[jpIndx+0], s_sphPrtPos[jpIndx+1], s_sphPrtPos[jpIndx+2]);
			vel = Vec3(s_sphPrtVel[jpIndx+0], s_sphPrtVel[jpIndx+1], s_sphPrtVel[jpIndx+2]);
		}
}

//SPH�@��SM�@�ŋ��߂����x�ƈʒu����`��� CPU
void IceObject::LinerInterPolationCPU(const int pIndx, const Vec3& pos, const Vec3& vel)
{
	int sphIndx = pIndx*4;
	double intrps = 1.0-m_fInterPolationCoefficience[pIndx];	//��ԌW��

	s_sphPrtVel[sphIndx+0] = vel[0] * m_fInterPolationCoefficience[pIndx] + s_sphPrtVel[sphIndx+0] * intrps;
	s_sphPrtVel[sphIndx+1] = vel[1] * m_fInterPolationCoefficience[pIndx] + s_sphPrtVel[sphIndx+1] * intrps;
	s_sphPrtVel[sphIndx+2] = vel[2] * m_fInterPolationCoefficience[pIndx] + s_sphPrtVel[sphIndx+2] * intrps;

	s_sphPrtPos[sphIndx+0] = pos[0] * m_fInterPolationCoefficience[pIndx] + s_sphPrtPos[sphIndx+0] * intrps;
	s_sphPrtPos[sphIndx+1] = pos[1] * m_fInterPolationCoefficience[pIndx] + s_sphPrtPos[sphIndx+1] * intrps;
	s_sphPrtPos[sphIndx+2] = pos[2] * m_fInterPolationCoefficience[pIndx] + s_sphPrtPos[sphIndx+2] * intrps;
}