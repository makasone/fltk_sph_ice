#include "IceObject.h"

float* IceObject::s_sphPrtPos;
float* IceObject::s_sphPrtVel;

float* IceObject::m_fInterPolationCoefficience;

//�f�o�C�X�|�C���^
float* IceObject::sd_sphPrtPos;
float* IceObject::sd_sphPrtVel;

float* IceObject::sd_sldPrtPos;	
float* IceObject::sd_sldPrtVel;

int IceObject::sm_particleNum;
int IceObject::sm_tetraNum;
int IceObject::sm_clusterNum;


IceObject::IceObject(float* pos, float* vel, int pMaxNum, int cMaxNum, int tMaxNum)
{
	//���̂Ƃ��날��܂�Ӗ��Ȃ��݂���
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

	for(int i = 0; i < sm_particleNum; ++i)
	{
		m_fInterPolationCoefficience[i] = 1.0f;
	}
}

//GPU�����ŗp����ϐ��̏�����
void IceObject::InitGPU()
{
	//�ő̍\���̏�����
	m_iceStrct->InitGPU();

	//TODO::�N���X�^��GPU�������������ɒu��

	//�e���q�̍ŏI�I�Ȉʒu�E���x�f�[�^
	cudaMalloc((void**)&sd_sldPrtPos,	sizeof(float) * MAXCLUSTER * SM_DIM);
	cudaMalloc((void**)&sd_sldPrtVel,	sizeof(float) * MAXCLUSTER * SM_DIM);

	//�ŏI�ʒu�E���x�����݂̃f�[�^�ŏ�����
	float* fPoses = new float[MAXCLUSTER * SM_DIM];
	float* fVeles = new float[MAXCLUSTER * SM_DIM];

	//s_pfPrtPos�Ȃǂ̓f�[�^�̒��g��DIM=4�ō���Ă���̂ŁC�������Ȃ��Ƃ����Ȃ�
	//TODO::���q�T�C�Y���傫���Ȃ�ƁC���������m�ۂł��Ȃ���������Ȃ��̂ɒ���
	int sphDIM = 4;
	for(int i = 0; i < MAXCLUSTER; ++i)
	{
		for(int j = 0; j < SM_DIM; ++j)
		{
			fPoses[i*SM_DIM+j] = s_sphPrtPos[i*sphDIM+j];
			fVeles[i*SM_DIM+j] = s_sphPrtVel[i*sphDIM+j];
		}
	}

	cudaMemcpy(sd_sldPrtPos, fPoses, sizeof(float) * MAXCLUSTER * SM_DIM, cudaMemcpyHostToDevice);
	cudaMemcpy(sd_sldPrtVel, fVeles, sizeof(float) * MAXCLUSTER * SM_DIM, cudaMemcpyHostToDevice);

	////�������̓]�������܂����������̊m�F
	////�ꎞ�z��̃��Z�b�g
	//for(int i = 0; i < MAXCLUSTER; ++i)
	//{
	//	for(int j = 0; j < SM_DIM; ++j)
	//	{
	//		fPoses[i*SM_DIM+j] = 0.0f;
	//		fVeles[i*SM_DIM+j] = 0.0f;
	//	}
	//}

	////�f�[�^��]��
	//cudaMemcpy(fPoses, d_FinalPos, sizeof(float) * MAXCLUSTER * SM_DIM, cudaMemcpyDeviceToHost);
	//cudaMemcpy(fVeles, d_FinalVel, sizeof(float) * MAXCLUSTER * SM_DIM, cudaMemcpyDeviceToHost);

	////�z�X�g���̃f�[�^��]���������ʂ��_���v
	//ofstream ofs( "DtoH_Test.txt" );
	//ofs << "DtoH_Test" << endl;
	//
	////�f�o�C�X���̃f�[�^��]��
	//for(int i = 0; i < MAXCLUSTER; i++)
	//{
	//	ofs << "particle" << i << " pos::(" << fPoses[i*SM_DIM+0] << ", " << fPoses[i*SM_DIM+1] << ", " << fPoses[i*SM_DIM+2] << ")" << endl;
	//	ofs << "sph     " << i << " pos::(" << s_pfPrtPos[i*sphDIM+0] << ", " << s_pfPrtPos[i*sphDIM+1] << ", " << s_pfPrtPos[i*sphDIM+2] << ")" << endl;
	//	ofs << "particle" << i << " vel::(" << fVeles[i*SM_DIM+0] << ", " << fVeles[i*SM_DIM+1] << ", " << fVeles[i*SM_DIM+2] << ")" << endl;
	//	ofs << "sph     " << i << " vel::(" << s_pfPrtVel[i*sphDIM+0] << ", " << s_pfPrtVel[i*sphDIM+1] << ", " << s_pfPrtVel[i*sphDIM+2] << ")" << endl;
	//}

	delete[] fPoses;
	delete[] fVeles;
}

//�ő̂̉^���v�Z
void IceObject::StepObjMove()
{
	//GPU��p�����N���X�^�̉^���v�Z
	//Ice_SM::UpdateGPU();

	////TODO::�������Ȃ����D
	//for(int i = 0; i < sm_particleNum; i++)
	//{
	//	if(GetPtoCNum(i) == 0){	continue;	}
	//	m_iceMove[i]->CopyDeviceToInstance(i);
	//}
}

//���̂ƌő̂̍ŏI�I�ȉ^���v�Z
void IceObject::StepInterPolation()
{
//CPU
	//for(int i = 0; i < sm_particleNum; ++i)
	//{
	//	if(GetParticleNum() <= i){	continue;	}	//�Z���݂̂̎����̂Ƃ��ɕK�v�ɂȂ�D
	//	if(GetPtoCNum(i) <= 0){		continue;	}
	
	//	Vec3 pos,vel;
	
	//	//�ő̉^���̍ŏI�ʒu�v�Z
	//	CalcAverageCPU(i, pos, vel);

	//	//�t�̂ƌő̂̕��
	//	LinerInterPolationCPU(i, pos, vel);	//���`���
	//}

//GPU
	//�ő̉^���̍ŏI�ʒu�v�Z
	float* smPrtPos = Ice_SM::GetDevicePosPointer();
	float* smPrtVel = Ice_SM::GetDeviceVelPointer();
	int* indxSet = Ice_SM::GetDeviceIndexSetPointer();

	sd_sphPrtPos = Ice_SM::GetDeviceSPHPosPointer();
	sd_sphPrtVel = Ice_SM::GetDeviceSPHVelPointer();

	int* PtoCIndx = IceStructure::GetDevicePtoCIndxPointer();
	int* PtoC = IceStructure::GetDevicePtoCPointer();
	int PNumMax = m_iceStrct->GetPNumMax();
	int PtoCMax = m_iceStrct->GetPtoCMax();
	int PtoCParamSize = 3;

	LaunchCalcAverageGPU(sd_sldPrtPos, sd_sldPrtVel, sd_sphPrtPos, sd_sphPrtVel, smPrtPos, smPrtVel, indxSet, PtoCIndx, PtoC, PNumMax, PtoCMax, PtoCParamSize);

	//�t�̂ƌő̂̕��
	LaunchInterPolationGPU();
}

//�e�N���X�^�̌v�Z���ʂ̕��ς��C�ő̂̍ŏI�I�ȉ^���v�Z���ʂƂ���
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

		//�N���X�^�̐��Ŋ���
		if(shapeNum != 0.0)
		{
			pos /= shapeNum;
			vel /= shapeNum;
		}		
		//�ǂ̃N���X�^�ɂ��܂܂�Ă��Ȃ��ꍇ�C�^����SPH�@�ɏ]��
		else
		{
			int jpIndx = pIndx*4;
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