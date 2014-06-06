#include <stdio.h>
#include <math.h>
#include <rx_cu_common.cuh>	//�g��Ȃ��H

#define SM_DIM 3

//�����n���ϐ������Ԉ���Ă��Ă��G���[���o�Ȃ��̂Œ���

void LaunchShapeMathcingGPU(float* curPos, float* vel, int* pIndxes, float dt, int prtNum);
__global__ void Update(float* curPos, float* vel, int* pIndxes, float dt, int prtNum);
__device__ void ExternalForce(float* curPos, float* vel, int* pIndxes, float dt, int prtNum);
__device__ void ProjectPos(float* curPos, float* vel, int* pIndxes, float dt, int prtNum);
__device__ void Integrate(float* curPos, float* vel, int* pIndxes, float dt, int prtNum);

float* d_Pos;	//�f�o�C�X���̗��q�ʒu
float* d_Vel;	//�f�o�C�X���̗��q���x

//�p�����[�^�̏�����
void InitParam()
{
}

//GPU����
void LaunchShapeMatchingGPU(float* curPos, float* vel, int* pIndxes, float dt, int prtNum)
{
	
	//printf("LaunchGPUKernel");

	dim3 grid(1, 1);
	dim3 block(1, 1, 1);

	//�^���v�Z
	Update <<< grid , block >>> (curPos, vel, pIndxes, dt, prtNum);
}


//GPU�̈ʒu�E���x�X�V
__global__
void Update(
	float* curPos, 
	float* vel,
	int* pIndxes,
	float dt,
	int prtNum)
{
	//printf("d_Integrate\n");	//�߂��Ⴍ����o��̂Œ���

	ExternalForce(curPos, vel, pIndxes, dt, prtNum);
	ProjectPos(curPos, vel, pIndxes, dt, prtNum);
	Integrate(curPos, vel, pIndxes, dt, prtNum);
}

__device__
	void ExternalForce(float* curPos, float* vel, int* pIndxes, float dt, int prtNum)
{
	// �d�͂̉e����t���C���x�𔽉f
	for(int i = 0; i < prtNum; ++i)
	{
		int pIndx = pIndxes[i]*4;
		int cIndx = i*SM_DIM;

		for(int j = 0; j < SM_DIM; j++)
		{
			int jpIndx = pIndx+j;
			int jcIndx = cIndx+j;

			//curPos[jcIndx] = d_Pos[jpIndx]+d_Vel[jpIndx]*dt;
		}
	}

	// ���E�ǂ̉e��
	//���������Ȃ�d���Ȃ邪�C����͂���݂���
	//double res = 0.9;	// �����W��
	//for(int i = 0; i < prtNum; ++i){
	//	//if(m_pFix[i]) continue;
	//	//Vec3 &p = m_vCurPos[i];
	//	//Vec3 &np = m_vNewPos[i];
	//	//Vec3 &v = m_vVel[i];
	//	//if(np[0] < m_v3Min[0] || np[0] > m_v3Max[0]){
	//	//	np[0] = p[0]-v[0]*dt*res;
	//	//	np[1] = p[1];
	//	//	np[2] = p[2];
	//	//}
	//	//if(np[1] < m_v3Min[1] || np[1] > m_v3Max[1]){
	//	//	np[1] = p[1]-v[1]*dt*res;
	//	//	np[0] = p[0] ;
	//	//	np[2] = p[2];
	//	//}
	//	//if(np[2] < m_v3Min[2] || np[2] > m_v3Max[2]){
	//	//	np[2] = p[2]-v[2]*dt*res;
	//	//	np[0] = p[0];
	//	//	np[1] = p[1];
	//	//}

	//	//clamp(curPos, i*SM_DIM);
	//}
}

__device__
	void ProjectPos(float* curPos, float* vel, int* pIndxes, float dt, int prtNum)
{
	if(prtNum <= 1) return;

	float3 cm = make_float3(0.0, 0.0, 0.0);
	float3 cm_org = make_float3(0.0, 0.0, 0.0);	// �d�S

	double mass = 0.0;	// ������

	// �d�S���W�̌v�Z
	for(int i = 0; i < prtNum;++i){
		//double m = m_pMass[i];
		double m = 1.0;
		//if(m_pFix[i]) m *= 300.0;	// �Œ�_�̎��ʂ�傫������
		
		int cIndx = i*SM_DIM;
		mass += m;

		cm.x += curPos[cIndx+0]*m;
		cm.y += curPos[cIndx+1]*m;
		cm.z += curPos[cIndx+2]*m;
	}

	cm.x /= mass;
	cm.y /= mass;
	cm.z /= mass;

	//cm_org = m_vec3OrgCm;

	//rxMatrix3 Apq(0.0), Aqq(0.0);
	//Vec3 p, q;

	//// Apq = ��mpq^T
	//// Aqq = ��mqq^T
	//for(int i = 0; i < prtNum; ++i)
	//{
	//	int cIndx = i*SM_DIM;

	//	for(int j = 0; j < SM_DIM; j++)
	//	{
	//		//p[j] = m_pNewPos[cIndx+j]-cm[j];

	//		//p[j] = newPos[cIndx+j]-cm[j];
	//		p[j] = curPos[cIndx+j]-cm[j];
	//		q[j] = m_pOrgPos[cIndx+j]-cm_org[j];
	//	}

	//	double m = m_pMass[i];

	//	Apq(0,0) += m*p[0]*q[0];
	//	Apq(0,1) += m*p[0]*q[1];
	//	Apq(0,2) += m*p[0]*q[2];
	//	Apq(1,0) += m*p[1]*q[0];
	//	Apq(1,1) += m*p[1]*q[1];
	//	Apq(1,2) += m*p[1]*q[2];
	//	Apq(2,0) += m*p[2]*q[0];
	//	Apq(2,1) += m*p[2]*q[1];
	//	Apq(2,2) += m*p[2]*q[2];

	//	Aqq(0,0) += m*q[0]*q[0];
	//	Aqq(0,1) += m*q[0]*q[1];
	//	Aqq(0,2) += m*q[0]*q[2];
	//	Aqq(1,0) += m*q[1]*q[0];
	//	Aqq(1,1) += m*q[1]*q[1];
	//	Aqq(1,2) += m*q[1]*q[2];
	//	Aqq(2,0) += m*q[2]*q[0];
	//	Aqq(2,1) += m*q[2]*q[1];
	//	Aqq(2,2) += m*q[2]*q[2];
	//}

	////Apq�̍s�񎮂����߁C���]���邩�𔻒�
	////�s����ȏꍇ�������̂Ł~
	////if( Apq.Determinant() < 0.0 && m_iNumVertices >= 4)
	////{
	//	//cout << "before det < 0" << endl;
	//	//�P�@�����𔽓]
	//	//Apq(0,2) = -Apq(0,2);
	//	//Apq(1,2) = -Apq(1,2);
	//	//Apq(2,2) = -Apq(2,2);

	//	//�Q�@a2��a3������
	//	//double tmp;
	//	//tmp = Apq(0,2);
	//	//Apq(0,2) = Apq(0,1);
	//	//Apq(0,1) = tmp;

	//	//tmp = Apq(1,2);
	//	//Apq(1,2) = Apq(1,1);
	//	//Apq(1,1) = tmp;

	//	//tmp = Apq(2,2);
	//	//Apq(2,2) = Apq(2,1);
	//	//Apq(2,1) = tmp;
	////}
	//

	//rxMatrix3 R, S;
	////PolarDecomposition(Apq, R, S, m_mtrxBeforeU);
	//PolarDecomposition(Apq, R, S);

	////double end1 = qc.End()/*/100*/;

	//if(m_bLinearDeformation)
	//{
	//	// Linear Deformations
	//	rxMatrix3 A;
	//	A = Apq*Aqq.Inverse();	// A = Apq*Aqq^-1

	//	// �̐ϕۑ��̂��߂Ɂ�(det(A))�Ŋ���
	//	if(m_bVolumeConservation){
	//		double det = fabs(A.Determinant());
	//		if(det > RX_FEQ_EPS){
	//			det = 1.0/sqrt(det);
	//			if(det > 2.0) det = 2.0;
	//			A *= det;
	//		}
	//	}

	//	//cout << "�v���J�n2" << endl;

	//	// �ڕW���W���v�Z���C���݂̒��_���W���ړ�
	//	for(int i = 0; i < m_iNumVertices; ++i){
	//		if(m_pFix[i]) continue;

	//		int cIndx = i*SM_DIM;

	//		// ��]�s��R�̑���̍s��RL=��A+(1-��)R���v�Z
	//		for(int j = 0; j < SM_DIM; j++)
	//		{
	//			q[j] = m_pOrgPos[cIndx+j]-cm_org[j];
	//		}

	//		Vec3 gp(R*q+cm);

	//		for(int j = 0; j < SM_DIM; j++)
	//		{
	//			int jcIndx = cIndx+j;

	//			m_pCurPos[jcIndx] += (gp[j]-m_pCurPos[jcIndx])*m_dAlphas[i];
	//		}
	//	}
	//}
}

/*!
 * ���x�ƈʒu�̍X�V
 *  - �V�����ʒu�ƌ��݂̈ʒu���W���瑬�x���Z�o
 * @param[in] dt �^�C���X�e�b�v��
 */
__device__
	void Integrate(float* curPos, float* vel, int* pIndxes, float dt, int prtNum)
{
	double dt1 = 1.0/dt;

	for(int i = 0; i < prtNum; ++i)
	{
		int pIndx = pIndxes[i]*4;

		for(int j = 0; j < SM_DIM; j++)
		{
			int cIndx = i*SM_DIM+j;

			//vel[cIndx] = (curPos[cIndx] - d_Pos[pIndx+j]) * dt1;/*+ m_v3Gravity * dt * 1.0*/;
		}
	}
	
}

void clamp(float* pos, int cIndx)
{
	//if(pos[cIndx+0] < m_v3Min[0]) pos[cIndx+0] = m_v3Min[0];
	//if(pos[cIndx+0] > m_v3Max[0]) pos[cIndx+0] = m_v3Max[0];
	//if(pos[cIndx+1] < m_v3Min[1]) pos[cIndx+1] = m_v3Min[1];
	//if(pos[cIndx+1] > m_v3Max[1]) pos[cIndx+1] = m_v3Max[1];
	//if(pos[cIndx+2] < m_v3Min[2]) pos[cIndx+2] = m_v3Min[2];
	//if(pos[cIndx+2] > m_v3Max[2]) pos[cIndx+2] = m_v3Max[2];
}