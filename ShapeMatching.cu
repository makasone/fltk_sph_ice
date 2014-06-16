#ifndef _GPU_SHAPE_MATCHING_H_
#define _GPU_SHAPE_MATCHING_H_

#include <math.h>	//�g����H
#include <stdio.h>

//#include <cstdio>
//#include <cmath>
//#include <cstdlib>
//
//#include "rx_utility.h"
//#include <rx_matrix.h>
//
//#include "rx_nnsearch.h"

//#include <rx_cu_common.cu>
#include <rx_cu_common.cuh>	//�搶����`�����֗��@�\���g����@�Ȃ����A�C���N���[�h����ƐԂ��Ȃ�Ȃ�

#define SM_DIM 3

//�����n���ϐ������Ԉ���Ă��Ă��G���[���o�Ȃ��̂Œ���

void LaunchShapeMathcingGPU(float* prtPos, float* prtVel, float* orgPos, float* curPos, float* vel, int* pIndxes, float dt, int prtNum);
__global__ void Update(float* prtPos, float* prtVel, float* orgPos, float* curPos, float* vel, int* pIndxes, float dt, int prtNum);
__device__ void ExternalForce(float* prtPos, float* prtVel, float* curPos, float* vel, int* pIndxes, float dt, int prtNum);
__device__ void ProjectPos(float* prtPos, float* prtVel, float* orgPos, float* curPos, float* vel, int* pIndxes, float dt, int prtNum);
__device__ void Integrate(float* prtPos, float* prtVel, float* curPos, float* vel, int* pIndxes, float dt, int prtNum);

__device__ void PolarDecomposition(matrix3x3 &A, matrix3x3 &R, matrix3x3 &S);

//�s�񉉎Z
__device__ void MakeIdentity(matrix3x3 &M);
__device__ matrix3x3 Transpose(const matrix3x3 &M);
__device__ matrix3x3 Inverse(const matrix3x3 &M);

__device__ matrix3x3 Multiple(matrix3x3 &M1, matrix3x3 &M2);
__device__ float3 Multiple(matrix3x3 &M1, float3& V);

//GPU����
void LaunchShapeMatchingGPU(float* prtPos, float* prtVel, float* orgPos, float* curPos, float* vel, int* pIndxes, float dt, int prtNum)
{
	
	//printf("LaunchGPUKernel");

	dim3 grid(1, 1);
	dim3 block(1, 1, 1);

	//�^���v�Z
	Update <<< grid , block >>> (prtPos, prtVel, orgPos, curPos, vel, pIndxes, dt, prtNum);
}


//GPU�̈ʒu�E���x�X�V
__global__
void Update(
	float* prtPos,
	float* prtVel, 
	float* orgPos,
	float* curPos, 
	float* vel,
	int* pIndxes,
	float dt,
	int prtNum)
{
	//printf("d_Integrate\n");	//�߂��Ⴍ����o��̂Œ���

	ExternalForce(prtPos, prtVel, curPos, vel, pIndxes, dt, prtNum);
	ProjectPos(prtPos, prtVel, orgPos, curPos, vel, pIndxes, dt, prtNum);
	Integrate(prtPos, prtVel, curPos, vel, pIndxes, dt, prtNum);
}

__device__
	void ExternalForce(float* prtPos, float* prtVel, float* curPos, float* vel, int* pIndxes, float dt, int prtNum)
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

			curPos[jcIndx] = prtPos[jpIndx]+prtVel[jpIndx]*dt;
		}
	}

	// ���E�ǂ̉e��
	//���������Ȃ�d���Ȃ邪�C����͂���݂���
	//float res = 0.9;	// �����W��
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
	void ProjectPos(float* prtPos, float* prtVel, float* orgPos, float* curPos, float* vel, int* pIndxes, float dt, int prtNum)
{
	if(prtNum <= 1) return;

	float3 cm = make_float3(0.0, 0.0, 0.0);
	float3 cm_org = make_float3(0.0, 0.0, 0.0);	// �d�S

	float mass = 0.0;	// ������

	// �d�S���W�̌v�Z
	for(int i = 0; i < prtNum;++i){
		//float m = m_pMass[i];
		float m = 1.0;
		//if(m_pFix[i]) m *= 300.0;	// �Œ�_�̎��ʂ�傫������
		
		int cIndx = i*SM_DIM;
		mass += m;

		cm.x += curPos[cIndx+0]*m;
		cm.y += curPos[cIndx+1]*m;
		cm.z += curPos[cIndx+2]*m;

		cm_org.x += orgPos[cIndx+0]*m;
		cm_org.y += orgPos[cIndx+1]*m;
		cm_org.z += orgPos[cIndx+2]*m;
	}

	cm.x /= mass;
	cm.y /= mass;
	cm.z /= mass;
	
	cm_org.x /= mass;
	cm_org.y /= mass;
	cm_org.z /= mass;

	matrix3x3 Apq, Aqq;
	float3 p, q;

	// Apq = ��mpq^T
	// Aqq = ��mqq^T
	for(int i = 0; i < prtNum; ++i)
	{
		int cIndx = i*SM_DIM;

		p.x = curPos[cIndx+0]-cm.x;
		p.y = curPos[cIndx+1]-cm.y;
		p.z = curPos[cIndx+2]-cm.z;

		q.x = orgPos[cIndx+0]-cm_org.x;
		q.y = orgPos[cIndx+1]-cm_org.y;
		q.z = orgPos[cIndx+2]-cm_org.z;

		float m = 1.0f;

		Apq.e[0].x += m*p.x*q.x;
		Apq.e[0].y += m*p.x*q.y;
		Apq.e[0].z += m*p.x*q.z;
		Apq.e[1].x += m*p.y*q.x;
		Apq.e[1].y += m*p.y*q.y;
		Apq.e[1].z += m*p.y*q.z;
		Apq.e[2].x += m*p.z*q.x;
		Apq.e[2].y += m*p.z*q.y;
		Apq.e[2].z += m*p.z*q.z;

		Aqq.e[0].x += m*q.x*q.x;
		Aqq.e[0].y += m*q.x*q.y;
		Aqq.e[0].z += m*q.x*q.z;
		Aqq.e[1].x += m*q.y*q.x;
		Aqq.e[1].y += m*q.y*q.y;
		Aqq.e[1].z += m*q.y*q.z;
		Aqq.e[2].x += m*q.z*q.x;
		Aqq.e[2].y += m*q.z*q.y;
		Aqq.e[2].z += m*q.z*q.z;
	}

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
	//	//float tmp;
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

	matrix3x3 R, S;
	////PolarDecomposition(Apq, R, S, m_mtrxBeforeU);
	PolarDecomposition(Apq, R, S);

	//if(m_bLinearDeformation)
	{
		// Linear Deformations
		matrix3x3 A;
		A = Multiple(Apq, Inverse(Aqq));	// A = Apq*Aqq^-1

		//// �̐ϕۑ��̂��߂Ɂ�(det(A))�Ŋ���
		//if(m_bVolumeConservation){
		//	float det = fabs(A.Determinant());
		//	if(det > RX_FEQ_EPS){
		//		det = 1.0/sqrt(det);
		//		if(det > 2.0) det = 2.0;
		//		A *= det;
		//	}
		//}

		// �ڕW���W���v�Z���C���݂̒��_���W���ړ�
		for(int i = 0; i < prtNum; ++i){
			//if(m_pFix[i]) continue;

			int cIndx = i*SM_DIM;

			// ��]�s��R�̑���̍s��RL=��A+(1-��)R���v�Z
			q.x = orgPos[cIndx+0]-cm_org.x;
			q.y = orgPos[cIndx+1]-cm_org.y;
			q.z = orgPos[cIndx+2]-cm_org.z;

			//Vec3 gp(R*q+cm);
			float3 Rq = Multiple(R, q);
			float3 gp;
			gp.x = Rq.x + cm.x;
			gp.y = Rq.y + cm.y;
			gp.z = Rq.z + cm.z;
			
			curPos[cIndx+0] += (gp.x-curPos[cIndx+0])*1.0f/*m_dAlphas[i]*/;
			curPos[cIndx+1] += (gp.y-curPos[cIndx+1])*1.0f/*m_dAlphas[i]*/;
			curPos[cIndx+2] += (gp.z-curPos[cIndx+2])*1.0f/*m_dAlphas[i]*/;
		}
	}
}

/*!
 * ���x�ƈʒu�̍X�V
 *  - �V�����ʒu�ƌ��݂̈ʒu���W���瑬�x���Z�o
 * @param[in] dt �^�C���X�e�b�v��
 */
__device__
	void Integrate(float* prtPos, float* prtVel, float* curPos, float* vel, int* pIndxes, float dt, int prtNum)
{
	float dt1 = 1.0/dt;

	for(int i = 0; i < prtNum; ++i)
	{
		int pIndx = pIndxes[i]*4;

		for(int j = 0; j < SM_DIM; j++)
		{
			int cIndx = i*SM_DIM+j;

			vel[cIndx] = (curPos[cIndx] - prtPos[pIndx+j]) * dt1;/*+ m_v3Gravity * dt * 1.0*/;
		}
	}
	
}

__device__  
	void clamp(float* pos, int cIndx)
{
	//if(pos[cIndx+0] < m_v3Min[0]) pos[cIndx+0] = m_v3Min[0];
	//if(pos[cIndx+0] > m_v3Max[0]) pos[cIndx+0] = m_v3Max[0];
	//if(pos[cIndx+1] < m_v3Min[1]) pos[cIndx+1] = m_v3Min[1];
	//if(pos[cIndx+1] > m_v3Max[1]) pos[cIndx+1] = m_v3Max[1];
	//if(pos[cIndx+2] < m_v3Min[2]) pos[cIndx+2] = m_v3Min[2];
	//if(pos[cIndx+2] > m_v3Max[2]) pos[cIndx+2] = m_v3Max[2];
}

/*!
 * Jacobi�@�ɂ��ŗL�l�̎Z�o
 * @param[inout] a ���Ώ̍s��D�v�Z��C�Ίp�v�f�ɌŗL�l������
 * @param[out] v �ŗL�x�N�g��(a�Ɠ����T�C�Y)
 * @param[in] n �s��̃T�C�Y(n�~n)
 * @param[in] eps �����덷
 * @param[in] iter_max �ő唽����
 * @return ������
 */
__device__ 
 int EigenJacobiMethod(float *a, float *v, int n, float eps = 1e-8, int iter_max = 100)
{
	float *bim, *bjm;
	float bii, bij, bjj, bji;
 
	bim = new float[n];
	bjm = new float[n];
 
	for(int i = 0; i < n; ++i){
		for(int j = 0; j < n; ++j){
			v[i*n+j] = (i == j) ? 1.0 : 0.0;
		}
	}
 
	int cnt = 0;
	for(;;){
		int i = -1, j = -1;
 
		float x = 0.0;
		for(int ia = 0; ia < n; ++ia){
			for(int ja = 0; ja < n; ++ja){
				int idx = ia*n+ja;
				if(ia != ja && fabs(a[idx]) > x){
					i = ia;
					j = ja;
					x = fabs(a[idx]);
				}
			}
		}

		if(i == -1 || j == -1) return 0;
 
		float aii = a[i*n+i];
		float ajj = a[j*n+j];
		float aij = a[i*n+j];
 
		float m_dAlpha, m_dBeta;
		m_dAlpha = (aii-ajj)/2.0;
		m_dBeta  = sqrt(m_dAlpha*m_dAlpha+aij*aij);
 
		float st, ct;
		ct = sqrt((1.0+fabs(m_dAlpha)/m_dBeta)/2.0);    // sin��
		st = (((aii-ajj) >= 0.0) ? 1.0 : -1.0)*aij/(2.0*m_dBeta*ct);    // cos��
 
		// A = PAP�̌v�Z
		for(int m = 0; m < n; ++m){
			if(m == i || m == j) continue;
 
			float aim = a[i*n+m];
			float ajm = a[j*n+m];
 
			bim[m] =  aim*ct+ajm*st;
			bjm[m] = -aim*st+ajm*ct;
		}
 
		bii = aii*ct*ct+2.0*aij*ct*st+ajj*st*st;
		bij = 0.0;
 
		bjj = aii*st*st-2.0*aij*ct*st+ajj*ct*ct;
		bji = 0.0;
 
		for(int m = 0; m < n; ++m){
			a[i*n+m] = a[m*n+i] = bim[m];
			a[j*n+m] = a[m*n+j] = bjm[m];
		}
		a[i*n+i] = bii;
		a[i*n+j] = bij;
		a[j*n+j] = bjj;
		a[j*n+i] = bji;
 
		// V = PV�̌v�Z
		for(int m = 0; m < n; ++m){
			float vmi = v[m*n+i];
			float vmj = v[m*n+j];
 
			bim[m] =  vmi*ct+vmj*st;
			bjm[m] = -vmi*st+vmj*ct;
		}
		for(int m = 0; m < n; ++m){
			v[m*n+i] = bim[m];
			v[m*n+j] = bjm[m];
		}
 
		float e = 0.0;
		for(int ja = 0; ja < n; ++ja){
			for(int ia = 0; ia < n; ++ia){
				if(ia != ja){
					e += fabs(a[ja*n+ia]);
				}
			}
		}
		if(e < eps) break;
 
		cnt++;
		if(cnt > iter_max) break;
	}
 
	delete [] bim;
	delete [] bjm;
 
	return cnt;
}


/*!
 * �ɕ����ŉ�]�s��ƑΏ̍s��ɕ��� A=RS
 * @param[in] A ���͍s��
 * @param[out] R ��]�s��(�����s�� R^-1 = R^T)
 * @param[out] S �Ώ̍s��
 */
__device__ 
	void PolarDecomposition(matrix3x3 &A, matrix3x3 &R, matrix3x3 &S)
{
	// S = (A^T A)^(1/2)�����߂�
	matrix3x3 ATA;
	ATA = Multiple(Transpose(A), A);	// (A^T A)�̌v�Z

	MakeIdentity(R);

	// (A^T A)���ŗL�l�������đΊp�s��ƒ����s������߂�
	// M^(1/2) = U^T M' U 
	//  M = (A^T A), M':�Ίp�s��̕����������������, U:�����s��

	//�s���float�ɕϊ�
	float* pATA = new float[9];
	float* pU = new float[9];

	for(int i = 0; i < 3; i++)
	{
		int indx = i*3;
		pATA[indx+0] = ATA.e[i].x;
		pATA[indx+1] = ATA.e[i].y;
		pATA[indx+2] = ATA.e[i].z;
	}

	EigenJacobiMethod(pATA, pU, 3);

	//float*���s��ɕϊ�	
	matrix3x3 U;

	for(int i = 0; i < 3; i++)
	{
		int indx = i*3;
		ATA.e[i].x = pATA[indx+0];
		ATA.e[i].y = pATA[indx+1];
		ATA.e[i].z = pATA[indx+2];

		U.e[i].x = pU[indx+0];
		U.e[i].y = pU[indx+1];
		U.e[i].z = pU[indx+2];
	}

	// �Ίp�s��̕��������Ƃ��āC�t�s��v�Z�̂��߂ɋt���ɂ��Ă���
	float l0 = (ATA.e[0].x <= 0.0) ? 0.0 : 1.0/sqrt(ATA.e[0].x);
	float l1 = (ATA.e[1].y <= 0.0) ? 0.0 : 1.0/sqrt(ATA.e[1].y);
	float l2 = (ATA.e[2].z <= 0.0) ? 0.0 : 1.0/sqrt(ATA.e[2].z);

	//// U^T M' U �̋t�s��v�Z
	matrix3x3 S1;
	S1.e[0].x = l0*U.e[0].x*U.e[0].x + l1*U.e[0].y*U.e[0].y + l2*U.e[0].z*U.e[0].z;
	S1.e[0].y = l0*U.e[0].x*U.e[1].x + l1*U.e[0].y*U.e[1].y + l2*U.e[0].z*U.e[1].z;
	S1.e[0].z = l0*U.e[0].x*U.e[2].x + l1*U.e[0].y*U.e[2].y + l2*U.e[0].z*U.e[2].z;
	S1.e[1].x = S1.e[0].x;
	S1.e[1].y = l0*U.e[1].x*U.e[1].x + l1*U.e[1].y*U.e[1].y + l2*U.e[1].z*U.e[1].z;
	S1.e[1].z = l0*U.e[1].x*U.e[2].x + l1*U.e[1].y*U.e[2].y + l2*U.e[1].y*U.e[2].z;
	S1.e[2].x = S1.e[0].z;
	S1.e[2].y = S1.e[1].z;
	S1.e[2].z = l0*U.e[2].x*U.e[2].x + l1*U.e[2].y*U.e[2].y + l2*U.e[2].z*U.e[2].z;

	R = Multiple(A, S1);	// R = A S^-1
	S = Multiple(Transpose(R), A); // S = R^-1 A = R^T A
}

//�z�񏉊���
__device__
	void MakeIdentity(matrix3x3 &M)
{
	M.e[0].x = 1.0;
	M.e[0].y = 0.0;
	M.e[0].z = 0.0;

	M.e[1].x = 0.0;
	M.e[1].y = 1.0;
	M.e[1].z = 0.0;

	M.e[2].x = 0.0;
	M.e[2].y = 0.0;
	M.e[2].z = 1.0;
}

//�]�u�s��
__device__
	matrix3x3 Transpose(const matrix3x3 &M)
{
	matrix3x3 T;

	T.e[0].x = M.e[0].x;
	T.e[0].y = M.e[1].x;
	T.e[0].z = M.e[2].x;

	T.e[1].x = M.e[0].y;
	T.e[1].y = M.e[1].y;
	T.e[1].z = M.e[2].y;

	T.e[2].x = M.e[0].z;
	T.e[2].y = M.e[1].z;
	T.e[2].z = M.e[2].z;

	return T;
}

//�t�s��
__device__
	matrix3x3 Inverse(const matrix3x3 &M)
{
	matrix3x3 I;

	I.e[0].x = M.e[0].x;
	I.e[0].y = M.e[1].x;
	I.e[0].z = M.e[2].x;

	I.e[1].x = M.e[0].y;
	I.e[1].y = M.e[1].y;
	I.e[1].z = M.e[2].y;

	I.e[2].x = M.e[0].z;
	I.e[2].y = M.e[1].z;
	I.e[2].z = M.e[2].z;

	float d = M.e[0].x*M.e[1].y*M.e[2].z- 
			 M.e[0].x*M.e[2].y*M.e[1].z+ 
			 M.e[1].x*M.e[2].y*M.e[0].z- 
			 M.e[1].x*M.e[0].y*M.e[2].z+ 
			 M.e[2].x*M.e[0].y*M.e[1].z- 
			 M.e[2].x*M.e[1].y*M.e[0].z;

	if(d == 0) d = 1;

	I.e[0].x =  (M.e[1].y*M.e[2].z-M.e[1].z*M.e[2].y)/d;
	I.e[0].y = -(M.e[0].y*M.e[2].z-M.e[0].z*M.e[2].y)/d;
	I.e[0].z =  (M.e[0].y*M.e[1].z-M.e[0].z*M.e[1].y)/d;
	I.e[1].x = -(M.e[1].x*M.e[2].z-M.e[1].z*M.e[2].x)/d;
	I.e[1].y =  (M.e[0].x*M.e[2].z-M.e[0].z*M.e[2].x)/d;
	I.e[1].z = -(M.e[0].x*M.e[1].z-M.e[0].z*M.e[1].x)/d;
	I.e[2].x =  (M.e[1].x*M.e[2].y-M.e[1].y*M.e[2].x)/d;
	I.e[2].y = -(M.e[0].x*M.e[2].y-M.e[0].y*M.e[2].x)/d;
	I.e[2].z =  (M.e[0].x*M.e[1].y-M.e[0].y*M.e[1].x)/d;

	return I;
}

__device__ 
	matrix3x3 Multiple(matrix3x3 &M1, matrix3x3 &M2)
{
	matrix3x3 M;

	M.e[0].x = M1.e[0].x * M2.e[0].x + M1.e[0].y * M2.e[1].x + M1.e[0].z * M2.e[2].x;
	M.e[0].y = M1.e[0].x * M2.e[0].y + M1.e[0].y * M2.e[1].y + M1.e[0].z * M2.e[2].y;
	M.e[0].z = M1.e[0].x * M2.e[0].z + M1.e[0].y * M2.e[1].z + M1.e[0].z * M2.e[2].z;
	 
	M.e[1].x = M1.e[1].x * M2.e[0].x + M1.e[1].y * M2.e[1].x + M1.e[1].z * M2.e[2].x;
	M.e[1].y = M1.e[1].x * M2.e[0].y + M1.e[1].y * M2.e[1].y + M1.e[1].z * M2.e[2].y;
	M.e[1].z = M1.e[1].x * M2.e[0].z + M1.e[1].y * M2.e[1].z + M1.e[1].z * M2.e[2].z;
	 
	M.e[2].x = M1.e[2].x * M2.e[0].x + M1.e[2].y * M2.e[1].x + M1.e[2].z * M2.e[2].x;
	M.e[2].y = M1.e[2].x * M2.e[0].y + M1.e[2].y * M2.e[1].y + M1.e[2].z * M2.e[2].y;
	M.e[2].z = M1.e[2].x * M2.e[0].z + M1.e[2].y * M2.e[1].z + M1.e[2].z * M2.e[2].z;

	return M;
}

__device__ 
	float3 Multiple(matrix3x3 &M1, float3 &V)
{
	float3 M;

	M.x = M1.e[0].x * V.x + M1.e[0].y * V.y + M1.e[0].z * V.z;
	M.y = M1.e[1].x * V.y + M1.e[1].y * V.y + M1.e[1].z * V.z;
	M.z = M1.e[2].x * V.z + M1.e[2].y * V.y + M1.e[2].z * V.z;

	return M;
}

#endif