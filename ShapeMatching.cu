//SM�@��GPU���� �^���v�Z
//�Ƃ肠��������Ă݂āC�f�[�^�\���͌�Œ���

#ifndef _GPU_SHAPE_MATCHING_H_
#define _GPU_SHAPE_MATCHING_H_

#include <math.h>
#include <stdio.h>
#include <iostream>
#include <cuda_runtime.h>

using namespace std;

//�ȉ��͂��܂��C���N���[�h�ł��Ȃ�����
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
void LaunchShapeMathcingGPU(int prtNum, float* prtPos, float* prtVel, float* orgPos, float* curPos, float* orgCm, float* curCm, float* vel, int* pIndxes, int* indxSet, float dt);
__global__ void Update(float* prtPos, float* prtVel, float* orgPos, float* curPos, float* orgCm, float* curCm, float* vel, int* pIndxes, int* indxSet, float dt, int prtNum, int side);
__device__ void ExternalForce(float* prtPos, float* prtVel, float* curPos, float* vel, int* pIndxes, int* indxSet, float dt, int prtNum, int side);
__device__ void ProjectPos(float* prtPos, float* prtVel, float* orgPos, float* curPos, float* orgCm, float* vel, int* pIndxes, int* indxSet, float dt, int prtNum, int side);
__device__ void Integrate(float* prtPos, float* prtVel, float* curPos, float* vel, int* pIndxes, int* indxSet, float dt, int prtNum, int side);


void LaunchShapeMathcingIterationGPU(int prtNum, float* prtPos, float* prtVel, float* sldPos, float* sldVel, float* orgPos, float* curPos, float* orgCm, float* curCm, float* vel, int* pIndxes, int* indxSet, float dt);
__global__ void UpdateIteration(float* prtPos, float* prtVel, float* sldPos, float* sldVel, float* orgPos, float* curPos, float* orgCm, float* curCm, float* vel, int* pIndxes, int* indxSet, float dt, int prtNum, int side);
__device__ void ExternalForceIteration(float* sldPos, float* sldVel, float* curPos, float* vel, int* pIndxes, int* indxSet, float dt, int prtNum, int side);
__device__ void ProjectPosIteration(float* prtPos, float* prtVel, float* orgPos, float* curPos, float* vel, int* pIndxes, int* indxSet, float dt, int prtNum, int side);

__device__ void PolarDecomposition(matrix3x3 &A, matrix3x3 &R, matrix3x3 &S);

//�s�񉉎Z
//TODO::�Ă��Ƃ��Ȃ̂ł��������g���₷��
__device__ void Clamp(float* curPos, int pIndx);
__device__ void MakeIdentity(matrix3x3 &M);
__device__ matrix3x3 Transpose(const matrix3x3 &M);
__device__ matrix3x3 Inverse(const matrix3x3 &M);

__device__ matrix3x3 Multiple(const matrix3x3 &M1, const matrix3x3 &M2);
__device__ float3 Multiple(const matrix3x3 &M1, const float3& V);

//�^���v�Z
void LaunchShapeMatchingGPU(int prtNum, float* prtPos, float* prtVel, float* orgPos, float* curPos, float* orgCm, float* curCm, float* vel, int* pIndxes, int* indxSet, float dt)
{
	//TODO::�����̂�O��ɂɃO���b�h�������Ă���̂ŁC�����镨�̂ɑΉ��ł���悤�ɂ���
	int side = pow( prtNum, 1.0/3.0 ) + 0.5;	//�����̂̂P�ӂ̒��_��

	dim3 grid(side, side);
	dim3 block(side, 1, 1);


	//�^���v�Z
	Update<<<grid ,block>>>(prtPos, prtVel, orgPos, curPos, orgCm, curCm, vel, pIndxes, indxSet, dt, prtNum, side);
	
	cudaThreadSynchronize();
}

//�^���v�Z
void LaunchShapeMatchingIterationGPU(int prtNum, float* prtPos, float* prtVel, float* sldPos, float* sldVel, float* orgPos, float* curPos, float* orgCm, float* curCm, float* vel, int* pIndxes, int* indxSet, float dt)
{
	//TODO::�����̂�O��ɂɃO���b�h�������Ă���̂ŁC�����镨�̂ɑΉ��ł���悤�ɂ���
	int side = pow( prtNum, 1.0/3.0 ) + 0.5;	//�����̂̂P�ӂ̒��_��

	dim3 grid(side, side);
	dim3 block(side, 1, 1);

	//�^���v�Z
	UpdateIteration<<<grid ,block>>>(prtPos, prtVel, sldPos, sldVel, orgPos, curPos, orgCm, curCm, vel, pIndxes, indxSet, dt, prtNum, side);

	cudaThreadSynchronize();
}

//GPU�̈ʒu�E���x�X�V
__global__
void Update(
	float* prtPos,
	float* prtVel, 
	float* orgPos,
	float* curPos,
	float* orgCm,
	float* curCm,
	float* vel,
	int* pIndxes,
	int* indxSet,
	float dt,
	int prtNum,
	int side)
{
	ExternalForce(prtPos, prtVel, curPos, vel, pIndxes, indxSet, dt, prtNum, side);
	ProjectPos(prtPos, prtVel, orgPos, curPos, orgCm, vel, pIndxes, indxSet, dt, prtNum, side);
	Integrate(prtPos, prtVel, curPos, vel, pIndxes, indxSet, dt, prtNum, side);
}

//GPU�̈ʒu�E���x�X�V
__global__
void UpdateIteration(
	float* prtPos,
	float* prtVel,
	float* sldPos,
	float* sldVel, 
	float* orgPos,
	float* curPos,
	float* orgCm,
	float* curCm,
	float* vel,
	int* pIndxes,
	int* indxSet,
	float dt,
	int prtNum,
	int side)
{
	ExternalForceIteration(sldPos, sldVel, curPos, vel, pIndxes, indxSet, dt, prtNum, side);
	ProjectPos(prtPos, prtVel, orgPos, curPos, orgCm, vel, pIndxes, indxSet, dt, prtNum, side);
}

__device__
	void ExternalForceIteration(float* sldPos, float* sldVel, float* curPos, float* vel, int* pIndxes, int* indxSet, float dt, int prtNum, int side)
{
	if(blockIdx.x > side) return;
	if(blockIdx.y > side) return;
	if(threadIdx.x > side) return;
	if(threadIdx.y > side) return;

	//�v�Z����N���X�^�̔���
	//int clusterIndx = blockIdx.x * side * side + blockIdx.y * side + threadIdx.x;
	int clusterIndx = threadIdx.x * side * side + threadIdx.y * side + blockIdx.x;

	int startIndx = indxSet[clusterIndx*2+0];
	int endIndx = indxSet[clusterIndx*2+1];
	//printf("startIndx = %d, endIndx = %d\n", startIndx, endIndx);

	// �d�͂̉e����t���C���x�𔽉f
	for(int i = startIndx; i < endIndx+1; ++i)
	{
		int pIndx = pIndxes[i]*3;
		int cIndx = i*SM_DIM;

		for(int j = 0; j < SM_DIM; ++j)
		{
			int jpIndx = pIndx+j;
			int jcIndx = cIndx+j;

			curPos[jcIndx] = sldPos[jpIndx];
		}
	}

	// ���E�ǂ̉e��
	//���������Ȃ�d���Ȃ邪�C����͂���݂���
	float res = 0.9;	// �����W��
	for(int i = startIndx; i < endIndx; ++i){
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

		Clamp(curPos, i*SM_DIM);
	}
}

__device__
	void ExternalForce(float* prtPos, float* prtVel, float* curPos, float* vel, int* pIndxes, int* indxSet, float dt, int prtNum, int side)
{
	if(blockIdx.x > side) return;
	if(blockIdx.y > side) return;
	if(threadIdx.x > side) return;
	if(threadIdx.y > side) return;

	//�v�Z����N���X�^�̔���
	int clusterIndx = blockIdx.x * side * side + blockIdx.y * side + threadIdx.x;
	//int clusterIndx = threadIdx.x * side * side + threadIdx.y * side + blockIdx.x;

	int startIndx	= indxSet[clusterIndx*2+0];
	int endIndx		= indxSet[clusterIndx*2+1];

	// �d�͂̉e����t���C���x�𔽉f
	for(int i = startIndx; i < endIndx+1; ++i)
	{
		int pIndx = pIndxes[i]*4;
		int cIndx = i*SM_DIM;

		curPos[cIndx+0] = prtPos[pIndx+0] + prtVel[pIndx+0]*dt;
		curPos[cIndx+1] = prtPos[pIndx+1] + prtVel[pIndx+1]*dt;
		curPos[cIndx+2] = prtPos[pIndx+2] + prtVel[pIndx+2]*dt;
	}

	// ���E�ǂ̉e��
	//���������Ȃ�d���Ȃ邪�C����͂���݂���
	float res = 0.9;	// �����W��
	for(int i = startIndx; i < endIndx; ++i){
		//if(m_pFix[i]) continue;
		//Vec3 &p = m_vCurPos[i];
		//Vec3 &np = m_vNewPos[i];
		//Vec3 &v = m_vVel[i];
		//if(np[0] < m_v3Min[0] || np[0] > m_v3Max[0]){
		//	np[0] = p[0]-v[0]*dt*res;
		//	np[1] = p[1];
		//	np[2] = p[2];
		//}
		//if(np[1] < m_v3Min[1] || np[1] > m_v3Max[1]){
		//	np[1] = p[1]-v[1]*dt*res;
		//	np[0] = p[0] ;
		//	np[2] = p[2];
		//}
		//if(np[2] < m_v3Min[2] || np[2] > m_v3Max[2]){
		//	np[2] = p[2]-v[2]*dt*res;
		//	np[0] = p[0];
		//	np[1] = p[1];
		//}

		Clamp(curPos, i*SM_DIM);
	}
}

__device__
	void ProjectPos(float* prtPos, float* prtVel, float* orgPos, float* curPos, float* orgCm, float* vel, int* pIndxes, int* indxSet, float dt, int prtNum, int side)
{
	if(blockIdx.x > side) return;
	if(blockIdx.y > side) return;
	if(threadIdx.x > side) return;
	if(threadIdx.y > side) return;

	//�v�Z����N���X�^�̔���
	int clusterIndx = blockIdx.x * side * side + blockIdx.y * side + threadIdx.x;
	//int clusterIndx = threadIdx.x * side * side + threadIdx.y * side + blockIdx.x;

	if(prtNum <= 1) return;
	if(clusterIndx > prtNum) return;

	int startIndx	 = indxSet[clusterIndx*2+0];
	int endIndx		 = indxSet[clusterIndx*2+1];
	
	float cm_x, cm_y, cm_z;
	float cm_org_x, cm_org_y, cm_org_z;

	float mass = 0.0f;	// ������

	// �d�S���W�̌v�Z
	for(int i = startIndx; i < endIndx+1; ++i)
	{
		//float m = m_pMass[i];
		float m = 1.0f;
		//if(m_pFix[i]) m *= 300.0;	// �Œ�_�̎��ʂ�傫������
		
		int cIndx = i*SM_DIM;
		mass += m;

		cm_x += curPos[cIndx+0]*m;
		cm_y += curPos[cIndx+1]*m;
		cm_z += curPos[cIndx+2]*m;
	}

	cm_x /= mass;
	cm_y /= mass;
	cm_z /= mass;

	//�O�v�Z���������f�[�^�𗘗p
	cm_org_x = orgCm[clusterIndx*SM_DIM+0];
	cm_org_y = orgCm[clusterIndx*SM_DIM+1];
	cm_org_z = orgCm[clusterIndx*SM_DIM+2];

	matrix3x3 Apq;
	matrix3x3 Aqq;
	
	float p_x, p_y, p_z;
	float q_x, q_y, q_z;

	// Apq = ��mpq^T
	// Aqq = ��mqq^T
	for(int i = startIndx; i < endIndx+1; ++i)
	{
		int cIndx = i*SM_DIM;

		p_x = curPos[cIndx+0]-cm_x;
		p_y = curPos[cIndx+1]-cm_y;
		p_z = curPos[cIndx+2]-cm_z;

		q_x = orgPos[cIndx+0]-cm_org_x;
		q_y = orgPos[cIndx+1]-cm_org_y;
		q_z = orgPos[cIndx+2]-cm_org_z;

		float m = 1.0f;

		Apq.e[0].x += m*p_x*q_x;
		Apq.e[0].y += m*p_x*q_y;
		Apq.e[0].z += m*p_x*q_z;
		Apq.e[1].x += m*p_y*q_x;
		Apq.e[1].y += m*p_y*q_y;
		Apq.e[1].z += m*p_y*q_z;
		Apq.e[2].x += m*p_z*q_x;
		Apq.e[2].y += m*p_z*q_y;
		Apq.e[2].z += m*p_z*q_z;

		Aqq.e[0].x += m*q_x*q_x;
		Aqq.e[0].y += m*q_x*q_y;
		Aqq.e[0].z += m*q_x*q_z;
		Aqq.e[1].x += m*q_y*q_x;
		Aqq.e[1].y += m*q_y*q_y;
		Aqq.e[1].z += m*q_y*q_z;
		Aqq.e[2].x += m*q_z*q_x;
		Aqq.e[2].y += m*q_z*q_y;
		Aqq.e[2].z += m*q_z*q_z;
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
		//matrix3x3 A;
		//A = Multiple(Apq, Inverse(Aqq));	// A = Apq*Aqq^-1

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
		for(int i = startIndx; i < endIndx+1; ++i)
		{
			//if(m_pFix[i]) continue;

			int cIndx = i*SM_DIM;

			// ��]�s��R�̑���̍s��RL=��A+(1-��)R���v�Z
			q_x = orgPos[cIndx+0] - cm_org_x;
			q_y = orgPos[cIndx+1] - cm_org_y;
			q_z = orgPos[cIndx+2] - cm_org_z;

			float Rq_x, Rq_y, Rq_z;
			Rq_x = R.e[0].x * q_x + R.e[0].y * q_y + R.e[0].z * q_z;
			Rq_y = R.e[1].x * q_y + R.e[1].y * q_y + R.e[1].z * q_z;
			Rq_z = R.e[2].x * q_z + R.e[2].y * q_y + R.e[2].z * q_z;
			
			curPos[cIndx+0] += (Rq_x + cm_x - curPos[cIndx+0]) * 1.0f/*m_dAlphas[i]*/;
			curPos[cIndx+1] += (Rq_y + cm_y - curPos[cIndx+1]) * 1.0f/*m_dAlphas[i]*/;
			curPos[cIndx+2] += (Rq_z + cm_z - curPos[cIndx+2]) * 1.0f/*m_dAlphas[i]*/;
		}
	}
}

/*!
 * ���x�ƈʒu�̍X�V
 *  - �V�����ʒu�ƌ��݂̈ʒu���W���瑬�x���Z�o
 * @param[in] dt �^�C���X�e�b�v��
 */
__device__
	void Integrate(float* prtPos, float* prtVel, float* curPos, float* vel, int* pIndxes, int* indxSet, float dt, int prtNum, int side)
{
	if(blockIdx.x > side) return;
	if(blockIdx.y > side) return;
	if(threadIdx.x > side) return;
	if(threadIdx.y > side) return;

	//�v�Z����N���X�^�̔���
	int clusterIndx = blockIdx.x * side * side + blockIdx.y * side + threadIdx.x;
	//int clusterIndx = threadIdx.x * side * side + threadIdx.y * side + blockIdx.x;

	int startIndx = indxSet[clusterIndx*2+0];
	int endIndx = indxSet[clusterIndx*2+1];

	float dt1 = 1.0f/dt;
	float gravity[3] = {0.0f, -9.81f, 0.0f};
	float param = 1.0f;
	gravity[0] *= dt * param;
	gravity[1] *= dt * param;
	gravity[2] *= dt * param;

	for(int i = startIndx; i < endIndx+1; ++i)
	{
		int pIndx = pIndxes[i]*4;
		int cIndx = i*SM_DIM;

		vel[cIndx+0] = (curPos[cIndx+0] - prtPos[pIndx+0]) * dt1 + gravity[0];
		vel[cIndx+1] = (curPos[cIndx+1] - prtPos[pIndx+1]) * dt1 + gravity[1];
		vel[cIndx+2] = (curPos[cIndx+2] - prtPos[pIndx+2]) * dt1 + gravity[2];
	}
}

__device__  
	void Clamp(float* pos, int cIndx)
{
	//if(pos[cIndx+0] < m_v3Min[0]) pos[cIndx+0] = m_v3Min[0];
	//if(pos[cIndx+0] > m_v3Max[0]) pos[cIndx+0] = m_v3Max[0];
	//if(pos[cIndx+1] < m_v3Min[1]) pos[cIndx+1] = m_v3Min[1];
	//if(pos[cIndx+1] > m_v3Max[1]) pos[cIndx+1] = m_v3Max[1];
	//if(pos[cIndx+2] < m_v3Min[2]) pos[cIndx+2] = m_v3Min[2];
	//if(pos[cIndx+2] > m_v3Max[2]) pos[cIndx+2] = m_v3Max[2];
	if(pos[cIndx+0] > 0.75f) pos[cIndx+0] = 0.75f;
	if(pos[cIndx+0] < -0.75f) pos[cIndx+0] = -0.75f;
	if(pos[cIndx+1] > 1.5f) pos[cIndx+1] = 1.5f;
	if(pos[cIndx+1] < -1.5f) pos[cIndx+1] = -1.5f;
	if(pos[cIndx+2] > 2.0f) pos[cIndx+2] = 2.0f;
	if(pos[cIndx+2] < -2.0f) pos[cIndx+2] = -2.0f;
}






/*!
 * Jacobi�@�ɂ��ŗL�l�̎Z�o+
 * @param[inout] a ���Ώ̍s��D�v�Z��C�Ίp�v�f�ɌŗL�l������
 * @param[out] v �ŗL�x�N�g��(a�Ɠ����T�C�Y)
 * @param[in] n �s��̃T�C�Y(n�~n)
 * @param[in] eps �����덷
 * @param[in] iter_max �ő唽����
 * @return ������
 */
__device__ 
// int d_EigenJacobiMethod(float *a, float *v, int n, float eps = 1e-8, int iter_max = 100)
 int d_EigenJacobiMethod(matrix3x3& a, matrix3x3& v, int n, float eps = 1e-8, int iter_max = 100)
{
	float bii, bij, bjj, bji;

	float bim[3] = {0.0f, 0.0f, 0.0f};
	float bjm[3] = {0.0f, 0.0f, 0.0f};

	MakeIdentity(v);

	int cnt = 0;

	for(;;){
		int i = -1, j = -1;
 
		float x = 0.0f;
		for(int ia = 0; ia < n; ++ia){
				//n=3�Ƃ��Ă���̂ŁC�������
				if(ia != 0 && fabs(a.e[ia].x) > x){
					i = ia;
					j = 0;
					x = fabs(a.e[ia].x);
				}

				if(ia != 1 && fabs(a.e[ia].y) > x){
					i = ia;
					j = 1;
					x = fabs(a.e[ia].y);
				}

				if(ia != 2 && fabs(a.e[ia].z) > x){
					i = ia;
					j = 2;
					x = fabs(a.e[ia].z);
				}
		}

		if(i == -1 || j == -1) return 0;
 
		float aii = 0.0f;
		float ajj = 0.0f;
		float aij = 0.0f;

		//(i,i)
		if(i == 0){	aii = a.e[0].x;	}
		if(i == 1){	aii = a.e[1].y;	}
		if(i == 2){	aii = a.e[2].z;	}

		//(j,j), (i,j)
		if(j == 0){	ajj = a.e[0].x;		aij = a.e[i].x;	}
		if(j == 1){	ajj = a.e[1].y;		aij = a.e[i].y;	}
		if(j == 2){	ajj = a.e[2].z;		aij = a.e[i].z;	}

		float m_dAlpha, m_dBeta;
		m_dAlpha = (aii-ajj)/2.0f;
		m_dBeta  = sqrt(m_dAlpha*m_dAlpha+aij*aij);
 
		float st, ct;
		ct = sqrt((1.0+fabs(m_dAlpha)/m_dBeta)/2.0f);    // sin��
		st = (((aii-ajj) >= 0.0f) ? 1.0f : -1.0f)*aij/(2.0f*m_dBeta*ct);    // cos��
 
		// A = PAP�̌v�Z
		for(int m = 0; m < n; ++m){
			if(m == i || m == j) continue;
			
			float aim = 0.0f;
			float ajm = 0.0f;

			if(m == 0){	aim = a.e[i].x; ajm = a.e[j].x;	}
			if(m == 1){	aim = a.e[i].y; ajm = a.e[j].y;	}
			if(m == 2){	aim = a.e[i].z; ajm = a.e[j].z;	}
 
			bim[m] =  aim*ct+ajm*st;
			bjm[m] = -aim*st+ajm*ct;
		}
 
		bii = aii*ct*ct+2.0*aij*ct*st+ajj*st*st;
		bij = 0.0f;
 
		bjj = aii*st*st-2.0*aij*ct*st+ajj*ct*ct;
		bji = 0.0f;
 
		for(int m = 0; m < n; ++m){

			if(m == 0){	a.e[i].x = bim[m];	}
			if(m == 1){	a.e[i].y = bim[m];	}
			if(m == 2){	a.e[i].z = bim[m];	}

			if(i == 0){	a.e[m].x = bim[m];	}
			if(i == 1){	a.e[m].y = bim[m];	}
			if(i == 2){	a.e[m].z = bim[m];	}

			if(m == 0){	a.e[j].x = bjm[m];	}
			if(m == 1){	a.e[j].y = bjm[m];	}
			if(m == 2){	a.e[j].z = bjm[m];	}

			if(j == 0){	a.e[m].x = bjm[m];	}
			if(j == 1){	a.e[m].y = bjm[m];	}
			if(j == 2){	a.e[m].z = bjm[m];	}
		}

		//(i,i), (j, i)
		if(i == 0){	a.e[0].x = bii;		a.e[j].x = bji;	}
		if(i == 1){	a.e[1].y = bii;		a.e[j].y = bji;	}
		if(i == 2){	a.e[2].z = bii;		a.e[j].z = bji;	}

		//(j,j), (i,j)
		if(j == 0){	a.e[0].x = bjj;		a.e[i].x = bij;	}
		if(j == 1){	a.e[1].y = bjj;		a.e[i].y = bij;	}
		if(j == 2){	a.e[2].z = bjj;		a.e[i].z = bij;	}
 
		// V = PV�̌v�Z
		for(int m = 0; m < n; ++m){
			float vmi = 0.0f;
			float vmj = 0.0f;

			if(i == 0){	vmi = v.e[m].x; }
			if(i == 1){	vmi = v.e[m].y; }
			if(i == 2){	vmi = v.e[m].z; }

			if(j == 0){	vmj = v.e[m].x; }
			if(j == 1){	vmj = v.e[m].y; }
			if(j == 2){	vmj = v.e[m].z; }

			bim[m] =  vmi*ct+vmj*st;
			bjm[m] = -vmi*st+vmj*ct;
		}
		for(int m = 0; m < n; ++m){

			if(i == 0){	v.e[m].x = bim[m]; }
			if(i == 1){	v.e[m].y = bim[m]; }
			if(i == 2){	v.e[m].z = bim[m]; }

			if(j == 0){	v.e[m].x = bjm[m]; }
			if(j == 1){	v.e[m].y = bjm[m]; }
			if(j == 2){	v.e[m].z = bjm[m]; }
		}
 
		float e = 0.0f;
		for(int ja = 0; ja < n; ++ja){
			for(int ia = 0; ia < n; ++ia){
				if(ia != ja){
					if(ia == 0){	e += fabs(a.e[ja].x); }
					if(ia == 1){	e += fabs(a.e[ja].y); }
					if(ia == 2){	e += fabs(a.e[ja].z); }
				}
			}
		}
		if(e < eps) break;
 
		cnt++;
		if(cnt > iter_max) break;
	}
 
	return cnt;


	//float bii, bij, bjj, bji;

	//for(int i = 0; i < n; ++i){
	//	for(int j = 0; j < n; ++j){
	//		v[i*n+j] = (i == j) ? 1.0f : 0.0f;
	//	}
	//}
 //
	//int cnt = 0;
	//for(;;){
	//	int i = -1, j = -1;
 //
	//	float x = 0.0;
	//	for(int ia = 0; ia < n; ++ia){
	//		for(int ja = 0; ja < n; ++ja){
	//			int idx = ia*n+ja;
	//			if(ia != ja && fabs(a[idx]) > x){
	//				i = ia;
	//				j = ja;
	//				x = fabs(a[idx]);
	//			}
	//		}
	//	}

	//	if(i == -1 || j == -1) return 0;
 //
	//	float aii = a[i*n+i];
	//	float ajj = a[j*n+j];
	//	float aij = a[i*n+j];
 //
	//	float m_dAlpha, m_dBeta;
	//	m_dAlpha = (aii-ajj)/2.0f;
	//	m_dBeta  = sqrt(m_dAlpha*m_dAlpha+aij*aij);
 //
	//	float st, ct;
	//	ct = sqrt((1.0+fabs(m_dAlpha)/m_dBeta)/2.0f);    // sin��
	//	st = (((aii-ajj) >= 0.0f) ? 1.0f : -1.0f)*aij/(2.0*m_dBeta*ct);    // cos��
 //
	//	// A = PAP�̌v�Z
	//	for(int m = 0; m < n; ++m){
	//		if(m == i || m == j) continue;
 //
	//		float aim = a[i*n+m];
	//		float ajm = a[j*n+m];
 //
	//		bim[m] =  aim*ct+ajm*st;
	//		bjm[m] = -aim*st+ajm*ct;
	//	}
 //
	//	bii = aii*ct*ct+2.0f*aij*ct*st+ajj*st*st;
	//	bij = 0.0f;
 //
	//	bjj = aii*st*st-2.0f*aij*ct*st+ajj*ct*ct;
	//	bji = 0.0f;
 //
	//	for(int m = 0; m < n; ++m){
	//		a[i*n+m] = a[m*n+i] = bim[m];
	//		a[j*n+m] = a[m*n+j] = bjm[m];
	//	}
	//	a[i*n+i] = bii;
	//	a[i*n+j] = bij;
	//	a[j*n+j] = bjj;
	//	a[j*n+i] = bji;
 //
	//	// V = PV�̌v�Z
	//	for(int m = 0; m < n; ++m){
	//		float vmi = v[m*n+i];
	//		float vmj = v[m*n+j];
 //
	//		bim[m] =  vmi*ct+vmj*st;
	//		bjm[m] = -vmi*st+vmj*ct;
	//	}
	//	for(int m = 0; m < n; ++m){
	//		v[m*n+i] = bim[m];
	//		v[m*n+j] = bjm[m];
	//	}
 //
	//	float e = 0.0f;
	//	for(int ja = 0; ja < n; ++ja){
	//		for(int ia = 0; ia < n; ++ia){
	//			if(ia != ja){
	//				e += fabs(a[ja*n+ia]);
	//			}
	//		}
	//	}
	//	if(e < eps) break;
 //
	//	cnt++;
	//	if(cnt > iter_max) break;
	//}
 //
	//delete [] bim;
	//delete [] bjm;
 //
	//return cnt;
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
	matrix3x3 ATA = Multiple(Transpose(A), A);	// (A^T A)�̌v�Z

	//Multiple���g��Ȃ��o�[�W����
	//matrix3x3 TrA = Transpose(A);
	//matrix3x3 ATA;

	//ATA.e[0].x = TrA.e[0].x * A.e[0].x + TrA.e[0].y * A.e[1].x + TrA.e[0].z * A.e[2].x;
	//ATA.e[0].y = TrA.e[0].x * A.e[0].y + TrA.e[0].y * A.e[1].y + TrA.e[0].z * A.e[2].y;
	//ATA.e[0].z = TrA.e[0].x * A.e[0].z + TrA.e[0].y * A.e[1].z + TrA.e[0].z * A.e[2].z;
	// 			 			  			 			  			 		 
	//ATA.e[1].x = TrA.e[1].x * A.e[0].x + TrA.e[1].y * A.e[1].x + TrA.e[1].z * A.e[2].x;
	//ATA.e[1].y = TrA.e[1].x * A.e[0].y + TrA.e[1].y * A.e[1].y + TrA.e[1].z * A.e[2].y;
	//ATA.e[1].z = TrA.e[1].x * A.e[0].z + TrA.e[1].y * A.e[1].z + TrA.e[1].z * A.e[2].z;
	// 			 			  			 			  			 		  
	//ATA.e[2].x = TrA.e[2].x * A.e[0].x + TrA.e[2].y * A.e[1].x + TrA.e[2].z * A.e[2].x;
	//ATA.e[2].y = TrA.e[2].x * A.e[0].y + TrA.e[2].y * A.e[1].y + TrA.e[2].z * A.e[2].y;
	//ATA.e[2].z = TrA.e[2].x * A.e[0].z + TrA.e[2].y * A.e[1].z + TrA.e[2].z * A.e[2].z;

	MakeIdentity(R);

	matrix3x3 U;

	// (A^T A)���ŗL�l�������đΊp�s��ƒ����s������߂�
	//  M^(1/2) = U^T M' U 
	//  M = (A^T A), M':�Ίp�s��̕����������������, U:�����s��

	d_EigenJacobiMethod(ATA, U, 3);

	// �Ίp�s��̕��������Ƃ��āC�t�s��v�Z�̂��߂ɋt���ɂ��Ă���
	float l0 = (ATA.e[0].x <= 0.0) ? 0.0f : 1.0f/sqrt(ATA.e[0].x);
	float l1 = (ATA.e[1].y <= 0.0) ? 0.0f : 1.0f/sqrt(ATA.e[1].y);
	float l2 = (ATA.e[2].z <= 0.0) ? 0.0f : 1.0f/sqrt(ATA.e[2].z);

	//// U^T M' U �̋t�s��v�Z
	matrix3x3 S1;
	S1.e[0].x = l0*U.e[0].x*U.e[0].x + l1*U.e[0].y*U.e[0].y + l2*U.e[0].z*U.e[0].z;
	S1.e[0].y = l0*U.e[0].x*U.e[1].x + l1*U.e[0].y*U.e[1].y + l2*U.e[0].z*U.e[1].z;
	S1.e[0].z = l0*U.e[0].x*U.e[2].x + l1*U.e[0].y*U.e[2].y + l2*U.e[0].z*U.e[2].z;
	S1.e[1].x = S1.e[0].y;
	S1.e[1].y = l0*U.e[1].x*U.e[1].x + l1*U.e[1].y*U.e[1].y + l2*U.e[1].z*U.e[1].z;
	S1.e[1].z = l0*U.e[1].x*U.e[2].x + l1*U.e[1].y*U.e[2].y + l2*U.e[1].z*U.e[2].z;
	S1.e[2].x = S1.e[0].z;
	S1.e[2].y = S1.e[1].z;
	S1.e[2].z = l0*U.e[2].x*U.e[2].x + l1*U.e[2].y*U.e[2].y + l2*U.e[2].z*U.e[2].z;

	R = Multiple(A, S1);			// R = A S^-1
	S = Multiple(Transpose(R), A);	// S = R^-1 A = R^T A
}

//�z�񏉊���
__device__
	void MakeIdentity(matrix3x3 &M)
{
	M.e[0].x = 1.0f;
	M.e[0].y = 0.0f;
	M.e[0].z = 0.0f;

	M.e[1].x = 0.0f;
	M.e[1].y = 1.0f;
	M.e[1].z = 0.0f;

	M.e[2].x = 0.0f;
	M.e[2].y = 0.0f;
	M.e[2].z = 1.0f;
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
	matrix3x3 Multiple(const matrix3x3 &M1, const matrix3x3 &M2)
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
	float3 Multiple(const matrix3x3 &M1, const float3 &V)
{
	float3 M;

	M.x = M1.e[0].x * V.x + M1.e[0].y * V.y + M1.e[0].z * V.z;
	M.y = M1.e[1].x * V.y + M1.e[1].y * V.y + M1.e[1].z * V.z;
	M.z = M1.e[2].x * V.z + M1.e[2].y * V.y + M1.e[2].z * V.z;

	return M;
}

#endif