//�e�N���X�^�̏�񂩂�ŏI�I�Ȍő̂̏����Z�o
//���͒P���ɕ���

#ifndef _GPU_ICE_STRUCTURE_H_
#define _GPU_ICE_STRUCTURE_H_

#include <iostream>
#include <cuda_runtime.h>
#include <rx_cu_common.cuh>

using namespace std;

void LaunchCalcAverageGPU
	(
	int prtNum,
	float* sldPrtPos, 
	float* sldPrtVel, 
	float* sphPrtPos, 
	float* sphPrtVel, 
	float* smPrtPos, 
	float* smPrtVel,
	int* smIndxSet,
	int* PtoCIndx,
	int* PtoC,
	int PNumMax,
	int PtoCMax,
	int PtoCParamSize
	);

__global__ void CalcAverage
	(
	float* sldPrtPos, 
	float* sldPrtVel, 
	float* sphPrtPos,
	float* sphPrtVel, 
	float* smPrtPos, 
	float* smPrtVel, 
	int* indxSet, 
	int* PtoCIndx, 
	int* PtoC, 
	int PNumMax, 
	int PtoCMax, 
	int PtoCParamSize,
	int side
	);

//__device__ int GetPtoCIndx(int pIndx);
//__device__ int GetPtoC(int pIndx, int lIndx, int oIndx);

void LaunchCalcAverageGPU(int prtNum, float* sldPrtPos, float* sldPrtVel, float* sphPrtPos, float* sphPrtVel, float* smPrtPos, float* smPrtVel, int* smIndxSet, int* PtoCIndx, int* PtoC, int PNumMax, int PtoCMax, int PtoCParamSize)
{
	int side = pow( prtNum, 1.0/3.0 ) + 0.5;	//�����̂̂P�ӂ̒��_��

	dim3 grid(side, side);
	dim3 block(side, 1, 1);

	//�^���v�Z
	CalcAverage<<<grid ,block>>>(sldPrtPos, sldPrtVel, sphPrtPos, sphPrtVel, smPrtPos, smPrtVel, smIndxSet, PtoCIndx, PtoC, PNumMax, PtoCMax, PtoCParamSize, side);
	
	cudaThreadSynchronize();
}

__global__
	void CalcAverage(float* sldPrtPos, float* sldPrtVel, float* sphPrtPos, float* sphPrtVel, float* smPrtPos, float* smPrtVel, int* indxSet, int* PtoCIndx, int* PtoC, int PNumMax, int PtoCMax, int PtoCParamSize, int side)
{
	//�v�Z���闱�q�̔���
	int pIndx = blockIdx.x * side * side + blockIdx.y * side + threadIdx.x;

	//���ꂼ��̃x�N�g�������������ς��Ƃ�
	float pos_x, pos_y, pos_z;
	float vel_x, vel_y, vel_z;

	float clusterNum = 0.0f;					//�N���X�^�̐�

	int pTocIndx = PtoCIndx[pIndx];

	for(int j = 0; j < pTocIndx; ++j)
	{
		//pIndx�Ԗڂ̗��q��������j�ڂ̃N���X�^
		int jcIndx = PtoC[(PtoCMax*PtoCParamSize) * pIndx + PtoCMax * 0 + j];
		int joIndx = PtoC[(PtoCMax*PtoCParamSize) * pIndx + PtoCMax * 1 + j];

		if(jcIndx == -1 || joIndx == -1){	continue;	}

		//�N���X�^jcIndx�Ԃ�joIndx�Ԗڂ̒��_
		int startIndx = indxSet[jcIndx*2+0];

		int smIndx = startIndx*3 + joIndx*3;

		pos_x += smPrtPos[smIndx+0];
		pos_y += smPrtPos[smIndx+1];
		pos_z += smPrtPos[smIndx+2];

		vel_x += smPrtVel[smIndx+0];
		vel_y += smPrtVel[smIndx+1];
		vel_z += smPrtVel[smIndx+2];

		clusterNum += 1.0f;
	}

	//�N���X�^�̐��Ŋ���
	if(clusterNum != 0.0f)
	{
		pos_x *= 1/clusterNum;
		pos_y *= 1/clusterNum;
		pos_z *= 1/clusterNum;

		vel_x *= 1/clusterNum;
		vel_y *= 1/clusterNum;
		vel_z *= 1/clusterNum;
	}		

	//�ő̂̍ŏI�I�ȉ^���v�Z����
	int sldIndx = pIndx*3;

	sldPrtPos[sldIndx+0] = pos_x;
	sldPrtPos[sldIndx+1] = pos_y;
	sldPrtPos[sldIndx+2] = pos_z;

	sldPrtVel[sldIndx+0] = vel_x;
	sldPrtVel[sldIndx+1] = vel_y;
	sldPrtVel[sldIndx+2] = vel_z;
}

//__device__ int GetPtoCIndx(int pIndx)
//{
//	return 0;
//}
//
//__device__ int GetPtoC(int pIndx, int lIndx, int oIndx)
//{
//	return 0;
//}


#endif