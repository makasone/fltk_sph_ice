//各クラスタの情報から最終的な固体の情報を算出
//今は単純に平均

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
	int side = pow( prtNum, 1.0/3.0 ) + 0.5;	//立方体の１辺の頂点数

	dim3 grid(side, side);
	dim3 block(side, 1, 1);

	//運動計算
	CalcAverage<<<grid ,block>>>(sldPrtPos, sldPrtVel, sphPrtPos, sphPrtVel, smPrtPos, smPrtVel, smIndxSet, PtoCIndx, PtoC, PNumMax, PtoCMax, PtoCParamSize, side);
	
	cudaThreadSynchronize();
}

__global__
	void CalcAverage(float* sldPrtPos, float* sldPrtVel, float* sphPrtPos, float* sphPrtVel, float* smPrtPos, float* smPrtVel, int* indxSet, int* PtoCIndx, int* PtoC, int PNumMax, int PtoCMax, int PtoCParamSize, int side)
{
	//計算する粒子の判定
	int pIndx = blockIdx.x * side * side + blockIdx.y * side + threadIdx.x;

	//それぞれのベクトルを合成し平均をとる
	float pos_x, pos_y, pos_z;
	float vel_x, vel_y, vel_z;

	float clusterNum = 0.0f;					//クラスタの数

	int pTocIndx = PtoCIndx[pIndx];

	for(int j = 0; j < pTocIndx; ++j)
	{
		//pIndx番目の粒子が属するj個目のクラスタ
		int jcIndx = PtoC[(PtoCMax*PtoCParamSize) * pIndx + PtoCMax * 0 + j];
		int joIndx = PtoC[(PtoCMax*PtoCParamSize) * pIndx + PtoCMax * 1 + j];

		if(jcIndx == -1 || joIndx == -1){	continue;	}

		//クラスタjcIndx番のjoIndx番目の頂点
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

	//クラスタの数で割る
	if(clusterNum != 0.0f)
	{
		pos_x *= 1/clusterNum;
		pos_y *= 1/clusterNum;
		pos_z *= 1/clusterNum;

		vel_x *= 1/clusterNum;
		vel_y *= 1/clusterNum;
		vel_z *= 1/clusterNum;
	}		

	//固体の最終的な運動計算結果
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