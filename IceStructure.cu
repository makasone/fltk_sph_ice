//各クラスタの情報から最終的な固体の情報を算出
//今は単純に平均

#ifndef _GPU_ICE_STRUCTURE_H_
#define _GPU_ICE_STRUCTURE_H_

#include <iostream>
#include <cuda_runtime.h>
#include <rx_cu_common.cuh>

using namespace std;

#define EDGE 17

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

__global__ void CalcAverage(float* sldPrtPos, float* sldPrtVel, float* sphPrtPos, float* sphPrtVel, float* smPrtPos, float* smPrtVel, int* indxSet, int* PtoCIndx, int* PtoC, int PNumMax, int PtoCMax, int PtoCParamSize);

//__device__ int GetPtoCIndx(int pIndx);
//__device__ int GetPtoC(int pIndx, int lIndx, int oIndx);

void LaunchCalcAverageGPU(int prtNum, float* sldPrtPos, float* sldPrtVel, float* sphPrtPos, float* sphPrtVel, float* smPrtPos, float* smPrtVel, int* smIndxSet, int* PtoCIndx, int* PtoC, int PNumMax, int PtoCMax, int PtoCParamSize)
{
	int n = pow( prtNum, 1.0/3.0 ) + 0.5;	//立方体の１辺の頂点数

	dim3 grid(n, n);
	dim3 block(n, 1, 1);

	//運動計算
	CalcAverage<<<grid ,block>>>(sldPrtPos, sldPrtVel, sphPrtPos, sphPrtVel, smPrtPos, smPrtVel, smIndxSet, PtoCIndx, PtoC, PNumMax, PtoCMax, PtoCParamSize);
}

__global__
	void CalcAverage(float* sldPrtPos, float* sldPrtVel, float* sphPrtPos, float* sphPrtVel, float* smPrtPos, float* smPrtVel, int* indxSet, int* PtoCIndx, int* PtoC, int PNumMax, int PtoCMax, int PtoCParamSize)
{
	//計算する粒子の判定
	int pIndx = blockIdx.x * EDGE * EDGE + blockIdx.y * EDGE + threadIdx.x;

	//それぞれのベクトルを合成し平均をとる
	float3 pos = make_float3(0.0f, 0.0f, 0.0f);
	float3 vel = make_float3(0.0f, 0.0f, 0.0f);
	float clusterNum = 0.0f;					//クラスタの数

	//TODO::開始添字，終了添字
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

		pos.x += smPrtPos[smIndx+0];
		pos.y += smPrtPos[smIndx+1];
		pos.z += smPrtPos[smIndx+2];

		vel.x += smPrtVel[smIndx+0];
		vel.y += smPrtVel[smIndx+1];
		vel.z += smPrtVel[smIndx+2];

		clusterNum += 1.0f;
	}

	//クラスタの数で割る
	if(clusterNum != 0.0f)
	{
		pos.x *= 1/clusterNum;
		pos.y *= 1/clusterNum;
		pos.z *= 1/clusterNum;

		vel.x *= 1/clusterNum;
		vel.y *= 1/clusterNum;
		vel.z *= 1/clusterNum;

		//pos.x /= clusterNum;
		//pos.y /= clusterNum;
		//pos.z /= clusterNum;

		//vel.x /= clusterNum;
		//vel.y /= clusterNum;
		//vel.z /= clusterNum;
	}		
	//どのクラスタにも含まれていない場合，運動はSPH法に従う
	else
	{
		int sphIndx = pIndx*4;

		pos.x = sphPrtPos[sphIndx+0];
		pos.y = sphPrtPos[sphIndx+1];
		pos.z = sphPrtPos[sphIndx+2];

		vel.x = sphPrtVel[sphIndx+0];
		vel.y = sphPrtVel[sphIndx+1];
		vel.z = sphPrtVel[sphIndx+2];
	}

	//固体の最終的な運動計算結果
	int sldIndx = pIndx*3;

	sldPrtPos[sldIndx+0] = pos.x;
	sldPrtPos[sldIndx+1] = pos.y;
	sldPrtPos[sldIndx+2] = pos.z;

	sldPrtVel[sldIndx+0] = vel.x;
	sldPrtVel[sldIndx+1] = vel.y;
	sldPrtVel[sldIndx+2] = vel.z;

	//適当に線形補間もどき
	sphPrtPos[pIndx*4+0] = pos.x;
	sphPrtPos[pIndx*4+1] = pos.y;
	sphPrtPos[pIndx*4+2] = pos.z;

	sphPrtVel[pIndx*4+0] = vel.x;
	sphPrtVel[pIndx*4+1] = vel.y;
	sphPrtVel[pIndx*4+2] = vel.z;
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