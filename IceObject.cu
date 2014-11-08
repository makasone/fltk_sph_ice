#ifndef _GPU_ICE_OBJECT_
#define _GPU_ICE_OBJECT_

#include <iostream>
#include <cuda_runtime.h>
#include <rx_cu_common.cuh>

//#include "Path_SM.cuh"

__global__ void UpdateSMFromPath(
	float* prtPos,
	float* prtVel,
	//-----------------------SM---------------------
	float* orgPos,
	float* curPos,
	float* orgCm,
	float* curCm,
	float* clstrApq,
	float* vel,
	int* pIndxes, 
	int* indxSet,
	//-----------------------Path-------------------
	short int* PRTtoPTH,
	short int* PTHandPrfxSet,
	float* prfxPos,
	float* prfxApq,
	//----------------------Struct-------------------
	int* CtoP,
	int* CtoPNum,
	int CtoPSizeY,
	int CtoPSizeZ,

	float dt,
	int prtNum,
	int side
);

__device__ void CalcCmSum
(
	float* curCm,
	float* prfxPos,
	short int* PRTtoPTH,
	short int* PTHandPrfxSet,
	int* CtoP,
	int* CtoPNum,
	int CtoPSizeY,
	int CtoPSizeZ,
	int prtNum
);

__device__ void CalcApqSum
(
	float* orgCmes,
	float* curCmes,
	float* prfxApq,
	float* clstrApq,
	short int* PRTtoPTH,
	short int* PTHandPrfxSet,
	int* CtoP,
	int* CtoPNumes,
	int CtoPSizeY,
	int CtoPSizeZ,
	int prtNum
);

void LauchUpdateSMFromPath(
	int prtNum, float* prtPos, float* prtVel, 
	float* orgPos, float* curPos, float* orgCm, float* curCm, float* clstrApq, float* vel, int* pIndxes, int* indxSet,
	short int* PRTtoPTH,
	short int* PTHandPrfxSet,
	float* prfxPos,
	float* prfxApq,
	int* CtoP,
	int* CtoPNum,
	int CtoPSizeY,
	int CtoPSizeZ,
	float dt)
{	//cout << __FUNCTION__ << endl;

	int side = pow(prtNum, 1.0/3.0) + 0.5;	//立方体の１辺の頂点数

	dim3 grid(side, side);
	dim3 block(side, 1, 1);

	//
	UpdateSMFromPath<<<grid ,block>>>(
		prtPos, prtVel,
		orgPos, curPos, orgCm, curCm, clstrApq, vel, pIndxes, indxSet,
		PRTtoPTH,
		PTHandPrfxSet,
		prfxPos,
		prfxApq,
		CtoP,
		CtoPNum,
	    CtoPSizeY,
	    CtoPSizeZ,
		dt, prtNum, side
	);

	cudaThreadSynchronize();
}


__global__
	void UpdateSMFromPath(
	float* prtPos, float* prtVel,
	float* orgPos, float* curPos, float* orgCm, float* curCm, float* clstrApq, float* vel, int* pIndxes, int* indxSet,
	short int* PRTtoPTH,
	short int* PTHandPrfxSet,
	float* prfxPos,
	float* prfxApq,
	int* CtoP,
	int* CtoPNum,
	int CtoPSizeY,
	int CtoPSizeZ,
	float dt, int prtNum, int side)
{
	//TODO: ここからPath_SM.cuhをインクルードして，CalcCmSumなんかを呼びだそうとしたが，うまくいかなかった
	
	CalcCmSum(curCm, prfxPos, PRTtoPTH, PTHandPrfxSet, CtoP, CtoPNum, CtoPSizeY, CtoPSizeZ, prtNum);
	CalcApqSum(orgCm, curCm, prfxApq, clstrApq, PRTtoPTH, PTHandPrfxSet, CtoP, CtoPNum, CtoPSizeY, CtoPSizeZ, prtNum);
}

//TODO: pathの数を１本として固定しているのに注意
__device__
void CalcCmSum
(
	float* curCm,
	float* prfxPos,
	short int* PRTtoPTH,
	short int* PTHandPrfxSet,
	int* CtoP,
	int* CtoPNumes,
	int CtoPSizeY,
	int CtoPSizeZ,
	int prtNum
)
{
	int cIndx = blockIdx.x * (blockDim.x * blockDim.y) + blockIdx.y * gridDim.x * (blockDim.x * blockDim.y) + threadIdx.x;
	int CtoPNum = CtoPNumes[cIndx];
	
	float3 cmSum;
	cmSum.x = 0.0f;	cmSum.y = 0.0f;	cmSum.z = 0.0f;

	int PTHandPrfxSet_S = (prtNum * 2) * cIndx + prtNum * 0;
	int PTHandPrfxSet_E = (prtNum * 2) * cIndx + prtNum * 1;

	int CtoPIndx = (CtoPSizeY * CtoPSizeZ) * cIndx + CtoPSizeY * 0;

	for(int iPrt = 0; iPrt < CtoPNum; iPrt++)
	{
		int start	= PTHandPrfxSet[PTHandPrfxSet_S + iPrt];
		int end		= PTHandPrfxSet[PTHandPrfxSet_E + iPrt];

		if(start == -1 || end == -1){	break;	}

		int prtIndx = CtoP[CtoPIndx + iPrt];
		int pthIndx = PRTtoPTH[prtNum * 0 + prtIndx];

		//prfxSumはデータの並びが普通のベクトル配列と違うのに注意
		int prfxEndIndx		= (end * 1 + pthIndx);
		int prfxStartIndx	= ((start-1) * 1 + pthIndx);

		if(start == 0)
		{
			cmSum.x += prfxPos[prtNum*0 + prfxEndIndx];
			cmSum.y += prfxPos[prtNum*1 + prfxEndIndx];
			cmSum.z += prfxPos[prtNum*2 + prfxEndIndx];
		}
		else
		{
			cmSum.x += prfxPos[prtNum*0 + prfxEndIndx] - prfxPos[prtNum*0 + prfxStartIndx];
			cmSum.y += prfxPos[prtNum*1 + prfxEndIndx] - prfxPos[prtNum*1 + prfxStartIndx];
			cmSum.z += prfxPos[prtNum*2 + prfxEndIndx] - prfxPos[prtNum*2 + prfxStartIndx];
		}
	}

	curCm[cIndx*3+0] = cmSum.x;
	curCm[cIndx*3+1] = cmSum.y;
	curCm[cIndx*3+2] = cmSum.z;
}

__device__
void CalcApqSum
(
	float* orgCmes,
	float* curCm,
	float* prfxApq,
	float* clstrApq,
	short int* PRTtoPTH,
	short int* PTHandPrfxSet,
	int* CtoP,
	int* CtoPNumes,
	int CtoPSizeY,
	int CtoPSizeZ,
	int prtNum
)
{
	int cIndx = blockIdx.x * (blockDim.x * blockDim.y) + blockIdx.y * gridDim.x * (blockDim.x * blockDim.y) + threadIdx.x;
	int CtoPNum = CtoPNumes[cIndx];

	float3 cmSum;
	cmSum.x = curCm[cIndx*3+0] / CtoPNum;
	cmSum.y = curCm[cIndx*3+1] / CtoPNum;
	cmSum.z = curCm[cIndx*3+2] / CtoPNum;

	float3 orgCm;
	orgCm.x = orgCmes[cIndx*3+0];
	orgCm.y = orgCmes[cIndx*3+1];
	orgCm.z = orgCmes[cIndx*3+2];

	//mtt0T
	matrix3x3 mtt0T;
	mtt0T.e[0].x = cmSum.x * orgCm.x;
	mtt0T.e[0].y = cmSum.x * orgCm.y;
	mtt0T.e[0].z = cmSum.x * orgCm.z;

	mtt0T.e[1].x = cmSum.y * orgCm.x;
	mtt0T.e[1].y = cmSum.y * orgCm.y;
	mtt0T.e[1].z = cmSum.y * orgCm.z;

	mtt0T.e[2].x = cmSum.z * orgCm.x;
	mtt0T.e[2].y = cmSum.z * orgCm.y;
	mtt0T.e[2].z = cmSum.z * orgCm.z;

	matrix3x3 ApqSum;
	ApqSum.e[0].x = -CtoPNum * mtt0T.e[0].x;
	ApqSum.e[0].y = -CtoPNum * mtt0T.e[0].y;
	ApqSum.e[0].z = -CtoPNum * mtt0T.e[0].z;
	
	ApqSum.e[1].x = -CtoPNum * mtt0T.e[1].x;
	ApqSum.e[1].y = -CtoPNum * mtt0T.e[1].y;
	ApqSum.e[1].z = -CtoPNum * mtt0T.e[1].z;
	
	ApqSum.e[2].x = -CtoPNum * mtt0T.e[2].x;
	ApqSum.e[2].y = -CtoPNum * mtt0T.e[2].y;
	ApqSum.e[2].z = -CtoPNum * mtt0T.e[2].z;

	int PTHandPrfx_S = (prtNum * 2) * cIndx + prtNum * 0;
	int PTHandPrfx_E = (prtNum * 2) * cIndx + prtNum * 1;

	int CtoPIndx = (CtoPSizeY * CtoPSizeZ) * cIndx + CtoPSizeY * 0;

	for(int iPrt = 0; iPrt < CtoPNum; iPrt++)
	{
		int start	= PTHandPrfxSet[PTHandPrfx_S + iPrt];
		int end		= PTHandPrfxSet[PTHandPrfx_E + iPrt];

		if(start == -1 || end == -1){	break;	}

		int prtIndx = CtoP[CtoPIndx + iPrt];
		int pthIndx = PRTtoPTH[prtNum * 0 + prtIndx];

		int prfxEndIndx		= (end * 1 + pthIndx);
		int prfxStartIndx	= ((start-1) * 1 + pthIndx);

		if(start == 0)
		{
			ApqSum.e[0].x += prfxApq[prtNum*0 + prfxEndIndx];
			ApqSum.e[0].y += prfxApq[prtNum*1 + prfxEndIndx];
			ApqSum.e[0].z += prfxApq[prtNum*2 + prfxEndIndx];
			
			ApqSum.e[1].x += prfxApq[prtNum*3 + prfxEndIndx];
			ApqSum.e[1].y += prfxApq[prtNum*4 + prfxEndIndx];
			ApqSum.e[1].z += prfxApq[prtNum*5 + prfxEndIndx];
			
			ApqSum.e[2].x += prfxApq[prtNum*6 + prfxEndIndx];
			ApqSum.e[2].y += prfxApq[prtNum*7 + prfxEndIndx];
			ApqSum.e[2].z += prfxApq[prtNum*8 + prfxEndIndx];
		}
		else
		{
			ApqSum.e[0].x += prfxApq[prtNum*0 + prfxEndIndx] - prfxApq[prtNum*0 + prfxStartIndx];
			ApqSum.e[0].y += prfxApq[prtNum*1 + prfxEndIndx] - prfxApq[prtNum*1 + prfxStartIndx];
			ApqSum.e[0].z += prfxApq[prtNum*2 + prfxEndIndx] - prfxApq[prtNum*2 + prfxStartIndx];
			
			ApqSum.e[1].x += prfxApq[prtNum*3 + prfxEndIndx] - prfxApq[prtNum*3 + prfxStartIndx];
			ApqSum.e[1].y += prfxApq[prtNum*4 + prfxEndIndx] - prfxApq[prtNum*4 + prfxStartIndx];
			ApqSum.e[1].z += prfxApq[prtNum*5 + prfxEndIndx] - prfxApq[prtNum*5 + prfxStartIndx];
	
			ApqSum.e[2].x += prfxApq[prtNum*6 + prfxEndIndx] - prfxApq[prtNum*6 + prfxStartIndx];
			ApqSum.e[2].y += prfxApq[prtNum*7 + prfxEndIndx] - prfxApq[prtNum*7 + prfxStartIndx];
			ApqSum.e[2].z += prfxApq[prtNum*8 + prfxEndIndx] - prfxApq[prtNum*8 + prfxStartIndx];
		}
	}

	int clstrIndx = cIndx * 9;
	clstrApq[clstrIndx+0] = ApqSum.e[0].x;
	clstrApq[clstrIndx+1] = ApqSum.e[0].y;
	clstrApq[clstrIndx+2] = ApqSum.e[0].z;
				
	clstrApq[clstrIndx+3] = ApqSum.e[1].x;
	clstrApq[clstrIndx+4] = ApqSum.e[1].y;
	clstrApq[clstrIndx+5] = ApqSum.e[1].z;
				
	clstrApq[clstrIndx+6] = ApqSum.e[2].x;
	clstrApq[clstrIndx+7] = ApqSum.e[2].y;
	clstrApq[clstrIndx+8] = ApqSum.e[2].z;
}

#endif