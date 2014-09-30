//液体と固体の運動計算結果を補間し，最終的な結果を計算
//HACK::今は単純に線形補間
//HACK::立方体しかできないのに注意

#ifndef _GPU_ICE_INTERPOLATION_H_
#define _GPU_ICE_INTERPOLATION_H_

#include <iostream>
#include <cuda_runtime.h>


using namespace std;


void LaunchInterPolationGPU(int prtNum, float* sldPrtPos, float* sldPrtVel, float* sphPrtPos, float* sphPrtVel);
__global__ void LinerInterPolation(float* sldPrtPos, float* sldPrtVel, float* sphPrtPos, float* sphPrtVel, int side);


void LaunchInterPolationGPU(int prtNum, float* sldPrtPos, float* sldPrtVel, float* sphPrtPos, float* sphPrtVel)
{	//cout << __FUNCTION__ << endl;

	int side = pow(prtNum, 1.0/3.0) + 0.5;	//立方体の１辺の頂点数

	dim3 grid(side, side);
	dim3 block(side, 1, 1);

	//線形補間
	LinerInterPolation<<<grid ,block>>>(sldPrtPos, sldPrtVel, sphPrtPos, sphPrtVel, side);

	cudaThreadSynchronize();
}

__global__
	void LinerInterPolation(float* sldPrtPos, float* sldPrtVel, float* sphPrtPos, float* sphPrtVel, int side)
{
	//計算する粒子の判定
	int pIndx = blockIdx.x * side * side + blockIdx.y * side + threadIdx.x;
	int sphIndx = pIndx * 4;
	int sldIndx = pIndx * 3;

	//線形補間
	sphPrtPos[sphIndx+0] = sldPrtPos[sldIndx+0];
	sphPrtPos[sphIndx+1] = sldPrtPos[sldIndx+1];
	sphPrtPos[sphIndx+2] = sldPrtPos[sldIndx+2];

	sphPrtVel[sphIndx+0] = sldPrtVel[sldIndx+0];
	sphPrtVel[sphIndx+1] = sldPrtVel[sldIndx+1];
	sphPrtVel[sphIndx+2] = sldPrtVel[sldIndx+2];
}





#endif