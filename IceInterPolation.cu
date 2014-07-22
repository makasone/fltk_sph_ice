//液体と固体の運動計算結果を補間し，最終的な結果を計算
//今は単純に線形補間

#ifndef _GPU_ICE_INTERPOLATION_H_
#define _GPU_ICE_INTERPOLATION_H_

#include <iostream>

using namespace std;

#define EDGE 17
//#define EDGE 27

void LaunchInterPolationGPU(int prtNum, float* sldPrtPos, float* sldPrtVel, float* sphPrtPos, float* sphPrtVel);

__global__ void LinerInterPolation(float* sldPrtPos, float* sldPrtVel, float* sphPrtPos, float* sphPrtVel);

void LaunchInterPolationGPU(int prtNum, float* sldPrtPos, float* sldPrtVel, float* sphPrtPos, float* sphPrtVel)
{	//cout << __FUNCTION__ << endl;

	int n = pow(prtNum, 1.0/3.0) + 0.5;	//立方体の１辺の頂点数

	dim3 grid(n, n);
	dim3 block(n, 1, 1);

	//線形補間
	LinerInterPolation<<<grid ,block>>>(sldPrtPos, sldPrtVel, sphPrtPos, sphPrtVel);

	cudaThreadSynchronize();
}

__global__
	void LinerInterPolation(float* sldPrtPos, float* sldPrtVel, float* sphPrtPos, float* sphPrtVel)
{
	//計算する粒子の判定
	int pIndx = blockIdx.x * EDGE * EDGE + blockIdx.y * EDGE + threadIdx.x;

	//線形補間
	sphPrtPos[pIndx*4+0] = sldPrtPos[pIndx*3+0];
	sphPrtPos[pIndx*4+1] = sldPrtPos[pIndx*3+1];
	sphPrtPos[pIndx*4+2] = sldPrtPos[pIndx*3+2];

	sphPrtVel[pIndx*4+0] = sldPrtVel[pIndx*3+0];
	sphPrtVel[pIndx*4+1] = sldPrtVel[pIndx*3+1];
	sphPrtVel[pIndx*4+2] = sldPrtVel[pIndx*3+2];
}





#endif