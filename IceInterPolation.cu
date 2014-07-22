//�t�̂ƌő̂̉^���v�Z���ʂ��Ԃ��C�ŏI�I�Ȍ��ʂ��v�Z
//���͒P���ɐ��`���

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

	int n = pow(prtNum, 1.0/3.0) + 0.5;	//�����̂̂P�ӂ̒��_��

	dim3 grid(n, n);
	dim3 block(n, 1, 1);

	//���`���
	LinerInterPolation<<<grid ,block>>>(sldPrtPos, sldPrtVel, sphPrtPos, sphPrtVel);

	cudaThreadSynchronize();
}

__global__
	void LinerInterPolation(float* sldPrtPos, float* sldPrtVel, float* sphPrtPos, float* sphPrtVel)
{
	//�v�Z���闱�q�̔���
	int pIndx = blockIdx.x * EDGE * EDGE + blockIdx.y * EDGE + threadIdx.x;

	//���`���
	sphPrtPos[pIndx*4+0] = sldPrtPos[pIndx*3+0];
	sphPrtPos[pIndx*4+1] = sldPrtPos[pIndx*3+1];
	sphPrtPos[pIndx*4+2] = sldPrtPos[pIndx*3+2];

	sphPrtVel[pIndx*4+0] = sldPrtVel[pIndx*3+0];
	sphPrtVel[pIndx*4+1] = sldPrtVel[pIndx*3+1];
	sphPrtVel[pIndx*4+2] = sldPrtVel[pIndx*3+2];
}





#endif