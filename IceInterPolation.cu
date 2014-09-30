//�t�̂ƌő̂̉^���v�Z���ʂ��Ԃ��C�ŏI�I�Ȍ��ʂ��v�Z
//HACK::���͒P���ɐ��`���
//HACK::�����̂����ł��Ȃ��̂ɒ���

#ifndef _GPU_ICE_INTERPOLATION_H_
#define _GPU_ICE_INTERPOLATION_H_

#include <iostream>
#include <cuda_runtime.h>


using namespace std;


void LaunchInterPolationGPU(int prtNum, float* sldPrtPos, float* sldPrtVel, float* sphPrtPos, float* sphPrtVel);
__global__ void LinerInterPolation(float* sldPrtPos, float* sldPrtVel, float* sphPrtPos, float* sphPrtVel, int side);


void LaunchInterPolationGPU(int prtNum, float* sldPrtPos, float* sldPrtVel, float* sphPrtPos, float* sphPrtVel)
{	//cout << __FUNCTION__ << endl;

	int side = pow(prtNum, 1.0/3.0) + 0.5;	//�����̂̂P�ӂ̒��_��

	dim3 grid(side, side);
	dim3 block(side, 1, 1);

	//���`���
	LinerInterPolation<<<grid ,block>>>(sldPrtPos, sldPrtVel, sphPrtPos, sphPrtVel, side);

	cudaThreadSynchronize();
}

__global__
	void LinerInterPolation(float* sldPrtPos, float* sldPrtVel, float* sphPrtPos, float* sphPrtVel, int side)
{
	//�v�Z���闱�q�̔���
	int pIndx = blockIdx.x * side * side + blockIdx.y * side + threadIdx.x;
	int sphIndx = pIndx * 4;
	int sldIndx = pIndx * 3;

	//���`���
	sphPrtPos[sphIndx+0] = sldPrtPos[sldIndx+0];
	sphPrtPos[sphIndx+1] = sldPrtPos[sldIndx+1];
	sphPrtPos[sphIndx+2] = sldPrtPos[sldIndx+2];

	sphPrtVel[sphIndx+0] = sldPrtVel[sldIndx+0];
	sphPrtVel[sphIndx+1] = sldPrtVel[sldIndx+1];
	sphPrtVel[sphIndx+2] = sldPrtVel[sldIndx+2];
}





#endif