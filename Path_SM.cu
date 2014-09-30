//
#ifndef _GPU_PATH_SM_H_
#define _GPU_PATH_SM_H_

#include <iostream>
#include <cuda_runtime.h>

using namespace std;


void LaunchUpdatePrefixSumGPU(
	int prtNum,
	int PosSizeX,
	int PosSizeY,
	int ApqSizeX,
	int ApqSizeY,
	int* md_2DiPTHtoPRT,
	int* md_2DiPRTtoPTH,
	float* md_2Df3PrfxPos,
	float* md_2Df9PrfxApq,
	int* md_3DiPTHandPrfxSet,
	float* md_f3OrgPos,
	float* md_f3OrgCm,
	const float* md_fPos,
	const float* md_fVel
	);

__global__ void UpdatePrefixSumGPU(
	int prtNum,
	int PosSizeX,
	int PosSizeY,
	int ApqSizeX,
	int ApqSizeY,
	int* md_2DiPTHtoPRT,
	int* md_2DiPRTtoPTH,
	float* md_2Df3PrfxPos,
	float* md_2Df9PrfxApq,
	int* md_3DiPTHandPrfxSet,
	float* md_f3OrgPos,
	float* md_f3OrgCm,
	const float* md_fPos,
	const float* md_fVel
	);

__device__ void UpdatePrefixSumPosGPU();
__device__ void UpdatePrefixSumApqGPU();

void LaunchUpdatePrefixSumGPU(
	int prtNum,
	int PosSizeX,
	int PosSizeY,
	int ApqSizeX,
	int ApqSizeY,
	int* md_2DiPTHtoPRT,
	int* md_2DiPRTtoPTH,
	float* md_2Df3PrfxPos,
	float* md_2Df9PrfxApq,
	int* md_3DiPTHandPrfxSet,
	float* md_f3OrgPos,
	float* md_f3OrgCm,
	const float* md_fPos,
	const float* md_fVel
	)
{	//cout << __FUNCTION__ << endl;

	int side = pow(prtNum, 1.0/3.0) + 0.5;	//—§•û‘Ì‚Ì‚P•Ó‚Ì’¸“_”

	dim3 grid(side, side);
	dim3 block(side, 1, 1);

	//üŒ`•âŠÔ
	UpdatePrefixSumGPU<<<grid ,block>>>
		(
		prtNum,
		PosSizeX, PosSizeY,
		ApqSizeX, ApqSizeY,
		md_2DiPTHtoPRT, md_2DiPRTtoPTH,
		md_2Df3PrfxPos, md_2Df9PrfxApq,
		md_3DiPTHandPrfxSet,
		md_f3OrgPos,	md_f3OrgCm,
		md_fPos,		md_fVel
		);

	cudaThreadSynchronize();
}

__global__
void UpdatePrefixSumGPU(
	int prtNum,
	int PosSizeX,
	int PosSizeY,
	int ApqSizeX,
	int ApqSizeY,
	int* md_2DiPTHtoPRT,
	int* md_2DiPRTtoPTH,
	float* md_2Df3PrfxPos,
	float* md_2Df9PrfxApq,
	int* md_3DiPTHandPrfxSet,
	float* md_f3OrgPos,
	float* md_f3OrgCm,
	const float* md_fPos,
	const float* md_fVel
	)
{
	UpdatePrefixSumPosGPU();
	UpdatePrefixSumApqGPU();
}

__device__
void UpdatePrefixSumPosGPU()
{
}

__device__
void UpdatePrefixSumApqGPU()
{
}

#endif