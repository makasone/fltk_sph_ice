#ifndef _GPU_PATH_SM_H_
#define _GPU_PATH_SM_H_

#include <iostream>
#include <cuda_runtime.h>
#include <thrust\scan.h>
//#include <cuda\cudpp-2.2\include\cudpp.h>

using namespace std;

#define TESTDATANUM 1024


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

void ThrustTest();
void ThrustTest2();

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

__device__ void UpdatePrefixSumPosGPU(
	int* PTHtoPRT,
	float* prfixPosOut,
	const float* pos,
	const float* vel,
	int num
);

__device__ void UpdatePrefixSumApqGPU(
	int* PTHtoPRT,
	float* prfixPosOut,
	float* orgPos,
	const float* pos,
	const float* vel,
	int num
);

//__device__ void CalcCmSum
//(
//
//);
//
//__device__ void CalcApqSum
//(
//
//);

__device__ void scan_A(int* inData, int* outData, int num);
__device__ void scan_B(int* inData, int* outData, int num);
__device__ void scan_C(int* inData, int* outData, int num);
__device__ void scan_D(int* inData, int* outData, int num);

#endif