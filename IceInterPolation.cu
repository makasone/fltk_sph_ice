//液体と固体の運動計算結果を補間し，最終的な結果を計算
//今は単純に線形補間

#ifndef _GPU_ICE_INTERPOLATION_H_
#define _GPU_ICE_INTERPOLATION_H_

#include <iostream>

using namespace std;

void LaunchInterPolationGPU();
__global__ void CalcAverage();		//実験として補間までいれちゃう
__device__ int GetPtoC(int pIndx, int lIndx, int oIndx);

void LaunchInterPolationGPU()
{	//cout << __FUNCTION__ << endl;
	return;
	dim3 grid(1, 1);
	dim3 block(729, 1, 1);

	//運動計算
//	CalcAverage<<<grid ,block>>>();
}






#endif