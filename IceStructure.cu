//各クラスタの情報から最終的な固体の情報を算出
//今は単純に平均

#ifndef _GPU_ICESTRUCTURE_H_
#define _GPU_ICESTRUCTURE_H_

void LaunchCalcAverage();
__global__ void CalcAverage();		//実験として補間までいれちゃう
__device__ int GetPtoC(int pIndx, int lIndx, int oIndx);

void LaunchCalcAverage()
{
	dim3 grid(1, 1);
	dim3 block(729, 1, 1);

	//運動計算
//	CalcAverage<<<grid ,block>>>();
}






#endif