//�t�̂ƌő̂̉^���v�Z���ʂ��Ԃ��C�ŏI�I�Ȍ��ʂ��v�Z
//���͒P���ɐ��`���

#ifndef _GPU_ICE_INTERPOLATION_H_
#define _GPU_ICE_INTERPOLATION_H_

#include <iostream>

using namespace std;

void LaunchInterPolationGPU();
__global__ void CalcAverage();		//�����Ƃ��ĕ�Ԃ܂ł��ꂿ�Ⴄ
__device__ int GetPtoC(int pIndx, int lIndx, int oIndx);

void LaunchInterPolationGPU()
{	//cout << __FUNCTION__ << endl;
	return;
	dim3 grid(1, 1);
	dim3 block(729, 1, 1);

	//�^���v�Z
//	CalcAverage<<<grid ,block>>>();
}






#endif