//�e�N���X�^�̏�񂩂�ŏI�I�Ȍő̂̏����Z�o
//���͒P���ɕ���

#ifndef _GPU_ICESTRUCTURE_H_
#define _GPU_ICESTRUCTURE_H_

void LaunchCalcAverage();
__global__ void CalcAverage();		//�����Ƃ��ĕ�Ԃ܂ł��ꂿ�Ⴄ
__device__ int GetPtoC(int pIndx, int lIndx, int oIndx);

void LaunchCalcAverage()
{
	dim3 grid(1, 1);
	dim3 block(729, 1, 1);

	//�^���v�Z
//	CalcAverage<<<grid ,block>>>();
}






#endif