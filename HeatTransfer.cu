#ifndef _GPU_HEATTRANSFER_CU_
#define _GPU_HEATTRANSFER_CU_

#include <math.h>
#include <stdio.h>
#include <iostream>
#include <cuda_runtime.h>

using namespace std;

#include "rx_cu_common.cuh"	//�搶����`�����֗��@�\���g����

void LaunchHeatTransferGPU
(
	float* heats,
	float* temps,
	float* dtemps,
	const int* surfParticles,
	float airTemp,
	float cffCnHt,
	int prtNum
	//const vector<vector<rxNeigh>>& neights,
	//const vector<int>& objNeight,
	//float floor,
	//float effRadius,
	//const float* pos,
	//const float* dens
);

__global__ void UpdateHeatTransfar(float* heats, const float* temps, const int* surfParticles, float airTemp, float ht);
__global__ void UpdateAirAndParticle(float* heats, const float* temps, const int* surfParticles, float airTemp, float ht);
__global__ void UpdatePhase(float* heats, float* temps, float* dtemps);

void LaunchHeatTransferGPU(float* heats, float* temps, float* dtemps, const int* surfParticles, float airTemp, float cffCnHt, int prtNum
	/*, const vector<vector<rxNeigh>>& neights, const vector<int>& objNeight, float floor, float effRadius, const float* pos, const float* dens*/
	)
{
	uint numThreads = min(THREAD_NUM, prtNum);
	uint numBlocks = (prtNum % numThreads != 0) ? (prtNum / numThreads + 1) : (prtNum / numThreads);

	//�\�ʗ��q�Ƌ�C�C���q�Ɨ��q�C�M���Ɨ��q
	//UpdateHeatTransfar<<<numBlocks, numThreads>>>(heats, temps, surfParticles, airTemp, ht);
	UpdateAirAndParticle<<<numBlocks, numThreads>>>(heats, temps, surfParticles, airTemp, cffCnHt);
	UpdatePhase<<<numBlocks, numThreads>>>(heats, temps, dtemps);

	cudaThreadSynchronize();
}

//�M����
__global__
void UpdateAirAndParticle(float* heats, const float* temps, const int* surParticles, float airTemp, float ht)
{
	int indx = __umul24(blockIdx.x,blockDim.x)+threadIdx.x;

//�\�ʗ��q�Ƌ�C
	if( surParticles[indx] == -1 ) return;

	double airNum = 20.0 - (double)surParticles[indx];
	if(airNum < 0) airNum = 0.0;

	double surfaceArea = airNum / 20.0;								//��C�ƐG��Ă���\�ʐρ@0�`1.0 20�͓K��

	if( surfaceArea < 0.0) surfaceArea = 0.0;
	if( surfaceArea > 1.0) surfaceArea = 1.0;

	double qHeat = ht * (airTemp - temps[indx]) * surfaceArea;		//�j���[�g���̗�p�@���̎�����M�ʂ��v�Z
	heats[indx] += qHeat;											//�M�ʂ����Z
}

//���x�ƔM�ʂ̏����@���M�E���M�̌v�Z�C���ω�����
__global__
void UpdatePhase(float* heats, float* temps, float* dtemps)
{
	int indx = __umul24(blockIdx.x,blockDim.x)+threadIdx.x;

	////���ԏ�Ԃւ̕ω����o
	//if( mPhase[i] == -2 && mTemps[i] > 250.0f ){				//�X�̏ꍇ
	//	mPhase[i] = -1;											//�X���ԏ��
	//	mTemps[i] = 250.0f;
	//	mHeats[i] = 0;
	//}
	//else if( mPhase[i] == 2 && mTemps[i] < 250.0f ){			//���̏ꍇ
	//	mPhase[i] = 1;											//�����ԏ��
	//	mTemps[i] = 250.0f;
	//	mHeats[i] = mLatentHeat;								//�Z����M
	//}

	//���M�E���M�̌v�Z
	if(false/*mPhase[i] == -1 || mPhase[i] == 1*/){					//�����Ԃ��X���Ԃ̏ꍇ
		////���M�v�Z
		////�ω����x��M�ʂɕϊ����Đ��M�v�Z�@�i���x�͕ω����Ȃ��j
		//float spcfHt = 2.1f + (2.1f * mHeats[i] / mLatentHeat);	//��M�@�ܗL�M�ʂŐ��ƕX�̔�M����
		//mHeats[i]	+= mTempsDelta[i] * spcfHt * 1.0f;			//���x�ω���M�ʂɊ��Z�@���ʂ͂P�ŌŒ�

		////���M�ω����猰�M�ω��֖߂锻��
		//if(mPhase[i] == -1 && mHeats[i] < -mLatentHeat/300.0f){	//���P�@<0�@�ɂ���ƁC�n���Ȃ��Ȃ�
		//	//�X���ԏ�ԁ��X
		//	mPhase[i] = -2;										//�X�ւ̑��ω�
		//	mTemps[i] = 249.0f;
		//	mHeats[i] = 0;										//�M�i�|�e���V�����G�l���M�[�j����o
		//}
		//else if( mPhase[i] == 1 && mHeats[i] > mLatentHeat ){
		//	//�����ԏ�ԁ���
		//	mPhase[i] = 2;										//���ւ̑��ω�
		//	mTemps[i] = 251.0f;
		//	mHeats[i] = 0;										//�M�i�|�e���V�����G�l���M�[�j���z��
		//}

		////���ω�����
		//if( mPhase[i] == -1 && mHeats[i] > mLatentHeat ){		//�ܗL�M�ʂ��Z����M������
		//	mPhase[i] = 2;										//���ւ̑��ω�
		//	mTemps[i] = 251.0f;
		//	mHeats[i] = 0;										//�M�i�|�e���V�����G�l���M�[�j���z��
		//	mPhaseChange[i] = 1;
		//}
		//else if( mPhase[i] == 1 && mHeats[i] < 0.0f ){			//�ܗL�M�ʂ��ÌŐ��M���g���؂�
		//	mPhase[i] = -2;										//�X�ւ̑��ω�
		//	mTemps[i] = 249.0f;
		//	mHeats[i] = 0;										//�M�i�|�e���V�����G�l���M�[�j����o
		//	mPhaseChange[i] = 1;
		//}
	}
	else{	
		//���M�v�Z
		//�ω��M�ʂ̓K�p
		float spcfHt = (temps[indx] > 250.0f)? 4.2f : 2.1f;		//��M�@���ƕX�Ŕ�M��ϓ��@��4.2�@�X2.1
		temps[indx] =	temps[indx] + heats[indx] / (spcfHt * 1.0f);	//�M�ʂ����x�Ɋ��Z�@���ʂ͌Œ肵�Ă���
		heats[indx] = 0.0f;										//�������@�t���[�����ƂɔM�ʁi���M�j�͒~�ς���Ȃ�

		//�ω����x��K�p
		temps[indx] += dtemps[indx];
		//if( temps[indx] > mTempMax) temps[indx] = mTempMax;
		//if( temps[indx] < mTempMin) temps[indx] = mTempMin;
	}
}

#endif