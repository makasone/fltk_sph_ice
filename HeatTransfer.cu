#ifndef _GPU_HEATTRANSFER_CU_
#define _GPU_HEATTRANSFER_CU_

#include <math.h>
#include <stdio.h>
#include <iostream>
#include <cuda_runtime.h>

using namespace std;

#include "rx_cu_common.cuh"	//�搶����`�����֗��@�\���g����

#define TEMP_MAX (500.0f)
#define TEMP_MIN (0.0f)
#define NEIGHT_MAX 30
#define PI (3.14159265359f)
#define LATENT_HEAT (165.0f)

void LaunchHeatTransferGPU
(
	float* heats,
	float* temps,
	float* dtemps,
	const int* surfParticles,
	const int* neightsParticles,
	const float* densty,
	int* phase,
	int* phaseChangef,
	const int* meltPrtIndx,
	float radius, 
	float airTemp,
	float cffCnHt,
	float cffCnTd,
	float dt,
	int prtNum
);

__global__ void UpdateHeatTransfar(float* heats, const float* temps, const int* surfParticles, float airTemp, float ht);
__global__ void UpdateAirAndParticle(float* heats, const float* temps, const int* surfParticles, float airTemp, float ht);
__global__ void UpdateParticleAndParticle(float dt, const int* neighs, const float* densty, float* temps, float* dtemps, float radius, float cffCnTd, float LWspline);
__global__ void UpdateMeltParticle(const int* meltPrtIndx, float* heats);
__global__ void UpdatePhase(float* heats, float* temps, float* dtemps, int* phase, int* phaseChangef);

void LaunchHeatTransferGPU(
	float* heats, float* temps, float* dtemps, const int* surfParticles, const int* neightsPrt, const float* densty, int* phase, int* phaseChangef,
	const int* meltPrtIndx, float radius, float airTemp, float cffCnHt, float cffCnTd, float dt, int prtNum)
{
	uint numThreads = min(THREAD_NUM, prtNum);
	uint numBlocks = (prtNum % numThreads != 0) ? (prtNum / numThreads + 1) : (prtNum / numThreads);

	float Wspline = 10.0f/(7.0f*PI*(float)pow((double)radius, (double)2.0));

	//�\�ʗ��q�Ƌ�C�C���q�Ɨ��q�C�M���Ɨ��q
	//UpdateHeatTransfar<<<numBlocks, numThreads>>>(heats, temps, surfParticles, airTemp, ht);
	UpdateAirAndParticle<<<numBlocks, numThreads>>>(heats, temps, surfParticles, airTemp, cffCnHt);
	UpdateParticleAndParticle<<<numBlocks, numThreads>>>(dt, neightsPrt, densty, temps, dtemps, radius, cffCnTd, Wspline);
	UpdateMeltParticle<<<numBlocks, numThreads>>>(meltPrtIndx, heats);
	UpdatePhase<<<numBlocks, numThreads>>>(heats, temps, dtemps, phase, phaseChangef);

	cudaThreadSynchronize();
}

__device__
float KernelSpline(const float d, const float r, const float Wspline)
{
	float q = d/r;
    if(q >= 0.0f && q < 1.0f){
        return Wspline*(1.0f-1.5f*q*q+0.75f*q*q*q);
    }
    else if(q >= 1.0f && q < 2.0f){
        return Wspline*0.25f*(2.0f-q)*(2.0f-q)*(2.0f-q);
    }
    else{
        return 0.0f;
    }
}

//�M����
__global__
void UpdateAirAndParticle(float* heats, const float* temps, const int* surParticles, float airTemp, float ht)
{
	int indx = __umul24(blockIdx.x,blockDim.x)+threadIdx.x;

//�\�ʗ��q�Ƌ�C
	if( surParticles[indx] <= 0) return;

	float airNum = 20.0f - (float)surParticles[indx];
	if(airNum < 0) airNum = 0.0f;

	float surfaceArea = airNum / 20.0f;								//��C�ƐG��Ă���\�ʐρ@0�`1.0 20�͓K��

	if( surfaceArea < 0.0f) surfaceArea = 0.0f;
	if( surfaceArea > 1.0f) surfaceArea = 1.0f;

	float qHeat = ht * (airTemp - temps[indx]) * surfaceArea;		//�j���[�g���̗�p�@���̎�����M�ʂ��v�Z
	heats[indx] += qHeat;											//�M�ʂ����Z
}

//���q�Ɨ��q
__global__
void UpdateParticleAndParticle(float dt, const int* neighs, const float* densty, float* temps, float* dtemps, float radius, float cffCnTd, float Wspline)
{
	int indx = __umul24(blockIdx.x,blockDim.x)+threadIdx.x;

	float tmp = 0.0f;
	float nowTemp = temps[indx];

	//�ߖT���q�̐������܂킷
	int nNum = neighs[indx * NEIGHT_MAX + NEIGHT_MAX - 1];
	for(unsigned j = 0; j < nNum; j++){
		int id = neighs[indx * NEIGHT_MAX + j];
		if(indx == id )	continue;							//�ߖT���q�Ɏ������܂܂�邱�Ƃ����邽��

		float d = densty[id];
		if( d < 0.05f ) d = 0.05f;							//���x������������or�O�̏ꍇ������̂Œ���

		float dis = /*mNeighborhoodsDis[i][j]*/1.0f;

		tmp += dt * 1.0f * (temps[id]-nowTemp) / d * KernelSpline(dis, radius, Wspline);
	}

	dtemps[indx] = tmp * cffCnTd;		//�M�g�U�W��
}

//����̗��q��Z��
//����Ȃ̂Ŗ��ʂ��炯
__global__
void UpdateMeltParticle(const int* meltPrtIndx, float* heats)
{
	int indx = __umul24(blockIdx.x,blockDim.x)+threadIdx.x;
	if(1000 <= indx || indx == 0)	return;		//1000�͓K��
	if(meltPrtIndx[indx] == -1)		return;

	int meltIndx = meltPrtIndx[indx];

	heats[meltIndx] += 300.0f;
}

//���x�ƔM�ʂ̏����@���M�E���M�̌v�Z�C���ω�����
__global__
void UpdatePhase(float* heats, float* temps, float* dtemps, int* phase, int* phaseChangef)
{
	int indx = __umul24(blockIdx.x,blockDim.x)+threadIdx.x;

	//���ԏ�Ԃւ̕ω����o
	if(phase[indx] == -2 && temps[indx] > 250.0f){				//�X�̏ꍇ
		phase[indx] = -1;										//�X���ԏ��
		temps[indx] = 250.0f;
		heats[indx] = 0;
	}
	else if(phase[indx] == 2 && temps[indx] < 250.0f){			//���̏ꍇ
		phase[indx] = 1;										//�����ԏ��
		temps[indx] = 250.0f;
		heats[indx] = LATENT_HEAT;								//�Z����M
	}

	//���M�E���M�̌v�Z
	if(phase[indx] == -1 || phase[indx] == 1){						//�����Ԃ��X���Ԃ̏ꍇ
		//���M�v�Z
		//�ω����x��M�ʂɕϊ����Đ��M�v�Z�@�i���x�͕ω����Ȃ��j
		float spcfHt = 2.1f + (2.1f * heats[indx] / LATENT_HEAT);	//��M�@�ܗL�M�ʂŐ��ƕX�̔�M����
		heats[indx] += dtemps[indx] * spcfHt * 1.0f;				//���x�ω���M�ʂɊ��Z�@���ʂ͂P�ŌŒ�

		//���M�ω����猰�M�ω��֖߂锻��
		if(phase[indx] == -1 && heats[indx] < -LATENT_HEAT/300.0f){	//���P�@<0�@�ɂ���ƁC�n���Ȃ��Ȃ�
			//�X���ԏ�ԁ��X
			phase[indx] = -2;										//�X�ւ̑��ω�
			temps[indx] = 249.0f;
			heats[indx] = 0;										//�M�i�|�e���V�����G�l���M�[�j����o
		}
		else if(phase[indx] == 1 && heats[indx] > LATENT_HEAT){
			//�����ԏ�ԁ���
			phase[indx] = 2;										//���ւ̑��ω�
			temps[indx] = 251.0f;
			heats[indx] = 0;										//�M�i�|�e���V�����G�l���M�[�j���z��
		}

		//���ω�����
		if(phase[indx] == -1 && heats[indx] > LATENT_HEAT){			//�ܗL�M�ʂ��Z����M������
			phase[indx] = 2;										//���ւ̑��ω�
			temps[indx] = 251.0f;
			heats[indx] = 0;										//�M�i�|�e���V�����G�l���M�[�j���z��
			phaseChangef[indx] = 1;
		}
		else if(phase[indx] == 1 && heats[indx] < 0.0f){			//�ܗL�M�ʂ��ÌŐ��M���g���؂�
			phase[indx] = -2;										//�X�ւ̑��ω�
			temps[indx] = 249.0f;
			heats[indx] = 0;										//�M�i�|�e���V�����G�l���M�[�j����o
			phaseChangef[indx] = 1;
		}
	}
	else{	
		//���M�v�Z
		//�ω��M�ʂ̓K�p
		float spcfHt = (temps[indx] > 250.0f)? 4.2f : 2.1f;			//��M�@���ƕX�Ŕ�M��ϓ��@��4.2�@�X2.1
		temps[indx] = temps[indx] + heats[indx] / (spcfHt * 1.0f);	//�M�ʂ����x�Ɋ��Z�@���ʂ͌Œ肵�Ă���
		heats[indx] = 0.0f;											//�������@�t���[�����ƂɔM�ʁi���M�j�͒~�ς���Ȃ�

		//�ω����x��K�p
		temps[indx] += dtemps[indx];
		if(temps[indx] > TEMP_MAX) temps[indx] = TEMP_MAX;
		if(temps[indx] < TEMP_MIN) temps[indx] = TEMP_MIN;
	}
}

#endif