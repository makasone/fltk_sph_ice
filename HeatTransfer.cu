#ifndef _GPU_HEATTRANSFER_CU_
#define _GPU_HEATTRANSFER_CU_

#include <math.h>
#include <stdio.h>
#include <iostream>
#include <cuda_runtime.h>

using namespace std;

#include "rx_cu_common.cuh"	//先生が定義した便利機能が使える

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

	//表面粒子と空気，粒子と粒子，熱源と粒子
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

//熱処理
__global__
void UpdateAirAndParticle(float* heats, const float* temps, const int* surParticles, float airTemp, float ht)
{
	int indx = __umul24(blockIdx.x,blockDim.x)+threadIdx.x;

//表面粒子と空気
	if( surParticles[indx] <= 0) return;

	float airNum = 20.0f - (float)surParticles[indx];
	if(airNum < 0) airNum = 0.0f;

	float surfaceArea = airNum / 20.0f;								//空気と触れている表面積　0〜1.0 20は適当

	if( surfaceArea < 0.0f) surfaceArea = 0.0f;
	if( surfaceArea > 1.0f) surfaceArea = 1.0f;

	float qHeat = ht * (airTemp - temps[indx]) * surfaceArea;		//ニュートンの冷却法則の式から熱量を計算
	heats[indx] += qHeat;											//熱量を加算
}

//粒子と粒子
__global__
void UpdateParticleAndParticle(float dt, const int* neighs, const float* densty, float* temps, float* dtemps, float radius, float cffCnTd, float Wspline)
{
	int indx = __umul24(blockIdx.x,blockDim.x)+threadIdx.x;

	float tmp = 0.0f;
	float nowTemp = temps[indx];

	//近傍粒子の数だけまわす
	int nNum = neighs[indx * NEIGHT_MAX + NEIGHT_MAX - 1];
	for(unsigned j = 0; j < nNum; j++){
		int id = neighs[indx * NEIGHT_MAX + j];
		if(indx == id )	continue;							//近傍粒子に自分が含まれることがあるため

		float d = densty[id];
		if( d < 0.05f ) d = 0.05f;							//密度が小さすぎるor０の場合があるので調整

		float dis = /*mNeighborhoodsDis[i][j]*/1.0f;

		tmp += dt * 1.0f * (temps[id]-nowTemp) / d * KernelSpline(dis, radius, Wspline);
	}

	dtemps[indx] = tmp * cffCnTd;		//熱拡散係数
}

//特定の粒子を融解
//試作なので無駄だらけ
__global__
void UpdateMeltParticle(const int* meltPrtIndx, float* heats)
{
	int indx = __umul24(blockIdx.x,blockDim.x)+threadIdx.x;
	if(1000 <= indx || indx == 0)	return;		//1000は適当
	if(meltPrtIndx[indx] == -1)		return;

	int meltIndx = meltPrtIndx[indx];

	heats[meltIndx] += 300.0f;
}

//温度と熱量の処理　顕熱・潜熱の計算，相変化判定
__global__
void UpdatePhase(float* heats, float* temps, float* dtemps, int* phase, int* phaseChangef)
{
	int indx = __umul24(blockIdx.x,blockDim.x)+threadIdx.x;

	//中間状態への変化検出
	if(phase[indx] == -2 && temps[indx] > 250.0f){				//氷の場合
		phase[indx] = -1;										//氷中間状態
		temps[indx] = 250.0f;
		heats[indx] = 0;
	}
	else if(phase[indx] == 2 && temps[indx] < 250.0f){			//水の場合
		phase[indx] = 1;										//水中間状態
		temps[indx] = 250.0f;
		heats[indx] = LATENT_HEAT;								//融解潜熱
	}

	//顕熱・潜熱の計算
	if(phase[indx] == -1 || phase[indx] == 1){						//水中間か氷中間の場合
		//潜熱計算
		//変化温度を熱量に変換して潜熱計算　（温度は変化しない）
		float spcfHt = 2.1f + (2.1f * heats[indx] / LATENT_HEAT);	//比熱　含有熱量で水と氷の比熱を補間
		heats[indx] += dtemps[indx] * spcfHt * 1.0f;				//温度変化を熱量に換算　質量は１で固定

		//潜熱変化から顕熱変化へ戻る判定
		if(phase[indx] == -1 && heats[indx] < -LATENT_HEAT/300.0f){	//改善　<0　にすると，溶けなくなる
			//氷中間状態→氷
			phase[indx] = -2;										//氷への相変化
			temps[indx] = 249.0f;
			heats[indx] = 0;										//熱（ポテンシャルエネルギー）を放出
		}
		else if(phase[indx] == 1 && heats[indx] > LATENT_HEAT){
			//水中間状態→水
			phase[indx] = 2;										//水への相変化
			temps[indx] = 251.0f;
			heats[indx] = 0;										//熱（ポテンシャルエネルギー）を吸収
		}

		//相変化判定
		if(phase[indx] == -1 && heats[indx] > LATENT_HEAT){			//含有熱量が融解潜熱を上回る
			phase[indx] = 2;										//水への相変化
			temps[indx] = 251.0f;
			heats[indx] = 0;										//熱（ポテンシャルエネルギー）を吸収
			phaseChangef[indx] = 1;
		}
		else if(phase[indx] == 1 && heats[indx] < 0.0f){			//含有熱量が凝固潜熱を使い切る
			phase[indx] = -2;										//氷への相変化
			temps[indx] = 249.0f;
			heats[indx] = 0;										//熱（ポテンシャルエネルギー）を放出
			phaseChangef[indx] = 1;
		}
	}
	else{	
		//顕熱計算
		//変化熱量の適用
		float spcfHt = (temps[indx] > 250.0f)? 4.2f : 2.1f;			//比熱　水と氷で比熱を変動　水4.2　氷2.1
		temps[indx] = temps[indx] + heats[indx] / (spcfHt * 1.0f);	//熱量を温度に換算　質量は固定している
		heats[indx] = 0.0f;											//初期化　フレームごとに熱量（顕熱）は蓄積されない

		//変化温度を適用
		temps[indx] += dtemps[indx];
		if(temps[indx] > TEMP_MAX) temps[indx] = TEMP_MAX;
		if(temps[indx] < TEMP_MIN) temps[indx] = TEMP_MIN;
	}
}

#endif