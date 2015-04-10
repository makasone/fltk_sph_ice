#ifndef _GPU_HEATTRANSFER_CU_
#define _GPU_HEATTRANSFER_CU_

#include <math.h>
#include <stdio.h>
#include <iostream>
#include <cuda_runtime.h>

using namespace std;

#include "rx_cu_common.cuh"	//先生が定義した便利機能が使える

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

	//表面粒子と空気，粒子と粒子，熱源と粒子
	//UpdateHeatTransfar<<<numBlocks, numThreads>>>(heats, temps, surfParticles, airTemp, ht);
	UpdateAirAndParticle<<<numBlocks, numThreads>>>(heats, temps, surfParticles, airTemp, cffCnHt);
	UpdatePhase<<<numBlocks, numThreads>>>(heats, temps, dtemps);

	cudaThreadSynchronize();
}

//熱処理
__global__
void UpdateAirAndParticle(float* heats, const float* temps, const int* surParticles, float airTemp, float ht)
{
	int indx = __umul24(blockIdx.x,blockDim.x)+threadIdx.x;

//表面粒子と空気
	if( surParticles[indx] == -1 ) return;

	double airNum = 20.0 - (double)surParticles[indx];
	if(airNum < 0) airNum = 0.0;

	double surfaceArea = airNum / 20.0;								//空気と触れている表面積　0〜1.0 20は適当

	if( surfaceArea < 0.0) surfaceArea = 0.0;
	if( surfaceArea > 1.0) surfaceArea = 1.0;

	double qHeat = ht * (airTemp - temps[indx]) * surfaceArea;		//ニュートンの冷却法則の式から熱量を計算
	heats[indx] += qHeat;											//熱量を加算
}

//温度と熱量の処理　顕熱・潜熱の計算，相変化判定
__global__
void UpdatePhase(float* heats, float* temps, float* dtemps)
{
	int indx = __umul24(blockIdx.x,blockDim.x)+threadIdx.x;

	////中間状態への変化検出
	//if( mPhase[i] == -2 && mTemps[i] > 250.0f ){				//氷の場合
	//	mPhase[i] = -1;											//氷中間状態
	//	mTemps[i] = 250.0f;
	//	mHeats[i] = 0;
	//}
	//else if( mPhase[i] == 2 && mTemps[i] < 250.0f ){			//水の場合
	//	mPhase[i] = 1;											//水中間状態
	//	mTemps[i] = 250.0f;
	//	mHeats[i] = mLatentHeat;								//融解潜熱
	//}

	//顕熱・潜熱の計算
	if(false/*mPhase[i] == -1 || mPhase[i] == 1*/){					//水中間か氷中間の場合
		////潜熱計算
		////変化温度を熱量に変換して潜熱計算　（温度は変化しない）
		//float spcfHt = 2.1f + (2.1f * mHeats[i] / mLatentHeat);	//比熱　含有熱量で水と氷の比熱を補間
		//mHeats[i]	+= mTempsDelta[i] * spcfHt * 1.0f;			//温度変化を熱量に換算　質量は１で固定

		////潜熱変化から顕熱変化へ戻る判定
		//if(mPhase[i] == -1 && mHeats[i] < -mLatentHeat/300.0f){	//改善　<0　にすると，溶けなくなる
		//	//氷中間状態→氷
		//	mPhase[i] = -2;										//氷への相変化
		//	mTemps[i] = 249.0f;
		//	mHeats[i] = 0;										//熱（ポテンシャルエネルギー）を放出
		//}
		//else if( mPhase[i] == 1 && mHeats[i] > mLatentHeat ){
		//	//水中間状態→水
		//	mPhase[i] = 2;										//水への相変化
		//	mTemps[i] = 251.0f;
		//	mHeats[i] = 0;										//熱（ポテンシャルエネルギー）を吸収
		//}

		////相変化判定
		//if( mPhase[i] == -1 && mHeats[i] > mLatentHeat ){		//含有熱量が融解潜熱を上回る
		//	mPhase[i] = 2;										//水への相変化
		//	mTemps[i] = 251.0f;
		//	mHeats[i] = 0;										//熱（ポテンシャルエネルギー）を吸収
		//	mPhaseChange[i] = 1;
		//}
		//else if( mPhase[i] == 1 && mHeats[i] < 0.0f ){			//含有熱量が凝固潜熱を使い切る
		//	mPhase[i] = -2;										//氷への相変化
		//	mTemps[i] = 249.0f;
		//	mHeats[i] = 0;										//熱（ポテンシャルエネルギー）を放出
		//	mPhaseChange[i] = 1;
		//}
	}
	else{	
		//顕熱計算
		//変化熱量の適用
		float spcfHt = (temps[indx] > 250.0f)? 4.2f : 2.1f;		//比熱　水と氷で比熱を変動　水4.2　氷2.1
		temps[indx] =	temps[indx] + heats[indx] / (spcfHt * 1.0f);	//熱量を温度に換算　質量は固定している
		heats[indx] = 0.0f;										//初期化　フレームごとに熱量（顕熱）は蓄積されない

		//変化温度を適用
		temps[indx] += dtemps[indx];
		//if( temps[indx] > mTempMax) temps[indx] = mTempMax;
		//if( temps[indx] < mTempMin) temps[indx] = mTempMin;
	}
}

#endif