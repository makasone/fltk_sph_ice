#ifndef _GPU_HEATTRANSFER_CU_
#define _GPU_HEATTRANSFER_CU_

#include <math.h>
#include <stdio.h>
#include <iostream>
#include <cuda_runtime.h>

using namespace std;

#include "rx_cu_common.cuh"	//æ¶‚ª’è‹`‚µ‚½•Ö—˜‹@”\‚ªg‚¦‚é

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

	//•\–Ê—±q‚Æ‹ó‹CC—±q‚Æ—±qC”MŒ¹‚Æ—±q
	//UpdateHeatTransfar<<<numBlocks, numThreads>>>(heats, temps, surfParticles, airTemp, ht);
	UpdateAirAndParticle<<<numBlocks, numThreads>>>(heats, temps, surfParticles, airTemp, cffCnHt);
	UpdatePhase<<<numBlocks, numThreads>>>(heats, temps, dtemps);

	cudaThreadSynchronize();
}

//”Mˆ—
__global__
void UpdateAirAndParticle(float* heats, const float* temps, const int* surParticles, float airTemp, float ht)
{
	int indx = __umul24(blockIdx.x,blockDim.x)+threadIdx.x;

//•\–Ê—±q‚Æ‹ó‹C
	if( surParticles[indx] == -1 ) return;

	double airNum = 20.0 - (double)surParticles[indx];
	if(airNum < 0) airNum = 0.0;

	double surfaceArea = airNum / 20.0;								//‹ó‹C‚ÆG‚ê‚Ä‚¢‚é•\–ÊÏ@0`1.0 20‚Í“K“–

	if( surfaceArea < 0.0) surfaceArea = 0.0;
	if( surfaceArea > 1.0) surfaceArea = 1.0;

	double qHeat = ht * (airTemp - temps[indx]) * surfaceArea;		//ƒjƒ…[ƒgƒ“‚Ì—â‹p–@‘¥‚Ì®‚©‚ç”M—Ê‚ğŒvZ
	heats[indx] += qHeat;											//”M—Ê‚ğ‰ÁZ
}

//‰·“x‚Æ”M—Ê‚Ìˆ—@Œ°”MEö”M‚ÌŒvZC‘Š•Ï‰»”»’è
__global__
void UpdatePhase(float* heats, float* temps, float* dtemps)
{
	int indx = __umul24(blockIdx.x,blockDim.x)+threadIdx.x;

	////’†ŠÔó‘Ô‚Ö‚Ì•Ï‰»ŒŸo
	//if( mPhase[i] == -2 && mTemps[i] > 250.0f ){				//•X‚Ìê‡
	//	mPhase[i] = -1;											//•X’†ŠÔó‘Ô
	//	mTemps[i] = 250.0f;
	//	mHeats[i] = 0;
	//}
	//else if( mPhase[i] == 2 && mTemps[i] < 250.0f ){			//…‚Ìê‡
	//	mPhase[i] = 1;											//…’†ŠÔó‘Ô
	//	mTemps[i] = 250.0f;
	//	mHeats[i] = mLatentHeat;								//—Z‰ğö”M
	//}

	//Œ°”MEö”M‚ÌŒvZ
	if(false/*mPhase[i] == -1 || mPhase[i] == 1*/){					//…’†ŠÔ‚©•X’†ŠÔ‚Ìê‡
		////ö”MŒvZ
		////•Ï‰»‰·“x‚ğ”M—Ê‚É•ÏŠ·‚µ‚Äö”MŒvZ@i‰·“x‚Í•Ï‰»‚µ‚È‚¢j
		//float spcfHt = 2.1f + (2.1f * mHeats[i] / mLatentHeat);	//”ä”M@ŠÜ—L”M—Ê‚Å…‚Æ•X‚Ì”ä”M‚ğ•âŠÔ
		//mHeats[i]	+= mTempsDelta[i] * spcfHt * 1.0f;			//‰·“x•Ï‰»‚ğ”M—Ê‚ÉŠ·Z@¿—Ê‚Í‚P‚ÅŒÅ’è

		////ö”M•Ï‰»‚©‚çŒ°”M•Ï‰»‚Ö–ß‚é”»’è
		//if(mPhase[i] == -1 && mHeats[i] < -mLatentHeat/300.0f){	//‰ü‘P@<0@‚É‚·‚é‚ÆC—n‚¯‚È‚­‚È‚é
		//	//•X’†ŠÔó‘Ô¨•X
		//	mPhase[i] = -2;										//•X‚Ö‚Ì‘Š•Ï‰»
		//	mTemps[i] = 249.0f;
		//	mHeats[i] = 0;										//”Miƒ|ƒeƒ“ƒVƒƒƒ‹ƒGƒlƒ‹ƒM[j‚ğ•úo
		//}
		//else if( mPhase[i] == 1 && mHeats[i] > mLatentHeat ){
		//	//…’†ŠÔó‘Ô¨…
		//	mPhase[i] = 2;										//…‚Ö‚Ì‘Š•Ï‰»
		//	mTemps[i] = 251.0f;
		//	mHeats[i] = 0;										//”Miƒ|ƒeƒ“ƒVƒƒƒ‹ƒGƒlƒ‹ƒM[j‚ğ‹zû
		//}

		////‘Š•Ï‰»”»’è
		//if( mPhase[i] == -1 && mHeats[i] > mLatentHeat ){		//ŠÜ—L”M—Ê‚ª—Z‰ğö”M‚ğã‰ñ‚é
		//	mPhase[i] = 2;										//…‚Ö‚Ì‘Š•Ï‰»
		//	mTemps[i] = 251.0f;
		//	mHeats[i] = 0;										//”Miƒ|ƒeƒ“ƒVƒƒƒ‹ƒGƒlƒ‹ƒM[j‚ğ‹zû
		//	mPhaseChange[i] = 1;
		//}
		//else if( mPhase[i] == 1 && mHeats[i] < 0.0f ){			//ŠÜ—L”M—Ê‚ª‹ÃŒÅö”M‚ğg‚¢Ø‚é
		//	mPhase[i] = -2;										//•X‚Ö‚Ì‘Š•Ï‰»
		//	mTemps[i] = 249.0f;
		//	mHeats[i] = 0;										//”Miƒ|ƒeƒ“ƒVƒƒƒ‹ƒGƒlƒ‹ƒM[j‚ğ•úo
		//	mPhaseChange[i] = 1;
		//}
	}
	else{	
		//Œ°”MŒvZ
		//•Ï‰»”M—Ê‚Ì“K—p
		float spcfHt = (temps[indx] > 250.0f)? 4.2f : 2.1f;		//”ä”M@…‚Æ•X‚Å”ä”M‚ğ•Ï“®@…4.2@•X2.1
		temps[indx] =	temps[indx] + heats[indx] / (spcfHt * 1.0f);	//”M—Ê‚ğ‰·“x‚ÉŠ·Z@¿—Ê‚ÍŒÅ’è‚µ‚Ä‚¢‚é
		heats[indx] = 0.0f;										//‰Šú‰»@ƒtƒŒ[ƒ€‚²‚Æ‚É”M—ÊiŒ°”Mj‚Í’~Ï‚³‚ê‚È‚¢

		//•Ï‰»‰·“x‚ğ“K—p
		temps[indx] += dtemps[indx];
		//if( temps[indx] > mTempMax) temps[indx] = mTempMax;
		//if( temps[indx] < mTempMin) temps[indx] = mTempMin;
	}
}

#endif