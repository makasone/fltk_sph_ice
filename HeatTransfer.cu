#ifndef _GPU_HEATTRANSFER_CU_
#define _GPU_HEATTRANSFER_CU_

#include <math.h>
#include <stdio.h>
#include <iostream>
#include <cuda_runtime.h>

using namespace std;

#include "rx_cu_common.cuh"	//æ¶‚ª’è‹`‚µ‚½•Ö—˜‹@”\‚ªg‚¦‚é

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

	//•\–Ê—±q‚Æ‹ó‹CC—±q‚Æ—±qC”MŒ¹‚Æ—±q
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

//”Mˆ—
__global__
void UpdateAirAndParticle(float* heats, const float* temps, const int* surParticles, float airTemp, float ht)
{
	int indx = __umul24(blockIdx.x,blockDim.x)+threadIdx.x;

//•\–Ê—±q‚Æ‹ó‹C
	if( surParticles[indx] <= 0) return;

	float airNum = 20.0f - (float)surParticles[indx];
	if(airNum < 0) airNum = 0.0f;

	float surfaceArea = airNum / 20.0f;								//‹ó‹C‚ÆG‚ê‚Ä‚¢‚é•\–ÊÏ@0`1.0 20‚Í“K“–

	if( surfaceArea < 0.0f) surfaceArea = 0.0f;
	if( surfaceArea > 1.0f) surfaceArea = 1.0f;

	float qHeat = ht * (airTemp - temps[indx]) * surfaceArea;		//ƒjƒ…[ƒgƒ“‚Ì—â‹p–@‘¥‚Ì®‚©‚ç”M—Ê‚ğŒvZ
	heats[indx] += qHeat;											//”M—Ê‚ğ‰ÁZ
}

//—±q‚Æ—±q
__global__
void UpdateParticleAndParticle(float dt, const int* neighs, const float* densty, float* temps, float* dtemps, float radius, float cffCnTd, float Wspline)
{
	int indx = __umul24(blockIdx.x,blockDim.x)+threadIdx.x;

	float tmp = 0.0f;
	float nowTemp = temps[indx];

	//‹ß–T—±q‚Ì”‚¾‚¯‚Ü‚í‚·
	int nNum = neighs[indx * NEIGHT_MAX + NEIGHT_MAX - 1];
	for(unsigned j = 0; j < nNum; j++){
		int id = neighs[indx * NEIGHT_MAX + j];
		if(indx == id )	continue;							//‹ß–T—±q‚É©•ª‚ªŠÜ‚Ü‚ê‚é‚±‚Æ‚ª‚ ‚é‚½‚ß

		float d = densty[id];
		if( d < 0.05f ) d = 0.05f;							//–§“x‚ª¬‚³‚·‚¬‚éor‚O‚Ìê‡‚ª‚ ‚é‚Ì‚Å’²®

		float dis = /*mNeighborhoodsDis[i][j]*/1.0f;

		tmp += dt * 1.0f * (temps[id]-nowTemp) / d * KernelSpline(dis, radius, Wspline);
	}

	dtemps[indx] = tmp * cffCnTd;		//”MŠgUŒW”
}

//“Á’è‚Ì—±q‚ğ—Z‰ğ
//ì‚È‚Ì‚Å–³‘Ê‚¾‚ç‚¯
__global__
void UpdateMeltParticle(const int* meltPrtIndx, float* heats)
{
	int indx = __umul24(blockIdx.x,blockDim.x)+threadIdx.x;
	if(1000 <= indx || indx == 0)	return;		//1000‚Í“K“–
	if(meltPrtIndx[indx] == -1)		return;

	int meltIndx = meltPrtIndx[indx];

	heats[meltIndx] += 300.0f;
}

//‰·“x‚Æ”M—Ê‚Ìˆ—@Œ°”MEö”M‚ÌŒvZC‘Š•Ï‰»”»’è
__global__
void UpdatePhase(float* heats, float* temps, float* dtemps, int* phase, int* phaseChangef)
{
	int indx = __umul24(blockIdx.x,blockDim.x)+threadIdx.x;

	//’†ŠÔó‘Ô‚Ö‚Ì•Ï‰»ŒŸo
	if(phase[indx] == -2 && temps[indx] > 250.0f){				//•X‚Ìê‡
		phase[indx] = -1;										//•X’†ŠÔó‘Ô
		temps[indx] = 250.0f;
		heats[indx] = 0;
	}
	else if(phase[indx] == 2 && temps[indx] < 250.0f){			//…‚Ìê‡
		phase[indx] = 1;										//…’†ŠÔó‘Ô
		temps[indx] = 250.0f;
		heats[indx] = LATENT_HEAT;								//—Z‰ğö”M
	}

	//Œ°”MEö”M‚ÌŒvZ
	if(phase[indx] == -1 || phase[indx] == 1){						//…’†ŠÔ‚©•X’†ŠÔ‚Ìê‡
		//ö”MŒvZ
		//•Ï‰»‰·“x‚ğ”M—Ê‚É•ÏŠ·‚µ‚Äö”MŒvZ@i‰·“x‚Í•Ï‰»‚µ‚È‚¢j
		float spcfHt = 2.1f + (2.1f * heats[indx] / LATENT_HEAT);	//”ä”M@ŠÜ—L”M—Ê‚Å…‚Æ•X‚Ì”ä”M‚ğ•âŠÔ
		heats[indx] += dtemps[indx] * spcfHt * 1.0f;				//‰·“x•Ï‰»‚ğ”M—Ê‚ÉŠ·Z@¿—Ê‚Í‚P‚ÅŒÅ’è

		//ö”M•Ï‰»‚©‚çŒ°”M•Ï‰»‚Ö–ß‚é”»’è
		if(phase[indx] == -1 && heats[indx] < -LATENT_HEAT/300.0f){	//‰ü‘P@<0@‚É‚·‚é‚ÆC—n‚¯‚È‚­‚È‚é
			//•X’†ŠÔó‘Ô¨•X
			phase[indx] = -2;										//•X‚Ö‚Ì‘Š•Ï‰»
			temps[indx] = 249.0f;
			heats[indx] = 0;										//”Miƒ|ƒeƒ“ƒVƒƒƒ‹ƒGƒlƒ‹ƒM[j‚ğ•úo
		}
		else if(phase[indx] == 1 && heats[indx] > LATENT_HEAT){
			//…’†ŠÔó‘Ô¨…
			phase[indx] = 2;										//…‚Ö‚Ì‘Š•Ï‰»
			temps[indx] = 251.0f;
			heats[indx] = 0;										//”Miƒ|ƒeƒ“ƒVƒƒƒ‹ƒGƒlƒ‹ƒM[j‚ğ‹zû
		}

		//‘Š•Ï‰»”»’è
		if(phase[indx] == -1 && heats[indx] > LATENT_HEAT){			//ŠÜ—L”M—Ê‚ª—Z‰ğö”M‚ğã‰ñ‚é
			phase[indx] = 2;										//…‚Ö‚Ì‘Š•Ï‰»
			temps[indx] = 251.0f;
			heats[indx] = 0;										//”Miƒ|ƒeƒ“ƒVƒƒƒ‹ƒGƒlƒ‹ƒM[j‚ğ‹zû
			phaseChangef[indx] = 1;
		}
		else if(phase[indx] == 1 && heats[indx] < 0.0f){			//ŠÜ—L”M—Ê‚ª‹ÃŒÅö”M‚ğg‚¢Ø‚é
			phase[indx] = -2;										//•X‚Ö‚Ì‘Š•Ï‰»
			temps[indx] = 249.0f;
			heats[indx] = 0;										//”Miƒ|ƒeƒ“ƒVƒƒƒ‹ƒGƒlƒ‹ƒM[j‚ğ•úo
			phaseChangef[indx] = 1;
		}
	}
	else{	
		//Œ°”MŒvZ
		//•Ï‰»”M—Ê‚Ì“K—p
		float spcfHt = (temps[indx] > 250.0f)? 4.2f : 2.1f;			//”ä”M@…‚Æ•X‚Å”ä”M‚ğ•Ï“®@…4.2@•X2.1
		temps[indx] = temps[indx] + heats[indx] / (spcfHt * 1.0f);	//”M—Ê‚ğ‰·“x‚ÉŠ·Z@¿—Ê‚ÍŒÅ’è‚µ‚Ä‚¢‚é
		heats[indx] = 0.0f;											//‰Šú‰»@ƒtƒŒ[ƒ€‚²‚Æ‚É”M—ÊiŒ°”Mj‚Í’~Ï‚³‚ê‚È‚¢

		//•Ï‰»‰·“x‚ğ“K—p
		temps[indx] += dtemps[indx];
		if(temps[indx] > TEMP_MAX) temps[indx] = TEMP_MAX;
		if(temps[indx] < TEMP_MIN) temps[indx] = TEMP_MIN;
	}
}

#endif