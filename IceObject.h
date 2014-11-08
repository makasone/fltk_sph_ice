//全相変化オブジェクトを管理するクラス

#ifndef _ICE_OBJECT_
#define _ICE_OBJECT_

#include "rx_sph.h"
#include "rx_utility.h"
#include "rx_matrix.h"

#include "Ice_SM.h"
#include "IceStructure.h"
#include "IceTetrahedra.h"
#include "HeatTransfar.h"
#include "Surf_SM.h"
#include "QueryCounter.h"

#include <time.h>
#include <cuda_runtime.h>
#include <cuda_gl_interop.h>

using namespace std;

//#define MK_USE_GPU
//#define USE_PATH
//#define USE_ITR
//#define USE_SELECTED

#define SELECTED 2

const int g_iterationNum = 1;

class IceObject
{
private:
	//TODO::液体運動計算クラスも置きたい
	//TODO::四面体処理も置く
	static float* s_sphPrtPos;						//sph粒子位置のホストポインタ
	static float* s_sphPrtVel;						//sph粒子速度のホストポインタ

//--------------------------------------GPU__------------------------------------------------------------
	static float* sd_sphPrtPos;						//sph粒子位置のデバイスポインタ
	static float* sd_sphPrtVel;						//sph粒子速度のデバイスポインタ

	static float* sd_sldPrtPos;						//総和計算による最終的な粒子位置のデバイスポインタ
	static float* sd_sldPrtVel;						//総和計算による最終的な粒子速度のデバイスポインタ
//--------------------------------------__GPU------------------------------------------------------------

	static int sm_particleNum;						//現在の粒子数
	static int sm_tetraNum;							//現在の四面体数
	static int sm_clusterNum;						//現在のクラスタ数
	static int sm_layerNum;							//探索レイヤー数
	
	static int sm_maxParticleNum;					//最大粒子数

	//構造管理クラス
	IceStructure* m_iceStrct;

	//固体運動計算クラス
	vector<Ice_SM*> m_iceMove;

	//高速計算用クラス
	//TODO::ポインタにしたら？
	Surf_SM m_SurfSm;

	//熱処理
	HeatTransfar* m_heatTransfer;

	
	static float* m_fInterPolationCoefficience;		//線形補間係数

public:
	IceObject(int pMaxNum, int cMaxNum, int tMaxNum, int prtNum, float* hp, float* hv, float* dp, float* dv, int layer, int maxParticleNum);
	~IceObject();

	void InitIceObj(int pMaxNum, int cMaxNum, int tMaxNum);
	void InitIceObjGPU();
	void InitTetra();
	void InitCluster(Vec3 boundarySpaceHigh, Vec3 boundarySpaceLow, float timeStep, int itr);
	void InitStrct();
	void InitPath();
	void InitHeatTransfer(float effectiveRadius, float timeStep, float tempMax, float tempMin, float latentHeat, float cffcntHt, float cffcntTd);

	static void InitInterPolation();
	void InitGPU();

	void SetSPHDevicePointer(float* pos, float* vel){	sd_sphPrtPos = pos; sd_sphPrtVel = vel;	}
	void SetSPHHostPointer(float* pos, float* vel){		s_sphPrtPos = pos;	s_sphPrtVel = vel;	}
	void SetSearchLayerNum(int layer){					sm_layerNum = layer;				}
	void SetMaxParticleNum(int particleNum){			sm_maxParticleNum = particleNum;	}
	void SetAirTemp(float temp){						m_heatTransfer->setAirTemp(temp);	}

	void SetClusterMoveInfo(int pIndx);
	void SetClusterStrctInfo(int cIndx, int *PtoCNum);

	static int GetParticleNum(){	return sm_particleNum;		}
	static int GetClusterNum(){		return sm_clusterNum;		}
	static int GetTetrahedraNum(){	return sm_tetraNum;			}
	static int GetLayerNum(){		return sm_layerNum;			}

	static int GetMaxClusterNum(){	return sm_maxParticleNum;	}

	float* GetTemps(){	return m_heatTransfer->getTemps();}

	Ice_SM* GetMoveObj(int cIndx){	return m_iceMove[cIndx];	}

	void StepObjMoveCPU();									//運動計算
	void StepObjMoveGPU();									//GPUによる運動計算

	void StepObjMoveUsePath();								//高速化手法を用いた運動計算
	void StepObjMoveGPUUsePath();

	void StepObjMoveIteration();							//反復処理を用いた運動計算
	void StepObjMoveIterationUsePath();						//高速化手法を用いた反復運動計算
	void StepObjMoveIterationWeighted();					//重み付け＋反復処理

	void StepObjMoveSelected();

	void StepObjCalcWidhIteration();						//固体の運動計算，総和計算，補間処理　GPU処理　反復処理あり

	void StepInterPolation();								//線形補間　いずれは処理が複雑になるのでクラスにしたい．
	void StepInterPolationForCluster();
	void StepInterPolationForClusterWeighted();
	void StepInterPolationSelected();
	void StepWeightedInterPolation();
	
	void StepHeatTransfer(const int* surfParticles, const vector<vector<rxNeigh>>& neights, float floor, float effRadius, const float* pos, const float* dens);		//熱処理
	void StepIceStructure();								//相変化処理
	void StepMelting();										//融解
	void StepFreezing();									//凝固

	void SearchMeltParticle(vector<unsigned>& pList);
	void SearchFreezeParticle(vector<unsigned>& pList);

	void CalcAverageCPU(int pIndx, Vec3& pos, Vec3& vel);
	void CalcAverageSelected(int pIndx, Vec3& pos, Vec3& vel);
	void CalcWeightedVector(int pIndx, Vec3& pos, Vec3& vel);

	void LinerInterPolationCPU(int pIndx, const Vec3& pos, const Vec3& vel);
	void LinerInterPolationForClusterCPU(const int pIndx, const Vec3& pos, const Vec3& vel);

	void WarmParticle(int pIndx, float temp, float heat){	m_heatTransfer->WarmParticle(pIndx, temp, heat);	}
	void MeltParticle(int pIndx){	m_heatTransfer->MeltParticle(pIndx);	}

	void ReConstructCluster(vector<unsigned>& particleList, vector<unsigned>& clusterList);

	void UpdateUnSelectedCluster(int pIndx, const Vec3& pos, const Vec3& vel, const vector<Vec3>& prePos);




	//デバッグ
	void DebugTetraInfo();
	void DebugClusterInfo();
	void DebugObjMoveUsePathWithGPU();
	void DebugDeformationAmount();
	void DebugDeformationAverage();

	//テスト
	void TestStepInterPolation();
	void TestUpdateSMFromPath();







	//--------------IceStructureと同じ動きをするために一時的に作った関数__----------------------------------
	//ちゃんと実装すれば全部消せる
	void SetParticleNum(int pNum){	sm_particleNum = pNum;	m_iceStrct->SetParticleNum(pNum);	}
	void SetTetraNum(int tNum){		sm_tetraNum = tNum;	m_iceStrct->SetTetraNum(tNum);		}
	void SetClusterNum(int cNum){	sm_clusterNum = cNum; m_iceStrct->SetClusterNum(cNum);	}

	void InitTetraInfo(){	m_iceStrct->InitTetraInfo();	}
	void InitClusterInfo(){	m_iceStrct->InitClusterInfo();	}

	int GetPtoTNum(int i){	return m_iceStrct->GetPtoTNum(i);	}
	int GetTtoPNum(int i){	return m_iceStrct->GetTtoPNum(i);	}
	int GetPtoCNum(int i){	return m_iceStrct->GetPtoCNum(i);	}
	int GetCtoPNum(int i){	return m_iceStrct->GetCtoPNum(i);	}

	int GetPtoTIndx(int pIndx){	return m_iceStrct->GetPtoTIndx(pIndx);	}
	int GetTtoPIndx(int tIndx){ return m_iceStrct->GetTtoPIndx(tIndx);	}
	int GetPtoCIndx(int i){	return m_iceStrct->GetPtoCIndx(i);	}
	int GetCtoPIndx(int i){	return m_iceStrct->GetCtoPIndx(i);	}

	int GetPtoC(int i, int j, int k){	return m_iceStrct->GetPtoC(i, j, k);	}
	int GetPtoT(int i, int j, int k){	return m_iceStrct->GetPtoT(i, j, k);	}
	int GetTtoP(int i, int j){			return m_iceStrct->GetTtoP(i, j);		}

	int GetPrtclNum(){	return m_iceStrct->GetParticleNum();	}
	int GetTetraNum(){	return m_iceStrct->GetTetraNum();	}

	int GetPtoCMax(){ return m_iceStrct->GetPtoCMax();	}

	int GetNTNum(int i){ return m_iceStrct->GetNTNum(i);	}

	int GetNeighborTetra(int i, int j, int k){	return m_iceStrct->GetNeighborTetra(i, j, k);	}

	void SetPtoTIndx(int pIndx){	m_iceStrct->SetPtoTIndx(pIndx, m_iceStrct->GetPtoTNum(pIndx));	}
	void SetTtoPIndx(int tIndx){	m_iceStrct->SetTtoPIndx(tIndx, m_iceStrct->GetTtoPNum(tIndx));	}
	void SetPtoCIndx(int i, int j){	m_iceStrct->SetPtoCIndx(i, j);	}
	void SetCtoPIndx(int i, int j){	m_iceStrct->SetCtoPIndx(i, j);	}

	void SetPtoT(int i, int j, int k, int l){	m_iceStrct->SetPtoT(i, j, k, l);	}
	void SetTtoP(int i, vector<int> j)	{		m_iceStrct->SetTtoP(i, j);			}
	void SetPtoC(int i, int j, int k, int l, int m){	m_iceStrct->SetPtoC(i, j, k, l, m);	}
	void SetCtoP(int cIndx, const vector<int>& pIndxList, int* pLayerList){	m_iceStrct->SetCtoP(cIndx, pIndxList, pLayerList);	}

	void SetNeighborTetra(int i, int layer){	m_iceStrct->SetNeighborTetra(i, layer);	}

	void CountPtoT(int pIndx){	m_iceStrct->CountPtoT(pIndx);	}
	void CountTtoP(int tIndx){	m_iceStrct->CountTtoP(tIndx);	}
	void CountPtoC(int pIndx){	m_iceStrct->CountPtoC(pIndx);	}
	void CountCtoP(int cIndx){	m_iceStrct->CountCtoP(cIndx);	}
	//--------------__IceStructureと同じ動きをするために一時的に作った関数----------------------------------

};


//GPU処理
//TODO: 各クラスタの運動計算結果と固体構造のデータが必要なので，ここにおいている
extern void LaunchCalcAverageGPU
(
	int prtNum,
	float* sldPrtPos,
	float* sldPrtVel,
	float* sphPrtPos,
	float* sphPrtVel,
	float* smPrtPos,
	float* smPrtVel,
	int* smIndxSet,
	int* PtoCIndx,
	int* PtoC,
	int PNumMax,
	int PtoCMax,
	int PtoCParamSize
);

//TODO: いずれはクラスにしたいが，とりあえずここにおいている
extern void LaunchInterPolationGPU
(
	int prtNum,
	float* sldPrtPos,
	float* sldPrtVel,
	float* sphPrtPos,
	float* sphPrtVel
);

//PrefixSumのデータを用いてSM法で使うデータを更新
//この処理は，Ice_SMとSurf_SMのデータを用いるので，ここに置く必要がある
extern void LauchUpdateSMFromPath
(
	int prtNum,
	float* prtPos,
	float* prtVel, 
	//------------------SM----------------------
	float* orgPos,
	float* curPos,
	float* orgCm,
	float* curCm,
	float* clstrApq,
	float* vel,
	int* pIndxes, 
	int* startEndSet,
	//-----------------Path---------------------
	short int* PRTtoPTH,
	short int* PTHandPrfxSet,
	float* prfxPos,
	float* prfxApq,
	//-----------------Struct-------------------
	int* CtoP,
	int* CtoPNum,
	int CtoPSizeY,
	int CtoPSizeZ,

	float dt
);


#endif