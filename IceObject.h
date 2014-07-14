//全相変化オブジェクトを管理するクラス

#ifndef _ICE_OBJECT_
#define _ICE_OBJECT_

#include "Ice_SM.h"
#include "IceStructure.h"
#include "QueryCounter.h"

#include <time.h>
#include <cuda_runtime.h>
#include <cuda_gl_interop.h>

using namespace std;

class IceObject
{
private:
	//TODO::液体運動計算クラスも置きたい
	//TODO::四面体処理も置く
	static float* s_sphPrtPos;				//位置のポインタ　読み込み専用
	static float* s_sphPrtVel;				//速度のポインタ　読み込み専用

//--------------------------------------GPU__------------------------------------------------------------
	static float* sd_sldPrtPos;						//粒子位置のデバイスポインタ
	static float* sd_sldPrtVel;						//粒子速度のデバイスポインタ

	static float* sd_ObjPrtPos;						//総和計算による最終的な粒子位置
	static float* sd_ObjPrtVel;						//速度
//--------------------------------------__GPU------------------------------------------------------------

	static int sm_particleNum;						//粒子数
	static int sm_tetraNum;							//四面体数
	static int sm_clusterNum;						//クラスタ数

	//固体運動計算クラス
	vector<Ice_SM*> m_iceMove;

	//構造管理クラス
	IceStructure* m_iceStrct;

	static float* m_fInterPolationCoefficience;			//線形補間係数

public:
	IceObject(float* pos, float* vel, int pMaxNum, int cMaxNum, int tMaxNum);
	~IceObject();

	static void SetSPHPointer(float* pos, float* vel){		s_sphPrtPos = pos;	s_sphPrtVel = vel;	};

	void InitIceObj(int pMaxNum, int cMaxNum, int tMaxNum);
	void InitCluster(Ice_SM* sm){	m_iceMove.push_back(sm);	}	//ポインタをコピーしているだけ　一時的な実装
	static void InitInterPolation();

	void StepObjMove();
	void StepInterPolation();								//線形補間　いずれは処理が複雑になるのでクラスにしたい．
	
	void CalcAverageCPU(const int pIndx, Vec3& pos, Vec3& vel);
	void LinerInterPolationCPU(const int pIndx, const Vec3& pos, const Vec3& vel);

	//--------------IceStructureと同じ動きをするために一時的に作った関数__----------------------------------
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

	int GetParticleNum(){	return m_iceStrct->GetParticleNum();	}

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

#endif