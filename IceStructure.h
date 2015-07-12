/*!
  @file IceStructure.h
	
  @brief 氷構造を操作するクラス
 
  @author Ryou Nakasone
  @date 2013-07
*/
//---------------------------------------------------------------------------
//このファイルのクラス・関数が二回以上コンパイルされるのを防ぐための処理（インクルードガード）
#ifndef IceStructure_H
#define IceStructure_H
//---------------------------------------------------------------------------

#include <iostream>
#include <vector>
#include <algorithm>
#include <math.h>
#include <cuda_runtime.h>

#include <UtilityScript\mk_Vector2D.h>
#include <UtilityScript\mk_Vector3D.h>
#include "QueryCounter.h"

#include <IceTetrahedra.h>
#include "Ice_SM.h"

using namespace std;

class IceStructure
{
public:	//TODO: 全てpublicにしない
	IceStructure();
	IceStructure(int pNum, int cNum, int tNum, int layer);

	~IceStructure(void);

	static int* GetDevicePtoCIndxPointer(){	return sd_piPtoCIndx;	}
	static int* GetDevicePtoCPointer(){		return sd_piPtoC;	}

	//初期化
	void InitTetraInfo();									//四面体情報のメモリ確保
	void InitClusterInfo();									//クラスタ情報のメモリ確保
	void InitSelectCluster(vector<Ice_SM*>& iceSM);
	void InitSelectClusterFromClusterSet(vector<Ice_SM*>& iceSM);

	void InitGPU();											//GPU処理で用いるデータの初期化

	//相変化
	//融解処理
	void StepObjMelt(vector<unsigned>& pList, vector<unsigned>& cList, vector<unsigned>& tList, vector<unsigned>& cLayerList, vector<unsigned>& tLayerList);

	void SearchMeltParticle(vector<int>& pList);	
	void SearchReconstruct_Cluster_Melt(const vector<unsigned>& pList, vector<unsigned>& cList, vector<unsigned>& lList);
	void SearchReconstruct_Tetra_Melt(const vector<unsigned>& pList, vector<unsigned>& tList, vector<unsigned>& lList);

	void UpdateInfo_Melt_PandT(const vector<unsigned>& pList);
	void UpdateInfo_Melt_PandC(const vector<unsigned>& pList, const vector<unsigned>& cList);

	void UpdateInfo_Delete_TandP(const vector<int>& tList, const vector<int>& deleteList);

	void SetInfo_Cluster(const vector<unsigned>& pList, const vector<unsigned>& cList, const vector<unsigned>& lList);
	void SetInfo_Tetra(const vector<unsigned>& pList, const vector<unsigned>& tList, const vector<unsigned>& lList);

	void CheckDeleteCluster(void);
	void CheckDeleteTetra(vector<int>& tList, vector<int>& lList);
	void CheckSameTetra(int tIndx, const vector<int>& searchList, vector<int>& deleteList);
	void CheckIncludeTetra(int tIndx, const vector<int>& searchList, vector<int>& deleteList);

	void RemoveReconstructTetra(vector<int>& tList, vector<int>& lList, vector<int>& deleteTList);

	//凝固処理
	void StepObjFreeze();
	
	void SearchFreezeParticle(vector<int>& pList);
	void SearchReconstructCluster_Freeze(const vector<int>& pList, vector<int>& cList, vector<int>& lList);
	void SearchReconstructTetra_Freeze(const vector<int>& pList, vector<int>& tList, vector<int>& lList);

	void SetFreezeTetraInfo(vector<int>& pList);
	void SetFreezeClusterInfo(const vector<int>& pList);

	//アクセッサ
	void SetParticleNum(int pNum){	m_iPNum = pNum; }		//現在の粒子数
	void SetClusterNum(int cNum){	m_iCNum = cNum; }		//現在のクラスタ数
	void SetTetraNum(int tNum){		m_iTNum = tNum;	}		//現在の四面体数

	int GetParticleNum(void){	return m_iPNum;	}
	int GetClusterNum(void){	return m_iCNum;	}
	int GetTetraNum(void){		return m_iTNum;	}

	int GetPNumMax(void){	return m_iPNumMax;	}
	int GetCNumMax(void){	return m_iCNumMax;	}

	int GetCtoPMax(void){	return m_iCtoPMax;	}
	int GetPtoCMax(void){	return m_iPtoCMax;	}

	void CountTetrahedra(int tIndx, vector<int>& pList);					//四面体情報の登録
	void CountClusterParticle(int cIndx, vector<int>& pList, int pNum);		//クラスタの情報登録

	void CountPtoC(int pIndx){	m_piPtoCNum[pIndx]++;	}
	void CountPtoT(int pIndx){	m_piPtoTNum[pIndx]++;	}
	void CountCtoP(int cIndx){	m_piCtoPNum[cIndx]++;	}
	void CountCtoT(int cIndx){	}
	void CountTtoP(int tIndx){	m_piTtoPNum[tIndx]++;	}
	void CountTtoC(int tIndx){	}

	void CountNT(int tIndx){	m_piNTNum[tIndx]++;		}

	void SetPtoCIndx(int pIndx, int indx){	m_piPtoCIndx[pIndx] = indx;	} 
	void SetCtoPIndx(int cIndx, int indx){	m_piCtoPIndx[cIndx] = indx;	}
	void SetPtoTIndx(int pIndx, int indx){	m_piPtoTIndx[pIndx] = indx; }
	void SetTtoPIndx(int tIndx, int indx){	m_piTtoPIndx[tIndx] = indx; }

	int	 GetPtoTNum(int pIndx){		return m_piPtoTNum[pIndx];	}
	int	 GetPtoCNum(int pIndx){		return m_piPtoCNum[pIndx];	}
	int  GetTtoPNum(int tIndx){		return m_piTtoPNum[tIndx];	}
	int  GetCtoPNum(int cIndx){		return m_piCtoPNum[cIndx];	}

	int  GetNTNum(int tIndx){		return m_piNTNum[tIndx];	}

	int  GetPtoTIndx(int pIndx){	return m_piPtoTIndx[pIndx];	}
	int  GetPtoCIndx(int pIndx){	return m_piPtoCIndx[pIndx];	}
	int  GetTtoPIndx(int tIndx){	return m_piTtoPIndx[tIndx];	}
	int  GetCtoPIndx(int cIndx){	return m_piCtoPIndx[cIndx];	}

	int  GetPtoTFreeIndx(int pIndx);
	int  GetPtoCFreeIndx(int pIndx);

	void SetPtoT(int pIndx, int lIndx, int tIndx, int oIndx);				//粒子が属する四面体の登録　　粒子番号，粒子内順序，四面体番号，四面体内順序
	void SetPtoC(int pIndx, int lIndx, int cIndx, int oIndx, int layer);	//粒子が属するクラスタの登録　粒子番号，粒子内順序，クラスタ番号，クラスタ内順序
	
	void SetCtoP(int cIndx, const vector<int>& pIndxList, int* pLayerList);	//クラスタに属する粒子の登録　クラスタ番号，粒子リスト
	void SetTtoP(int tIndx, vector<int>& pIndxList);						//四面体に属する粒子の登録　　四面体番号，粒子リスト
	
	void SetTetraInfo(int tIndx, int* PtoTNum);								//初期化に用いる四面体情報の登録
	//void SetClusterInfo(int cIndx, int* PtoCNum, int 

	void SetNeighborTetra(int tIndx, int layer);
	void SetNeighborTetraFromLayer(int tIndx, int searchLayer, int deleteLayer);

	void SetSelectRadius(float radius){	m_selectRadius = radius;	}
	float GetSelectRadius(){	return m_selectRadius;	}

	//削除
	void DeleteTtoP(int tIndx, int lIndx);
	void DeletePtoC(int pIndx, int lIndx);
	void DeletePtoT(int pIndx, int lIndx);

	int  GetTtoP(int tIndx, int lIndx);

	int GetPtoC(int pIndx, int lIndx, int oIndx);
	int GetPtoT(int pIndx, int lIndx, int oIndx);
	const int& GetCtoP(const int& cIndx, const int& lIndx, const int& oIndx);

	int GetNeighborTetra(int tIndx, int lIndx, int oIndx);

	//初期化
	void ClearPtoT(int pIndx);
	void ClearPtoC(int pIndx);
	void ClearCtoP(int cIndx);
	void ClearTtoP(int tIndx);

	void ClearNeighborTetra(int tIndx);
	void ClearNeighborTetraFromLayer(int tIndx, int layer);

	//判定
	int  CheckNeighborTetra(int tIndx, int checkTIndx);

	//運動計算クラスタの選択
	void UpdateSelectCluster(const vector<unsigned>& prtList, vector<unsigned>& neighborClusters, const vector<Ice_SM*>& iceSM);
	void UpdateSelectClusterFronSet(const vector<unsigned>& prtList, vector<unsigned>& neighborClusters, const vector<Ice_SM*>& iceSM);
	void UpdateNeighborOfSelectCluster(vector<Ice_SM*>& iceSM);
	void ResetSelectCluster(vector<Ice_SM*>& iceSM);

	void UpdateMotionCalcCluster(unsigned cIndx, short unsigned num);
	inline short unsigned GetMotionCalcCluster(unsigned cIndx) const{	return m_psuSelectClusterIndx[cIndx];	}

	//デバッグ
	void DebugPtoT(int pIndx);
	void DebugPtoC(int pIndx);
	void DebugCtoP(int cIndx);
	void DebugTtoP(int tIndx);
	void DebugNeighborTetra(int tIndx);
	void DebugStepObjMelt(vector<unsigned>& pList, vector<unsigned>& cList);

//テスト
	void TestStepObjMelt(vector<unsigned>& pList, vector<unsigned>& cList, vector<unsigned>& tList, vector<unsigned>& cLayerList, vector<unsigned>& tLayerList);

	//フラグ　未使用
	void SetPFlag(int indx, bool state){	m_pbPFlag[indx] = state;	}
	void SetCFlag(int indx, bool state){	m_pbCFlag[indx] = state;	}
	void SetTFlag(int indx, bool state){	m_pbTFlag[indx] = state;	}

	bool GetPFlag(int indx){	return m_pbPFlag[indx];	}
	bool GetCFlag(int indx){	return m_pbCFlag[indx];	}
	bool GetTFlag(int indx){	return m_pbTFlag[indx];	}

	void ResetPFlag(int endNum);
	void ResetCFlag(int endNum);
	void ResetTFlag(int endNum);

	//---------------------------------------_GPU------------------------------------------
	int* GetDeviceCtoPNumPointer(){	return sd_piCtoPNum;	}
	int* GetDeviceCtoPointer(){		return sd_piCtoP;		}

	//---------------------------------------GPU_------------------------------------------

protected:

	int m_iPNumMax;								//最大粒子数
	int m_iPNum;								//現在の粒子数
	
	int m_iPtoCMax;								//粒子がクラスタに所属する最大数
	int m_iCtoPMax;								//クラスタが含む粒子の最大数

	int m_iTNumMax;
	int m_iTNum;

	int m_iPtoTMax;
	int m_iTtoPMax;

	int m_iCNumMax;								//最大クラスタ数
	int m_iCNum;								//現在のクラスタ数

	int m_iNeighborMax;							//近傍粒子の最大数
	int m_iLayer;

	float m_selectRadius;						//運動計算クラスタを選択する際の半径

	//とりあえず，擬似多次元配列にポインタは使わない
	//粒子→
	mk_Vector3D<int> m_mk3DiPtoC;				//粒子→クラスタ
	mk_Vector3D<int> m_mk3DiPtoT;				//粒子→四面体

	int*   m_piPtoCNum;							//粒子→クラスタの個数
	int*   m_piPtoTNum;							//粒子→四面体の個数

	int*   m_piPtoCIndx;
	int*   m_piPtoTIndx;

	//クラスタ→
	mk_Vector3D<int> m_mk3DiCtoP;				//クラスタ→粒子　クラスタに所属する粒子を返す　粒子の接続情報
	//mk_Vector3D<int> m_mk3DiCtoT;				//クラスタ→四面体

	int*   m_piCtoPNum;							//クラスタ→粒子の個数
	int*   m_piCtoTNum;							//クラスタ→四面体の個数

	int*   m_piCtoPIndx;
	int*   m_piCtoTIndx;

	//四面体→
	mk_Vector2D<int> m_mk2DiTtoP;				//四面体→粒子
	//mk_Vector3D<int> m_mk3DiTtoC;				//四面体→クラスタ

	int*   m_piTtoPNum;							//四面体→粒子の個数
	int*   m_piTtoCNum;							//四面体→クラスタの個数

	int*   m_piTtoPIndx;
	int*   m_piTtoCIndx;

	//近傍四面体
	mk_Vector3D<int> m_mk3DiNeighborTetra;		//近傍四面体

	int*   m_piNTNum;							//各近傍四面体の個数

	short unsigned* m_psuSelectClusterIndx;		//運動計算する粒子のフラグ　運動計算するなら1，しないなら0
												//本当はboolでいいが，今後のことも考えてshort unsignedにしておいた

//--------------------------------------GPU__------------------------------------------------------------
	//粒子→
	static int* sd_piPtoC;
	static int* sd_piPtoT;

	static int* sd_piPtoCNum;
	static int* sd_piPtoTNum;							//粒子→四面体の個数

	static int* sd_piPtoCIndx;
	static int* sd_piPtoTIndx;

	//クラスタ→
	static int* sd_piCtoP;

	static int* sd_piCtoPNum;
	//static int* sd_piCtoTNum;
	
	static int* sd_piCtoPIndx;
	//static int*   sd_piCtoTIndx;
	
	//四面体→
	static int* sd_piTtoP;

	static int* sd_piTtoPNum;
	static int* sd_piTtoCNum;

	static int* sd_piTtoPIndx;
	static int* sd_piTtoCIndx;

	//近傍四面体
	static int* sd_piNeighborTetra;

	static int* sd_piNeighborTetraTNum;
//--------------------------------------__GPU------------------------------------------------------------

	//探索用フラグ　未使用
	bool* m_pbPFlag;							//粒子
	bool* m_pbCFlag;							//クラスタ
	bool* m_pbTFlag;							//四面体

	//未使用
	vector<vector<int>>	m_vviEdgeToCluster;					//辺→クラスタ	ある辺がどのクラスタに所属しているか　必ず１個以上
	vector<vector<int>>	m_vviEdgeToParticle;				//辺→粒子　　　ある辺がどの粒子を接続しているか　　　必ず２個

	vector<vector<int>>	m_vviClusterToEdge;					//クラスタ→辺　あるクラスタにどの辺が所属しているか　必ず１個か３個
	vector<vector<int>>	m_vviParticleToEdge;				//粒子→辺　　　ある粒子はどの辺で接続されているか

};
#endif