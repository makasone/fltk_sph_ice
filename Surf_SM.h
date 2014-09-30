/*!
  @file Sur_SM.h

  @brief ShapeMatching法を基にした相変化シミュレーション
  @ref R. Diziol et al., "Robust Real-Time Deformation of Incompressible Suface Meshes", SIGGRAPH2005. 
 
  @author Ryo Nakasone
  @date 2013-10
*/

#ifndef _SURF_SM_H
#define _SURF_SM_H

#include <iostream>
#include <vector>
#include <algorithm>
#include <cuda_runtime.h>

#include <Ice_SM.h>
#include <IceStructure.h>

#include "rx_utility.h"
#include "rx_matrix.h"

#include <UtilityScript\mk_Vector2D.h>
#include <UtilityScript\mk_Vector3D.h>

using namespace std;

//表面粒子のみを用いたSM法のためのクラス．
//たくさんのパスを管理して，高速化のためのパラメータを計算して返す．

extern void LaunchUpdatePrefixSumGPU
(
	int prtNum,
	int PosSizeX,
	int PosSizeY,
	int ApqSizeX,
	int ApqSizeY,
	int* md_2DiPTHtoPRT,
	int* md_2DiPRTtoPTH,
	float* md_2Df3PrfxPos,
	float* md_2Df9PrfxApq,
	int* md_3DiPTHandPrfxSet,
	float* md_f3OrgPos,
	float* md_f3OrgCm,
	const float* md_fPos,
	const float* md_fVel
);

class Surf_SM
{
private:
	mk_Vector2D<Vec3>		m_mk2Dvec3_PrfxPos;		//パスごとの，現在の位置情報のためのprefixSum
	mk_Vector2D<rxMatrix3>	m_mk2Dmat3_PrfxApq;		//パスごとの，変形行列のためのprefixSum　（完全な変形行列Apqではないよ　詳しくは論文を参照）
	
	mk_Vector2D<int>		m_mk2DiPTHtoPRT;		//パス→粒子　不規則配列になる
	mk_Vector2D<int>		m_mk2DiPRTtoPTH;		//粒子→パス　0:パス番号，1:パス内番号　粒子は１つのパスにしか属さない

	mk_Vector3D<int>		m_mk3DiPTHandPrfxSet;	//各クラス多における，パスとprefixSumの番地セット[0]：始点　[1]：終点　prefixSum番地のみでいい　path番号は粒子から経由して取得できる

	vector<Vec3> m_vvec3OrgPos;						//粒子の初期位置
	vector<Vec3> m_vvec3OrgCm;						//クラスタの初期重心

	const float* m_fPos;							//粒子の位置へのポインタ
	const float* m_fVel;							//粒子の速度へのポインタ

	int m_iPrtclNum;

	IceStructure* m_strct;
	vector<Ice_SM*> m_iceSM;

//---------------------------------------------GPU__---------------------------------------------------------
	float* md_2Df3PrfxPos;							//パスごとの，現在の位置情報のためのprefixSum
	float* md_2Df9PrfxApq;							//パスごとの，変形行列のためのprefixSum　（完全な変形行列Apqではないよ　詳しくは論文を参照）
	
	int* md_2DiPTHtoPRT;							//パス→粒子　不規則配列になる
	int* md_2DiPRTtoPTH;							//粒子→パス　0:パス番号，1:パス内番号　粒子は１つのパスにしか属さない

	int* md_3DiPTHandPrfxSet;						//各クラス多における，パスとprefixSumの番地セット[0]：始点　[1]：終点　prefixSum番地のみでいい　path番号は粒子から経由して取得できる

	float* md_f3OrgPos;								//クラスタ内の粒子の初期位置		Ice_SMのデバイスポインタをコピーして使う
	float* md_f3OrgCm;								//クラスタの初期重心				Ice_SMのデバイスポインタをコピーして使う

	const float* md_fPos;							//粒子の位置へのポインタ
	const float* md_fVel;							//粒子の速度へのポインタ
//---------------------------------------------__GPU---------------------------------------------------------

public:
	void InitPath(const float* pos, const float* vel, const vector<Ice_SM*> iceSM, IceStructure* strct, int pthSize, int prtNum);	//パス作成
	void InitPathGPU();
	
	void InitPathPrfxIndxSet(const vector<Ice_SM*> iceSM, IceStructure* strct);				//どのパスのどの部分が必要なのか，をクラスタごとに計算
	void InitOrgPos(int prtNum);
	void InitOrgCm();

	void UpdatePrefixSum();
	void UpdatePrefixSumItr();
	void UpdatePrefixSumGPU();

	void UpdatePrefixSumPos();		//重心　PrefixSumの計算
	void UpdatePrefixSumApq();		//変形行列　PrefixSumの計算

	void UpdatePrefixSumPosItr();
	void UpdatePrefixSumApqItr();

	Vec3 CalcCmSum(const int& cIndx);											//prefixSumからクラスタの重心を計算して返す
	const Vec3 CalcCmFromPrfxSm(const int& path, const int& start, const int& end);

	rxMatrix3 CalcApqSum(const int& cIndx);
	const rxMatrix3 CalcApqFromPrfxSm(const int& path, const int& start, const int& end);


//デバッグ
	void DebugPathDataPos();
	void DebugPathDataApq();
	void DebugPathPrfxIndxSet();
};

#endif