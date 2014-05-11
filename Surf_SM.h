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

#include <Ice_SM.h>

class IceStructure;

#include "rx_utility.h"
#include "rx_matrix.h"

#include <UtilityScript\mk_Vector2D.h>
#include <UtilityScript\mk_Vector3D.h>

using namespace std;

//表面粒子のみを用いたＳＭ法のためのクラス．
//たくさんのパスを管理して，高速化のためのパラメータを計算して返す．

class Surf_SM
{
private:
	mk_Vector2D<Vec3>		m_mk2Dvec3_PrfxPos;		//パスごとの，現在の位置情報のためのprefixSum
	mk_Vector2D<rxMatrix3>	m_mk2Dmat3_PrfxApq;		//パスごとの，変形行列のためのprefixSum　（完全な変形行列Apqではないよ　詳しくは論文を参照）
	
	mk_Vector2D<int>		m_mk2DiPTHtoPRT;		//パス→粒子
	mk_Vector2D<int>		m_mk2DiPRTtoPTH;		//粒子→パス　パス番号，パス内番号　粒子は１つのパスにしか属さない

	mk_Vector3D<int>		m_mk3DiPTHandPrfxSet;	//各クラス多における，パスとprefixSumの番地セット[0]：始点　[1]：終点　prefixSum番地のみでいい　path番号は粒子から経由して取得できる

	vector<Vec3> m_vvec3OrgPos;						//粒子の初期位置	初期情報だけ保存しておく
	//float* m_pos;									//粒子の現在地へのポインタ
	//vector<Vec3> m_vvec3OrgCm;					//クラスタの初期重心	融解・凝固時には更新

	IceStructure* m_strct;

public:
	void InitPath(const float* pos, const float* vel, const vector<Ice_SM*> iceSM, IceStructure* strct, int pthSize);	//パス作成
	void InitOrgPos(const float* pos, int pNum);
	void InitPathPrfxIndxSet(const vector<Ice_SM*> iceSM, IceStructure* strct);	//どのパスのどの部分が必要なのか，をクラスタごとに計算
	
	void SetPathDataApq();				

	void UpdatePrefixSum(const float* p, const float* vel);
	void UpdatePrefixSumPos(const float* pos, const float* vel);		//重心　PrefixSumの計算
	void UpdatePrefixSumApq(const float* pos);		//変形行列　PrefixSumの計算
	
	Vec3 ClacCmSum(int cIndx, const float* pos);						//prefixSumからクラスタの重心を計算して返す
	Vec3 CalcCmFromPrfxSm(int path, int start, int end);

	Vec3		GetPos();							//ある範囲における重心ベクトルの総和を返す
	rxMatrix3	GetApq();							//ある範囲における変形行列の総和を返す
	int			GetPath(int path, int indx);		//パスに所属する粒子番号を返す　パス番号，順番

	//デバッグ
	void DebugPathDataPos();
	void DebugPathDataApq();
	void DebugPathPrfxIndxSet();
};

#endif