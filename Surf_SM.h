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
	mk_Vector2D<Vec3>		m_vvvec_PrfxPos;		//パスごとの，現在の位置情報のためのprefixSum
	mk_Vector2D<rxMatrix3>	m_vvmat_PrfxApq;		//パスごとの，変形行列のためのprefixSum　（完全な変形行列Apqではないよ　詳しくは論文を参照）
	
	mk_Vector2D<int>	m_mk2DiPTHoPRT;				//パス→粒子
	vector<int>			m_viPRTtoPTH;				//粒子→パス　粒子は１つのパスにしか属さない

public:
	void MakePath(const float* pos, int num, int size);						//パス作成

	void CalcPrefixSum();							//PrefixSumの計算

	Vec3 GetPos();									//ある範囲における位置ベクトルの総和を返す
	rxMatrix3 GetApq();								//ある範囲における変形行列の総和を返す

	int GetPath(int path, int indx);				//パスに所属する粒子番号を返す　パス番号，順番
};

#endif