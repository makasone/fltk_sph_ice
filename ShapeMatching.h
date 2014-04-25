/*!
  @file rx_shape_matching.h
	
  @brief Shape Matching法による弾性変形
  @ref M. Muller et al., "Meshless deformations based on shape matching", SIGGRAPH2005. 
 
  @author Makoto Fujisawa
  @date 2013-07
*/

#ifndef _RX_SHAPE_MATCHING_H_
#define _RX_SHAPE_MATCHING_H_


//-----------------------------------------------------------------------------
// インクルードファイル
//-----------------------------------------------------------------------------
#include <cstdio>
#include <cmath>
#include <cstdlib>

#include <vector>
#include <string>

#include "rx_utility.h"
#include "rx_matrix.h"

#include "rx_nnsearch.h"

using namespace std;

//! 衝突判定用関数
typedef void (*CollisionFunc)(Vec3&, Vec3&, Vec3&, int);

//-----------------------------------------------------------------------------
// Shape Matchingクラスの宣言
//-----------------------------------------------------------------------------
class rxShapeMatching
{
protected:
	// 形状データ
	int m_iNumVertices;
	vector<Vec3> m_vOrgPos;			//!< オリジナルの頂点位置
	vector<Vec3> m_vCurPos;			//!< 現在の頂点位置
	vector<Vec3> m_vNewPos;			//!< 次のステップの頂点位置
	vector<Vec3> m_vGoalPos;		//!< 目標頂点位置
	vector<double> m_vMass;			//!< 頂点質量(変形時の重み)
	vector<Vec3> m_vVel;			//!< 頂点速度
	vector<bool> m_vFix;			//!< 頂点固定フラグ

	// シミュレーションパラメータ
	double m_dDt;					//!< タイムステップ幅
	Vec3 m_v3Min, m_v3Max;			//!< シミュレーション空間の大きさ
	Vec3 m_v3Gravity;				//!< 重力加速度ベクトル
	
	double m_dAlpha;				//!< stiffnessパラメータ[0,1] (速度計算に使用)
	double m_dBeta;					//!< deformationパラメータ[0,1]

	bool m_bLinearDeformation;		//!< Linear/Quadratic deformation切り替えフラグ
	bool m_bVolumeConservation;		//!< 変形時の体積保存性(√det(A)で割るかどうか)

	int m_iObjectNo;				//!< オブジェクト番号

	CollisionFunc m_fpCollision;

public:
	//! コンストラクタとデストラクタ
	rxShapeMatching(int obj);
	~rxShapeMatching();

	void Clear();
	void AddVertex(const Vec3 &pos, double mass);

	void Update();

	// アクセスメソッド
	void SetTimeStep(double dt){ m_dDt = dt; }
	void SetSimulationSpace(Vec3 minp, Vec3 maxp){ m_v3Min = minp; m_v3Max = maxp; }
	void SetStiffness(double alpha, double beta){ m_dAlpha = alpha; m_dBeta = beta; }

	void SetCurrentPos(int oIndx, const Vec3 &pos){ m_vCurPos[oIndx] = pos; }
	void SetOriginalPos(int oIndx, Vec3 pos){ m_vOrgPos[oIndx] = pos; }
	void SetNewPos(int oIndx, Vec3 pos){ m_vNewPos[oIndx] = pos; }
	void SetGoalPos(int oIndx, Vec3 pos){ m_vGoalPos[oIndx] = pos; }

	void SetVelocity(int oIndx, const Vec3 &vel){ m_vVel[oIndx] = vel; }

	void SetCollisionFunc(CollisionFunc func){ m_fpCollision = func; }

	int GetNumVertices() const { return m_iNumVertices; }
	Vec3 GetVertexPos(int i){ return m_vCurPos[i]; }
	Vec3 GetNewPos(int i){ return m_vNewPos[i]; }
	Vec3 GetOriginalPos(int i){ return m_vOrgPos[i]; }
	Vec3 GetGoalPos(int i) { return m_vGoalPos[i]; }
	Vec3 GetVertexVel(int i){ return m_vVel[i]; }
	double GetMass(int i){ return m_vMass[i]; }

	void FixVertex(int i, const Vec3 &pos);
	void UnFixVertex(int i);
	bool IsFixed(int i) { return m_vFix[i]; }

protected:
	//! 頂点位置の初期化
	void initialize(void);

	// Shape Matching法の計算
	void calExternalForces(double dt);
	void calCollision(double dt);
	void shapeMatching(double dt);
	void integrate(double dt);

	void clamp(Vec3 &pos) const
	{
		if(pos[0] < m_v3Min[0]) pos[0] = m_v3Min[0];
		if(pos[0] > m_v3Max[0]) pos[0] = m_v3Max[0];
		if(pos[1] < m_v3Min[1]) pos[1] = m_v3Min[1];
		if(pos[1] > m_v3Max[1]) pos[1] = m_v3Max[1];
		if(pos[2] < m_v3Min[2]) pos[2] = m_v3Min[2];
		if(pos[2] > m_v3Max[2]) pos[2] = m_v3Max[2];
	}

};




/*!
 * Jacobi法による固有値の算出
 * @param[inout] a 実対称行列．計算後，対角要素に固有値が入る
 * @param[out] v 固有ベクトル(aと同じサイズ)
 * @param[in] n 行列のサイズ(n×n)
 * @param[in] eps 収束誤差
 * @param[in] iter_max 最大反復回数
 * @return 反復回数
 */
inline int EigenJacobiMethod(double *a, double *v, int n, double eps = 1e-8, int iter_max = 100)
{
	double *bim, *bjm;
	double bii, bij, bjj, bji;
 
	bim = new double[n];
	bjm = new double[n];
 
	for(int i = 0; i < n; ++i){
		for(int j = 0; j < n; ++j){
			v[i*n+j] = (i == j) ? 1.0 : 0.0;
		}
	}
 
	int cnt = 0;
	for(;;){
		int i = -1, j = -1;
 
		double x = 0.0;
		for(int ia = 0; ia < n; ++ia){
			for(int ja = 0; ja < n; ++ja){
				int idx = ia*n+ja;
				if(ia != ja && fabs(a[idx]) > x){
					i = ia;
					j = ja;
					x = fabs(a[idx]);
				}
			}
		}

		if(i == -1 || j == -1) return 0;
 
		double aii = a[i*n+i];
		double ajj = a[j*n+j];
		double aij = a[i*n+j];
 
		double m_dAlpha, m_dBeta;
		m_dAlpha = (aii-ajj)/2.0;
		m_dBeta  = sqrt(m_dAlpha*m_dAlpha+aij*aij);
 
		double st, ct;
		ct = sqrt((1.0+fabs(m_dAlpha)/m_dBeta)/2.0);    // sinθ
		st = (((aii-ajj) >= 0.0) ? 1.0 : -1.0)*aij/(2.0*m_dBeta*ct);    // cosθ
 
		// A = PAPの計算
		for(int m = 0; m < n; ++m){
			if(m == i || m == j) continue;
 
			double aim = a[i*n+m];
			double ajm = a[j*n+m];
 
			bim[m] =  aim*ct+ajm*st;
			bjm[m] = -aim*st+ajm*ct;
		}
 
		bii = aii*ct*ct+2.0*aij*ct*st+ajj*st*st;
		bij = 0.0;
 
		bjj = aii*st*st-2.0*aij*ct*st+ajj*ct*ct;
		bji = 0.0;
 
		for(int m = 0; m < n; ++m){
			a[i*n+m] = a[m*n+i] = bim[m];
			a[j*n+m] = a[m*n+j] = bjm[m];
		}
		a[i*n+i] = bii;
		a[i*n+j] = bij;
		a[j*n+j] = bjj;
		a[j*n+i] = bji;
 
		// V = PVの計算
		for(int m = 0; m < n; ++m){
			double vmi = v[m*n+i];
			double vmj = v[m*n+j];
 
			bim[m] =  vmi*ct+vmj*st;
			bjm[m] = -vmi*st+vmj*ct;
		}
		for(int m = 0; m < n; ++m){
			v[m*n+i] = bim[m];
			v[m*n+j] = bjm[m];
		}
 
		double e = 0.0;
		for(int ja = 0; ja < n; ++ja){
			for(int ia = 0; ia < n; ++ia){
				if(ia != ja){
					e += fabs(a[ja*n+ia]);
				}
			}
		}
		if(e < eps) break;
 
		cnt++;
		if(cnt > iter_max) break;
	}
 
	delete [] bim;
	delete [] bjm;
 
	return cnt;
}


/*!
 * 極分解で回転行列と対称行列に分解 A=RS
 * @param[in] A 入力行列
 * @param[out] R 回転行列(直交行列 R^-1 = R^T)
 * @param[out] S 対称行列
 */
inline void PolarDecomposition(const rxMatrix3 &A, rxMatrix3 &R, rxMatrix3 &S)
{
	R.makeIdentity();

	// S = (A^T A)^(1/2)を求める
	rxMatrix3 ATA;
	// (A^T A)の計算
	ATA = A.Transpose()*A;

	rxMatrix3 U;
	R.makeIdentity();
	U.makeIdentity();

	// (A^T A)を固有値分解して対角行列と直交行列を求める
	// M^(1/2) = U^T M' U 
	//  M = (A^T A), M':対角行列の平方根を取ったもの, U:直交行列, 
	EigenJacobiMethod(&ATA, &U, 3);

	// 対角行列の平方根をとって，逆行列計算のために逆数にしておく
	real l0 = (ATA(0,0) <= 0.0) ? 0.0 : 1.0/sqrt(ATA(0,0));
	real l1 = (ATA(1,1) <= 0.0) ? 0.0 : 1.0/sqrt(ATA(1,1));
	real l2 = (ATA(2,2) <= 0.0) ? 0.0 : 1.0/sqrt(ATA(2,2));

	// U^T M' U の逆行列計算
	rxMatrix3 S1;
	S1(0,0) = l0*U(0,0)*U(0,0) + l1*U(0,1)*U(0,1) + l2*U(0,2)*U(0,2);
	S1(0,1) = l0*U(0,0)*U(1,0) + l1*U(0,1)*U(1,1) + l2*U(0,2)*U(1,2);
	S1(0,2) = l0*U(0,0)*U(2,0) + l1*U(0,1)*U(2,1) + l2*U(0,2)*U(2,2);
	S1(1,0) = S1(0,1);
	S1(1,1) = l0*U(1,0)*U(1,0) + l1*U(1,1)*U(1,1) + l2*U(1,2)*U(1,2);
	S1(1,2) = l0*U(1,0)*U(2,0) + l1*U(1,1)*U(2,1) + l2*U(1,2)*U(2,2);
	S1(2,0) = S1(0,2);
	S1(2,1) = S1(1,2);
	S1(2,2) = l0*U(2,0)*U(2,0) + l1*U(2,1)*U(2,1) + l2*U(2,2)*U(2,2);

	R = A*S1;	// R = A S^-1
	S = R.Transpose()*A; // S = R^-1 A = R^T A
}









#endif
