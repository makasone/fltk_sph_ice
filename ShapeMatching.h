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

#include <cuda_runtime.h>
#include <cuda_gl_interop.h>

using namespace std;

//! 衝突判定用関数
typedef void (*CollisionFunc)(Vec3&, Vec3&, Vec3&, int);


//GPU処理
extern void LaunchShapeMatchingGPU(
	/*unsigned int num_particles,
	float (*pos)[2],
	float time,
	float dt*/);

#define MAXPARTICLE 100
#define SM_DIM 3


//-----------------------------------------------------------------------------
// Shape Matchingクラスの宣言
//-----------------------------------------------------------------------------
class rxShapeMatching
{
protected:
	// 形状データ
	double* m_pOrgPos;				//!< オリジナルの頂点位置
	double* m_pCurPos;				//!< 現在の頂点位置
	double* m_pNewPos;				//!< 次のステップの頂点位置
	double* m_pGoalPos;				//!< 目標頂点位置
	double* m_pMass;				//!< 頂点質量(変形時の重み)
	double* m_pVel;					//!< 頂点速度

	bool* m_pFix;					//!< 頂点固定フラグ

	int m_iNumVertices;

	// シミュレーションパラメータ
	double m_dDt;					//!< タイムステップ幅
	Vec3 m_v3Min, m_v3Max;			//!< シミュレーション空間の大きさ
	Vec3 m_v3Gravity;				//!< 重力加速度ベクトル
	
	//追加
	rxMatrix3 m_mtrxBeforeU;		//前フレームでの直交行列　warm startという高速化手法
	rxMatrix3 m_mtrxBeforeATA;

	double m_dAlpha;				//!< stiffnessパラメータ[0,1] (速度計算に使用)
	double m_dBeta;					//!< deformationパラメータ[0,1]

	bool m_bLinearDeformation;		//!< Linear/Quadratic deformation切り替えフラグ
	bool m_bVolumeConservation;		//!< 変形時の体積保存性(√det(A)で割るかどうか)

	int m_iObjectNo;				//!< オブジェクト番号

	CollisionFunc m_fpCollision;

	//GPU
	//デバイス側へのポインタ
	double* d_OrgPos;
	double* d_CurPos;
	double* d_NewPos;
	double* d_GoalPos;
	double* d_Mass;
	double* d_Vel;

	bool* d_Fix;

public:
	//! コンストラクタとデストラクタ
	rxShapeMatching(int obj);
	~rxShapeMatching();

	void InitGPU();

	void Clear();
	void AddVertex(const Vec3 &pos, double mass);

	void Update();

	// アクセスメソッド
	void SetTimeStep(double dt){ m_dDt = dt; }
	void SetSimulationSpace(Vec3 minp, Vec3 maxp){ m_v3Min = minp; m_v3Max = maxp; }
	void SetStiffness(double alpha, double beta){ m_dAlpha = alpha; m_dBeta = beta; }

	void SetCurrentPos(int oIndx, const Vec3 &pos)
	{ 
		m_pCurPos[oIndx*SM_DIM+0] = pos[0];
		m_pCurPos[oIndx*SM_DIM+1] = pos[1];
		m_pCurPos[oIndx*SM_DIM+2] = pos[2];
	}

	void SetOriginalPos(int oIndx, Vec3 pos)
	{ 
		m_pOrgPos[oIndx*SM_DIM+0] = pos[0];
		m_pOrgPos[oIndx*SM_DIM+1] = pos[1];
		m_pOrgPos[oIndx*SM_DIM+2] = pos[2];
	}

	void SetNewPos(int oIndx, Vec3 pos)
	{
		m_pNewPos[oIndx*SM_DIM+0] = pos[0];
		m_pNewPos[oIndx*SM_DIM+1] = pos[1];
		m_pNewPos[oIndx*SM_DIM+2] = pos[2];
	}

	void SetGoalPos(int oIndx, Vec3 pos)
	{
		m_pGoalPos[oIndx*SM_DIM+0] = pos[0];
		m_pGoalPos[oIndx*SM_DIM+1] = pos[1];
		m_pGoalPos[oIndx*SM_DIM+2] = pos[2];
	}

	void SetVelocity(int oIndx, const Vec3 &vel)
	{
		m_pVel[oIndx*SM_DIM+0] = vel[0];
		m_pVel[oIndx*SM_DIM+1] = vel[1];
		m_pVel[oIndx*SM_DIM+2] = vel[2];
	}

	void SetCollisionFunc(CollisionFunc func){ m_fpCollision = func; }

	int GetNumVertices() const { return m_iNumVertices; }

	const Vec3 GetVertexPos(int i){ return Vec3(m_pCurPos[i*SM_DIM+0], m_pCurPos[i*SM_DIM+1], m_pCurPos[i*SM_DIM+2]); }
	const Vec3 GetNewPos(int i){ return Vec3(m_pNewPos[i*SM_DIM+0], m_pNewPos[i*SM_DIM+1], m_pNewPos[i*SM_DIM+2]); }
	const Vec3 GetOrgPos(int i){ return Vec3(m_pOrgPos[i*SM_DIM+0], m_pOrgPos[i*SM_DIM+1], m_pOrgPos[i*SM_DIM+2]); }
	const Vec3 GetGoalPos(int i) { return Vec3(m_pGoalPos[i*SM_DIM+0], m_pGoalPos[i*SM_DIM+1], m_pGoalPos[i*SM_DIM+2]); }
	const Vec3 GetVertexVel(int i){ return Vec3(m_pVel[i*SM_DIM+0], m_pVel[i*SM_DIM+1], m_pVel[i*SM_DIM+2]); }
	double GetMass(int i){ return m_pMass[i]; }

	void FixVertex(int i, const Vec3 &pos);
	void UnFixVertex(int i);
	bool IsFixed(int i) { return m_pFix[i]; }

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

	void clamp(double* pos, int cIndx) const
	{
		if(pos[cIndx+0] < m_v3Min[0]) pos[cIndx+0] = m_v3Min[0];
		if(pos[cIndx+0] > m_v3Max[0]) pos[cIndx+0] = m_v3Max[0];
		if(pos[cIndx+1] < m_v3Min[1]) pos[cIndx+1] = m_v3Min[1];
		if(pos[cIndx+1] > m_v3Max[1]) pos[cIndx+1] = m_v3Max[1];
		if(pos[cIndx+2] < m_v3Min[2]) pos[cIndx+2] = m_v3Min[2];
		if(pos[cIndx+2] > m_v3Max[2]) pos[cIndx+2] = m_v3Max[2];
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
 
	//for(int i = 0; i < n; ++i){
	//	for(int j = 0; j < n; ++j){
	//		v[i*n+j] = (i == j) ? 1.0 : 0.0;
	//	}
	//}
 
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

/*!
 * 極分解で回転行列と対称行列に分解 A=RS
 * warm startあり	1フレーム目は単位行列　効果があったかは不明
 * @param[in] A 入力行列
 * @param[out] R 回転行列(直交行列 R^-1 = R^T)
 * @param[out] S 対称行列
 */
inline void PolarDecomposition(const rxMatrix3 &A, rxMatrix3 &R, rxMatrix3 &S, rxMatrix3& bfrU)
{
	// S = (A^T A)^(1/2)を求める
	//warm start	1フレーム目は単位行列
	// (A^T A)の計算
	rxMatrix3 ATA(bfrU.Transpose()*A.Transpose()*A*bfrU);

	R.makeIdentity();

	// (A^T A)を固有値分解して対角行列と直交行列を求める
	// M^(1/2) = U^T M' U 
	//  M = (A^T A), M':対角行列の平方根を取ったもの, U:直交行列, 
	EigenJacobiMethod(&ATA, &bfrU, 3);

	// 対角行列の平方根をとって，逆行列計算のために逆数にしておく
	real l0 = (ATA(0,0) <= 0.0) ? 0.0 : 1.0/sqrt(ATA(0,0));
	real l1 = (ATA(1,1) <= 0.0) ? 0.0 : 1.0/sqrt(ATA(1,1));
	real l2 = (ATA(2,2) <= 0.0) ? 0.0 : 1.0/sqrt(ATA(2,2));

	// U^T M' U の逆行列計算
	rxMatrix3 S1;
	S1(0,0) = l0*bfrU(0,0)*bfrU(0,0) + l1*bfrU(0,1)*bfrU(0,1) + l2*bfrU(0,2)*bfrU(0,2);
	S1(0,1) = l0*bfrU(0,0)*bfrU(1,0) + l1*bfrU(0,1)*bfrU(1,1) + l2*bfrU(0,2)*bfrU(1,2);
	S1(0,2) = l0*bfrU(0,0)*bfrU(2,0) + l1*bfrU(0,1)*bfrU(2,1) + l2*bfrU(0,2)*bfrU(2,2);
	S1(1,0) = S1(0,1);
	S1(1,1) = l0*bfrU(1,0)*bfrU(1,0) + l1*bfrU(1,1)*bfrU(1,1) + l2*bfrU(1,2)*bfrU(1,2);
	S1(1,2) = l0*bfrU(1,0)*bfrU(2,0) + l1*bfrU(1,1)*bfrU(2,1) + l2*bfrU(1,2)*bfrU(2,2);
	S1(2,0) = S1(0,2);
	S1(2,1) = S1(1,2);
	S1(2,2) = l0*bfrU(2,0)*bfrU(2,0) + l1*bfrU(2,1)*bfrU(2,1) + l2*bfrU(2,2)*bfrU(2,2);

	R = A*S1;	// R = A S^-1
	S = R.Transpose()*A; // S = R^-1 A = R^T A
}

#endif
