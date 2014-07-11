//SM法のGPU実装 運動計算
//とりあえずやってみて，データ構造は後で直す

#ifndef _GPU_SHAPE_MATCHING_H_
#define _GPU_SHAPE_MATCHING_H_

#include <math.h>
#include <stdio.h>
#include <cuda_runtime.h>

//以下はうまくインクルードできなかった
//#include <cstdio>
//#include <cmath>
//#include <cstdlib>
//
//#include "rx_utility.h"
//#include <rx_matrix.h>
//
//#include "rx_nnsearch.h"

//#include <rx_cu_common.cu>
#include <rx_cu_common.cuh>	//先生が定義した便利機能が使える　なぜか、インクルードすると赤くならない

#define SM_DIM 3

//引き渡す変数名が間違っていてもエラーが出ないので注意

void LaunchShapeMathcingGPU(float* prtPos, float* prtVel, float* orgPos, float* curPos, float* vel, int* pIndxes, int* indxSet, float dt, int prtNum);
__global__ void Update(float* prtPos, float* prtVel, float* orgPos, float* curPos, float* vel, int* pIndxes, int* indxSet, float dt, int prtNum);
__device__ void ExternalForce(float* prtPos, float* prtVel, float* curPos, float* vel, int* pIndxes, int* indxSet, float dt, int prtNum);
__device__ void ProjectPos(float* prtPos, float* prtVel, float* orgPos, float* curPos, float* vel, int* pIndxes, int* indxSet, float dt, int prtNum);
__device__ void Integrate(float* prtPos, float* prtVel, float* curPos, float* vel, int* pIndxes, int* indxSet, float dt, int prtNum);

__device__ void PolarDecomposition(matrix3x3 &A, matrix3x3 &R, matrix3x3 &S);

//行列演算
//TODO::てきとうなのでもう少し使いやすく
__device__ void MakeIdentity(matrix3x3 &M);
__device__ matrix3x3 Transpose(const matrix3x3 &M);
__device__ matrix3x3 Inverse(const matrix3x3 &M);

__device__ matrix3x3 Multiple(matrix3x3 &M1, matrix3x3 &M2);
__device__ float3 Multiple(matrix3x3 &M1, float3& V);

//GPU処理
void LaunchShapeMatchingGPU(float* prtPos, float* prtVel, float* orgPos, float* curPos, float* vel, int* pIndxes, int* indxSet, float dt, int prtNum)
{
	//printf("LaunchGPUKernel");

	dim3 grid(1, 1);
	dim3 block(729, 1, 1);

	//運動計算
	Update <<< grid , block >>> (prtPos, prtVel, orgPos, curPos, vel, pIndxes, indxSet, dt, prtNum);
}


//GPUの位置・速度更新
__global__
void Update(
	float* prtPos,
	float* prtVel, 
	float* orgPos,
	float* curPos, 
	float* vel,
	int* pIndxes,
	int* indxSet,
	float dt,
	int prtNum)
{
	ExternalForce(prtPos, prtVel, curPos, vel, pIndxes, indxSet, dt, prtNum);
	ProjectPos(prtPos, prtVel, orgPos, curPos, vel, pIndxes, indxSet, dt, prtNum);
	Integrate(prtPos, prtVel, curPos, vel, pIndxes, indxSet, dt, prtNum);
}

__device__
	void ExternalForce(float* prtPos, float* prtVel, float* curPos, float* vel, int* pIndxes, int* indxSet, float dt, int prtNum)
{
	//計算するクラスタの判定
	int clusterIndx = threadIdx.x;

	int startIndx = indxSet[clusterIndx*2+0];
	int endIndx = indxSet[clusterIndx*2+1];
	//printf("startIndx = %d, endIndx = %d\n", startIndx, endIndx);

	// 重力の影響を付加，速度を反映
	for(int i = startIndx; i < endIndx+1; ++i)
	{
		int pIndx = pIndxes[i]*4;
		int cIndx = i*SM_DIM;

		for(int j = 0; j < SM_DIM; ++j)
		{
			int jpIndx = pIndx+j;
			int jcIndx = cIndx+j;

			curPos[jcIndx] = prtPos[jpIndx]+prtVel[jpIndx]*dt;
		}
		
		//printf("prtVel(%f, %f, %f)\n", prtVel[pIndx], prtVel[pIndx+1], prtVel[pIndx+2]);
		//printf("pIndxes[i] = %d, prtPos(%f, %f, %f)\n", pIndxes[i], prtPos[pIndx], prtPos[pIndx+1], prtPos[pIndx+2]);
		//printf("gpu:: i = %d, cIndx = %d, curPos(%f, %f, %f)\n", i, cIndx, curPos[cIndx], curPos[cIndx+1], curPos[cIndx+2]);
	}



	// 境界壁の影響
	//処理がかなり重くなるが，安定はするみたい
	//float res = 0.9;	// 反発係数
	//for(int i = 0; i < prtNum; ++i){
	//	//if(m_pFix[i]) continue;
	//	//Vec3 &p = m_vCurPos[i];
	//	//Vec3 &np = m_vNewPos[i];
	//	//Vec3 &v = m_vVel[i];
	//	//if(np[0] < m_v3Min[0] || np[0] > m_v3Max[0]){
	//	//	np[0] = p[0]-v[0]*dt*res;
	//	//	np[1] = p[1];
	//	//	np[2] = p[2];
	//	//}
	//	//if(np[1] < m_v3Min[1] || np[1] > m_v3Max[1]){
	//	//	np[1] = p[1]-v[1]*dt*res;
	//	//	np[0] = p[0] ;
	//	//	np[2] = p[2];
	//	//}
	//	//if(np[2] < m_v3Min[2] || np[2] > m_v3Max[2]){
	//	//	np[2] = p[2]-v[2]*dt*res;
	//	//	np[0] = p[0];
	//	//	np[1] = p[1];
	//	//}

	//	//clamp(curPos, i*SM_DIM);
	//}
}

__device__
	void ProjectPos(float* prtPos, float* prtVel, float* orgPos, float* curPos, float* vel, int* pIndxes, int* indxSet, float dt, int prtNum)
{
	//計算するクラスタの判定
	int clusterIndx = threadIdx.x;

	int startIndx = indxSet[clusterIndx*2+0];
	int endIndx = indxSet[clusterIndx*2+1];

	if(prtNum <= 1) return;

	float3 cm = make_float3(0.0, 0.0, 0.0);
	float3 cm_org = make_float3(0.0, 0.0, 0.0);	// 重心

	float mass = 0.0f;	// 総質量

	// 重心座標の計算
	for(int i = startIndx; i < endIndx+1; ++i)
	{
		//float m = m_pMass[i];
		float m = 1.0f;
		//if(m_pFix[i]) m *= 300.0;	// 固定点の質量を大きくする
		
		int cIndx = i*SM_DIM;
		mass += m;

		cm.x += curPos[cIndx+0]*m;
		cm.y += curPos[cIndx+1]*m;
		cm.z += curPos[cIndx+2]*m;

		cm_org.x += orgPos[cIndx+0]*m;
		cm_org.y += orgPos[cIndx+1]*m;
		cm_org.z += orgPos[cIndx+2]*m;
	}

	cm.x /= mass;
	cm.y /= mass;
	cm.z /= mass;
	
	cm_org.x /= mass;
	cm_org.y /= mass;
	cm_org.z /= mass;

	matrix3x3 Apq, Aqq;
	float3 p, q;

	// Apq = Σmpq^T
	// Aqq = Σmqq^T
	for(int i = startIndx; i < endIndx+1; ++i)
	{
		int cIndx = i*SM_DIM;

		p.x = curPos[cIndx+0]-cm.x;
		p.y = curPos[cIndx+1]-cm.y;
		p.z = curPos[cIndx+2]-cm.z;

		q.x = orgPos[cIndx+0]-cm_org.x;
		q.y = orgPos[cIndx+1]-cm_org.y;
		q.z = orgPos[cIndx+2]-cm_org.z;

		float m = 1.0f;

		Apq.e[0].x += m*p.x*q.x;
		Apq.e[0].y += m*p.x*q.y;
		Apq.e[0].z += m*p.x*q.z;
		Apq.e[1].x += m*p.y*q.x;
		Apq.e[1].y += m*p.y*q.y;
		Apq.e[1].z += m*p.y*q.z;
		Apq.e[2].x += m*p.z*q.x;
		Apq.e[2].y += m*p.z*q.y;
		Apq.e[2].z += m*p.z*q.z;

		Aqq.e[0].x += m*q.x*q.x;
		Aqq.e[0].y += m*q.x*q.y;
		Aqq.e[0].z += m*q.x*q.z;
		Aqq.e[1].x += m*q.y*q.x;
		Aqq.e[1].y += m*q.y*q.y;
		Aqq.e[1].z += m*q.y*q.z;
		Aqq.e[2].x += m*q.z*q.x;
		Aqq.e[2].y += m*q.z*q.y;
		Aqq.e[2].z += m*q.z*q.z;
	}

	////Apqの行列式を求め，反転するかを判定
	////不安定な場合が多いので×
	////if( Apq.Determinant() < 0.0 && m_iNumVertices >= 4)
	////{
	//	//cout << "before det < 0" << endl;
	//	//１　符号を反転
	//	//Apq(0,2) = -Apq(0,2);
	//	//Apq(1,2) = -Apq(1,2);
	//	//Apq(2,2) = -Apq(2,2);

	//	//２　a2とa3を交換
	//	//float tmp;
	//	//tmp = Apq(0,2);
	//	//Apq(0,2) = Apq(0,1);
	//	//Apq(0,1) = tmp;

	//	//tmp = Apq(1,2);
	//	//Apq(1,2) = Apq(1,1);
	//	//Apq(1,1) = tmp;

	//	//tmp = Apq(2,2);
	//	//Apq(2,2) = Apq(2,1);
	//	//Apq(2,1) = tmp;
	////}
	//

	matrix3x3 R, S;
	////PolarDecomposition(Apq, R, S, m_mtrxBeforeU);
	PolarDecomposition(Apq, R, S);

	//if(m_bLinearDeformation)
	{
		// Linear Deformations
		matrix3x3 A;
		A = Multiple(Apq, Inverse(Aqq));	// A = Apq*Aqq^-1

		//// 体積保存のために√(det(A))で割る
		//if(m_bVolumeConservation){
		//	float det = fabs(A.Determinant());
		//	if(det > RX_FEQ_EPS){
		//		det = 1.0/sqrt(det);
		//		if(det > 2.0) det = 2.0;
		//		A *= det;
		//	}
		//}

		// 目標座標を計算し，現在の頂点座標を移動
		for(int i = startIndx; i < endIndx+1; ++i)
		{
			//if(m_pFix[i]) continue;

			int cIndx = i*SM_DIM;

			// 回転行列Rの代わりの行列RL=βA+(1-β)Rを計算
			q.x = orgPos[cIndx+0]-cm_org.x;
			q.y = orgPos[cIndx+1]-cm_org.y;
			q.z = orgPos[cIndx+2]-cm_org.z;

			//Vec3 gp(R*q+cm);
			float3 Rq = Multiple(R, q);
			float3 gp;
			gp.x = Rq.x + cm.x;
			gp.y = Rq.y + cm.y;
			gp.z = Rq.z + cm.z;
			
			curPos[cIndx+0] += (gp.x-curPos[cIndx+0])*1.0f/*m_dAlphas[i]*/;
			curPos[cIndx+1] += (gp.y-curPos[cIndx+1])*1.0f/*m_dAlphas[i]*/;
			curPos[cIndx+2] += (gp.z-curPos[cIndx+2])*1.0f/*m_dAlphas[i]*/;
		}
	}
}

/*!
 * 速度と位置の更新
 *  - 新しい位置と現在の位置座標から速度を算出
 * @param[in] dt タイムステップ幅
 */
__device__
	void Integrate(float* prtPos, float* prtVel, float* curPos, float* vel, int* pIndxes, int* indxSet, float dt, int prtNum)
{
	//計算するクラスタの判定
	int clusterIndx = threadIdx.x;

	int startIndx = indxSet[clusterIndx*2+0];
	int endIndx = indxSet[clusterIndx*2+1];

	float dt1 = 1.0f/dt;
	float gravity[3] = {0.0f, -9.81f, 0.0f};

	for(int i = startIndx; i < endIndx+1; ++i)
	{
		int pIndx = pIndxes[i]*4;

		for(int j = 0; j < SM_DIM; j++)
		{
			int cIndx = i*SM_DIM+j;

			vel[cIndx] = (curPos[cIndx] - prtPos[pIndx+j]) * dt1 + gravity[j] * dt * 0.1f;
		}
	}
}

__device__  
	void clamp(float* pos, int cIndx)
{
	//if(pos[cIndx+0] < m_v3Min[0]) pos[cIndx+0] = m_v3Min[0];
	//if(pos[cIndx+0] > m_v3Max[0]) pos[cIndx+0] = m_v3Max[0];
	//if(pos[cIndx+1] < m_v3Min[1]) pos[cIndx+1] = m_v3Min[1];
	//if(pos[cIndx+1] > m_v3Max[1]) pos[cIndx+1] = m_v3Max[1];
	//if(pos[cIndx+2] < m_v3Min[2]) pos[cIndx+2] = m_v3Min[2];
	//if(pos[cIndx+2] > m_v3Max[2]) pos[cIndx+2] = m_v3Max[2];
}

/*!
 * Jacobi法による固有値の算出
 * @param[inout] a 実対称行列．計算後，対角要素に固有値が入る
 * @param[out] v 固有ベクトル(aと同じサイズ)
 * @param[in] n 行列のサイズ(n×n)
 * @param[in] eps 収束誤差
 * @param[in] iter_max 最大反復回数
 * @return 反復回数
 */
__device__ 
 int d_EigenJacobiMethod(double *a, double *v, int n, double eps = 1e-8, int iter_max = 100)
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
 
		float m_dAlpha, m_dBeta;
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
__device__ 
	void PolarDecomposition(matrix3x3 &A, matrix3x3 &R, matrix3x3 &S)
{
	// S = (A^T A)^(1/2)を求める
	matrix3x3 ATA;
	ATA = Multiple(Transpose(A), A);	// (A^T A)の計算

	MakeIdentity(R);

	// (A^T A)を固有値分解して対角行列と直交行列を求める
	// M^(1/2) = U^T M' U 
	//  M = (A^T A), M':対角行列の平方根を取ったもの, U:直交行列

	//行列をfloatに変換
	double* pATA = new double[9];
	double* pU = new double[9];

	for(int i = 0; i < 3; i++)
	{
		int indx = i*3;
		pATA[indx+0] = ATA.e[i].x;
		pATA[indx+1] = ATA.e[i].y;
		pATA[indx+2] = ATA.e[i].z;
	}

	d_EigenJacobiMethod(pATA, pU, 3);

	//float*を行列に変換	
	matrix3x3 U;

	for(int i = 0; i < 3; i++)
	{
		int indx = i*3;
		ATA.e[i].x = pATA[indx+0];
		ATA.e[i].y = pATA[indx+1];
		ATA.e[i].z = pATA[indx+2];

		U.e[i].x = pU[indx+0];
		U.e[i].y = pU[indx+1];
		U.e[i].z = pU[indx+2];
	}

	// 対角行列の平方根をとって，逆行列計算のために逆数にしておく
	double l0 = (ATA.e[0].x <= 0.0) ? 0.0 : 1.0/sqrt(ATA.e[0].x);
	double l1 = (ATA.e[1].y <= 0.0) ? 0.0 : 1.0/sqrt(ATA.e[1].y);
	double l2 = (ATA.e[2].z <= 0.0) ? 0.0 : 1.0/sqrt(ATA.e[2].z);

	//// U^T M' U の逆行列計算
	matrix3x3 S1;
	S1.e[0].x = l0*U.e[0].x*U.e[0].x + l1*U.e[0].y*U.e[0].y + l2*U.e[0].z*U.e[0].z;
	S1.e[0].y = l0*U.e[0].x*U.e[1].x + l1*U.e[0].y*U.e[1].y + l2*U.e[0].z*U.e[1].z;
	S1.e[0].z = l0*U.e[0].x*U.e[2].x + l1*U.e[0].y*U.e[2].y + l2*U.e[0].z*U.e[2].z;
	S1.e[1].x = S1.e[0].y;
	S1.e[1].y = l0*U.e[1].x*U.e[1].x + l1*U.e[1].y*U.e[1].y + l2*U.e[1].z*U.e[1].z;
	S1.e[1].z = l0*U.e[1].x*U.e[2].x + l1*U.e[1].y*U.e[2].y + l2*U.e[1].z*U.e[2].z;
	S1.e[2].x = S1.e[0].z;
	S1.e[2].y = S1.e[1].z;
	S1.e[2].z = l0*U.e[2].x*U.e[2].x + l1*U.e[2].y*U.e[2].y + l2*U.e[2].z*U.e[2].z;

	R = Multiple(A, S1);	// R = A S^-1
	S = Multiple(Transpose(R), A); // S = R^-1 A = R^T A

	delete[] pATA;
	delete[] pU;
}

//配列初期化
__device__
	void MakeIdentity(matrix3x3 &M)
{
	M.e[0].x = 1.0;
	M.e[0].y = 0.0;
	M.e[0].z = 0.0;

	M.e[1].x = 0.0;
	M.e[1].y = 1.0;
	M.e[1].z = 0.0;

	M.e[2].x = 0.0;
	M.e[2].y = 0.0;
	M.e[2].z = 1.0;
}

//転置行列
__device__
	matrix3x3 Transpose(const matrix3x3 &M)
{
	matrix3x3 T;

	T.e[0].x = M.e[0].x;
	T.e[0].y = M.e[1].x;
	T.e[0].z = M.e[2].x;

	T.e[1].x = M.e[0].y;
	T.e[1].y = M.e[1].y;
	T.e[1].z = M.e[2].y;

	T.e[2].x = M.e[0].z;
	T.e[2].y = M.e[1].z;
	T.e[2].z = M.e[2].z;

	return T;
}

//逆行列
__device__
	matrix3x3 Inverse(const matrix3x3 &M)
{
	matrix3x3 I;

	I.e[0].x = M.e[0].x;
	I.e[0].y = M.e[1].x;
	I.e[0].z = M.e[2].x;

	I.e[1].x = M.e[0].y;
	I.e[1].y = M.e[1].y;
	I.e[1].z = M.e[2].y;

	I.e[2].x = M.e[0].z;
	I.e[2].y = M.e[1].z;
	I.e[2].z = M.e[2].z;

	float d = M.e[0].x*M.e[1].y*M.e[2].z- 
			 M.e[0].x*M.e[2].y*M.e[1].z+ 
			 M.e[1].x*M.e[2].y*M.e[0].z- 
			 M.e[1].x*M.e[0].y*M.e[2].z+ 
			 M.e[2].x*M.e[0].y*M.e[1].z- 
			 M.e[2].x*M.e[1].y*M.e[0].z;

	if(d == 0) d = 1;

	I.e[0].x =  (M.e[1].y*M.e[2].z-M.e[1].z*M.e[2].y)/d;
	I.e[0].y = -(M.e[0].y*M.e[2].z-M.e[0].z*M.e[2].y)/d;
	I.e[0].z =  (M.e[0].y*M.e[1].z-M.e[0].z*M.e[1].y)/d;
	I.e[1].x = -(M.e[1].x*M.e[2].z-M.e[1].z*M.e[2].x)/d;
	I.e[1].y =  (M.e[0].x*M.e[2].z-M.e[0].z*M.e[2].x)/d;
	I.e[1].z = -(M.e[0].x*M.e[1].z-M.e[0].z*M.e[1].x)/d;
	I.e[2].x =  (M.e[1].x*M.e[2].y-M.e[1].y*M.e[2].x)/d;
	I.e[2].y = -(M.e[0].x*M.e[2].y-M.e[0].y*M.e[2].x)/d;
	I.e[2].z =  (M.e[0].x*M.e[1].y-M.e[0].y*M.e[1].x)/d;

	return I;
}

__device__ 
	matrix3x3 Multiple(matrix3x3 &M1, matrix3x3 &M2)
{
	matrix3x3 M;

	M.e[0].x = M1.e[0].x * M2.e[0].x + M1.e[0].y * M2.e[1].x + M1.e[0].z * M2.e[2].x;
	M.e[0].y = M1.e[0].x * M2.e[0].y + M1.e[0].y * M2.e[1].y + M1.e[0].z * M2.e[2].y;
	M.e[0].z = M1.e[0].x * M2.e[0].z + M1.e[0].y * M2.e[1].z + M1.e[0].z * M2.e[2].z;
	 
	M.e[1].x = M1.e[1].x * M2.e[0].x + M1.e[1].y * M2.e[1].x + M1.e[1].z * M2.e[2].x;
	M.e[1].y = M1.e[1].x * M2.e[0].y + M1.e[1].y * M2.e[1].y + M1.e[1].z * M2.e[2].y;
	M.e[1].z = M1.e[1].x * M2.e[0].z + M1.e[1].y * M2.e[1].z + M1.e[1].z * M2.e[2].z;
	 
	M.e[2].x = M1.e[2].x * M2.e[0].x + M1.e[2].y * M2.e[1].x + M1.e[2].z * M2.e[2].x;
	M.e[2].y = M1.e[2].x * M2.e[0].y + M1.e[2].y * M2.e[1].y + M1.e[2].z * M2.e[2].y;
	M.e[2].z = M1.e[2].x * M2.e[0].z + M1.e[2].y * M2.e[1].z + M1.e[2].z * M2.e[2].z;

	return M;
}

__device__ 
	float3 Multiple(matrix3x3 &M1, float3 &V)
{
	float3 M;

	M.x = M1.e[0].x * V.x + M1.e[0].y * V.y + M1.e[0].z * V.z;
	M.y = M1.e[1].x * V.y + M1.e[1].y * V.y + M1.e[1].z * V.z;
	M.z = M1.e[2].x * V.z + M1.e[2].y * V.y + M1.e[2].z * V.z;

	return M;
}

#endif