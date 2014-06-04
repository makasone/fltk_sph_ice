/*! 
  @file rx_shape_matching.cpp
	
  @brief Shape Matching法による弾性変形
 
  @author Makoto Fujisawa
  @date 2013-07
*/


//-----------------------------------------------------------------------------
// インクルードファイル
//-----------------------------------------------------------------------------
#include "ShapeMatching.h"


//-----------------------------------------------------------------------------
// Shape Matchingクラスの実装
//-----------------------------------------------------------------------------
/*!
 * デフォルトコンストラクタ
 */
rxShapeMatching::rxShapeMatching(int obj)
{
	//パラメータの初期化
	m_iObjectNo = obj;

	m_dDt = 0.01;

	m_v3Gravity = Vec3(0.0, -9.81/**0.005*/, 0.0);
	//m_v3Gravity = Vec3(0.0, 0.0, 0.0);

	m_v3Min = Vec3(-1.0);
	m_v3Max = Vec3(1.0);

	m_dAlpha = 1.0;	
	m_dBeta = 1.0;	// これが小さいと硬い材質になる
	
	m_bLinearDeformation = true;
	m_bVolumeConservation = true;

	m_fpCollision = 0;

	//メモリ確保
	//クラスタが保存できる最大粒子数をMAXPARTICLEで定義する．
	m_pOrgPos = new double[MAXPARTICLE*SM_DIM];
	m_pCurPos = new double[MAXPARTICLE*SM_DIM];
	m_pNewPos = new double[MAXPARTICLE*SM_DIM];
	m_pGoalPos = new double[MAXPARTICLE*SM_DIM];
	m_pVel = new double[MAXPARTICLE*SM_DIM];

	m_pMass = new double[MAXPARTICLE];
	m_pFix = new bool[MAXPARTICLE];

	Clear();
}

/*!
 * デストラクタ
 */
rxShapeMatching::~rxShapeMatching()
{
}

/*
 *	GPUの初期化
 */
void rxShapeMatching::InitGPU()
{
	//デバイス側のメモリを確保
	//クラスタが保存できる最大粒子数をMAXPARTICLEで定義する．
	cudaMalloc((void**)&d_OrgPos, sizeof(double) * MAXPARTICLE * SM_DIM);

	//ホスト側のデータをデバイス側へ転送
	//cudaMemcpy(d_OrgPos, m_pOrgPos, sizeof(float) * MAXPARTICLE * SM_DIM, cudaMemcpyHostToDevice);
}


/*!
 * 頂点位置の初期化
 */
void rxShapeMatching::initialize(void)
{
	for(int i = 0; i < m_iNumVertices; ++i){

		for(int j = 0; j < SM_DIM; j++)
		{
			m_pCurPos[i*SM_DIM+j] = m_pOrgPos[i*SM_DIM+j];
			m_pNewPos[i*SM_DIM+j] = m_pOrgPos[i*SM_DIM+j];
			m_pGoalPos[i*SM_DIM+j] = m_pOrgPos[i*SM_DIM+j];
			m_pVel[i*SM_DIM+j] = 0.0;
		}

		m_pFix[i] = false;
	}
}

/*!
 * 全頂点を初期化
 */
void rxShapeMatching::Clear()
{
	m_iNumVertices = 0;

	for(int i = 0; i < MAXPARTICLE; i++)
	{
		for(int j = 0; j < SM_DIM; j++)
		{
			m_pOrgPos[i*SM_DIM+j] = 0.0;
			m_pCurPos[i*SM_DIM+j] = 0.0;
			m_pNewPos[i*SM_DIM+j] = 0.0;
			m_pGoalPos[i*SM_DIM+j] = 0.0;
			m_pVel[i*SM_DIM+j] = 0.0;
		}

		m_pMass[i] = 0.0;
		m_pFix[i] = false;
	}
}


/*!
 * 頂点を追加
 * @param[in] pos 頂点位置
 * @param[out] mass 頂点質量
 */
void rxShapeMatching::AddVertex(const Vec3 &pos, double mass)
{
	for(int i = 0; i < SM_DIM; i++)
	{
		m_pOrgPos[m_iNumVertices*SM_DIM+i] = pos[i];
		m_pCurPos[m_iNumVertices*SM_DIM+i] = pos[i];
		m_pNewPos[m_iNumVertices*SM_DIM+i] = pos[i];
		m_pGoalPos[m_iNumVertices*SM_DIM+i] = pos[i];
		m_pVel[m_iNumVertices*SM_DIM+i] = 0.0;
	}

	//m_pFixは特にすること無し
	m_pMass[m_iNumVertices] = mass;

	m_iNumVertices++;

	initialize();	//これ必要？
}

/*!
 * 外力
 *  - 重力と境界壁からの力の影響
 * @param[in] dt タイムステップ幅
 */
void rxShapeMatching::calExternalForces(double dt)
{
	// 重力の影響を付加
	for(int i = 0; i < m_iNumVertices; ++i){
		if(m_pFix[i]) continue;
		for(int j = 0; j < SM_DIM; j++)
		{
			m_pVel[i*SM_DIM+j] += m_v3Gravity[j]*dt;
			m_pNewPos[i*SM_DIM+j] = m_pCurPos[i*SM_DIM+j]+m_pVel[i*SM_DIM+j]*dt;
			m_pGoalPos[i*SM_DIM+j] = m_pOrgPos[i*SM_DIM+j];
		}
	}

	// 境界壁の影響
	double res = 0.5;	// 反発係数
	for(int i = 0; i < m_iNumVertices; ++i){
		if(m_pFix[i]) continue;

		if(m_pNewPos[i*SM_DIM+0] < m_v3Min[0] || m_pNewPos[i*SM_DIM+0] > m_v3Max[0])
		{
			m_pNewPos[i*SM_DIM+0] = m_pCurPos[i*SM_DIM+0]-m_pVel[i*SM_DIM+0]*dt*res;
			m_pNewPos[i*SM_DIM+1] = m_pCurPos[i*SM_DIM+1];
			m_pNewPos[i*SM_DIM+2] = m_pCurPos[i*SM_DIM+2];
		}
		if(m_pNewPos[i*SM_DIM+1] < m_v3Min[1] || m_pNewPos[i*SM_DIM+1] > m_v3Max[1])
		{
			m_pNewPos[i*SM_DIM+1] = m_pCurPos[i*SM_DIM+1]-m_pVel[i*SM_DIM+1]*dt*res;
			m_pNewPos[i*SM_DIM+0] = m_pCurPos[i*SM_DIM+0];
			m_pNewPos[i*SM_DIM+2] = m_pCurPos[i*SM_DIM+2];
		}
		if(m_pNewPos[i*SM_DIM+2] < m_v3Min[2] || m_pNewPos[i*SM_DIM+2] > m_v3Max[2])
		{
			m_pNewPos[i*SM_DIM+2] = m_pCurPos[i*SM_DIM+2]-m_pVel[i*SM_DIM+2]*dt*res;
			m_pNewPos[i*SM_DIM+0] = m_pCurPos[i*SM_DIM+0];
			m_pNewPos[i*SM_DIM+1] = m_pCurPos[i*SM_DIM+1];
		}

		Vec3 np(m_pNewPos[i*SM_DIM+0], m_pNewPos[i*SM_DIM+1], m_pNewPos[i*SM_DIM+2]);
		clamp(np);

		m_pNewPos[i*SM_DIM+0] = np[0];
		m_pNewPos[i*SM_DIM+1] = np[1];
		m_pNewPos[i*SM_DIM+2] = np[2];
	}
}

void rxShapeMatching::calCollision(double dt)
{
	// 他のオブジェクトとの衝突
	//if(m_fpCollision != 0){
	//	for(int i = 0; i < m_iNumVertices; ++i){
	//		if(m_pFix[i]) continue;
	//		Vec3 &p = m_pCurPos[i];
	//		Vec3 &np = m_pNewPos[i];
	//		Vec3 &v = m_pVel[i];
	//		m_fpCollision(p, np, v, m_iObjectNo);
	//	}
	//}
}

/*!
 * Shape Matching法
 *  - 目標位置を計算して，m_pNewPosをその位置に近づける
 * @param[in] dt タイムステップ幅
 */
void rxShapeMatching::shapeMatching(double dt)
{
	if(m_iNumVertices <= 1) return;

	Vec3 cm(0.0), cm_org(0.0);	// 重心
	double mass = 0.0f;	// 総質量

	// 重心座標の計算
	for(int i = 0; i < m_iNumVertices;++i){
		double m = m_pMass[i];
		if(m_pFix[i]) m *= 300.0;	// 固定点の質量を大きくする
		mass += m;

		Vec3 np(m_pNewPos[i*SM_DIM+0], m_pNewPos[i*SM_DIM+1], m_pNewPos[i*SM_DIM+3]);
		Vec3 op(m_pOrgPos[i*SM_DIM+0], m_pOrgPos[i*SM_DIM+1], m_pOrgPos[i*SM_DIM+3]);

		cm += np*m;
		cm_org += op*m;
	}

	cm /= mass;
	cm_org /= mass;

	rxMatrix3 Apq(0.0), Aqq(0.0);
	Vec3 p, q;

	// Apq = Σmpq^T
	// Aqq = Σmqq^T
	for(int i = 0; i < m_iNumVertices; ++i){

		Vec3 np(m_pNewPos[i*SM_DIM+0], m_pNewPos[i*SM_DIM+1], m_pNewPos[i*SM_DIM+3]);
		Vec3 op(m_pOrgPos[i*SM_DIM+0], m_pOrgPos[i*SM_DIM+1], m_pOrgPos[i*SM_DIM+3]);

		p = np-cm;
		q = op-cm_org;
		double m = m_pMass[i];

		Apq(0,0) += m*p[0]*q[0];
		Apq(0,1) += m*p[0]*q[1];
		Apq(0,2) += m*p[0]*q[2];
		Apq(1,0) += m*p[1]*q[0];
		Apq(1,1) += m*p[1]*q[1];
		Apq(1,2) += m*p[1]*q[2];
		Apq(2,0) += m*p[2]*q[0];
		Apq(2,1) += m*p[2]*q[1];
		Apq(2,2) += m*p[2]*q[2];

		Aqq(0,0) += m*q[0]*q[0];
		Aqq(0,1) += m*q[0]*q[1];
		Aqq(0,2) += m*q[0]*q[2];
		Aqq(1,0) += m*q[1]*q[0];
		Aqq(1,1) += m*q[1]*q[1];
		Aqq(1,2) += m*q[1]*q[2];
		Aqq(2,0) += m*q[2]*q[0];
		Aqq(2,1) += m*q[2]*q[1];
		Aqq(2,2) += m*q[2]*q[2];
	}

	rxMatrix3 R, S;
	PolarDecomposition(Apq, R, S);

	if(m_bLinearDeformation){
		// Linear Deformations
		rxMatrix3 A;
		A = Apq*Aqq.Inverse();	// A = Apq*Aqq^-1

		// 体積保存のために√(det(A))で割る
		if(m_bVolumeConservation){
			double det = fabs(A.Determinant());
			if(det > RX_FEQ_EPS){
				det = 1.0/sqrt(det);
				if(det > 2.0) det = 2.0;
				A *= det;
			}
		}

		// 回転行列Rの代わりの行列RL=βA+(1-β)Rを計算
		rxMatrix3 RL = m_dBeta*A+(1.0-m_dBeta)*R;

		// 目標座標を計算し，現在の頂点座標を移動
		for(int i = 0; i < m_iNumVertices; ++i){
			if(m_pFix[i]) continue;

			Vec3 np(m_pNewPos[i*SM_DIM+0], m_pNewPos[i*SM_DIM+1], m_pNewPos[i*SM_DIM+3]);
			Vec3 op(m_pOrgPos[i*SM_DIM+0], m_pOrgPos[i*SM_DIM+1], m_pOrgPos[i*SM_DIM+3]);

			q = op-cm_org;
			Vec3 gp = RL*q+cm;

			for(int j = 0; j < SM_DIM; j++)
			{
				m_pGoalPos[i*SM_DIM+j] = gp[j];
				m_pNewPos[i*SM_DIM+j] += (m_pGoalPos[i*SM_DIM+j]-np[j])*m_dAlpha;
			}
		}
	}
	else{
		//// Quadratic Deformations
		//double Atpq[3][9];
		//for(int j = 0; j < 9; ++j){
		//	Atpq[0][j] = 0.0;
		//	Atpq[1][j] = 0.0;
		//	Atpq[2][j] = 0.0;
		//}
		//rxMatrixN<double,9> Atqq;
		//Atqq.SetValue(0.0);

		//for(int i = 0; i < m_iNumVertices; ++i){
		//	p = m_pNewPos[i]-cm;
		//	q = m_pOrgPos[i]-cm_org;

		//	// q~の計算
		//	double qt[9];
		//	qt[0] = q[0];      qt[1] = q[1];      qt[2] = q[2];
		//	qt[3] = q[0]*q[0]; qt[4] = q[1]*q[1]; qt[5] = q[2]*q[2];
		//	qt[6] = q[0]*q[1]; qt[7] = q[1]*q[2]; qt[8] = q[2]*q[0];

		//	// A~pq = Σmpq~ の計算
		//	double m = m_vMass[i];
		//	for(int j = 0; j < 9; ++j){
		//		Atpq[0][j] += m*p[0]*qt[j];
		//		Atpq[1][j] += m*p[1]*qt[j];
		//		Atpq[2][j] += m*p[2]*qt[j];
		//	}

		//	// A~qq = Σmq~q~ の計算
		//	for(int j = 0; j < 9; ++j){
		//		for(int k = 0; k < 9; ++k){
		//			Atqq(j,k) += m*qt[j]*qt[k];
		//		}
		//	}
		//}

		//// A~qqの逆行列算出
		//Atqq.Invert();

		//double At[3][9];
		//for(int i = 0; i < 3; ++i){
		//	for(int j = 0; j < 9; j++){
		//		At[i][j] = 0.0f;
		//		for(int k = 0; k < 9; k++){
		//			At[i][j] += Atpq[i][k]*Atqq(k,j);
		//		}

		//		// βA~+(1-β)R~
		//		At[i][j] *= m_dBeta;
		//		if(j < 3){
		//			At[i][j] += (1.0f-m_dBeta)*R(i,j);
		//		}
		//	}
		//}

		// // a00a11a22+a10a21a02+a20a01a12-a00a21a12-a20a11a02-a10a01a22
		//double det = At[0][0]*At[1][1]*At[2][2]+At[1][0]*At[2][1]*At[0][2]+At[2][0]*At[0][1]*At[1][2]
		//			-At[0][0]*At[2][1]*At[1][2]-At[2][0]*At[1][1]*At[0][2]-At[1][0]*At[0][1]*At[2][2];

		//// 体積保存のために√(det(A))で割る
		//if(m_bVolumeConservation){
		//	if(det != 0.0f){
		//		det = 1.0f/sqrt(fabs(det));
		//		if(det > 2.0f) det = 2.0f;
		//		for(int i = 0; i < 3; ++i){
		//			for(int j = 0; j < 3; ++j){
		//				At[i][j] *= det;
		//			}
		//		}
		//	}
		//}


		//// 目標座標を計算し，現在の頂点座標を移動
		//for(int i = 0; i < m_iNumVertices; ++i){
		//	if(m_pFix[i]) continue;
		//	q = m_pOrgPos[i]-cm_org;

		//	for(int k = 0; k < 3; ++k){
		//		m_pGoalPos[i][k] = At[k][0]*q[0]+At[k][1]*q[1]+At[k][2]*q[2]
		//						  +At[k][3]*q[0]*q[0]+At[k][4]*q[1]*q[1]+At[k][5]*q[2]*q[2]+
		//						  +At[k][6]*q[0]*q[1]+At[k][7]*q[1]*q[2]+At[k][8]*q[2]*q[0];
		//	}

		//	m_pGoalPos[i] += cm;
		//	m_pNewPos[i] += (m_pGoalPos[i]-m_pNewPos[i])*m_dAlpha;
		//}

	}
}

/*!
 * 速度と位置の更新
 *  - 新しい位置と現在の位置座標から速度を算出
 * @param[in] dt タイムステップ幅
 */
void rxShapeMatching::integrate(double dt)
{
	double dt1 = 1.0f/dt;
	for(int i = 0; i < m_iNumVertices; ++i){
		for(int j = 0; j < SM_DIM; j++)
		{
			m_pVel[i*SM_DIM+j] = (m_pNewPos[i*SM_DIM+j]-m_pCurPos[i*SM_DIM+j])*dt1;
			m_pCurPos[i*SM_DIM+j] = m_pNewPos[i*SM_DIM+j];
		}
	}
}

/*!
 * シミュレーションステップを進める
 */
void rxShapeMatching::Update()
{	
	calExternalForces(m_dDt);
	calCollision(m_dDt);
	shapeMatching(m_dDt);
	integrate(m_dDt);
}

/*!
 * 固定頂点を設定
 * @param[in] i 頂点インデックス
 * @param[in] pos 固定位置
 */
void rxShapeMatching::FixVertex(int i, const Vec3 &pos)
{
	m_pNewPos[i*SM_DIM+0] = pos[0];
	m_pNewPos[i*SM_DIM+1] = pos[1];
	m_pNewPos[i*SM_DIM+2] = pos[2];

	m_pFix[i] = true;
}

/*!
 * 頂点の固定を解除
 * @param[in] i 頂点インデックス
 */
void rxShapeMatching::UnFixVertex(int i)
{
	m_pFix[i] = false;
}


