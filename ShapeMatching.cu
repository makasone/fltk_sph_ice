#include <stdio.h>
#include <math.h>
#include <rx_cu_common.cuh>	//使わない？

#define SM_DIM 3

//引き渡す変数名が間違っていてもエラーが出ないので注意

void LaunchShapeMathcingGPU(float* curPos, float* vel, int* pIndxes, float dt, int prtNum);
__global__ void Update(float* curPos, float* vel, int* pIndxes, float dt, int prtNum);
__device__ void ExternalForce(float* curPos, float* vel, int* pIndxes, float dt, int prtNum);
__device__ void ProjectPos(float* curPos, float* vel, int* pIndxes, float dt, int prtNum);
__device__ void Integrate(float* curPos, float* vel, int* pIndxes, float dt, int prtNum);

float* d_Pos;	//デバイス側の粒子位置
float* d_Vel;	//デバイス側の粒子速度

//パラメータの初期化
void InitParam()
{
}

//GPU処理
void LaunchShapeMatchingGPU(float* curPos, float* vel, int* pIndxes, float dt, int prtNum)
{
	
	//printf("LaunchGPUKernel");

	dim3 grid(1, 1);
	dim3 block(1, 1, 1);

	//運動計算
	Update <<< grid , block >>> (curPos, vel, pIndxes, dt, prtNum);
}


//GPUの位置・速度更新
__global__
void Update(
	float* curPos, 
	float* vel,
	int* pIndxes,
	float dt,
	int prtNum)
{
	//printf("d_Integrate\n");	//めちゃくちゃ出るので注意

	ExternalForce(curPos, vel, pIndxes, dt, prtNum);
	ProjectPos(curPos, vel, pIndxes, dt, prtNum);
	Integrate(curPos, vel, pIndxes, dt, prtNum);
}

__device__
	void ExternalForce(float* curPos, float* vel, int* pIndxes, float dt, int prtNum)
{
	// 重力の影響を付加，速度を反映
	for(int i = 0; i < prtNum; ++i)
	{
		int pIndx = pIndxes[i]*4;
		int cIndx = i*SM_DIM;

		for(int j = 0; j < SM_DIM; j++)
		{
			int jpIndx = pIndx+j;
			int jcIndx = cIndx+j;

			//curPos[jcIndx] = d_Pos[jpIndx]+d_Vel[jpIndx]*dt;
		}
	}

	// 境界壁の影響
	//処理がかなり重くなるが，安定はするみたい
	//double res = 0.9;	// 反発係数
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
	void ProjectPos(float* curPos, float* vel, int* pIndxes, float dt, int prtNum)
{
	if(prtNum <= 1) return;

	float3 cm = make_float3(0.0, 0.0, 0.0);
	float3 cm_org = make_float3(0.0, 0.0, 0.0);	// 重心

	double mass = 0.0;	// 総質量

	// 重心座標の計算
	for(int i = 0; i < prtNum;++i){
		//double m = m_pMass[i];
		double m = 1.0;
		//if(m_pFix[i]) m *= 300.0;	// 固定点の質量を大きくする
		
		int cIndx = i*SM_DIM;
		mass += m;

		cm.x += curPos[cIndx+0]*m;
		cm.y += curPos[cIndx+1]*m;
		cm.z += curPos[cIndx+2]*m;
	}

	cm.x /= mass;
	cm.y /= mass;
	cm.z /= mass;

	//cm_org = m_vec3OrgCm;

	//rxMatrix3 Apq(0.0), Aqq(0.0);
	//Vec3 p, q;

	//// Apq = Σmpq^T
	//// Aqq = Σmqq^T
	//for(int i = 0; i < prtNum; ++i)
	//{
	//	int cIndx = i*SM_DIM;

	//	for(int j = 0; j < SM_DIM; j++)
	//	{
	//		//p[j] = m_pNewPos[cIndx+j]-cm[j];

	//		//p[j] = newPos[cIndx+j]-cm[j];
	//		p[j] = curPos[cIndx+j]-cm[j];
	//		q[j] = m_pOrgPos[cIndx+j]-cm_org[j];
	//	}

	//	double m = m_pMass[i];

	//	Apq(0,0) += m*p[0]*q[0];
	//	Apq(0,1) += m*p[0]*q[1];
	//	Apq(0,2) += m*p[0]*q[2];
	//	Apq(1,0) += m*p[1]*q[0];
	//	Apq(1,1) += m*p[1]*q[1];
	//	Apq(1,2) += m*p[1]*q[2];
	//	Apq(2,0) += m*p[2]*q[0];
	//	Apq(2,1) += m*p[2]*q[1];
	//	Apq(2,2) += m*p[2]*q[2];

	//	Aqq(0,0) += m*q[0]*q[0];
	//	Aqq(0,1) += m*q[0]*q[1];
	//	Aqq(0,2) += m*q[0]*q[2];
	//	Aqq(1,0) += m*q[1]*q[0];
	//	Aqq(1,1) += m*q[1]*q[1];
	//	Aqq(1,2) += m*q[1]*q[2];
	//	Aqq(2,0) += m*q[2]*q[0];
	//	Aqq(2,1) += m*q[2]*q[1];
	//	Aqq(2,2) += m*q[2]*q[2];
	//}

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
	//	//double tmp;
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

	//rxMatrix3 R, S;
	////PolarDecomposition(Apq, R, S, m_mtrxBeforeU);
	//PolarDecomposition(Apq, R, S);

	////double end1 = qc.End()/*/100*/;

	//if(m_bLinearDeformation)
	//{
	//	// Linear Deformations
	//	rxMatrix3 A;
	//	A = Apq*Aqq.Inverse();	// A = Apq*Aqq^-1

	//	// 体積保存のために√(det(A))で割る
	//	if(m_bVolumeConservation){
	//		double det = fabs(A.Determinant());
	//		if(det > RX_FEQ_EPS){
	//			det = 1.0/sqrt(det);
	//			if(det > 2.0) det = 2.0;
	//			A *= det;
	//		}
	//	}

	//	//cout << "計測開始2" << endl;

	//	// 目標座標を計算し，現在の頂点座標を移動
	//	for(int i = 0; i < m_iNumVertices; ++i){
	//		if(m_pFix[i]) continue;

	//		int cIndx = i*SM_DIM;

	//		// 回転行列Rの代わりの行列RL=βA+(1-β)Rを計算
	//		for(int j = 0; j < SM_DIM; j++)
	//		{
	//			q[j] = m_pOrgPos[cIndx+j]-cm_org[j];
	//		}

	//		Vec3 gp(R*q+cm);

	//		for(int j = 0; j < SM_DIM; j++)
	//		{
	//			int jcIndx = cIndx+j;

	//			m_pCurPos[jcIndx] += (gp[j]-m_pCurPos[jcIndx])*m_dAlphas[i];
	//		}
	//	}
	//}
}

/*!
 * 速度と位置の更新
 *  - 新しい位置と現在の位置座標から速度を算出
 * @param[in] dt タイムステップ幅
 */
__device__
	void Integrate(float* curPos, float* vel, int* pIndxes, float dt, int prtNum)
{
	double dt1 = 1.0/dt;

	for(int i = 0; i < prtNum; ++i)
	{
		int pIndx = pIndxes[i]*4;

		for(int j = 0; j < SM_DIM; j++)
		{
			int cIndx = i*SM_DIM+j;

			//vel[cIndx] = (curPos[cIndx] - d_Pos[pIndx+j]) * dt1;/*+ m_v3Gravity * dt * 1.0*/;
		}
	}
	
}

void clamp(float* pos, int cIndx)
{
	//if(pos[cIndx+0] < m_v3Min[0]) pos[cIndx+0] = m_v3Min[0];
	//if(pos[cIndx+0] > m_v3Max[0]) pos[cIndx+0] = m_v3Max[0];
	//if(pos[cIndx+1] < m_v3Min[1]) pos[cIndx+1] = m_v3Min[1];
	//if(pos[cIndx+1] > m_v3Max[1]) pos[cIndx+1] = m_v3Max[1];
	//if(pos[cIndx+2] < m_v3Min[2]) pos[cIndx+2] = m_v3Min[2];
	//if(pos[cIndx+2] > m_v3Max[2]) pos[cIndx+2] = m_v3Max[2];
}