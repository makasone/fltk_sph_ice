/*!
  @file Ice_SM.h
	
  @brief ShapeMatching法を基にした相変化シミュレーション
  @ref M. Muller et al., "Meshless deformations based on shape matching", SIGGRAPH2005. 
 
  @author Ryo Nakasone
  @date 2013-10
*/

#ifndef _ICE_SM_
#define _ICE_SM_

#include "ShapeMatching.h"
#include "QueryCounter.h"

#include <time.h>

using namespace std;

//GPU処理
extern void LaunchShapeMatchingGPU(
	float* prtPos,
	float* prtVel, 
	float* orgPos,
	float* curPos,
	float* vel,
	int* pIndxes, 
	float dt,
	int prtNum
	);


class Ice_SM : public rxShapeMatching
{
protected:

	unsigned m_iIndxNum;				//配列で実装したため、穴あきに対応するための最大添字番号

	float* m_fpAlphas;					//!< stiffnessパラメータ[0,1] (速度計算に使用)
	float* m_fpBetas;					//!< deformationパラメータ[0,1]	未使用
	
	int* m_ipLayeres;

	Vec3 m_vec3OrgCm;					//初期のクラスタの重心
	Vec3 m_vec3NowCm;					//現在のクラスタの重心
	//vector<Vec3> m_vvec3OrgQ;			//初期の位置-重心

	Vec3 m_vec3DisCm;					//重心位置の変位

	rxMatrix3	m_mtrx3Apq;				//変形行列Apq
	rxMatrix3	m_mtrx3AqqInv;			//変形行列Aqqの逆行列	前計算可能

	vector<int> m_iLinearDeformation;	//!< Linear/Quadratic deformation切り替えフラグ　未使用
	vector<int> m_iVolumeConservation;	//!< 変形時の体積保存性(√det(A)で割るかどうか)　未使用

	static const float* s_pfPrtPos;		//読み込み専用
	static const float* s_pfPrtVel;		//読み込み専用

	static float* sd_PrtPos;		//デバイスポインタ
	static float* sd_PrtVel;		//デバイスポインタ

//--------------------------------------GPU------------------------------------------------------------
	//デバイス側へのポインタ
	static float* d_OrgPos;
	static float* d_CurPos;
	static float* d_NewPos;
	static float* d_GoalPos;
	static float* d_Mass;
	static float* d_Vel;

	static bool* d_Fix;

	static int* d_PIndxes;

	static int* d_IndxSet;					//クラスタのデータの開始添字と終了添字を保存

	static int s_vertSum;					//全クラスタに含まれる粒子の総数

//--------------------------------------GPU------------------------------------------------------------

		//bounds壁にぶつかったときの反発力
		//allowFlip反転を許すかのフラグ

public:
	Ice_SM(int obj);
	~Ice_SM();

	static void SetParticlePosAndVel(const float* pos, const float* vel)
	{
		s_pfPrtPos = pos;	s_pfPrtVel = vel;
	}

	static void Ice_SM::InitGPU(const vector<Ice_SM*>& sm, float* d_pos, float* d_vel);

	void InitGPU_Instance();

	void AddVertex(const Vec3 &pos, double mass, int pIndx);
	
	void Update();
	static void UpdateGPU();

	void CopyDeviceToInstance();

	void ShapeMatching(float* newPos, double dt);
	void ShapeMatchingSolid(float* newPos, double dt);

	void calExternalForces(float* newPos, double dt);
	void integrate(float* newPos, double dt);

	void SetAlphas(int indx, float alpha){	m_fpAlphas[indx] = alpha;	}
	void SetBetas (int indx, float beta){	m_fpBetas[indx] = beta;		}

	void SetLayer(int indx, int layer){	m_ipLayeres[indx] = layer;	}

	void SetNowCm(Vec3& nowCm){	m_vec3NowCm = nowCm;	}
	void SetApq(rxMatrix3& Apq){	m_mtrx3Apq = Apq;	}

	void SetLinerFalg(int indx, int flag){	m_iLinearDeformation[indx] = flag; }
	void SetVolumeFlag(int indx, int flag){	m_iVolumeConservation[indx] = flag;}

	float GetAlphas(int indx){	return m_fpAlphas[indx];	}
	float GetBetas (int indx){	return m_fpBetas[indx];	}

	Vec3 GetCm(void){		return m_vec3NowCm;	}
	Vec3 GetOrgCm(void){	return m_vec3OrgCm;	}

	rxMatrix3 GetApq(void){	return m_mtrx3Apq;	}

	int GetLayer(int indx){	return m_ipLayeres[indx];	}

	void Remove(int indx);
	void Clear();

	//いずれ使う？
	int GetLinerFlags(int indx){	return m_iLinearDeformation[indx]; }
	int GetVolumeFlags(int indx){	return m_iVolumeConservation[indx]; }

	bool CheckHole(int oIndx);
	bool CheckIndx(int pIndx);
	int	 SearchIndx(int pIndx);

	const Vec3& GetDisVec(){	return m_vec3DisCm;	}
	void CalcDisplaceMentVectorCm();

	//デバッグ
	void DebugIndx(void);
	void DebugLayer(void);

};

#endif