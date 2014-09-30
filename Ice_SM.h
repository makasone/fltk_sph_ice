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
extern void LaunchShapeMatchingGPU
(
	int prtNum,
	float* prtPos,
	float* prtVel, 
	float* orgPos,
	float* curPos,
	float* orgCm,
	float* curCm,
	float* vel,
	int* pIndxes, 
	int* d_IndxSet,
	float dt
);

extern void LaunchShapeMatchingIterationGPU
(
	int prtNum,
	float* prtPos, 
	float* prtVel,
	float* sldPos,
	float* sldVel, 
	float* orgPos,
	float* curPos,
	float* orgCm,
	float* curCm,
	float* vel,
	int* pIndxes, 
	int* d_IndxSet,
	float dt
);

class Ice_SM : public rxShapeMatching
{
protected:

	unsigned m_iIndxNum;						//配列で実装したため、穴あきに対応するための最大添字番号

	float* m_fpAlphas;							//!< stiffnessパラメータ[0,1] (速度計算に使用)
	float* m_fpBetas;							//!< deformationパラメータ[0,1]	未使用
	
	int* m_ipLayeres;

	static int s_iIterationNum;					//反復回数

	Vec3 m_vec3OrgCm;							//初期のクラスタの重心
	Vec3 m_vec3NowCm;							//現在のクラスタの重心
	//vector<Vec3> m_vvec3OrgQ;					//初期の位置-重心

	Vec3 m_vec3DisCm;							//重心位置の変位

	rxMatrix3	m_mtrx3Apq;						//変形行列Apq
	rxMatrix3	m_mtrx3AqqInv;					//変形行列Aqqの逆行列	前計算可能

	vector<int> m_iLinearDeformation;			//!< Linear/Quadratic deformation切り替えフラグ　未使用
	vector<int> m_iVolumeConservation;			//!< 変形時の体積保存性(√det(A)で割るかどうか)　未使用

	static const float* s_pfPrtPos;				//位置のホストポインタ　読み込み専用
	static const float* s_pfPrtVel;				//速度のホストポインタ　読み込み専用

	static float* s_pfSldPos;					//クラスタの最終的な位置
	static float* s_pfSldVel;					//クラスタの最終的な速度

//--------------------------------------GPU------------------------------------------------------------
	static float* sd_PrtPos;					//最終粒子位置のデバイスポインタ
	static float* sd_PrtVel;					//最終粒子速度のデバイスポインタ

	static float* d_OrgPos;
	static float* d_CurPos;
	
	static float* d_OrgCm;						//初期重心
	static float* d_CurCm;						//現在の重心　まだ使わない

	static float* d_Mass;
	static float* d_Vel;

	static bool* d_Fix;

	static int* d_PIndxes;
	static int* d_IndxSet;						//クラスタのデータの開始添字と終了添字を保存

	static int s_vertNum;						//全てのクラスタが対象とする粒子数
	static int s_vertSum;						//全クラスタに含まれる粒子の総数

//--------------------------------------GPU------------------------------------------------------------

		//bounds壁にぶつかったときの反発力
		//allowFlip反転を許すかのフラグ

public:
	Ice_SM(int obj);
	~Ice_SM();
	
	static void Ice_SM::InitGPU(const vector<Ice_SM*>& sm, float* d_pos, float* d_vel, int prtNum);

	static void SetPrtPointerPosAndVel(const float* pos, const float* vel){	s_pfPrtPos = pos;	s_pfPrtVel = vel;	}
	static void SetDevicePosPointer(float* d_pos){	sd_PrtPos = d_pos;	}
	
	static float* GetSldPosPointer(){	return s_pfSldPos;	}
	static float* GetSldVelPointer(){	return s_pfSldVel;	}

	static float* GetDeviceSPHPosPointer(){	return sd_PrtPos;	}
	static float* GetDeviceSPHVelPointer(){	return sd_PrtVel;	}

	static float* GetDevicePosPointer(){	return d_CurPos;	}
	static float* GetDeviceVelPointer(){	return d_Vel;		}
	static float* GetOrgPosPointer(){		return d_OrgPos;	}
	static float* GetOrgCmPointer(){		return d_OrgCm;		}

	static int* GetDeviceIndexSetPointer(){	return d_IndxSet;	}
	static int	GetVertexNum(){				return s_vertNum;	}

	void InitGPU_Instance();
	static void InitFinalParamPointer(int vrtxNum);

	void AddVertex(const Vec3 &pos, double mass, int pIndx);
	
	void UpdateCPU();
	void UpdateUsePathCPU();

	static void UpdateGPU();
	static void UpdateUsePathGPU();
	static void UpdateIterationGPU(float* sldPos, float* sldVel);

	static void CalcAverage();
	void CopyDeviceToInstance(int num);

	void ShapeMatchingUsePath();
	void ShapeMatchingSolid();
	void ShapeMatchingIteration();

	void calExternalForces();
	void calExternalForcesIteration();
	
	void integrate(double dt);
	void integrateIteration();

	void SetAlphas(int indx, float alpha){	m_fpAlphas[indx] = alpha;	}
	void SetBetas (int indx, float beta){	m_fpBetas[indx] = beta;		}

	void SetLayer(int indx, int layer){	m_ipLayeres[indx] = layer;	}

	static void SetIterationNum(int itr){	s_iIterationNum = itr;	}

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
	static int GetIteration(){	return s_iIterationNum;	}

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