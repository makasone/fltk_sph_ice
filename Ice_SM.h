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

#include <math.h>
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

extern void LaunchShapeMatchingUsePathGPU
(
	int prtNum,
	float* prtPos,
	float* prtVel, 
	float* orgPos,
	float* curPos,
	float* orgCm,
	float* curCm,
	float* curApq,
	float* vel,
	int* pIndxes, 
	int* d_IndxSet,
	float dt
);

class Ice_SM : public rxShapeMatching
{
public:
	//領域分割で指定するラベル　XYZ
	//TODO: いずれは2進数で扱えるように
	enum PrtErea{
		TOP_LEFT__FORE = 122,
		TOP_RIGHT_FORE = 222,
		TOP_LEFT__BACK = 221,
		TOP_RIGHT_BACK = 121,
		BOT_LEFT__FORE = 112,
		BOT_RIGHT_FORE = 212,
		BOT_LEFT__BACK = 111,
		BOT_RIGHT_BACK = 211,

		NUM
	};

	//ビットでなくintで領域分割ラベルを指定しているので，アクセスしやすい構造体を用意
	struct EreaData{
		int x;
		int y;
		int z;

		EreaData(){};
		EreaData(int data){
			EreaData temp = IntToEreaData(data);
			x = temp.x; y = temp.y; z = temp.z;
		}
		EreaData(int X, int Y, int Z){
			if(X < 0 || Y < 0 || Z < 0){
				cout << __FUNCTION__ << " ERROR" << endl;
			}
			x = X; y = Y; z = Z;
		}
	};

protected:

	float* m_pPrePos;							//前フレームの頂点位置

	float* m_fpAlphas;							//!< stiffnessパラメータ[0,1] (速度計算に使用)
	float* m_fpBetas;							//!< deformationパラメータ[0,1]	未使用
	
	int* m_ipLayeres;
	EreaData* m_ipErea;

	vector<int> m_ivNeighborFeatureCluster;		//周辺に存在する特徴クラスタ（運動計算をするクラスタ）の近傍情報を格納
												//自身が特徴クラスタの場合，少しい遠い特徴クラスタの情報を持つことになる
												

	static int s_iIterationNum;					//反復回数
	static float s_fItrStiffness;				//変形量の閾値

	Vec3 m_vec3OrgCm;							//初期のクラスタの重心
	Vec3 m_vec3NowCm;							//現在のクラスタの重心
	Vec3 m_vec3PreCm;							//前フレームのクラスタの重心
	//vector<Vec3> m_vvec3OrgQ;					//初期の位置-重心

	Vec3 m_vec3DisCm;							//重心位置の変位
	float m_fDefAmount;							//変形量
	float m_fDefPriority;						//変形優先度合い　変形量を周りのクラスタと比べた時，どのくらい大きいかの指標
												//平均を上回っている場合プラス，下回っている場合マイナスとなる

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

	static float* d_Apq;

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
	Ice_SM();
	Ice_SM(int obj);
	Ice_SM(int obj, int prtNum);
	Ice_SM(const Ice_SM& copy);

	~Ice_SM();
	
	//演算子のオーバーロード
	Ice_SM& operator=(const Ice_SM& copy);

	void ReleaseMemory();
	void Copy(const Ice_SM& copy);

	//アクセッサ---------------------------------------------------------
	static void Ice_SM::InitGPU(const vector<Ice_SM*>& sm, float* d_pos, float* d_vel, int prtNum, int maxprtNum);

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
	static float* GetCurCmPointer(){		return d_CurCm;		}
	static float* GetDeviceApqPointer(){	return d_Apq;		}

	static int* GetDeviceIndexesPointer(){	return d_PIndxes;	}
	static int* GetDeviceIndexSetPointer(){	return d_IndxSet;	}
	//-------------------------------------------------------------------

	void InitGPU_Instance();
	static void InitFinalParamPointer(int vrtxNum);
	static void ResetFinalParamPointer(unsigned clusterNum);
	
	void UpdateCPU();
	void UpdateUsePathCPU();

	void UpdatePrePos(int pIndx);

	static void UpdateGPU();
	static void UpdateUsePathGPU();
	static void UpdateIterationGPU(float* sldPos, float* sldVel);

	static void CalcAverage();
	void CopyDeviceToInstance(int num);

	void ShapeMatchingUsePath();
	void ShapeMatchingSolid();
	void ShapeMatchingIteration();
	void ShapeMatchingSelected(int selected);

	void calExternalForces();
	void calExternalForcesIteration();

	void CalcMass();
	void CalcOrgCm();

	void ClassifyAllOrgParticle();
	void ClassifyOrgParticle(int indx);

	void integrate(double dt);
	void integrateIteration();

	int AddAnotherClusterVertex(const Vec3& orgPos, const Vec3& curPos, const Vec3& vel, double mass, int pIndx, double alpha, double beta, int layer);

	void AddVertex(const Vec3 &pos, double mass, int pIndx);
	void AddVertex(const Vec3 &pos, const Vec3& vel, double mass, int pIndx);
	void Remove(int indx);
	void Clear();

	static unsigned CalcEreaDistance(const EreaData& ereaA, const EreaData& ereaB);
	static EreaData IntToEreaData(int erea);
	static int EreaDataToInt(EreaData erea);

//アクセッサ
	void SetAlphas(int indx, float alpha){	m_fpAlphas[indx] = alpha;	}
	void SetBetas (int indx, float beta){	m_fpBetas[indx] = beta;		}

	void SetLayer(int indx, int layer){	m_ipLayeres[indx] = layer;	}

	static void SetIterationNum(int itr){	s_iIterationNum = itr;	}
	static void SetItrStiffness(float stiff){	s_fItrStiffness = stiff;	}

	void SetNowCm(Vec3 nowCm){	m_vec3NowCm = nowCm;	}
	void SetPreCm(Vec3 preCm){	m_vec3PreCm = preCm;	}
	void SetApq(rxMatrix3& Apq){	m_mtrx3Apq = Apq;	}
	void SetPrePos(int pIndx, const Vec3& nowPos);

	void SetDefAmount(float amount){	m_fDefAmount = amount;	}
	void SetDefPriority(float priority){	m_fDefPriority = priority;	}

	void SetLinerFalg(int indx, int flag){	m_iLinearDeformation[indx] = flag; }
	void SetVolumeFlag(int indx, int flag){	m_iVolumeConservation[indx] = flag;}

	float GetAlphas(int indx) const {	return m_fpAlphas[indx];	}
	float GetBetas (int indx) const {	return m_fpBetas[indx];	}

	Vec3 GetCm(void){		return m_vec3NowCm;	}
	const Vec3& GetOrgCm(void) const {	return m_vec3OrgCm;	}
	Vec3 GetPreCm(void){	return m_vec3PreCm;	}
	Vec3 GetPrePos(int pIndx) const {	return Vec3(m_pPrePos[pIndx*SM_DIM+0],m_pPrePos[pIndx*SM_DIM+1],m_pPrePos[pIndx*SM_DIM+2]);	}

	rxMatrix3 GetApq() const {	return m_mtrx3Apq;	}

	float GetDefAmount() const {	return m_fDefAmount;	}
	float GetDefPriority(){	return m_fDefPriority;	}

	int GetLayer(int indx) const {	return m_ipLayeres[indx];	}
	const EreaData& erea(int indx) const { return m_ipErea[indx];	}

	static int GetIteration(){		return s_iIterationNum;	}
	static float GetItrStiffness(){	return s_fItrStiffness;	}

	int neighborFeatureCluster(int indx) const {	return m_ivNeighborFeatureCluster[indx];	}
	int neighborFeatureClusterNum() const {	return m_ivNeighborFeatureCluster.size();	}
	void SetNeighborFeatureCluster(int indx, int cIndx){	m_ivNeighborFeatureCluster[indx] = cIndx;	}
	void AddNeighborFeatureCluster(int cIndx) { m_ivNeighborFeatureCluster.push_back(cIndx);	}
	void RemoveNeibhborFeatureCluster(int cIndx)
	{ m_ivNeighborFeatureCluster.erase(remove(m_ivNeighborFeatureCluster.begin(), m_ivNeighborFeatureCluster.end(), cIndx), m_ivNeighborFeatureCluster.end());	}

	void ClearNeighborFeaturceCluster(){
			m_ivNeighborFeatureCluster.clear();
			m_ivNeighborFeatureCluster.shrink_to_fit();
	}

	void OrganizeNeighborFeatureCluster(){
		sort(m_ivNeighborFeatureCluster.begin(), m_ivNeighborFeatureCluster.end());
		m_ivNeighborFeatureCluster.erase( 
			unique(m_ivNeighborFeatureCluster.begin(), m_ivNeighborFeatureCluster.end()),
			m_ivNeighborFeatureCluster.end()
		);
		m_ivNeighborFeatureCluster.shrink_to_fit();
	}

//いずれ使う？
	int GetLinerFlags(int indx){	return m_iLinearDeformation[indx]; }
	int GetVolumeFlags(int indx){	return m_iVolumeConservation[indx]; }

	bool CheckIndx(int pIndx) const;
	int	 SearchIndx(int pIndx) const;

	const Vec3& GetDisVec(){	return m_vec3DisCm;	}
	void CalcDisplaceMentVectorCm();

//デバッグ
	void DebugIndx(void);
	void DebugLayer(void);

	void DebugClusterInfo();
	void DebugNeighborFeatureClusterInfo();
};

#endif