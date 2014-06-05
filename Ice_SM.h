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

//

class Ice_SM : public rxShapeMatching
{
protected:
	vector<int> m_iParticleIndxes;		//!< クラスタに所属する粒子の番号

	vector<double> m_dAlphas;			//!< stiffnessパラメータ[0,1] (速度計算に使用)
	vector<double> m_dBetas;			//!< deformationパラメータ[0,1]

	Vec3 m_vec3OrgCm;					//初期のクラスタの重心
	Vec3 m_vec3NowCm;					//現在のクラスタの重心
	//vector<Vec3> m_vvec3OrgQ;			//初期の位置-重心

	Vec3 m_vec3DisCm;					//重心位置の変位

	rxMatrix3	m_mtrx3Apq;				//変形行列Apq
	rxMatrix3	m_mtrx3AqqInv;			//変形行列Aqqの逆行列	前計算可能

	vector<int> m_iLayeres;

	vector<int> m_iLinearDeformation;	//!< Linear/Quadratic deformation切り替えフラグ
	vector<int> m_iVolumeConservation;	//!< 変形時の体積保存性(√det(A)で割るかどうか)

	static const float* s_pfPrtPos;		//読み込み専用
	static const float* s_pfPrtVel;		//読み込み専用

	//ここに辺と面を持たせる
		//bounds壁にぶつかったときの反発力
		//allowFlip反転を許すかのフラグ　いらない？

public:
	Ice_SM(int obj);
	~Ice_SM();

	static void SetParticlePosAndVel(const float* pos, const float* vel)
	{
		s_pfPrtPos = pos;	s_pfPrtVel = vel;
	}

	void AddVertex(const Vec3 &pos, double mass, int pIndx);
	
	void Update();
	void ShapeMatching(double dt);
	void calExternalForces(float* newPos, double dt);
	void integrate(float* newPos, double dt);

	void ShapeMatchingSolid(float* newPos, double dt);
	void calExternalForcesSolid(double dt);
	void integrateSolid(double);


	void SetAlphas(int indx, int alpha){ m_dAlphas[indx] = alpha;	}
	void SetBetas (int indx, int beta){	m_dBetas[indx] = beta;		}

	void SetLayer(int indx, int layer){	m_iLayeres[indx] = layer;	}

	void SetNowCm(Vec3 nowCm){	m_vec3NowCm = nowCm;	}
	void SetApq(rxMatrix3 Apq){	m_mtrx3Apq = Apq;	}

	void SetLinerFalg(int indx, int flag){	m_iLinearDeformation[indx] = flag; }
	void SetVolumeFlag(int indx, int flag){	m_iVolumeConservation[indx] = flag;}

	int GetParticleIndx(int indx){ return m_iParticleIndxes[indx]; }
	const vector<int>& GetVertexIndxList(){ return m_iParticleIndxes;	}

	double GetAlphas(int indx){	return m_dAlphas[indx];	}
	double GetBetas (int indx){	return m_dBetas[indx];	}

	Vec3 GetCm(void){		return m_vec3NowCm;	}
	Vec3 GetOrgCm(void){	return m_vec3OrgCm;	}

	rxMatrix3 GetApq(void){	return m_mtrx3Apq;	}

	int GetLayer(int indx){	return m_iLayeres[indx];	}

	void Remove(int indx);
	void Clear();

	//いずれ使う？
	int GetLinerFlags(int indx){	return m_iLinearDeformation[indx]; }
	int GetVolumeFlags(int indx){	return m_iVolumeConservation[indx]; }

	bool CheckIndx(int pIndx);
	int	 SearchIndx(int pIndx);

	const Vec3& GetDisVec(){	return m_vec3DisCm;	}
	void CalcDisplaceMentVectorCm();

	//デバッグ
	void DebugIndx(void);
	void DebugLayer(void);
};

#endif