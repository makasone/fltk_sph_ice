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

using namespace std;

//

class Ice_SM : public rxShapeMatching
{
protected:
	vector<int> m_iParticleIndxes;		//!< クラスタに所属する粒子の番号

	vector<double> m_dAlphas;			//!< stiffnessパラメータ[0,1] (速度計算に使用)
	vector<double> m_dBetas;			//!< deformationパラメータ[0,1]

	Vec3 m_vec3Cm;						//クラスタの重心
	rxMatrix3	m_mtrx3Apq;				//変形行列

	vector<int> m_iLayeres;

	vector<int> m_iLinearDeformation;	//!< Linear/Quadratic deformation切り替えフラグ
	vector<int> m_iVolumeConservation;	//!< 変形時の体積保存性(√det(A)で割るかどうか)

	//ここに辺と面を持たせる
		//bounds壁にぶつかったときの反発力
		//allowFlip反転を許すかのフラグ　いらない？

public:
	Ice_SM(int obj);
	~Ice_SM();

	void AddVertex(const Vec3 &pos, double mass, int pIndx);

	void ShapeMatching(double dt);
	void Update();

	void SetAlphas(int indx, int alpha){ m_dAlphas[indx] = alpha;	}
	void SetBetas (int indx, int beta){	m_dBetas[indx] = beta;		}

	void SetLayer(int indx, int layer){	m_iLayeres[indx] = layer;	}

	void SetLinerFalg(int indx, int flag){	m_iLinearDeformation[indx] = flag; }
	void SetVolumeFlag(int indx, int flag){	m_iVolumeConservation[indx] = flag;}

	int GetParticleIndx(int indx){ return m_iParticleIndxes[indx]; }
	const vector<int>& GetVertexIndxList(){ return m_iParticleIndxes;	}

	double GetAlphas(int indx){	return m_dAlphas[indx];	}
	double GetBetas (int indx){	return m_dBetas[indx];	}

	Vec3 GetCm(void){		return m_vec3Cm;	}

	rxMatrix3 GetApq(void){	return m_mtrx3Apq;	}

	int GetLayer(int indx){	return m_iLayeres[indx];	}

	void Remove(int indx);
	void Clear();

	//いずれ使う？
	int GetLinerFlags(int indx){	return m_iLinearDeformation[indx]; }
	int GetVolumeFlags(int indx){	return m_iVolumeConservation[indx]; }

	bool CheckIndx(int pIndx);



	//デバッグ
	void DebugIndx(void);
	void DebugLayer(void);
};

#endif