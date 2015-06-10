/*!
  @file Ice_OrientedParticle.h
	
  @brief Oriented Particle法による弾性体変形                                                                                                                                                                                                                                                                                                                                                                                                                                                        
  @ref M. Muller et al., "Solid Simulation with Oriented Particles", SIGGRAPH2011. 
 
  @author Ryo Nakasone
  @date 2015-4
*/

#ifndef _ORIENTED_PARTICLE_H_
#define _ORIENTED_PARTICLE_H_

#include <math.h>

#include "rx_utility.h"
#include "rx_matrix.h"

#include "rx_nnsearch.h"

class rxShapeMatching;

using namespace std;

struct mk_Quaternion{
	float	w, x, y, z;

	mk_Quaternion(){
		mk_Quaternion(1.0f, 0.0f, 0.0f, 0.0f);
	}

	mk_Quaternion(float xx, float yy, float zz, float ww){
		x = xx; y = yy; z = zz; w = ww;
	}
};

class OrientedParticle
{
private:
	Vec3 m_vec3AngularVel;								//角速度
	Vec3 m_vec3ElipsoidRadius;							//楕円の各軸の半径

	mk_Quaternion m_QuatOrgOrientation;					//初期姿勢
	mk_Quaternion m_QuatCurOrientation;					//現在姿勢
	mk_Quaternion m_QuatPrdOrientation;					//推定姿勢

	rxMatrix3 m_mtrx3PrtA_elipsoid;						//楕円の変形勾配の線形近似
	rxMatrix3 m_mtrx3PrdMomentMtrx;						//モーメントマトリックス

	rxMatrix3 m_mtrx3Rotation;							//回転行列

	Vec3 m_vec3OrgPos;									//初期位置
	Vec3 m_vec3PrdPos;									//推定位置

	Vec3 m_vec3Force;									//力

	float m_fMass;										//質量
	int m_iId;											//粒子番号

	rxShapeMatching* m_smCluster;						//粒子に対応するクラスタ

public:
//コンストラクタ
	OrientedParticle();
	OrientedParticle(int id, float mass, Vec3 orgPos);

	void Init();

	~OrientedParticle();

	void Integrate();
	void Integrate_Itr();

	void Update();
	void Update_Itr();
	void Update_ItrEnd();

//アクセッサ
	int Id() const {		return m_iId;	}
	float Mass() const {	return m_fMass;	}
	Vec3 OrgPos() const {	return m_vec3OrgPos;	}
	Vec3& PrdPos(){	return m_vec3PrdPos;	}

	rxMatrix3 A_elipsoid() const {		return m_mtrx3PrtA_elipsoid;	}
	rxMatrix3 MomentMatrix() const {	return m_mtrx3PrdMomentMtrx;	}

	void Rotation(rxMatrix3 rotation){	m_mtrx3Rotation = rotation;	}

	void CurOrientation(mk_Quaternion orientation){	m_QuatCurOrientation = orientation;	}

	Vec3 AngularVel() const {	return m_vec3AngularVel;	}
	void AngularVel(Vec3 angl){	m_vec3AngularVel = angl;	}

	void Cluster(rxShapeMatching* cluster){	m_smCluster = cluster;	}

private:
	
	void IntegrateEstimatedPos();
	void IntegrateEstimatedOrientation();
	void IntegrateA_elipsoid();
	void IntegrateMomentMatrix();

	void IntegrateA_elipsoid_Itr();
	void IntegrateMomentMatrix_Itr();


	void UpdatePosAndVel();
	void UpdateOrientation();
	void UpdateAngularVel();
};

#endif