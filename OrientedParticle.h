/*!
  @file OrientedParticle.h
	
  @brief Oriented Particle法による弾性体変形                                                                                                                                                                                                                                                                                                                                                                                                                                                        
  @ref M. Muller et al., "Solid Simulation with Oriented Particles", SIGGRAPH2011. 
 
  @author Ryo Nakasone
  @date 2015-4
*/

#ifndef _ORIENTED_PARTICLE_H_
#define _ORIENTED_PARTICLE_H_

#include <math.h>
#include <cmath>
#include <iostream>

#include "IceStructure.h"
#include "mk_Quaternion.h"

#include "rx_utility.h"
#include "rx_matrix.h"

#include "rx_nnsearch.h"

class OrientedParticleBaseElasticObject;

using namespace std;

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
	rxMatrix3 m_mtrx3Symmetric;							//回転以外の成分を含む行列　拡大縮小行列？

	Vec3 m_vec3OrgPos;									//初期位置
	Vec3 m_vec3PrdPos;									//推定位置

	Vec3 m_vec3Vel;										//速度

	Vec3 m_vec3Axis;									//回転軸

	Vec3 m_vec3Force;									//力

	float m_fMass;										//質量
	int m_iId;											//粒子番号

	OrientedParticleBaseElasticObject* m_smCluster;		//粒子に対応するクラスタ

public:
//コンストラクタ
	OrientedParticle();
	OrientedParticle(int id, float mass, Vec3 orgPos);

	~OrientedParticle();

	void Init();

	void Integrate();
	void Integrate_Itr();
	void Integrate_NotSampled();
	void Integrate_NotSampled2(const IceStructure* iceStrct);
	void Integrate_NotSampled_Itr();
	void Integrate_NotSampled2_Itr(const IceStructure* iceStrct);

	void Update();
	void Update_Itr();
	void Update_Itr_NotSampled(const IceStructure* iceStrct);
	void Update_ItrEnd();

	void UpdatePosAndVel();
	void UpdateOrientation();
	void UpdateAngularVel();

//アクセッサ
	int Id() const {		return m_iId;	}
	float Mass() const {	return m_fMass;	}

	inline const Vec3& OrgPos() const {	return m_vec3OrgPos;	}

	inline const Vec3 PrdPos() const {	return m_vec3PrdPos;	}
	inline void PrdPos(Vec3 pos){	m_vec3PrdPos = pos;		}

	inline const Vec3 Velocity() const {	return m_vec3Vel;	}
	inline void Velocity(Vec3 vel){	m_vec3Vel = vel;	}

	inline const rxMatrix3& A_elipsoid() const {		return m_mtrx3PrtA_elipsoid;	}
	inline const rxMatrix3& MomentMatrix() const {	return m_mtrx3PrdMomentMtrx;	}

	inline void ElipsoidRadius(Vec3 radius){	m_vec3ElipsoidRadius = radius;	}

	inline void Rotation(rxMatrix3 rotation){	m_mtrx3Rotation = rotation;	}
	inline const rxMatrix3& Rotation() const { return m_mtrx3Rotation;	}

	inline void Symmetric(rxMatrix3 sym){	m_mtrx3Symmetric = sym;	}
	inline rxMatrix3 Symmetric(){		return m_mtrx3Symmetric;	}

	inline void CurOrientation(mk_Quaternion orientation){	m_QuatCurOrientation = orientation;	}
	inline const mk_Quaternion& CurOrientation() const {	return m_QuatCurOrientation;	}

	inline void PrdOrientation(mk_Quaternion orientation){	m_QuatPrdOrientation = orientation;	}
	inline const mk_Quaternion& PrdOrientation() const {	return m_QuatPrdOrientation;	}

	inline const Vec3& Axis() const {	return m_vec3Axis;	}

	inline const Vec3& AngularVel() const {	return m_vec3AngularVel;	}
	inline void AngularVel(Vec3 angl){	m_vec3AngularVel = angl;	}

	inline void Cluster(OrientedParticleBaseElasticObject* cluster){	m_smCluster = cluster;	}

private:
	void IntegrateEstimatedPos();
	void IntegrateEstimatedOrientation();
	void InterpolateOrientation(const IceStructure* iceStrct);
	void IntegrateA_elipsoid();
	void IntegrateMomentMatrix();

	void IntegrateA_elipsoid_Itr();
	void IntegrateMomentMatrix_Itr();
};

#endif