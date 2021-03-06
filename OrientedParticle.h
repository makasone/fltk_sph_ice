/*!
  @file OrientedParticle.h
	
  @brief Oriented Particle@Ιζιe«ΜΟ`                                                                                                                                                                                                                                                                                                                                                                                                                                                        
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
	Vec3 m_vec3AngularVel;								//p¬x
	Vec3 m_vec3ElipsoidRadius;							//Θ~Μe²ΜΌa

	mk_Quaternion m_QuatOrgOrientation;					//ϊp¨
	mk_Quaternion m_QuatCurOrientation;					//»έp¨
	mk_Quaternion m_QuatPrdOrientation;					//θp¨

	rxMatrix3 m_mtrx3PrtA_elipsoid;						//Θ~ΜΟ`ωzΜό`ί
	rxMatrix3 m_mtrx3PrdMomentMtrx;						//[g}gbNX

	rxMatrix3 m_mtrx3Rotation;							//ρ]sρ
	rxMatrix3 m_mtrx3Symmetric;							//ρ]ΘOΜ¬ͺπάήsρ@gεk¬sρH

	Vec3 m_vec3OrgPos;									//ϊΚu
	Vec3 m_vec3PrdPos;									//θΚu

	Vec3 m_vec3Vel;										//¬x

	Vec3 m_vec3Axis;									//ρ]²

	Vec3 m_vec3Force;									//Ν

	float m_fMass;										//ΏΚ
	int m_iId;											//±qΤ

	OrientedParticleBaseElasticObject* m_smCluster;		//±qΙΞ·ιNX^

public:
//RXgN^
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

//ANZbT
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