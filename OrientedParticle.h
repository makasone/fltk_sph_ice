/*!
  @file Ice_OrientedParticle.h
	
  @brief Oriented Particle�@�ɂ��e���̕ό`                                                                                                                                                                                                                                                                                                                                                                                                                                                        
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
	Vec3 m_vec3AngularVel;								//�p���x
	Vec3 m_vec3ElipsoidRadius;							//�ȉ~�̊e���̔��a

	mk_Quaternion m_QuatOrgOrientation;					//�����p��
	mk_Quaternion m_QuatCurOrientation;					//���ݎp��
	mk_Quaternion m_QuatPrdOrientation;					//����p��

	rxMatrix3 m_mtrx3PrtA_elipsoid;						//�ȉ~�̕ό`���z�̐��`�ߎ�
	rxMatrix3 m_mtrx3PrdMomentMtrx;						//���[�����g�}�g���b�N�X

	rxMatrix3 m_mtrx3Rotation;							//��]�s��

	Vec3 m_vec3OrgPos;									//�����ʒu
	Vec3 m_vec3PrdPos;									//����ʒu

	Vec3 m_vec3Force;									//��

	float m_fMass;										//����
	int m_iId;											//���q�ԍ�

	rxShapeMatching* m_smCluster;						//���q�ɑΉ�����N���X�^

public:
//�R���X�g���N�^
	OrientedParticle();
	OrientedParticle(int id, float mass, Vec3 orgPos);

	void Init();

	~OrientedParticle();

	void Integrate();
	void Integrate_Itr();

	void Update();
	void Update_Itr();
	void Update_ItrEnd();

//�A�N�Z�b�T
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