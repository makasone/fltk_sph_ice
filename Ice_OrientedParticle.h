/*!
  @file Ice_OrientedParticle.h
	
  @brief Oriented Particle�@�ɂ��e���̕ό`                                                                                                                                                                                                                                                                                                                                                                                                                                                        
  @ref M. Muller et al., "Solid Simulation with Oriented Particles", SIGGRAPH2011. 
 
  @author Ryo Nakasone
  @date 2015-4
*/

#ifndef _ICE_ORIENTED_
#define _ICE_ORIENTED_

#include <math.h>
#include <time.h>

#include "IceStructure.h"
#include "mk_Quaternion.h"
#include "ShapeMatching.h"
#include "OrientedParticle.h"

#include "QueryCounter.h"

using namespace std;

class OrientedParticleBaseElasticObject : public rxShapeMatching
{
private:
	Vec3 m_vec3OrgCm;									//�����̃N���X�^�̏d�S
	Vec3 m_vec3NowCm;									//���݂̃N���X�^�̏d�S

	float m_fDefAmount;									//�ό`��

	rxMatrix3	m_mtrx3Apq;								//�ό`�s��Apq
	rxMatrix3	m_mtrx3AqqInv;							//�ό`�s��Aqq�̋t�s��	�O�v�Z�\

	static const float* s_pfPrtPos;						//�ʒu�̃z�X�g�|�C���^�@�ǂݍ��ݐ�p
	static const float* s_pfPrtVel;						//���x�̃z�X�g�|�C���^�@�ǂݍ��ݐ�p

	static float* s_pfClstrPos;							//�e�N���X�^�̍ŏI�I�Ȉʒu
	static float* s_pfClstrVel;							//�e�N���X�^�̍ŏI�I�ȑ��x

	vector<OrientedParticle*> m_vOrientedPrtes;			//�ȉ~���q
	//IceStructure* m_pIceStructure;						//���̂Ɋւ���\���̏��@��������i�Ȃ̂Ŏg���Â炢

public:
	//�R���X�g���N�^
	OrientedParticleBaseElasticObject();
	OrientedParticleBaseElasticObject(int obj);
	OrientedParticleBaseElasticObject(const OrientedParticleBaseElasticObject& copy);

	~OrientedParticleBaseElasticObject();
	
	//���Z�q�̃I�[�o�[���[�h
	OrientedParticleBaseElasticObject& operator=(const OrientedParticleBaseElasticObject& copy);

	void Copy(const OrientedParticleBaseElasticObject& copy);

	static void AllocateStaticMemory(const float* prtPos, const float* prtVel, int prtNum);

	void ReleaseMemory();
	
	static float* GetClstrPosPointer() {	return s_pfClstrPos;	}
	static float* GetClstrVelPointer() {	return s_pfClstrVel;	}

	rxMatrix3 Apq() const { return m_mtrx3Apq;	}

	float DefAmount() const { return m_fDefAmount;	}

	OrientedParticle* Particle(int i) const{ return m_vOrientedPrtes.at(i);	}

	void AddParticle(OrientedParticle* orientedPrt);
	void AddParticle(const Vec3 &pos, double mass, int pIndx);

	void InitOrientation();

	void RemoveParticle(int pIndx);
	void ClearParticle();

	void CalcOrgCm();

	void UpdateCluster();
	void UpdateCluster_Sampling(const IceStructure* ice_struct);

	static void CopyPrtToClstrPos(unsigned prtNum);

private:
	void CalcNowCm();
	void CalcClusterMomentMatrix(rxMatrix3& Apq, rxMatrix3& Aqq);

	void CalcForce(float dt);
	void CalcVelocity(float dt);

	rxMatrix3 InterpolateRotation(const IceStructure* ice_struct);
	rxMatrix3 InterpolateApq(const IceStructure* ice_struct);

	void DampVelocity(float dt);	
	void ProjectConstraint(float dt);
	void DistanceConstraint(float dt);
	void ApplyEstimatedPosition(float dt);
	void CorrectVelocity(float dt);

	void ShapeMatchingNormal();
};

#endif