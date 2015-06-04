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

#include "ShapeMatching.h"
#include "OrientedParticle.h"

#include "QueryCounter.h"

using namespace std;

class OrientedParticleBaseElasticObject : public rxShapeMatching
{
public:
	struct mk_Quaternion{
		float	x, y, z, w;

		mk_Quaternion(){
			mk_Quaternion(1.0f, 0.0f, 0.0f, 0.0f);
		}

		mk_Quaternion(float xx, float yy, float zz, float ww){
			x = xx; y = yy; z = zz; w = ww;
		}
	};

private:
	Vec3 m_vec3OrgCm;									//�����̃N���X�^�̏d�S
	Vec3 m_vec3NowCm;									//���݂̃N���X�^�̏d�S

	float m_fDefAmount;									//�ό`��

	rxMatrix3	m_mtrx3Apq;								//�ό`�s��Apq
	rxMatrix3	m_mtrx3AqqInv;							//�ό`�s��Aqq�̋t�s��	�O�v�Z�\

	Vec3 m_vec3AngularVel;								//�p���x

	Vec3 m_vec3ElipsoidRadius;							//�ȉ~�̊e���̔��a

	static const float* s_pfPrtPos;						//�ʒu�̃z�X�g�|�C���^�@�ǂݍ��ݐ�p
	static const float* s_pfPrtVel;						//���x�̃z�X�g�|�C���^�@�ǂݍ��ݐ�p

	static float* s_pfClstrPos;							//�e�N���X�^�̍ŏI�I�Ȉʒu
	static float* s_pfClstrVel;							//�e�N���X�^�̍ŏI�I�ȑ��x

	static vector<mk_Quaternion> s_vQuatOrgOrientation;	//�e���q�̏����̎p��
	static vector<mk_Quaternion> s_vQuatCurOrientation;	//�e���q�̌��݂̎p��
	static vector<mk_Quaternion> s_vQuatPrdOrientation;	//�e���q�̈ꎞ�I�Ȏp��
	static vector<Vec3> s_vvec3AngularVel;				//�e���q�̊p���x
	static vector<rxMatrix3> s_vmtrx3PrtAe;				//�ȉ~�p�[�e�B�N���̃��[�����g�}�g���b�N�X?
	static vector<rxMatrix3> s_vmtrx3PrtMomentMtrx;		//�ȉ~�p�[�e�B�N���̃��[�����g�}�g���b�N�X?
	static vector<Vec3> s_vvec3X0;						//�������W
	static vector<Vec3> s_vvec3Xp;						//�ꎞ�I�Ȉʒu
	static vector<Vec3> s_vvec3Force;					//���q�ɓ�����
	static vector<rxMatrix3> s_vmtrx3Rotation;			//��]�s��

	vector<OrientedParticle*> m_vOrientedPrtes;			//�ȉ~���q

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

	static vector<Vec3>& Xpes(){ return s_vvec3Xp;	}

	rxMatrix3 Apq() const { return m_mtrx3Apq;	}

	float DefAmount() const { return m_fDefAmount;	}

	void AddParticle(const Vec3 &pos, double mass, int pIndx);
	void RemoveParticle(int pIndx);
	void ClearParticle();

	void CalcOrgCm();

	void UpdateCluster();
	static void IntegrateParticle();
	static void IntegrateParticleItr();
	static void UpdateParticle();
	static void UpdateParticleVel();
	static void UpdateParticleItr();

	static void CopyPrtToClstrPos();

private:
	void CalcNowCm();
	void CalcClusterMomentMatrix(rxMatrix3& Apq, rxMatrix3& Aqq);

	void CalcForce(float dt);
	void CalcVelocity(float dt);
	void DampVelocity(float dt);	
	void ProjectConstraint(float dt);
	void ApplyEstimatedPosition(float dt);
	void CorrectVelocity(float dt);

	static void CalcAe();
	static void CalcMomentMatrix();
	static void CalcEstimatedPos();
	static void CalcEstimatedPosItr();
	static void CalcEstimatedOrientation();

	void ShapeMatchingNormal();
};

#endif