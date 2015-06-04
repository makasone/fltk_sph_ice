/*!
  @file Ice_OrientedParticle.h
	
  @brief Oriented Particle�@�ɂ��e���̕ό`                                                                                                                                                                                                                                                                                                                                                                                                                                                        
  @ref M. Muller et al., "Solid Simulation with Oriented Particles", SIGGRAPH2011. 
 
  @author Ryo Nakasone
  @date 2015-4
*/

#ifndef _ORIENTED_PARTICLE_
#define _ORIENTED_PARTICLE_

#include <math.h>

using namespace std;

class OrientedParticle
{
//public:
//	struct mk_Quaternion{
//		float	x, y, z, w;
//
//		mk_Quaternion(){
//			mk_Quaternion(1.0f, 0.0f, 0.0f, 0.0f);
//		}
//
//		mk_Quaternion(float xx, float yy, float zz, float ww){
//			x = xx; y = yy; z = zz; w = ww;
//		}
//	};
//
//private:
//	Vec3 m_vec3AngularVel;								//�p���x
//	Vec3 m_vec3ElipsoidRadius;							//�ȉ~�̊e���̔��a
//
//	mk_Quaternion m_QuatOrgOrientation;					//�����p��
//	mk_Quaternion m_QuatCurOrientation;					//���ݎp��
//	mk_QUaternion m_QuatPrdOrientation;					//����p��
//
//	rxMatrix3 m_mtrx3PrtAe;								//�ȉ~�̃��[�����g�}�g���b�N�X A_elipsoid
//	rxMatrix3 m_mtrx3PrdMomentMtrx;						//���胂�[�����g�}�g���b�N�X
//
//	rxMatrix3 m_mtrx3Rotation;							//��]�s��
//
//	Vec3 m_vec3OrgPos;									//�����ʒu
//	Vec3 m_vec3PrdPos;									//����ʒu
//
//	Vec3 m_vec3Force;									//��

public:
	//�R���X�g���N�^
	OrientedParticle();
	~OrientedParticle();
};

#endif