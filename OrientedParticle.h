/*!
  @file Ice_OrientedParticle.h
	
  @brief Oriented Particle法による弾性体変形                                                                                                                                                                                                                                                                                                                                                                                                                                                        
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
//	Vec3 m_vec3AngularVel;								//角速度
//	Vec3 m_vec3ElipsoidRadius;							//楕円の各軸の半径
//
//	mk_Quaternion m_QuatOrgOrientation;					//初期姿勢
//	mk_Quaternion m_QuatCurOrientation;					//現在姿勢
//	mk_QUaternion m_QuatPrdOrientation;					//推定姿勢
//
//	rxMatrix3 m_mtrx3PrtAe;								//楕円のモーメントマトリックス A_elipsoid
//	rxMatrix3 m_mtrx3PrdMomentMtrx;						//推定モーメントマトリックス
//
//	rxMatrix3 m_mtrx3Rotation;							//回転行列
//
//	Vec3 m_vec3OrgPos;									//初期位置
//	Vec3 m_vec3PrdPos;									//推定位置
//
//	Vec3 m_vec3Force;									//力

public:
	//コンストラクタ
	OrientedParticle();
	~OrientedParticle();
};

#endif