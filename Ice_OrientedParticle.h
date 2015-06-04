/*!
  @file Ice_OrientedParticle.h
	
  @brief Oriented Particle法による弾性体変形                                                                                                                                                                                                                                                                                                                                                                                                                                                        
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
	Vec3 m_vec3OrgCm;									//初期のクラスタの重心
	Vec3 m_vec3NowCm;									//現在のクラスタの重心

	float m_fDefAmount;									//変形量

	rxMatrix3	m_mtrx3Apq;								//変形行列Apq
	rxMatrix3	m_mtrx3AqqInv;							//変形行列Aqqの逆行列	前計算可能

	Vec3 m_vec3AngularVel;								//角速度

	Vec3 m_vec3ElipsoidRadius;							//楕円の各軸の半径

	static const float* s_pfPrtPos;						//位置のホストポインタ　読み込み専用
	static const float* s_pfPrtVel;						//速度のホストポインタ　読み込み専用

	static float* s_pfClstrPos;							//各クラスタの最終的な位置
	static float* s_pfClstrVel;							//各クラスタの最終的な速度

	static vector<mk_Quaternion> s_vQuatOrgOrientation;	//各粒子の初期の姿勢
	static vector<mk_Quaternion> s_vQuatCurOrientation;	//各粒子の現在の姿勢
	static vector<mk_Quaternion> s_vQuatPrdOrientation;	//各粒子の一時的な姿勢
	static vector<Vec3> s_vvec3AngularVel;				//各粒子の角速度
	static vector<rxMatrix3> s_vmtrx3PrtAe;				//楕円パーティクルのモーメントマトリックス?
	static vector<rxMatrix3> s_vmtrx3PrtMomentMtrx;		//楕円パーティクルのモーメントマトリックス?
	static vector<Vec3> s_vvec3X0;						//初期座標
	static vector<Vec3> s_vvec3Xp;						//一時的な位置
	static vector<Vec3> s_vvec3Force;					//粒子に働く力
	static vector<rxMatrix3> s_vmtrx3Rotation;			//回転行列

	vector<OrientedParticle*> m_vOrientedPrtes;			//楕円粒子

public:
	//コンストラクタ
	OrientedParticleBaseElasticObject();
	OrientedParticleBaseElasticObject(int obj);
	OrientedParticleBaseElasticObject(const OrientedParticleBaseElasticObject& copy);

	~OrientedParticleBaseElasticObject();
	
	//演算子のオーバーロード
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