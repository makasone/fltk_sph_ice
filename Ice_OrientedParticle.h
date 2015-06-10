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
private:
	Vec3 m_vec3OrgCm;									//初期のクラスタの重心
	Vec3 m_vec3NowCm;									//現在のクラスタの重心

	float m_fDefAmount;									//変形量

	rxMatrix3	m_mtrx3Apq;								//変形行列Apq
	rxMatrix3	m_mtrx3AqqInv;							//変形行列Aqqの逆行列	前計算可能

	static const float* s_pfPrtPos;						//位置のホストポインタ　読み込み専用
	static const float* s_pfPrtVel;						//速度のホストポインタ　読み込み専用

	static float* s_pfClstrPos;							//各クラスタの最終的な位置
	static float* s_pfClstrVel;							//各クラスタの最終的な速度

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

	rxMatrix3 Apq() const { return m_mtrx3Apq;	}

	float DefAmount() const { return m_fDefAmount;	}

	void AddParticle(OrientedParticle* orientedPrt);
	void AddParticle(const Vec3 &pos, double mass, int pIndx);

	void RemoveParticle(int pIndx);
	void ClearParticle();

	void CalcOrgCm();

	void UpdateCluster();

	static void CopyPrtToClstrPos(unsigned prtNum);

private:
	void CalcNowCm();
	void CalcClusterMomentMatrix(rxMatrix3& Apq, rxMatrix3& Aqq);

	void CalcForce(float dt);
	void CalcVelocity(float dt);
	void DampVelocity(float dt);	
	void ProjectConstraint(float dt);
	void ApplyEstimatedPosition(float dt);
	void CorrectVelocity(float dt);

	void ShapeMatchingNormal();
};

#endif