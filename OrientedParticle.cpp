#ifndef _ORIENTED_PARTICLE_CP_
#define _ORIENTED_PARTICLE_CP_

//quaternionを使うためのライブラリ
//TODO: なぜかここに置かないとエラーが出る　他のヘッダをincludeする前に置かないといけないみたい？
#include <Eigen\Dense>
#include <Eigen\Geometry>
using namespace Eigen;

#include <math.h>

#include "OrientedParticle.h"
#include "Ice_OrientedParticle.h"
#include "ShapeMatching.h"

typedef OrientedParticleBaseElasticObject OrientedCluster;

using namespace std;

const unsigned X = 0;
const unsigned Y = 1;
const unsigned Z = 2;



//-------------------------------------------------------------------------------------------------------------------
//うまくEigenをインクルードできなかったために作った変換関数
Quaternionf ConvertQuaternion(mk_Quaternion mk_q)
{
	return Quaternionf(mk_q.w, mk_q.x, mk_q.y, mk_q.z);
}

//なぜか，参照で渡さないと "error C2719: 'q': __declspec(align('16')) の仮引数は配置されません。" というエラーが出る
mk_Quaternion ConvertQuaternion(const Quaternionf& q)
{
	return mk_Quaternion(q.x(), q.y(), q.z(), q.w());
}
//--------------------------------------------------------------------------------------------------------------------




OrientedParticle::OrientedParticle()
{

}

OrientedParticle::OrientedParticle(int id, float mass, Vec3 orgPos)
{
	m_iId = id;
	m_fMass = mass;
	m_vec3OrgPos = orgPos;
	m_vec3PrdPos = orgPos;

	Init();
}

OrientedParticle::~OrientedParticle()
{

}

void OrientedParticle::Init()
{
	m_vec3AngularVel = Vec3(0.0);
	m_vec3ElipsoidRadius = Vec3(1.0);

	m_QuatOrgOrientation = ConvertQuaternion(Quaternionf::Identity());
	m_QuatCurOrientation = ConvertQuaternion(Quaternionf::Identity());
	m_QuatPrdOrientation = ConvertQuaternion(Quaternionf::Identity());

	m_mtrx3PrtA_elipsoid = rxMatrix3::Identity();
	m_mtrx3PrdMomentMtrx = rxMatrix3::Identity();
	m_mtrx3Rotation = rxMatrix3::Identity();

	m_smCluster = 0;
}

void OrientedParticle::Integrate()
{
	IntegrateEstimatedPos();
	IntegrateEstimatedOrientation();
	IntegrateA_elipsoid();
	IntegrateMomentMatrix();
}

void OrientedParticle::Integrate_Itr()
{
	IntegrateA_elipsoid_Itr();
	IntegrateMomentMatrix_Itr();
}

void OrientedParticle::Update()
{
	UpdatePosAndVel();
	UpdateOrientation();
	UpdateAngularVel();
}

void OrientedParticle::Update_Itr()
{
	//最終的な速度のためにclstrPosとVelはそのままじゃないといけない
	//姿勢のみ更新
	UpdateOrientation();
}

void OrientedParticle::Update_ItrEnd()
{
	UpdatePosAndVel();
	UpdateAngularVel();
}

//前フレームの速度から現フレームの位置を推定
void OrientedParticle::IntegrateEstimatedPos()
{
	Vec3 gravity = m_smCluster->gravity();
	float dt = m_smCluster->dt();

	//TODO: 本当はClstrPosなんていらないはず　いずれなくす
	float* clstrPos = OrientedCluster::GetClstrPosPointer();
	float* clstrVel = OrientedCluster::GetClstrVelPointer();

	int pIndx = m_iId * SM_DIM;
	m_vec3PrdPos[X] = clstrPos[pIndx+X] + (clstrVel[pIndx+X] + (gravity[X] * dt)) * dt;
	m_vec3PrdPos[Y] = clstrPos[pIndx+Y] + (clstrVel[pIndx+Y] + (gravity[Y] * dt)) * dt;
	m_vec3PrdPos[Z] = clstrPos[pIndx+Z] + (clstrVel[pIndx+Z] + (gravity[Z] * dt)) * dt;
}

//姿勢の推定
void OrientedParticle::IntegrateEstimatedOrientation()
{
	Vec3 angularVel = m_vec3AngularVel;
	float len = length(angularVel);
	Quaternionf qCurrent = ConvertQuaternion(m_QuatCurOrientation);
	Quaternionf qPredict = Quaternionf::Identity();
	
	if(len <= 0.0001f){
		qPredict = qCurrent;
	}
	else{
		Vec3 dir = angularVel / len;
		float ang = len * 0.01f;
		
		//クォータニオン同士の掛け算で回転
		float halfAng = ang * 0.5f;
		Vec3 vec = sin(halfAng) * dir;
		Quaternionf dq(cos(halfAng), vec[0], vec[1], vec[2]);
		
		qPredict = dq * qCurrent;
	}

	m_QuatPrdOrientation = ConvertQuaternion(qPredict);
	//m_QuatPrdOrientation = ConvertQuaternion(qCurrent); //角速度を無視する場合　ただし姿勢が一定になってしまう　反復処理で使う？
}

//楕円の変形勾配の線形近似
//論文の式(9)
void OrientedParticle::IntegrateA_elipsoid()
{
	//各パーティクルの回転を計算
	Quaternionf q0 = ConvertQuaternion(m_QuatOrgOrientation);
	Quaternionf qp = ConvertQuaternion(m_QuatPrdOrientation);

	Quaternionf r = qp * q0.inverse();

	Matrix<float, 3, 3, RowMajor> tmp_r = r.matrix();
	rxMatrix3 R(
		tmp_r(0, 0), tmp_r(0, 1), tmp_r(0, 2), 
		tmp_r(1, 0), tmp_r(1, 1), tmp_r(1, 2), 
		tmp_r(2, 0), tmp_r(2, 1), tmp_r(2, 2)
	);

	float m = 1.0f * 0.2f;

	rxMatrix3 tmp = rxMatrix3::Identity();
	tmp(0, 0) = m * 1.0f;
	tmp(1, 1) = m * 1.0f;
	tmp(2, 2) = m * 1.0f;

	rxMatrix3 Ae = tmp * R;
	m_mtrx3PrtA_elipsoid = Ae;
}

//最小二乗法で求めた正規方程式（の一部）高速化のためにばらしている
void OrientedParticle::IntegrateMomentMatrix()
{
	Vec3 orgPos = m_vec3OrgPos;
	Vec3 curPos = m_vec3PrdPos;

	rxMatrix3 momentMtrx;

	momentMtrx(0,0) = curPos[X] * orgPos[X];
	momentMtrx(0,1)	= curPos[X] * orgPos[Y];
	momentMtrx(0,2)	= curPos[X] * orgPos[Z];
	momentMtrx(1,0)	= curPos[Y] * orgPos[X];
	momentMtrx(1,1)	= curPos[Y] * orgPos[Y];
	momentMtrx(1,2)	= curPos[Y] * orgPos[Z];
	momentMtrx(2,0)	= curPos[Z] * orgPos[X];
	momentMtrx(2,1)	= curPos[Z] * orgPos[Y];
	momentMtrx(2,2)	= curPos[Z] * orgPos[Z];

	m_mtrx3PrdMomentMtrx = momentMtrx;
}

//反復処理--------------------------------------------------------------------------------------------------------------------------
//楕円の変形勾配の線形近似
//論文の式(9)
void OrientedParticle::IntegrateA_elipsoid_Itr()
{
	//各パーティクルの回転を計算
	Quaternionf q0 = ConvertQuaternion(m_QuatOrgOrientation);
	Quaternionf qp = ConvertQuaternion(m_QuatPrdOrientation);

	Quaternionf r = qp * q0.inverse();

	Matrix<float, 3, 3, RowMajor> tmp_r = r.matrix();
	rxMatrix3 R(
		tmp_r(0, 0), tmp_r(0, 1), tmp_r(0, 2), 
		tmp_r(1, 0), tmp_r(1, 1), tmp_r(1, 2), 
		tmp_r(2, 0), tmp_r(2, 1), tmp_r(2, 2)
	);

	float m = 1.0f * 0.2f;

	rxMatrix3 tmp = rxMatrix3::Identity();
	tmp(0, 0) = m * 1.0f;
	tmp(1, 1) = m * 1.0f;
	tmp(2, 2) = m * 1.0f;

	rxMatrix3 Ae = tmp * R;
	m_mtrx3PrtA_elipsoid = Ae;
}

//最小二乗法で求めた正規方程式（の一部）高速化のためにばらしている
void OrientedParticle::IntegrateMomentMatrix_Itr()
{
	Vec3 orgPos = m_vec3OrgPos;
	Vec3 curPos = m_vec3PrdPos;	//Prdで良いのか？

	rxMatrix3 momentMtrx;

	momentMtrx(0,0) = curPos[X] * orgPos[X];
	momentMtrx(0,1)	= curPos[X] * orgPos[Y];
	momentMtrx(0,2)	= curPos[X] * orgPos[Z];
	momentMtrx(1,0)	= curPos[Y] * orgPos[X];
	momentMtrx(1,1)	= curPos[Y] * orgPos[Y];
	momentMtrx(1,2)	= curPos[Y] * orgPos[Z];
	momentMtrx(2,0)	= curPos[Z] * orgPos[X];
	momentMtrx(2,1)	= curPos[Z] * orgPos[Y];
	momentMtrx(2,2)	= curPos[Z] * orgPos[Z];

	m_mtrx3PrdMomentMtrx = momentMtrx;
}
//--------------------------------------------------------------------------------------------------------------------------反復処理

//速度算出，位置更新
void OrientedParticle::UpdatePosAndVel()
{
	float dt1 = 1.0f / m_smCluster->dt();

	//TODO: 本当はClstrPosなんていらないはず　いずれなくす
	float* clstrPos = OrientedCluster::GetClstrPosPointer();
	float* clstrVel = OrientedCluster::GetClstrVelPointer();

	int cIndx = m_iId * SM_DIM;
	clstrVel[cIndx+X] = (m_vec3PrdPos[X] - clstrPos[cIndx+X]) * dt1;
	clstrVel[cIndx+Y] = (m_vec3PrdPos[Y] - clstrPos[cIndx+Y]) * dt1;
	clstrVel[cIndx+Z] = (m_vec3PrdPos[Z] - clstrPos[cIndx+Z]) * dt1;

	clstrPos[cIndx+X] = m_vec3PrdPos[X];
	clstrPos[cIndx+Y] = m_vec3PrdPos[Y];
	clstrPos[cIndx+Z] = m_vec3PrdPos[Z];
}

//姿勢の更新
void OrientedParticle::UpdateOrientation()
{
	rxMatrix3 rotateMtrx = m_mtrx3Rotation;

	Matrix<float, 3, 3, RowMajor> tmp_r;
	tmp_r <<	rotateMtrx(0, 0), rotateMtrx(0, 1), rotateMtrx(0, 2), 
				rotateMtrx(1, 0), rotateMtrx(1, 1), rotateMtrx(1, 2), 
				rotateMtrx(2, 0), rotateMtrx(2, 1), rotateMtrx(2, 2);

	Quaternionf rotate(tmp_r);
	m_QuatPrdOrientation = ConvertQuaternion(rotate * ConvertQuaternion(m_QuatOrgOrientation));
}

//角速度の更新
void OrientedParticle::UpdateAngularVel()
{
	//角速度
	Quaternionf qp = ConvertQuaternion(m_QuatPrdOrientation);
	Quaternionf q_inv = ConvertQuaternion(m_QuatCurOrientation).inverse();

	Quaternionf r = qp * q_inv;

	if(r.w() < 0.0f){
		r.x() = -r.x();
		r.y() = -r.y();
		r.z() = -r.z();
		r.w() = -r.w();
	}

	//axis
	Vec3 axis;
	Vec3 vec(r.x(), r.y(), r.z());
	float len = norm(vec);
	if(len <= 0.0001f){
		axis = Vec3(1.0f, 0.0f, 0.0f);
	}
	else{
		axis = vec/len;
	}

	//angle
	float angle;
	Quaternionf normlized = r;
	normlized.normalize();

	float w = abs(normlized.w());
	if(w < 1.0f){
		angle = abs(acos(w) * 2.0f);
	}
	else{
		angle = 0.0f;
	}

	//angular vel
	float dt1 = 1.0f / m_smCluster->dt();

	if(angle < 0.0001f){
		m_vec3AngularVel = Vec3(0.0f);
	}
	else{
		m_vec3AngularVel = axis * angle * dt1;
	}

	Quaternionf qp_norm = ConvertQuaternion(m_QuatPrdOrientation);
	m_QuatCurOrientation = ConvertQuaternion(qp_norm.normalized());
}


#endif