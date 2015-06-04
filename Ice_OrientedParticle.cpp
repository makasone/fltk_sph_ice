/*!
  @file Ice_OrientedParticle.h
	
  @brief Oriented Particle法による弾性体変形                                                                                                                                                                                                                                                                                                                                                                                                                                                        
  @ref M. Muller et al., "Solid Simulation with Oriented Particles", SIGGRAPH2011. 
  参考にしたページ，プログラム http://yuki-koyama.hatenablog.com/entry/2015/01/30/150419
 
  @author Ryo Nakasone
  @date 2015-4
*/

//TODO: なぜかここに置かないとエラーが出る
#include <Eigen\Dense>
#include <Eigen\Geometry>
using namespace Eigen;

#include "Ice_OrientedParticle.h"

typedef OrientedParticleBaseElasticObject OrientedPrtObj;

const float* OrientedPrtObj::s_pfPrtPos = 0;
const float* OrientedPrtObj::s_pfPrtVel = 0;

float* OrientedPrtObj::s_pfClstrPos = 0;
float* OrientedPrtObj::s_pfClstrVel = 0;

vector<OrientedPrtObj::mk_Quaternion> OrientedPrtObj::s_vQuatOrgOrientation;
vector<OrientedPrtObj::mk_Quaternion> OrientedPrtObj::s_vQuatCurOrientation;
vector<OrientedPrtObj::mk_Quaternion> OrientedPrtObj::s_vQuatPrdOrientation;
vector<Vec3> OrientedPrtObj::s_vvec3AngularVel;
vector<rxMatrix3> OrientedPrtObj::s_vmtrx3PrtAe;		
vector<rxMatrix3> OrientedPrtObj::s_vmtrx3PrtMomentMtrx;
vector<Vec3> OrientedPrtObj::s_vvec3X0;
vector<Vec3> OrientedPrtObj::s_vvec3Xp;	
vector<Vec3> OrientedPrtObj::s_vvec3Force;
vector<rxMatrix3> OrientedPrtObj::s_vmtrx3Rotation;

//-------------------------------------------------------------------------------------------------------------------
//うまくEigenをインクルードできなかったために作った変換関数
Quaternionf ConvertQuaternion(OrientedPrtObj::mk_Quaternion mk_q)
{
	return Quaternionf(mk_q.w, mk_q.x, mk_q.y, mk_q.z);
}

//なぜか，参照で渡さないと "error C2719: 'q': __declspec(align('16')) の仮引数は配置されません。" というエラーが出る
OrientedPrtObj::mk_Quaternion ConvertQuaternion(const Quaternionf& q)
{
	return OrientedPrtObj::mk_Quaternion(q.x(), q.y(), q.z(), q.w());
}
//--------------------------------------------------------------------------------------------------------------------

//OrientedParticleBaseElasticObject::OrientedPrtObj() : rxShapeMatching()
//{
//
//}

/*!
 * コンストラクタ
 */
OrientedParticleBaseElasticObject::OrientedParticleBaseElasticObject(int obj) : rxShapeMatching(obj)
{
	m_mtrx3Apq = rxMatrix3(0.0);

	m_fDefAmount = 0.0f;

	m_vec3AngularVel = Vec3(0.0f);

	m_vec3ElipsoidRadius = Vec3(1.0f);

	//AllocateMemory(prtNum);
}

/*!
 * コピーコンストラクタ
 */
OrientedParticleBaseElasticObject::OrientedParticleBaseElasticObject(const OrientedPrtObj& copy) : rxShapeMatching(copy)
{
	Copy(copy);
}

/*!
 * デストラクタ
 */
OrientedParticleBaseElasticObject::~OrientedParticleBaseElasticObject()
{
	ReleaseMemory();
}

//代入演算子でコピー
OrientedPrtObj& OrientedPrtObj::operator=(const OrientedPrtObj& copy)
{
	if(this != &copy){
		rxShapeMatching::ReleaseMemory();
		ReleaseMemory();

		rxShapeMatching::Copy(copy);
		Copy(copy);
	}

	return *this;
}

//メモリ解放
void OrientedPrtObj::ReleaseMemory()
{

}

//データのコピー
void OrientedPrtObj::Copy(const OrientedPrtObj& copy)
{
	m_mtrx3Apq = copy.Apq();
	m_iIndxNum = copy.GetIndxNum();
	m_fDefAmount = copy.DefAmount();
}

void OrientedPrtObj::AllocateStaticMemory(const float* prtPos, const float* prtVel, int vrtxNum)
{
	s_pfPrtPos = prtPos;
	s_pfPrtVel = prtVel;

	s_pfClstrPos = new float[vrtxNum*SM_DIM];
	s_pfClstrVel = new float[vrtxNum*SM_DIM];

	if(s_pfClstrPos == 0 || s_pfClstrVel == 0){
		cout << __FUNCTION__ << " Error!" << endl;
		return;
	}

	for(int i = 0; i < vrtxNum; ++i){
		int pIndx = i*4;
		int cIndx = i*SM_DIM;

		for(int j = 0; j < SM_DIM; j++){
			s_pfClstrPos[cIndx+j] = s_pfPrtPos[pIndx+j];
			s_pfClstrVel[cIndx+j] = s_pfPrtVel[pIndx+j];
		}

		//角速度
		s_vvec3AngularVel.push_back(Vec3(0.0));
	
		//姿勢
		s_vQuatCurOrientation.push_back(ConvertQuaternion(Quaternionf::Identity()));
		s_vQuatOrgOrientation.push_back(ConvertQuaternion(Quaternionf::Identity()));
		s_vQuatPrdOrientation.push_back(ConvertQuaternion(Quaternionf::Identity()));
		
		//行列
		s_vmtrx3PrtAe.push_back(rxMatrix3::Identity());
		s_vmtrx3PrtMomentMtrx.push_back(rxMatrix3::Identity());
		s_vmtrx3Rotation.push_back(rxMatrix3::Identity());

		//
		s_vvec3Force.push_back(Vec3(0.0f));
		s_vvec3X0.push_back(Vec3(s_pfClstrPos[cIndx+0], s_pfClstrPos[cIndx+1], s_pfClstrPos[cIndx+2]));
		s_vvec3Xp.push_back(Vec3(s_pfClstrPos[cIndx+0], s_pfClstrPos[cIndx+1], s_pfClstrPos[cIndx+2]));
	}
}

void OrientedPrtObj::AddParticle(const Vec3 &pos, double mass, int pIndx)
{
	rxShapeMatching::AddVertex(pos, mass, pIndx);

	//最大添字番号の更新
	if(m_iNumVertices > m_iIndxNum){	m_iIndxNum = m_iNumVertices;	}

	//重心位置の更新
	CalcOrgCm();
}

void OrientedPrtObj::RemoveParticle(int indx)
{
	rxShapeMatching::Remove(indx);
}

void OrientedPrtObj::ClearParticle()
{
	int indx = m_iIndxNum;

	rxShapeMatching::Clear();
}

//パーティクル情報の更新
void OrientedPrtObj::IntegrateParticle()
{
	CalcEstimatedPos();
	CalcEstimatedOrientation();
	CalcAe();
	CalcMomentMatrix();
}

void OrientedPrtObj::IntegrateParticleItr()
{
	CalcEstimatedPosItr();
	CalcEstimatedOrientation();
	CalcAe();
	CalcMomentMatrix();
}

void OrientedPrtObj::UpdateParticle()
{
	for(unsigned pIndx = 0; pIndx < s_vQuatOrgOrientation.size(); pIndx++){
		//姿勢
		rxMatrix3 rotateMtrx = s_vmtrx3Rotation[pIndx];
		Matrix<float, 3, 3, RowMajor> tmp_r;
		tmp_r <<	rotateMtrx(0, 0), rotateMtrx(0, 1), rotateMtrx(0, 2), 
					rotateMtrx(1, 0), rotateMtrx(1, 1), rotateMtrx(1, 2), 
					rotateMtrx(2, 0), rotateMtrx(2, 1), rotateMtrx(2, 2);

		Quaternionf rotate(tmp_r);
		s_vQuatPrdOrientation[pIndx] = ConvertQuaternion(rotate * ConvertQuaternion(s_vQuatOrgOrientation[pIndx]));

		//速度算出，位置更新
		float dt1 = 1.0f/0.01f;
		for(int j = 0; j < SM_DIM; j++){
			int cIndx = pIndx*SM_DIM+j;
			s_pfClstrVel[cIndx] = (s_vvec3Xp[pIndx][j] - s_pfClstrPos[cIndx]) * dt1;
			s_pfClstrPos[cIndx] = s_vvec3Xp[pIndx][j];
		}

		//角速度
		Quaternionf qp = ConvertQuaternion(s_vQuatPrdOrientation[pIndx]);
		Quaternionf q_inv = ConvertQuaternion(s_vQuatCurOrientation[pIndx]).inverse();

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
		if(angle < 0.0001f){
			s_vvec3AngularVel[pIndx] = Vec3(0.0f);
		}
		else{
			s_vvec3AngularVel[pIndx] = axis * angle * dt1;
		}

		Quaternionf qp_norm = ConvertQuaternion(s_vQuatPrdOrientation[pIndx]);
		s_vQuatCurOrientation[pIndx] = ConvertQuaternion(qp_norm.normalized());
	}
}

void OrientedPrtObj::UpdateParticleItr()
{
	for(unsigned pIndx = 0; pIndx < s_vQuatOrgOrientation.size(); pIndx++){
		//速度算出，位置更新
		float dt1 = 1.0f/0.01f;
		for(int j = 0; j < SM_DIM; j++){
			int cIndx = pIndx*SM_DIM+j;
			s_pfClstrPos[cIndx] = s_vvec3Xp[pIndx][j];
		}

		//姿勢
		rxMatrix3 rotateMtrx = s_vmtrx3Rotation[pIndx];
		Matrix<float, 3, 3, RowMajor> tmp_r;
		tmp_r <<	rotateMtrx(0, 0), rotateMtrx(0, 1), rotateMtrx(0, 2), 
					rotateMtrx(1, 0), rotateMtrx(1, 1), rotateMtrx(1, 2), 
					rotateMtrx(2, 0), rotateMtrx(2, 1), rotateMtrx(2, 2);

		Quaternionf rotate(tmp_r);
		s_vQuatPrdOrientation[pIndx] = ConvertQuaternion(rotate * ConvertQuaternion(s_vQuatOrgOrientation[pIndx]));

		//角速度
		Quaternionf qp = ConvertQuaternion(s_vQuatPrdOrientation[pIndx]);
		Quaternionf q_inv = ConvertQuaternion(s_vQuatCurOrientation[pIndx]).inverse();

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
		if(angle < 0.0001f){
			s_vvec3AngularVel[pIndx] = Vec3(0.0f);
		}
		else{
			s_vvec3AngularVel[pIndx] = axis * angle * dt1;
		}

		//姿勢
		Quaternionf qp_norm = ConvertQuaternion(s_vQuatPrdOrientation[pIndx]);
		s_vQuatCurOrientation[pIndx] = ConvertQuaternion(qp_norm.normalized());
	}
}

void OrientedPrtObj::UpdateParticleVel()
{
	for(unsigned pIndx = 0; pIndx < s_vQuatOrgOrientation.size(); pIndx++){
		//速度算出，位置更新
		float dt1 = 1.0f/0.01f;
		for(int j = 0; j < SM_DIM; j++){
			int cIndx = pIndx*SM_DIM+j;
			s_pfClstrVel[cIndx] = (s_pfPrtPos[pIndx*4+j] - s_pfClstrPos[cIndx]) * dt1;
		}
	}
}

void OrientedPrtObj::UpdateCluster()
{
	float dt = m_dDt;

	////CalcForce(dt);
	////CalcVelocity(dt);
	////DampVelocity(dt);
	ProjectConstraint(dt);
	////ApplyEstimatedPosition(dt);
	////CorrectVelocity(dt);

	////従来のシェイプマッチング
	//ShapeMatchingNormal();
}

//初期重心位置
void OrientedPrtObj::CalcOrgCm()
{
	//重心の更新
	float massSum = 0.0f;	// 総質量
	m_vec3OrgCm = Vec3(0.0);

	// 重心座標の計算
	for(int i = 0; i < m_iIndxNum;++i){
		if(CheckHole(i)){	continue;	}

		float m = m_pMass[i];
		int pIndx = i*SM_DIM;

		massSum += m;
		m_vec3OrgCm += Vec3(m_pOrgPos[pIndx+0], m_pOrgPos[pIndx+1], m_pOrgPos[pIndx+2]) * m;
	}

	m_vec3OrgCm /= massSum;
}

//現在の重心位置
void OrientedPrtObj::CalcNowCm()
{
	//重心の更新
	float massSum = 0.0f;	// 総質量
	m_vec3NowCm = Vec3(0.0);

	// 重心座標の計算
	for(int i = 0; i < m_iIndxNum;++i){
		if(CheckHole(i)){	continue;	}

		float m = m_pMass[i];
		if(m_pFix[i]) m *= 300.0;	// 固定点の質量を大きくする
		
		int pIndx = i*SM_DIM;

		massSum += m;
		m_vec3NowCm += Vec3(m_pCurPos[pIndx+0], m_pCurPos[pIndx+1], m_pCurPos[pIndx+2]) * m;
	}

	m_vec3NowCm /= massSum;
}

//楕円の回転を考慮したMomentMatrix
//論文の式(9)
void OrientedPrtObj::CalcAe()
{
	//各パーティクルの回転を計算
	for(unsigned i = 0; i < s_vQuatOrgOrientation.size(); i++){
		//case1
		Quaternionf q0 = ConvertQuaternion(s_vQuatOrgOrientation[i]);
		Quaternionf qp = ConvertQuaternion(s_vQuatPrdOrientation[i]);

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
		s_vmtrx3PrtAe[i] = Ae;
	}
}

//頂点位置からのMomentMatrix
void OrientedPrtObj::CalcMomentMatrix()
{
	for(unsigned i = 0; i < s_vmtrx3PrtMomentMtrx.size(); i++){
		Vec3 orgPos = s_vvec3X0[i];
		Vec3 curPos = s_vvec3Xp[i];

		rxMatrix3 momentMtrx;

		momentMtrx(0,0) = curPos[0] * orgPos[0];
		momentMtrx(0,1)	= curPos[0] * orgPos[1];
		momentMtrx(0,2)	= curPos[0] * orgPos[2];
		momentMtrx(1,0)	= curPos[1] * orgPos[0];
		momentMtrx(1,1)	= curPos[1] * orgPos[1];
		momentMtrx(1,2)	= curPos[1] * orgPos[2];
		momentMtrx(2,0)	= curPos[2] * orgPos[0];
		momentMtrx(2,1)	= curPos[2] * orgPos[1];
		momentMtrx(2,2)	= curPos[2] * orgPos[2];

		s_vmtrx3PrtMomentMtrx[i] = momentMtrx;
	}
}

//シェイプマッチングのテスト
void OrientedPrtObj::ShapeMatchingNormal()
{
	double dt = m_dDt;

	//最終位置・速度をクラスタの各粒子に反映　各クラスタの各粒子は統一される
	for(int i = 0; i < m_iIndxNum; ++i)
	{
		if( CheckHole(i) ){	continue;	}
		int pIndx = m_iPIndxes[i]*4;
		int cIndx = i*SM_DIM;

		for(int j = 0; j < SM_DIM; j++)
		{
			int jpIndx = pIndx+j;
			int jcIndx = cIndx+j;
			
			m_pCurPos[jcIndx] = s_pfPrtPos[jpIndx] + (s_pfPrtVel[jpIndx] + (m_v3Gravity[j] * dt)) *dt;
		}
	}

	// 境界壁の影響
	//処理がかなり重くなるが，安定はするみたい
	double res = 0.9;	// 反発係数
	for(int i = 0; i < m_iIndxNum; ++i)
	{
		if( CheckHole(i) ){	continue;	}
		if(m_pFix[i]) continue;

		int pIndx =  m_iPIndxes[i]*4;
		int cIndx = i*SM_DIM;

		if(m_pCurPos[cIndx+0] < m_v3Min[0] || m_pCurPos[cIndx+0] > m_v3Max[0]){
			m_pCurPos[cIndx+0] = s_pfPrtPos[pIndx+0] - (s_pfPrtVel[pIndx+0] + (m_v3Gravity[0] * dt))*dt*res;
			m_pCurPos[cIndx+1] = s_pfPrtPos[pIndx+1];
			m_pCurPos[cIndx+2] = s_pfPrtPos[pIndx+2];
		}

		if(m_pCurPos[cIndx+1] < m_v3Min[1] || m_pCurPos[cIndx+1] > m_v3Max[1]){
			m_pCurPos[cIndx+1] = s_pfPrtPos[pIndx+1] - (s_pfPrtVel[pIndx+1] + (m_v3Gravity[1] * dt))*dt*res;
			m_pCurPos[cIndx+0] = s_pfPrtPos[pIndx+0] ;
			m_pCurPos[cIndx+2] = s_pfPrtPos[pIndx+2];
		}

		if(m_pCurPos[cIndx+2] < m_v3Min[2] || m_pCurPos[cIndx+2] > m_v3Max[2]){
			m_pCurPos[cIndx+2] = s_pfPrtPos[pIndx+2] - (s_pfPrtVel[pIndx+2] + (m_v3Gravity[2] * dt))*dt*res;
			m_pCurPos[cIndx+0] = s_pfPrtPos[pIndx+0];
			m_pCurPos[cIndx+1] = s_pfPrtPos[pIndx+1];
		}

		clamp(m_pCurPos, cIndx);
	}

	if(m_iNumVertices <= 1) return;

	
	Vec3 p, q;
	Vec3 cm(0.0), cm_org(0.0);	// 重心
	double mass = 0.0;	// 総質量
	rxMatrix3 R, S;

	//前フレームと現フレームの位置の差の大きさから質量を決定
	//CalcMass();

	// 重心座標の計算
	for(int i = 0; i < m_iIndxNum; ++i)
	{
		if( CheckHole(i) ){	continue;	}

		double m = m_pMass[i];
		if(m_pFix[i]) m *= 300.0;	// 固定点の質量を大きくする

		int cIndx = i*SM_DIM;
		mass += m;

		for(int j = 0; j < SM_DIM; j++)
		{
			cm[j] += m_pCurPos[cIndx+j]*m;
		}
	}

	cm /= mass;
	m_vec3NowCm = cm;
	cm_org = m_vec3OrgCm;

	// Apq = Σmpq^T
	// Aqq = Σmqq^T
	rxMatrix3 Apq(0.0), Aqq(0.0);
	for(int i = 0; i < m_iIndxNum; ++i)
	{
		if( CheckHole(i) ){	continue;	}

		int cIndx = i*SM_DIM;

		for(int j = 0; j < SM_DIM; j++)
		{
			p[j] = m_pCurPos[cIndx+j]-cm[j];
			q[j] = m_pOrgPos[cIndx+j]-cm_org[j];
		}

		double m = m_pMass[i];
		p *= m;

		Apq(0,0) += p[0] * q[0];
		Apq(0,1) += p[0] * q[1];
		Apq(0,2) += p[0] * q[2];
		Apq(1,0) += p[1] * q[0];
		Apq(1,1) += p[1] * q[1];
		Apq(1,2) += p[1] * q[2];
		Apq(2,0) += p[2] * q[0];
		Apq(2,1) += p[2] * q[1];
		Apq(2,2) += p[2] * q[2];
	}


	PolarDecomposition(Apq, R, S);

	if(m_bLinearDeformation){
		m_fDefAmount = 0.0f;

		// 目標座標を計算し，現在の頂点座標を移動
		for(int i = 0; i < m_iIndxNum; ++i)
		{
			if( CheckHole(i) ){	continue;	}

			if(m_pFix[i]) continue;

			int cIndx = i*SM_DIM;

			// 回転行列Rの代わりの行列RL=βA+(1-β)Rを計算
			for(int j = 0; j < SM_DIM; j++)
			{
				q[j] = m_pOrgPos[cIndx+j]-cm_org[j];
			}

			Vec3 gp(R*q+cm);

			for(int j = 0; j < SM_DIM; j++)
			{
				float defAmount = (gp[j]-m_pCurPos[cIndx+j]);
				m_pCurPos[cIndx+j] += defAmount;

				m_fDefAmount += abs(defAmount);
			}
		}
	}

	double dt1 = 1.0/dt;
	for(int i = 0; i < m_iIndxNum; ++i)
	{
		if( CheckHole(i) ){	continue;	}

		int pIndx = m_iPIndxes[i]*4;

		for(int j = 0; j < SM_DIM; j++)
		{
			int cIndx = i*SM_DIM+j;
			m_pVel[cIndx] = (m_pCurPos[cIndx] - s_pfPrtPos[pIndx+j]) * dt1/* + m_v3Gravity[j] * dt*/;
		}
	}
}

//各粒子に加わる力を更新　主には重力など
void OrientedPrtObj::CalcForce(float dt)
{
	//for(vector<Vec3>::iterator it = s_vvec3Force.begin(); it != s_vvec3Force.end(); it++){
	//	*it = Vec3(0.0);

	//	//重力を適用
	//	(*it) += m_v3Gravity;
	//}
}

void OrientedPrtObj::CalcVelocity(float dt)
{
}

void OrientedPrtObj::DampVelocity(float dt)
{
}

//マウスによるドラッグを反映させるために，無理やり値を更新
void OrientedPrtObj::CopyPrtToClstrPos()
{
	for(unsigned i = 0; i < s_vQuatPrdOrientation.size(); i++){
		int pIndx = i * SM_DIM;

		for(unsigned j = 0; j < SM_DIM; j++){
			s_pfClstrPos[pIndx+j] = s_pfPrtPos[i*4+j];
			s_pfClstrVel[pIndx+j] = s_pfPrtVel[i*4+j];
		}
	}
}

//位置更新，推定
void OrientedPrtObj::CalcEstimatedPos()
{
	Vec3 gravity(0.0f, -9.81f, 0.0f);

	for(unsigned i = 0; i < s_vQuatPrdOrientation.size(); i++){
		int pIndx = i * SM_DIM;

		for(unsigned j = 0; j < SM_DIM; j++){
			s_vvec3Xp[i][j] = s_pfClstrPos[pIndx+j] + (s_pfClstrVel[pIndx+j] + (gravity[j] * 0.01f)) * 0.01f;
		}
	}
}

//反復用の位置更新
void OrientedPrtObj::CalcEstimatedPosItr()
{
	for(unsigned i = 0; i < s_vQuatPrdOrientation.size(); i++){
		int pIndx = i * SM_DIM;

		for(unsigned j = 0; j < SM_DIM; j++){
			s_vvec3Xp[i][j] = s_pfClstrPos[pIndx+j];
		}
	}
}

//姿勢推定
void OrientedPrtObj::CalcEstimatedOrientation()
{
	for(unsigned i = 0; i < s_vQuatPrdOrientation.size(); i++){
		Vec3 angularVel = s_vvec3AngularVel[i];
		float len = length(angularVel);
		Quaternionf qCurrent = ConvertQuaternion(s_vQuatCurOrientation[i]);
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

		s_vQuatPrdOrientation[i] = ConvertQuaternion(qPredict);
		//s_vQuatPrdOrientation[i] = ConvertQuaternion(qCurrent); //角速度を無視する場合
	}
}

void OrientedPrtObj::CalcClusterMomentMatrix(rxMatrix3& Apq, rxMatrix3& Aqq)
{
	rxMatrix3 Asum = rxMatrix3(0.0f);
	for(int i = 0; i < m_iIndxNum; ++i){
		if( CheckHole(i) ){	continue;	}

		int prtIndx = m_iPIndxes[i];
		Asum += s_vmtrx3PrtAe[prtIndx];	//回転を考慮するとおかしくなる
		Asum += s_vmtrx3PrtMomentMtrx[prtIndx];
	}

	rxMatrix3 MccT;

	MccT(0,0) = 1.0f * m_vec3NowCm[0] * m_vec3OrgCm[0];
	MccT(0,1) = 1.0f * m_vec3NowCm[0] * m_vec3OrgCm[1];
	MccT(0,2) = 1.0f * m_vec3NowCm[0] * m_vec3OrgCm[2];
	MccT(1,0) = 1.0f * m_vec3NowCm[1] * m_vec3OrgCm[0];
	MccT(1,1) = 1.0f * m_vec3NowCm[1] * m_vec3OrgCm[1];
	MccT(1,2) = 1.0f * m_vec3NowCm[1] * m_vec3OrgCm[2];
	MccT(2,0) = 1.0f * m_vec3NowCm[2] * m_vec3OrgCm[0];
	MccT(2,1) = 1.0f * m_vec3NowCm[2] * m_vec3OrgCm[1];
	MccT(2,2) = 1.0f * m_vec3NowCm[2] * m_vec3OrgCm[2];

	float massSum = (float)m_iIndxNum;

	Apq = Asum - massSum * MccT;
}

void OrientedPrtObj::ProjectConstraint(float dt)
{
	//最終位置・速度をクラスタの各粒子に反映　各クラスタの各粒子は統一される
	for(int i = 0; i < m_iIndxNum; ++i){
		if( CheckHole(i) ){	continue;	}

		int pIndx = m_iPIndxes[i];
		int cIndx = i * SM_DIM;
		for(int j = 0; j < SM_DIM; j++){
			m_pCurPos[cIndx+j] = s_vvec3Xp[pIndx][j];
		}
	}

	// 境界壁の影響
	//処理がかなり重くなるが，安定はするみたい
	double res = 0.9;	// 反発係数
	for(int i = 0; i < m_iIndxNum; ++i){
		if( CheckHole(i) ){	continue;	}
		if(m_pFix[i]) continue;

		int pIndx = m_iPIndxes[i];
		int cIndx = i*SM_DIM;

		if(m_pCurPos[cIndx+0] < m_v3Min[0] || m_pCurPos[cIndx+0] > m_v3Max[0]){
			m_pCurPos[cIndx+0] = s_vvec3Xp[pIndx][0] - (s_pfClstrVel[cIndx+0] + (m_v3Gravity[0] * dt))*dt*res;
			m_pCurPos[cIndx+1] = s_vvec3Xp[pIndx][1];
			m_pCurPos[cIndx+2] = s_vvec3Xp[pIndx][2];
		}

		if(m_pCurPos[cIndx+1] < m_v3Min[1] || m_pCurPos[cIndx+1] > m_v3Max[1]){
			m_pCurPos[cIndx+1] = s_vvec3Xp[pIndx][1] - (s_pfClstrVel[cIndx+1] + (m_v3Gravity[1] * dt))*dt*res;
			m_pCurPos[cIndx+0] = s_vvec3Xp[pIndx][0] ;
			m_pCurPos[cIndx+2] = s_vvec3Xp[pIndx][2];
		}

		if(m_pCurPos[cIndx+2] < m_v3Min[2] || m_pCurPos[cIndx+2] > m_v3Max[2]){
			m_pCurPos[cIndx+2] = s_vvec3Xp[pIndx][2] - (s_pfClstrVel[cIndx+2] + (m_v3Gravity[2] * dt))*dt*res;
			m_pCurPos[cIndx+0] = s_vvec3Xp[pIndx][0];
			m_pCurPos[cIndx+1] = s_vvec3Xp[pIndx][1];
		}

		clamp(m_pCurPos, cIndx);
	}

	//弾性体の制約　シェイプマッチング計算
	CalcNowCm();							//クラスタの重心を更新
	
	rxMatrix3 Apq(0.0), Aqq(0.0);
	CalcClusterMomentMatrix(Apq, Aqq);		//パーティクルの姿勢を考慮したモーメントマトリックス

	rxMatrix3 R, S;
	PolarDecomposition(Apq, R, S);
	if(R.Determinant() < 0.0f){
		R = -1.0f * R;
	}

	// 目標座標を計算し，現在の頂点座標を移動
	m_fDefAmount = 0.0f;

	for(int i = 0; i < m_iIndxNum; ++i){
		if( CheckHole(i) ){	continue;	}
		if(m_pFix[i]) continue;

		int cIndx = i*SM_DIM;

		Vec3 q;
		for(int j = 0; j < SM_DIM; j++){
			q[j] = m_pOrgPos[cIndx+j]-m_vec3OrgCm[j];
		}

		Vec3 gp(R*q+m_vec3NowCm);

		for(int j = 0; j < SM_DIM; j++){
			float defAmount = (gp[j]-m_pCurPos[cIndx+j]);
			m_pCurPos[cIndx+j] += defAmount;

			m_fDefAmount += abs(defAmount);
		}
	}

	s_vmtrx3Rotation[m_iPIndxes[0]] = R;
}

void OrientedPrtObj::ApplyEstimatedPosition(float dt)
{
}

void OrientedPrtObj::CorrectVelocity(float dt)
{
}