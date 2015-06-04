/*!
  @file Ice_OrientedParticle.h
	
  @brief Oriented Particle�@�ɂ��e���̕ό`                                                                                                                                                                                                                                                                                                                                                                                                                                                        
  @ref M. Muller et al., "Solid Simulation with Oriented Particles", SIGGRAPH2011. 
  �Q�l�ɂ����y�[�W�C�v���O���� http://yuki-koyama.hatenablog.com/entry/2015/01/30/150419
 
  @author Ryo Nakasone
  @date 2015-4
*/

//TODO: �Ȃ��������ɒu���Ȃ��ƃG���[���o��
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
//���܂�Eigen���C���N���[�h�ł��Ȃ��������߂ɍ�����ϊ��֐�
Quaternionf ConvertQuaternion(OrientedPrtObj::mk_Quaternion mk_q)
{
	return Quaternionf(mk_q.w, mk_q.x, mk_q.y, mk_q.z);
}

//�Ȃ����C�Q�Ƃœn���Ȃ��� "error C2719: 'q': __declspec(align('16')) �̉������͔z�u����܂���B" �Ƃ����G���[���o��
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
 * �R���X�g���N�^
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
 * �R�s�[�R���X�g���N�^
 */
OrientedParticleBaseElasticObject::OrientedParticleBaseElasticObject(const OrientedPrtObj& copy) : rxShapeMatching(copy)
{
	Copy(copy);
}

/*!
 * �f�X�g���N�^
 */
OrientedParticleBaseElasticObject::~OrientedParticleBaseElasticObject()
{
	ReleaseMemory();
}

//������Z�q�ŃR�s�[
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

//���������
void OrientedPrtObj::ReleaseMemory()
{

}

//�f�[�^�̃R�s�[
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

		//�p���x
		s_vvec3AngularVel.push_back(Vec3(0.0));
	
		//�p��
		s_vQuatCurOrientation.push_back(ConvertQuaternion(Quaternionf::Identity()));
		s_vQuatOrgOrientation.push_back(ConvertQuaternion(Quaternionf::Identity()));
		s_vQuatPrdOrientation.push_back(ConvertQuaternion(Quaternionf::Identity()));
		
		//�s��
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

	//�ő�Y���ԍ��̍X�V
	if(m_iNumVertices > m_iIndxNum){	m_iIndxNum = m_iNumVertices;	}

	//�d�S�ʒu�̍X�V
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

//�p�[�e�B�N�����̍X�V
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
		//�p��
		rxMatrix3 rotateMtrx = s_vmtrx3Rotation[pIndx];
		Matrix<float, 3, 3, RowMajor> tmp_r;
		tmp_r <<	rotateMtrx(0, 0), rotateMtrx(0, 1), rotateMtrx(0, 2), 
					rotateMtrx(1, 0), rotateMtrx(1, 1), rotateMtrx(1, 2), 
					rotateMtrx(2, 0), rotateMtrx(2, 1), rotateMtrx(2, 2);

		Quaternionf rotate(tmp_r);
		s_vQuatPrdOrientation[pIndx] = ConvertQuaternion(rotate * ConvertQuaternion(s_vQuatOrgOrientation[pIndx]));

		//���x�Z�o�C�ʒu�X�V
		float dt1 = 1.0f/0.01f;
		for(int j = 0; j < SM_DIM; j++){
			int cIndx = pIndx*SM_DIM+j;
			s_pfClstrVel[cIndx] = (s_vvec3Xp[pIndx][j] - s_pfClstrPos[cIndx]) * dt1;
			s_pfClstrPos[cIndx] = s_vvec3Xp[pIndx][j];
		}

		//�p���x
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
		//���x�Z�o�C�ʒu�X�V
		float dt1 = 1.0f/0.01f;
		for(int j = 0; j < SM_DIM; j++){
			int cIndx = pIndx*SM_DIM+j;
			s_pfClstrPos[cIndx] = s_vvec3Xp[pIndx][j];
		}

		//�p��
		rxMatrix3 rotateMtrx = s_vmtrx3Rotation[pIndx];
		Matrix<float, 3, 3, RowMajor> tmp_r;
		tmp_r <<	rotateMtrx(0, 0), rotateMtrx(0, 1), rotateMtrx(0, 2), 
					rotateMtrx(1, 0), rotateMtrx(1, 1), rotateMtrx(1, 2), 
					rotateMtrx(2, 0), rotateMtrx(2, 1), rotateMtrx(2, 2);

		Quaternionf rotate(tmp_r);
		s_vQuatPrdOrientation[pIndx] = ConvertQuaternion(rotate * ConvertQuaternion(s_vQuatOrgOrientation[pIndx]));

		//�p���x
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

		//�p��
		Quaternionf qp_norm = ConvertQuaternion(s_vQuatPrdOrientation[pIndx]);
		s_vQuatCurOrientation[pIndx] = ConvertQuaternion(qp_norm.normalized());
	}
}

void OrientedPrtObj::UpdateParticleVel()
{
	for(unsigned pIndx = 0; pIndx < s_vQuatOrgOrientation.size(); pIndx++){
		//���x�Z�o�C�ʒu�X�V
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

	////�]���̃V�F�C�v�}�b�`���O
	//ShapeMatchingNormal();
}

//�����d�S�ʒu
void OrientedPrtObj::CalcOrgCm()
{
	//�d�S�̍X�V
	float massSum = 0.0f;	// ������
	m_vec3OrgCm = Vec3(0.0);

	// �d�S���W�̌v�Z
	for(int i = 0; i < m_iIndxNum;++i){
		if(CheckHole(i)){	continue;	}

		float m = m_pMass[i];
		int pIndx = i*SM_DIM;

		massSum += m;
		m_vec3OrgCm += Vec3(m_pOrgPos[pIndx+0], m_pOrgPos[pIndx+1], m_pOrgPos[pIndx+2]) * m;
	}

	m_vec3OrgCm /= massSum;
}

//���݂̏d�S�ʒu
void OrientedPrtObj::CalcNowCm()
{
	//�d�S�̍X�V
	float massSum = 0.0f;	// ������
	m_vec3NowCm = Vec3(0.0);

	// �d�S���W�̌v�Z
	for(int i = 0; i < m_iIndxNum;++i){
		if(CheckHole(i)){	continue;	}

		float m = m_pMass[i];
		if(m_pFix[i]) m *= 300.0;	// �Œ�_�̎��ʂ�傫������
		
		int pIndx = i*SM_DIM;

		massSum += m;
		m_vec3NowCm += Vec3(m_pCurPos[pIndx+0], m_pCurPos[pIndx+1], m_pCurPos[pIndx+2]) * m;
	}

	m_vec3NowCm /= massSum;
}

//�ȉ~�̉�]���l������MomentMatrix
//�_���̎�(9)
void OrientedPrtObj::CalcAe()
{
	//�e�p�[�e�B�N���̉�]���v�Z
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

//���_�ʒu�����MomentMatrix
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

//�V�F�C�v�}�b�`���O�̃e�X�g
void OrientedPrtObj::ShapeMatchingNormal()
{
	double dt = m_dDt;

	//�ŏI�ʒu�E���x���N���X�^�̊e���q�ɔ��f�@�e�N���X�^�̊e���q�͓��ꂳ���
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

	// ���E�ǂ̉e��
	//���������Ȃ�d���Ȃ邪�C����͂���݂���
	double res = 0.9;	// �����W��
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
	Vec3 cm(0.0), cm_org(0.0);	// �d�S
	double mass = 0.0;	// ������
	rxMatrix3 R, S;

	//�O�t���[���ƌ��t���[���̈ʒu�̍��̑傫�����玿�ʂ�����
	//CalcMass();

	// �d�S���W�̌v�Z
	for(int i = 0; i < m_iIndxNum; ++i)
	{
		if( CheckHole(i) ){	continue;	}

		double m = m_pMass[i];
		if(m_pFix[i]) m *= 300.0;	// �Œ�_�̎��ʂ�傫������

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

	// Apq = ��mpq^T
	// Aqq = ��mqq^T
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

		// �ڕW���W���v�Z���C���݂̒��_���W���ړ�
		for(int i = 0; i < m_iIndxNum; ++i)
		{
			if( CheckHole(i) ){	continue;	}

			if(m_pFix[i]) continue;

			int cIndx = i*SM_DIM;

			// ��]�s��R�̑���̍s��RL=��A+(1-��)R���v�Z
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

//�e���q�ɉ����͂��X�V�@��ɂ͏d�͂Ȃ�
void OrientedPrtObj::CalcForce(float dt)
{
	//for(vector<Vec3>::iterator it = s_vvec3Force.begin(); it != s_vvec3Force.end(); it++){
	//	*it = Vec3(0.0);

	//	//�d�͂�K�p
	//	(*it) += m_v3Gravity;
	//}
}

void OrientedPrtObj::CalcVelocity(float dt)
{
}

void OrientedPrtObj::DampVelocity(float dt)
{
}

//�}�E�X�ɂ��h���b�O�𔽉f�����邽�߂ɁC�������l���X�V
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

//�ʒu�X�V�C����
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

//�����p�̈ʒu�X�V
void OrientedPrtObj::CalcEstimatedPosItr()
{
	for(unsigned i = 0; i < s_vQuatPrdOrientation.size(); i++){
		int pIndx = i * SM_DIM;

		for(unsigned j = 0; j < SM_DIM; j++){
			s_vvec3Xp[i][j] = s_pfClstrPos[pIndx+j];
		}
	}
}

//�p������
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
			
			//�N�H�[�^�j�I�����m�̊|���Z�ŉ�]
			float halfAng = ang * 0.5f;
			Vec3 vec = sin(halfAng) * dir;
			Quaternionf dq(cos(halfAng), vec[0], vec[1], vec[2]);
			
			qPredict = dq * qCurrent;
		}

		s_vQuatPrdOrientation[i] = ConvertQuaternion(qPredict);
		//s_vQuatPrdOrientation[i] = ConvertQuaternion(qCurrent); //�p���x�𖳎�����ꍇ
	}
}

void OrientedPrtObj::CalcClusterMomentMatrix(rxMatrix3& Apq, rxMatrix3& Aqq)
{
	rxMatrix3 Asum = rxMatrix3(0.0f);
	for(int i = 0; i < m_iIndxNum; ++i){
		if( CheckHole(i) ){	continue;	}

		int prtIndx = m_iPIndxes[i];
		Asum += s_vmtrx3PrtAe[prtIndx];	//��]���l������Ƃ��������Ȃ�
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
	//�ŏI�ʒu�E���x���N���X�^�̊e���q�ɔ��f�@�e�N���X�^�̊e���q�͓��ꂳ���
	for(int i = 0; i < m_iIndxNum; ++i){
		if( CheckHole(i) ){	continue;	}

		int pIndx = m_iPIndxes[i];
		int cIndx = i * SM_DIM;
		for(int j = 0; j < SM_DIM; j++){
			m_pCurPos[cIndx+j] = s_vvec3Xp[pIndx][j];
		}
	}

	// ���E�ǂ̉e��
	//���������Ȃ�d���Ȃ邪�C����͂���݂���
	double res = 0.9;	// �����W��
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

	//�e���̂̐���@�V�F�C�v�}�b�`���O�v�Z
	CalcNowCm();							//�N���X�^�̏d�S���X�V
	
	rxMatrix3 Apq(0.0), Aqq(0.0);
	CalcClusterMomentMatrix(Apq, Aqq);		//�p�[�e�B�N���̎p�����l���������[�����g�}�g���b�N�X

	rxMatrix3 R, S;
	PolarDecomposition(Apq, R, S);
	if(R.Determinant() < 0.0f){
		R = -1.0f * R;
	}

	// �ڕW���W���v�Z���C���݂̒��_���W���ړ�
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