/*!
  @file Ice_OrientedParticle.h
	
  @brief Oriented Particle�@�ɂ��e���̕ό`                                                                                                                                                                                                                                                                                                                                                                                                                                                        
  @ref M. Muller et al., "Solid Simulation with Oriented Particles", SIGGRAPH2011. 
  �Q�l�ɂ����y�[�W�C�v���O���� http://yuki-koyama.hatenablog.com/entry/2015/01/30/150419
 
  @author Ryo Nakasone
  @date 2015-4
*/
//���l�v�Z���C�u����
//TODO: �Ȃ��������ɒu���Ȃ��ƃG���[���o��@���̃w�b�_��include����O�ɒu���Ȃ��Ƃ����Ȃ��݂����H
#include <Eigen\Core>
#include <Eigen\Dense>
#include <Eigen\Geometry>
using namespace Eigen;

#include "Ice_OrientedParticle.h"

typedef OrientedParticleBaseElasticObject OrientedPrtObj;

const float* OrientedPrtObj::s_pfPrtPos = 0;
const float* OrientedPrtObj::s_pfPrtVel = 0;

float* OrientedPrtObj::s_pfClstrPos = 0;
float* OrientedPrtObj::s_pfClstrVel = 0;

const unsigned X = 0;
const unsigned Y = 1;
const unsigned Z = 2;

/*!
 * �R���X�g���N�^
 */
OrientedParticleBaseElasticObject::OrientedParticleBaseElasticObject(int obj) : rxShapeMatching(obj)
{
	m_mtrx3Apq = rxMatrix3(0.0);

	m_fDefAmount = 0.0f;

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
	}
}

void OrientedPrtObj::AddParticle(OrientedParticle* orientedPrt)
{
	Vec3 pos = orientedPrt->OrgPos();
	double mass = orientedPrt->Mass();
	int id = orientedPrt->Id();

	m_vOrientedPrtes.push_back(orientedPrt);

	AddParticle(pos, mass, id);
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

//�听�����͂��s���C�N���X�^���\�����闱�q�z�u���珉���p��������
void OrientedPrtObj::InitOrientation()
{
	CalcOrgCm();
	Vector3f ave = Vector3f(m_vec3OrgCm[X], m_vec3OrgCm[Y], m_vec3OrgCm[Z]);
	Matrix3f covarianceMtrx = Matrix3f::Zero();

	for(int i = 0; i < m_iIndxNum; i++){
		if(CheckHole(i)){	continue;	}

		int pIndx = i * SM_DIM;
		Vector3f pos = Vector3f(m_pOrgPos[pIndx+X], m_pOrgPos[pIndx+Y], m_pOrgPos[pIndx+Z]);
		Vector3f tmp = pos - ave;

		covarianceMtrx.row(0) += tmp(0) * tmp;
		covarianceMtrx.row(1) += tmp(1) * tmp;
		covarianceMtrx.row(2) += tmp(2) * tmp;
	}

	//EigenSolver���ƃR���p�C���œ����G���[���o��
	//EigenSolver<Matrix3f> es;
	//es.compute(covarianceMtrx, true);

	SelfAdjointEigenSolver<Matrix3f> es;
	es.compute(covarianceMtrx);

	//compute�̌��ʂ̓\�[�g����Ă��Ȃ��̂ŁC�ŗL�l�̑傫�����Ƀ\�[�g
	//eigen��Vector��stl�Ƒ����������݂���
	//�Ȃ���vector<pair>�ɓ�����Vector3d���o�͂��悤�Ƃ���ƃR���p�C�����ʂ�Ȃ�
	vector<pair<float, int>> eigenValVec;
	eigenValVec.push_back(pair<float, int>(es.eigenvalues()[0], 0));
	eigenValVec.push_back(pair<float, int>(es.eigenvalues()[1], 1));
	eigenValVec.push_back(pair<float, int>(es.eigenvalues()[2], 2));

	sort(eigenValVec.begin(), eigenValVec.end());

	//�ŗL�l�C�ŗL�x�N�g���őȉ~�̑傫���Ǝp��������
	//�\�[�g�������ʂ���C�ŗL�l�̏��������ɌŗL�x�N�g�������o��
	//�ȉ~�̊e���̌a�͌ŗL�l���狁�߂�

	//�{����eigen��vector3d���g���������COrientedParticle�ł͐搶��Vec3���g���Ă���̂Ŏd���Ȃ�
	//�ꉞ�C(1.0, 1.0, 1.0)�̎��ɏ����ł����킹�邽�߂ɐ��K�����ā�3�{���Ă݂�
	Vec3 radius = Vec3(0.0);
	Vec3 norm_eigenVal = Vec3(eigenValVec[X].first, eigenValVec[Y].first, eigenValVec[Z].first);
	radius = Unit(norm_eigenVal) * sqrt(3.0f);

	m_vOrientedPrtes[0]->ElipsoidRadius(radius);

	//�ő�ŗL�l�̕����������p���ƂȂ�
	int indx = eigenValVec[Z].second;
	Vector3f orientedVec = es.eigenvectors().row(indx);

	//�P	�����̃x�N�g������p�����`������@���킩��Ȃ��D�s��ɒ�����H
	//�Q	�d���Ȃ��̂ŁC���W�n����]�����ƍl����D
	//		(1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)�̐��K������ꂪ
	//		�ŗL�x�N�g���ƂȂ�悤�ȉ�]�s����ŏ����@��p���ċ��߂�



	////Test
	//Matrix3d testMatrix;
	//testMatrix <<	1.0,0.446,-0.56,

	//				0.446,1.0,-0.239,

	//				-0.56,0.239,1.0;

	//SelfAdjointEigenSolver<Matrix3d> es;
	//es.compute(testMatrix);

	//cout << "Input matrix:" << endl

 //      << testMatrix << endl

 //      << endl

 //      << "Eigenvalues:" << endl

 //      << es.eigenvalues() << endl

 //      << endl

 //      << "Eigenvectors:"<< endl

 //      << es.eigenvectors() << endl;

	//vector<pair<float, int>> eigenValVec;
	//eigenValVec.push_back(pair<float, int>(es.eigenvalues()[1], 1));
	//eigenValVec.push_back(pair<float, int>(es.eigenvalues()[2], 2));
	//eigenValVec.push_back(pair<float, int>(es.eigenvalues()[0], 0));

	//sort(eigenValVec.begin(), eigenValVec.end());

	////�e�X�g
	//Vector3f testVector = Vector3f::Zero();
	//testVector << 1, 2, 3;
	//Matrix3f testMatrix = Matrix3f::Zero();

	//testMatrix.row(0) += testVector(0) * testVector;
	//testMatrix.row(1) += testVector(1) * testVector;
	//testMatrix.row(2) += testVector(2) * testVector;

	//cout << "Test" << endl;
	//cout << testMatrix << endl;
}

void OrientedPrtObj::UpdateCluster()
{
	float dt = m_dDt;

	////CalcForce(dt);
	////CalcVelocity(dt);
	////DampVelocity(dt);
	ProjectConstraint(dt);
	DistanceConstraint(dt);
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
void OrientedPrtObj::CopyPrtToClstrPos(unsigned prtNum)
{
	int cIndx = 0;
	int pIndx = 0;
	int num = prtNum;

	//#pragma omp parallel for private(cIndx, pIndx)
	for(int i = 0; i < num; i++){
		cIndx = i * SM_DIM;
		pIndx = i * 4;

		s_pfClstrPos[cIndx+X] = s_pfPrtPos[pIndx+X];
		s_pfClstrVel[cIndx+X] = s_pfPrtVel[pIndx+X];

		s_pfClstrPos[cIndx+Y] = s_pfPrtPos[pIndx+Y];
		s_pfClstrVel[cIndx+Y] = s_pfPrtVel[pIndx+Y];

		s_pfClstrPos[cIndx+Z] = s_pfPrtPos[pIndx+Z];
		s_pfClstrVel[cIndx+Z] = s_pfPrtVel[pIndx+Z];
	}
}

void OrientedPrtObj::CalcClusterMomentMatrix(rxMatrix3& Apq, rxMatrix3& Aqq)
{
	rxMatrix3 Asum = rxMatrix3(0.0f);
	for(int i = 0; i < m_iIndxNum; ++i){
		if( CheckHole(i) ){	continue;	}

		OrientedParticle* prt = m_vOrientedPrtes[i];
		Asum += prt->A_elipsoid();
		Asum += prt->MomentMatrix();
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

	//Aqq
	Vec3 q(0.0f);

	for(int i = 0; i < m_iIndxNum; ++i){
		if( CheckHole(i) ){	continue;	}
		Vec3 op(m_pOrgPos[i*SM_DIM+0], m_pOrgPos[i*SM_DIM+1], m_pOrgPos[i*SM_DIM+2]);

		q = op - m_vec3OrgCm;
		double m = m_pMass[i];

		Aqq(0,0) += m*q[0]*q[0];
		Aqq(0,1) += m*q[0]*q[1];
		Aqq(0,2) += m*q[0]*q[2];
		Aqq(1,0) += m*q[1]*q[0];
		Aqq(1,1) += m*q[1]*q[1];
		Aqq(1,2) += m*q[1]*q[2];
		Aqq(2,0) += m*q[2]*q[0];
		Aqq(2,1) += m*q[2]*q[1];
		Aqq(2,2) += m*q[2]*q[2];
	}
}

void OrientedPrtObj::ProjectConstraint(float dt)
{
	//�X�V�ʒu���N���X�^�̊e���q�ɔ��f�@�e�N���X�^�̊e���q�͓��ꂳ���
	for(int i = 0; i < m_iIndxNum; ++i){
		if( CheckHole(i) ){	continue;	}

		int pIndx = m_iPIndxes[i];
		int cIndx = i * SM_DIM;
		Vec3 prdPos = m_vOrientedPrtes[i]->PrdPos();

		for(int j = 0; j < SM_DIM; j++){
			m_pCurPos[cIndx+j] = prdPos[j];
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
		Vec3 prdPos = m_vOrientedPrtes[i]->PrdPos();

		if(m_pCurPos[cIndx+0] < m_v3Min[0] || m_pCurPos[cIndx+0] > m_v3Max[0]){
			m_pCurPos[cIndx+0] = prdPos[0] - (s_pfClstrVel[cIndx+0] + (m_v3Gravity[0] * dt))*dt*res;
			m_pCurPos[cIndx+1] = prdPos[1];
			m_pCurPos[cIndx+2] = prdPos[2];
		}

		if(m_pCurPos[cIndx+1] < m_v3Min[1] || m_pCurPos[cIndx+1] > m_v3Max[1]){
			m_pCurPos[cIndx+1] = prdPos[1] - (s_pfClstrVel[cIndx+1] + (m_v3Gravity[1] * dt))*dt*res;
			m_pCurPos[cIndx+0] = prdPos[0] ;
			m_pCurPos[cIndx+2] = prdPos[2];
		}

		if(m_pCurPos[cIndx+2] < m_v3Min[2] || m_pCurPos[cIndx+2] > m_v3Max[2]){
			m_pCurPos[cIndx+2] = prdPos[2] - (s_pfClstrVel[cIndx+2] + (m_v3Gravity[2] * dt))*dt*res;
			m_pCurPos[cIndx+0] = prdPos[0];
			m_pCurPos[cIndx+1] = prdPos[1];
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

	rxMatrix3 A(0.0);
	A = Apq*Aqq.Inverse();	// A = Apq*Aqq^-1

	// �ڕW���W���v�Z���C���݂̒��_���W���ړ�
	m_fDefAmount = 0.0f;

	float alpha = 0.1f;
	float beta = 0.2f;

	rxMatrix3 RL = beta*A+(1.0f-beta)*R;

	for(int i = 0; i < m_iIndxNum; ++i){
		if( CheckHole(i) ){	continue;	}
		if(m_pFix[i]) continue;

		int cIndx = i*SM_DIM;

		Vec3 q;
		for(int j = 0; j < SM_DIM; j++){
			q[j] = m_pOrgPos[cIndx+j]-m_vec3OrgCm[j];
		}

		Vec3 gp(R*q+m_vec3NowCm);
		//Vec3 gp(RL*q+m_vec3NowCm);

		for(int j = 0; j < SM_DIM; j++){
			//float defAmount = (gp[j]-m_pCurPos[cIndx+j]) * alpha;
			float defAmount = (gp[j]-m_pCurPos[cIndx+j]);
			m_pCurPos[cIndx+j] += defAmount;

			m_fDefAmount += abs(defAmount);
		}
	}

	//�p���̂��߂ɉ�]�s���ۑ�
	m_vOrientedPrtes[0]->Rotation(R);
}

//��������
void OrientedPrtObj::DistanceConstraint(float dt)
{

}

void OrientedPrtObj::ApplyEstimatedPosition(float dt)
{
}

void OrientedPrtObj::CorrectVelocity(float dt)
{
}