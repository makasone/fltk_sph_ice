/*!
  @file Ice_OrientedParticle.h
	
  @brief Oriented Particle�@�ɂ��e���̕ό`                                                                                                                                                                                                                                                                                                                                                                                                                                                        
  @ref M. Muller et al., "Solid Simulation with Oriented Particles", SIGGRAPH2011. 
  �Q�l�ɂ����y�[�W�C�v���O���� http://yuki-koyama.hatenablog.com/entry/2015/01/30/150419
 
  @author Ryo Nakasone
  @date 2015-4
*/
#include "Ice_OrientedParticle.h"

typedef OrientedParticleBaseElasticObject OrientedPrtObj;

const float* OrientedPrtObj::s_pfPrtPos = 0;
const float* OrientedPrtObj::s_pfPrtVel = 0;

float* OrientedPrtObj::s_pfClstrPos = 0;
float* OrientedPrtObj::s_pfClstrVel = 0;

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
	for(unsigned i = 0; i < prtNum; i++){
		int pIndx = i * SM_DIM;

		for(unsigned j = 0; j < SM_DIM; j++){
			s_pfClstrPos[pIndx+j] = s_pfPrtPos[i*4+j];
			s_pfClstrVel[pIndx+j] = s_pfPrtVel[i*4+j];
		}
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

	//�p���̂��߂ɉ�]�s���ۑ�
	m_vOrientedPrtes[0]->Rotation(R);
}

void OrientedPrtObj::ApplyEstimatedPosition(float dt)
{
}

void OrientedPrtObj::CorrectVelocity(float dt)
{
}