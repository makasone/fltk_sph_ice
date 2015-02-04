/*! 
  @file rx_shape_matching.cpp
	
  @brief Shape Matching�@�ɂ��e���ό`
 
  @author Makoto Fujisawa
  @date 2013-07
*/


//-----------------------------------------------------------------------------
// �C���N���[�h�t�@�C��
//-----------------------------------------------------------------------------
#include "ShapeMatching.h"


//-----------------------------------------------------------------------------
// Shape Matching�N���X�̎���
//-----------------------------------------------------------------------------
/*!
 * �f�t�H���g�R���X�g���N�^
 */
rxShapeMatching::rxShapeMatching()
{
	m_pOrgPos = 0;
	m_pCurPos = 0;
	//m_pNewPos = new float[MAXPARTICLE*SM_DIM];
	//m_pGoalPos = new float[MAXPARTICLE*SM_DIM];
	
	m_pVel = 0;

	m_pMass = 0;
	m_pFix = 0;
}

rxShapeMatching::rxShapeMatching(int obj)
{
	//�p�����[�^�̏�����
	m_iObjectNo = obj;

	m_dDt = 0.01;

	m_v3Gravity = Vec3(0.0, -9.81/**0.005*/, 0.0);
	//m_v3Gravity = Vec3(0.0, 0.0, 0.0);

	m_v3Min = Vec3(-1.0);
	m_v3Max = Vec3(1.0);

	m_dAlpha = 1.0;	
	m_dBeta = 1.0;	// ���ꂪ�������ƍd���ގ��ɂȂ�
	
	m_bLinearDeformation = true;
	m_bVolumeConservation = true;

	m_fpCollision = 0;

	//�������m��
	//�N���X�^���ۑ��ł���ő嗱�q����MAXPARTICLE�Œ�`����D
	//���������s���ꍇ�̓R�����g�A�E�g
	m_pOrgPos = new float[MAXPARTICLE*SM_DIM];
	m_pCurPos = new float[MAXPARTICLE*SM_DIM];
	//m_pNewPos = new float[MAXPARTICLE*SM_DIM];
	//m_pGoalPos = new float[MAXPARTICLE*SM_DIM];
	
	m_pVel = new float[MAXPARTICLE*SM_DIM];

	m_pMass = new float[MAXPARTICLE];
	m_pFix = new bool[MAXPARTICLE];

	m_iPIndxes.resize(MAXPARTICLE, MAXINT);

	Clear();
}

/*!
 * �R�s�[�R���X�g���N�^
 */
rxShapeMatching::rxShapeMatching(const rxShapeMatching& copy)
{
	cout << __FUNCTION__ << endl;
	Copy(copy);
}

/*!
 * �f�X�g���N�^
 */
rxShapeMatching::~rxShapeMatching()
{
	ReleaseMemory();
}

//������Z�q�ŃR�s�[
rxShapeMatching& rxShapeMatching::operator=(const rxShapeMatching& copy)
{
	if(this != &copy){
		ReleaseMemory();

		Copy(copy);
	}

	return *this;
}


/*!
 * ���_�ʒu�̏�����
 */
void rxShapeMatching::initialize(void)
{
	for(int i = 0; i < m_iNumVertices; ++i){

		for(int j = 0; j < SM_DIM; j++)
		{
			m_pCurPos[i*SM_DIM+j] = m_pOrgPos[i*SM_DIM+j];
			//m_pNewPos[i*SM_DIM+j] = m_pOrgPos[i*SM_DIM+j];
			//m_pGoalPos[i*SM_DIM+j] = m_pOrgPos[i*SM_DIM+j];
			m_pVel[i*SM_DIM+j] = 0.0;
		}

		m_pFix[i] = false;
	}
}

/*!
 * �S���_��������
 */
void rxShapeMatching::Clear()
{
	m_iNumVertices = 0;

	for(int i = 0; i < MAXPARTICLE; i++)
	{
		for(int j = 0; j < SM_DIM; j++)
		{
			m_pOrgPos[i*SM_DIM+j] = 0.0;
			m_pCurPos[i*SM_DIM+j] = 0.0;
			//m_pNewPos[i*SM_DIM+j] = 0.0;
			//m_pGoalPos[i*SM_DIM+j] = 0.0;
			m_pVel[i*SM_DIM+j] = 0.0;
		}

		m_pMass[i] = 0.0;
		m_pFix[i] = false;

		m_iPIndxes[i] = MAXINT;
	}
}

//�������̊J��
void rxShapeMatching::ReleaseMemory()
{
	if(m_pOrgPos != 0){
		delete[] m_pOrgPos;
		m_pOrgPos = 0;
	}

	if(m_pCurPos != 0){
		delete[] m_pCurPos;
		m_pCurPos = 0;
	}

	if(m_pVel != 0){
		delete[] m_pVel;
		m_pVel = 0;
	}

	if(m_pMass != 0){
		delete[] m_pMass;
		m_pMass = 0;
	}

	if(m_pFix != 0){
		delete[] m_pFix;
		m_pFix = 0;
	}
}

//�f�[�^�R�s�[
void rxShapeMatching::Copy(const rxShapeMatching& copy)
{
	//�p�����[�^�̏�����
	m_iObjectNo = copy.objNo();
	m_dDt = copy.dt();
	m_v3Gravity = copy.gravity();

	m_v3Min = copy.areaMin();
	m_v3Max = copy.areaMax();

	m_dAlpha = copy.alpha();
	m_dBeta = copy.beta();
	
	m_bLinearDeformation = copy.isLiner();
	m_bVolumeConservation = copy.isVolumeConserve();

	m_fpCollision = 0;

	//�������m��
	//�N���X�^���ۑ��ł���ő嗱�q����MAXPARTICLE�Œ�`����D
	//���������s���ꍇ�̓R�����g�A�E�g
	m_pOrgPos = new float[MAXPARTICLE*SM_DIM];
	m_pCurPos = new float[MAXPARTICLE*SM_DIM];
	//m_pNewPos = new float[MAXPARTICLE*SM_DIM];
	//m_pGoalPos = new float[MAXPARTICLE*SM_DIM];
	
	m_pVel = new float[MAXPARTICLE*SM_DIM];

	m_pMass = new float[MAXPARTICLE];
	m_pFix = new bool[MAXPARTICLE];

	m_iPIndxes.resize(MAXPARTICLE);

	Clear();

	//�e�����R�s�[���ď�����
	for(int i = 0; i < MAXPARTICLE; i++){
		for(int j = 0; j < SM_DIM; j++){
			m_pOrgPos[i*SM_DIM+j] = copy.GetOrgPos(i)[j];
			m_pCurPos[i*SM_DIM+j] = copy.GetVertexPos(i)[j];

			m_pVel[i*SM_DIM+j] = copy.GetVertexVel(i)[j];
		}

		m_pMass[i] = copy.GetMass(i);
		m_pFix[i] = copy.IsFixed(i);

		m_iPIndxes[i] = copy.GetParticleIndx(i);
	}
}

/*!
 * ���_��ǉ�
 * @param[in] pos ���_�ʒu
 * @param[out] mass ���_����
 */
void rxShapeMatching::AddVertex(const Vec3 &pos, double mass, int pIndx)
{
	AddVertex(pos, Vec3(0.0), mass, pIndx);
}

/*!
 * ���_��ǉ�
 */
void rxShapeMatching::AddVertex(const Vec3 &pos, const Vec3& vel, double mass, int pIndx)
{
	AddVertex(pos, pos, vel, mass, pIndx);
}

/*!
 * ���_��ǉ�
 * @param[in] pos ���_�ʒu
 * @param[out] mass ���_����
 */
int rxShapeMatching::AddVertex(const Vec3& orgPos, const Vec3 &curPos, const Vec3& vel, double mass, int pIndx)
{
	//�J���Ă���ꏊ��T��
	int freeIndx = MAXINT;
	for(int i = 0; i < MAXPARTICLE; i++){
		if(m_iPIndxes[i] == MAXINT){
			freeIndx = i;
			break;
		}
	}

	if(freeIndx == MAXINT){
		cerr << "Cant Add" << endl;
		return -1;
	}

	for(int i = 0; i < SM_DIM; i++){
		m_pOrgPos[freeIndx*SM_DIM+i] = orgPos[i];
		m_pCurPos[freeIndx*SM_DIM+i] = curPos[i];

		m_pVel[freeIndx*SM_DIM+i] = vel[i];
	}

	//m_pFix�͓��ɂ��邱�Ɩ���
	m_pMass[freeIndx] = mass;

	m_iPIndxes[freeIndx] = pIndx;
	//if(m_iPIndxes.size() >= 2)
	//{
	//	sort(m_iPIndxes.begin(), m_iPIndxes.end());	//�\�[�g
	//}

	m_iNumVertices++;

	return freeIndx;
}

//���_�폜
void rxShapeMatching::Remove(int indx)
{
	if(m_iNumVertices <= 0) return;
	m_iNumVertices--;

	m_iPIndxes[indx] = MAXINT;

	for(int i = 0; i < SM_DIM; i++)
	{
		m_pOrgPos[indx*SM_DIM+i] = 0.0;
		m_pCurPos[indx*SM_DIM+i] = 0.0;
		m_pVel[indx*SM_DIM+i] = 0.0;
	}

	m_pMass[indx] = 0.0;
	m_pFix[indx] = false;
}

/*!
 * �O��
 *  - �d�͂Ƌ��E�ǂ���̗͂̉e��
 * @param[in] dt �^�C���X�e�b�v��
 */
void rxShapeMatching::calExternalForces(double dt)
{
	// �d�͂̉e����t��
	for(int i = 0; i < m_iNumVertices; ++i){
		if(m_pFix[i]) continue;
		for(int j = 0; j < SM_DIM; j++)
		{
			m_pVel[i*SM_DIM+j] += m_v3Gravity[j]*dt;
			m_pNewPos[i*SM_DIM+j] = m_pCurPos[i*SM_DIM+j]+m_pVel[i*SM_DIM+j]*dt;
			m_pGoalPos[i*SM_DIM+j] = m_pOrgPos[i*SM_DIM+j];
		}
	}

	// ���E�ǂ̉e��
	double res = 0.5;	// �����W��
	for(int i = 0; i < m_iNumVertices; ++i){
		if(m_pFix[i]) continue;

		if(m_pNewPos[i*SM_DIM+0] < m_v3Min[0] || m_pNewPos[i*SM_DIM+0] > m_v3Max[0])
		{
			m_pNewPos[i*SM_DIM+0] = m_pCurPos[i*SM_DIM+0]-m_pVel[i*SM_DIM+0]*dt*res;
			m_pNewPos[i*SM_DIM+1] = m_pCurPos[i*SM_DIM+1];
			m_pNewPos[i*SM_DIM+2] = m_pCurPos[i*SM_DIM+2];
		}
		if(m_pNewPos[i*SM_DIM+1] < m_v3Min[1] || m_pNewPos[i*SM_DIM+1] > m_v3Max[1])
		{
			m_pNewPos[i*SM_DIM+1] = m_pCurPos[i*SM_DIM+1]-m_pVel[i*SM_DIM+1]*dt*res;
			m_pNewPos[i*SM_DIM+0] = m_pCurPos[i*SM_DIM+0];
			m_pNewPos[i*SM_DIM+2] = m_pCurPos[i*SM_DIM+2];
		}
		if(m_pNewPos[i*SM_DIM+2] < m_v3Min[2] || m_pNewPos[i*SM_DIM+2] > m_v3Max[2])
		{
			m_pNewPos[i*SM_DIM+2] = m_pCurPos[i*SM_DIM+2]-m_pVel[i*SM_DIM+2]*dt*res;
			m_pNewPos[i*SM_DIM+0] = m_pCurPos[i*SM_DIM+0];
			m_pNewPos[i*SM_DIM+1] = m_pCurPos[i*SM_DIM+1];
		}

		Vec3 np(m_pNewPos[i*SM_DIM+0], m_pNewPos[i*SM_DIM+1], m_pNewPos[i*SM_DIM+2]);
		clamp(np);

		m_pNewPos[i*SM_DIM+0] = np[0];
		m_pNewPos[i*SM_DIM+1] = np[1];
		m_pNewPos[i*SM_DIM+2] = np[2];
	}
}

void rxShapeMatching::calCollision(double dt)
{
	// ���̃I�u�W�F�N�g�Ƃ̏Փ�
	//if(m_fpCollision != 0){
	//	for(int i = 0; i < m_iNumVertices; ++i){
	//		if(m_pFix[i]) continue;
	//		Vec3 &p = m_pCurPos[i];
	//		Vec3 &np = m_pNewPos[i];
	//		Vec3 &v = m_pVel[i];
	//		m_fpCollision(p, np, v, m_iObjectNo);
	//	}
	//}
}

/*!
 * Shape Matching�@
 *  - �ڕW�ʒu���v�Z���āCm_pNewPos�����̈ʒu�ɋ߂Â���
 * @param[in] dt �^�C���X�e�b�v��
 */
void rxShapeMatching::shapeMatching(double dt)
{
	if(m_iNumVertices <= 1) return;

	Vec3 cm(0.0), cm_org(0.0);	// �d�S
	double mass = 0.0f;	// ������

	// �d�S���W�̌v�Z
	for(int i = 0; i < m_iNumVertices;++i){
		double m = m_pMass[i];
		if(m_pFix[i]) m *= 300.0;	// �Œ�_�̎��ʂ�傫������
		mass += m;

		Vec3 np(m_pNewPos[i*SM_DIM+0], m_pNewPos[i*SM_DIM+1], m_pNewPos[i*SM_DIM+3]);
		Vec3 op(m_pOrgPos[i*SM_DIM+0], m_pOrgPos[i*SM_DIM+1], m_pOrgPos[i*SM_DIM+3]);

		cm += np*m;
		cm_org += op*m;
	}

	cm /= mass;
	cm_org /= mass;

	rxMatrix3 Apq(0.0), Aqq(0.0);
	Vec3 p, q;

	// Apq = ��mpq^T
	// Aqq = ��mqq^T
	for(int i = 0; i < m_iNumVertices; ++i){

		Vec3 np(m_pNewPos[i*SM_DIM+0], m_pNewPos[i*SM_DIM+1], m_pNewPos[i*SM_DIM+3]);
		Vec3 op(m_pOrgPos[i*SM_DIM+0], m_pOrgPos[i*SM_DIM+1], m_pOrgPos[i*SM_DIM+3]);

		p = np-cm;
		q = op-cm_org;
		double m = m_pMass[i];

		Apq(0,0) += m*p[0]*q[0];
		Apq(0,1) += m*p[0]*q[1];
		Apq(0,2) += m*p[0]*q[2];
		Apq(1,0) += m*p[1]*q[0];
		Apq(1,1) += m*p[1]*q[1];
		Apq(1,2) += m*p[1]*q[2];
		Apq(2,0) += m*p[2]*q[0];
		Apq(2,1) += m*p[2]*q[1];
		Apq(2,2) += m*p[2]*q[2];

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

	rxMatrix3 R, S;
	PolarDecomposition(Apq, R, S);

	if(m_bLinearDeformation){
		// Linear Deformations
		rxMatrix3 A;
		A = Apq*Aqq.Inverse();	// A = Apq*Aqq^-1

		// �̐ϕۑ��̂��߂Ɂ�(det(A))�Ŋ���
		if(m_bVolumeConservation){
			double det = fabs(A.Determinant());
			if(det > RX_FEQ_EPS){
				det = 1.0/sqrt(det);
				if(det > 2.0) det = 2.0;
				A *= det;
			}
		}

		// ��]�s��R�̑���̍s��RL=��A+(1-��)R���v�Z
		rxMatrix3 RL = m_dBeta*A+(1.0-m_dBeta)*R;

		// �ڕW���W���v�Z���C���݂̒��_���W���ړ�
		for(int i = 0; i < m_iNumVertices; ++i){
			if(m_pFix[i]) continue;

			Vec3 np(m_pNewPos[i*SM_DIM+0], m_pNewPos[i*SM_DIM+1], m_pNewPos[i*SM_DIM+3]);
			Vec3 op(m_pOrgPos[i*SM_DIM+0], m_pOrgPos[i*SM_DIM+1], m_pOrgPos[i*SM_DIM+3]);

			q = op-cm_org;
			Vec3 gp = RL*q+cm;

			for(int j = 0; j < SM_DIM; j++)
			{
				m_pGoalPos[i*SM_DIM+j] = gp[j];
				m_pNewPos[i*SM_DIM+j] += (m_pGoalPos[i*SM_DIM+j]-np[j])*m_dAlpha;
			}
		}
	}
	else{
		//// Quadratic Deformations
		//double Atpq[3][9];
		//for(int j = 0; j < 9; ++j){
		//	Atpq[0][j] = 0.0;
		//	Atpq[1][j] = 0.0;
		//	Atpq[2][j] = 0.0;
		//}
		//rxMatrixN<double,9> Atqq;
		//Atqq.SetValue(0.0);

		//for(int i = 0; i < m_iNumVertices; ++i){
		//	p = m_pNewPos[i]-cm;
		//	q = m_pOrgPos[i]-cm_org;

		//	// q~�̌v�Z
		//	double qt[9];
		//	qt[0] = q[0];      qt[1] = q[1];      qt[2] = q[2];
		//	qt[3] = q[0]*q[0]; qt[4] = q[1]*q[1]; qt[5] = q[2]*q[2];
		//	qt[6] = q[0]*q[1]; qt[7] = q[1]*q[2]; qt[8] = q[2]*q[0];

		//	// A~pq = ��mpq~ �̌v�Z
		//	double m = m_vMass[i];
		//	for(int j = 0; j < 9; ++j){
		//		Atpq[0][j] += m*p[0]*qt[j];
		//		Atpq[1][j] += m*p[1]*qt[j];
		//		Atpq[2][j] += m*p[2]*qt[j];
		//	}

		//	// A~qq = ��mq~q~ �̌v�Z
		//	for(int j = 0; j < 9; ++j){
		//		for(int k = 0; k < 9; ++k){
		//			Atqq(j,k) += m*qt[j]*qt[k];
		//		}
		//	}
		//}

		//// A~qq�̋t�s��Z�o
		//Atqq.Invert();

		//double At[3][9];
		//for(int i = 0; i < 3; ++i){
		//	for(int j = 0; j < 9; j++){
		//		At[i][j] = 0.0f;
		//		for(int k = 0; k < 9; k++){
		//			At[i][j] += Atpq[i][k]*Atqq(k,j);
		//		}

		//		// ��A~+(1-��)R~
		//		At[i][j] *= m_dBeta;
		//		if(j < 3){
		//			At[i][j] += (1.0f-m_dBeta)*R(i,j);
		//		}
		//	}
		//}

		// // a00a11a22+a10a21a02+a20a01a12-a00a21a12-a20a11a02-a10a01a22
		//double det = At[0][0]*At[1][1]*At[2][2]+At[1][0]*At[2][1]*At[0][2]+At[2][0]*At[0][1]*At[1][2]
		//			-At[0][0]*At[2][1]*At[1][2]-At[2][0]*At[1][1]*At[0][2]-At[1][0]*At[0][1]*At[2][2];

		//// �̐ϕۑ��̂��߂Ɂ�(det(A))�Ŋ���
		//if(m_bVolumeConservation){
		//	if(det != 0.0f){
		//		det = 1.0f/sqrt(fabs(det));
		//		if(det > 2.0f) det = 2.0f;
		//		for(int i = 0; i < 3; ++i){
		//			for(int j = 0; j < 3; ++j){
		//				At[i][j] *= det;
		//			}
		//		}
		//	}
		//}


		//// �ڕW���W���v�Z���C���݂̒��_���W���ړ�
		//for(int i = 0; i < m_iNumVertices; ++i){
		//	if(m_pFix[i]) continue;
		//	q = m_pOrgPos[i]-cm_org;

		//	for(int k = 0; k < 3; ++k){
		//		m_pGoalPos[i][k] = At[k][0]*q[0]+At[k][1]*q[1]+At[k][2]*q[2]
		//						  +At[k][3]*q[0]*q[0]+At[k][4]*q[1]*q[1]+At[k][5]*q[2]*q[2]+
		//						  +At[k][6]*q[0]*q[1]+At[k][7]*q[1]*q[2]+At[k][8]*q[2]*q[0];
		//	}

		//	m_pGoalPos[i] += cm;
		//	m_pNewPos[i] += (m_pGoalPos[i]-m_pNewPos[i])*m_dAlpha;
		//}

	}
}

/*!
 * ���x�ƈʒu�̍X�V
 *  - �V�����ʒu�ƌ��݂̈ʒu���W���瑬�x���Z�o
 * @param[in] dt �^�C���X�e�b�v��
 */
void rxShapeMatching::integrate(double dt)
{
	double dt1 = 1.0f/dt;
	for(int i = 0; i < m_iNumVertices; ++i){
		for(int j = 0; j < SM_DIM; j++)
		{
			m_pVel[i*SM_DIM+j] = (m_pNewPos[i*SM_DIM+j]-m_pCurPos[i*SM_DIM+j])*dt1;
			m_pCurPos[i*SM_DIM+j] = m_pNewPos[i*SM_DIM+j];
		}
	}
}

/*!
 * �V�~�����[�V�����X�e�b�v��i�߂�
 */
void rxShapeMatching::Update()
{	
	calExternalForces(m_dDt);
	calCollision(m_dDt);
	shapeMatching(m_dDt);
	integrate(m_dDt);
}

/*!
 * �Œ蒸�_��ݒ�
 * @param[in] i ���_�C���f�b�N�X
 * @param[in] pos �Œ�ʒu
 */
void rxShapeMatching::FixVertex(int i, const Vec3 &pos)
{
	m_pCurPos[i*SM_DIM+0] = pos[0];
	m_pCurPos[i*SM_DIM+1] = pos[1];
	m_pCurPos[i*SM_DIM+2] = pos[2];

	m_pFix[i] = true;
}

/*!
 * ���_�̌Œ������
 * @param[in] i ���_�C���f�b�N�X
 */
void rxShapeMatching::UnFixVertex(int i)
{
	m_pFix[i] = false;
}


