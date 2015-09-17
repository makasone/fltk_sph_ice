/*!
  @file ElasticObject_OP.h
	
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

#include "ElasticObject_OP.h"

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

	//�P	�����̃x�N�g������p�����`������@���킩��Ȃ��D�Q����Ȃ��`�ł���݂��������ǁc�@�s��ɒ�����H
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

void OrientedPrtObj::UpdateCluster_OP()
{
	float dt = m_dDt;

	////CalcForce(dt);
	////CalcVelocity(dt);
	////DampVelocity(dt);
	ProjectConstraint(dt);
	//DistanceConstraint(dt);
	////ApplyEstimatedPosition(dt);
	////CorrectVelocity(dt);

	////�]���̃V�F�C�v�}�b�`���O
	//ShapeMatchingNormal();
}

//���܂������Ȃ�����
void OrientedPrtObj::UpdateCluster_Sampling(const IceStructure* ice_struct)
{
	float dt = m_dDt;

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
	
	////�ߖT�N���X�^�̉�]�s����Ԃ��Ď��g�̉�]�s������߂�
	//rxMatrix3 R_intrp = InterpolateRotation(ice_struct);

	//rxMatrix3 Apq(0.0), Aqq(0.0);
	//CalcClusterMomentMatrix(Apq, Aqq);		//�p�[�e�B�N���̎p�����l���������[�����g�}�g���b�N�X

	//rxMatrix3 R_polor, S;
	//PolarDecomposition(Apq, R_polor, S);

	////if(m_iObjectNo < 5){
	////	cout << m_iObjectNo << endl;
	////	cout << "R_intrp:\n" << R_intrp << endl;
	////	cout << "R_polor:\n" << R_polor << endl;
	////}

	//rxMatrix3 R = R_intrp;

	//�ߖT�N���X�^�̕ό`���z�e���\������Apq���ߎ�
	rxMatrix3 Apq = InterpolateApq(ice_struct);
	rxMatrix3 R = Apq;

	//rxMatrix3 R_polor, S;
	//PolarDecomposition(Apq, R_polor, S);
	//rxMatrix3 R = R_polor;

	//��Ԃ�����]�s���p���ĉ^���v�Z
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

//�N���X�^�Ɋ܂܂�Ă��闱�q�̎p���E��]�s�񂩂玩�g�̉�]�s�����
//�p�����Ԃ��邾���ł悢��������Ȃ����C�ꉞ��]�s��ŋ��߂�
//���܂������Ȃ�����
rxMatrix3 OrientedPrtObj::InterpolateRotation(const IceStructure* ice_struct)
{
	float weight = 0.5f;	//�܂��͕���
	mk_Quaternion tmp_q = m_vOrientedPrtes[0]->CurOrientation();
	Quaternionf myQuat = Quaternionf(tmp_q.w, tmp_q.x, tmp_q.y, tmp_q.z);
	//Quaternionf myQuat = Quaternionf::Identity();
	rxMatrix3 R = rxMatrix3::Identity();

	//i = 0�͎����Ȃ̂ŏ���
	for(int i = 1; i < m_iIndxNum; i++){
		if(CheckHole(i)){	continue;	}
		int pIndx = m_iPIndxes[i];
		if(ice_struct->GetMotionCalcCluster(pIndx) == 0){	continue;	}

		rxMatrix3 rotateMtrx = m_vOrientedPrtes[i]->Rotation();

		//��]�s�񁨃N�H�[�^�j�I��
		Matrix<float, 3, 3, RowMajor> tmp_r1;
		tmp_r1 <<	rotateMtrx(0, 0), rotateMtrx(0, 1), rotateMtrx(0, 2), 
					rotateMtrx(1, 0), rotateMtrx(1, 1), rotateMtrx(1, 2), 
					rotateMtrx(2, 0), rotateMtrx(2, 1), rotateMtrx(2, 2);

		Quaternionf rotateQuat(tmp_r1);

		//�N�H�[�^�j�I�����m�ŋ��ʐ��`���
		Quaternionf rotateSlerp = myQuat.slerp(weight, rotateQuat);

		//�N�H�[�^�j�I������]�s��
		Matrix<float, 3, 3, RowMajor> tmp_r2 = rotateSlerp.matrix();
		rxMatrix3 R_Slerp(
			tmp_r2(0, 0), tmp_r2(0, 1), tmp_r2(0, 2), 
			tmp_r2(1, 0), tmp_r2(1, 1), tmp_r2(1, 2), 
			tmp_r2(2, 0), tmp_r2(2, 1), tmp_r2(2, 2)
		);

		//���όv�Z
		R *= R_Slerp;
	}

	return R;
}

//�ߖT�N���X�^�̕ό`���z�e���\����p���Ă��̃N���X�^�̕ό`���z�e���\�����ߎ�
//���܂������Ȃ�����
rxMatrix3 OrientedPrtObj::InterpolateApq(const IceStructure* ice_struct)
{
	float weight = 0.5f;	//�Ă��Ƃ��ɕ���

	//��]����
	rxMatrix3 R = rxMatrix3::Identity();

	mk_Quaternion tmp_q = m_vOrientedPrtes[0]->CurOrientation();
	//Quaternionf myQuat = Quaternionf(tmp_q.w, tmp_q.x, tmp_q.y, tmp_q.z);

	Quaternionf myQuat = Quaternionf::Identity();

	//i = 0�͎����Ȃ̂ŏ���
	for(int i = 1; i < m_iIndxNum; i++){
		if(CheckHole(i)){	continue;	}
		int pIndx = m_iPIndxes[i];
		if(ice_struct->GetMotionCalcCluster(pIndx) == 0){	continue;	}

		rxMatrix3 rotateMtrx = m_vOrientedPrtes[i]->Rotation();

		//��]�s�񁨃N�H�[�^�j�I��
		Matrix<float, 3, 3, RowMajor> tmp_r1;
		tmp_r1 <<	rotateMtrx(0, 0), rotateMtrx(0, 1), rotateMtrx(0, 2), 
					rotateMtrx(1, 0), rotateMtrx(1, 1), rotateMtrx(1, 2), 
					rotateMtrx(2, 0), rotateMtrx(2, 1), rotateMtrx(2, 2);
		
		Quaternionf rotateQuat(tmp_r1);

		//�N�H�[�^�j�I�����m�ŋ��ʐ��`���
		Quaternionf rotateSlerp = myQuat.slerp(weight, rotateQuat);

		//�N�H�[�^�j�I������]�s��
		Matrix<float, 3, 3, RowMajor> tmp_r2 = rotateSlerp.matrix();
		rxMatrix3 R_Slerp(
			tmp_r2(0, 0), tmp_r2(0, 1), tmp_r2(0, 2), 
			tmp_r2(1, 0), tmp_r2(1, 1), tmp_r2(1, 2), 
			tmp_r2(2, 0), tmp_r2(2, 1), tmp_r2(2, 2)
		);

		//���όv�Z
		R *= R_Slerp;
	}

	//��]�ȊO�̐���
	Matrix<float, 3, 3, RowMajor> S_temp;
	S_temp <<	0.0f, 0.0f, 0.0f,
				0.0f, 0.0f, 0.0f,
				0.0f, 0.0f, 0.0f;

	for(int i = 1; i < m_iIndxNum; i++){
		if(CheckHole(i)){	continue;	}
		int pIndx = m_iPIndxes[i];
		if(ice_struct->GetMotionCalcCluster(pIndx) == 0){	continue;	}

		rxMatrix3 sym = m_vOrientedPrtes[i]->Symmetric();
		
		Matrix<float, 3, 3, RowMajor> eigen_Sim;
		eigen_Sim <<	sym(0, 0), sym(0, 1), sym(0, 2), 
						sym(1, 0), sym(1, 1), sym(1, 2), 
						sym(2, 0), sym(2, 1), sym(2, 2);

		//�w���s����Ă���ł����̂��H
		Matrix<float, 3, 3, RowMajor> log_sim = eigen_Sim.array().log();
		S_temp += weight * log_sim;
	}
	//�ΐ��s����Ă���ł����̂��H
	S_temp = S_temp.array().exp();

	rxMatrix3 S(
			S_temp(0, 0), S_temp(0, 1), S_temp(0, 2), 
			S_temp(1, 0), S_temp(1, 1), S_temp(1, 2), 
			S_temp(2, 0), S_temp(2, 1), S_temp(2, 2)
	);

	rxMatrix3 Apq =R * S;

	return Apq;
}


//�V�F�C�v�}�b�`���O�@
void OrientedPrtObj::UpdateCluster_SM()
{
	double dt = m_dDt;

	//�ŏI�ʒu�E���x���N���X�^�̊e���q�ɔ��f�@�e�N���X�^�̊e���q�͓��ꂳ���
	for(int i = 0; i < m_iIndxNum; ++i)
	{
		if( CheckHole(i) ){	continue;	}
		int pIndx = m_iPIndxes[i]*4;
		int cIndx = i*SM_DIM;
		int indx = m_iPIndxes[i]*3;


		Vec3 prdPos = m_vOrientedPrtes[i]->PrdPos();
		Vec3 vel = m_vOrientedPrtes[i]->Velocity();

		for(int j = 0; j < SM_DIM; j++)
		{
			int jpIndx = pIndx+j;
			int jcIndx = cIndx+j;

			m_pCurPos[jcIndx] = s_pfPrtPos[jpIndx] + (s_pfPrtVel[jpIndx] + (m_v3Gravity[j] * dt)) *dt;
			//m_pCurPos[jcIndx] = s_pfPrtPos[jpIndx] + (s_pfClstrVel[indx+j] + (m_v3Gravity[j] * dt)) *dt;	//�����Ȃ�
			//m_pCurPos[jcIndx] = s_pfPrtPos[jpIndx] + (vel[j] + (m_v3Gravity[j] * dt)) *dt;	//�Ȃ����ł��Ȃ�
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
		Vec3 vel = m_vOrientedPrtes[i]->Velocity();

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

	m_vOrientedPrtes[0]->Rotation(R);

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

void OrientedPrtObj::UpdateCluster_SM_Itr()
{
	double dt = m_dDt;

	//// ���E�ǂ̉e��
	////���������Ȃ�d���Ȃ邪�C����͂���݂���
	//double res = 0.9;	// �����W��
	//for(int i = 0; i < m_iIndxNum; ++i)
	//{
	//	if( CheckHole(i) ){	continue;	}
	//	if(m_pFix[i]) continue;

	//	int pIndx =  m_iPIndxes[i]*4;
	//	int cIndx = i*SM_DIM;

	//	if(m_pCurPos[cIndx+0] < m_v3Min[0] || m_pCurPos[cIndx+0] > m_v3Max[0]){
	//		m_pCurPos[cIndx+0] = s_pfPrtPos[pIndx+0] - (s_pfPrtVel[pIndx+0] + (m_v3Gravity[0] * dt))*dt*res;
	//		m_pCurPos[cIndx+1] = s_pfPrtPos[pIndx+1];
	//		m_pCurPos[cIndx+2] = s_pfPrtPos[pIndx+2];
	//	}

	//	if(m_pCurPos[cIndx+1] < m_v3Min[1] || m_pCurPos[cIndx+1] > m_v3Max[1]){
	//		m_pCurPos[cIndx+1] = s_pfPrtPos[pIndx+1] - (s_pfPrtVel[pIndx+1] + (m_v3Gravity[1] * dt))*dt*res;
	//		m_pCurPos[cIndx+0] = s_pfPrtPos[pIndx+0] ;
	//		m_pCurPos[cIndx+2] = s_pfPrtPos[pIndx+2];
	//	}

	//	if(m_pCurPos[cIndx+2] < m_v3Min[2] || m_pCurPos[cIndx+2] > m_v3Max[2]){
	//		m_pCurPos[cIndx+2] = s_pfPrtPos[pIndx+2] - (s_pfPrtVel[pIndx+2] + (m_v3Gravity[2] * dt))*dt*res;
	//		m_pCurPos[cIndx+0] = s_pfPrtPos[pIndx+0];
	//		m_pCurPos[cIndx+1] = s_pfPrtPos[pIndx+1];
	//	}

	//	clamp(m_pCurPos, cIndx);
	//}

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

	if(m_iNumVertices <= 1) return;

	Vec3 p, q;
	Vec3 cm(0.0), cm_org(0.0);	// �d�S
	double mass = 0.0;	// ������
	rxMatrix3 R, S;

	CalcNowCm();							//�N���X�^�̏d�S���X�V

	if(m_iNumVertices <= 1) return;

	rxMatrix3 Apq(0.0), Aqq(0.0);

	cm = m_vec3NowCm;
	cm_org = m_vec3OrgCm;

	// Apq = ��mpq^T
	// Aqq = ��mqq^T
	for(int i = 0; i < m_iIndxNum; ++i)
	{
		if( CheckHole(i) ){	continue;	}

		int pIndx = m_iPIndxes[i] * SM_DIM;
		int cIndx = i*SM_DIM;

		Vec3 prdPos = m_vOrientedPrtes[i]->PrdPos();

		for(int j = 0; j < SM_DIM; j++)
		{
			p[j] = /*prdPos[j]*/m_pCurPos[cIndx+j]-cm[j];
			q[j] = m_pOrgPos[cIndx+j]-cm_org[j];
		}

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

		//Aqq(0,0) += m*q[0]*q[0];
		//Aqq(0,1) += m*q[0]*q[1];
		//Aqq(0,2) += m*q[0]*q[2];
		//Aqq(1,0) += m*q[1]*q[0];
		//Aqq(1,1) += m*q[1]*q[1];
		//Aqq(1,2) += m*q[1]*q[2];
		//Aqq(2,0) += m*q[2]*q[0];
		//Aqq(2,1) += m*q[2]*q[1];
		//Aqq(2,2) += m*q[2]*q[2];
	}

	//Apq�̍s�񎮂����߁C���]���邩�𔻒�
	//TODO:�s����ȏꍇ�������̂Ł~
	if(Apq.Determinant() < 0.0 && m_iNumVertices >= 10){
		//�P�@�����𔽓]
		Apq(0,2) = -Apq(0,2);
		Apq(1,2) = -Apq(1,2);
		Apq(2,2) = -Apq(2,2);
	}

	//PolarDecomposition(Apq, R, S, m_mtrxBeforeU); //warm start
	PolarDecomposition(Apq, R, S);

	if(m_bLinearDeformation)
	{
		//// Linear Deformations
		//rxMatrix3 A;
		//A = Apq*Aqq.Inverse();	// A = Apq*Aqq^-1

		//// �̐ϕۑ��̂��߂Ɂ�(det(A))�Ŋ���
		//if(m_bVolumeConservation){
		//	double det = fabs(A.Determinant());
		//	if(det > RX_FEQ_EPS){
		//		det = 1.0/sqrt(det);
		//		if(det > 2.0) det = 2.0;
		//		A *= det;
		//	}
		//}

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

	//���E����
	for(int i = 0; i < m_iIndxNum; ++i)
	{
		clamp(m_pCurPos, i*SM_DIM);
	}

	m_vOrientedPrtes[0]->Rotation(R);
}

//��������ɂ��^���v�Z
void OrientedPrtObj::UpdateCluster_DC()
{
	double dt = m_dDt;

	//�ŏI�ʒu�E���x���N���X�^�̊e���q�ɔ��f�@�e�N���X�^�̊e���q�͓��ꂳ���
	for(int i = 0; i < m_iIndxNum; ++i)
	{
		if( CheckHole(i) ){	continue;	}
		int pIndx = m_iPIndxes[i]*4;
		int cIndx = i*SM_DIM;
		int indx = m_iPIndxes[i]*3;


		Vec3 prdPos = m_vOrientedPrtes[i]->PrdPos();
		Vec3 vel = m_vOrientedPrtes[i]->Velocity();

		for(int j = 0; j < SM_DIM; j++)
		{
			int jpIndx = pIndx+j;
			int jcIndx = cIndx+j;

			m_pCurPos[jcIndx] = s_pfPrtPos[jpIndx] + (s_pfPrtVel[jpIndx] + (m_v3Gravity[j] * dt)) *dt;
			//m_pCurPos[jcIndx] = s_pfPrtPos[jpIndx] + (s_pfClstrVel[indx+j] + (m_v3Gravity[j] * dt)) *dt;	//�����Ȃ�
			//m_pCurPos[jcIndx] = s_pfPrtPos[jpIndx] + (vel[j] + (m_v3Gravity[j] * dt)) *dt;	//�Ȃ����ł��Ȃ�
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
		Vec3 vel = m_vOrientedPrtes[i]->Velocity();

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

	
	Vec3 cm(0.0), cm_org(0.0);	// �d�S
	double mass = 0.0;	// ������
	
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

	//��������
	Vec3 mainPrtPos_cur = GetVertexPos(0);
	Vec3 mainPrtPos_org = GetOrgPos(0);

	//�d�S�Ƃ̋�������@�d�S�ʒu�͌Œ�
	Vec3 nowDx_center = m_vec3NowCm - mainPrtPos_cur;
	float nowDist_center = norm(nowDx_center);

	Vec3 orgDx_center = m_vec3OrgCm - mainPrtPos_org;
	float orgDist_center = norm(orgDx_center);

	Vec3 dirVec_center = nowDist_center > 0.0f ? nowDx_center / nowDist_center : Vec3(0.0f);

	mainPrtPos_cur += dirVec_center * (orgDist_center - nowDist_center);
	SetCurrentPos(0, mainPrtPos_cur);

	m_fDefAmount = 0.0f;

	//��������@�N���X�^�ɑΉ����闱�q�ƁC���̑��̗��q
	for(int i = 1; i < m_iIndxNum; i++){
		if(CheckHole(i)){	continue;	}
		if(m_pFix[i]){		continue;	}

		Vec3 vertexPos = GetVertexPos(i);

		Vec3 nowDx = vertexPos - mainPrtPos_cur;
		float nowDist = norm(nowDx);

		Vec3 orgDx = GetOrgPos(i) - mainPrtPos_org;
		float orgDist = norm(orgDx);

		Vec3 dirVec = nowDist > 0.0f ? nowDx / nowDist : Vec3(0.0f);
		Vec3 dx = 1.0f/2.0f * dirVec * (orgDist - nowDist);

		mainPrtPos_cur -= dx;
		SetCurrentPos(i, vertexPos + dx);

		m_fDefAmount += abs(orgDist - nowDist);
	}

	SetCurrentPos(0, mainPrtPos_cur);

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

//���������@��������ɂ��^���v�Z
void OrientedPrtObj::UpdateCluster_DC_Itr()
{
	double dt = m_dDt;

	//// ���E�ǂ̉e��
	////���������Ȃ�d���Ȃ邪�C����͂���݂���
	//double res = 0.9;	// �����W��
	//for(int i = 0; i < m_iIndxNum; ++i)
	//{
	//	if( CheckHole(i) ){	continue;	}
	//	if(m_pFix[i]) continue;

	//	int pIndx =  m_iPIndxes[i]*4;
	//	int cIndx = i*SM_DIM;

	//	if(m_pCurPos[cIndx+0] < m_v3Min[0] || m_pCurPos[cIndx+0] > m_v3Max[0]){
	//		m_pCurPos[cIndx+0] = s_pfPrtPos[pIndx+0] - (s_pfPrtVel[pIndx+0] + (m_v3Gravity[0] * dt))*dt*res;
	//		m_pCurPos[cIndx+1] = s_pfPrtPos[pIndx+1];
	//		m_pCurPos[cIndx+2] = s_pfPrtPos[pIndx+2];
	//	}

	//	if(m_pCurPos[cIndx+1] < m_v3Min[1] || m_pCurPos[cIndx+1] > m_v3Max[1]){
	//		m_pCurPos[cIndx+1] = s_pfPrtPos[pIndx+1] - (s_pfPrtVel[pIndx+1] + (m_v3Gravity[1] * dt))*dt*res;
	//		m_pCurPos[cIndx+0] = s_pfPrtPos[pIndx+0] ;
	//		m_pCurPos[cIndx+2] = s_pfPrtPos[pIndx+2];
	//	}

	//	if(m_pCurPos[cIndx+2] < m_v3Min[2] || m_pCurPos[cIndx+2] > m_v3Max[2]){
	//		m_pCurPos[cIndx+2] = s_pfPrtPos[pIndx+2] - (s_pfPrtVel[pIndx+2] + (m_v3Gravity[2] * dt))*dt*res;
	//		m_pCurPos[cIndx+0] = s_pfPrtPos[pIndx+0];
	//		m_pCurPos[cIndx+1] = s_pfPrtPos[pIndx+1];
	//	}

	//	clamp(m_pCurPos, cIndx);
	//}

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

	if(m_iNumVertices <= 1) return;

	Vec3 cm(0.0), cm_org(0.0);	// �d�S
	double mass = 0.0;	// ������
	
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

	//��������
	Vec3 mainPrtPos_cur = GetVertexPos(0);
	Vec3 mainPrtPos_org = GetOrgPos(0);

	//�d�S�Ƃ̋�������@�d�S�ʒu�͌Œ�
	Vec3 nowDx_center = m_vec3NowCm - mainPrtPos_cur;
	float nowDist_center = norm(nowDx_center);

	Vec3 orgDx_center = m_vec3OrgCm - mainPrtPos_org;
	float orgDist_center = norm(orgDx_center);

	Vec3 dirVec_center = nowDist_center > 0.0f ? 1.0f / nowDist_center * nowDx_center : Vec3(0.0f);

	Vec3 fixPos = GetVertexPos(0) + dirVec_center * (orgDist_center - nowDist_center);
	SetCurrentPos(0, fixPos);

	mainPrtPos_cur = GetVertexPos(0);
	m_fDefAmount = 0.0f;

	//��������@�N���X�^�ɑΉ����闱�q�ƁC���̑��̗��q
	for(int i = 1; i < m_iIndxNum; i++){
		if(CheckHole(i)){	continue;	}
		if(m_pFix[i]){		continue;	}

		Vec3 nowDx = GetVertexPos(i) - mainPrtPos_cur;
		float nowDist = norm(nowDx);

		Vec3 orgDx = GetOrgPos(i) - mainPrtPos_org;
		float orgDist = norm(orgDx);

		Vec3 dirVec = nowDist > 0.0f ? 1.0f / nowDist * nowDx : Vec3(0.0f);
		Vec3 dx = 1.0f/2.0f * dirVec * (orgDist - nowDist);

		SetCurrentPos(0, mainPrtPos_cur - dx);
		SetCurrentPos(i, GetVertexPos(i) + dx);

		mainPrtPos_cur = GetVertexPos(0);

		m_fDefAmount += norm2(dx);
	}

	//���E����
	for(int i = 0; i < m_iIndxNum; ++i)
	{
		clamp(m_pCurPos, i*SM_DIM);
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
			m_pCurPos[cIndx+0] = prdPos[0] - (s_pfClstrVel[pIndx * SM_DIM+0] + (m_v3Gravity[0] * dt))*dt*res;
			m_pCurPos[cIndx+1] = prdPos[1];
			m_pCurPos[cIndx+2] = prdPos[2];
		}

		if(m_pCurPos[cIndx+1] < m_v3Min[1] || m_pCurPos[cIndx+1] > m_v3Max[1]){
			m_pCurPos[cIndx+1] = prdPos[1] - (s_pfClstrVel[pIndx * SM_DIM+1] + (m_v3Gravity[1] * dt))*dt*res;
			m_pCurPos[cIndx+0] = prdPos[0] ;
			m_pCurPos[cIndx+2] = prdPos[2];
		}

		if(m_pCurPos[cIndx+2] < m_v3Min[2] || m_pCurPos[cIndx+2] > m_v3Max[2]){
			m_pCurPos[cIndx+2] = prdPos[2] - (s_pfClstrVel[pIndx * SM_DIM+2] + (m_v3Gravity[2] * dt))*dt*res;
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

	//��Ԃ̂��߂ɉ�]�ȊO�̕������܂ލs���ۑ�
	m_vOrientedPrtes[0]->Symmetric(S);
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

//�N���X�^���ɂ��闱�q�Ԃ̓��ς̑��a�@�N�H�[�^�j�I���ł̓���
//�S�ē����p���Ȃ�1�C�S�Ĕ��΂Ȃ�-1
float OrientedPrtObj::QuatInnerDotSum() const
{
	int prtNum = 1;
	float dot = 0.0f;
	mk_Quaternion baseQuat = m_vOrientedPrtes[0]->CurOrientation();

	for(unsigned i = 1; i < m_iIndxNum; i++){
		if( CheckHole(i) ){	continue;	}
		if(m_pFix[i]){		continue;	}

		mk_Quaternion quat = m_vOrientedPrtes[i]->CurOrientation();

		dot += baseQuat.x * quat.x + baseQuat.y * quat.y + baseQuat.z * quat.z + baseQuat.w * quat.w;

		prtNum++;
	}

	return dot/(float)prtNum;
}

float OrientedPrtObj::DirectionalVecInnerDotSum() const
{
	int prtNum = 1;
	float dot = 0.0f;
	Vec3 baseAxis = m_vOrientedPrtes[0]->Axis();

	for(unsigned i = 1; i < m_iIndxNum; i++){
		if( CheckHole(i) ){	continue;	}
		if(m_pFix[i]){		continue;	}

		Vec3 axis = m_vOrientedPrtes[i]->Axis();

		dot += baseAxis[0] * axis[0] + baseAxis[1] * axis[1] + baseAxis[2] * axis[2];

		prtNum++;
	}

	return dot/(float)prtNum;
}

float OrientedPrtObj::FinalDefAmount() const
{
	float defAmount = 0.0f;

	//�ŏI�ʒu�̏d�S
	Vec3 finalCm(0.0f);
	for(unsigned i = 0; i < m_iIndxNum; i++){
		if( CheckHole(i) ){	continue;	}
		if(m_pFix[i]){		continue;	}

		//int pIndx = GetParticleIndx(i) * SM_DIM;
		//finalCm += Vec3(s_pfClstrPos[pIndx + X], s_pfClstrPos[pIndx + Y], s_pfClstrPos[pIndx + Z]);

		int pIndx = GetParticleIndx(i) * 4;
		finalCm += Vec3(s_pfPrtPos[pIndx + X], s_pfPrtPos[pIndx + Y], s_pfPrtPos[pIndx + Z]);
	}

	finalCm /= GetNumVertices();

	for(unsigned i = 1; i < m_iIndxNum; i++){
		if( CheckHole(i) ){	continue;	}
		if(m_pFix[i]){		continue;	}

		//�����ʒu���Ǐ����W�n�ɒ����āC��]�s����|����
		rxMatrix3 R = m_vOrientedPrtes[i]->Rotation();
		Vec3 orgPos = R * (GetOrgPos(i) - m_vec3OrgCm)/* + finalCm*/;

		//���݈ʒu���Ǐ����W�n�ɕϊ�
		//int pIndx = GetParticleIndx(i) * SM_DIM;
		//Vec3 finalPos = Vec3(s_pfClstrPos[pIndx + X], s_pfClstrPos[pIndx + Y], s_pfClstrPos[pIndx + Z])/* - finalCm*/;

		int pIndx = GetParticleIndx(i) * 4;
		Vec3 finalPos = Vec3(s_pfPrtPos[pIndx + X], s_pfPrtPos[pIndx + Y], s_pfPrtPos[pIndx + Z]) - finalCm;

		//�Ǐ����W�n�ɂ�����ړ��ʂ�ό`�ʂƂ���
		defAmount += norm2(orgPos-finalPos);
	}

	return defAmount / GetNumVertices();
}

float OrientedPrtObj::VolumeChange() const
{
	float defAmount = 0.0f;

	//�ŏI�ʒu�̏d�S
	Vec3 finalCm(0.0f);
	for(unsigned i = 0; i < m_iIndxNum; i++){
		if( CheckHole(i) ){	continue;	}
		if(m_pFix[i]){		continue;	}

		int pIndx = GetParticleIndx(i) * 4;
		finalCm += Vec3(s_pfPrtPos[pIndx + X], s_pfPrtPos[pIndx + Y], s_pfPrtPos[pIndx + Z]);
	}

	finalCm /= GetNumVertices();

	rxMatrix3 Apq(0.0);
	Vec3 p(0.0f);
	Vec3 q(0.0f);

	for(int i = 0; i < m_iIndxNum; ++i)
	{
		if( CheckHole(i) ){	continue;	}

		int cIndx = i*SM_DIM;
		int pIndx = GetParticleIndx(i)*4;

		for(int j = 0; j < SM_DIM; j++)
		{
			p[j] = s_pfPrtPos[pIndx+j]-finalCm[j];
			q[j] = m_pOrgPos[cIndx+j]-m_vec3OrgCm[j];
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

	defAmount = Apq.Determinant();

	return defAmount;
}

bool OrientedPrtObj::IsThereOutlier() const
{
	float defAmount = 0.0f;
	int outlierNum = 0;

	//�ŏI�ʒu�̏d�S
	Vec3 finalCm(0.0f);
	for(unsigned i = 0; i < m_iIndxNum; i++){
		if( CheckHole(i) ){	continue;	}
		if(m_pFix[i]){		continue;	}

		int pIndx = GetParticleIndx(i) * 4;
		finalCm += Vec3(s_pfPrtPos[pIndx + X], s_pfPrtPos[pIndx + Y], s_pfPrtPos[pIndx + Z]);
	}

	finalCm /= GetNumVertices();

	for(unsigned i = 1; i < m_iIndxNum; i++){
		if( CheckHole(i) ){	continue;	}
		if(m_pFix[i]){		continue;	}

		//�����ʒu���Ǐ����W�n�ɒ����āC��]�s����|����
		rxMatrix3 R = m_vOrientedPrtes[i]->Rotation();
		Vec3 orgPos = R * (GetOrgPos(i) - m_vec3OrgCm);

		//���݈ʒu���Ǐ����W�n�ɕϊ�
		int pIndx = GetParticleIndx(i) * 4;
		Vec3 finalPos = Vec3(s_pfPrtPos[pIndx + X], s_pfPrtPos[pIndx + Y], s_pfPrtPos[pIndx + Z]) - finalCm;

		//�傫���ړ��������q�����邩�ǂ���
		if(norm2(orgPos-finalPos) > 0.007f){
			outlierNum++;
		}
	}

	//return outlierNum > 0;
	return outlierNum > (GetNumVertices()/3);
}