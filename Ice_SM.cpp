/*!
  @file Ice_SM.h
	
  @brief ShapeMatching�@����ɂ������ω��V�~�����[�V����
  @ref M. Muller et al., "Meshless deformations based on shape matching", SIGGRAPH2005. 
 
  @author Ryo Nakasone
  @date 2013-10
*/

#include "Ice_SM.h"

const float* Ice_SM::s_pfPrtPos;
const float* Ice_SM::s_pfPrtVel;

//�f�o�C�X�|�C���^
float* Ice_SM::sd_PrtPos;	
cudaGraphicsResource* Ice_SM::sd_PrtPosVbo;
float* Ice_SM::sd_PrtVel;

float* Ice_SM::d_OrgPos;
float* Ice_SM::d_CurPos;
float* Ice_SM::d_NewPos;
float* Ice_SM::d_GoalPos;
float* Ice_SM::d_Mass;
float* Ice_SM::d_Vel;

bool* Ice_SM::d_Fix;

int* Ice_SM::d_PIndxes;

int* Ice_SM::d_IndxSet;					//�N���X�^�̃f�[�^�̊J�n�Y���ƏI���Y����ۑ�

int Ice_SM::s_vertSum;					//�S�N���X�^�Ɋ܂܂�闱�q�̑���



Ice_SM::Ice_SM(int obj) : rxShapeMatching(obj)
{
	m_mtrx3Apq = rxMatrix3(0.0);
	m_mtrxBeforeU.makeIdentity();

	m_iIndxNum = 0;

	m_fpAlphas = new float[MAXPARTICLE];
	m_fpBetas =	new float[MAXPARTICLE];

	m_ipLayeres = new int[MAXPARTICLE];
}

/*!
 * �f�X�g���N�^
 */
Ice_SM::~Ice_SM()
{
}

void Ice_SM::InitGPU(const vector<Ice_SM*>& ice_sm, float* d_pos, cudaGraphicsResource* d_pos_vbo, float* d_vel)
{
	//�f�o�C�X�|�C���^�̃A�h���X��ۑ�
	sd_PrtPos = d_pos;
	sd_PrtPosVbo = d_pos_vbo;
	sd_PrtVel = d_vel;

	//�f�o�C�X���̃��������m��
	//�i�ő�N���X�^���j�~�i�N���X�^���ۑ��ł���ő嗱�q���j�@�Ń��������m�ہD
	cudaMalloc((void**)&d_OrgPos,	sizeof(float) * MAXCLUSTER * MAXPARTICLE * SM_DIM);
	cudaMalloc((void**)&d_CurPos,	sizeof(float) * MAXCLUSTER * MAXPARTICLE * SM_DIM);
	cudaMalloc((void**)&d_Vel,		sizeof(float) * MAXCLUSTER * MAXPARTICLE * SM_DIM);

	cudaMalloc((void**)&d_Mass,		sizeof(float) * MAXCLUSTER * MAXPARTICLE);
	cudaMalloc((void**)&d_Fix,		sizeof(bool)  * MAXCLUSTER * MAXPARTICLE);
	cudaMalloc((void**)&d_PIndxes,	sizeof(int)   * MAXCLUSTER * MAXPARTICLE);

	cudaMalloc((void**)&d_IndxSet,	sizeof(int)   * MAXCLUSTER * 2);

	//CPU�̃f�[�^���P�����z��ɃR�s�[
	//�z�X�g���̃f�[�^���f�o�C�X���֓]�����ď�����
		//���������m�ۂł��Ȃ��̂ŁA�������ăR�s�[����
	float* oPoses = new float[MAXPARTICLE * SM_DIM];
	float* cPoses = new float[MAXPARTICLE * SM_DIM];
	float* veles  = new float[MAXPARTICLE * SM_DIM];

	float* mass = new float[MAXPARTICLE];
	bool* fix = new bool[MAXPARTICLE];
	int* indx = new int[MAXPARTICLE];

	int* set = new int[2];

	s_vertSum = 0;

	for(int i = 0; i < MAXCLUSTER; i++)
	{
		Ice_SM* sm = ice_sm[i];

		for(int j = 0; j < MAXPARTICLE; j++)
		{
			//���݈ʒu�C�����ʒu�C���x
			Vec3 oPos = sm->GetVertexPos(j);
			Vec3 cPos = sm->GetOrgPos(j);
			Vec3 vel  = sm->GetVertexVel(j);

			for(int k = 0; k < SM_DIM; k++)
			{
				oPoses[j*SM_DIM + k] = oPos[k];
				cPoses[j*SM_DIM + k] = cPos[k];
				veles [j*SM_DIM + k] = vel[k];
			}

			//���ʁC�I���t���O�C���q�ԍ�
			mass[j] = sm->GetMass(j);
			fix[j] = false;
			indx[j] = sm->GetParticleIndx(j);
		}

		//�z��̂ǂ�����ǂ��܂ł����܂��Ă��邩
		set[0] = i*MAXPARTICLE;
		set[1] = i*MAXPARTICLE + sm->GetNumVertices()-1;

		//�����q��
		s_vertSum += sm->GetNumVertices();

		//������
		int vecSize = i * MAXPARTICLE * SM_DIM;

		cudaMemcpy(d_OrgPos+vecSize,	oPoses,	sizeof(float) * MAXPARTICLE * SM_DIM, cudaMemcpyHostToDevice);
		cudaMemcpy(d_CurPos+vecSize,	cPoses,	sizeof(float) * MAXPARTICLE * SM_DIM, cudaMemcpyHostToDevice);
		cudaMemcpy(d_Vel   +vecSize,	veles,	sizeof(float) * MAXPARTICLE * SM_DIM, cudaMemcpyHostToDevice);

		int smSize = i * MAXPARTICLE;

		cudaMemcpy(d_Mass   +smSize,	mass,	sizeof(float) * MAXPARTICLE, cudaMemcpyHostToDevice);
		cudaMemcpy(d_Fix    +smSize,	fix,	sizeof(bool)  * MAXPARTICLE, cudaMemcpyHostToDevice);
		cudaMemcpy(d_PIndxes+smSize,	indx,	sizeof(int)   * MAXPARTICLE, cudaMemcpyHostToDevice);

		cudaMemcpy(d_IndxSet+i*2,		set,	sizeof(int) * 2, cudaMemcpyHostToDevice);
	}

	////���Z�b�g
	//for(int i = 0; i < MAXCLUSTER; i++)
	//{
	//	Ice_SM* sm = ice_sm[i];

	//	for(int j = 0; j < MAXPARTICLE; j++)
	//	{
	//		for(int k = 0; k < SM_DIM; k++)
	//		{
	//			oPoses[j*SM_DIM + k] = 0.0f;
	//			cPoses[j*SM_DIM + k] = 0.0f;
	//			veles [j*SM_DIM + k] = 0.0f;
	//		}

	//		mass[j] = 0.0f;
	//		fix [j] = false;
	//		indx[j] = 0;
	//	}
	//}

	////�z�X�g���̃f�[�^��]���������ʂ��_���v
	//ofstream ofs( "DtoH_Test.txt" );
	//ofs << "DtoH_Test" << endl;

	////�f�o�C�X���̃f�[�^���z�X�g�֕������ē]��
	//for(int i = 0; i < MAXCLUSTER; i++)
	//{
	//	ofs << "�N���X�^ " << i << endl;

	//	int vecSize = i * MAXPARTICLE * SM_DIM;
	//	Ice_SM* sm = ice_sm[i];

	//	cudaMemcpy(oPoses,	d_OrgPos+vecSize,	sizeof(float) * MAXPARTICLE * SM_DIM, cudaMemcpyDeviceToHost);
	//	cudaMemcpy(cPoses,	d_CurPos+vecSize,	sizeof(float) * MAXPARTICLE * SM_DIM, cudaMemcpyDeviceToHost);
	//	cudaMemcpy(veles,	d_Vel   +vecSize,	sizeof(float) * MAXPARTICLE * SM_DIM, cudaMemcpyDeviceToHost);

	//	int smSize = i * MAXPARTICLE;

	//	cudaMemcpy(mass,	d_Mass   +smSize,	sizeof(float) * MAXPARTICLE, cudaMemcpyDeviceToHost);
	//	cudaMemcpy(fix,		d_Fix    +smSize,	sizeof(bool)  * MAXPARTICLE, cudaMemcpyDeviceToHost);
	//	cudaMemcpy(indx,	d_PIndxes+smSize,	sizeof(int)   * MAXPARTICLE, cudaMemcpyDeviceToHost);

	//	cudaMemcpy(set,		d_IndxSet+i*2,		sizeof(int)   * 2, cudaMemcpyDeviceToHost);

	//	//�]���ł��Ă��邩�̊m�F
	//	for(int j = 0; j < MAXPARTICLE; j++)
	//	{
	//		//ofs << j << " :: "; 
	//		
	//		//�����ʒu
	//		//ofs << "oPoses = ( "; 
	//		//
	//		//for(int k = 0; k < SM_DIM; k++)
	//		//{
	//		//	ofs << oPoses[j*SM_DIM + k] << " ";
	//		//}
	//		//ofs << "), ";

	//		////mass[j] = j;
	//		////fix[j] = false;
	//		//if(indx[j] == -1){	ofs << "end" << endl; break;		}
	//		//else{				ofs << "indx = " << indx[j] << ",";	}

	//		//ofs << endl;
	//	}

	//	//�J�n�ʒu�ƏI���ʒu
	//	ofs << "set = ( " << set[0] << " ," << set[1] << ")" << endl;

	//	//�����P�ł���΂悢
	//	//cout << "start = " << set[0] << ", end = " << set[1] << ", Size = " << set[1]-set[0] << ", sm = " << sm->GetNumVertices();
	//	//cout << endl;
	//}

	//�]������̃f�[�^���ׂĂ݂�
	//for(int i = 0; i < MAXCLUSTER;i++)
	//{
	//	ice_sm[i]->CopyDeviceToInstance(i);
	//}

	delete[] oPoses;
	delete[] cPoses;
	delete[] veles;
	delete[] mass;
	delete[] fix;
	delete[] indx;
	delete[] set;
}

void Ice_SM::InitGPU_Instance()
{


}

void Ice_SM::AddVertex(const Vec3 &pos, double mass, int pIndx)
{
	rxShapeMatching::AddVertex(pos, mass, pIndx);

	//�ő�Y���ԍ��̍X�V
	if(m_iNumVertices > m_iIndxNum)
	{
		m_iIndxNum = m_iNumVertices;
	}

	//m_iLinearDeformation.push_back(0);
	//m_iVolumeConservation.push_back(0);

	//�d�S�̍X�V
	double massSum = 0.0;	// ������
	m_vec3OrgCm = Vec3(0.0);

	// �d�S���W�̌v�Z
	for(int i = 0; i < m_iIndxNum;++i)
	{
		if( CheckHole(i) ){	continue;	}

		double m = m_pMass[i];
		//if(m_pFix[i]) m *= 300.0;	// �Œ�_�̎��ʂ�傫������
		massSum += m;
		m_vec3OrgCm += Vec3(m_pOrgPos[i*SM_DIM+0], m_pOrgPos[i*SM_DIM+1], m_pOrgPos[i*SM_DIM+2]) * m;
	}
	m_vec3OrgCm /= massSum;

	//�ό`�s��Aqq�̍X�V
	Vec3 p, q;
	for(int i = 0; i < m_iIndxNum; ++i)
	{
		if( CheckHole(i) ){	continue;	}

		q = Vec3(m_pOrgPos[i*SM_DIM+0], m_pOrgPos[i*SM_DIM+1], m_pOrgPos[i*SM_DIM+2])-m_vec3OrgCm;
		double m = m_pMass[i];

		m_mtrx3AqqInv(0,0) += m*q[0]*q[0];
		m_mtrx3AqqInv(0,1) += m*q[0]*q[1];
		m_mtrx3AqqInv(0,2) += m*q[0]*q[2];
		m_mtrx3AqqInv(1,0) += m*q[1]*q[0];
		m_mtrx3AqqInv(1,1) += m*q[1]*q[1];
		m_mtrx3AqqInv(1,2) += m*q[1]*q[2];
		m_mtrx3AqqInv(2,0) += m*q[2]*q[0];
		m_mtrx3AqqInv(2,1) += m*q[2]*q[1];
		m_mtrx3AqqInv(2,2) += m*q[2]*q[2];
	}

	m_mtrx3AqqInv = m_mtrx3AqqInv.Inverse();

	////Q�̒ǉ�
	//m_vvec3OrgQ.resize(m_iNumVertices);
	//for(int i = 0; i < m_iNumVertices;++i)
	//{
	//	m_vvec3OrgQ[i] = m_vOrgPos[i]-m_vec3OrgCm;
	//}
}

void Ice_SM::Remove(int indx)
{
	if(m_iNumVertices <= 0) return;
	m_iNumVertices--;

	m_iPIndxes[indx] = -1;

	for(int i = 0; i < SM_DIM; i++)
	{
		m_pOrgPos[indx*SM_DIM+i] = 0.0;
		m_pCurPos[indx*SM_DIM+i] = 0.0;
		//m_pNewPos[indx*SM_DIM+i] = 0.0;
		//m_pGoalPos[indx*SM_DIM+i] = 0.0;
		m_pVel[indx*SM_DIM+i] = 0.0;
	}

	m_pMass[indx] = 0.0;
	m_pFix[indx] = false;

	m_ipLayeres[indx] = -1;

	m_fpAlphas[indx] = -1;
	m_fpBetas[indx] = -1;

	//m_iLinearDeformation.erase	( m_iLinearDeformation	.begin() + indx);
	//m_iVolumeConservation.erase	( m_iVolumeConservation	.begin() + indx);
}

void Ice_SM::Clear()
{
	rxShapeMatching::Clear();

	m_iIndxNum = 0;

	for(int i = 0; i < MAXPARTICLE; i++)
	{
		m_fpAlphas[i] = 0.0f;
		m_fpBetas[i] = 0.0f;

		m_ipLayeres[i] = -1;
	}

	//m_iLinearDeformation	.clear();
	//m_iVolumeConservation	.clear();
}

bool Ice_SM::CheckIndx(int pIndx)
{//	cout << __FUNCTION__ << endl;

	//vector<int>::iterator check = find( m_iPIndxes.begin(), m_iPIndxes.end(), pIndx);
	//if( check == m_iPIndxes.end() )	return false;
	//return true;

	for(int i = 0; i < m_iIndxNum; i++)
	{
		if(m_iPIndxes[i] == pIndx)
		{
			return true;
		}
	}

	return false;
}

int	Ice_SM::SearchIndx(int pIndx)
{
	//vector<int>::iterator check = find( m_iPIndxes.begin(), m_iPIndxes.end(), pIndx);
	//if( check == m_iPIndxes.end() )	return -1;
	//return m_iPIndxes.size() - (m_iPIndxes.end() - check);

	for(int i = 0; i < m_iIndxNum; i++)
	{
		if(m_iPIndxes[i] == pIndx)
		{
			return i;
		}
	}

	return -1;
}

/*
 *	�Y���������ǂ����̃`�F�b�N
 */
bool Ice_SM::CheckHole(int oIndx)
{
	return (m_iPIndxes[oIndx] == -1);
}

/*!
 * Shape Matching�@
 *  - �ڕW�ʒu���v�Z���āCm_vNewPos�����̈ʒu�ɋ߂Â���
 *  - �e���q�Ƀ��ƃ������������o�[�W����
 * @param[in] dt �^�C���X�e�b�v��
 */
void Ice_SM::ShapeMatching(float* newPos, double dt)
{
	if(m_iNumVertices <= 1) return;

	QueryCounter qc;

	double mass = 0.0;	// ������

	// �d�S���W�̌v�Z
	for(int i = 0; i < m_iIndxNum; ++i){
		//if( CheckHole(i) ){	continue;	}
		//if(m_pFix[i]){	mass += m_pMass[i]*300;	}	// �Œ�_�̎��ʂ�傫������
		mass += m_pMass[i];
	}

	Vec3 cm(1 / mass * m_vec3NowCm), cm_org(m_vec3OrgCm);	// �d�S
	Vec3 p(0.0), q(0.0);

	//Apq�̍s�񎮂����߁C���]���邩�𔻒�
	//�s����ȏꍇ�������̂Ł~
	//if( Apq.Determinant() < 0.0 && m_iNumVertices >= 4)
	//{
		//cout << "before det < 0" << endl;
		//�P�@�����𔽓]
		//Apq(0,2) = -Apq(0,2);
		//Apq(1,2) = -Apq(1,2);
		//Apq(2,2) = -Apq(2,2);

		//�Q�@a2��a3������
		//double tmp;
		//tmp = Apq(0,2);
		//Apq(0,2) = Apq(0,1);
		//Apq(0,1) = tmp;

		//tmp = Apq(1,2);
		//Apq(1,2) = Apq(1,1);
		//Apq(1,1) = tmp;

		//tmp = Apq(2,2);
		//Apq(2,2) = Apq(2,1);
		//Apq(2,1) = tmp;
	//}
	
	//cout << "�v���J�n1" << endl;
	//qc.Start();

	rxMatrix3 R, S;
	//PolarDecomposition(m_mtrx3Apq, R, S, m_mtrxBeforeU);
	PolarDecomposition(m_mtrx3Apq, R, S);

	//double end1 = qc.End()/*/100*/;

	if(m_bLinearDeformation)
	{
		// Linear Deformations
		rxMatrix3 A(m_mtrx3Apq * m_mtrx3AqqInv);	// A = Apq*Aqq^-1
		//A = Apq*Aqq.Inverse();	//A = Apq*Aqq^-1

		//�̐ϕۑ��̂��߂Ɂ�(det(A))�Ŋ���
		if(m_bVolumeConservation){
			double det = fabs(A.Determinant());
			if(det > RX_FEQ_EPS){
				det = 1.0/sqrt(det);
				if(det > 2.0) det = 2.0;
				A *= det;
			}
		}

		//cout << "�v���J�n2" << endl;
		//qc.Start();

		// �ڕW���W���v�Z���C���݂̒��_���W���ړ�
		for(int i = 0; i <  m_iIndxNum; ++i)
		{
			if( CheckHole(i) ){	continue;	}
			if(m_pFix[i]) continue;
			// ��]�s��R�̑���̍s��RL=��A+(1-��)R���v�Z
			//rxMatrix3 RL = m_dBetas[i]*A+(1.0-m_dBetas[i])*R;
			//q = m_vOrgPos[i]-cm_org;
			//m_vGoalPos[i] = RL*q+cm;
			//m_vNewPos[i] += (m_vGoalPos[i]-m_vNewPos[i])*m_dAlphas[i];

			//Vec3 np(m_pNewPos[i*SM_DIM+0], m_pNewPos[i*SM_DIM+1], m_pNewPos[i*SM_DIM+2]);
			int cIndx = i*SM_DIM;
			Vec3 op(m_pOrgPos[cIndx+0], m_pOrgPos[cIndx+1], m_pOrgPos[cIndx+2]);

			q = op-cm_org;
			Vec3 gp(R*q+cm);
			//m_vGoalPos[i] = R*m_vvec3OrgQ[i]+cm;

			for(int j = 0; j < SM_DIM; j++)
			{
				int jcIndx = cIndx+j;
				//m_pGoalPos[i*SM_DIM+j] = gp[j];
				//m_pNewPos[i*SM_DIM+j] += (gp[j]-np[j])*m_fpAlphas[i];
				m_pCurPos[jcIndx] += (gp[j]-m_pCurPos[jcIndx])*m_fpAlphas[i];
			}
		}
		//double end2 = qc.End()/*/100*/;
		//cout << "�v���I��1::" << end1 << endl;
		//cout << "�v���I��2::" << end2 << endl;
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
		//	p = m_vNewPos[i]-cm;
		//	q = m_vOrgPos[i]-cm_org;

		//	// q~�̌v�Z
		//	double qt[9];
		//	qt[0] = q[0];      qt[1] = q[1];      qt[2] = q[2];
		//	qt[3] = q[0]*q[0]; qt[4] = q[1]*q[1]; qt[5] = q[2]*q[2];
		//	qt[6] = q[0]*q[1]; qt[7] = q[1]*q[2]; qt[8] = q[2]*q[0];

		//	// A~pq = ��mpq~ �̌v�Z
		//	double m = m_pMass[i];
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
		//		At[i][j] *= m_dBetas[i];
		//		if(j < 3){
		//			At[i][j] += (1.0f-m_dBetas[i])*R(i,j);
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
		//	q = m_vOrgPos[i]-cm_org;

		//	for(int k = 0; k < 3; ++k){
		//		m_vGoalPos[i][k] = At[k][0]*q[0]+At[k][1]*q[1]+At[k][2]*q[2]
		//						  +At[k][3]*q[0]*q[0]+At[k][4]*q[1]*q[1]+At[k][5]*q[2]*q[2]+
		//						  +At[k][6]*q[0]*q[1]+At[k][7]*q[1]*q[2]+At[k][8]*q[2]*q[0];
		//	}

		//	m_vGoalPos[i] += cm;
		//	m_vNewPos[i] += (m_vGoalPos[i]-m_vNewPos[i])*m_dAlphas[i];
		//}
	}

}

/*!
 * �O��
 *  - �d�͂Ƌ��E�ǂ���̗͂̉e��
 * @param[in] dt �^�C���X�e�b�v��
 */
void Ice_SM::calExternalForces(float* newPos, double dt)
{
	// �d�͂̉e����t���C���x�𔽉f
	for(int i = 0; i < m_iIndxNum; ++i)
	{
		if( CheckHole(i) ){	continue;	}
		int pIndx = m_iPIndxes[i]*4;
		int cIndx = i*SM_DIM;

		for(int j = 0; j < SM_DIM; j++)
		{
			int jpIndx = pIndx+j;
			int jcIndx = cIndx+j;

			//m_pNewPos[jcIndx] = s_pfPrtPos[jpIndx]+s_pfPrtVel[jpIndx]*dt;
			//m_pGoalPos[jcIndx] = m_pOrgPos[jcIndx];
			
			//newPos[jcIndx] = s_pfPrtPos[jpIndx]+s_pfPrtVel[jpIndx]*dt;

			m_pCurPos[jcIndx] = s_pfPrtPos[jpIndx]+s_pfPrtVel[jpIndx]*dt;
		}
	}

	// ���E�ǂ̉e��
	//���������Ȃ�d���Ȃ邪�C����͂���݂���
	double res = 0.9;	// �����W��
	for(int i = 0; i < m_iIndxNum; ++i)
	{
		//if( CheckHole(i) ){	continue;	}
		////if(m_pFix[i]) continue;
		//Vec3 &p = m_pCurPos[i];

		//Vec3 &v = m_pVel[i];
		//if(np[0] < m_v3Min[0] || np[0] > m_v3Max[0]){
		//	np[0] = p[0]-v[0]*dt*res;
		//	np[1] = p[1];
		//	np[2] = p[2];
		//}
		//if(np[1] < m_v3Min[1] || np[1] > m_v3Max[1]){
		//	np[1] = p[1]-v[1]*dt*res;
		//	np[0] = p[0] ;
		//	np[2] = p[2];
		//}
		//if(np[2] < m_v3Min[2] || np[2] > m_v3Max[2]){
		//	np[2] = p[2]-v[2]*dt*res;
		//	np[0] = p[0];
		//	np[1] = p[1];
		//}

		//clamp(newPos, i*SM_DIM);
		clamp(m_pCurPos, i*SM_DIM);
	}
}


//----------------------------------------------�\���b�h��---------------------------------------------
/*!
 * Shape Matching�@
 *  - �ڕW�ʒu���v�Z���āCm_vNewPos�����̈ʒu�ɋ߂Â���
 *  - �e���q�Ƀ��ƃ������������o�[�W����
 * @param[in] dt �^�C���X�e�b�v��
 */
void Ice_SM::ShapeMatchingSolid(float* newPos, double dt)
{
	if(m_iNumVertices <= 1) return;

	QueryCounter qc;

	Vec3 cm(0.0), cm_org(0.0);	// �d�S
	double mass = 0.0;	// ������

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
			//cm[j] += m_pNewPos[cIndx+j]*m;
			//cm[j] += newPos[cIndx+j]*m;
			cm[j] += m_pCurPos[cIndx+j]*m;
		}
	}

	cm /= mass;
	cm_org = m_vec3OrgCm;

	//if(m_iObjectNo < 1)
	//{
	//	cout << "cpu, number::" << m_iObjectNo << " cm = " << cm << ", cm_org = " << cm_org << endl;
	//}

	rxMatrix3 Apq(0.0), Aqq(0.0);
	Vec3 p, q;

	// Apq = ��mpq^T
	// Aqq = ��mqq^T
	for(int i = 0; i < m_iIndxNum; ++i)
	{
		if( CheckHole(i) ){	continue;	}

		int cIndx = i*SM_DIM;

		for(int j = 0; j < SM_DIM; j++)
		{
			//p[j] = m_pNewPos[cIndx+j]-cm[j];

			//p[j] = newPos[cIndx+j]-cm[j];
			p[j] = m_pCurPos[cIndx+j]-cm[j];
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

	//if(m_iObjectNo < 1)
	//{
	//	cout << "cpu, number::" << m_iObjectNo << " Apq = " << Apq << ", Aqq = " << Apq << endl;
	//}

	//Apq�̍s�񎮂����߁C���]���邩�𔻒�
	//�s����ȏꍇ�������̂Ł~
	//if( Apq.Determinant() < 0.0 && m_iNumVertices >= 4)
	//{
		//cout << "before det < 0" << endl;
		//�P�@�����𔽓]
		//Apq(0,2) = -Apq(0,2);
		//Apq(1,2) = -Apq(1,2);
		//Apq(2,2) = -Apq(2,2);

		//�Q�@a2��a3������
		//double tmp;
		//tmp = Apq(0,2);
		//Apq(0,2) = Apq(0,1);
		//Apq(0,1) = tmp;

		//tmp = Apq(1,2);
		//Apq(1,2) = Apq(1,1);
		//Apq(1,1) = tmp;

		//tmp = Apq(2,2);
		//Apq(2,2) = Apq(2,1);
		//Apq(2,1) = tmp;
	//}
	
	//qc.Start();

	rxMatrix3 R, S;
	//PolarDecomposition(Apq, R, S, m_mtrxBeforeU);
	PolarDecomposition(Apq, R, S);

	//if(m_iObjectNo < 1)
	//{
	//	cout << "cpu, number::" << m_iObjectNo << " R = " << R << ", S = " << S << endl;
	//}

	//double end1 = qc.End()/*/100*/;

	if(m_bLinearDeformation)
	{
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

		//cout << "�v���J�n2" << endl;
		//qc.Start();

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
				int jcIndx = cIndx+j;

				//m_pGoalPos[jcIndx] = gp[j];
				//m_pNewPos[jcIndx] += (gp[j]-m_pNewPos[jcIndx])*m_dAlphas[i];
				//newPos[jcIndx] += (gp[j]-newPos[jcIndx])*m_dAlphas[i];
				m_pCurPos[jcIndx] += (gp[j]-m_pCurPos[jcIndx])*m_fpAlphas[i];
			}
		}

		//double end2 = qc.End()/*/100*/;
		//cout << "�v���I��1::" << end1 << endl;
		//cout << "�v���I��2::" << end2 << endl;
	}
}

//----------------------------------------------�\���b�h��---------------------------------------------

/*!
 * ���x�ƈʒu�̍X�V
 *  - �V�����ʒu�ƌ��݂̈ʒu���W���瑬�x���Z�o
 * @param[in] dt �^�C���X�e�b�v��
 */
void Ice_SM::integrate(float* newPos, double dt)
{
	double dt1 = 1.0/dt;
	for(int i = 0; i < m_iIndxNum; ++i)
	{
		if( CheckHole(i) ){	continue;	}
		int pIndx = m_iPIndxes[i]*4;

		for(int j = 0; j < SM_DIM; j++)
		{
			int cIndx = i*SM_DIM+j;
			//m_pVel[cIndx] = (m_pNewPos[cIndx] - s_pfPrtPos[pIndx+j]) * dt1;/*+ m_v3Gravity * dt * 1.0*/;
			//m_pCurPos[cIndx] = m_pNewPos[cIndx];
			//m_pVel[cIndx] = (newPos[cIndx] - s_pfPrtPos[pIndx+j]) * dt1;/*+ m_v3Gravity * dt * 1.0*/;
			//m_pCurPos[cIndx] = newPos[cIndx];
			m_pVel[cIndx] = (m_pCurPos[cIndx] - s_pfPrtPos[pIndx+j]) * dt1 + m_v3Gravity[j] * dt * 0.1f;
		}
	}
}

void Ice_SM::CalcDisplaceMentVectorCm()
{
	Vec3 preVec(0.0);
	Vec3 nowVec(0.0);
	double mass = 0.0;

	for(int i = 0; i < m_iIndxNum; i++)
	{
		if( CheckHole(i) ){	continue;	}

		int pIndx = m_iPIndxes[i] * 4;
		int cIndx = i*SM_DIM;

		preVec += Vec3(s_pfPrtPos[pIndx+0], s_pfPrtPos[pIndx+1], s_pfPrtPos[pIndx+2]);
		nowVec += Vec3(m_pGoalPos[cIndx+0], m_pGoalPos[cIndx+1], m_pGoalPos[cIndx+2]);
		mass += m_pMass[i];
	}

	preVec /= mass;
	nowVec /= mass;

	m_vec3DisCm = nowVec-preVec;
}

/*!
 * �V�~�����[�V�����X�e�b�v��i�߂�
 */
void Ice_SM::Update()
{//	cout << "Ice_update" << endl;
	//calCollision(m_dDt);
	QueryCounter qc;

	float* newPos = new float[m_iNumVertices*SM_DIM];

	qc.Start();
	calExternalForces(newPos, m_dDt);
	double end1 = qc.End()/*/100*/;

	qc.Start();
	//ShapeMatching(newPos, m_dDt);		//�p�X��p����������
	ShapeMatchingSolid(newPos, m_dDt);	//���ʂ̌v�Z
	double end2 = qc.End()/*/100*/;

	qc.Start();
	integrate(newPos, m_dDt);
	double end3 = qc.End()/*/100*/;

	//cout << "�v���I��1::" << end1 << endl;
	//cout << "�v���I��2::" << end2 << endl;
	//cout << "�v���I��3::" << end3 << endl;

	//CalcDisplaceMentVectorCm();	//�d�S�ړ���

	delete newPos;
}

/*
 *	GPU��p���ăV�~�����[�V������i�߂� OpenMP�͎g���Ȃ��̂ɒ���
 */
void Ice_SM::UpdateGPU()
{
	//�z�X�g���̃f�[�^���f�o�C�X���֓]��
	//cudaMemcpy(d_CurPos, m_pCurPos, sizeof(float) * m_iNumVertices * SM_DIM, cudaMemcpyHostToDevice);
	//cudaMemcpy(d_Vel,		m_pVel,	sizeof(float) * m_iNumVertices * SM_DIM, cudaMemcpyHostToDevice);

	//GPU����
	LaunchShapeMatchingGPU(sd_PrtPos, sd_PrtPosVbo, sd_PrtVel, d_OrgPos, d_CurPos, d_Vel, d_PIndxes, d_IndxSet, 0.02, s_vertSum);
}

//GPU�̌v�Z���ʂ��e�C���X�^���X�փR�s�[����
void Ice_SM::CopyDeviceToInstance(int num)
{	//cout << __FUNCTION__ << " num = " << num << endl;

	//�f�o�C�X���̃f�[�^���z�X�g���֓]��
	float* cPoses = new float[MAXPARTICLE * SM_DIM];
	float* veles  = new float[MAXPARTICLE * SM_DIM];

	int vecSize = num * MAXPARTICLE * SM_DIM;

	cudaMemcpy(cPoses,	d_CurPos+vecSize,	sizeof(float) * MAXPARTICLE * SM_DIM, cudaMemcpyDeviceToHost);
	cudaMemcpy(veles,	d_Vel   +vecSize,	sizeof(float) * MAXPARTICLE * SM_DIM, cudaMemcpyDeviceToHost);

	for(int j = 0; j < MAXPARTICLE; j++)
	{
		if(GetParticleIndx(j) == -1){	continue;	}

		//�ʒu�C���x
		Vec3 cPos = Vec3(cPoses[j*SM_DIM + 0], cPoses[j*SM_DIM + 1], cPoses[j*SM_DIM + 2]);
		Vec3 vel  = Vec3(veles [j*SM_DIM + 0], veles [j*SM_DIM + 1], veles [j*SM_DIM + 2]);

		//if(num < 2)
		//{
		//	////CPU�ł̉^���v�Z���ʂƔ�r
		//	Vec3 cpuPos = GetVertexPos(j);
		//	Vec3 cpuVel = GetVertexVel(j);

		//	cout << "num = " << num << endl;
		//	cout << "posDiv = " << norm2(cpuPos)-norm2(cPos) << endl;
		//	cout << "velDiv = " << norm2(cpuVel)-norm2(vel) << endl;
		//}

		SetCurrentPos(j, cPos);
		SetVelocity(j, vel);
	}

	delete[] cPoses;
	delete[] veles;
}

//----------------------------------------�f�o�b�O--------------------------------------
void Ice_SM::DebugIndx()
{	cout << "DebugIndx Indx = ";
	for(int i = 0; i < m_iIndxNum; i++)
	{
		cout << " " << m_iPIndxes[i];
	}
	cout << endl;
}

void Ice_SM::DebugLayer()
{	cout << "DebugLayer layer =";
	for(int i = 0; i < m_iIndxNum; i++)
	{
		cout << " " << m_ipLayeres[i];
	}
	cout << endl;
}