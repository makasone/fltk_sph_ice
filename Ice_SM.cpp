/*!
  @file Ice_SM.h
	
  @brief ShapeMatching�@(GPU����)
  @ref M. Muller et al., "Meshless deformations based on shape matching", SIGGRAPH2005. 
 
  @author Ryo Nakasone
  @date 2014-7
*/

#include "Ice_SM.h"

const float* Ice_SM::s_pfPrtPos;
const float* Ice_SM::s_pfPrtVel;

float* Ice_SM::s_pfSldPos;
float* Ice_SM::s_pfSldVel;

//�f�o�C�X�|�C���^
float* Ice_SM::sd_PrtPos;	
float* Ice_SM::sd_PrtVel;

float* Ice_SM::d_OrgPos;
float* Ice_SM::d_CurPos;

float* Ice_SM::d_OrgCm;
float* Ice_SM::d_CurCm;

float* Ice_SM::d_Apq;

float* Ice_SM::d_Mass;
float* Ice_SM::d_Vel;

bool* Ice_SM::d_Fix;

int* Ice_SM::d_PIndxes;
int* Ice_SM::d_IndxSet;					//�N���X�^�̃f�[�^�̊J�n�Y���ƏI���Y����ۑ�

int Ice_SM::s_vertNum;
int Ice_SM::s_vertSum;					//�S�N���X�^�Ɋ܂܂�闱�q�̑���

int Ice_SM::s_iIterationNum;


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

void Ice_SM::InitFinalParamPointer(int vrtxNum)
{
	s_pfSldPos = new float[vrtxNum*SM_DIM];
	s_pfSldVel = new float[vrtxNum*SM_DIM];

	for(int i = 0; i < vrtxNum; ++i)
	{
		int pIndx = i*4;
		int cIndx = i*SM_DIM;

		for(int j = 0; j < SM_DIM; j++)
		{
			s_pfSldPos[cIndx+j] = s_pfPrtPos[pIndx+j];
			s_pfSldVel[cIndx+j] = s_pfPrtVel[pIndx+j];
		}
	}
}

//GPU�v�Z�̂��߂̏�����
void Ice_SM::InitGPU(const vector<Ice_SM*>& ice_sm, float* d_pos, float* d_vel, int prtNum)
{
	//�f�o�C�X�|�C���^�̃A�h���X��ۑ�
	sd_PrtPos = d_pos;
	sd_PrtVel = d_vel;

	s_vertNum = prtNum;

	//�f�o�C�X���̃��������m��
	//�i�ő�N���X�^���j�~�i�N���X�^���ۑ��ł���ő嗱�q���j�@�Ń��������m�ہD
	cudaMalloc((void**)&d_OrgPos,	sizeof(float) * MAXCLUSTER * MAXPARTICLE * SM_DIM);
	cudaMalloc((void**)&d_CurPos,	sizeof(float) * MAXCLUSTER * MAXPARTICLE * SM_DIM);
	cudaMalloc((void**)&d_Vel,		sizeof(float) * MAXCLUSTER * MAXPARTICLE * SM_DIM);

	cudaMalloc((void**)&d_OrgCm,	sizeof(float) * MAXCLUSTER * SM_DIM);
	cudaMalloc((void**)&d_CurCm,	sizeof(float) * MAXCLUSTER * SM_DIM);

	cudaMalloc((void**)&d_Apq,		sizeof(float) * MAXCLUSTER * 9);	//TODO: ���������Ă��Ȃ��̂ɒ���

	cudaMalloc((void**)&d_Mass,		sizeof(float) * MAXCLUSTER * MAXPARTICLE);
	cudaMalloc((void**)&d_Fix,		sizeof(bool)  * MAXCLUSTER * MAXPARTICLE);
	cudaMalloc((void**)&d_PIndxes,	sizeof(int)   * MAXCLUSTER * MAXPARTICLE);

	cudaMalloc((void**)&d_IndxSet,	sizeof(int)   * MAXCLUSTER * 2);

	//CPU�̃f�[�^���P�����z��ɃR�s�[
	//�z�X�g���̃f�[�^���f�o�C�X���֓]�����ď�����
		//��C�ɑ�ʂ̃������͊m�ۂł��Ȃ��̂ŁA�������ăR�s�[����
	float* oPoses = new float[MAXPARTICLE * SM_DIM];
	float* cPoses = new float[MAXPARTICLE * SM_DIM];
	float* veles  = new float[MAXPARTICLE * SM_DIM];

	float* orgCm = new float[MAXCLUSTER * SM_DIM];
	float* curCm = new float[MAXCLUSTER * SM_DIM];

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

		Vec3 orgCmVec = sm->GetOrgCm();
		Vec3 curCmVec = sm->GetCm();
		
		for(int j = 0; j < SM_DIM; j++)
		{
			orgCm[i*SM_DIM+j] = orgCmVec[j];
			curCm[i*SM_DIM+j] = curCmVec[j];
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

	//������
	cudaMemcpy(d_OrgCm,	orgCm,	sizeof(float) * MAXCLUSTER * SM_DIM, cudaMemcpyHostToDevice);
	cudaMemcpy(d_CurCm,	curCm,	sizeof(float) * MAXCLUSTER * SM_DIM, cudaMemcpyHostToDevice);

//�f�o�b�O
	////�ꎞ�z��̃��Z�b�g
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

	delete[] orgCm;
	delete[] curCm;

	delete[] mass;
	delete[] fix;
	delete[] indx;
	delete[] set;
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
void Ice_SM::ShapeMatchingUsePath()
{
	if(m_iNumVertices <= 1) return;

	double mass = 0.0;	// ������

	// �d�S���W�̌v�Z
	for(int i = 0; i < m_iIndxNum; ++i){
		//if( CheckHole(i) ){	continue;	}
		//if(m_pFix[i]){	mass += m_pMass[i]*300;	}	// �Œ�_�̎��ʂ�傫������
		mass += m_pMass[i];
	}

	Vec3 cm(1 / mass * m_vec3NowCm);	// ���݂̏d�S
	Vec3 q(0.0);

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

	rxMatrix3 R, S;
	//PolarDecomposition(m_mtrx3Apq, R, S, m_mtrxBeforeU);
	PolarDecomposition(m_mtrx3Apq, R, S);

	if(m_bLinearDeformation)
	{
		// Linear Deformations
		rxMatrix3 A(m_mtrx3Apq * m_mtrx3AqqInv);	// A = Apq*Aqq^-1

		//�̐ϕۑ��̂��߂Ɂ�(det(A))�Ŋ���
		if(m_bVolumeConservation){
			double det = fabs(A.Determinant());
			if(det > RX_FEQ_EPS){
				det = 1.0/sqrt(det);
				if(det > 2.0) det = 2.0;
				A *= det;
			}
		}

		// �ڕW���W���v�Z���C���݂̒��_���W���ړ�
		for(int i = 0; i <  m_iIndxNum; ++i)
		{
			if( CheckHole(i) ){	continue;	}
			if(m_pFix[i]) continue;

			// ��]�s��R�̑���̍s��RL=��A+(1-��)R���v�Z
			int pIndx = m_iPIndxes[i] * SM_DIM;
			int cIndx = i*SM_DIM;

			for(int j = 0; j < SM_DIM; j++)
			{
				q[j] = m_pOrgPos[cIndx+j]-m_vec3OrgCm[j];
			}

			Vec3 gp(R*q+cm);

			for(int j = 0; j < SM_DIM; j++)
			{
				//m_pCurPos[cIndx+j] += (gp[j]-s_pfSldPos[pIndx+j])  * m_fpAlphas[i];
				m_pCurPos[cIndx+j] += (gp[j]-m_pCurPos[cIndx+j]) * m_fpAlphas[i];
			}
		}
	}

	//���E�����@���ꂪ�Ȃ��ƌ��݈ʒu�����E�𖳎�����̂ŁC���q�����蔲����
	for(int i = 0; i < m_iIndxNum; ++i)
	{
		clamp(m_pCurPos, i*SM_DIM);
	}
}

/*!
 * �O��
 *  - �d�͂Ƌ��E�ǂ���̗͂̉e��
 * @param[in] dt �^�C���X�e�b�v��
 */
void Ice_SM::calExternalForces()
{
	double dt = m_dDt;

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
			
			m_pCurPos[jcIndx] = s_pfPrtPos[jpIndx] + s_pfPrtVel[jpIndx]*dt;
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
void Ice_SM::ShapeMatchingSolid()
{
	if(m_iNumVertices <= 1) return;

	rxMatrix3 Apq(0.0), Aqq(0.0);
	Vec3 p, q;
	Vec3 cm(0.0), cm_org(0.0);	// �d�S
	double mass = 0.0;	// ������
	rxMatrix3 R, S;

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
	cm_org = m_vec3OrgCm;

	// Apq = ��mpq^T
	// Aqq = ��mqq^T
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

	////Apq�̍s�񎮂����߁C���]���邩�𔻒�
	////�s����ȏꍇ�������̂Ł~
	//if( Apq.Determinant() < 0.0 && m_iNumVertices >= 4)
	//{
	//	//cout << "before det < 0" << endl;
	//	//�P�@�����𔽓]
	//	Apq(0,2) = -Apq(0,2);
	//	Apq(1,2) = -Apq(1,2);
	//	Apq(2,2) = -Apq(2,2);

	//	//�Q�@a2��a3������
	//	//double tmp;
	//	//tmp = Apq(0,2);
	//	//Apq(0,2) = Apq(0,1);
	//	//Apq(0,1) = tmp;

	//	//tmp = Apq(1,2);
	//	//Apq(1,2) = Apq(1,1);
	//	//Apq(1,1) = tmp;

	//	//tmp = Apq(2,2);
	//	//Apq(2,2) = Apq(2,1);
	//	//Apq(2,1) = tmp;
	//}

	//PolarDecomposition(Apq, R, S, m_mtrxBeforeU);
	PolarDecomposition(Apq, R, S);

	if(m_bLinearDeformation)
	{
		// Linear Deformations
		rxMatrix3 A(Apq*Aqq.Inverse());

		// �̐ϕۑ��̂��߂Ɂ�(det(A))�Ŋ���
		if(m_bVolumeConservation){
			double det = fabs(A.Determinant());
			if(det > RX_FEQ_EPS){
				det = 1.0/sqrt(det);
				if(det > 2.0) det = 2.0;
				A *= det;
			}
		}

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
				m_pCurPos[cIndx+j] += (gp[j]-m_pCurPos[cIndx+j]) * m_fpAlphas[i];
			}
		}
	}
}

//----------------------------------------------�\���b�h��---------------------------------------------

/*!
 * ���x�ƈʒu�̍X�V
 *  - �V�����ʒu�ƌ��݂̈ʒu���W���瑬�x���Z�o
 * @param[in] dt �^�C���X�e�b�v��
 */
void Ice_SM::integrate(double dt)
{
	double dt1 = 1.0/dt;
	for(int i = 0; i < m_iIndxNum; ++i)
	{
		if( CheckHole(i) ){	continue;	}
		int pIndx = m_iPIndxes[i]*4;

		for(int j = 0; j < SM_DIM; j++)
		{
			int cIndx = i*SM_DIM+j;
			m_pVel[cIndx] = (m_pCurPos[cIndx] - s_pfPrtPos[pIndx+j]) * dt1 + m_v3Gravity[j] * dt * 0.1f;	//TODO::0.1f������Ƃ��܂������Č�����
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
 * �ŏI�ʒu���e�N���X�^�ɔ��f
 * @param[in] dt �^�C���X�e�b�v��
 */
void Ice_SM::calExternalForcesIteration()
{
	// �N���X�^���̗��q�̈ʒu���C���q�̈ʒu�ōX�V
	for(int i = 0; i < m_iIndxNum; ++i)
	{
		if( CheckHole(i) ){	continue;	}
		int pIndx = m_iPIndxes[i]*SM_DIM;
		int cIndx = i*SM_DIM;

		for(int j = 0; j < SM_DIM; j++)
		{
			int jpIndx = pIndx+j;
			int jcIndx = cIndx+j;

			m_pCurPos[jcIndx] = s_pfSldPos[jpIndx];
		}
	}

	// ���E�ǂ̉e��
	//���������Ȃ�d���Ȃ邪�C���肷��
	double res = 0.9;	// �����W��
	for(int i = 0; i < m_iIndxNum; ++i)
	{
	//	//if( CheckHole(i) ){	continue;	}
	//	////if(m_pFix[i]) continue;
	//	//Vec3 &p = m_pCurPos[i];

	//	//Vec3 &v = m_pVel[i];
	//	//if(np[0] < m_v3Min[0] || np[0] > m_v3Max[0]){
	//	//	np[0] = p[0]-v[0]*dt*res;
	//	//	np[1] = p[1];
	//	//	np[2] = p[2];
	//	//}
	//	//if(np[1] < m_v3Min[1] || np[1] > m_v3Max[1]){
	//	//	np[1] = p[1]-v[1]*dt*res;
	//	//	np[0] = p[0] ;
	//	//	np[2] = p[2];
	//	//}
	//	//if(np[2] < m_v3Min[2] || np[2] > m_v3Max[2]){
	//	//	np[2] = p[2]-v[2]*dt*res;
	//	//	np[0] = p[0];
	//	//	np[1] = p[1];
	//	//}

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
void Ice_SM::ShapeMatchingIteration()
{
	if(m_iNumVertices <= 1) return;

	rxMatrix3 Apq(0.0), Aqq(0.0);
	Vec3 p, q;
	Vec3 cm(0.0), cm_org(0.0);	// �d�S
	double mass = 0.0;	// ������
	rxMatrix3 R, S;

	// �d�S���W�̌v�Z
	for(int i = 0; i < m_iIndxNum; ++i)
	{
		if( CheckHole(i) ){	continue;	}

		double m = m_pMass[i];

		if(m_pFix[i]) m *= 300.0;	// �Œ�_�̎��ʂ�傫������
		mass += m;

		int pIndx = m_iPIndxes[i] * SM_DIM;

		for(int j = 0; j < SM_DIM; j++)
		{
			cm[j] += s_pfSldPos[pIndx+j]*m;
		}
	}

	cm /= mass;
	cm_org = m_vec3OrgCm;

	// Apq = ��mpq^T
	// Aqq = ��mqq^T
	for(int i = 0; i < m_iIndxNum; ++i)
	{
		if( CheckHole(i) ){	continue;	}

		int pIndx = m_iPIndxes[i] * SM_DIM;
		int cIndx = i*SM_DIM;

		for(int j = 0; j < SM_DIM; j++)
		{
			p[j] = s_pfSldPos[pIndx+j]-cm[j];
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

	//Apq�̍s�񎮂����߁C���]���邩�𔻒�
	//�s����ȏꍇ�������̂Ł~
	if( Apq.Determinant() < 0.0 && m_iNumVertices >= 4)
	{
		//cout << "before det < 0" << endl;
		//�P�@�����𔽓]
		Apq(0,2) = -Apq(0,2);
		Apq(1,2) = -Apq(1,2);
		Apq(2,2) = -Apq(2,2);

		////�Q�@a2��a3������
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
	}

	//PolarDecomposition(Apq, R, S, m_mtrxBeforeU);
	PolarDecomposition(Apq, R, S);

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

			int pIndx = m_iPIndxes[i] * SM_DIM;
			int cIndx = i*SM_DIM;

			// ��]�s��R�̑���̍s��RL=��A+(1-��)R���v�Z
			for(int j = 0; j < SM_DIM; j++)
			{
				q[j] = m_pOrgPos[cIndx+j]-cm_org[j];
			}

			Vec3 gp(R*q+cm);

			for(int j = 0; j < SM_DIM; j++)
			{
				int jpIndx = pIndx+j;
				int jcIndx = cIndx+j;

				m_pCurPos[jcIndx] += (gp[j]-s_pfSldPos[jpIndx])*m_fpAlphas[i];
			}
		}
	}

	//���E����
	for(int i = 0; i < m_iIndxNum; ++i)
	{
		clamp(m_pCurPos, i*SM_DIM);
	}
}

/*!
 * ���x�ƈʒu�̍X�V
 *  - �V�����ʒu�ƌ��݂̈ʒu���W���瑬�x���Z�o
 * @param[in] dt �^�C���X�e�b�v��
 */
void Ice_SM::integrateIteration()
{
	double dt1 = 1.0 / m_dDt;
	for(int i = 0; i < m_iIndxNum; ++i)
	{
		if( CheckHole(i) ){	continue;	}
		int pIndx = m_iPIndxes[i]*4;

		for(int j = 0; j < SM_DIM; j++)
		{
			int cIndx = i*SM_DIM+j;
			m_pVel[cIndx] = (m_pCurPos[cIndx] - s_pfPrtPos[pIndx+j]) * dt1 + m_v3Gravity[j] * m_dDt * 0.1f;	//TODO::�p�����[�^0.1f���g��Ȃ��Ă������悤��
		}
	}
}

/*!
 *�@CPU�ŉ^���v�Z
 */
void Ice_SM::UpdateCPU()
{
	calExternalForces();		//���݈ʒu�̌v�Z
	ShapeMatchingSolid();		//���݈ʒu�̍X�V ���ʂ̌v�Z
	integrate(m_dDt);			//���x�̌v�Z

	//CalcDisplaceMentVectorCm();	//�d�S�ړ���
}

/*!
 *�@�p�X��p������������@�@CPU�ŉ^���v�Z
 */
void Ice_SM::UpdateUsePathCPU()
{
	calExternalForces();			//���݈ʒu�̌v�Z
	ShapeMatchingUsePath();			//���݈ʒu�̍X�V �p�X��p����������
	integrate(m_dDt);				//���x�̌v�Z

	//CalcDisplaceMentVectorCm();	//�d�S�ړ���
}

/*
 *	GPU��p�����^���v�Z
 */
void Ice_SM::UpdateGPU()
{
	LaunchShapeMatchingGPU(s_vertNum, sd_PrtPos, sd_PrtVel, d_OrgPos, d_CurPos, d_OrgCm, d_CurCm, d_Vel, d_PIndxes, d_IndxSet, 0.01);
}

void Ice_SM::UpdateUsePathGPU()
{
	LaunchShapeMatchingUsePathGPU(s_vertNum, sd_PrtPos, sd_PrtVel, d_OrgPos, d_CurPos, d_OrgCm, d_CurCm, d_Apq, d_Vel, d_PIndxes, d_IndxSet, 0.01);
}

/*
 *	GPU��p�����^���v�Z�@������������
 */
void Ice_SM::UpdateIterationGPU(float* sldPos, float* sldVel)
{
	LaunchShapeMatchingIterationGPU(s_vertNum, sd_PrtPos, sd_PrtVel, sldPos, sldVel, d_OrgPos, d_CurPos, d_OrgCm, d_CurCm, d_Vel, d_PIndxes, d_IndxSet, 0.01);
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