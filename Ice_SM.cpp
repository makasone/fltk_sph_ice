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
float Ice_SM::s_fItrStiffness;

Ice_SM::Ice_SM() : rxShapeMatching()
{
	m_pPrePos = 0;

	m_fpAlphas = 0;
	m_fpBetas =	0;
	
	m_ipLayeres = 0;
	m_ipErea = 0;
}

Ice_SM::Ice_SM(int obj) : rxShapeMatching(obj)
{
	m_mtrx3Apq = rxMatrix3(0.0);
	m_mtrxBeforeU.makeIdentity();

	m_fDefAmount = 0.0f;

	m_pPrePos = new float[maxNum*SM_DIM];

	m_fpAlphas = new float[maxNum];
	m_fpBetas =	new float[maxNum];

	m_ipLayeres = new int[maxNum];
	m_ipErea = new EreaData[maxNum];
}

Ice_SM::Ice_SM(int obj, int prtNum) : rxShapeMatching(obj, prtNum)
{
	m_mtrx3Apq = rxMatrix3(0.0);
	m_mtrxBeforeU.makeIdentity();

	m_fDefAmount = 0.0f;

	m_pPrePos = new float[maxNum*SM_DIM];

	m_fpAlphas = new float[maxNum];
	m_fpBetas =	new float[maxNum];

	m_ipLayeres = new int[maxNum];
	m_ipErea = new EreaData[maxNum];
}

/*!
 * �R�s�[�R���X�g���N�^
 */
Ice_SM::Ice_SM(const Ice_SM& copy) : rxShapeMatching(copy)
{
	Copy(copy);
}

/*!
 * �f�X�g���N�^
 */
Ice_SM::~Ice_SM()
{
	ReleaseMemory();
}

//������Z�q�ŃR�s�[
Ice_SM& Ice_SM::operator=(const Ice_SM& copy)
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
void Ice_SM::ReleaseMemory()
{
	if(m_pPrePos != 0){
		delete[] m_pPrePos;
		m_pPrePos = 0;
	}

	if(m_fpAlphas != 0){
		delete[] m_fpAlphas;
		m_fpAlphas = 0;
	}

	if(m_fpBetas != 0){
		delete[] m_fpBetas;
		m_fpBetas = 0;
	}

	if(m_ipLayeres != 0){
		delete[] m_ipLayeres;
		m_ipLayeres = 0;
	}

	if(m_ipErea != 0){
		delete[] m_ipErea;
		m_ipErea = 0;
	}
}

//�f�[�^�̃R�s�[
void Ice_SM::Copy(const Ice_SM& copy)
{
	m_mtrx3Apq = copy.GetApq();
	m_mtrxBeforeU.makeIdentity();	//TODO: �R�s�[���ĂȂ��@���Ȃ��Ă����H

	m_iIndxNum = copy.GetIndxNum();
	m_fDefAmount = copy.GetDefAmount();

	m_pPrePos = new float[maxNum*SM_DIM];

	m_fpAlphas = new float[maxNum];
	m_fpBetas =	new float[maxNum];

	m_ipLayeres = new int[maxNum];

	m_ipErea = new EreaData[maxNum];

	//�e���̃R�s�[
	for(int i = 0; i < maxNum; i++){
		for(int j = 0; j < SM_DIM; j++){
			m_pPrePos[i*SM_DIM+j] = copy.GetPrePos(i)[j];
		}

		m_fpAlphas[i] = copy.GetAlphas(i);
		m_fpBetas[i] = copy.GetBetas(i);

		m_ipLayeres[i] = copy.GetLayer(i);
		m_ipErea[i] = copy.erea(i);
	}

	//�ߖT�N���X�^���̃R�s�[
	m_ivNeighborFeatureCluster.clear();
	m_ivNeighborFeatureCluster.shrink_to_fit();
	int size = copy.neighborFeatureClusterNum();

	for(int i = 0; i < size; i++){
		AddNeighborFeatureCluster(copy.neighborFeatureCluster(i));
	}
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
void Ice_SM::InitGPU(const vector<Ice_SM*>& ice_sm, float* d_pos, float* d_vel, int prtNum, int MAXCLUSTER)
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
	if(m_iNumVertices > m_iIndxNum){	m_iIndxNum = m_iNumVertices;	}

	//�d�S�ʒu�̍X�V
	CalcOrgCm();

	//�e���q�̏�������̈���̍X�V
	ClassifyAllOrgParticle();

	//�g���ĂȂ�
	//m_iLinearDeformation.push_back(0);
	//m_iVolumeConservation.push_back(0);

	////�ό`�s��Aqq�̍X�V
	//Vec3 p, q;
	//for(int i = 0; i < m_iIndxNum; ++i)
	//{
	//	if( CheckHole(i) ){	continue;	}

	//	q = Vec3(m_pOrgPos[i*SM_DIM+0], m_pOrgPos[i*SM_DIM+1], m_pOrgPos[i*SM_DIM+2])-m_vec3OrgCm;
	//	double m = m_pMass[i];

	//	m_mtrx3AqqInv(0,0) += m*q[0]*q[0];
	//	m_mtrx3AqqInv(0,1) += m*q[0]*q[1];
	//	m_mtrx3AqqInv(0,2) += m*q[0]*q[2];
	//	m_mtrx3AqqInv(1,0) += m*q[1]*q[0];
	//	m_mtrx3AqqInv(1,1) += m*q[1]*q[1];
	//	m_mtrx3AqqInv(1,2) += m*q[1]*q[2];
	//	m_mtrx3AqqInv(2,0) += m*q[2]*q[0];
	//	m_mtrx3AqqInv(2,1) += m*q[2]*q[1];
	//	m_mtrx3AqqInv(2,2) += m*q[2]*q[2];
	//}

	//m_mtrx3AqqInv = m_mtrx3AqqInv.Inverse();

	////Q�̒ǉ�
	//m_vvec3OrgQ.resize(m_iNumVertices);
	//for(int i = 0; i < m_iNumVertices;++i)
	//{
	//	m_vvec3OrgQ[i] = m_vOrgPos[i]-m_vec3OrgCm;
	//}

	////�O�t���[���̈ʒu�X�V
	//for(int i = 0; i < m_iIndxNum; ++i)
	//{
	//	m_pPrePos[i*SM_DIM+0] = m_pCurPos[i*SM_DIM+0];
	//	m_pPrePos[i*SM_DIM+1] = m_pCurPos[i*SM_DIM+1];
	//	m_pPrePos[i*SM_DIM+2] = m_pCurPos[i*SM_DIM+2];
	//}

	//m_vec3PreCm = m_vec3OrgCm;
}

void Ice_SM::AddVertex(const Vec3 &pos, const Vec3& vel, double mass, int pIndx)
{
	rxShapeMatching::AddVertex(pos, vel, mass, pIndx);

	//�ő�Y���ԍ��̍X�V
	if(m_iNumVertices > m_iIndxNum){	m_iIndxNum = m_iNumVertices;	}
}

int Ice_SM::AddAnotherClusterVertex(const Vec3& orgPos, const Vec3& curPos, const Vec3& vel, double mass, int pIndx, double alpha, double beta, int layer)
{
	int indx = rxShapeMatching::AddVertex(orgPos, curPos, vel, mass, pIndx);

	if(indx != -1){
		SetAlphas(indx, alpha);
		SetBetas(indx, beta);
		SetLayer(indx, layer);
	}

	//�ő�Y���ԍ��̍X�V
	if(m_iNumVertices > m_iIndxNum){	m_iIndxNum = m_iNumVertices;	}

	return indx;
}

void Ice_SM::Remove(int indx)
{
	rxShapeMatching::Remove(indx);

	m_ipLayeres[indx] = -1;
	m_ipErea[indx] = EreaData();

	m_fpAlphas[indx] = -1;
	m_fpBetas[indx] = -1;

	//m_iLinearDeformation.erase	( m_iLinearDeformation	.begin() + indx);
	//m_iVolumeConservation.erase	( m_iVolumeConservation	.begin() + indx);
}

void Ice_SM::Clear()
{
	int indx = m_iIndxNum;

	rxShapeMatching::Clear();

	for(int i = 0; i < m_iIndxNum; i++)
	{
		m_fpAlphas[i] = 0.0f;
		m_fpBetas[i] = 0.0f;

		m_ipLayeres[i] = -1;
		m_ipErea[i] = EreaData();
	}
	
	ClearNeighborFeaturceCluster();

	//m_iLinearDeformation	.clear();
	//m_iVolumeConservation	.clear();
}

//TODO: ����Ȃ��K���Ȏ����@�Ȃ�ׂ��g��Ȃ�
bool Ice_SM::CheckIndx(int pIndx) const 
{//	cout << __FUNCTION__ << endl;

	//vector<int>::iterator check = find( m_iPIndxes.begin(), m_iPIndxes.end(), pIndx);
	//if( check == m_iPIndxes.end() )	return false;
	//return true;

	for(int i = 0; i < m_iIndxNum; i++)
	{
		if(m_iPIndxes[i] == pIndx){	return true;	}
	}

	return false;
}

//TODO: ����Ȃ��K���Ȏ����@�Ȃ�ׂ��g��Ȃ�
int	Ice_SM::SearchIndx(int pIndx) const
{
	////vector find�o�[�W�����@�߂��Ⴍ����d��
	//vector<int>::iterator begin = m_iPIndxes.begin();
	//vector<int>::iterator end = m_iPIndxes.end();
	//vector<int>::iterator check = find(begin, end, pIndx);

	//if(check == end)	return MAXINT;
	//return m_iPIndxes.size() - (end - check);

	////vector binary_search�o�[�W����
	//�\�[�g�ł��Ȃ��̂ł���
	//vector<int>::iterator begin = m_iPIndxes.begin();
	//vector<int>::iterator end = m_iPIndxes.end();

	//if(!binary_search(begin, end, pIndx)) return MAXINT;
	//return *( std::lower_bound(begin, end, pIndx) );

	//�z��o�[�W����
	for(int i = 0; i < m_iIndxNum; i++)
	{
		if(m_iPIndxes[i] == pIndx){	return i;	}
	}

	return MAXINT;

	////�z��{algorizm �ȉ��̕��@�ł͂ł��Ȃ�
	//int* last = m_iPIndxes + (int)m_iIndxNum;
	//int* result = find(m_iPIndxes, last, pIndx);

	//if(last == result) return MAXINT;
	//return *result;
}

void Ice_SM::ResetFinalParamPointer(unsigned clusterNum)
{
	for(unsigned cIndx = 0; cIndx < clusterNum; cIndx++)
	{
		for(unsigned dim = 0; dim < SM_DIM; dim++)
		{
			s_pfSldPos[cIndx*SM_DIM+dim] = 0.0f;
			s_pfSldVel[cIndx*SM_DIM+dim] = 0.0f;
		}
	}
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
		if( CheckHole(i) ){	continue;	}
		//if(m_pFix[i]){	mass += m_pMass[i]*300;	}	// �Œ�_�̎��ʂ�傫������
		mass += m_pMass[i];
	}

	Vec3 cm(1 / mass * m_vec3NowCm);	// ���݂̏d�S
	Vec3 q(0.0);

	//Apq�̍s�񎮂����߁C���]���邩�𔻒�
	//�s����ȏꍇ�������̂Ł~
	if( m_mtrx3Apq.Determinant() < 0.0 && m_iNumVertices >= 4)
	{
		//cout << "before det < 0" << endl;
		//�P�@�����𔽓]
		m_mtrx3Apq(0,2) = -m_mtrx3Apq(0,2);
		m_mtrx3Apq(1,2) = -m_mtrx3Apq(1,2);
		m_mtrx3Apq(2,2) = -m_mtrx3Apq(2,2);

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
	}

	rxMatrix3 R, S;
	//PolarDecomposition(m_mtrx3Apq, R, S, m_mtrxBeforeU);
	PolarDecomposition(m_mtrx3Apq, R, S);

	if(m_bLinearDeformation)
	{
		//// Linear Deformations
		//rxMatrix3 A(m_mtrx3Apq * m_mtrx3AqqInv);	// A = Apq*Aqq^-1

		////�̐ϕۑ��̂��߂Ɂ�(det(A))�Ŋ���
		//if(m_bVolumeConservation){
		//	double det = fabs(A.Determinant());
		//	if(det > RX_FEQ_EPS){
		//		det = 1.0/sqrt(det);
		//		if(det > 2.0) det = 2.0;
		//		A *= det;
		//	}
		//}

		//// �ڕW���W���v�Z���C���݂̒��_���W���ړ�
		//for(int i = 0; i <  m_iIndxNum; ++i)
		//{
		//	if( CheckHole(i) ){	continue;	}
		//	if(m_pFix[i]) continue;

		//	// ��]�s��R�̑���̍s��RL=��A+(1-��)R���v�Z
		//	int pIndx = m_iPIndxes[i] * SM_DIM;
		//	int cIndx = i*SM_DIM;

		//	for(int j = 0; j < SM_DIM; j++)
		//	{
		//		q[j] = m_pOrgPos[cIndx+j]-m_vec3OrgCm[j];
		//	}

		//	Vec3 gp(R*q+cm);

		//	for(int j = 0; j < SM_DIM; j++)
		//	{
		//		m_pCurPos[cIndx+j] += (gp[j]-m_pCurPos[cIndx+j]) * m_fpAlphas[i];
		//	}
		//}

		m_fDefAmount = 0.0;

		// �ڕW���W���v�Z���C���݂̒��_���W���ړ�
		for(int i = 0; i < m_iIndxNum; ++i)
		{
			if( CheckHole(i) ){	continue;	}
			if(m_pFix[i]) continue;

			int cIndx = i*SM_DIM;

			// ��]�s��R�̑���̍s��RL=��A+(1-��)R���v�Z
			for(int j = 0; j < SM_DIM; j++)
			{
				q[j] = m_pOrgPos[cIndx+j]-m_vec3OrgCm[j];
			}

			Vec3 gp(R*q+cm);

			for(int j = 0; j < SM_DIM; j++)
			{
				float defAmount = (gp[j]-m_pCurPos[cIndx+j]) * m_fpAlphas[i];
				m_pCurPos[cIndx+j] += defAmount;
				m_fDefAmount += abs(defAmount);
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
		//// �̐ϕۑ��̂��߂Ɂ�(det(A))�Ŋ���
		////if(m_bVolumeConservation){
			////beta == 0��O��ɂ��Ă���̂ŁCA���g��Ȃ��Ƃ��Ă���
			//rxMatrix3 A(Apq*Aqq.Inverse());
			//double det = fabs(A.Determinant());

			//if(det > RX_FEQ_EPS){
			//	//det = 1.0/sqrt(det);
			//	det = pow(det, 1.0/3.0);	//�_���I�ɂ͂��������������̂ł́H
			//	if(det > 2.0) det = 2.0;
			//	A *= det;
			//}
		////}

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
				float defAmount = (gp[j]-m_pCurPos[cIndx+j]) * m_fpAlphas[i];
				m_pCurPos[cIndx+j] += defAmount;

				m_fDefAmount += abs(defAmount);
			}
		}
	}
}

void Ice_SM::ShapeMatchingSelected(int selected)
{
	if(m_iNumVertices <= 1) return;

	rxMatrix3 Apq(0.0), Aqq(0.0);
	Vec3 p, q;
	Vec3 cm(0.0), cm_org(0.0);	// �d�S
	double mass = 0.0;			// ������
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
	m_vec3NowCm = cm;
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
	}

	PolarDecomposition(Apq, R, S);

	m_fDefAmount = 0.0f;

	if(m_bLinearDeformation)
	{
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
				float defAmount = (gp[j]-m_pCurPos[cIndx+j]) * m_fpAlphas[i];

				m_pCurPos[cIndx+j] += defAmount;
				m_fDefAmount += abs(defAmount);
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
			m_pVel[cIndx] = (m_pCurPos[cIndx] - s_pfPrtPos[pIndx+j]) * dt1/* + m_v3Gravity[j] * dt*/;
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
 * �����p�@�ő̂̍ŏI�ʒu���e���q�ɔ��f�@���x�͔��f���Ȃ�
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
	//���������Ȃ�d���Ȃ邪�C����͂���݂���
	double res = 0.9;	// �����W��
	double dt = m_dDt;
	for(int i = 0; i < m_iIndxNum; ++i)
	{
		if( CheckHole(i) ){	continue;	}
		if(m_pFix[i]) continue;

		int pIndx =  m_iPIndxes[i]*SM_DIM;
		int cIndx = i*SM_DIM;

		//if(m_pCurPos[cIndx+0] < m_v3Min[0] || m_pCurPos[cIndx+0] > m_v3Max[0]){
		//	m_pCurPos[cIndx+0] = s_pfSldPos[pIndx+0]-s_pfSldVel[pIndx+0]*dt*res;
		//	m_pCurPos[cIndx+1] = s_pfSldPos[pIndx+1];
		//	m_pCurPos[cIndx+2] = s_pfSldPos[pIndx+2];
		//}
		//if(m_pCurPos[cIndx+1] < m_v3Min[1] || m_pCurPos[cIndx+1] > m_v3Max[1]){
		//	m_pCurPos[cIndx+1] = s_pfSldPos[pIndx+1]-s_pfSldVel[pIndx+1]*dt*res;
		//	m_pCurPos[cIndx+0] = s_pfSldPos[pIndx+0] ;
		//	m_pCurPos[cIndx+2] = s_pfSldPos[pIndx+2];
		//}
		//if(m_pCurPos[cIndx+2] < m_v3Min[2] || m_pCurPos[cIndx+2] > m_v3Max[2]){
		//	m_pCurPos[cIndx+2] = s_pfSldPos[pIndx+2]-s_pfSldVel[pIndx+2]*dt*res;
		//	m_pCurPos[cIndx+0] = s_pfSldPos[pIndx+0];
		//	m_pCurPos[cIndx+1] = s_pfSldPos[pIndx+1];
		//}

		clamp(m_pCurPos, cIndx);
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
	m_vec3NowCm = cm;
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

		//m_fDefAmount = 0.0f;

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
				float defAmount = (gp[j]-m_pCurPos[cIndx+j]) * m_fpAlphas[i];
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
		int sldIndx = m_iPIndxes[i]*SM_DIM;

		for(int j = 0; j < SM_DIM; j++)
		{
			int cIndx = i*SM_DIM+j;
			m_pVel[cIndx] = (m_pCurPos[cIndx] - s_pfPrtPos[pIndx+j]) * dt1;
			//m_pVel[cIndx] = (s_pfSldPos[sldIndx+j] - s_pfPrtPos[pIndx+j]) * dt1;
		}
	}
}

//�O�X�e�b�v�ƌ��X�e�b�v�̈ʒu���玿�ʂ�����
void Ice_SM::CalcMass()
{
	float distSum = 0.0f;
	for(int i = 0; i < m_iIndxNum; i++)
	{
		float dist =	abs(m_pCurPos[i*SM_DIM+0] - m_pPrePos[i*SM_DIM+0])
					+	abs(m_pCurPos[i*SM_DIM+1] - m_pPrePos[i*SM_DIM+1])
					+	abs(m_pCurPos[i*SM_DIM+2] - m_pPrePos[i*SM_DIM+2]);
		
		distSum += dist*dist;
	}

	float average = 1.0f/(float)m_iIndxNum;
	float massSum = 0.0f;
	
	//���ړ��ʂƕό`�ʂ�臒l�Ƃ���
	if(m_fDefAmount > 0.1f)
	{
		for(int i = 0; i < m_iIndxNum; i++)
		{
			float dist =	abs(m_pCurPos[i*SM_DIM+0] - m_pPrePos[i*SM_DIM+0])
						+	abs(m_pCurPos[i*SM_DIM+1] - m_pPrePos[i*SM_DIM+1])
						+	abs(m_pCurPos[i*SM_DIM+2] - m_pPrePos[i*SM_DIM+2]);
			
			m_pMass[i] = 1.0f + dist * dist/distSum - average;
			massSum += m_pMass[i];

			//cout << "mass = " << m_pMass[i] << endl;
		}

		//cout << "massSum = " << massSum << " , vertexNum = " << m_iIndxNum << endl;
	}
	else
	{
		for(int i = 0; i < m_iIndxNum; i++)
		{
			m_pMass[i] = 1.0f;
		}
	}
}

//���݂̏d�S���v�Z
void Ice_SM::CalcOrgCm()
{
	//�d�S�̍X�V
	double massSum = 0.0;	// ������
	m_vec3OrgCm = Vec3(0.0);

	// �d�S���W�̌v�Z
	for(int i = 0; i < m_iIndxNum;++i)
	{
		if( CheckHole(i) ){	continue;	}

		double m = m_pMass[i];
		int pIndx = i*SM_DIM;

		massSum += m;
		m_vec3OrgCm += Vec3(m_pOrgPos[pIndx+0], m_pOrgPos[pIndx+1], m_pOrgPos[pIndx+2]) * m;
	}

	m_vec3OrgCm /= massSum;
}

//���ׂĂ�org���q���e�̈�ɕ���
void Ice_SM::ClassifyAllOrgParticle()
{
	for(int i = 0; i < m_iIndxNum;++i)
	{
		if( CheckHole(i) ){	continue;	}

		ClassifyOrgParticle(i);
	}
}

//���q���e�̈�ɕ��ށ@�܂��͂W���� ���Z�ʑ̂̏ꍇ����n~3�ŕ�������������
void Ice_SM::ClassifyOrgParticle(int indx)
{
	double X = GetOrgPos(indx)[0] - GetOrgCm()[0];
	double Y = GetOrgPos(indx)[1] - GetOrgCm()[1];
	double Z = GetOrgPos(indx)[2] - GetOrgCm()[2];

	//�W�p�^�[���ɕ���
	if(X >= 0.0 && Y >= 0.0 && Z >= 0.0){
		m_ipErea[indx] = EreaData(TOP_LEFT__FORE);
	}
	else if(X < 0.0 && Y >= 0.0 && Z >= 0.0){
		m_ipErea[indx] = EreaData(TOP_RIGHT_FORE);
	}
	else if(X >= 0.0 && Y >= 0.0 && Z < 0.0){
		m_ipErea[indx] = EreaData(TOP_LEFT__BACK);
	}
	else if(X < 0.0 && Y >= 0.0 && Z < 0.0){
		m_ipErea[indx] = EreaData(TOP_RIGHT_BACK);
	}

	else if(X >= 0.0 && Y < 0.0 && Z >= 0.0){
		m_ipErea[indx] = EreaData(BOT_LEFT__FORE);
	}
	else if(X < 0.0 && Y < 0.0 && Z >= 0.0){
		m_ipErea[indx] = EreaData(BOT_RIGHT_FORE);
	}
	else if(X >= 0.0 && Y < 0.0 && Z < 0.0){
		m_ipErea[indx] = EreaData(BOT_LEFT__BACK);
	}
	else if(X < 0.0 && Y < 0.0 && Z < 0.0){
		m_ipErea[indx] = EreaData(BOT_RIGHT_BACK);
	}

	return;
}

//�̈�Ԃ̋������v�Z
unsigned Ice_SM::CalcEreaDistance(const EreaData& ereaA, const EreaData& ereaB)
{
	return abs(ereaA.x - ereaB.x) + abs(ereaA.y - ereaB.y) + abs(ereaA.z - ereaB.z);
}

//�̈�͑S��3���̕ϐ��Ŏw��
Ice_SM::EreaData Ice_SM::IntToEreaData(int erea)
{
	if(erea >= 1000){
		cout << __FUNCTION__ << " Error ����͂R���܂�" << endl;
	}

	return EreaData(erea/100, erea%100 /10, erea%10); 
}

int Ice_SM::EreaDataToInt(EreaData erea)
{
	return erea.x * 100 + erea.y * 10 + erea.z;
}

void Ice_SM::UpdatePrePos(int pIndx)
{
	pIndx = pIndx*SM_DIM;

	m_pPrePos[pIndx+0] = m_pCurPos[pIndx+0];
	m_pPrePos[pIndx+1] = m_pCurPos[pIndx+1];
	m_pPrePos[pIndx+2] = m_pCurPos[pIndx+2];
}

void Ice_SM::SetPrePos(int pIndx, const Vec3& nowPos)
{
	pIndx = pIndx*SM_DIM;
	m_pPrePos[pIndx+0] = nowPos[0];
	m_pPrePos[pIndx+1] = nowPos[1];
	m_pPrePos[pIndx+2] = nowPos[2];
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
		if(GetParticleIndx(j) == MAXINT){	continue;	}

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

void Ice_SM::DebugClusterInfo()
{
	for(unsigned i = 0; i < GetIndxNum(); i++ )
	{
		if(CheckHole(i)){
			continue;
		}

		int pIndx = GetParticleIndx(i);
		int ereaNum = Ice_SM::EreaDataToInt(erea(i));

		cout << "	i = "  << i 
			<< " pIndx = " << pIndx
			<< " layer = " << GetLayer(i)
			<< " erea   = " << ereaNum;
		cout << endl;
	}
}

void Ice_SM::DebugNeighborFeatureClusterInfo()
{
	cout << "	NeighborClusterInfo:size=" << neighborFeatureClusterNum() << ",";
	for(unsigned i = 0; i < neighborFeatureClusterNum(); i++){
		cout << " " << neighborFeatureCluster(i);
	}
	cout << endl;
}