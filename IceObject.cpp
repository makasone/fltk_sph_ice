#include "IceObject.h"

float* IceObject::s_sphPrtPos;
float* IceObject::s_sphPrtVel;

float* IceObject::m_fInterPolationCoefficience;

//�f�o�C�X�|�C���^
float* IceObject::sd_sphPrtPos;
float* IceObject::sd_sphPrtVel;

float* IceObject::sd_sldPrtPos;	
float* IceObject::sd_sldPrtVel;

int IceObject::sm_particleNum;
int IceObject::sm_tetraNum;
int IceObject::sm_clusterNum;
int IceObject::sm_layerNum;
int IceObject::sm_maxParticleNum;

IceObject::IceObject(int pMaxNum, int cMaxNum, int tMaxNum, int prtNum, float* hp, float* hv, float* dp, float* dv, int layer, int maxParticleNum)
{
	InitIceObj(pMaxNum, cMaxNum, tMaxNum);

	SetParticleNum(prtNum);

	SetSPHHostPointer(hp, hv);						//CPU�����ŗp����|�C���^�̓o�^
	SetSPHDevicePointer(dp, dv);					//GPU�����ŗp����|�C���^�̓o�^

	SetSearchLayerNum(layer);						//�T�����C���[��
	SetMaxParticleNum(maxParticleNum);				//�ő嗱�q��

	InitTetra();									//�l�ʑ̂̏�����
}

IceObject::~IceObject()
{
}

//���ꂼ��̃N���X�E�ϐ��̏�����
void IceObject::InitIceObj(int pMaxNum, int cMaxNum, int tMaxNum)
{	cout << __FUNCTION__ << endl;
	//���̂̍\���̏�����
	m_iceStrct = new IceStructure(pMaxNum, cMaxNum, tMaxNum, sm_layerNum);

	//�^���v�Z���s��SM�N���X�^�̏�����

	//��ԏ����̂��߂̃p�����[�^�̏�����
	InitInterPolation();
}

void IceObject::InitSelectCluster()
{
	////�K���ɉ����̔{���őI��
	////TODO:: �������X���C�_�[�ŉςɂ���Ɩʔ�������
	//int SELECT = 7;

	//for(int cIndx = 0; cIndx < sm_particleNum; cIndx++)
	//{
	//	if(cIndx%SELECT != 1){		continue;	}
	//	m_iceStrct->UpdateMotionCalcCluster(cIndx, 1);
	//}

	//�ߖT�N���X�^���ŃN���X�^��I��
	//�N���X�^�̑S�̏W�����炠��N���X�^��I�сC���̃N���X�^�̋ߖT�Ɋ܂܂��N���X�^����菜��
	//������J��Ԃ��ăN���X�^�W�����炷�ׂẴN���X�^���I�΂��ƏI��

	//�N���X�^�W���̏�����
	vector<unsigned> clusters;
	for(int cIndx = 0; cIndx < sm_clusterNum; cIndx++)
	{
		clusters.push_back(cIndx);
	}

	//�I��
	while(clusters.size() != 0)
	{
		unsigned cIndx = *clusters.begin();

		m_iceStrct->UpdateMotionCalcCluster(cIndx, 1);

		//�ߖT�N���X�^����菜��
		for(int indx = 0; indx < m_iceSM[cIndx]->GetNumVertices(); indx++)
		{
			unsigned icIndx = (unsigned)m_iceSM[cIndx]->GetParticleIndx(indx);

			//stl�ō폜�́Cerase��remove��g�ݍ��킹�čs��
			clusters.erase(remove(clusters.begin(), clusters.end(), icIndx), clusters.end());  
		}

		clusters.erase(remove(clusters.begin(), clusters.end(), cIndx), clusters.end());  
	}
}

void IceObject::InitTetra()
{	cout << __FUNCTION__ << endl;

	IceTetrahedra &tetra = IceTetrahedra::GetInstance();		//�V���O���g���Ȏl�ʑ̊Ǘ��N���X��p�ӂ��Ă݂�
	tetra.InitTetra(s_sphPrtPos, sm_particleNum);

cout << __FUNCTION__ << " check1" << endl;

	m_iceStrct->SetTetraNum(tetra.GetTetraListSize());			//���l�ʑ̐���o�^

cout << __FUNCTION__ << " check2" << endl;

	//�e�l�ʑ̂Ɋ܂܂�闱�q���̃J�E���g
	for(unsigned i = 0; i < tetra.GetTetraListSize(); i++)
	{
		m_iceStrct->CountTetrahedra(i, tetra.GetTetraList(i));
	}

cout << __FUNCTION__ << " check3" << endl;

	//�������m��
	m_iceStrct->InitTetraInfo();

cout << __FUNCTION__ << " check4" << endl;

	//���q���������Ă���N���X�^���̔z����R�s�[
	int *PtoTNum = new int[sm_particleNum];

	for(int i = 0; i < sm_particleNum; i++)
	{
		PtoTNum[i] = m_iceStrct->GetPtoTNum(i);
	}

cout << __FUNCTION__ << " check5" << endl;

	//�l�ʑ̃f�[�^�o�^
	for(unsigned i = 0; i < tetra.GetTetraListSize(); i++)
	{
		m_iceStrct->SetTetraInfo(i, PtoTNum);
	}
	delete[] PtoTNum;

cout << __FUNCTION__ << " check6" << endl;

	//�ߖT�l�ʑ̃f�[�^�o�^
	for(unsigned i = 0; i < tetra.GetTetraListSize(); i++)
	{
		m_iceStrct->SetNeighborTetra(i, sm_layerNum);
	}

cout << __FUNCTION__ << " check7" << endl;

	//�f�o�b�O
	//m_iceObj->DebugTetraInfo();
}

//�N���X�^�̏�����
void IceObject::InitCluster(Vec3 boundarySpaceHigh, Vec3 boundarySpaceLow, float timeStep, int itr)
{	cout << __FUNCTION__ << endl;

	//�ϐ��̏�����
	for(vector<Ice_SM*>::iterator it = m_iceSM.begin(); it != m_iceSM.end(); ++it)
	{
		if(*it) delete *it;
	}
	
	sm_clusterNum = 0;	//�N���X�^���̏�����

	Ice_SM::SetIterationNum(itr);
	Ice_SM::SetPrtPointerPosAndVel(s_sphPrtPos, s_sphPrtVel);

	cout << __FUNCTION__ << ", check1" << endl;

	//�l�ʑ̃��X�g�����ɁC���q���ɃN���X�^�쐬
	for(int i = 0; i < sm_particleNum; ++i)
	{
		//�N���X�^������
		m_iceSM.push_back(new Ice_SM(sm_clusterNum));
		m_iceSM[sm_clusterNum]->SetSimulationSpace(boundarySpaceLow, boundarySpaceHigh);
		m_iceSM[sm_clusterNum]->SetTimeStep(timeStep);
		m_iceSM[sm_clusterNum]->SetCollisionFunc(0);
		m_iceSM[sm_clusterNum]->SetStiffness(1.0, 1.0);

		//�l�ʑ̃��X�g�����ɁC�N���X�^�֗��q��o�^
		SetClusterMoveInfo(i);

		sm_clusterNum++;
	}

	cout << __FUNCTION__ << ", check2" << endl;

	Ice_SM::InitFinalParamPointer(sm_clusterNum);

	////MakeClusterFromNeight();
	////MakeOneCluster();
	
	//TODO::���͂𐶂ނ��߂ɗ��q���ʂ�������

//�f�o�b�O
	//�N���X�^�Ɋ܂܂�闱�q
	//for(int i = 0; i < ICENUM; i++)
	//{
	//	cout << "pIndx = " << endl;

	//	for(int j = 0; j < m_sm_cluster[i]->GetNumVertices(); j++)
	//	{
	//		cout << " j = " << m_sm_cluster[i]->GetParticleIndx(j);
	//	}
	//	cout << endl;
	//}

	//DebugClusterInfo();
}

//�I�u�W�F�N�g�̍\���̏�����
void IceObject::InitStrct()
{	cout << __FUNCTION__ << endl;

	//�J�E���g
	for(int i = 0; i < sm_particleNum; i++)
	{
		int pNum = m_iceSM[i]->GetNumVertices();
		vector<int> pList;

		for(int j = 0; j < pNum; j++)
		{
			pList.push_back(m_iceSM[i]->GetParticleIndx(j));
		}
		m_iceStrct->CountClusterParticle(i, pList, pNum);
	}
	cout << __FUNCTION__ << " �J�E���g" << endl;

	//�������m��
	InitClusterInfo();
	cout << __FUNCTION__ << " �������m��" << endl;

	//���q���������Ă���N���X�^���̔z����R�s�[
	int *PtoCNum = new int[sm_particleNum];
	cout << __FUNCTION__ << " �z����R�s�[�@�������m��" << endl;

	for(int i = 0; i < sm_particleNum; i++)
	{
		PtoCNum[i] = GetPtoCNum(i);
	}
	cout << __FUNCTION__ << " �z����R�s�[" << endl;

	//�N���X�^�Ɨ��q�̊֘A���̓o�^
	for(int i = 0; i < sm_particleNum; i++)
	{
		SetClusterStrctInfo(i, PtoCNum);	//�J�E���g��O�ōs���Ă��邽�߁C��������g��
	}
	cout << __FUNCTION__ << " �N���X�^�Ɨ��q�̊֘A���̓o�^�o�^" << endl;

	delete[] PtoCNum;

//�f�o�b�O
	//for(int i = 0; i < m_iClusteresNum; i++)
	//{
	//	m_ice->DebugCtoP(i);
	//}

	//for(int i = 0; i < m_pPS->GetNumParticles(); i++)
	//{
	//	m_ice->DebugPtoC(i);
	//}
}

//�M����������
void IceObject::InitHeatTransfer(float effectiveRadius, float timeStep, float tempMax, float tempMin, float latentHeat, float cffcntHt, float cffcntTd)
{//	cout << __FUNCTION__ << endl;
	m_heatTransfer = new HeatTransfar(sm_particleNum*2);		//�ŏ��ɍő吔���m�ۂ��Ă����āC�g���͍̂쐬���ꂽ�p�[�e�B�N���܂łƂ���
	m_heatTransfer->setCarnelConstant(effectiveRadius);			//�J�[�l���֐��̒萔�̂��߂̏���
	m_heatTransfer->setNumVertices(sm_particleNum);				//�p�[�e�B�N���̐����擾

	m_heatTransfer->setTimeStep(timeStep);
	m_heatTransfer->setTempMax(tempMax);
	m_heatTransfer->setTempMin(tempMin);
	m_heatTransfer->setLatentHeat(latentHeat);
	m_heatTransfer->setCffCntHt(cffcntHt);
	m_heatTransfer->setCffCntTd(cffcntTd);
}

//�^���n�₻�ꂼ��̎�@�̂��߂̏�����
void IceObject::InitMoveMethod()
{	cout << __FUNCTION__ << endl;

	//NULL�ŏ�����
	m_iceClsuterMove = NULL;
	m_iceJudeMove = NULL;
	m_iceInterPolation = NULL;
	m_iceInterPolationJudge = NULL;
	m_iceCalcMethod = NULL;

	//�S�ăm�[�}��
	ChangeMode_ClusterMove_Normal();
	ChangeMode_JudgeMove_Normal();
	ChangeMode_InterPolation_Normal();
	ChangeMode_IntrpJudge_Normal();
	ChangeMode_CalcMethod_Normal();

	m_SurfSm = NULL;

	IceObject::InitInterPolation();					//���`��Ԃ̏�����

#ifdef MK_USE_GPU
	InitGPU();										//�\���̏�񂪊���������GPU��������
#endif
}

//�������ɗp����p�X�̏�����
void IceObject::InitPath()
{
	m_SurfSm = new Surf_SM(s_sphPrtPos, s_sphPrtVel, m_iceSM, m_iceStrct, sm_clusterNum, sm_particleNum);
}

//���[�h�ؑ�
//�v�Z��@
void IceObject::ChangeMode_CalcMethod_Normal()
{	cout << __FUNCTION__ << endl;

	if(m_iceCalcMethod != NULL) delete m_iceCalcMethod;
	m_iceCalcMethod = new Ice_CalcMethod_Normal(m_iceClsuterMove);
}

void IceObject::ChangeMode_CalcMethod_Iteration()
{	cout << __FUNCTION__ << endl;

	if(m_iceCalcMethod != NULL) delete m_iceCalcMethod;
	m_iceCalcMethod = new Ice_CalcMethod_Iteration(m_iceSM, m_iceClsuterMove, m_iceInterPolation);
}


//�^���v�Z�N���X�^�I��
void IceObject::ChangeMode_JudgeMove_Normal()
{	cout << __FUNCTION__ << endl;

	ResetSelectCluster();		//�I��I�^���v�Z�̏�����

	if(m_iceJudeMove != NULL) delete m_iceJudeMove;
	m_iceJudeMove = new Ice_JudgeMove_Normal(m_iceSM);

	m_iceClsuterMove->SetJudgeMove(m_iceJudeMove);
}

void IceObject::ChangeMode_JudgeMove_Spears()
{	cout << __FUNCTION__ << endl;

	UpdateSelectCluster();		//�I��I�^���v�Z�̍X�V

	if(m_iceJudeMove != NULL) delete m_iceJudeMove;
	m_iceJudeMove = new Ice_JudgeMove_Spears(m_iceSM, m_iceStrct);

	m_iceClsuterMove->SetJudgeMove(m_iceJudeMove);
}

//�^���v�Z
void IceObject::ChangeMode_ClusterMove_Normal()
{	cout << __FUNCTION__ << endl;

	if(m_iceClsuterMove != NULL) delete m_iceClsuterMove;
	m_iceClsuterMove = new Ice_ClusterMove_Normal(m_iceSM);

	m_iceClsuterMove->SetJudgeMove(m_iceJudeMove);

	if(m_iceCalcMethod != NULL)	m_iceCalcMethod->SetObjMove(m_iceClsuterMove);
}

void IceObject::ChangeMode_ClusterMove_Path()
{	cout << __FUNCTION__ << endl;

	if(m_SurfSm == NULL){	InitPath();	}	//�����ȉ^���v�Z�̂��߂̃p�X������

	if(m_iceClsuterMove != NULL) delete m_iceClsuterMove;
	m_iceClsuterMove = new Ice_ClusterMove_FastPath(m_iceSM, m_SurfSm);

	m_iceClsuterMove->SetJudgeMove(m_iceJudeMove);

	if(m_iceCalcMethod != NULL)	m_iceCalcMethod->SetObjMove(m_iceClsuterMove);
}

//�ŏI���ʕ�ԃN���X�^�I��
void IceObject::ChangeMode_IntrpJudge_Normal()
{	cout << __FUNCTION__ << endl;

	if(m_iceInterPolationJudge != NULL) delete m_iceInterPolationJudge;
	m_iceInterPolationJudge = new Ice_InterPolationJudge_Normal(m_iceSM);

	m_iceInterPolation->SetIntrpJudge(m_iceInterPolationJudge);
}

void IceObject::ChangeMode_IntrpJudge_Spears()
{	cout << __FUNCTION__ << endl;

	if(m_iceInterPolationJudge != NULL) delete m_iceInterPolationJudge;
	m_iceInterPolationJudge = new Ice_InterPolationJudge_Spears(m_iceSM, m_iceStrct);

	m_iceInterPolation->SetIntrpJudge(m_iceInterPolationJudge);
}

//�ŏI���ʕ��
void IceObject::ChangeMode_InterPolation_Normal()
{	cout << __FUNCTION__ << endl;

	if(m_iceInterPolation != NULL) delete m_iceInterPolation;
	m_iceInterPolation = new Ice_InterPolation_Normal(m_iceSM, m_iceStrct);

	m_iceInterPolation->SetIntrpJudge(m_iceInterPolationJudge);
}

void IceObject::ChangeMode_InterPolation_Weight()
{	cout << __FUNCTION__ << endl;

	if(m_iceInterPolation != NULL) delete m_iceInterPolation;
	m_iceInterPolation = new Ice_InterPolation_Weight(m_iceSM, m_iceStrct);

	m_iceInterPolation->SetIntrpJudge(m_iceInterPolationJudge);
}

//�N���X�^�̉^���v�Z���̓o�^
void IceObject::SetClusterMoveInfo(int pIndx)
{
	//���闱�q���܂܂��l�ʑ̂��ߖT�l�ʑ̂Ƃ��C�e�ߖT�l�ʑ̂Ɋ܂܂�闱�q�ŃN���X�^���쐬
	for(int i = 0; i < GetPtoTIndx(pIndx); i++)
	{
		int itIndx = GetPtoT(pIndx, i, 0);
		int ioIndx = GetPtoT(pIndx, i, 1);

		if(itIndx == -1 || ioIndx == -1){ continue;	}

		//�l�ʑ̂̑S�Ă̗��q���N���X�^�ɓo�^
		for(int j = 0; j < GetTtoPIndx(itIndx); j++)
		{
			int jpIndx = GetTtoP(itIndx, j);

			if(jpIndx == -1){	continue;	}
			if(m_iceSM[pIndx]->CheckIndx(jpIndx)){	continue;	}

			int pNum = m_iceSM[pIndx]->GetNumVertices();
			float mass = 1.0f;

			m_iceSM[pIndx]->AddVertex( Vec3(s_sphPrtPos[jpIndx*4+0], s_sphPrtPos[jpIndx*4+1], s_sphPrtPos[jpIndx*4+2] ), mass, jpIndx);
			m_iceSM[pIndx]->SetAlphas(pNum, 1.0);
			m_iceSM[pIndx]->SetBetas (pNum, 0.0);
			m_iceSM[pIndx]->SetLayer (pNum, 0);
		}
	}

	if(sm_layerNum == 0){	return;	}

	//�ߖT�l�ʑ̂�layer���ǂ�C���q��ǉ����Ă���
	//TODO::�s����ɂȂ�Ȃ�Clayer�������ق�beta��������
	for(int i = 0; i < GetPtoTIndx(pIndx); i++)
	{
		int itIndx = GetPtoT(pIndx, i, 0);
		int ioIndx = GetPtoT(pIndx, i, 1);

		if(itIndx == -1 || ioIndx == -1){	continue;	}

		for(int j = 0; j < GetNTNum(itIndx); j++)
		{
			int jtIndx = GetNeighborTetra(itIndx, j, 0);
			int jlIndx = GetNeighborTetra(itIndx, j, 1);

			if(jtIndx == -1 || jlIndx == -1){	continue;	}

			//�l�ʑ̂̑S�Ă̗��q���N���X�^�ɓo�^
			for(int k = 0; k < GetTtoPIndx(jtIndx); k++)
			{
				int kpIndx = GetTtoP(jtIndx, k);

				if(kpIndx == -1){	continue;	}
				if(m_iceSM[pIndx]->CheckIndx(kpIndx)){	/*cout << "contineue kpIndx = " << kpIndx << endl;*/ continue;	}

				int pNum = m_iceSM[pIndx]->GetNumVertices();
				float mass = 1.0f;

				m_iceSM[pIndx]->AddVertex( Vec3(s_sphPrtPos[kpIndx*4+0], s_sphPrtPos[kpIndx*4+1], s_sphPrtPos[kpIndx*4+2] ), mass, kpIndx);
				m_iceSM[pIndx]->SetAlphas(pNum, 1.0);
				m_iceSM[pIndx]->SetBetas (pNum, 0.0);
				m_iceSM[pIndx]->SetLayer (pNum, jlIndx);
				//cout << "pIndx = " << pIndx << " GetNumVertices = " << m_iceSM[pIndx]->GetNumVertices() << endl;
			}
		}
	}
}

//���q���N���X�^��N���X�^�����q��o�^
void IceObject::SetClusterStrctInfo(int cIndx, int* PtoCNum)
{
	//���q�������Ă���l�ʑ̂̔ԍ���o�^���邽�߂̏���
	//pCountList�ɂ́CcIndx�Ԗڂ̃N���X�^�Ɋ܂܂��e���q���C���ꂼ�ꂢ���̃N���X�^�ɑ����邩�����߂ĕۑ�����
	int* pCountList = new int[m_iceSM[cIndx]->GetNumVertices()];
	int* pLayerList = new int[m_iceSM[cIndx]->GetNumVertices()];			//���q�̏������C���[
	
	for(int i = 0; i < m_iceSM[cIndx]->GetNumVertices(); i++)
	{
		int pIndx = m_iceSM[cIndx]->GetParticleIndx(i);

		//TODO::��������z�肵�Ă��Ȃ����߁C����ď㏑�����Ă��܂��Ă���
		//TODO::-1��T�����ď㏑������̂ɐ؂�ւ���D
		for(int j = 0; j < GetPtoCMax(); j++)
		{
			//�z���ŋ󂢂Ă���X�y�[�X��T���i-1�łȂ��ꏊ��T���j
			if(GetPtoC(pIndx, j, 0) != -1 || GetPtoC(pIndx, j, 1) != -1){	continue;	}
			
			if(j >= GetPtoCIndx(pIndx))
			{
				SetPtoCIndx(pIndx, j+1);		//���݂̖���Indx���傫���Ȃ�X�V
			}

			pCountList[i] = j;					//�z���ł̏ꏊ��ۑ�
			break;
		}

		pLayerList[i] = m_iceSM[cIndx]->GetLayer(i);					//���q�����w�ڂ̋ߖT�Ȃ̂����擾
	}

	//���q�ƃN���X�^�̏��o�^
	vector<int> pIndxList;

	for(int i = 0; i < GetCtoPNum(cIndx); i++)
	{
		int pIndx = m_iceSM[cIndx]->GetParticleIndx(i);
			
		SetPtoC(pIndx, pCountList[i], cIndx, i, pLayerList[i]);		//���q���������Ă���N���X�^��o�^

		pIndxList.push_back(pIndx);									//�N���X�^�Ɋ܂܂�闱�q�W��
	}

	SetCtoP(cIndx, pIndxList, pLayerList);
	
	delete[] pCountList;
	delete[] pLayerList;
}

//��ԏ����̂��߂̃p�����[�^�̏�����
void IceObject::InitInterPolation()
{	cout << __FUNCTION__ << endl;
	m_fInterPolationCoefficience = new float[sm_particleNum];	//���`��ԌW��

	for(int i = 0; i < sm_particleNum; ++i)
	{
		m_fInterPolationCoefficience[i] = 1.0f;
	}
}

//GPU�����ŗp����ϐ��̏�����
void IceObject::InitGPU()
{
	//�ő̍\���̏�����
	m_iceStrct->InitGPU();

	//�N���X�^�Ɋւ���GPU�̏�����
	//���ׂẴN���X�^�ɃA�N�Z�X�ł���悤�Ƀ|�C���^��n���Ă���
	Ice_SM::InitGPU(m_iceSM, sd_sphPrtPos, sd_sphPrtVel, sm_particleNum, sm_maxParticleNum);

#if defined(USE_PATH)
{
	m_SurfSm->InitPathGPU();		//�������ɗp����p�X�̏�����
}
#endif
	InitIceObjGPU();
}

void IceObject::InitIceObjGPU()
{
	int MAXCLUSTER = GetMaxClusterNum();

	//�e���q�̍ŏI�I�Ȉʒu�E���x�f�[�^
	cudaMalloc((void**)&sd_sldPrtPos,	sizeof(float) * MAXCLUSTER * SM_DIM);
	cudaMalloc((void**)&sd_sldPrtVel,	sizeof(float) * MAXCLUSTER * SM_DIM);

	//�ŏI�ʒu�E���x�����݂̃f�[�^�ŏ�����
	float* fPoses = new float[MAXCLUSTER * SM_DIM];
	float* fVeles = new float[MAXCLUSTER * SM_DIM];

	//s_pfPrtPos�Ȃǂ̓f�[�^�̒��g��DIM=4�ō���Ă���̂ŁC�������Ȃ��Ƃ����Ȃ�
	//TODO::���q�T�C�Y���傫���Ȃ�ƁC���������m�ۂł��Ȃ���������Ȃ��̂ɒ���
	int sphDIM = 4;
	for(int i = 0; i < MAXCLUSTER; ++i)
	{
		for(int j = 0; j < SM_DIM; ++j)
		{
			fPoses[i*SM_DIM+j] = s_sphPrtPos[i*sphDIM+j];
			fVeles[i*SM_DIM+j] = s_sphPrtVel[i*sphDIM+j];
		}
	}

	cudaMemcpy(sd_sldPrtPos, fPoses, sizeof(float) * MAXCLUSTER * SM_DIM, cudaMemcpyHostToDevice);
	cudaMemcpy(sd_sldPrtVel, fVeles, sizeof(float) * MAXCLUSTER * SM_DIM, cudaMemcpyHostToDevice);

	////�������̓]�������܂����������̊m�F
	////�ꎞ�z��̃��Z�b�g
	//for(int i = 0; i < MAXCLUSTER; ++i)
	//{
	//	for(int j = 0; j < SM_DIM; ++j)
	//	{
	//		fPoses[i*SM_DIM+j] = 0.0f;
	//		fVeles[i*SM_DIM+j] = 0.0f;
	//	}
	//}

	////�f�[�^��]��
	//cudaMemcpy(fPoses, d_FinalPos, sizeof(float) * MAXCLUSTER * SM_DIM, cudaMemcpyDeviceToHost);
	//cudaMemcpy(fVeles, d_FinalVel, sizeof(float) * MAXCLUSTER * SM_DIM, cudaMemcpyDeviceToHost);

	////�z�X�g���̃f�[�^��]���������ʂ��_���v
	//ofstream ofs( "DtoH_Test.txt" );
	//ofs << "DtoH_Test" << endl;
	//
	////�f�o�C�X���̃f�[�^��]��
	//for(int i = 0; i < MAXCLUSTER; i++)
	//{
	//	ofs << "particle" << i << " pos::(" << fPoses[i*SM_DIM+0] << ", " << fPoses[i*SM_DIM+1] << ", " << fPoses[i*SM_DIM+2] << ")" << endl;
	//	ofs << "sph     " << i << " pos::(" << s_pfPrtPos[i*sphDIM+0] << ", " << s_pfPrtPos[i*sphDIM+1] << ", " << s_pfPrtPos[i*sphDIM+2] << ")" << endl;
	//	ofs << "particle" << i << " vel::(" << fVeles[i*SM_DIM+0] << ", " << fVeles[i*SM_DIM+1] << ", " << fVeles[i*SM_DIM+2] << ")" << endl;
	//	ofs << "sph     " << i << " vel::(" << s_pfPrtVel[i*sphDIM+0] << ", " << s_pfPrtVel[i*sphDIM+1] << ", " << s_pfPrtVel[i*sphDIM+2] << ")" << endl;
	//}

	delete[] fPoses;
	delete[] fVeles;
}

//�ő̂̉^���v�Z�C�v�Z���ʂ̓���
//�������e�͕ʃN���X�ɋL�q
void IceObject::StepObjMoveCPU()
{
#ifdef USE_ITR	//�p�X

	#ifdef USE_PATH
		StepObjMoveIterationUsePath();
	#else
		#ifdef USE_SELECTED
			StepObjMoveIterationSelected();
		#else
			StepObjMoveIteration();
			//m_iceObj->StepObjMoveIterationWeighted();	//�d�ݕt���{����
		#endif
	#endif

#else

	m_iceCalcMethod->StepObjMove();				RXTIMER("StepObjMove");

#endif

	//�t�̂ƌő̂̉^���̐��`���
	m_iceInterPolation->StepInterPolation();	RXTIMER("StepInterPolation");

//�f�o�b�O
	//DebugNowMoveMethod();
}

//�ő̂̉^���v�Z
void IceObject::StepObjMoveNormal()
{
	#pragma omp parallel
	{
	#pragma omp for
		for(int i = 0; i < sm_clusterNum; ++i)
		{	
			if(m_iceSM[i]->GetNumVertices() == 0){		continue;	}
			m_iceSM[i]->UpdateCPU();				//�^���v�Z
		}
	}//#pragma omp parallel
}

void IceObject::StepObjMoveSelected()
{
	#pragma omp parallel
	{
	#pragma omp for
		for(int i = 0; i < sm_clusterNum; ++i)
		{
			if(m_iceStrct->GetMotionCalcCluster(i) == 0){	continue;	}
			if(m_iceSM[i]->GetNumVertices() == 0){		continue;	}

			m_iceSM[i]->UpdateCPU();
		}
	}//#pragma omp parallel
}

//GPU
void IceObject::StepObjMoveGPU()
{
	#ifdef USE_ITR
	//TODO: ������
	//StepObjCalcWidhIteration();	RXTIMER("StepObjCalcWidhIteration");
#else

	#ifdef USE_PATH
		StepObjMoveGPUUsePath();	RXTIMER("StepObjMoveGPUUsePath");
	#else
		StepObjMoveGPUNormal();		RXTIMER("StepObjMoveGPU");			//�N���X�^�̉^���v�Z
	#endif
		
	StepInterPolationNormal();		RXTIMER("StepInterpolationGPU");	//���a�v�Z�C���`���

#endif
}

void IceObject::StepObjMoveGPUNormal()
{
	Ice_SM::SetDevicePosPointer(sd_sphPrtPos);	//VBO���g���Ă���̂ŁC���ꂪ�Ȃ��ƃG���[���o��̂ɒ���
	Ice_SM::UpdateGPU();
}

//�p�X��p�����������@�ő̂̉^���v�Z
void IceObject::StepObjMoveUsePath()
{
	//prefixSum�̍X�V
	m_SurfSm->UpdatePrefixSum();

	//�N���X�^�̃p�����[�^�X�V
	#pragma omp parallel for
	for(int i = 0; i < sm_clusterNum; ++i)
	{	
		if(GetPtoCNum(i) == 0){	continue;	}
		m_iceSM[i]->SetNowCm(m_SurfSm->CalcCmSum(i));	//�d�S�̍X�V
		m_iceSM[i]->SetApq(m_SurfSm->CalcApqSum(i));		//�ό`�s��̍X�V
	}

	//�^���v�Z
	#pragma omp parallel for
	for(int i = 0; i < sm_clusterNum; ++i)
	{	
		if(GetPtoCNum(i) == 0){	continue; }
		m_iceSM[i]->UpdateUsePathCPU();
	}
}

void IceObject::StepObjMoveGPUUsePath()
{
	//prefixSum�̍X�V
	m_SurfSm->SetDevicePointer(sd_sphPrtPos, sd_sphPrtVel);	//VBO���g���Ă���̂ŁC���ꂪ�Ȃ��ƃG���[���o��̂ɒ���
	m_SurfSm->UpdatePrefixSumGPU();
	RXTIMER("UpdatePrefixSumGPU");

	//�N���X�^�̃p�����[�^�X�V
	LauchUpdateSMFromPath
	(
		sm_clusterNum,
		sd_sphPrtPos,
		sd_sphPrtVel,
		//SM-------------------------------------------------
		Ice_SM::GetOrgPosPointer(),
		Ice_SM::GetDevicePosPointer(),
		Ice_SM::GetOrgCmPointer(),
		Ice_SM::GetCurCmPointer(),
		Ice_SM::GetDeviceApqPointer(),
		Ice_SM::GetDeviceVelPointer(),
		Ice_SM::GetDeviceIndexesPointer(),
		Ice_SM::GetDeviceIndexSetPointer(),
		//Path-----------------------------------------------
		m_SurfSm->GetDevicePRTtoPTHPointer(),
		m_SurfSm->GetDevicePTHandPrfxSetPointer(),
		m_SurfSm->GetDecvicePrfxPos(),
		m_SurfSm->GetDecvicePrfxApq(),
		//IceStruct------------------------------------------
		m_iceStrct->GetDeviceCtoPointer(),
		m_iceStrct->GetDeviceCtoPNumPointer(),
		m_iceStrct->GetCtoPMax(),
		2,
		0.01		//TODO: 0.02�ł́H���낢��΂�΂���ۂ�
	);

	RXTIMER("LauchUpdateSMFromPath");

	//�^���v�Z
	Ice_SM::SetDevicePosPointer(sd_sphPrtPos);				//VBO���g���Ă���̂ŁC���ꂪ�Ȃ��ƃG���[���o��̂ɒ���
	Ice_SM::UpdateUsePathGPU();
	RXTIMER("UpdateUsePathGPU");

//�f�o�b�O
	//DebugObjMoveUsePathWithGPU();
}

//������p�����ő̂̉^���v�Z
void IceObject::StepObjMoveIteration()
{
	#pragma omp parallel for
	for(int i = 0; i < sm_clusterNum; ++i)
	{	
		if(m_iceSM[i]->GetNumVertices() == 0){	continue;	}
		//GetMoveObj(i)->calExternalForces();				//���x�𔽉f���Ĉʒu�X�V
		//GetMoveObj(i)->ShapeMatchingSolid();			//���x�𔽉f������Ԃ�SM�@
		GetMoveObj(i)->UpdateCPU();
	}

	StepInterPolationForCluster();						//�e�N���X�^�̑��a�v�Z�C���`��Ԃ����Čő̂̍ŏI�ʒu������

	//�^���v�Z�𔽕�
	for(int itr = 1; itr < Ice_SM::GetIteration(); ++itr)
	{
		#pragma omp parallel for
		for(int i = 0; i < sm_clusterNum; ++i)
		{	
			if(m_iceSM[i]->GetNumVertices() == 0){	continue;	}
			GetMoveObj(i)->ShapeMatchingIteration();		//���݂̗��q�ʒu��p����SM�@
		}	

		StepInterPolationForCluster();							//�e�N���X�^�̑��a�v�Z�C���`��Ԃ����Čő̂̍ŏI�ʒu������
	}

	//���x�Z�o
	#pragma omp parallel for
	for(int i = 0; i < sm_clusterNum; ++i)
	{	
		if(m_iceSM[i]->GetNumVertices() == 0){	continue;	}
		GetMoveObj(i)->integrateIteration();
	}
}

//�p�X��p���������^���v�Z
void IceObject::StepObjMoveIterationUsePath()
{
	m_SurfSm->UpdatePrefixSum();								//prefixSum�̍X�V

	//�e�N���X�^�̃f�[�^�X�V
	#pragma omp parallel for
	for(int i = 0; i < sm_clusterNum; ++i)
	{	
		if(GetPtoCNum(i) == 0){	continue;	}
		GetMoveObj(i)->SetNowCm(m_SurfSm->CalcCmSum(i));	//�d�S�̍X�V
		GetMoveObj(i)->SetApq(m_SurfSm->CalcApqSum(i));		//�ό`�s��̍X�V
		GetMoveObj(i)->ShapeMatchingUsePath();				//���݂̈ʒu��SM�@
	}

	StepInterPolationForCluster();							//�e�N���X�^�̑��a�v�Z�C���`��Ԃ����Čő̂̍ŏI�ʒu������

	//�^���v�Z�𔽕�
	for(int itr = 1; itr < Ice_SM::GetIteration(); ++itr)
	{
		//����
		{
			m_SurfSm->UpdatePrefixSumItr();							//prefixSum�̍X�V
		}

		//�e�N���X�^�̃f�[�^�X�V
		#pragma omp parallel for
		for(int i = 0; i < sm_clusterNum; ++i)
		{	
			if(GetPtoCNum(i) == 0){	continue;	}
			GetMoveObj(i)->SetNowCm(m_SurfSm->CalcCmSum(i));		//�d�S�̍X�V
			GetMoveObj(i)->SetApq(m_SurfSm->CalcApqSum(i));		//�ό`�s��̍X�V
			GetMoveObj(i)->ShapeMatchingUsePath();				//���݂̈ʒu��SM�@
		}

		StepInterPolationForCluster();							//�e�N���X�^�̑��a�v�Z�C���`��Ԃ����Čő̂̍ŏI�ʒu������
	}

	//���x�Z�o
	#pragma omp parallel for
	for(int i = 0; i < sm_clusterNum; ++i)
	{	
		if(GetPtoCNum(i) == 0){	continue;	}
		GetMoveObj(i)->integrateIteration();
	}
}

//�ő̂̉^���v�Z�C���a�v�Z�C��ԏ����@GPU�����@������������
void IceObject::StepObjCalcWidhIteration()
{
	for(int itr = 0; itr < g_iterationNum; ++itr)
	{
		if(itr == 0)
		{
			Ice_SM::UpdateGPU();
		}
		else
		{
			Ice_SM::UpdateIterationGPU(sd_sldPrtPos, sd_sldPrtVel);
		}

		StepInterPolationForCluster();			//���a�v�Z�C���`���
	}

	StepInterPolationNormal();
}

//�ő̂̉^���v�Z�@�d�ݕt���{����
void IceObject::StepObjMoveIterationWeighted()
{
	//�^���v�Z�𔽕�
	for(int itr = 0; itr < Ice_SM::GetIteration(); ++itr)
	{
		//����
		if(itr == 0)
		{
			#pragma omp parallel for
			for(int i = 0; i < sm_clusterNum; ++i)
			{	
				if(GetPtoCNum(i) == 0){	continue;	}
				GetMoveObj(i)->calExternalForces();				//���x�𔽉f���Ĉʒu�X�V
				GetMoveObj(i)->ShapeMatchingSolid();			//���x�𔽉f������Ԃ�SM�@
			}
		}
		//����
		else
		{
			#pragma omp parallel for
			for(int i = 0; i < sm_clusterNum; ++i)
			{	
				if(GetPtoCNum(i) == 0){	continue;	}
				GetMoveObj(i)->ShapeMatchingIteration();		//���݂̗��q�ʒu��p����SM�@
			}
		}

		StepInterPolationForClusterWeighted();					//�e�N���X�^�̑��a�v�Z�C���`��Ԃ����Čő̂̍ŏI�ʒu������
	}

	//���x�Z�o
	#pragma omp parallel for
	for(int i = 0; i < sm_clusterNum; ++i)
	{	
		if(GetPtoCNum(i) == 0){	continue;	}
		GetMoveObj(i)->integrateIteration();
	}
}

//�X�p�[�X�{����
void IceObject::StepObjMoveIterationSelected()
{
	//�^���v�Z�𔽕�
	for(int itr = 0; itr < Ice_SM::GetIteration(); ++itr)
	{
		//����
		if(itr == 0)
		{
			#pragma omp parallel for
			for(int i = 0; i < sm_clusterNum; ++i)
			{	
				if(GetPtoCNum(i) == 0){	continue;	}
				if(m_iceStrct->GetMotionCalcCluster(i) == 0){	continue;	}

				GetMoveObj(i)->calExternalForces();				//���x�𔽉f���Ĉʒu�X�V
				GetMoveObj(i)->ShapeMatchingSolid();			//���x�𔽉f������Ԃ�SM�@
			}
		}
		//����
		else
		{
			#pragma omp parallel for
			for(int i = 0; i < sm_clusterNum; ++i)
			{	
				if(GetPtoCNum(i) == 0){	continue;	}

				if(m_iceStrct->GetMotionCalcCluster(i) == 0){	continue;	}
				GetMoveObj(i)->ShapeMatchingIteration();		//���݂̗��q�ʒu��p����SM�@
			}
		}

		StepInterPolationSelectedForCluster();					//�e�N���X�^�̑��a�v�Z�C���`��Ԃ����Čő̂̍ŏI�ʒu������
	}

	//���x�Z�o
	#pragma omp parallel for
	for(int i = 0; i < sm_clusterNum; ++i)
	{	
		if(GetPtoCNum(i) == 0){	continue;	}
		if(m_iceStrct->GetMotionCalcCluster(i) == 0){	continue;	}

		GetMoveObj(i)->integrateIteration();
	}
}

void IceObject::StepInterPolationNormal()
{
#if defined(MK_USE_GPU)

	//�ő̉^���̍ŏI�ʒu�v�Z
	float* smPrtPos = Ice_SM::GetDevicePosPointer();
	float* smPrtVel = Ice_SM::GetDeviceVelPointer();
	int* indxSet = Ice_SM::GetDeviceIndexSetPointer();

	sd_sphPrtPos = Ice_SM::GetDeviceSPHPosPointer();
	sd_sphPrtVel = Ice_SM::GetDeviceSPHVelPointer();

	int* PtoCIndx = IceStructure::GetDevicePtoCIndxPointer();
	int* PtoC = IceStructure::GetDevicePtoCPointer();
	int PNumMax = m_iceStrct->GetPNumMax();
	int PtoCMax = m_iceStrct->GetPtoCMax();
	int PtoCParamSize = 3;

	//�ő̉^���̍ŏI�ʒu�v�Z
	LaunchCalcAverageGPU(sm_particleNum, sd_sldPrtPos, sd_sldPrtVel, sd_sphPrtPos, sd_sphPrtVel, smPrtPos, smPrtVel, indxSet, PtoCIndx, PtoC, PNumMax, PtoCMax, PtoCParamSize);

	//�t�̂ƌő̂̕��
	LaunchInterPolationGPU(sm_particleNum, sd_sldPrtPos, sd_sldPrtVel, sd_sphPrtPos, sd_sphPrtVel);

#else

	Vec3 pos,vel;

	#pragma omp parallel for private(pos, vel)
	for(int i = 0; i < sm_particleNum; ++i)
	{
		//if(GetParticleNum() <= i){	continue;	}	//�Z���݂̂̎����̂Ƃ��ɕK�v�ɂȂ�D
		//if(GetPtoCNum(i) <= 0){		continue;	}
	
		//�ő̉^���̍ŏI�ʒu�v�Z
		CalcAverageCPU(i, pos, vel);

		//�t�̂ƌő̂̕��
		LinerInterPolationCPU(i, pos, vel);			//���`���
	}

//�f�o�b�O
	//DebugDeformationAmount();
	//DebugDeformationAverage();

#endif
}

//���q��I��I�ɉ^���v�Z������p�^�[���̃e�X�g
void IceObject::StepInterPolationSelected()
{
	//Vec3 pos,vel;

	//vector<Vec3> preSphPos;
	//preSphPos.resize(sm_particleNum);
	//for(int i = 0; i < sm_particleNum; i++)
	//{
	//	Vec3 pos;
	//	pos[0] = s_sphPrtPos[i*4+0];
	//	pos[1] = s_sphPrtPos[i*4+1];
	//	pos[2] = s_sphPrtPos[i*4+2];

	//	preSphPos.push_back(pos);
	//}

	//#pragma omp parallel for private(pos, vel)
	//for(int i = 0; i < sm_particleNum; ++i)
	//{
	//	if(GetParticleNum() <= i){	continue;	}	//�Z���݂̂̎����̂Ƃ��ɕK�v�ɂȂ�D
	//	if(GetPtoCNum(i) <= 0){		continue;	}
	//
	//	//�ő̉^���̍ŏI�ʒu�v�Z
	//	CalcAverageSelected(i, pos, vel);

	//	//�t�̂ƌő̂̕��
	//	LinerInterPolationCPU(i, pos, vel);			//���`���
	//}

	vector<unsigned> addParticleNum(sm_particleNum, 0);
	vector<float> deformationSum(sm_particleNum, 0.0f);

	float* sldPos = Ice_SM::GetSldPosPointer();
	float* sldVel = Ice_SM::GetSldVelPointer();

	//sld�̏�����
	for(int cIndx = 0; cIndx < sm_particleNum; cIndx++)
	{
		for(int dim = 0; dim < SM_DIM; dim++)
		{
			sldPos[cIndx*SM_DIM+dim] = 0.0f;
			sldVel[cIndx*SM_DIM+dim] = 0.0f;
		}
	}

	for(int cIndx = 0; cIndx < sm_particleNum; cIndx++)
	{
		//��ԗ��q�N���X�^�Ɋ܂܂��f�[�^���g�킸�C�v�Z���q�N���X�^�Ɋ܂܂��f�[�^�݂̂��g��
		if(m_iceStrct->GetMotionCalcCluster(cIndx) == 0){	continue;	}

		//�N���X�^���܂�ł���e���q�̈ʒu�E���x���擾
		//�e�ʒu�E���x��static�ȍŏI�ʒu�ɑ���
		for(int oIndx = 0; oIndx < GetMoveObj(cIndx)->GetIndxNum(); oIndx++)
		{
			int pIndx = GetMoveObj(cIndx)->GetParticleIndx(oIndx);

			if(pIndx == MAXINT){	continue;	}

			float defAmount = GetMoveObj(cIndx)->GetDefAmount();
			//defAmount = pow(defAmount, 2.0f);

			Vec3 pos = GetMoveObj(cIndx)->GetVertexPos(oIndx) * defAmount;
			Vec3 vel = GetMoveObj(cIndx)->GetVertexVel(oIndx) * defAmount;

			for(int dim = 0; dim < SM_DIM; dim++)
			{
				sldPos[pIndx*SM_DIM+dim] += pos[dim];
				sldVel[pIndx*SM_DIM+dim] += vel[dim];
			}

			//�����������J�E���g
			//addParticleNum[pIndx] += 1;

			//���ό`�ʂ��J�E���g
			deformationSum[pIndx] += defAmount;
		}
	}

	//���������ŕ��ς��C���q�ʒu�ɔ��f
	for(int i = 0; i < sm_particleNum; i++)
	{
		int sphIndx = i*4;
		int smIndx = i*SM_DIM;

		//float clusterNum = addParticleNum[i];
		float clusterNum = deformationSum[i];

		//if(GetPtoCNum(i) == 0){	continue;	}
		if(clusterNum == 0){	continue;	}

		////�f�o�b�O
		//if(clusterNum == 0)
		//{	
		//	for(int dim = 0; dim < SM_DIM; dim++)
		//	{
		//		s_sphPrtPos[sphIndx+dim] = 0.0f;
		//		s_sphPrtVel[sphIndx+dim] = 0.0f;
		//	}
		//	continue;	
		//}

		//���q�ʒu�ɔ��f
		for(int dim = 0; dim < SM_DIM; dim++)
		{
			s_sphPrtPos[sphIndx+dim] = sldPos[smIndx+dim]/clusterNum;
			s_sphPrtVel[sphIndx+dim] = sldVel[smIndx+dim]/clusterNum;
		}
	}

	//#pragma omp parallel for private(pos, vel)
	//for(int i = 0; i < sm_particleNum; ++i)
	//{
	//	if(GetParticleNum() <= i){	continue;	}	//�Z���݂̂̎����̂Ƃ��ɕK�v�ɂȂ�D
	//	if(GetPtoCNum(i) <= 0){		continue;	}

	//	//�ό`�ʂ̔��f
	//	UpdateUnSelectedCluster(i, pos, vel, preSphPos);
	//}


//�f�o�b�O
	//DebugDeformationAmount();
	//DebugDeformationAverage();
}

void IceObject::UpdateUnSelectedCluster(int cIndx, const Vec3& pos, const Vec3& vel, const vector<Vec3>& preSphPos)
{
	if(m_iceStrct->GetMotionCalcCluster(cIndx) == 0){	return;	}

	float defAmount = 0.0f;

	for(int i = 0; i < m_iceSM[cIndx]->GetVertexNum(); ++i)
	{
		int pIndx = m_iceSM[cIndx]->GetParticleIndx(i);
		
		if(m_iceStrct->GetMotionCalcCluster(pIndx) == 0){	continue;	}

		Vec3 nowPos = m_iceSM[cIndx]->GetVertexPos(i);
		Vec3 prePos = preSphPos[pIndx];

		defAmount += abs(nowPos[0]-prePos[0]);
		defAmount += abs(nowPos[1]-prePos[1]);
		defAmount += abs(nowPos[2]-prePos[2]);
	}

	m_iceSM[cIndx]->SetDefAmount(defAmount);
}

void IceObject::StepWeightedInterPolation()
{
	Vec3 pos,vel;

	#pragma omp parallel for private(pos, vel)
	for(int i = 0; i < sm_particleNum; ++i)
	{
		if(GetParticleNum() <= i){	continue;	}	//�Z���݂̂̎����̂Ƃ��ɕK�v�ɂȂ�D
		if(GetPtoCNum(i) <= 0){		continue;	}
	
		//�ő̉^���̍ŏI�ʒu�v�Z
		CalcWeightedVector(i, pos, vel);

		//�t�̂ƌő̂̕��
		LinerInterPolationCPU(i, pos, vel);			//���`���
	}

//�f�o�b�O
	//DebugDeformationAmount();
}

void IceObject::StepInterPolationForClusterWeighted()
{
	Vec3 pos,vel;

	#pragma omp parallel for private(pos, vel)
	for(int i = 0; i < sm_particleNum; ++i)
	{
		if(GetParticleNum() <= i){	continue;	}		//�Z���݂̂̎����̂Ƃ��ɕK�v�ɂȂ�D
		if(GetPtoCNum(i) <= 0){		continue;	}

		//�ő̉^���̍ŏI�ʒu�v�Z
		//CalcAverageCPU(i, pos, vel);
		CalcWeightedVector(i, pos, vel);

		//�t�̂ƌő̂̕��
		LinerInterPolationForClusterCPU(i, pos, vel);	//���`���
	}
}

//����������p�������̂ƌő̂̍ŏI�I�ȉ^���v�Z
void IceObject::StepInterPolationForCluster()
{
#if defined(MK_USE_GPU)
	//�ő̉^���̍ŏI�ʒu�v�Z
	float* smPrtPos = Ice_SM::GetDevicePosPointer();
	float* smPrtVel = Ice_SM::GetDeviceVelPointer();
	int* indxSet = Ice_SM::GetDeviceIndexSetPointer();

	sd_sphPrtPos = Ice_SM::GetDeviceSPHPosPointer();
	sd_sphPrtVel = Ice_SM::GetDeviceSPHVelPointer();

	int* PtoCIndx = IceStructure::GetDevicePtoCIndxPointer();
	int* PtoC = IceStructure::GetDevicePtoCPointer();
	int PNumMax = m_iceStrct->GetPNumMax();
	int PtoCMax = m_iceStrct->GetPtoCMax();
	int PtoCParamSize = 3;

	//�ő̉^���̍ŏI�ʒu�v�Z
	LaunchCalcAverageGPU(sm_particleNum, sd_sldPrtPos, sd_sldPrtVel, sd_sphPrtPos, sd_sphPrtVel, smPrtPos, smPrtVel, indxSet, PtoCIndx, PtoC, PNumMax, PtoCMax, PtoCParamSize);

#else
	Vec3 pos,vel;

	#pragma omp parallel for private(pos, vel)
	for(int i = 0; i < sm_particleNum; ++i)
	{
		if(GetPtoCNum(i) <= 0){		continue;	}
		if(m_iceStrct->GetMotionCalcCluster(i) == 0){	continue;	}

		//�ő̉^���̍ŏI�ʒu�v�Z
		CalcAverageCPU(i, pos, vel);

		//�t�̂ƌő̂̕��
		LinerInterPolationForClusterCPU(i, pos, vel);	//���`���
	}

#endif
}

void IceObject::StepInterPolationSelectedForCluster()
{
	vector<unsigned> addParticleNum(sm_particleNum, 0);
	vector<float> deformationSum(sm_particleNum, 0.0f);

	float* sldPos = Ice_SM::GetSldPosPointer();
	float* sldVel = Ice_SM::GetSldVelPointer();

	//sld�̏�����
	for(int cIndx = 0; cIndx < sm_particleNum; cIndx++)
	{
		for(int dim = 0; dim < SM_DIM; dim++)
		{
			sldPos[cIndx*SM_DIM+dim] = 0.0f;
			sldVel[cIndx*SM_DIM+dim] = 0.0f;
		}
	}

	for(int cIndx = 0; cIndx < sm_particleNum; cIndx++)
	{
		//��ԗ��q�N���X�^�Ɋ܂܂��f�[�^���g�킸�C�v�Z���q�N���X�^�Ɋ܂܂��f�[�^�݂̂��g��
		if(m_iceStrct->GetMotionCalcCluster(cIndx) == 0){	continue;	}

		//�N���X�^���܂�ł���e���q�̈ʒu�E���x���擾
		//�e�ʒu�E���x��static�ȍŏI�ʒu�ɑ���
		for(int oIndx = 0; oIndx < GetMoveObj(cIndx)->GetIndxNum(); oIndx++)
		{
			int pIndx = GetMoveObj(cIndx)->GetParticleIndx(oIndx);

			if(pIndx == MAXINT){	continue;	}

			//float defAmount = 1.0f;
			float defAmount = GetMoveObj(cIndx)->GetDefAmount();
			//defAmount = pow(defAmount, 2.0f);

			Vec3 pos = GetMoveObj(cIndx)->GetVertexPos(oIndx) * defAmount;
			Vec3 vel = GetMoveObj(cIndx)->GetVertexVel(oIndx) * defAmount;

			for(int dim = 0; dim < SM_DIM; dim++)
			{
				sldPos[pIndx*SM_DIM+dim] += pos[dim];
				sldVel[pIndx*SM_DIM+dim] += vel[dim];
			}

			////�����������J�E���g
			//addParticleNum[pIndx] += 1;

			//���ό`�ʂ��J�E���g
			deformationSum[pIndx] += defAmount;
		}
	}

	//���������ŕ��ς��C���q�ʒu�ɔ��f
	for(int i = 0; i < sm_particleNum; i++)
	{
		int smIndx = i*SM_DIM;

		//float clusterNum = addParticleNum[i];
		float clusterNum = deformationSum[i];

		//if(GetPtoCNum(i) == 0){	continue;	}
		if(clusterNum == 0){	continue;	}

		////�f�o�b�O
		//if(clusterNum == 0)
		//{	
		//	for(int dim = 0; dim < SM_DIM; dim++)
		//	{
		//		s_sphPrtPos[sphIndx+dim] = 0.0f;
		//		s_sphPrtVel[sphIndx+dim] = 0.0f;
		//	}
		//	continue;	
		//}

		//���q�ʒu�ɔ��f
		for(int dim = 0; dim < SM_DIM; dim++)
		{
			sldPos[smIndx+dim] = sldPos[smIndx+dim]/clusterNum;
			sldVel[smIndx+dim] = sldVel[smIndx+dim]/clusterNum;
		}
	}
}

//�e�N���X�^�̌v�Z���ʂ̕��ς��C�ő̂̍ŏI�I�ȉ^���v�Z���ʂƂ���
void IceObject::CalcAverageCPU(int pIndx, Vec3& pos, Vec3& vel)
{
	//���ꂼ��̃x�N�g�������������ς��Ƃ�
	pos = Vec3(0.0, 0.0, 0.0);
	vel = Vec3(0.0, 0.0, 0.0);
	double shapeNum = 0.0;		//�N���X�^�̐�
	int cIndx, oIndx;

	for(int i = 0; i < GetPtoCIndx(pIndx); ++i)
	{
		cIndx = GetPtoC(pIndx, i, 0);
		oIndx = GetPtoC(pIndx, i, 1);

		if(cIndx == -1 || oIndx == -1){	continue;	}

		pos += m_iceSM[cIndx]->GetVertexPos(oIndx);
		vel += m_iceSM[cIndx]->GetVertexVel(oIndx);

		shapeNum += 1.0;
	}

	//�N���X�^�̐��Ŋ���
	if(shapeNum > 0.0)
	{
		pos /= shapeNum;
		vel /= shapeNum;
	}		
	//�ǂ̃N���X�^�ɂ��܂܂�Ă��Ȃ��ꍇ�C�^����SPH�@�ɏ]��
	else
	{
		int jpIndx = pIndx*4;
		pos = Vec3(s_sphPrtPos[jpIndx+0], s_sphPrtPos[jpIndx+1], s_sphPrtPos[jpIndx+2]);
		vel = Vec3(s_sphPrtVel[jpIndx+0], s_sphPrtVel[jpIndx+1], s_sphPrtVel[jpIndx+2]);
	}
}

void IceObject::CalcAverageSelected(int pIndx, Vec3& pos, Vec3& vel)
{
	//���ꂼ��̃x�N�g�������������ς��Ƃ�
	pos = Vec3(0.0, 0.0, 0.0);
	vel = Vec3(0.0, 0.0, 0.0);
	double shapeNum = 0.0;		//�N���X�^�̐�
	int cIndx, oIndx;

	for(int i = 0; i < GetPtoCIndx(pIndx); ++i)
	{
		cIndx = GetPtoC(pIndx, i, 0);
		oIndx = GetPtoC(pIndx, i, 1);

		if(cIndx == -1 || oIndx == -1){	continue;	}

		//��ԗ��q�N���X�^�Ɋ܂܂��f�[�^���g�킸�C�v�Z���q�N���X�^�Ɋ܂܂��f�[�^�݂̂��g��
		if(m_iceStrct->GetMotionCalcCluster(cIndx) == 0){	continue;	}

		pos += m_iceSM[cIndx]->GetVertexPos(oIndx);
		vel += m_iceSM[cIndx]->GetVertexVel(oIndx);

		shapeNum += 1.0;
	}

	//�N���X�^�̐��Ŋ���
	if(shapeNum != 0.0)
	{
		pos /= shapeNum;
		vel /= shapeNum;
	}
	//�ǂ̃N���X�^�ɂ��܂܂�Ă��Ȃ��ꍇ�C�^����SPH�@�ɏ]��
	else
	{
		int jpIndx = pIndx*4;
		//pos = Vec3(s_sphPrtPos[jpIndx+0], s_sphPrtPos[jpIndx+1], s_sphPrtPos[jpIndx+2]);
		//vel = Vec3(s_sphPrtVel[jpIndx+0], s_sphPrtVel[jpIndx+1], s_sphPrtVel[jpIndx+2]);
		
		//cout << "pIndx = " << pIndx << endl;
		pos = Vec3(0.0f, 0.0f, 0.0f);
		vel = Vec3(0.0f, 0.0f, 0.0f);
	}
}

void IceObject::CalcWeightedVector(int pIndx, Vec3& pos, Vec3& vel)
{
	//���ꂼ��̃x�N�g�������������ς��Ƃ�
	pos = Vec3(0.0, 0.0, 0.0);
	vel = Vec3(0.0, 0.0, 0.0);

	int cIndx, oIndx;
	int clusterNum = GetPtoCIndx(pIndx);	//�e�X�g�Ƃ������Ƃ�
	int average = 1.0f/(float)clusterNum;

	//�ʒu�C���x�x�N�g���Z�o
	float deformAmountSum = 0.0f;

	for(int i = 0; i < GetPtoCIndx(pIndx); ++i)
	{
		cIndx = GetPtoC(pIndx, i, 0);
		oIndx = GetPtoC(pIndx, i, 1);

		if(cIndx == -1 || oIndx == -1){	continue;	}

		float defAmount = m_iceSM[cIndx]->GetDefAmount();
		//defAmount *= defAmount;				//����
		//defAmount *= defAmount * defAmount;	//�O���

		pos += m_iceSM[cIndx]->GetVertexPos(oIndx) * defAmount;
		vel += m_iceSM[cIndx]->GetVertexVel(oIndx) * defAmount;

		deformAmountSum += defAmount;
	}

	pos /= deformAmountSum;
	vel /= deformAmountSum;
	

	//for(int i = 0; i < GetPtoCIndx(pIndx); ++i)
	//{
	//	cIndx = GetPtoC(pIndx, i, 0);
	//	oIndx = GetPtoC(pIndx, i, 1);

	//	if(cIndx == -1 || oIndx == -1){	continue;	}

	//	float defAmount = m_iceSM[cIndx]->GetDefAmount()/deformAmountSum;
	//	//defAmount *= defAmount;		//����

	//	m_iceSM[cIndx]->SetDefPriority(defAmount*10.0f/* - average*/);	//���ςƂǂ̂��炢�������邩��\��
	//}

	//TODO::�t�o�[�W�����@���܂������Ȃ��@������deformAmountSum��0������
	//float deformAmountSum = 0.0f;
	//for(int i = 0; i < GetPtoCIndx(pIndx); ++i)
	//{
	//	cIndx = GetPtoC(pIndx, i, 0);
	//	oIndx = GetPtoC(pIndx, i, 1);

	//	if(cIndx == -1 || oIndx == -1){	continue;	}

	//	deformAmountSum += m_iceSM[cIndx]->GetDefAmount();
	//}

	////�ʒu�C���x�x�N�g���Z�o
	//for(int i = 0; i < GetPtoCIndx(pIndx); ++i)
	//{
	//	cIndx = GetPtoC(pIndx, i, 0);
	//	oIndx = GetPtoC(pIndx, i, 1);

	//	if(cIndx == -1 || oIndx == -1){	continue;	}

	//	float defAmount = m_iceSM[cIndx]->GetDefAmount();
	//	pos += m_iceSM[cIndx]->GetVertexPos(oIndx) * (deformAmountSum - defAmount)/deformAmountSum;
	//	vel += m_iceSM[cIndx]->GetVertexVel(oIndx) * (deformAmountSum - defAmount)/deformAmountSum;
	//}
}

//SPH�@��SM�@�ŋ��߂����x�ƈʒu����`��ԁ@CPU
void IceObject::LinerInterPolationCPU(int pIndx, const Vec3& pos, const Vec3& vel)
{
	int sphIndx = pIndx*4;
	double intrps = 1.0 - m_fInterPolationCoefficience[pIndx];	//��ԌW��

	for(int i = 0; i < 3; i++)
	{
		s_sphPrtVel[sphIndx+i] = vel[i] * m_fInterPolationCoefficience[pIndx] + s_sphPrtVel[sphIndx+i] * intrps;
		s_sphPrtPos[sphIndx+i] = pos[i] * m_fInterPolationCoefficience[pIndx] + s_sphPrtPos[sphIndx+i] * intrps;
	}
}

//SPH�@��SM�@�ŋ��߂����x�ƈʒu����`��ԁ@���f��̓N���X�^ CPU
void IceObject::LinerInterPolationForClusterCPU(int pIndx, const Vec3& pos, const Vec3& vel)
{
	int sphIndx = pIndx*4;
	double intrps = 1.0-m_fInterPolationCoefficience[pIndx];	//��ԌW��
	float* p = Ice_SM::GetSldPosPointer();
	float* v = Ice_SM::GetSldVelPointer();

	v[pIndx*3+0] = vel[0] * m_fInterPolationCoefficience[pIndx] + s_sphPrtVel[sphIndx+0] * intrps;
	v[pIndx*3+1] = vel[1] * m_fInterPolationCoefficience[pIndx] + s_sphPrtVel[sphIndx+1] * intrps;
	v[pIndx*3+2] = vel[2] * m_fInterPolationCoefficience[pIndx] + s_sphPrtVel[sphIndx+2] * intrps;

	p[pIndx*3+0] = pos[0] * m_fInterPolationCoefficience[pIndx] + s_sphPrtPos[sphIndx+0] * intrps;
	p[pIndx*3+1] = pos[1] * m_fInterPolationCoefficience[pIndx] + s_sphPrtPos[sphIndx+1] * intrps;
	p[pIndx*3+2] = pos[2] * m_fInterPolationCoefficience[pIndx] + s_sphPrtPos[sphIndx+2] * intrps;
}




//�M����
void IceObject::StepHeatTransfer(const int* surfParticles, const vector<vector<rxNeigh>>& neights, float floor, float effRadius, const float* pos, const float* dens)
{
	vector<int> ids;														//�\�ʗ��q�̓Y����
	vector<float> dists;													//�\�ʗ��q�̋���

	//������
	m_heatTransfer->resetNeighborhoodsId();
	m_heatTransfer->resetNeighborhoodsDist();

	//�ߖT���q�̓Y�����Ƌ����̐ݒ�
	for( int i = 0; i < sm_particleNum; i++ )
	{
		//�\�ʗ��q������@�P�Ȃ�Ε\�ʗ��q�C�O�Ȃ�������q�@���̍ۂ̐��l�͋ߖT���q������\��
		//�ߖT���q�����ŕ\�ʐς��ߎ�
		if( surfParticles[i] == 1 )
		{
			//���ɋ߂��C���͂̍������q�́C�\�ʂł͂Ȃ���ʂƂ���
			//���q���ɂ���ăp�����[�^��ς��Ȃ��Ƃ����Ȃ��D
			//�\�ʗ��q����ɕs�������̂ŏC�����K�v�D
			//0.75�ŉ��Q�i�Ƃ��
			if(pos[i*4+1] < floor+effRadius*0.2)				//1331�@���P�i
			{
				if(dens[i] < 950.0)
				{
					m_heatTransfer->setSurfaceParticleNums(i, (int)( neights[i].size() ));		//�\�ʈ���
				}
				else
				{
					m_heatTransfer->setSurfaceParticleNums(i, -1);								//��ʈ���
				}
			}
			else
			{
				m_heatTransfer->setSurfaceParticleNums(i, (int)( neights[i].size() ));
			}
		}
		else
		{
			m_heatTransfer->setSurfaceParticleNums(i, -1);
		}
		
		//������
		ids.clear();
		dists.clear();

		for( unsigned j = 0; j < neights[i].size(); j++)
		{
			if( i == (int)( neights[i][j].Idx ) ) continue;							//�������g���Ȃ�
			ids.push_back( (int)( neights[i][j].Idx ) );
			dists.push_back( (float)( neights[i][j].Dist ) );
		}

		m_heatTransfer->AddNeighborhoodsId( ids );
		m_heatTransfer->AddNeighborhoodsDist( dists );
	}

	//�M�����v�Z
	m_heatTransfer->heatAirAndParticle(); 		 										//�M�����@��C�Ɨ��q
	m_heatTransfer->heatParticleAndParticle(dens, effRadius);	//�M�����@���q��
	m_heatTransfer->calcTempAndHeat();													//�M�ʂ̉��x�ϊ��C���x�̔M�ʕϊ�

//�f�o�b�O
	//for(int i = 0; i < sm_particleNum; i++)
	//{
	//	cout << m_heatTransfer->getTemps()[i] << endl;
	//}
}

//�`�M�����C���ω�������ő̃f�[�^�̍X�V
void IceObject::StepIceStructure()
{
	StepMelting();			//�Z��
	//StepFreezing();		//�Ì�
}

void IceObject::StepMelting()
{
	vector<unsigned> prtList;			//�Z���������q�̏W��
	vector<unsigned> clusterList;		//�Ē�`����N���X�^�̏W��
	vector<unsigned> cLayerList;		//�Ē�`����N���X�^�̃��C���[
	vector<unsigned> tetraList;			//�Ē�`����l�ʑ̂̏W��
	vector<unsigned> tLayerList;		//�Ē�`����l�ʑ̂̃��C���[

	SearchMeltParticle(prtList);		//RXTIMER("SearchMeltParticle");	//�Z�𗱎q�T��

	//�������v�Z���Ԃ̖w�ǂ��߂Ă���
	//m_iceStrct->StepObjMelt
	//	(	prtList,
	//		clusterList,
	//		tetraList,
	//		cLayerList,
	//		tLayerList);				//�Z������
	//RXTIMER("m_iceStrct->StepObjMelt");

	//m_iceStrct->TestStepObjMelt
	//	(	prtList,
	//		clusterList,
	//		tetraList,
	//		cLayerList,
	//		tLayerList);				//�Z������
	
	QueryCounter counter1;
	QueryCounter counter2;

	counter1.Start();

	ReConstructCluster(prtList, clusterList);	//RXTIMER("ReConstructCluster");	//�N���X�^�̍č\�z

	double end1 = counter1.End();

	//�^���v�Z����N���X�^�̍Ē�`
	if(prtList.size() == 0) return;

	counter2.Start();

	UpdateSelectCluster();

	double end2 = counter2.End();

	cout << "Time" << endl;
	cout << "ReConstructCluster:	" << end1 << endl;
	cout << "UpdateSelectCLuster:	" << end2 << endl;
}

void IceObject::StepFreezing()
{
	vector<unsigned> vuFreezePrtList;			//�Ìł������q�̏W��

	SearchMeltParticle(vuFreezePrtList);		//�Ìŗ��q�T��
	m_iceStrct->StepObjFreeze();
}

void IceObject::SearchMeltParticle(vector<unsigned>& pList)
{
	for(int i = 0; i < sm_particleNum; ++i)
	{	
		if(m_iceSM[i]->GetNumVertices() == 0)			continue;	//�N���X�^�������Ă��Ȃ�
		if( m_heatTransfer->getPhaseChange(i) != 1 )	continue;	//���]�ڂ̏����𖞂����Ă���
		if( m_heatTransfer->getPhase(i) != 2 )			continue;	//���ւƑ��]�ڂ��Ă���
		//if( m_iceStrct->GetPtoCNum(i) == 0 )			continue;	//�N���X�^�Ɋ܂܂�Ă��Ȃ�
		//if( m_iceStrct->GetPtoTNum(i) == 0 )			continue;	//�l�ʑ̂Ɋ܂܂�Ă��Ȃ�

		//if(pList.size() > 500){	break;	}							//�Z�𗱎q���̐���

		m_fInterPolationCoefficience[i] = 0.0f;						//���`��Ԃ��Ȃ�
		m_heatTransfer->setPhaseChange(i, 0);						//���]�ڂ��I��������Ƃ�`����
		pList.push_back(i);											//�Z�𗱎q�̋L�^
	}

//�f�o�b�O
	if(pList.size() == 0) return;

	cout << __FUNCTION__ << endl;

	//for(unsigned i = 0; i < pList.size(); i++)
	//{
	//	cout << pList[i] << endl;
	//}

	cout << "pList.size() = " << pList.size() << endl;
}

void IceObject::SearchFreezeParticle(vector<unsigned>& pList)
{
//	for(int i = 0; i < m_pPS->GetNumParticles(); i++)
//	{	
//		if(m_ht->getPhaseChange(i) != 1)				continue;	//���]�ڂ̏����𖞂����Ă��Ȃ��ꍇ�͖߂�
//		if(m_ht->getPhase(i) != -2)						continue;	//�X�ւƑ��]�ڂ��Ă��Ȃ��ꍇ�͖߂�
//		if(m_ice->GetParticleNum() <= i)				continue;	//�Z���݂̂̎����̂Ƃ��ɕK�v�ɂȂ�D
//		if(m_ice->GetPtoCNum(i) != 0)					continue;	//�N���X�^�Ɋ܂܂�Ă���
//		if(m_ice->GetPtoTNum(i) != 0)					continue;	//�N���X�^�Ɋ܂܂�Ă���
////		if(pList.size() > 1){	break;	}							//�Ìŗ��q���̐���
//		
//		pList.push_back(i);											//�Ìŗ��q�̋L�^
//	}
}

void IceObject::ReConstructCluster(vector<unsigned>& pList, vector<unsigned>& cList)
{
	int pListSize = pList.size();
	int cListSize = cList.size();

	if(pListSize == 0 /*|| cListSize == 0*/){	return; }

	////���Ƃœ񕪒T�����邽�߂�pList��sort����
	//std::sort(pList.begin(), pList.end());

	vector<bool> clusterFlags(sm_particleNum, true);
	for(unsigned i = 0; i < pList.size(); i++)
	{
		clusterFlags[pList[i]] = false;
	}

	//�N���X�^���痱�q����菜��
	for(int i = 0; i < pListSize; i++)
	{
		int pIndx = pList[i];

		//�Z���N���X�^��layer0���擾
		vector<int> connectPIndx;
		for(int oIndx = 0; oIndx < GetMoveObj(pIndx)->GetIndxNum(); oIndx++)
		{
			int layer = GetMoveObj(pIndx)->GetLayer(oIndx);

			if(layer == -1) continue;
			if(layer != 0)	break;

			int ipIndx = GetMoveObj(pIndx)->GetParticleIndx(oIndx);

			connectPIndx.push_back(ipIndx);
		}

		//���̑��̃N���X�^����Z�𗱎q����菜��
		int oIndx = 0;
		int ccpIndx = 0;
		int coIndx = 0;
		int ccIndx = 0;
		int connectPIndxSize = connectPIndx.size();
		Ice_SM* moveObj;
		vector<unsigned>::iterator begin = pList.begin();
		vector<unsigned>::iterator end = pList.end();

		#pragma omp parallel
		{
		#pragma omp for private(oIndx, ccpIndx, coIndx, moveObj, ccIndx, begin, end)
			for(int cIndx = 0; cIndx < sm_particleNum; cIndx++)
			{
				//if(std::binary_search(begin, end, (unsigned)cIndx)){	continue;	}	//�d��
				if(clusterFlags[cIndx] == false)	continue;

				moveObj = GetMoveObj(cIndx);
				oIndx = moveObj->SearchIndx(pIndx);

				if(oIndx == MAXINT) continue;	//���݂��Ȃ�

				moveObj->Remove(oIndx);

				//�Z�𗱎q��layer0�ƂȂ闱�q��S�Ď�菜��
				//������������layer0�͎c��
				for(ccIndx = 0; ccIndx < connectPIndxSize; ccIndx++)
				{
					ccpIndx = connectPIndx[ccIndx];
					coIndx = moveObj->SearchIndx(ccpIndx);
	
					if(coIndx == MAXINT)				continue;	//���݂��Ȃ�
					if(moveObj->GetLayer(coIndx) == 0)	continue;	//������layer0�͎c��

					moveObj->Remove(coIndx);
				}

				moveObj->CalcCm();
			}
		}//end #pragma omp parallel
	}

	#pragma omp parallel
	{
	#pragma omp for
		for(int pIndx = 0; pIndx < pListSize; pIndx++)
		{
			GetMoveObj(pList[pIndx])->Clear();
		}
	}

	//vector<int> checkTList;
	//int j = 0, k = 0, l = 0;
	//int icIndx = 0;
	//int jtIndx = 0, joIndx = 0;
	//int kpIndx = 0, ktIndx = 0, klIndx = 0;
	//int lpIndx = 0;
	//int pNum = 0;

	////�Ē�`�N���X�^�̏�����
	//#pragma omp parallel
	//{
	//#pragma omp for
	//	for(int i = 0; i < cListSize; ++i)
	//	{
	//		if(std::find(pList.begin(), pList.end(), cList[i]) != pList.end()){	continue;	}
	//		GetMoveObj(cList[i])->Clear();
	//	}
	//}//end #pragma omp parallel

	////�N���X�^�̍Ē�`
	//#pragma omp parallel
	//{
	//#pragma omp for private(checkTList, j, k, l, icIndx, jtIndx, joIndx, kpIndx, ktIndx, klIndx, lpIndx, pNum)
	//	for(int i = 0; i < cListSize; ++i)
	//	{
	//		checkTList.clear();
	//		icIndx = cList[i];
	//		if(std::find(pList.begin(), pList.end(), icIndx) != pList.end()){	continue;	}

	//		//�N���X�^���Ē�`����ہC��{�ƂȂ闱�q��������l�ʑ̂��珉�����q�𓾂�D
	//		//�N���X�^�ԍ��������q�ԍ��Ȃ̂ɒ���
	//		//�ȉ����֐��ɂ���ƁC�G���[���o�Ă��܂������Ȃ�
	//		for(j = 0; j < GetPtoTIndx(icIndx); j++)
	//		{
	//			jtIndx = GetPtoT(icIndx, j, 0);
	//			joIndx = GetPtoT(icIndx, j, 1);

	//			if(jtIndx == -1 || joIndx == -1){ continue;	}
	//			if(std::find(checkTList.begin(), checkTList.end(), jtIndx) != checkTList.end())
	//			{	continue;	}
	//			else
	//			{	checkTList.push_back(jtIndx);	}

	//			//�l�ʑ̂̑S�Ă̗��q���N���X�^�ɓo�^
	//			for(k = 0; k < GetTtoPIndx(jtIndx); k++)
	//			{
	//				kpIndx = GetTtoP(jtIndx, k);

	//				if(kpIndx == -1){	continue;	}
	//				if(GetMoveObj(icIndx)->CheckIndx(kpIndx)){	continue;	}

	//				pNum = GetMoveObj(icIndx)->GetNumVertices();

	//				GetMoveObj(icIndx)->AddVertex(Vec3(s_sphPrtPos[kpIndx*4+0], s_sphPrtPos[kpIndx*4+1], s_sphPrtPos[kpIndx*4+2]), 1.0, kpIndx);
	//				//GetMoveObj(icIndx)->SetVelocity(pNum, Vec3(s_sphPrtVel[kpIndx*4+0], s_sphPrtVel[kpIndx*4+1], s_sphPrtVel[kpIndx*4+2]));
	//				GetMoveObj(icIndx)->SetAlphas(pNum, 1.0);
	//				GetMoveObj(icIndx)->SetBetas (pNum, 0.0);
	//				GetMoveObj(icIndx)->SetLayer (pNum, 0);
	//			}
	//		}

	//		//�ߖT�l�ʑ̂�layer���ǂ�C���q��ǉ����Ă���
	//		//TODO::�s����ɂȂ�Ȃ�Clayer�������ق�beta��������
	//		for(j = 0; j < GetPtoTIndx(icIndx); ++j)
	//		{
	//			jtIndx = GetPtoT(icIndx, j, 0);
	//			joIndx = GetPtoT(icIndx, j, 1);

	//			if(jtIndx == -1 || joIndx == -1){ continue;	}

	//			for(k = 0; k < GetNTNum(jtIndx); k++)
	//			{
	//				ktIndx = GetNeighborTetra(jtIndx, k, 0);
	//				klIndx = GetNeighborTetra(jtIndx, k, 1);

	//				if(ktIndx == -1 || klIndx == -1){	continue;	}
	//				if(std::find(checkTList.begin(), checkTList.end(), ktIndx) != checkTList.end())
	//				{	continue;	}
	//				else
	//				{	checkTList.push_back(ktIndx);	}

	//				//�l�ʑ̂̑S�Ă̗��q���N���X�^�ɓo�^
	//				for(l = 0; l < GetTtoPIndx(ktIndx); l++)
	//				{
	//					lpIndx = GetTtoP(ktIndx, l);

	//					if(lpIndx == -1){	continue;	}
	//					if(GetMoveObj(icIndx)->CheckIndx(lpIndx)){	/*cout << "contineue kpIndx = " << kpIndx << endl;*/ continue;	}

	//					pNum = GetMoveObj(icIndx)->GetNumVertices();

	//					GetMoveObj(icIndx)->AddVertex( Vec3(s_sphPrtPos[lpIndx*4+0], s_sphPrtPos[lpIndx*4+1], s_sphPrtPos[lpIndx*4+2] ), 1.0, lpIndx);
	//					//GetMoveObj(icIndx)->SetVelocity(pNum, Vec3(s_sphPrtVel[lpIndx*4+0], s_sphPrtVel[lpIndx*4+1], s_sphPrtVel[lpIndx*4+2]));
	//					GetMoveObj(icIndx)->SetAlphas(pNum, 1.0);
	//					GetMoveObj(icIndx)->SetBetas (pNum, 0.0);
	//					GetMoveObj(icIndx)->SetLayer (pNum, klIndx);
	//				}
	//			}
	//		}
	//	}
	//}//end #pragma omp parallel


	////�Z���N���X�^�����Z�b�g
	//#pragma omp parallel
	//{
	//#pragma omp for
	//	for(int i = 0; i < pListSize; ++i)
	//	{
	//		GetMoveObj(pList[i])->Clear();
	//	}
	//}//end #pragma omp parallel


	//////�N���X�^�̌ő̏���������
	////for(int i = 0; i < (int)cList.size(); ++i)
	////{
	////	m_iceStrct->ClearCtoP(cList[i]);
	////}

	//////�N���X�^�̌ő̏����ēx�o�^
	////for(unsigned i = 0; i < cList.size(); ++i)
	////{
	////	if(std::find(pList.begin(), pList.end(), cList[i]) != pList.end())
	////	{
	////		continue;	
	////	}

	////	int cIndx = cList[i];

	////	vector<int> pIndxList;
	////	int pNum = GetMoveObj(cIndx)->GetNumVertices();
	////	int* pLayerList = new int[pNum];									//���q�̏������C���[
	////
	////	//layer�̂��߂̃R�s�[
	////	for(int i = 0; i < pNum; i++)
	////	{
	////		pLayerList[i] = GetMoveObj(cIndx)->GetLayer(i);		//���q�����w�ڂ̋ߖT�Ȃ̂����擾
	////	}
	////
	////	for(int i = 0; i < pNum; i++)
	////	{
	////		int pIndx = GetMoveObj(cIndx)->GetParticleIndx(i);
	////		int freeIndx = m_iceStrct->GetPtoCFreeIndx(pIndx);
	////
	////		m_iceStrct->SetPtoC(pIndx, freeIndx, cIndx, i, pLayerList[i]);		//���q���������Ă���l�ʑ̂�o�^
	////		pIndxList.push_back(pIndx);
	////	}
	////
	////	m_iceStrct->SetCtoP(cIndx, pIndxList, pLayerList);						//�l�ʑ̂��܂�ł��闱�q��o�^
	////
	////	delete[] pLayerList;
	////}
}

void IceObject::ResetSelectCluster()
{
	//�S�ẴN���X�^��I������
	for(int cIndx = 0; cIndx < sm_clusterNum; cIndx++)
	{
		m_iceStrct->UpdateMotionCalcCluster(cIndx, 1);		
	}
}

void IceObject::UpdateSelectCluster()
{
	//�I���N���X�^�̃��Z�b�g
	for(int cIndx = 0; cIndx < sm_clusterNum; cIndx++)
	{
		m_iceStrct->UpdateMotionCalcCluster(cIndx, 0);		
	}

	//�N���X�^�W���̏�����
	vector<unsigned> clusters;
	for(int cIndx = 0; cIndx < sm_clusterNum; cIndx++)
	{
		if(m_iceSM[cIndx]->GetNumVertices() <= 0){	continue;	}
		clusters.push_back(cIndx);
	}

	//�I��
	while(clusters.size() != 0)
	{
		unsigned cIndx = *clusters.begin();

		m_iceStrct->UpdateMotionCalcCluster(cIndx, 1);

		//�ߖT�N���X�^����菜��
		for(int indx = 0; indx < m_iceSM[cIndx]->GetNumVertices(); indx++)
		{
			unsigned icIndx = (unsigned)m_iceSM[cIndx]->GetParticleIndx(indx);

			//stl�ō폜�́Cerase��remove��g�ݍ��킹�čs��
			clusters.erase(remove(clusters.begin(), clusters.end(), icIndx), clusters.end());  
		}

		clusters.erase(remove(clusters.begin(), clusters.end(), cIndx), clusters.end());  
	}

//�f�o�b�O
	//DebugUpdateSelectCluster();
}

//---------------------------------------------�f�o�b�O------------------------------

//�l�ʑ�
void IceObject::DebugTetraInfo()
{	cout << __FUNCTION__ << endl;

	unsigned tetraSize = IceTetrahedra::GetInstance().GetTetraListSize();

	//�l�ʑ́����q
	for(unsigned i = 0; i < tetraSize; i++ )
	{
		m_iceStrct->DebugTtoP(i);
	}

	//���q���l�ʑ�
	for(int i = 0; i < sm_particleNum; i++)
	{
		m_iceStrct->DebugPtoT(i);
	}

	//�ߖT�l�ʑ�
	for(unsigned i = 0; i < tetraSize; i++ )
	{
		m_iceStrct->DebugNeighborTetra(i);
	}
}

//�N���X�^
void IceObject::DebugClusterInfo()
{	cout << __FUNCTION__ << endl;

	//�ܗL���q���̑����C���ρC�ŏ��l�C�ő�l
	int sum = 0;
	int max = m_iceSM[0]->GetNumVertices();
	int min = m_iceSM[0]->GetNumVertices();

	for(int i = 0; i < sm_clusterNum; i++)
	{
		int num =  m_iceSM[i]->GetNumVertices();
		
		if(max < num){	max = num;	}

		if(min > num){	min = num;	}

		sum += num;
	}

	cout << "sum = " << sum << " ave = " << (float)(sum/sm_clusterNum) << endl;
	cout << "max = " << max << " min = " << min << endl;
}

void IceObject::DebugDeformationAmount()
{	//cout << __FUNCTION__ << endl;
	
	int deformClusterNum = 0;
	for(int i = 0; i < sm_clusterNum; ++i)
	{	
		if(GetPtoCNum(i) == 0){	continue;	}

		float defAmount = m_iceSM[i]->GetDefAmount();
		if(defAmount < 0.5f){	continue;	}

		//cout << "i = " << i << "::" <<  defAmount << endl;

		deformClusterNum++;
	}

	//cout << "Deformation Amount deformClusterNum::" << deformClusterNum << endl;

	ofstream ofs( "DeformationAmount.txt", ios::app);
	ofs << deformClusterNum << endl;
}

void IceObject::DebugDeformationAverage()
{
	float deformClusterAverage = 0.0f;
	for(int i = 0; i < sm_clusterNum; ++i)
	{	
		if(GetPtoCNum(i) == 0){	continue;	}

		float defAmount = m_iceSM[i]->GetDefAmount();

		deformClusterAverage += defAmount;
	}

	//cout << "Deformation Average deformClusterNum::" << deformClusterNum << endl;

	ofstream ofs( "DeformationAverage.txt", ios::app);
	ofs << deformClusterAverage/(float)sm_clusterNum << endl;
}

void IceObject::DebugObjMoveUsePathWithGPU()
{
//	QueryCounter counter;
//	QueryCounter counter2;	
//	QueryCounter counter3;
//
//counter.Start();
//
//	//prefixSum�̍X�V
//	m_SurfSm.SetDevicePointer(sd_sphPrtPos, sd_sphPrtVel);	//VBO���g���Ă���̂ŁC���ꂪ�Ȃ��ƃG���[���o��̂ɒ���
//	m_SurfSm.UpdatePrefixSumGPU();
//
//	double prefixSumTime = counter.End()/1000.0;
////cout << "UpdatePrefixSumGPU�F	" << prefixSumTime << endl;
//
//counter2.Start();
//
//	//�N���X�^�̃p�����[�^�X�V
//	//TestUpdateSMFromPath();
//	LauchUpdateSMFromPath
//	(
//		sm_clusterNum,
//		sd_sphPrtPos,
//		sd_sphPrtVel,
//		//SM-------------------------------------------------
//		Ice_SM::GetOrgPosPointer(),
//		Ice_SM::GetDevicePosPointer(),
//		Ice_SM::GetOrgCmPointer(),
//		Ice_SM::GetCurCmPointer(),
//		Ice_SM::GetDeviceApqPointer(),
//		Ice_SM::GetDeviceVelPointer(),
//		Ice_SM::GetDeviceIndexesPointer(),
//		Ice_SM::GetDeviceIndexSetPointer(),
//		//Path-----------------------------------------------
//		m_SurfSm.GetDevicePRTtoPTHPointer(),
//		m_SurfSm.GetDevicePTHandPrfxSetPointer(),
//		m_SurfSm.GetDecvicePrfxPos(),
//		m_SurfSm.GetDecvicePrfxApq(),
//		//IceStruct------------------------------------------
//		m_iceStrct->GetDeviceCtoPointer(),
//		m_iceStrct->GetDeviceCtoPNumPointer(),
//		m_iceStrct->GetCtoPMax(),
//		2,
//		0.01		//TODO: 0.02�ł́H���낢��΂�΂���ۂ�
//	);
//
//	double updatePosApqTime = counter2.End()/1000.0;
////cout << "LauchUpdateSMFromPath�F	" << updatePosApqTime << endl;
//
//counter3.Start();
//
//	//�^���v�Z
//	Ice_SM::SetDevicePosPointer(sd_sphPrtPos);				//VBO���g���Ă���̂ŁC���ꂪ�Ȃ��ƃG���[���o��̂ɒ���
//	Ice_SM::UpdateUsePathGPU();
//
//	double updateSMTime = counter3.End()/1000.0;
////cout << "UpdateUsePathGPU�F	" << updateSMTime << endl;
//
////cout << "���v�F		" << prefixSumTime+updatePosApqTime+updateSMTime << endl;
//
//
////�z�X�g���̃f�[�^��]���������ʂ��_���v
//	//�Â��t�@�C�����폜
//	
//
//	//prefixSum
//	ofstream ofs( "LOG_prefixSumTime.txt", ios::app);
//	ofs << prefixSumTime << endl;
//
//	//updatePosApq
//	ofstream ofs1( "LOG_updateTime.txt", ios::app);
//	ofs1 << updatePosApqTime << endl;
//
//	//ShapeMatching
//	ofstream ofs2( "LOG_SMTime.txt", ios::app);
//	ofs2 << updateSMTime << endl;
}

void IceObject::DebugUpdateSelectCluster()
{	cout << __FUNCTION__ << endl;
	for(int cIndx = 0; cIndx < sm_particleNum; cIndx++)
	{
		if(m_iceStrct->GetMotionCalcCluster(cIndx) == 0){	continue;	}
		cout << "cIndx = " << cIndx << " , " << m_iceStrct->GetMotionCalcCluster(cIndx) << endl;
	}
}

//�e��@�̃N���X����\��
void IceObject::DebugNowMoveMethod()
{	cout << __FUNCTION__ << endl;

	cout << "CalcMethod	" << typeid(*m_iceCalcMethod).name() << endl;				//�v�Z���@
	cout << "JudgeMove	" << typeid(*m_iceJudeMove).name() << endl;					//�^���v�Z�Ώۂ𔻒�
	cout << "ClsuterMove	" << typeid(*m_iceClsuterMove).name() << endl;			//�^���v�Z���@
	cout << "IntrpoJudge	" << typeid(*m_iceInterPolationJudge).name() << endl;	//�ŏI�������ʂɗp����Ώۂ𔻒�
	cout << "InterPolation	" << typeid(*m_iceInterPolation).name() << endl;		//�ŏI�������ʂ����߂�
}

//---------------------------------------------------�e�X�g---------------------------------------------
void IceObject::TestUpdateSMFromPath()
{	cout << __FUNCTION__ << endl;

//GPU��
	LauchUpdateSMFromPath
	(
		sm_clusterNum,
		sd_sphPrtPos,
		sd_sphPrtVel,
		//SM-------------------------------------------------
		Ice_SM::GetOrgPosPointer(),
		Ice_SM::GetDevicePosPointer(),
		Ice_SM::GetOrgCmPointer(),
		Ice_SM::GetCurCmPointer(),
		Ice_SM::GetDeviceApqPointer(),
		Ice_SM::GetDeviceVelPointer(),
		Ice_SM::GetDeviceIndexesPointer(),
		Ice_SM::GetDeviceIndexSetPointer(),
		//Path-----------------------------------------------
		m_SurfSm->GetDevicePRTtoPTHPointer(),
		m_SurfSm->GetDevicePTHandPrfxSetPointer(),
		m_SurfSm->GetDecvicePrfxPos(),
		m_SurfSm->GetDecvicePrfxApq(),
		//IceStruct------------------------------------------
		m_iceStrct->GetDeviceCtoPointer(),
		m_iceStrct->GetDeviceCtoPNumPointer(),
		m_iceStrct->GetCtoPMax(),
		2,
		0.01		//TODO: 0.02�ł́H���낢��΂�΂���ۂ�
	);

//CPU��
	//prefixSum�̍X�V
	m_SurfSm->UpdatePrefixSum();

	//�N���X�^�̃p�����[�^�X�V
	#pragma omp parallel for
	for(int i = 0; i < sm_clusterNum; ++i)
	{	
		if(GetPtoCNum(i) == 0){	continue;	}
		m_iceSM[i]->SetNowCm(m_SurfSm->CalcCmSum(i));		//�d�S�̍X�V
		m_iceSM[i]->SetApq(m_SurfSm->CalcApqSum(i));		//�ό`�s��̍X�V
	}

	////CurCm�̔�r
	//float* dCurCm = new float[MAXCLUSTER * SM_DIM];
	//cudaMemcpy(dCurCm, Ice_SM::GetCurCmPointer(), sizeof(float) * MAXCLUSTER * SM_DIM, cudaMemcpyDeviceToHost);

	//for(int iCrstr = 0; iCrstr < MAXCLUSTER; iCrstr++)
	//{
	//	int cIndx = iCrstr * 3;
	//	cout << "GPU:[" << iCrstr << "] = " << endl;
	//	cout <<	"(" << dCurCm[cIndx+0] << ", " << dCurCm[cIndx+1] << ", " << dCurCm[cIndx+2] << ")" << endl;
	//	cout << endl;

	//	Vec3 cm = m_iceSM[iCrstr]->GetCm();
	//	cout << "CPU:[" << iCrstr << "] = " << endl;
	//	cout <<	"(" << cm[0] << ", " << cm[1] << ", " << cm[2] << ")" << endl;
	//	cout << endl;
	//}

	//delete[] dCurCm;

	////CurApq�̔�r
	//float* dCurApq = new float[MAXCLUSTER * 9];
	//cudaMemcpy(dCurApq, Ice_SM::GetDeviceApqPointer(), sizeof(float) * MAXCLUSTER * 9, cudaMemcpyDeviceToHost);

	//for(int iCrstr = 0; iCrstr < MAXCLUSTER; iCrstr++)
	//{
	//	int cIndx = iCrstr * 9;
	//	cout << "GPU:[" << iCrstr << "] = " << endl;
	//	cout <<	"(" << dCurApq[cIndx+0] << ", " << dCurApq[cIndx+1] << ", " << dCurApq[cIndx+2] << ")" << endl;
	//	cout <<	"(" << dCurApq[cIndx+3] << ", " << dCurApq[cIndx+4] << ", " << dCurApq[cIndx+5] << ")" << endl;
	//	cout <<	"(" << dCurApq[cIndx+6] << ", " << dCurApq[cIndx+7] << ", " << dCurApq[cIndx+8] << ")" << endl;
	//	cout << endl;

	//	rxMatrix3 apq = m_iceSM[iCrstr]->GetApq();
	//	cout << "CPU:[" << iCrstr << "] = " << endl;
	//	cout <<	"(" << apq(0, 0) << ", " << apq(0, 1) << ", " << apq(0, 2) << ")" << endl;
	//	cout <<	"(" << apq(1, 0) << ", " << apq(1, 1) << ", " << apq(1, 2) << ")" << endl;
	//	cout <<	"(" << apq(2, 0) << ", " << apq(2, 1) << ", " << apq(2, 2) << ")" << endl;
	//	cout << endl;
	//}

	//delete[] dCurApq;

}

void IceObject::TestStepInterPolation()
{
	Vec3 pos,vel;

	#pragma omp parallel for private(pos, vel)
	for(int i = 0; i < sm_particleNum; ++i)
	{
		if(GetParticleNum() <= i){	continue;	}	//�Z���݂̂̎����̂Ƃ��ɕK�v�ɂȂ�D
		if(GetPtoCNum(i) <= 0){		continue;	}
	
		//�ő̉^���̍ŏI�ʒu�v�Z
		CalcAverageCPU(i, pos, vel);

		//�t�̂ƌő̂̕��
		LinerInterPolationCPU(i, pos, vel);			//���`���
	}
}