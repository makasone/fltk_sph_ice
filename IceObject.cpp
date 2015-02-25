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

//�e�����a�����ɉ^���v�Z����N���X�^��I���@�Ȃ�ׂ��a�ɂȂ�悤�ɑI��
void IceObject::InitSelectCluster(float radius)
{
	m_iceStrct->SetSelectRadius(radius);
	m_iceStrct->InitSelectCluster(m_iceSM);

	//DebugUpdateSelectCluster();
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
void IceObject::InitCluster(Vec3 boundarySpaceHigh, Vec3 boundarySpaceLow, float timeStep, int itr, const vector<vector<rxNeigh>>& neights)
{	cout << __FUNCTION__ << endl;

	//�ϐ��̏�����
	for(vector<Ice_SM*>::iterator it = m_iceSM.begin(); it != m_iceSM.end(); ++it)
	{
		if(*it) delete *it;
	}
	
	sm_clusterNum = 0;	//�N���X�^���̏�����

	Ice_SM::SetIterationNum(itr);
	Ice_SM::SetItrStiffness(1.0f);	//TODO::�K��
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
		m_iceSM[sm_clusterNum]->SetStiffness(1.0, 0.0);

		//�l�ʑ̃��X�g�����ɁC�N���X�^�֗��q��o�^
#ifdef USE_NEIGHT
		SetClusterMoveInfoFromNeight(i, neights);
#else
		SetClusterMoveInfo(i);
#endif
		sm_clusterNum++;
	}

cout << __FUNCTION__ << ", check2" << endl;

	//���ʂ̏C��
	UpdateParticleMass_Normal();		//����
	//UpdateParticleMass_Average();		//���ρ@LSM�̕��@
	//UpdateParticleMass_Direction();	//�����x�N�g���Ƃ̗ގ��x�ŏd�ݕt���@Ijiri��̕��@

	////MakeOneCluster();
	
	Ice_SM::InitFinalParamPointer(sm_clusterNum);

//�f�o�b�O
	//DebugClusterInfo();
	//DebugNeights(neights);
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

	m_iceStrct->InitSelectClusterFromClusterSet(m_iceSM);

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
	m_fpStepObjMove = NULL;
	m_iceClsuterMove = NULL;
	m_iceJudeMove = NULL;
	m_iceConvolution = NULL;
	m_iceConvoJudge = NULL;
	m_iceCalcMethod = NULL;

	//�S�ăm�[�}��
	ChangeMode_Normal();
	ChangeMode_ClusterMove_Normal();
	ChangeMode_JudgeMove_Normal();
	ChangeMode_Convolution_Normal();
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
void IceObject::ChangeMode_Debug()
{	FUNCTION_NAME

	m_fpStepObjMove = &IceObject::StepObjMoveCPUDebug;

	//�f�[�^�t�H���_�̒��g��S�č폜
	string path = RESULT_DATA_PATH;
	DebugClearDirectry(path);
}

void IceObject::ChangeMode_Normal()
{	FUNCTION_NAME
	
	m_fpStepObjMove = &IceObject::StepObjMoveCPUNormal;
}

//�v�Z��@
void IceObject::ChangeMode_CalcMethod_Normal()
{	FUNCTION_NAME

	if(m_iceCalcMethod != NULL) delete m_iceCalcMethod;
	m_iceCalcMethod = new Ice_CalcMethod_Normal(m_iceClsuterMove);
}

void IceObject::ChangeMode_CalcMethod_Itr_Num()
{	FUNCTION_NAME

	if(m_iceCalcMethod != NULL) delete m_iceCalcMethod;
	m_iceCalcMethod = new Ice_CalcMethod_Iteration(m_iceSM, m_iceClsuterMove, m_iceConvolution);
}

void IceObject::ChangeMode_CalcMethod_Itr_Stiff()
{	FUNCTION_NAME;
	
	if(m_iceCalcMethod != NULL) delete m_iceCalcMethod;
	m_iceCalcMethod = new Ice_CalcMethod_Itr_Stiffness(m_iceSM, m_iceClsuterMove, m_iceConvolution);
}

void IceObject::ChangeMode_CalcMethod_Itr_Expand()
{	FUNCTION_NAME;

	if(m_iceCalcMethod != NULL) delete m_iceCalcMethod;
	m_iceCalcMethod = new Ice_CalcMethod_Itr_Expand(m_iceSM, m_iceClsuterMove, m_iceConvolution);
}

void IceObject::ChangeMode_CalcMethod_Itr_Exp_Stiff()
{	FUNCTION_NAME;

	if(m_iceCalcMethod != NULL) delete m_iceCalcMethod;
	m_iceCalcMethod = new Ice_CalcMethod_Itr_Exp_Stiff(m_iceSM, m_iceClsuterMove, m_iceConvolution);
}


//�^���v�Z�N���X�^�I��
void IceObject::ChangeMode_JudgeMove_Normal()
{	FUNCTION_NAME

	ResetSelectCluster();		//�I��I�^���v�Z�̏�����

	if(m_iceJudeMove != NULL) delete m_iceJudeMove;
	m_iceJudeMove = new Ice_JudgeMove_Normal(m_iceSM);

	m_iceClsuterMove->SetJudgeMove(m_iceJudeMove);
}

void IceObject::ChangeMode_JudgeMove_Spears()
{	FUNCTION_NAME

	if(m_iceJudeMove != NULL) delete m_iceJudeMove;
	m_iceJudeMove = new Ice_JudgeMove_Spears(m_iceSM, m_iceStrct);

	m_iceClsuterMove->SetJudgeMove(m_iceJudeMove);
}

//�^���v�Z
void IceObject::ChangeMode_ClusterMove_Normal()
{	FUNCTION_NAME

	if(m_iceClsuterMove != NULL) delete m_iceClsuterMove;
	m_iceClsuterMove = new Ice_ClusterMove_Normal(m_iceSM);

	m_iceClsuterMove->SetJudgeMove(m_iceJudeMove);

	if(m_iceCalcMethod != NULL)	m_iceCalcMethod->SetObjMove(m_iceClsuterMove);
}

void IceObject::ChangeMode_ClusterMove_Path()
{	FUNCTION_NAME

	if(m_SurfSm == NULL){	InitPath();	}	//�����ȉ^���v�Z�̂��߂̃p�X������

	if(m_iceClsuterMove != NULL) delete m_iceClsuterMove;
	m_iceClsuterMove = new Ice_ClusterMove_FastPath(m_iceSM, m_SurfSm);

	m_iceClsuterMove->SetJudgeMove(m_iceJudeMove);

	if(m_iceCalcMethod != NULL)	m_iceCalcMethod->SetObjMove(m_iceClsuterMove);
}

//�ŏI���ʕ�ԃN���X�^�I��
void IceObject::ChangeMode_IntrpJudge_Normal()
{	cout << __FUNCTION__ << endl;

	if(m_iceConvoJudge != NULL) delete m_iceConvoJudge;
	m_iceConvoJudge = new Ice_ConvoJudge_Normal(m_iceSM);

	m_iceConvolution->SetConvoJudge(m_iceConvoJudge);
}

void IceObject::ChangeMode_IntrpJudge_Spears()
{	cout << __FUNCTION__ << endl;

	if(m_iceConvoJudge != NULL) delete m_iceConvoJudge;
	m_iceConvoJudge = new Ice_ConvoJudge_Spears(m_iceSM, m_iceStrct);

	m_iceConvolution->SetConvoJudge(m_iceConvoJudge);
}

//�ŏI���ʕ��
void IceObject::ChangeMode_Convolution_Normal()
{	cout << __FUNCTION__ << endl;

	if(m_iceConvolution != NULL) delete m_iceConvolution;
	m_iceConvolution = new Ice_Convolution_Normal(m_iceSM, m_iceStrct);

	m_iceConvolution->SetConvoJudge(m_iceConvoJudge);
}

//�ό`�ʂŏd�݂����ĕ��
void IceObject::ChangeMode_Convolution_Weight()
{	cout << __FUNCTION__ << endl;

	if(m_iceConvolution != NULL) delete m_iceConvolution;
	m_iceConvolution = new Ice_Convolution_Weight(m_iceSM, m_iceStrct);

	m_iceConvolution->SetConvoJudge(m_iceConvoJudge);
}

//�����x�N�g���Ƃ̗ގ��x�ŏd�݂�t���ĕ��
void IceObject::ChangeMode_Convolution_Anisotropic()
{
	FUNCTION_NAME;

	if(m_iceConvolution != NULL) delete m_iceConvolution;
	m_iceConvolution = new Ice_Convolution_Anisotropic(m_iceSM, m_iceStrct);

	m_iceConvolution->SetConvoJudge(m_iceConvoJudge);
}

//�N���X�^�̉^���v�Z���̓o�^
void IceObject::SetClusterMoveInfo(int pIndx)
{
	//�ŏ��̈�ڂ͕K��pIndx�Ƃ���
	int pNum = m_iceSM[pIndx]->GetNumVertices();
	float mass = 1.0f;

	m_iceSM[pIndx]->AddVertex( Vec3(s_sphPrtPos[pIndx*4+0], s_sphPrtPos[pIndx*4+1], s_sphPrtPos[pIndx*4+2] ), mass, pIndx);
	m_iceSM[pIndx]->SetAlphas(pNum, 1.0);
	m_iceSM[pIndx]->SetBetas (pNum, 0.0);
	m_iceSM[pIndx]->SetLayer (pNum, 0);

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
			int sphIndx = jpIndx*4;

			m_iceSM[pIndx]->AddVertex( Vec3(s_sphPrtPos[sphIndx+0], s_sphPrtPos[sphIndx+1], s_sphPrtPos[sphIndx+2] ), mass, jpIndx);
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

void IceObject::SetClusterMoveInfoFromNeight(int pIndx, const vector<vector<rxNeigh>>& neights)
{
	//�ŏ��̈�ڂ͕K��pIndx�Ƃ���
	int pNum = m_iceSM[pIndx]->GetNumVertices();
	float mass = 1.0f;

	m_iceSM[pIndx]->AddVertex( Vec3(s_sphPrtPos[pIndx*4+0], s_sphPrtPos[pIndx*4+1], s_sphPrtPos[pIndx*4+2]), mass, pIndx);
	m_iceSM[pIndx]->SetAlphas(pNum, 1.0);
	m_iceSM[pIndx]->SetBetas (pNum, 0.0);
	m_iceSM[pIndx]->SetLayer (pNum, 0);

	//�ߖT���q���N���X�^�ɒǉ�
	for(unsigned id_np = 0; id_np < neights[pIndx].size(); id_np++)
	{
		int np_pIndx = neights[pIndx][id_np].Idx;
		int pNum = m_iceSM[pIndx]->GetNumVertices();
		double mass = 1.0;

		if(m_iceSM[pIndx]->CheckIndx(np_pIndx)){	continue;	}

		m_iceSM[pIndx]->AddVertex( Vec3(s_sphPrtPos[np_pIndx*4+0], s_sphPrtPos[np_pIndx*4+1], s_sphPrtPos[np_pIndx*4+2]), mass, np_pIndx);
		m_iceSM[pIndx]->SetAlphas(pNum, 1.0);
		m_iceSM[pIndx]->SetBetas (pNum, 0.0);
		m_iceSM[pIndx]->SetLayer (pNum, 0);
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

//�f�o�b�O���̕`��
void IceObject::Display()
{
	//�e�X�g
	//���݂�CalcMethod�̃��[�h�m�F
	const char* calcMethodName1 = typeid(Ice_CalcMethod_Itr_Stiffness).name();
	const char* calcMethodName2 = typeid(Ice_CalcMethod_Itr_Exp_Stiff).name();
	const char* nowName = typeid(*(m_iceCalcMethod)).name();

	//����iteration_stiffness���[�h���̊m�F�@strcmp�͕����񂪈�v�����0�ɂȂ�
	if(strcmp(calcMethodName1, nowName) != 0
	&& strcmp(calcMethodName2, nowName) != 0) return ;

	//rigid�I�u�W�F�N�g�̈ʒu��`�悵�Ă݂�
	if(strcmp(calcMethodName1, nowName) == 0){
		((Ice_CalcMethod_Itr_Stiffness*)m_iceCalcMethod)->DebugStiffness();
	}
	else if(strcmp(calcMethodName2, nowName) == 0){
		((Ice_CalcMethod_Itr_Exp_Stiff*)m_iceCalcMethod)->DebugStiffness();
	}
}

//�ő̂̉^���v�Z�C�v�Z���ʂ̓���
//�������e�͕ʃN���X�ɋL�q
void IceObject::StepObjMoveCPU()
{
	//CPUNormal��CPUDebug
	(this->*m_fpStepObjMove)();
}

void IceObject::StepObjMoveCPUNormal()
{
	//�ő̂̉^���v�Z
	m_iceCalcMethod->StepObjMove();			RXTIMER("StepObjMove");

	//�ő̂̍ŏI�ʒu����
	m_iceConvolution->StepConvolution();	RXTIMER("StepConvolution");

	//�t�̂ƌő̂̐��`���
	StepInterPolation();					RXTIMER("StepInterPolation");
	//StepInterPolationKenjya();
}

void IceObject::StepObjMoveCPUDebug()
{
	//�ő̂̉^���v�Z
	m_iceCalcMethod->StepObjMoveDebug();

	//�ő̂̍ŏI�ʒu����
	m_iceConvolution->StepConvolutionDebug();

	//�t�̂ƌő̂̐��`���
	StepInterPolation();

	//�f�[�^�쐬
	//TODO::��ɃN���X������
	DebugDeformationAmount();
	DebugDeformationAverage();
}

//�t�̂ƌő̂̐��`���
void IceObject::StepInterPolation()
{
	for(int pIndx = 0; pIndx < sm_particleNum; pIndx++)
	{
		if(m_iceSM[pIndx]->GetNumVertices() == 0){			continue;	}

		int sphIndx = pIndx*4;
		int sldIndx = pIndx*3;
		double intrps = 1.0 - m_fInterPolationCoefficience[pIndx];	//��ԌW��
	
		float* sldVel = Ice_SM::GetSldVelPointer();
		float* sldPos = Ice_SM::GetSldPosPointer();
	
		for(int i = 0; i < 3; i++)
		{
			s_sphPrtVel[sphIndx+i] = sldVel[sldIndx+i] * m_fInterPolationCoefficience[pIndx] + s_sphPrtVel[sphIndx+i] * intrps;
			s_sphPrtPos[sphIndx+i] = sldPos[sldIndx+i] * m_fInterPolationCoefficience[pIndx] + s_sphPrtPos[sphIndx+i] * intrps;
		}
	}
}

void IceObject::StepInterPolationKenjya()
{
	//���ό`��
	float defSum = 0.0f;
	for(int pIndx = 0; pIndx < m_iceSM.size(); pIndx++){
		defSum += m_iceSM[pIndx]->GetDefAmount();
	}

	for(int pIndx = 0; pIndx < sm_particleNum; pIndx++)
	{
		if(m_iceSM[pIndx]->GetNumVertices() == 0){			continue;	}

		int sphIndx = pIndx*4;
		int sldIndx = pIndx*3;

		float defAmount = m_iceSM[pIndx]->GetDefAmount() / defSum;
		if(defAmount > 1.0f) defAmount = 1.0f;
		if(defAmount < 0.1f) defAmount = 0.1f;
		double intrps = 1.0 - defAmount;	//��ԌW��
	
		float* sldVel = Ice_SM::GetSldVelPointer();
		float* sldPos = Ice_SM::GetSldPosPointer();
	
		for(int i = 0; i < 3; i++)
		{
			s_sphPrtVel[sphIndx+i] = sldVel[sldIndx+i] * defAmount + s_sphPrtVel[sphIndx+i] * intrps;
			s_sphPrtPos[sphIndx+i] = sldPos[sldIndx+i] * defAmount + s_sphPrtPos[sphIndx+i] * intrps;
		}
	}
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

#endif
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
#endif
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
	
	//�^���v�Z����N���X�^�̍Ē�`
	if(prtList.size() == 0) return;

QueryCounter counter1;
QueryCounter counter2;

counter1.Start();
	vector<unsigned> neighborList;
	//CashNeighborList(prtList, neighborList);
	ReConstructCluster(prtList, clusterList);	//RXTIMER("ReConstructCluster");	//�N���X�^�̍č\�z
double end1 = counter1.End();

counter2.Start();
	m_iceStrct->UpdateSelectCluster(prtList, neighborList, m_iceSM);
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
		for(unsigned oIndx = 0; oIndx < GetMoveObj(pIndx)->GetIndxNum(); oIndx++)
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

#ifndef USE_NEIGHT
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
#endif
				moveObj->CalcOrgCm();
				moveObj->ClassifyAllOrgParticle();
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
	for(int cIndx = 0; cIndx < sm_particleNum; cIndx++)
	{
		m_iceStrct->UpdateMotionCalcCluster(cIndx, 1);		
	}
}

void IceObject::CashNeighborList(const vector<unsigned>& prtList, vector<unsigned>& neighborList)
{
	for(unsigned lIndx = 0; lIndx < prtList.size(); lIndx++)
	{
		unsigned cIndx = prtList[lIndx];

		for(int indx = 0; indx < (int)m_iceSM[cIndx]->GetIndxNum(); indx++)
		{
			int pIndx = m_iceSM[cIndx]->GetParticleIndx(indx);
			if(pIndx == MAXINT) continue;

			if(find(neighborList.begin(), neighborList.end(), (unsigned)pIndx) != neighborList.end()) continue;

			neighborList.push_back((unsigned)pIndx);
		}
	}
}

//���ʏC���@�S��1.0f
void IceObject::UpdateParticleMass_Normal()
{
	unsigned pNum = IceObject::GetParticleNum();
	float mass = 1.0f;

	//�e���q�̊e�N���X�^�ł̎��ʂ��C��
	for(int cIndx = 0; cIndx < pNum; cIndx++)
	{
		for(int oIndx = 0; oIndx < m_iceSM[cIndx]->GetIndxNum(); oIndx++)
		{
			int pIndx = m_iceSM[cIndx]->GetParticleIndx(oIndx);
			if(pIndx == MAXINT){	continue;	}

			m_iceSM[cIndx]->SetMass(oIndx, mass);
		}
	}

	//�e�N���X�^�̏d�S���X�V
	for(int pIndx = 0; pIndx < sm_particleNum; pIndx++)
	{
		m_iceSM[pIndx]->CalcOrgCm();
		m_iceSM[pIndx]->ClassifyAllOrgParticle();
	}
}

//���ʏC���@�N���X�^���ŕ���
void IceObject::UpdateParticleMass_Average()
{
	unsigned pNum = IceObject::GetParticleNum();
	vector<unsigned> addParticleNum(pNum, 0);

	//�e���q�������̃N���X�^�Ɋ܂܂�Ă��邩�C�̃J�E���g
	for(int cIndx = 0; cIndx < pNum; cIndx++)
	{
		for(int oIndx = 0; oIndx < m_iceSM[cIndx]->GetIndxNum(); oIndx++)
		{
			int pIndx = m_iceSM[cIndx]->GetParticleIndx(oIndx);
			if(pIndx == MAXINT){	continue;	}

			//���q���̃J�E���g
			addParticleNum[pIndx] += 1;
		}
	}

	float mass = 1.0f;

	//�e���q�̊e�N���X�^�ł̎��ʂ��C��
	for(int cIndx = 0; cIndx < pNum; cIndx++)
	{
		for(int oIndx = 0; oIndx < m_iceSM[cIndx]->GetIndxNum(); oIndx++)
		{
			int pIndx = m_iceSM[cIndx]->GetParticleIndx(oIndx);
			if(pIndx == MAXINT){	continue;	}

			m_iceSM[cIndx]->SetMass(oIndx, mass/(float)addParticleNum[pIndx]);
		}
	}

	//�e�N���X�^�̏d�S���X�V
	for(int pIndx = 0; pIndx < sm_particleNum; pIndx++)
	{
		m_iceSM[pIndx]->CalcOrgCm();
		m_iceSM[pIndx]->ClassifyAllOrgParticle();
	}
}

//���ʏC���@�����x�N�g���Ƃ̋ߎ��x��
void IceObject::UpdateParticleMass_Direction()
{
	//�e�N���X�^�̏d�S���X�V�@����1.0f�ŏd�S���v�Z
	for(int pIndx = 0; pIndx < sm_particleNum; pIndx++){	m_iceSM[pIndx]->CalcOrgCm();	}

	unsigned pNum = IceObject::GetParticleNum();
	vector<unsigned> addParticleNum(pNum, 0);

	//�e���q�������̃N���X�^�Ɋ܂܂�Ă��邩�C�̃J�E���g
	for(int cIndx = 0; cIndx < pNum; cIndx++)
	{
		for(int oIndx = 0; oIndx < m_iceSM[cIndx]->GetIndxNum(); oIndx++)
		{
			int pIndx = m_iceSM[cIndx]->GetParticleIndx(oIndx);
			if(pIndx == MAXINT){	continue;	}

			//���q���̃J�E���g
			addParticleNum[pIndx] += 1;
		}
	}

	Vec3 dirVec = Vec3(0.0f, -1.0f, 1.0f);
	dirVec = Unit(dirVec);

	//�e���q�̊e�N���X�^�ł̎��ʂ��C��
	for(int cIndx = 0; cIndx < pNum; cIndx++)
	{
		Vec3 orgCm = m_iceSM[cIndx]->GetOrgCm();

		for(int oIndx = 0; oIndx < m_iceSM[cIndx]->GetIndxNum(); oIndx++)
		{
			int pIndx = m_iceSM[cIndx]->GetParticleIndx(oIndx);
			if(pIndx == MAXINT){	continue;	}

			//ijiri��2012�̎�@
			float mass = 1.0f;
			Vec3 prtPos = m_iceSM[cIndx]->GetOrgPos(oIndx);
			Vec3 prtVec = Unit(prtPos-orgCm);
			float weight = mass * (dot(dirVec, prtVec) * dot(dirVec, prtVec) + 0.01f);

			m_iceSM[cIndx]->SetMass(oIndx, weight);
		}
	}

	//�e�N���X�^�̏d�S���X�V�@�C�����ʂŏd�S���v�Z
	for(int pIndx = 0; pIndx < sm_particleNum; pIndx++){	m_iceSM[pIndx]->CalcOrgCm();	}
}

//�����ʒu��ۑ�
void IceObject::SaveInitPos()
{
	for(int pIndx = 0; pIndx < sm_particleNum; pIndx++)
	{
		int sphIndx = pIndx*4;
		m_VecInitPos.push_back(Vec3(s_sphPrtPos[sphIndx+0], s_sphPrtPos[sphIndx+1], s_sphPrtPos[sphIndx+2]));
	}
}

//�ʒu�Ƒ��x������������
void IceObject::ResetPosAndVel()
{
//���q
	for(int pIndx = 0; pIndx < sm_particleNum; pIndx++)
	{
		int sphIndx = pIndx*4;

		s_sphPrtPos[sphIndx+0] = m_VecInitPos[pIndx][0];
		s_sphPrtPos[sphIndx+1] = m_VecInitPos[pIndx][1];
		s_sphPrtPos[sphIndx+2] = m_VecInitPos[pIndx][2];

		s_sphPrtVel[sphIndx+0] = 0.0f;
		s_sphPrtVel[sphIndx+1] = 0.0f;
		s_sphPrtVel[sphIndx+2] = 0.0f;
	}

//�N���X�^
	for(int pIndx = 0; pIndx < sm_particleNum; pIndx++)
	{
		int smIndx = pIndx*3;
		float* p = Ice_SM::GetSldPosPointer();
		float* v = Ice_SM::GetSldVelPointer();
	
		v[smIndx+0] = 0.0f;
		v[smIndx+1] = 0.0f;
		v[smIndx+2] = 0.0f;
	
		p[smIndx+0] = m_VecInitPos[pIndx][0];
		p[smIndx+1] = m_VecInitPos[pIndx][1];
		p[smIndx+2] = m_VecInitPos[pIndx][2];
	}
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

//�ߖT���q
void IceObject::DebugNeights(const vector<vector<rxNeigh>>& neights)
{	cout << __FUNCTION__ << endl;
	
	string result = RESULT_DATA_PATH;
	result += "DebugNeights.txt";
	ofstream ofs(result, ios::app);

	for(unsigned pIndx = 0; pIndx < neights.size(); pIndx++){
		if(neights[pIndx].size() == 0){
			continue;
		}
		ofs << "Particle:" << pIndx << endl;
		for(unsigned id_np = 0; id_np < neights[pIndx].size(); id_np++){
			int np_pIndx = neights[pIndx][id_np].Idx;
			ofs << "	" << np_pIndx << endl;
		}
	}
}

//����臒l�𒴂������q���J�E���g
void IceObject::DebugDeformationAmount()
{	//cout << __FUNCTION__ << endl;

	int deformClusterNum = 0;
	for(int i = 0; i < sm_clusterNum; ++i)
	{	
		if(m_iceSM[i]->GetNumVertices() == 0){	continue;	}

		float defAmount = m_iceSM[i]->GetDefAmount();
		if(defAmount < 1.0f){	continue;	}

		deformClusterNum++;
	}

	string result = RESULT_DATA_PATH;
	result += "DeformationClusterNum.txt";
	ofstream ofs(result, ios::app);
	ofs << deformClusterNum << endl;
}

//���ϕό`��
void IceObject::DebugDeformationAverage()
{
	float deformClusterAverage = 0.0f;
	for(int i = 0; i < sm_clusterNum; ++i)
	{	
		if(m_iceSM[i]->GetNumVertices() == 0){	continue;	}

		float defAmount = m_iceSM[i]->GetDefAmount();

		deformClusterAverage += defAmount;
	}

	string result = RESULT_DATA_PATH;
	result += "DeformationAverage.txt";
	ofstream ofs(result, ios::app);
	ofs << deformClusterAverage/(float)sm_clusterNum << endl;
}

//�f�B���N�g���������ׂč폜
void IceObject::DebugClearDirectry(string path)
{	cout << __FUNCTION__ << " path = " << path << endl;

	const boost::filesystem::path dir(path);
	boost::filesystem::directory_iterator end;

	//�C�e���[�^���ō폜�̎d�����킩��Ȃ������̂ŁC�S�t�@�C�����擾���Ă���폜���Ă���
	vector<boost::filesystem::path> files;

	//�f�B���N�g�����̃t�@�C���擾
	for(boost::filesystem::directory_iterator p(dir); p != end; ++p)
	{
		if(boost::filesystem::is_directory(*p))
		{  
			cout << "	" << p->path() << " is directry" << std::endl;	//�f�B���N�g���̎�
		}
		else
		{
			cout << "	" << p->path() << " is file" <<  endl;
		}

		files.push_back(p->path());
	}

	cout << endl;

	//�t�@�C���폜
	vector<boost::filesystem::path>::iterator deleteFile;
	for(deleteFile = files.begin(); deleteFile != files.end(); deleteFile++)
	{
		if(boost::filesystem::is_directory(*deleteFile))
		{
			cout << "	" <<  *deleteFile << " �̓f�B���N�g���ł�" << endl;
			continue;
		}

		if(boost::filesystem::remove(*deleteFile))
		{
			cout << "	" << *deleteFile << " ���폜���܂���" << endl;
		}
		else
		{
			cout << "	" << "�G���[ " <<  *deleteFile << " ���폜�ł��܂���" << endl;
		}
	}
}

//�^���v�Z���Ă��Ȃ��N���X�^�̕ό`�ʂ��C�O�t���[���Ƃ̍������狁�߂�
//�����ʂ��s���S�������̂Ŏg��Ȃ�
//�������ǁC�S�ẴN���X�^���^���v�Z�����Ă���C���ʂ̎Z�o�Ɏg���g��Ȃ��𔻒f����C�Ƃ������@�ɂ���
void IceObject::UpdateUnSelectedClusterDefAmount()
{
	float* sldPos = Ice_SM::GetSldPosPointer();

	for(int cIndx = 0; cIndx < sm_particleNum; cIndx++)
	{
		if(m_iceStrct->GetMotionCalcCluster(cIndx) != 0){	continue;	}
	
		float defAmount = 0.0f;
	
		//���݂̏d�S���v�Z
		Vec3 nowCm(0.0);
		for(int i = 0; i < m_iceSM[cIndx]->GetIndxNum(); ++i)
		{
			int pIndx = m_iceSM[cIndx]->GetParticleIndx(i);			
			if(pIndx == MAXINT){	continue;	}

			nowCm += Vec3(sldPos[pIndx*3+0], sldPos[pIndx*3+1], sldPos[pIndx*3+2]);
		}

		nowCm /= (float)m_iceSM[cIndx]->GetNumVertices();

		//�ό`�ʂ����߂�
		for(int i = 0; i < m_iceSM[cIndx]->GetIndxNum(); ++i)
		{
			int pIndx = m_iceSM[cIndx]->GetParticleIndx(i);			
			if(pIndx == MAXINT){	continue;	}
	
			//�O�t���[���̑��΍��W
			Vec3 prePos = m_iceSM[cIndx]->GetPrePos(i) - m_iceSM[cIndx]->GetPreCm();

			//���݂̑��΍��W
			Vec3 nowPos = Vec3(sldPos[pIndx*3+0], sldPos[pIndx*3+1], sldPos[pIndx*3+2]) - nowCm;

			defAmount += abs(nowPos[0]-prePos[0]);
			defAmount += abs(nowPos[1]-prePos[1]);
			defAmount += abs(nowPos[2]-prePos[2]);

			//�O�t���[���̈ʒu�x�N�g�����X�V
			m_iceSM[cIndx]->SetPrePos(i, Vec3(sldPos[pIndx*3+0], sldPos[pIndx*3+1], sldPos[pIndx*3+2]));
		}
	
		m_iceSM[cIndx]->SetDefAmount(defAmount);
		m_iceSM[cIndx]->SetPreCm(nowCm);

		cout << cIndx << ":" << defAmount << endl;
	}
}

void IceObject::MakeDirectory(string pathName)
{	cout << __FUNCTION__ << " " << pathName << endl;

	const boost::filesystem::path path(pathName);

    boost::system::error_code error;
    const bool result = boost::filesystem::create_directory(path, error);
    if (!result || error) {
        std::cout << "�f�B���N�g���̍쐬�Ɏ��s" << std::endl;
    }
}

//�t�@�C�����Ɏg���镶������C���ݎ�������쐬
string IceObject::MakeDirNameTimeStamp()
{	
	FUNCTION_NAME;

	//���݂̓���
	struct tm;
	unsigned int size = 26;
	char timeStr[80];
	time_t localtime;
	time(&localtime);
	ctime_s(timeStr, size, &localtime);
	
	//�f�B���N�g�������쐬
	string dirName(timeStr);
	replace(dirName.begin(), dirName.end(), ' ', '_');	//�X�y�[�X���A���_�[�X�R�A�ɒu�� string�̊֐��łȂ��̂ɒ���
	replace(dirName.begin(), dirName.end(), ':', '_');	//�R�������A���_�[�X�R�A�ɒu��
	struct pred { bool operator()(char c) const { return /*c == ' ' ||*/ c == '\t' || c == '\r' || c == '\n'; } };
	dirName.erase(std::remove_if(dirName.begin(), dirName.end(), pred()), dirName.end());	//�^�u�C���^�[���C���s���폜

	return dirName;
}

void IceObject::DebugMakeGraph()
{	

	FUNCTION_NAME;

	//�f�o�b�O���[�h���̊m�F
	if(m_fpStepObjMove != &IceObject::StepObjMoveCPUDebug)	return;

	//���ݎ�������f�B���N�g�����ɗp���镶������쐬
	string dirName = MakeDirNameTimeStamp();
	dirName.insert(0, RESULT_IMAGE_PATH);			//�f�B���N�g�����쐬����p�X

	//�f�B���N�g���쐬
	MakeDirectory(dirName);

	//�ό`�N���X�^��
	DebugMakeGraph_DefClusterNum(dirName);

	//���ϕό`��
	DebugMakeGraph_DefAmountAverage(dirName);

	//�����񐔁@TODO::���ꂢ�ȕ��@����Ȃ��c
	DebugMakeGraph_IterationNum(dirName);

	//���݂̃��[�h���e�L�X�g�ŕۑ�
	ofstream ofs(dirName+"/nowMode.txt");

	string nowMode = DebugNowMoveMethod();
	ofs << dirName << endl;
	ofs << nowMode << endl;

////�T���v��
//	vector<double> x,y;
//	x.resize(100);
//	y.resize(100);
//	//double x[100], y[100];
//
//	for (int i = 0; i<100; i++) {
//		x[i] = 0.0628*i;
//		y[i] = sin(x[i]);
//	}
//
//	vector<double> x2,y2;
//	x2.resize(100);
//	y2.resize(100);
//	//double x2[100], y2[100];
//
//	for (int i = 0; i<100; i++) {
//		x2[i] = 0.0628*i;
//		y2[i] = -sin(x[i]);
//	}
//
//	CGnuplot gp;
//	gp.SetTitle("Trigonometric Functions");
//	gp.Command("set key left bottom");
//	gp.SetLabel("x", "y");
//	gp.SetPlotType(CGnuplot::PLOT_TYPE_LINES_POINTS);
//	gp.SetXRange(-1.0, 1.0);
//	gp.SetYRange(-1.0, 1.0);
//
//	MultiPlot plotData;
//	plotData.push_back(svv("sin(x)", x, y));
//	plotData.push_back(svv("cos(x)", x2, y2));
//	gp.Multiplot(plotData);
//
//	gp.DumpToPng("result/images/test");
//	//���܂��t�@�C�������Ȃ��@�Ō�̍s�����������Ȃ�
//	gp.DumpToFile("result/data/test.txt");

	////�����ɏI������ƃt�@�C��������Ȃ������̂ł�����Ƒ҂�
	//int sleepTime = 3;	//�b
	//Sleep(sleepTime*1000);
}

//�ό`�N���X�^��
//�t�@�C����ǂݍ���Ńf�[�^�쐬
//TODO::������̓t�@�C���𒼐ړn���ăO���t�����
void IceObject::DebugMakeGraph_DefClusterNum(string dirName)
{
	FUNCTION_NAME;

	vector<int> defNum_x, defNum_y;

	//�f�[�^�t�@�C�����J��
	string filePath = RESULT_DATA_PATH;
	filePath += "DeformationClusterNum.txt";
	ifstream ifs(filePath);

	//�t�@�C���̑��݊m�F
	if(ifs.fail()) {	cerr << filePath << " is do not exist.\n";	return;	}

	//�ϐ��̗p�ӁC������
	int count = 0;
	string str;

	//������̓ǂݍ���
	while( getline(ifs, str) )
	{
		int data = 0;

		sscanf(str.data(), "%d", &data);

		defNum_x.push_back(count++);
		defNum_y.push_back(data);
	}

	ifs.close();

	//�g�p�����f�[�^�t�@�C�����ړ�
	string newFileName = dirName;
	newFileName += "/DeformationClusterNum.txt";
	const boost::filesystem::path path(filePath);
	const boost::filesystem::path dest(newFileName);

    try {
        boost::filesystem::rename(path, dest);
    }
    catch (boost::filesystem::filesystem_error& ex) {
        std::cout << ex.what() << std::endl;
        throw;
    }

	//�O���t�쐬
	CGnuplot gp;
	gp.SetTitle("Deformation Cluster Num Graph");
	gp.Command("set key left bottom");
	gp.SetLabel("TimeStep", "ClusterNum");
	gp.SetPlotType(CGnuplot::PLOT_TYPE_LINES_POINTS);
	gp.SetXRange(0, (int)defNum_x.size());
	gp.SetYRange(0, sm_particleNum);
	gp.Plot(defNum_x, defNum_y);

	string defNumResultPath = dirName;
	defNumResultPath += "/DeformationClusterNum";
	gp.DumpToPng(defNumResultPath.c_str());

	//�����ɏI������ƃt�@�C��������Ȃ������̂ł�����Ƒ҂�
	int sleepTime = 2;	//�b
	Sleep(sleepTime*1000);
}

//���ϕό`�ʃN���X�^��
void IceObject::DebugMakeGraph_DefAmountAverage(string dirName)
{
	FUNCTION_NAME;

	vector<float> defNum_x;
	vector<float> defNum_y;

	//�t�@�C����ǂݍ���
	string filePath = RESULT_DATA_PATH;
	filePath += "DeformationAverage.txt";
	ifstream ifs(filePath);

	//�t�@�C���̑��݊m�F
	if(ifs.fail()) {	cerr << filePath << " is do not exist.\n";	return;	}

	//�ϐ��̗p�ӁC������
	int count = 0;
	float max = 0;
	string str;

	//������̓ǂݍ���
	while( getline(ifs, str) )
	{
		float data = 0;

		sscanf(str.data(), "%f", &data);

		if(data > max){	max = data;	}

		defNum_x.push_back(count++);
		defNum_y.push_back(data);
	}

	ifs.close();

	//�g�p�����f�[�^�t�@�C�����ړ�
	string newFileName = dirName;
	newFileName += "/DeformationAverage.txt";
	const boost::filesystem::path path(filePath);
	const boost::filesystem::path dest(newFileName);

    try {
        boost::filesystem::rename(path, dest);
    }
    catch (boost::filesystem::filesystem_error& ex) {
        std::cout << ex.what() << std::endl;
        throw;
    }

	//�O���t�쐬
	CGnuplot gp;
	gp.SetTitle("Deformation Average Graph");
	gp.Command("set key left bottom");
	gp.SetLabel("TimeStep", "Average Deformation Amount");
	gp.SetPlotType(CGnuplot::PLOT_TYPE_LINES_POINTS);
	gp.SetXRange(0, (int)defNum_x.size());
	//gp.SetYRange(0, max*2.0);	//�c���͍ő�l�̓�{
	gp.SetYRange(0, 0.5);			//��r���₷���悤�ɌŒ�
	gp.Plot(defNum_x, defNum_y);

	string defNumResultPath = dirName;
	defNumResultPath += "/DefotmationAverage";
	gp.DumpToPng(defNumResultPath.c_str());

	//�����ɏI������ƃt�@�C��������Ȃ������̂ł�����Ƒ҂�
	int sleepTime = 2;	//�b
	Sleep(sleepTime*1000);
}

void IceObject::DebugMakeGraph_IterationNum(string dirName)
{
	//���[�h�m�F
	//TODO: ����܂��Y��Ȃ�������Ȃ�
	const char* calcMethodName1 = typeid(Ice_CalcMethod_Itr_Stiffness).name();
	const char* calcMethodName2 = typeid(Ice_CalcMethod_Itr_Exp_Stiff).name();
	const char* nowName = typeid(*(m_iceCalcMethod)).name();

	//����iteration_stiffness���[�h���̊m�F�@strcmp�͕����񂪈�v�����0�ɂȂ�
	if(strcmp(calcMethodName1, nowName) != 0
	&& strcmp(calcMethodName2, nowName) != 0) return ;

	FUNCTION_NAME;

	vector<float> defNum_x;
	vector<float> defNum_y;

	//�t�@�C����ǂݍ���
	string filePath = RESULT_DATA_PATH;
	filePath += "Itr_Stiffness_CountNum.txt";
	ifstream ifs(filePath);

	//�t�@�C���̑��݊m�F
	if(ifs.fail()) {	cerr << filePath << " is do not exist.\n";	return;	}

	//�ϐ��̗p�ӁC������
	int count = 0;
	float max = 0;
	string str;

	//������̓ǂݍ���
	while( getline(ifs, str) )
	{
		float data = 0;

		sscanf(str.data(), "%f", &data);

		if(data > max){	max = data;	}

		defNum_x.push_back(count++);
		defNum_y.push_back(data);
	}

	ifs.close();

	//�g�p�����f�[�^�t�@�C�����ړ�
	string newFileName = dirName;
	newFileName += "/Itr_Stiffness_CountNum.txt";
	const boost::filesystem::path path(filePath);
	const boost::filesystem::path dest(newFileName);

    try {
        boost::filesystem::rename(path, dest);
    }
    catch (boost::filesystem::filesystem_error& ex) {
        std::cout << ex.what() << std::endl;
        throw;
    }

	//�O���t�쐬
	CGnuplot gp;
	gp.SetTitle("Iteration Num Graph");
	gp.Command("set key left bottom");
	gp.SetLabel("TimeStep", "Iteration Num");
	gp.SetPlotType(CGnuplot::PLOT_TYPE_LINES_POINTS);
	gp.SetXRange(0, (int)defNum_x.size());
	//gp.SetYRange(0, max*2.0);	//�c���͍ő�l�̓�{
	gp.SetYRange(0, 250);			//��r���₷���悤�ɌŒ�
	gp.Plot(defNum_x, defNum_y);

	string defNumResultPath = dirName;
	defNumResultPath += "/Itr_Stiffness_CountNum";
	gp.DumpToPng(defNumResultPath.c_str());

	//�����ɏI������ƃt�@�C��������Ȃ������̂ł�����Ƒ҂�
	int sleepTime = 2;	//�b
	Sleep(sleepTime*1000);
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
string IceObject::DebugNowMoveMethod()
{	cout << __FUNCTION__ << endl;

	string nowMode;
	string end = "\n";
	nowMode = string("CalcMethod	") +	typeid(*m_iceCalcMethod).name() +	end		//�v�Z���@
			+ string("JudgeMove	") +		typeid(*m_iceJudeMove).name() +		end		//�^���v�Z�Ώۂ𔻒�
			+ string("ClsuterMove	") +	typeid(*m_iceClsuterMove).name() +	end		//�^���v�Z���@
			+ string("ConvopoJudge	") +	typeid(*m_iceConvoJudge).name() +	end		//�ŏI�������ʂɗp����Ώۂ𔻒�
			+ string("Convolution	" ) +	typeid(*m_iceConvolution).name() +	end;	//�ŏI�������ʂ����߂�

	cout << nowMode << endl;

	return nowMode;
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

//�ő̂̏㕔���Œ�
void IceObject::TestFixUpperPos()
{
	int side = 17;

	//���q
	for(int pIndx = 0; pIndx < sm_particleNum; pIndx++)
	{
		if(pIndx%(side*side) < side * (side-1)) continue;
		//if( pIndx%(side*side) % side >= 10) continue;

		int sphIndx = pIndx*4;

		s_sphPrtPos[sphIndx+0] = m_VecInitPos[pIndx][0];
		s_sphPrtPos[sphIndx+1] = m_VecInitPos[pIndx][1];
		s_sphPrtPos[sphIndx+2] = m_VecInitPos[pIndx][2];

		s_sphPrtVel[sphIndx+0] = 0.0f;
		s_sphPrtVel[sphIndx+1] = 0.0f;
		s_sphPrtVel[sphIndx+2] = 0.0f;
	}

//�N���X�^
	for(int pIndx = 0; pIndx < sm_particleNum; pIndx++)
	{
		if(pIndx%(side*side) < side * (side-1)) continue;
		//if( pIndx%(side*side) % side >= 10) continue;

		int smIndx = pIndx*3;
		float* p = Ice_SM::GetSldPosPointer();
		float* v = Ice_SM::GetSldVelPointer();
	
		v[smIndx+0] = 0.0f;
		v[smIndx+1] = 0.0f;
		v[smIndx+2] = 0.0f;
	
		p[smIndx+0] = m_VecInitPos[pIndx][0];
		p[smIndx+1] = m_VecInitPos[pIndx][1];
		p[smIndx+2] = m_VecInitPos[pIndx][2];
	}
}

//�ő̂̑��ʂ��Œ�
void IceObject::TestFixSidePos()
{
	int side = 17;

	//���q
	for(int pIndx = 0; pIndx < sm_particleNum/side; pIndx++)
	{
		int sphIndx = pIndx*4;

		s_sphPrtPos[sphIndx+0] = m_VecInitPos[pIndx][0];
		s_sphPrtPos[sphIndx+1] = m_VecInitPos[pIndx][1];
		s_sphPrtPos[sphIndx+2] = m_VecInitPos[pIndx][2];

		s_sphPrtVel[sphIndx+0] = 0.0f;
		s_sphPrtVel[sphIndx+1] = 0.0f;
		s_sphPrtVel[sphIndx+2] = 0.0f;
	}

//�N���X�^
	for(int pIndx = 0; pIndx < sm_particleNum/side; pIndx++)
	{
		int smIndx = pIndx*3;
		float* p = Ice_SM::GetSldPosPointer();
		float* v = Ice_SM::GetSldVelPointer();
	
		v[smIndx+0] = 0.0f;
		v[smIndx+1] = 0.0f;
		v[smIndx+2] = 0.0f;
	
		p[smIndx+0] = m_VecInitPos[pIndx][0];
		p[smIndx+1] = m_VecInitPos[pIndx][1];
		p[smIndx+2] = m_VecInitPos[pIndx][2];
	}
}

//�t�@�C����ǂݍ��񂾃p�����[�^�ŃV�~�����[�V���������s
void IceObject::TestSimulationFromFile(string fileName)
{
	cout << __FUNCTION__ << " " << fileName << endl;

	//�t�@�C���I�[�v��
	//�f�[�^�ǂݍ���
	//�f�[�^���f
}