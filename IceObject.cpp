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
	for(vector<Ice_SM*>::iterator it = m_iceMove.begin(); it != m_iceMove.end(); ++it)
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
		m_iceMove.push_back(new Ice_SM(sm_clusterNum));
		m_iceMove[sm_clusterNum]->SetSimulationSpace(boundarySpaceLow, boundarySpaceHigh);
		m_iceMove[sm_clusterNum]->SetTimeStep(timeStep);
		m_iceMove[sm_clusterNum]->SetCollisionFunc(0);
		m_iceMove[sm_clusterNum]->SetStiffness(1.0, 1.0);

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
		int pNum = m_iceMove[i]->GetNumVertices();
		vector<int> pList;

		for(int j = 0; j < pNum; j++)
		{
			pList.push_back(m_iceMove[i]->GetParticleIndx(j));
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

//�������ɗp����p�X�̏�����
void IceObject::InitPath()
{
	m_SurfSm.InitPath(s_sphPrtPos, s_sphPrtVel, m_iceMove, m_iceStrct, sm_clusterNum, sm_particleNum);
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
			if(m_iceMove[pIndx]->CheckIndx(jpIndx)){	continue;	}

			int pNum = m_iceMove[pIndx]->GetNumVertices();
			float mass = 1.0f;

			m_iceMove[pIndx]->AddVertex( Vec3(s_sphPrtPos[jpIndx*4+0], s_sphPrtPos[jpIndx*4+1], s_sphPrtPos[jpIndx*4+2] ), mass, jpIndx);
			m_iceMove[pIndx]->SetAlphas(pNum, 1.0);
			m_iceMove[pIndx]->SetBetas (pNum, 0.0);
			m_iceMove[pIndx]->SetLayer (pNum, 0);
		}
	}

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
				if(m_iceMove[pIndx]->CheckIndx(kpIndx)){	/*cout << "contineue kpIndx = " << kpIndx << endl;*/ continue;	}

				int pNum = m_iceMove[pIndx]->GetNumVertices();
				float mass = 1.0f;

				m_iceMove[pIndx]->AddVertex( Vec3(s_sphPrtPos[kpIndx*4+0], s_sphPrtPos[kpIndx*4+1], s_sphPrtPos[kpIndx*4+2] ), mass, kpIndx);
				m_iceMove[pIndx]->SetAlphas(pNum, 1.0);
				m_iceMove[pIndx]->SetBetas (pNum, 0.0);
				m_iceMove[pIndx]->SetLayer (pNum, jlIndx);
				//cout << "pIndx = " << pIndx << " GetNumVertices = " << m_iceMove[pIndx]->GetNumVertices() << endl;
			}
		}
	}
}

//���q���N���X�^��N���X�^�����q��o�^
void IceObject::SetClusterStrctInfo(int cIndx, int* PtoCNum)
{
	//���q�������Ă���l�ʑ̂̔ԍ���o�^���邽�߂̏���
	//pCountList�ɂ́CcIndx�Ԗڂ̃N���X�^�Ɋ܂܂��e���q���C���ꂼ�ꂢ���̃N���X�^�ɑ����邩�����߂ĕۑ�����
	int* pCountList = new int[m_iceMove[cIndx]->GetNumVertices()];
	int* pLayerList = new int[m_iceMove[cIndx]->GetNumVertices()];			//���q�̏������C���[
	
	for(int i = 0; i < m_iceMove[cIndx]->GetNumVertices(); i++)
	{
		int pIndx = m_iceMove[cIndx]->GetParticleIndx(i);

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

		pLayerList[i] = m_iceMove[cIndx]->GetLayer(i);					//���q�����w�ڂ̋ߖT�Ȃ̂����擾
	}

	//���q�ƃN���X�^�̏��o�^
	vector<int> pIndxList;

	for(int i = 0; i < GetCtoPNum(cIndx); i++)
	{
		int pIndx = m_iceMove[cIndx]->GetParticleIndx(i);
			
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
	Ice_SM::InitGPU(m_iceMove, sd_sphPrtPos, sd_sphPrtVel, sm_particleNum, sm_maxParticleNum);

#if defined(USE_PATH)
{
	m_SurfSm.InitPathGPU();		//�������ɗp����p�X�̏�����
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

//�ő̂̉^���v�Z
void IceObject::StepObjMoveCPU()
{
	#pragma omp parallel
	{
	#pragma omp for
		for(int i = 0; i < sm_clusterNum; ++i)
		{	
			if(GetPtoCNum(i) == 0){	continue;	}
			m_iceMove[i]->UpdateCPU();				//�^���v�Z
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
			if(GetPtoCNum(i) == 0){	continue;	}

			if(i%SELECTED != 1){	continue;	}	//�e�X�g
			m_iceMove[i]->calExternalForces();					//���݈ʒu�̌v�Z
			m_iceMove[i]->ShapeMatchingSelected(SELECTED);		//���݈ʒu�̍X�V ���ʂ̌v�Z
			m_iceMove[i]->integrate(0.02f);						//���x�̌v�Z
		}
	}//#pragma omp parallel
}

//GPU
void IceObject::StepObjMoveGPU()
{
	Ice_SM::SetDevicePosPointer(sd_sphPrtPos);	//VBO���g���Ă���̂ŁC���ꂪ�Ȃ��ƃG���[���o��̂ɒ���
	Ice_SM::UpdateGPU();
}

//�p�X��p�����������@�ő̂̉^���v�Z
void IceObject::StepObjMoveUsePath()
{
	//prefixSum�̍X�V
	m_SurfSm.UpdatePrefixSum();

	//�N���X�^�̃p�����[�^�X�V
	#pragma omp parallel for
	for(int i = 0; i < sm_clusterNum; ++i)
	{	
		if(GetPtoCNum(i) == 0){	continue;	}
		m_iceMove[i]->SetNowCm(m_SurfSm.CalcCmSum(i));		//�d�S�̍X�V
		m_iceMove[i]->SetApq(m_SurfSm.CalcApqSum(i));		//�ό`�s��̍X�V
	}

	//�^���v�Z
	#pragma omp parallel for
	for(int i = 0; i < sm_clusterNum; ++i)
	{	
		if(GetPtoCNum(i) == 0){	continue; }
		m_iceMove[i]->UpdateUsePathCPU();
	}
}

void IceObject::StepObjMoveGPUUsePath()
{
	//prefixSum�̍X�V
	m_SurfSm.SetDevicePointer(sd_sphPrtPos, sd_sphPrtVel);	//VBO���g���Ă���̂ŁC���ꂪ�Ȃ��ƃG���[���o��̂ɒ���
	m_SurfSm.UpdatePrefixSumGPU();
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
		m_SurfSm.GetDevicePRTtoPTHPointer(),
		m_SurfSm.GetDevicePTHandPrfxSetPointer(),
		m_SurfSm.GetDecvicePrfxPos(),
		m_SurfSm.GetDecvicePrfxApq(),
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

//�p�X��p���������^���v�Z
void IceObject::StepObjMoveIterationUsePath()
{
	//�^���v�Z�𔽕�
	for(int itr = 0; itr < Ice_SM::GetIteration(); ++itr)
	{
		//����
		if(itr == 0)
		{
			m_SurfSm.UpdatePrefixSum();								//prefixSum�̍X�V
		}
		//����
		else
		{
			m_SurfSm.UpdatePrefixSumItr();							//prefixSum�̍X�V
		}

		//�e�N���X�^�̃f�[�^�X�V
		#pragma omp parallel for
		for(int i = 0; i < sm_clusterNum; ++i)
		{	
			if(GetPtoCNum(i) == 0){	continue;	}
			GetMoveObj(i)->SetNowCm(m_SurfSm.CalcCmSum(i));		//�d�S�̍X�V
			GetMoveObj(i)->SetApq(m_SurfSm.CalcApqSum(i));		//�ό`�s��̍X�V
		}

		//�^���v�Z
		#pragma omp parallel for
		for(int i = 0; i < sm_clusterNum; ++i)
		{	
			if(GetPtoCNum(i) == 0){	continue;	}
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

	StepInterPolation();
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





//���̂ƌő̂̍ŏI�I�ȉ^���v�Z
void IceObject::StepInterPolation()
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
		if(GetParticleNum() <= i){	continue;	}	//�Z���݂̂̎����̂Ƃ��ɕK�v�ɂȂ�D
		if(GetPtoCNum(i) <= 0){		continue;	}
	
		//�ő̉^���̍ŏI�ʒu�v�Z
		CalcAverageCPU(i, pos, vel);

		//�t�̂ƌő̂̕��
		LinerInterPolationCPU(i, pos, vel);			//���`���
	}

//�f�o�b�O
	//DebugDeformationAmount();
	DebugDeformationAverage();

#endif
}

//���q��I��I�ɉ^���v�Z������p�^�[���̃e�X�g
void IceObject::StepInterPolationSelected()
{
	Vec3 pos,vel;

	vector<Vec3> preSphPos;
	preSphPos.resize(sm_particleNum);
	for(int i = 0; i < sm_particleNum; i++)
	{
		Vec3 pos;
		pos[0] = s_sphPrtPos[i*4+0];
		pos[1] = s_sphPrtPos[i*4+1];
		pos[2] = s_sphPrtPos[i*4+2];

		preSphPos.push_back(pos);
	}

	#pragma omp parallel for private(pos, vel)
	for(int i = 0; i < sm_particleNum; ++i)
	{
		if(GetParticleNum() <= i){	continue;	}	//�Z���݂̂̎����̂Ƃ��ɕK�v�ɂȂ�D
		if(GetPtoCNum(i) <= 0){		continue;	}
	
		//�ő̉^���̍ŏI�ʒu�v�Z
		CalcAverageSelected(i, pos, vel);

		//�t�̂ƌő̂̕��
		LinerInterPolationCPU(i, pos, vel);			//���`���
	}

	#pragma omp parallel for private(pos, vel)
	for(int i = 0; i < sm_particleNum; ++i)
	{
		if(GetParticleNum() <= i){	continue;	}	//�Z���݂̂̎����̂Ƃ��ɕK�v�ɂȂ�D
		if(GetPtoCNum(i) <= 0){		continue;	}

		//�ό`�ʂ̔��f
		UpdateUnSelectedCluster(i, pos, vel, preSphPos);
	}


//�f�o�b�O
	//DebugDeformationAmount();
	DebugDeformationAverage();
}

void IceObject::UpdateUnSelectedCluster(int cIndx, const Vec3& pos, const Vec3& vel, const vector<Vec3>& preSphPos)
{
	if(cIndx%SELECTED == 1)	return;	//�^���v�Z���Ă��Ȃ��N���X�^���Ώ�

	float defAmount = 0.0f;

	for(int i = 0; i < m_iceMove[cIndx]->GetVertexNum(); ++i)
	{
		int pIndx = m_iceMove[cIndx]->GetParticleIndx(i);
		
		if(pIndx%SELECTED == 1)	continue;

		Vec3 nowPos = m_iceMove[cIndx]->GetVertexPos(i);
		Vec3 prePos = preSphPos[pIndx];

		defAmount += abs(nowPos[0]-prePos[0]);
		defAmount += abs(nowPos[1]-prePos[1]);
		defAmount += abs(nowPos[2]-prePos[2]);
	}

	m_iceMove[cIndx]->SetDefAmount(defAmount);
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
		if(GetParticleNum() <= i){	continue;	}		//�Z���݂̂̎����̂Ƃ��ɕK�v�ɂȂ�D
		if(GetPtoCNum(i) <= 0){		continue;	}

		//�ő̉^���̍ŏI�ʒu�v�Z
		CalcAverageCPU(i, pos, vel);

		//�t�̂ƌő̂̕��
		LinerInterPolationForClusterCPU(i, pos, vel);	//���`���
	}

#endif
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

		pos += m_iceMove[cIndx]->GetVertexPos(oIndx);
		vel += m_iceMove[cIndx]->GetVertexVel(oIndx);

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
		if(cIndx%SELECTED != 1){		continue;	}

		pos += m_iceMove[cIndx]->GetVertexPos(oIndx);
		vel += m_iceMove[cIndx]->GetVertexVel(oIndx);

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
		pos = Vec3(s_sphPrtPos[jpIndx+0], s_sphPrtPos[jpIndx+1], s_sphPrtPos[jpIndx+2]);
		vel = Vec3(s_sphPrtVel[jpIndx+0], s_sphPrtVel[jpIndx+1], s_sphPrtVel[jpIndx+2]);
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

		float defAmount = m_iceMove[cIndx]->GetDefAmount();
		//defAmount *= defAmount;				//����
		//defAmount *= defAmount * defAmount;	//�O���

		pos += m_iceMove[cIndx]->GetVertexPos(oIndx) * defAmount;
		vel += m_iceMove[cIndx]->GetVertexVel(oIndx) * defAmount;

		deformAmountSum += defAmount;
	}

	pos /= deformAmountSum;
	vel /= deformAmountSum;
	

	//for(int i = 0; i < GetPtoCIndx(pIndx); ++i)
	//{
	//	cIndx = GetPtoC(pIndx, i, 0);
	//	oIndx = GetPtoC(pIndx, i, 1);

	//	if(cIndx == -1 || oIndx == -1){	continue;	}

	//	float defAmount = m_iceMove[cIndx]->GetDefAmount()/deformAmountSum;
	//	//defAmount *= defAmount;		//����

	//	m_iceMove[cIndx]->SetDefPriority(defAmount*10.0f/* - average*/);	//���ςƂǂ̂��炢�������邩��\��
	//}

	//TODO::�t�o�[�W�����@���܂������Ȃ��@������deformAmountSum��0������
	//float deformAmountSum = 0.0f;
	//for(int i = 0; i < GetPtoCIndx(pIndx); ++i)
	//{
	//	cIndx = GetPtoC(pIndx, i, 0);
	//	oIndx = GetPtoC(pIndx, i, 1);

	//	if(cIndx == -1 || oIndx == -1){	continue;	}

	//	deformAmountSum += m_iceMove[cIndx]->GetDefAmount();
	//}

	////�ʒu�C���x�x�N�g���Z�o
	//for(int i = 0; i < GetPtoCIndx(pIndx); ++i)
	//{
	//	cIndx = GetPtoC(pIndx, i, 0);
	//	oIndx = GetPtoC(pIndx, i, 1);

	//	if(cIndx == -1 || oIndx == -1){	continue;	}

	//	float defAmount = m_iceMove[cIndx]->GetDefAmount();
	//	pos += m_iceMove[cIndx]->GetVertexPos(oIndx) * (deformAmountSum - defAmount)/deformAmountSum;
	//	vel += m_iceMove[cIndx]->GetVertexVel(oIndx) * (deformAmountSum - defAmount)/deformAmountSum;
	//}
}

//SPH�@��SM�@�ŋ��߂����x�ƈʒu����`��ԁ@CPU
void IceObject::LinerInterPolationCPU(int pIndx, const Vec3& pos, const Vec3& vel)
{
	int sphIndx = pIndx*4;
	double intrps = 1.0-m_fInterPolationCoefficience[pIndx];	//��ԌW��

	s_sphPrtVel[sphIndx+0] = vel[0] * m_fInterPolationCoefficience[pIndx] + s_sphPrtVel[sphIndx+0] * intrps;
	s_sphPrtVel[sphIndx+1] = vel[1] * m_fInterPolationCoefficience[pIndx] + s_sphPrtVel[sphIndx+1] * intrps;
	s_sphPrtVel[sphIndx+2] = vel[2] * m_fInterPolationCoefficience[pIndx] + s_sphPrtVel[sphIndx+2] * intrps;

	s_sphPrtPos[sphIndx+0] = pos[0] * m_fInterPolationCoefficience[pIndx] + s_sphPrtPos[sphIndx+0] * intrps;
	s_sphPrtPos[sphIndx+1] = pos[1] * m_fInterPolationCoefficience[pIndx] + s_sphPrtPos[sphIndx+1] * intrps;
	s_sphPrtPos[sphIndx+2] = pos[2] * m_fInterPolationCoefficience[pIndx] + s_sphPrtPos[sphIndx+2] * intrps;
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

	SearchMeltParticle(prtList);		//�Z�𗱎q�T��

	m_iceStrct->StepObjMelt(prtList,clusterList, tetraList, cLayerList, tLayerList);	//�Z������
	
	ReConstructCluster(prtList, clusterList);	//TODO::�N���X�^�̍č\�z
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
		if( m_heatTransfer->getPhaseChange(i) != 1 )	continue;	//���]�ڂ̏����𖞂����Ă���
		if( m_heatTransfer->getPhase(i) != 2 )			continue;	//���ւƑ��]�ڂ��Ă���
		if( m_iceStrct->GetParticleNum() <= i)			continue;	//�Z���݂̂̎����̂Ƃ��ɕK�v�ɂȂ�D
		if( m_iceStrct->GetPtoCNum(i) == 0 )			continue;	//�N���X�^�Ɋ܂܂�Ă��Ȃ�
		if( m_iceStrct->GetPtoTNum(i) == 0 )			continue;	//�l�ʑ̂Ɋ܂܂�Ă��Ȃ�

		if(pList.size() > 50){	break;	}							//�Z�𗱎q���̐���

		m_fInterPolationCoefficience[i] = 0.0f;						//���`��Ԃ��Ȃ�
		m_heatTransfer->setPhaseChange(i, 0);						//���]�ڂ��I��������Ƃ�`����
		pList.push_back(i);											//�Z�𗱎q�̋L�^
	}

//�f�o�b�O
	//cout << __FUNCTION__ << endl;

	//for(unsigned i = 0; i < pList.size(); i++)
	//{
	//	cout << pList[i] << endl;
	//}

	//cout << "pList.size() = " << pList.size() << endl;
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

	if(pListSize == 0 || cListSize == 0){	return; }

	vector<int> checkTList;
	int j = 0, k = 0, l = 0;
	int icIndx = 0;
	int jtIndx = 0, joIndx = 0;
	int kpIndx = 0, ktIndx = 0, klIndx = 0;
	int lpIndx = 0;
	int pNum = 0;

	//�Ē�`�N���X�^�̏�����
	#pragma omp parallel
	{
	#pragma omp for
		for(int i = 0; i < cListSize; ++i)
		{
			if(std::find(pList.begin(), pList.end(), cList[i]) != pList.end()){	continue;	}
			GetMoveObj(cList[i])->Clear();
		}
	}//end #pragma omp parallel

	//�N���X�^�̍Ē�`
	#pragma omp parallel
	{
	#pragma omp for private(checkTList, j, k, l, icIndx, jtIndx, joIndx, kpIndx, ktIndx, klIndx, lpIndx, pNum)
		for(int i = 0; i < cListSize; ++i)
		{
			checkTList.clear();
			icIndx = cList[i];
			if(std::find(pList.begin(), pList.end(), icIndx) != pList.end()){	continue;	}

			//�N���X�^���Ē�`����ہC��{�ƂȂ闱�q��������l�ʑ̂��珉�����q�𓾂�D
			//�N���X�^�ԍ��������q�ԍ��Ȃ̂ɒ���
			//�ȉ����֐��ɂ���ƁC�G���[���o�Ă��܂������Ȃ�
			for(j = 0; j < GetPtoTIndx(icIndx); j++)
			{
				jtIndx = GetPtoT(icIndx, j, 0);
				joIndx = GetPtoT(icIndx, j, 1);

				if(jtIndx == -1 || joIndx == -1){ continue;	}
				if(std::find(checkTList.begin(), checkTList.end(), jtIndx) != checkTList.end())
				{	continue;	}
				else
				{	checkTList.push_back(jtIndx);	}

				//�l�ʑ̂̑S�Ă̗��q���N���X�^�ɓo�^
				for(k = 0; k < GetTtoPIndx(jtIndx); k++)
				{
					kpIndx = GetTtoP(jtIndx, k);

					if(kpIndx == -1){	continue;	}
					if(GetMoveObj(icIndx)->CheckIndx(kpIndx)){	continue;	}

					pNum = GetMoveObj(icIndx)->GetNumVertices();

					GetMoveObj(icIndx)->AddVertex(Vec3(s_sphPrtPos[kpIndx*4+0], s_sphPrtPos[kpIndx*4+1], s_sphPrtPos[kpIndx*4+2]), 1.0, kpIndx);
					//GetMoveObj(icIndx)->SetVelocity(pNum, Vec3(s_sphPrtVel[kpIndx*4+0], s_sphPrtVel[kpIndx*4+1], s_sphPrtVel[kpIndx*4+2]));
					GetMoveObj(icIndx)->SetAlphas(pNum, 1.0);
					GetMoveObj(icIndx)->SetBetas (pNum, 0.0);
					GetMoveObj(icIndx)->SetLayer (pNum, 0);
				}
			}

			//�ߖT�l�ʑ̂�layer���ǂ�C���q��ǉ����Ă���
			//TODO::�s����ɂȂ�Ȃ�Clayer�������ق�beta��������
			for(j = 0; j < GetPtoTIndx(icIndx); ++j)
			{
				jtIndx = GetPtoT(icIndx, j, 0);
				joIndx = GetPtoT(icIndx, j, 1);

				if(jtIndx == -1 || joIndx == -1){ continue;	}

				for(k = 0; k < GetNTNum(jtIndx); k++)
				{
					ktIndx = GetNeighborTetra(jtIndx, k, 0);
					klIndx = GetNeighborTetra(jtIndx, k, 1);

					if(ktIndx == -1 || klIndx == -1){	continue;	}
					if(std::find(checkTList.begin(), checkTList.end(), ktIndx) != checkTList.end())
					{	continue;	}
					else
					{	checkTList.push_back(ktIndx);	}

					//�l�ʑ̂̑S�Ă̗��q���N���X�^�ɓo�^
					for(l = 0; l < GetTtoPIndx(ktIndx); l++)
					{
						lpIndx = GetTtoP(ktIndx, l);

						if(lpIndx == -1){	continue;	}
						if(GetMoveObj(icIndx)->CheckIndx(lpIndx)){	/*cout << "contineue kpIndx = " << kpIndx << endl;*/ continue;	}

						pNum = GetMoveObj(icIndx)->GetNumVertices();

						GetMoveObj(icIndx)->AddVertex( Vec3(s_sphPrtPos[lpIndx*4+0], s_sphPrtPos[lpIndx*4+1], s_sphPrtPos[lpIndx*4+2] ), 1.0, lpIndx);
						//GetMoveObj(icIndx)->SetVelocity(pNum, Vec3(s_sphPrtVel[lpIndx*4+0], s_sphPrtVel[lpIndx*4+1], s_sphPrtVel[lpIndx*4+2]));
						GetMoveObj(icIndx)->SetAlphas(pNum, 1.0);
						GetMoveObj(icIndx)->SetBetas (pNum, 0.0);
						GetMoveObj(icIndx)->SetLayer (pNum, klIndx);
					}
				}
			}
		}
	}//end #pragma omp parallel
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
	int max = m_iceMove[0]->GetNumVertices();
	int min = m_iceMove[0]->GetNumVertices();

	for(int i = 0; i < sm_clusterNum; i++)
	{
		int num =  m_iceMove[i]->GetNumVertices();
		
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

		float defAmount = m_iceMove[i]->GetDefAmount();
		if(defAmount < 0.5f){	continue;	}

		//cout << "i = " << i << "::" <<  defAmount << endl;

#ifdef USE_SELECTED
		deformClusterNum += SELECTED;
#else
		deformClusterNum++;
#endif
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

		float defAmount = m_iceMove[i]->GetDefAmount();

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
		m_SurfSm.GetDevicePRTtoPTHPointer(),
		m_SurfSm.GetDevicePTHandPrfxSetPointer(),
		m_SurfSm.GetDecvicePrfxPos(),
		m_SurfSm.GetDecvicePrfxApq(),
		//IceStruct------------------------------------------
		m_iceStrct->GetDeviceCtoPointer(),
		m_iceStrct->GetDeviceCtoPNumPointer(),
		m_iceStrct->GetCtoPMax(),
		2,
		0.01		//TODO: 0.02�ł́H���낢��΂�΂���ۂ�
	);

//CPU��
	//prefixSum�̍X�V
	m_SurfSm.UpdatePrefixSum();

	//�N���X�^�̃p�����[�^�X�V
	#pragma omp parallel for
	for(int i = 0; i < sm_clusterNum; ++i)
	{	
		if(GetPtoCNum(i) == 0){	continue;	}
		m_iceMove[i]->SetNowCm(m_SurfSm.CalcCmSum(i));		//�d�S�̍X�V
		m_iceMove[i]->SetApq(m_SurfSm.CalcApqSum(i));		//�ό`�s��̍X�V
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

	//	Vec3 cm = m_iceMove[iCrstr]->GetCm();
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

	//	rxMatrix3 apq = m_iceMove[iCrstr]->GetApq();
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