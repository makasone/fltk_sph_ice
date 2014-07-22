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

IceObject::IceObject()
{
}

IceObject::~IceObject()
{
}

//���ꂼ��̃N���X�E�ϐ��̏�����
void IceObject::InitIceObj(int pMaxNum, int cMaxNum, int tMaxNum)
{	cout << __FUNCTION__ << endl;
	//���̂̍\���̏�����
	m_iceStrct = new IceStructure(pMaxNum, cMaxNum, tMaxNum);

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

//�N���X�^�̏�����	//�ꎞ�I�Ȏ���
void IceObject::InitCluster(Ice_SM* sm)
{
	m_iceMove.push_back(sm);	//�|�C���^���R�s�[���Ă��邾��
}

//�N���X�^�̏�����
void IceObject::InitCluster(Vec3 boundarySpaceHigh, Vec3 boundarySpaceLow, float timeStep)
{
	//�ȉ��͂�����ƈڍs�����ꍇ�̃R�[�h
	//�ϐ��̏�����
	for(vector<Ice_SM*>::iterator it = m_iceMove.begin(); it != m_iceMove.end(); ++it)
	{
		if(*it) delete *it;
	}
	
	sm_clusterNum = 0;	//�N���X�^���̏�����

	Ice_SM::SetParticlePosAndVel(s_sphPrtPos, s_sphPrtVel);	//TOD::�|�C���^���R�s�[����̂ł͂Ȃ��ォ��v�b�V�������ق�������

	////�K�w�\���ɂ�鍄���̊m��
	////�����P
	////MakeClusterHigh();

	//�l�ʑ̃��X�g�����ɁC���q���ɃN���X�^�쐬
	for(int i = 0; i < sm_particleNum; ++i)
	{
		//�N���X�^������
		m_iceMove.push_back(new Ice_SM(sm_clusterNum));
		m_iceMove[sm_clusterNum]->SetSimulationSpace(boundarySpaceLow, boundarySpaceHigh);
		m_iceMove[sm_clusterNum]->SetTimeStep(timeStep);
		m_iceMove[sm_clusterNum]->SetCollisionFunc(0);
		m_iceMove[sm_clusterNum]->SetStiffness(1.0, 1.0);

		//�l�ʑ̃��X�g�����ɁC���q���ɃN���X�^�쐬
		SetClusterMoveInfo(i);

		sm_clusterNum++;
	}

	////MakeClusterFromNeight();
	////MakeOneCluster();

	//�N���X�^�Ɋւ���GPU�̏�����
	//���ׂẴN���X�^�ɃA�N�Z�X�ł���悤�Ƀ|�C���^��n���Ă���
	Ice_SM::InitGPU(m_iceMove, sd_sphPrtPos, sd_sphPrtVel, sm_particleNum);
	Ice_SM::InitFinalParamPointer(sm_clusterNum);

	//�N���X�^�̃f�[�^��GPU�֏�����
	for(int i = 0; i < sm_particleNum; i++)
	{
		m_iceMove[i]->InitGPU_Instance();
	}

	//�\�ʗ��q�ɏd�݂�t����
	//((RXSPH*)m_pPS)->DetectSurfaceParticles();								//�\�ʗ��q���o
	//int *surfaceParticles = (int *)( ((RXSPH*)m_pPS)->GetArraySurf() );		//�\�ʗ��q
	//float weightedMass = 300.0f;

	////�\�ʗ��q�̎��ʂ���
	//for(int pIndx = 0; pIndx < ICENUM; pIndx++)
	//{
	//	if(surfaceParticles[pIndx] != 1){ continue;	}
	//	else{	cout << "surfaceParticles[" << pIndx << "] = " << surfaceParticles[pIndx] << endl;	}

	//	for(int pNum = 0; pNum < m_ice->GetPtoCIndx(pIndx); ++pNum)
	//	{
	//		int cIndx = m_ice->GetPtoT(pIndx, pNum, 0);
	//		int oIndx = m_ice->GetPtoT(pIndx, pNum, 1);

	//		if(cIndx == -1 || oIndx == -1){ continue;	}

	//		m_sm_cluster[cIndx]->SetMass(oIndx, weightedMass);
	//	}
	//}

	//�\�ʗ��q�̃N���X�^���ׂĂ̗��q
	//for(int cIndx = 0; cIndx < ICENUM; ++cIndx)
	//{
	//	if(surfaceParticles[cIndx] != 1){ continue;	}
	//	else{	cout << "surfaceParticles[" << cIndx << "] = " << surfaceParticles[cIndx] << endl;	}

	//	for(int pIndx = 0; pIndx < /*m_ice->GetCtoPIndx(cIndx)*/m_sm_cluster[cIndx]->GetNumVertices(); ++pIndx)
	//	{
	//		m_sm_cluster[cIndx]->SetMass(pIndx, weightedMass);
	//	}

	//	//	for(int pNum = 0; pNum < m_ice->GetPtoCIndx(pIndx); ++pNum)
	//	//	{
	//	//		int cIndx = m_ice->GetPtoT(pIndx, pNum, 0);
	//	//		int oIndx = m_ice->GetPtoT(pIndx, pNum, 1);

	//	//		if(cIndx == -1 || oIndx == -1){ continue;	}

	//	//		m_sm_cluster[cIndx]->SetMass(oIndx, weightedMass);
	//	//	}
	//}

	//�S�̂̔����̗��q
	//�Ȃ������܂������Ȃ��H
	//for(int cIndx = 0; cIndx < ICENUM; cIndx++)
	//{
	//	//for(int pNum = 0; pNum < m_sm_cluster[cIndx]->GetNumVertices(); pNum++)
	//	//{
	//	//	int pIndx = m_sm_cluster[cIndx]->GetParticleIndx(pNum);
	//	//	
	//	//	if(surfaceParticles[pIndx] != 1){ continue;	}
	//	//	else{	cout << "surfaceParticles[" << pIndx << "] = " << surfaceParticles[pIndx] << endl;	}
	//	//
	//	//	m_sm_cluster[cIndx]->SetMass(pNum, weightedMass);
	//	//}
	//	
	//	for(int pNum = 0; pNum < m_sm_cluster[cIndx]->GetNumVertices(); pNum++)
	//	{
	//		int pIndx = m_sm_cluster[cIndx]->GetParticleIndx(pNum);
	//		
	//		//if(surfaceParticles[pIndx] != 1){ continue;	}
	//		//else{	cout << "surfaceParticles[" << pIndx << "] = " << surfaceParticles[pIndx] << endl;	}
	//	
	//		if(pIndx > ICENUM/2 || pIndx < 0){ continue;	}
	//		else{	/*cout << "pInx = " << pIndx << endl;*/	}

	//		m_sm_cluster[cIndx]->SetMass(pNum, weightedMass);
	//	}
	//}
	

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

	DebugClusterInfo();
}

//�I�u�W�F�N�g�̍\���̏�����
void IceObject::InitStrct()
{
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
	cout << __FUNCTION__ << " check0" << endl;

	//�������m��
	InitClusterInfo();
	cout << __FUNCTION__ << " check1" << endl;

	//���q���������Ă���N���X�^���̔z����R�s�[
	int *PtoCNum = new int[sm_particleNum];
	cout << __FUNCTION__ << " check2" << endl;
	for(int i = 0; i < sm_particleNum; i++)
	{
		PtoCNum[i] = GetPtoCNum(i);
	}
	cout << __FUNCTION__ << " check3" << endl;

	//�N���X�^�Ɨ��q�̊֘A���̓o�^
	for(int i = 0; i < sm_particleNum; i++)
	{
		SetClusterStrctInfo(i, PtoCNum);	//�J�E���g��O�ōs���Ă��邽�߁C��������g��
	}
	cout << __FUNCTION__ << " check4" << endl;

	delete[] PtoCNum;

	//TODO::�N���X�^�Ǝl�ʑ̂̊֘A���̓o�^
	cout << __FUNCTION__ << " check5" << endl;
	//m_ice->InitPath(p, v, m_sm_cluster, ICENUM);			//�������̂��߂̃p�X�쐬


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


//�N���X�^�̉^���v�Z���̓o�^
void IceObject::SetClusterMoveInfo(int pIndx)
{
	//((RXSPH*)m_pPS)->DetectSurfaceParticles();								//�\�ʗ��q���o
	//int *surfaceParticles = (int *)( ((RXSPH*)m_pPS)->GetArraySurf() );		//�\�ʗ��q
	//float weightedMass = 300.0f;
	//int layer = 1;

	//���闱�q���܂܂��l�ʑ̂��ߖT�N���X�^�Ƃ��C�e�l�ʑ̂Ɋ܂܂�闱�q�ŃN���X�^���쐬
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

			//if(jpIndx < ICENUM/2)					//���q�̔����ɏd�ݕt��
			//if(surfaceParticles[jpIndx] == 1)
			//if( (jpIndx < 13*13*layer) || (ICENUM-13*13*layer < jpIndx) )
			//if( ()									//x�w
			//||										//y�w
			//)
			//{
			//	//cout << "jpIndx == " << jpIndx << endl;
			//	mass = weightedMass;
			//}

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

				//if(kpIndx < ICENUM/2)					//���q�̔����ɏd�ݕt��
				//if(surfaceParticles[kpIndx] == 1)
				//if( (kpIndx < 13*13*layer) || (ICENUM-13*13*layer < kpIndx) ) 
				//{
				//	//cout << "kpIndx == " << kpIndx << endl;
				//	mass = weightedMass;
				//}

				m_iceMove[pIndx]->AddVertex( Vec3(s_sphPrtPos[kpIndx*4+0], s_sphPrtPos[kpIndx*4+1], s_sphPrtPos[kpIndx*4+2] ), mass, kpIndx);
				m_iceMove[pIndx]->SetAlphas(pNum, 1.0);
				m_iceMove[pIndx]->SetBetas (pNum, 0.0);
				m_iceMove[pIndx]->SetLayer (pNum, jlIndx);
				//cout << "pIndx = " << pIndx << " vNum = " << m_sm_cluster[pIndx]->GetNumVertices() << endl;
			}
		}
	}
}

//
void IceObject::SetClusterStrctInfo(int cIndx, int* PtoCNum)
{
	//���q�������Ă���l�ʑ̂̔ԍ���o�^���邽�߂̏���
	//pCountList�ɂ́CcIndx�Ԗڂ̃N���X�^�Ɋ܂܂��e���q���C���ꂼ�ꂢ���̃N���X�^�ɑ����邩�����߂ĕۑ�����
	int* pCountList = new int[m_iceMove[cIndx]->GetNumVertices()];
	int* pLayerList = new int[m_iceMove[cIndx]->GetNumVertices()];			//���q�̏������C���[
	
	for(int i = 0; i < m_iceMove[cIndx]->GetNumVertices(); i++)
	{
		int pIndx = m_iceMove[cIndx]->GetParticleIndx(i);

		//��������z�肵�Ă��Ȃ����߁C����ď㏑�����Ă��܂��Ă���
		//-1��T�����ď㏑������̂ɐ؂�ւ���D
		for(int j = 0; j < GetPtoCMax(); j++)
		{
			if(GetPtoC(pIndx, j, 0) != -1 || GetPtoC(pIndx, j, 1) != -1){	continue;	}
			
			if(j >= GetPtoCIndx(pIndx))
			{
				SetPtoCIndx(pIndx, j+1);		//���݂�Indx���傫���Ȃ�X�V
			}

			pCountList[i] = j;
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
		pIndxList.push_back(pIndx);
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

	//TODO::�N���X�^��GPU�������������ɒu��

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
void IceObject::StepObjMove()
{

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
	for(int i = 0; i < sm_particleNum; ++i)
	{
		if(GetParticleNum() <= i){	continue;	}	//�Z���݂̂̎����̂Ƃ��ɕK�v�ɂȂ�D
		if(GetPtoCNum(i) <= 0){		continue;	}
	
		Vec3 pos,vel;
	
		//�ő̉^���̍ŏI�ʒu�v�Z
		CalcAverageCPU(i, pos, vel);

		//�t�̂ƌő̂̕��
		LinerInterPolationCPU(i, pos, vel);	//���`���
	}
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
#else
	for(int i = 0; i < sm_particleNum; ++i)
	{
		if(GetParticleNum() <= i){	continue;	}	//�Z���݂̂̎����̂Ƃ��ɕK�v�ɂȂ�D
		if(GetPtoCNum(i) <= 0){		continue;	}
	
		Vec3 pos,vel;
	
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

		for(int j = 0; j < GetPtoCIndx(pIndx); ++j)
		{
			int jcIndx = GetPtoC(pIndx, j, 0);
			int joIndx = GetPtoC(pIndx, j, 1);

			if(jcIndx == -1 || joIndx == -1){	continue;	}

			pos += m_iceMove[jcIndx]->GetVertexPos(joIndx);
			vel += m_iceMove[jcIndx]->GetVertexVel(joIndx);

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

//�f�o�b�O

//�l�ʑ�
void IceObject::DebugTetraInfo()
{
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
{
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