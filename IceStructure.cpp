//-----------------------------------------------------------------------------
// �C���N���[�h�t�@�C��
//-----------------------------------------------------------------------------
#include "IceStructure.h"

//���q��
int* IceStructure::sd_piPtoC;
int* IceStructure::sd_piPtoT;

int* IceStructure::sd_piPtoCNum;
int* IceStructure::sd_piPtoTNum;

int* IceStructure::sd_piPtoCIndx;
int* IceStructure::sd_piPtoTIndx;

//�N���X�^��
int* IceStructure::sd_piCtoP;

int* IceStructure::sd_piCtoPNum;

int* IceStructure::sd_piCtoPIndx;

//�l�ʑ́�
int* IceStructure::sd_piTtoP;

int* IceStructure::sd_piTtoPNum;
int* IceStructure::sd_piTtoCNum;

int* IceStructure::sd_piTtoPIndx;
int* IceStructure::sd_piTtoCIndx;

int* IceStructure::sd_piNeighborTetra;

int* IceStructure::sd_piNeighborTetraTNum;

/*!
 * @param[in] pNumMax�@�ő嗱�q��
 * @param[in] cNumMax�@�ő�N���X�^��
 * @param[in] tNumMax�@�ő�l�ʑ̐�
 */
IceStructure::IceStructure()
{
}

IceStructure::IceStructure(int pNumMax, int cNumMax, int tNumMax, int layer)
{	cout << __FUNCTION__ << endl;

	//�ő吔�̓o�^
	m_iPNumMax = pNumMax;
	m_iCNumMax = cNumMax;
	m_iTNumMax = tNumMax;

	//���q���̏�����
	m_iPtoCMax = m_iCNumMax*0.1f;	//1331 layer2 0.4 layer3 0.75
									//2197 layer2 0.4 layer3 0.4 layer4 0.5
									//4913 layer2 0.3
									//4913 layer1 0.1
									//12167 layer1 0.01
									//19683 layer1 0.01

	m_iPtoTMax = m_iTNumMax*0.1f;	//1331 layer2 0.3 layer3 0.5
									//2197 layer2 0.3 layer3 0.3
									//4913 layer2 0.3
									//12167 layer1 0.01
									//19683 layer1 0.01

cout << __FUNCTION__ << " check1" << endl;

	m_piPtoCNum = new int[m_iPNumMax];
	m_piPtoTNum = new int[m_iPNumMax];

	m_piPtoCIndx = new int[m_iPNumMax];
	m_piPtoTIndx = new int[m_iPNumMax];

	for(int i = 0; i < m_iPNumMax; i++)
	{
		m_piPtoCNum[i] = 0;
		m_piPtoTNum[i] = 0;
		
		m_piPtoCIndx[i] = 0;
		m_piPtoTIndx[i] = 0;
	}

	//�N���X�^���̏�����
	//CtoTMax�͗��q���Ɠ������̂Œ�`���Ȃ�
	//���������e�X�g����Ƃ���CtoT�̓R�����g��
	m_iCtoPMax = m_iPNumMax*0.1f;	//1331 layer2 0.5 layer3 0.75
									//2197 layer2 0.5 layre3 0.5
									//4913 layer2
									//4913 layer1 0.1
									//12167 layer1 0.01
									//19683 layer1 0.01
	//m_iCtoPMax = m_iPNumMax;		//�P��N���X�^

	m_piCtoPNum = new int[m_iCNumMax];
	m_piCtoPIndx = new int[m_iCNumMax];

cout << __FUNCTION__ << " check2" << endl;

	for(int i = 0; i < m_iCNumMax; i++)
	{
		m_piCtoPNum[i] = 0;
		//m_piCtoTNum[i] = 0;
		
		m_piCtoPIndx[i] = 0;
		//m_piCtoTIndx[i] = 0;
	}

	//�l�ʑ̏��̏�����
	//TtoPMax�͍ő�S�ŌŒ�
	//TtoCMax�͕K�v�Ȃ�
	m_piTtoPNum = new int[m_iTNumMax];
	m_piTtoCNum = new int[m_iTNumMax];

	m_piTtoPIndx = new int[m_iTNumMax];
	m_piTtoCIndx = new int[m_iTNumMax];

cout << __FUNCTION__ << " check3" << endl;

	for(int i = 0; i < m_iTNumMax; i++)
	{
		m_piTtoPNum[i] = 0;
		m_piTtoCNum[i] = 0;
		
		m_piTtoPIndx[i] = 0;
		m_piTtoCIndx[i] = 0;
	}

	//�ߖT�l�ʑ�
	int ntnSize = m_iTNumMax * 1.0f;
	m_piNTNum = new int[ntnSize];

cout << __FUNCTION__ << " check3.5" << endl;

	//m_iNeighborMax = m_iTNumMax*0.1f;		//1331  layer2 0.3 layer3 0.75
											//2197  layer2 0.3 layre3 0.3 layer4 0.4
											//4913  layer2 0.1
											//6859  layer1 0.01
											//9261  layer1 0.005
											//12167 layer1

	m_iNeighborMax = 150;					//layer1�Ȃ���v

	cout << "m_iTNumMax = " << m_iTNumMax << ", m_iNeighborMax = " << m_iNeighborMax << endl;

	m_mk3DiNeighborTetra.SetSize(m_iTNumMax, m_iNeighborMax, 2);

cout << __FUNCTION__ << " check4" << endl;

	for(int i = 0; i < m_iTNumMax; i++)
	{
		m_piNTNum[i] = 0;

		for(int j = 0; j < m_iNeighborMax; j++)
		{
			for(int k = 0; k < 2; k++)
			{
				m_mk3DiNeighborTetra(i, j, k) = -1;
			}
		}
	}

	m_iLayer = layer;

	//�t���O
	//m_pbPFlag = new bool[m_iPNumMax];
	//m_pbCFlag = new bool[m_iCNumMax];
	m_pbTFlag = new bool[m_iTNumMax];

	//ResetPFlag(m_iPNumMax);
	//ResetCFlag(m_iCNumMax);
	ResetTFlag(m_iTNumMax);

	//�I��I�^���v�Z
	m_psuSelectClusterIndx = new short unsigned[m_iPNumMax];
	
	for(int i = 0; i < m_iPNumMax; i++)
	{
		UpdateMotionCalcCluster(i, 1);
	}

cout << __FUNCTION__ << " check5" << endl;
}

IceStructure::~IceStructure(void)
{
}

//------------------------------------------�t���O�Ǘ�-------------------------------------------------
/*!
 * �T���p���q�t���O�̏�����
 */
void IceStructure::ResetPFlag(int endNum)
{
	for(int i = 0; i < endNum; i++)
	{
		m_pbPFlag[i] = false;
	}
}

/*!
 * �T���p�N���X�^�t���O�̏�����
 */
void IceStructure::ResetCFlag(int endNum)
{
	for(int i = 0; i < endNum; i++)
	{
		m_pbCFlag[i] = false;
	}
}

/*!
 * �T���p�l�ʑ̃t���O�̏�����
 */
void IceStructure::ResetTFlag(int endNum)
{
	for(int i = 0; i < endNum; i++)
	{
		m_pbTFlag[i] = false;
	}
}
//------------------------------------------�t���O�Ǘ�-------------------------------------------------

//------------------------------------------������-------------------------------------------------
/*!
 * �v�Z���@�l�ʑ̏��̗̈�m�ہ@���q�x�[�X����
 */
void IceStructure::InitTetraInfo()
{
	//���q���l�ʑ�
	m_mk3DiPtoT.SetSize(m_iPNumMax, m_iPtoTMax, 2);

	for(int i = 0; i < m_iPNumMax; i++)
	{		
		for(int j = 0; j < m_iPtoTMax; j++)
		{
			for(int k = 0; k < 2; k++)
			{
				m_mk3DiPtoT(i, j, k) = -1;
			}
		}
	}

	//�l�ʑ́����q
	m_mk2DiTtoP.SetSize(m_iTNumMax, 4);
	
	for(int i = 0; i < m_iTNumMax; i++)
	{
		for(int j = 0; j < 4; j++)
		{
			m_mk2DiTtoP(i, j) = -1;
		}
	}

	//�z��̓Y�����̏�����
	for(int i = 0; i < m_iPNumMax; i++)
	{
		m_piPtoTIndx[i] = m_piPtoTNum[i];
	}

	for(int i = 0; i < m_iTNumMax; i++)
	{
		m_piTtoPIndx[i] = m_piTtoPNum[i];
	}
}

/*!
 * �v�Z���@�N���X�^���̗̈�m�ہ@���q�x�[�X����
 */
void IceStructure::InitClusterInfo()
{	cout << __FUNCTION__ << endl;

	//���q���N���X�^
	m_mk3DiPtoC.SetSize(m_iPNumMax, m_iPtoCMax, 3);
	cout << __FUNCTION__ << ", check1" << endl;

	for(int i = 0; i < m_iPNumMax; i++)
	{
		for(int j = 0; j < m_iPtoCMax; j++)
		{
			//[0] = �N���X�^�ԍ��C[1] = �N���X�^���ł̔ԍ��C[2] = layer�ԍ�
			for(int k = 0; k < 3; k++)
			{
				m_mk3DiPtoC(i, j, k) = -1;
			}
		}
	}

	//�N���X�^�����q
	m_mk3DiCtoP.SetSize(m_iCNumMax, m_iCtoPMax, 2);
	cout << __FUNCTION__ << ", check2" << endl;

	for(int i = 0; i < m_iCNumMax; i++)
	{
		for(int j = 0; j < m_iCtoPMax; j++)
		{
			//[0] = ���q�ԍ��C[1] = layer�ԍ�
			for(int k = 0; k < 2; k++)
			{
				m_mk3DiCtoP(i, j, k) = -1;
			}
		}
	}

	//�z��̓Y����Indx�̏�����
	for(int i = 0; i < m_iPNumMax; i++)
	{
		m_piPtoCIndx[i] = m_piPtoCNum[i];
	}

	for(int i = 0; i < m_iCNumMax; i++)
	{
		m_piCtoPIndx[i] = m_piCtoPNum[i];
	}

	cout << __FUNCTION__ << ", check3" << endl;
}

//GPU�����ŗp����f�[�^�̏�����
//���̂Ƃ���̓R�s�[�ł܂��Ȃ�
void IceStructure::InitGPU()
{	cout << __FUNCTION__ << endl;

	//�f�o�C�X���̃��������m��
		//���q��
	cudaMalloc((void**)&sd_piPtoT,		sizeof(int) * m_iPNumMax * m_iPtoTMax * 2);
	cudaMalloc((void**)&sd_piPtoC,		sizeof(int) * m_iPNumMax * m_iPtoCMax * 3);

	cudaMalloc((void**)&sd_piPtoCNum,	sizeof(int) * m_iPNumMax);
	cudaMalloc((void**)&sd_piPtoTNum,	sizeof(int) * m_iPNumMax);

	cudaMalloc((void**)&sd_piPtoCIndx,	sizeof(int) * m_iPNumMax);
	cudaMalloc((void**)&sd_piPtoTIndx,	sizeof(int) * m_iPNumMax);

		//�N���X�^��
	cudaMalloc((void**)&sd_piCtoP,		sizeof(int) * m_iCNumMax * m_iCtoPMax * 2);

	cudaMalloc((void**)&sd_piCtoPNum,	sizeof(int) * m_iCNumMax);
	cudaMalloc((void**)&sd_piCtoPIndx,	sizeof(int) * m_iCNumMax);

		//�l�ʑ́�
	cudaMalloc((void**)&sd_piTtoP,	sizeof(int) * m_iTNumMax * 4);

	cudaMalloc((void**)&sd_piTtoPNum,	sizeof(int) * m_iTNumMax);
	cudaMalloc((void**)&sd_piTtoCNum,	sizeof(int) * m_iTNumMax);

	cudaMalloc((void**)&sd_piTtoPIndx,	sizeof(int) * m_iTNumMax);
	cudaMalloc((void**)&sd_piTtoCIndx,	sizeof(int) * m_iTNumMax);

		//�ߖT�l�ʑ�
	cudaMalloc((void**)&sd_piNeighborTetra,		sizeof(int) * m_iTNumMax * m_iNeighborMax * 2);
	cudaMalloc((void**)&sd_piNeighborTetraTNum,	sizeof(int) * m_iTNumMax);

	//������
	cudaMemcpy(sd_piPtoCNum, m_piPtoCNum, sizeof(int) * m_iPNumMax, cudaMemcpyHostToDevice);
	cudaMemcpy(sd_piPtoTNum, m_piPtoTNum, sizeof(int) * m_iPNumMax, cudaMemcpyHostToDevice);

	cudaMemcpy(sd_piPtoCIndx, m_piPtoCIndx, sizeof(int) * m_iPNumMax, cudaMemcpyHostToDevice);
	cudaMemcpy(sd_piPtoTIndx, m_piPtoTIndx, sizeof(int) * m_iPNumMax, cudaMemcpyHostToDevice);

	cudaMemcpy(sd_piCtoPNum,	m_piCtoPNum,	sizeof(int) * m_iCNumMax, cudaMemcpyHostToDevice);
	cudaMemcpy(sd_piCtoPIndx,	m_piCtoPIndx,	sizeof(int) * m_iCNumMax, cudaMemcpyHostToDevice);

	cudaMemcpy(sd_piTtoPNum,	m_piTtoPNum,	sizeof(int) * m_iTNumMax, cudaMemcpyHostToDevice);
	cudaMemcpy(sd_piTtoCNum,	m_piTtoCNum,	sizeof(int) * m_iTNumMax, cudaMemcpyHostToDevice);

	cudaMemcpy(sd_piTtoPIndx,	m_piTtoPIndx,	sizeof(int) * m_iTNumMax, cudaMemcpyHostToDevice);
	cudaMemcpy(sd_piTtoCIndx,	m_piTtoCIndx,	sizeof(int) * m_iTNumMax, cudaMemcpyHostToDevice);

		//Vector�ŊǗ����Ă���f�[�^��z��ɂ��邽�߂�data()���g���Ă���
	cudaMemcpy(sd_piNeighborTetra,		m_mk3DiNeighborTetra.Get().data(),	sizeof(int) * m_iTNumMax * m_iNeighborMax * 2,	cudaMemcpyHostToDevice);
	cudaMemcpy(sd_piNeighborTetraTNum,	m_piNTNum,							sizeof(int) * m_iTNumMax,						cudaMemcpyHostToDevice);

	cudaMemcpy(sd_piPtoT,	m_mk3DiPtoT.Get().data(),	sizeof(int) * m_iPNumMax * m_iPtoTMax * 2,	cudaMemcpyHostToDevice);
	cudaMemcpy(sd_piPtoC,	m_mk3DiPtoC.Get().data(),	sizeof(int) * m_iPNumMax * m_iPtoCMax * 3,	cudaMemcpyHostToDevice);

	cudaMemcpy(sd_piTtoP,	m_mk2DiTtoP.Get().data(),	sizeof(int) * m_iTNumMax * 4,				cudaMemcpyHostToDevice);

	cudaMemcpy(sd_piCtoP,	m_mk3DiCtoP.Get().data(),	sizeof(int) * m_iCNumMax * m_iCtoPMax * 2,	cudaMemcpyHostToDevice);

////�f�o�b�O
//	int* testA = new int[m_iPNumMax];
//	int* testB = new int[m_iPNumMax];
//	int* testC = new int[m_iCNumMax];
//	int* testD = new int[m_iCNumMax];
//
//	//�f�o�C�X����z�X�g�ւ̃R�s�[
//	cudaMemcpy(testA, sd_piPtoCIndx, sizeof(int) * m_iPNumMax, cudaMemcpyDeviceToHost);
//	cudaMemcpy(testB, sd_piPtoTIndx, sizeof(int) * m_iPNumMax, cudaMemcpyDeviceToHost);
//	cudaMemcpy(testC, sd_piCtoPNum,	 sizeof(int) * m_iCNumMax, cudaMemcpyDeviceToHost);
//	cudaMemcpy(testD, sd_piCtoPIndx, sizeof(int) * m_iCNumMax, cudaMemcpyDeviceToHost);
//
	////�z�X�g���̃f�[�^��]���������ʂ��_���v
	//ofstream ofs( "DtoH_Test.txt" );
	//ofs << "DtoH_Test" << endl;

	////for(int i = 0; i < m_iPNumMax; i++)
	//for(int i = 0; i < m_iTNumMax * m_iNeighborMax * 2; i++)
	//{
	//	////�������Ȃ�0�ɂȂ�͂�
	//	//ofs << i << "      host-device:: "
	//	//	<< abs(m_piPtoCIndx[i]-testA[i]) << ", "
	//	//	<< abs(m_piPtoTIndx[i]-testB[i]) << ", "
	//	//	<< abs(m_piCtoPNum[i] -testC[i]) << ", "
	//	//	<< abs(m_piCtoPIndx[i]-testD[i]) << ", "
	//	//	<< endl;

	//	//if(m_mk3DiNeighborTetra.Get()[i] == -1 && m_mk3DiNeighborTetra.Get().data()[i]){	continue;	}

	//	//ofs << i << "        data()-Get() = " << m_mk3DiNeighborTetra.Get().data()[i] - m_mk3DiNeighborTetra.Get()[i] << endl;
	//	//ofs << i << " Get() = " << m_mk3DiNeighborTetra.Get()[i] << endl;
	//	//ofs << i << " Get() = " << m_mk3DiNeighborTetra.Get().data()[i] << endl;
	//}
//
//	delete[] testA;
//	delete[] testB;
//	delete[] testC;
//	delete[] testD;
}

//------------------------------------------__������-------------------------------------------------

//------------------------------------------�Q���ω�-------------------------------------------------
void IceStructure::StepObjMelt(
	vector<unsigned>& pList,
	vector<unsigned>& cList,
	vector<unsigned>& tList,
	vector<unsigned>& cLayerList,
	vector<unsigned>& tLayerList)
{
	QueryCounter counter;
	counter.Start();
		double end = counter.End();
	SearchReconstruct_Tetra_Melt(pList, tList, tLayerList);			//�Ē�`�l�ʑ̂̒T��
	SearchReconstruct_Cluster_Melt(pList, cList, cLayerList);		//�Ē�`�N���X�^�̒T��
	UpdateInfo_Melt_PandT(pList);									//���q�E�l�ʑ̏��̍X�V
	UpdateInfo_Melt_PandC(pList, cList);							//���q�E�N���X�^���̍X�V

	////CheckDeleteCluster();													//����C��܊֌W�ɂ���N���X�^���폜
	////CheckDeleteTetra(viTetraList, viTLayerList);							//����C��܊֌W�ɂ���l�ʑ̂��폜

	SetInfo_Tetra(pList, tList, tLayerList);						//���q�E�ߖT�l�ʑ̏��̍Ē�`
	//RXTIMER("SetInfo_Tetra");

//�f�o�b�O
	//DebugStepObjMelt(pList, cList);
}

void IceStructure::StepObjFreeze()
{
	vector<int> viParticleList;														//�Ìł������q�W��
	vector<int> viClusterList;														//�Ē�`����N���X�^�̏W��
	vector<int> viCLayerList;														//�Ē�`����N���X�^�̃��C���[
	vector<int> viTetraList;														//�Ē�`����l�ʑ̂̏W��
	vector<int> viTLayerList;														//�Ē�`����l�ʑ̂̃��C���[

	//SearchFreezeParticle(viParticleList);											//�Ìŗ��q�̒T��
	//SetFreezeTetraInfo(viParticleList);												//�Ìŗ��q�Ɋւ���l�ʑ̂̍쐬
	//SetFreezeClusterInfo(viParticleList);											//�Ìŗ��q�Ɋւ���N���X�^�̍쐬
	//SearchReconstructTetra_Freeze(viParticleList, viTetraList, viTLayerList);		//�Ē�`�l�ʑ̂̒T��
	//SearchReconstructCluster_Freeze(viParticleList, viClusterList, viCLayerList);	//�Ē�`�N���X�^�̒T��

	////CheckDeleteCluster();															//����C��܊֌W�ɂ���N���X�^���폜
	////CheckDeleteTetra(viTetraList, viTLayerList);									//����C��܊֌W�ɂ���l�ʑ̂��폜

	//SetTetraInfo(viParticleList, viTetraList, viTLayerList);						//���q�E�ߖT�l�ʑ̏��̍Ē�`
	//SetClusterInfo(viParticleList, viClusterList, viCLayerList);					//���q�E�N���X�^���̍Ē�`

	////�f�o�b�O
	//if(viParticleList.size() == 0 || viClusterList.size() == 0){	return;	}
	//cout << "Debug " << __FUNCTION__ << "viParticleList.size = " << viParticleList.size() << " " << endl;

	//for(unsigned i = 0; i < viParticleList.size(); i++)
	//{
	//	cout << " " << viParticleList[i];
	//}
	//cout << endl;

	//cout << "viClusterList.size =  " << viClusterList.size() << " ";
	//for(unsigned i = 0; i < viClusterList.size(); i++)
	//{
	//	cout << " " << viClusterList[i];
	//}
	//cout << endl;

	//cout << "viCLayerList:: ";
	//for(unsigned i = 0; i < viCLayerList.size(); i++)
	//{
	//	cout << " " << viCLayerList[i];
	//}
	//cout << endl;

	//cout << "viTetraList.size = " << viTetraList.size() << " ";
	//for(unsigned i = 0; i < viTetraList.size(); i++)
	//{
	//	cout << " " << viTetraList[i];
	//}
	//cout << endl;

	//cout << "viTLayerList:: "3
	//for(unsigned i = 0; i < viTLayerList.size(); i++)
	//{
	//	cout << " " << viTLayerList[i];
	//}
	//cout << endl;

	//�N���X�^�����q
	//for(int i = 0; i < m_iClusteresNum; i++){	m_ice->DebugCtoP(i);	}

	//���q���N���X�^
	//for(int i = 0; i < m_pPS->GetNumParticles(); i++){	m_ice->DebugPtoC(i);	}

	//SM�N���X�^�Ɋ܂܂�闱�q�͋@�\�Ŋm�F�ł���
	//for(int i = 0; i < ICENUM; i++)
	//{
	//	cout << "pIndx = " << endl;

	//	for(int j = 0; j < m_sm_cluster[i]->GetNumVertices(); j++)
	//	{
	//		cout << " j = " << m_sm_cluster[i]->GetParticleIndx(j);
	//	}
	//	cout << endl;
	//}

	//�l�ʑ́����q�͋@�\�Ŋm�F�ł���

	//���q���l�ʑ�
	//for(int i = 0; i < m_pPS->GetNumParticles(); i++){	m_ice->DebugPtoT(i);	}

	//�ߖT�l�ʑ�
	//for(unsigned i = 0; i < m_vviTetraList.size(); i++ ){	m_ice->DebugNeighborTetra(i);	}
}

void IceStructure::SearchReconstruct_Tetra_Melt(const vector<unsigned>& pList, vector<unsigned>& tList, vector<unsigned>& lList)
{
	unsigned pListSize = pList.size();
	if(pListSize == 0){	return;}

	ResetTFlag(m_iTNum);	//�l�ʑ̒T���t���O�̏�����

	//�P�@���q���܂܂�Ă����C����l�ʑ�
	//�Q�@�P�̎l�ʑ̂̋ߖT�l�ʑ́@�߂��Ⴍ����d��
	//�܂�̓N���X�^���\�������l�ʑ́@TODO::�o��������Ȃ̂ŁC�����̌v�Z�R�X�g��������ΏC���\
	
	//�P �Z�𗱎q���܂܂�Ă����l�ʑ�
	for(unsigned i = 0; i < pListSize; i++)
	{
		int ipIndx = pList[i];

		for(int j = 0; j < GetPtoTIndx(ipIndx); j++)
		{
			if(GetPtoT(ipIndx, j, 0) == -1 || GetPtoT(ipIndx, j, 1) == -1){	continue;	}

			//�l�ʑ̒T���t���O������ɒT���������ǂ����𔻒�
			if( GetTFlag(GetPtoT(ipIndx, j, 0)) )	{	continue;								}
			else									{	SetTFlag(GetPtoT(ipIndx, j, 0), true);	}

			tList.push_back(GetPtoT(ipIndx, j, 0));
			lList.push_back(1);								//0��1���̔��f�͂ł��Ȃ��̂�1�ɍ��킹��D
		}
	}

	//return;	//�Q���������Ȃ��Ȃ�߂��Ⴍ���ᑁ���Ȃ�

	//�Q�@�P�̎l�ʑ̂̋ߖT�l�ʑ�
	int tetraNum = tList.size();
	for(int i = 0; i < tetraNum; i++)
	{
		int itIndx = tList[i];

		for(int j = 0; j < GetNTNum(itIndx); j++)
		{
			int jtIndx = GetNeighborTetra(itIndx, j, 0);
			int jlIndx = GetNeighborTetra(itIndx, j, 1);

			if(GetTFlag(jtIndx))
			{
				vector<unsigned>::iterator check = std::find(tList.begin(), tList.end(), jtIndx);
				
				//TODO::���Ɋ܂܂�Ă���̂Ȃ�Clayer���ׂď������ق���D�悷��
				int layerIndx = check - tList.begin();
				if(lList[layerIndx] > jlIndx)
				{
					lList[layerIndx] = jlIndx;
				}
				continue;
			}
			else
			{
				SetTFlag(jtIndx, true);	
			}

			tList.push_back(jtIndx);
			lList.push_back(jlIndx);
		}
	}
}

void IceStructure::SearchReconstruct_Cluster_Melt(const vector<unsigned>& pList, vector<unsigned>& cList, vector<unsigned>& lList)
{
	unsigned pListSize = pList.size();
	if(pListSize == 0){	return;}	

	//�Ē�`�N���X�^�́C�Z�𗱎q���������Ă����N���X�^
	for(unsigned i = 0;  i < pListSize; i++)
	{
		unsigned ipIndx = pList[i];

		for(unsigned j = 0, ctopIndx = GetCtoPIndx(ipIndx); j < ctopIndx; j++)
		{
			unsigned jcIndx = GetPtoC(ipIndx, j, 0);
			unsigned joIndx = GetPtoC(ipIndx, j, 1);
			unsigned jlIndx = GetPtoC(ipIndx, j, 2);

			if(jcIndx == -1 || joIndx == -1){	continue;	}

			vector<unsigned>::iterator check = std::find(cList.begin(), cList.end(), jcIndx);

			//TODO::���ɍĒ�`�N���X�^�Ƃ��Ď擾����Ă���Ȃ�Clayer���ׂď������ق���D�悷��
			if(check != cList.end())
			{
				int layerIndx = check - cList.begin();
				if(lList[layerIndx] > jlIndx)
				{
					lList[layerIndx] = jlIndx;
				}
				continue;
			}

			cList.push_back(jcIndx);
			lList.push_back(jlIndx);
		}
	}
}

void IceStructure::UpdateInfo_Melt_PandT(const vector<unsigned>& pList)
{
	unsigned pListSize = pList.size();
	if(pListSize == 0){	return;}	


	for(unsigned i = 0; i < pListSize; i++)
	{
		int ipIndx = pList[i];

		for(int j = 0; j < GetPtoTIndx(ipIndx); j++)
		{
			int tIndx = GetPtoT(ipIndx, j, 0);
			int oIndx = GetPtoT(ipIndx, j, 1);	//���̏ꍇ�͂��̂܂ܓY����
			
			if(tIndx == -1 || oIndx == -1){ continue;	}

			DeleteTtoP(tIndx, oIndx);
		}

		ClearPtoT(ipIndx);
	}
}

void IceStructure::UpdateInfo_Melt_PandC(const vector<unsigned>& pList, const vector<unsigned>& cList)
{
	int pListSize = pList.size();
	int cListSize = cList.size();

	if(pListSize == 0 || cListSize == 0){	return; }

	//���񏈗��ŗp����ϐ����܂Ƃ߂Ē�`
	int j= 0, k = 0;
	int icIndx = 0;
	int jpIndx = 0;

	//�Z�𗱎q���N���X�^�����C�Z���N���X�^�Ɋ܂܂�Ă������q�����菜��
	#pragma omp parallel
	{
	#pragma omp for private(j,k,icIndx,jpIndx)
		for(int i = 0; i < pListSize; i++)
		{
			icIndx = pList[i];		//�e���q�ɃN���X�^���p�ӂ���Ă���̂�pIndx->cIndx�Ƃ��Ă���
	
			for(j = 0; j < GetCtoPIndx(icIndx); j++)
			{
				jpIndx = GetCtoP(icIndx, j, 0);		//�Z���N���X�^�Ɋ܂܂�Ă������q
													//���̗��q����C�Z�������N���X�^�̏�����菜��
				
				if(jpIndx == -1){	continue;	}

				for(k = 0; k < GetPtoCIndx(jpIndx); k++)
				{					
					if(GetPtoC(jpIndx, k, 0) == -1
					|| GetPtoC(jpIndx, k, 1) == -1
					|| GetPtoC(jpIndx, k, 0) != icIndx)
					{
						continue;
					}
	
					#pragma omp critical (DeletePtoC)	//TODO�F�F��ɃJ�E���g�����ق������񉻂ł��Ă悢
					{
						DeletePtoC(jpIndx, k);
					}

					break;			//�����N���X�^�ɕ����������邱�Ƃ͖����̂ŁCbreak
				}
			}
		}
	}//end #pragma omp parallel

	//�Z�𗱎q���N���X�^�����폜
	#pragma omp parallel
	{
	#pragma omp for
		for(int i = 0; i < pListSize; i++)
		{
			ClearPtoC(pList[i]);
		}
	}//end #pragma omp parallel

	//�Z���N���X�^�����q�����폜
	#pragma omp parallel
	{
	#pragma omp for
		for(int i = 0; i < pListSize; i++)
		{
			ClearCtoP(pList[i]);
			//m_iceObj->GetMoveObj(pList[i])->Clear();
		}
	}//end #pragma omp parallel

	//�Ē�`�N���X�^�Ɋ܂܂�闱�q����A�Ē�`�N���X�^�̏�������
	#pragma omp parallel
	{
	#pragma omp for private(j,k,icIndx,jpIndx)
		for(int i = 0; i < cListSize; i++)
		{
			icIndx = cList[i];
	
			for(j = 0; j < GetCtoPIndx(icIndx); j++)
			{
				jpIndx = GetCtoP(icIndx, j, 0);
				
				if(jpIndx == -1){	continue;	}

				for(k = 0; k < GetPtoCIndx(jpIndx); k++)
				{
					if(GetPtoC(jpIndx, k, 0) == -1
					|| GetPtoC(jpIndx, k, 1) == -1
					|| GetPtoC(jpIndx, k, 0) != icIndx)
					{
						continue;
					}

					#pragma omp critical (DeletePtoC)	//TODO�F�F��ɃJ�E���g�����ق������񉻂ł��Ă悢
					{
						DeletePtoC(jpIndx, k);
					}

					break;			//�����N���X�^�ɕ����������邱�Ƃ͖����̂ŁCbreak
				}
			}
		}
	}//end #pragma omp parallel
}

void IceStructure::SetInfo_Tetra(const vector<unsigned>& pList, const vector<unsigned>& tList, const vector<unsigned>& lList)
{
	int pListSize = pList.size();
	int tListSize = tList.size();

	if(pListSize == 0 || tListSize == 0){	return; }

	int itIndx = 0;
	int ilayer = 0;

	#pragma omp parallel
	{
	#pragma omp for private(itIndx,ilayer)
		//�ߖT�l�ʑ̂̍Ē�`
		for(int i = 0; i < tListSize; i++)
		{
			itIndx = tList[i];
			ilayer = lList[i];

			ClearNeighborTetraFromLayer(itIndx, ilayer);
			SetNeighborTetraFromLayer(itIndx, m_iLayer, ilayer);	//���������ɏd��
		}
	}//end #pragma omp parallel
}


//------------------------------------------���ω��Q-------------------------------------------------




//-------------------------------------------�擾----------------------------------------
/*!
 * ���q���l�ʑ́@�z��̋󂫏ꏊ��T���ēY������Ԃ�
 */
int IceStructure::GetPtoTFreeIndx(int pIndx)
{
	int freeIndx = -1;

	for(int i = 0; i < m_iPtoTMax; i++)
	{
		if(GetPtoT(pIndx, i, 0) != -1 || GetPtoT(pIndx, i, 1) != -1){	continue;	}
		freeIndx = i;	break;
	}

	if(freeIndx == m_iPtoTMax || freeIndx == -1)
	{	
		cout << __FUNCTION__ << " Error::�z��ɋ󂫂�����܂��� " << endl;
		freeIndx = 0;
	}

	return freeIndx;
}

/*!
 * ���q���N���X�^�@�z��̋󂫏ꏊ��T���ēY������Ԃ�
 */
int IceStructure::GetPtoCFreeIndx(int pIndx)
{
	int freeIndx = -1;

	for(int i = 0; i < m_iPtoCMax; i++)
	{
		//if(GetPtoC(pIndx, i)[0] != -1 || GetPtoC(pIndx, i)[1] != -1){	continue;	}
		if(GetPtoC(pIndx, i, 0) != -1 || GetPtoC(pIndx, i, 1) != -1 || GetPtoC(pIndx, i, 2) != -1)
		{
			continue;
		}

		freeIndx = i;	break;
	}

	if(freeIndx == m_iPtoCMax || freeIndx == -1)
	{	
		cout << __FUNCTION__ << " Error::�z��ɋ󂫂�����܂��� " << endl;
		freeIndx = 0;
	}

	return freeIndx;
}

//-------------------------------------------�擾----------------------------------------

//-------------------------------------------��������----------------------------------------
/*!
 * �o�^�����@���q���l�ʑ́@TODO::��������-1�������Ė��߂Ă������ق�������
 * @param[in] pIndx�@���q�ԍ�
 * @param[in] lIndx�@���q��������l�Ԗڂ̃N���X�^
 * @param[in] cIndx�@�N���X�^�ԍ�
 * @param[in] oIndx�@�N���X�^���ł̏����ԍ�
 */
void IceStructure::SetPtoT(int pIndx, int lIndx, int tIndx, int oIndx)
{
	//�G���[�`�F�b�N
	if(m_iPtoTMax <= lIndx)
	{
		cout << __FUNCTION__ << " Error::���q��������l�ʑ̏����i�[�ł��Ȃ��Ȃ�܂����D������������܂���" << endl;
		cout << m_iPtoTMax << "<" << lIndx << endl;
		return;
	}
	else if(lIndx < 0)
	{
		cout << __FUNCTION__ << " Error::�Y�������s���Ȓl�ł��D" << endl;
		cout << "lIndx = " << lIndx << endl;
	}

	m_mk3DiPtoT(pIndx, lIndx, 0) = tIndx;
	m_mk3DiPtoT(pIndx, lIndx, 1) = oIndx;
}

/*!
 * �o�^�����@�l�ʑ́����q
 * @param[in] tIndx�@�@�@�l�ʑ̔ԍ�
 * @param[in] pIndxList�@�������q�z��
 */
void IceStructure::SetTtoP(int tIndx, vector<int>& pIndxList)
{
	//�G���[�`�F�b�N
	if(4 < pIndxList.size())
	{
		cout << __FUNCTION__ << " Error::�l�ʑ̂��܂ޗ��q�̏����i�[�ł��Ȃ��Ȃ�܂����D������������܂���" << endl;
		cout << 4 << "<" << pIndxList.size() << endl;
		return;
	}

	for(int i = 0; i < pIndxList.size(); i++)
	{
		m_mk2DiTtoP(tIndx, i) = pIndxList[i];
	}
}

/*!
 * �o�^�����@���q���N���X�^
 * @param[in] pIndx�@���q�ԍ�
 * @param[in] lIndx�@���q��������l�Ԗڂ̃N���X�^
 * @param[in] cIndx�@�N���X�^�ԍ�
 * @param[in] oIndx�@�N���X�^���ł̏����ԍ�
 * @param[in] layer�@�N���X�^���ł̑w��
 */
void IceStructure::SetPtoC(int pIndx, int lIndx, int cIndx, int oIndx, int layer)
{
	//�G���[�`�F�b�N
	if(m_iPtoCMax <= lIndx)
	{
		cout << __FUNCTION__ << " Error::���q��������N���X�^�����i�[�ł��Ȃ��Ȃ�܂����D������������܂���" << endl;
		cout << m_iPtoCMax << "<" << lIndx << endl;
		return;
	}
	
	m_mk3DiPtoC(pIndx, lIndx, 0) = cIndx;
	m_mk3DiPtoC(pIndx, lIndx, 1) = oIndx;
	m_mk3DiPtoC(pIndx, lIndx, 2) = layer;
}

/*!
 * �o�^�����@�N���X�^�����q
 * @param[in] cIndx�@�@�@�N���X�^�ԍ�
 * @param[in] pIndxList�@�������q�z��
 */
void IceStructure::SetCtoP(int cIndx, const vector<int>& pIndxList, int* pLayerList)
{
	//�G���[�`�F�b�N
	if(m_iCtoPMax <= pIndxList.size())
	{
		cout << __FUNCTION__ << " Error::�N���X�^���܂ޗ��q�̏����i�[�ł��Ȃ��Ȃ�܂����D������������܂���" << endl;
		cout << m_iCtoPMax << "<" << pIndxList.size() << endl;
		return;
	}

	for(int i = 0; i < pIndxList.size(); i++)
	{
		m_mk3DiCtoP(cIndx, i, 0) = pIndxList[i];
		m_mk3DiCtoP(cIndx, i, 1) = pLayerList[i];
	}
}

/*!
 * �o�^�����@�����̋ߖT�l�ʑ�
 * @param[in] tIndx�@�@�l�ʑ̔ԍ�
 * @param[in] PtoTNum�@
 */
void IceStructure::SetTetraInfo(int tIndx, int* PtoTNum)
{
	IceTetrahedra &tetra = IceTetrahedra::GetInstance();		//����܂肱���ł͌Ăт����Ȃ�����

	//���q�������Ă���l�ʑ̂̔ԍ���o�^���邽�߂̏���
	//pCountList�ɂ́CtIndx�Ԗڂ̎l�ʑ̂Ɋ܂܂��e���q���C���ꂼ�ꂢ���̎l�ʑ̂ɑ����邩�����߂ĕۑ�����
	int* pCountList = new int[tetra.GetTetraList(tIndx).size()];
	//cout << __FUNCTION__ << "::check0" << endl;
	for(int j = 0; j < tetra.GetTetraList(tIndx).size(); j++)
	{
		int pIndx = tetra.GetTetraList(tIndx)[j];
		pCountList[j] = GetPtoTNum(pIndx)-PtoTNum[pIndx];
		PtoTNum[pIndx]--;
	}
	//cout << __FUNCTION__ << "::check1" << endl;
	//���q�Ǝl�ʑ̂̏��o�^
	vector<int>& pIndxList = tetra.GetTetraList(tIndx);

	for(int i = 0; i < GetTtoPNum(tIndx); i++)
	{
		SetPtoT(pIndxList[i], pCountList[i], tIndx, i);
	}
	//cout << __FUNCTION__ << "::check2" << endl;

	SetTtoP(tIndx, pIndxList);
	//cout << __FUNCTION__ << "::check3" << endl;

	delete[] pCountList;
}


/*!
 * �o�^�����@�����̋ߖT�l�ʑ�
 * @param[in] tIndx�@�l�ʑ̔ԍ�
 * @param[in] layer�@�T���K�w
 */
void IceStructure::SetNeighborTetra(int tIndx, int layer)
{
	/*	
		����l�ʑ�A�Ɋ܂܂�Ă��闱�q���C���̕����̎l�ʑ̂Ɋ܂܂�Ă���ꍇ�C
		���̕����̎l�ʑ͎̂l�ʑ�A�̋ߖT�l�ʑ̂ƂȂ�
	*/
	vector<int> pIndxList;

	//layer=1�w��
	for(int i = 0; i < GetTtoPIndx(tIndx); i++)
	{
		int ipIndx = GetTtoP(tIndx, i);
		if(ipIndx == -1){	continue;	}

		if(find(pIndxList.begin(), pIndxList.end(), ipIndx) != pIndxList.end())
		{	
			continue;
		}
		pIndxList.push_back(ipIndx);

		//�T��
		for(int j = 0; j < GetPtoTIndx(ipIndx); j++)
		{
			if(GetPtoT(ipIndx, j, 0) == -1
			|| GetPtoT(ipIndx, j, 1) == -1
			|| GetPtoT(ipIndx, j, 0) == tIndx)
			{	
				continue;	
			}

			//�����N���X�^�����Ɋ܂�ł��Ȃ����̃`�F�b�N
			if(CheckNeighborTetra( tIndx, GetPtoT(ipIndx, j, 0) ) != -1){	continue;	}

			//�ߖT�N���X�^�̓o�^�{�J�E���g
			if(GetNTNum(tIndx) >= m_iNeighborMax)
			{
				cout << __FUNCTION__ << "�ߖT���q�̃�����������Ȃ��Ȃ�܂���" << endl;
				cout << GetNTNum(tIndx) << " < " << m_iNeighborMax << endl;
			}
			else
			{
				m_mk3DiNeighborTetra(tIndx, GetNTNum(tIndx), 0) = GetPtoT(ipIndx, j, 0);
				m_mk3DiNeighborTetra(tIndx, GetNTNum(tIndx), 1) = 1;
				CountNT(tIndx);
			}
		}
	}

	//layer�w�ڂ����ǂ��ċߖT�l�ʑ̂��擾����
	int nowSize = 0, nowIndx = 0;
	int d_addNum = 0;
	for(int i = 2; i <= layer; i++)
	{
		nowSize = GetNTNum(tIndx);
		d_addNum = GetNTNum(tIndx);

		//�T������N���X�^��nowSize��nowIndx�Ő������Ă���
		for(int j = nowIndx; j < nowSize; j++)
		{
			int jtIndx = GetNeighborTetra(tIndx, j, 0);			//�ߖT�l�ʑ̂̂ЂƂ�
			
			//�ߖT�l�ʑ̂Ɋ܂܂�闱�q�����̎l�ʑ̂ɂ��܂܂�Ă���ꍇ�C���̎l�ʑ̂��ߖT�Ƃ��ēo�^
			for(int k = 0; k < GetTtoPIndx(jtIndx); k++)
			{
				int kpIndx = GetTtoP(jtIndx, k);				//�ߖT�l�ʑ̂Ɋ܂܂�Ă��闱�q�̂ЂƂ�
				if(kpIndx == -1){ continue;	}

				if(find(pIndxList.begin(), pIndxList.end(), kpIndx) != pIndxList.end())
				{
					continue;
				}
				pIndxList.push_back(kpIndx);

				for(int l = 0; l < GetPtoTIndx(kpIndx); l++)
				{
					if(GetPtoT(kpIndx, l, 0) == -1 
					|| GetPtoT(kpIndx, l, 1) == -1
					|| GetPtoT(kpIndx, l, 0) == tIndx)
					{
						continue;
					}

					//�����l�ʑ̂����Ɋ܂�ł��Ȃ����̃`�F�b�N
					if(CheckNeighborTetra( tIndx, GetPtoT(kpIndx, l, 0) ) != -1)
					{	//cout << "check4" << endl;
						continue;
					}

					//�ߖT�N���X�^�̓o�^�{�J�E���g
					if(GetNTNum(tIndx) >= m_iNeighborMax)
					{
						cout << __FUNCTION__ << " �ߖT�l�ʑ̂̃�����������Ȃ��Ȃ�܂��� " << GetNTNum(tIndx) << endl;
						cout << GetNTNum(tIndx) << " < " << m_iNeighborMax << endl;
					}
					else
					{
//						if(GetNTNum(tIndx) > 200 ) return;	//�ߖT�l�ʑ̐��̐���
						m_mk3DiNeighborTetra(tIndx, GetNTNum(tIndx), 0) = GetPtoT(kpIndx, l, 0);
						m_mk3DiNeighborTetra(tIndx, GetNTNum(tIndx), 1) = i;
						CountNT(tIndx);
					}
				}
			}
		}
		nowIndx = nowSize;										//���̃��[�v�J�n���̃X�^�[�g�ԍ����X�V
		//cout << "tIndx = " << tIndx << " addNum = " << GetNTNum(tIndx)-d_addNum << endl;
	}
}

/*!
 * �o�^�����@�����̋ߖT�l�ʑ�
 * @param[in] tIndx�@�l�ʑ̔ԍ�
 * @param[in] searchLayer�@�T���I���K�w searchLayer >= 1
 * @param[in] deleteLayer�@�T���J�n�K�w deleteLayer >= 1
 */
void IceStructure::SetNeighborTetraFromLayer(int tIndx, int searchLayer, int deleteLayer)
{//	cout << __FUNCTION__ << endl;
	/*	
		����l�ʑ�A�Ɋ܂܂�Ă��闱�q���C���̕����̎l�ʑ̂Ɋ܂܂�Ă���ꍇ�C
		���̕����̎l�ʑ͎̂l�ʑ�A�̋ߖT�l�ʑ̂ƂȂ�
	*/
	vector<int> pIndxList;		//�T���ςݗ��q��ۑ��@����̂������ł��Ȃ葁���Ȃ�D

	//layer=1�w��
	for(int i = 0; i < GetTtoPIndx(tIndx); i++)
	{
		int ipIndx = GetTtoP(tIndx, i);
		if(ipIndx == -1){	continue;	}

		if(find(pIndxList.begin(), pIndxList.end(), ipIndx) != pIndxList.end()){	continue;		}
		pIndxList.push_back(ipIndx);

		if(deleteLayer != 1){ continue;	}		//1�w�ڂ̂ݍs��

		//�T��
		for(int j = 0; j < GetPtoTIndx(ipIndx); j++)
		{
			if(GetPtoT(ipIndx, j, 0) == -1
			|| GetPtoT(ipIndx, j, 1) == -1
			|| GetPtoT(ipIndx, j, 0) == tIndx)
			{
				continue;
			}

			//�����N���X�^�����Ɋ܂�ł��Ȃ����̃`�F�b�N
			if(CheckNeighborTetra(tIndx, GetPtoT(ipIndx, j, 0)) != -1){	continue;	}

			//�ߖT�N���X�^�̓o�^�{�J�E���g
			if(GetNTNum(tIndx) >= m_iNeighborMax)
			{
				cout << __FUNCTION__ << "�ߖT���q�̃�����������Ȃ��Ȃ�܂���" << endl;
				cout << GetNTNum(tIndx) << " < " << m_iNeighborMax << endl;
			}
			else
			{
				m_mk3DiNeighborTetra(tIndx, GetNTNum(tIndx), 0) = GetPtoT(ipIndx, j, 0);
				m_mk3DiNeighborTetra(tIndx, GetNTNum(tIndx), 1) = 1;
				CountNT(tIndx);
			}
		}
	}

	if(deleteLayer == 1)
	{
		deleteLayer += 1;			//�S�ď���������ꍇ��+1����
	}
		
	//layer�w�ڂ����ǂ��ċߖT�l�ʑ̂��擾����
	int nowSize = 0, nowIndx = 0;

	for(int i = 2; i <= searchLayer; i++)
	{
		nowSize = GetNTNum(tIndx);

		//�T������N���X�^��nowSize��nowIndx�Ő������Ă���
		for(int j = nowIndx; j < nowSize; j++)
		{
			int jtIndx = GetNeighborTetra(tIndx, j, 0);				//�ߖT�l�ʑ̂̂ЂƂ�
			
			//�ߖT�l�ʑ̂Ɋ܂܂�闱�q�����̎l�ʑ̂ɂ��܂܂�Ă���ꍇ�C���̎l�ʑ̂��ߖT�Ƃ��ēo�^
			//TODO::�l�ʑ́����q���ߖT�l�ʑ́@�ł͂Ȃ��C�l�ʑ́��ߖT�l�ʑ́@�Ƃ���
			//TODO::A��B�̋ߖT�ł���Ȃ�CB��A�̋ߖT�ł���@�𗘗p����
			//TODO::�Z���������q���������Ă����l�ʑ̂��番�����邩���t�Z����D
			//TODO::�č\�����Ȃ��l�ʑ̂𗘗p����
			for(int k = 0; k < GetTtoPIndx(jtIndx); k++)
			{
				int kpIndx = GetTtoP(jtIndx, k);					//�ߖT�l�ʑ̂Ɋ܂܂�Ă��闱�q�̂ЂƂ�
				if(kpIndx == -1){	continue;	}
				if(find(pIndxList.begin(), pIndxList.end(), kpIndx) != pIndxList.end()){	continue;}
				pIndxList.push_back(kpIndx);

				if(i < deleteLayer){	continue;	}				//�ےT������K�v�Ȃ��Ȃ�C���q��ǉ����������Ŗ߂�

				//���q���������Ă���ߖT�l�ʑ̂�T��
				for(int l = 0; l < GetPtoTIndx(kpIndx); l++)
				{
					if(GetPtoT(kpIndx, l, 0) == -1
					|| GetPtoT(kpIndx, l, 1) == -1
					|| GetPtoT(kpIndx, l, 0) == tIndx)
					{
						continue;
					}

					if(CheckNeighborTetra( tIndx, GetPtoT(kpIndx, l, 0) ) != -1){	continue;	}	//�����l�ʑ̂����Ɋ܂�ł��Ȃ����̃`�F�b�N

					//�ߖT�N���X�^�̓o�^�{�J�E���g
					if(GetNTNum(tIndx) >= m_iNeighborMax)
					{
						cout << __FUNCTION__ << " �ߖT�l�ʑ̂̃�����������Ȃ��Ȃ�܂��� " << endl;
						cout << GetNTNum(tIndx) << " < " << m_iNeighborMax << endl;
					}
					else
					{
//						if(GetNTNum(tIndx) > 200 ) return;	//�ߖT�l�ʑ̐��̐���
						m_mk3DiNeighborTetra(tIndx, GetNTNum(tIndx), 0) = GetPtoT(kpIndx, l, 0);
						m_mk3DiNeighborTetra(tIndx, GetNTNum(tIndx), 1) = i;
						CountNT(tIndx);
					}
				}
			}
		}
		nowIndx = nowSize;											//���̃��[�v�J�n���̃X�^�[�g�ԍ����X�V
	}
}

//-------------------------------------------��������----------------------------------------

/*!
 * �擾�����@���q���l�ʑ�
 * @param[in] pIndx�@���q�ԍ�
 * @param[in] lIndx�@���q���ԍ�
 * @param[in] oIndx	0->�N���X�^�̔ԍ�, 1->�N���X�^���ł̗��q�ԍ�, 2->�K�w
 */

int IceStructure::GetPtoT(int pIndx, int lIndx, int oIndx)
{
	return m_mk3DiPtoT(pIndx, lIndx, oIndx);
}

/*!
 * �擾�����@�l�ʑ́����q
 * @param[in] tIndx�@�l�ʑ̔ԍ�
 * @param[in] lIndx�@�l�ʑ̓��ԍ�
 */
int IceStructure::GetTtoP(int tIndx, int lIndx)
{
	return m_mk2DiTtoP(tIndx, lIndx);
}

/*!
 * �擾�����@���q���N���X�^
 * @param[in] pIndx	���q�ԍ�
 * @param[in] lIndx	���q���ԍ�
 * @param[in] oIndx	0->�N���X�^�̔ԍ�, 1->�N���X�^���ł̗��q�ԍ�, 2->�K�w
 */
int IceStructure::GetPtoC(int pIndx, int lIndx, int oIndx)
{
	return m_mk3DiPtoC(pIndx, lIndx, oIndx);
}

/*!
 * �擾�����@�N���X�^�����q
 * @param[in] cIndx �N���X�^�ԍ�
 * @param[in] lIndx �N���X�^�����q�ԍ�
 * @param[in] oIndx 
 */
const int& IceStructure::GetCtoP(const int& cIndx, const int& lIndx, const int& oIndx)
{
	return m_mk3DiCtoP(cIndx, lIndx, oIndx);
}

/*!
 * �擾�����@�l�ʑ́��ߖT�l�ʑ�
 * @param[in] tIndx�@�l�ʑ̔ԍ�
 * @param[in] lIndx�@�l�ʑ̓��ԍ�
 * @param[in] oIndx 
 */
int IceStructure::GetNeighborTetra(int tIndx, int lIndx, int oIndx)
{
	return m_mk3DiNeighborTetra(tIndx, lIndx, oIndx);
}

/*!
 * �l�ʑ̂Ɋ܂܂�Ă��闱�q���̃J�E���g�C���q����������l�ʑ̐��̃J�E���g
 * @param[in] tIndx�@�l�ʑ̔ԍ�
 * @param[in] pList�@�l�ʑ̂Ɋ܂܂�闱�q���X�g
 */
void IceStructure::CountTetrahedra(int tIndx, vector<int>& pList)
{
	for(unsigned i = 0; i < pList.size(); i++)
	{
		int pIndx = pList[i];		

		CountPtoT(pIndx);
		CountTtoP(tIndx);

		//Indx�̍X�V
		if(GetPtoTNum(pIndx) >= GetPtoTIndx(pIndx))
		{
			SetPtoTIndx(pIndx, GetPtoTNum(pIndx));
		}
	}
	
	//Indx�̍X�V
	if(GetTtoPNum(tIndx) >= GetTtoPIndx(tIndx))
	{
		SetTtoPIndx(tIndx, GetTtoPNum(tIndx));
	}
}

/*!
 * �N���X�^�Ɋ܂܂�Ă��闱�q���̃J�E���g�C���q����������l�ʑ̐��̃J�E���g
 * @param[in] cIndx�@���q�ԍ�
 * @param[in] pList�@�N���X�^�Ɋ܂܂�闱�q���X�g
 */
void IceStructure::CountClusterParticle(int cIndx, vector<int>& pList, int pNum)
{
	for(int j = 0; j < pNum; j++)
	{
		int jpIndx = pList[j];
		CountPtoC(jpIndx);										//���q���ڑ��N���X�^�ɏ���������̃J�E���g
		CountCtoP(cIndx);										//�ڑ��N���X�^�����q���܂ތ��̃J�E���g
	
		//Indx�̍X�V
		if(GetPtoCNum(jpIndx) >= GetPtoCIndx(jpIndx))
		{
			SetPtoCIndx(jpIndx, GetPtoCNum(jpIndx));
		}
	}

	//Indx�̍X�V
	if(GetCtoPNum(cIndx) >= GetCtoPIndx(cIndx))
	{
		SetCtoPIndx(cIndx, GetCtoPNum(cIndx));
	}
}

/*!
 * �폜�����@�l�ʑ́����q
 * @param[in] tIndx�@�l�ʑ̔ԍ�
 * @param[in] lIndx�@�l�ʑ̓��ԍ�
 */
void IceStructure::DeleteTtoP(int tIndx, int lIndx)
{//	cout << __FUNCTION__ << " check1" << endl;
	//����ClIndx�͓Y����
	if(4 <= lIndx)
	{
		cout << __FUNCTION__ << " Error:: 4 < lIndx �l�ʑ̔z��ւ̃A�N�Z�X�G���[" << endl;
		cout << 4 << "<" << lIndx << endl;
		return;
	}

	if(m_piTtoPNum[tIndx] > 0)
	{
		m_piTtoPNum[tIndx]--;
	}
	else
	{
		cout << __FUNCTION__ << " Error::m_piTtoPNum[" << tIndx << "] < 0" << endl;
		return;
	}

	m_mk2DiTtoP(tIndx, lIndx) = -1;
}

/*!
 * �폜�����@���q���N���X�^
 * @param[in] pIndx�@���q�ԍ�
 * @param[in] lIndx�@���q���ԍ�
 */
void IceStructure::DeletePtoC(int pIndx, int lIndx)
{//	cout << __FUNCTION__ << " check1" << endl;
	//����ClIndx�͓Y����
	if(m_iPtoCMax <= lIndx)
	{
		cout << __FUNCTION__ << " Error:: m_iPtoCMax < lIndx" << endl;
		cout << m_iPtoCMax << "<" << lIndx << endl;
		return;
	}

	if(m_piPtoCNum[pIndx] > 0)
	{
		m_piPtoCNum[pIndx]--;
	}
	else
	{
		cout << __FUNCTION__ << " Error::m_piPtoCNum[" << pIndx << "] < 0" << endl;
		return;
	}

	//m_pppiPtoC[pIndx][lIndx][0] = -1;
	//m_pppiPtoC[pIndx][lIndx][1] = -1;
	//m_pppiPtoC[pIndx][lIndx][2] = -1;

	for(int i = 0; i < 3; i++)
	{
		m_mk3DiPtoC(pIndx, lIndx, i) = -1;
	}
}

/*!
 * �폜�����@���q���l�ʑ�
 * @param[in] pIndx�@���q�ԍ�
 * @param[in] lIndx�@���q���ԍ�
 */
void IceStructure::DeletePtoT(int pIndx, int lIndx)
{//	cout << __FUNCTION__ << " check1" << endl;
	//����ClIndx�͓Y����
	if(m_iPtoTMax <= lIndx)
	{
		cout << __FUNCTION__ << " Error:: m_iPtoTMax < lIndx" << endl;
		cout << m_iPtoTMax << "<" << lIndx << endl;
		return;
	}

	if(m_piPtoTNum[pIndx] > 0)
	{
		m_piPtoTNum[pIndx]--;
	}
	else
	{
		cout << __FUNCTION__ << " Error::m_piPtoTNum[" << pIndx << "] < 0, lIndx = " << lIndx << endl;
		return;
	}

	for(int i = 0; i < 2; i++)
	{
		m_mk3DiPtoT(pIndx, lIndx, i) = -1;
	}
}

/*!
 * �����������@���q���l�ʑ�
 */
void IceStructure::ClearPtoT(int pIndx)
{
	for(int i = 0; i < m_iPtoTMax; i++)
	{
		for(int j = 0; j < 2; j++)
		{
			m_mk3DiPtoT(pIndx, i, j) = -1;
		}
	}

	m_piPtoTIndx[pIndx] = 0;
	m_piPtoTNum[pIndx] = 0;
}

/*!
 * �����������@���q���N���X�^
 */
void IceStructure::ClearPtoC(int pIndx)
{
	for(int i = 0; i < m_iPtoCMax; i++)
	{
		for(int j = 0; j < 3; j++)
		{
			m_mk3DiPtoC(pIndx, i, j) = -1;
		}
	}

	m_piPtoCIndx[pIndx] = 0;
	m_piPtoCNum[pIndx] = 0;
}

/*!
 * �����������@�N���X�^�����q
 */
void IceStructure::ClearCtoP(int cIndx)
{
	for(int i = 0; i < m_piCtoPIndx[cIndx]; i++)
	{
		m_mk3DiCtoP(cIndx, i, 0) = -1;
		m_mk3DiCtoP(cIndx, i, 1) = -1;
	}

	m_piCtoPIndx[cIndx] = 0;
	m_piCtoPNum[cIndx] = 0;
}

/*!
 * �����������@�l�ʑ́����q
 */
void IceStructure::ClearTtoP(int tIndx)
{
	for(int i = 0; i < m_piTtoPIndx[tIndx]; i++)
	{
		m_mk2DiTtoP(tIndx, i) = -1;
	}

	m_piTtoPIndx[tIndx] = 0;
	m_piTtoPNum[tIndx] = 0;

	ClearNeighborTetra(tIndx);		//�ߖT�l�ʑ̂�������
}

/*!
 * �����������@�ߖT�l�ʑ�
 */
void IceStructure::ClearNeighborTetra(int tIndx)
{//	cout << __FUNCTION__ << endl;
	for(int i = 0; i < m_piNTNum[tIndx]; i++)
	{
		m_mk3DiNeighborTetra(tIndx, i, 0) = -1;
		m_mk3DiNeighborTetra(tIndx, i, 1) = -1;
	}

	m_piNTNum[tIndx] = 0;
}

/*!
 * �����I�����������@�ߖT�l�ʑ�
 *@paramm[in] layer �폜����w layer >= 1
 */
void IceStructure::ClearNeighborTetraFromLayer(int tIndx, int layer)
{
	for(int i = m_piNTNum[tIndx]-1; 0 <=i ; i--)
	{
		if(m_mk3DiNeighborTetra(tIndx, i, 1) >= layer)
		{
			m_mk3DiNeighborTetra(tIndx, i, 0) = -1;
			m_mk3DiNeighborTetra(tIndx, i, 1) = -1;
			m_piNTNum[tIndx]--;
		}
		else
		{
			break;
		}
	}
}

/*!
 * ���菈���@����l�ʑ̂��ߖT�l�ʑ̂Ƃ��Ċ܂܂�Ă��邩�̃`�F�b�N
 * @param[in] tIndx�@�@�@�l�ʑ̔ԍ�
 * @param[in] checkTIndx�@�m�F����l�ʑ̔ԍ�
 */
int IceStructure::CheckNeighborTetra(int tIndx, int checkTIndx)
{
	int findIndx = -1;
	
	for(int k = 0; k < GetNTNum(tIndx); k++)
	{
		if(checkTIndx == GetNeighborTetra(tIndx, k, 0))
		{
			findIndx = k; break;
		}
	}

	return findIndx;
}

//-------------------------------------�X�V----------------------------------------

//-------------------------------------�I��I�^���v�Z------------------------------------------

void IceStructure::UpdateMotionCalcCluster(unsigned cIndx, short unsigned num)
{
	//cout << __FUNCTION__ << ", cIndx = " << cIndx << ", num = " << num << endl;
	m_psuSelectClusterIndx[cIndx] = num;
}

short unsigned IceStructure::GetMotionCalcCluster(unsigned cIndx)
{
	return m_psuSelectClusterIndx[cIndx];
}


//-------------------------------------�f�o�b�O----------------------------------------
void IceStructure::DebugPtoT(int pIndx)
{	cout << __FUNCTION__ << " pIndx = " << pIndx;
	cout << " num=" << GetPtoTNum(pIndx) << " Indx=" << GetPtoTIndx(pIndx);
	
	for(int i = 0; i < GetPtoTIndx(pIndx); i++)
	{
		cout << " c=" << GetPtoT(pIndx, i, 0) << " o=" << GetPtoT(pIndx, i, 1);
	}
	cout << endl;
}

void IceStructure::DebugPtoC(int pIndx)
{	cout << __FUNCTION__ << " pIndx = " << pIndx;
	cout << " num=" << GetPtoCNum(pIndx) << " Indx=" << GetPtoCIndx(pIndx);
	
	for(int i = 0; i < GetPtoCIndx(pIndx); i++)
	{
		cout << " c=" << GetPtoC(pIndx, i, 0) << " o=" << GetPtoC(pIndx, i, 1);
	}
	cout << endl;

}

void IceStructure::DebugCtoP(int cIndx)
{	cout << __FUNCTION__ << " cIndx = " << cIndx;
	cout << " num=" << GetCtoPNum(cIndx) << " Indx=" << GetCtoPIndx(cIndx);

	for(int i = 0; i < GetCtoPIndx(cIndx); i++)
	{
		cout <<" p=" << GetCtoP(cIndx, i, 0) << " l=" << GetCtoP(cIndx, i, 1);
	}
	cout << endl;
}

void IceStructure::DebugTtoP(int tIndx)
{	cout << __FUNCTION__ << " tIndx=" << tIndx;
	cout << " num=" << GetTtoPNum(tIndx) << " Indx=" << GetTtoPIndx(tIndx);

	for(int i = 0; i < GetTtoPIndx(tIndx); i++)
	{
		cout << " " << GetTtoP(tIndx, i);
	}
	cout << endl;
}

void IceStructure::DebugNeighborTetra(int tIndx)
{	cout << __FUNCTION__ << " tIndx=" << tIndx;
	cout << " num = " << GetNTNum(tIndx);
	for(unsigned j = 0; j < GetNTNum(tIndx); j++)
	{
		cout << " NC=" << GetNeighborTetra(tIndx, j, 0) << " Ly=" << GetNeighborTetra(tIndx, j, 1);
	}
	cout << endl;
}

void IceStructure::DebugStepObjMelt(vector<unsigned>& pList, vector<unsigned>& cList)
{	cout << __FUNCTION__ << endl;
	if(pList.size() == 0){	return;	}
	
	cout << "pList.size = " << pList.size() << " ";
	sort(pList.begin(), pList.end());

	for(unsigned i = 0; i < pList.size(); i++)
	{
		cout << " " << pList[i];
	}
	cout << endl;

	cout << "cList.size =  " << cList.size() << " ";
	sort(cList.begin(), cList.end());

	for(unsigned i = 0; i < cList.size(); i++)
	{
		cout << " " << cList[i];
	}
	cout << endl;

	//cout << "viCLayerList:: ";
	//for(unsigned i = 0; i < viCLayerList.size(); i++)
	//{
	//	cout << " " << viCLayerList[i];
	//}
	//cout << endl;

	//cout << "viTetraList.size = " << viTetraList.size() << " ";
	//for(unsigned i = 0; i < viTetraList.size(); i++)
	//{
	//	cout << " " << viTetraList[i];
	//}
	//cout << endl;

	//cout << "viTLayerList:: ";
	//for(unsigned i = 0; i < viTLayerList.size(); i++)
	//{
	//	cout << " " << viTLayerList[i];
	//}
	//cout << endl;

	////�N���X�^�����q
	//for(int i = 0; i < m_iClusteresNum; i++){	m_ice->DebugCtoP(i);	}

	//���q���N���X�^
	//for(int i = 0; i < m_pPS->GetNumParticles(); i++){	m_ice->DebugPtoC(i);	}

	//SM�N���X�^�Ɋ܂܂�闱�q�͋@�\�Ŋm�F�ł���
	//for(int i = 0; i < ICENUM; i++)
	//{
	//	cout << "pIndx = " << endl;

	//	for(int j = 0; j < m_sm_cluster[i]->GetNumVertices(); j++)
	//	{
	//		cout << " j = " << m_sm_cluster[i]->GetParticleIndx(j);
	//	}
	//	cout << endl;
	//}

	//�l�ʑ́����q�͋@�\�Ŋm�F�ł���

	//���q���l�ʑ�
	//for(int i = 0; i < m_pPS->GetNumParticles(); i++){	m_ice->DebugPtoT(i);	}

	//�ߖT�l�ʑ�
	//for(unsigned i = 0; i < m_vviTetraList.size(); i++ ){	m_ice->DebugNeighborTetra(i);	}
}

//�e�X�g
void IceStructure::TestStepObjMelt(
	vector<unsigned>& pList,
	vector<unsigned>& cList,
	vector<unsigned>& tList,
	vector<unsigned>& cLayerList,
	vector<unsigned>& tLayerList)
{
	QueryCounter counter1;
	QueryCounter counter2;
	QueryCounter counter3;
	QueryCounter counter4;
	QueryCounter counter5;

counter1.Start();
	SearchReconstruct_Tetra_Melt(pList, tList, tLayerList);			//�Ē�`�l�ʑ̂̒T��
double end1 = counter1.End();

counter2.Start();
	SearchReconstruct_Cluster_Melt(pList, cList, cLayerList);		//�Ē�`�N���X�^�̒T��
double end2 = counter2.End();

counter3.Start();
	UpdateInfo_Melt_PandT(pList);									//���q�E�l�ʑ̏��̍X�V
double end3 = counter3.End();

counter4.Start();
	UpdateInfo_Melt_PandC(pList, cList);							//���q�E�N���X�^���̍X�V
double end4 = counter4.End();

	////CheckDeleteCluster();													//����C��܊֌W�ɂ���N���X�^���폜
	////CheckDeleteTetra(viTetraList, viTLayerList);							//����C��܊֌W�ɂ���l�ʑ̂��폜

counter5.Start();
	SetInfo_Tetra(pList, tList, tLayerList);						//���q�E�ߖT�l�ʑ̏��̍Ē�`
double end5 = counter5.End();

	if(pList.size() == 0) return;
	cout << "SearchReconstruct_Tetra_Melt	:" << end1 << endl;
	cout << "SearchReconstruct_Cluster_Melt	:" << end2 << endl;
	cout << "UpdateInfo_Melt_PandT		:" << end3 << endl;
	cout << "UpdateInfo_Melt_PandC		:" << end4 << endl;
	cout << "SetInfo_Tetra			:" << end5 << endl;

//�f�o�b�O
	//DebugStepObjMelt(pList, cList);
}