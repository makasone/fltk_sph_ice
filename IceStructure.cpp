//-----------------------------------------------------------------------------
// �C���N���[�h�t�@�C��
//-----------------------------------------------------------------------------
#include "IceStructure.h"

/*!
 * @param[in] pNumMax�@�ő嗱�q��
 * @param[in] cNumMax�@�ő�N���X�^��
 * @param[in] tNumMax�@�ő�l�ʑ̐�
 */
IceStructure::IceStructure(int pNumMax, int cNumMax, int tNumMax)
{
	//�ő吔�̓o�^
	m_iPNumMax = pNumMax;
	m_iCNumMax = cNumMax;
	m_iTNumMax = tNumMax;

	//���q���̏�����
	m_iPtoCMax = m_iCNumMax*0.4;	//1331 layer2 0.4 layer3 0.75
									//2197 layer2 0.4 layer3 0.4 layer4 0.5
	m_iPtoTMax = m_iTNumMax*0.4;	//1331 layer2 0.3 layer3 0.5
									//2197 layer2 0.3 layer3 0.3

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
	m_iCtoPMax = m_iPNumMax*0.5;	//1331 layer2 0.5 layer3 0.75
									//2197 layer2 0.5 layre3 0.5
	m_piCtoPNum = new int[m_iCNumMax];
	m_piCtoTNum = new int[m_iCNumMax];

	m_piCtoPIndx = new int[m_iCNumMax];
	m_piCtoTIndx = new int[m_iCNumMax];

	for(int i = 0; i < m_iCNumMax; i++)
	{
		m_piCtoPNum[i] = 0;
		m_piCtoTNum[i] = 0;
		
		m_piCtoPIndx[i] = 0;
		m_piCtoTIndx[i] = 0;
	}

	//�l�ʑ̏��̏�����
	//TtoPMax�͍ő�S�ŌŒ�
	//TtoCMax�͕K�v�Ȃ�
	m_piTtoPNum = new int[m_iTNumMax];
	m_piTtoCNum = new int[m_iTNumMax];

	m_piTtoPIndx = new int[m_iTNumMax];
	m_piTtoCIndx = new int[m_iTNumMax];

	for(int i = 0; i < m_iTNumMax; i++)
	{
		m_piTtoPNum[i] = 0;
		m_piTtoCNum[i] = 0;
		
		m_piTtoPIndx[i] = 0;
		m_piTtoCIndx[i] = 0;
	}

	//�ߖT�l�ʑ�
	m_pppiNeighborTetra = new int**[m_iTNumMax];
	m_piNTNum = new int[m_iTNumMax];

	m_iNeighborMax = m_iTNumMax*0.3;		//1331 layer2 0.3 layer3 0.75
											//2197 layer2 0.3 layre3 0.3 layer4 0.4
											//3375 layer2
	for(int i = 0; i < m_iTNumMax; i++)
	{
		m_piNTNum[i] = 0;
		m_pppiNeighborTetra[i] = new int*[m_iNeighborMax];

		for(int j = 0; j < m_iNeighborMax; j++)
		{
			m_pppiNeighborTetra[i][j] = new int[2];
			m_pppiNeighborTetra[i][j][0] = -1;
			m_pppiNeighborTetra[i][j][1] = -1;
		}
	}

	//�t���O
	m_pbPFlag = new bool[m_iPNumMax];
	m_pbCFlag = new bool[m_iCNumMax];
	m_pbTFlag = new bool[m_iTNumMax];

	ResetPFlag(m_iPNumMax);
	ResetCFlag(m_iCNumMax);
	ResetTFlag(m_iTNumMax);
}


/*!
 * @param[in] pNum�@�ő嗱�q��
 * @param[in] cNum�@�ő�N���X�^��
 */
IceStructure::IceStructure(int pNum, int cNum)
{
	m_iPNumMax = pNum;
	m_iCNumMax = cNum;

	//���q
//	m_iPtoCMax = m_iCNumMax*0.5;
	m_iPtoCMax = m_iCNumMax*0.65;
	m_piPtoCNum_Connect = new int[m_iPNumMax];
	m_piPtoCNum_Calc	= new int[m_iPNumMax];

	m_piPtoCIndx_Connect	= new int[m_iPNumMax];
	m_piPtoCIndx_Calc		= new int[m_iPNumMax];

	//�������@���q���N���X�^�ɏ���������@�z��̓Y����
	for(int i = 0; i < m_iPNumMax; i++)
	{
		m_piPtoCNum_Connect[i] = 0;
		m_piPtoCNum_Calc[i] = 0;

		m_piPtoCIndx_Connect[i] = 0;
		m_piPtoCIndx_Calc[i] = 0;
	}

	//�N���X�^
	m_iCtoPMax = m_iPNumMax*0.75;
	m_piCtoPNum_Connect = new int[m_iCNumMax];
	m_piCtoPNum_Calc	= new int[m_iCNumMax];

	m_piCtoPIndx_Connect	= new int[m_iCNumMax];
	m_piCtoPIndx_Calc		= new int[m_iCNumMax];

	//�������@�N���X�^�Ɋ܂܂�闱�q�̌��@�z��̓Y����
	for(int i = 0; i < m_iCNumMax; i++)
	{
		m_piCtoPNum_Connect[i] = 0;
		m_piCtoPNum_Calc[i] = 0;

		m_piCtoPIndx_Connect[i] = 0;
		m_piCtoPIndx_Calc[i] = 0;
	}

	//�ߖT�N���X�^
	m_ppiNeighborCluster = new int**[m_iCNumMax];				//�ڑ��N���X�^�̋ߖT�ƂȂ�N���X�^��ۑ�
	m_ppiNCNum = new int[m_iCNumMax];							//�ڑ��N���X�^�̋ߖT�ƂȂ�N���X�^�̌�
	
//	m_iNeighborMax = m_iCNumMax * 0.35;
	m_iNeighborMax = m_iCNumMax * 0.5;							//�ߖT�N���X�^�̌�

	for(int i = 0; i < m_iCNumMax; i++)
	{
		m_ppiNCNum[i] = 0;
		m_ppiNeighborCluster[i] = new int*[m_iNeighborMax];
		for(int j = 0; j < m_iNeighborMax; j++)
		{
			m_ppiNeighborCluster[i][j] = new int[2];
			m_ppiNeighborCluster[i][j][0] = -1;
			m_ppiNeighborCluster[i][j][1] = -1;
		}
	}

	//�t���O
	m_pbPFlag = new bool[m_iPNumMax];
	m_pbCFlag = new bool[m_iCNumMax];

	ResetPFlag(m_iPNumMax);
	ResetCFlag(m_iCNumMax);
}

IceStructure::~IceStructure(void)
{
}

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

/*!
 * �ڑ����@���q�ƃN���X�^�̗̈�m��
 */
void IceStructure::InitStateConnect()
{//	cout << __FUNCTION__ << endl;
 
	//���q���N���X�^
	m_pppiPtoC_Connect = new int**[m_iPNumMax];

	for(int i = 0; i < m_iPNumMax; i++)
	{
		m_pppiPtoC_Connect[i] = new int*[m_iPtoCMax/4];

		for(int j = 0; j < m_iPtoCMax/4; j++)
		{
			m_pppiPtoC_Connect[i][j] = new int[2];				//[0] = �N���X�^�ԍ��@[1] = �N���X�^���ł̔ԍ�
			m_pppiPtoC_Connect[i][j][0] = -1;
			m_pppiPtoC_Connect[i][j][1] = -1;
		}
	}

	//�N���X�^�����q
	m_ppiCtoP_Connect = new int*[m_iCNumMax];

	for(int i = 0; i < m_iCNumMax; i++)
	{
		m_ppiCtoP_Connect[i] = new int[4];						//�S�ŌŒ�

		for(int j = 0; j < 4; j++)
		{
			m_ppiCtoP_Connect[i][j] = -1;
		}
	}

	//�Y�����̏�����
	for(int i = 0; i < m_iPNumMax; i++)
	{
		m_piPtoCIndx_Connect[i] = m_piPtoCNum_Connect[i];		//�Y�����̏�����
	}

	for(int i = 0; i < m_iCNumMax; i++)
	{
		m_piCtoPIndx_Connect[i] = m_piCtoPNum_Connect[i];		//�Y�����̏�����
	}
}


/*!
 * �v�Z���@���q�ƃN���X�^�̏�����
 */
void IceStructure::InitStateCalc()
{//	cout << __FUNCTION__ << endl;
	//���q
	m_pppiPtoC_Calc = new int**[m_iPNumMax];
	int j = 0;

	for(int i = 0; i < m_iPNumMax; i++)
	{
		m_pppiPtoC_Calc[i] = new int*[m_iPtoCMax];

		for(j = 0; j < m_iPtoCMax; j++)
		{
			m_pppiPtoC_Calc[i][j] = new int[3];
			m_pppiPtoC_Calc[i][j][0] = -1;				//[0] = �N���X�^�ԍ�
			m_pppiPtoC_Calc[i][j][1] = -1;				//[1] = �N���X�^���ł̔ԍ�
			m_pppiPtoC_Calc[i][j][2] = -1;				//[2] = layer�ԍ�
		}
	}

//	cout << __FUNCTION__ << "Particle end" << endl;

	//�N���X�^
	m_ppiCtoP_Calc = new int**[m_iCNumMax];

	for(int i = 0; i < m_iCNumMax; i++)
	{
		m_ppiCtoP_Calc[i] = new int*[m_iCtoPMax];

		for(j = 0; j < m_iCtoPMax; j++)
		{
			m_ppiCtoP_Calc[i][j] = new int[2];
			m_ppiCtoP_Calc[i][j][0] = -1;
			m_ppiCtoP_Calc[i][j][1] = -1;
		}
	}

	//�Y�����̏�����
	for(int i = 0; i < m_iPNumMax; i++)
	{
		m_piPtoCIndx_Calc[i] = m_piPtoCNum_Calc[i];		//�Y�����̏�����
	}

	for(int i = 0; i < m_iCNumMax; i++)
	{
		m_piCtoPIndx_Calc[i] = m_piCtoPNum_Calc[i];		//�Y�����̏�����
	}
}

/*!
 * �v�Z���@�l�ʑ̏��̗̈�m�ہ@���q�x�[�X����
 */
void IceStructure::InitTetraInfo()
{
	//���q���l�ʑ�
	m_pppiPtoT = new int**[m_iPNumMax];

	for(int i = 0; i < m_iPNumMax; i++)
	{
		m_pppiPtoT[i] = new int*[m_iPtoTMax];
		
		for(int j = 0; j < m_iPtoTMax; j++)
		{
			m_pppiPtoT[i][j] = new int[2];
			m_pppiPtoT[i][j][0] = -1;
			m_pppiPtoT[i][j][1] = -1;
		}
	}

	//�l�ʑ́����q
	m_ppiTtoP = new int*[m_iTNumMax];

	for(int i = 0; i < m_iTNumMax; i++)
	{
		m_ppiTtoP[i] = new int[4];
		
		for(int j = 0; j < 4; j++)
		{
			m_ppiTtoP[i][j] = -1;
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
{
	//���q���N���X�^
	m_pppiPtoC = new int**[m_iPNumMax];

	for(int i = 0; i < m_iPNumMax; i++)
	{
		m_pppiPtoC[i] = new int*[m_iPtoCMax];
		
		for(int j = 0; j < m_iPtoCMax; j++)
		{
			m_pppiPtoC[i][j] = new int[3];
			m_pppiPtoC[i][j][0] = -1;				//[0] = �N���X�^�ԍ�
			m_pppiPtoC[i][j][1] = -1;				//[1] = �N���X�^���ł̔ԍ�
			m_pppiPtoC[i][j][2] = -1;				//[2] = layer�ԍ�
		}
	}

	//�N���X�^�����q
	m_pppiCtoP = new int**[m_iCNumMax];

	for(int i = 0; i < m_iCNumMax; i++)
	{
		m_pppiCtoP[i] = new int*[m_iCtoPMax];

		for(int j = 0; j < m_iCtoPMax; j++)
		{
			m_pppiCtoP[i][j] = new int[2];
			m_pppiCtoP[i][j][0] = -1;				//[0] = ���q�ԍ�
			m_pppiCtoP[i][j][1] = -1;				//[1] = layer�ԍ�
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
}

//-------------------------------------���q�x�[�X����----------------------------------------
/*!
 * ���q���l�ʑ́@�z��̋󂫏ꏊ��T���ēY������Ԃ�
 */
int IceStructure::GetPtoTFreeIndx(int pIndx)
{
	int freeIndx = -1;

	for(int i = 0; i < m_iPtoTMax; i++)
	{
		if(GetPtoT(pIndx, i)[0] != -1 || GetPtoT(pIndx, i)[1] != -1){	continue;	}
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
		if(GetPtoC(pIndx, i)[0] != -1 || GetPtoC(pIndx, i)[1] != -1){	continue;	}
		freeIndx = i;	break;
	}

	if(freeIndx == m_iPtoCMax || freeIndx == -1)
	{	
		cout << __FUNCTION__ << " Error::�z��ɋ󂫂�����܂��� " << endl;
		freeIndx = 0;
	}

	return freeIndx;
}

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

	m_pppiPtoT[pIndx][lIndx][0] = tIndx;
	m_pppiPtoT[pIndx][lIndx][1] = oIndx;
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
		m_ppiTtoP[tIndx][i] = pIndxList[i];
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

	m_pppiPtoC[pIndx][lIndx][0] = cIndx;
	m_pppiPtoC[pIndx][lIndx][1] = oIndx;
	m_pppiPtoC[pIndx][lIndx][2] = layer;
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
		m_pppiCtoP[cIndx][i][0] = pIndxList[i];
		m_pppiCtoP[cIndx][i][1] = pLayerList[i];
	}
}

/*!
 * �o�^�����@�����̋ߖT�l�ʑ�
 * @param[in] tIndx�@�l�ʑ̔ԍ�
 * @param[in] layer�@�T���K�w
 */
void IceStructure::SetNeighborTetra(int tIndx, int layer)
{	//	cout << __FUNCTION__ << endl;
	
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
			int* jtSet = GetPtoT(ipIndx, j);

			if(jtSet[0] == -1 || jtSet[1] == -1 || jtSet[0] == tIndx){	continue;	}

			//�����N���X�^�����Ɋ܂�ł��Ȃ����̃`�F�b�N
			if(CheckNeighborTetra(tIndx, jtSet[0]) != -1){	continue;	}

			//�ߖT�N���X�^�̓o�^�{�J�E���g
			if(GetNTNum(tIndx) >= m_iNeighborMax)
			{
				cout << __FUNCTION__ << "�ߖT���q�̃�����������Ȃ��Ȃ�܂���" << endl;
				cout << GetNTNum(tIndx) << " < " << m_iNeighborMax << endl;
			}
			else
			{
				m_pppiNeighborTetra[tIndx][GetNTNum(tIndx)][0] = jtSet[0];
				m_pppiNeighborTetra[tIndx][GetNTNum(tIndx)][1] = 1;
				CountNT(tIndx);
			}
		}
	}

	//layer�w�ڂ����ǂ��ċߖT�l�ʑ̂��擾����
	int nowSize = 0, nowIndx = 0;

	for(int i = 2; i <= layer; i++)
	{
		nowSize = GetNTNum(tIndx);

		//�T������N���X�^��nowSize��nowIndx�Ő������Ă���
		for(int j = nowIndx; j < nowSize; j++)
		{
			int jtIndx = GetNeighborTetra(tIndx, j)[0];			//�ߖT�l�ʑ̂̂ЂƂ�
			
			//�ߖT�l�ʑ̂Ɋ܂܂�闱�q�����̎l�ʑ̂ɂ��܂܂�Ă���ꍇ�C���̎l�ʑ̂��ߖT�Ƃ��ēo�^
			for(int k = 0; k < GetTtoPIndx(jtIndx); k++)
			{
				int kpIndx = GetTtoP(jtIndx, k);				//�ߖT�l�ʑ̂Ɋ܂܂�Ă��闱�q�̂ЂƂ�
				if(kpIndx == -1){	continue;	}

				if(find(pIndxList.begin(), pIndxList.end(), kpIndx) != pIndxList.end())
				{	
					continue;
				}
				pIndxList.push_back(kpIndx);

				for(int l = 0; l < GetPtoTIndx(kpIndx); l++)
				{
					int* ltSet = GetPtoT(kpIndx, l);
					if(ltSet[0] == -1 || ltSet[1] == -1 || ltSet[0] == tIndx){	continue;	}

					//�����l�ʑ̂����Ɋ܂�ł��Ȃ����̃`�F�b�N
					if(CheckNeighborTetra(tIndx, ltSet[0]) != -1){	continue;	}

					//�ߖT�N���X�^�̓o�^�{�J�E���g
					if(GetNTNum(tIndx) >= m_iNeighborMax)
					{
						cout << __FUNCTION__ << " �ߖT�l�ʑ̂̃�����������Ȃ��Ȃ�܂��� " << GetNTNum(tIndx) << endl;
						cout << GetNTNum(tIndx) << " < " << m_iNeighborMax << endl;
					}
					else
					{
//						if(GetNTNum(tIndx) > 200 ) return;	//�ߖT�l�ʑ̐��̐���
						m_pppiNeighborTetra[tIndx][GetNTNum(tIndx)][0] = ltSet[0];
						m_pppiNeighborTetra[tIndx][GetNTNum(tIndx)][1] = i;
						CountNT(tIndx);
					}
				}
			}
		}
		nowIndx = nowSize;										//���̃��[�v�J�n���̃X�^�[�g�ԍ����X�V
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
			int* jtSet = GetPtoT(ipIndx, j);

			if(jtSet[0] == -1 || jtSet[1] == -1 || jtSet[0] == tIndx){	continue;	}

			//�����N���X�^�����Ɋ܂�ł��Ȃ����̃`�F�b�N
			if(CheckNeighborTetra(tIndx, jtSet[0]) != -1){	continue;	}

			//�ߖT�N���X�^�̓o�^�{�J�E���g
			if(GetNTNum(tIndx) >= m_iNeighborMax)
			{
				cout << __FUNCTION__ << "�ߖT���q�̃�����������Ȃ��Ȃ�܂���" << endl;
				cout << GetNTNum(tIndx) << " < " << m_iNeighborMax << endl;
			}
			else
			{
				m_pppiNeighborTetra[tIndx][GetNTNum(tIndx)][0] = jtSet[0];
				m_pppiNeighborTetra[tIndx][GetNTNum(tIndx)][1] = 1;
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
			int jtIndx = GetNeighborTetra(tIndx, j)[0];				//�ߖT�l�ʑ̂̂ЂƂ�
			
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
					int* ltSet = GetPtoT(kpIndx, l);
					if(ltSet[0] == -1 || ltSet[1] == -1 || ltSet[0] == tIndx){	continue;	}
					if(CheckNeighborTetra(tIndx, ltSet[0]) != -1){	continue;	}	//�����l�ʑ̂����Ɋ܂�ł��Ȃ����̃`�F�b�N

					//�ߖT�N���X�^�̓o�^�{�J�E���g
					if(GetNTNum(tIndx) >= m_iNeighborMax)
					{
						cout << __FUNCTION__ << " �ߖT�l�ʑ̂̃�����������Ȃ��Ȃ�܂��� " << endl;
						cout << GetNTNum(tIndx) << " < " << m_iNeighborMax << endl;
					}
					else
					{
//						if(GetNTNum(tIndx) > 200 ) return;	//�ߖT�l�ʑ̐��̐���
						m_pppiNeighborTetra[tIndx][GetNTNum(tIndx)][0] = ltSet[0];
						m_pppiNeighborTetra[tIndx][GetNTNum(tIndx)][1] = i;
						CountNT(tIndx);
					}
				}
			}
		}
		nowIndx = nowSize;											//���̃��[�v�J�n���̃X�^�[�g�ԍ����X�V
	}
}


/*!
 * �擾�����@���q���l�ʑ�
 * @param[in] pIndx�@���q�ԍ�
 * @param[in] lIndx�@���q���ԍ�
 */
int* IceStructure::GetPtoT(int pIndx, int lIndx)
{//	cout << __FUNCTION__ << endl;
	return m_pppiPtoT[pIndx][lIndx];
}

/*!
 * �擾�����@�l�ʑ́����q
 * @param[in] tIndx�@�l�ʑ̔ԍ�
 * @param[in] lIndx�@�l�ʑ̓��ԍ�
 */
int IceStructure::GetTtoP(int tIndx, int lIndx)
{//	cout << __FUNCTION__ << endl;
	return m_ppiTtoP[tIndx][lIndx];
}

/*!
 * �擾�����@���q���N���X�^
 * @param[in] pIndx�@���q�ԍ�
 * @param[in] lIndx�@���q���ԍ�
 */
int* IceStructure::GetPtoC(int pIndx, int lIndx)
{
	return m_pppiPtoC[pIndx][lIndx];
}

/*!
 * �擾�����@�N���X�^�����q
 * @param[in] cIndx�@���q�ԍ�
 * @param[in] lIndx�@���q���ԍ�
 */
int* IceStructure::GetCtoP(int cIndx, int lIndx)
{
	return m_pppiCtoP[cIndx][lIndx];
}

/*!
 * �擾�����@�l�ʑ́��ߖT�l�ʑ�
 * @param[in] tIndx�@�l�ʑ̔ԍ�
 * @param[in] lIndx�@�l�ʑ̓��ԍ�
 */
int* IceStructure::GetNeighborTetra(int tIndx, int lIndx)
{//	cout << __FUNCTION__ << endl;
	return m_pppiNeighborTetra[tIndx][lIndx];
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

	m_ppiTtoP[tIndx][lIndx] = -1;
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

	m_pppiPtoC[pIndx][lIndx][0] = -1;
	m_pppiPtoC[pIndx][lIndx][1] = -1;
	m_pppiPtoC[pIndx][lIndx][2] = -1;
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

	m_pppiPtoT[pIndx][lIndx][0] = -1;
	m_pppiPtoT[pIndx][lIndx][1] = -1;
}

/*!
 * �����������@���q���l�ʑ�
 */
void IceStructure::ClearPtoT(int pIndx)
{
	for(int i = 0; i < m_iPtoTMax; i++)
	{
		m_pppiPtoT[pIndx][i][0] = -1;
		m_pppiPtoT[pIndx][i][1] = -1;
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
		m_pppiPtoC[pIndx][i][0] = -1;
		m_pppiPtoC[pIndx][i][1] = -1;
		m_pppiPtoC[pIndx][i][2] = -1;
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
		m_pppiCtoP[cIndx][i][0] = -1;
		m_pppiCtoP[cIndx][i][1] = -1;
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
		m_ppiTtoP[tIndx][i] = -1;
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
		m_pppiNeighborTetra[tIndx][i][0] = -1;
		m_pppiNeighborTetra[tIndx][i][1] = -1;
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
		if(m_pppiNeighborTetra[tIndx][i][1] >= layer)
		{
			m_pppiNeighborTetra[tIndx][i][0] = -1;
			m_pppiNeighborTetra[tIndx][i][1] = -1;
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
		if(checkTIndx == GetNeighborTetra(tIndx, k)[0])
		{
			findIndx = k; break;
		}
	}

	return findIndx;
}

//-------------------------------------���q�x�[�X����----------------------------------------


//---------------------------------------�ڑ����-------------------------------------------
/*!
 * �ǉ������@���q���N���X�^
 */
void IceStructure::AddPtoC_Connect()
{

}

void IceStructure::AddPtoC_Connect(int cIndx, int oIndx)
{

}


/*!
 * �ǉ������@�N���X�^�����q
 */
void IceStructure::AddCtoP_Connect(vector<int> pIndxes)
{
}

/*!
 * �폜�����@���q���N���X�^
 */
void IceStructure::DeletePtoC_Connect(int pIndx, int coIndx)
{//	cout << __FUNCTION__ << " check1" << endl;
	if(m_pppiPtoC_Connect[pIndx][coIndx][0] == -1 || m_pppiPtoC_Connect[pIndx][coIndx][1] == -1)
	{	cout << __FUNCTION__ << " Error::-1 " << endl;
		return;
	}

	m_piPtoCNum_Connect[pIndx]--;

	if(m_piPtoCNum_Connect[pIndx] < 0)
	{
		cout << __FUNCTION__ << " Error::minus m_piPtoCNum_Connect[" << pIndx << "] = " << m_piPtoCNum_Connect[pIndx] << endl;
		m_piPtoCNum_Connect[pIndx] = 0;
		return;
	}

	//A�@-1��������@�@�����@���ʂȃ������������Ă����C����̏ꍇ�ɂ���Ă͋��������������Ȃ�
	//�Ƃ肠����-1�ɂ��������D����͕ʂ̋@��ɁD
	m_pppiPtoC_Connect[pIndx][coIndx][0] = -1;
	m_pppiPtoC_Connect[pIndx][coIndx][1] = -1;
//	cout << __FUNCTION__ << " check2" << endl;

	//B�@���l�����@�����I�@�R�s�[���������d���Ȃ�
}

/*!
 * �폜�����@�N���X�^�����q
 */
void IceStructure::DeleteCtoP_Connect(int cIndx, int oIndx)
{//	cout << __FUNCTION__ << " check1" << endl;
	//oIndx�͓Y�����ł͂Ȃ����Ȃ̂ŁC�T������oIndx�ڂ̂������������D
	int indx = 0;
	for(int i = 0; i < 4; i++)
	{
		if(m_ppiCtoP_Connect[cIndx][i] == -1){	continue;	}
		if(indx == oIndx){ indx = i;	break;	}				//oIndx�ڂ̏ꍇ�CIndx��ۑ�����break
		indx++;
	}

	if(indx == 4 || m_ppiCtoP_Connect[cIndx][indx] == -1){	return;	}		//break�Ŕ����Ȃ������ꍇ�̃`�F�b�N


	if(m_piCtoPNum_Connect[cIndx] <= 0)
	{
		cout << __FUNCTION__ << " Error::minus m_piCtoPNum_Connect[cIndx] = " << m_piCtoPNum_Connect[cIndx] << endl;
		return;
	}
	else
	{
		m_piCtoPNum_Connect[cIndx]--;
		if(m_piCtoPNum_Connect[cIndx]==0){		m_piCtoPIndx_Connect[cIndx] = 0;	}
	}

	//A�@-1��������@�@�����@���ʂȃ������������Ă����C����̏ꍇ�ɂ���Ă͋��������������Ȃ�
	m_ppiCtoP_Connect[cIndx][indx] = -1;
}

/*!
 * �ύX�����@���q���N���X�^
 * @param[in] pIndx�@���q�ԍ�
 * @param[in] lIndx�@���q��������l�Ԗڂ̃N���X�^
 * @param[in] cIndx�@�N���X�^�ԍ�
 * @param[in] oIndx�@�N���X�^���ł̏����ԍ�
 */
void IceStructure::SetPtoC_Connect(int pIndx, int lIndx, int cIndx, int oIndx)
{//	cout << "SetParticleToCluster" << endl;
	//�G���[�`�F�b�N
	if(m_iPtoCMax/4 < lIndx)
	{
		cout << __FUNCTION__ << " Error::���q��������ڑ��N���X�^�����i�[�ł��Ȃ��Ȃ�܂����D������������܂���" << endl;
		cout << m_iPtoCMax/4 << "<" << lIndx << endl;
		return;
	}

	m_pppiPtoC_Connect[pIndx][lIndx][0] = cIndx;
	m_pppiPtoC_Connect[pIndx][lIndx][1] = oIndx;
}

/*!
 * �ύX�����@�N���X�^�����q
 */
void IceStructure::SetCtoP_Connect(int cIndx, vector<int> pIndxList)
{
	if(4 < pIndxList.size())
	{
		cout << __FUNCTION__ << " Error::�v�Z�N���X�^���܂ޗ��q�̏����i�[�ł��Ȃ��Ȃ�܂����D������������܂���" << endl;
		cout << 4 << "<" << pIndxList.size() << endl;
		return;
	}

	for(int i = 0; i < pIndxList.size(); i++)
	{
		m_ppiCtoP_Connect[cIndx][i] = pIndxList[i];
	}
}

void IceStructure::SetPtoCIndx_Connect(int pIndx, int pNum)
{
	//�G���[�`�F�b�N
	if(pNum >= m_iPtoCMax/4)
	{
		cout << __FUNCTION__ << " pIndx = " << pIndx << " pNum = " << pNum << endl;
		cout << "�z��ɋ󂫂��Ȃ��Ȃ�܂����D" << endl;

		m_piPtoCIndx_Calc[pIndx] = 0;
		return;
	}
	m_piPtoCIndx_Connect[pIndx] = pNum;
}

void IceStructure::SetCtoPIndx_Connect(int cIndx, int cNum)
{
	//�G���[�`�F�b�N
	if(cNum > 4)
	{
		cout << __FUNCTION__ << " cIndx = " << cIndx << " cNum = " << cNum << endl;
		cout << "�ڑ��N���X�^���ێ��ł��闱�q�������E�ł��D�z��ɋ󂫂��Ȃ��Ȃ�܂����D" << endl;
		return;
	}
	m_piCtoPIndx_Connect[cIndx] = cNum;
}

/*!
 * �擾�����@���q���N���X�^
 */
int* IceStructure::GetPtoC_Connect(int pIndx, int lIndx)
{
	return m_pppiPtoC_Connect[pIndx][lIndx];
}


/*!
 * �擾�����@�N���X�^�����q
 */
int IceStructure::GetCtoP_Connect(int cIndx, int lIndx)
{//	cout << __FUNCTION__ << endl;
	return m_ppiCtoP_Connect[cIndx][lIndx];
}

/*!
 * �X�V�����@���q���N���X�^
 */
void IceStructure::ResetOrderPToC(int cIndx, int oIndx)
{//	cout << __FUNCTION__ << endl;
	for(int i = 0; i < m_iPNum; i++)
	{
		for(int j = 0; j < m_piPtoCIndx_Connect[i]; j++)
		{
			int* coSet = GetPtoC_Connect(i, j);
			int jcIndx = coSet[0];
			int joIndx = coSet[1];

			if(jcIndx == -1 || joIndx == -1){	continue;	}
			if(jcIndx != cIndx){				continue;	}
			if(joIndx < oIndx){					continue;	}

			coSet[1]--;
		}
	}
//	cout << __FUNCTION__ << "End" << endl;
}

/*!
 * �����������@���q���N���X�^
 */
void IceStructure::ClearPtoC_Connect(int pIndx)
{	
	for(int i = 0; i < m_iPtoCMax/4; i++)
	{
		m_pppiPtoC_Connect[pIndx][i][0] = -1;
		m_pppiPtoC_Connect[pIndx][i][1] = -1;
	}

	m_piPtoCNum_Connect[pIndx] = 0;
	m_piPtoCIndx_Connect[pIndx] = 0;
}

/*!
 * �����������@�N���X�^�����q
 */
void IceStructure::ClearCtoP_Connect(int cIndx)
{
	for(int i = 0; i < 4; i++)
	{
		m_ppiCtoP_Connect[cIndx][i] = -1;
	}

	m_piCtoPIndx_Connect[cIndx] = 0;
	m_piCtoPNum_Connect[cIndx] = 0;
}

//------------------------------------�v�Z�����N���X�^-------------------------------------
/*!
 * �ǉ������@���q���N���X�^
 */
void IceStructure::AddPtoC_Calc()
{

}

void IceStructure::AddPtoC_Calc(int cIndx, int oIndx)
{

}

/*!
 * �ǉ������@�N���X�^�����q
 */
void IceStructure::AddCtoP_Calc(vector<int> pIndxes)
{

}

/*!
 * �폜�����@���q���N���X�^
 */
void IceStructure::DeletePtoC_Calc(int pIndx, int oIndx)
{
	if(m_pppiPtoC_Calc[pIndx][oIndx][0] == -1 || m_pppiPtoC_Calc[pIndx][oIndx][1] == -1)
	{
		return;
	}

	m_piPtoCNum_Calc[pIndx]--;

	if(m_piPtoCNum_Calc[pIndx] < 0)
	{
		cout << __FUNCTION__ << " Error::minus m_piPtoCNum_Calc[pIndx] = " << m_piPtoCNum_Calc[pIndx] << endl;
		m_piPtoCNum_Calc[pIndx] = 0;
		return;
	}

	m_pppiPtoC_Calc[pIndx][oIndx][0] = -1;
	m_pppiPtoC_Calc[pIndx][oIndx][1] = -1;
	m_pppiPtoC_Calc[pIndx][oIndx][2] = -1;
}

/*!
 * �폜�����@�N���X�^�����q
 */
void IceStructure::DeleteCtoP_Calc(int cIndx, int oIndx)
{
	//oIndx�͓Y�����ł͂Ȃ����Ȃ̂ŁC�T������oIndx�ڂ̂������������D
	int indx = -1;
	for(int i = 0; i < GetCtoPIndx_Calc(cIndx); i++)
	{
		if(m_ppiCtoP_Calc[cIndx][i][0] == -1 || m_ppiCtoP_Calc[cIndx][i][1] == -1){	continue;	}
		indx++;
		if(indx == oIndx){ indx = i;	break;	}				//oIndx�ڂ̏ꍇ�CIndx��ۑ�����break
	}

	if(m_ppiCtoP_Calc[cIndx][indx][0] == -1 || m_ppiCtoP_Calc[cIndx][indx][1] == -1){	return;	}

	m_piCtoPNum_Calc[cIndx]--;
	if(m_piCtoPNum_Calc[cIndx] < 0)
	{
		cout << __FUNCTION__ << " Error::minus m_piCtoPNum_Calc[cIndx] = " << m_piCtoPNum_Calc[cIndx] << endl;
		return;
	}

	//A�@-1��������@�@�����@���ʂȃ������������Ă����C����̏ꍇ�ɂ���Ă͋��������������Ȃ�
	//���̕ϐ��͎��ۂ͏z�͂��Ȃ�
	m_ppiCtoP_Calc[cIndx][indx][0] = -1;
	m_ppiCtoP_Calc[cIndx][indx][1] = -1;
}

/*!
 * �ύX�����@���q���N���X�^
 * @param[in] pIndx�@���q�ԍ�
 * @param[in] lIndx�@���q��������l�Ԗڂ̃N���X�^
 * @param[in] cIndx�@�N���X�^�ԍ�
 * @param[in] oIndx�@�N���X�^���ł̏����ԍ�
 * @param[in] layer  ���q�̑�����w�ԍ�
 */
void IceStructure::SetPtoC_Calc(int pIndx, int lIndx, int cIndx, int oIndx, int layer)
{//	cout << "SetParticleToCluster" << endl;
	if(m_iPtoCMax < lIndx)
	{	
		cout << __FUNCTION__ << " Error::���q��������v�Z�N���X�^�����i�[�ł��Ȃ��Ȃ�܂����D������������܂���" << endl;
		cout << m_iPtoCMax << "<" << lIndx << endl;
		cout << "pIndx = " << pIndx << " cIndx = " << cIndx << " oIndx = " << oIndx << " layer = " << layer << endl;
		return;
	}

	m_pppiPtoC_Calc[pIndx][lIndx][0] = cIndx;
	m_pppiPtoC_Calc[pIndx][lIndx][1] = oIndx;
	m_pppiPtoC_Calc[pIndx][lIndx][2] = layer;
}

/*!
 * �ύX�����@�N���X�^�����q
 */
void IceStructure::SetCtoP_Calc(int cIndx, vector<int> pIndxList, int* pLayerList)
{
	if(m_iCtoPMax < pIndxList.size())
	{
		cout << __FUNCTION__ << " Error::�v�Z�N���X�^���܂ޗ��q�̏����i�[�ł��Ȃ��Ȃ�܂����D������������܂���" << endl;
		cout << m_iCtoPMax << "<" << pIndxList.size() << endl;
		return;
	}

	for(int i = 0; i < pIndxList.size(); i++)
	{
		m_ppiCtoP_Calc[cIndx][i][0] = pIndxList[i];
		m_ppiCtoP_Calc[cIndx][i][1] = pLayerList[i];
	}
}

void IceStructure::SetPtoCIndx_Calc(int pIndx, int pNum)
{
	//�G���[�`�F�b�N
	if(pNum >= m_iPtoCMax)
	{
		cout << __FUNCTION__ << " pIndx = " << pIndx << " pNum = " << pNum << endl;
		cout << "�z��ɋ󂫂��Ȃ��Ȃ�܂����D" << endl;

		m_piPtoCIndx_Calc[pIndx] = 0;
		return;
	}
	m_piPtoCIndx_Calc[pIndx] = pNum;
}

void IceStructure::SetCtoPIndx_Calc(int cIndx, int cNum)
{
	//�G���[�`�F�b�N
	if(cNum >= m_iCtoPMax)
	{
		cout << __FUNCTION__ << " cIndx = " << cIndx << " cNum = " << cNum << endl;
		cout << "�v�Z�N���X�^���ێ��ł��闱�q�������E�ł��D�z��ɋ󂫂��Ȃ��Ȃ�܂����D" << endl;
		return;
	}
	m_piCtoPIndx_Calc[cIndx] = cNum;
}

/*!
 * �����������@���q���N���X�^
 */
void IceStructure::ClearPtoC_Calc(int pIndx)
{	
	for(int i = 0; i < m_iPtoCMax; i++)
	{
		m_pppiPtoC_Calc[pIndx][i][0] = -1;
		m_pppiPtoC_Calc[pIndx][i][1] = -1;
		m_pppiPtoC_Calc[pIndx][i][1] = -1;
	}

	m_piPtoCNum_Calc[pIndx] = 0;
	m_piPtoCIndx_Calc[pIndx] = 0;
}

/*!
 * �����������@�N���X�^�����q
 */
void IceStructure::ClearCtoP_Calc(int cIndx)
{
	for(int i = 0; i < m_piCtoPIndx_Calc[cIndx]; i++)
	{
		m_ppiCtoP_Calc[cIndx][i][0] = -1;
		m_ppiCtoP_Calc[cIndx][i][1] = -1;
	}

	m_piCtoPIndx_Calc[cIndx] = 0;
	m_piCtoPNum_Calc[cIndx] = 0;
}

/*!
 * �擾�����@���q���N���X�^
 */
int* IceStructure::GetPtoC_Calc(int pIndx, int lIndx)
{//	cout << __FUNCTION__ << endl;
//	cout << "m_iPtoCMax = " << m_iPtoCMax << " lIndx = " << lIndx << endl;

	if(m_iPtoCMax <= lIndx)
	{	
		cout << __FUNCTION__ << " Error::���q��������v�Z�N���X�^�����i�[�ł��Ȃ��Ȃ�܂����D������������܂���" << endl;
		cout << m_iPtoCMax << "<" << lIndx << " pIndx = " << pIndx << endl;
		return m_pppiPtoC_Calc[pIndx][0];
	}
	return m_pppiPtoC_Calc[pIndx][lIndx];
}

/*!
 * �擾�����@�N���X�^�����q
 */
int* IceStructure::GetCtoP_Calc(int cIndx, int lIndx)
{
	return m_ppiCtoP_Calc[cIndx][lIndx];
}

/*!
 * �X�V�����@���q���N���X�^
 */
void IceStructure::ResetOrderPToCCalc(int pIndx, int cIndx, int oIndx)
{

}

//-----------------------------------�ߖT�N���X�^-------------------------------------------
/*!
 * �ǉ�����
 */
void IceStructure::AddNeighborCluster(int cIndx, int layer)
{
}

/*!
 * �ύX�����@��������p�@�ڑ��N���X�^�̋ߖT�N���X�^��o�^
 */
void IceStructure::SetNeighborCluster(int cIndx, int layer)
{//	cout << __FUNCTION__ << endl;
	//�����N���X�^�Ɋ܂܂�闱�q�����̃N���X�^�ɂ��܂܂�Ă���ꍇ�C���̃N���X�^���ߖT�Ƃ��ēo�^
	int j = 0, k = 0;

	int ipIndx = 0, kpIndx = 0;
	int jcIndx = 0;

	int* jcIndxList;
	bool bFind = false;

	for(int i = 0; i < GetCtoPIndx_Connect(cIndx); i++)
	{
		ipIndx = GetCtoP_Connect(cIndx, i);
		if(ipIndx == -1){	continue;	}

		//�T��
		for(j = 0; j < GetPtoCIndx_Connect(ipIndx); j++)
		{
			jcIndxList = GetPtoC_Connect(ipIndx, j);
			if(jcIndxList[0] == -1 || jcIndxList[1] == -1){	continue;	}
			if(jcIndxList[0] == cIndx){	continue;	}

			//�����N���X�^�����Ɋ܂�ł��Ȃ����̃`�F�b�N
			bFind = false;
			for(k = 0; k < m_ppiNCNum[cIndx]; k++)
			{
				if(jcIndxList[0] == m_ppiNeighborCluster[cIndx][k][0]){ bFind = true; break;}
			}
			if( bFind == true ){	continue;	}

			//�ߖT�N���X�^�̓o�^�{�J�E���g
			if(m_ppiNCNum[cIndx] >= m_iNeighborMax){ cout << "Init �ߖT���q�̃�����������Ȃ��Ȃ�܂���" << endl;	}
			else
			{
				m_ppiNeighborCluster[cIndx][m_ppiNCNum[cIndx]][0] = jcIndxList[0];
				m_ppiNeighborCluster[cIndx][m_ppiNCNum[cIndx]][1] = 0;
				m_ppiNCNum[cIndx]++;
			}
		}
	}

//	cout << __FUNCTION__ << "Check1" << endl;
	//���񉻂̖��c
	int nowSize = 0, nowIndx = 0;
	int l = 0, m = 0;
	int* lcIndxList;
	
	//layer�w�����ǂ��ċߖT�N���X�^���擾����
	for(int i = 0; i < layer; i++)
	{
		nowSize = m_ppiNCNum[cIndx];

		//�T������N���X�^��nowSize��nowIndx�Ő������Ă���
		for(j = nowIndx; j < nowSize; j++)
		{
			jcIndx = m_ppiNeighborCluster[cIndx][j][0];				//�ߖT�N���X�^�̂ЂƂ�
			//�ߖT�N���X�^�Ɋ܂܂�闱�q�����̃N���X�^�ɂ��܂܂�Ă���ꍇ�C���̃N���X�^���ߖT�Ƃ��ēo�^
			for(k = 0; k < GetCtoPIndx_Connect(jcIndx); k++)
			{		
				kpIndx = GetCtoP_Connect(jcIndx, k);				//�ߖT�N���X�^�Ɋ܂܂�Ă��闱�q�̂ЂƂ�
				if(kpIndx == -1){	continue;	}

				for(l = 0; l < GetPtoCIndx_Connect(kpIndx); l++)
				{
					lcIndxList = GetPtoC_Connect(kpIndx, l);
					if(lcIndxList[0] == -1 || lcIndxList[1] == -1){	continue;	}
					if(lcIndxList[0] == cIndx){	continue;	}

					//�����N���X�^�����Ɋ܂�ł��Ȃ����̃`�F�b�N
					bFind = false;
					for(m = 0; m < m_ppiNCNum[cIndx]; m++)
					{
						if(lcIndxList[0] == m_ppiNeighborCluster[cIndx][m][0]){ bFind = true; break;}
					}
					if(bFind == true){	continue;	}

					//�ߖT�N���X�^�̓o�^�{�J�E���g
					if(m_ppiNCNum[cIndx] >= m_iNeighborMax){ cout << __FUNCTION__ << " �ߖT���q�̃�����������Ȃ��Ȃ�܂��� " << m_ppiNCNum[cIndx] << endl;	}
					else
					{
//						if(m_ppiNCNum[cIndx] > 200 ) return;	//�ߖT�N���X�^���̐���
						m_ppiNeighborCluster[cIndx][m_ppiNCNum[cIndx]][0] = lcIndxList[0];
						m_ppiNeighborCluster[cIndx][m_ppiNCNum[cIndx]][1] = i+1;
						m_ppiNCNum[cIndx]++;
					}
				}
			}
		}
		nowIndx = nowSize;											//���̃��[�v�J�n���̃X�^�[�g�ԍ����X�V
	}

//	cout << __FUNCTION__ << " m_ppiNCNum[cIndx] = " << m_ppiNCNum[cIndx] << endl;
//	cout << __FUNCTION__ << "Check2" << endl;
}


/*!
 * �ύX�����@���݂̋ߖT�N���X�^�̏�����ɁC�ڑ��N���X�^�̋ߖT�N���X�^��o�^
 * @param[in] cIndx�@�N���X�^�ԍ�
 * @param[in] layer�@�T���I�����C���[��
 * @param[in] jlayer �T���J�n���C���[��
 */

void IceStructure::SetNeighborClusterFromCluster(int cIndx, int layer, int jlayer)
{//	cout << __FUNCTION__ << endl;

	//�ߖT�N���X�^���̏�����
	//jlayer�ȏ���������C�T���J�n�ʒu�����߂�
	//��Cjlayer=2�Ȃ�C���߂�layer=1�ƂȂ�ꏊ��startNum�ƂȂ�@�܂��Clayer>=2�ƂȂ�ꏊ�������������
	int startNum = 0;

	for(int i = m_ppiNCNum[cIndx]-1; 0 <= i; i--)
	{
		if(m_ppiNeighborCluster[cIndx][i][1] < jlayer-1)
		{
			startNum = i+1;
			break;
		}
		else if(m_ppiNeighborCluster[cIndx][i][1] >= jlayer)
		{
			m_ppiNeighborCluster[cIndx][i][0] = -1;
			m_ppiNeighborCluster[cIndx][i][1] = -1;
		}
	}

	m_ppiNCNum[cIndx] = startNum;

	//�ߖT�N���X�^�̎擾
	//jlayer�w����J�n�@jlayer-1�w�܂ł̃N���X�^�͎擾�ς�
	int j = 0, k = 0, l = 0, m = 0;
	int ipIndx = 0, kpIndx = 0;
	int jcIndx = 0;
	int nowSize = 0, nowIndx = 0;

	int* lcIndxList;

	bool bFind = false;

	//layer�w�����ǂ��ċߖT�N���X�^���擾����
	for(int i = 0; i < layer; i++)
	{
		nowSize = startNum;

		//�T������N���X�^��nowSize��nowIndx�Ő������Ă���
		for(j = nowIndx; j < nowSize; j++)
		{
			jcIndx = m_ppiNeighborCluster[cIndx][j][0];				//�ߖT�N���X�^�̂ЂƂ�
			//�ߖT�N���X�^�Ɋ܂܂�闱�q�����̃N���X�^�ɂ��܂܂�Ă���ꍇ�C���̃N���X�^���ߖT�Ƃ��ēo�^
			for(k = 0; k < GetCtoPIndx_Connect(jcIndx); k++)
			{		
				kpIndx = GetCtoP_Connect(jcIndx, k);				//�ߖT�N���X�^�Ɋ܂܂�Ă��闱�q�̂ЂƂ�
				if(kpIndx == -1){	continue;	}

				for(l = 0; l < GetPtoCIndx_Connect(kpIndx); l++)
				{
					lcIndxList = GetPtoC_Connect(kpIndx, l);
					if(lcIndxList[0] == -1 || lcIndxList[1] == -1){	continue;	}
					if(lcIndxList[0] == cIndx){	continue;	}

					//�����N���X�^�����Ɋ܂�ł��Ȃ����̃`�F�b�N
					bFind = false;
					for(m = 0; m < m_ppiNCNum[cIndx]; m++)
					{
						if(lcIndxList[0] == m_ppiNeighborCluster[cIndx][m][0]){ bFind = true; break;}
					}
					if(bFind == true){	continue;	}

					//�ߖT�N���X�^�̓o�^�{�J�E���g
					if(m_ppiNCNum[cIndx] >= m_iNeighborMax){ cout << __FUNCTION__ << " �ߖT���q�̃�����������Ȃ��Ȃ�܂��� " << m_ppiNCNum[cIndx] << endl;	}
					else
					{
						m_ppiNeighborCluster[cIndx][m_ppiNCNum[cIndx]][0] = lcIndxList[0];
						m_ppiNeighborCluster[cIndx][m_ppiNCNum[cIndx]][1] = i+1;
//						if(m_ppiNCNum[cIndx] > 200 ) return;	//�ߖT�N���X�^���̐���
						m_ppiNCNum[cIndx]++;
					}
				}
			}
		}
		nowIndx = nowSize;											//���̃��[�v�J�n���̃X�^�[�g�ԍ����X�V
	}

//	cout << __FUNCTION__ << " m_ppiNCNum[cIndx] = " << m_ppiNCNum[cIndx] << endl;
}

/*!
 * ����������
 */
void IceStructure::ClearNeighborCluster(int cIndx)
{//	cout << __FUNCTION__ << "check1" << endl;
	for(int i = 0; i < m_ppiNCNum[cIndx]; i++)
	{
		m_ppiNeighborCluster[cIndx][i][0] = -1;
		m_ppiNeighborCluster[cIndx][i][1] = -1;
	}

	m_ppiNCNum[cIndx] = 0;
//	cout << __FUNCTION__ << "check1" << endl;
}

/*!
 * �擾����
 */
int* IceStructure::GetNeighborCluster(int cIndx, int lIndx)
{//	cout << __FUNCTION__ << endl;
	return m_ppiNeighborCluster[cIndx][lIndx];
}

//--------------------------------�ڑ����ƌv�Z�����̊֘A-------------------------------
/*!
 * ����������
 */
void IceStructure::MakeCalcToConnect()
{
	m_pppiCalcToConnect = new int**[m_iCNumMax];
	int j = 0;

	for(int i = 0; i < m_iCNumMax; i++)
	{
//		cout << "i = " << i << " maxCnum = " << m_iCNumMax << endl;
		m_pppiCalcToConnect[i] = new int*[m_iCtoPMax];

		for(j = 0; j < m_iCtoPMax; j++)
		{
			m_pppiCalcToConnect[i][j] = new int[2];
			m_pppiCalcToConnect[i][j][0] = -1;
			m_pppiCalcToConnect[i][j][1] = -1;
		}
	}
}

/*!
 * �ǉ�����
 */
void IceStructure::AddCalcToConnect(int caIndx, int coIndx, int oIndx)
{

}

/*!
 * �폜����
 */
void IceStructure::DeleteCalcToConnect(int caIndx, int poIndx)
{

}

/*!
 * �o�^����
 * @param[in] caIndx�@�v�Z�N���X�^
 * @param[in] ocaIndx�@�v�Z�N���X�^�ɏ������闱�q�̏����ԍ�
 * @param[in] coIndx�@���q���������Ă���ڑ��ԍ�
 * @param[in] oIndx�@�@���q���ڑ��N���X�^�ɏ������Ă��鏇���ԍ�
 */
void IceStructure::SetCalcToConnect(int caIndx, int ocaIndx, int coIndx, int oIndx)
{//	cout << __FUNCTION__ << endl;
	if(m_iCtoPMax <= ocaIndx)
	{
		cout << __FUNCTION__ << " ������������܂���" << endl;
		cout << m_iCtoPMax << "<" << ocaIndx << endl;
		return;
	}

	m_pppiCalcToConnect[caIndx][ocaIndx][0] = coIndx;
	m_pppiCalcToConnect[caIndx][ocaIndx][1] = oIndx;
}

/*!
 * �擾����
 */
int* IceStructure::GetCalcToConnect(int caIndx, int ocaIndx)
{//	cout << __FUNCTION__ << " caIndx = " << caIndx << " ocaIndx = " << ocaIndx;
	if(m_iCtoPMax <= ocaIndx)
	{
		cout << __FUNCTION__ << " �s���A�N�Z�X�ł�" << endl;
		cout << m_iCtoPMax << "<" << ocaIndx << endl;
		return m_pppiCalcToConnect[0][0];
	}
//	cout << __FUNCTION__ << "End" << endl;
	return m_pppiCalcToConnect[caIndx][ocaIndx];
}

/*!
 * �X�V����
 */
void IceStructure::ResetOrderCaToCo(int caIndx, int poIndx)
{

}

/*!
 * ����������
 */
void IceStructure::ClearCalcToConnect(int caIndx)
{//	cout << __FUNCTION__ << "check1" << endl;
	for(int i = 0; i < m_iCtoPMax; i++)
	{
		m_pppiCalcToConnect[caIndx][i][0] = -1;
		m_pppiCalcToConnect[caIndx][i][1] = -1;
	}
//	cout << __FUNCTION__ << "check2" << endl;
}

//----------------------------------�f�o�b�O�@�\---------------------------------------------
void IceStructure::DebugPtoC_Connect(int pIndx)
{	cout << "DebugPtoC_Connect pIndx=" << pIndx;
	cout << " num=" << GetPtoCNum_Connect(pIndx) << " Indx=" << GetPtoCIndx_Connect(pIndx);
	
	for(int i = 0; i < GetPtoCIndx_Connect(pIndx); i++)
	{
			int* coSet = GetPtoC_Connect(pIndx, i);
			cout <<" c=" << coSet[0] << " o=" << coSet[1];
	}
	cout << endl;
}

void IceStructure::DebugCtoP_Connect(int cIndx)
{	cout << "DebugCtoP_Connect cIndx=" << cIndx;
	cout << " num=" << GetCtoPNum_Connect(cIndx) << " Indx=" << GetCtoPIndx_Connect(cIndx);

	for(int i = 0; i < GetCtoPIndx_Connect(cIndx); i++)
	{
		cout << " " << GetCtoP_Connect(cIndx, i);
	}
	cout << endl;
}

void IceStructure::DebugPtoC_Calc(int pIndx)
{	cout << "DebugPtoC_Calc	 pIndx=" << pIndx;
	cout << " num=" << GetPtoCNum_Calc(pIndx) << " Indx=" << GetPtoCIndx_Calc(pIndx);

	for(int i = 0; i < GetPtoCIndx_Calc(pIndx); i++)
	{
		int* coSet = GetPtoC_Calc(pIndx, i);
		cout << " c=" << coSet[0] << " o=" << coSet[1] << " L=" << coSet[2];
	}
	cout << endl;
}

void IceStructure::DebugCtoP_Calc(int cIndx)
{	cout << "DebugCtoP_Calc	cIndx=" << cIndx;
	cout << " num=" << GetCtoPNum_Calc(cIndx) << " Indx=" << GetCtoPIndx_Calc(cIndx);

	for(int i = 0; i < GetCtoPIndx_Calc(cIndx); i++)
	{
		int* coSet = GetCtoP_Calc(cIndx, i);
		cout << " p=" << coSet[0] << " L=" << coSet[1];
	}
	cout << endl;
}

void IceStructure::DebugCalcToConnect(int cIndx)
{	cout << "DebugCalcToConnect	 cIndx=" << cIndx;

	for(unsigned j = 0; j < GetCtoPIndx_Calc(cIndx); j++ )
	{
		cout << " co=" << m_pppiCalcToConnect[cIndx][j][0] << " o=" << m_pppiCalcToConnect[cIndx][j][1];
	}
	cout << endl;
}

void IceStructure::DebugNeighborCluster(int cIndx)
{	cout << "DebugNeighborCluster cIndx=" << cIndx;
	cout << " num = " << GetNCNum(cIndx);
	for(unsigned j = 0; j < GetNCNum(cIndx); j++)
	{
		cout << " NC=" << GetNeighborCluster(cIndx, j)[0] << " Ly=" << GetNeighborCluster(cIndx, j)[1];
	}
	cout << endl;
}

//--------------------------------------------���q�x�[�X-------------------------------------
void IceStructure::DebugPtoT(int pIndx)
{	cout << __FUNCTION__ << " pIndx = " << pIndx;
	cout << " num=" << GetPtoTNum(pIndx) << " Indx=" << GetPtoTIndx(pIndx);
	
	for(int i = 0; i < GetPtoTIndx(pIndx); i++)
	{
			int* coSet = GetPtoT(pIndx, i);
			cout <<" c=" << coSet[0] << " o=" << coSet[1];
	}
	cout << endl;
}

void IceStructure::DebugPtoC(int pIndx)
{	cout << __FUNCTION__ << " pIndx = " << pIndx;
	cout << " num=" << GetPtoCNum(pIndx) << " Indx=" << GetPtoCIndx(pIndx);
	
	for(int i = 0; i < GetPtoCIndx(pIndx); i++)
	{
			int* coSet = GetPtoC(pIndx, i);
			cout <<" c=" << coSet[0] << " o=" << coSet[1];
	}
	cout << endl;

}

void IceStructure::DebugCtoP(int cIndx)
{	cout << __FUNCTION__ << " cIndx = " << cIndx;
	cout << " num=" << GetCtoPNum(cIndx) << " Indx=" << GetCtoPIndx(cIndx);

	for(int i = 0; i < GetCtoPIndx(cIndx); i++)
	{
		int* plSet = GetCtoP(cIndx, i);
		cout <<" p=" << plSet[0] << " l=" << plSet[1];
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
		cout << " NC=" << GetNeighborTetra(tIndx, j)[0] << " Ly=" << GetNeighborTetra(tIndx, j)[1];
	}
	cout << endl;
}