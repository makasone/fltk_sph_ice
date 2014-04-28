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
	m_mk3DiPtoC.SetSize(m_iPNumMax, m_iPtoCMax, 3);

	for(int i = 0; i < m_iPNumMax; i++)
	{
		for(int j = 0; j < m_iPtoCMax; j++)
		{
			//[0] = �N���X�^�ԍ�
			//[1] = �N���X�^���ł̔ԍ�
			//[2] = layer�ԍ�
			for(int k = 0; k < 3; k++)
			{
				m_mk3DiPtoC(i, j, k) = -1;
			}
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
				m_pppiNeighborTetra[tIndx][GetNTNum(tIndx)][0] = GetPtoT(ipIndx, j, 0);
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
					if(GetPtoT(kpIndx, j, 0) == -1 
					|| GetPtoT(kpIndx, j, 1) == -1
					|| GetPtoT(kpIndx, j, 0) == tIndx)
					{
						continue;
					}

					//�����l�ʑ̂����Ɋ܂�ł��Ȃ����̃`�F�b�N
					if(CheckNeighborTetra( tIndx, GetPtoT(kpIndx, j, 0) ) != -1)
					{
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
						m_pppiNeighborTetra[tIndx][GetNTNum(tIndx)][0] = GetPtoT(kpIndx, j, 0);
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
				m_pppiNeighborTetra[tIndx][GetNTNum(tIndx)][0] = GetPtoT(ipIndx, j, 0);
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
						m_pppiNeighborTetra[tIndx][GetNTNum(tIndx)][0] = GetPtoT(kpIndx, l, 0);
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
	return m_ppiTtoP[tIndx][lIndx];
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

//--------------------------------------------���q�x�[�X-------------------------------------
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