/*!
  @file IceStructure.h
	
  @brief �_�C�ӁC�ʂ�p���ĕX�\���𑀍삷��N���X
 
  @author Ryou Nakasone
  @date 2013-07
*/
//---------------------------------------------------------------------------
//���̃t�@�C���̃N���X�E�֐������ȏ�R���p�C�������̂�h�����߂̏����i�C���N���[�h�K�[�h�j
#ifndef IceStructure_H
#define IceStructure_H
//---------------------------------------------------------------------------

#include <iostream>
#include <vector>
#include <algorithm>

//#include "rx_matrix.h"
#include <UtilityScript\mk_Vector2D.h>
#include <UtilityScript\mk_Vector3D.h>

using namespace std;

class IceStructure
{
public:
	IceStructure(int pNum, int cNum, int tNum);
	~IceStructure(void);

	void SetParticleNum(int pNum){	m_iPNum = pNum; }		//���݂̗��q��
	void SetClusterNum(int cNum){	m_iCNum = cNum; }		//���݂̃N���X�^��
	void SetTetraNum(int tNum){		m_iTNum = tNum;	}		//���݂̎l�ʑ̐�

	int GetParticleNum(void){		return m_iPNum;	}
	int GetClusterNum(void){		return m_iCNum;	}
	int GetTetraNum(void){			return m_iTNum;	}

	void InitTetraInfo();											//�l�ʑ̏��̃������m��
	void InitClusterInfo();											//�N���X�^���̃������m��

	int  GetPtoCMax(void){	return m_iPtoCMax;	}

	void CountPtoC(int pIndx){	m_piPtoCNum[pIndx]++;	}
	void CountPtoT(int pIndx){	m_piPtoTNum[pIndx]++;	}
	void CountCtoP(int cIndx){	m_piCtoPNum[cIndx]++;	}
	void CountCtoT(int cIndx){	}
	void CountTtoP(int tIndx){	m_piTtoPNum[tIndx]++;	}
	void CountTtoC(int tIndx){	}

	void CountNT(int tIndx){	m_piNTNum[tIndx]++;		}

	void SetPtoCIndx(int pIndx, int indx){	m_piPtoCIndx[pIndx] = indx;	} 
	void SetCtoPIndx(int cIndx, int indx){	m_piCtoPIndx[cIndx] = indx;	}
	void SetPtoTIndx(int pIndx, int indx){	m_piPtoTIndx[pIndx] = indx; }
	void SetTtoPIndx(int tIndx, int indx){	m_piTtoPIndx[tIndx] = indx; }

	int	 GetPtoTNum(int pIndx){		return m_piPtoTNum[pIndx];	}
	int	 GetPtoCNum(int pIndx){		return m_piPtoCNum[pIndx];	}
	int  GetTtoPNum(int tIndx){		return m_piTtoPNum[tIndx];	}
	int  GetCtoPNum(int cIndx){		return m_piCtoPNum[cIndx];	}

	int  GetNTNum(int tIndx){		return m_piNTNum[tIndx];	}

	int  GetPtoTIndx(int pIndx){	return m_piPtoTIndx[pIndx];	}
	int  GetPtoCIndx(int pIndx){	return m_piPtoCIndx[pIndx];	}
	int  GetTtoPIndx(int tIndx){	return m_piTtoPIndx[tIndx];	}
	int  GetCtoPIndx(int cIndx){	return m_piCtoPIndx[cIndx];	}

	int  GetPtoTFreeIndx(int pIndx);
	int  GetPtoCFreeIndx(int pIndx);

	void SetPtoT(int pIndx, int lIndx, int tIndx, int oIndx);				//���q��������l�ʑ̂̓o�^�@�@���q�ԍ��C���q�������C�l�ʑ̔ԍ��C�l�ʑ̓�����
	void SetPtoC(int pIndx, int lIndx, int cIndx, int oIndx, int layer);	//���q��������N���X�^�̓o�^�@���q�ԍ��C���q�������C�N���X�^�ԍ��C�N���X�^������
	
	void SetCtoP(int cIndx, const vector<int>& pIndxList, int* pLayerList);	//�N���X�^�ɑ����闱�q�̓o�^�@�N���X�^�ԍ��C���q���X�g
	void SetTtoP(int tIndx, vector<int>& pIndxList);						//�l�ʑ̂ɑ����闱�q�̓o�^�@�@�l�ʑ̔ԍ��C���q���X�g
	
	void SetNeighborTetra(int tIndx, int layer);
	void SetNeighborTetraFromLayer(int tIndx, int searchLayer, int deleteLayer);

	void DeleteTtoP(int tIndx, int lIndx);
	void DeletePtoC(int pIndx, int lIndx);
	void DeletePtoT(int pIndx, int lIndx);


	int  GetTtoP(int tIndx, int lIndx);

	//�C��
	int GetPtoC(int pIndx, int lIndx, int oIndx);
	int GetPtoT(int pIndx, int lIndx, int oIndx);
	int GetCtoP(int cIndx, int lIndx, int oIndx);

	int* GetNeighborTetra(int tIndx, int lIndx);

	void ClearPtoT(int pIndx);
	void ClearPtoC(int pIndx);
	void ClearCtoP(int cIndx);
	void ClearTtoP(int tIndx);

	void ClearNeighborTetra(int tIndx);
	void ClearNeighborTetraFromLayer(int tIndx, int layer);

	int  CheckNeighborTetra(int tIndx, int checkTIndx);

	//�f�o�b�O
	void DebugPtoT(int pIndx);
	void DebugPtoC(int pIndx);
	void DebugCtoP(int cIndx);
	void DebugTtoP(int tIndx);
	void DebugNeighborTetra(int tIndx);

	//�t���O�@���g�p
	void SetPFlag(int indx, bool state){	m_pbPFlag[indx] = state;	}
	void SetCFlag(int indx, bool state){	m_pbCFlag[indx] = state;	}
	void SetTFlag(int indx, bool state){	m_pbTFlag[indx] = state;	}

	bool GetPFlag(int indx){	return m_pbPFlag[indx];	}
	bool GetCFlag(int indx){	return m_pbCFlag[indx];	}
	bool GetTFlag(int indx){	return m_pbTFlag[indx];	}

	void ResetPFlag(int endNum);
	void ResetCFlag(int endNum);
	void ResetTFlag(int endNum);

protected:

	int m_iPNumMax;								//�ő嗱�q��
	int m_iPNum;								//���݂̗��q��
	
	int m_iPtoCMax;								//���q���N���X�^�ɏ�������ő吔�@connect��calc�̔����ł���
	int m_iCtoPMax;								//�N���X�^���܂ޗ��q�̍ő吔�@�@�@connect�͍ő�S�ŌŒ�	

	int m_iTNumMax;
	int m_iTNum;

	int m_iPtoTMax;
	int m_iTtoPMax;

	int m_iCNumMax;								//�ő�N���X�^��
	int m_iCNum;								//���݂̃N���X�^��

	int m_iNeighborMax;							//�ߖT���q�̍ő吔

	//�Ƃ肠�����C�|�C���^�͎g��Ȃ�
	//���q��
	mk_Vector3D<int> m_mk3DiPtoC;				//���q���N���X�^
	mk_Vector3D<int> m_mk3DiPtoT;				//���q���l�ʑ�

	int*   m_piPtoCNum;							//���q���N���X�^�̌�
	int*   m_piPtoTNum;							//���q���l�ʑ̂̌�

	int*   m_piPtoCIndx;
	int*   m_piPtoTIndx;

	//�N���X�^��
	mk_Vector3D<int> m_mk3DiCtoP;				//�N���X�^�����q�@�N���X�^�ɏ������闱�q��Ԃ��@���q�̐ڑ����
	mk_Vector3D<int> m_mk3DiCtoT;				//�N���X�^���l�ʑ�				

	int*   m_piCtoPNum;							//�N���X�^�����q�̌�
	int*   m_piCtoTNum;							//�N���X�^���l�ʑ̂̌�

	int*   m_piCtoPIndx;
	int*   m_piCtoTIndx;

	//�l�ʑ́�
	int**  m_ppiTtoP;							//�l�ʑ́����q
	int*** m_pppiTtoC;							//�l�ʑ́��N���X�^

	int*   m_piTtoPNum;							//�l�ʑ́����q�̌�
	int*   m_piTtoCNum;							//�l�ʑ́��N���X�^�̌�

	int*   m_piTtoPIndx;
	int*   m_piTtoCIndx;

	//�ߖT�l�ʑ�
	int*** m_pppiNeighborTetra;					//�ߖT�l�ʑ�
	int*   m_piNTNum;							//�e�ߖT�l�ʑ̂̌�

	//�T���p�t���O�@���g�p
	bool* m_pbPFlag;							//���q
	bool* m_pbCFlag;							//�N���X�^
	bool* m_pbTFlag;							//�l�ʑ�

	vector<vector<int>>	m_vviEdgeToCluster;					//�Ӂ��N���X�^	����ӂ��ǂ̃N���X�^�ɏ������Ă��邩�@�K���P�ȏ�
	vector<vector<int>>	m_vviEdgeToParticle;				//�Ӂ����q�@�@�@����ӂ��ǂ̗��q��ڑ����Ă��邩�@�@�@�K���Q��

	vector<vector<int>>	m_vviClusterToEdge;					//�N���X�^���Ӂ@����N���X�^�ɂǂ̕ӂ��������Ă��邩�@�K���P���R��
	vector<vector<int>>	m_vviParticleToEdge;				//���q���Ӂ@�@�@���闱�q�͂ǂ̕ӂŐڑ�����Ă��邩

};
#endif