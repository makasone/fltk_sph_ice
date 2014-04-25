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

#include <stdio.h>
#include <iostream>
#include <vector>
#include <map>
#include <algorithm>
#include <omp.h>
#include "rx_matrix.h"


using namespace std;

class IceStructure
{
public:
	IceStructure(int pNumMax, int cNumMax, int tNumMax);				//�l�ʑ̃x�[�X
	IceStructure(int pNum, int cNum);
	~IceStructure(void);

	void SetParticleNum(int pNum){	m_iPNum = pNum; }		//���݂̗��q��
	void SetClusterNum(int cNum){	m_iCNum = cNum; }		//���݂̃N���X�^��
	void SetTetraNum(int tNum){		m_iTNum = tNum;	}		//���݂̎l�ʑ̐�

	int GetParticleNum(void){		return m_iPNum;	}
	int GetClusterNum(void){		return m_iCNum;	}
	int GetTetraNum(void){			return m_iTNum;	}

//-------------------------------------���q�x�[�X����----------------------------------------
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

	int* GetPtoT(int pIndx, int lIndx);
	int* GetPtoC(int pIndx, int lIndx);
	int* GetCtoP(int cIndx, int lIndx);
	int  GetTtoP(int tIndx, int lIndx);

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

//-------------------------------------���q�x�[�X����----------------------------------------


//-------------------------------------�l�ʑ̃x�[�X����-------------------------------------
	void InitStateConnect();								//�ڑ����̃T�C�Y�m��
	void InitStateCalc();									//�v�Z���̃T�C�Y�m��

	//�ڑ����p�A�N�Z�b�T
	void AddPtoC_Connect();
	void AddPtoC_Connect(int cIndx, int oIndx);
	void AddCtoP_Connect(vector<int> pIndxes);

	void DeletePtoC_Connect(int pIndx, int coIndx);
	void DeleteCtoP_Connect(int cIndx, int oIndx);

	void SetPtoC_Connect(int pIndx, int lIndx, int cIndx, int oIndx);			//���q�ԍ��C�N���X�^�ԍ��C�����ԍ��@�ڑ����̍쐬
	void SetCtoP_Connect(int cIndx, vector<int> pIndxList);						//����������N���X�^�ԍ��C���q�ԍ��̃��X�g

	void SetPtoCNum_Connect(int pIndx, int pNum){	m_piPtoCNum_Connect[pIndx] = pNum;	}
	void SetCtoPNum_Connect(int cIndx, int cNum){	m_piCtoPNum_Connect[cIndx] = cNum;	}

	void SetPtoCIndx_Connect(int pIndx, int pNum);
	void SetCtoPIndx_Connect(int cIndx, int cNum);

	void CountPtoC_Connect(int pIndx){ m_piPtoCNum_Connect[pIndx]++; }			//���q���������Ă���ڑ��N���X�^�̃J�E���^
	void CountCtoP_Connect(int cIndx){ m_piCtoPNum_Connect[cIndx]++; }			//�ڑ��N���X�^�ɏ������Ă��闱�q�̃J�E���^

	int* GetPtoC_Connect(int pIndx, int lIndx);
	int GetCtoP_Connect(int cIndx, int lIndx);

	int GetPtoCNum_Connect(int pIndx){ return m_piPtoCNum_Connect[pIndx]; }		//SM�@�ƘA�g���邽�߂ɂ��d�v�ɂȂ�
	int GetCtoPNum_Connect(int cIndx){ return m_piCtoPNum_Connect[cIndx]; }		//SM�@�ƘA�g���邽�߂ɂ��d�v�ɂȂ�

	int GetPtoCIndx_Connect(int pIndx){ return m_piPtoCIndx_Connect[pIndx];	}
	int GetCtoPIndx_Connect(int cIndx){ return m_piCtoPIndx_Connect[cIndx];	}

	void ClearPtoC_Connect(int pIndx);
	void ClearCtoP_Connect(int cIndx);

	void ResetOrderPToC(int cIndx, int oIndx);

	//�v�Z�����N���X�^�p�A�N�Z�b�T
	void AddPtoC_Calc();
	void AddPtoC_Calc(int cIndx, int oIndx);
	void AddCtoP_Calc(vector<int> pIndxes);

	void DeletePtoC_Calc(int pIndx, int coIndx);
	void DeleteCtoP_Calc(int cIndx, int oIndx);

	void SetPtoC_Calc(int pIndx, int lIndx, int cIndx, int oIndx, int layer);	//���q�ԍ��C�N���X�^�ԍ��C�����ԍ��@�v�Z�����N���X�^�̏��
	void SetCtoP_Calc(int cIndx, vector<int> pIndxList, int* pLayerList);		//����������N���X�^�ԍ��C���q�ԍ��̃��X�g

	void SetPtoCNum_Calc(int pIndx, int pNum){	m_piPtoCNum_Calc[pIndx] = pNum;	}
	void SetCtoPNum_Calc(int cIndx, int cNum){	m_piCtoPNum_Calc[cIndx] = cNum;	}

	void SetPtoCIndx_Calc(int pIndx, int pNum);
	void SetCtoPIndx_Calc(int cIndx, int cNum);

	void CountPtoC_Calc(int pIndx){ m_piPtoCNum_Calc[pIndx]++; }				//���q���������Ă���v�Z�N���X�^�̃J�E���^
	void CountCtoP_Calc(int cIndx){ m_piCtoPNum_Calc[cIndx]++; }				//�v�Z�N���X�^�ɏ������Ă��闱�q�̃J�E���^

	int* GetPtoC_Calc(int pIndx, int lIndx);
	int* GetCtoP_Calc(int cIndx, int lIndx);

	int GetPtoCNum_Calc(int pIndx){	return m_piPtoCNum_Calc[pIndx]; }
	int GetCtoPNum_Calc(int cIndx){ return m_piCtoPNum_Calc[cIndx]; }

	int GetPtoCIndx_Calc(int pIndx){ return m_piPtoCIndx_Calc[pIndx];	}
	int GetCtoPIndx_Calc(int cIndx){ return m_piCtoPIndx_Calc[cIndx];	}

	void ClearPtoC_Calc(int pIndx);
	void ClearCtoP_Calc(int cIndx);

	void ResetOrderPToCCalc(int pIndx, int cIndx, int oIndx);

	//�ߖT�N���X�^�p�A�N�Z�b�T
	void AddNeighborCluster(int cIndx, int layer);							//�ߖT�N���X�^�̒ǉ��@����ߖT�N���X�^���X�g�̃f�[�^��ǉ�
	void SetNeighborCluster(int cIndx, int layer);							//�ߖT�N���X�^�̏C���@����ߖT�N���X�^���X�g�̃f�[�^����������
	void SetNeighborClusterFromCluster(int cIndx, int layer, int jlayer);	//�ߖT�N���X�^�̏C���@�ߖT�N���X�^�̎����𗘗p������������
	int* GetNeighborCluster(int cIndx, int lIndx);							//�ߖT�N���X�^�ԍ��@�N���X�^�̋ߖT�N���X�^�̃��X�g��Ԃ�
	int GetNCNum(int cIndx){ return m_ppiNCNum[cIndx]; }					//�ߖT�N���X�^�̐�
	void ClearNeighborCluster(int cIndx);

	////�v�Z�����N���X�^���ڑ����N���X�^
	void MakeCalcToConnect(void);
	void AddCalcToConnect(int caIndx, int coIndx, int oIndx);
	void DeleteCalcToConnect(int caIndx, int poIndx);
	int* GetCalcToConnect(int cIndx, int lIndx);
	void SetCalcToConnect(int caIndx, int ocaIndx, int coIndx, int oIndx);
	void ResetOrderCaToCo(int cIndx, int poIndx);
	void ClearCalcToConnect(int caIndx);

	//�f�o�b�O
	void DebugPtoC_Connect(int pIndx);
	void DebugCtoP_Connect(int cIndx);

	void DebugPtoC_Calc(int pIndx);
	void DebugCtoP_Calc(int cIndx);

	void DebugCalcToConnect(int cIndx);

	void DebugNeighborCluster(int cIndx);
//-------------------------------------�l�ʑ̃x�[�X����-------------------------------------

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

	int m_iNeighborMax;							//�ߖT���q�̍ő吔
	
	int m_iPtoCMax;								//���q���N���X�^�ɏ�������ő吔�@connect��calc�̔����ł���
	int m_iCtoPMax;								//�N���X�^���܂ޗ��q�̍ő吔�@�@�@connect�͍ő�S�ŌŒ�	
//-------------------------------------���q�x�[�X����----------------------------------------
	int m_iTNumMax;
	int m_iTNum;

	int m_iPtoTMax;
	int m_iTtoPMax;

	//���q��
	int*** m_pppiPtoC;							//���q���N���X�^�@���Ԗڂ̃N���X�^���ŉ��ԖڂȂ̂��𔻒�@0�Ԃ���n�܂�̂ɒ���
	int*** m_pppiPtoT;							//���q���l�ʑ�

	int*   m_piPtoCNum;							//���q���N���X�^�̌�
	int*   m_piPtoTNum;							//���q���l�ʑ̂̌�

	int*   m_piPtoCIndx;
	int*   m_piPtoTIndx;

	//�N���X�^��
	int***  m_pppiCtoP;							//�N���X�^�����q�@�N���X�^�ɏ������闱�q��Ԃ��@���q�̐ڑ����
	int*** m_pppiCtoT;							//�N���X�^���l�ʑ�

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
//-------------------------------------���q�x�[�X����----------------------------------------

//-------------------------------------�l�ʑ̃x�[�X����-------------------------------------
	int m_iCNumMax;								//�ő�N���X�^��
	int m_iCNum;								//���݂̃N���X�^��

	//�ڑ����
	int*** m_pppiPtoC_Connect;					//���q���N���X�^�@���Ԗڂ̃N���X�^���ŉ��ԖڂȂ̂��𔻒�@0�Ԃ���n�܂�̂ɒ���
	int**  m_ppiCtoP_Connect;					//�N���X�^�����q�@�N���X�^�ɏ������闱�q��Ԃ��@���q�̐ڑ����

	int* m_piPtoCNum_Connect;					//���q���ڑ��N���X�^�ɏ������Ă����
	int* m_piCtoPNum_Connect;					//�N���X�^�����q���܂�ł����

	int* m_piPtoCIndx_Connect;					//���q���ڑ��N���X�^�ɏ������Ă��邱�Ƃ�ۑ�����z��́C���݂̓Y�����ԍ�
	int* m_piCtoPIndx_Connect;					//�ڑ��N���X�^�����q���܂�ł��邱�Ƃ�ۑ�����z��́C���݂̓Y�����ԍ�

	//�v�Z�����N���X�^
	int*** m_pppiPtoC_Calc;						//���q���N���X�^�@���Ԗڂ̃N���X�^���ŉ��ԖڂȂ̂��𔻒�@0�Ԃ���n�܂�̂ɒ���
	int*** m_ppiCtoP_Calc;						//�N���X�^�����q�@�N���X�^�ɏ������闱�q��Ԃ��@���q�̐ڑ����

	int* m_piPtoCNum_Calc;						//���q���v�Z�N���X�^�ɏ������Ă����
	int* m_piCtoPNum_Calc;						//�N���X�^�����q���܂�ł����

	int* m_piPtoCIndx_Calc;						//���q���v�Z�N���X�^�ɏ������Ă��邱�Ƃ�ۑ�����z��́C���݂̓Y�����ԍ�
	int* m_piCtoPIndx_Calc;						//�v�Z�N���X�^�����q���܂�ł��邱�Ƃ�ۑ�����z��́C���݂̓Y�����ԍ�

	//�ߖT�N���X�^
	int*** m_ppiNeighborCluster;				//�ߖT�N���X�^�̑g�@���ꂼ��ɗ��q�̃��X�g��p��
	int*  m_ppiNCNum;							//�ڑ��N���X�^�̋ߖT�ƂȂ�N���X�^�̌�

	//�v�Z�����N���X�^���ڑ����N���X�^
	int*** m_pppiCalcToConnect;					//�v�Z�����N���X�^���ڑ����N���X�^
//-------------------------------------�l�ʑ̃x�[�X����-------------------------------------

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