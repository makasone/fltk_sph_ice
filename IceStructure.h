/*!
  @file IceStructure.h
	
  @brief �X�\���𑀍삷��N���X
 
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
#include <math.h>
#include <cuda_runtime.h>

#include <UtilityScript\mk_Vector2D.h>
#include <UtilityScript\mk_Vector3D.h>
#include "QueryCounter.h"

#include <IceTetrahedra.h>
#include "Ice_SM.h"

using namespace std;

class IceStructure
{
public:	//TODO: �S��public�ɂ��Ȃ�
	IceStructure();
	IceStructure(int pNum, int cNum, int tNum, int layer);

	~IceStructure(void);

	static int* GetDevicePtoCIndxPointer(){	return sd_piPtoCIndx;	}
	static int* GetDevicePtoCPointer(){		return sd_piPtoC;	}

	//������
	void InitTetraInfo();									//�l�ʑ̏��̃������m��
	void InitClusterInfo();									//�N���X�^���̃������m��
	void InitSelectCluster(vector<Ice_SM*>& iceSM);
	void InitSelectClusterFromClusterSet(vector<Ice_SM*>& iceSM);

	void InitGPU();											//GPU�����ŗp����f�[�^�̏�����

	//���ω�
	//�Z������
	void StepObjMelt(vector<unsigned>& pList, vector<unsigned>& cList, vector<unsigned>& tList, vector<unsigned>& cLayerList, vector<unsigned>& tLayerList);

	void SearchMeltParticle(vector<int>& pList);	
	void SearchReconstruct_Cluster_Melt(const vector<unsigned>& pList, vector<unsigned>& cList, vector<unsigned>& lList);
	void SearchReconstruct_Tetra_Melt(const vector<unsigned>& pList, vector<unsigned>& tList, vector<unsigned>& lList);

	void UpdateInfo_Melt_PandT(const vector<unsigned>& pList);
	void UpdateInfo_Melt_PandC(const vector<unsigned>& pList, const vector<unsigned>& cList);

	void UpdateInfo_Delete_TandP(const vector<int>& tList, const vector<int>& deleteList);

	void SetInfo_Cluster(const vector<unsigned>& pList, const vector<unsigned>& cList, const vector<unsigned>& lList);
	void SetInfo_Tetra(const vector<unsigned>& pList, const vector<unsigned>& tList, const vector<unsigned>& lList);

	void CheckDeleteCluster(void);
	void CheckDeleteTetra(vector<int>& tList, vector<int>& lList);
	void CheckSameTetra(int tIndx, const vector<int>& searchList, vector<int>& deleteList);
	void CheckIncludeTetra(int tIndx, const vector<int>& searchList, vector<int>& deleteList);

	void RemoveReconstructTetra(vector<int>& tList, vector<int>& lList, vector<int>& deleteTList);

	//�Ìŏ���
	void StepObjFreeze();
	
	void SearchFreezeParticle(vector<int>& pList);
	void SearchReconstructCluster_Freeze(const vector<int>& pList, vector<int>& cList, vector<int>& lList);
	void SearchReconstructTetra_Freeze(const vector<int>& pList, vector<int>& tList, vector<int>& lList);

	void SetFreezeTetraInfo(vector<int>& pList);
	void SetFreezeClusterInfo(const vector<int>& pList);

	//�A�N�Z�b�T
	void SetParticleNum(int pNum){	m_iPNum = pNum; }		//���݂̗��q��
	void SetClusterNum(int cNum){	m_iCNum = cNum; }		//���݂̃N���X�^��
	void SetTetraNum(int tNum){		m_iTNum = tNum;	}		//���݂̎l�ʑ̐�

	int GetParticleNum(void){	return m_iPNum;	}
	int GetClusterNum(void){	return m_iCNum;	}
	int GetTetraNum(void){		return m_iTNum;	}

	int GetPNumMax(void){	return m_iPNumMax;	}
	int GetCNumMax(void){	return m_iCNumMax;	}

	int GetCtoPMax(void){	return m_iCtoPMax;	}
	int GetPtoCMax(void){	return m_iPtoCMax;	}

	void CountTetrahedra(int tIndx, vector<int>& pList);					//�l�ʑ̏��̓o�^
	void CountClusterParticle(int cIndx, vector<int>& pList, int pNum);		//�N���X�^�̏��o�^

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
	
	void SetTetraInfo(int tIndx, int* PtoTNum);								//�������ɗp����l�ʑ̏��̓o�^
	//void SetClusterInfo(int cIndx, int* PtoCNum, int 

	void SetNeighborTetra(int tIndx, int layer);
	void SetNeighborTetraFromLayer(int tIndx, int searchLayer, int deleteLayer);

	void SetSelectRadius(float radius){	m_selectRadius = radius;	}
	float GetSelectRadius(){	return m_selectRadius;	}

	//�폜
	void DeleteTtoP(int tIndx, int lIndx);
	void DeletePtoC(int pIndx, int lIndx);
	void DeletePtoT(int pIndx, int lIndx);

	int  GetTtoP(int tIndx, int lIndx);

	int GetPtoC(int pIndx, int lIndx, int oIndx);
	int GetPtoT(int pIndx, int lIndx, int oIndx);
	const int& GetCtoP(const int& cIndx, const int& lIndx, const int& oIndx);

	int GetNeighborTetra(int tIndx, int lIndx, int oIndx);

	//������
	void ClearPtoT(int pIndx);
	void ClearPtoC(int pIndx);
	void ClearCtoP(int cIndx);
	void ClearTtoP(int tIndx);

	void ClearNeighborTetra(int tIndx);
	void ClearNeighborTetraFromLayer(int tIndx, int layer);

	//����
	int  CheckNeighborTetra(int tIndx, int checkTIndx);

	//�^���v�Z�N���X�^�̑I��
	void UpdateSelectCluster(const vector<unsigned>& prtList, vector<unsigned>& neighborClusters, const vector<Ice_SM*>& iceSM);
	void UpdateSelectClusterFronSet(const vector<unsigned>& prtList, vector<unsigned>& neighborClusters, const vector<Ice_SM*>& iceSM);
	void UpdateNeighborOfSelectCluster(vector<Ice_SM*>& iceSM);
	void ResetSelectCluster(vector<Ice_SM*>& iceSM);

	void UpdateMotionCalcCluster(unsigned cIndx, short unsigned num);
	inline short unsigned GetMotionCalcCluster(unsigned cIndx) const{	return m_psuSelectClusterIndx[cIndx];	}

	//�f�o�b�O
	void DebugPtoT(int pIndx);
	void DebugPtoC(int pIndx);
	void DebugCtoP(int cIndx);
	void DebugTtoP(int tIndx);
	void DebugNeighborTetra(int tIndx);
	void DebugStepObjMelt(vector<unsigned>& pList, vector<unsigned>& cList);

//�e�X�g
	void TestStepObjMelt(vector<unsigned>& pList, vector<unsigned>& cList, vector<unsigned>& tList, vector<unsigned>& cLayerList, vector<unsigned>& tLayerList);

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

	//---------------------------------------_GPU------------------------------------------
	int* GetDeviceCtoPNumPointer(){	return sd_piCtoPNum;	}
	int* GetDeviceCtoPointer(){		return sd_piCtoP;		}

	//---------------------------------------GPU_------------------------------------------

protected:

	int m_iPNumMax;								//�ő嗱�q��
	int m_iPNum;								//���݂̗��q��
	
	int m_iPtoCMax;								//���q���N���X�^�ɏ�������ő吔
	int m_iCtoPMax;								//�N���X�^���܂ޗ��q�̍ő吔

	int m_iTNumMax;
	int m_iTNum;

	int m_iPtoTMax;
	int m_iTtoPMax;

	int m_iCNumMax;								//�ő�N���X�^��
	int m_iCNum;								//���݂̃N���X�^��

	int m_iNeighborMax;							//�ߖT���q�̍ő吔
	int m_iLayer;

	float m_selectRadius;						//�^���v�Z�N���X�^��I������ۂ̔��a

	//�Ƃ肠�����C�[���������z��Ƀ|�C���^�͎g��Ȃ�
	//���q��
	mk_Vector3D<int> m_mk3DiPtoC;				//���q���N���X�^
	mk_Vector3D<int> m_mk3DiPtoT;				//���q���l�ʑ�

	int*   m_piPtoCNum;							//���q���N���X�^�̌�
	int*   m_piPtoTNum;							//���q���l�ʑ̂̌�

	int*   m_piPtoCIndx;
	int*   m_piPtoTIndx;

	//�N���X�^��
	mk_Vector3D<int> m_mk3DiCtoP;				//�N���X�^�����q�@�N���X�^�ɏ������闱�q��Ԃ��@���q�̐ڑ����
	//mk_Vector3D<int> m_mk3DiCtoT;				//�N���X�^���l�ʑ�

	int*   m_piCtoPNum;							//�N���X�^�����q�̌�
	int*   m_piCtoTNum;							//�N���X�^���l�ʑ̂̌�

	int*   m_piCtoPIndx;
	int*   m_piCtoTIndx;

	//�l�ʑ́�
	mk_Vector2D<int> m_mk2DiTtoP;				//�l�ʑ́����q
	//mk_Vector3D<int> m_mk3DiTtoC;				//�l�ʑ́��N���X�^

	int*   m_piTtoPNum;							//�l�ʑ́����q�̌�
	int*   m_piTtoCNum;							//�l�ʑ́��N���X�^�̌�

	int*   m_piTtoPIndx;
	int*   m_piTtoCIndx;

	//�ߖT�l�ʑ�
	mk_Vector3D<int> m_mk3DiNeighborTetra;		//�ߖT�l�ʑ�

	int*   m_piNTNum;							//�e�ߖT�l�ʑ̂̌�

	short unsigned* m_psuSelectClusterIndx;		//�^���v�Z���闱�q�̃t���O�@�^���v�Z����Ȃ�1�C���Ȃ��Ȃ�0
												//�{����bool�ł������C����̂��Ƃ��l����short unsigned�ɂ��Ă�����

//--------------------------------------GPU__------------------------------------------------------------
	//���q��
	static int* sd_piPtoC;
	static int* sd_piPtoT;

	static int* sd_piPtoCNum;
	static int* sd_piPtoTNum;							//���q���l�ʑ̂̌�

	static int* sd_piPtoCIndx;
	static int* sd_piPtoTIndx;

	//�N���X�^��
	static int* sd_piCtoP;

	static int* sd_piCtoPNum;
	//static int* sd_piCtoTNum;
	
	static int* sd_piCtoPIndx;
	//static int*   sd_piCtoTIndx;
	
	//�l�ʑ́�
	static int* sd_piTtoP;

	static int* sd_piTtoPNum;
	static int* sd_piTtoCNum;

	static int* sd_piTtoPIndx;
	static int* sd_piTtoCIndx;

	//�ߖT�l�ʑ�
	static int* sd_piNeighborTetra;

	static int* sd_piNeighborTetraTNum;
//--------------------------------------__GPU------------------------------------------------------------

	//�T���p�t���O�@���g�p
	bool* m_pbPFlag;							//���q
	bool* m_pbCFlag;							//�N���X�^
	bool* m_pbTFlag;							//�l�ʑ�

	//���g�p
	vector<vector<int>>	m_vviEdgeToCluster;					//�Ӂ��N���X�^	����ӂ��ǂ̃N���X�^�ɏ������Ă��邩�@�K���P�ȏ�
	vector<vector<int>>	m_vviEdgeToParticle;				//�Ӂ����q�@�@�@����ӂ��ǂ̗��q��ڑ����Ă��邩�@�@�@�K���Q��

	vector<vector<int>>	m_vviClusterToEdge;					//�N���X�^���Ӂ@����N���X�^�ɂǂ̕ӂ��������Ă��邩�@�K���P���R��
	vector<vector<int>>	m_vviParticleToEdge;				//���q���Ӂ@�@�@���闱�q�͂ǂ̕ӂŐڑ�����Ă��邩

};
#endif