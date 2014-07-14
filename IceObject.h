//�S���ω��I�u�W�F�N�g���Ǘ�����N���X

#ifndef _ICE_OBJECT_
#define _ICE_OBJECT_

#include "Ice_SM.h"
#include "IceStructure.h"
#include "QueryCounter.h"

#include <time.h>
#include <cuda_runtime.h>
#include <cuda_gl_interop.h>

using namespace std;

class IceObject
{
private:
	//TODO::�t�̉^���v�Z�N���X���u������
	//TODO::�l�ʑ̏������u��
	static float* s_sphPrtPos;				//�ʒu�̃|�C���^�@�ǂݍ��ݐ�p
	static float* s_sphPrtVel;				//���x�̃|�C���^�@�ǂݍ��ݐ�p

//--------------------------------------GPU__------------------------------------------------------------
	static float* sd_sldPrtPos;						//���q�ʒu�̃f�o�C�X�|�C���^
	static float* sd_sldPrtVel;						//���q���x�̃f�o�C�X�|�C���^

	static float* sd_ObjPrtPos;						//���a�v�Z�ɂ��ŏI�I�ȗ��q�ʒu
	static float* sd_ObjPrtVel;						//���x
//--------------------------------------__GPU------------------------------------------------------------

	static int sm_particleNum;						//���q��
	static int sm_tetraNum;							//�l�ʑ̐�
	static int sm_clusterNum;						//�N���X�^��

	//�ő̉^���v�Z�N���X
	vector<Ice_SM*> m_iceMove;

	//�\���Ǘ��N���X
	IceStructure* m_iceStrct;

	static float* m_fInterPolationCoefficience;			//���`��ԌW��

public:
	IceObject(float* pos, float* vel, int pMaxNum, int cMaxNum, int tMaxNum);
	~IceObject();

	static void SetSPHPointer(float* pos, float* vel){		s_sphPrtPos = pos;	s_sphPrtVel = vel;	};

	void InitIceObj(int pMaxNum, int cMaxNum, int tMaxNum);
	void InitCluster(Ice_SM* sm){	m_iceMove.push_back(sm);	}	//�|�C���^���R�s�[���Ă��邾���@�ꎞ�I�Ȏ���
	static void InitInterPolation();

	void StepObjMove();
	void StepInterPolation();								//���`��ԁ@������͏��������G�ɂȂ�̂ŃN���X�ɂ������D
	
	void CalcAverageCPU(const int pIndx, Vec3& pos, Vec3& vel);
	void LinerInterPolationCPU(const int pIndx, const Vec3& pos, const Vec3& vel);

	//--------------IceStructure�Ɠ������������邽�߂Ɉꎞ�I�ɍ�����֐�__----------------------------------
	void SetParticleNum(int pNum){	sm_particleNum = pNum;	m_iceStrct->SetParticleNum(pNum);	}
	void SetTetraNum(int tNum){		sm_tetraNum = tNum;	m_iceStrct->SetTetraNum(tNum);		}
	void SetClusterNum(int cNum){	sm_clusterNum = cNum; m_iceStrct->SetClusterNum(cNum);	}

	void InitTetraInfo(){	m_iceStrct->InitTetraInfo();	}
	void InitClusterInfo(){	m_iceStrct->InitClusterInfo();	}

	int GetPtoTNum(int i){	return m_iceStrct->GetPtoTNum(i);	}
	int GetTtoPNum(int i){	return m_iceStrct->GetTtoPNum(i);	}
	int GetPtoCNum(int i){	return m_iceStrct->GetPtoCNum(i);	}
	int GetCtoPNum(int i){	return m_iceStrct->GetCtoPNum(i);	}

	int GetPtoTIndx(int pIndx){	return m_iceStrct->GetPtoTIndx(pIndx);	}
	int GetTtoPIndx(int tIndx){ return m_iceStrct->GetTtoPIndx(tIndx);	}
	int GetPtoCIndx(int i){	return m_iceStrct->GetPtoCIndx(i);	}
	int GetCtoPIndx(int i){	return m_iceStrct->GetCtoPIndx(i);	}

	int GetPtoC(int i, int j, int k){	return m_iceStrct->GetPtoC(i, j, k);	}
	int GetPtoT(int i, int j, int k){	return m_iceStrct->GetPtoT(i, j, k);	}
	int GetTtoP(int i, int j){			return m_iceStrct->GetTtoP(i, j);		}

	int GetParticleNum(){	return m_iceStrct->GetParticleNum();	}

	int GetPtoCMax(){ return m_iceStrct->GetPtoCMax();	}

	int GetNTNum(int i){ return m_iceStrct->GetNTNum(i);	}

	int GetNeighborTetra(int i, int j, int k){	return m_iceStrct->GetNeighborTetra(i, j, k);	}

	void SetPtoTIndx(int pIndx){	m_iceStrct->SetPtoTIndx(pIndx, m_iceStrct->GetPtoTNum(pIndx));	}
	void SetTtoPIndx(int tIndx){	m_iceStrct->SetTtoPIndx(tIndx, m_iceStrct->GetTtoPNum(tIndx));	}
	void SetPtoCIndx(int i, int j){	m_iceStrct->SetPtoCIndx(i, j);	}
	void SetCtoPIndx(int i, int j){	m_iceStrct->SetCtoPIndx(i, j);	}

	void SetPtoT(int i, int j, int k, int l){	m_iceStrct->SetPtoT(i, j, k, l);	}
	void SetTtoP(int i, vector<int> j)	{		m_iceStrct->SetTtoP(i, j);			}
	void SetPtoC(int i, int j, int k, int l, int m){	m_iceStrct->SetPtoC(i, j, k, l, m);	}
	void SetCtoP(int cIndx, const vector<int>& pIndxList, int* pLayerList){	m_iceStrct->SetCtoP(cIndx, pIndxList, pLayerList);	}

	void SetNeighborTetra(int i, int layer){	m_iceStrct->SetNeighborTetra(i, layer);	}

	void CountPtoT(int pIndx){	m_iceStrct->CountPtoT(pIndx);	}
	void CountTtoP(int tIndx){	m_iceStrct->CountTtoP(tIndx);	}
	void CountPtoC(int pIndx){	m_iceStrct->CountPtoC(pIndx);	}
	void CountCtoP(int cIndx){	m_iceStrct->CountCtoP(cIndx);	}
	//--------------__IceStructure�Ɠ������������邽�߂Ɉꎞ�I�ɍ�����֐�----------------------------------

};

#endif