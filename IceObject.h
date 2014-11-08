//�S���ω��I�u�W�F�N�g���Ǘ�����N���X

#ifndef _ICE_OBJECT_
#define _ICE_OBJECT_

#include "rx_sph.h"
#include "rx_utility.h"
#include "rx_matrix.h"

#include "Ice_SM.h"
#include "IceStructure.h"
#include "IceTetrahedra.h"
#include "HeatTransfar.h"
#include "Surf_SM.h"
#include "QueryCounter.h"

#include <time.h>
#include <cuda_runtime.h>
#include <cuda_gl_interop.h>

using namespace std;

//#define MK_USE_GPU
//#define USE_PATH
//#define USE_ITR
//#define USE_SELECTED

#define SELECTED 2

const int g_iterationNum = 1;

class IceObject
{
private:
	//TODO::�t�̉^���v�Z�N���X���u������
	//TODO::�l�ʑ̏������u��
	static float* s_sphPrtPos;						//sph���q�ʒu�̃z�X�g�|�C���^
	static float* s_sphPrtVel;						//sph���q���x�̃z�X�g�|�C���^

//--------------------------------------GPU__------------------------------------------------------------
	static float* sd_sphPrtPos;						//sph���q�ʒu�̃f�o�C�X�|�C���^
	static float* sd_sphPrtVel;						//sph���q���x�̃f�o�C�X�|�C���^

	static float* sd_sldPrtPos;						//���a�v�Z�ɂ��ŏI�I�ȗ��q�ʒu�̃f�o�C�X�|�C���^
	static float* sd_sldPrtVel;						//���a�v�Z�ɂ��ŏI�I�ȗ��q���x�̃f�o�C�X�|�C���^
//--------------------------------------__GPU------------------------------------------------------------

	static int sm_particleNum;						//���݂̗��q��
	static int sm_tetraNum;							//���݂̎l�ʑ̐�
	static int sm_clusterNum;						//���݂̃N���X�^��
	static int sm_layerNum;							//�T�����C���[��
	
	static int sm_maxParticleNum;					//�ő嗱�q��

	//�\���Ǘ��N���X
	IceStructure* m_iceStrct;

	//�ő̉^���v�Z�N���X
	vector<Ice_SM*> m_iceMove;

	//�����v�Z�p�N���X
	//TODO::�|�C���^�ɂ�����H
	Surf_SM m_SurfSm;

	//�M����
	HeatTransfar* m_heatTransfer;

	
	static float* m_fInterPolationCoefficience;		//���`��ԌW��

public:
	IceObject(int pMaxNum, int cMaxNum, int tMaxNum, int prtNum, float* hp, float* hv, float* dp, float* dv, int layer, int maxParticleNum);
	~IceObject();

	void InitIceObj(int pMaxNum, int cMaxNum, int tMaxNum);
	void InitIceObjGPU();
	void InitTetra();
	void InitCluster(Vec3 boundarySpaceHigh, Vec3 boundarySpaceLow, float timeStep, int itr);
	void InitStrct();
	void InitPath();
	void InitHeatTransfer(float effectiveRadius, float timeStep, float tempMax, float tempMin, float latentHeat, float cffcntHt, float cffcntTd);

	static void InitInterPolation();
	void InitGPU();

	void SetSPHDevicePointer(float* pos, float* vel){	sd_sphPrtPos = pos; sd_sphPrtVel = vel;	}
	void SetSPHHostPointer(float* pos, float* vel){		s_sphPrtPos = pos;	s_sphPrtVel = vel;	}
	void SetSearchLayerNum(int layer){					sm_layerNum = layer;				}
	void SetMaxParticleNum(int particleNum){			sm_maxParticleNum = particleNum;	}
	void SetAirTemp(float temp){						m_heatTransfer->setAirTemp(temp);	}

	void SetClusterMoveInfo(int pIndx);
	void SetClusterStrctInfo(int cIndx, int *PtoCNum);

	static int GetParticleNum(){	return sm_particleNum;		}
	static int GetClusterNum(){		return sm_clusterNum;		}
	static int GetTetrahedraNum(){	return sm_tetraNum;			}
	static int GetLayerNum(){		return sm_layerNum;			}

	static int GetMaxClusterNum(){	return sm_maxParticleNum;	}

	float* GetTemps(){	return m_heatTransfer->getTemps();}

	Ice_SM* GetMoveObj(int cIndx){	return m_iceMove[cIndx];	}

	void StepObjMoveCPU();									//�^���v�Z
	void StepObjMoveGPU();									//GPU�ɂ��^���v�Z

	void StepObjMoveUsePath();								//��������@��p�����^���v�Z
	void StepObjMoveGPUUsePath();

	void StepObjMoveIteration();							//����������p�����^���v�Z
	void StepObjMoveIterationUsePath();						//��������@��p���������^���v�Z
	void StepObjMoveIterationWeighted();					//�d�ݕt���{��������

	void StepObjMoveSelected();

	void StepObjCalcWidhIteration();						//�ő̂̉^���v�Z�C���a�v�Z�C��ԏ����@GPU�����@������������

	void StepInterPolation();								//���`��ԁ@������͏��������G�ɂȂ�̂ŃN���X�ɂ������D
	void StepInterPolationForCluster();
	void StepInterPolationForClusterWeighted();
	void StepInterPolationSelected();
	void StepWeightedInterPolation();
	
	void StepHeatTransfer(const int* surfParticles, const vector<vector<rxNeigh>>& neights, float floor, float effRadius, const float* pos, const float* dens);		//�M����
	void StepIceStructure();								//���ω�����
	void StepMelting();										//�Z��
	void StepFreezing();									//�Ì�

	void SearchMeltParticle(vector<unsigned>& pList);
	void SearchFreezeParticle(vector<unsigned>& pList);

	void CalcAverageCPU(int pIndx, Vec3& pos, Vec3& vel);
	void CalcAverageSelected(int pIndx, Vec3& pos, Vec3& vel);
	void CalcWeightedVector(int pIndx, Vec3& pos, Vec3& vel);

	void LinerInterPolationCPU(int pIndx, const Vec3& pos, const Vec3& vel);
	void LinerInterPolationForClusterCPU(const int pIndx, const Vec3& pos, const Vec3& vel);

	void WarmParticle(int pIndx, float temp, float heat){	m_heatTransfer->WarmParticle(pIndx, temp, heat);	}
	void MeltParticle(int pIndx){	m_heatTransfer->MeltParticle(pIndx);	}

	void ReConstructCluster(vector<unsigned>& particleList, vector<unsigned>& clusterList);

	void UpdateUnSelectedCluster(int pIndx, const Vec3& pos, const Vec3& vel, const vector<Vec3>& prePos);




	//�f�o�b�O
	void DebugTetraInfo();
	void DebugClusterInfo();
	void DebugObjMoveUsePathWithGPU();
	void DebugDeformationAmount();
	void DebugDeformationAverage();

	//�e�X�g
	void TestStepInterPolation();
	void TestUpdateSMFromPath();







	//--------------IceStructure�Ɠ������������邽�߂Ɉꎞ�I�ɍ�����֐�__----------------------------------
	//�����Ǝ�������ΑS��������
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

	int GetPrtclNum(){	return m_iceStrct->GetParticleNum();	}
	int GetTetraNum(){	return m_iceStrct->GetTetraNum();	}

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


//GPU����
//TODO: �e�N���X�^�̉^���v�Z���ʂƌő̍\���̃f�[�^���K�v�Ȃ̂ŁC�����ɂ����Ă���
extern void LaunchCalcAverageGPU
(
	int prtNum,
	float* sldPrtPos,
	float* sldPrtVel,
	float* sphPrtPos,
	float* sphPrtVel,
	float* smPrtPos,
	float* smPrtVel,
	int* smIndxSet,
	int* PtoCIndx,
	int* PtoC,
	int PNumMax,
	int PtoCMax,
	int PtoCParamSize
);

//TODO: ������̓N���X�ɂ��������C�Ƃ肠���������ɂ����Ă���
extern void LaunchInterPolationGPU
(
	int prtNum,
	float* sldPrtPos,
	float* sldPrtVel,
	float* sphPrtPos,
	float* sphPrtVel
);

//PrefixSum�̃f�[�^��p����SM�@�Ŏg���f�[�^���X�V
//���̏����́CIce_SM��Surf_SM�̃f�[�^��p����̂ŁC�����ɒu���K�v������
extern void LauchUpdateSMFromPath
(
	int prtNum,
	float* prtPos,
	float* prtVel, 
	//------------------SM----------------------
	float* orgPos,
	float* curPos,
	float* orgCm,
	float* curCm,
	float* clstrApq,
	float* vel,
	int* pIndxes, 
	int* startEndSet,
	//-----------------Path---------------------
	short int* PRTtoPTH,
	short int* PTHandPrfxSet,
	float* prfxPos,
	float* prfxApq,
	//-----------------Struct-------------------
	int* CtoP,
	int* CtoPNum,
	int CtoPSizeY,
	int CtoPSizeZ,

	float dt
);


#endif