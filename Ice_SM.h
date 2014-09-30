/*!
  @file Ice_SM.h
	
  @brief ShapeMatching�@����ɂ������ω��V�~�����[�V����
  @ref M. Muller et al., "Meshless deformations based on shape matching", SIGGRAPH2005. 
 
  @author Ryo Nakasone
  @date 2013-10
*/

#ifndef _ICE_SM_
#define _ICE_SM_

#include "ShapeMatching.h"
#include "QueryCounter.h"

#include <time.h>

using namespace std;

//GPU����
extern void LaunchShapeMatchingGPU
(
	int prtNum,
	float* prtPos,
	float* prtVel, 
	float* orgPos,
	float* curPos,
	float* orgCm,
	float* curCm,
	float* vel,
	int* pIndxes, 
	int* d_IndxSet,
	float dt
);

extern void LaunchShapeMatchingIterationGPU
(
	int prtNum,
	float* prtPos, 
	float* prtVel,
	float* sldPos,
	float* sldVel, 
	float* orgPos,
	float* curPos,
	float* orgCm,
	float* curCm,
	float* vel,
	int* pIndxes, 
	int* d_IndxSet,
	float dt
);

class Ice_SM : public rxShapeMatching
{
protected:

	unsigned m_iIndxNum;						//�z��Ŏ����������߁A�������ɑΉ����邽�߂̍ő�Y���ԍ�

	float* m_fpAlphas;							//!< stiffness�p�����[�^[0,1] (���x�v�Z�Ɏg�p)
	float* m_fpBetas;							//!< deformation�p�����[�^[0,1]	���g�p
	
	int* m_ipLayeres;

	static int s_iIterationNum;					//������

	Vec3 m_vec3OrgCm;							//�����̃N���X�^�̏d�S
	Vec3 m_vec3NowCm;							//���݂̃N���X�^�̏d�S
	//vector<Vec3> m_vvec3OrgQ;					//�����̈ʒu-�d�S

	Vec3 m_vec3DisCm;							//�d�S�ʒu�̕ψ�

	rxMatrix3	m_mtrx3Apq;						//�ό`�s��Apq
	rxMatrix3	m_mtrx3AqqInv;					//�ό`�s��Aqq�̋t�s��	�O�v�Z�\

	vector<int> m_iLinearDeformation;			//!< Linear/Quadratic deformation�؂�ւ��t���O�@���g�p
	vector<int> m_iVolumeConservation;			//!< �ό`���̑̐ϕۑ���(��det(A)�Ŋ��邩�ǂ���)�@���g�p

	static const float* s_pfPrtPos;				//�ʒu�̃z�X�g�|�C���^�@�ǂݍ��ݐ�p
	static const float* s_pfPrtVel;				//���x�̃z�X�g�|�C���^�@�ǂݍ��ݐ�p

	static float* s_pfSldPos;					//�N���X�^�̍ŏI�I�Ȉʒu
	static float* s_pfSldVel;					//�N���X�^�̍ŏI�I�ȑ��x

//--------------------------------------GPU------------------------------------------------------------
	static float* sd_PrtPos;					//�ŏI���q�ʒu�̃f�o�C�X�|�C���^
	static float* sd_PrtVel;					//�ŏI���q���x�̃f�o�C�X�|�C���^

	static float* d_OrgPos;
	static float* d_CurPos;
	
	static float* d_OrgCm;						//�����d�S
	static float* d_CurCm;						//���݂̏d�S�@�܂��g��Ȃ�

	static float* d_Mass;
	static float* d_Vel;

	static bool* d_Fix;

	static int* d_PIndxes;
	static int* d_IndxSet;						//�N���X�^�̃f�[�^�̊J�n�Y���ƏI���Y����ۑ�

	static int s_vertNum;						//�S�ẴN���X�^���ΏۂƂ��闱�q��
	static int s_vertSum;						//�S�N���X�^�Ɋ܂܂�闱�q�̑���

//--------------------------------------GPU------------------------------------------------------------

		//bounds�ǂɂԂ������Ƃ��̔�����
		//allowFlip���]���������̃t���O

public:
	Ice_SM(int obj);
	~Ice_SM();
	
	static void Ice_SM::InitGPU(const vector<Ice_SM*>& sm, float* d_pos, float* d_vel, int prtNum);

	static void SetPrtPointerPosAndVel(const float* pos, const float* vel){	s_pfPrtPos = pos;	s_pfPrtVel = vel;	}
	static void SetDevicePosPointer(float* d_pos){	sd_PrtPos = d_pos;	}
	
	static float* GetSldPosPointer(){	return s_pfSldPos;	}
	static float* GetSldVelPointer(){	return s_pfSldVel;	}

	static float* GetDeviceSPHPosPointer(){	return sd_PrtPos;	}
	static float* GetDeviceSPHVelPointer(){	return sd_PrtVel;	}

	static float* GetDevicePosPointer(){	return d_CurPos;	}
	static float* GetDeviceVelPointer(){	return d_Vel;		}
	static float* GetOrgPosPointer(){		return d_OrgPos;	}
	static float* GetOrgCmPointer(){		return d_OrgCm;		}

	static int* GetDeviceIndexSetPointer(){	return d_IndxSet;	}
	static int	GetVertexNum(){				return s_vertNum;	}

	void InitGPU_Instance();
	static void InitFinalParamPointer(int vrtxNum);

	void AddVertex(const Vec3 &pos, double mass, int pIndx);
	
	void UpdateCPU();
	void UpdateUsePathCPU();

	static void UpdateGPU();
	static void UpdateUsePathGPU();
	static void UpdateIterationGPU(float* sldPos, float* sldVel);

	static void CalcAverage();
	void CopyDeviceToInstance(int num);

	void ShapeMatchingUsePath();
	void ShapeMatchingSolid();
	void ShapeMatchingIteration();

	void calExternalForces();
	void calExternalForcesIteration();
	
	void integrate(double dt);
	void integrateIteration();

	void SetAlphas(int indx, float alpha){	m_fpAlphas[indx] = alpha;	}
	void SetBetas (int indx, float beta){	m_fpBetas[indx] = beta;		}

	void SetLayer(int indx, int layer){	m_ipLayeres[indx] = layer;	}

	static void SetIterationNum(int itr){	s_iIterationNum = itr;	}

	void SetNowCm(Vec3& nowCm){	m_vec3NowCm = nowCm;	}
	void SetApq(rxMatrix3& Apq){	m_mtrx3Apq = Apq;	}

	void SetLinerFalg(int indx, int flag){	m_iLinearDeformation[indx] = flag; }
	void SetVolumeFlag(int indx, int flag){	m_iVolumeConservation[indx] = flag;}

	float GetAlphas(int indx){	return m_fpAlphas[indx];	}
	float GetBetas (int indx){	return m_fpBetas[indx];	}

	Vec3 GetCm(void){		return m_vec3NowCm;	}
	Vec3 GetOrgCm(void){	return m_vec3OrgCm;	}

	rxMatrix3 GetApq(void){	return m_mtrx3Apq;	}

	int GetLayer(int indx){	return m_ipLayeres[indx];	}
	static int GetIteration(){	return s_iIterationNum;	}

	void Remove(int indx);
	void Clear();

	//������g���H
	int GetLinerFlags(int indx){	return m_iLinearDeformation[indx]; }
	int GetVolumeFlags(int indx){	return m_iVolumeConservation[indx]; }

	bool CheckHole(int oIndx);
	bool CheckIndx(int pIndx);
	int	 SearchIndx(int pIndx);

	const Vec3& GetDisVec(){	return m_vec3DisCm;	}
	void CalcDisplaceMentVectorCm();

	//�f�o�b�O
	void DebugIndx(void);
	void DebugLayer(void);

};

#endif