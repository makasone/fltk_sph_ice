/*!
  @file Sur_SM.h

  @brief ShapeMatching�@����ɂ������ω��V�~�����[�V����
  @ref R. Diziol et al., "Robust Real-Time Deformation of Incompressible Suface Meshes", SIGGRAPH2005. 
 
  @author Ryo Nakasone
  @date 2013-10
*/

#ifndef _SURF_SM_H
#define _SURF_SM_H

#include <iostream>
#include <vector>
#include <algorithm>
#include <cuda_runtime.h>

#include <Ice_SM.h>
#include <IceStructure.h>

#include "rx_utility.h"
#include "rx_matrix.h"

#include <UtilityScript\mk_Vector2D.h>
#include <UtilityScript\mk_Vector3D.h>
#include "QueryCounter.h"

using namespace std;

//�\�ʗ��q�݂̂�p����SM�@�̂��߂̃N���X�D
//��������̃p�X���Ǘ����āC�������̂��߂̃p�����[�^���v�Z���ĕԂ��D

extern void LaunchUpdatePrefixSumGPU
(
	int prtNum,
	int PosSizeX,
	int PosSizeY,
	int ApqSizeX,
	int ApqSizeY,
	int* md_2DiPTHtoPRT,
	int* md_2DiPRTtoPTH,
	float* md_2Df3PrfxPos,
	float* md_2Df9PrfxApq,
	int* md_3DiPTHandPrfxSet,
	float* md_f3OrgPos,
	float* md_f3OrgCm,
	unsigned* dgroupPos,
	unsigned* dgroupApq,
	const float* md_fPos,
	const float* md_fVel
);

class Surf_SM
{
private:
	mk_Vector2D<Vec3>		m_mk2Dvec3_PrfxPos;		//�p�X���Ƃ́C���݂̈ʒu���̂��߂�prefixSum
	mk_Vector2D<rxMatrix3>	m_mk2Dmat3_PrfxApq;		//�p�X���Ƃ́C�ό`�s��̂��߂�prefixSum�@�i���S�ȕό`�s��Apq�ł͂Ȃ���@�ڂ����͘_�����Q�Ɓj
	
	mk_Vector2D<int>		m_mk2DiPTHtoPRT;		//�p�X�����q�@�s�K���z��ɂȂ�
	mk_Vector2D<int>		m_mk2DiPRTtoPTH;		//���q���p�X�@0:�p�X�ԍ��C1:�p�X���ԍ��@���q�͂P�̃p�X�ɂ��������Ȃ�

	mk_Vector3D<int>		m_mk3DiPTHandPrfxSet;	//�e�N���X�^�ɂ�����C�p�X��prefixSum�̔Ԓn�Z�b�g[0]�F�n�_�@[1]�F�I�_�@prefixSum�Ԓn�݂̂ł����@path�ԍ��͗��q����o�R���Ď擾�ł���

	vector<Vec3> m_vvec3OrgPos;						//���q�̏����ʒu
	vector<Vec3> m_vvec3OrgCm;						//�N���X�^�̏����d�S

	const float* m_fPos;							//���q�̈ʒu�ւ̃|�C���^
	const float* m_fVel;							//���q�̑��x�ւ̃|�C���^

	int m_iPrtclNum;

	IceStructure* m_strct;
	vector<Ice_SM*> m_iceSM;

//---------------------------------------------GPU__---------------------------------------------------------
	float* md_2Df3PrfxPos;							//�p�X���Ƃ́C���݂̈ʒu���̂��߂�prefixSum
	float* md_2Df9PrfxApq;							//�p�X���Ƃ́C�ό`�s��̂��߂�prefixSum�@�i���S�ȕό`�s��Apq�ł͂Ȃ���@�ڂ����͘_�����Q�Ɓj
	
	int* md_2DiPTHtoPRT;							//�p�X�����q�@�s�K���z��ɂȂ�
	int* md_2DiPRTtoPTH;							//���q���p�X�@0:�p�X�ԍ��C1:�p�X���ԍ��@���q�͂P�̃p�X�ɂ��������Ȃ�

	int* md_3DiPTHandPrfxSet;						//�e�N���X���ɂ�����C�p�X��prefixSum�̔Ԓn�Z�b�g[0]�F�n�_�@[1]�F�I�_�@prefixSum�Ԓn�݂̂ł����@path�ԍ��͗��q����o�R���Ď擾�ł���

	float* md_f3OrgPos;								//���q�̏����ʒu

	float* md_f3ClusterOrgPos;						//�N���X�^���̗��q�̏����ʒu		Ice_SM�̃f�o�C�X�|�C���^���R�s�[���Ďg��
	float* md_f3ClusterOrgCm;						//�N���X�^�̏����d�S				Ice_SM�̃f�o�C�X�|�C���^���R�s�[���Ďg��

	unsigned* md_uPosGroup;							//prefixSum�̃O���[�v���Ɏg�����߂̃t���O
	unsigned* md_uApqGroup;							//prefixSum�̃O���[�v���Ɏg�����߂̃t���O

	const float* md_fPos;							//���q�̈ʒu�ւ̃|�C���^
	const float* md_fVel;							//���q�̑��x�ւ̃|�C���^
//---------------------------------------------__GPU---------------------------------------------------------

public:
	void InitPath(const float* pos, const float* vel, const vector<Ice_SM*> iceSM, IceStructure* strct, int pthSize, int prtNum);	//�p�X�쐬
	void InitPathGPU();
	
	void InitPathPrfxIndxSet(const vector<Ice_SM*> iceSM, IceStructure* strct);				//�ǂ̃p�X�̂ǂ̕������K�v�Ȃ̂��C���N���X�^���ƂɌv�Z
	void InitOrgPos(int prtNum);
	void InitOrgCm();

	void UpdatePrefixSum();
	void UpdatePrefixSumItr();
	void UpdatePrefixSumGPU();

	void UpdatePrefixSumPos();		//�d�S�@PrefixSum�̌v�Z
	void UpdatePrefixSumApq();		//�ό`�s��@PrefixSum�̌v�Z

	void UpdatePrefixSumPosItr();
	void UpdatePrefixSumApqItr();

	Vec3 CalcCmSum(const int& cIndx);											//prefixSum����N���X�^�̏d�S���v�Z���ĕԂ�
	const Vec3 CalcCmFromPrfxSm(const int& path, const int& start, const int& end);

	rxMatrix3 CalcApqSum(const int& cIndx);
	const rxMatrix3 CalcApqFromPrfxSm(const int& path, const int& start, const int& end);

//�A�N�Z�b�T
	int* GetDevicePRTtoPTHPointer(){		return md_2DiPRTtoPTH;		}
	int* GetDevicePTHandPrfxSetPointer(){	return md_3DiPTHandPrfxSet;	}
	float* GetDecvicePrfxPos(){	return md_2Df3PrfxPos;	}
	float* GetDecvicePrfxApq(){	return md_2Df9PrfxApq;	}

	void SetDevicePointer(const float* dPos, const float* dVel){	md_fPos = dPos; md_fVel = dVel;	}


//�f�o�b�O
	void DebugInit();
	void DebugPathDataPos();
	void DebugPathDataApq();
	void DebugPathPrfxIndxSet();

//�e�X�g
	void TestUpdatePrefixSum();

};

#endif