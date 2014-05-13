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

#include <Ice_SM.h>

class IceStructure;

#include "rx_utility.h"
#include "rx_matrix.h"

#include <UtilityScript\mk_Vector2D.h>
#include <UtilityScript\mk_Vector3D.h>

using namespace std;

//�\�ʗ��q�݂̂�p�����r�l�@�̂��߂̃N���X�D
//��������̃p�X���Ǘ����āC�������̂��߂̃p�����[�^���v�Z���ĕԂ��D

class Surf_SM
{
private:
	mk_Vector2D<Vec3>		m_mk2Dvec3_PrfxPos;		//�p�X���Ƃ́C���݂̈ʒu���̂��߂�prefixSum
	mk_Vector2D<rxMatrix3>	m_mk2Dmat3_PrfxApq;		//�p�X���Ƃ́C�ό`�s��̂��߂�prefixSum�@�i���S�ȕό`�s��Apq�ł͂Ȃ���@�ڂ����͘_�����Q�Ɓj
	
	mk_Vector2D<int>		m_mk2DiPTHtoPRT;		//�p�X�����q
	mk_Vector2D<int>		m_mk2DiPRTtoPTH;		//���q���p�X�@�p�X�ԍ��C�p�X���ԍ��@���q�͂P�̃p�X�ɂ��������Ȃ�

	mk_Vector3D<int>		m_mk3DiPTHandPrfxSet;	//�e�N���X���ɂ�����C�p�X��prefixSum�̔Ԓn�Z�b�g[0]�F�n�_�@[1]�F�I�_�@prefixSum�Ԓn�݂̂ł����@path�ԍ��͗��q����o�R���Ď擾�ł���

	vector<Vec3> m_vvec3OrgPos;						//���q�̏����ʒu	������񂾂��ۑ����Ă���
	vector<Vec3> m_vvec3OrgCm;						//�N���X�^�̏����d�S

	const float* m_fPos;									//���q�̈ʒu�ւ̃|�C���^
	const float* m_fVel;									//���q�̑��x�ւ̃|�C���^

	IceStructure* m_strct;
	vector<Ice_SM*> m_iceSM;

public:
	void InitPath(const float* pos, const float* vel, const vector<Ice_SM*> iceSM, IceStructure* strct, int pthSize);	//�p�X�쐬
	void InitPathPrfxIndxSet(const vector<Ice_SM*> iceSM, IceStructure* strct);	//�ǂ̃p�X�̂ǂ̕������K�v�Ȃ̂��C���N���X�^���ƂɌv�Z
	void InitOrgPos(int prtNum);
	void InitOrgCm();

	void UpdatePrefixSum();
	void UpdatePrefixSumPos();		//�d�S�@PrefixSum�̌v�Z
	void UpdatePrefixSumApq();		//�ό`�s��@PrefixSum�̌v�Z
	
	const Vec3 CalcCmSum(int cIndx);								//prefixSum����N���X�^�̏d�S���v�Z���ĕԂ�
	const Vec3 CalcCmFromPrfxSm(int path, int start, int end);

	const rxMatrix3 CalcApqSum(int cIndx);
	const rxMatrix3 CalcApqFromPrfxSm(int path, int start, int end);

	//�f�o�b�O
	void DebugPathDataPos();
	void DebugPathDataApq();
	void DebugPathPrfxIndxSet();
};

#endif