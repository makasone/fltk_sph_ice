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
	//float* m_pos;									//���q�̌��ݒn�ւ̃|�C���^
	//vector<Vec3> m_vvec3OrgCm;					//�N���X�^�̏����d�S	�Z���E�ÌŎ��ɂ͍X�V

	IceStructure* m_strct;

public:
	void InitPath(const float* pos, const float* vel, const vector<Ice_SM*> iceSM, IceStructure* strct, int pthSize);	//�p�X�쐬
	void InitOrgPos(const float* pos, int pNum);
	void InitPathPrfxIndxSet(const vector<Ice_SM*> iceSM, IceStructure* strct);	//�ǂ̃p�X�̂ǂ̕������K�v�Ȃ̂��C���N���X�^���ƂɌv�Z
	
	void SetPathDataApq();				

	void UpdatePrefixSum(const float* p, const float* vel);
	void UpdatePrefixSumPos(const float* pos, const float* vel);		//�d�S�@PrefixSum�̌v�Z
	void UpdatePrefixSumApq(const float* pos);		//�ό`�s��@PrefixSum�̌v�Z
	
	Vec3 ClacCmSum(int cIndx, const float* pos);						//prefixSum����N���X�^�̏d�S���v�Z���ĕԂ�
	Vec3 CalcCmFromPrfxSm(int path, int start, int end);

	Vec3		GetPos();							//����͈͂ɂ�����d�S�x�N�g���̑��a��Ԃ�
	rxMatrix3	GetApq();							//����͈͂ɂ�����ό`�s��̑��a��Ԃ�
	int			GetPath(int path, int indx);		//�p�X�ɏ������闱�q�ԍ���Ԃ��@�p�X�ԍ��C����

	//�f�o�b�O
	void DebugPathDataPos();
	void DebugPathDataApq();
	void DebugPathPrfxIndxSet();
};

#endif