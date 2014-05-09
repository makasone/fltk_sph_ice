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

	vector<Vec3> m_vvec3OrgPos;						//���q�̏����ʒu	������񂾂��ۑ����Ă���
	//vector<Vec3> m_vvec3OrgCm;						//�N���X�^�̏����d�S	�Z���E�ÌŎ��ɂ͍X�V

public:
	void MakePath(const float* pos, int prtNum, int pthSize);	//�p�X�쐬

	void CalcPrefixSumPos(const float* pos);	//�d�S�@PrefixSum�̌v�Z
	void CalcPrefixSumApq(const float* pos);			//�ό`�s��@PrefixSum�̌v�Z
	
	void SetPathDataApq();

	void InitOrgPos(const float* pos, int pNum);

	Vec3 GetPos();									//����͈͂ɂ�����d�S�x�N�g���̑��a��Ԃ�
	rxMatrix3 GetApq();								//����͈͂ɂ�����ό`�s��̑��a��Ԃ�

	int GetPath(int path, int indx);				//�p�X�ɏ������闱�q�ԍ���Ԃ��@�p�X�ԍ��C����

	//�f�o�b�O
	void DebugPathDataPos();
	void DebugPathDataApq();
};

#endif