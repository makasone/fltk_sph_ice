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
	mk_Vector2D<Vec3>		m_vvvec_PrfxPos;		//�p�X���Ƃ́C���݂̈ʒu���̂��߂�prefixSum
	mk_Vector2D<rxMatrix3>	m_vvmat_PrfxApq;		//�p�X���Ƃ́C�ό`�s��̂��߂�prefixSum�@�i���S�ȕό`�s��Apq�ł͂Ȃ���@�ڂ����͘_�����Q�Ɓj
	
	mk_Vector2D<int>	m_mk2DiPTHoPRT;				//�p�X�����q
	vector<int>			m_viPRTtoPTH;				//���q���p�X�@���q�͂P�̃p�X�ɂ��������Ȃ�

public:
	void MakePath(const float* pos, int num, int size);						//�p�X�쐬

	void CalcPrefixSum();							//PrefixSum�̌v�Z

	Vec3 GetPos();									//����͈͂ɂ�����ʒu�x�N�g���̑��a��Ԃ�
	rxMatrix3 GetApq();								//����͈͂ɂ�����ό`�s��̑��a��Ԃ�

	int GetPath(int path, int indx);				//�p�X�ɏ������闱�q�ԍ���Ԃ��@�p�X�ԍ��C����
};

#endif