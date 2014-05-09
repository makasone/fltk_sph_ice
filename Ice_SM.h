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

using namespace std;

//

class Ice_SM : public rxShapeMatching
{
protected:
	vector<int> m_iParticleIndxes;		//!< �N���X�^�ɏ������闱�q�̔ԍ�

	vector<double> m_dAlphas;			//!< stiffness�p�����[�^[0,1] (���x�v�Z�Ɏg�p)
	vector<double> m_dBetas;			//!< deformation�p�����[�^[0,1]

	Vec3 m_vec3Cm;						//�N���X�^�̏d�S
	rxMatrix3	m_mtrx3Apq;				//�ό`�s��

	vector<int> m_iLayeres;

	vector<int> m_iLinearDeformation;	//!< Linear/Quadratic deformation�؂�ւ��t���O
	vector<int> m_iVolumeConservation;	//!< �ό`���̑̐ϕۑ���(��det(A)�Ŋ��邩�ǂ���)

	//�����ɕӂƖʂ���������
		//bounds�ǂɂԂ������Ƃ��̔�����
		//allowFlip���]���������̃t���O�@����Ȃ��H

public:
	Ice_SM(int obj);
	~Ice_SM();

	void AddVertex(const Vec3 &pos, double mass, int pIndx);

	void ShapeMatching(double dt);
	void Update();

	void SetAlphas(int indx, int alpha){ m_dAlphas[indx] = alpha;	}
	void SetBetas (int indx, int beta){	m_dBetas[indx] = beta;		}

	void SetLayer(int indx, int layer){	m_iLayeres[indx] = layer;	}

	void SetLinerFalg(int indx, int flag){	m_iLinearDeformation[indx] = flag; }
	void SetVolumeFlag(int indx, int flag){	m_iVolumeConservation[indx] = flag;}

	int GetParticleIndx(int indx){ return m_iParticleIndxes[indx]; }
	const vector<int>& GetVertexIndxList(){ return m_iParticleIndxes;	}

	double GetAlphas(int indx){	return m_dAlphas[indx];	}
	double GetBetas (int indx){	return m_dBetas[indx];	}

	Vec3 GetCm(void){		return m_vec3Cm;	}

	rxMatrix3 GetApq(void){	return m_mtrx3Apq;	}

	int GetLayer(int indx){	return m_iLayeres[indx];	}

	void Remove(int indx);
	void Clear();

	//������g���H
	int GetLinerFlags(int indx){	return m_iLinearDeformation[indx]; }
	int GetVolumeFlags(int indx){	return m_iVolumeConservation[indx]; }

	bool CheckIndx(int pIndx);



	//�f�o�b�O
	void DebugIndx(void);
	void DebugLayer(void);
};

#endif