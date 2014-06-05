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

//

class Ice_SM : public rxShapeMatching
{
protected:
	vector<int> m_iParticleIndxes;		//!< �N���X�^�ɏ������闱�q�̔ԍ�

	vector<double> m_dAlphas;			//!< stiffness�p�����[�^[0,1] (���x�v�Z�Ɏg�p)
	vector<double> m_dBetas;			//!< deformation�p�����[�^[0,1]

	Vec3 m_vec3OrgCm;					//�����̃N���X�^�̏d�S
	Vec3 m_vec3NowCm;					//���݂̃N���X�^�̏d�S
	//vector<Vec3> m_vvec3OrgQ;			//�����̈ʒu-�d�S

	Vec3 m_vec3DisCm;					//�d�S�ʒu�̕ψ�

	rxMatrix3	m_mtrx3Apq;				//�ό`�s��Apq
	rxMatrix3	m_mtrx3AqqInv;			//�ό`�s��Aqq�̋t�s��	�O�v�Z�\

	vector<int> m_iLayeres;

	vector<int> m_iLinearDeformation;	//!< Linear/Quadratic deformation�؂�ւ��t���O
	vector<int> m_iVolumeConservation;	//!< �ό`���̑̐ϕۑ���(��det(A)�Ŋ��邩�ǂ���)

	static const float* s_pfPrtPos;		//�ǂݍ��ݐ�p
	static const float* s_pfPrtVel;		//�ǂݍ��ݐ�p

	//�����ɕӂƖʂ���������
		//bounds�ǂɂԂ������Ƃ��̔�����
		//allowFlip���]���������̃t���O�@����Ȃ��H

public:
	Ice_SM(int obj);
	~Ice_SM();

	static void SetParticlePosAndVel(const float* pos, const float* vel)
	{
		s_pfPrtPos = pos;	s_pfPrtVel = vel;
	}

	void AddVertex(const Vec3 &pos, double mass, int pIndx);
	
	void Update();
	void ShapeMatching(double dt);
	void calExternalForces(float* newPos, double dt);
	void integrate(float* newPos, double dt);

	void ShapeMatchingSolid(float* newPos, double dt);
	void calExternalForcesSolid(double dt);
	void integrateSolid(double);


	void SetAlphas(int indx, int alpha){ m_dAlphas[indx] = alpha;	}
	void SetBetas (int indx, int beta){	m_dBetas[indx] = beta;		}

	void SetLayer(int indx, int layer){	m_iLayeres[indx] = layer;	}

	void SetNowCm(Vec3 nowCm){	m_vec3NowCm = nowCm;	}
	void SetApq(rxMatrix3 Apq){	m_mtrx3Apq = Apq;	}

	void SetLinerFalg(int indx, int flag){	m_iLinearDeformation[indx] = flag; }
	void SetVolumeFlag(int indx, int flag){	m_iVolumeConservation[indx] = flag;}

	int GetParticleIndx(int indx){ return m_iParticleIndxes[indx]; }
	const vector<int>& GetVertexIndxList(){ return m_iParticleIndxes;	}

	double GetAlphas(int indx){	return m_dAlphas[indx];	}
	double GetBetas (int indx){	return m_dBetas[indx];	}

	Vec3 GetCm(void){		return m_vec3NowCm;	}
	Vec3 GetOrgCm(void){	return m_vec3OrgCm;	}

	rxMatrix3 GetApq(void){	return m_mtrx3Apq;	}

	int GetLayer(int indx){	return m_iLayeres[indx];	}

	void Remove(int indx);
	void Clear();

	//������g���H
	int GetLinerFlags(int indx){	return m_iLinearDeformation[indx]; }
	int GetVolumeFlags(int indx){	return m_iVolumeConservation[indx]; }

	bool CheckIndx(int pIndx);
	int	 SearchIndx(int pIndx);

	const Vec3& GetDisVec(){	return m_vec3DisCm;	}
	void CalcDisplaceMentVectorCm();

	//�f�o�b�O
	void DebugIndx(void);
	void DebugLayer(void);
};

#endif