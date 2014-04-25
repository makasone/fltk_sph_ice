/*!
  @file rx_particle_on_surf.h
	
  @brief �A�֐��\�ʂɃp�[�e�B�N����z�u
 
  @author Makoto Fujisawa
  @date   2013-06
*/


#ifndef _RX_PARTICLE_ON_SURFACE_H_
#define _RX_PARTICLE_ON_SURFACE_H_


//-----------------------------------------------------------------------------
// �C���N���[�h�t�@�C��
//-----------------------------------------------------------------------------
// C�W��
#include <cmath>

// STL
#include <vector>
#include <string>

#include <iostream>

#include "rx_utility.h"
#include "rx_nnsearch.h"


//-----------------------------------------------------------------------------
// ��`
//-----------------------------------------------------------------------------
typedef unsigned int uint;

#ifndef RXREAL
	#define RXREAL float
#endif
#ifndef DIM
	#define DIM 4
#endif





//-----------------------------------------------------------------------------
// �A�֐����l�ʂɃp�[�e�B�N����z�u����N���X
//  -�Q�l: A. Witkin and P. Heckbert, "Using particles to sample and control implicit surfaces", SIGGRAPH1994.
//-----------------------------------------------------------------------------
class rxParticleOnSurf
{
	struct rxSurfParticle
	{
		Vec3 Vel;
		Vec3 Ep;
		RXREAL Sigma;
		int Flag;

		rxSurfParticle() : Vel(Vec3(0.0)), Ep(Vec3(0.0)), Sigma(0.0), Flag(0) {}
	};

protected:
	uint m_uNumParticles;

	RXREAL m_fParticleRadius;
	RXREAL m_fEffectiveRadius;		//!< �L�����a
	RXREAL m_fKernelRadius;			//!< �J�[�l���̉e���͈�

	vector<RXREAL> m_vPos;			//!< �ߖT�T�����[�`���ɓn�����߂Ɉʒu�����ʊǗ�
	vector<rxSurfParticle> m_vPs;	//!< �p�[�e�B�N�����

	// Repulsion�p�����[�^
	RXREAL m_fAlpha;	//!< repulsion amplitude
	RXREAL m_fEh;		//!< desired energy
	RXREAL m_fRho;		//!< �p�[�e�B�N����energy��m_fEh�ɕۂ��߂̌W��
	RXREAL m_fPhi;		//!< �p�[�e�B�N���ʒu���Ȗʏ�ɕۂ��߂̌W��
	RXREAL m_fBeta;		//!< zero-divide�h�~�p
	RXREAL m_fSh;		//!< desired repulsion radius (���[�U�w��)
	RXREAL m_fSmax;		//!< maximum repulsion radius

	RXREAL m_fGamma;	//!< equilibrium speed (�Ђɑ΂���{��)
	RXREAL m_fNu;		//!< �p�[�e�B�N������̂��߂̌W��
	RXREAL m_fDelta;	//!< �p�[�e�B�N���폜�̂��߂̌W��

	void *m_pFuncPtr;
	Vec4 (*m_fpFunc)(void*, double, double, double);

	rxNNGrid *m_pNNGrid;			//!< �����O���b�h�ɂ��ߖT�T��
	vector< vector<rxNeigh> > m_vNeighs;	//!< �ߖT�p�[�e�B�N��

public:
	//! �R���X�g���N�^
	rxParticleOnSurf();

	//! �f�X�g���N�^
	~rxParticleOnSurf(){}

	// �p�[�e�B�N��������&�j��
	void Initialize(const vector<Vec3> &vrts, double rad, Vec3 minp, Vec3 maxp, 
					Vec4 (*func)(void*, double, double, double), void* func_ptr);
	void Finalize(void);

	//! �p�[�e�B�N����
	int	GetNumParticles() const { return m_uNumParticles; }

	//! �p�[�e�B�N�����a
	float GetParticleRadius(void){ return m_fParticleRadius; }

	//! �p�[�e�B�N���f�[�^�̎擾
	RXREAL* GetPositionArray(void){ return &m_vPos[0]; }

	//! �p�[�e�B�N���ʒu�̍X�V
	void Update(double dt, int &num_iter, RXREAL &eps);

protected:
	//! �Ђ̍X�V
	void updateSigma(double dt);

	//! �����͂̌v�Z
	void applyRepulsion(double dt);
	void applyRepulsion2(double dt);

	//! �\�ʂɉ������ړ�
	void applyFloating(double dt, double &v_avg);

	//!< �p�[�e�B�N������/�j���̔���
	void testFissionDeath(void);

public:
	// �ߖT�擾
	void GetNearestNeighbors(Vec3 pos, vector<rxNeigh> &neighs, RXREAL h = -1.0);
	void GetNearestNeighbors(int idx, RXREAL *p, vector<rxNeigh> &neighs, RXREAL h = -1.0);

	// �����Z���Ƀp�[�e�B�N�����i�[
	void SetParticlesToCell(void);
	void SetParticlesToCell(RXREAL *prts, int n, RXREAL h);

};


#endif // _RX_PARTICLE_ON_SURFACE_H_

