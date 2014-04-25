/*!
  @file rx_particle_on_surf.cpp
	
  @brief �A�֐��\�ʂɃp�[�e�B�N����z�u

  @ref A. Witkin and P. Heckbert, "Using particles to sample and control implicit surfaces", SIGGRAPH1994.
 
  @author Makoto Fujisawa
  @date   2013-06
*/


//-----------------------------------------------------------------------------
// �C���N���[�h�t�@�C��
//-----------------------------------------------------------------------------
#include "rx_particle_on_surf.h"


//-----------------------------------------------------------------------------
// rxParticleOnSurf�N���X�̎���
//-----------------------------------------------------------------------------
/*!
 * �R���X�g���N�^
 */
rxParticleOnSurf::rxParticleOnSurf()
{
	m_uNumParticles = 0;

	// �ߖT�T���Z��
	m_pNNGrid = new rxNNGrid(DIM);

	m_fAlpha = 6.0;
	m_fEh = 0.8*m_fAlpha;
	m_fRho = 15.0;
	m_fPhi = 15.0;
	m_fBeta = 10.0;

	m_fGamma = 4.0;
	m_fNu = 0.2;
	m_fDelta = 0.7;

	m_pFuncPtr = 0;
}

/*!
 * �p�[�e�B�N��
 * @param[in] vrts �p�[�e�B�N�������ʒu
 * @param[in] rad �p�[�e�B�N�����a
 */
void rxParticleOnSurf::Initialize(const vector<Vec3> &vrts, double rad, Vec3 minp, Vec3 maxp, 
								   Vec4 (*func)(void*, double, double, double), void* func_ptr)
{
	m_fParticleRadius = rad;
	m_fEffectiveRadius = 3.0*m_fParticleRadius;
	m_fSh = m_fParticleRadius;
	m_fSmax = 1.5*m_fSh;

	m_uNumParticles = (uint)vrts.size();

	m_fpFunc = func;
	m_pFuncPtr = func_ptr;

	// �������m��
	m_vPos.resize(m_uNumParticles*DIM, 0.0);
	m_vPs.resize(m_uNumParticles);


	// �����p�[�e�B�N���ʒu
	for(uint i = 0; i < m_uNumParticles; ++i){
		for(int j = 0; j < 3; ++j){
			m_vPos[DIM*i+j]= vrts[i][j];
		}
		m_vPs[i].Sigma = m_fSh;
	}

	// �����Z���ݒ�
	m_pNNGrid->Setup(minp, maxp, m_fEffectiveRadius, m_uNumParticles);
	m_vNeighs.resize(m_uNumParticles);

}

/*!
 * �m�ۂ����������̉��
 */
void rxParticleOnSurf::Finalize(void)
{
	m_vPos.clear();
	m_vPs.clear();
}

/*!
 * repulsion radius �� �̍X�V
 *  - "Using particles to sample and control implicit surfaces"�̎�(10)�`(13)
 * @param[in] dt �X�V�p���ԃX�e�b�v��
 */
void rxParticleOnSurf::updateSigma(double dt)
{
	RXREAL h = m_fEffectiveRadius;

	for(uint i = 0; i < m_uNumParticles; ++i){
		Vec3 pos0(m_vPos[DIM*i+0], m_vPos[DIM*i+1], m_vPos[DIM*i+2]);
		RXREAL si = m_vPs[i].Sigma;

		double D = 0.0, Ds = 0.0;
		for(vector<rxNeigh>::iterator itr = m_vNeighs[i].begin() ; itr != m_vNeighs[i].end(); ++itr){
			int j = itr->Idx;
			if(j < 0 || i == j) continue;

			Vec3 pos1(m_vPos[DIM*j+0], m_vPos[DIM*j+1], m_vPos[DIM*j+2]);

			Vec3 rij = pos1-pos0;

			RXREAL r = norm(rij);

			if(r <= h){
				double Eij = m_fAlpha*exp(-r*r/(2.0*si*si));
				D += Eij;			// ���݂̔����G�l���M(��(10)�̏�)
				Ds += r*r*Eij;		// �G�l���MD�̃Ђɂ�����(��(13))
			}
		}
		Ds /= si*si*si;

		double Dv = -m_fRho*(D-m_fEh);// �^�[�Q�b�g�G�l���M�ɋ߂Â��邽�߂̐��`�t�B�[�h�o�b�N(��(10))
		double sv = Dv/(Ds+m_fBeta);	// �Ђ̕ω���(��(12))

		// �Ђ̍X�V
		m_vPs[i].Sigma += sv*dt;
	}
}

/*!
 * �����ɂ�鑬�x�̌v�Z
 *  - �K���I�ȃp�[�e�B�N���̒ǉ�/�폜�̂��߂̃Ђ̍X�V���܂ރo�[�W����
 *  - "Using particles to sample and control implicit surfaces"��4.3�߁C��(9)
 * @param[in] dt �X�V�p���ԃX�e�b�v��
 */
void rxParticleOnSurf::applyRepulsion2(double dt)
{
	RXREAL h = m_fEffectiveRadius;

	// �ߖT�T���Z���փp�[�e�B�N�����i�[ & �ߖT�T��
	SetParticlesToCell();

	// repulsion radius �̍X�V
	updateSigma(dt);

	for(uint i = 0; i < m_uNumParticles; ++i){
		Vec3 pos0(m_vPos[DIM*i+0], m_vPos[DIM*i+1], m_vPos[DIM*i+2]);
		RXREAL si = m_vPs[i].Sigma;

		Vec3 Ep(0.0);
		for(vector<rxNeigh>::iterator itr = m_vNeighs[i].begin() ; itr != m_vNeighs[i].end(); ++itr){
			int j = itr->Idx;
			if(j < 0 || i == j) continue;

			Vec3 pos1(m_vPos[DIM*j+0], m_vPos[DIM*j+1], m_vPos[DIM*j+2]);

			Vec3 rij = pos1-pos0;

			RXREAL r = norm(rij);

			if(r <= h){
				double Eij = m_fAlpha*exp(-r*r/(2.0*si*si));
	
				RXREAL sj = m_vPs[j].Sigma;
				double Eji = m_fAlpha*exp(-r*r/(2.0*sj*sj));
				Vec3 rji = pos0-pos1;

				Ep += (rij/(si*si))*Eij-(rji/(sj*sj))*Eji;	// ��(9)
			}
		}
		Ep *= si*si;

		m_vPs[i].Ep = Ep;
	}
}

/*!
 * �����ɂ�鑬�x�̌v�Z
 *  - �Ђ̕ύX���܂܂Ȃ��V���v���ȃo�[�W����
 *  - "Using particles to sample and control implicit surfaces"��4.1��
 * @param[in] dt �X�V�p���ԃX�e�b�v��
 */
void rxParticleOnSurf::applyRepulsion(double dt)
{
	RXREAL h = m_fEffectiveRadius;

	// �ߖT�T���Z���փp�[�e�B�N�����i�[ & �ߖT�T��
	SetParticlesToCell();

	for(uint i = 0; i < m_uNumParticles; ++i){
		Vec3 pos0(m_vPos[DIM*i+0], m_vPos[DIM*i+1], m_vPos[DIM*i+2]);
		RXREAL si = m_vPs[i].Sigma;

		Vec3 Ep(0.0);
		for(vector<rxNeigh>::iterator itr = m_vNeighs[i].begin() ; itr != m_vNeighs[i].end(); ++itr){
			int j = itr->Idx;
			if(j < 0 || i == j) continue;

			Vec3 pos1(m_vPos[DIM*j+0], m_vPos[DIM*j+1], m_vPos[DIM*j+2]);

			Vec3 rij = pos1-pos0;

			RXREAL r = norm(rij);

			if(r <= h){
				double Eij = m_fAlpha*exp(-r*r/(2.0*si*si));	// �����G�l���M
				Ep += rij*Eij;
			}
		}

		m_vPs[i].Ep = Ep;
	}
}

/*!
 * �p�[�e�B�N���ʒu�X�V
 *  - �p�[�e�B�N�����A�֐��Ȗʏ�ɍڂ�悤�ɔ����ɂ�鑬�x���C�����āC�ʒu���X�V
 * @param[in] dt �X�V�p���ԃX�e�b�v��
 * @param[out] v_avg �ړ��ʂ�2�敽�ϒl
 */
void rxParticleOnSurf::applyFloating(double dt, double &v_avg)
{
	v_avg = 0.0;
	if(!m_uNumParticles) return;

	RXREAL h = m_fEffectiveRadius;
	for(uint i = 0; i < m_uNumParticles; ++i){
		Vec3 p(m_vPos[DIM*i+0], m_vPos[DIM*i+1], m_vPos[DIM*i+2]);	// ���݂̃p�[�e�B�N�����W
		Vec3 ve = m_vPs[i].Ep;	// �����ɂ�鑬�x

		Vec4 fv = m_fpFunc(m_pFuncPtr, p[0], p[1], p[2]);	// �p�[�e�B�N�����W�ɂ�����A�֐��l�Ƃ��̌��z���擾
		Vec3 fx(fv[0], fv[1], fv[2]);	// ���z
		RXREAL f = fv[3];				// �A�֐��l

		// �p�[�e�B�N���ڗ����x
		Vec3 v = ve-((dot(fx, ve)+m_fPhi*f)/(dot(fx, fx)))*fx;

		m_vPos[DIM*i+0] -= v[0]*dt;
		m_vPos[DIM*i+1] -= v[1]*dt;
		m_vPos[DIM*i+2] -= v[2]*dt;

		m_vPs[i].Vel = v;

		v_avg += norm2(v*dt);
	}
	v_avg /= m_uNumParticles;
}

/*!
 * �p�[�e�B�N���̕���ƍ폜����
 * @param[in] dt �X�V�p���ԃX�e�b�v��
 */
void rxParticleOnSurf::testFissionDeath(void)
{
	if(!m_uNumParticles) return;

	// �ߖT�T���Z���փp�[�e�B�N�����i�[ & �ߖT�T��
	SetParticlesToCell();

	int num_remove = 0, num_fission = 0;
	RXREAL h = m_fEffectiveRadius;
	for(uint i = 0; i < m_uNumParticles; ++i){
		Vec3 p(m_vPos[DIM*i+0], m_vPos[DIM*i+1], m_vPos[DIM*i+2]);	// �p�[�e�B�N�����W
		Vec3 v = m_vPs[i].Vel;	// �p�[�e�B�N�����x
		RXREAL si = m_vPs[i].Sigma;	// �������a��

		RXREAL vp = norm(v);

		// �p�[�e�B�N�������t��Ԃɋ߂����ǂ���
		if(vp < m_fGamma*si){
			// �����G�l���MD�̌v�Z
			double D = 0.0;
			for(vector<rxNeigh>::iterator itr = m_vNeighs[i].begin() ; itr != m_vNeighs[i].end(); ++itr){
				int j = itr->Idx;
				if(j < 0 || i == j) continue;

				Vec3 pos1(m_vPos[DIM*j+0], m_vPos[DIM*j+1], m_vPos[DIM*j+2]);
				Vec3 rij = pos1-p;
				RXREAL r = norm(rij);
				if(r <= h){
					double Eij = m_fAlpha*exp(-r*r/(2.0*si*si));
					D += Eij;			// �����G�l���M
				}
			}
			RXREAL R = RXFunc::Frand();	// [0,1]�̗���

			// �Ђ��傫������ or �G�l���M�[���K�؂ŃЂ��K��l�ȏ� �� ����(�p�[�e�B�N���ǉ�)
			if((si > m_fSmax) || (D > m_fNu*m_fEh && si > m_fSh)){
				m_vPs[i].Flag = 1;	// �ǉ��t���O��ON
				num_fission++;
			}
			// �Ђ����������� and �������g�����e�X�g��ʂ��� �� �폜
			else if((si < m_fDelta*m_fSh) && (R > si/(m_fDelta*m_fSh))){
				m_vPs[i].Flag = 2;	// �폜�t���O��ON
				num_remove++;
			}

		}

		// �\�ʂ��痣�ꂷ�����p�[�e�B�N�����폜
		Vec4 fv = m_fpFunc(m_pFuncPtr, p[0], p[1], p[2]);	// �p�[�e�B�N�����W�ɂ�����A�֐��l�Ƃ��̌��z���擾
		RXREAL f = fv[3];				// �A�֐��l
		if(fabs(f) > 2.0*m_fSmax){
			m_vPs[i].Flag = 2;	// �폜�t���O��ON
			num_remove++;
		}
	}

	// �p�[�e�B�N���폜
	if(num_remove){
		int cnt = 0;
		vector<rxSurfParticle>::iterator itr = m_vPs.begin();
		vector<RXREAL>::iterator jtr = m_vPos.begin();
		while(itr != m_vPs.end()){
			if(itr->Flag == 2){
				itr = m_vPs.erase(itr);
				jtr = m_vPos.erase(jtr, jtr+DIM);
				cnt++;
			}
			else{
				++itr;
				jtr += DIM;
			}
		}

		//cout << cnt << " particles are removed." << endl;
	}
}

/*!
 * �p�[�e�B�N���ʒu�X�V
 * @param[in] dt �X�V�p���ԃX�e�b�v��
 * @param[inout] num_iter �ő唽����+���ۂ̔�����
 * @param[inout] eps �C������p�ړ����e��+�덷
 */
void rxParticleOnSurf::Update(double dt, int &num_iter, RXREAL &eps)
{
	double v_avg;
	int k;
	for(k = 0; k < num_iter; ++k){
		// �����ɂ�鑬�x�̌v�Z
		applyRepulsion2(dt);

		// �p�[�e�B�N���ʒu�̍X�V
		applyFloating(dt, v_avg);

		// �p�[�e�B�N���ǉ�/�폜
		testFissionDeath();

		v_avg = sqrt(v_avg);
		if(v_avg < eps) break;
	}

	eps = v_avg;
	num_iter = k;
}


//-----------------------------------------------------------------------------
// �ߖT�T��
//-----------------------------------------------------------------------------
/*!
 * �ߖT���q�T��
 * @param[in] idx �T�����S�p�[�e�B�N���C���f�b�N�X
 * @param[in] prts �p�[�e�B�N���ʒu
 * @param[out] neighs �T�����ʊi�[����ߖT���R���e�i
 * @param[in] h �L�����a
 */
void rxParticleOnSurf::GetNearestNeighbors(int idx, RXREAL *prts, vector<rxNeigh> &neighs, RXREAL h)
{
	if(idx < 0 || idx >= (int)m_uNumParticles) return;

	Vec3 pos;
	pos[0] = prts[DIM*idx+0];
	pos[1] = prts[DIM*idx+1];
	pos[2] = prts[DIM*idx+2];

	if(h < 0.0) h = m_fEffectiveRadius;

	m_pNNGrid->GetNN(pos, prts, m_uNumParticles, neighs, h);
}

/*!
 * �ߖT���q�T��
 * @param[in] idx �T�����S�p�[�e�B�N���C���f�b�N�X
 * @param[out] neighs �T�����ʊi�[����ߖT���R���e�i
 * @param[in] h �L�����a
 */
void rxParticleOnSurf::GetNearestNeighbors(Vec3 pos, vector<rxNeigh> &neighs, RXREAL h)
{
	if(h < 0.0) h = m_fEffectiveRadius;

	m_pNNGrid->GetNN(pos, &m_vPos[0], m_uNumParticles, neighs, h);
}


/*!
 * �S�p�[�e�B�N���𕪊��Z���Ɋi�[
 */
void rxParticleOnSurf::SetParticlesToCell(RXREAL *prts, int n, RXREAL h)
{
	// �����Z���ɗ��q��o�^
	m_pNNGrid->SetObjectToCell(prts, n);

	// �ߖT���q�T��
	if(h < 0.0) h = m_fEffectiveRadius;
	for(int i = 0; i < (int)m_uNumParticles; i++){
		m_vNeighs[i].clear();
		GetNearestNeighbors(i, prts, m_vNeighs[i], h);
	}
}
void rxParticleOnSurf::SetParticlesToCell(void)
{
	SetParticlesToCell(&m_vPos[0], m_uNumParticles, m_fEffectiveRadius);
}