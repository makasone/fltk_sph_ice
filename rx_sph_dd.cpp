/*!
  @file rx_sph_dd.cpp
	
  @brief SPH�@(Double Density)�̎���
	- S.Clavet et al., "Particle-based Viscoelastic Fluid Simulation", SCA2005, 2005. 
	- http://www.iro.umontreal.ca/labs/infographie/papers/Clavet-2005-PVFS/index.html
 
  @author Makoto Fujisawa
  @date   2008-10,2012-11
*/
// FILE --rx_sph_dd.cpp--


//-----------------------------------------------------------------------------
// �C���N���[�h�t�@�C��
//-----------------------------------------------------------------------------
#include "rx_sph.h"

#include "rx_cu_funcs.cuh"
#include "rx_cu_common.cuh"

#include "rx_pcube.h"


//-----------------------------------------------------------------------------
// �O���[�o���ϐ�
//-----------------------------------------------------------------------------
double g_fPresK = 2.0;


//-----------------------------------------------------------------------------
// rxSPH�N���X�̎���
//-----------------------------------------------------------------------------
/*!
 * �R���X�g���N�^
 * @param[in] use_opengl VBO�g�p�t���O
 */
rxDDSPH::rxDDSPH(bool use_opengl) : 
	rxParticleSystemBase(use_opengl), 
	m_hNrm(0), 
	m_hFrc(0),
	m_hVelOld(0),
	m_hDens(0), 
	m_hPres(0), 
	m_hSurf(0), 
	m_pBoundary(0)
{
	m_v3Gravity = Vec3(0.0, -9.80665*0.005, 0.0);

	m_fBuoyancy = 0.0f;

	// �Փ˔���p
	m_fDamping = 0.0;
	m_fRestitution = 0.0;

	// double density
	m_fK = (RXREAL)g_fPresK;
	m_fKnear = 2.0f*m_fK;
	m_fViscC = 2.0f;
	m_fViscBeta = 1.0f;
	m_fElasKspring = 0.3f;
	m_fPlasAlpha = 0.3f;
	m_fPlasGamma = 0.1f;	// [0, 0.2]���炢

	// �ߖT�T���Z��
	m_pNNGrid = new rxNNGrid(DIM);

	m_uNumParticles = 0;


	m_iColorType = RX_RAMP;
}

/*!
 * �f�X�g���N�^
 */
rxDDSPH::~rxDDSPH()
{
	Finalize();
}


/*!
 * �V�~�����[�V�����̏�����
 * @param[in] max_particles �ő�p�[�e�B�N����
 * @param[in] boundary_ext ���E�̑傫��(�e�ӂ̒�����1/2)
 * @param[in] dens �������x
 * @param[in] mass �p�[�e�B�N���̎���
 * @param[in] kernel_particle �L�����ah�ȓ��̃p�[�e�B�N����
 */
void rxDDSPH::Initialize(const rxSPHEnviroment &env)
{
	// MARK:Initialize
	RXCOUT << "[rxDDSPH::Initialize]" << endl;

	m_fInitDens        = env.dens;
	m_fMass            = env.mass;
	m_iKernelParticles = env.kernel_particles;

	RXREAL volume = env.max_particles*m_fMass/m_fInitDens;

	m_fEffectiveRadius = pow(((m_iKernelParticles*volume)/(env.max_particles*RX_PI)), 1.0/2.0);
	m_fParticleRadius  = 0.5f*m_fEffectiveRadius;

	//m_fInitDens = 10.0;
	m_fInitDensNear = 2.0*m_fInitDens;
	m_fInitMaxDens = m_fInitDens;

	//m_fMaxDens = 10.0;
	m_fSurfDens = 0.82;

	
	RXCOUT << "particle : " << endl;
	RXCOUT << " n_max = " << env.max_particles << endl;
	RXCOUT << " h = " << m_fEffectiveRadius << endl;
	RXCOUT << " r = " << m_fParticleRadius << endl;
	RXCOUT << " dens = " << m_fInitDens << endl;
	RXCOUT << " mass = " << m_fMass << endl;
	RXCOUT << " kernel_particles = " << m_iKernelParticles << endl;
	//RXCOUT << " volume = " << volume << endl;

	RXREAL h = m_fEffectiveRadius;
	RXREAL r = m_fParticleRadius;

	//
	// ���E�ݒ�
	//
	m_pBoundary = new rxSolidBox(env.boundary_cen-env.boundary_ext, env.boundary_cen+env.boundary_ext, -1);
	//m_pBoundary = new rxSolidSphere(Vec3(0.0, 0.1, 0.0), 0.25, -1);
	
	// �V�~�����[�V�������̑傫��
	m_v3EnvMin = m_pBoundary->GetMin();
	m_v3EnvMax = m_pBoundary->GetMax();
	RXCOUT << "simlation range : " << m_v3EnvMin << " - " << m_v3EnvMax << endl;

	Vec3 world_size = m_v3EnvMax-m_v3EnvMin;
	Vec3 world_origin = m_v3EnvMin;
	
	double expansion = 0.01;
	world_origin -= 0.5*expansion*world_size;
	world_size *= (1.0+expansion); // �V�~�����[�V�������S�̂𕢂��悤�ɐݒ�

	m_v3EnvMin = world_origin;
	m_v3EnvMax = world_origin+world_size;
	

	// �J�[�l���֐��̒萔
	m_fWpoly6	=  315.0/(64.0*RX_PI*pow(h, (RXREAL)9.0));
	m_fGWpoly6	= -945.0/(32.0*RX_PI*pow(h, (RXREAL)9.0));
	m_fLWpoly6	= -945.0/(32.0*RX_PI*pow(h, (RXREAL)9.0));

	m_uNumParticles = 0;
	m_bCalDens = true;

	Allocate(env.max_particles);
}

/*!
 * �������̊m��
 *  - �ő�p�[�e�B�N�����Ŋm��
 * @param[in] max_particles �ő�p�[�e�B�N����
 */
void rxDDSPH::Allocate(int max_particles)
{
	// MARK:Allocate
	assert(!m_bInitialized);

	//m_uNumParticles = max_particles;
	m_uMaxParticles = max_particles;
	
	unsigned int size  = m_uMaxParticles*DIM;

	unsigned int size1 = m_uMaxParticles*2;
	unsigned int mem_size  = sizeof(RXREAL)*size;
	unsigned int mem_size1 = sizeof(RXREAL)*size1;

	//
	// �������m��
	//
	m_hPos = new RXREAL[size];
	m_hVel = new RXREAL[size];
	m_hNrm = new RXREAL[size];
	m_hFrc = new RXREAL[size];
	m_hVelOld = new RXREAL[size];
	memset(m_hPos, 0, mem_size);
	memset(m_hVel, 0, mem_size);
	memset(m_hNrm, 0, mem_size);
	memset(m_hFrc, 0, mem_size);
	memset(m_hVelOld, 0, mem_size);

	m_hDens = new RXREAL[size1];
	m_hPres = new RXREAL[size1];
	memset(m_hDens, 0, mem_size1);
	memset(m_hPres, 0, mem_size1);

	m_hPredictPos = new RXREAL[size];
	memset(m_hPredictPos, 0, mem_size);
	m_hDist = new RXREAL[size1];
	memset(m_hDist, 0, mem_size1);

	m_hSurf = new uint[m_uMaxParticles];
	memset(m_hSurf, 0, sizeof(RXREAL)*m_uMaxParticles);

	m_hTmp = new RXREAL[m_uMaxParticles];
	memset(m_hTmp, 0, sizeof(RXREAL)*m_uMaxParticles);

	if(m_bUseOpenGL){
		m_posVBO = createVBO(mem_size);	
		m_colorVBO = createVBO(m_uMaxParticles*4*sizeof(RXREAL));

		SetColorVBO(RX_RAMP);
	}

	// �����Z���ݒ�
	m_pNNGrid->Setup(m_v3EnvMin, m_v3EnvMax, m_fEffectiveRadius, m_uMaxParticles);
	m_vNeighs.resize(m_uMaxParticles);

	m_bInitialized = true;
}

/*!
 * �m�ۂ����������̉��
 */
void rxDDSPH::Finalize(void)
{
	assert(m_bInitialized);

	// ���������
	if(m_hPos) delete [] m_hPos;
	if(m_hVel) delete [] m_hVel;
	if(m_hNrm) delete [] m_hNrm;
	if(m_hFrc) delete [] m_hFrc;
	if(m_hVelOld) delete [] m_hVelOld;

	if(m_hDens) delete [] m_hDens;
	if(m_hPres) delete [] m_hPres;

	if(m_hPredictPos) delete [] m_hPredictPos;
	if(m_hDist) delete [] m_hDist;

	if(m_hSurf) delete [] m_hSurf;

	if(m_hTmp) delete [] m_hTmp;

	if(m_bUseOpenGL){
		glDeleteBuffers(1, (const GLuint*)&m_posVBO);
		glDeleteBuffers(1, (const GLuint*)&m_colorVBO);
	}

	if(m_pNNGrid) delete m_pNNGrid;
	m_vNeighs.clear();

	if(m_pBoundary) delete m_pBoundary;

	int num_solid = (int)m_vSolids.size();
	for(int i = 0; i < num_solid; ++i){
		delete m_vSolids[i];
	}
	m_vSolids.clear();	

	m_uNumParticles = 0;
	m_uMaxParticles = 0;
}




/*!
 * �p�[�e�B�N���f�[�^�̎擾
 * @return �p�[�e�B�N���f�[�^�̃f�o�C�X�������|�C���^
 */
RXREAL* rxDDSPH::GetParticle(void)
{
	return m_hPos;
}

/*!
 * �p�[�e�B�N���f�[�^�̎擾
 * @return �p�[�e�B�N���f�[�^�̃f�o�C�X�������|�C���^
 */
RXREAL* rxDDSPH::GetParticleDevice(void)
{
	return 0;
}

/*!
 * �J���[�l�pVBO�̕ҏW
 * @param[in] type �F�̃x�[�X�Ƃ��镨���l
 */
void rxDDSPH::SetColorVBO(int type)
{
	// MRK:SetColorVBO

	switch(type){
	case RX_DENSITY:
		SetColorVBOFromArray(m_hDens, 2, false, 2.0*m_fInitMaxDens);
		break;

	case RX_PRESSURE:
		SetColorVBOFromArray(m_hPres, 2);
		break;

	case RX_TEST:
		SetColorVBOFromArray(m_hTmp, 1, false, m_fTmpMax);
		break;
		
	case RX_SURFACE:	// �\�ʃp�[�e�B�N�������F��ς��ĕ`��
		if(m_bUseOpenGL){
			// �J���[�o�b�t�@�ɒl��ݒ�
			glBindBufferARB(GL_ARRAY_BUFFER, m_colorVBO);
			RXREAL *data = (RXREAL*)glMapBufferARB(GL_ARRAY_BUFFER, GL_WRITE_ONLY);
			RXREAL *ptr = data;
			for(uint i = 0; i < m_uNumParticles; ++i){
				uint s = m_hSurf[i];

				*ptr++ = (RXREAL)s;
				*ptr++ = 0.0f;
				*ptr++ = 0.0f;
				*ptr++ = 1.0f;
			}
			glUnmapBufferARB(GL_ARRAY_BUFFER);
		}
		break;

	case RX_RAMP:
		if(m_bUseOpenGL){
			// �J���[�o�b�t�@�ɒl��ݒ�
			glBindBufferARB(GL_ARRAY_BUFFER, m_colorVBO);
			RXREAL *data = (RXREAL*)glMapBufferARB(GL_ARRAY_BUFFER, GL_WRITE_ONLY);
			RXREAL *ptr = data;
			for(uint i = 0; i < m_uNumParticles; ++i){
				RXREAL t = i/(RXREAL)m_uNumParticles;
#if 0
				*ptr++ = rand()/(RXREAL)RAND_MAX;
				*ptr++ = rand()/(RXREAL)RAND_MAX;
				*ptr++ = rand()/(RXREAL)RAND_MAX;
#else
				RX_COLOR_RAMP(t, ptr);
				ptr += 3;
#endif
				*ptr++ = 1.0f;
			}
			glUnmapBufferARB(GL_ARRAY_BUFFER);
		}
		break;

	default:
		break;
	}

}


/*!
 * SPH��1�X�e�b�v�i�߂�
 * @param[in] dt ���ԃX�e�b�v��
 * @param[in] step ���݂̃X�e�b�v��
 * @retval ture  �v�Z����
 * @retval false �ő�X�e�b�v���𒴂��Ă��܂�
 */
bool rxDDSPH::Update(RXREAL dt, int step)
{
	// MARK:Update : ���ԍX�V
	assert(m_bInitialized);

	Vec3 v_half, acc_f;

	static bool init = false;

	// �����z�u���珉�����x���v�Z
	if(m_bCalDens){
		// �ߖT���q�T��
		SetParticlesToCell();

		calDoubleDensity(m_hPos, m_fInitDens, m_fInitDensNear, m_fInitMaxDens);
		if(m_fInitMaxDens > RX_FEQ_EPS){
			m_bCalDens = false;
		}
	}


	for(uint j = 0; j < m_solverIterations; ++j){
		// �ߖT���q�T��
		SetParticlesToCell();

		// �O�͍�(�d�́C�S���Ȃ�)
		calExternalForceToVel(dt);

		// �\���ʒu�X�V(prediction-relaxation scheme)
		for(uint i = 0; i < m_uNumParticles; ++i){
			for(uint j = 0; j < DIM; ++j){
				m_hPredictPos[DIM*i+j] = m_hPos[DIM*i+j]+dt*m_hVel[DIM*i+j];
			}
		}

		// �\���ʒu�ł̖��x�v�Z
		calDoubleDensity(m_hPredictPos);

		// ���q�ԃo�l�̒ǉ��ƍ폜
		adjustSprings(dt);

		// ���q�ԃo�l�ɂ��Elasticity,Plasticity
		applySpringDisplacements(dt);

		// �񈳏k��(double density relaxation)
		calDoubleDensityRelaxation(m_hPredictPos, dt);

		// �ʒu�C���x�X�V
		for(uint i = 0; i < m_uNumParticles; ++i){
			Vec3 x0, x1, v;
			for(int k = 0; k < 3; ++k){
				int idx = DIM*i+k;
				x0[k] = m_hPos[idx];
				x1[k] = m_hPredictPos[idx];
				v[k] = m_hVel[idx];
			}

			// ���E�Ƃ̏Փ˔���
			calCollisionSolid(x0, x1, v, dt);

			// �V�������x�ƈʒu�ōX�V
			for(int k = 0; k < 3; ++k){
				int idx = DIM*i+k;
				m_hPredictPos[idx] = x1[k];
				m_hVel[idx] = (m_hPredictPos[idx]-m_hPos[idx])/dt;
				m_hPos[idx] = m_hPredictPos[idx];
			}
		}

		RXTIMER("sph");

		init = false;
	}

	if(m_bCalNormal){
		CalNormal();
	}


	SetArrayVBO(RX_POSITION, m_hPos, 0, m_uNumParticles);

	SetColorVBO(m_iColorType);

	RXTIMER("color(vbo)");

	return true;
}





//-----------------------------------------------------------------------------
// Double Density
//-----------------------------------------------------------------------------
/*!
 * ���x���v�Z(�������x�v�Z�p)
 * @param[in] ppos �p�[�e�B�N�����W
 * @param[out] avg_dens ���ϖ��x
 * @param[out] avg_dens_near ���ϖ��x(near)
 * @param[out] max_dens �ő喧�x
 */
void rxDDSPH::calDoubleDensity(RXREAL *ppos, double &avg_dens, double &avg_dens_near, double &max_dens)
{
	avg_dens = 0.0;
	avg_dens_near = 0.0;
	max_dens = 0.0;
	double min_dens = RX_FEQ_INF;

	RXREAL h = m_fEffectiveRadius;
	RXREAL r0 = m_fInitDens;

	RXREAL k = m_fK;
	RXREAL knear = m_fKnear;

	for(uint i = 0; i < m_uNumParticles; ++i){
		Vec3 pos0 = Vec3(ppos[DIM*i+0], ppos[DIM*i+1], ppos[DIM*i+2]);

		m_hDens[2*i+0] = 0.0;	// Density
		m_hDens[2*i+1] = 0.0;	// Near Density

		// �ߖT���q
		for(vector<rxNeigh>::iterator itr = m_vNeighs[i].begin() ; itr != m_vNeighs[i].end(); ++itr){
			int j = itr->Idx;
			if(j < 0 || i == j) continue;

			Vec3 pos1 = Vec3(ppos[DIM*j+0], ppos[DIM*j+1], ppos[DIM*j+2]);
			Vec3 rij = pos1-pos0;

			RXREAL r = norm(rij);//sqrt(itr->Dist2);
			RXREAL q = r/h;

			if(q <= 1.0 && q >= 0.0){
				// density
				m_hDens[2*i]   += m_fMass*(1.0-q)*(1.0-q);

				// near density
				m_hDens[2*i+1] += m_fMass*(1.0-q)*(1.0-q)*(1.0-q);
			}
		}

		// ���ϖ��x�̌v�Z
		avg_dens += m_hDens[2*i];
		avg_dens_near += m_hDens[2*i+1];

		if(m_hDens[2*i] > max_dens) max_dens = m_hDens[2*i];
		if(m_hDens[2*i] < min_dens) min_dens = m_hDens[2*i];
	}

	if(m_uNumParticles){
		avg_dens /= (double)m_uNumParticles;
		avg_dens_near /= (double)m_uNumParticles;
	}

	RXCOUT << "minimum density : " << min_dens << endl;
	RXCOUT << "maximum density : " << max_dens << endl;
}

/*!
 * ���x�ƈ��͂��v�Z
 * @param[in] ppos �p�[�e�B�N�����W
 */
void rxDDSPH::calDoubleDensity(RXREAL *ppos)
{
	// MRK:calDensity
	RXREAL h = m_fEffectiveRadius;
	RXREAL r0 = m_fInitDens;

	RXREAL k = m_fK;
	RXREAL knear = m_fKnear;

	for(uint i = 0; i < m_uNumParticles; ++i){
		Vec3 pos0 = Vec3(ppos[DIM*i+0], ppos[DIM*i+1], ppos[DIM*i+2]);

		m_hDens[2*i+0] = 0.0;	// Density
		m_hDens[2*i+1] = 0.0;	// Near Density

		// �ߖT���q
		for(vector<rxNeigh>::iterator itr = m_vNeighs[i].begin() ; itr != m_vNeighs[i].end(); ++itr){
			int j = itr->Idx;
			if(j < 0 || i == j) continue;

			Vec3 pos1 = Vec3(ppos[DIM*j+0], ppos[DIM*j+1], ppos[DIM*j+2]);
			Vec3 rij = pos1-pos0;

			RXREAL r = norm(rij);//sqrt(itr->Dist2);
			RXREAL q = r/h;

			if(q <= 1.0 && q >= 0.0){
				// density
				m_hDens[2*i]   += m_fMass*(1.0-q)*(1.0-q);

				// near density
				m_hDens[2*i+1] += m_fMass*(1.0-q)*(1.0-q)*(1.0-q);
			}
		}

		m_hPres[2*i]   = k*(m_hDens[2*i]-r0);	// pressure
		m_hPres[2*i+1] = knear*m_hDens[2*i+1];	// near-pressure
	}
}

/*!
 * double density relaxation�ɂ��񈳏k�������̓K�p
 * @param[in] ppos �p�[�e�B�N�����W
 * @param[in] dt �^�C���X�e�b�v��
 */
void rxDDSPH::calDoubleDensityRelaxation(RXREAL *ppos, RXREAL dt)
{
	RXREAL h = m_fEffectiveRadius;

	Vec3 avg_dx = 0.0;
	for(uint i = 0; i < m_uNumParticles; ++i){
		Vec3 pos0 = Vec3(ppos[DIM*i+0], ppos[DIM*i+1], ppos[DIM*i+2]);

		RXREAL prsi  = m_hPres[2*i];
		RXREAL prsni = m_hPres[2*i+1];

		Vec3 dx(0.0);
		for(vector<rxNeigh>::iterator itr = m_vNeighs[i].begin() ; itr != m_vNeighs[i].end(); ++itr){
			int j = itr->Idx;
			if(j < 0 || i == j) continue;

			Vec3 pos1 = Vec3(ppos[DIM*j+0], ppos[DIM*j+1], ppos[DIM*j+2]);
			Vec3 rij = pos1-pos0;

			RXREAL r = norm(rij);//sqrt(itr->Dist2);
			RXREAL q = r/h;

			if(q <= 1.0 && q > 0.0){
				rij = Unit(rij);

				RXREAL prsj  = m_hPres[2*j];
				RXREAL prsnj = m_hPres[2*j+1];

				RXREAL prs, prsn;
				prs  = 0.5*(prsi+prsj);
				prsn = 0.5*(prsni+prsnj);

				Vec3 D = m_fMass*(0.5*prs*(1.0-q)+0.5*prsn*(1.0-q)*(1.0-q))*rij/m_hDens[2*j];

				dx -= 0.5*dt*dt*D;
				//da += 0.5*D*m_params.Mass;
			}
		}

		m_hPredictPos[DIM*i]   += dx[0];
		m_hPredictPos[DIM*i+1] += dx[1];
	}

}

/*!
 * �O�͍��C�S�����̌v�Z
 * @param[in] dt �^�C���X�e�b�v��
 */
void rxDDSPH::calExternalForceToVel(RXREAL dt)
{
	RXREAL h  = m_fEffectiveRadius;
	RXREAL r0 = m_fInitDens;

	RXREAL c  = m_fViscC;
	RXREAL B  = m_fViscBeta;

	// ���x��̕ۑ�
	for(uint i = 0; i < m_uNumParticles; ++i){
		for(uint j = 0; j < DIM; ++j){
			int idx = DIM*i+j;
			m_hVelOld[idx] = m_hVel[idx];
		}
	}

	for(uint i = 0; i < m_uNumParticles; ++i){
		Vec3 pos0, vel0;
		pos0 = Vec3(m_hPos[DIM*i+0], m_hPos[DIM*i+1], m_hPos[DIM*i+2]);
		vel0 = Vec3(m_hVelOld[DIM*i+0], m_hVelOld[DIM*i+1], m_hVelOld[DIM*i+2]);

		// �S�����̌v�Z
		Vec3 rij, vji;
		Vec3 pfrc(0.0), vfrc(0.0);
		for(vector<rxNeigh>::iterator itr = m_vNeighs[i].begin() ; itr != m_vNeighs[i].end(); ++itr){
			int j = itr->Idx;
			if(j < 0) continue;

			if((int)i < j){
				RXREAL r = sqrt(itr->Dist2);
				RXREAL q = r/h;

				if(q < 1.0 && q > RX_FEQ_EPS){
					Vec3 pos1, vel1;
					pos1 = Vec3(m_hPos[DIM*j+0], m_hPos[DIM*j+1], m_hPos[DIM*j+2]);
					vel1 = Vec3(m_hVelOld[DIM*j+0], m_hVelOld[DIM*j+1], m_hVelOld[DIM*j+2]);

					rij = Unit(pos1-pos0);
					RXREAL u = (RXREAL)dot(vel0-vel1, rij);

					if(u > 0.0){
						// �S��
						Vec3 I = dt*(1.0-q)*(c*u+B*u*u)*rij;
						m_hVel[DIM*i+0] -= 0.5*I[0];
						m_hVel[DIM*i+1] -= 0.5*I[1];
						m_hVel[DIM*i+2] -= 0.5*I[2];
						m_hVel[DIM*j+0] += 0.5*I[0];
						m_hVel[DIM*j+1] += 0.5*I[1];
						m_hVel[DIM*j+2] += 0.5*I[2];
					}
				}
			}
		}

		// �O�͍��̌v�Z
		Vec3 g;
		g[0] = m_v3Gravity[0];
		g[1] = m_v3Gravity[1];
		g[2] = m_v3Gravity[2];
		Vec3 v = dt*g;

		m_hVel[DIM*i+0] += v[0];
		m_hVel[DIM*i+1] += v[1];
		m_hVel[DIM*i+2] += v[2];
	}
}


/*!
 * ���q�ԃo�l�̒ǉ��ƍ폜
 * @param[in] dt �^�C���X�e�b�v��
 */
void rxDDSPH::adjustSprings(RXREAL dt)
{
	RXREAL h = m_fEffectiveRadius;
	RXREAL r0 = m_fInitDens;

	RXREAL gamma = m_fPlasGamma; // fraction of L : [0, 0.2]
	RXREAL alpha = m_fPlasAlpha; // plasticity constant

	map<int, rxSPHSpring>::iterator itr;
	for(itr = m_mapSpring.begin(); itr != m_mapSpring.end(); ++itr){
		itr->second.enable = false;
	}

	// ���q�ԃo�l�̒ǉ��ƒ���
	for(uint i = 0; i < m_uNumParticles; ++i){
		for(vector<rxNeigh>::iterator itr = m_vNeighs[i].begin() ; itr != m_vNeighs[i].end(); ++itr){
			int j = itr->Idx;
			if(j < 0) continue;

			if((int)i < j){
				RXREAL r = sqrt(itr->Dist2);
				RXREAL q = r/h;

				int idx = j*m_uNumParticles+i;
				if(!m_mapSpring[idx].enable){
					// ���q�ԃo�l��V���ɒǉ�
					m_mapSpring[idx].pi = i;
					m_mapSpring[idx].pj = j;
					m_mapSpring[idx].L = r;
				}

				m_mapSpring[idx].enable = true;
				m_mapSpring[idx].r = r;

				RXREAL d = gamma*m_mapSpring[idx].L;

				if(r > m_mapSpring[idx].L+d){
					m_mapSpring[idx].L += dt*alpha*(r-m_mapSpring[idx].L-d);
				}
				else if(r < m_mapSpring[idx].L-d){
					m_mapSpring[idx].L -= dt*alpha*(m_mapSpring[idx].L-d-r);
				}
			}
		}
	}

	// ���q�ԃo�l�̍폜
	map<int, rxSPHSpring>::iterator iter = m_mapSpring.begin();
	while(iter != m_mapSpring.end()){
		if(!iter->second.enable){
			m_mapSpring.erase(iter++);
		}
		else{
			++iter;
		}
	}
}


/*!
 * Elasticity,Plasticity
 * @param[in] dt �^�C���X�e�b�v��
 */
void rxDDSPH::applySpringDisplacements(RXREAL dt)
{
	RXREAL h = m_fEffectiveRadius;
	RXREAL r0 = m_fInitDens;

	RXREAL ks = m_fElasKspring; // �o�l�萔

	map<int, rxSPHSpring>::iterator itr;
	for(itr = m_mapSpring.begin(); itr != m_mapSpring.end(); ++itr){
		if(itr->second.enable){
			int i = itr->second.pi;
			int j = itr->second.pj;

			Vec3 pos0, pos1;
			pos0 = Vec3(m_hPos[DIM*i+0], m_hPos[DIM*i+1], m_hPos[DIM*i+2]);
			pos1 = Vec3(m_hPos[DIM*j+0], m_hPos[DIM*j+1], m_hPos[DIM*j+2]);

			Vec3 rij = Unit(pos1-pos0);
			Vec3 D = dt*dt*ks*(1.0-itr->second.L/h)*(itr->second.L-itr->second.r)*rij;

			pos0 -= 0.5*D;
			pos1 += 0.5*D;

			for(int k = 0; k < DIM; ++k){
				m_hPos[DIM*i+k] = pos0[k];
				m_hPos[DIM*j+k] = pos1[k];
			}
		}
	}
}


/*!
 * �ő̃I�u�W�F�N�g�Ƃ̏Փ˔���C�Փˉ���
 * @param[in] pos0 �O�X�e�b�v�̈ʒu
 * @param[inout] pos1 �V�����ʒu
 * @param[inout] vel ���x
 * @param[in] dt �^�C���X�e�b�v��
 * @return �Փ˃I�u�W�F�N�g�̐�
 */
int rxDDSPH::calCollisionSolid(Vec3 &pos0, Vec3 &pos1, Vec3 &vel, RXREAL dt)
{
	int c = 0;
	rxCollisionInfo coli;

	// �ő̃I�u�W�F�N�g�Ƃ̏Փˏ���
	for(vector<rxSolid*>::iterator i = m_vSolids.begin(); i != m_vSolids.end(); ++i){
		if((*i)->GetDistanceR(pos1, m_fParticleRadius, coli)){
			RXREAL res = m_fDamping;
			res = (res > 0) ? (res*fabs(coli.Penetration())/(dt*norm(vel))) : 0.0f;
			vel -= (1+res)*dot(vel, coli.Normal())*coli.Normal();
			pos1 = coli.Contact();
		}
	}

	// �V�~�����[�V������ԋ��E�Ƃ̏Փˏ���
	if(m_pBoundary->GetDistanceR(pos1, m_fParticleRadius, coli)){
		RXREAL res = m_fDamping;
		res = (res > 0) ? (res*fabs(coli.Penetration())/(dt*norm(vel))) : 0.0f;
		vel -= (1+res)*dot(vel, coli.Normal())*coli.Normal();
		pos1 = coli.Contact()+coli.Normal()*RXFunc::Rand(0.0, 0.001);
	}

	return c;
}



//-----------------------------------------------------------------------------
// �\�ʃp�[�e�B�N���̌��o
//-----------------------------------------------------------------------------
/*!
 * �\�ʃp�[�e�B�N�����o
 *  - B. Solenthaler, Y. Zhang and R. Pajarola, "Efficient Refinement of Dynamic Point Data",
 *    Proceedings Eurographics/IEEE VGTC Symposium on Point-Based Graphics, 2007.
 *  - 3.1 Surface Particle Detection �� ��(2) 
 *  - d_i,cm ��臒l�ȏ�ŕ\�ʃp�[�e�B�N���Ɣ���Ə�����Ă��邪�C
 *    d_i,cm �͋ߖT�p�[�e�B�N���̏d�S�ւ̃x�N�g���Ŏ��ۂɂ͂��̃x�N�g���̒��� |d_i,cm| ���g���Ĕ���
 */
void rxDDSPH::DetectSurfaceParticles(void)
{
	RXREAL h = m_fEffectiveRadius;
	RXREAL r = m_fParticleRadius;

	//m_hPos = GetArrayVBO(RX_POSITION);

	// �ߖT���q�T��
	SetParticlesToCell();

	for(uint i = 0; i < m_uNumParticles; ++i){
		double d = CalDistToNormalizedMassCenter(i);
		int nn_num = (int)m_vNeighs[i].size();

		// �f�o�b�O�p
		//m_hTmp[i] = (RXREAL)d;
		//m_fTmpMax = r;

		//m_hSurf[i] = 0;
		if(nn_num <= 3){	// �ߖT�p�[�e�B�N������3�ȉ��Ȃ�Ε\��
			m_hSurf[i] = 1;
		}
		else{				// 3���傫���ꍇ�͋ߖT�d�S�܂ł̋����Ŕ��f
			if(m_hSurf[i]){
				// �O�X�e�b�v�ŕ\�ʂ������珬����臒l�Ŕ��f
				if(d < g_fSurfThr[0]*r) m_hSurf[i] = 0;
			}
			else{
				if(d > g_fSurfThr[1]*r || d < 0.0) m_hSurf[i] = 1;
			}
		}

		// ���x���g���Ĕ��f����ꍇ
		//if(m_hDens[2*i] < 0.7*m_fInitMaxDens){
		//	m_hSurf[i] = 1;
		//}
	}
}

/*!
 * �ߖT�p�[�e�B�N���̐��K���d�S�܂ł̋������v�Z
 * @param[in] i �p�[�e�B�N���C���f�b�N�X
 */
double rxDDSPH::CalDistToNormalizedMassCenter(const int i)
{
	Vec3 sum_pos(0.0);
	double sum_mass = 0.0;

	Vec3 pos0 = Vec3(m_hPos[DIM*i+0], m_hPos[DIM*i+1], m_hPos[DIM*i+2]);

	for(vector<rxNeigh>::iterator itr = m_vNeighs[i].begin() ; itr != m_vNeighs[i].end(); ++itr){
		int j = itr->Idx;
		if(j < 0 && i == j) continue;

		//RXREAL r = sqrt(itr->Dist2);
		Vec3 pos1 = Vec3(m_hPos[DIM*j+0], m_hPos[DIM*j+1], m_hPos[DIM*j+2]);

		sum_pos  += (pos0-pos1)*m_fMass;
		sum_mass += m_fMass;
	}

	double dis = DBL_MIN;
	if(sum_mass != 0.0){
		sum_pos /= sum_mass;
		dis     = norm(sum_pos);
	}

	return dis;
}

/*!
 * �\�ʃp�[�e�B�N�����̎擾
 *  �p�[�e�B�N�����Ɠ����傫����uint�z��ŕ\�ʃp�[�e�B�N���Ȃ��1, �����łȂ����0���i�[����Ă���
 */
uint* rxDDSPH::GetArraySurf(void)
{
	return m_hSurf;
}

/*!
 * �w�肵�����W�̋ߖT�̕\�ʃp�[�e�B�N�������擾����
 * @param[in] pos �T�����S���W
 * @param[in] h �T�����a(0�ȉ��̒l��ݒ肵����L�����a��p����)
 * @param[out] sp �\�ʃp�[�e�B�N��
 * @return ���������p�[�e�B�N����
 */
int rxDDSPH::GetSurfaceParticles(const Vec3 pos, RXREAL h, vector<rxSurfaceParticle> &sps)
{
	int num = 0;
	sps.clear();

	DetectSurfaceParticles();

	if(pos[0] < m_v3EnvMin[0]) return 0;
	if(pos[0] > m_v3EnvMax[0]) return 0;
	if(pos[1] < m_v3EnvMin[1]) return 0;
	if(pos[1] > m_v3EnvMax[1]) return 0;
	if(pos[2] < m_v3EnvMin[2]) return 0;
	if(pos[2] > m_v3EnvMax[2]) return 0;

	if(h <= 0.0) h = m_fEffectiveRadius;

	vector<rxNeigh> ne;
	m_pNNGrid->GetNN(pos, m_hPos, m_uNumParticles, ne, h);

	// �ߖT���q
	for(vector<rxNeigh>::iterator itr = ne.begin(); itr != ne.end(); ++itr){
		int j = itr->Idx;
		if(j < 0) continue;

		int surf = m_hSurf[j];

		if(surf){
			rxSurfaceParticle sp;
			sp.pos = Vec3(m_hPos[DIM*j+0], m_hPos[DIM*j+1], m_hPos[DIM*j+2]);
			sp.vel = Vec3(m_hVel[DIM*j+0], m_hVel[DIM*j+1], m_hVel[DIM*j+2]);
			sp.nrm = Vec3(m_hNrm[DIM*j+0], m_hNrm[DIM*j+1], m_hNrm[DIM*j+2]);
			sp.idx = j;
			sp.d = sqrt(itr->Dist2);
			sps.push_back(sp);
			num++;
		}
	}

	return num;
}



//-----------------------------------------------------------------------------
// �@���v�Z
//-----------------------------------------------------------------------------
/*!
 * �@�����v�Z
 *  - �z�쏃��C"���q�@�ɂ�闬�̃V�~�����[�V�����̍��i�������_�����O"�C�É���w���Ƙ_���C2007.
 *  - p.24, 3.4.3 �@���̏C�� �� ��(3.8) 
 *  - ��(3.8) �ł� ��w(v-r)(v-r) �ƂȂ��Ă��邪 ��w(v-r) �̊ԈႢ
 */
void rxDDSPH::CalNormal(void)
{
	RXREAL h = m_fEffectiveRadius;

	// �ߖT���q�T��
	SetParticlesToCell();

	//CalNormalFromDensity();

	for(uint i = 0; i < m_uNumParticles; ++i){
		Vec3 pos0 = Vec3(m_hPos[DIM*i+0], m_hPos[DIM*i+1], m_hPos[DIM*i+2]);
		Vec3 nrm(0.0);

		// �ߖT���q
		for(vector<rxNeigh>::iterator itr = m_vNeighs[i].begin() ; itr != m_vNeighs[i].end(); ++itr){
			int j = itr->Idx;
			if(j < 0 || i == j) continue;

			Vec3 pos1 = Vec3(m_hPos[DIM*j+0], m_hPos[DIM*j+1], m_hPos[DIM*j+2]);
			Vec3 rij = pos0-pos1;

			RXREAL rr = sqrt(itr->Dist2)/h;

			RXREAL w = (rr < 1.0 ? 1.0/rr-1.0 : 0.0);

			nrm += w*rij;
		}

		//normalize(nrm);
		//nrm *= -1;

		m_hNrm[DIM*i+0] = nrm[0];
		m_hNrm[DIM*i+1] = nrm[1];
	}
}

/*!
 * ���x���z����@�����v�Z
 */
void rxDDSPH::CalNormalFromDensity(void)
{
	RXREAL h = m_fEffectiveRadius;

	for(uint i = 0; i < m_uNumParticles; ++i){
		Vec3 pos0;
		pos0[0] = m_hPos[DIM*i+0];
		pos0[1] = m_hPos[DIM*i+1];
		pos0[2] = m_hPos[DIM*i+2];

		Vec3 nrm(0.0);

		// �ߖT���q
		for(vector<rxNeigh>::iterator itr = m_vNeighs[i].begin() ; itr != m_vNeighs[i].end(); ++itr){
			int j = itr->Idx;
			if(j < 0 || i == j) continue;

			Vec3 pos1;
			pos1[0] = m_hPos[DIM*j+0];
			pos1[1] = m_hPos[DIM*j+1];
			pos1[2] = m_hPos[DIM*j+2];

			Vec3 rij = pos0-pos1;

			RXREAL r = sqrt(itr->Dist2);
			RXREAL q = h*h-r*r;

			nrm += (m_fMass/m_hDens[j])*m_fGWpoly6*q*q*rij;
		}

		//normalize(nrm);
		//nrm *= -1;

		m_hNrm[DIM*i+0] = nrm[0];
		m_hNrm[DIM*i+1] = nrm[1];
	}
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
void rxDDSPH::GetNearestNeighbors(int idx, RXREAL *prts, vector<rxNeigh> &neighs, RXREAL h)
{
	if(idx < 0 || idx >= (int)m_uNumParticles) return;

	Vec3 pos;
	pos[0] = prts[DIM*idx+0];
	pos[1] = prts[DIM*idx+1];
	pos[2] = prts[DIM*idx+2];

	if(h < 0.0) h = m_fEffectiveRadius;

	m_pNNGrid->GetNN(pos, prts, m_uNumParticles, neighs, h);
	//m_pNNGrid->GetNN_Direct(pos, prts, m_uNumParticles, neighs, h);	// �O���b�h���g��Ȃ���������
}

/*!
 * �ߖT���q�T��
 * @param[in] idx �T�����S�p�[�e�B�N���C���f�b�N�X
 * @param[out] neighs �T�����ʊi�[����ߖT���R���e�i
 * @param[in] h �L�����a
 */
void rxDDSPH::GetNearestNeighbors(Vec3 pos, vector<rxNeigh> &neighs, RXREAL h)
{
	if(h < 0.0) h = m_fEffectiveRadius;

	m_pNNGrid->GetNN(pos, m_hPos, m_uNumParticles, neighs, h);
	//m_pNNGrid->GetNN_Direct(pos, prts, m_uNumParticles, neighs, h);	// �O���b�h���g��Ȃ���������
}


/*!
 * �S�p�[�e�B�N���𕪊��Z���Ɋi�[
 */
void rxDDSPH::SetParticlesToCell(RXREAL *prts, int n, RXREAL h)
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
void rxDDSPH::SetParticlesToCell(void)
{
	SetParticlesToCell(m_hPos, m_uNumParticles, m_fEffectiveRadius);
}

/*!
 * �����Z���Ɋi�[���ꂽ�|���S�������擾
 * @param[in] gi,gj,gk �Ώە����Z��
 * @param[out] polys �|���S��
 * @return �i�[�|���S����
 */
int rxDDSPH::GetPolygonsInCell(int gi, int gj, int gk, vector<int> &polys)
{
	return m_pNNGrid->GetPolygonsInCell(gi, gj, gk, polys);
}
/*!
 * �����Z�����̃|���S���̗L���𒲂ׂ�
 * @param[in] gi,gj,gk �Ώە����Z��
 * @return �|���S�����i�[����Ă����true
 */
bool rxDDSPH::IsPolygonsInCell(int gi, int gj, int gk)
{
	return m_pNNGrid->IsPolygonsInCell(gi, gj, gk);
}

/*!
 * �|���S���𕪊��Z���Ɋi�[
 */
void rxDDSPH::SetPolygonsToCell(void)
{
	m_pNNGrid->SetPolygonsToCell(m_hVrts, m_iNumVrts, m_hTris, m_iNumTris);
}


//-----------------------------------------------------------------------------
// OpenGL�`��
//-----------------------------------------------------------------------------
/*!
 * �T���p�Z���̕`��
 * @param[in] i,j,k �O���b�h��̃C���f�b�N�X
 */
void rxDDSPH::DrawCell(int i, int j, int k)
{
	if(m_pNNGrid) m_pNNGrid->DrawCell(i, j, k);
}

/*!
 * �T���p�O���b�h�̕`��
 * @param[in] col �p�[�e�B�N�����܂܂��Z���̐F
 * @param[in] col2 �|���S�����܂܂��Z���̐F
 * @param[in] sel �����_���ɑI�����ꂽ�Z���̂ݕ`��(1�ŐV�����Z����I���C2�ł��łɑI������Ă���Z����`��C0�ł��ׂẴZ����`��)
 */
void rxDDSPH::DrawCells(Vec3 col, Vec3 col2, int sel)
{
	if(m_pNNGrid) m_pNNGrid->DrawCells(col, col2, sel, m_hPos);
}


/*!
 * �ő̏�Q���̕`��
 */
void rxDDSPH::DrawObstacles(void)
{
	for(vector<rxSolid*>::iterator i = m_vSolids.begin(); i != m_vSolids.end(); ++i){
		(*i)->Draw();
	}
}



//-----------------------------------------------------------------------------
// �ő̃I�u�W�F�N�g
//-----------------------------------------------------------------------------
/*!
 * �O�p�`�|���S���ɂ���Q��
 * @param[in] vrts ���_
 * @param[in] tris ���b�V��
 */
void rxDDSPH::SetPolygonObstacle(const vector<Vec3> &vrts, const vector<Vec3> &nrms, const vector< vector<int> > &tris, Vec3 vel)
{
	int vn = (int)vrts.size();
	int n = (int)tris.size();

	if(m_hVrts) delete [] m_hVrts;
	if(m_hTris) delete [] m_hTris;
	m_hVrts = new RXREAL[vn*3];
	m_hTris = new int[n*3];

	for(int i = 0; i < vn; ++i){
		for(int j = 0; j < 3; ++j){
			m_hVrts[3*i+j] = vrts[i][j];
		}
	}

	for(int i = 0; i < n; ++i){
		for(int j = 0; j < 3; ++j){
			m_hTris[3*i+j] = tris[i][j];
		}
	}

	m_iNumVrts = vn;
	m_iNumTris = n;
	RXCOUT << "the number of triangles : " << m_iNumTris << endl;

	m_pNNGrid->SetPolygonsToCell(m_hVrts, m_iNumVrts, m_hTris, m_iNumTris);

	//setPolysToCell(m_hVrts, m_iNumVrts, m_hTris, m_iNumTris);
}

/*!
 * �{�b�N�X�^��Q��
 * @param[in] cen �{�b�N�X���S���W
 * @param[in] ext �{�b�N�X�̑傫��(�ӂ̒�����1/2)
 * @param[in] ang �{�b�N�X�̊p�x(�I�C���[�p)
 * @param[in] flg �L��/�����t���O
 */
void rxDDSPH::SetBoxObstacle(Vec3 cen, Vec3 ext, Vec3 ang, Vec3 vel, int flg)
{
	rxSolidBox *box = new rxSolidBox(cen-ext, cen+ext, 1);
	double m[16];
	EulerToMatrix(m, ang[0], ang[1], ang[2]);
	box->SetMatrix(m);

	m_vSolids.push_back(box);

}

/*!
 * ���^��Q��
 * @param[in] cen ���̒��S���W
 * @param[in] rad ���̂̔��a
 * @param[in] flg �L��/�����t���O
 */
void rxDDSPH::SetSphereObstacle(Vec3 cen, double rad, Vec3 vel, int flg)
{
	rxSolidSphere *sphere = new rxSolidSphere(cen, rad, 1);
	m_vSolids.push_back(sphere);
}

/*!
 * ���^��Q���𓮂���
 * @param[in] b ���̔ԍ�
 * @param[in] disp �ړ���
 */
void rxDDSPH::MoveSphereObstacle(int b, Vec3 disp)
{
	if(m_vSolids.empty()) return;
	m_vSolids[b]->SetPosition(m_vSolids[b]->GetPosition()+disp);
}

/*!
 * ���^��Q���̈ʒu���擾
 * @param[in] b ���̔ԍ�
 */
Vec3 rxDDSPH::GetSphereObstaclePos(int b)
{
	if(m_vSolids.empty()) return Vec3(0.0);
	return m_vSolids[b]->GetPosition();
}

/*!
 * VBO����z�X�g�������փf�[�^��]���C�擾
 * @param[in] type �f�[�^�̎��
 * @return �z�X�g��������̃f�[�^
 */
RXREAL* rxDDSPH::GetArrayVBO(rxParticleArray type, bool d2h, int num)
{
	assert(m_bInitialized);
 
	if(num == -1) num = m_uNumParticles;
	RXREAL* hdata = 0;

	unsigned int vbo = 0;

	switch(type){
	default:
	case RX_POSITION:
		hdata = m_hPos;
		vbo = m_posVBO;
		break;

	case RX_VELOCITY:
		hdata = m_hVel;
		break;

	case RX_DENSITY:
		hdata = m_hDens;
		break;

	case RX_PRESSURE:
		hdata = m_hPres;
		break;

	case RX_NORMAL:
		hdata = m_hNrm;
		break;

	case RX_FORCE:
		hdata = m_hFrc;
		break;
	}

	return hdata;
}

/*!
 * �z�X�g����������VBO�������փf�[�^��]��
 * @param[in] type �f�[�^�̎��
 * @param[in] data �z�X�g��������̃f�[�^
 * @param[in] start �f�[�^�̊J�n�C���f�b�N�X
 * @param[in] count �ǉ���
 */
void rxDDSPH::SetArrayVBO(rxParticleArray type, const RXREAL* data, int start, int count)
{
	assert(m_bInitialized);
 
	switch(type){
	default:
	case RX_POSITION:
		{
			if(m_bUseOpenGL){
				glBindBuffer(GL_ARRAY_BUFFER, m_posVBO);
				glBufferSubData(GL_ARRAY_BUFFER, start*2*sizeof(RXREAL), count*2*sizeof(RXREAL), data);
				glBindBuffer(GL_ARRAY_BUFFER, 0);
			}
		}
		break;

	case RX_VELOCITY:
		//CuCopyArrayToDevice(m_dVel, data, start*DIM*sizeof(RXREAL), count*DIM*sizeof(RXREAL));
		break;

	case RX_NORMAL:
		//CuCopyArrayToDevice(m_dNrm, data, start*DIM*sizeof(RXREAL), count*DIM*sizeof(RXREAL));
		break;
	}	   
}


//-----------------------------------------------------------------------------
// �A�֐��l
//-----------------------------------------------------------------------------
double rxDDSPH::GetImplicit(double x, double y, double z)
{
	return CalColorField(x, y, z);
}

/*!
 * �p�[�e�B�N������O���b�h�̉A�֐��l���v�Z
 * @param[in] n �O���b�h��
 * @param[in] minp �O���b�h�̍ŏ����W
 * @param[in] d �O���b�h��
 * @param[out] hF �A�֐��l(nx�~ny�~nz�̔z��)
 */
void rxDDSPH::CalImplicitField(int n[3], Vec3 minp, Vec3 d, RXREAL *hF)
{
	int slice0 = n[0];
	int slice1 = n[0]*n[1];

	for(int k = 0; k < n[2]; ++k){
		for(int j = 0; j < n[1]; ++j){
			for(int i = 0; i < n[0]; ++i){
				int idx = k*slice1+j*slice0+i;
				Vec3 pos = minp+Vec3(i, j, k)*d;
				hF[idx] = GetImplicit(pos[0], pos[1], pos[2]);
			}
		}
	}
}

/*!
 * �p�[�e�B�N������O���b�h�̉A�֐��l���v�Z
 * @param[in] n �O���b�h��
 * @param[in] minp �O���b�h�̍ŏ����W
 * @param[in] d �O���b�h��
 * @param[out] hF �A�֐��l(nx�~ny�~nz�̔z��)
 */
void rxDDSPH::CalImplicitFieldDevice(int n[3], Vec3 minp, Vec3 d, RXREAL *dF)
{
	RXREAL *hF = new RXREAL[n[0]*n[1]*n[2]];

	CalImplicitField(n, minp, d, hF);
	//CuCopyArrayToDevice(dF, hF, 0, n[0]*n[1]*n[2]*sizeof(RXREAL));

	delete [] hF;
}


/*!
 * �J���[�t�B�[���h�l�v�Z
 *  - Poly6�J�[�l��
 * @param[in] x,y,z �v�Z�ʒu
 * @return �J���[�t�B�[���h�l
 */
double rxDDSPH::CalColorField(double x, double y, double z)
{
	// MRK:CalColorField
	RXREAL c = 0.0;
	Vec3 pos(x, y, z);

	if(pos[0] < m_v3EnvMin[0]) return c;
	if(pos[0] > m_v3EnvMax[0]) return c;
	if(pos[1] < m_v3EnvMin[1]) return c;
	if(pos[1] > m_v3EnvMax[1]) return c;
	if(pos[2] < m_v3EnvMin[2]) return c;
	if(pos[2] > m_v3EnvMax[2]) return c;

	RXREAL h = m_fEffectiveRadius;

	vector<rxNeigh> ne;
	m_pNNGrid->GetNN(pos, m_hPos, m_uNumParticles, ne, h);

	// �ߖT���q
	for(vector<rxNeigh>::iterator itr = ne.begin(); itr != ne.end(); ++itr){
		int j = itr->Idx;
		if(j < 0) continue;

		Vec3 pos1;
		pos1[0] = m_hPos[DIM*j+0];
		pos1[1] = m_hPos[DIM*j+1];
		pos1[2] = m_hPos[DIM*j+2];

		RXREAL r = sqrt(itr->Dist2);

		RXREAL q = h*h-r*r;
		c += m_fMass*m_fWpoly6*q*q*q;

		//RXREAL q = r/h;
		//if(q <= 1.0 && q >= 0.0){
		//	c += m_fMass*(1.0-q)*(1.0-q);
		//}
	}

	return c;
}


//-----------------------------------------------------------------------------
// MARK:�V�~���f�[�^�̏o��
//-----------------------------------------------------------------------------
/*!
 * �V�~�����[�V�����ݒ�(�p�[�e�B�N�����C�͈́C���x�C���ʂȂ�)
 * @param[in] fn �o�̓t�@�C����
 */
void rxDDSPH::OutputSetting(string fn)
{
	ofstream fout;
	fout.open(fn.c_str());
	if(!fout){
		RXCOUT << fn << " couldn't open." << endl;
		return;
	}

	fout << m_uNumParticles << endl;
	fout << m_v3EnvMin[0] << " " << m_v3EnvMin[1] << " " << m_v3EnvMin[2] << endl;
	fout << m_v3EnvMax[0] << " " << m_v3EnvMax[1] << " " << m_v3EnvMax[2] << endl;
	fout << m_fInitDens << endl;
	fout << m_fMass << endl;
	fout << m_iKernelParticles << endl;

	fout.close();
}


