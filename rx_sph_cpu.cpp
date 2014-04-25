/*!
  @file rx_sph_cpu.cpp
	
  @brief SPH�@(GPU)�̎���
 
  @author Makoto Fujisawa
  @date   2008-10,2011-06
*/
// FILE --rx_sph_cpu.cpp--

//-----------------------------------------------------------------------------
// �C���N���[�h�t�@�C��
//-----------------------------------------------------------------------------
#include "rx_sph.h"

#include "rx_cu_funcs.cuh"
#include "rx_cu_common.cuh"

#include "rx_pcube.h"




//-----------------------------------------------------------------------------
// rxSPH�N���X�̎���
//-----------------------------------------------------------------------------
/*!
 * �R���X�g���N�^
 * @param[in] use_opengl VBO�g�p�t���O
 */
rxSPH::rxSPH(bool use_opengl) : 
	rxParticleSystemBase(use_opengl), 
	m_hNrm(0), 
	m_hFrc(0),
	m_hDens(0), 
	m_hPres(0), 
	m_hSurf(0), 
	m_hUpPos(0), 
	m_hPosW(0), 
	m_hEigen(0), 
	m_hRMatrix(0), 
	m_hG(0), 
	m_hVrts(0), 
	m_hTris(0), 
	m_pBoundary(0)
{
	m_v3Gravity = Vec3(0.0, -9.82, 0.0);

	m_fViscosity = 0.01;
	m_fGasStiffness = 3.0;
	m_fBuoyancy = 0.0f;

	m_fDamping = 0.0;
	m_fRestitution = 0.0;

	m_bUseVorticity = false;
	m_bUseWaveletTurb = false;
	m_bGridVelocity = false;
	m_bCalNormal = false;
	m_bUpsampling = false;
	m_bSubParticle = false;
	m_bCalAnisotropic = false;

	// �ߖT�T���Z��
	m_pNNGrid = new rxNNGrid(DIM);
	m_pNNGridB = new rxNNGrid(DIM);	// ���E�p�[�e�B�N���p

	m_uNumParticles = 0;
	m_uNumBParticles = 0;
	m_iNumTris = 0;

	m_iColorType = RX_RAMP;
}

/*!
 * �f�X�g���N�^
 */
rxSPH::~rxSPH()
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
void rxSPH::Initialize(const rxSPHEnviroment &env)
{
	// MARK:Initialize
	RXCOUT << "[rxSPH::Initialize]" << endl;

	m_fRestDens        = env.dens;
	m_fMass            = env.mass;
	m_iKernelParticles = env.kernel_particles;

	RXREAL volume = m_iKernelParticles*m_fMass/m_fRestDens;

	m_fEffectiveRadius = pow(((3.0*volume)/(4.0*RX_PI)), 1.0/3.0);
	m_fParticleRadius = pow((RX_PI/(6.0*m_iKernelParticles)), 1.0/3.0)*m_fEffectiveRadius;
	//m_fParticleRadius = 0.5f*m_fEffectiveRadius;

	m_fViscosity = env.viscosity;
	m_fGasStiffness = env.gas_k;

	RXREAL h = m_fEffectiveRadius;
	RXREAL r = m_fParticleRadius;

	// �J�[�l���֐��̒萔
	m_fAw = KernelCoefPoly6(h, 3, 1);
	m_fAg = KernelCoefSpiky(h, 3, 2);
	m_fAl = KernelCoefVisc(h, 3, 3);

	// �J�[�l���֐�
	m_fpW  = KernelPoly6;
	m_fpGW = KernelSpikyG<Vec3>;
	m_fpLW = KernelViscL;

	// �������x�̌v�Z
	m_fRestDens = calRestDensity(h);

	RXCOUT << "particle : " << endl;
	RXCOUT << " n_max = " << env.max_particles << endl;
	RXCOUT << " h = " << m_fEffectiveRadius << endl;
	RXCOUT << " r = " << m_fParticleRadius << endl;
	RXCOUT << " dens = " << m_fRestDens << endl;
	RXCOUT << " mass = " << m_fMass << endl;
	RXCOUT << " kernel_particles = " << m_iKernelParticles << endl;
	RXCOUT << " volume = " << volume << endl;
	RXCOUT << " viscosity = " << m_fViscosity << endl;


	//
	// ���E�ݒ�
	//
	m_pBoundary = new rxSolidBox(env.boundary_cen-env.boundary_ext, env.boundary_cen+env.boundary_ext, -1);
	//m_pBoundary = new rxSolidSphere(Vec3(0.0, 0.1, 0.0), 0.25, -1);

#ifdef RX_USE_BOUNDARY_PARTICLE
	// ���E�p�[�e�B�N������
	m_uNumBParticles = m_pBoundary->GenerateParticlesOnSurf(0.85*m_fParticleRadius, &m_hPosB);
#endif

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

	m_uNumParticles = 0;

	Allocate(env.max_particles);
}

/*!
 * �������̊m��
 *  - �ő�p�[�e�B�N�����Ŋm��
 * @param[in] max_particles �ő�p�[�e�B�N����
 */
void rxSPH::Allocate(int max_particles)
{
	// MARK:Allocate
	assert(!m_bInitialized);

	//m_uNumParticles = max_particles;
	m_uMaxParticles = max_particles;
	
	unsigned int size  = m_uMaxParticles*DIM;
	unsigned int size1 = m_uMaxParticles;
	unsigned int mem_size  = sizeof(RXREAL)*size;
	unsigned int mem_size1 = sizeof(RXREAL)*size1;

	//
	// �������m��
	//
	m_hPos = new RXREAL[size];
	m_hVel = new RXREAL[size];
	m_hNrm = new RXREAL[size];
	m_hFrc = new RXREAL[size];
	memset(m_hPos, 0, mem_size);
	memset(m_hVel, 0, mem_size);
	memset(m_hNrm, 0, mem_size);
	memset(m_hFrc, 0, mem_size);

	m_hDens = new RXREAL[size1];
	m_hPres = new RXREAL[size1];
	memset(m_hDens, 0, mem_size1);
	memset(m_hPres, 0, mem_size1);

	m_hSurf = new uint[m_uMaxParticles];
	memset(m_hSurf, 0, sizeof(uint)*m_uMaxParticles);

	m_hAttr = new int[size1];
	memset(m_hAttr, 0, size1*sizeof(uint));
	m_hTmp = new RXREAL[m_uMaxParticles];
	memset(m_hTmp, 0, sizeof(RXREAL)*m_uMaxParticles);

	m_hDebugVec = new RXREAL[size];
	memset(m_hDebugVec, 0, mem_size);

	// Anisotropic kernel
	m_hUpPos = new RXREAL[size];
	m_hPosW = new RXREAL[size];
	m_hEigen = new RXREAL[m_uMaxParticles*3];
	m_hRMatrix = new RXREAL[m_uMaxParticles*9];
	m_hG = new RXREAL[m_uMaxParticles*9];

	m_vNeighs.resize(m_uMaxParticles);

	if(m_bUseOpenGL){
		m_posVBO = createVBO(mem_size);	
		m_colorVBO = createVBO(m_uMaxParticles*4*sizeof(RXREAL));

		SetColorVBO(RX_RAMP);
	}

	// �����Z���ݒ�
	m_pNNGrid->Setup(m_v3EnvMin, m_v3EnvMax, m_fEffectiveRadius, m_uMaxParticles);
	m_vNeighs.resize(m_uMaxParticles);

	if(m_uNumBParticles){
		Vec3 minp = m_pBoundary->GetMin()-Vec3(4.0*m_fParticleRadius);
		Vec3 maxp = m_pBoundary->GetMax()+Vec3(4.0*m_fParticleRadius);
		m_pNNGridB->Setup(minp, maxp, m_fEffectiveRadius, m_uNumBParticles);

		// �����Z���ɗ��q��o�^
		m_pNNGridB->SetObjectToCell(m_hPosB, m_uNumBParticles);

		//rxNNGrid::rxCell cell = m_pNNGridB->GetCellData();
		//Dump<uint>("_hash_cpu.txt", cell.hGridParticleHash, m_uNumBParticles, 1);
		//Dump<uint>("_index_cpu.txt", cell.hSortedIndex, m_uNumBParticles, 1);


		// ���E�p�[�e�B�N���̑̐�
		m_hVolB = new RXREAL[m_uNumBParticles];
		memset(m_hVolB, 0, sizeof(RXREAL)*m_uNumBParticles);
		calBoundaryVolumes(m_hPosB, m_hVolB, m_fMass, m_uNumBParticles, m_fEffectiveRadius);

		//Dump<RXREAL>("_volb_cpu.txt", m_hVolB, m_uNumBParticles, 1);
	}

	m_bInitialized = true;
}

/*!
 * �m�ۂ����������̉��
 */
void rxSPH::Finalize(void)
{
	assert(m_bInitialized);

	// ���������
	if(m_hPos) delete [] m_hPos;
	if(m_hVel) delete [] m_hVel;
	if(m_hNrm) delete [] m_hNrm;
	if(m_hFrc) delete [] m_hFrc;

	if(m_hDens) delete [] m_hDens;
	if(m_hPres) delete [] m_hPres;

	if(m_hPosB) delete [] m_hPosB;
	if(m_hVolB) delete [] m_hVolB;

	if(m_hSurf) delete [] m_hSurf;
	if(m_hAttr) delete [] m_hAttr;
	if(m_hTmp) delete [] m_hTmp;

	if(m_hDebugVec) delete [] m_hDebugVec;

	// Anisotoropic kernel
	if(m_hUpPos) delete [] m_hUpPos;
	if(m_hPosW) delete [] m_hPosW;
	if(m_hEigen) delete [] m_hEigen;
	if(m_hRMatrix) delete [] m_hRMatrix;
	if(m_hG) delete [] m_hG;

	m_vNeighs.clear();

	if(m_bUseOpenGL){
		glDeleteBuffers(1, (const GLuint*)&m_posVBO);
		glDeleteBuffers(1, (const GLuint*)&m_colorVBO);
	}

	if(m_pNNGrid) delete m_pNNGrid;
	m_vNeighs.clear();
	if(m_pNNGridB) delete m_pNNGridB;

	if(m_hVrts) delete [] m_hVrts;
	if(m_hTris) delete [] m_hTris;

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
 * SPH��1�X�e�b�v�i�߂�
 * @param[in] dt ���ԃX�e�b�v��
 * @retval ture  �v�Z����
 * @retval false �ő�X�e�b�v���𒴂��Ă��܂�
 */
bool rxSPH::Update(RXREAL dt, int step)
{
	// �����p�[�e�B�N����ǉ�
	if(!m_vInletLines.empty()){
		int start = (m_iInletStart == -1 ? 0 : m_iInletStart);
		int num = 0;
		vector<rxInletLine>::iterator itr = m_vInletLines.begin();
		for(; itr != m_vInletLines.end(); ++itr){
			rxInletLine iline = *itr;
			if(iline.span > 0 && step%(iline.span) == 0){
				int count = addParticles(m_iInletStart, iline);
				num += count;
			}
		}
		SetArrayVBO(RX_POSITION, m_hPos, start, num);
		SetArrayVBO(RX_VELOCITY, m_hVel, start, num);
	}


	assert(m_bInitialized);


	RXREAL h = m_fEffectiveRadius;
	int min_iter = 2;
	int max_iter = 20;
	RXREAL eta = 0.03;

	static bool init = true;

	// �ߖT���q�T��
	SetParticlesToCell();

	// ���x�v�Z
	calDensity(m_hPos, m_hDens, h);

	//if(!init){	// �^�C���X�e�b�v���̏C��
	//	dt = calTimeStep(dt, eta, m_hFrc, m_hVel, m_hDens);
	//}

	calForce(m_hPos, m_hVel, m_hDens, m_hPres, m_hFrc, h);

//	RXTIMER("force calculation");

	// �ʒu�E���x�̍X�V
	integrate(m_hPos, m_hVel, m_hDens, m_hFrc, m_hPos, m_hVel, dt);

//	RXTIMER("update position");


	SetArrayVBO(RX_POSITION, m_hPos, 0, m_uNumParticles);

	SetColorVBO(m_iColorType);

	RXTIMER("color(vbo)");

	init = false;
	return true;
}

/*!
 * �p�[�e�B�N���f�[�^�̎擾
 * @return �p�[�e�B�N���f�[�^�̃f�o�C�X�������|�C���^
 */
RXREAL* rxSPH::GetParticle(void)
{
	return m_hPos;
}

/*!
 * �p�[�e�B�N���f�[�^�̎擾
 * @return �p�[�e�B�N���f�[�^�̃f�o�C�X�������|�C���^
 */
RXREAL* rxSPH::GetParticleDevice(void)
{
	if(!m_uNumParticles) return 0;

	RXREAL *dPos = 0;
	CuAllocateArray((void**)&dPos, m_uNumParticles*4*sizeof(RXREAL));

	CuCopyArrayToDevice(dPos, m_hPos, 0, m_uNumParticles*4*sizeof(RXREAL));

	return dPos;
}


/*!
 * �J���[�l�pVBO�̕ҏW
 * @param[in] type �F�̃x�[�X�Ƃ��镨���l
 */
void rxSPH::SetColorVBO(int type)
{
	switch(type){
	case RX_DENSITY:
		SetColorVBOFromArray(m_hDens, 1, false, m_fRestDens*2.0);
		break;

	case RX_PRESSURE:
		SetColorVBOFromArray(m_hPres, 1);
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

	case RX_CONSTANT:
		if(m_bUseOpenGL){
			// �J���[�o�b�t�@�ɒl��ݒ�
			glBindBufferARB(GL_ARRAY_BUFFER, m_colorVBO);
			RXREAL *data = (RXREAL*)glMapBufferARB(GL_ARRAY_BUFFER, GL_WRITE_ONLY);
			RXREAL *ptr = data;
			for(uint i = 0; i < m_uNumParticles; ++i){
				RXREAL t = i/(RXREAL)m_uNumParticles;
				*ptr++ = 0.15f;
				*ptr++ = 0.15f;
				*ptr++ = 0.95f;
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
 * ���x���v�Z
 */
void rxSPH::calDensity(const RXREAL *ppos, RXREAL *pdens, RXREAL h)
{
	for(uint i = 0; i < m_uNumParticles; ++i){
		Vec3 pos0;
		pos0[0] = ppos[DIM*i+0];
		pos0[1] = ppos[DIM*i+1];
		pos0[2] = ppos[DIM*i+2];

		pdens[i] = 0.0;

		// �ߖT���q
		for(vector<rxNeigh>::iterator itr = m_vNeighs[i].begin() ; itr != m_vNeighs[i].end(); ++itr){
			int j = itr->Idx;
			if(j < 0) continue;

			Vec3 pos1;
			pos1[0] = ppos[DIM*j+0];
			pos1[1] = ppos[DIM*j+1];
			pos1[2] = ppos[DIM*j+2];

			RXREAL r = sqrt(itr->Dist2);

			pdens[i] += m_fMass*m_fpW(r, h, m_fAw);
		}
	}


#ifdef RX_USE_BOUNDARY_PARTICLE
	// ���E�p�[�e�B�N���̉e�����v�Z
	for(uint i = 0; i < m_uNumParticles; ++i){
		Vec3 pos0;
		pos0[0] = ppos[DIM*i+0];
		pos0[1] = ppos[DIM*i+1];
		pos0[2] = ppos[DIM*i+2];

		// �ߖT���q
		vector<rxNeigh> neigh;
		GetNearestNeighborsB(pos0, neigh, h);

		RXREAL brho = 0.0;
		for(vector<rxNeigh>::iterator itr = neigh.begin() ; itr != neigh.end(); ++itr){
			int j = itr->Idx;
			if(j < 0) continue;

			Vec3 pos1;
			pos1[0] = m_hPosB[DIM*j+0];
			pos1[1] = m_hPosB[DIM*j+1];
			pos1[2] = m_hPosB[DIM*j+2];

			RXREAL r = norm(pos1-pos0);

			brho += m_fRestDens*m_hVolB[j]*m_fpW(r, h, m_fAw);
		}

		pdens[i] += brho;
	}
#endif
}


/*!
 * ���͍��C�S�����̌v�Z
 */
void rxSPH::calForce(const RXREAL *ppos, const RXREAL *pvel, const RXREAL *pdens, RXREAL *ppres, RXREAL *pfrc, RXREAL h)
{
	RXREAL r0 = m_fRestDens;
	
	//RXREAL B = 1119000;
	for(uint i = 0; i < m_uNumParticles; ++i){
		pfrc[4*i+0] = 0.0;
		pfrc[4*i+1] = 0.0;
		pfrc[4*i+2] = 0.0;
		pfrc[4*i+3] = 0.0;

		ppres[i] = m_fGasStiffness*(pdens[i]-r0);
		//ppres[i] = B*(pow((double)(pdens[i]/r0), 7.0)-1.0);
	}


	for(uint i = 0; i < m_uNumParticles; ++i){
		Vec3 pos0, vel0;
		pos0 = Vec3(ppos[4*i+0], ppos[4*i+1], ppos[4*i+2]);
		vel0 = Vec3(pvel[4*i+0], pvel[4*i+1], pvel[4*i+2]);

		RXREAL prsi = ppres[i]/(pdens[i]*pdens[i]);

		Vec3 Fp(0.0), Fev(0.0);
		for(vector<rxNeigh>::iterator itr = m_vNeighs[i].begin() ; itr != m_vNeighs[i].end(); ++itr){
			int j = itr->Idx;
			if(j < 0 || i == j) continue;

			Vec3 pos1, vel1;
			pos1 = Vec3(ppos[4*j+0], ppos[4*j+1], ppos[4*j+2]);
			vel1 = Vec3(pvel[4*j+0], pvel[4*j+1], pvel[4*j+2]);

			Vec3 rij = pos0-pos1;
			Vec3 vji = vel1-vel0;

			RXREAL r = norm(rij);//sqrt(itr->Dist2);

			RXREAL prsj = ppres[j]/(pdens[j]*pdens[j]);

			// ����
			Fp += m_fMass*(prsi+prsj)*m_fpGW(r, h, m_fAg, rij);

			// �S��
			Fev += m_fMass*(vji/m_hDens[j])*m_fpLW(r, h, m_fAl, 3);
		}

		Vec3 force(0.0);
		force += -Fp;
		force += m_fViscosity*Fev;

		// �O�͍�
		force += m_v3Gravity;

		pfrc[4*i+0] += force[0];
		pfrc[4*i+1] += force[1];
		pfrc[4*i+2] += force[2];
	}

#ifdef RX_USE_BOUNDARY_PARTICLE
	// ���E�p�[�e�B�N���̉e�����v�Z
	for(uint i = 0; i < m_uNumParticles; ++i){
		Vec3 pos0;
		pos0[0] = ppos[DIM*i+0];
		pos0[1] = ppos[DIM*i+1];
		pos0[2] = ppos[DIM*i+2];

		RXREAL prsi = ppres[i]/(pdens[i]*pdens[i]);

		// �ߖT���q
		vector<rxNeigh> neigh;
		GetNearestNeighborsB(pos0, neigh, h);

		Vec3 bp(0.0);
		for(vector<rxNeigh>::iterator itr = neigh.begin() ; itr != neigh.end(); ++itr){
			int j = itr->Idx;
			if(j < 0) continue;

			Vec3 pos1;
			pos1[0] = m_hPosB[DIM*j+0];
			pos1[1] = m_hPosB[DIM*j+1];
			pos1[2] = m_hPosB[DIM*j+2];

			Vec3 rij = pos0-pos1;
			RXREAL r = norm(rij);

			bp += -m_fRestDens*m_hVolB[j]*prsi*m_fpGW(r, h, m_fAg, rij)/m_fMass;
		}


		m_hDebugVec[DIM*i+0] = bp[0];
		m_hDebugVec[DIM*i+1] = bp[1];
		m_hDebugVec[DIM*i+2] = bp[2];

		pfrc[DIM*i+0] += bp[0];
		pfrc[DIM*i+1] += bp[1];
		pfrc[DIM*i+2] += bp[2];
	}
#endif
}


/*!
 * ���E�p�[�e�B�N���̑̐ς��v�Z
 *  - "Versatile Rigid-Fluid Coupling for Incompressible SPH", 2.2 ��(3)�̏�
 * @param[in] bpos ���E�p�[�e�B�N���̈ʒu
 * @param[out] bvol ���E�p�[�e�B�N���̑̐�
 * @param[in] mass �p�[�e�B�N������
 * @param[in] h �L�����a
 */
void rxSPH::calBoundaryVolumes(const RXREAL *bpos, RXREAL *bvol, RXREAL mass, uint n, RXREAL h)
{
	for(uint i = 0; i < n; ++i){
		Vec3 pos0;
		pos0[0] = bpos[DIM*i+0];
		pos0[1] = bpos[DIM*i+1];
		pos0[2] = bpos[DIM*i+2];

		// �ߖT���q
		vector<rxNeigh> neigh;
		GetNearestNeighborsB(pos0, neigh, h);

		RXREAL mw = 0.0;
		for(vector<rxNeigh>::iterator itr = neigh.begin() ; itr != neigh.end(); ++itr){
			int j = itr->Idx;
			if(j < 0) continue;

			Vec3 pos1;
			pos1[0] = bpos[DIM*j+0];
			pos1[1] = bpos[DIM*j+1];
			pos1[2] = bpos[DIM*j+2];

			RXREAL r = norm(pos1-pos0);

			mw += mass*m_fpW(r, h, m_fAw);
		}

		bvol[i] = mass/mw;
	}
}



/*!
 * rest density�̌v�Z
 *  - �ߖT�Ƀp�[�e�B�N�����~���l�߂��Ă���Ƃ��Ė��x���v�Z����
 * @param[in] h �L�����a
 * @return rest density
 */
RXREAL rxSPH::calRestDensity(RXREAL h)
{
	RXREAL r0 = 0.0;
	RXREAL l = 2*GetParticleRadius();
	int n = (int)ceil(m_fKernelRadius/l)+1;
	for(int x = -n; x <= n; ++x){
		for(int y = -n; y <= n; ++y){
			for(int z = -n; z <= n; ++z){
				Vec3 rij = Vec3(x*l, y*l, z*l);
				r0 += m_fMass*m_fpW(norm(rij), h, m_fAw);
			}
		}
	}
	return r0;
}


/*!
 * ���ԃX�e�b�v���̏C��
 *  - Ihmsen et al., "Boundary Handling and Adaptive Time-stepping for PCISPH", Proc. VRIPHYS, pp.79-88, 2010.
 * @param[in] dt ���݂̃^�C�v�X�e�b�v
 * @param[in] eta_avg ���x�ϓ��̕���(���[�U�w��)
 * @param[in] pfrc ���͍��ɂ��͏�(���ۂɂ͉����x)
 * @param[in] pvel �p�[�e�B�N�����x
 * @param[in] pdens �p�[�e�B�N�����x
 * @return �C�����ꂽ�^�C�v�X�e�b�v��
 */
RXREAL rxSPH::calTimeStep(RXREAL &dt, RXREAL eta_avg, const RXREAL *pfrc, const RXREAL *pvel, const RXREAL *pdens)
{
	RXREAL h = m_fEffectiveRadius;
	RXREAL r0 = m_fRestDens;
	RXREAL new_dt = dt;

	// �ő�́C�ő呬�x�C���x�΍��̕��ςƍő���Z�o
	RXREAL ft_max = 0.0;
	RXREAL vt_max = 0.0;
	RXREAL rerr_max = 0.0, rerr_avg = 0.0;
	for(uint i = 0; i < m_uNumParticles; ++i){
		Vec3 f, v;
		for(int k = 0; k < 3; ++k){
			f[k] = pfrc[DIM*i+k];
			v[k] = pvel[DIM*i+k];
		}
		RXREAL ft = norm(f);
		RXREAL vt = norm(v);

		if(ft > ft_max) ft_max = ft;
		if(vt > vt_max) vt_max = vt;

		RXREAL rerr = (pdens[i]-r0)/r0;
		if(rerr > rerr_max) rerr_max = rerr;
		rerr_avg += rerr;
	}
	rerr_avg /= (RXREAL)m_uNumParticles;

	ft_max = sqrt(h/ft_max);
	vt_max = h/vt_max;

	int inc = 0;

	if(0.19*ft_max > dt && rerr_max < 4.5*eta_avg && rerr_avg < 0.9*eta_avg && 0.39*vt_max > dt){
		inc = 1;
	}
	else{
		if(0.2*ft_max < dt && rerr_max > 5.5*eta_avg && rerr_avg >= eta_avg && 0.4*vt_max <= dt){
			inc = -1;
		}
	}

	new_dt += 0.002*inc*dt;

	return new_dt;
}


/*!
 * �@�����v�Z
 */
void rxSPH::calNormal(void)
{
	RXREAL h = m_fEffectiveRadius;

	for(uint i = 0; i < m_uNumParticles; ++i){
		Vec3 pos0;
		pos0[0] = m_hPos[4*i+0];
		pos0[1] = m_hPos[4*i+1];
		pos0[2] = m_hPos[4*i+2];

		Vec3 nrm(0.0);

		// �ߖT���q
		for(vector<rxNeigh>::iterator itr = m_vNeighs[i].begin() ; itr != m_vNeighs[i].end(); ++itr){
			int j = itr->Idx;
			if(j < 0 || i == j) continue;

			Vec3 pos1;
			pos1[0] = m_hPos[4*j+0];
			pos1[1] = m_hPos[4*j+1];
			pos1[2] = m_hPos[4*j+2];

			Vec3 rij = pos0-pos1;

			RXREAL r = sqrt(itr->Dist2);

			nrm += m_fMass*m_fpW(r, h, m_fAw)*rij/m_hDens[j];
		}

		//normalize(nrm);
		//nrm *= -1;

		m_hNrm[4*i+0] = nrm[0];
		m_hNrm[4*i+1] = nrm[1];
		m_hNrm[4*i+2] = nrm[2];
		m_hNrm[4*i+3] = 0.0;
	}

}


/*!
 * ���C/�����ƎO�p�`�̌���
 * @param[in] P0,P1 ���C/�����̒[�_or���C��̓_
 * @param[in] V0,V1,V2 �O�p�`�̒��_���W
 * @param[out] I ��_���W
 * @retval 1 ��_I�Ō��� 
 * @retval 0 ��_�Ȃ�
 * @retval 2 �O�p�`�̕��ʓ�
 * @retval -1 �O�p�`��"degenerate"�ł���(�ʐς�0�C�܂�C�������_�ɂȂ��Ă���)
 */
static int intersectSegmentTriangle(Vec3 P0, Vec3 P1,			// Segment
									Vec3 V0, Vec3 V1, Vec3 V2,	// Triangle
									Vec3 &I, Vec3 &n, float rp)			// Intersection point (return)
{
	// �O�p�`�̃G�b�W�x�N�g���Ɩ@��
	Vec3 u = V1-V0;
	Vec3 v = V2-V0;
	n = Unit(cross(u, v));
	if(RXFunc::IsZeroVec(n)){
		return -1;	// �O�p�`��"degenerate"�ł���(�ʐς�0)
	}

	// ����
	Vec3 dir = P1-P0;
	double a = dot(n, P0-V0);
	double b = dot(n, dir);
	if(fabs(b) < 1e-10){	// �����ƎO�p�`���ʂ����s
		//if(fabs(a+b) < rp){
		//	I = P1+n*(rp-fabs(a+b));
		//	return 1;
		//}
		//else if(fabs(a) < rp){
		//	I = P1+n*(rp-fabs(a));
		//	return 1;
		//}
		//else{
			if(a == 0){
				return 2;	// ���������ʏ�
			}
			else{
				return 0;	// ��_�Ȃ�
			}
		//}
	}

	// ��_�v�Z

	// 2�[�_�����ꂼ��قȂ�ʂɂ��邩�ǂ����𔻒�
	float r = -a/b;
	Vec3 offset = Vec3(0.0);
	float dn = 0;
	float sign_n = 1;
	//if(r < 0.0 || fabs(a) > fabs(b) || b > 0){
	//	return 0;
	//}
	if(a < 0){
		return 0;
	}

	if(r < 0.0){
		//if(fabs(a) < rp){
		//	r = 0.0;
		//	offset = -fabs(a)*n;

		//	dn = fabs(a);
		//}
		//else{
			return 0;
		//}
	}
	else{
		if(fabs(a) > fabs(b)){
			//if(fabs(a+b) < rp){
			//	r = 1.0;
			//	offset = -fabs(a+b)*n;

			//	dn = fabs(a+b);
			//}
			//else{
				return 0;
			//}
		}
		else{
			if(b > 0){
				return 0;
			}
		}
	}

	// �����ƕ��ʂ̌�_
	I = P0+r*dir;//+offset;

	// ��_���O�p�`���ɂ��邩�ǂ����̔���
	double uu, uv, vv, wu, wv, D;
	uu = dot(u, u);
	uv = dot(u, v);
	vv = dot(v, v);
	Vec3 w = I-V0;
	wu = dot(w, u);
	wv = dot(w, v);
	D = uv*uv-uu*vv;

	double s, t;
	s = (uv*wv-vv*wu)/D;
	if(s < 0.0 || s > 1.0){
		return 0;
	}
	
	t = (uv*wu-uu*wv)/D;
	if(t < 0.0 || (s+t) > 1.0){
		return 0;
	}

	//I -= offset;

	// ���a�������߂�
	//Vec3 back = -dir*((rp-dn)/(dot(-dir, sign_n*n)));
	//I += back;

	return 1;
}

/*!
 * �|���S���I�u�W�F�N�g�Ƃ̏Փ˔���C�Փˉ���
 * @param[in] grid_hash ��������O���b�h�̃n�b�V��
 * @param[in] pos0 �O�X�e�b�v�̈ʒu
 * @param[inout] pos1 �V�����ʒu
 * @param[inout] vel ���x
 * @param[in] dt �^�C���X�e�b�v��
 * @return �Փ˃I�u�W�F�N�g�̐�
 */
int rxSPH::calCollisionPolygon(uint grid_hash, Vec3 &pos0, Vec3 &pos1, Vec3 &vel, RXREAL dt)
{
	vector<int> polys_in_cell;

	int c = 0;
	if(m_pNNGrid->GetPolygonsInCell(grid_hash, polys_in_cell)){
		for(int j = 0; j < (int)polys_in_cell.size(); ++j){
			int pidx = polys_in_cell[j];

			int vidx[3];
			vidx[0] = m_hTris[3*pidx+0];
			vidx[1] = m_hTris[3*pidx+1];
			vidx[2] = m_hTris[3*pidx+2];

			Vec3 vrts[3];
			vrts[0] = Vec3(m_hVrts[3*vidx[0]], m_hVrts[3*vidx[0]+1], m_hVrts[3*vidx[0]+2]);
			vrts[1] = Vec3(m_hVrts[3*vidx[1]], m_hVrts[3*vidx[1]+1], m_hVrts[3*vidx[1]+2]);
			vrts[2] = Vec3(m_hVrts[3*vidx[2]], m_hVrts[3*vidx[2]+1], m_hVrts[3*vidx[2]+2]);

			Vec3 cp, n;
			if(intersectSegmentTriangle(pos0, pos1, vrts[0], vrts[1], vrts[2], cp, n, m_fParticleRadius) == 1){
				double d = length(pos1-cp);
				n = Unit(n);

				RXREAL res = m_fRestitution;
				res = (res > 0) ? (res*fabs(d)/(dt*length(vel))) : 0.0f;
				Vec3 vr = -(1.0+res)*n*dot(n, vel);

				double l = norm(pos1-pos0);
				pos1 = cp+vr*(dt*d/l);
				vel += vr;

				c++;
			}
		}
	}

	return c;
}


/*!
 * �ő̃I�u�W�F�N�g�Ƃ̏Փ˔���C�Փˉ���
 * @param[in] pos0 �O�X�e�b�v�̈ʒu
 * @param[inout] pos1 �V�����ʒu
 * @param[inout] vel ���x
 * @param[in] dt �^�C���X�e�b�v��
 * @return �Փ˃I�u�W�F�N�g�̐�
 */
int rxSPH::calCollisionSolid(Vec3 &pos0, Vec3 &pos1, Vec3 &vel, RXREAL dt)
{
	int c = 0;
	rxCollisionInfo coli;

	// �ő̃I�u�W�F�N�g�Ƃ̏Փˏ���
	for(vector<rxSolid*>::iterator i = m_vSolids.begin(); i != m_vSolids.end(); ++i){
		if((*i)->GetDistanceR(pos1, m_fParticleRadius, coli)){
			RXREAL res = m_fRestitution;
			res = (res > 0) ? (res*fabs(coli.Penetration())/(dt*norm(vel))) : 0.0f;
			//vel -= (1+res)*dot(vel, coli.Normal())*coli.Normal();
			pos1 = coli.Contact();
		}
	}

	// �V�~�����[�V������ԋ��E�Ƃ̏Փˏ���
	if(m_pBoundary->GetDistanceR(pos1, m_fParticleRadius, coli)){
		RXREAL res = m_fRestitution;
		res = (res > 0) ? (res*fabs(coli.Penetration())/(dt*norm(vel))) : 0.0f;
		//vel -= (1+res)*dot(vel, coli.Normal())*coli.Normal();
		pos1 = coli.Contact();
	}

	return c;
}

/*!
 * �ʒu�E���x�̍X�V
 * @param[in] pos �p�[�e�B�N���ʒu
 * @param[in] vel �p�[�e�B�N�����x
 * @param[in] frc �p�[�e�B�N���ɂ������
 * @param[out] pos_new �X�V�p�[�e�B�N���ʒu
 * @param[out] vel_new �X�V�p�[�e�B�N�����x
 * @param[in] dt �^�C���X�e�b�v��
 */
void rxSPH::integrate(const RXREAL *pos, const RXREAL *vel, const RXREAL *dens, const RXREAL *acc, 
						  RXREAL *pos_new, RXREAL *vel_new, RXREAL dt)
{
	for(uint i = 0; i < m_uNumParticles; ++i){
		Vec3 x, x_old, v, a, v_old;
		for(int k = 0; k < 3; ++k){
			x[k] = pos[DIM*i+k];
			v[k] = vel[DIM*i+k];
			a[k] = acc[DIM*i+k];
		}
		x_old = x;

		// �V�������x�ƈʒu
		v += dt*a;
		x += dt*v;

		// �|���S���I�u�W�F�N�g�Ƃ̌�������
		if(m_iNumTris != 0){
			uint grid_hash0 = m_pNNGrid->CalGridHash(x_old);
			calCollisionPolygon(grid_hash0, x_old, x, v, dt);

			uint grid_hash1 = m_pNNGrid->CalGridHash(x);
			if(grid_hash1 != grid_hash0){
				calCollisionPolygon(grid_hash1, x_old, x, v, dt);
			}
		}

		// ���E�Ƃ̏Փ˔���
		calCollisionSolid(x_old, x, v, dt);

		// �V�������x�ƈʒu�ōX�V
		for(int k = 0; k < 3; ++k){
			pos_new[DIM*i+k] = x[k];
			vel_new[DIM*i+k] = v[k];
		}

	}
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
void rxSPH::DetectSurfaceParticles(void)
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
double rxSPH::CalDistToNormalizedMassCenter(const int i)
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
uint* rxSPH::GetArraySurf(void)
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
int rxSPH::GetSurfaceParticles(const Vec3 pos, RXREAL h, vector<rxSurfaceParticle> &sps)
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
void rxSPH::CalNormal(void)
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
void rxSPH::CalNormalFromDensity(void)
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

			nrm += m_fMass*m_fpW(r, h, m_fAw)*rij/m_hDens[j];
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
void rxSPH::GetNearestNeighbors(int idx, RXREAL *prts, vector<rxNeigh> &neighs, RXREAL h)
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
void rxSPH::GetNearestNeighbors(Vec3 pos, vector<rxNeigh> &neighs, RXREAL h)
{
	if(h < 0.0) h = m_fEffectiveRadius;

	m_pNNGrid->GetNN(pos, m_hPos, m_uNumParticles, neighs, h);
	//m_pNNGrid->GetNN_Direct(pos, prts, m_uNumParticles, neighs, h);	// �O���b�h���g��Ȃ���������
}

/*!
 * �ߖT���E���q�T��
 * @param[in] idx �T�����S�p�[�e�B�N���C���f�b�N�X
 * @param[out] neighs �T�����ʊi�[����ߖT���R���e�i
 * @param[in] h �T�����a
 */
void rxSPH::GetNearestNeighborsB(Vec3 pos, vector<rxNeigh> &neighs, RXREAL h)
{
	if(h < 0.0) h = m_fEffectiveRadius;

	m_pNNGridB->GetNN(pos, m_hPosB, m_uNumBParticles, neighs, h);
}

/*!
 * �S�p�[�e�B�N���𕪊��Z���Ɋi�[
 */
void rxSPH::SetParticlesToCell(RXREAL *prts, int n, RXREAL h)
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
void rxSPH::SetParticlesToCell(void)
{
	SetParticlesToCell(m_hPos, m_uNumParticles, m_fEffectiveRadius);
}


/*!
 * �����Z���Ɋi�[���ꂽ�|���S�������擾
 * @param[in] gi,gj,gk �Ώە����Z��
 * @param[out] polys �|���S��
 * @return �i�[�|���S����
 */
int rxSPH::GetPolygonsInCell(int gi, int gj, int gk, vector<int> &polys)
{
	return m_pNNGrid->GetPolygonsInCell(gi, gj, gk, polys);
}
/*!
 * �����Z�����̃|���S���̗L���𒲂ׂ�
 * @param[in] gi,gj,gk �Ώە����Z��
 * @return �|���S�����i�[����Ă����true
 */
bool rxSPH::IsPolygonsInCell(int gi, int gj, int gk)
{
	return m_pNNGrid->IsPolygonsInCell(gi, gj, gk);
}


/*!
 * �|���S���𕪊��Z���Ɋi�[
 */
void rxSPH::SetPolygonsToCell(void)
{
	m_pNNGrid->SetPolygonsToCell(m_hVrts, m_iNumVrts, m_hTris, m_iNumTris);
}



/*!
 * �T���p�Z���̕`��
 * @param[in] i,j,k �O���b�h��̃C���f�b�N�X
 */
void rxSPH::DrawCell(int i, int j, int k)
{
	if(m_pNNGrid) m_pNNGrid->DrawCell(i, j, k);
}

/*!
 * �T���p�O���b�h�̕`��
 * @param[in] col �p�[�e�B�N�����܂܂��Z���̐F
 * @param[in] col2 �|���S�����܂܂��Z���̐F
 * @param[in] sel �����_���ɑI�����ꂽ�Z���̂ݕ`��(1�ŐV�����Z����I���C2�ł��łɑI������Ă���Z����`��C0�ł��ׂẴZ����`��)
 */
void rxSPH::DrawCells(Vec3 col, Vec3 col2, int sel)
{
	if(m_pNNGrid) m_pNNGrid->DrawCells(col, col2, sel, m_hPos);

}

/*!
 * �ő̏�Q���̕`��
 */
void rxSPH::DrawObstacles(void)
{
	for(vector<rxSolid*>::iterator i = m_vSolids.begin(); i != m_vSolids.end(); ++i){
		(*i)->Draw();
	}
}




/*!
 * �O�p�`�|���S���ɂ���Q��
 * @param[in] vrts ���_
 * @param[in] tris ���b�V��
 */
void rxSPH::SetPolygonObstacle(const vector<Vec3> &vrts, const vector<Vec3> &nrms, const vector< vector<int> > &tris, Vec3 vel)
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
void rxSPH::SetBoxObstacle(Vec3 cen, Vec3 ext, Vec3 ang, Vec3 vel, int flg)
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
void rxSPH::SetSphereObstacle(Vec3 cen, double rad, Vec3 vel, int flg)
{
	rxSolidSphere *sphere = new rxSolidSphere(cen, rad, 1);
	m_vSolids.push_back(sphere);
}

/*!
 * ���^��Q���𓮂���
 * @param[in] b ���̔ԍ�
 * @param[in] disp �ړ���
 */
void rxSPH::MoveSphereObstacle(int b, Vec3 disp)
{
	if(m_vSolids.empty()) return;
	m_vSolids[b]->SetPosition(m_vSolids[b]->GetPosition()+disp);
}

/*!
 * ���^��Q���̈ʒu���擾
 * @param[in] b ���̔ԍ�
 */
Vec3 rxSPH::GetSphereObstaclePos(int b)
{
	if(m_vSolids.empty()) return Vec3(0.0);
	return m_vSolids[b]->GetPosition();
}


/*!
 * VBO����z�X�g�������փf�[�^��]���C�擾
 * @param[in] type �f�[�^�̎��
 * @return �z�X�g��������̃f�[�^
 */
RXREAL* rxSPH::GetArrayVBO(rxParticleArray type, bool d2h, int num)
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

	case RX_BOUNDARY_PARTICLE:
		hdata = m_hPosB;
		break;

	case RX_UPDATED_POSITION:
		hdata = m_hUpPos;
		break;

	case RX_EIGEN_VALUE:
		hdata = m_hEigen;
		break;

	case RX_ROTATION_MATRIX:
		hdata = m_hRMatrix;
		break;

	case RX_TRANSFORMATION:
		hdata = m_hG;
		break;

	case RX_DEBUG_VECTOR:
		hdata = m_hDebugVec;
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
void rxSPH::SetArrayVBO(rxParticleArray type, const RXREAL* data, int start, int count)
{
	assert(m_bInitialized);
 
	switch(type){
	default:
	case RX_POSITION:
		{
			if(m_bUseOpenGL){
				glBindBuffer(GL_ARRAY_BUFFER, m_posVBO);
				glBufferSubData(GL_ARRAY_BUFFER, start*4*sizeof(RXREAL), count*4*sizeof(RXREAL), data);
				glBindBuffer(GL_ARRAY_BUFFER, 0);
			}
		}
		break;

	case RX_VELOCITY:
		//CuCopyArrayToDevice(m_dVel, data, start*DIM*sizeof(RXREAL), count*DIM*sizeof(RXREAL));
		//{
		//	if(m_bUseOpenGL){
		//		glBindBuffer(GL_ARRAY_BUFFER, m_pos);
		//		glBufferSubData(GL_ARRAY_BUFFER, start*4*sizeof(RXREAL), count*4*sizeof(RXREAL), data);
		//		glBindBuffer(GL_ARRAY_BUFFER, 0);
		//	}
		//}
		break;

	case RX_NORMAL:
		//CuCopyArrayToDevice(m_dNrm, data, start*DIM*sizeof(RXREAL), count*DIM*sizeof(RXREAL));
		break;
	}	   
}
void rxSPH::SetArrayVBO(rxParticleArray type, const int* data, int start, int count)
{
	assert(m_bInitialized);
 
	switch(type){
	default:
	case RX_ATTRIBUTE:
		break;
	}	   
}




//-----------------------------------------------------------------------------
// MARK:�A�֐��l
//-----------------------------------------------------------------------------
double rxSPH::GetImplicit(double x, double y, double z)
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
void rxSPH::CalImplicitField(int n[3], Vec3 minp, Vec3 d, RXREAL *hF)
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
 * @param[in] pnx,pny,pnz �O���b�h���̎w�� nx=2^pnx
 * @param[in] minp �O���b�h�̍ŏ����W
 * @param[in] d �O���b�h��
 * @param[out] hF �A�֐��l(nx�~ny�~nz�̔z��)
 */
void rxSPH::CalImplicitFieldDevice(int n[3], Vec3 minp, Vec3 d, RXREAL *dF)
{
	RXREAL *hF = new RXREAL[n[0]*n[1]*n[2]];

	CalImplicitField(n, minp, d, hF);
	CuCopyArrayToDevice(dF, hF, 0, n[0]*n[1]*n[2]*sizeof(RXREAL));

	delete [] hF;
}

/*!
 * �J���[�t�B�[���h�l�v�Z
 * @param[in] pos �v�Z�ʒu
 * @return �J���[�t�B�[���h�l
 */
double rxSPH::CalColorField(double x, double y, double z)
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

		c += m_fMass*m_fpW(r, h, m_fAw);
	}

	return c;
}


//�ǉ��@�X�p
double rxSPH::GetImplicitSolid(double x, double y, double z)
{
	return CalColorFieldSolid(x, y, z);
}

/*!�ǉ��@�X�p
 * �p�[�e�B�N������O���b�h�̉A�֐��l���v�Z
 * @param[in] n �O���b�h��
 * @param[in] minp �O���b�h�̍ŏ����W
 * @param[in] d �O���b�h��
 * @param[out] hF �A�֐��l(nx�~ny�~nz�̔z��)
 */
void rxSPH::CalImplicitFieldSolid(int n[3], Vec3 minp, Vec3 d, RXREAL *hF)
{
	int slice0 = n[0];
	int slice1 = n[0]*n[1];

	for(int k = 0; k < n[2]; ++k){
		for(int j = 0; j < n[1]; ++j){
			for(int i = 0; i < n[0]; ++i){
				int idx = k*slice1+j*slice0+i;
				Vec3 pos = minp+Vec3(i, j, k)*d;
				hF[idx] = GetImplicitSolid(pos[0], pos[1], pos[2]);
			}
		}
	}
}

/*!�ǉ��@�X�p
 * �p�[�e�B�N������O���b�h�̉A�֐��l���v�Z
 * @param[in] pnx,pny,pnz �O���b�h���̎w�� nx=2^pnx
 * @param[in] minp �O���b�h�̍ŏ����W
 * @param[in] d �O���b�h��
 * @param[out] hF �A�֐��l(nx�~ny�~nz�̔z��)
 */
void rxSPH::CalImplicitFieldDeviceSolid(int n[3], Vec3 minp, Vec3 d, RXREAL *dF, float *m_fIntrps)
{
	RXREAL *hF = new RXREAL[n[0]*n[1]*n[2]];

	CalImplicitFieldSolid(n, minp, d, hF);
	CuCopyArrayToDevice(dF, hF, 0, n[0]*n[1]*n[2]*sizeof(RXREAL));

	delete [] hF;
}

/*!
 * �ǉ��@�X�p�J���[�t�B�[���h�l�v�Z
 * @param[in] pos �v�Z�ʒu
 * @return �J���[�t�B�[���h�l
 */
double rxSPH::CalColorFieldSolid(double x, double y, double z)
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
	//�O���b�h���_�̋ߖT���q�T���@GetNN()��GPU�ɂ͂Ȃ��̂ŁC�ȉ��Œ�`
	//RXREAL h = m_params.EffectiveRadius;

	//for(uint l = 0; l < m_uNumParticles; ++l){
	//	Vec3 pos = Vec3(m_hPos[DIM*l+0], m_hPos[DIM*l+1], m_hPos[DIM*l+2]);
	//	m_vNeighs[l].clear();

	//	// �����Z���C���f�b�N�X�̎Z�o
	//	int x = (pos[0]-m_v3EnvMin[0])/m_params.CellWidth.x;
	//	int y = (pos[1]-m_v3EnvMin[1])/m_params.CellWidth.y;
	//	int z = (pos[2]-m_v3EnvMin[2])/m_params.CellWidth.z;

	//	int numArdGrid = (int)(h/m_params.CellWidth.x+1);
	//	for(int k = -numArdGrid; k <= numArdGrid; ++k){
	//		for(int j = -numArdGrid; j <= numArdGrid; ++j){
	//			for(int i = -numArdGrid; i <= numArdGrid; ++i){
	//				int i1 = x+i;
	//				int j1 = y+j;
	//				int k1 = z+k;
	//				if(i1 < 0 || i1 >= m_params.GridSize.x || j1 < 0 || j1 >= m_params.GridSize.y || k1 < 0 || k1 >= m_params.GridSize.z){
	//					continue;
	//				}

	//				getNeighborsInCell(pos, m_hPos, i1, j1, k1, m_vNeighs[l], h);
	//			}
	//		}
	//	}
	//}

	// �ߖT���q
	for(vector<rxNeigh>::iterator itr = ne.begin(); itr != ne.end(); ++itr){
		int j = itr->Idx;
		if(j < 0) continue;

		Vec3 pos1;
		pos1[0] = m_hPos[DIM*j+0];
		pos1[1] = m_hPos[DIM*j+1];
		pos1[2] = m_hPos[DIM*j+2];

		RXREAL r = sqrt(itr->Dist2);

		c += m_fMass*m_fpW(r, h, m_fAw);
	}

	return c;
}

//-----------------------------------------------------------------------------
// �ٕ����J�[�l��
//-----------------------------------------------------------------------------

/*!
 * Anisotropic kernel�̌v�Z
 *  - J. Yu and G. Turk, Reconstructing Surfaces of Particle-Based Fluids Using Anisotropic Kernels, SCA2010. 
 */
void rxSPH::CalAnisotropicKernel(void)
{
	// MARK:CalAnisotropicKernel
	if(m_uNumParticles == 0) return;

	// VBO���烁������
	RXREAL r0 = m_fRestDens;
	RXREAL h0 = m_fEffectiveRadius;
	RXREAL h = 2.0*h0;

	// 2�{�̒T�����a�ŋߖT�T��
	SetParticlesToCell(m_hPos, m_uNumParticles, h);

	RXREAL lambda = 0.9;

	// �X�V�ʒu�̌v�Z
	for(uint i = 0; i < m_uNumParticles; ++i){
		Vec3 pos0;
		pos0[0] = m_hPos[4*i+0];
		pos0[1] = m_hPos[4*i+1];
		pos0[2] = m_hPos[4*i+2];

		// �ߖT���q
		Vec3 posw(0.0);
		RXREAL sumw = 0.0f;
		for(vector<rxNeigh>::iterator itr = m_vNeighs[i].begin() ; itr != m_vNeighs[i].end(); ++itr){
			int j = itr->Idx;
			if(j < 0) continue;

			Vec3 pos1;
			pos1[0] = m_hPos[4*j+0];
			pos1[1] = m_hPos[4*j+1];
			pos1[2] = m_hPos[4*j+2];

			RXREAL r = sqrt(itr->Dist2);
			
			if(r < h){
				RXREAL q = 1-r/h;
				RXREAL wij = q*q*q;
				posw += pos1*wij;
				sumw += wij;
			}
		}

		m_hPosW[DIM*i+0] = posw[0]/sumw;
		m_hPosW[DIM*i+1] = posw[1]/sumw;
		m_hPosW[DIM*i+2] = posw[2]/sumw;
		m_hPosW[DIM*i+3] = 0.0f;

		m_hUpPos[DIM*i+0] = (1-lambda)*m_hPos[DIM*i+0]+lambda*m_hPosW[DIM*i+0];
		m_hUpPos[DIM*i+1] = (1-lambda)*m_hPos[DIM*i+1]+lambda*m_hPosW[DIM*i+1];
		m_hUpPos[DIM*i+2] = (1-lambda)*m_hPos[DIM*i+2]+lambda*m_hPosW[DIM*i+2];

		//m_hUpPos[DIM*i+0] = m_hPos[DIM*i+0];
		//m_hUpPos[DIM*i+1] = m_hPos[DIM*i+1];
		//m_hUpPos[DIM*i+2] = m_hPos[DIM*i+2];
	}

	SetParticlesToCell(m_hUpPos, m_uNumParticles, h);

	// covariance matrix�̌v�Z
	for(uint i = 0; i < m_uNumParticles; ++i){
		Vec3 xi;
		xi[0] = m_hUpPos[4*i+0];
		xi[1] = m_hUpPos[4*i+1];
		xi[2] = m_hUpPos[4*i+2];

		Vec3 xiw(0.0);
		RXREAL sumw = 0.0f;
		for(vector<rxNeigh>::iterator itr = m_vNeighs[i].begin() ; itr != m_vNeighs[i].end(); ++itr){
			int j = itr->Idx;
			if(j < 0) continue;

			Vec3 xj;
			xj[0] = m_hUpPos[4*j+0];
			xj[1] = m_hUpPos[4*j+1];
			xj[2] = m_hUpPos[4*j+2];

			RXREAL r = sqrt(itr->Dist2);
			
			if(r < h){
				RXREAL q = 1-r/h;
				RXREAL wij = q*q*q;
				xiw += xj*wij;
				sumw += wij;
			}
		}

		xiw /= sumw;

		matrix3x3 c;
		for(int k = 0; k < 3; ++k){
			c.e[k].x = 0.0;
			c.e[k].y = 0.0;
			c.e[k].z = 0.0;
		}

		int n = 0;
		sumw = 0.0f;
		for(vector<rxNeigh>::iterator itr = m_vNeighs[i].begin() ; itr != m_vNeighs[i].end(); ++itr){
			int j = itr->Idx;
			if(j < 0) continue;

			Vec3 xj;
			xj[0] = m_hUpPos[4*j+0];
			xj[1] = m_hUpPos[4*j+1];
			xj[2] = m_hUpPos[4*j+2];

			RXREAL r = sqrt(itr->Dist2);

			if(r < h){
				RXREAL q = 1-r/h;
				RXREAL wij = q*q*q;

				Vec3 dxj = xj-xiw;

				for(int k = 0; k < 3; ++k){
					c.e[k].x += wij*dxj[k]*dxj[0];
					c.e[k].y += wij*dxj[k]*dxj[1];
					c.e[k].z += wij*dxj[k]*dxj[2];
				}

				sumw += wij;

				n++;
			}
		}

		for(int k = 0; k < 3; ++k){
			c.e[k].x /= sumw;
			c.e[k].y /= sumw;
			c.e[k].z /= sumw;
		}

		// ���ْl����
		float w[3], u[9], v[9];
		for(int k = 0; k < 3; ++k){
			u[k*3+0] = c.e[k].x;
			u[k*3+1] = c.e[k].y;
			u[k*3+2] = c.e[k].z;
		}

		RxSVDecomp3(w, u, v, 1.0e-10);
		
		// �ŗL�l��
		Vec3 sigma;
		for(int j = 0; j < 3; ++j){
			sigma[j] = w[j];
		}

		// �ŗL�x�N�g��(��]�s��R)
		rxMatrix3 R;
		for(int j = 0; j < 3; ++j){
			for(int k = 0; k < 3; ++k){
				R(k, j) = u[k*3+j];
			}
		}

		// �f�o�b�O�p�ɒl��Ҕ�
		for(int j = 0; j < 3; ++j){
			m_hEigen[3*i+j] = sigma[j];
		}
		for(int j = 0; j < 3; ++j){
			for(int k = 0; k < 3; ++k){
				m_hRMatrix[9*i+3*j+k] = R(j, k);
			}
		}

		
		int ne = 10;//m_params.KernelParticles*0.8;
		RXREAL ks = 1400;
		RXREAL kn = 0.5;
		RXREAL kr = 4.0;
		if(n > ne){
			for(int j = 1; j < 3; ++j){
				sigma[j] = RX_MAX(sigma[j], sigma[0]/kr);
			}
			sigma *= ks;
		}
		else{
			for(int j = 0; j < 3; ++j){
				sigma[j] = kn*1.0;
			}
		}


		// �J�[�l���ό`�s��G
		rxMatrix3 G;
		for(int j = 0; j < 3; ++j){
			for(int k = 0; k < 3; ++k){
				RXREAL x = 0;
				for(int l = 0; l < 3; ++l){
					//x += m_hRMatrix[9*i+3*j+l]*m_hRMatrix[9*i+3*k+l]/m_hEigen[3*i+l];
					x += R(j, l)*R(k, l)/sigma[l];
				}
				
				G(j, k) = x;
				//m_hG[9*i+3*j+k] = x;
			}
		}

		double max_diag = -1.0e10;
		for(int j = 0; j < 3; ++j){
			for(int k = 0; k < 3; ++k){
				if(G(j, k) > max_diag) max_diag = G(j, k);
			}
		}

		for(int j = 0; j < 3; ++j){
			for(int k = 0; k < 3; ++k){
				G(j, k) /= max_diag;
			}
		}
				
		G = G.Inverse();


		for(int j = 0; j < 3; ++j){
			for(int k = 0; k < 3; ++k){
				m_hG[9*i+3*j+k] = G(j, k);
			}
		}

	}
}



//-----------------------------------------------------------------------------
// MARK:�V�~���f�[�^�̏o��
//-----------------------------------------------------------------------------
/*!
 * �V�~�����[�V�����ݒ�(�p�[�e�B�N�����C�͈́C���x�C���ʂȂ�)
 * @param[in] fn �o�̓t�@�C����
 */
void rxSPH::OutputSetting(string fn)
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
	fout << m_fRestDens << endl;
	fout << m_fMass << endl;
	fout << m_iKernelParticles << endl;

	fout.close();
}

