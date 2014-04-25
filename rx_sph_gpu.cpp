/*!
  @file rx_sph_gpu.cpp
	
  @brief SPH�@(GPU)�̎���
 
  @author Makoto Fujisawa
  @date   2008-10,2011-06
*/
// FILE --rx_sph_gpu.cpp--


//-----------------------------------------------------------------------------
// �C���N���[�h�t�@�C��
//-----------------------------------------------------------------------------
#include "rx_sph.h"

#include "rx_cu_funcs.cuh"
#include <cuda_runtime.h>

#include "rx_pcube.h"



//-----------------------------------------------------------------------------
// �O���[�o���ϐ�
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
// rxSPH_GPU�N���X�̎���
//-----------------------------------------------------------------------------
/*!
 * �R���X�g���N�^
 * @param[in] use_opengl VBO�g�p�t���O
 */
rxSPH_GPU::rxSPH_GPU(bool use_opengl) : 
	rxParticleSystemBase(use_opengl), 
	m_hNrm(0), 
	m_hFrc(0),
	m_hDens(0), 
	m_hPres(0), 
	m_hSurf(0), 
	m_hVwt(0), 
	m_hEt(0), 
	m_hTurb(0), 
	m_hUpPos(0), 
	m_hPosW(0), 
	m_hCMatrix(0), 
	m_hEigen(0), 
	m_hRMatrix(0), 
	m_hG(0), 
	m_dPos(0),
	m_dVel(0),
	m_dNrm(0), 
	m_dFrc(0), 
	m_dDens(0), 
	m_dPres(0), 
	m_dAttr(0), 
	m_dEt(0), 
	m_dTurb(0), 
	m_dUpPos(0), 
	m_dPosW(0), 
	m_dCMatrix(0), 
	m_dEigen(0), 
	m_dRMatrix(0), 
	m_dG(0), 
	m_dPosB(0), 
	m_dVolB(0),
	m_dVrts(0), 
	m_dTris(0)
{
	m_params.Gravity = make_float3(m_v3Gravity[0], m_v3Gravity[1], m_v3Gravity[2]);

	m_params.InitDensity = 1.0f;
	m_params.Dt = 0.01f;

	m_params.Pressure = 101325.0f;

	m_params.Tension = 0.0728f;
	m_params.Viscosity = 0.01f;
	m_params.GasStiffness = 3.0f;

	m_params.Restitution = 0.0f;

	m_params.VorticityConfinement = 1.0f;
	m_params.Threshold = 1.0f;
	m_params.Buoyancy = 0.0f;

	m_params.BoxNum = 0;
	m_params.SphereNum = 0;

	m_fRestitution = 0.0f;

	m_hNoiseTile[0] = 0; 
	m_hNoiseTile[1] = 0;
	m_hNoiseTile[2] = 0;

	m_bUseVorticity = false;
	m_bUseWaveletTurb = false;
	m_bGridVelocity = false;
	m_bCalNormal = false;
	m_bUpsampling = false;
	m_bSubParticle = false;
	m_bCalAnisotropic = false;

	m_fEigenMax = 1.0;

	m_dCellData.dSortedPolyIdx = 0;
	m_dCellData.dGridPolyHash = 0;
	m_dCellData.dPolyCellStart = 0;
	m_dCellData.dPolyCellEnd = 0;
	m_dCellData.uNumPolyHash = 0;

	m_dCellDataB.dSortedPolyIdx = 0;
	m_dCellDataB.dGridPolyHash = 0;
	m_dCellDataB.dPolyCellStart = 0;
	m_dCellDataB.dPolyCellEnd = 0;
	m_dCellDataB.uNumPolyHash = 0;

	m_hSubPosPack = 0;
	m_hSubRadPack = 0;
	m_hSubRatPack = 0;
	m_hSubOcc = 0;
	m_hSubOccScan = 0;

	m_fTime = 0.0f;
	m_uNumParticles = 0;
	m_iNumVrts = 0;
	m_iNumTris = 0;

	m_uNumBParticles = 0;

	m_bSubPacked = false;

	m_iColorType = RX_RAMP;
}

/*!
 * �f�X�g���N�^
 */
rxSPH_GPU::~rxSPH_GPU()
{
	Finalize();
	CuClearData();
}



/*!
 * �V�~�����[�V�����̏�����
 * @param[in] max_particles �ő�p�[�e�B�N����
 * @param[in] boundary_ext ���E�̑傫��(�e�ӂ̒�����1/2)
 * @param[in] dens �������x
 * @param[in] mass �p�[�e�B�N���̎���
 * @param[in] kernel_particle �L�����ah�ȓ��̃p�[�e�B�N����
 */
void rxSPH_GPU::Initialize(const rxSPHEnviroment &env)
{
	// MARK:Initialize
	RXCOUT << "[rxSPH_GPU::Initialize]" << endl;

	m_params.Density = env.dens;
	m_params.Mass    = env.mass;
	m_params.KernelParticles = env.kernel_particles;

	m_params.Volume = env.kernel_particles*m_params.Mass/m_params.Density;

	m_params.EffectiveRadius = pow(((3.0*m_params.Volume)/(4.0*RX_PI)), 1.0/3.0);
	m_fParticleRadius = pow((RX_PI/(6.0*m_params.KernelParticles)), 1.0/3.0)*m_params.EffectiveRadius;
	//m_fParticleRadius = 0.5f*m_params.EffectiveRadius;
	m_params.ParticleRadius = m_fParticleRadius;

	m_params.Viscosity = env.viscosity;
	m_params.GasStiffness = env.gas_k;

	RXREAL h = m_params.EffectiveRadius;
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
	m_params.Density = calRestDensity(h);

	RXCOUT << "particle : " << endl;
	RXCOUT << " n_max = " << env.max_particles << endl;
	RXCOUT << " h = " << m_params.EffectiveRadius << endl;
	RXCOUT << " r = " << m_params.ParticleRadius << endl;
	RXCOUT << " dens = " << m_params.Density << endl;
	RXCOUT << " mass = " << m_params.Mass << endl;
	RXCOUT << " kernel_particles = " << m_params.KernelParticles << endl;
	RXCOUT << " volume = " << m_params.Volume << endl;
	RXCOUT << " viscosity = " << m_params.Viscosity << endl;

	//
	// ���E�ݒ�
	//
	// ���E�͈̔�
	// �V�~�����[�V�������̑傫��
	m_v3EnvMin = env.boundary_cen-env.boundary_ext;
	m_v3EnvMax = env.boundary_cen+env.boundary_ext;
	RXCOUT << "simlation range : " << m_v3EnvMin << " - " << m_v3EnvMax << endl;

	m_params.BoundaryMin = MAKE_FLOAT3V(m_v3EnvMin);
	m_params.BoundaryMax = MAKE_FLOAT3V(m_v3EnvMax);

	Vec3 world_size = m_v3EnvMax-m_v3EnvMin;
	Vec3 world_origin = m_v3EnvMin;

	double expansion = 0.01;
	world_origin -= 0.5*expansion*world_size;
	world_size *= (1.0+expansion); // �V�~�����[�V�������S�̂𕢂��悤�ɐݒ�

	m_v3EnvMin = world_origin;
	m_v3EnvMax = world_origin+world_size;

	// �p�[�e�B�N���폜�̈�
	m_iDeleteRegion = env.use_delete_region;
	if(m_iDeleteRegion){
		for(int i = 0; i < m_iDeleteRegion; ++i){
			m_vDeleteRegion.push_back(env.delete_region[i][0]);
			m_vDeleteRegion.push_back(env.delete_region[i][1]); 

			cout << "region for deleting particles : " << m_vDeleteRegion[2*i] << " - " << m_vDeleteRegion[2*i+1] << endl;
		}
	}

	// �����Z���ݒ�
	double cell_width;
	setupCells(m_dCellData, m_gridSize, cell_width, m_v3EnvMin, m_v3EnvMax, h);
	m_params.CellWidth = MAKE_FLOAT3(cell_width, cell_width, cell_width);
	
	m_gridSortBits = 24;	// �O���b�h�����������͂��̒l�𑝂₷

	m_params.GridSize = m_gridSize;
	m_params.NumCells = m_dCellData.uNumCells;
	m_params.NumBodies = m_uNumParticles;

	m_params.WorldOrigin = MAKE_FLOAT3(world_origin[0], world_origin[1], world_origin[2]);
	m_params.WorldMax = MAKE_FLOAT3(m_v3EnvMax[0], m_v3EnvMax[1], m_v3EnvMax[2] );

	RXCOUT << "grid for nn search : " << endl;
	RXCOUT << "  size   : " << m_params.GridSize << endl;
	RXCOUT << "  num    : " << m_params.NumCells << endl;
	RXCOUT << "  origin : " << m_params.WorldOrigin << endl;
	RXCOUT << "  width  : " << m_params.CellWidth << endl;

	m_uNumArdGrid = (int)(m_params.EffectiveRadius/m_params.CellWidth.x)+1;
	RXCOUT << "  numArdGrid : " << m_uNumArdGrid << endl;

	// �E�F�[�u���b�g�ϊ��X�P�[��
	RXCOUT << "wavelet turbulence : " << endl;
	RXCOUT << " scale = " << g_fWaveletScale << endl;

	// GPU�ɓn�����߂̕����Z�����
	m_dCellData.uNumArdGrid = m_uNumArdGrid;
	m_dCellData.uNumCells = m_dCellData.uNumCells;
	m_dCellData.uNumParticles = m_uNumParticles;

	// �J�[�l���֐��̒萔
	m_params.Wpoly6  =  315.0/(64.0*RX_PI*pow(h, (RXREAL)9.0));
	m_params.GWpoly6 = -945.0/(32.0*RX_PI*pow(h, (RXREAL)9.0));
	m_params.LWpoly6 = -945.0/(32.0*RX_PI*pow(h, (RXREAL)9.0));

	m_params.Wspiky  =  15.0/(RX_PI*pow(h, (RXREAL)6.0));
	m_params.GWspiky = -45.0/(RX_PI*pow(h, (RXREAL)6.0));
	m_params.LWspiky = -90.0/(RX_PI*pow(h, (RXREAL)6.0));

	m_params.Wvisc   = 15.0/(2.0*RX_PI*pow(h, (RXREAL)3.0));
	m_params.GWvisc  = 15.0/(2.0*RX_PI*pow(h, (RXREAL)3.0));
	m_params.LWvisc  = 45.0/(RX_PI*pow(h, (RXREAL)6.0));

	m_uNumParticles = 0;

	Allocate(env.max_particles);
}

inline uint calUintPow(uint x, uint y)
{
	uint x_y = 1;
	for(uint i=0; i < y;i++) x_y *= x;
	return x_y;
}

/*!
 * �������̊m��
 *  - �ő�p�[�e�B�N�����Ŋm��
 * @param[in] max_particles �ő�p�[�e�B�N����
 */
void rxSPH_GPU::Allocate(int maxParticles)
{
	// MARK:Allocate
	assert(!m_bInitialized);

	//m_uNumParticles = maxParticles;
	m_uMaxParticles = maxParticles;
	m_uMaxSubLevel	= MAX_SUB_LEVEL;//maxSubLevel;

	unsigned int size  = m_uMaxParticles*DIM;
	unsigned int size1 = m_uMaxParticles;
	unsigned int mem_size  = sizeof(RXREAL)*size;
	unsigned int mem_size1 = sizeof(RXREAL)*size1;

	unsigned int size_sub	= m_uMaxParticles*( calUintPow(2, m_uMaxSubLevel+1) - 1);//�S�ẴT�u���q
	unsigned int size_sub_m	= m_uMaxParticles*( calUintPow(2, m_uMaxSubLevel)   - 1);//���x����(Max-1)�܂ł̃T�u���q
	unsigned int size_sub_e = (m_uMaxSubLevel+1) * m_uMaxParticles;

	m_uMaxSubParticles = size_sub;

	CuDeviceProp();

	//
	// CPU���������m��
	//
	// GPU�Ƃ̃f�[�^�]���������f�[�^�̓y�[�W���b�N�������Ɋm��
	cudaMallocHost((void**)&m_hPos, mem_size);
	cudaMallocHost((void**)&m_hVel, mem_size);
	cudaMallocHost((void**)&m_hAttr, sizeof(int)*size1);
	cudaMallocHost((void**)&m_hEt, mem_size1);
	cudaMallocHost((void**)&m_hTurb, mem_size);

	// �ʏ�̃������m��
	m_hNrm = new RXREAL[size];
	m_hFrc = new RXREAL[size];
	m_hDens = new RXREAL[size1];
	m_hPres = new RXREAL[size1];
	memset(m_hNrm, 0, mem_size);
	memset(m_hFrc, 0, mem_size);
	memset(m_hDens, 0, mem_size1);
	memset(m_hPres, 0, mem_size1);

	//�ǉ��F�ߖT���q
	m_vNeighs.resize(m_uMaxParticles);

	//�ǉ��F�����O���b�h�\���̂̔z��m��
	m_hSortedIndex = new uint[m_uMaxParticles];
	m_hGridParticleHash = new uint[m_uMaxParticles];
	m_hCellStart = new uint[m_dCellData.uNumCells];
	m_hCellEnd = new uint[m_dCellData.uNumCells];

	m_hSurf = new uint[m_uMaxParticles];
	memset(m_hSurf, 0, sizeof(uint)*m_uMaxParticles);
	m_hTmp = new RXREAL[m_uMaxParticles];
	memset(m_hTmp, 0, sizeof(RXREAL)*m_uMaxParticles);

	// �E�F�[�u���b�g����
	m_hVwt = new RXREAL[size];
	memset(m_hVwt, 0, mem_size);

	// SPS����
	m_hSubPos	= new RXREAL[size_sub*4];
	m_hSubAxis	= new RXREAL[size_sub_m*4];//!< �T�u���q�̉�]��(�P�ʃx�N�g��)
	m_hSubChild	= new RXREAL[size_sub_m*4];//!< �T�u���q�̎q1�ւ̒P�ʃx�N�g��
	m_hSubEt	= new RXREAL[size_sub_m];
	memset(m_hSubPos,	0, size_sub*4*sizeof(RXREAL));
	memset(m_hSubAxis,	0, size_sub_m*4*sizeof(RXREAL));
	memset(m_hSubChild,	0, size_sub_m*4*sizeof(RXREAL));
	memset(m_hSubEt,	0, size_sub_m*sizeof(RXREAL));

	// Anisotropic kernel
	m_hUpPos   = new RXREAL[size];
	m_hPosW    = new RXREAL[size];
	m_hCMatrix = new RXREAL[m_uMaxParticles*9];
	m_hEigen   = new RXREAL[m_uMaxParticles*3];
	m_hRMatrix = new RXREAL[m_uMaxParticles*9];
	m_hG       = new RXREAL[m_uMaxParticles*9];


	//
	// GPU���������m��
	//
	if(m_bUseOpenGL){
		cout << "use opengl" << endl;
		m_posVBO = createVBO(mem_size);	
		CuRegisterGLBufferObject(m_posVBO, &m_pPosResource);

		m_subposVBO	= createVBO(size_sub*4*sizeof(RXREAL));
		CuRegisterGLBufferObject(m_subposVBO, &m_pSubPosResource);

		m_colorVBO = createVBO(sizeof(RXREAL)*4*m_uMaxParticles);
		SetColorVBO(RX_RAMP);
	}
	else{
		CuAllocateArray((void**)&m_dPos, mem_size);
	}

	CuAllocateArray((void**)&m_dVel,    mem_size);
	CuAllocateArray((void**)&m_dNrm,    mem_size);
	CuAllocateArray((void**)&m_dFrc,    mem_size);
	CuAllocateArray((void**)&m_dDens,   mem_size1);
	CuAllocateArray((void**)&m_dPres,   mem_size1);
	CuAllocateArray((void**)&m_dAttr,   size1*sizeof(int));

	// �E�F�[�u���b�g����
	CuAllocateArray((void**)&m_dEt,     mem_size1);
	CuAllocateArray((void**)&m_dTurb,   mem_size);

	// SPS����
	CuAllocateArray((void**)&m_dSubPos,		size_sub*4*sizeof(RXREAL));
	CuAllocateArray((void**)&m_dSubAxis,	size_sub_m*4*sizeof(RXREAL));
	CuAllocateArray((void**)&m_dSubChild,	size_sub_m*4*sizeof(RXREAL));
	CuAllocateArray((void**)&m_dSubEt,		size_sub_m*sizeof(RXREAL));
	CuAllocateArray((void**)&m_dSubRand,	size_sub*sizeof(uint));

	// Anisotropic kernel
	CuAllocateArray((void**)&m_dUpPos,      mem_size);
	CuAllocateArray((void**)&m_dPosW,       mem_size);
	CuAllocateArray((void**)&m_dCMatrix,	m_uMaxParticles*9*sizeof(RXREAL));
	CuAllocateArray((void**)&m_dEigen,		m_uMaxParticles*3*sizeof(RXREAL));
	CuAllocateArray((void**)&m_dRMatrix,	m_uMaxParticles*9*sizeof(RXREAL));
	CuAllocateArray((void**)&m_dG,          m_uMaxParticles*9*sizeof(RXREAL));

	// �p�[�e�B�N���O���b�h�f�[�^
	CuAllocateArray((void**)&m_dCellData.dSortedPos, mem_size);
	CuAllocateArray((void**)&m_dCellData.dSortedVel, mem_size);
	CuAllocateArray((void**)&m_dCellData.dGridParticleHash,  m_uMaxParticles*sizeof(uint));
	CuAllocateArray((void**)&m_dCellData.dSortedIndex, m_uMaxParticles*sizeof(uint));
	CuAllocateArray((void**)&m_dCellData.dCellStart, m_dCellData.uNumCells*sizeof(uint));
	CuAllocateArray((void**)&m_dCellData.dCellEnd, m_dCellData.uNumCells*sizeof(uint));

	// �T�u�p�[�e�B�N���O���b�h�f�[�^
	m_dSubCellData.uSubNumAllParticles	= size_sub;
	m_dSubCellData.uSubNumCells			= m_dCellData.uNumCells;
	m_dSubCellData.uSubNumArdGrid		= m_dCellData.uNumArdGrid;
	m_dSubCellData.uSubMaxLevel			= MAX_SUB_LEVEL;
	m_dSubCellData.fEtcri				= g_fEtCri;
	CuAllocateArray((void**)&m_dSubCellData.dSubUnsortPos, size_sub*4*sizeof(RXREAL));
	CuAllocateArray((void**)&m_dSubCellData.dSubUnsortRad, size_sub*sizeof(RXREAL));
	CuAllocateArray((void**)&m_dSubCellData.dSubUnsortRat, size_sub*sizeof(RXREAL));
	CuAllocateArray((void**)&m_dSubCellData.dSubSortedPos, size_sub*4*sizeof(RXREAL));
	CuAllocateArray((void**)&m_dSubCellData.dSubSortedRad, size_sub*sizeof(RXREAL));
	CuAllocateArray((void**)&m_dSubCellData.dSubSortedRat, size_sub*sizeof(RXREAL));
	CuAllocateArray((void**)&m_dSubCellData.dSubOcc, size_sub*sizeof(uint));
	CuAllocateArray((void**)&m_dSubCellData.dSubOccScan, size_sub*sizeof(uint));
	CuAllocateArray((void**)&m_dSubCellData.dSubSortedIndex, size_sub*sizeof(uint));
	CuAllocateArray((void**)&m_dSubCellData.dSubGridParticleHash, size_sub*sizeof(uint));
	CuAllocateArray((void**)&m_dSubCellData.dSubCellStart, m_dSubCellData.uSubNumCells*sizeof(uint));
	CuAllocateArray((void**)&m_dSubCellData.dSubCellEnd, m_dSubCellData.uSubNumCells*sizeof(uint));

	CuSetParameters(&m_params);
	m_bInitialized = true;

#ifdef RX_USE_BOUNDARY_PARTICLE
	InitBoundary();
#endif
}

void rxSPH_GPU::InitBoundary(void)
{
	rxSolid* boundary = new rxSolidBox(m_v3EnvMin, m_v3EnvMax, -1);

	// ���E�p�[�e�B�N������
	m_uNumBParticles = boundary->GenerateParticlesOnSurf(3.0*m_params.ParticleRadius, &m_hPosB);

	delete boundary;

	if(m_uNumBParticles){
		unsigned int sizeb  = m_uNumBParticles*DIM;
		unsigned int mem_sizeb  = sizeof(RXREAL)*sizeb;

		Vec3 minp = m_v3EnvMin-Vec3(4.0*m_fParticleRadius);
		Vec3 maxp = m_v3EnvMax+Vec3(4.0*m_fParticleRadius);

		// �����Z���ݒ�
		double cell_width;
		setupCells(m_dCellDataB, m_gridSizeB, cell_width, minp, maxp, m_params.EffectiveRadius);
		m_params.CellWidthB = MAKE_FLOAT3(cell_width, cell_width, cell_width);

		m_params.GridSizeB = m_gridSizeB;
		m_params.NumCellsB = m_dCellDataB.uNumCells;
		m_params.NumBodiesB = m_uNumBParticles;

		m_params.WorldOriginB = MAKE_FLOAT3(minp[0], minp[1], minp[2]);
		m_params.WorldMaxB = MAKE_FLOAT3(maxp[0], maxp[1], maxp[2] );

		RXCOUT << "grid for nn search (boundary) : " << endl;
		RXCOUT << "  size   : " << m_params.GridSizeB << endl;
		RXCOUT << "  num    : " << m_params.NumCellsB << endl;
		RXCOUT << "  origin : " << m_params.WorldOriginB << endl;
		RXCOUT << "  width  : " << m_params.CellWidthB << endl;

		// GPU�ɓn�����߂̕����Z�����
		m_dCellDataB.uNumArdGrid = m_uNumArdGrid;
		m_dCellDataB.uNumParticles = m_uNumBParticles;

		// ���E�p�[�e�B�N���O���b�h�f�[�^
		CuAllocateArray((void**)&m_dCellDataB.dSortedPos, mem_sizeb);
		CuAllocateArray((void**)&m_dCellDataB.dSortedVel, mem_sizeb);
		CuAllocateArray((void**)&m_dCellDataB.dGridParticleHash,  m_uNumBParticles*sizeof(uint));
		CuAllocateArray((void**)&m_dCellDataB.dSortedIndex, m_uNumBParticles*sizeof(uint));
		CuAllocateArray((void**)&m_dCellDataB.dCellStart, m_dCellDataB.uNumCells*sizeof(uint));
		CuAllocateArray((void**)&m_dCellDataB.dCellEnd, m_dCellDataB.uNumCells*sizeof(uint));

		// ���E�p�[�e�B�N���̈ʒu�����f�o�C�X�ɓ]��
		CuAllocateArray((void**)&m_dPosB, mem_sizeb);
		//SetArrayVBO(RX_BOUNDARY_PARTICLE, m_hPosB, 0, m_uNumBParticles);
		CuCopyArrayToDevice(m_dPosB, m_hPosB, 0, m_uNumBParticles*DIM*sizeof(RXREAL));

		
		CuSetParameters(&m_params);	//for(int i = 0; i < m_uNumBParticles*DIM; ++i) m_hPosB[i] = 0.0;

		//CuCopyArrayFromDevice(m_hPosB, m_dPosB, 0, m_uNumBParticles*DIM*sizeof(RXREAL));
		//Dump<RXREAL>("_posb_gpu.txt", m_hPosB, m_uNumBParticles, DIM);

		// �ߖT�T���������p�O���b�h�f�[�^�̍쐬
		// �����Z���̃n�b�V�����v�Z
		CuCalcHashB(m_dCellDataB.dGridParticleHash, m_dCellDataB.dSortedIndex, m_dPosB, 
					m_params.WorldOriginB, m_params.CellWidthB, m_params.GridSizeB, m_uNumBParticles);


		// �n�b�V���Ɋ�Â��p�[�e�B�N�����\�[�g
		CuSort(m_dCellDataB.dGridParticleHash, m_dCellDataB.dSortedIndex, m_uNumBParticles);

		// �p�[�e�B�N���z����\�[�g���ꂽ���Ԃɕ��ёւ��C
		// �e�Z���̎n�܂�ƏI���̃C���f�b�N�X������
		CuReorderDataAndFindCellStartB(m_dCellDataB, m_dPosB);

		//uint *tmp = new uint[m_uNumBParticles];
		//CuCopyArrayFromDevice(tmp, m_dCellDataB.dGridParticleHash, 0, m_uNumBParticles*sizeof(uint));
		//Dump<uint>("_hash_gpu.txt", tmp, m_uNumBParticles, 1);
		//CuCopyArrayFromDevice(tmp, m_dCellDataB.dSortedIndex, 0, m_uNumBParticles*sizeof(uint));
		//Dump<uint>("_index_gpu.txt", tmp, m_uNumBParticles, 1);
		//delete [] tmp;


		// ���E�p�[�e�B�N���̑̐�
		m_hVolB = new RXREAL[m_uNumBParticles];
		memset(m_hVolB, 0, sizeof(RXREAL)*m_uNumBParticles);
		CuAllocateArray((void**)&m_dVolB, sizeof(RXREAL)*m_uNumBParticles);

		CuSphBoundaryVolume(m_dVolB, m_params.Mass, m_dCellDataB);

		CuCopyArrayFromDevice(m_hVolB, m_dVolB, 0, m_uNumBParticles*sizeof(RXREAL));
		//m_hVolB = GetArrayVBO(RX_BOUNDARY_PARTICLE_VOL, true, m_uNumBParticles);

		//Dump<RXREAL>("_volb_gpu.txt", m_hVolB, m_uNumBParticles, 1);
	}
}

/*!
 * �m�ۂ����������̉��
 */
void rxSPH_GPU::Finalize(void)
{
	assert(m_bInitialized);

	// �y�[�W���b�N���������
	cudaFreeHost(m_hPos);
	cudaFreeHost(m_hVel);
	cudaFreeHost(m_hAttr);
	cudaFreeHost(m_hEt);
	cudaFreeHost(m_hTurb);

	// �ʏ탁�������
	if(m_hNrm) delete [] m_hNrm;
	if(m_hFrc) delete [] m_hFrc;
	if(m_hDens) delete [] m_hDens;
	if(m_hPres) delete [] m_hPres;

	if(m_hPosB) delete [] m_hPosB;
	if(m_hVolB) delete [] m_hVolB;

	if(m_hSurf) delete [] m_hSurf;
	if(m_hTmp) delete [] m_hTmp;

	// �E�F�[�u���b�g����
	if(m_hVwt) delete [] m_hVwt;

	// SPS����
	if(m_hSubPos) delete [] m_hSubPos;
	if(m_hSubAxis) delete [] m_hSubAxis;
	if(m_hSubChild) delete [] m_hSubChild;
	if(m_hSubEt) delete [] m_hSubEt;

	if(m_hSubPosPack) delete [] m_hSubPosPack;
	if(m_hSubRadPack) delete [] m_hSubRadPack;
	if(m_hSubRatPack) delete [] m_hSubRatPack;
	if(m_hSubOcc) delete [] m_hSubOcc;
	if(m_hSubOccScan) delete [] m_hSubOccScan;

	// Anisotoropic kernel
	if(m_hUpPos) delete [] m_hUpPos;
	if(m_hPosW) delete [] m_hPosW;
	if(m_hCMatrix) delete [] m_hCMatrix;
	if(m_hEigen) delete [] m_hEigen;
	if(m_hRMatrix) delete [] m_hRMatrix;
	if(m_hG) delete [] m_hG;

	m_vNeighs.clear();


	// GPU���������
	CuFreeArray(m_dVel);
	CuFreeArray(m_dNrm);
	CuFreeArray(m_dFrc);
	CuFreeArray(m_dDens);
	CuFreeArray(m_dPres);
	CuFreeArray(m_dAttr);

	CuFreeArray(m_dPosB);
	CuFreeArray(m_dVolB);

	CuFreeArray(m_dEt);
	CuFreeArray(m_dTurb);

	//�T�u���q
	CuFreeArray(m_dSubPos);
	CuFreeArray(m_dSubAxis);
	CuFreeArray(m_dSubChild);
	CuFreeArray(m_dSubEt);
	CuFreeArray(m_dSubRand);

	CuFreeArray(m_dUpPos);
	CuFreeArray(m_dPosW);
	CuFreeArray(m_dCMatrix);
	CuFreeArray(m_dEigen);
	CuFreeArray(m_dRMatrix);
	CuFreeArray(m_dG);

	CuFreeArray(m_dCellData.dSortedPos);
	CuFreeArray(m_dCellData.dSortedVel);
	CuFreeArray(m_dCellData.dGridParticleHash);
	CuFreeArray(m_dCellData.dSortedIndex);
	CuFreeArray(m_dCellData.dCellStart);
	CuFreeArray(m_dCellData.dCellEnd);

	CuFreeArray(m_dCellDataB.dSortedPos);
	CuFreeArray(m_dCellDataB.dSortedVel);
	CuFreeArray(m_dCellDataB.dGridParticleHash);
	CuFreeArray(m_dCellDataB.dSortedIndex);
	CuFreeArray(m_dCellDataB.dCellStart);
	CuFreeArray(m_dCellDataB.dCellEnd);

	CuFreeArray(m_dSubCellData.dSubUnsortPos);
	CuFreeArray(m_dSubCellData.dSubUnsortRad);
	CuFreeArray(m_dSubCellData.dSubUnsortRat);
	CuFreeArray(m_dSubCellData.dSubSortedPos);
	CuFreeArray(m_dSubCellData.dSubSortedRad);
	CuFreeArray(m_dSubCellData.dSubSortedRat);
	CuFreeArray(m_dSubCellData.dSubOcc);
	CuFreeArray(m_dSubCellData.dSubOccScan);
	CuFreeArray(m_dSubCellData.dSubGridParticleHash);
	CuFreeArray(m_dSubCellData.dSubSortedIndex);
	CuFreeArray(m_dSubCellData.dSubCellStart);
	CuFreeArray(m_dSubCellData.dSubCellEnd);

	if(m_dVrts) CuFreeArray(m_dVrts);
	if(m_dTris) CuFreeArray(m_dTris);

	if(m_dCellData.dSortedPolyIdx) CuFreeArray(m_dCellData.dSortedPolyIdx);
	if(m_dCellData.dGridPolyHash)  CuFreeArray(m_dCellData.dGridPolyHash);
	if(m_dCellData.dPolyCellStart) CuFreeArray(m_dCellData.dPolyCellStart);
	if(m_dCellData.dPolyCellEnd)   CuFreeArray(m_dCellData.dPolyCellEnd);

	if(m_dCellDataB.dSortedPolyIdx) CuFreeArray(m_dCellDataB.dSortedPolyIdx);
	if(m_dCellDataB.dGridPolyHash)  CuFreeArray(m_dCellDataB.dGridPolyHash);
	if(m_dCellDataB.dPolyCellStart) CuFreeArray(m_dCellDataB.dPolyCellStart);
	if(m_dCellDataB.dPolyCellEnd)   CuFreeArray(m_dCellDataB.dPolyCellEnd);

	//�ǉ��F
	if(m_hSortedIndex) delete [] m_hSortedIndex;
	if(m_hGridParticleHash) delete [] m_hGridParticleHash;
	if(m_hCellStart) delete [] m_hCellStart;
	if(m_hCellEnd) delete [] m_hCellEnd;


	if(m_bUseOpenGL){
		CuUnregisterGLBufferObject(m_pPosResource);
		glDeleteBuffers(1, (const GLuint*)&m_posVBO);
		CuUnregisterGLBufferObject(m_pSubPosResource);
		glDeleteBuffers(1, (const GLuint*)&m_subposVBO);

		glDeleteBuffers(1, (const GLuint*)&m_colorVBO);
	}
	else{
		CuFreeArray(m_dPos);
	}

	m_uNumParticles = 0;
	m_uMaxParticles = 0;

	m_uNumBParticles = 0;
	//CuSPHClean();
}



#if MAX_BOX_NUM
void rxSPH_GPU::ToggleSolidFlg(int t)
{
	if(!m_params.BoxNum) return;

	RX_TOGGLE(m_params.BoxFlg[m_params.BoxNum-1], t);
	CuSetParameters(&m_params);
}
#endif


/*!
 * GPU���̃p�����[�^�t�@�C�����X�V
 */
void rxSPH_GPU::UpdateParams(void)
{
	m_params.ParticleRadius = m_fParticleRadius;
	m_params.Gravity = make_float3(m_v3Gravity[0], m_v3Gravity[1], m_v3Gravity[2]);
	m_params.GlobalDamping = m_fDamping;

	m_params.Restitution = m_fRestitution;
	CuSetParameters(&m_params);
}

/*!
 * SPH��1�X�e�b�v�i�߂�
 * @param[in] dt ���ԃX�e�b�v��
 * @retval ture  �v�Z����
 * @retval false �ő�X�e�b�v���𒴂��Ă��܂�
 */
bool rxSPH_GPU::Update(RXREAL dt, int step)
{
	//RXTIMER_RESET;

	// �����p�[�e�B�N����ǉ�
	if(!m_vInletLines.empty()){
		int start = (m_iInletStart == -1 ? 0 : m_iInletStart);
		int num = 0;
		int attr = 0;
		vector<rxInletLine>::iterator itr = m_vInletLines.begin();
		for(; itr != m_vInletLines.end(); ++itr){
			rxInletLine iline = *itr;
			if(iline.span > 0 && step%(iline.span) == 0){
				int count = addParticles(m_iInletStart, iline, attr);
				//AddSubParticles(m_iInletStart, count);
				num += count;
			}
			attr++;
		}

		//�o�O���������̂Ő搶���C��
		if(num){
			SetArrayVBO(RX_POSITION, &m_hPos[DIM*start], start, num);
			SetArrayVBO(RX_VELOCITY, &m_hVel[DIM*start], start, num);
			SetArrayVBO(RX_ATTRIBUTE, &m_hAttr[DIM*start], start, num);
		}
	}

	if(!m_uNumParticles) return false;

	assert(m_bInitialized);

	static bool init = true;

	m_params.Dt = dt;
	UpdateParams();

	// GPU�p�ϐ��Ƀp�[�e�B�N������ݒ�
	m_params.NumBodies = m_uNumParticles;
	m_dCellData.uNumParticles = m_uNumParticles;

	// �p�[�e�B�N�����W�z����}�b�s���O
	RXREAL *dPos, *dSubPos;
	if(m_bUseOpenGL){
		dPos = (RXREAL*)CuMapGLBufferObject(&m_pPosResource);
		dSubPos = (RXREAL*)CuMapGLBufferObject(&m_pSubPosResource);
	}
	else{
		dPos = (RXREAL*)m_dPos;
		dSubPos = (RXREAL*)m_dSubPos;
	}

	//if(init){
	//	InitBoundary();
	//}

	for(uint j = 0; j < m_solverIterations; ++j){
		// �ߖT�T���������p�O���b�h�f�[�^�̍쐬
		// �����Z���̃n�b�V�����v�Z
		CuCalcHash(m_dCellData.dGridParticleHash, m_dCellData.dSortedIndex, dPos, m_dAttr, m_uNumParticles);

		// �n�b�V���Ɋ�Â��p�[�e�B�N�����\�[�g
		CuSort(m_dCellData.dGridParticleHash, m_dCellData.dSortedIndex, m_uNumParticles);

		// �p�[�e�B�N���z����\�[�g���ꂽ���Ԃɕ��ёւ��C
		// �e�Z���̎n�܂�ƏI���̃C���f�b�N�X������
		CuReorderDataAndFindCellStart(m_dCellData, dPos, m_dVel);

		//RXTIMER("cell construction");

		// ���x�v�Z
		CuSphDensity(m_dDens, m_dPres, m_dCellData);

#ifdef RX_USE_BOUNDARY_PARTICLE
		// ���E�p�[�e�B�N���ɂ�閧�x
		CuSphBoundaryDensity(m_dDens, m_dPres, dPos, m_dVolB, m_dCellDataB, m_uNumParticles);
#endif

		if(m_bCalNormal){
			// �@���v�Z
			CuSphNormal(m_dNrm, m_dDens, m_dCellData);

			m_hNrm = GetArrayVBO(RX_NORMAL);
		}

		// �p�[�e�B�N���ɂ������(�O�́C���́C�S���g�U�͂Ȃ�)���v�Z
		CuSphForces(m_dDens, m_dPres, m_dFrc, m_dCellData, dt);

#ifdef RX_USE_BOUNDARY_PARTICLE
		// ���E�p�[�e�B�N���ɂ���
		CuSphBoundaryForces(m_dDens, m_dPres, dPos, m_dVolB, m_dFrc, m_dCellDataB, m_uNumParticles);
#endif

//		RXTIMER("force calculation");

		// �p�[�e�B�N��CWT
		RXREAL scale = g_fWaveletScale*m_params.EffectiveRadius;
		if(m_bUseWaveletTurb){
			// �p�[�e�B�N������̃E�F�[�u���b�g��͂ɂ�藐���G�l���M���v�Z
			CalWaveletTurbulence(scale, dPos, dt);
			//RXTIMER("wavelet turbulence");
		}
		else if(m_iColorType == RX_ENERGY_SPECTRUM){
			// �p�[�e�B�N������̃E�F�[�u���b�g��͂ɂ�藐���G�l���M���v�Z(�`��p)
			CalParticleCWT(scale);
			RXTIMER("energy spectrum");
		}

		// �ʒu�C���x�̍X�V
		if(m_iNumTris == 0){
			CuSphIntegrate(dPos, m_dVel, m_dFrc, m_dDens, m_dAttr, dt, m_uNumParticles);
		}
		else{	// �|���S���ɂ��ő̋��E�L��
			CuSphIntegrateWithPolygon(dPos, m_dVel, m_dFrc, m_dDens, m_dAttr, 
									  m_dVrts, m_dTris, m_iNumTris, dt, m_dCellData);
		}


		// �p�[�e�B�N���̈ʒu�����Ďw��͈͂ɗ����炻�̃p�[�e�B�N���̑�����-1�ɂ���D
		if(m_iDeleteRegion){
			for(int k = 0; k < m_iDeleteRegion; ++k){
				float minp[3], maxp[3], farpoint[3] = {0, 0, 1000000};
				minp[0] = m_vDeleteRegion[2*k+0][0];
				minp[1] = m_vDeleteRegion[2*k+0][1];
				minp[2] = m_vDeleteRegion[2*k+0][2];
				maxp[0] = m_vDeleteRegion[2*k+1][0];
				maxp[1] = m_vDeleteRegion[2*k+1][1];
				maxp[2] = m_vDeleteRegion[2*k+1][2];
				CuSphCheckDelete(dPos, m_dVel, m_dAttr, minp, maxp, farpoint, m_uNumParticles);
				//CuSphCheckDeleteX(dPos, m_dVel, m_dAttr, minp[0], farpoint, m_uNumParticles);
			}
		}


//		RXTIMER("update position");

		// �T�u�p�[�e�B�N���ʒu�̍X�V
		if(m_bSubParticle){
			UpdateSubParticle(dPos ,dSubPos, scale, dt);		
			RXTIMER("update position(sp)");
		}

		init = false;
	}




	if(m_bUseOpenGL){
		CuUnmapGLBufferObject(m_pPosResource);
		CuUnmapGLBufferObject(m_pSubPosResource);
	}

	m_fTime += dt;

	// VBO���烁������
	//m_hPos = GetArrayVBO(RX_POSITION);
	//m_hVel = GetArrayVBO(RX_VELOCITY);

	SetColorVBO(m_iColorType);

	//�ǉ��FGPU����CPU�Ƀf�[�^���R�s�[
	searchNeighbors();		

	//RXTIMER("color(vbo)");

	m_bSubPacked = false;
	m_bCalAnisotropic = false;

	return true;
}

/*!
 * CPU�p�̋ߖT�T���@�ǉ��FGPU����CPU�Ƀf�[�^���R�s�[���C�ߖT�T�����s���Ă���D
 */
void rxSPH_GPU::searchNeighbors(void)
{
	//// �ߖT�T���������p�O���b�h�f�[�^�̍쐬
	//// �����Z���̃n�b�V�����v�Z
	//CuCalcHash(m_dCellData.dGridParticleHash, m_dCellData.dSortedIndex, m_dPos, m_dAttr, m_uNumParticles);

	//// �n�b�V���Ɋ�Â��p�[�e�B�N�����\�[�g
	//CuSort(m_dCellData.dGridParticleHash, m_dCellData.dSortedIndex, m_uNumParticles);

	//// �p�[�e�B�N���z����\�[�g���ꂽ���Ԃɕ��ёւ��C
	//// �e�Z���̎n�܂�ƏI���̃C���f�b�N�X������
	//CuReorderDataAndFindCellStart(m_dCellData, m_dPos, m_dVel);

	// GPU -> CPU
	CuCopyArrayFromDevice(m_hSortedIndex, m_dCellData.dSortedIndex, 0, m_uNumParticles*sizeof(uint));
	CuCopyArrayFromDevice(m_hGridParticleHash, m_dCellData.dGridParticleHash, 0, m_uNumParticles*sizeof(uint));
	CuCopyArrayFromDevice(m_hCellStart, m_dCellData.dCellStart, 0, m_params.NumCells*sizeof(uint));
	CuCopyArrayFromDevice(m_hCellEnd, m_dCellData.dCellEnd, 0, m_params.NumCells*sizeof(uint));

	m_hPos = GetArrayVBO(RX_POSITION);
	//CuCopyArrayFromDevice(m_hPos, m_dPos, 0, m_uNumParticles*DIM*sizeof(RXREAL));

	RXREAL h = m_params.EffectiveRadius;

	for(uint l = 0; l < m_uNumParticles; ++l){
		Vec3 pos = Vec3(m_hPos[DIM*l+0], m_hPos[DIM*l+1], m_hPos[DIM*l+2]);
		m_vNeighs[l].clear();

		// �����Z���C���f�b�N�X�̎Z�o
		int x = (pos[0]-m_v3EnvMin[0])/m_params.CellWidth.x;
		int y = (pos[1]-m_v3EnvMin[1])/m_params.CellWidth.y;
		int z = (pos[2]-m_v3EnvMin[2])/m_params.CellWidth.z;

		int numArdGrid = (int)(h/m_params.CellWidth.x+1);
		for(int k = -numArdGrid; k <= numArdGrid; ++k){
			for(int j = -numArdGrid; j <= numArdGrid; ++j){
				for(int i = -numArdGrid; i <= numArdGrid; ++i){
					int i1 = x+i;
					int j1 = y+j;
					int k1 = z+k;
					if(i1 < 0 || (unsigned)i1 >= m_params.GridSize.x || j1 < 0 || (unsigned)j1 >= m_params.GridSize.y || k1 < 0 || (unsigned)k1 >= m_params.GridSize.z){
						continue;
					}

					getNeighborsInCell(pos, m_hPos, i1, j1, k1, m_vNeighs[l], h);
				}
			}
		}
	}
}

/*!
 * �����Z�����̗��q����ߖT�����o
 * @param[in] pos �T�����S
 * @param[in] p �p�[�e�B�N���ʒu
 * @param[in] gi,gj,gk �Ώە����Z��
 * @param[out] neighs �T�����ʊi�[����ߖT���R���e�i
 * @param[in] h �L�����a
 */
void rxSPH_GPU::getNeighborsInCell(Vec3 pos, RXREAL *p, int gi, int gj, int gk, vector<rxNeigh> &neighs, RXREAL h)
{
	RXREAL h2 = h*h;

	uint grid_hash = gk*m_params.GridSize.y*m_params.GridSize.x+gj*m_params.GridSize.x+gi;

	uint start_index = m_hCellStart[grid_hash];
	if(start_index != 0xffffffff){	// �Z������łȂ����̃`�F�b�N
		uint end_index = m_hCellEnd[grid_hash];
		for(uint j = start_index; j < end_index; ++j){
			uint idx = m_hSortedIndex[j];

			Vec3 xij;
			xij[0] = pos[0]-p[DIM*idx+0];
			xij[1] = pos[1]-p[DIM*idx+1];
			xij[2] = pos[2]-p[DIM*idx+2];

			rxNeigh neigh;
			neigh.Dist2 = norm2(xij);

			if(neigh.Dist2 <= h2){
				neigh.Idx = idx;
				//neigh.Dist = sqrt(neigh.Dist2);

				neighs.push_back(neigh);
			}
		}
	}
}


/*!
 * rest density�̌v�Z
 *  - �ߖT�Ƀp�[�e�B�N�����~���l�߂��Ă���Ƃ��Ė��x���v�Z����
 * @param[in] h �L�����a
 * @return rest density
 */
RXREAL rxSPH_GPU::calRestDensity(RXREAL h)
{
	double a = KernelCoefPoly6(h, 3, 1);
	RXREAL r0 = 0.0;
	RXREAL l = 2*GetParticleRadius();
	int n = (int)ceil(m_params.EffectiveRadius/l)+1;
	for(int x = -n; x <= n; ++x){
		for(int y = -n; y <= n; ++y){
			for(int z = -n; z <= n; ++z){
				Vec3 rij = Vec3(x*l, y*l, z*l);
				r0 += m_params.Mass*KernelPoly6(norm(rij), h, a);
			}
		}
	}
	return r0;
}


/*!
 * �p�[�e�B�N���f�[�^�̎擾
 * @return �p�[�e�B�N���f�[�^�̃f�o�C�X�������|�C���^
 */
RXREAL* rxSPH_GPU::GetParticle(void)
{
	return m_hPos;
}

/*!
 * �p�[�e�B�N���f�[�^�̎擾
 * @return �p�[�e�B�N���f�[�^�̃f�o�C�X�������|�C���^
 */
RXREAL* rxSPH_GPU::GetParticleDevice(void)
{
	return (m_bUseOpenGL ? (RXREAL*)CuMapGLBufferObject(&m_pPosResource) : (RXREAL*)m_dPos);
}
/*!
 * �p�[�e�B�N���f�[�^�̃A���}�b�v(VBO)
 */
void rxSPH_GPU::UnmapParticle(void)
{
	if(m_bUseOpenGL) CuUnmapGLBufferObject(m_pPosResource);
}


/*!
 * �p�[�e�B�N�����x�̘A���E�F�[�u���b�g�ϊ�
 * @param[in] scale �E�F�[�u���b�g�X�P�[��
 */
void rxSPH_GPU::CalParticleCWT(RXREAL scale)
{
	RXREAL coef_et = g_fCoefEt*KOL2;
	RXREAL et_max  = 0.5;

	CuSphES(m_dEt, scale, coef_et, et_max, m_dCellData, m_params.CellWidth);
}

/*!
 * �p�[�e�B�N���E�F�[�u���b�g������̌v�Z
 * @param[in] scale �E�F�[�u���b�g�X�P�[��
 */
void rxSPH_GPU::CalWaveletTurbulence(RXREAL scale, RXREAL *dPos, RXREAL dt)
{
	static bool initialized = false;

	RXREAL coef_et = g_fCoefEt*KOL2;
	RXREAL et_max  = g_fMaxEt;

	// �p�[�e�B�N�����x�̃E�F�[�u���b�g�ϊ��ɂ�藐���G�l���M���v�Z
	CuSphES(m_dEt, scale, coef_et, et_max, m_dCellData, m_params.CellWidth);

	RXTIMER("energy spectrum");

	RXREAL r = m_params.ParticleRadius;
	Vec3 dim  = GetDim();
	Vec3 minp = GetMin();

	// �m�C�Y�^�C���쐬
	if(!initialized){
		int tnx = (int)(dim[0]/r)+1;
		int tny = (int)(dim[1]/r)+1;
		int tnz = (int)(dim[2]/r)+1;

		m_hNoiseTile[0] = rxWaveletNoise::GenerateNoiseTile3r(tnx, tny, tnz);
		CuSetWaveletTile3DBT(m_hNoiseTile[0], tnx, tny, tnz, 0);

		m_fNoiseTileWidth = dim[0]/tnx;
		m_iNoiseTileN[0] = tnx;
		m_iNoiseTileN[1] = tny;
		m_iNoiseTileN[2] = tnz;

		RXCOUT << "wavelet turbulence : " << endl;
		RXCOUT << " scale = " << scale << " (" << g_fWaveletScale << " x h)" << endl;
		RXCOUT << " coef_et = " << g_fCoefEt << endl;
		RXCOUT << " Etmax = " << et_max << endl;
		RXCOUT << " tile_n = (" << tnx << ", " << tny << ", " << tnz << ")" << endl;

		//if(m_bSubParticle){
		//	RXCOUT << " sp_et = " << g_fCoefTurb << endl;
		//	RXCOUT << " sp_et_mesh = " << g_fCoefTurbForMesh << endl;
		//	RXCOUT << " sp_et_cri = " << g_fEtCri << endl;
		//}

		RXTIMER_RESET;
	}

	RXREAL s = (scale)/m_fNoiseTileWidth;
	if(!initialized) RXCOUT << "s = " << s << endl;

	int minb = log((double)1.0/s)/log(2.0);
	
	// Wavelet turbulence�̒ǉ�
	int maxb = 0;//(int)(log((double)(m_v3EnvMax[0]-m_v3EnvMin[0])/(r*2))/log(2.0));
	int first = minb;
	int nband = maxb-minb;

	if(!initialized){
		RXCOUT << " band of turbulence :  " << minb << " - " << maxb << endl;
	}

	RXREAL pdim[3], pmin[3];
	for(int i = 0; i < 3; ++i){
		pdim[i] = dim[i]/(double)(m_iNoiseTileN[i]);
		pmin[i] = minp[i];
	}

	CuAddWTurb3D(m_dTurb, m_dFrc, m_dEt, m_dDens, 0, first, nband, dPos, pmin, pdim, m_uNumParticles, dt);

	RXTIMER("wavelet turbulence");

	initialized = true;
}


//-----------------------------------------------------------------------------
// MARK:�T�u�p�[�e�B�N��
//-----------------------------------------------------------------------------

/*!
 * �T�u�p�[�e�B�N���f�[�^�̎擾
 * @return �T�u�p�[�e�B�N���f�[�^�̃f�o�C�X�������|�C���^
 */
RXREAL* rxSPH_GPU::GetSubParticleDev(void)
{
	return (m_bUseOpenGL ? (RXREAL*)CuMapGLBufferObject(&m_pSubPosResource) : (RXREAL*)m_dSubPos);
}
void rxSPH_GPU::UnmapSubParticle(void)
{
	if(m_bUseOpenGL) CuUnmapGLBufferObject(m_pSubPosResource);
}

void rxSPH_GPU::InitSubParticles(const double &et_cri)
{
	if(!m_bSubParticle) return;

	RXCOUT << "InitSubParticles" << endl;

	SetEtcri(et_cri);
	
	m_dSubCellData.uNumParticles = m_uNumParticles;
	m_dSubCellData.uMaxParticles = m_uMaxParticles;
	for(uint level = 0; level <= m_dSubCellData.uSubMaxLevel; ++level){
		m_dSubCellData.uSubHeadIndex[level]			= m_dSubCellData.uMaxParticles * (calUintPow(2,level)-1);
		//m_dSubCellData.uSubNumLevel[level]			=;
		m_dSubCellData.uSubNumEach[level]			= calUintPow(2,level);
		m_dSubCellData.fSubRad[level]				= m_params.ParticleRadius * pow(2.0f,-(float)level/3.0f);
		m_dSubCellData.fSubEffectiveFactor			= m_params.EffectiveRadius / m_params.ParticleRadius;
		m_dSubCellData.fSubEffectiveRadius[level]	= m_dSubCellData.fSubEffectiveFactor * m_dSubCellData.fSubRad[level];
	}
	m_uNumSubParticles = m_uNumParticles * ( calUintPow(2, m_uMaxSubLevel+1) - 1);

	AddSubParticles(0, m_uNumParticles);
}

void rxSPH_GPU::AddSubParticles(int start, int count)
{
	if(!m_bSubParticle) return;

	RXCOUT << "AddSubParticles (" << count << "/" << m_uNumParticles << ") - " << m_params.ParticleRadius << endl;

	if(!count) return;

	m_dSubCellData.uNumParticles = m_uNumParticles;
	m_uNumSubParticles = m_uNumParticles * ( calUintPow(2, m_uMaxSubLevel+1) - 1);


	RXREAL *dPos,*dSubPos;
	if(m_bUseOpenGL){
		dPos = (RXREAL*)CuMapGLBufferObject(&m_pPosResource);
		dSubPos = (RXREAL*)CuMapGLBufferObject(&m_pSubPosResource);
	}
	else{
		dPos = (RXREAL*)m_dPos;
		dSubPos = (RXREAL*)m_dSubPos;
	}

	CuSubAdd(dSubPos, dPos, m_dSubAxis, m_dSubChild, m_dSubRand, 
			 m_params.ParticleRadius, m_uMaxSubLevel, m_uNumParticles, m_uMaxParticles, start, count);

	if(m_bUseOpenGL){
		CuUnmapGLBufferObject(m_pPosResource);
		CuUnmapGLBufferObject(m_pSubPosResource);
	}

}

/*!
 * �T�u�p�[�e�B�N���̈ʒu���X�V
 */
void rxSPH_GPU::UpdateSubParticle(RXREAL *dPos, RXREAL *dSubPos, RXREAL scale,RXREAL dt)
{
	float radius = m_params.ParticleRadius;

	CuSubUpdate(dPos, m_dEt, dSubPos, m_dSubChild, m_dSubAxis, m_dSubEt, m_dSubRand,
				g_fCoefTurb, radius, 0.50, m_uMaxSubLevel, m_uNumParticles, m_uMaxParticles, scale, dt);
}


/*!
 * �T�u�p�[�e�B�N���f�[�^�̎擾
 * @return �T�u�p�[�e�B�N���f�[�^�̃f�o�C�X�������|�C���^
 */
RXREAL* rxSPH_GPU::GetAllSubParticleDev(void)
{
	return (m_bUseOpenGL ? (RXREAL*)CuMapGLBufferObject(&m_pSubPosResource) : (RXREAL*)m_dSubPos);
}
void    rxSPH_GPU::UnmapAllSubParticle(void)
{
	if(m_bUseOpenGL) CuUnmapGLBufferObject(m_pSubPosResource);
}

/*!
 * �l�߂��T�u�p�[�e�B�N���f�[�^�̎擾
 * @return �l�߂��T�u�p�[�e�B�N���f�[�^�̃f�o�C�X�������|�C���^
 */
RXREAL* rxSPH_GPU::GetSubParticlePosDev(void)
{
	return (RXREAL*)m_dSubCellData.dSubUnsortPos;
}
RXREAL* rxSPH_GPU::GetSubParticleRadDev(void)
{
	return m_dSubCellData.dSubUnsortRad;
}
RXREAL* rxSPH_GPU::GetSubParticleRatDev(void)
{
	return m_dSubCellData.dSubUnsortRat;
}

// �T�u�p�[�e�B�N����
int rxSPH_GPU::GetNumAllSubParticles(void)
{
	return m_dSubCellData.uSubNumAllParticles;
}
int rxSPH_GPU::GetNumValidSubParticles(void)
{
	return m_dSubCellData.uSubNumValidParticles;
}


/*!
 * �l�߂��T�u�p�[�e�B�N���f�[�^�̎擾
 * @return �l�߂��T�u�p�[�e�B�N���f�[�^�̃f�o�C�X�������|�C���^
 */
RXREAL* rxSPH_GPU::GetSubParticlePos(void)
{
	if(!m_hSubPosPack){
		m_hSubPosPack = new RXREAL[4*GetNumAllSubParticles()];
	}
	int n = GetNumValidSubParticles();
	CuCopyArrayFromDevice(m_hSubPosPack, m_dSubCellData.dSubUnsortPos, 0, 4*n*sizeof(RXREAL));
	return m_hSubPosPack;
}
RXREAL* rxSPH_GPU::GetSubParticleRad(void)
{
	if(!m_hSubRadPack){
		m_hSubRadPack = new RXREAL[GetNumAllSubParticles()];
	}
	int n = GetNumValidSubParticles();
	CuCopyArrayFromDevice(m_hSubRadPack, m_dSubCellData.dSubUnsortRad, 0, n*sizeof(RXREAL));
	return m_hSubRadPack;
}
RXREAL* rxSPH_GPU::GetSubParticleRat(void)
{
	if(!m_hSubRatPack){
		m_hSubRatPack = new RXREAL[GetNumAllSubParticles()];
	}
	int n = GetNumValidSubParticles();
	CuCopyArrayFromDevice(m_hSubRatPack, m_dSubCellData.dSubUnsortRat, 0, n*sizeof(RXREAL));
	return m_hSubRadPack;
}
unsigned int* rxSPH_GPU::GetSubParticleOcc(void)
{
	if(!m_hSubOcc){
		m_hSubOcc = new unsigned int[GetNumAllSubParticles()];
	}
	int n = GetNumValidSubParticles();
	CuCopyArrayFromDevice(m_hSubOcc, m_dSubCellData.dSubOcc, 0, n*sizeof(unsigned int));
	return m_hSubOcc;
}
unsigned int* rxSPH_GPU::GetSubParticleOccScan(void)
{
	if(!m_hSubOccScan){
		m_hSubOccScan = new unsigned int[GetNumAllSubParticles()];
	}
	int n = GetNumValidSubParticles();
	CuCopyArrayFromDevice(m_hSubOccScan, m_dSubCellData.dSubOccScan, 0, n*sizeof(unsigned int));
	return m_hSubOccScan;
}


/*!
 * Et�Ɋ�Â��Ċe�T�u�p�[�e�B�N���̍����䗦�Ɣ��a���v�Z���ċl�߂�
 */
void rxSPH_GPU::CalRadiusAndRatio(bool force)
{
	if(m_bSubPacked && !force) return;

	RXREAL *dSubPos;
	if(m_bUseOpenGL){
		dSubPos = (RXREAL*)CuMapGLBufferObject(&m_pSubPosResource);
	}
	else{
		dSubPos = (RXREAL*)m_dSubPos;
	}

	RXREAL et_sub_scale = g_fCoefTurbForMesh;
	m_dSubCellData.fEtcri = g_fEtCri;

	// Et�Ɋ�Â��Ċe�T�u�p�[�e�B�N���̍����䗦�Ɣ��a���v�Z
	CuSubSetUnsortArray(m_dEt, dSubPos, m_dSubCellData, et_sub_scale);
	m_bSubPacked = true;

	if(m_bUseOpenGL){
		CuUnmapGLBufferObject(m_pSubPosResource);
	}
}

/*!
 * �T�u�p�[�e�B�N������O���b�h�̉A�֐��l���v�Z
 * @param[in] n[3] �O���b�h��
 * @param[in] minp �O���b�h�̍ŏ����W
 * @param[in] d �O���b�h��
 * @param[out] hF �A�֐��l(nx�~ny�~nz�̔z��)
 * @param[in] turb ���b�V���ό`�ɂ�闐���������邩�ǂ����̃t���O
 */
void rxSPH_GPU::CalSubImplicitFieldDevice(int n[3], Vec3 minp, Vec3 d, RXREAL *dF, int turb)
{
	// HACK:CalSubImplicitFieldDevice
	if(!m_dSubCellData.uNumParticles) return;

	RXREAL *dPos, *dSubPos;
	RXREAL et_scale = g_fWaveletScale*m_params.EffectiveRadius;
	
	if(m_bUseOpenGL){
		dPos = (RXREAL*)CuMapGLBufferObject(&m_pPosResource);
		dSubPos = (RXREAL*)CuMapGLBufferObject(&m_pSubPosResource);
	}
	else{
		dPos = (RXREAL*)m_dPos;
		dSubPos = (RXREAL*)m_dSubPos;
	}

	RXREAL et_sub_scale = g_fCoefTurbForMesh;//*g_fCoefTurbForMesh;
	m_dSubCellData.fEtcri = g_fEtCri;

	// Et�Ɋ�Â��Ċe�T�u�p�[�e�B�N���̍����䗦�Ɣ��a���v�Z
	CuSubSetUnsortArray(m_dEt, dSubPos, m_dSubCellData, et_sub_scale);
	m_bSubPacked = true;

	// �����Z���̃n�b�V�����v�Z
	CuSubCalcHash(m_dSubCellData);
	CuSubCheckRatio(m_dSubCellData);//�n�b�V���ɍő�l������
	
	// �n�b�V���Ɋ�Â��p�[�e�B�N�����\�[�g
	CuSort(m_dSubCellData.dSubGridParticleHash, m_dSubCellData.dSubSortedIndex, m_dSubCellData.uSubNumValidParticles);

	// �p�[�e�B�N���z����\�[�g���ꂽ���Ԃɕ��ёւ��C
	// �e�Z���̎n�܂�ƏI���̃C���f�b�N�X������
	CuSubReorderDataAndFindCellStart(m_dSubCellData);

	CuSubSphGridDensity(dF, m_dSubCellData, n[0], n[1], n[2], minp[0], minp[1], minp[2], d[0], d[1], d[2]);

	if(m_bUseOpenGL){
		CuUnmapGLBufferObject(m_pPosResource);
		CuUnmapGLBufferObject(m_pSubPosResource);
	}
	RXTIMER("SUB:Density");
}


void rxSPH_GPU::AddSphereField(int n[3], Vec3 minp, Vec3 d, RXREAL *hF, 
							   const Vec3 &pos, const double &eRad, const double& ratio)
{
	uint nMin[3], nMax[3];
	Vec3 relPos = pos-minp;
	double rr = eRad*eRad;
	double inv_rr = 1.0/rr;

	for(uint i = 0; i < 3; ++i){
		nMin[i] = uint((relPos[i]-eRad)/d[i])+1;
		nMin[i] = ((nMin[i] >= 0) ? nMin[i]: 0);

		nMax[i] = uint((relPos[i]+eRad)/d[i]);
		nMax[i] = ((nMax[i] <= (uint)n[i])? nMax[i]: (uint)n[i]);
	}
	
	for(uint k = nMin[2]; k <= nMax[2]; ++k){
		for(uint j = nMin[1]; j <= nMax[1]; ++j){
			for(uint i = nMin[0]; i <= nMax[0]; ++i){
				Vec3 ijkRelPos = Vec3(d[0]*(double)i, d[1]*(double)j, d[2]*(double)k);

				double dd = norm2(relPos-ijkRelPos);

				if(dd < rr){
					double tmp = 1.0-dd*inv_rr;
					hF[(i+j*(n[0]+1)+k*(n[0]+1)*(n[1]+1))] += ratio*tmp*tmp*tmp;
				}
			}
		}
	}
}

/*!
 * �T�u�p�[�e�B�N������O���b�h�̉A�֐��l���v�Z
 * @param[in] n[3] �O���b�h��
 * @param[in] minp �O���b�h�̍ŏ����W
 * @param[in] d �O���b�h��
 * @param[out] hF �A�֐��l(nx�~ny�~nz�̔z��)
 * @param[out]	scalarField3D �O���b�h�f�[�^
 * @param[in]	sfSetting	�O���b�h�̐ݒ�
 * @param[in]	et_cri		��G�l���M�[�X�y�N�g��
 */
void rxSPH_GPU::CalSubImplicitField(int n[3], Vec3 minp, Vec3 d, RXREAL *hF, const double &et_cri)
{
	if(hF){
		delete [] hF;
		hF = NULL;
	}
	hF = new RXREAL[(n[0]+1)*(n[1]+1)*(n[2]+1)];

	//�f�[�^���[�h
	unsigned int sizeSub	= m_uMaxParticles * ( calUintPow(2,m_uMaxSubLevel+1) - 1);
	unsigned int sizeSub4	= sizeSub *4;
	unsigned int sizeSubM	= m_uMaxParticles * ( calUintPow(2,m_uMaxSubLevel)   - 1);
	//unsigned int mem_size_sub	= sizeof(RXREAL)*sizeSub;
	unsigned int mem_size_subM	= sizeof(RXREAL)*sizeSubM;
	unsigned int mem_size_sub4	= sizeof(RXREAL)*sizeSub4;
	
	CuCopyArrayFromDevice(m_hSubEt,  m_dSubEt, 0, mem_size_subM);
	float *subEt	= m_hSubEt;
	float4 *subPos;
	if(m_bUseOpenGL){
		m_hSubPos = GetArrayVBO(RX_SUBPOSITION);
		subPos = (float4*)m_hSubPos;
	}
	else{
		CuCopyArrayFromDevice(m_hSubPos, m_dSubPos, 0, mem_size_sub4);
		subPos = (float4*)m_hSubPos;
	}

	//�l�ݒ�A���O�v�Z
	double eRad = m_params.EffectiveRadius;

	double *lWpoly6	= new double [m_uMaxSubLevel+1];
	double *lMass	= new double [m_uMaxSubLevel+1];
	double *leRad	= new double [m_uMaxSubLevel+1];
	double *lRatio	= new double [m_uMaxSubLevel+1];
	
	for(uint level=0; level<=m_uMaxSubLevel; level++){
		leRad[level]	= eRad * pow(2.0,-(double)level/3.0);
		lWpoly6[level]	= 315.0 / (64.0*RX_PI*pow(leRad[level], 3.0));
		lMass[level]	= m_params.Mass * pow(0.5,(double)level);
		lRatio[level]	= lMass[level] * lWpoly6[level];
	}


	const double inv_log_2_m5d9 = -9.0/(5.0*log(2.0));
	const double log_et_cri = log(et_cri);
	const double radius	= m_params.ParticleRadius;
	//const double rt = 1.1;

	//
	for(uint index=0;index<m_uNumParticles;index++){

		double et_sub = subEt[index];
		double etLevel = (log_et_cri - log(et_sub)) * inv_log_2_m5d9;

		//���x��0
		if(etLevel<=0.0){
			
			float4 posf4 = subPos[index];
			Vec3 pos = Vec3((double)posf4.x, (double)posf4.y, (double)posf4.z);
			//RXCOUT << "Pos[" << index << "] = (" << pos[0] << "," << pos[1] << "," << pos[2] << ")" << endl;
			AddSphereField(n, minp, d, hF, pos, leRad[0], lRatio[0]);
			//scalarField3D.AddColor(pos, leRad[0], lRatio[0]);
		}
		//���x��0�`Max
		else if(etLevel < (double)m_uMaxSubLevel){

			uint uLevel = (int)etLevel;
			uint dLevel = uLevel+1;

			double uRatio	= (double)dLevel - etLevel;
			double dRatio	= 1.0 - uRatio;

			//uint uSubIndexHead	= m_uNumParticles * (uintPow(2,uLevel) - 1) + uintPow(2,uLevel) * index;
			//uint uSubIndexBack	= uSubIndexHead + uintPow(2,uLevel) - 1;
			//uint dSubIndexHead	= m_uNumParticles * (uintPow(2,dLevel) - 1) + uintPow(2,dLevel) * index;
			//uint dSubIndexBack	= dSubIndexHead + uintPow(2,dLevel) - 1;
			//
			////uLevel
			//for(uint uSubIndex = uSubIndexHead; uSubIndex <= uSubIndexBack; uSubIndex++){
			//	float4 posf4 = subPos[uSubIndex];
			//	Vec3 pos = Vec3((double)posf4.x, (double)posf4.y, (double)posf4.z);
			//	scalarField3D.AddColor(pos, leRad[uLevel], uRatio*lRatio[uLevel]);
			//}
			////dLevel
			//for(uint dSubIndex = dSubIndexHead; dSubIndex <= dSubIndexBack; dSubIndex++){
			//	float4 posf4 = subPos[dSubIndex];
			//	Vec3 pos = Vec3((double)posf4.x, (double)posf4.y, (double)posf4.z);
			//	scalarField3D.AddColor(pos, leRad[dLevel], dRatio*lRatio[dLevel]);
			//}

			uint dSubIndex = m_uNumParticles * ( calUintPow(2,dLevel) - 1) + index;
			uint dBlockNum = calUintPow(2,dLevel);
			uint uSubIndex = m_uNumParticles * ( calUintPow(2,uLevel) - 1) + index;
			uint uBlockNum = calUintPow(2,uLevel);

			//uLevel
			for(uint uBlockIndex = 0 ; uBlockIndex < uBlockNum; uBlockIndex++){
				float4 posf4 = subPos[uSubIndex];
				Vec3 pos = Vec3((double)posf4.x, (double)posf4.y, (double)posf4.z);
				AddSphereField(n, minp, d, hF, pos, leRad[uLevel], uRatio*lRatio[uLevel]);
				//scalarField3D.AddColor(pos, leRad[uLevel], uRatio*lRatio[uLevel]);
				uSubIndex += m_uNumParticles;
			}
			//dLevel
			for(uint dBlockIndex = 0 ; dBlockIndex < dBlockNum; dBlockIndex++){
				float4 posf4 = subPos[dSubIndex];
				Vec3 pos = Vec3((double)posf4.x, (double)posf4.y, (double)posf4.z);
				AddSphereField(n, minp, d, hF, pos, leRad[dLevel], dRatio*lRatio[dLevel]);
				//scalarField3D.AddColor(pos, leRad[dLevel], dRatio*lRatio[dLevel]);
				uSubIndex += m_uNumParticles;
			}

		}
		//���x��Max
		else{
			//uint subIndexHead	= m_uNumParticles * (uintPow(2,m_uMaxSubLevel) - 1) + uintPow(2,m_uMaxSubLevel) * index;
			//uint subIndexBack	= subIndexHead + uintPow(2,m_uMaxSubLevel) - 1;

			//for(uint subIndex = subIndexHead; subIndex <= subIndexBack; subIndex++){
			//	float4 posf4 = subPos[subIndex];
			//	Vec3 pos = Vec3((double)posf4.x, (double)posf4.y, (double)posf4.z);
			//	scalarField3D.AddColor(pos, leRad[m_uMaxSubLevel], lRatio[m_uMaxSubLevel]);
			//}

			uint mSubIndex = m_uNumParticles * ( calUintPow(2,m_uMaxSubLevel) - 1) + index;
			uint mBlockNum = calUintPow(2,m_uMaxSubLevel);

			for(uint mBlockIndex = 0 ; mBlockIndex < mBlockNum; mBlockIndex++){
				float4 posf4 = subPos[mSubIndex];
				Vec3 pos = Vec3((double)posf4.x, (double)posf4.y, (double)posf4.z);
				AddSphereField(n, minp, d, hF, pos, leRad[m_uMaxSubLevel], lRatio[m_uMaxSubLevel]);
				//scalarField3D.AddColor(pos, leRad[m_uMaxSubLevel], lRatio[m_uMaxSubLevel]);
				mSubIndex += m_uNumParticles;
			}
			
		}
	}

	delete[] lWpoly6;
	delete[] lMass;
	delete[] leRad;
	delete[] lRatio;

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
void rxSPH_GPU::DetectSurfaceParticles(void)
{
	RXREAL h = m_params.EffectiveRadius;
	RXREAL r = m_params.ParticleRadius;

	//m_hPos = GetArrayVBO(RX_POSITION);

	// �ߖT���q�T��
	SetParticlesToCell();

	for(uint i = 0; i < m_uNumParticles; ++i){
		double d = CalDistToNormalizedMassCenter(i);
		int nn_num = (int)m_vNeighs[i].size();
//		cout << "nn_num = " << nn_num << endl;
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
double rxSPH_GPU::CalDistToNormalizedMassCenter(const int i)
{
	Vec3 sum_pos(0.0);
	double sum_mass = 0.0;

	Vec3 pos0 = Vec3(m_hPos[DIM*i+0], m_hPos[DIM*i+1], m_hPos[DIM*i+2]);

	for(vector<rxNeigh>::iterator itr = m_vNeighs[i].begin() ; itr != m_vNeighs[i].end(); ++itr){
		int j = itr->Idx;
		if(j < 0 && i == j) continue;

		//RXREAL r = sqrt(itr->Dist2);
		Vec3 pos1 = Vec3(m_hPos[DIM*j+0], m_hPos[DIM*j+1], m_hPos[DIM*j+2]);

		sum_pos  += (pos0-pos1)*m_params.Mass;
		sum_mass += m_params.Mass;
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
uint* rxSPH_GPU::GetArraySurf(void)
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
int rxSPH_GPU::GetSurfaceParticles(const Vec3 pos, RXREAL h, vector<rxSurfaceParticle> &sps)
{
	int num = 0;
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
void rxSPH_GPU::CalNormal(void)
{
}

/*!
 * ���x���z����@�����v�Z
 */
void rxSPH_GPU::CalNormalFromDensity(void)
{
	// �@���v�Z
	CuSphNormal(m_dNrm, m_dDens, m_dCellData);

	m_hNrm = GetArrayVBO(RX_NORMAL);
}



//-----------------------------------------------------------------------------
// Anisotropic Kernel
//-----------------------------------------------------------------------------
/*!
 * Anisotropic kernel�̌v�Z
 *  - J. Yu and G. Turk, Reconstructing Surfaces of Particle-Based Fluids Using Anisotropic Kernels, SCA2010. 
 */
void rxSPH_GPU::CalAnisotropicKernel(void)
{
	if(m_uNumParticles == 0) return;

	RXREAL r0 = m_params.Density;
	RXREAL h0 = m_params.EffectiveRadius;
	RXREAL h = 2.0*h0;

	//SetParticlesToCell();

	UpdateParams();
	m_params.NumBodies = m_uNumParticles;
	m_dCellData.uNumParticles = m_uNumParticles;

	RXREAL *dPos;
	if(m_bUseOpenGL){
		dPos = (RXREAL*)CuMapGLBufferObject(&m_pPosResource);
	}
	else{
		dPos = (RXREAL*)m_dPos;
	}

	// �ߖT�T���������p�O���b�h�f�[�^�̍쐬
	// �����Z���̃n�b�V�����v�Z
	CuCalcHash(m_dCellData.dGridParticleHash, m_dCellData.dSortedIndex, dPos, m_dAttr, m_uNumParticles);

	// �n�b�V���Ɋ�Â��p�[�e�B�N�����\�[�g
	CuSort(m_dCellData.dGridParticleHash, m_dCellData.dSortedIndex, m_uNumParticles);

	// �p�[�e�B�N���z����\�[�g���ꂽ���Ԃɕ��ёւ��C
	// �e�Z���̎n�܂�ƏI���̃C���f�b�N�X������
	CuReorderDataAndFindCellStart(m_dCellData, dPos, m_dVel);

	if(m_bUseOpenGL){
		CuUnmapGLBufferObject(m_pPosResource);
	}


	// �������J�[�l�����S(��6)�Əd�ݕt�����ψʒu(��10)�̌v�Z
	RXREAL lambda = 0.9;
	CuSphCalUpdatedPosition(m_dUpPos, m_dPosW, lambda, h, m_dCellData);

	// �������ʒu�ŋߖT�T���Z�����X�V
	setObjectToCell(m_dUpPos);

	// �������ʒu�ł̏d�ݕt�����ψʒu�̍Čv�Z��covariance matrix�̌v�Z(��9)
	CuSphCalCovarianceMatrix(m_dPosW, m_dCMatrix, h, m_dCellData);

	RXTIMER("aniso : cmatrix");

	// ���ْl�����ɂ����ْl,���كx�N�g�����v�Z
	CuSphSVDecomposition(m_dCMatrix, m_dPosW, m_dEigen, m_dRMatrix, m_uNumParticles);

	RXTIMER("aniso : svd");
	
	// �ό`�s��G�̌v�Z
	CuSphCalTransformMatrix(m_dEigen, m_dRMatrix, m_dG, m_uNumParticles);

	//CuCopyArrayFromDevice(m_hG, m_dG, 0, m_uNumParticles*9*sizeof(RXREAL));

	CuCopyArrayFromDevice(m_hEigen, m_dEigen, 0, m_uNumParticles*3*sizeof(RXREAL));
	m_fEigenMax = 0.0;
	for(uint i = 0; i < m_uNumParticles; ++i){
		for(int j = 0; j < 3; ++j){
			if(m_hEigen[3*i+j] > m_fEigenMax) m_fEigenMax = m_hEigen[3*i+j];
		}
	}
	//RXCOUT << "eigen_max = " << m_fEigenMax << endl;

	m_bCalAnisotropic = true;

	RXTIMER("aniso : G");

}


/*!
 * �J���[�l�pVBO�̕ҏW
 * @param[in] type �F�̃x�[�X�Ƃ��镨���l
 */
void rxSPH_GPU::SetColorVBO(int type)
{
	// MRK:SetColorVBO

	switch(type){
	case RX_DENSITY:
		CuCopyArrayFromDevice(m_hDens, m_dDens, 0, sizeof(RXREAL)*m_uNumParticles);
		SetColorVBOFromArray(m_hDens, 1, false, m_params.Density*3);
		break;

	case RX_PRESSURE:
		CuCopyArrayFromDevice(m_hPres, m_dPres, 0, sizeof(RXREAL)*m_uNumParticles);
		SetColorVBOFromArray(m_hPres, 1);
		break;

	case RX_ENERGY_SPECTRUM:
		CuCopyArrayFromDevice(m_hEt, m_dEt, 0, sizeof(RXREAL)*m_uNumParticles);
		SetColorVBOFromArray(m_hEt, 1);
		break;

	case RX_CONSTANT:
		if(m_bUseOpenGL){
			// �J���[�o�b�t�@�ɒl��ݒ�
			glBindBufferARB(GL_ARRAY_BUFFER, m_colorVBO);
			RXREAL *data = (RXREAL*)glMapBufferARB(GL_ARRAY_BUFFER, GL_WRITE_ONLY);
			RXREAL *ptr = data;
			for(uint i = 0; i < m_uNumParticles; ++i){
				RXREAL t = i/(RXREAL)m_uNumParticles;
				RX_COLOR_RAMP<RXREAL>(m_hAttr[i]/3.0, ptr);
				ptr += 3;

				//*ptr++ = 0.15f;
				//*ptr++ = 0.15f;
				//*ptr++ = 0.95f;
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


void rxSPH_GPU::CalMaxDensity(int k)
{
	CuCopyArrayFromDevice(m_hDens, m_dDens, 0, sizeof(RXREAL)*2*m_uNumParticles);

	RXREAL max_dens = 0.0;
	for(uint i = 0; i < m_uNumParticles; ++i){
		if(m_hDens[2*i+k] > max_dens) max_dens = m_hDens[2*i+k];
	}

	RXCOUT << "Density  : " << max_dens << endl;
}


/*!
 * �O���b�h�n�b�V���l�̌v�Z
 * @param[in] x,y,z �O���b�h�ʒu
 * @return �O���b�h�n�b�V���l
 */
uint rxSPH_GPU::calGridHash(uint x, uint y, uint z)
{
	x = (x < 0 ? 0 : (x >= m_params.GridSize.x ? m_params.GridSize.x-1 : x));
	y = (y < 0 ? 0 : (y >= m_params.GridSize.y ? m_params.GridSize.y-1 : y));
	z = (z < 0 ? 0 : (z >= m_params.GridSize.z ? m_params.GridSize.z-1 : z));
	return z*m_params.GridSize.y*m_params.GridSize.x+y*m_params.GridSize.x+x;
}
/*!
 * �O���b�h�n�b�V���l�̌v�Z
 * @param[in] pos �p�[�e�B�N�����W
 * @return �O���b�h�n�b�V���l
 */
uint rxSPH_GPU::calGridHash(Vec3 pos)
{
	pos -= m_v3EnvMin;

	// �����Z���C���f�b�N�X�̎Z�o
	uint x = pos[0]/m_params.CellWidth.x;
	uint y = pos[1]/m_params.CellWidth.y;
	uint z = pos[2]/m_params.CellWidth.z;
	return calGridHash(x, y, z);
}

/*!
 * ��ԕ����@�̏���
 * @param[out] cell �����O���b�h�f�[�^
 * @param[out] gridsize �e���̃O���b�h��
 * @param[out] cell_width �Z����
 * @param[in] vMin ���̍ŏ����W
 * @param[in] vMax ���̍ő���W
 * @param[in] h �e�����a
 */
void rxSPH_GPU::setupCells(rxParticleCell &cell, uint3 &gridsize, double &cell_width, Vec3 vMin, Vec3 vMax, double h)
{
	if(h < RX_EPS) return;

	Vec3 world_size = vMax-vMin;
	Vec3 world_origin = vMin;

	double max_axis = RXFunc::Max3(world_size);

	int d = (int)(log(max_axis/h)/log(2.0)+0.5);
	int n = (int)(pow(2.0, (double)d)+0.5);
	cell_width = max_axis/n;

	//d = (int)(log(world_size[0]/cell_width)/log(2.0)+0.5);
	//gridsize.x = (int)(pow(2.0, (double)d)+0.5);
	//d = (int)(log(world_size[1]/cell_width)/log(2.0)+0.5);
	//gridsize.y = (int)(pow(2.0, (double)d)+0.5);;
	//d = (int)(log(world_size[2]/cell_width)/log(2.0)+0.5);
	//gridsize.z = (int)(pow(2.0, (double)d)+0.5);;

	gridsize.x = (int)(world_size[0]/cell_width)+1;
	gridsize.y = (int)(world_size[1]/cell_width)+1;
	gridsize.z = (int)(world_size[2]/cell_width)+1;

	cell.uNumCells = gridsize.x*gridsize.y*gridsize.z;
}



/*!
 * �S�p�[�e�B�N���𕪊��Z���Ɋi�[
 *  - �e�p�[�e�B�N���̑�����O���b�h�n�b�V�����v�Z���Ċi�[����
 * @param[in] p �S�p�[�e�B�N���̍��W���i�[�����z��
 */
void rxSPH_GPU::setObjectToCell(RXREAL *p)
{
	// �ߖT�T���������p�O���b�h�f�[�^�̍쐬
	// �����Z���̃n�b�V�����v�Z
	CuCalcHash(m_dCellData.dGridParticleHash, m_dCellData.dSortedIndex, p, m_dAttr, m_uNumParticles);

	// �n�b�V���Ɋ�Â��p�[�e�B�N�����\�[�g
	CuSort(m_dCellData.dGridParticleHash, m_dCellData.dSortedIndex, m_uNumParticles);

	// �p�[�e�B�N���z����\�[�g���ꂽ���Ԃɕ��ёւ��C
	// �e�Z���̎n�܂�ƏI���̃C���f�b�N�X������
	CuReorderDataAndFindCellStart(m_dCellData, p, m_dVel);
}

/*!
 * �S�p�[�e�B�N���𕪊��Z���Ɋi�[
 */
void rxSPH_GPU::SetParticlesToCell(RXREAL *prts, int n, RXREAL h)
{
}
void rxSPH_GPU::SetParticlesToCell(void)
{
	SetParticlesToCell(m_hPos, m_uNumParticles, m_params.EffectiveRadius);
}

/*!
 * �|���S���𕪊��Z���Ɋi�[
 * @param[in] vrts �|���S�����_
 * @param[in] nv ���_��
 * @param[in] tris ���b�V��
 * @param[in] nt ���b�V����
 */
void rxSPH_GPU::setPolysToCell(RXREAL *vrts, int nv, int* tris, int nt)
{
	// MRK:setPolysToCell
	uint *hPolyCellStart = new uint[m_dCellData.uNumCells];
	uint *hPolyCellEnd = new uint[m_dCellData.uNumCells];

	int mem_size2 = m_dCellData.uNumCells*sizeof(uint);
	memset(hPolyCellStart, 0xffffffff, mem_size2);
	memset(hPolyCellEnd,   0,          mem_size2);

	int num_hash = 0;

	// �e�p�[�e�B�N���̃O���b�h�n�b�V���̌v�Z
	vector<uint> tri_hash, tri_idx;
	vector<Vec3> tri_vrts, tri_vrts_c;
	tri_vrts.resize(3);
	tri_vrts_c.resize(3);
	for(int i = 0; i < nt; i++){
		for(int j = 0; j < 3; ++j){
			Vec3 pos;
			pos[0] = vrts[3*tris[3*i+j]+0];
			pos[1] = vrts[3*tris[3*i+j]+1];
			pos[2] = vrts[3*tris[3*i+j]+2];
			tri_vrts[j] = pos;
		}

		Vec3 nrm = Unit(cross(tri_vrts[1]-tri_vrts[0], tri_vrts[2]-tri_vrts[0]));

		// �|���S����BBox
		Vec3 bmin, bmax;
		bmin = tri_vrts[0];
		bmax = tri_vrts[0];
		for(int j = 1; j < 3; ++j){
			for(int k = 0; k < 3; ++k){
				if(tri_vrts[j][k] < bmin[k]) bmin[k] = tri_vrts[j][k];
				if(tri_vrts[j][k] > bmax[k]) bmax[k] = tri_vrts[j][k];
			}
		}

		// BBox�Əd�Ȃ�Z��
		bmin -= m_v3EnvMin;
		bmax -= m_v3EnvMin;

		// �����Z���C���f�b�N�X�̎Z�o
		int bmin_gidx[3], bmax_gidx[3];
		bmin_gidx[0] = bmin[0]/m_params.CellWidth.x;
		bmin_gidx[1] = bmin[1]/m_params.CellWidth.y;
		bmin_gidx[2] = bmin[2]/m_params.CellWidth.z;
		bmax_gidx[0] = bmax[0]/m_params.CellWidth.x;
		bmax_gidx[1] = bmax[1]/m_params.CellWidth.y;
		bmax_gidx[2] = bmax[2]/m_params.CellWidth.z;

		bmin_gidx[0] = RX_CLAMP(bmin_gidx[0], 0, (int)m_params.GridSize.x-1);
		bmin_gidx[1] = RX_CLAMP(bmin_gidx[1], 0, (int)m_params.GridSize.y-1);
		bmin_gidx[2] = RX_CLAMP(bmin_gidx[2], 0, (int)m_params.GridSize.z-1);
		bmax_gidx[0] = RX_CLAMP(bmax_gidx[0], 0, (int)m_params.GridSize.x-1);
		bmax_gidx[1] = RX_CLAMP(bmax_gidx[1], 0, (int)m_params.GridSize.y-1);
		bmax_gidx[2] = RX_CLAMP(bmax_gidx[2], 0, (int)m_params.GridSize.z-1);

		// �e�Z���Ƀ|���S�����܂܂�邩���`�F�b�N
		Vec3 len = Vec3(m_params.CellWidth.x, m_params.CellWidth.y, m_params.CellWidth.z);
		Vec3 cen(0.0);
		for(int x = bmin_gidx[0]; x <= bmax_gidx[0]; ++x){
			for(int y = bmin_gidx[1]; y <= bmax_gidx[1]; ++y){
				for(int z = bmin_gidx[2]; z <= bmax_gidx[2]; ++z){
					cen = m_v3EnvMin+Vec3(x+0.5, y+0.5, z+0.5)*len;

					for(int j = 0; j < 3; ++j){
						tri_vrts_c[j] = (tri_vrts[j]-cen)/len;
					}

					if(RXFunc::polygon_intersects_cube(tri_vrts_c, nrm, 0)){
						// �n�b�V���l�v�Z
						uint hash = calGridHash(x, y, z);

						tri_idx.push_back((uint)i);
						tri_hash.push_back(hash);

						num_hash++;
					}
				}
			}
		}
	}

	RXCOUT << "polygon hash : " << num_hash << endl;

	m_dCellData.uNumPolyHash = (uint)num_hash;

	if(num_hash){
		int mem_size1 = m_dCellData.uNumPolyHash*sizeof(uint);
		uint *hSortedPolyIdx = new uint[m_dCellData.uNumPolyHash];
		uint *hGridPolyHash  = new uint[m_dCellData.uNumPolyHash];
		memcpy(hSortedPolyIdx, &tri_idx[0], mem_size1);
		memcpy(hGridPolyHash, &tri_hash[0], mem_size1);

		// �O���b�h�n�b�V���Ń\�[�g
		if(m_dCellData.dSortedPolyIdx) CuFreeArray(m_dCellData.dSortedPolyIdx);
		if(m_dCellData.dGridPolyHash) CuFreeArray(m_dCellData.dGridPolyHash);
		CuAllocateArray((void**)&m_dCellData.dSortedPolyIdx, mem_size1);
		CuAllocateArray((void**)&m_dCellData.dGridPolyHash,  mem_size1);

		CuCopyArrayToDevice(m_dCellData.dSortedPolyIdx, &tri_idx[0],  0, mem_size1);
		CuCopyArrayToDevice(m_dCellData.dGridPolyHash,  &tri_hash[0], 0, mem_size1);

		CuSort(m_dCellData.dGridPolyHash, m_dCellData.dSortedPolyIdx, m_dCellData.uNumPolyHash);

		CuCopyArrayFromDevice(hSortedPolyIdx, m_dCellData.dSortedPolyIdx, 0, mem_size1);
		CuCopyArrayFromDevice(hGridPolyHash,  m_dCellData.dGridPolyHash,  0, mem_size1);

		// �p�[�e�B�N���z����\�[�g���ꂽ���Ԃɕ��ёւ��C
		// �e�Z���̎n�܂�ƏI���̃C���f�b�N�X������
		for(uint i = 0; i < m_dCellData.uNumPolyHash; ++i){
			uint hash = hGridPolyHash[i];

			if(i == 0){
				hPolyCellStart[hash] = i;
			}
			else{
				uint prev_hash = hGridPolyHash[i-1];

				if(i == 0 || hash != prev_hash){
					hPolyCellStart[hash] = i;
					if(i > 0){
						hPolyCellEnd[prev_hash] = i;
					}
				}

				if(i == m_uNumParticles-1){
					hPolyCellEnd[hash] = i+1;
				}
			}
		}

		if(m_dCellData.dPolyCellStart) CuFreeArray(m_dCellData.dPolyCellStart);
		if(m_dCellData.dPolyCellEnd) CuFreeArray(m_dCellData.dPolyCellEnd);

		CuAllocateArray((void**)&m_dCellData.dPolyCellStart, mem_size2);
		CuAllocateArray((void**)&m_dCellData.dPolyCellEnd,   mem_size2);

		CuCopyArrayToDevice(m_dCellData.dPolyCellStart,  hPolyCellStart, 0, mem_size2);
		CuCopyArrayToDevice(m_dCellData.dPolyCellEnd,    hPolyCellEnd,   0, mem_size2);

		delete [] hGridPolyHash;
		delete [] hSortedPolyIdx;
	}

	delete [] hPolyCellStart;
	delete [] hPolyCellEnd;
}

/*!
 * �|���S���𕪊��Z���Ɋi�[
 */
void rxSPH_GPU::SetPolygonsToCell(void)
{
	setPolysToCell(&m_vVrts[0], m_iNumVrts, &m_vTris[0], m_iNumTris);
}


/*!
 * �T���p�Z���̕`��
 * @param[in] i,j,k �O���b�h��̃C���f�b�N�X
 */
void rxSPH_GPU::DrawCell(int i, int j, int k)
{
	glPushMatrix();
	glTranslated(m_v3EnvMin[0], m_v3EnvMin[1], m_v3EnvMin[2]);
	glTranslatef((i+0.5)*m_params.CellWidth.x, (j+0.5)*m_params.CellWidth.y, (k+0.5)*m_params.CellWidth.z);
	glutWireCube(m_params.CellWidth.x);
	glPopMatrix();
}


/*!
 * �T���p�O���b�h�̕`��
 * @param[in] col �p�[�e�B�N�����܂܂��Z���̐F
 * @param[in] col2 �|���S�����܂܂��Z���̐F
 * @param[in] sel �����_���ɑI�����ꂽ�Z���̂ݕ`��(1�ŐV�����Z����I���C2�ł��łɑI������Ă���Z����`��C0�ł��ׂẴZ����`��)
 */
void rxSPH_GPU::DrawCells(Vec3 col, Vec3 col2, int sel)
{
	glPushMatrix();

	uint *hCellStart     = new uint[m_dCellData.uNumCells];
	uint *hPolyCellStart = new uint[m_dCellData.uNumCells];
	CuCopyArrayFromDevice(hCellStart,     m_dCellData.dCellStart,     0, m_dCellData.uNumCells*sizeof(uint));
	CuCopyArrayFromDevice(hPolyCellStart, m_dCellData.dPolyCellStart, 0, m_dCellData.uNumCells*sizeof(uint));

	if(sel){
		uint *hCellEnd     = new uint[m_dCellData.uNumCells];
		uint *hSortedIndex = new uint[m_uNumParticles];
		CuCopyArrayFromDevice(hCellEnd,     m_dCellData.dCellEnd,     0, m_dCellData.uNumCells*sizeof(uint));
		CuCopyArrayFromDevice(hSortedIndex, m_dCellData.dSortedIndex, 0, m_uNumParticles*sizeof(uint));

		// �����_���ɑI�񂾃Z���Ƃ��̒��̃p�[�e�B�N���̂ݕ`��
		static int grid_hash = 0;
		static uint start_index = 0xffffffff;
		if(sel == 1){
			do{
				grid_hash = RXFunc::Nrand(m_dCellData.uNumCells-1);
				start_index = hCellStart[grid_hash];
			}while(start_index == 0xffffffff);
		}

		uint w = grid_hash%(m_params.GridSize.x*m_params.GridSize.y);
		DrawCell(w%m_params.GridSize.x, w/m_params.GridSize.x, grid_hash/(m_params.GridSize.x*m_params.GridSize.y));

		glColor3d(1.0, 0.0, 0.0);
		glPointSize(10.0);
		glBegin(GL_POINTS);

		int c = 0;
		uint end_index = hCellEnd[grid_hash];
		for(uint j = start_index; j < end_index; ++j){
			uint idx = hSortedIndex[j];
			Vec3 pos;
			pos[0] = m_hPos[4*idx+0];
			pos[1] = m_hPos[4*idx+1];
			pos[2] = m_hPos[4*idx+2];
			
			glVertex3dv(pos);

			c++;
		}
		glEnd();
		cout << "cell(" << grid_hash << ") : " << c << endl;

		delete [] hCellEnd;
		delete [] hSortedIndex;
	}
	else{
		int cnt = 0;
		// �p�[�e�B�N�� or �|���S�����܂ޑS�Z���̕`��
		RXFOR3(0, (int)m_params.GridSize.x, 0, (int)m_params.GridSize.y, 0, (int)m_params.GridSize.z){
			bool disp = false;
			uint grid_hash = calGridHash(i, j, k);
			uint start_index = hCellStart[grid_hash];
			uint start_index_poly = 0xffffffff;
		
			if(m_dCellData.uNumPolyHash) start_index_poly = hPolyCellStart[grid_hash];

			if(start_index != 0xffffffff){
				glColor3dv(col2.data);
				disp = true;
			}
			if(start_index_poly != 0xffffffff){
				glColor3dv(col.data);
				disp = true;
			}

			cnt++;

			if(disp){
				DrawCell(i, j, k);
			}
		}

		cout << cnt << endl;
	}

	delete [] hCellStart;
	delete [] hPolyCellStart;

	glPopMatrix();
}

/*!
 * �ő̏�Q���̕`��
 */
void rxSPH_GPU::DrawObstacles(void)
{
#if MAX_BOX_NUM
	for(int i = 0; i < (int)m_params.BoxNum; ++i){
		if(!m_params.BoxFlg[i]) continue;

		Vec3 bcen = Vec3(m_params.BoxCen[i].x, m_params.BoxCen[i].y, m_params.BoxCen[i].z);
		Vec3 bext = Vec3(m_params.BoxExt[i].x, m_params.BoxExt[i].y, m_params.BoxExt[i].z);
		float bmat[16];
		GetGLMatrix(m_params.BoxRot[i], bmat);

		glPushMatrix();
		glTranslated(bcen[0], bcen[1], bcen[2]);
		glMultMatrixf(bmat);
		glScalef(2.0*bext[0], 2.0*bext[1], 2.0*bext[2]);
		//glRotated(brot, 0, 0, 1);
		glutWireCube(1.0);
		glPopMatrix();
	}
#endif

#if MAX_SPHERE_NUM
	for(int i = 0; i < (int)m_params.SphereNum; ++i){
		if(!m_params.SphereFlg[i]) continue;

		Vec3 scen  = Vec3(m_params.SphereCen[i].x, m_params.SphereCen[i].y, m_params.SphereCen[i].z);
		float srad = m_params.SphereRad[i];

		glPushMatrix();
		glTranslated(scen[0], scen[1], scen[2]);
		glutWireSphere(srad, 32, 32);
		glPopMatrix();
	}
#endif

}





/*!
 * �O�p�`�|���S���ɂ���Q��
 * @param[in] vrts ���_
 * @param[in] tris ���b�V��
 */
void rxSPH_GPU::SetPolygonObstacle(const vector<Vec3> &vrts, const vector<Vec3> &nrms, const vector< vector<int> > &tris, Vec3 vel)
{
	int vn = (int)vrts.size();
	int n = (int)tris.size();
	m_iNumVrts += vn;
	m_iNumTris += n;

	//if(nrms.empty()){
		for(int i = 0; i < vn; ++i){
			for(int j = 0; j < 3; ++j){
				m_vVrts.push_back(vrts[i][j]);
			}
		}
	//}
	//else{
	//	for(int i = 0; i < vn; ++i){
	//		Vec3 v = vrts[i]+nrms[i]*m_fParticleRadius;
	//		for(int j = 0; j < 3; ++j){
	//			m_vVrts.push_back(v[j]);
	//		}
	//	}
	//}

	for(int i = 0; i < n; ++i){
		for(int j = 0; j < 3; ++j){
			m_vTris.push_back(tris[i][j]);
		}
	}
	//int p = m_params.PolyNum;
	//m_params.PolyVel[p] = MAKE_FLOAT3(vel[0], vel[1], vel[2]);
	//m_params.PolyNum++;

	// GPU�������̊m�ۂƓ]��
	if(m_dVrts) CuFreeArray(m_dVrts);
	if(m_dTris) CuFreeArray(m_dTris);
	m_dVrts = 0;
	m_dTris = 0;

	CuAllocateArray((void**)&m_dVrts, m_iNumVrts*3*sizeof(RXREAL));
	CuAllocateArray((void**)&m_dTris, m_iNumTris*3*sizeof(int));

	CuCopyArrayToDevice(m_dVrts, &m_vVrts[0], 0, m_iNumVrts*3*sizeof(RXREAL));
	CuCopyArrayToDevice(m_dTris, &m_vTris[0], 0, m_iNumTris*3*sizeof(int));

	RXCOUT << "the number of triangles : " << m_iNumTris << endl;

	//SetParticlesToCell();
	setPolysToCell(&m_vVrts[0], m_iNumVrts, &m_vTris[0], m_iNumTris);

}

/*!
 * �{�b�N�X�^��Q��
 * @param[in] cen �{�b�N�X���S���W
 * @param[in] ext �{�b�N�X�̑傫��(�ӂ̒�����1/2)
 * @param[in] ang �{�b�N�X�̊p�x(�I�C���[�p)
 * @param[in] flg �L��/�����t���O
 */
void rxSPH_GPU::SetBoxObstacle(Vec3 cen, Vec3 ext, Vec3 ang, Vec3 vel, int flg)
{
	// ��Q��
	int b = m_params.BoxNum;

	if(b < MAX_BOX_NUM){
		m_params.BoxCen[b] = MAKE_FLOAT3(cen[0], cen[1], cen[2]);
		m_params.BoxExt[b] = MAKE_FLOAT3(ext[0], ext[1], ext[2]);
		m_params.BoxRot[b] = EulerToMatrix(ang[0], ang[1], ang[2]);
		m_params.BoxInvRot[b] = Inverse(m_params.BoxRot[b]);
		m_params.BoxFlg[b] = flg;
		b++;

		m_params.BoxNum = b;
	}
}

/*!
 * ���^��Q��
 * @param[in] cen ���̒��S���W
 * @param[in] rad ���̂̔��a
 * @param[in] flg �L��/�����t���O
 */
void rxSPH_GPU::SetSphereObstacle(Vec3 cen, double rad, Vec3 vel, int flg)
{
	// ��Q��
	int b = m_params.SphereNum;

	if(b < MAX_SPHERE_NUM){
		m_params.SphereCen[b] = MAKE_FLOAT3(cen[0], cen[1], cen[2]);
		m_params.SphereRad[b] = (RXREAL)rad;
		m_params.SphereFlg[b] = flg;
		b++;

		m_params.SphereNum = b;
	}
}

/*!
 * ���^��Q���𓮂���
 * @param[in] b ���̔ԍ�
 * @param[in] disp �ړ���
 */
void rxSPH_GPU::MoveSphereObstacle(int b, Vec3 disp)
{
	// ��Q��
	if(b >= (int)m_params.SphereNum) return;

	if(b < MAX_SPHERE_NUM){
		m_params.SphereCen[b].x += disp[0];
		m_params.SphereCen[b].y += disp[1];
		m_params.SphereCen[b].z += disp[2];
	}
}

/*!
 * ���^��Q���̈ʒu���擾
 * @param[in] b ���̔ԍ�
 */
Vec3 rxSPH_GPU::GetSphereObstaclePos(int b)
{
	if(b >= (int)m_params.SphereNum || (int)m_params.SphereNum == 0) return Vec3(0.0);
	if(b < 0) b = (int)m_params.SphereNum-1;
	return Vec3(m_params.SphereCen[b].x, m_params.SphereCen[b].y, m_params.SphereCen[b].z);
}

/*!
 * VBO,�f�o�C�X����������z�X�g�������փf�[�^��]���C�擾
 * @param[in] type �f�[�^�̎��
 * @param[in] d2h �f�o�C�X����z�X�g�ւ̃f�[�^�]���t���O
 * @return �z�X�g��������̃f�[�^
 */
RXREAL* rxSPH_GPU::GetArrayVBO(rxParticleArray type, bool d2h, int num)
{
	assert(m_bInitialized);

	if(num == -1) num = m_uNumParticles;
 
	RXREAL* hdata = 0;
	RXREAL* ddata = 0;

	cudaGraphicsResource **graphics_resource = 0;
	int d = DIM;

	switch(type){
	default:
	case RX_POSITION:
		hdata = m_hPos;
		ddata = m_dPos;
		graphics_resource = &m_pPosResource;
		break;

	case RX_VELOCITY:
		hdata = m_hVel;
		ddata = m_dVel;
		break;

	case RX_DENSITY:
		hdata = m_hDens;
		ddata = m_dDens;
		d = 1;
		break;

	case RX_PRESSURE:
		hdata = m_hPres;
		ddata = m_dPres;
		d = 1;
		break;

	case RX_NORMAL:
		hdata = m_hNrm;
		ddata = m_dNrm;
		break;

	case RX_FORCE:
		hdata = m_hFrc;
		ddata = m_dFrc;
		break;

	case RX_BOUNDARY_PARTICLE:
		hdata = m_hPosB;
		ddata = m_dPosB;
		break;

	case RX_BOUNDARY_PARTICLE_VOL:
		hdata = m_hVolB;
		ddata = m_dVolB;
		d = 1;
		break;

	case RX_TURB_VELOCITY:
		hdata = m_hTurb;
		ddata = m_dTurb;
		d2h = true;
		break;

	case RX_UPDATED_POSITION:
		hdata = m_hUpPos;
		ddata = m_dUpPos;
		d2h = true;
		break;

	case RX_EIGEN_VALUE:
		hdata = m_hEigen;
		ddata = m_dEigen;
		d = 3;
		d2h = true;
		break;

	case RX_ROTATION_MATRIX:
		hdata = m_hRMatrix;
		ddata = m_dRMatrix;
		d = 9;
		d2h = true;
		break;

	case RX_TRANSFORMATION:
		hdata = m_hG;
		ddata = m_dG;
		d = 9;
		d2h = true;
		break;
	}

	if(d2h) CuCopyArrayFromDevice(hdata, ddata, graphics_resource, num*d*sizeof(RXREAL));

	return hdata;
}

/*!
 * �z�X�g����������VBO,�f�o�C�X�������փf�[�^��]��
 * @param[in] type �f�[�^�̎��
 * @param[in] data �z�X�g��������̃f�[�^
 * @param[in] start �f�[�^�̊J�n�C���f�b�N�X
 * @param[in] count �ǉ���
 */
void rxSPH_GPU::SetArrayVBO(rxParticleArray type, const RXREAL* data, int start, int count)
{
	assert(m_bInitialized);
 
	switch(type){
	case RX_POSITION:
		{
			if(m_bUseOpenGL){
				CuUnregisterGLBufferObject(m_pPosResource);
				glBindBuffer(GL_ARRAY_BUFFER, m_posVBO);
				glBufferSubData(GL_ARRAY_BUFFER, start*4*sizeof(RXREAL), count*4*sizeof(RXREAL), data);
				glBindBuffer(GL_ARRAY_BUFFER, 0);
				CuRegisterGLBufferObject(m_posVBO, &m_pPosResource);
			}
		}
		break;

	case RX_VELOCITY:
		CuCopyArrayToDevice(m_dVel, data, start*DIM*sizeof(RXREAL), count*DIM*sizeof(RXREAL));
		CuCopyArrayToDevice(m_dTurb, data, start*DIM*sizeof(RXREAL), count*DIM*sizeof(RXREAL));
		break;

	case RX_NORMAL:
		CuCopyArrayToDevice(m_dNrm, data, start*DIM*sizeof(RXREAL), count*DIM*sizeof(RXREAL));
		break;

	case RX_BOUNDARY_PARTICLE:
		CuCopyArrayToDevice(m_dPosB, data, start*DIM*sizeof(RXREAL), count*DIM*sizeof(RXREAL));
		break;

	default:
		break;
	}	   
}
void rxSPH_GPU::SetArrayVBO(rxParticleArray type, const int* data, int start, int count)
{
	assert(m_bInitialized);
 
	switch(type){
	case RX_ATTRIBUTE:
		CuCopyArrayToDevice(m_dAttr, data, start*sizeof(int), count*sizeof(int));
		break;
	default:
		break;
	}	   
}



//-----------------------------------------------------------------------------
// MARK:�A�֐��l
//-----------------------------------------------------------------------------
double rxSPH_GPU::GetImplicit(double x, double y, double z)
{
	return 0.0;
}

/*!
 * �p�[�e�B�N������O���b�h�̉A�֐��l���v�Z
 * @param[in] pnx,pny,pnz �O���b�h���̎w�� nx=2^pnx
 * @param[in] minp �O���b�h�̍ŏ����W
 * @param[in] d �O���b�h��
 * @param[out] hF �A�֐��l(nx�~ny�~nz�̔z��)
 */
void rxSPH_GPU::CalImplicitField(int n[3], Vec3 minp, Vec3 d, RXREAL *hF)
{
	unsigned int memSize = sizeof(RXREAL)*n[0]*n[1]*n[2];
	//RXCOUT << memSize/sizeof(RXREAL) << endl;
	
	RXREAL *dF = 0;
	CuAllocateArray((void**)&dF, memSize);

	CuSphGridDensity(dF, m_dCellData, n[0], n[1], n[2], minp[0], minp[1], minp[2], d[0], d[1], d[2]);

	CuCopyArrayFromDevice(hF, dF, 0, memSize);

	if(dF) CuFreeArray(dF);
}


/*!
 * �p�[�e�B�N������O���b�h�̉A�֐��l���v�Z
 * @param[in] n[3] �O���b�h��
 * @param[in] minp �O���b�h�̍ŏ����W
 * @param[in] d �O���b�h��
 * @param[out] hF �A�֐��l(nx�~ny�~nz�̔z��)
 */
void rxSPH_GPU::CalImplicitFieldDevice(int n[3], Vec3 minp, Vec3 d, RXREAL *dF)
{
	RXREAL *dPos;
	if(m_bUseOpenGL){
		dPos = (RXREAL*)CuMapGLBufferObject(&m_pPosResource);
	}
	else{
		dPos = (RXREAL*)m_dPos;
	}

	// �����Z���̃n�b�V�����v�Z
	CuCalcHash(m_dCellData.dGridParticleHash, m_dCellData.dSortedIndex, dPos, m_dAttr, m_uNumParticles);

	// �n�b�V���Ɋ�Â��p�[�e�B�N�����\�[�g
	CuSort(m_dCellData.dGridParticleHash, m_dCellData.dSortedIndex, m_uNumParticles);

	// �p�[�e�B�N���z����\�[�g���ꂽ���Ԃɕ��ёւ��C
	// �e�Z���̎n�܂�ƏI���̃C���f�b�N�X������
	CuReorderDataAndFindCellStart(m_dCellData, dPos, m_dVel);

	if(m_bCalAnisotropic){
		// �p�[�e�B�N�����x��p�����{�����[���f�[�^(�ٕ����J�[�l��)
		CuSphGridDensityAniso(dF, m_dG, m_fEigenMax, m_dCellData, n[0], n[1], n[2], minp[0], minp[1], minp[2], d[0], d[1], d[2]);
	}
	else{
		// �p�[�e�B�N�����x��p�����{�����[���f�[�^
		CuSphGridDensity(dF, m_dCellData, n[0], n[1], n[2], minp[0], minp[1], minp[2], d[0], d[1], d[2]);
	}

	if(m_bUseOpenGL){
		CuUnmapGLBufferObject(m_pPosResource);
	}
}

//�ǉ��@�X�p
double rxSPH_GPU::GetImplicitSolid(double x, double y, double z, float* fIceCheck)
{
	return CalColorFieldSolid(x, y, z, fIceCheck);
}

/*!�ǉ��@�X�p
 * �p�[�e�B�N������O���b�h�̉A�֐��l���v�Z
 * @param[in] n �O���b�h��
 * @param[in] minp �O���b�h�̍ŏ����W
 * @param[in] d �O���b�h��
 * @param[in] fIceCheck �X���ǂ����̃t���O�@0.0f�����@����ȊO���X
 * @param[out] hF �A�֐��l(nx�~ny�~nz�̔z��)
 */
void rxSPH_GPU::CalImplicitFieldSolid(int n[3], Vec3 minp, Vec3 d, RXREAL *hF, float* fIceCheck)
{
	int slice0 = n[0];
	int slice1 = n[0]*n[1];

	for(int k = 0; k < n[2]; ++k){
		for(int j = 0; j < n[1]; ++j){
			for(int i = 0; i < n[0]; ++i){
				int idx = k*slice1+j*slice0+i;
				Vec3 pos = minp+Vec3(i, j, k)*d;
				hF[idx] = GetImplicitSolid(pos[0], pos[1], pos[2], fIceCheck);
			}
		}
	}
}

/*!�ǉ��@�X�p
 * �p�[�e�B�N������O���b�h�̉A�֐��l���v�Z
 * @param[in] pnx,pny,pnz �O���b�h���̎w�� nx=2^pnx
 * @param[in] minp �O���b�h�̍ŏ����W
 * @param[in] d �O���b�h��
 * @param[in] fIceCheck �X���ǂ����̔���z��@0.0f�����@����ȊO���X
 * @param[out] dF �A�֐��l(nx�~ny�~nz�̔z��)
 */
void rxSPH_GPU::CalImplicitFieldDeviceSolid(int n[3], Vec3 minp, Vec3 d, RXREAL *dF, float* dfIceCheck)
{
//	RXREAL *hF = new RXREAL[n[0]*n[1]*n[2]];

	// CPU��
//	CalImplicitFieldSolid(n, minp, d, hF, fIceCheck);
//	CuCopyArrayToDevice(dF, hF, 0, n[0]*n[1]*n[2]*sizeof(RXREAL));

	// GPU��
	CuIceMeshMake(dF, m_dCellData, n[0], n[1], n[2], minp[0], minp[1], minp[2], d[0], d[1], d[2], dfIceCheck);

//	delete [] hF;
}

/*!
 * �ǉ��@�X�p�J���[�t�B�[���h�l�v�Z
 * @param[in] pos �v�Z�ʒu
 * @return �J���[�t�B�[���h�l
 */
double rxSPH_GPU::CalColorFieldSolid(double x, double y, double z, float* fIceCheck)
{
	RXREAL c = 0.0;
	Vec3 pos(x, y, z);

	if(pos[0] < m_v3EnvMin[0]) return c;
	if(pos[0] > m_v3EnvMax[0]) return c;
	if(pos[1] < m_v3EnvMin[1]) return c;
	if(pos[1] > m_v3EnvMax[1]) return c;
	if(pos[2] < m_v3EnvMin[2]) return c;
	if(pos[2] > m_v3EnvMax[2]) return c;

	RXREAL h = m_params.EffectiveRadius;

	vector<rxNeigh> ne;

	//m_pNNGrid->GetNN(pos, m_hPos, m_uNumParticles, ne, h);
	//GetNN�́Cpos�𒆐S�ɁC���q�z��m_hPos����ɂ����ߖT�����C�e�����ah�ŁCne�ɓ����Ƃ�������
	//�e�O���b�h�ɂǂ̗��q�������Ă���̂���m�肽���D
	//

	//�O���b�h���_�̋ߖT���q�T���@GetNN()��GPU�ɂ͂Ȃ��̂ŁC�ȉ��Œ�`

	// �����Z���C���f�b�N�X�̎Z�o
	int cx = (pos[0]-m_v3EnvMin[0])/m_params.CellWidth.x;
	int cy = (pos[1]-m_v3EnvMin[1])/m_params.CellWidth.y;
	int cz = (pos[2]-m_v3EnvMin[2])/m_params.CellWidth.z;

	int numArdGrid = (int)(h/m_params.CellWidth.x+1);
	for(int k = -numArdGrid; k <= numArdGrid; ++k){
		for(int j = -numArdGrid; j <= numArdGrid; ++j){
			for(int i = -numArdGrid; i <= numArdGrid; ++i){
				int i1 = cx+i;
				int j1 = cy+j;
				int k1 = cz+k;
				if(i1 < 0 || (unsigned)i1 >= m_params.GridSize.x || j1 < 0 || (unsigned)j1 >= m_params.GridSize.y || k1 < 0 || (unsigned)k1 >= m_params.GridSize.z){
					continue;
				}

				getNeighborsInCell(pos, m_hPos, i1, j1, k1, ne, h);
			}
		}
	}

	// �ߖT���q
	for(vector<rxNeigh>::iterator itr = ne.begin(); itr != ne.end(); ++itr){
		int j = itr->Idx;
		if(j < 0) continue;

		Vec3 pos1;
		pos1[0] = m_hPos[DIM*j+0];
		pos1[1] = m_hPos[DIM*j+1];
		pos1[2] = m_hPos[DIM*j+2];

//		RXREAL r = length(pos-pos1);
		RXREAL r = sqrt(itr->Dist2);

		//�e�X�g
//		if(r < 0.95*h){
//			c += m_params.Mass*m_fpW(r, 0.95*h, m_fAw);	//h�Cm_fAw�͒萔
//		}

		if( fIceCheck[j] < 0.0f ) continue;		//�X�łȂ��Ȃ烁�b�V�������Ȃ�
		c += m_params.Mass*m_fpW(r, h, m_fAw);	//h�Cm_fAw�͒萔
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
void rxSPH_GPU::OutputSetting(string fn)
{
	ofstream fout;
	fout.open(fn.c_str());
	if(!fout){
		RXCOUT << fn << " couldn't open." << endl;
		return;
	}

	fout << m_uNumParticles << endl;
	fout << m_params.BoundaryMin.x << " " << m_params.BoundaryMin.y << " " << m_params.BoundaryMin.z << endl;
	fout << m_params.BoundaryMax.x << " " << m_params.BoundaryMax.y << " " << m_params.BoundaryMax.z << endl;
	fout << m_params.Density << endl;
	fout << m_params.Mass    << endl;
	fout << m_params.KernelParticles << endl;

	fout.close();
}
