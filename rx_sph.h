/*!
  @file rx_sph.h
	
  @brief SPH�@
 
  @author Makoto Fujisawa
  @date 2008-10,2011-06
*/
// FILE --rx_sph.h--

#ifndef _RX_SPH_H_
#define _RX_SPH_H_


//-----------------------------------------------------------------------------
// �C���N���[�h�t�@�C��
//-----------------------------------------------------------------------------
#include "rx_sph_commons.h"

#include "rx_ps.h"			// �p�[�e�B�N���V�X�e�����N���X
#include "rx_nnsearch.h"	// �O���b�h�����ɂ��ߖT�T��

#include "rx_sph_solid.h"
#include "rx_wavelet_noise.h"

#include "rx_kernel.h"

#include "rx_cu_common.cuh"



//-----------------------------------------------------------------------------
// ��`
//-----------------------------------------------------------------------------
//#define RX_USE_BOUNDARY_PARTICLE	// ���E�p�[�e�B�N���̗L��

//#define GL_REAL GL_DOUBLE
#define GL_REAL GL_FLOAT

// ���Ԍv��
class rxTimerAvg;
extern rxTimerAvg g_Time;

#define RX_USE_TIMER

#ifdef RX_USE_TIMER
#define RXTIMER_CLEAR g_Time.ClearTime()
#define RXTIMER_RESET g_Time.ResetTime()
#define RXTIMER(x) g_Time.Split(x)
#define RXTIMER_PRINT g_Time.Print()
#define RXTIMER_STRING(x) g_Time.PrintToString(x)
#else
#define RXTIMER_CLEAR
#define RXTIMER_RESET
#define RXTIMER(x) 
#define RXTIMER_PRINT
#define RXTIMER_STRING(x)
#endif

const RXREAL KOL2 = (RXREAL)0.561231024;	// 2^(-5/6)

// �O���[�o���ϐ��̐錾
extern double g_fCoefEt;
extern double g_fMaxEt;
extern double g_fWaveletScale;
extern double g_fMaxEnergySpectrum;
extern double g_fMaxWaveletTurb;
extern double g_fESScale;
extern double g_fVCCoef;
extern double g_fCoefTurb;
extern double g_fCoefTurbForMesh;
extern double g_fEtCri;

extern double g_fNoiseScl;
extern double g_fNoiseEthr;
extern double g_fNoiseMag;


inline RXREAL3 MAKE_RXREAL3(RXREAL x, RXREAL y, RXREAL z)
{
	return make_float3(x, y, z);
}
inline RXREAL2 MAKE_RXREAL2(RXREAL x, RXREAL y)
{
	return make_float2(x, y);
}
inline RXREAL3 MAKE_FLOAT3V(Vec3 x)
{
	return make_float3((FLOAT)x[0], (FLOAT)x[1], (FLOAT)x[2]);
}



//! SPH�V�[���̃p�����[�^
struct rxSPHEnviroment
{
	#define MAX_DELETE_REGIONS 64

	int max_particles;			//!< �ő�p�[�e�B�N����
	Vec3 boundary_cen;			//!< ���E�̒��S
	Vec3 boundary_ext;			//!< ���E�̑傫��(�e�ӂ̒�����1/2)
	RXREAL dens;				//!< �������x
	RXREAL mass;				//!< �p�[�e�B�N���̎���
	RXREAL kernel_particles;	//!< �L�����ah�ȓ��̃p�[�e�B�N����
	RXREAL dt;					//!< ���ԃX�e�b�v��
	RXREAL viscosity;			//!< ���S���W��
	RXREAL gas_k;				//!< �K�X�萔

	int use_inlet;				//!< �������E�����̗L��
	RXREAL et_cri;				//!< �����`���p�̌W��

	RXREAL epsilon;				//!< CFM�̊ɘa�W��
	RXREAL eta;					//!< ���x�ϓ����e��
	int min_iter;				//!< ���R�r�����ŏ���
	int max_iter;				//!< ���R�r�����ő吔

	int use_ap;					//!< �l�H����ON/OFF (0 or 1)
	RXREAL ap_k;				//!< �l�H���͂̂��߂̌W��k (�{��)
	RXREAL ap_n;				//!< �l�H���͂̂��߂̌W��n (n��)
	RXREAL ap_q;				//!< �l�H���͌v�Z���̊�J�[�l���l�v�Z�p�W��(�L�����ah�ɑ΂���W��, [0,1])

	int use_delete_region;		//!< �p�[�e�B�N���폜�̈�̗L��(��)
	Vec3 delete_region[MAX_DELETE_REGIONS][2];	//!< �폜�̈�͈̔�(�ŏ��C�ő���W)

	// �\�ʃ��b�V��
	Vec3 mesh_boundary_cen;		//!< ���b�V���������E�̒��S
	Vec3 mesh_boundary_ext;		//!< ���b�V���������E�̑傫��(�e�ӂ̒�����1/2)
	int mesh_vertex_store;		//!< ���_������|���S������\������Ƃ��̌W��
	int mesh_max_n;				//!< MC�@�p�O���b�h�̍ő啪����

	//�ǉ��F�F�t�@�C������ǂݍ��񂾃p�����[�^
	//�M����
	float htTimeStep;
	float tempMax;
	float tempMin;

	float latentHeat;

	double cffcntHt;
	double cffcntTd;

	//SM�@
	float smTimeStep;

	//�X�\��
	int layer;

	rxSPHEnviroment()
	{
		max_particles = 50000;
		boundary_cen = Vec3(0.0);
		boundary_ext = Vec3(2.0, 0.8, 0.8);
		dens = (RXREAL)998.29;
		mass = (RXREAL)0.04;
		kernel_particles = (RXREAL)20.0;
		dt = 0.01;
		viscosity = 1.0e-3;
		gas_k = 3.0;

		mesh_vertex_store = 10;
		use_inlet = 0;
		mesh_max_n = 128;
		et_cri = 1.0;

		epsilon = 0.001;
		eta = 0.05;
		min_iter = 2;
		max_iter = 10;

		use_ap = true;
		ap_k = 0.1;
		ap_n = 4.0;
		ap_q = 0.2;

		use_delete_region = 0;
		for(int i = 0; i < MAX_DELETE_REGIONS; ++i){
			delete_region[i][0] = Vec3(0.0);
			delete_region[i][1] = Vec3(0.0);
		}
	}
};



//! �\�ʃp�[�e�B�N��
struct rxSurfaceParticle
{
	Vec3 pos;					//!< ���S���W
	Vec3 nrm;					//!< �@��
	Vec3 vel;					//!< ���x
	RXREAL d;					//!< �T�����S����̋���
	int idx;					//!< �p�[�e�B�N���C���f�b�N�X
};

extern double g_fSurfThr[2];


//-----------------------------------------------------------------------------
// MARK:rxSPH�N���X�̐錾
//-----------------------------------------------------------------------------
class rxSPH : public rxParticleSystemBase
{
private:
	// �p�[�e�B�N��
	RXREAL *m_hNrm;					//!< �p�[�e�B�N���@��
	RXREAL *m_hFrc;					//!< �p�[�e�B�N���ɂ������
	RXREAL *m_hDens;				//!< �p�[�e�B�N�����x
	RXREAL *m_hPres;				//!< �p�[�e�B�N������

	// �\�ʐ����p(Anisotropic kernel)
	RXREAL *m_hUpPos;				//!< �������p�[�e�B�N���ʒu
	RXREAL *m_hPosW;				//!< �d�ݕt�����ύ��W
	RXREAL *m_hCMatrix;				//!< �����U�s��
	RXREAL *m_hEigen;				//!< �����U�s��̓��ْl
	RXREAL *m_hRMatrix;				//!< ��]�s��(�����U�s��̓��كx�N�g��)
	RXREAL *m_hG;					//!< �ό`�s��

	uint *m_hSurf;					//!< �\�ʃp�[�e�B�N��

	// ���E�E�ő�
	rxSolid *m_pBoundary;			//!< �V�~�����[�V������Ԃ̋��E
	vector<rxSolid*> m_vSolids;		//!< �ő̕���
	RXREAL *m_hVrts;				//!< �ő̃|���S���̒��_
	int m_iNumVrts;					//!< �ő̃|���S���̒��_��
	int *m_hTris;					//!< �ő̃|���S��
	int m_iNumTris;					//!< �ő̃|���S���̐�
	RXREAL *m_hSVels;				//!< �ő̃|���S�����x

	// ��ԕ����i�q�֘A
	rxNNGrid *m_pNNGrid;			//!< �����O���b�h�ɂ��ߖT�T��
	vector< vector<rxNeigh> > m_vNeighs;	//!< �ߖT�p�[�e�B�N��

	rxNNGrid *m_pNNGridB;			//!< ���E�p�[�e�B�N���p�����O���b�h

	// ���q�p�����[�^
	uint m_iKernelParticles;		//!< �J�[�l�����̃p�[�e�B�N����
	RXREAL m_fRestDens, m_fMass;	//!< ���x�C����
	RXREAL m_fEffectiveRadius;		//!< �L�����a
	RXREAL m_fKernelRadius;			//!< �J�[�l���̉e���͈�

	// �V�~�����[�V�����p�����[�^
	RXREAL m_fGasStiffness;			//!< �K�X�萔
	RXREAL m_fViscosity;			//!< �S���W��
	RXREAL m_fBuoyancy;				//!< ����

	// �J�[�l���֐��̌v�Z�̍ۂɗp������萔�W��
	double m_fAw;
	double m_fAg;
	double m_fAl;

	// �J�[�l���֐�
	double (*m_fpW)(double, double, double);
	Vec3 (*m_fpGW)(double, double, double, Vec3);
	double (*m_fpLW)(double, double, double, double);

protected:
	rxSPH(){}

public:
	//! �R���X�g���N�^
	rxSPH(bool use_opengl);

	//! �f�X�g���N�^
	virtual ~rxSPH();

	// �p�[�e�B�N�����a
	float GetEffectiveRadius(){ return m_fEffectiveRadius; }

	// �ߖT�p�[�e�B�N��
	uint* GetNeighborList(const int &i, int &n);

public:
	//
	// ���z�֐�
	//
	// �p�[�e�B�N���f�[�^
	virtual RXREAL* GetParticle(void);
	virtual RXREAL* GetParticleDevice(void);
	virtual void UnmapParticle(void){}

	// �V�~�����[�V�����X�e�b�v
	virtual bool Update(RXREAL dt, int step = 0);

	// �V�[���̐ݒ�
	virtual void SetPolygonObstacle(const vector<Vec3> &vrts, const vector<Vec3> &nrms, const vector< vector<int> > &tris, Vec3 vel);
	virtual void SetBoxObstacle(Vec3 cen, Vec3 ext, Vec3 ang, Vec3 vel, int flg);
	virtual void SetSphereObstacle(Vec3 cen, double rad, Vec3 vel, int flg);
	virtual void MoveSphereObstacle(int b, Vec3 disp);
	virtual Vec3 GetSphereObstaclePos(int b = -1);

	// �z�X�g<->VBO�ԓ]��
	virtual RXREAL* GetArrayVBO(rxParticleArray type, bool d2h = true, int num = -1);
	virtual void SetArrayVBO(rxParticleArray type, const RXREAL* data, int start, int count);
	virtual void SetArrayVBO(rxParticleArray type, const int* data, int start, int count);
	virtual void SetColorVBO(int type);

	// �A�֐��l�v�Z
	virtual double GetImplicit(double x, double y, double z);
	virtual void CalImplicitField(int n[3], Vec3 minp, Vec3 d, RXREAL *hF);
	virtual void CalImplicitFieldDevice(int n[3], Vec3 minp, Vec3 d, RXREAL *dF);
	//�ǉ�
	double GetImplicitSolid(double x, double y, double z);
	void CalImplicitFieldSolid(int n[3], Vec3 minp, Vec3 d, RXREAL *hF);
	void CalImplicitFieldDeviceSolid(int n[3], Vec3 minp, Vec3 d, RXREAL *dF, float *m_fIntrps);

	// SPH���o��
	virtual void OutputSetting(string fn);

	// �`��֐�
	virtual void DrawCell(int i, int j, int k);
	virtual void DrawCells(Vec3 col, Vec3 col2, int sel = 0);
	virtual void DrawObstacles(void);


public:
	// SPH������
	void Initialize(const rxSPHEnviroment &env);
	void Allocate(int max_particles);
	void Finalize(void);

	// �ߖT�擾
	void GetNearestNeighbors(Vec3 pos, vector<rxNeigh> &neighs, RXREAL h = -1.0);
	void GetNearestNeighbors(int idx, RXREAL *p, vector<rxNeigh> &neighs, RXREAL h = -1.0);

	void GetNearestNeighborsB(Vec3 pos, vector<rxNeigh> &neighs, RXREAL h = -1.0);

	// �����Z���Ƀp�[�e�B�N�����i�[
	virtual void SetParticlesToCell(void);
	virtual void SetParticlesToCell(RXREAL *prts, int n, RXREAL h);

	// ���^�{�[���ɂ��A�֐��l
	double CalColorField(double x, double y, double z);
	//�ǉ�
	double CalColorFieldSolid(double x, double y, double z);

	// �����Z���Ƀ|���S�����i�[
	virtual void SetPolygonsToCell(void);

	int  GetPolygonsInCell(int gi, int gj, int gk, vector<int> &polys);
	bool IsPolygonsInCell(int gi, int gj, int gk);

	
	// �\�ʃp�[�e�B�N�����o
	void DetectSurfaceParticles(void);					// �\�ʃp�[�e�B�N���̌��o
	double CalDistToNormalizedMassCenter(const int i);	// �ߖT�p�[�e�B�N���̏d�S�܂ł̋����v�Z
	uint* GetArraySurf(void);							// �\�ʃp�[�e�B�N�����̎擾

	// �\�ʃp�[�e�B�N�����̎擾
	int GetSurfaceParticles(const Vec3 pos, RXREAL h, vector<rxSurfaceParticle> &sp);

	// �@���v�Z
	void CalNormalFromDensity(void);
	void CalNormal(void);

	//�ǉ��F�F�ߖT���q�̎擾
	vector<vector<rxNeigh>>& GetNeights(void){ return m_vNeighs; }

	//
	// Anisotropic kernel
	//
	virtual void CalAnisotropicKernel(void);

	//
	// SPS����
	//
	int	GetMaxSubLevel() const { return 0; }

	// �T�u�p�[�e�B�N���f�[�^(�f�o�C�X������)
	RXREAL* GetSubParticleDev(void){ return 0; }
	void    UnmapSubParticle(void){}
	RXREAL* GetAllSubParticleDev(void){ return 0; }
	void    UnmapAllSubParticle(void){}
	RXREAL* GetSubParticlePosDev(void){ return 0; }
	RXREAL* GetSubParticleRadDev(void){ return 0; }
	RXREAL* GetSubParticleRatDev(void){ return 0; }

	int GetNumAllSubParticles(void){ return 0; }
	int GetNumValidSubParticles(void){ return 0; }

	// �T�u�p�[�e�B�N���f�[�^(�z�X�g������)
	RXREAL* GetSubParticlePos(void){ return 0; }
	RXREAL* GetSubParticleRad(void){ return 0; }
	RXREAL* GetSubParticleRat(void){ return 0; }
	unsigned int* GetSubParticleOcc(void){ return 0; }
	unsigned int* GetSubParticleOccScan(void){ return 0; }

	RXREAL GetSubRadius(int level){ return 0.0f; }

	//�T�u���q
	void SetEtcri(const double &et_cri){}
	void InitSubParticles(const double &et_cri){}
	void AddSubParticles(int start, int count){}
	void UpdateSubParticle(RXREAL *dPos, RXREAL *dSubPos, RXREAL scale, RXREAL dt){}
	void CalRadiusAndRatio(bool force = false){}

	void CalSubImplicitFieldDevice(int n[3], Vec3 minp, Vec3 d, RXREAL *dF, int turb){}
	void CalSubImplicitField(int n[3], Vec3 minp, Vec3 d, RXREAL *hF, const double &et_cri){}

	void AddSphereField(int n[3], Vec3 minp, Vec3 d, RXREAL *hF, 
						const Vec3 &pos, const double &eRad, const double& ratio){}



protected:
	// CPU�ɂ��SPH�v�Z
	void calDensity(const RXREAL *pos, RXREAL *dens, RXREAL h);
	void calNormal(void);
	void calForce(const RXREAL *ppos, const RXREAL *pvel, const RXREAL *pdens, RXREAL *ppres, RXREAL *pfrc, RXREAL h);

	// rest density�̌v�Z
	RXREAL calRestDensity(RXREAL h);

	// �X�̋��E�p�[�e�B�N���̑̐ς��v�Z
	void calBoundaryVolumes(const RXREAL *bpos, RXREAL *bvol, RXREAL mass, uint n, RXREAL h);

	// ���ԃX�e�b�v���̏C��
	RXREAL calTimeStep(RXREAL &dt, RXREAL eta_avg, const RXREAL *pfrc, const RXREAL *pvel, const RXREAL *pdens);

	// �ʒu�Ƒ��x�̍X�V
	void integrate(const RXREAL *pos, const RXREAL *vel, const RXREAL *dens, const RXREAL *acc, 
				   RXREAL *pos_new, RXREAL *vel_new, RXREAL dt);

	// �Փ˔���
	int calCollisionPolygon(uint grid_hash, Vec3 &pos0, Vec3 &pos1, Vec3 &vel, RXREAL dt);
	int calCollisionSolid(Vec3 &pos0, Vec3 &pos1, Vec3 &vel, RXREAL dt);
};


//-----------------------------------------------------------------------------
// MARK:rxSPH_GPU�N���X�̐錾
//-----------------------------------------------------------------------------
class rxSPH_GPU : public rxParticleSystemBase
{
private:
	//
	// �����o�ϐ�(GPU�ϐ�)
	//
	RXREAL *m_dPos;			//!< �p�[�e�B�N���ʒu
	RXREAL *m_dVel;			//!< �p�[�e�B�N�����x
	RXREAL *m_dNrm;			//!< �p�[�e�B�N���@��
	RXREAL *m_dFrc;			//!< �p�[�e�B�N���ɂ������
	RXREAL *m_dDens;		//!< �p�[�e�B�N�����x
	RXREAL *m_dPres;		//!< �p�[�e�B�N������

	RXREAL *m_dPosB;		//!< ���E�p�[�e�B�N��
	RXREAL *m_dVolB;		//!< ���E�p�[�e�B�N���̑̐�

	int *m_dAttr;			//!< �p�[�e�B�N������

	cudaGraphicsResource *m_pPosResource;	//!< OpenGL(��VBO)-CUDA�Ԃ̃f�[�^�]�����������߂̃n���h��

	// �E�F�[�u���b�g����
	RXREAL *m_dEt;			//!< ���x�G�l���M�[�X�y�N�g��
	RXREAL *m_dTurb;		//!< �E�F�[�u���b�g�����ɂ�鑬�x��

	//�T�u���q GPU data
	RXREAL *m_dSubPos;		//!< �T�u���q�̐�΍��W float4(0.0,x,y,z)
	RXREAL *m_dSubChild;	//!< �T�u���q�̎q1�ւ̒P�ʃx�N�g��
	RXREAL *m_dSubAxis;		//!< �T�u���q�̉�]��(�P�ʃx�N�g��)
	RXREAL *m_dSubEt;		//!< �T�u���q�̃G�l���M�[�X�y�N�g��
	uint   *m_dSubRand;		//!< �����e�[�u��

	uint m_subposVBO;		//!< �T�u���q���WVBO
	cudaGraphicsResource *m_pSubPosResource;	//!< OpenGL(��VBO)-CUDA�Ԃ̃f�[�^�]�����������߂̃n���h��

	// SPS����
	RXREAL *m_dUpPos;		//!< �������p�[�e�B�N���ʒu
	RXREAL *m_dPosW;		//!< �d�ݕt�����ύ��W
	RXREAL *m_dCMatrix;		//!< �����U�s��
	RXREAL *m_dEigen;		//!< �����U�s��̓��ْl
	RXREAL *m_dRMatrix;		//!< ��]�s��(�����U�s��̓��كx�N�g��)
	RXREAL *m_dG;			//!< �ό`�s��

	// �\�ʃ��b�V��
	RXREAL *m_dVrts;		//!< �ő̃��b�V�����_
	int    *m_dTris;		//!< �ő̃��b�V��

	// �V�~�����[�V�����p�����[�^
	rxSimParams m_params;	//!< �V�~�����[�V�����p�����[�^(GPU�ւ̃f�[�^�n���p)
	uint3 m_gridSize;		//!< �ߖT�T���O���b�h�̊e���̕�����
	uint m_numGridCells;	//!< �ߖT�T���O���b�h��������
	
	// ��ԕ���(GPU)
	rxParticleCell m_dCellData;	//!< �ߖT�T���O���b�h
	rxSubParticleCell m_dSubCellData;	//!< �T�u���q�p�ߖT�T���O���b�h
	uint m_gridSortBits;		//!< �n�b�V���l�ɂ���\�[�g���̊����

	rxParticleCell m_dCellDataB;	//!< ���E�p�[�e�B�N���p�ߖT�T���O���b�h
	uint3 m_gridSizeB;		//!< �ߖT�T���O���b�h�̊e���̕�����

	//
	// �����o�ϐ�(CPU�ϐ�)
	//
	RXREAL *m_hNrm;			//!< �p�[�e�B�N���@��
	RXREAL *m_hFrc;			//!< �p�[�e�B�N���ɂ������
	RXREAL *m_hDens;		//!< �p�[�e�B�N�����x
	RXREAL *m_hPres;		//!< �p�[�e�B�N������

	uint *m_hSurf;			//!< �\�ʃp�[�e�B�N��

	// �\�ʐ����p(Anisotropic kernel)
	RXREAL *m_hUpPos;		//!< �������p�[�e�B�N���ʒu
	RXREAL *m_hPosW;		//!< �d�ݕt�����ύ��W
	RXREAL *m_hCMatrix;		//!< �����U�s��
	RXREAL *m_hEigen;		//!< �����U�s��̓��ْl
	RXREAL *m_hRMatrix;		//!< ��]�s��(�����U�s��̓��كx�N�g��)
	RXREAL *m_hG;			//!< �ό`�s��
	RXREAL  m_fEigenMax;	//!< ���ْl�̍ő�l(�T�����a�g���ɗp����)

	// �E�F�[�u���b�g����
	RXREAL *m_hVwt;			//!< �p�[�e�B�N�����x�̃E�F�[�u���b�g�ϊ�
	RXREAL *m_hEt;			//!< ���x�̃G�l���M�[�X�y�N�g����
	RXREAL *m_hTurb;		//!< �E�F�[�u���b�g�����ɂ�鑬�x��
	RXREAL *m_hNoiseTile[3];//!< �E�F�[�u���b�g�m�C�Y�^�C��
	RXREAL m_fNoiseTileWidth;	//!< �m�C�Y�^�C���̕�
	uint   m_iNoiseTileN[3];	//!< �m�C�Y�^�C���̉𑜓x

	// SPS����
	RXREAL *m_hSubPos;		//!< �T�u���q�̐�΍��W
	RXREAL *m_hSubChild;	//!< �T�u���q�̎q1�ւ̒P�ʃx�N�g��
	RXREAL *m_hSubAxis;		//!< �T�u���q�̉�]��(�P�ʃx�N�g��)
	RXREAL *m_hSubEt;		//!< �T�u���q�̃G�l���M�[�X�y�N�g��

	uint m_uNumSubParticles;//!< �T�u���q��
	uint m_uMaxSubParticles;//!< �T�u���q�ő吔
	uint m_uMaxSubLevel;	//!< �T�u���q�ő僌�x��

	RXREAL *m_hSubPosPack;
	RXREAL *m_hSubRadPack;
	RXREAL *m_hSubRatPack;
	unsigned int *m_hSubOcc;
	unsigned int *m_hSubOccScan;
	bool m_bSubPacked;		//!< �T�u�p�[�e�B�N����Pack������true

	// �\�ʃ��b�V��
	vector<RXREAL> m_vVrts;	//!< �ő̃��b�V�����_
	int m_iNumVrts;			//!< �ő̃��b�V�����_��
	vector<int> m_vTris;	//!< �ő̃��b�V��
	int m_iNumTris;			//!< �ő̃��b�V����

	bool m_bCalNormal;		//!< �@���v�Z�t���O

	//�ǉ��F�ߖT���q��GPU����̃R�s�[�f�[�^
	uint* m_hSortedIndex;		//!< �n�b�V���l�Ń\�[�g�����p�[�e�B�N���C���f�b�N�X
	uint* m_hGridParticleHash;	//!< �e�p�[�e�B�N���̃O���b�h�n�b�V���l
	uint* m_hCellStart;			//!< �\�[�g���X�g���̊e�Z���̃X�^�[�g�C���f�b�N�X
	uint* m_hCellEnd;			//!< �\�[�g���X�g���̊e�Z���̃G���h�C���f�b�N�X
	uint  m_uNumCells;			//!< ���Z����
	vector< vector<rxNeigh> > m_vNeighs;	//!< �ߖT�p�[�e�B�N��

	// �J�[�l���֐��̌v�Z�̍ۂɗp������萔�W��
	double m_fAw;
	double m_fAg;
	double m_fAl;

	// �J�[�l���֐�
	double (*m_fpW)(double, double, double);
	Vec3 (*m_fpGW)(double, double, double, Vec3);
	double (*m_fpLW)(double, double, double, double);

protected:
	rxSPH_GPU(){}

public:
	//! �R���X�g���N�^
	rxSPH_GPU(bool use_opengl);

	//! �f�X�g���N�^
	~rxSPH_GPU();

	// �p�[�e�B�N�����a
	//�ǉ��F
	void SetEffectiveRadius(float r){	m_params.EffectiveRadius = r;	}
	float GetEffectiveRadius(){ return m_params.EffectiveRadius; }
	
	// �ߖT�p�[�e�B�N��
	uint* GetNeighborList(const int &i, int &n);

	// �V�~�����[�V�����p�����[�^
	rxSimParams GetParams(void){ return m_params; }
	void UpdateParams(void);

	// �t���O�ؑ�
	void ToggleNormalCalc(int t = -1){ RX_TOGGLE(m_bCalNormal, t); }			//!< �p�[�e�B�N���@���̌v�Z
	bool IsNormalCalc(void) const { return m_bCalNormal; }
#if MAX_BOX_NUM
	void ToggleSolidFlg(int t = -1);
#endif

public:
	//
	// ���z�֐�
	//
	// �p�[�e�B�N���f�[�^
	virtual RXREAL* GetParticle(void);
	virtual RXREAL* GetParticleDevice(void);
	virtual void UnmapParticle(void);

	// �V�~�����[�V�����X�e�b�v
	virtual bool Update(RXREAL dt, int step = 0);

	// �V�[���̐ݒ�
	virtual void SetPolygonObstacle(const vector<Vec3> &vrts, const vector<Vec3> &nrms, const vector< vector<int> > &tris, Vec3 vel);
	virtual void SetBoxObstacle(Vec3 cen, Vec3 ext, Vec3 ang, Vec3 vel, int flg);
	virtual void SetSphereObstacle(Vec3 cen, double rad, Vec3 vel, int flg);
	virtual void MoveSphereObstacle(int b, Vec3 disp);
	virtual Vec3 GetSphereObstaclePos(int b = -1);

	// �z�X�g<->VBO�ԓ]��
	virtual RXREAL* GetArrayVBO(rxParticleArray type, bool d2h = true, int num = -1);
	virtual void SetArrayVBO(rxParticleArray type, const RXREAL* data, int start, int count);
	virtual void SetArrayVBO(rxParticleArray type, const int* data, int start, int count);
	virtual void SetColorVBO(int type);


	// �A�֐��l�v�Z
	virtual double GetImplicit(double x, double y, double z);
	virtual void CalImplicitField(int n[3], Vec3 minp, Vec3 d, RXREAL *hF);
	virtual void CalImplicitFieldDevice(int n[3], Vec3 minp, Vec3 d, RXREAL *dF);

	//�ǉ�
	double GetImplicitSolid(double x, double y, double z, float* fIceCheck);
	void CalImplicitFieldSolid(int n[3], Vec3 minp, Vec3 d, RXREAL *hF, float* fIceCheck);
	void CalImplicitFieldDeviceSolid(int n[3], Vec3 minp, Vec3 d, RXREAL *dF, float* fIceCheck);

	// SPH���o��
	virtual void OutputSetting(string fn);

	// �`��֐�
	virtual void DrawCell(int i, int j, int k);
	virtual void DrawCells(Vec3 col, Vec3 col2, int sel = 0);
	virtual void DrawObstacles(void);

protected:
	void setObjectToCell(RXREAL *p);


public:
	// SPH������
	void Initialize(const rxSPHEnviroment &env);
	void Allocate(int max_particles);
	void Finalize(void);

	// �����Z���Ƀp�[�e�B�N�����i�[
	virtual void SetParticlesToCell(void);
	virtual void SetParticlesToCell(RXREAL *prts, int n, RXREAL h);
	virtual void SetPolygonsToCell(void);

	int  GetPolygonsInCell(int gi, int gj, int gk, vector<int> &polys){ return 0; }
	bool IsPolygonsInCell(int gi, int gj, int gk){ return false; }
	
	//�ǉ��F�F�J���[�t�B�[���h�l�v�Z�H
	double CalColorFieldSolid(double x, double y, double z, float* fIceCheck);
	
	void CalMaxDensity(int k);

	// �\�ʃp�[�e�B�N�����o
	void DetectSurfaceParticles(void);					// �\�ʃp�[�e�B�N���̌��o
	double CalDistToNormalizedMassCenter(const int i);	// �ߖT�p�[�e�B�N���̏d�S�܂ł̋����v�Z
	uint* GetArraySurf(void);							// �\�ʃp�[�e�B�N�����̎擾

	// �\�ʃp�[�e�B�N�����̎擾
	int GetSurfaceParticles(const Vec3 pos, RXREAL h, vector<rxSurfaceParticle> &sp);

	// �@���v�Z
	void CalNormalFromDensity(void);
	void CalNormal(void);

	// ���E�p�[�e�B�N���̏�����
	void InitBoundary(void);

	//�ǉ��F�F�ߖT���q�擾
	vector<vector<rxNeigh>>& GetNeights(void){ return m_vNeighs; }

	//
	// Anisotropic kernel
	//
	virtual void CalAnisotropicKernel(void);

	//
	// �E�F�[�u���b�g����
	//
	void CalParticleCWT(RXREAL scale);
	void CalWaveletTurbulence(RXREAL scale, RXREAL *dPos, RXREAL dt);

	//
	// SPS����
	//
	int	GetMaxSubLevel() const { return m_uMaxSubLevel; }

	// �T�u�p�[�e�B�N���f�[�^(�f�o�C�X������)
	RXREAL* GetSubParticleDev(void);
	void    UnmapSubParticle(void);
	RXREAL* GetAllSubParticleDev(void);
	void    UnmapAllSubParticle(void);
	RXREAL* GetSubParticlePosDev(void);
	RXREAL* GetSubParticleRadDev(void);
	RXREAL* GetSubParticleRatDev(void);

	int GetNumAllSubParticles(void);
	int GetNumValidSubParticles(void);

	// �T�u�p�[�e�B�N���f�[�^(�z�X�g������)
	RXREAL* GetSubParticlePos(void);
	RXREAL* GetSubParticleRad(void);
	RXREAL* GetSubParticleRat(void);
	unsigned int* GetSubParticleOcc(void);
	unsigned int* GetSubParticleOccScan(void);

	RXREAL GetSubRadius(int level){ return m_dSubCellData.fSubRad[level]; }

	//�T�u���q
	void SetEtcri(const double &et_cri){ m_dSubCellData.fEtcri = et_cri;}
	void InitSubParticles(const double &et_cri);
	void AddSubParticles(int start, int count);
	void UpdateSubParticle(RXREAL *dPos, RXREAL *dSubPos, RXREAL scale, RXREAL dt);
	void CalRadiusAndRatio(bool force = false);

	void CalSubImplicitFieldDevice(int n[3], Vec3 minp, Vec3 d, RXREAL *dF, int turb);
	void CalSubImplicitField(int n[3], Vec3 minp, Vec3 d, RXREAL *hF, const double &et_cri);

	void AddSphereField(int n[3], Vec3 minp, Vec3 d, RXREAL *hF, 
						const Vec3 &pos, const double &eRad, const double& ratio);

protected:
	// rest density�̌v�Z
	RXREAL calRestDensity(RXREAL h);

	// �����Z���̏����ݒ�
	void setupCells(rxParticleCell &cell, uint3 &gridsize, double &cell_width, Vec3 vMin, Vec3 vMax, double h);

	// �O���b�h�n�b�V���̌v�Z
	uint calGridHash(uint x, uint y, uint z);
	uint calGridHash(Vec3 pos);

	// �|���S���𕪊��Z���Ɋi�[
	void setPolysToCell(RXREAL *vrts, int nv, int* tris, int nt);

	// CPU�p�̋ߖT���q�T��
	void searchNeighbors(void);
	void getNeighborsInCell(Vec3 pos, RXREAL *p, int gi, int gj, int gk, vector<rxNeigh> &neighs, RXREAL h);
};




//-----------------------------------------------------------------------------
// MARK:rxDDSPH�N���X�̐錾
//  - Double Density Relaxation �ɂ��SPH
//  - S.Clavet et al., "Particle-based Viscoelastic Fluid Simulation", SCA2005, 2005. 
//  - http://www.iro.umontreal.ca/labs/infographie/papers/Clavet-2005-PVFS/index.html
//-----------------------------------------------------------------------------
class rxDDSPH : public rxParticleSystemBase
{
	//! ���q�ԃo�l
	class rxSPHSpring
	{
	public:
		bool enable;
		double L;
		double r;
		int pi, pj;

		rxSPHSpring()
		{
			enable = false;
			L = r = 0.0;
		}

		~rxSPHSpring()
		{
		}
	};

protected:
	// �p�[�e�B�N��
	RXREAL *m_hNrm;					//!< �p�[�e�B�N���@��
	RXREAL *m_hFrc;					//!< �p�[�e�B�N���ɂ������
	RXREAL *m_hDens;				//!< �p�[�e�B�N�����x
	RXREAL *m_hPres;				//!< �p�[�e�B�N������
	RXREAL *m_hVelOld;				//!< �O�X�e�b�v���x

	RXREAL *m_hPredictPos;			//!< �\���ʒu
	RXREAL *m_hDist;				//!< �p�[�e�B�N���ԋ���

	uint *m_hSurf;					//!< �\�ʃp�[�e�B�N��


	// ���E�E�ő�
	rxSolid *m_pBoundary;			//!< �V�~�����[�V������Ԃ̋��E
	vector<rxSolid*> m_vSolids;		//!< �ő̕���
	RXREAL *m_hVrts;				//!< �ő̃|���S���̒��_
	int m_iNumVrts;					//!< �ő̃|���S���̒��_��
	int *m_hTris;					//!< �ő̃|���S��
	int m_iNumTris;					//!< �ő̃|���S���̐�

	// �ߖT���q�T��
	rxNNGrid *m_pNNGrid;			//!< �����O���b�h�ɂ��ߖT�T��
	vector< vector<rxNeigh> > m_vNeighs;	//!< �ߖT�p�[�e�B�N��

	// ���q�p�����[�^
	uint m_iKernelParticles;		//!< �J�[�l�����̃p�[�e�B�N����
	RXREAL m_fMass;					//!< ����
	RXREAL m_fEffectiveRadius;		//!< �L�����a

	// �V�~�����[�V�����p�����[�^
	RXREAL m_fBuoyancy;				//!< ����

	// Double Density�p�p�����[�^
	RXREAL m_fK;					//!< �K�X�萔
	RXREAL m_fKnear;				//!< near density�p�K�X�萔
	RXREAL m_fViscC;
	RXREAL m_fViscBeta;
	RXREAL m_fElasKspring;
	RXREAL m_fPlasAlpha;
	RXREAL m_fPlasGamma;

	// �S�e���͗p���q�ԃX�v�����O
	map<int, rxSPHSpring> m_mapSpring;

	// ���x
	double m_fInitDens;				//!< �������ϖ��x
	double m_fInitDensNear;			//!< �������ϖ��x(near)
	double m_fInitMaxDens;			//!< �����ő喧�x

	double m_fMaxDens;				//!< �ő喧�x
	double m_fSurfDens;				//!< �\�ʗ��q�Ƃ��閧�x�̔䗦

	bool m_bCalDens;				//!< rest density �v�Z�t���O

	// �J�[�l���֐��̌v�Z�̍ۂɗp������萔�W��
	double m_fWpoly6;				//!< Pory6�J�[�l���̒萔�W��
	double m_fGWpoly6;				//!< Pory6�J�[�l���̌��z�̒萔�W��
	double m_fLWpoly6;				//!< Pory6�J�[�l���̃��v���V�A���̒萔�W��

protected:
	rxDDSPH(){}

public:
	//! �R���X�g���N�^
	rxDDSPH(bool use_opengl);

	//! �f�X�g���N�^
	virtual ~rxDDSPH();

	// �p�[�e�B�N�����a
	float GetEffectiveRadius(){ return m_fEffectiveRadius; }

	// �ߖT�p�[�e�B�N��
	uint* GetNeighborList(const int &i, int &n);

public:
	//
	// ���z�֐�
	//
	// �p�[�e�B�N���f�[�^
	virtual RXREAL* GetParticle(void);
	virtual RXREAL* GetParticleDevice(void);
	virtual void UnmapParticle(void){}

	// �V�~�����[�V�����X�e�b�v
	virtual bool Update(RXREAL dt, int step = 0);

	// �V�[���̐ݒ�
	virtual void SetPolygonObstacle(const vector<Vec3> &vrts, const vector<Vec3> &nrms, const vector< vector<int> > &tris, Vec3 vel);
	virtual void SetBoxObstacle(Vec3 cen, Vec3 ext, Vec3 ang, Vec3 vel, int flg);
	virtual void SetSphereObstacle(Vec3 cen, double rad, Vec3 vel, int flg);
	virtual void MoveSphereObstacle(int b, Vec3 disp);
	virtual Vec3 GetSphereObstaclePos(int b = -1);

	// �z�X�g<->VBO�ԓ]��
	virtual RXREAL* GetArrayVBO(rxParticleArray type, bool d2h = true, int num = -1);
	virtual void SetArrayVBO(rxParticleArray type, const RXREAL* data, int start, int count);
	virtual void SetArrayVBO(rxParticleArray type, const int* data, int start, int count){}
	virtual void SetColorVBO(int type);
	
	// �A�֐��l�v�Z
	virtual double GetImplicit(double x, double y, double z);
	virtual void CalImplicitField(int n[3], Vec3 minp, Vec3 d, RXREAL *hF);
	virtual void CalImplicitFieldDevice(int n[3], Vec3 minp, Vec3 d, RXREAL *dF);

	// SPH���o��
	virtual void OutputSetting(string fn);

	// �`��֐�
	virtual void DrawCell(int i, int j, int k);
	virtual void DrawCells(Vec3 col, Vec3 col2, int sel = 0);
	virtual void DrawObstacles(void);

	
public:
	// SPH������
	void Initialize(const rxSPHEnviroment &env);
	void Allocate(int max_particles);
	void Finalize(void);

	// �ߖT�擾
	void GetNearestNeighbors(Vec3 pos, vector<rxNeigh> &neighs, RXREAL h = -1.0);
	void GetNearestNeighbors(int idx, RXREAL *p, vector<rxNeigh> &neighs, RXREAL h = -1.0);

	// �����Z���Ƀp�[�e�B�N�����i�[
	virtual void SetParticlesToCell(void);
	virtual void SetParticlesToCell(RXREAL *prts, int n, RXREAL h);

	// ���^�{�[���ɂ��A�֐��l
	double CalColorField(double x, double y, double z);

	// �����Z���Ƀ|���S�����i�[
	virtual void SetPolygonsToCell(void);

	int  GetPolygonsInCell(int gi, int gj, int gk, vector<int> &polys);
	bool IsPolygonsInCell(int gi, int gj, int gk);

	// �\�ʃp�[�e�B�N�����o
	void DetectSurfaceParticles(void);					// �\�ʃp�[�e�B�N���̌��o
	double CalDistToNormalizedMassCenter(const int i);	// �ߖT�p�[�e�B�N���̏d�S�܂ł̋����v�Z
	uint* GetArraySurf(void);							// �\�ʃp�[�e�B�N�����̎擾

	// �\�ʃp�[�e�B�N�����̎擾
	int GetSurfaceParticles(const Vec3 pos, RXREAL h, vector<rxSurfaceParticle> &sp);

	// �@���v�Z
	void CalNormalFromDensity(void);
	void CalNormal(void);

	//�ǉ��F�F�ߖT���q�̎擾
	vector<vector<rxNeigh>>& GetNeights(void){ return m_vNeighs; }

protected:
	// Double density relaxation method
	void calExternalForceToVel(RXREAL dt);
	void calDoubleDensity(RXREAL *ppos);
	void calDoubleDensity(RXREAL *ppos, double &avg_dens, double &avg_dens_near, double &max_dens);
	void calDoubleDensityRelaxation(RXREAL *ppos, RXREAL dt);
	void adjustSprings(RXREAL dt);
	void applySpringDisplacements(RXREAL dt);
	
	// �Փ˔���
	int calCollisionPolygon(uint grid_hash, Vec3 &pos0, Vec3 &pos1, Vec3 &vel, RXREAL dt);
	int calCollisionSolid(Vec3 &pos0, Vec3 &pos1, Vec3 &vel, RXREAL dt);

public:
	//
	// Anisotropic kernel - DD�p�͖�����
	//
	virtual void CalAnisotropicKernel(void){}

	//
	// SPS���� - CPU�p�͖�����
	//
	int	GetMaxSubLevel() const { return 0; }

	// �T�u�p�[�e�B�N���f�[�^(�f�o�C�X������)
	RXREAL* GetSubParticleDev(void){ return 0; }
	void    UnmapSubParticle(void){}
	RXREAL* GetAllSubParticleDev(void){ return 0; }
	void    UnmapAllSubParticle(void){}
	RXREAL* GetSubParticlePosDev(void){ return 0; }
	RXREAL* GetSubParticleRadDev(void){ return 0; }
	RXREAL* GetSubParticleRatDev(void){ return 0; }

	int GetNumAllSubParticles(void){ return 0; }
	int GetNumValidSubParticles(void){ return 0; }

	// �T�u�p�[�e�B�N���f�[�^(�z�X�g������)
	RXREAL* GetSubParticlePos(void){ return 0; }
	RXREAL* GetSubParticleRad(void){ return 0; }
	RXREAL* GetSubParticleRat(void){ return 0; }
	unsigned int* GetSubParticleOcc(void){ return 0; }
	unsigned int* GetSubParticleOccScan(void){ return 0; }

	RXREAL GetSubRadius(int level){ return 0.0f; }

	//�T�u���q
	void SetEtcri(const double &et_cri){}
	void InitSubParticles(const double &et_cri){}
	void AddSubParticles(int start, int count){}
	void UpdateSubParticle(RXREAL *dPos, RXREAL *dSubPos, RXREAL scale, RXREAL dt){}
	void CalRadiusAndRatio(bool force = false){}

	void CalSubImplicitFieldDevice(int n[3], Vec3 minp, Vec3 d, RXREAL *dF, int turb){}
	void CalSubImplicitField(int n[3], Vec3 minp, Vec3 d, RXREAL *hF, const double &et_cri){}

	void AddSphereField(int n[3], Vec3 minp, Vec3 d, RXREAL *hF, 
						const Vec3 &pos, const double &eRad, const double& ratio){}
};




//-----------------------------------------------------------------------------
// MARK:rxPBDSPH�N���X�̐錾
//  - Miles Macklin and Matthias Muller, "Position Based Fluids", Proc. SIGGRAPH 2013, 2013. 
//  - http://blog.mmacklin.com/publications/
//-----------------------------------------------------------------------------
class rxPBDSPH : public rxParticleSystemBase
{
private:
	// �p�[�e�B�N��
	RXREAL *m_hNrm;					//!< �p�[�e�B�N���@��
	RXREAL *m_hFrc;					//!< �p�[�e�B�N���ɂ������
	RXREAL *m_hDens;				//!< �p�[�e�B�N�����x
	RXREAL *m_hPres;				//!< �p�[�e�B�N������

	RXREAL *m_hS;					//!< Scaling factor for CFM
	RXREAL *m_hDp;					//!< �ʒu�C����

	RXREAL *m_hPredictPos;			//!< �\���ʒu
	RXREAL *m_hPredictVel;			//!< �\�����x

	RXREAL *m_hSb;					//!< ���E�p�[�e�B�N����Scaling factor

	// �\�ʐ����p(Anisotropic kernel)
	RXREAL *m_hUpPos;				//!< �������p�[�e�B�N���ʒu
	RXREAL *m_hPosW;				//!< �d�ݕt�����ύ��W
	RXREAL *m_hCMatrix;				//!< �����U�s��
	RXREAL *m_hEigen;				//!< �����U�s��̓��ْl
	RXREAL *m_hRMatrix;				//!< ��]�s��(�����U�s��̓��كx�N�g��)
	RXREAL *m_hG;					//!< �ό`�s��

	uint *m_hSurf;					//!< �\�ʃp�[�e�B�N��

	// ���E�E�ő�
	rxSolid *m_pBoundary;			//!< �V�~�����[�V������Ԃ̋��E
	vector<rxSolid*> m_vSolids;		//!< �ő̕���
	RXREAL *m_hVrts;				//!< �ő̃|���S���̒��_
	int m_iNumVrts;					//!< �ő̃|���S���̒��_��
	int *m_hTris;					//!< �ő̃|���S��
	int m_iNumTris;					//!< �ő̃|���S���̐�
	RXREAL *m_hSVels;				//!< �ő̃|���S�����x

	// ��ԕ����i�q�֘A
	rxNNGrid *m_pNNGrid;			//!< �����O���b�h�ɂ��ߖT�T��
	vector< vector<rxNeigh> > m_vNeighs;	//!< �ߖT�p�[�e�B�N��

	rxNNGrid *m_pNNGridB;			//!< ���E�p�[�e�B�N���p�����O���b�h


	// ���q�p�����[�^
	uint m_iKernelParticles;		//!< �J�[�l�����̃p�[�e�B�N����
	RXREAL m_fRestDens, m_fMass;	//!< ���x�C����
	RXREAL m_fEffectiveRadius;		//!< �L�����a
	RXREAL m_fKernelRadius;			//!< �J�[�l���̉e���͈�

	// �V�~�����[�V�����p�����[�^
	RXREAL m_fGasStiffness;			//!< �K�X�萔
	RXREAL m_fViscosity;			//!< �S���W��
	RXREAL m_fBuoyancy;				//!< ����

	RXREAL m_fEpsilon;				//!< CFM�̊ɘa�W��
	RXREAL m_fEta;					//!< ���x�ϓ���
	int m_iMinIterations;			//!< ���R�r�����ŏ�������
	int m_iMaxIterations;			//!< ���R�r�����ő唽����

	bool m_bArtificialPressure;		//!< �N���X�^�����O��h�����߂�Artificial Pressure����ǉ�����t���O
	RXREAL m_fApK;					//!< �l�H���͂̂��߂̌W��k
	RXREAL m_fApN;					//!< �l�H���͂̂��߂̌W��n
	RXREAL m_fApQ;					//!< �l�H���͌v�Z���̊�J�[�l���l�v�Z�p�W��(�L�����ah�ɑ΂���W��, [0,1])


	// �J�[�l���֐��̌v�Z�̍ۂɗp������萔�W��
	double m_fAw;
	double m_fAg;
	double m_fAl;

	// �J�[�l���֐�
	double (*m_fpW)(double, double, double);
	Vec3 (*m_fpGW)(double, double, double, Vec3);
	double (*m_fpLW)(double, double, double, double);

protected:
	rxPBDSPH(){}

public:
	//! �R���X�g���N�^
	rxPBDSPH(bool use_opengl);

	//! �f�X�g���N�^
	virtual ~rxPBDSPH();

	// �p�[�e�B�N�����a
	float GetEffectiveRadius(){ return m_fEffectiveRadius; }

	// �ߖT�p�[�e�B�N��
	uint* GetNeighborList(const int &i, int &n);

public:
	//
	// ���z�֐�
	//
	// �p�[�e�B�N���f�[�^
	virtual RXREAL* GetParticle(void);
	virtual RXREAL* GetParticleDevice(void);
	virtual void UnmapParticle(void){}

	// �V�~�����[�V�����X�e�b�v
	virtual bool Update(RXREAL dt, int step = 0);

	// �V�[���̐ݒ�
	virtual void SetPolygonObstacle(const vector<Vec3> &vrts, const vector<Vec3> &nrms, const vector< vector<int> > &tris, Vec3 vel);
	virtual void SetBoxObstacle(Vec3 cen, Vec3 ext, Vec3 ang, Vec3 vel, int flg);
	virtual void SetSphereObstacle(Vec3 cen, double rad, Vec3 vel, int flg);
	virtual void MoveSphereObstacle(int b, Vec3 disp);
	virtual Vec3 GetSphereObstaclePos(int b = -1);

	// �z�X�g<->VBO�ԓ]��
	virtual RXREAL* GetArrayVBO(rxParticleArray type, bool d2h = true, int num = -1);
	virtual void SetArrayVBO(rxParticleArray type, const RXREAL* data, int start, int count);
	virtual void SetArrayVBO(rxParticleArray type, const int* data, int start, int count);
	virtual void SetColorVBO(int type);

	// �A�֐��l�v�Z
	virtual double GetImplicit(double x, double y, double z);
	virtual void CalImplicitField(int n[3], Vec3 minp, Vec3 d, RXREAL *hF);
	virtual void CalImplicitFieldDevice(int n[3], Vec3 minp, Vec3 d, RXREAL *dF);

	// SPH���o��
	virtual void OutputSetting(string fn);

	// �`��֐�
	virtual void DrawCell(int i, int j, int k);
	virtual void DrawCells(Vec3 col, Vec3 col2, int sel = 0);
	virtual void DrawObstacles(void);


public:
	// SPH������
	void Initialize(const rxSPHEnviroment &env);
	void Allocate(int max_particles);
	void Finalize(void);

	// �ߖT�擾
	void GetNearestNeighbors(Vec3 pos, vector<rxNeigh> &neighs, RXREAL h = -1.0);
	void GetNearestNeighbors(int idx, RXREAL *p, vector<rxNeigh> &neighs, RXREAL h = -1.0);

	void GetNearestNeighborsB(Vec3 pos, vector<rxNeigh> &neighs, RXREAL h = -1.0);

	// �����Z���Ƀp�[�e�B�N�����i�[
	virtual void SetParticlesToCell(void);
	virtual void SetParticlesToCell(RXREAL *prts, int n, RXREAL h);

	// ���^�{�[���ɂ��A�֐��l
	double CalColorField(double x, double y, double z);

	// �����Z���Ƀ|���S�����i�[
	virtual void SetPolygonsToCell(void);

	int  GetPolygonsInCell(int gi, int gj, int gk, vector<int> &polys);
	bool IsPolygonsInCell(int gi, int gj, int gk);

	
	// �\�ʃp�[�e�B�N�����o
	void DetectSurfaceParticles(void);					// �\�ʃp�[�e�B�N���̌��o
	double CalDistToNormalizedMassCenter(const int i);	// �ߖT�p�[�e�B�N���̏d�S�܂ł̋����v�Z
	uint* GetArraySurf(void);							// �\�ʃp�[�e�B�N�����̎擾

	// �\�ʃp�[�e�B�N�����̎擾
	int GetSurfaceParticles(const Vec3 pos, RXREAL h, vector<rxSurfaceParticle> &sp);

	// �@���v�Z
	void CalNormalFromDensity(void);
	void CalNormal(void);

	//�ǉ��F�F�ߖT���q�̎擾
	vector<vector<rxNeigh>>& GetNeights(void){ return m_vNeighs; }

	// �l�����͍�
	bool& GetArtificialPressure(void){ return m_bArtificialPressure; }

	//
	// Anisotropic kernel
	//
	virtual void CalAnisotropicKernel(void);

	//
	// SPS����
	//
	int	GetMaxSubLevel() const { return 0; }

	// �T�u�p�[�e�B�N���f�[�^(�f�o�C�X������)
	RXREAL* GetSubParticleDev(void){ return 0; }
	void    UnmapSubParticle(void){}
	RXREAL* GetAllSubParticleDev(void){ return 0; }
	void    UnmapAllSubParticle(void){}
	RXREAL* GetSubParticlePosDev(void){ return 0; }
	RXREAL* GetSubParticleRadDev(void){ return 0; }
	RXREAL* GetSubParticleRatDev(void){ return 0; }

	int GetNumAllSubParticles(void){ return 0; }
	int GetNumValidSubParticles(void){ return 0; }

	// �T�u�p�[�e�B�N���f�[�^(�z�X�g������)
	RXREAL* GetSubParticlePos(void){ return 0; }
	RXREAL* GetSubParticleRad(void){ return 0; }
	RXREAL* GetSubParticleRat(void){ return 0; }
	unsigned int* GetSubParticleOcc(void){ return 0; }
	unsigned int* GetSubParticleOccScan(void){ return 0; }

	RXREAL GetSubRadius(int level){ return 0.0f; }

	//�T�u���q
	void SetEtcri(const double &et_cri){}
	void InitSubParticles(const double &et_cri){}
	void AddSubParticles(int start, int count){}
	void UpdateSubParticle(RXREAL *dPos, RXREAL *dSubPos, RXREAL scale, RXREAL dt){}
	void CalRadiusAndRatio(bool force = false){}

	void CalSubImplicitFieldDevice(int n[3], Vec3 minp, Vec3 d, RXREAL *dF, int turb){}
	void CalSubImplicitField(int n[3], Vec3 minp, Vec3 d, RXREAL *hF, const double &et_cri){}

	void AddSphereField(int n[3], Vec3 minp, Vec3 d, RXREAL *hF, 
						const Vec3 &pos, const double &eRad, const double& ratio){}



protected:
	// CPU�ɂ��SPH�v�Z
	void calDensity(const RXREAL *pos, RXREAL *dens, RXREAL h);
	void calNormal(void);
	void calForceExtAndVisc(const RXREAL *pos, const RXREAL *vel, const RXREAL *dens, RXREAL *frc, RXREAL h);

	// Scaling factor�̌v�Z
	void calScalingFactor(const RXREAL *ppos, RXREAL *pdens, RXREAL *pscl, RXREAL h, RXREAL dt);

	// Scaling factor�̌v�Z
	void calPositionCorrection(const RXREAL *ppos, const RXREAL *pscl, RXREAL *pdp, RXREAL h, RXREAL dt);

	// rest density�̌v�Z
	RXREAL calRestDensity(RXREAL h);

	// �X�̋��E�p�[�e�B�N���̑̐ς��v�Z
	void calBoundaryVolumes(const RXREAL *bpos, RXREAL *bvol, RXREAL mass, uint n, RXREAL h);

	// ���ԃX�e�b�v���̏C��
	RXREAL calTimeStep(RXREAL &dt, RXREAL eta_avg, const RXREAL *pfrc, const RXREAL *pvel, const RXREAL *pdens);
	// �ʒu�Ƒ��x�̍X�V(Leap-Frog)
	void integrate(const RXREAL *pos, const RXREAL *vel, const RXREAL *dens, const RXREAL *acc, 
				   RXREAL *pos_new, RXREAL *vel_new, RXREAL dt);

	// �Փ˔���
	int calCollisionPolygon(uint grid_hash, Vec3 &pos0, Vec3 &pos1, Vec3 &vel, RXREAL dt);
	int calCollisionSolid(Vec3 &pos0, Vec3 &pos1, Vec3 &vel, RXREAL dt);
};




//-----------------------------------------------------------------------------
// MARK:rxPBDSPH_GPU�N���X�̐錾
//  - Miles Macklin and Matthias Muller, "Position Based Fluids", Proc. SIGGRAPH 2013, 2013. 
//  - http://blog.mmacklin.com/publications/
//-----------------------------------------------------------------------------
class rxPBDSPH_GPU : public rxParticleSystemBase
{
private:
	//
	// �����o�ϐ�(GPU�ϐ�)
	//
	RXREAL *m_dPos;			//!< �p�[�e�B�N���ʒu
	RXREAL *m_dVel;			//!< �p�[�e�B�N�����x
	RXREAL *m_dNrm;			//!< �p�[�e�B�N���@��
	RXREAL *m_dFrc;			//!< �p�[�e�B�N���ɂ������
	RXREAL *m_dDens;		//!< �p�[�e�B�N�����x
	RXREAL *m_dPres;		//!< �p�[�e�B�N������

	RXREAL *m_dPosB;		//!< ���E�p�[�e�B�N��
	RXREAL *m_dVolB;		//!< ���E�p�[�e�B�N���̑̐�

	RXREAL *m_dS;			//!< Scaling factor for CFM
	RXREAL *m_dDp;			//!< �ʒu�C����
	RXREAL *m_dPredictPos;	//!< �\���ʒu
	RXREAL *m_dPredictVel;	//!< �\�����x

	RXREAL *m_dSb;			//!< ���E�p�[�e�B�N����Scaling factor

	RXREAL *m_dErr;			//!< ���x�ϓ��l	
	RXREAL *m_dErrScan;		//!< ���x�ϓ��l��Scan����

	int *m_dAttr;			//!< �p�[�e�B�N������

	cudaGraphicsResource *m_pPosResource;	//!< OpenGL(��VBO)-CUDA�Ԃ̃f�[�^�]�����������߂̃n���h��

	// �E�F�[�u���b�g����
	RXREAL *m_dEt;			//!< ���x�G�l���M�[�X�y�N�g��
	RXREAL *m_dTurb;		//!< �E�F�[�u���b�g�����ɂ�鑬�x��

	//�T�u���q GPU data
	RXREAL *m_dSubPos;		//!< �T�u���q�̐�΍��W float4(0.0,x,y,z)
	RXREAL *m_dSubChild;	//!< �T�u���q�̎q1�ւ̒P�ʃx�N�g��
	RXREAL *m_dSubAxis;		//!< �T�u���q�̉�]��(�P�ʃx�N�g��)
	RXREAL *m_dSubEt;		//!< �T�u���q�̃G�l���M�[�X�y�N�g��
	uint   *m_dSubRand;		//!< �����e�[�u��

	uint m_subposVBO;		//!< �T�u���q���WVBO
	cudaGraphicsResource *m_pSubPosResource;	//!< OpenGL(��VBO)-CUDA�Ԃ̃f�[�^�]�����������߂̃n���h��

	// SPS����
	RXREAL *m_dUpPos;		//!< �������p�[�e�B�N���ʒu
	RXREAL *m_dPosW;		//!< �d�ݕt�����ύ��W
	RXREAL *m_dCMatrix;		//!< �����U�s��
	RXREAL *m_dEigen;		//!< �����U�s��̓��ْl
	RXREAL *m_dRMatrix;		//!< ��]�s��(�����U�s��̓��كx�N�g��)
	RXREAL *m_dG;			//!< �ό`�s��

	// �\�ʃ��b�V��
	RXREAL *m_dVrts;		//!< �ő̃��b�V�����_
	int    *m_dTris;		//!< �ő̃��b�V��

	// �V�~�����[�V�����p�����[�^
	rxSimParams m_params;	//!< �V�~�����[�V�����p�����[�^(GPU�ւ̃f�[�^�n���p)
	uint3 m_gridSize;		//!< �ߖT�T���O���b�h�̊e���̕�����
	uint m_numGridCells;	//!< �ߖT�T���O���b�h��������
	
	// ��ԕ���(GPU)
	rxParticleCell m_dCellData;			//!< �ߖT�T���O���b�h
	rxSubParticleCell m_dSubCellData;	//!< �T�u���q�p�ߖT�T���O���b�h
	uint m_gridSortBits;				//!< �n�b�V���l�ɂ���\�[�g���̊����

	rxParticleCell m_dCellDataB;		//!< ���E�p�[�e�B�N���p�ߖT�T���O���b�h
	uint3 m_gridSizeB;					//!< ���E�p�[�e�B�N���p�ߖT�T���O���b�h�̊e���̕�����

	//
	// �����o�ϐ�(CPU�ϐ�)
	//
	RXREAL *m_hNrm;			//!< �p�[�e�B�N���@��
	RXREAL *m_hFrc;			//!< �p�[�e�B�N���ɂ������
	RXREAL *m_hDens;		//!< �p�[�e�B�N�����x
	RXREAL *m_hPres;		//!< �p�[�e�B�N������

	RXREAL *m_hS;			//!< Scaling factor for CFM
	RXREAL *m_hDp;			//!< �ʒu�C����
	RXREAL *m_hPredictPos;	//!< �\���ʒu
	RXREAL *m_hPredictVel;	//!< �\�����x

	RXREAL *m_hSb;			//!< ���E�p�[�e�B�N����Scaling factor

	uint *m_hSurf;			//!< �\�ʃp�[�e�B�N��

	// �\�ʐ����p(Anisotropic kernel)
	RXREAL *m_hUpPos;		//!< �������p�[�e�B�N���ʒu
	RXREAL *m_hPosW;		//!< �d�ݕt�����ύ��W
	RXREAL *m_hCMatrix;		//!< �����U�s��
	RXREAL *m_hEigen;		//!< �����U�s��̓��ْl
	RXREAL *m_hRMatrix;		//!< ��]�s��(�����U�s��̓��كx�N�g��)
	RXREAL *m_hG;			//!< �ό`�s��
	RXREAL  m_fEigenMax;	//!< ���ْl�̍ő�l(�T�����a�g���ɗp����)

	// �E�F�[�u���b�g����
	RXREAL *m_hVwt;			//!< �p�[�e�B�N�����x�̃E�F�[�u���b�g�ϊ�
	RXREAL *m_hEt;			//!< ���x�̃G�l���M�[�X�y�N�g����
	RXREAL *m_hTurb;		//!< �E�F�[�u���b�g�����ɂ�鑬�x��
	RXREAL *m_hNoiseTile[3];//!< �E�F�[�u���b�g�m�C�Y�^�C��
	RXREAL m_fNoiseTileWidth;	//!< �m�C�Y�^�C���̕�
	uint   m_iNoiseTileN[3];	//!< �m�C�Y�^�C���̉𑜓x

	// SPS����
	RXREAL *m_hSubPos;		//!< �T�u���q�̐�΍��W
	RXREAL *m_hSubChild;	//!< �T�u���q�̎q1�ւ̒P�ʃx�N�g��
	RXREAL *m_hSubAxis;		//!< �T�u���q�̉�]��(�P�ʃx�N�g��)
	RXREAL *m_hSubEt;		//!< �T�u���q�̃G�l���M�[�X�y�N�g��

	uint m_uNumSubParticles;//!< �T�u���q��
	uint m_uMaxSubParticles;//!< �T�u���q�ő吔
	uint m_uMaxSubLevel;	//!< �T�u���q�ő僌�x��

	RXREAL *m_hSubPosPack;
	RXREAL *m_hSubRadPack;
	RXREAL *m_hSubRatPack;
	unsigned int *m_hSubOcc;
	unsigned int *m_hSubOccScan;
	bool m_bSubPacked;		//!< �T�u�p�[�e�B�N����Pack������true

	// �\�ʃ��b�V��
	vector<RXREAL> m_vVrts;	//!< �ő̃��b�V�����_
	int m_iNumVrts;			//!< �ő̃��b�V�����_��
	vector<int> m_vTris;	//!< �ő̃��b�V��
	int m_iNumTris;			//!< �ő̃��b�V����

	bool m_bCalNormal;		//!< �@���v�Z�t���O

	RXREAL m_fEpsilon;				//!< CFM�̊ɘa�W��
	RXREAL m_fEta;					//!< ���x�ϓ���
	int m_iMinIterations;			//!< ���R�r�����ŏ�������
	int m_iMaxIterations;			//!< ���R�r�����ő唽����

	bool m_bArtificialPressure;		//!< �N���X�^�����O��h�����߂�Artificial Pressure����ǉ�����t���O

	//�ǉ��F�ߖT���q��GPU����̃R�s�[�f�[�^
	uint* m_hSortedIndex;		//!< �n�b�V���l�Ń\�[�g�����p�[�e�B�N���C���f�b�N�X
	uint* m_hGridParticleHash;	//!< �e�p�[�e�B�N���̃O���b�h�n�b�V���l
	uint* m_hCellStart;			//!< �\�[�g���X�g���̊e�Z���̃X�^�[�g�C���f�b�N�X
	uint* m_hCellEnd;			//!< �\�[�g���X�g���̊e�Z���̃G���h�C���f�b�N�X
	uint  m_uNumCells;			//!< ���Z����
	vector< vector<rxNeigh> > m_vNeighs;	//!< �ߖT�p�[�e�B�N��

protected:
	rxPBDSPH_GPU(){}

public:
	//! �R���X�g���N�^
	rxPBDSPH_GPU(bool use_opengl);

	//! �f�X�g���N�^
	~rxPBDSPH_GPU();

	// �p�[�e�B�N�����a
	float GetEffectiveRadius(){ return m_params.EffectiveRadius; }

	// �ߖT�p�[�e�B�N��
	uint* GetNeighborList(const int &i, int &n);

	// �V�~�����[�V�����p�����[�^
	rxSimParams GetParams(void){ return m_params; }
	void UpdateParams(void);

	// �t���O�ؑ�
	void ToggleNormalCalc(int t = -1){ RX_TOGGLE(m_bCalNormal, t); }			//!< �p�[�e�B�N���@���̌v�Z
	bool IsNormalCalc(void) const { return m_bCalNormal; }
#if MAX_BOX_NUM
	void ToggleSolidFlg(int t = -1);
#endif

public:
	//
	// ���z�֐�
	//
	// �p�[�e�B�N���f�[�^
	virtual RXREAL* GetParticle(void);
	virtual RXREAL* GetParticleDevice(void);
	virtual void UnmapParticle(void);

	// �V�~�����[�V�����X�e�b�v
	virtual bool Update(RXREAL dt, int step = 0);

	// �V�[���̐ݒ�
	virtual void SetPolygonObstacle(const vector<Vec3> &vrts, const vector<Vec3> &nrms, const vector< vector<int> > &tris, Vec3 vel);
	virtual void SetBoxObstacle(Vec3 cen, Vec3 ext, Vec3 ang, Vec3 vel, int flg);
	virtual void SetSphereObstacle(Vec3 cen, double rad, Vec3 vel, int flg);
	virtual void MoveSphereObstacle(int b, Vec3 disp);
	virtual Vec3 GetSphereObstaclePos(int b = -1);

	// �z�X�g<->VBO�ԓ]��
	virtual RXREAL* GetArrayVBO(rxParticleArray type, bool d2h = true, int num = -1);
	virtual void SetArrayVBO(rxParticleArray type, const RXREAL* data, int start, int count);
	virtual void SetArrayVBO(rxParticleArray type, const int* data, int start, int count);
	virtual void SetColorVBO(int type);


	// �A�֐��l�v�Z
	virtual double GetImplicit(double x, double y, double z);
	virtual void CalImplicitField(int n[3], Vec3 minp, Vec3 d, RXREAL *hF);
	virtual void CalImplicitFieldDevice(int n[3], Vec3 minp, Vec3 d, RXREAL *dF);

	// SPH���o��
	virtual void OutputSetting(string fn);

	// �`��֐�
	virtual void DrawCell(int i, int j, int k);
	virtual void DrawCells(Vec3 col, Vec3 col2, int sel = 0);
	virtual void DrawObstacles(void);

protected:
	void setObjectToCell(RXREAL *p, RXREAL *v);


public:
	// SPH������
	void Initialize(const rxSPHEnviroment &env);
	void Allocate(int max_particles);
	void Finalize(void);

	// �����Z���Ƀp�[�e�B�N�����i�[
	virtual void SetParticlesToCell(void);
	virtual void SetParticlesToCell(RXREAL *prts, int n, RXREAL h);
	virtual void SetPolygonsToCell(void);

	int  GetPolygonsInCell(int gi, int gj, int gk, vector<int> &polys){ return 0; }
	bool IsPolygonsInCell(int gi, int gj, int gk){ return false; }

	void CalMaxDensity(int k);

	// �\�ʃp�[�e�B�N�����o
	void DetectSurfaceParticles(void);					// �\�ʃp�[�e�B�N���̌��o
	double CalDistToNormalizedMassCenter(const int i);	// �ߖT�p�[�e�B�N���̏d�S�܂ł̋����v�Z
	uint* GetArraySurf(void);							// �\�ʃp�[�e�B�N�����̎擾

	// �\�ʃp�[�e�B�N�����̎擾
	int GetSurfaceParticles(const Vec3 pos, RXREAL h, vector<rxSurfaceParticle> &sp);

	// �@���v�Z
	void CalNormalFromDensity(void);
	void CalNormal(void);

	// �l�����͍�
	bool& GetArtificialPressure(void){ return m_bArtificialPressure; }

	// ���E�p�[�e�B�N���̏�����
	void InitBoundary(void);

	//�ǉ��F�F�ߖT���q�擾
	vector<vector<rxNeigh>>& GetNeights(void){ return m_vNeighs; }

	//
	// Anisotropic kernel
	//
	virtual void CalAnisotropicKernel(void);

	//
	// �E�F�[�u���b�g����
	//
	void CalParticleCWT(RXREAL scale);
	void CalWaveletTurbulence(RXREAL scale, RXREAL *dPos, RXREAL dt);

	//
	// SPS����
	//
	int	GetMaxSubLevel() const { return m_uMaxSubLevel; }

	// �T�u�p�[�e�B�N���f�[�^(�f�o�C�X������)
	RXREAL* GetSubParticleDev(void);
	void    UnmapSubParticle(void);
	RXREAL* GetAllSubParticleDev(void);
	void    UnmapAllSubParticle(void);
	RXREAL* GetSubParticlePosDev(void);
	RXREAL* GetSubParticleRadDev(void);
	RXREAL* GetSubParticleRatDev(void);

	int GetNumAllSubParticles(void);
	int GetNumValidSubParticles(void);

	// �T�u�p�[�e�B�N���f�[�^(�z�X�g������)
	RXREAL* GetSubParticlePos(void);
	RXREAL* GetSubParticleRad(void);
	RXREAL* GetSubParticleRat(void);
	unsigned int* GetSubParticleOcc(void);
	unsigned int* GetSubParticleOccScan(void);

	RXREAL GetSubRadius(int level){ return m_dSubCellData.fSubRad[level]; }

	//�T�u���q
	void SetEtcri(const double &et_cri){ m_dSubCellData.fEtcri = et_cri;}
	void InitSubParticles(const double &et_cri);
	void AddSubParticles(int start, int count);
	void UpdateSubParticle(RXREAL *dPos, RXREAL *dSubPos, RXREAL scale, RXREAL dt);
	void CalRadiusAndRatio(bool force = false);

	void CalSubImplicitFieldDevice(int n[3], Vec3 minp, Vec3 d, RXREAL *dF, int turb);
	void CalSubImplicitField(int n[3], Vec3 minp, Vec3 d, RXREAL *hF, const double &et_cri);

	void AddSphereField(int n[3], Vec3 minp, Vec3 d, RXREAL *hF, 
						const Vec3 &pos, const double &eRad, const double& ratio);

protected:
	// rest density�̌v�Z
	RXREAL calRestDensity(RXREAL h);

	// �����Z���̏����ݒ�
	void setupCells(rxParticleCell &cell, uint3 &gridsize, double &cell_width, Vec3 vMin, Vec3 vMax, double h);

	// �O���b�h�n�b�V���̌v�Z
	uint calGridHash(uint x, uint y, uint z);
	uint calGridHash(Vec3 pos);

	// �|���S���𕪊��Z���Ɋi�[
	void setPolysToCell(RXREAL *vrts, int nv, int* tris, int nt);

	// CPU�p�̋ߖT���q�T���@�ǉ��F�FGPU��CPU�̃f�[�^�]��
	void searchNeighbors(void);
	void getNeighborsInCell(Vec3 pos, RXREAL *p, int gi, int gj, int gk, vector<rxNeigh> &neighs, RXREAL h);
};




#endif	// _SPH_H_

