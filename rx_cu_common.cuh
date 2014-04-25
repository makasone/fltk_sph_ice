/*! 
  @file rx_cu_common.cuh
	
  @brief CUDA���ʃw�b�_
 
  @author Makoto Fujisawa
  @date 2009-08, 2011-06
*/
// FILE --rx_cu_common.cuh--

#ifndef _RX_CU_COMMON_CUH_
#define _RX_CU_COMMON_CUH_


//-----------------------------------------------------------------------------
// �C���N���[�h�t�@�C��
//-----------------------------------------------------------------------------
#include "vector_types.h"
#include "vector_functions.h"


//-----------------------------------------------------------------------------
// ��`
//-----------------------------------------------------------------------------


typedef unsigned int uint;
typedef unsigned char uchar;

#define FLOAT float
#define FLOAT3 float3

#define MAKE_FLOAT3 make_float3

#define RX_CUMC_USE_GEOMETRY

#define RX_USE_ATOMIC_FUNC // �vCompute capability 1.1�ȏ�(-arch sm_11)


// �e�N�X�`���������̎g�p�t���O
#ifndef __DEVICE_EMULATION__
#define USE_TEX 0
#endif

#if USE_TEX
#define FETCH(t, i) tex1Dfetch(t##Tex, i)
#else
#define FETCH(t, i) t[i]
#endif

#if USE_TEX
#define FETCHC(t, i) tex1Dfetch(t##Tex, i)
#else
#define FETCHC(t, i) cell.t[i]
#endif


// 1�u���b�N������̃X���b�h��(/����)
#define BLOCK_SIZE 16

// 1�u���b�N������̍ő�X���b�h��
#define THREAD_NUM 256

// �T���v���{�����[����p����ꍇ1, �A�֐���p����ꍇ��0
#define SAMPLE_VOLUME 1

// Shared Memory�̎g�p�t���O
#define USE_SHARED 1

// Shared Memory��p����ꍇ�̃X���b�h���̐���
#define NTHREADS 32

#define SKIP_EMPTY_VOXELS 1

// CWT
#define M_PI 3.141592653589793238462643383279502884197
#define M_2PI 6.283185307179586476925286766559005768394		// 2*PI
#define M_SQRTPI 1.772453850905516027298167483341145182797	// sqrt(PI)
#define M_SQRT2PI 2.506628274631000502415765284811045253006	// sqrt(2*PI)

#define MEXICAN_HAT_FC (1.0/M_PI)
#define MEXICAN_HAT_C 0.8673250705840776f // c =  2/(sqrt(3)*pi^(1/4))
#define MEXICAN_HAT_R 5.0f

#define RX_MAX_FILTER_SIZE 10
#define RX_BINOMIALS_SIZE (RX_MAX_FILTER_SIZE+1)*(RX_MAX_FILTER_SIZE+1)

//�T�u�p�[�e�B�N���Ɋւ���萔
#define MAX_SUB_LEVEL (3)
#define MAX_SUB_NUM (15) //2^(MAX_SUB_LEVEL+1)-1
#define RATIO_AXIS_ANGLE (0.50) //���̂Ԃ�Ɋւ�����萔
#define POW_2_M1D3 (0.793700526) //pow(2.0,-(1.0/3.0))
#define POW_2_M5D3 (0.314980262) //pow(2.0,-(5.0/3.0))
#define POW_2_M5D6 (0.561231024) //pow(2.0,-(5.0/6.0))
#define POW_2_M5D9 (0.680395000) //pow(2.0,-(5.0/9.0))
#define INV_LOG_2_M5D9 (-5.979470571) // 1/log(2^(-5/9));
//#define FLT_MAX         3.402823466e+38F        // max value 


// �����Z���k�c�C�X�^�[
#define      DCMT_SEED 4172
#define  MT_RNG_PERIOD 607

typedef struct{
	unsigned int matrix_a;
	unsigned int mask_b;
	unsigned int mask_c;
	unsigned int seed;
} mt_struct_stripped;

#define   MT_RNG_COUNT 4096
#define          MT_MM 9
#define          MT_NN 19
#define       MT_WMASK 0xFFFFFFFFU
#define       MT_UMASK 0xFFFFFFFEU
#define       MT_LMASK 0x1U
#define      MT_SHIFT0 12
#define      MT_SHIFTB 7
#define      MT_SHIFTC 15
#define      MT_SHIFT1 18

// �s��
struct matrix3x3
{
	float3 e[3];
};

struct matrix4x4
{
	float4 e[4];
};


#define MAX_POLY_NUM 10
#define MAX_BOX_NUM 10
#define MAX_SPHERE_NUM 10


//-----------------------------------------------------------------------------
//! SPH�V�~�����[�V�����p�����[�^
//-----------------------------------------------------------------------------
struct rxSimParams
{
	FLOAT3 Gravity;
	FLOAT GlobalDamping;
	FLOAT ParticleRadius;

	float3 BoundaryMax;
	float3 BoundaryMin;

	uint3 GridSize;
	uint NumCells;
	FLOAT3 WorldOrigin;
	FLOAT3 WorldMax;
	FLOAT3 CellWidth;

	uint NumBodies;
	uint MaxParticlesPerCell;

	uint3 GridSizeB;
	uint NumCellsB;
	FLOAT3 WorldOriginB;
	FLOAT3 WorldMaxB;
	FLOAT3 CellWidthB;
	uint NumBodiesB;


	FLOAT EffectiveRadius;
	FLOAT Mass;			// �p�[�e�B�N������[kg]
	FLOAT VorticityConfinement;

	FLOAT Buoyancy;

	FLOAT Density;		// ���x[kg/m^3]
	FLOAT Pressure;     // [Pa = N/m^2 = kg/m.s^2]

	FLOAT Tension;		// [N/m = kg/s^2]
	FLOAT Viscosity;	// [Pa.s = N.s/m^2 = kg/m.s]
	FLOAT GasStiffness;	// [J = N.m = kg.m^2/s^2]  // used for DC96 symmetric pressure force

	FLOAT Volume;
	FLOAT KernelParticles;
	FLOAT Restitution;

	FLOAT Threshold;

	FLOAT InitDensity;

	FLOAT Dt;

	FLOAT Wpoly6;		//!< Pory6�J�[�l���̒萔�W��
	FLOAT GWpoly6;		//!< Pory6�J�[�l���̌��z�̒萔�W��
	FLOAT LWpoly6;		//!< Pory6�J�[�l���̃��v���V�A���̒萔�W��
	FLOAT Wspiky;		//!< Spiky�J�[�l���̒萔�W��
	FLOAT GWspiky;		//!< Spiky�J�[�l���̌��z�̒萔�W��
	FLOAT LWspiky;		//!< Spiky�J�[�l���̃��v���V�A���̒萔�W��
	FLOAT Wvisc;		//!< Viscosity�J�[�l���̒萔�W��
	FLOAT GWvisc;		//!< Viscosity�J�[�l���̌��z�̒萔�W��
	FLOAT LWvisc;		//!< Viscosity�J�[�l���̃��v���V�A���̒萔�W��

	FLOAT Wd2;		//!< 
	FLOAT Wd3;		//!< 
	FLOAT GWd2;		//!< 
	FLOAT GWd3;		//!< 
	FLOAT Wd1;		//!< 

	FLOAT NoiseScl;		//!< �m�C�Y�t�����̃X�P�[��
	FLOAT NoiseEthr;	//!< �m�C�Y�t�����̃G�l���M�[�X�y�N�g��臒l
	FLOAT NoiseMag;		//!< �m�C�Y�t�����̃G�l���M�[�X�y�N�g���W��

	int AP;
	FLOAT AP_K;
	FLOAT AP_N;
	FLOAT AP_Q;
	FLOAT AP_WQ;

	uint   PolyNum;
	FLOAT3 PolyVel[MAX_POLY_NUM];

	uint   BoxNum;
#if MAX_BOX_NUM
	FLOAT3 BoxCen[MAX_BOX_NUM];
	FLOAT3 BoxExt[MAX_BOX_NUM];
	//FLOAT3 BoxVel[MAX_BOX_NUM];
	matrix3x3 BoxRot[MAX_BOX_NUM];
	matrix3x3 BoxInvRot[MAX_BOX_NUM];
	uint   BoxFlg[MAX_BOX_NUM];
#endif

	uint SphereNum;
#if MAX_SPHERE_NUM
	FLOAT3 SphereCen[MAX_SPHERE_NUM];
	FLOAT  SphereRad[MAX_SPHERE_NUM];
	//FLOAT3 SphereVel[MAX_SPHERE_NUM];
	uint   SphereFlg[MAX_SPHERE_NUM];
#endif
};



struct rxParticleCell
{
	float4* dSortedPos;			//!< �\�[�g�ς݃p�[�e�B�N�����W
	float4* dSortedVel;			//!< �\�[�g�ς݃p�[�e�B�N�����x

	uint* dSortedIndex;			//!< �\�[�g�ς݃p�[�e�B�N���C���f�b�N�X
	uint* dGridParticleHash;	//!< �e�p�[�e�B�N���̃O���b�h�n�b�V���l(�\�[�g�p�L�[)
	uint* dCellStart;			//!< �\�[�g���X�g���̊e�Z���̃X�^�[�g�C���f�b�N�X
	uint* dCellEnd;				//!< �\�[�g���X�g���̊e�Z���̃G���h�C���f�b�N�X
	uint  uNumParticles;		//!< ���p�[�e�B�N����
	uint  uNumCells;			//!< ���Z����
	uint  uNumArdGrid;			//!< �ߖT�T�����Q�ƃO���b�h�͈�

	uint* dSortedPolyIdx;		//!< �\�[�g�ς݃|���S���C���f�b�N�X
	uint* dGridPolyHash;		//!< �|���S���̃O���b�h�n�b�V���l(�\�[�g�p�L�[)
	uint* dPolyCellStart;		//!< �\�[�g���X�g���̊e�Z���̃X�^�[�g�C���f�b�N�X
	uint* dPolyCellEnd;			//!< �\�[�g���X�g���̊e�Z���̃G���h�C���f�b�N�X
	uint  uNumPolyHash;
};


//-----------------------------------------------------------------------------
//! �T�u�p�[�e�B�N��
//-----------------------------------------------------------------------------
struct rxSubParticleCell
{
	float4* dSubUnsortPos;	//!< �T�u�p�[�e�B�N�����W
	float*	dSubUnsortRad;	//!< �T�u�p�[�e�B�N�����a
	float*	dSubUnsortRat;	//!< �T�u�p�[�e�B�N���e���W��
	float4* dSubSortedPos;	//!< �T�u�p�[�e�B�N�����W(�O���b�h�n�b�V���Ń\�[�g)
	float*	dSubSortedRad;	//!< �T�u�p�[�e�B�N�����a(�O���b�h�n�b�V���Ń\�[�g)
	float*	dSubSortedRat;	//!< �T�u�p�[�e�B�N���e���W��(�O���b�h�n�b�V���Ń\�[�g)

	uint*   dSubOcc;		//!< �T�u�p�[�e�B�N���L��/����
	uint*   dSubOccScan;	//!< �T�u�p�[�e�B�N���L��/������Scan

	uint*	dSubSortedIndex;
	uint*	dSubGridParticleHash;	//!< �e�T�u�p�[�e�B�N���̃O���b�h�n�b�V���l
	uint*	dSubCellStart;
	uint*	dSubCellEnd;

	uint	uSubNumAllParticles;	//!< �\�ʐ����ɕK�v�ȃT�u�p�[�e�B�N����
	uint	uSubNumMCParticles;		//!< �\�ʐ����ɕK�v�ȃT�u�p�[�e�B�N����
	uint    uSubNumValidParticles;	//!< �����_�����O���ɗL���ȃT�u�p�[�e�B�N����
	uint	uSubNumCells;
	uint	uSubNumArdGrid;			//!< �ߖT�T�����Q�ƃO���b�h�͈�
	uint	uSubMaxLevel;
	uint	uNumParticles;		//!< ���x��0�̃p�[�e�B�N����
	uint	uMaxParticles;		//!< ���x��0�̍ő�p�[�e�B�N����
	uint	uSubHeadIndex[MAX_SUB_LEVEL+1];
	//uint	uSubNumLevel[MAX_SUB_LEVEL+1];
	uint	uSubNumEach[MAX_SUB_LEVEL+1];
	float	fSubRad[MAX_SUB_LEVEL+1];
	float	fSubEffectiveRadius[MAX_SUB_LEVEL+1];
	float	fSubEffectiveFactor;
	float	fEtcri;

	//���x���ʕ\�ʐ����p
	/*float*	dSubSelectRat[MAX_SUB_LEVEL+1];
	float4* dSubSelectPos[MAX_SUB_LEVEL+1];
	uint*	dSubSelectSortedIndex[MAX_SUB_LEVEL+1];
	uint*	dSubSelectCellStart[MAX_SUB_LEVEL+1];
	uint*	dSubSelectCellEnd[MAX_SUB_LEVEL+1];
	uint	uSubSelectNum[MAX_SUB_LEVEL+1];*/
};



//-----------------------------------------------------------------------------
// SSM
//-----------------------------------------------------------------------------
struct rxSsmParams
{
	float PrtRad;		//!< �p�[�e�B�N���̔��a
	float4 PMV[4];		//!< �������e�s��P�ƃ��f���r���[�s��MV���|�����s��
	float3 Tr;			//!< ���a�̓��e�ϊ�
	uint W, H;			//!< �`��̈�𑜓x
	float Spacing;		//!< �f�v�X�}�b�v�̃T���v�����O�Ԋu
	float Zmax;			//!< �֊s�ƂȂ�f�v�X����臒l
	int Nfilter;		//!< �f�v�X�l�������̃t�B���^�T�C�Y
	int Niters;			//!< �֊s�������̔�����
	int Ngx, Ngy;		//!< ���b�V�������p�O���b�h�̉𑜓x
};

//! �O���b�h�G�b�W
struct rxSSEdgeG
{
	float3 x0, x1;		//!< �[�_���W�ƃf�v�X�l
	float depth;		//!< �G�b�W�f�v�X�l
	int front_vertex;	//!< �G�b�W���_�̃C���f�b�N�X
	float dx;			//!< �f�v�X�l���������[�_����G�b�W���_�܂ł̋���
	int silhouette;
};


//! ���b�V�������p�O���b�h
struct rxSSGridG
{
	int i, j;
	int node_vrts[4];	//!< �m�[�h���_�C���f�b�N�X
	int num_nv;			//!< �m�[�h���_��
	int edge_vrts[4];	//!< �G�b�W���_(front vertex)
	int num_ev;			//!< �G�b�W���_��(front vertex)
	int back_vrts[6];	//!< �G�b�W���_(back vertex, back-2 vertex)
	int num_bv;			//!< �G�b�W���_��(back vertex)

	float node_depth[4];	//!< �m�[�h�̃f�v�X�l
	int vrot;

	int table_index0;	//!< �f�o�b�O�p:���b�V�����̂��߂̃C���f�b�N�X�l
	int table_index1;	//!< �f�o�b�O�p:���b�V�����̂��߂̃C���f�b�N�X�l
	int mesh_num;		//!< �f�o�b�O�p:���b�V����
	int mesh[6];		//!< �f�o�b�O�p:���b�V���C���f�b�N�X
	int back2;			//!< �f�o�b�O�p
	int v[14];
};



//-----------------------------------------------------------------------------
//! �p�b�N���邽�߂̃f�[�^�\�� 
//-----------------------------------------------------------------------------
template<class T> 
struct rxVPack
{
	T* dPos;
	T* dCompactedPos;
	uint* dOcc;
	uint* dOccScan;

	rxVPack()
	{
		dPos = 0;
		dCompactedPos = 0;
		dOcc = 0;
		dOccScan = 0;
	}
};

typedef rxVPack<float> rxVPackf;
typedef rxVPack<rxSSEdgeG> rxVPacke;

struct rxVrtAdd
{
	int num;
	int layer;
	int edge[2];
	float3 vrts[2];
};



#endif // _RX_CU_COMMON_CUH_
