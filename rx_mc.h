/*!
  @file rx_mc.h
	
  @brief �A�֐��\�ʂ���̃|���S������(MC�@)
	
	http://local.wasp.uwa.edu.au/~pbourke/geometry/polygonise/
 
  @author Raghavendra Chandrashekara (basesd on source code
			provided by Paul Bourke and Cory Gene Bloyd)
  @date   2010-03
*/


#ifndef _RX_MC_MESH_H_
#define _RX_MC_MESH_H_


//-----------------------------------------------------------------------------
// �C���N���[�h�t�@�C��
//-----------------------------------------------------------------------------
// C�W��
#include <cstdlib>

// OpenGL
#include <GL/glew.h>
#include <GL/glut.h>

// STL
#include <map>
#include <vector>
#include <string>

#include <iostream>

#include "rx_utility.h"
#include "rx_mesh.h"

#include "rx_cu_common.cuh"


//-----------------------------------------------------------------------------
// ���O���
//-----------------------------------------------------------------------------
using namespace std;


//-----------------------------------------------------------------------------
// ��`
//-----------------------------------------------------------------------------
typedef unsigned int uint;

#ifndef RXREAL
	#define RXREAL float
#endif


struct RxVertexID
{
	uint newID;
	double x, y, z;
};

typedef std::map<uint, RxVertexID> ID2VertexID;

struct RxTriangle
{
	uint pointID[3];
};

typedef std::vector<RxTriangle> RxTriangleVector;

struct RxScalarField
{
	uint iNum[3];
	Vec3 fWidth;
	Vec3 fMin;
};



//-----------------------------------------------------------------------------
// rxMCMeshCPU�N���X
//-----------------------------------------------------------------------------
class rxMCMeshCPU
{
public:
	// �R���X�g���N�^
	rxMCMeshCPU();

	// �f�X�g���N�^
	~rxMCMeshCPU();
	
	//! �A�֐�����O�p�`���b�V���𐶐�
	bool CreateMesh(RXREAL (*func)(double, double, double), Vec3 min_p, double h, int n[3], RXREAL threshold, 
					vector<Vec3> &vrts, vector<Vec3> &nrms, vector<rxFace> &face);

	//! �T���v���{�����[������O�p�`���b�V���𐶐�
	bool CreateMeshV(RXREAL *field, Vec3 min_p, double h, int n[3], RXREAL threshold, 
					 vector<Vec3> &vrts, vector<Vec3> &nrms, vector<rxFace> &face);

	//! �T���v���{�����[�����瓙�l�ʃ��b�V������
	void GenerateSurface(const RxScalarField sf, RXREAL *field, RXREAL threshold, 
						 vector<Vec3> &vrts, vector<Vec3> &nrms, vector<int> &tris);

	//! �֐�����T���v���{�����[�����쐬���ē��l�ʃ��b�V������
	void GenerateSurfaceV(const RxScalarField sf, RXREAL (*func)(double, double, double), RXREAL threshold, 
						  vector<Vec3> &vrts, vector<Vec3> &nrms, vector<int> &tris);

	//! �֐����瓙�l�ʃ��b�V������
	void GenerateSurfaceF(const RxScalarField sf, RXREAL (*func)(double, double, double), RXREAL threshold, 
						  vector<Vec3> &vrts, vector<Vec3> &nrms, vector<int> &tris);

	//! ���l�ʍ쐬������������true��Ԃ�
	bool IsSurfaceValid() const { return m_bValidSurface; }

	//! �쐬�����l�ʃ��b�V���̔j��
	void DeleteSurface();

	//! ���b�V�����ɗp�����O���b�h�̑傫��(���b�V���쐬���Ă��Ȃ��ꍇ�͕Ԓl��-1)
	int GetVolumeLengths(double& fVolLengthX, double& fVolLengthY, double& fVolLengthZ);

	// �쐬�������b�V���̏��
	uint GetNumVertices(void) const { return m_nVertices; }
	uint GetNumTriangles(void) const { return m_nTriangles; }
	uint GetNumNormals(void) const { return m_nNormals; }

protected:
	// MARK:�����o�ϐ�
	uint m_nVertices;	//!< ���l�ʃ��b�V���̒��_��
	uint m_nNormals;	//!< ���l�ʃ��b�V���̒��_�@����(�쐬����Ă���� �@����=���_��)
	uint m_nTriangles;	//!< ���l�ʃ��b�V���̎O�p�`�|���S����

	ID2VertexID m_i2pt3idVertices;			//!< ���l�ʂ��`�����钸�_�̃��X�g
	RxTriangleVector m_trivecTriangles;		//!< �O�p�`�|���S�����`�����钸�_�̃��X�g

	RxScalarField m_Grid;					//!< �����O���b�h���

	// �A�֐��l(�X�J���[�l)�擾�p�ϐ�(�ǂ��炩�̂ݗp����)
	const RXREAL* m_ptScalarField;				//!< �X�J���[�l��ێ�����T���v���{�����[��
	RXREAL (*m_fpScalarFunc)(double, double, double);	//!< �X�J���[�l��Ԃ��֐��|�C���^

	RXREAL m_tIsoLevel;							//!< 臒l

	bool m_bValidSurface;					//!< ���b�V�����������̉�


	// ���b�V���\�z�p�̃e�[�u��
	static const uint m_edgeTable[256];
	static const int m_triTable[256][16];



	// MARK:protected�����o�֐�

	//! �G�b�WID
	uint GetEdgeID(uint nX, uint nY, uint nZ, uint nEdgeNo);

	//! ���_ID
	uint GetVertexID(uint nX, uint nY, uint nZ);

	// �G�b�W��̓��l�_���v�Z
	RxVertexID CalculateIntersection(uint nX, uint nY, uint nZ, uint nEdgeNo);
	RxVertexID CalculateIntersectionF(uint nX, uint nY, uint nZ, uint nEdgeNo);

	//! �O���b�h�G�b�W���[�̉A�֐��l������^��Ԃœ��l�_���v�Z
	RxVertexID Interpolate(double fX1, double fY1, double fZ1, double fX2, double fY2, double fZ2, RXREAL tVal1, RXREAL tVal2);

	//! ���_�C���b�V���􉽏����o�͌`���Ŋi�[
	void RenameVerticesAndTriangles(vector<Vec3> &vrts, uint &nvrts, vector<int> &tris, uint &ntris);

	//! ���_�@���v�Z
	void CalculateNormals(const vector<Vec3> &vrts, uint nvrts, const vector<int> &tris, uint ntris, 
						  vector<Vec3> &nrms, uint &nnrms);

};


//-----------------------------------------------------------------------------
// rxMCMeshGPU�N���X
//-----------------------------------------------------------------------------
class rxMCMeshGPU
{
protected:
	// MC�@�p
	uint3 m_u3GridSize;				//!< �O���b�h��(nx,ny,nz)
	uint3 m_u3GridSizeMask;			//!< �O���b�h/�C���f�b�N�X�ϊ����̃}�X�N
	uint3 m_u3GridSizeShift;		//!< �O���b�h/�C���f�b�N�X�ϊ����̃V�t�g��

	float3 m_f3VoxelMin;			//!< �O���b�h�ŏ��ʒu
	float3 m_f3VoxelMax;			//!< �O���b�h�ő�ʒu
	float3 m_f3VoxelH;				//!< �O���b�h��
	uint m_uNumVoxels;				//!< ���O���b�h��
	uint m_uMaxVerts;				//!< �ő咸�_��
	uint m_uNumActiveVoxels;		//!< ���b�V�������݂���{�N�Z����
	uint m_uNumVrts;				//!< �����_��
	uint m_uNumTris;				//!< �����b�V����

	// �f�o�C�X������
	float *g_dfVolume;				//!< �A�֐��f�[�^���i�[����O���b�h
	float *g_dfNoise;				//!< �m�C�Y�l���i�[����O���b�h(�`�掞�̐F�����肷��̂Ɏg�p)
	uint *m_duVoxelVerts;			//!< �O���b�h�Ɋ܂܂�郁�b�V�����_��
	uint *m_duVoxelVertsScan;		//!< �O���b�h�Ɋ܂܂�郁�b�V�����_��(Scan)
	uint *m_duCompactedVoxelArray;	//!< ���b�V�����܂ރO���b�h���

	uint *m_duVoxelOccupied;		//!< �|���S���������ɑ��݂���{�N�Z���̃��X�g
	uint *m_duVoxelOccupiedScan;	//!< �|���S���������ɑ��݂���{�N�Z���̃��X�g(prefix scan)

#ifdef RX_CUMC_USE_GEOMETRY
	// �􉽏��𐶐�����Ƃ��ɕK�v�ȕϐ�
	uint3 m_u3EdgeSize[3];			//!< �G�b�W��(nx,ny,nz)
	uint m_uNumEdges[4];			//!< ���G�b�W��
#endif

#ifdef RX_CUMC_USE_GEOMETRY
	// �􉽏��𐶐�����Ƃ��ɕK�v�ȕϐ�
	uint *m_duVoxelCubeIdx;			//!< �O���b�h8���_�̉A�֐��l��臒l�ȏォ�ǂ������e�r�b�g�Ɋi�[�����ϐ�

	uint *m_duEdgeOccupied;			//!< �G�b�W�Ƀ��b�V�����_���܂ނ��ǂ���(x�����Cy����, z�����̏�)
	uint *m_duEdgeOccupiedScan;		//!< �G�b�W�Ƀ��b�V�����_���܂ނ��ǂ���(Scan)
	float *m_dfEdgeVrts;			//!< �G�b�W���ƂɎZ�o�������_���
	float *m_dfCompactedEdgeVrts;	//!< ���Ԃ��߂����_���
	uint *m_duIndexArray;			//!< �|���S���̊􉽏��
	float *m_dfVertexNormals;		//!< ���_�@��
	uint *m_duVoxelTriNum;			//!< �O���b�h���Ƃ̎O�p�`���b�V����
	uint *m_duVoxelTriNumScan;		//!< �O���b�h���Ƃ̎O�p�`���b�V����(Scan)
#else
	// �􉽏���K�v�Ƃ��Ȃ��Ƃ��̂ݗp����
	float4 *m_df4Vrts;				//!< �|���S�����_���W
	float4 *m_df4Nrms;				//!< �|���S�����_�@��
#endif

	// �z�X�g������
	float4 *m_f4VertPos;			//!< ���_���W
#ifdef RX_CUMC_USE_GEOMETRY
	uint3 *m_u3TriIdx;				//!< ���b�V���C���f�b�N�X
	uint *m_uScan;					//!< �f�o�b�O�p
#else
	float4 *m_f4VertNrm;			//!< ���_�@��
#endif

	int m_iVertexStore;
	bool m_bSet;


	// �A�֐��l(�X�J���[�l)�擾�p�ϐ�(�ǂ��炩�̂ݗp����)
	const float* m_ptScalarField;				//!< �X�J���[�l��ێ�����T���v���{�����[��
	float (*m_fpScalarFunc)(double, double, double);	//!< �X�J���[�l��Ԃ��֐��|�C���^

	float m_tIsoLevel;							//!< 臒l

	bool m_bValidSurface;					//!< ���b�V�����������̉�

public:
	// �R���X�g���N�^
	rxMCMeshGPU();

	// �f�X�g���N�^
	~rxMCMeshGPU();
	
	//! �A�֐�����O�p�`���b�V���𐶐�
	bool CreateMesh(float (*func)(double, double, double), Vec3 min_p, double h, int n[3], float threshold, string method, 
					vector<Vec3> &vrts, vector<Vec3> &nrms, vector<rxFace> &face);

	//! �T���v���{�����[������O�p�`���b�V���𐶐�
	bool CreateMeshV(Vec3 minp, double h, int n[3], float threshold, uint &nvrts, uint &ntris);
	
	//! �z��̊m��
	bool Set(Vec3 vMin, Vec3 vH, int n[3], unsigned int vm);

	//! �z��̍폜
	void Clean(void);

	//! FBO�Ƀf�[�^��ݒ�
	bool SetDataToFBO(GLuint uVrtVBO, GLuint uNrmVBO, GLuint uTriVBO);

	//! �z�X�g���z��Ƀf�[�^��ݒ�
	bool SetDataToArray(vector<Vec3> &vrts, vector<Vec3> &nrms, vector<rxFace> &face);

	//! �T���v���{�����[���Z�b�g
	void   SetSampleVolumeFromHost(float *hVolume);
	float* GetSampleVolumeDevice(void);

	void   SetSampleNoiseFromHost(float *hVolume);
	float* GetSampleNoiseDevice(void);

	// �ő咸�_��
	uint GetMaxVrts(void){ return m_uMaxVerts; }

	// �ő咸�_���v�Z�p�W��
	void SetVertexStore(int vs){ m_iVertexStore = vs; }

	
#ifdef RX_CUMC_USE_GEOMETRY
	//! ���_�f�[�^(�f�o�C�X������)
	float* GetVrtDev(void){ return (float*)m_dfCompactedEdgeVrts; }

	//! ���b�V���f�[�^(�f�o�C�X������)
	uint* GetIdxDev(void){ return (uint*)m_duIndexArray; }

	//! �@���f�[�^(�f�o�C�X������)
	float* GetNrmDev(void){ return (float*)m_dfVertexNormals; }

	//! ���_�f�[�^(�z�X�g������)
	float GetVertex(int idx, int dim)
	{
		if(dim == 0){
			return m_f4VertPos[idx].x;
		}
		else if(dim == 1){
			return m_f4VertPos[idx].y;
		}
		else{
			return m_f4VertPos[idx].z;
		}
	}
	void GetVertex2(int idx, float *x, float *y, float *z)
	{
		*x = m_f4VertPos[idx].x;
		*y = m_f4VertPos[idx].y;
		*z = m_f4VertPos[idx].z;
	}

	//! ���b�V���f�[�^(�z�X�g������)
	void GetTriIdx(int idx, unsigned int *vidx0, unsigned int *vidx1, unsigned int *vidx2)
	{
		*vidx0 = m_u3TriIdx[idx].x;
		*vidx1 = m_u3TriIdx[idx].y;
		*vidx2 = m_u3TriIdx[idx].z;
	}
#else // #ifdef RX_CUMC_USE_GEOMETRY

	//! ���_�f�[�^(�f�o�C�X������)
	float* GetVrtDev(void)
	{
		return (float*)m_df4Vrts;
	}

	//! �@���f�[�^(�f�o�C�X������)
	float* GetNrmDev(void)
	{
		return (float*)m_df4Nrms;
	}

#endif // #ifdef RX_CUMC_USE_GEOMETRY
};





/*!
 * ���b�V������(�T���v���{�����[���쐬)
 * @param[in] sf �����O���b�h���
 * @param[in] func �A�֐��l�擾�֐��|�C���^
 * @param[in] threshold 臒l
 * @param[out] vrts ���b�V�����_
 * @param[out] nrms ���b�V�����_�@��
 * @param[out] tris ���b�V���􉽏��(���_�ڑ����)
 */
static void GenerateValueArray(RXREAL **field, RXREAL (*func)(void*, double, double, double), void* func_ptr, const RxScalarField sf)
{
	int nx, ny, nz;
	nx = sf.iNum[0]+1;
	ny = sf.iNum[1]+1;
	nz = sf.iNum[2]+1;

	Vec3 minp = sf.fMin;
	Vec3 d = sf.fWidth;

	if(*field) delete [] *field;
	*field = new RXREAL[nx*ny*nz];
	for(int k = 0; k < nz; ++k){
		for(int j = 0; j < ny; ++j){
			for(int i = 0; i < nx; ++i){
				int idx = k*nx*ny+j*nx+i;
				Vec3 pos = minp+Vec3(i, j, k)*d;

				RXREAL val = func(func_ptr, pos[0], pos[1], pos[2]);
				(*field)[idx] = val;
			}
		}
	}
}



#endif // _RX_MC_MESH_H_

