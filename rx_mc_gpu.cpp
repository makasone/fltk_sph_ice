/*!
  @file rx_mc_cpu.cpp
	
  @brief �A�֐��\�ʂ���̃|���S������(MC�@)
	
	http://local.wasp.uwa.edu.au/~pbourke/geometry/polygonise/
 
  @author Raghavendra Chandrashekara (basesd on source code
			provided by Paul Bourke and Cory Gene Bloyd)
  @date   2010-03
*/


//-----------------------------------------------------------------------------
// �C���N���[�h�t�@�C��
//-----------------------------------------------------------------------------
#include <cmath>
#include "helper_math.h"

#include "rx_mc.h"

#include "rx_cu_funcs.cuh"
#include <cuda_runtime.h>


//-----------------------------------------------------------------------------
// rxMCMeshGPU�̎���
//-----------------------------------------------------------------------------
/*!
 * �R���X�g���N�^
 */
rxMCMeshGPU::rxMCMeshGPU()
{
	// MARK:�R���X�g���N�^
	m_uNumVoxels = 0;				//!< ���O���b�h��
	m_uMaxVerts = 0;				//!< �ő咸�_��
	m_uNumActiveVoxels = 0;			//!< ���b�V�������݂���{�N�Z����
	m_uNumVrts = 0;			//!< �����_��

	// �f�o�C�X������
	g_dfVolume = 0;					//!< �A�֐��f�[�^���i�[����O���b�h
	g_dfNoise = 0;					//!< �m�C�Y�l���i�[����O���b�h(�`�掞�̐F�����肷��̂Ɏg�p)
	m_duVoxelVerts = 0;				//!< �O���b�h�Ɋ܂܂�郁�b�V�����_��
	m_duVoxelVertsScan = 0;			//!< �O���b�h�Ɋ܂܂�郁�b�V�����_��(Scan)
	m_duCompactedVoxelArray = 0;	//!< ���b�V�����܂ރO���b�h���

	m_duVoxelOccupied = 0;			//!< �|���S���������ɑ��݂���{�N�Z���̃��X�g
	m_duVoxelOccupiedScan = 0;		//!< �|���S���������ɑ��݂���{�N�Z���̃��X�g(prefix scan)

#ifdef RX_CUMC_USE_GEOMETRY
	// �􉽏��𐶐�����Ƃ��ɕK�v�ȕϐ�
	m_duVoxelCubeIdx = 0;			//!< �O���b�h8���_�̉A�֐��l��臒l�ȏォ�ǂ������e�r�b�g�Ɋi�[�����ϐ�

	m_duEdgeOccupied = 0;			//!< �G�b�W�Ƀ��b�V�����_���܂ނ��ǂ���(x�����Cy����, z�����̏�)
	m_duEdgeOccupiedScan = 0;		//!< �G�b�W�Ƀ��b�V�����_���܂ނ��ǂ���(Scan)
	m_dfEdgeVrts = 0;				//!< �G�b�W���ƂɎZ�o�������_���
	m_dfCompactedEdgeVrts = 0;		//!< ���Ԃ��߂����_���
	m_duIndexArray = 0;			//!< �|���S���̊􉽏��
	m_dfVertexNormals = 0;			//!< ���_�@��
	m_duVoxelTriNum = 0;			//!< �O���b�h���Ƃ̎O�p�`���b�V����
	m_duVoxelTriNumScan = 0;		//!< �O���b�h���Ƃ̎O�p�`���b�V����(Scan)
#else
	// �􉽏���K�v�Ƃ��Ȃ��Ƃ��̂ݗp����
	m_df4Vrts = 0;					//!< �|���S�����_���W
	m_df4Nrms = 0;					//!< �|���S�����_�@��
#endif

	// �z�X�g������
	m_f4VertPos = 0;				//!< ���_���W
#ifdef RX_CUMC_USE_GEOMETRY
	m_u3TriIdx = 0;					//!< ���b�V���C���f�b�N�X
	m_uScan = 0;					//!< �f�o�b�O�p
#else
	m_f4VertNrm = 0;				//!< ���_�@��
#endif

	m_uMaxVerts = 0;
	m_iVertexStore = 4;
	m_bSet = false;

	m_ptScalarField = NULL;
	m_fpScalarFunc = 0;
	m_tIsoLevel = 0;
	m_bValidSurface = false;

	CuMCInit();
}

/*!
 * �f�X�g���N�^
 */
rxMCMeshGPU::~rxMCMeshGPU()
{
	Clean();
}



/*!
 * ���b�V�������{�N�Z���������̐ݒ�ƃf�o�C�X�������C�z�X�g�������̊m��
 * @param[out] max_verts �ő咸�_��
 * @param[in]  n �{�N�Z���������̏搔(2^n)
 */
bool rxMCMeshGPU::Set(Vec3 vMin, Vec3 vH, int n[3], unsigned int vm)
{
	float minp[3], maxp[3];
	for(int i = 0; i < 3; ++i){
		minp[i] = (float)(vMin[i]);
		maxp[i] = (float)(vMin[i]+vH[i]*n[i]);
	}

	// �O���b�h��
	m_u3GridSize = make_uint3(n[0], n[1], n[2]);

	// ���O���b�h��
	m_uNumVoxels = m_u3GridSize.x*m_u3GridSize.y*m_u3GridSize.z;

	// �O���b�h�S�̂̑傫��
	m_f3VoxelMin = make_float3(minp[0], minp[1], minp[2]);
	m_f3VoxelMax = make_float3(maxp[0], maxp[1], maxp[2]);
	float3 rng = m_f3VoxelMax-m_f3VoxelMin;
	m_f3VoxelH = make_float3(rng.x/m_u3GridSize.x, rng.y/m_u3GridSize.y, rng.z/m_u3GridSize.z);


#ifdef RX_CUMC_USE_GEOMETRY
	// �G�b�W��
	m_u3EdgeSize[0] = make_uint3(m_u3GridSize.x, m_u3GridSize.y+1, m_u3GridSize.z+1);
	m_u3EdgeSize[1] = make_uint3(m_u3GridSize.x+1, m_u3GridSize.y, m_u3GridSize.z+1);
	m_u3EdgeSize[2] = make_uint3(m_u3GridSize.x+1, m_u3GridSize.y+1, m_u3GridSize.z);
	m_uNumEdges[0] = m_u3GridSize.x*(m_u3GridSize.y+1)*(m_u3GridSize.z+1);
	m_uNumEdges[1] = (m_u3GridSize.x+1)*m_u3GridSize.y*(m_u3GridSize.z+1);
	m_uNumEdges[2] = (m_u3GridSize.x+1)*(m_u3GridSize.y+1)*m_u3GridSize.z;
	m_uNumEdges[3] = m_uNumEdges[0]+m_uNumEdges[1]+m_uNumEdges[2];
	//printf("mc edge num : %d, %d, %d - %d\n", m_uNumEdges[0], m_uNumEdges[1], m_uNumEdges[2], m_uNumEdges[3]);
#endif

	// �ő咸�_��
	m_uMaxVerts = m_u3GridSize.x*m_u3GridSize.y*vm;
	m_iVertexStore = vm;
	cout << "maximum vertex num : " << m_uMaxVerts << ", vertex store : " << m_iVertexStore << endl;

	m_uNumTris = 0;
	m_uNumVrts = 0;

	// �A�֐��l���i�[����{�����[��
	int size = m_u3GridSize.x*m_u3GridSize.y*m_u3GridSize.z*sizeof(float);
	CuAllocateArray((void**)&g_dfVolume, size);
	CuSetArrayValue((void*)g_dfVolume, 0, size);
	CuAllocateArray((void**)&g_dfNoise, size);
	CuSetArrayValue((void*)g_dfNoise, 0, size);

	// �e�[�u��
	CuInitMCTable();

	// �f�o�C�X������
	unsigned int memSize = sizeof(uint)*m_uNumVoxels;
	CuAllocateArray((void**)&m_duVoxelVerts,      memSize);
	CuAllocateArray((void**)&m_duVoxelVertsScan,  memSize);

#if SKIP_EMPTY_VOXELS
	CuAllocateArray((void**)&m_duVoxelOccupied, memSize);
	CuAllocateArray((void**)&m_duVoxelOccupiedScan, memSize);
	CuAllocateArray((void**)&m_duCompactedVoxelArray, memSize);
#endif

#ifdef RX_CUMC_USE_GEOMETRY

	CuAllocateArray((void**)&m_duVoxelCubeIdx,    memSize);

	CuAllocateArray((void**)&m_duVoxelTriNum,     memSize);
	CuAllocateArray((void**)&m_duVoxelTriNumScan, memSize);

	CuSetArrayValue((void*)m_duVoxelCubeIdx,    0, memSize);
	CuSetArrayValue((void*)m_duVoxelTriNum,     0, memSize);
	CuSetArrayValue((void*)m_duVoxelTriNumScan, 0, memSize);

	memSize = sizeof(uint)*m_uNumEdges[3];
	CuAllocateArray((void**)&m_duEdgeOccupied,	    memSize);
	CuAllocateArray((void**)&m_duEdgeOccupiedScan, memSize);

	CuSetArrayValue((void*)m_duEdgeOccupied,     0, memSize);
	CuSetArrayValue((void*)m_duEdgeOccupiedScan, 0, memSize);

	memSize = sizeof(float)*4*m_uNumEdges[3];
	CuAllocateArray((void**)&m_dfEdgeVrts,          memSize);
	CuAllocateArray((void**)&m_dfCompactedEdgeVrts, memSize);

	CuSetArrayValue((void*)m_dfEdgeVrts,          0, memSize);
	CuSetArrayValue((void*)m_dfCompactedEdgeVrts, 0, memSize);

	memSize = sizeof(float)*3*m_uNumEdges[3];
	CuAllocateArray((void**)&m_dfVertexNormals,     memSize);
	CuSetArrayValue((void*)m_dfVertexNormals,     0, memSize);
	
	memSize = sizeof(uint)*3*m_uMaxVerts*3;
	CuAllocateArray((void**)&m_duIndexArray, memSize);
	CuSetArrayValue((void*)m_duIndexArray, 0, memSize);

#else

	memSize = sizeof(float4)*m_uMaxVerts;
	CuAllocateArray((void**)&m_df4Vrts, memSize);
	CuAllocateArray((void**)&m_df4Nrms, memSize);

#endif

	// �z�X�g������
	m_f4VertPos = (float4*)malloc(sizeof(float4)*m_uMaxVerts);
	memset(m_f4VertPos, 0, sizeof(float4)*m_uMaxVerts);

#ifdef RX_CUMC_USE_GEOMETRY
	m_u3TriIdx = (uint3*)malloc(sizeof(uint3)*m_uMaxVerts*3);
	memset(m_u3TriIdx, 0, sizeof(uint3)*m_uMaxVerts*3);

	m_uScan = (uint*)malloc(sizeof(uint)*m_uNumEdges[3]);
	memset(m_uScan, 0, sizeof(uint)*m_uNumEdges[3]);
#else
	m_f4VertNrm = (float4*)malloc(sizeof(float4)*m_uMaxVerts);
	memset(m_f4VertNrm, 0, sizeof(float4)*m_uMaxVerts);
#endif

	m_bSet = true;
	return true;
}

/*!
 * �m�ۂ����z��̍폜
 */
void rxMCMeshGPU::Clean(void)
{
	if(m_bSet){
		CuFreeArray(g_dfVolume);
		CuFreeArray(g_dfNoise);

		CuFreeArray(m_duVoxelVerts);
		CuFreeArray(m_duVoxelVertsScan);

#if SKIP_EMPTY_VOXELS
		CuFreeArray(m_duVoxelOccupied);
		CuFreeArray(m_duVoxelOccupiedScan);
		CuFreeArray(m_duCompactedVoxelArray);
#endif

		CuCleanMCTable();

#ifdef RX_CUMC_USE_GEOMETRY
		CuFreeArray(m_duVoxelCubeIdx);
		CuFreeArray(m_duVoxelTriNum);
		CuFreeArray(m_duVoxelTriNumScan);
		CuFreeArray(m_duIndexArray);

		CuFreeArray(m_duEdgeOccupied);
		CuFreeArray(m_duEdgeOccupiedScan);
		CuFreeArray(m_dfEdgeVrts);
		CuFreeArray(m_dfCompactedEdgeVrts);

		CuFreeArray(m_dfVertexNormals);

		if(m_u3TriIdx != 0) free(m_u3TriIdx);
		if(m_uScan != 0) free(m_uScan);
#else
		CuFreeArray(m_df4Vrts);
		CuFreeArray(m_df4Nrms);

		if(m_f4VertNrm != 0) free(m_f4VertNrm);
#endif

		if(m_f4VertPos != 0) free(m_f4VertPos);

		m_bSet = false;
	}
}

/*!
 * �A�֐�����O�p�`���b�V���𐶐�
 * @param[in] func �A�֐��l�擾�p�֐��|�C���^
 * @param[in] min_p �O���b�h�̍ŏ����W
 * @param[in] h �O���b�h�̕�
 * @param[in] n[3] �O���b�h��(x,y,z)
 * @param[in] threshold �������l(�A�֐��l�����̒l�̂Ƃ�������b�V����)
 * @param[in] method ���b�V���������@("mc", "rmt", "bloomenthal")
 * @param[out] vrts ���_���W
 * @param[out] nrms ���_�@��
 * @param[out] tris ���b�V��
 * @retval true  ���b�V����������
 * @retval false ���b�V���������s
 */
bool rxMCMeshGPU::CreateMesh(float (*func)(double, double, double), Vec3 minp, double h, int n[3], float threshold, string method, 
									vector<Vec3> &vrts, vector<Vec3> &nrms, vector<rxFace> &face)
{
	if(func == NULL) return false;

	// �{�����[���f�[�^�쐬
	float *field = new float[n[0]*n[1]*n[2]];
	for(int k = 0; k < n[2]; ++k){
		for(int j = 0; j < n[1]; ++j){
			for(int i = 0; i < n[0]; ++i){
				int idx = k*n[0]*n[1]+j*n[0]+i;
				Vec3 pos = minp+Vec3(i, j, k)*h;

				float val = func(pos[0], pos[1], pos[2]);
				field[idx] = val;
			}
		}
	}

	// �{�����[���f�[�^���f�o�C�X�������ɃR�s�[
	int size = m_u3GridSize.x*m_u3GridSize.y*m_u3GridSize.z*sizeof(float);
	CuCopyArrayToDevice((void*)g_dfVolume, (void*)field, 0, size);

	// ���b�V������
	uint nvrts, ntris;
	CreateMeshV(minp, h, n, threshold, nvrts, ntris);

	// �f�[�^��z��Ɋi�[
	SetDataToArray(vrts, nrms, face);

	return true;
}


/*!
 * �f�o�C�X�{�����[���f�[�^����O�p�`���b�V���𐶐�
 * @param[in] field �T���v���{�����[��
 * @param[in] min_p �O���b�h�̍ŏ����W
 * @param[in] h �O���b�h�̕�
 * @param[in] n[3] �O���b�h��(x,y,z)
 * @param[in] threshold �������l(�A�֐��l�����̒l�̂Ƃ�������b�V����)
 * @param[in] method ���b�V���������@("mc", "rmt", "bloomenthal")
 * @param[out] vrts ���_���W
 * @param[out] nrms ���_�@��
 * @param[out] tris ���b�V��
 * @retval true  ���b�V����������
 * @retval false ���b�V���������s
 */
bool rxMCMeshGPU::CreateMeshV(Vec3 minp, double h, int n[3], float threshold, uint &nvrts, uint &ntris)
{
	m_uNumVrts = 0;

	// �{�N�Z�����̒��_���C���b�V�����̎Z�o
	CuMCCalTriNum(g_dfVolume, m_duVoxelCubeIdx, m_duVoxelVerts, m_duVoxelVertsScan, 
				  m_duVoxelTriNum, m_duVoxelTriNumScan, m_duVoxelOccupied, m_duVoxelOccupiedScan, m_duCompactedVoxelArray, 
				  m_u3GridSize, m_uNumVoxels, m_f3VoxelH, threshold, 
				  m_uNumActiveVoxels, m_uNumVrts, m_uNumTris);

	// �G�b�W���_�̎Z�o
	CuMCCalEdgeVrts(g_dfVolume, m_dfEdgeVrts, m_dfCompactedEdgeVrts, 
					m_duEdgeOccupied, m_duEdgeOccupiedScan, m_u3EdgeSize, m_uNumEdges, 
					m_u3GridSize, m_uNumVoxels, m_f3VoxelH, m_f3VoxelMin, threshold, 
					m_uNumVrts);

	if(m_uNumVrts){
		// ���b�V������
		CuMCCalTri(m_duIndexArray, m_duVoxelCubeIdx, m_duVoxelTriNumScan, m_duCompactedVoxelArray, 
				   m_duEdgeOccupiedScan, m_u3EdgeSize, m_uNumEdges, 
				   m_u3GridSize, m_uNumVoxels, m_f3VoxelH, threshold, 
				   m_uNumActiveVoxels, m_uNumVrts, m_uNumTris);

		// ���_�@���v�Z
		CuMCCalNrm(m_dfVertexNormals, m_duIndexArray, m_dfCompactedEdgeVrts, m_uNumVrts, m_uNumTris);
	}
	else{
		m_uNumTris = 0;
	}


	nvrts = m_uNumVrts;
	ntris = m_uNumTris;

	return true;
}

/*!
 * FBO�Ƀf�[�^��ݒ�
 * @param[in] uVrtVBO ���_FBO
 * @param[in] uNrmVBO �@��FBO
 * @param[in] uTriVBO ���b�V��FBO
 */
bool rxMCMeshGPU::SetDataToFBO(GLuint uVrtVBO, GLuint uNrmVBO, GLuint uTriVBO)
{
	int memsize = 0;

	// ���_���̎擾
	memsize = sizeof(float)*m_uNumVrts*4;

	glBindBuffer(GL_ARRAY_BUFFER, uVrtVBO);
	float *vrt_ptr = (float*)glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);

	CuCopyArrayFromDevice(vrt_ptr, GetVrtDev(), 0, memsize);

	glUnmapBuffer(GL_ARRAY_BUFFER);
	glBindBuffer(GL_ARRAY_BUFFER, 0);

	// �@�����̎擾
	memsize = sizeof(float)*m_uNumVrts*4;

	glBindBuffer(GL_ARRAY_BUFFER, uNrmVBO);
	float *nrm_ptr = (float*)glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);

	CuCopyArrayFromDevice(nrm_ptr, GetNrmDev(), 0, memsize);

	glUnmapBuffer(GL_ARRAY_BUFFER);
	glBindBuffer(GL_ARRAY_BUFFER, 0);

	// �ڑ����̎擾
	memsize = sizeof(uint)*m_uNumTris*3;

	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, uTriVBO);
	uint *tri_ptr = (uint*)glMapBuffer(GL_ELEMENT_ARRAY_BUFFER, GL_WRITE_ONLY);

	CuCopyArrayFromDevice(tri_ptr, GetIdxDev(), 0, memsize);

	glUnmapBuffer(GL_ELEMENT_ARRAY_BUFFER);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

	return true;
}

/*!
 * �z�X�g���z��Ƀf�[�^��ݒ�
 * @param[out] vrts ���_���W
 * @param[out] nrms ���_�@��
 * @param[out] tris ���b�V��
 */
bool rxMCMeshGPU::SetDataToArray(vector<Vec3> &vrts, vector<Vec3> &nrms, vector<rxFace> &face)
{
	// ���_���̎擾
	float *vrt_ptr = new float[m_uNumVrts*4];
	CuCopyArrayFromDevice(vrt_ptr, GetVrtDev(), 0, m_uNumVrts*4*sizeof(float));

	vrts.resize(m_uNumVrts);
	for(uint i = 0; i < m_uNumVrts; ++i){
		vrts[i][0] = vrt_ptr[4*i];
		vrts[i][1] = vrt_ptr[4*i+1];
		vrts[i][2] = vrt_ptr[4*i+2];
	}
	
	// �@�����̎擾
	CuCopyArrayFromDevice(vrt_ptr, GetNrmDev(), 0, m_uNumVrts*4*sizeof(float));

	nrms.resize(m_uNumVrts);
	for(uint i = 0; i < m_uNumVrts; ++i){
		nrms[i][0] = vrt_ptr[4*i];
		nrms[i][1] = vrt_ptr[4*i+1];
		nrms[i][2] = vrt_ptr[4*i+2];
	}

	delete [] vrt_ptr;

	// �ڑ����̎擾
	uint *tri_ptr = new uint[m_uNumTris*3];
	CuCopyArrayFromDevice(tri_ptr, GetIdxDev(), 0, m_uNumTris*3*sizeof(uint));

	face.resize(m_uNumTris);
	for(uint i = 0; i < m_uNumTris; ++i){
		face[i].vert_idx.resize(3);
		face[i][0] = tri_ptr[3*i];
		face[i][1] = tri_ptr[3*i+1];
		face[i][2] = tri_ptr[3*i+2];
	}

	delete [] tri_ptr;

	return true;
}

/*!
 * �T���v���{�����[����ݒ�
 * @param[in] hVolume �T���v���{�����[��
 */
void rxMCMeshGPU::SetSampleVolumeFromHost(float *hVolume)
{
	int size = m_u3GridSize.x*m_u3GridSize.y*m_u3GridSize.z*sizeof(float);
	CuCopyArrayToDevice((void*)g_dfVolume, (void*)hVolume, 0, size);
}
float* rxMCMeshGPU::GetSampleVolumeDevice(void)
{
	return g_dfVolume;
}

/*!
 * �m�C�Y�t���T���v���{�����[����ݒ�
 * @param[in] 
 * @return 
 */
void rxMCMeshGPU::SetSampleNoiseFromHost(float *hVolume)
{
	int size = m_u3GridSize.x*m_u3GridSize.y*m_u3GridSize.z*sizeof(float);
	CuCopyArrayToDevice((void*)g_dfNoise, (void*)hVolume, 0, size);
}
float* rxMCMeshGPU::GetSampleNoiseDevice(void)
{
	return g_dfNoise;
}


