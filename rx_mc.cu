/*! 
  @file rx_mc.cu
	
  @brief CUDA�ɂ�郁�b�V������(MC�@)

  �{�����[���f�[�^����MC(Marching Cube)�@��p����臒l�Ɋ�Â��\�ʂ𒊏o����D
  ����ɏ������邽�߂�Scan(Prefix Sum)�𗘗p���C
  Scan�̍��������ɂ�CUDPP��p����D
  
  MC�@�Ɋւ�����
  http://local.wasp.uwa.edu.au/~pbourke/geometry/polygonise/
  http://en.wikipedia.org/wiki/Marching_cubes
  
  MC�@�̃e�[�u��
  http://local.wasp.uwa.edu.au/~pbourke/geometry/polygonise/
  
  CUDPP(CUDA Data Parallel Primitives Library)�Ɋւ�����
  http://www.gpgpu.org/developer/cudpp
 
  @author Makoto Fujisawa
  @date 2009-08, 2011-06
*/
// FILE --rx_mc.cu--



//-----------------------------------------------------------------------------
// �C���N���[�h�t�@�C��
//-----------------------------------------------------------------------------
#include <cstdio>
#include <GL/glew.h>
#include <GL/glut.h>


#include "rx_cuda_utils.h"

#include "rx_mc_kernel.cu"
#include "rx_mc_tables.h"

#include "rx_cu_funcs.cuh"



//-----------------------------------------------------------------------------
// MARK:�O���[�o���ϐ�

// MC�@�̃e�[�u��
uint* g_puNumVrtsTable = 0;		//!< �{�N�Z�����̒��_���e�[�u��
uint* g_puEdgeTable = 0;		//!< �{�N�Z�����̃G�b�W�e�[�u��
uint* g_puTriTable = 0;			//!< �{�N�Z�����̃��b�V���\���e�[�u��


//-----------------------------------------------------------------------------
// CUDA�֐�
//-----------------------------------------------------------------------------
extern "C"
{



/*!
 * CUDA������
 */
bool CuMCInit(void)
{
	//int argc = 0;
	//char** argv = NULL;

	//if(cutCheckCmdLineFlag(argc, (const char**)argv, "device")){
	//	// �R�}���h���C���Ŏw�肳�ꂽCUDA�f�o�C�X���g�p
	//	cutilDeviceInit(argc, argv);
	//}
	//else{
		// �R�}���h���C���Ŏw�肳��Ă��Ȃ��ꍇ��Gflops/s���ł������f�o�C�X���g�p����
		cudaSetDevice(gpuGetMaxGflopsDeviceId());
	//}

	return true;
}

void DumpBuffer(uint *d_buffer, int nelements)
{
	uint bytes = nelements*sizeof(uint);
	uint *h_buffer = (uint *) malloc(bytes);
	RX_CUCHECK(cudaMemcpy(h_buffer, d_buffer, bytes, cudaMemcpyDeviceToHost) );
	for(int i=0; i<nelements; i++) {
		printf("%d: %u\n", i, h_buffer[i]);
	}
	printf("\n");
	free(h_buffer);
}

/*!
 * MC�@�̃e�[�u��
 */
void CuInitMCTable(void)
{
	cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc(32, 0, 0, 0, cudaChannelFormatKindUnsigned);

	RX_CUCHECK(cudaMalloc((void**) &g_puEdgeTable, 256*sizeof(uint)));
	RX_CUCHECK(cudaMemcpy(g_puEdgeTable, edgeTable, 256*sizeof(uint), cudaMemcpyHostToDevice) );
	RX_CUCHECK(cudaBindTexture(0, edgeTex, g_puEdgeTable, channelDesc) );

	RX_CUCHECK(cudaMalloc((void**) &g_puTriTable, 256*16*sizeof(uint)));
	RX_CUCHECK(cudaMemcpy(g_puTriTable, triTable, 256*16*sizeof(uint), cudaMemcpyHostToDevice) );
	RX_CUCHECK(cudaBindTexture(0, triTex, g_puTriTable, channelDesc) );

	RX_CUCHECK(cudaMalloc((void**) &g_puNumVrtsTable, 256*sizeof(uint)));
	RX_CUCHECK(cudaMemcpy(g_puNumVrtsTable, numVertsTable, 256*sizeof(uint), cudaMemcpyHostToDevice) );
	RX_CUCHECK(cudaBindTexture(0, numVertsTex, g_puNumVrtsTable, channelDesc) );
}

/*!
 * MC�@�̃e�[�u���̔j��
 */
void CuCleanMCTable(void)
{
	RX_CUCHECK(cudaFree(g_puEdgeTable));
	RX_CUCHECK(cudaFree(g_puTriTable));
	RX_CUCHECK(cudaFree(g_puNumVrtsTable));
}



#define DEBUG_BUFFERS 0


#ifdef RX_CUMC_USE_GEOMETRY

/*!
 * �e�{�N�Z����8���_�̓��O�𔻒肵�C�e�[�u�����璸�_���C�|���S���������߂�
 * @param[in] dVolume �A�֐��f�[�^���i�[����O���b�h
 * @param[out] dVoxBit �O���b�h8���_�̉A�֐��l��臒l�ȏォ�ǂ������e�r�b�g�Ɋi�[�����ϐ�
 * @param[out] dVoxVNum �O���b�h�Ɋ܂܂�郁�b�V�����_��
 * @param[out] dVoxVNumScan �O���b�h�Ɋ܂܂�郁�b�V�����_����Scan
 * @param[out] dVoxTNum �O���b�h���Ƃ̎O�p�`���b�V����
 * @param[out] dVoxTNumScan �O���b�h���Ƃ̎O�p�`���b�V������Scan
 * @param[out] dVoxOcc �O���b�h���Ƃ̃��b�V�����ݏ��(���b�V���������1,�����łȂ����0)
 * @param[out] dVoxOccScan �O���b�h���Ƃ̃��b�V�����ݏ���Scan
 * @param[out] dCompactedVox ���b�V�����܂ރO���b�h���
 * @param[in] grid_size  �O���b�h��(nx,ny,nz)
 * @param[in] num_voxels ���O���b�h��
 * @param[in] grid_width �O���b�h��
 * @param[in] threshold  ���b�V����臒l
 * @param[out] num_active_voxels ���b�V�����܂ރ{�N�Z���̐�
 * @param[out] nvrts ���_��
 * @param[out] ntris ���b�V����
*/
void CuMCCalTriNum(float *dVolume, uint *dVoxBit, uint *dVoxVNum, uint *dVoxVNumScan, 
				   uint *dVoxTNum, uint *dVoxTNumScan, uint *dVoxOcc, uint *dVoxOccScan, uint *dCompactedVox, 
				   uint3 grid_size, uint num_voxels, float3 grid_width, float threshold, 
				   uint &num_active_voxels, uint &nvrts, uint &ntris)
{
	uint lval, lsval;

	// 1�X���b�h/�{�N�Z��
	int threads = THREAD_NUM;
	dim3 grid((num_voxels+threads-1)/threads, 1, 1);
	if(grid.x > 65535){
		grid.y = (grid.x+32768-1)/32768;
		grid.x = 32768;
	}

	// �e�{�N�Z����8���_�̓��O�𔻒肵�C�e�[�u�����璸�_���C�|���S���������߂�
	ClassifyVoxel2<<<grid, threads>>>(dVoxBit, dVoxVNum, dVoxTNum, dVoxOcc, 
									  dVolume, grid_size, num_voxels, grid_width, threshold);
	RX_CUERROR("ClassifyVoxel2 failed");

	num_active_voxels = num_voxels;

#if SKIP_EMPTY_VOXELS
	// ��̃O���b�h���X�L�b�v����ꍇ

	// �O���b�h���̃��b�V���L���z���Scan
	CuScan(dVoxOccScan, dVoxOcc, num_voxels);

	// Exclusive scan (�Ō�̗v�f��0�Ԗڂ���n-2�Ԗڂ܂ł̍��v�ɂȂ��Ă���)�Ȃ̂ŁC
	// Scan�O�z��̍Ō�(n-1�Ԗ�)�̗v�f�ƍ��v���邱�ƂŃO���b�h�����v�Z
	RX_CUCHECK(cudaMemcpy((void*)&lval, (void*)(dVoxOcc+num_voxels-1), sizeof(uint), cudaMemcpyDeviceToHost));
	RX_CUCHECK(cudaMemcpy((void*)&lsval, (void*)(dVoxOccScan+num_voxels-1), sizeof(uint), cudaMemcpyDeviceToHost));
	num_active_voxels = lval+lsval;

	if(!num_active_voxels){
		nvrts = 0; ntris = 0;
		return;
	}

	CompactVoxels<<<grid, threads>>>(dCompactedVox, dVoxOcc, dVoxOccScan, num_voxels);
	RX_CUERROR("CompactVoxels failed");

#endif // #if SKIP_EMPTY_VOXELS

	// �o�[�e�b�N�X�J�E���g�pScan(prefix sum)�쐬
	CuScan(dVoxVNumScan, dVoxVNum, num_voxels);

	// �O�p�`���b�V��������Scan(prefix sum)
	CuScan(dVoxTNumScan, dVoxTNum, num_voxels);

	// Exclusive scan (�Ō�̗v�f��0�Ԗڂ���n-2�Ԗڂ܂ł̍��v�ɂȂ��Ă���)�Ȃ̂ŁC
	// Scan�O�z��̍Ō�(n-1�Ԗ�)�̗v�f�ƍ��v���邱�ƂŃ|���S�������v�Z
	RX_CUCHECK(cudaMemcpy((void*)&lval, (void*)(dVoxVNum+num_voxels-1), sizeof(uint), cudaMemcpyDeviceToHost));
	RX_CUCHECK(cudaMemcpy((void*)&lsval, (void*)(dVoxVNumScan+num_voxels-1), sizeof(uint), cudaMemcpyDeviceToHost));
	nvrts = lval+lsval;

	// Exclusive scan (�Ō�̗v�f��0�Ԗڂ���n-2�Ԗڂ܂ł̍��v�ɂȂ��Ă���)�Ȃ̂ŁC
	// Scan�O�z��̍Ō�(n-1�Ԗ�)�̗v�f�ƍ��v���邱�ƂŃ|���S�������v�Z
	RX_CUCHECK(cudaMemcpy((void*)&lval, (void*)(dVoxTNum+num_voxels-1), sizeof(uint), cudaMemcpyDeviceToHost));
	RX_CUCHECK(cudaMemcpy((void*)&lsval, (void*)(dVoxTNumScan+num_voxels-1), sizeof(uint), cudaMemcpyDeviceToHost));
	ntris = lval+lsval;
}
	

/*!
 * �G�b�W���Ƃɒ��_���W���v�Z
 * @param[in] dVolume �A�֐��f�[�^���i�[����O���b�h
 * @param[in] dNoise �m�C�Y�f�[�^���i�[����O���b�h
 * @param[out] dEdgeVrts �G�b�W���ƂɎZ�o�������_���
 * @param[out] dCompactedEdgeVrts ���Ԃ��߂����_���
 * @param[out] dEdgeOcc �G�b�W�Ƀ��b�V�����_���܂ނ��ǂ���(x�����Cy����, z�����̏�)
 * @param[out] dEdgeOccScan �G�b�W�Ƀ��b�V�����_���܂ނ��ǂ�����Scan
 * @param[in] edge_size[3] �G�b�W��(nx,ny,nz)
 * @param[in] num_edge[4] ���G�b�W��
 * @param[in] grid_size  �O���b�h��(nx,ny,nz)
 * @param[in] num_voxels ���O���b�h��
 * @param[in] grid_width �O���b�h��
 * @param[in] threshold  ���b�V����臒l
 * @param[inout] nvrts ���_��
 */
void CuMCCalEdgeVrts(float *dVolume, float *dEdgeVrts, float *dCompactedEdgeVrts, 
					 uint *dEdgeOcc, uint *dEdgeOccScan, uint3 edge_size[3], uint num_edge[4], 
					 uint3 grid_size, uint num_voxels, float3 grid_width, float3 grid_min, float threshold, 
					 uint &nvrts)
{
	uint lval, lsval;
	//
	// �G�b�W���Ƃɒ��_���W���v�Z
	//
	uint3 dir[3];
	dir[0] = make_uint3(1, 0, 0);
	dir[1] = make_uint3(0, 1, 0);
	dir[2] = make_uint3(0, 0, 1);

	uint cpos = 0;
	int threads = THREAD_NUM;
	dim3 grid;
	for(int i = 0; i < 3; ++i){
		// 1�X���b�h/�G�b�W
		grid = dim3((num_edge[i]+threads-1)/threads, 1, 1);
		if(grid.x > 65535){
			grid.y = (grid.x+32768-1)/32768;
			grid.x = 32768;
		}
		CalVertexEdge<<<grid, threads>>>(((float4*)dEdgeVrts)+cpos, dEdgeOcc+cpos, 
										  dVolume, dir[i], edge_size[i], grid_size, 
										  num_voxels, num_edge[i], grid_width, grid_min, threshold);

		cpos += num_edge[i];
	}
	RX_CUERROR("CalVertexEdge failed");
	RX_CUCHECK(cudaThreadSynchronize());


	// ���_�����l�߂�
	CuScan(dEdgeOccScan, dEdgeOcc, num_edge[3]);

	RX_CUCHECK(cudaMemcpy((void*)&lval, (void*)(dEdgeOcc+num_edge[3]-1), sizeof(uint), cudaMemcpyDeviceToHost));
	RX_CUCHECK(cudaMemcpy((void*)&lsval, (void*)(dEdgeOccScan+num_edge[3]-1), sizeof(uint), cudaMemcpyDeviceToHost));
	nvrts = lval + lsval;

	if(nvrts == 0){
		return;
	}

	grid = dim3((num_edge[3]+threads-1)/threads, 1, 1);
	if(grid.x > 65535){
		grid.y = (grid.x+32768-1)/32768;
		grid.x = 32768;
	}

	// compact edge vertex array
	CompactEdges<<<grid, threads>>>((float4*)dCompactedEdgeVrts, dEdgeOcc, dEdgeOccScan, 
									 (float4*)dEdgeVrts, num_edge[3]);
	RX_CUERROR("CompactEdges failed");
	RX_CUCHECK(cudaThreadSynchronize());

	//cudaMemcpy(m_f4VertPos, m_dfCompactedEdgeVrts, sizeof(float4)*m_uNumTotalVerts, cudaMemcpyDeviceToHost);
}


/*!
 * �ʑ����쐬
 * @param[out] dTris �|���S���̊􉽏��
 * @param[in] dVoxBit �O���b�h8���_�̉A�֐��l��臒l�ȏォ�ǂ������e�r�b�g�Ɋi�[�����ϐ�
 * @param[in] dVoxTNumScan �O���b�h���Ƃ̎O�p�`���b�V������Scan
 * @param[in] dCompactedVox ���b�V�����܂ރO���b�h���
 * @param[in] dEdgeOccScan �G�b�W�Ƀ��b�V�����_���܂ނ��ǂ�����Scan
 * @param[in] edge_size[3] �G�b�W��(nx,ny,nz)
 * @param[in] num_edge[4] ���G�b�W��
 * @param[in] grid_size  �O���b�h��(nx,ny,nz)
 * @param[in] num_voxels ���O���b�h��
 * @param[in] grid_width �O���b�h��
 * @param[in] threshold  ���b�V����臒l
 * @param[in] num_active_voxels ���b�V�����܂ރ{�N�Z���̐�
 * @param[in] nvrts ���_��
 * @param[in] ntris ���b�V����
*/
void CuMCCalTri(uint *dTris, uint *dVoxBit, uint *dVoxTNumScan, uint *dCompactedVox, 
				uint *dEdgeOccScan, uint3 edge_size[3], uint num_edge[4], 
				uint3 grid_size, uint num_voxels, float3 grid_width, float threshold, 
				uint num_active_voxels, uint nvrts, uint ntris)
{
	// 1�X���b�h/�A�N�e�B�u�{�N�Z��
	int threads = NTHREADS;
	dim3 grid((num_active_voxels+threads-1)/threads, 1, 1);
	if(grid.x > 65535){
		grid.y = (grid.x+32768-1)/32768;
		grid.x = 32768;
	}

	// �ʑ����쐬
	uint3 numEdge = make_uint3(num_edge[0], num_edge[1], num_edge[2]);
	GenerateTriangles3<<<grid, threads>>>((uint3*)dTris, dVoxTNumScan, dEdgeOccScan, dVoxBit, 
										  edge_size[0], edge_size[1], edge_size[2], numEdge, dCompactedVox, grid_size, 
										  num_active_voxels, grid_width, threshold, nvrts, ntris);
	RX_CUERROR("GenerateTriangles3 failed");
	RX_CUCHECK(cudaThreadSynchronize());

	//cudaMemcpy(m_u3TriIdx, dTris, sizeof(uint3)*m_uMaxVerts*3, cudaMemcpyDeviceToHost);
}

/*!
 * �@�����쐬
 * @param[out] dNrms ���_�@��
 * @param[in] dTris �|���S���̊􉽏��
 * @param[in] dCompactedEdgeVrts ���Ԃ��߂����_���
 * @param[in] nvrts ���_��
 * @param[in] ntris ���b�V����
*/
void CuMCCalNrm(float *dNrms, uint *dTris, float *dCompactedEdgeVrts, uint nvrts, uint ntris)
{
	RX_CUCHECK(cudaMemset((void*)dNrms, 0, sizeof(float3)*nvrts));

	// 1�X���b�h/���b�V��
	int threads = NTHREADS;
	dim3 grid((ntris+threads-1)/threads, 1, 1);
	if(grid.x > 65535){
		grid.y = (grid.x+32768-1)/32768;
		grid.x = 32768;
	}

	// ���b�V���@���̌v�Z�ƒ��_�ւ̒~��
	CalVertexNormalKernel<<<grid, threads>>>((float4*)dCompactedEdgeVrts, (uint3*)dTris, (float3*)dNrms, nvrts, ntris);
	RX_CUCHECK(cudaThreadSynchronize());


	// 1�X���b�h/���_
	grid = dim3((nvrts+threads-1)/threads, 1, 1);
	if(grid.x > 65535){
		grid.y = (grid.x+32768-1)/32768;
		grid.x = 32768;
	}

	// ���_�@���̐��K��
	NormalizeKernel<<<grid, threads>>>((float3*)dNrms, nvrts);
	RX_CUCHECK(cudaThreadSynchronize());
}


#else

void CuMCCreateMesh(GLuint pvbo, GLuint nvbo, float threshold, unsigned int &nvrts, unsigned int &ntris)
{
	// MARK:CuMCCreateMesh
	int threads = 128;
	dim3 grid(m_uNumVoxels/threads, 1, 1);
	// get around maximum grid size of 65535 in each dimension
	if(grid.x > 65535){
		grid.y = grid.x/32768;
		grid.x = 32768;
	}

	//
	// �e�{�N�Z����8���_�̓��O�𔻒�
	//
	classifyVoxel<<<grid, threads>>>(m_duVoxelVerts, m_duVoxelOccupied, g_dfVolume, 
									 m_u3GridSize, m_uNumVoxels, m_f3VoxelH, threshold);
	cutilCheckMsg("classifyVoxel failed");


#if SKIP_EMPTY_VOXELS
	// �{�N�Z����L�z���scan
	ThrustScanWrapper(m_duVoxelOccupiedScan, m_duVoxelOccupied, m_uNumVoxels);

	// read back values to calculate total number of non-empty voxels
	// since we are using an exclusive scan, the total is the last value of
	// the scan result plus the last value in the input array
	{
		uint lval, lsval;
		RX_CUCHECK(cudaMemcpy((void *) &lval, 
					   (void *) (m_duVoxelOccupied + m_uNumVoxels-1), 
					   sizeof(uint), cudaMemcpyDeviceToHost));
		RX_CUCHECK(cudaMemcpy((void *) &lsval, 
					   (void *) (m_duVoxelOccupiedScan + m_uNumVoxels-1), 
					   sizeof(uint), cudaMemcpyDeviceToHost));
		m_uNumActiveVoxels = lval + lsval;
	}

	if (m_uNumActiveVoxels==0) {
		// return if there are no full voxels
		m_uNumTotalVerts = 0;
		return;
	}

	// compact voxel index array
	compactVoxels<<<grid, threads>>>(m_duCompactedVoxelArray, m_duVoxelOccupied, m_duVoxelOccupiedScan, m_uNumVoxels);
	cutilCheckMsg("compactVoxels failed");

#endif // SKIP_EMPTY_VOXELS

	// scan voxel vertex count array
	ThrustScanWrapper(m_duVoxelVertsScan, m_duVoxelVerts, m_uNumVoxels);

	// readback total number of vertices
	{
		uint lval, lsval;
		RX_CUCHECK(cudaMemcpy((void *) &lval, 
					   (void *) (m_duVoxelVerts + m_uNumVoxels-1), 
					   sizeof(uint), cudaMemcpyDeviceToHost));
		RX_CUCHECK(cudaMemcpy((void *) &lsval, 
					   (void *) (m_duVoxelVertsScan + m_uNumVoxels-1), 
					   sizeof(uint), cudaMemcpyDeviceToHost));
		m_uNumTotalVerts = lval+lsval;

		nvrts = m_uNumTotalVerts;
		ntris = nvrts/3;
	}


	//
	// �O�p�`���b�V���쐬
	//

	// ���_�Ɩ@���o�b�t�@
	float4 *d_pos = 0, *d_nrm = 0;
	if(pvbo){
		RX_CUCHECK(cudaGLMapBufferObject((void**)&d_pos, pvbo));
		RX_CUCHECK(cudaGLMapBufferObject((void**)&d_nrm, nvbo));
	}
	else{
		d_pos = m_df4Vrts;
		d_nrm = m_df4Nrms;
	}

#if SKIP_EMPTY_VOXELS
	dim3 grid2((int) ceil(m_uNumActiveVoxels / (float) NTHREADS), 1, 1);
#else
	dim3 grid2((int) ceil(m_uNumVoxels / (float) NTHREADS), 1, 1);
#endif

	while(grid2.x > 65535) {
		grid2.x/=2;
		grid2.y*=2;
	}

	generateTriangles2<<<grid2, NTHREADS>>>(d_pos, d_nrm, 
						   m_duCompactedVoxelArray, m_duVoxelVertsScan, g_dfVolume, 
						   m_u3GridSize, m_f3VoxelH, m_f3VoxelMin, threshold, m_uNumActiveVoxels, m_uMaxVerts);
	cutilCheckMsg("generateTriangles2 failed");

	if(pvbo){
		RX_CUCHECK(cudaGLUnmapBufferObject(nvbo));
		RX_CUCHECK(cudaGLUnmapBufferObject(pvbo));
	}


}

#endif



}   // extern "C"
