/*! 
  @file rx_ssm.cu
	
  @brief SSM�@�ɂ�郁�b�V������
		- M. Muller et al., Screen space meshes, SCA2007, 2007. 

  @author Makoto Fujisawa
  @date 2011-05
*/
// FILE --rx_ssm.cu--



//-----------------------------------------------------------------------------
// �C���N���[�h�t�@�C��
//-----------------------------------------------------------------------------
#include <cstdio>
#include <GL/glew.h>
#include <GL/glut.h>

#include "rx_ssm_kernel.cu"
#include "rx_ssm_tables.h"

#include "rx_cu_funcs.cuh"


//-----------------------------------------------------------------------------
// MARK:�O���[�o���ϐ�
//-----------------------------------------------------------------------------
// MS�@�̃e�[�u��
uint* g_puSSMeshTable = 0;
uint* g_puSSEdgeTable = 0;
uint* g_puSSNodeTable = 0;
uint* g_puSSVRotTable = 0;


//-----------------------------------------------------------------------------
// CUDA�֐�
//-----------------------------------------------------------------------------
extern "C"
{


/*!
 * ���_���̏�����
 * @param[in] pack
 */
void CuInitVPackf(rxVPackf &pack, int size)
{
	if(pack.dPos) RX_CUCHECK(cudaFree(pack.dPos));
	RX_CUCHECK(cudaMalloc((void**)&pack.dPos, sizeof(float)*3*size));
	if(pack.dOcc) RX_CUCHECK(cudaFree(pack.dOcc));
	RX_CUCHECK(cudaMalloc((void**)&pack.dOcc, sizeof(uint)*size));
	if(pack.dOccScan) RX_CUCHECK(cudaFree(pack.dOccScan));
	RX_CUCHECK(cudaMalloc((void**)&pack.dOccScan, sizeof(uint)*size));
}

/*!
 * ���_���̏�����
 * @param[in] pack
 */
void CuInitVPacke(rxVPacke &pack, int size)
{
	if(pack.dPos) RX_CUCHECK(cudaFree(pack.dPos));
	RX_CUCHECK(cudaMalloc((void**)&pack.dPos, sizeof(rxSSEdgeG)*3*size));
	if(pack.dOcc) RX_CUCHECK(cudaFree(pack.dOcc));
	RX_CUCHECK(cudaMalloc((void**)&pack.dOcc, sizeof(uint)*size));
	if(pack.dOccScan) RX_CUCHECK(cudaFree(pack.dOccScan));
	RX_CUCHECK(cudaMalloc((void**)&pack.dOccScan, sizeof(uint)*size));
}

/*!
 * ���_���̏���
 * @param[in] pack
 */
void CuCleanVPackf(rxVPackf &pack)
{
	if(pack.dPos) RX_CUCHECK(cudaFree(pack.dPos));
	if(pack.dCompactedPos) RX_CUCHECK(cudaFree(pack.dCompactedPos));
	if(pack.dOcc) RX_CUCHECK(cudaFree(pack.dOcc));
	if(pack.dOccScan) RX_CUCHECK(cudaFree(pack.dOccScan));
}

/*!
 * ���_���̏���
 * @param[in] pack
 */
void CuCleanVPacke(rxVPacke &pack)
{
	if(pack.dPos) RX_CUCHECK(cudaFree(pack.dPos));
	if(pack.dCompactedPos) RX_CUCHECK(cudaFree(pack.dCompactedPos));
	if(pack.dOcc) RX_CUCHECK(cudaFree(pack.dOcc));
	if(pack.dOccScan) RX_CUCHECK(cudaFree(pack.dOccScan));
}

/*!
 * MS�@�̃e�[�u��
 */
void CuInitTable(void)
{
	cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc(32, 0, 0, 0, cudaChannelFormatKindSigned);

	RX_CUCHECK(cudaMalloc((void**) &g_puSSMeshTable, 256*19*sizeof(int)));
	RX_CUCHECK(cudaMemcpy(g_puSSMeshTable, g_MeshTable, 256*19*sizeof(int), cudaMemcpyHostToDevice) );
	RX_CUCHECK(cudaBindTexture(0, g_texSSMeshTable, g_puSSMeshTable, channelDesc) );

	RX_CUCHECK(cudaMalloc((void**) &g_puSSEdgeTable, 16*4*sizeof(int)));
	RX_CUCHECK(cudaMemcpy(g_puSSEdgeTable, g_EdgeTable, 16*4*sizeof(int), cudaMemcpyHostToDevice) );
	RX_CUCHECK(cudaBindTexture(0, g_texSSEdgeTable, g_puSSEdgeTable, channelDesc) );

	RX_CUCHECK(cudaMalloc((void**) &g_puSSNodeTable, 4*8*sizeof(int)));
	RX_CUCHECK(cudaMemcpy(g_puSSNodeTable, g_NodeTable, 4*8*sizeof(int), cudaMemcpyHostToDevice) );
	RX_CUCHECK(cudaBindTexture(0, g_texSSNodeTable, g_puSSNodeTable, channelDesc) );

	RX_CUCHECK(cudaMalloc((void**) &g_puSSVRotTable, 4*14*sizeof(int)));
	RX_CUCHECK(cudaMemcpy(g_puSSVRotTable, g_VrtRotTable, 4*14*sizeof(int), cudaMemcpyHostToDevice) );
	RX_CUCHECK(cudaBindTexture(0, g_texSSVRotTable, g_puSSVRotTable, channelDesc) );
}

/*!
 * MS�@�̃e�[�u���̔j��
 */
void CuCleanTable(void)
{
	RX_CUCHECK(cudaFree(g_puSSMeshTable));
	RX_CUCHECK(cudaFree(g_puSSEdgeTable));
	RX_CUCHECK(cudaFree(g_puSSNodeTable));
	RX_CUCHECK(cudaFree(g_puSSVRotTable));
}


//-----------------------------------------------------------------------------
// MARK:SSM
//-----------------------------------------------------------------------------

/*!
 * �p�����[�^���R���X�^���g�������ɓ]��
 */
void CuSetSSMParameters(rxSsmParams *hostParams)
{
	RX_CUCHECK( cudaMemcpyToSymbol(g_cSSMParams, hostParams, sizeof(rxSsmParams)) );
}
/*!
 * float�^�̔z���������
 * @param[in] dArray �l����������float�^�̔z��
 * @param[in] val �������l
 * @param[in] n �z��̃T�C�Y
 */
void CuInitFloatArray(float* dArray, float val, int n)
{
	// 1�X���b�h/�v�f
	uint grid, block;	// �O���b�h���u���b�N���C�u���b�N���X���b�h��
	block = THREAD_NUM;
	grid = DivCeil(n, block);

	// �J�[�l�����s
	initFloatArray<<<grid, block>>>(dArray, val, n);

	RX_CUERROR("initDepthMap kernel execution failed");	// �J�[�l�����s�G���[�`�F�b�N
	RX_CUCHECK(cudaThreadSynchronize());		// �S�ẴX���b�h���I���̂�҂�
}

/*!
 * �p�[�e�B�N�����a�̔z���������
 * @param[in] drad �o�͔��a�z��(�T�C�Y=n*3)
 * @param[in] dval ���͔��a�z��(�T�C�Y=n)
 * @param[in] n �z��̃T�C�Y
 */
void CuInitRadArray(float* drad, float* dval, int n)
{
	// 1�X���b�h/�v�f
	uint grid, block;	// �O���b�h���u���b�N���C�u���b�N���X���b�h��
	block = THREAD_NUM;
	grid = DivCeil(n, block);

	// �J�[�l�����s
	initRadArray<<<grid, block>>>((float3*)drad, dval, n);

	RX_CUERROR("initDepthMap kernel execution failed");	// �J�[�l�����s�G���[�`�F�b�N
	RX_CUCHECK(cudaThreadSynchronize());		// �S�ẴX���b�h���I���̂�҂�
}

/*!
 * �f�v�X�}�b�v�����ŏ�����
 * @param[in] dDMap �f�v�X�}�b�v((nx+1)x(ny+1))
 * @param[in] nx,ny �O���b�h�𑜓x
 */
void CuInitDepthMap(float* dDMap, int nx, int ny)
{
	int n = (nx+1)*(ny+1);

	// 1�X���b�h/�s�N�Z��
	uint grid, block;	// �O���b�h���u���b�N���C�u���b�N���X���b�h��
	block = THREAD_NUM;
	grid = DivCeil(n, block);

	// �J�[�l�����s
	initDepthMap<<<grid, block>>>(dDMap, n);

	RX_CUERROR("initDepthMap kernel execution failed");	// �J�[�l�����s�G���[�`�F�b�N
	RX_CUCHECK(cudaThreadSynchronize());		// �S�ẴX���b�h���I���̂�҂�
}

/*!
 * �p�[�e�B�N������f�v�X�}�b�v�̌v�Z
 * @param[in] dDMap �f�v�X�}�b�v((nx+1)x(ny+1)xfloat)
 * @param[inout] dPrtPos �p�[�e�B�N�����W(pnumxfloat3)�C�X�N���[���X�y�[�X�ł̍��W��Ԃ�
 * @param[out] dPrtRad �X�N���[���X�y�[�X�ł̃p�[�e�B�N�����a(pnumxfloat3)
 * @param[in] pnum �p�[�e�B�N����
 * @param[in] pdim 1�p�[�e�B�N���̃������T�C�Y(3 or 4)
 * @param[in] tr �p�[�e�B�N�����a�v�Z�p�W��(�v�f��3)
 * @param[in] pmv �������e�ƃ��f���r���[�s����|�����s��(�v�f��16=4x4)
 * @param[in] W,H �X�N���[���̉𑜓x
 * @param[in] dw,dh �O���b�h��
 * @param[in] nx,ny �O���b�h�𑜓x
 */
void CuCreateDepthMap(float* dDMap, float* dPrtPos, float* dPrtRad, int pnum, int pdim, float* tr, float* pmv, 
					  int W, int H, float dw, float dh, int nx, int ny)
{
	float3 tr3;
	tr3.x = tr[0];
	tr3.y = tr[1];
	tr3.z = tr[2];

	matrix4x4 mat_pmv;
	for(int i = 0; i < 4; ++i){
		mat_pmv.e[i].x = pmv[4*i+0];
		mat_pmv.e[i].y = pmv[4*i+1];
		mat_pmv.e[i].z = pmv[4*i+2];
		mat_pmv.e[i].w = pmv[4*i+3];
	}

	// 1�X���b�h/�p�[�e�B�N��
	uint grid, block;	// �O���b�h���u���b�N���C�u���b�N���X���b�h��
	block = THREAD_NUM;
	grid = DivCeil(pnum, block);

	// �J�[�l�����s
	if(pdim == 3){
		calDepthMap<<<grid, block>>>(dDMap, (float3*)dPrtPos, (float3*)dPrtRad, pnum, tr3, mat_pmv, 
									 (float)W, (float)H, dw, dh, nx, ny);
	}
	else{
		calDepthMap<<<grid, block>>>(dDMap, (float4*)dPrtPos, (float3*)dPrtRad, pnum, tr3, mat_pmv, 
									 (float)W, (float)H, dw, dh, nx, ny);
	}

	RX_CUERROR("calDepthMap kernel execution failed");	// �J�[�l�����s�G���[�`�F�b�N
	RX_CUCHECK(cudaThreadSynchronize());		// �S�ẴX���b�h���I���̂�҂�
}

/*!
 * �f�v�X�}�b�v�ɕ��������{��
 * @param[in] dDMap �f�v�X�}�b�v((nx)x(ny))
 * @param[in] nx,ny �}�b�v�𑜓x
 * @param[in] n_filter �t�B���^��
 * @param[in] binomials �񍀉��Z�W��(n_filter = 0, 1, 2,..,RX_MAX_FILTER_SIZE �̌W�������ԂɊi�[����Ă���)
 * @param[in] zmax �֊s�ƂȂ�f�v�X����臒l
 */
void CuDepthSmoothing(float* dDMap, int nx, int ny, int n_filter, float *binomials, float zmax)
{
	int b = (n_filter+1)*(n_filter+1);
	cudaMemcpyToSymbol(g_cBinomials, binomials, sizeof(float)*b);

	// 1�X���b�h/�s�N�Z��(�e�X���b�h���}�b�v�̃s�N�Z���ɑΉ�����悤��2�����I�ɔz�u)
	dim3 block(BLOCK_SIZE, BLOCK_SIZE);
	dim3 grid((nx+block.x-1)/block.x, (ny+block.y-1)/block.y);

	// �u���b�N���Ƃ̃V�F�A�[�h�������T�C�Y(�ő�16KB)
	int sbytes = (block.x+2*n_filter)*(block.y+2*n_filter)*sizeof(float);

	// �J�[�l�����s(x����������)
	smoothDepthMapX<<< grid, block, sbytes >>>(dDMap, nx, ny, n_filter, block.x+2*n_filter, zmax);

	RX_CUERROR("smoothDepthMapX kernel execution failed");	// �J�[�l�����s�G���[�`�F�b�N
	RX_CUCHECK(cudaThreadSynchronize());		// �S�ẴX���b�h���I���̂�҂�

	// �J�[�l�����s(y����������)
	smoothDepthMapY<<< grid, block, sbytes >>>(dDMap, nx, ny, n_filter, block.x+2*n_filter, zmax);

	RX_CUERROR("smoothDepthMapY kernel execution failed");	// �J�[�l�����s�G���[�`�F�b�N
	RX_CUCHECK(cudaThreadSynchronize());		// �S�ẴX���b�h���I���̂�҂�
}


/*!
 * �֊s�G�b�W�̌��o��front edge vertex�̌v�Z
 * @param[in] dDMap �f�v�X�}�b�v(�T�C�Y = (nx+1)x(ny+1))
 * @param[in] nx,ny �}�b�v�𑜓x
 * @param[in] dw,dh �O���b�h��
 * @param[in] zmax �֊s�ƂȂ�f�v�X����臒l
 * @param[out] dNodeVrts �m�[�h���_(�T�C�Y = (nx+1)x(ny+1))
 * @param[out] dMGrid ���b�V�������p�O���b�h(�T�C�Y = (nx)x(ny))
 * @param[out] num_node_vertex �m�[�h���_��
 */
void CuCalNodeVertex(float* dDMap, int nx, int ny, float dw, float dh, float zmax, 
					 rxVPackf &dNodeVrts, rxSSGridG* dMGrid, int &num_node_vertex)
{
	uint block1, grid1;
	dim3 block2, grid2;
	unsigned int size = (nx+1)*(ny+1);
	uint last_val, last_scan_val;

	RX_CUCHECK(cudaMemset((void*)dNodeVrts.dOcc, 0, size*sizeof(uint)));

	//
	// �m�[�h���_
	//
	// 1�X���b�h/�m�[�h
	block2 = dim3(BLOCK_SIZE, BLOCK_SIZE);
	grid2 = dim3(((nx+1)+block2.x-1)/block2.x, ((ny+1)+block2.y-1)/block2.y);
	
	// �J�[�l�����s
	calNodeVertex<<< grid2, block2 >>>(dDMap, nx+1, ny+1, dw, dh, zmax, (float3*)dNodeVrts.dPos, dNodeVrts.dOcc);

	RX_CUERROR("calNodeVertex kernel execution failed");	// �J�[�l�����s�G���[�`�F�b�N
	RX_CUCHECK(cudaThreadSynchronize());		// �S�ẴX���b�h���I���̂�҂�


	//
	// �m�[�h���_�̃p�b�L���O
	//
	// �m�[�h���_������Ȃ�1,�����łȂ��Ȃ�0���i�[���ꂽ�z���Scan
	CuScan(dNodeVrts.dOccScan, dNodeVrts.dOcc, size);

	// Exclusive scan (�Ō�̗v�f��0�Ԗڂ���n-2�Ԗڂ܂ł̍��v�ɂȂ��Ă���)�Ȃ̂ŁC
	// Scan�O�z��̍Ō�(n-1�Ԗ�)�̗v�f�ƍ��v���邱�Ƃŗ֊s�G�b�W�����v�Z
	RX_CUCHECK(cudaMemcpy((void*)&last_val, (void*)(dNodeVrts.dOcc+size-1), sizeof(uint), cudaMemcpyDeviceToHost));
	RX_CUCHECK(cudaMemcpy((void *) &last_scan_val, (void*)(dNodeVrts.dOccScan+size-1), sizeof(uint), cudaMemcpyDeviceToHost));
	num_node_vertex = last_val+last_scan_val;

	// �l�߂��֊s�G�b�W�����i�[����̈�̊m��
	if(dNodeVrts.dCompactedPos) RX_CUCHECK(cudaFree(dNodeVrts.dCompactedPos));
	RX_CUCHECK(cudaMalloc((void**)&dNodeVrts.dCompactedPos, num_node_vertex*3*sizeof(float)));

	// �֊s�G�b�W���l�߂�
	block1 = THREAD_NUM;
	grid1 = DivCeil(size, block1);

	// �J�[�l�����s
	compactNodeVertex<<<grid1, block1>>>((float3*)(dNodeVrts.dCompactedPos), dNodeVrts.dOcc, dNodeVrts.dOccScan, (float3*)(dNodeVrts.dPos), size);
	
	RX_CUERROR("compactNodeVertex kernel execution failed");	// �J�[�l�����s�G���[�`�F�b�N
	RX_CUCHECK(cudaThreadSynchronize());		// �S�ẴX���b�h���I���̂�҂�


	//
	// �m�[�h���_�C���f�b�N�X�����b�V���O���b�h�Ɋi�[
	//
	// 1�X���b�h/�O���b�h
	block2 = dim3(BLOCK_SIZE, BLOCK_SIZE);
	grid2 = dim3(((nx)+block2.x-1)/block2.x, ((ny)+block2.y-1)/block2.y);
		
	// �J�[�l�����s
	storeNodeVertex<<< grid2, block2 >>>(dMGrid, nx, ny, (float3*)(dNodeVrts.dCompactedPos), dNodeVrts.dOcc, dNodeVrts.dOccScan);

	RX_CUERROR("storeNodeVertex kernel execution failed");	// �J�[�l�����s�G���[�`�F�b�N
	RX_CUCHECK(cudaThreadSynchronize());		// �S�ẴX���b�h���I���̂�҂�
}



/*!
 * �֊s�G�b�W�̌��o��front edge vertex�̌v�Z
 * @param[in] dDMap �f�v�X�}�b�v((nx+1)x(ny+1))
 * @param[in] nx,ny �}�b�v�𑜓x
 * @param[in] dw,dh �O���b�h��
 * @param[in] zmax �֊s�ƂȂ�f�v�X����臒l
 * @param[in] dPrtPos �p�[�e�B�N�����W(pnumxfloat3)�C�X�N���[���X�y�[�X�ł̍��W��Ԃ�
 * @param[in] dPrtRad �X�N���[���X�y�[�X�ł̃p�[�e�B�N�����a(pnumxfloat3)
 * @param[in] pnum �p�[�e�B�N����
 * @param[in] pdim 1�p�[�e�B�N���̃������T�C�Y(3 or 4)
 * @param[in] W,H �X�N���[���̉𑜓x
 * @param[out] dEdge �G�b�W���(nx*(ny+1)+(nx+1)*ny)
 * @param[out] dCompactedEdge �֊s�G�b�W���(�G�b�W��񂩂�֊s�G�b�W�݂̂��l�߂���)
 * @param[out] dEdgeSil �֊s�G�b�W�Ȃ��1,�����łȂ����0���i�[����z��(nx*(ny+1)+(nx+1)*ny)
 * @param[out] dEdgeSilScan dEdgeSil��Prefix sum (nx*(ny+1)+(nx+1)*ny)
 * @param[out] dEdgeVrts �G�b�W���_(�T�C�Y = (nx*(ny+1)+(nx+1)*ny))
 * @param[out] dBackVrts �w�ʃG�b�W���_(�T�C�Y = (nx*(ny+1)+(nx+1)*ny))
 * @param[out] dMGrid ���b�V�������p�O���b�h(�T�C�Y = (nx)x(ny))
 * @param[in] num_nv �m�[�h���_��
 * @param[out] num_edge �֊s�G�b�W��
 * @param[out] num_front_edge �G�b�W���_��(�O��)
 * @param[out] num_back_edge �G�b�W���_��(�w��)
 */
void CuDetectSilhouetteEdgeVertex(float* dDMap, int nx, int ny, float dw, float dh, float zmax, 
								  float* dPrtPos, float* dPrtRad, int pnum, int pdim, int W, int H, 
								  rxVPacke &dEdge, rxVPackf &dEdgeVrts, rxVPackf &dBackVrts, rxSSGridG* dMGrid, 
								  int num_nv, int &num_edge, int &num_front_edge, int &num_back_edge)
{
	uint block1, grid1;
	dim3 block2, grid2;
	unsigned int size = nx*(ny+1)+(nx+1)*ny;
	int offset = (nx)*(ny+1);
	uint last_val, last_scan_val;

	//
	// �֊s�G�b�W���o
	//
	RX_CUCHECK(cudaMemset((void*)dEdge.dOcc, 0, size*sizeof(uint)));

	// x�����֊s�G�b�W���o
	// 1�X���b�h/�G�b�W
	block2 = dim3(BLOCK_SIZE, BLOCK_SIZE);
	grid2 = dim3(((nx)+block2.x-1)/block2.x, ((ny+1)+block2.y-1)/block2.y);
	
	// �J�[�l�����s
	detectSilhouetteEdgeX<<< grid2, block2 >>>(dDMap, nx, ny+1, dw, dh, zmax, dEdge.dPos, dEdge.dOcc);

	RX_CUERROR("detectSilhouetteEdgeX kernel execution failed");	// �J�[�l�����s�G���[�`�F�b�N
	RX_CUCHECK(cudaThreadSynchronize());		// �S�ẴX���b�h���I���̂�҂�
	
	// y�����֊s�G�b�W���o
	// 1�X���b�h/�G�b�W
	block2 = dim3(BLOCK_SIZE, BLOCK_SIZE);
	grid2 = dim3(((nx+1)+block2.x-1)/block2.x, ((ny)+block2.y-1)/block2.y);
		
	// �J�[�l�����s
	detectSilhouetteEdgeY<<< grid2, block2 >>>(dDMap, nx+1, ny, dw, dh, zmax, dEdge.dPos+offset, dEdge.dOcc+offset);

	RX_CUERROR("detectSilhouetteEdgeY kernel execution failed");	// �J�[�l�����s�G���[�`�F�b�N
	RX_CUCHECK(cudaThreadSynchronize());		// �S�ẴX���b�h���I���̂�҂�
	
	//
	// �֊s�G�b�W�̉ۂ�Prefix Sum�쐬���āC�G�b�W�����l�߂�
	//
	// �֊s�G�b�W�Ȃ�1,�����łȂ��Ȃ�0���i�[���ꂽ�z���Scan
	CuScan(dEdge.dOccScan, dEdge.dOcc, size);

	// Exclusive scan (�Ō�̗v�f��0�Ԗڂ���n-2�Ԗڂ܂ł̍��v�ɂȂ��Ă���)�Ȃ̂ŁC
	// Scan�O�z��̍Ō�(n-1�Ԗ�)�̗v�f�ƍ��v���邱�Ƃŗ֊s�G�b�W�����v�Z
	RX_CUCHECK(cudaMemcpy((void*)&last_val, (void*)(dEdge.dOcc+size-1), sizeof(uint), cudaMemcpyDeviceToHost));
	RX_CUCHECK(cudaMemcpy((void*)&last_scan_val, (void*)(dEdge.dOccScan+size-1), sizeof(uint), cudaMemcpyDeviceToHost));
	num_edge = last_val+last_scan_val;

	// �l�߂��֊s�G�b�W�����i�[����̈�̊m��
	if(dEdge.dCompactedPos) RX_CUCHECK(cudaFree(dEdge.dCompactedPos));
	RX_CUCHECK(cudaMalloc((void**)&dEdge.dCompactedPos, num_edge*sizeof(rxSSEdgeG)));

	// �֊s�G�b�W���l�߂�
	block1 = THREAD_NUM;
	grid1 = DivCeil(size, block1);

	// �J�[�l�����s
	compactSilhouetteEdges<<<grid1, block1>>>(dEdge.dCompactedPos, dEdge.dOcc, dEdge.dOccScan, dEdge.dPos, size);
	
	RX_CUERROR("CompactSilhouetteEdges kernel execution failed");	// �J�[�l�����s�G���[�`�F�b�N
	RX_CUCHECK(cudaThreadSynchronize());		// �S�ẴX���b�h���I���̂�҂�



	//
	// �֊s�G�b�W��front vertex���Z�o
	//
	// �G�b�W���_�L�����z��̏�����
	RX_CUCHECK(cudaMemset((void*)dEdgeVrts.dOcc, 0, size*sizeof(uint)));

	// 1�X���b�h/�p�[�e�B�N��
	block1 = THREAD_NUM;
	grid1 = DivCeil(pnum, block1);

	// �J�[�l�����s
	if(pdim == 3){
		calFrontEdgeVertex<<<grid1, block1>>>((float3*)dPrtPos, (float3*)dPrtRad, pnum, 
											  dEdge.dPos, dEdge.dOcc, offset, nx, ny, dw, dh, W, H, 
											  (float3*)dEdgeVrts.dPos, dEdgeVrts.dOcc);
	}
	else{
		calFrontEdgeVertex<<<grid1, block1>>>((float4*)dPrtPos, (float3*)dPrtRad, pnum, 
											  dEdge.dPos, dEdge.dOcc, offset, nx, ny, dw, dh, W, H, 
											  (float3*)dEdgeVrts.dPos, dEdgeVrts.dOcc);
	}

	RX_CUERROR("calFrontEdgeVertex kernel execution failed");	// �J�[�l�����s�G���[�`�F�b�N
	RX_CUCHECK(cudaThreadSynchronize());		// �S�ẴX���b�h���I���̂�҂�
	
	//
	// �G�b�W���_�L������Prefix Sum���쐬���āC�G�b�W���_���l�߂�
	//
	// �֊s�G�b�W�Ȃ�1,�����łȂ��Ȃ�0���i�[���ꂽ�z���Scan
	CuScan(dEdgeVrts.dOccScan, dEdgeVrts.dOcc, size);

	// �֊s�G�b�W�����v�Z
	RX_CUCHECK(cudaMemcpy((void*)&last_val, (void*)(dEdgeVrts.dOcc+size-1), sizeof(uint), cudaMemcpyDeviceToHost));
	RX_CUCHECK(cudaMemcpy((void*)&last_scan_val, (void*)(dEdgeVrts.dOccScan+size-1), sizeof(uint), cudaMemcpyDeviceToHost));
	num_front_edge = last_val+last_scan_val;

	// �l�߂��֊s�G�b�W�����i�[����̈�̊m��
	if(dEdgeVrts.dCompactedPos) RX_CUCHECK(cudaFree(dEdgeVrts.dCompactedPos));
	RX_CUCHECK(cudaMalloc((void**)&dEdgeVrts.dCompactedPos, num_front_edge*sizeof(float)*3));

	// �֊s�G�b�W���l�߂�
	block1 = THREAD_NUM;
	grid1 = DivCeil(size, block1);

	// �J�[�l�����s
	compactEdgeVertices<<<grid1, block1>>>((float3*)dEdgeVrts.dCompactedPos, dEdgeVrts.dOcc, dEdgeVrts.dOccScan, (float3*)dEdgeVrts.dPos, dEdge.dPos, size);
	
	RX_CUERROR("compactEdgeVertices kernel execution failed");	// �J�[�l�����s�G���[�`�F�b�N
	RX_CUCHECK(cudaThreadSynchronize());		// �S�ẴX���b�h���I���̂�҂�



	//
	// �֊s�G�b�W��back vertex���Z�o
	//
	// �G�b�W���_�L�����z��̏�����
	RX_CUCHECK(cudaMemset((void*)dBackVrts.dOcc, 0, size*sizeof(uint)));

	// x�����G�b�W
	// 1�X���b�h/�G�b�W
	block2 = dim3(BLOCK_SIZE, BLOCK_SIZE);
	grid2 = dim3(((nx)+block2.x-1)/block2.x, ((ny+1)+block2.y-1)/block2.y);
	
	// �J�[�l�����s
	calBackEdgeVertexX<<<grid2, block2>>>(dDMap, dEdge.dPos, dEdge.dOcc, nx, ny+1, dw, dh, W, H, 
										  (float3*)dEdgeVrts.dPos, (float3*)dBackVrts.dPos, dBackVrts.dOcc);

	RX_CUERROR("calBackEdgeVertexX kernel execution failed");	// �J�[�l�����s�G���[�`�F�b�N
	RX_CUCHECK(cudaThreadSynchronize());		// �S�ẴX���b�h���I���̂�҂�

	offset = (nx)*(ny+1);

	// y�����G�b�W
	// 1�X���b�h/�G�b�W
	block2 = dim3(BLOCK_SIZE, BLOCK_SIZE);
	grid2 = dim3(((nx+1)+block2.x-1)/block2.x, ((ny)+block2.y-1)/block2.y);
	
	// �J�[�l�����s
	calBackEdgeVertexY<<<grid2, block2>>>(dDMap, dEdge.dPos, dEdge.dOcc, nx+1, ny, dw, dh, W, H, 
										  (float3*)dEdgeVrts.dPos, (float3*)dBackVrts.dPos, dBackVrts.dOcc, offset);

	RX_CUERROR("calBackEdgeVertexY kernel execution failed");	// �J�[�l�����s�G���[�`�F�b�N
	RX_CUCHECK(cudaThreadSynchronize());		// �S�ẴX���b�h���I���̂�҂�

	//
	// �G�b�W���_�L������Prefix Sum���쐬���āC�G�b�W���_���l�߂�
	//
	// �֊s�G�b�W�Ȃ�1,�����łȂ��Ȃ�0���i�[���ꂽ�z���Scan
	CuScan(dBackVrts.dOccScan, dBackVrts.dOcc, size);

	// �֊s�G�b�W�����v�Z
	RX_CUCHECK(cudaMemcpy((void*)&last_val, (void*)(dBackVrts.dOcc+size-1), sizeof(uint), cudaMemcpyDeviceToHost));
	RX_CUCHECK(cudaMemcpy((void*)&last_scan_val, (void*)(dBackVrts.dOccScan+size-1), sizeof(uint), cudaMemcpyDeviceToHost));
	num_back_edge = last_val+last_scan_val;

	if(num_back_edge){
		// �l�߂��֊s�G�b�W�����i�[����̈�̊m��
		if(dBackVrts.dCompactedPos) RX_CUCHECK(cudaFree(dBackVrts.dCompactedPos));
		RX_CUCHECK(cudaMalloc((void**)&dBackVrts.dCompactedPos, num_back_edge*sizeof(float)*3));

		// �֊s�G�b�W���l�߂�
		block1 = THREAD_NUM;
		grid1 = DivCeil(size, block1);

		// �J�[�l�����s
		compactBackEdgeVertices<<<grid1, block1>>>((float3*)dBackVrts.dCompactedPos, dBackVrts.dOcc, dBackVrts.dOccScan, (float3*)dBackVrts.dPos, size);
	
		RX_CUERROR("compactEdgeVertices kernel execution failed");	// �J�[�l�����s�G���[�`�F�b�N
		RX_CUCHECK(cudaThreadSynchronize());		// �S�ẴX���b�h���I���̂�҂�
	}


	//
	// �G�b�W���_�C���f�b�N�X�����b�V���O���b�h�Ɋi�[
	//
	// 1�X���b�h/�O���b�h
	block2 = dim3(BLOCK_SIZE, BLOCK_SIZE);
	grid2 = dim3(((nx)+block2.x-1)/block2.x, ((ny)+block2.y-1)/block2.y);
		
	// �J�[�l�����s
	storeEdgeVertex<<< grid2, block2 >>>(dMGrid, nx, ny, offset, num_nv, num_nv+num_front_edge, 
										 (float3*)(dEdgeVrts.dCompactedPos), dEdgeVrts.dOcc, dEdgeVrts.dOccScan, 
										 (float3*)(dBackVrts.dCompactedPos), dBackVrts.dOcc, dBackVrts.dOccScan);

	RX_CUERROR("storeEdgeVertex kernel execution failed");	// �J�[�l�����s�G���[�`�F�b�N
	RX_CUCHECK(cudaThreadSynchronize());		// �S�ẴX���b�h���I���̂�҂�
}



/*!
 * ���b�V������
 * @param[in] dDMap �f�v�X�}�b�v(�T�C�Y = (nx+1)x(ny+1))
 * @param[in] nx,ny �}�b�v�𑜓x
 * @param[in] dw,dh �O���b�h��
 * @param[in] zmax �֊s�ƂȂ�f�v�X����臒l
 * @param[out] dTriNum �e�O���b�h�̃��b�V����
 * @param[out] dTriNumScan �e�O���b�h�̃��b�V������Scan
 * @param[out] dBack2Vrts �Ŕw�ʃG�b�W���_(�T�C�Y = (nx*(ny+1)+(nx+1)*ny))
 * @param[in] dVrts ���b�V�����_��
 * @param[in] num_vrts ���b�V�����_��
 * @param[out] num_back2_vrts �Ŕw�ʃG�b�W���_��
 * @param[out] dTriArray ���b�V��
 * @param[out] num_mesh ���b�V����
 */
void CuCalMesh(rxSSGridG* dMGrid, int nx, int ny, float dw, float dh, float zmax, uint* dTriNum, uint* dTriNumScan, 
			   rxVPackf &dBack2Vrts, float* dVrts, int num_vrts, int &num_back2_vrts, uint* &dTriArray, int &num_mesh)
{
	uint block1, grid1;
	dim3 block2, grid2;
	unsigned int size = nx*ny;
	unsigned int size_e = nx*(ny+1)+(nx+1)*ny;
	int offset = (nx)*(ny+1);
	uint last_val, last_scan_val;


	// �G�b�W���_�L�����z��̏�����
	RX_CUCHECK(cudaMemset((void*)dBack2Vrts.dOcc, 0, size_e*sizeof(uint)));

	// ���b�V�����z��̏�����
	RX_CUCHECK(cudaMemset((void*)dTriNum, 0, size*sizeof(uint)));

	//
	// �e�O���b�h�Ń��b�V������
	//
	// 1�X���b�h/�O���b�h
	block2 = dim3(BLOCK_SIZE, BLOCK_SIZE);
	grid2 = dim3(((nx)+block2.x-1)/block2.x, ((ny)+block2.y-1)/block2.y);

	// �J�[�l�����s
	calGridMesh<<< grid2, block2 >>>(dMGrid, nx, ny, zmax, (float3*)dVrts, num_vrts, dTriNum, 
									 (float3*)dBack2Vrts.dPos, dBack2Vrts.dOcc, offset);

	RX_CUERROR("calGridMesh kernel execution failed");	// �J�[�l�����s�G���[�`�F�b�N
	RX_CUCHECK(cudaThreadSynchronize());		// �S�ẴX���b�h���I���̂�҂�
	
	//
	// �G�b�W���_�L������Prefix Sum���쐬���āC�G�b�W���_���l�߂�
	//
	// �G�b�W���_������Ȃ�1,�����łȂ��Ȃ�0���i�[���ꂽ�z���Scan
	CuScan(dBack2Vrts.dOccScan, dBack2Vrts.dOcc, size_e);

	// �֊s�G�b�W�����v�Z
	RX_CUCHECK(cudaMemcpy((void*)&last_val, (void*)(dBack2Vrts.dOcc+size_e-1), sizeof(uint), cudaMemcpyDeviceToHost));
	RX_CUCHECK(cudaMemcpy((void*)&last_scan_val, (void*)(dBack2Vrts.dOccScan+size_e-1), sizeof(uint), cudaMemcpyDeviceToHost));
	num_back2_vrts = last_val+last_scan_val;

	if(num_back2_vrts){
		// �l�߂��֊s�G�b�W�����i�[����̈�̊m��
		if(dBack2Vrts.dCompactedPos) RX_CUCHECK(cudaFree(dBack2Vrts.dCompactedPos));
		RX_CUCHECK(cudaMalloc((void**)&dBack2Vrts.dCompactedPos, num_back2_vrts*sizeof(float)*3));

		// �֊s�G�b�W���l�߂�
		block1 = THREAD_NUM;
		grid1 = DivCeil(size_e, block1);

		// �J�[�l�����s
		compactBackEdgeVertices<<<grid1, block1>>>((float3*)dBack2Vrts.dCompactedPos, dBack2Vrts.dOcc, dBack2Vrts.dOccScan, (float3*)dBack2Vrts.dPos, size_e);
	
		RX_CUERROR("compactBackEdgeVertices kernel execution failed");	// �J�[�l�����s�G���[�`�F�b�N
		RX_CUCHECK(cudaThreadSynchronize());		// �S�ẴX���b�h���I���̂�҂�
	}

	// ���b�V�����z���Scan
	CuScan(dTriNumScan, dTriNum, size);
	
	// �����b�V�������v�Z
	RX_CUCHECK(cudaMemcpy((void*)&last_val, (void*)(dTriNum+size-1), sizeof(uint), cudaMemcpyDeviceToHost));
	RX_CUCHECK(cudaMemcpy((void*)&last_scan_val, (void*)(dTriNumScan+size-1), sizeof(uint), cudaMemcpyDeviceToHost));
	num_mesh = last_val+last_scan_val;

	// ���b�V�����i�[����̈�̊m��
	if(dTriArray) RX_CUCHECK(cudaFree(dTriArray));
	RX_CUCHECK(cudaMalloc((void**)&dTriArray, num_mesh*sizeof(uint)*3));

	// 1�X���b�h/�O���b�h
	block2 = dim3(BLOCK_SIZE, BLOCK_SIZE);
	grid2 = dim3(((nx)+block2.x-1)/block2.x, ((ny)+block2.y-1)/block2.y);

	//printf("num_vrts = %d\n", num_vrts);

	// �J�[�l�����s
	genGridMesh<<< grid2, block2 >>>(dMGrid, nx, ny, zmax, (float3*)dVrts, num_vrts, dTriArray, dTriNumScan, 
									 (float3*)dBack2Vrts.dPos, dBack2Vrts.dOcc, dBack2Vrts.dOccScan, offset, num_vrts);

	RX_CUERROR("genGridMesh kernel execution failed");	// �J�[�l�����s�G���[�`�F�b�N
	RX_CUCHECK(cudaThreadSynchronize());		// �S�ẴX���b�h���I���̂�҂�
}

/*!
 * �֊s������
 * @param[inout] dVrts ���b�V�����_��
 * @param[in] num_vrts ���b�V�����_����
 * @param[in] num_node_vrts �m�[�h���_��
 * @param[in] dTriArray ���b�V��
 * @param[in] num_mesh ���b�V����
 * @param[in] n_iter ������������
 */
void CuSilhouetteSmoothing(float* dVrts, int num_vrts, int num_node_vrts, uint* dTriArray, int num_mesh, int n_iter)
{
	uint block1, grid1;

	float4 *dAvgPos;
	RX_CUCHECK(cudaMalloc(&dAvgPos, (num_vrts-num_node_vrts)*sizeof(float4)));


	for(int l = 0; l < n_iter; ++l){
		RX_CUCHECK(cudaMemset(dAvgPos, 0, (num_vrts-num_node_vrts)*sizeof(float4)));

		// 1�X���b�h/���b�V��
		block1 = THREAD_NUM;
		grid1 = DivCeil(num_mesh, block1);

		// �J�[�l�����s(���ψʒu�̎Z�o)
		smoothSilhouette<<< grid1, block1 >>>((float3*)dVrts, num_node_vrts, dTriArray, num_mesh, dAvgPos);

		RX_CUERROR("smoothSilhouette kernel execution failed");	// �J�[�l�����s�G���[�`�F�b�N
		RX_CUCHECK(cudaThreadSynchronize());		// �S�ẴX���b�h���I���̂�҂�


		// 1�X���b�h/���_
		block1 = THREAD_NUM;
		grid1 = DivCeil(num_vrts-num_node_vrts, block1);

		// �J�[�l�����s(���ψʒu�̎Z�o)
		smoothSilhouette2<<< grid1, block1 >>>((float3*)dVrts, num_vrts, num_node_vrts, dAvgPos);

		RX_CUERROR("smoothSilhouette2 kernel execution failed");	// �J�[�l�����s�G���[�`�F�b�N
		RX_CUCHECK(cudaThreadSynchronize());		// �S�ẴX���b�h���I���̂�҂�
	}
}

/*!
 * �X�N���[���X�y�[�X���璸�_�������3D��Ԃɖ߂�
 * @param[in] dSSVrts �X�N���[���X�y�[�X���b�V�����_��
 * @param[in] num_vrts ���b�V�����_����
 * @param[in] mvq �������e�ƃ��f���r���[�t�s����|�����s��(�v�f��16=4x4)
 * @param[in] W,H �X�N���[���̉𑜓x
 * @param[out] dVrts 3D���b�V�����_��
 */
void CuTransformBack(float* dSSVrts, int num_vrts, float* mvq, float* q, int W, int H, float* dVrts)
{
	uint block1, grid1;

	float4 vec_q;
	vec_q.x = q[0];
	vec_q.y = q[1];
	vec_q.z = q[2];
	vec_q.w = q[3];

	matrix4x4 mat_mvq;
	for(int i = 0; i < 4; ++i){
		mat_mvq.e[i].x = mvq[4*i+0];
		mat_mvq.e[i].y = mvq[4*i+1];
		mat_mvq.e[i].z = mvq[4*i+2];
		mat_mvq.e[i].w = mvq[4*i+3];
	}

	// 1�X���b�h/���_
	block1 = THREAD_NUM;
	grid1 = DivCeil(num_vrts, block1);

	// �J�[�l�����s(���ψʒu�̎Z�o)
	transfromBack<<< grid1, block1 >>>((float3*)dSSVrts, (float3*)dVrts, num_vrts, mat_mvq, vec_q, W, H);

	RX_CUERROR("transfromBack kernel execution failed");	// �J�[�l�����s�G���[�`�F�b�N
	RX_CUCHECK(cudaThreadSynchronize());		// �S�ẴX���b�h���I���̂�҂�

}

/*!
 * ���_�@���v�Z
 * @param[in] dVrts ���b�V�����_��
 * @param[in] num_vrts ���b�V�����_����
 * @param[in] dTriArray ���b�V��
 * @param[in] num_mesh ���b�V����
 * @param[out] dNrms ���_�@��
 */
void CuCalVertexNormal(float* dVrts, int num_vrts, uint* dTriArray, int num_mesh, float* dNrms)
{
	uint block1, grid1;

	// �@���z��̏�����
	RX_CUCHECK(cudaMemset(dNrms, 0, num_vrts*sizeof(float3)));

	// 1�X���b�h/���b�V��
	block1 = THREAD_NUM;
	grid1 = DivCeil(num_mesh, block1);

	// �J�[�l�����s(���ψʒu�̎Z�o)
	sumFaceNormal<<< grid1, block1 >>>((float3*)dVrts, num_vrts, dTriArray, num_mesh, dNrms);

	RX_CUERROR("sumFaceNormal kernel execution failed");	// �J�[�l�����s�G���[�`�F�b�N
	RX_CUCHECK(cudaThreadSynchronize());		// �S�ẴX���b�h���I���̂�҂�


	// 1�X���b�h/���_
	block1 = THREAD_NUM;
	grid1 = DivCeil(num_vrts, block1);

	// �J�[�l�����s(���ψʒu�̎Z�o)
	normalizeNormal<<< grid1, block1 >>>((float3*)dNrms, num_vrts);

	RX_CUERROR("normalizeNormal kernel execution failed");	// �J�[�l�����s�G���[�`�F�b�N
	RX_CUCHECK(cudaThreadSynchronize());		// �S�ẴX���b�h���I���̂�҂�
}


}   // extern "C"
