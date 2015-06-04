/*! 
  @file rx_sph.cu
	
  @brief CUDA�ɂ��SPH

  @author Makoto Fujisawa
  @date 2009-08, 2011-06
*/
// FILE --rx_sph.cu--


//-----------------------------------------------------------------------------
// �C���N���[�h�t�@�C��
//-----------------------------------------------------------------------------
#include <cstdio>
#include <GL/glew.h>
#include <GL/glut.h>

#include "rx_sph_kernel.cu"
#include "rx_pbdsph_kernel.cu"

#include "rx_turb_kernel.cu"
#include "rx_aniso_kernel.cu"


//#include "rx_cu_funcs.cuh"
#include <thrust/device_vector.h>
#include <thrust/scan.h>


//-----------------------------------------------------------------------------
// MARK:�O���[�o���ϐ�
//-----------------------------------------------------------------------------
cudaArray *g_caNoiseTile = 0;
float *g_dNoiseTile[3] = {0, 0, 0};
uint g_udNoiseTileSize = 0;
uint g_uNoiseTileNum[3*3] = {0, 0, 0,  0, 0, 0,  0, 0, 0};


//-----------------------------------------------------------------------------
// CUDA�֐�
//-----------------------------------------------------------------------------
extern "C"
{
	
void CuSetParameters(rxSimParams *hostParams)
{
	// copy parameters to constant memory
	RX_CUCHECK( cudaMemcpyToSymbol(params, hostParams, sizeof(rxSimParams)) );
}

void CuClearData(void)
{
}

/*!
 * thrust::exclusive_scan�̌Ăяo��
 * @param[out] dScanData scan��̃f�[�^
 * @param[in] dData ���f�[�^
 * @param[in] num �f�[�^��
 */
void CuScanf(float* dScanData, float* dData, unsigned int num)
{
	thrust::exclusive_scan(thrust::device_ptr<float>(dData), 
						   thrust::device_ptr<float>(dData+num),
						   thrust::device_ptr<float>(dScanData));
}


/*!
 * �����Z���̃n�b�V�����v�Z
 * @param[in] 
 * @return 
 */
void CuCalcHash(uint* dGridParticleHash, uint* dSortedIndex, float* dPos, int *attr, int nprts)
{
	uint numThreads, numBlocks;
	computeGridSize(nprts, THREAD_NUM, numBlocks, numThreads);

	// �J�[�l�����s
	calcHashD2<<< numBlocks, numThreads >>>(dGridParticleHash,
										   dSortedIndex,
										   (float4*)dPos,
										   attr, 
										   nprts);
	//calcHashD<<< numBlocks, numThreads >>>(dGridParticleHash,
	//									   dSortedIndex,
	//									   (float4*)dPos,
	//									   nprts);
	
	RX_CUERROR("Kernel execution failed");	// �J�[�l���G���[�`�F�b�N
}


/*!
 * �p�[�e�B�N���z����\�[�g���ꂽ���Ԃɕ��ёւ��C
 * �e�Z���̎n�܂�ƏI���̃C���f�b�N�X������
 * @param[in] cell �p�[�e�B�N���O���b�h�f�[�^
 * @param[in] oldPos �p�[�e�B�N���ʒu
 * @param[in] oldVel �p�[�e�B�N�����x
 */
void CuReorderDataAndFindCellStart(rxParticleCell cell, float* oldPos, float* oldVel)
{
	uint numThreads, numBlocks;
	computeGridSize(cell.uNumParticles, THREAD_NUM, numBlocks, numThreads);

	RX_CUCHECK(cudaMemset(cell.dCellStart, 0xffffffff, cell.uNumCells*sizeof(uint)));

#if USE_TEX
	RX_CUCHECK(cudaBindTexture(0, dSortedPosTex, oldPos, cell.uNumParticles*sizeof(float4)));
	RX_CUCHECK(cudaBindTexture(0, dSortedVelTex, oldVel, cell.uNumParticles*sizeof(float4)));
#endif

	uint smemSize = sizeof(uint)*(numThreads+1);

	// �J�[�l�����s
	reorderDataAndFindCellStartD<<< numBlocks, numThreads, smemSize>>>(cell, (float4*)oldPos, (float4*)oldVel);

	RX_CUERROR("Kernel execution failed: CuReorderDataAndFindCellStartD");
	RX_CUCHECK(cudaThreadSynchronize());		// �S�ẴX���b�h���I���̂�҂�

#if USE_TEX
	RX_CUCHECK(cudaUnbindTexture(dSortedPosTex));
	RX_CUCHECK(cudaUnbindTexture(dSortedVelTex));
#endif
}


/*!
 * �����Z���̃n�b�V�����v�Z
 * @param[in] 
 * @return 
 */
void CuCalcHashB(uint* dGridParticleHash, uint* dSortedIndex, float* dPos, 
				 float3 world_origin, float3 cell_width, uint3 grid_size, int nprts)
{
	uint numThreads, numBlocks;
	computeGridSize(nprts, THREAD_NUM, numBlocks, numThreads);

	// �J�[�l�����s
	calcHashB<<< numBlocks, numThreads >>>(dGridParticleHash,
										   dSortedIndex,
										   (float4*)dPos,
										   world_origin, 
										   cell_width, 
										   grid_size, 
										   nprts);
	
	RX_CUERROR("Kernel execution failed : calcHashB");	// �J�[�l���G���[�`�F�b�N
}

/*!
 * �p�[�e�B�N���z����\�[�g���ꂽ���Ԃɕ��ёւ��C
 * �e�Z���̎n�܂�ƏI���̃C���f�b�N�X������
 * @param[in] cell �p�[�e�B�N���O���b�h�f�[�^
 * @param[in] oldPos �p�[�e�B�N���ʒu
 */
void CuReorderDataAndFindCellStartB(rxParticleCell cell, float* oldPos)
{
	uint numThreads, numBlocks;
	computeGridSize(cell.uNumParticles, THREAD_NUM, numBlocks, numThreads);

	RX_CUCHECK(cudaMemset(cell.dCellStart, 0xffffffff, cell.uNumCells*sizeof(uint)));

	uint smemSize = sizeof(uint)*(numThreads+1);

	// �J�[�l�����s
	reorderDataAndFindCellStartB<<< numBlocks, numThreads, smemSize>>>(cell, (float4*)oldPos);

	RX_CUERROR("Kernel execution failed: CuReorderDataAndFindCellStartB");
	RX_CUCHECK(cudaThreadSynchronize());		// �S�ẴX���b�h���I���̂�҂�
}


//-----------------------------------------------------------------------------
// MARK:3D SPH
//-----------------------------------------------------------------------------



/*!
 * �p�[�e�B�N�����x�̌v�Z(�J�[�l���Ăяo��)
 * @param[out] dDens �p�[�e�B�N�����x
 * @param[out] dPres �p�[�e�B�N������
 * @param[in]  cell �p�[�e�B�N���O���b�h�f�[�^
 */
void CuSphDensity(float* dDens, float* dPres, rxParticleCell cell)
{
	// MRK:CuSphDensity2D
#if USE_TEX
	RX_CUCHECK(cudaBindTexture(0, dSortedPosTex, cell.dSortedPos, cell.uNumParticles*sizeof(float4)));
	RX_CUCHECK(cudaBindTexture(0, dCellStartTex, cell.dCellStart, cell.uNumCells*sizeof(uint)));
	RX_CUCHECK(cudaBindTexture(0, dCellEndTex, cell.dCellEnd, cell.uNumCells*sizeof(uint)));	
#endif
	//RX_CUCHECK(cudaMemset((void*)dNewDens, 0, sizeof(float2)*cell.uNumParticles));

	// 1�X���b�h/�p�[�e�B�N��
	uint numThreads, numBlocks;
	computeGridSize(cell.uNumParticles, THREAD_NUM, numBlocks, numThreads);

	// �J�[�l�����s
	sphCalDensity<<< numBlocks, numThreads >>>(dDens, dPres, cell);

	RX_CUERROR("sphCalDensity kernel execution failed");	// �J�[�l�����s�G���[�`�F�b�N
	RX_CUCHECK(cudaThreadSynchronize());		// �S�ẴX���b�h���I���̂�҂�

#if USE_TEX
	RX_CUCHECK(cudaUnbindTexture(dSortedPosTex));
	RX_CUCHECK(cudaUnbindTexture(dCellStartTex));
	RX_CUCHECK(cudaUnbindTexture(dCellEndTex));
#endif
}

/*!
 * �p�[�e�B�N���@���̌v�Z
 * @param[out] dNewDens �p�[�e�B�N�����x
 * @param[out] dNewPres �p�[�e�B�N������
 * @param[in]  cell �p�[�e�B�N���O���b�h�f�[�^
 */
void CuSphNormal(float* dNrms, float* dDens, rxParticleCell cell)
{
	// MRK:CuSphNormal

	// 1�X���b�h/�p�[�e�B�N��
	uint numThreads, numBlocks;
	computeGridSize(cell.uNumParticles, THREAD_NUM, numBlocks, numThreads);

	// �J�[�l�����s
	sphCalNormal<<< numBlocks, numThreads >>>((float4*)dNrms, dDens, cell);

	RX_CUERROR("sphCalNormal kernel execution failed");	// �J�[�l�����s�G���[�`�F�b�N
	RX_CUCHECK(cudaThreadSynchronize());		// �S�ẴX���b�h���I���̂�҂�

}

/*!
 * �p�[�e�B�N���ɂ�����͂̌v�Z(�J�[�l���Ăяo��)
 * @param[in] dDens �p�[�e�B�N�����x
 * @param[in] dPres �p�[�e�B�N������
 * @param[out] dFrc �p�[�e�B�N���ɂ������
 * @param[in] cell �p�[�e�B�N���O���b�h�f�[�^
 * @param[in] dt ���ԃX�e�b�v��
 */
void CuSphForces(float* dDens, float* dPres, float* dFrc, rxParticleCell cell, float dt)
{
#if USE_TEX
	RX_CUCHECK(cudaBindTexture(0, dSortedPosTex, cell.dSortedPos, cell.uNumParticles*sizeof(float4)));
	RX_CUCHECK(cudaBindTexture(0, dSortedVelTex, cell.dSortedVel, cell.uNumParticles*sizeof(float4)));
	RX_CUCHECK(cudaBindTexture(0, dCellStartTex, cell.dCellStart, cell.uNumCells*sizeof(uint)));
	RX_CUCHECK(cudaBindTexture(0, dCellEndTex, cell.dCellEnd, cell.uNumCells*sizeof(uint)));	
#endif

	// 1�X���b�h/�p�[�e�B�N��
	uint numThreads, numBlocks;
	computeGridSize(cell.uNumParticles, THREAD_NUM, numBlocks, numThreads);

	// �J�[�l�����s
	sphCalForces<<< numBlocks, numThreads >>>(dDens, dPres, (float4*)dFrc, cell);

	RX_CUERROR("calForcesSPH kernel execution failed");	// �J�[�l�����s�G���[�`�F�b�N
	RX_CUCHECK(cudaThreadSynchronize());		// �S�ẴX���b�h���I���̂�҂�

#if USE_TEX
	RX_CUCHECK(cudaUnbindTexture(dSortedPosTex));
	RX_CUCHECK(cudaUnbindTexture(dSortedVelTex));
	RX_CUCHECK(cudaUnbindTexture(dCellStartTex));
	RX_CUCHECK(cudaUnbindTexture(dCellEndTex));
#endif
}

/*!
 * �p�[�e�B�N���ʒu�C���x�̍X�V
 * @param[inout] pos �p�[�e�B�N���ʒu
 * @param[inout] vel �p�[�e�B�N�����x
 * @param[inout] velOld �O�X�e�b�v�̃p�[�e�B�N�����x
 * @param[in] frc �p�[�e�B�N���ɂ������
 * @param[in] dens �p�[�e�B�N�����x
 * @param[in] dt ���ԃX�e�b�v��
 * @param[in] nprts �p�[�e�B�N����
 */
void CuSphIntegrate(float* pos, float* vel, float* frc, float* dens, int* attr, 
					float dt, uint nprts)
{
	// 1�X���b�h/�p�[�e�B�N��
	uint numThreads, numBlocks;
	computeGridSize(nprts, THREAD_NUM, numBlocks, numThreads);

	// �J�[�l�����s
	sphIntegrate<<< numBlocks, numThreads >>>((float4*)pos, (float4*)vel, (float4*)frc, dens, attr, 
											  dt, nprts);
	
	RX_CUERROR("sphIntegrate kernel execution failed");	// �J�[�l�����s�G���[�`�F�b�N
	RX_CUCHECK(cudaThreadSynchronize());		// �S�ẴX���b�h���I���̂�҂�
}


/*!
 * �p�[�e�B�N���ʒu�C���x�̍X�V
 * @param[inout] pos �p�[�e�B�N���ʒu
 * @param[inout] vel �p�[�e�B�N�����x
 * @param[inout] velOld �O�X�e�b�v�̃p�[�e�B�N�����x
 * @param[in] frc �p�[�e�B�N���ɂ������
 * @param[in] dens �p�[�e�B�N�����x
 * @param[in] dt ���ԃX�e�b�v��
 * @param[in] nprts �p�[�e�B�N����
 */
void CuSphIntegrateWithPolygon(float* pos, float* vel, float* frc, float* dens, int* attr, 
							   float* vrts, int* tris, int tri_num, float dt, rxParticleCell cell)
{
	// 1�X���b�h/�p�[�e�B�N��
	uint numThreads, numBlocks;
	computeGridSize(cell.uNumParticles, THREAD_NUM, numBlocks, numThreads);

	// �J�[�l�����s
	sphIntegrateWithPolygon<<< numBlocks, numThreads >>>((float4*)pos, (float4*)vel, (float4*)frc, dens, attr, 
											   (float3*)vrts, (int3*)tris, tri_num, dt, cell);
	
	RX_CUERROR("sphIntegrate kernel execution failed");	// �J�[�l�����s�G���[�`�F�b�N
	RX_CUCHECK(cudaThreadSynchronize());		// �S�ẴX���b�h���I���̂�҂�
}


/*!
 * �p�[�e�B�N���ʒu���`�F�b�N���č폜�̈���Ȃ�Α�����-1�ɂ��āC�͈͊O�ɔ�΂�
 * @param[inout] pos �p�[�e�B�N���ʒu
 * @param[inout] vel �p�[�e�B�N�����x
 * @param[in] dens �p�[�e�B�N�����x
 * @param[in] dt ���ԃX�e�b�v��
 * @param[in] nprts �p�[�e�B�N����
 */
void CuSphCheckDelete(float* pos, float* vel, int* attr, float minp[3], float maxp[3], float farpoint[3], uint nprts)
{
	// 1�X���b�h/�p�[�e�B�N��
	uint numThreads, numBlocks;
	computeGridSize(nprts, THREAD_NUM, numBlocks, numThreads);

	float3 minp3, maxp3, farpoint3;
	minp3.x = minp[0];
	minp3.y = minp[1];
	minp3.z = minp[2];
	maxp3.x = maxp[0];
	maxp3.y = maxp[1];
	maxp3.z = maxp[2];
	farpoint3.x = farpoint[0];
	farpoint3.y = farpoint[1];
	farpoint3.z = farpoint[2];


	// �J�[�l�����s
	checkDelete<<< numBlocks, numThreads >>>((float4*)pos, (float4*)vel, attr, minp3, maxp3, farpoint3, nprts);
	
	RX_CUERROR("sphIntegrate kernel execution failed");	// �J�[�l�����s�G���[�`�F�b�N
	RX_CUCHECK(cudaThreadSynchronize());		// �S�ẴX���b�h���I���̂�҂�
}

/*!
 * �p�[�e�B�N����x���W�l���`�F�b�N���Ă���l�ȏ�Ȃ�Α�����-1�ɂ��āC�͈͊O�ɔ�΂�
 * @param[inout] pos �p�[�e�B�N���ʒu
 * @param[inout] vel �p�[�e�B�N�����x
 * @param[in] dens �p�[�e�B�N�����x
 * @param[in] dt ���ԃX�e�b�v��
 * @param[in] nprts �p�[�e�B�N����
 */
void CuSphCheckDeleteX(float* pos, float* vel, int* attr, float xmax, float farpoint[3], uint nprts)
{
	// 1�X���b�h/�p�[�e�B�N��
	uint numThreads, numBlocks;
	computeGridSize(nprts, THREAD_NUM, numBlocks, numThreads);

	float3 farpoint3;
	farpoint3.x = farpoint[0];
	farpoint3.y = farpoint[1];
	farpoint3.z = farpoint[2];

	// �J�[�l�����s
	checkDeleteX<<< numBlocks, numThreads >>>((float4*)pos, (float4*)vel, attr, xmax, farpoint3, nprts);
	
	RX_CUERROR("sphIntegrate kernel execution failed");	// �J�[�l�����s�G���[�`�F�b�N
	RX_CUCHECK(cudaThreadSynchronize());		// �S�ẴX���b�h���I���̂�҂�
}

/*!
 * �O���b�h��̖��x���Z�o
 * @param[out] dGridD �O���b�h��̖��x�l
 * @param[in] cell �p�[�e�B�N���O���b�h�f�[�^
 * @param[in] nx,ny �O���b�h��
 * @param[in] x0,y0 �O���b�h�ŏ����W
 * @param[in] dx,dy �O���b�h��
 */
void CuSphGridDensity(float *dGridD, rxParticleCell cell, 
					  int nx, int ny, int nz, float x0, float y0, float z0, float dx, float dy, float dz)
{
#if USE_TEX
	RX_CUCHECK(cudaBindTexture(0, dSortedPosTex, cell.dSortedPos, cell.uNumParticles*sizeof(float4)));
	RX_CUCHECK(cudaBindTexture(0, dCellStartTex, cell.dCellStart, cell.uNumCells*sizeof(uint)));
	RX_CUCHECK(cudaBindTexture(0, dCellEndTex, cell.dCellEnd, cell.uNumCells*sizeof(uint)));	
#endif

	uint3  gnum = make_uint3(nx, ny, nz);
	float3 gmin = make_float3(x0, y0, z0);
	float3 glen = make_float3(dx, dy, dz);

	int numcell = gnum.x*gnum.y*gnum.z;

	int threads = 128;
	dim3 grid((numcell+threads-1)/threads, 1, 1);
	if(grid.x > 65535){
		grid.y = (grid.x+32768-1)/32768;
		grid.x = 32768;
	}

	// �J�[�l�����s
	sphCalDensityInGrid<<<grid, threads>>>(dGridD, cell, gnum, gmin, glen);

	RX_CUERROR("Kernel execution failed");	// �J�[�l�����s�G���[�`�F�b�N
	RX_CUCHECK(cudaThreadSynchronize());		// �S�ẴX���b�h���I���̂�҂�

#if USE_TEX
	RX_CUCHECK(cudaUnbindTexture(dSortedPosTex));
	RX_CUCHECK(cudaUnbindTexture(dCellStartTex));
	RX_CUCHECK(cudaUnbindTexture(dCellEndTex));
#endif
}

/*!�ǉ��@�X�p
 * �O���b�h��̖��x���Z�o
 * @param[out] dGridD �O���b�h��̖��x�l
 * @param[in] cell �p�[�e�B�N���O���b�h�f�[�^
 * @param[in] nx,ny �O���b�h��
 * @param[in] x0,y0 �O���b�h�ŏ����W
 * @param[in] dx,dy �O���b�h��
 */
void CuIceMeshMake(float *dGridD, rxParticleCell cell, 
					  int nx, int ny, int nz, float x0, float y0, float z0, float dx, float dy, float dz, float *bIceFlag)
{
#if USE_TEX
	RX_CUCHECK(cudaBindTexture(0, dSortedPosTex, cell.dSortedPos, cell.uNumParticles*sizeof(float4)));
	RX_CUCHECK(cudaBindTexture(0, dCellStartTex, cell.dCellStart, cell.uNumCells*sizeof(uint)));
	RX_CUCHECK(cudaBindTexture(0, dCellEndTex, cell.dCellEnd, cell.uNumCells*sizeof(uint)));	
#endif

	uint3  gnum = make_uint3(nx, ny, nz);
	float3 gmin = make_float3(x0, y0, z0);
	float3 glen = make_float3(dx, dy, dz);

	int numcell = gnum.x*gnum.y*gnum.z;

	int threads = 128;
	dim3 grid((numcell+threads-1)/threads, 1, 1);
	if(grid.x > 65535){
		grid.y = (grid.x+32768-1)/32768;
		grid.x = 32768;
	}

	// �J�[�l�����s
	sphCalDensityInGridIceMesh<<<grid, threads>>>(dGridD, cell, gnum, gmin, glen, bIceFlag);

	RX_CUERROR("Kernel execution failed");	// �J�[�l�����s�G���[�`�F�b�N
	RX_CUCHECK(cudaThreadSynchronize());		// �S�ẴX���b�h���I���̂�҂�

#if USE_TEX
	RX_CUCHECK(cudaUnbindTexture(dSortedPosTex));
	RX_CUCHECK(cudaUnbindTexture(dCellStartTex));
	RX_CUCHECK(cudaUnbindTexture(dCellEndTex));
#endif
}


//�ǉ��F�\�ʃp�[�e�B�N���̌��o
#define NEIGHT_MAX 30	//�ߖT�ő嗱�q��
void CuSphDetectSurfaceParticles(int* neights, int* surface, rxParticleCell celldata, float* pos, float radius)
{
	uint numThreads, numBlocks;
	computeGridSize(celldata.uNumParticles, THREAD_NUM, numBlocks, numThreads);

	detectSurfaceParticles<<<numBlocks, numThreads>>>(neights, surface, celldata, NEIGHT_MAX, (float4*)pos, radius);

	RX_CUERROR("Kernel execution failed");	// �J�[�l�����s�G���[�`�F�b�N
	cudaThreadSynchronize();
}

//�ǉ��F�ߖT�p�[�e�B�N���̌��o
void CuSphDetectNeighborParticles(int* neights, rxParticleCell celldata, float radius, int prtNum)
{
	uint numThreads, numBlocks;
	computeGridSize(celldata.uNumParticles, THREAD_NUM, numBlocks, numThreads);

	detectNeighborParticles<<<numBlocks, numThreads>>>(neights, celldata, radius, NEIGHT_MAX);

	RX_CUERROR("Kernel execution failed");	// �J�[�l�����s�G���[�`�F�b�N
	cudaThreadSynchronize();
}

/*!
 * ���E�p�[�e�B�N���̑̐ς��v�Z
 *  - "Versatile Rigid-Fluid Coupling for Incompressible SPH", 2.2 ��(3)�̏�
 * @param[out] dVolB ���E�p�[�e�B�N���̑̐�
 * @param[in]  mass �p�[�e�B�N������
 * @param[in]  cell �p�[�e�B�N���O���b�h�f�[�^
 */
void CuSphBoundaryVolume(float* dVolB, float mass, rxParticleCell cell)
{
	// 1�X���b�h/�p�[�e�B�N��
	uint numThreads, numBlocks;
	computeGridSize(cell.uNumParticles, THREAD_NUM, numBlocks, numThreads);

	// �J�[�l�����s
	sphCalBoundaryVolume<<< numBlocks, numThreads >>>(dVolB, cell);

	RX_CUERROR("kernel execution failed : sphCalBoundaryVolume");	// �J�[�l�����s�G���[�`�F�b�N
	RX_CUCHECK(cudaThreadSynchronize());		// �S�ẴX���b�h���I���̂�҂�
}

/*!
 * �p�[�e�B�N�����x�̌v�Z(�J�[�l���Ăяo��)
 * @param[out] dDens �p�[�e�B�N�����x
 * @param[out] dPres �p�[�e�B�N������
 * @param[in]  cell �p�[�e�B�N���O���b�h�f�[�^
 */
void CuSphBoundaryDensity(float* dDens, float* dPres, float* dPos, float* dVolB, rxParticleCell bcell, uint pnum)
{
	// 1�X���b�h/�p�[�e�B�N��
	uint numThreads, numBlocks;
	computeGridSize(pnum, THREAD_NUM, numBlocks, numThreads);

	// �J�[�l�����s
	sphCalBoundaryDensity<<< numBlocks, numThreads >>>(dDens, dPres, (float4*)dPos, dVolB, bcell, pnum);

	RX_CUERROR("kernel execution failed : sphCalBoundaryDensity");	// �J�[�l�����s�G���[�`�F�b�N
	RX_CUCHECK(cudaThreadSynchronize());		// �S�ẴX���b�h���I���̂�҂�
}

/*!
 * �p�[�e�B�N���ɂ�����͂̌v�Z(�J�[�l���Ăяo��)
 * @param[in] dDens �p�[�e�B�N�����x
 * @param[in] dPres �p�[�e�B�N������
 * @param[out] dFrc �p�[�e�B�N���ɂ������
 * @param[in] cell �p�[�e�B�N���O���b�h�f�[�^
 * @param[in] dt ���ԃX�e�b�v��
 */
void CuSphBoundaryForces(float* dDens, float* dPres, float* dPos, float* dVolB, float* dFrc, rxParticleCell bcell, uint pnum)
{
	// 1�X���b�h/�p�[�e�B�N��
	uint numThreads, numBlocks;
	computeGridSize(pnum, THREAD_NUM, numBlocks, numThreads);

	// �J�[�l�����s
	sphCalBoundaryForce<<< numBlocks, numThreads >>>(dDens, dPres, (float4*)dPos, dVolB, (float4*)dFrc, bcell, pnum);

	RX_CUERROR("kernel execution failed : sphCalBoundaryForce");	// �J�[�l�����s�G���[�`�F�b�N
	RX_CUCHECK(cudaThreadSynchronize());		// �S�ẴX���b�h���I���̂�҂�
}

//-----------------------------------------------------------------------------
// PBDSPH
//-----------------------------------------------------------------------------

/*!
 * �p�[�e�B�N�����x�̌v�Z(�J�[�l���Ăяo��)
 * @param[out] dDens �p�[�e�B�N�����x
 * @param[out] dPres �p�[�e�B�N������
 * @param[in]  cell �p�[�e�B�N���O���b�h�f�[�^
 */
void CuPbdSphDensity(float* dDens, rxParticleCell cell)
{
#if USE_TEX
	RX_CUCHECK(cudaBindTexture(0, dSortedPosTex, cell.dSortedPos, cell.uNumParticles*sizeof(float4)));
	RX_CUCHECK(cudaBindTexture(0, dCellStartTex, cell.dCellStart, cell.uNumCells*sizeof(uint)));
	RX_CUCHECK(cudaBindTexture(0, dCellEndTex, cell.dCellEnd, cell.uNumCells*sizeof(uint)));	
#endif
	//RX_CUCHECK(cudaMemset((void*)dNewDens, 0, sizeof(float2)*cell.uNumParticles));

	// 1�X���b�h/�p�[�e�B�N��
	uint numThreads, numBlocks;
	computeGridSize(cell.uNumParticles, THREAD_NUM, numBlocks, numThreads);

	// �J�[�l�����s
	pbdsphCalDensity<<< numBlocks, numThreads >>>(dDens, cell);

	RX_CUERROR("pbdsphCalDensity kernel execution failed");	// �J�[�l�����s�G���[�`�F�b�N
	RX_CUCHECK(cudaThreadSynchronize());		// �S�ẴX���b�h���I���̂�҂�

#if USE_TEX
	RX_CUCHECK(cudaUnbindTexture(dSortedPosTex));
	RX_CUCHECK(cudaUnbindTexture(dCellStartTex));
	RX_CUCHECK(cudaUnbindTexture(dCellEndTex));
#endif
}

/*!
 * �p�[�e�B�N���ɂ�����͂̌v�Z(�J�[�l���Ăяo��)
 * @param[in] dDens �p�[�e�B�N�����x
 * @param[out] dFrc �p�[�e�B�N���ɂ������
 * @param[in] cell �p�[�e�B�N���O���b�h�f�[�^
 * @param[in] dt ���ԃX�e�b�v��
 */
void CuPbdSphExternalForces(float* dDens, float* dFrc, rxParticleCell cell, float dt)
{
#if USE_TEX
	RX_CUCHECK(cudaBindTexture(0, dSortedPosTex, cell.dSortedPos, cell.uNumParticles*sizeof(float4)));
	RX_CUCHECK(cudaBindTexture(0, dSortedVelTex, cell.dSortedVel, cell.uNumParticles*sizeof(float4)));
	RX_CUCHECK(cudaBindTexture(0, dCellStartTex, cell.dCellStart, cell.uNumCells*sizeof(uint)));
	RX_CUCHECK(cudaBindTexture(0, dCellEndTex, cell.dCellEnd, cell.uNumCells*sizeof(uint)));	
#endif

	// 1�X���b�h/�p�[�e�B�N��
	uint numThreads, numBlocks;
	computeGridSize(cell.uNumParticles, THREAD_NUM, numBlocks, numThreads);

	// �J�[�l�����s
	pbdsphCalExternalForces<<< numBlocks, numThreads >>>(dDens, (float4*)dFrc, cell);

	RX_CUERROR("calForcesSPH kernel execution failed");	// �J�[�l�����s�G���[�`�F�b�N
	RX_CUCHECK(cudaThreadSynchronize());		// �S�ẴX���b�h���I���̂�҂�

#if USE_TEX
	RX_CUCHECK(cudaUnbindTexture(dSortedPosTex));
	RX_CUCHECK(cudaUnbindTexture(dSortedVelTex));
	RX_CUCHECK(cudaUnbindTexture(dCellStartTex));
	RX_CUCHECK(cudaUnbindTexture(dCellEndTex));
#endif
}

/*!
 * �X�P�[�����O�t�@�N�^�̌v�Z
 * @param[in] dPos �p�[�e�B�N�����S���W
 * @param[out] dDens �p�[�e�B�N�����x
 * @param[out] dScl �X�P�[�����O�t�@�N�^
 * @param[in] eps �ɘa�W��
 * @param[in] cell �p�[�e�B�N���O���b�h�f�[�^
 */
void CuPbdSphScalingFactor(float* dPos, float* dDens, float* dScl, float eps, rxParticleCell cell)
{
#if USE_TEX
	RX_CUCHECK(cudaBindTexture(0, dSortedPosTex, cell.dSortedPos, cell.uNumParticles*sizeof(float4)));
	RX_CUCHECK(cudaBindTexture(0, dCellStartTex, cell.dCellStart, cell.uNumCells*sizeof(uint)));
	RX_CUCHECK(cudaBindTexture(0, dCellEndTex, cell.dCellEnd, cell.uNumCells*sizeof(uint)));	
#endif
	//RX_CUCHECK(cudaMemset((void*)dNewDens, 0, sizeof(float2)*cell.uNumParticles));

	// 1�X���b�h/�p�[�e�B�N��
	uint numThreads, numBlocks;
	computeGridSize(cell.uNumParticles, THREAD_NUM, numBlocks, numThreads);

	// �J�[�l�����s
	pbdsphCalScalingFactor<<< numBlocks, numThreads >>>((float4*)dPos, dDens, dScl, eps, cell);

	RX_CUERROR("pbdsphCalScalingFactor kernel execution failed");	// �J�[�l�����s�G���[�`�F�b�N
	RX_CUCHECK(cudaThreadSynchronize());		// �S�ẴX���b�h���I���̂�҂�

#if USE_TEX
	RX_CUCHECK(cudaUnbindTexture(dSortedPosTex));
	RX_CUCHECK(cudaUnbindTexture(dCellStartTex));
	RX_CUCHECK(cudaUnbindTexture(dCellEndTex));
#endif
}

/*!
 * ���ϖ��x�ϓ��̌v�Z
 *  - ���ׂẴp�[�e�B�N�����x�̏������x�Ƃ̍����J�[�l���Ōv�Z���CPrefix Sum (Scan)�ł��̍��v�����߂�
 * @param[out] dErrScan �ϓ��l��Scan���ʂ��i�[����z��
 * @param[out] dErr �p�[�e�B�N�����x�ϓ��l
 * @param[in] dDens �p�[�e�B�N�����x
 * @param[in] rest_dens �������x
 * @param[in] nprts �p�[�e�B�N����
 * @return ���ϖ��x�ϓ�
 */
float CuPbdSphCalDensityFluctuation(float* dErrScan, float* dErr, float* dDens, float rest_dens, uint nprts)
{
	// 1�X���b�h/�p�[�e�B�N��
	uint numThreads, numBlocks;
	computeGridSize(nprts, THREAD_NUM, numBlocks, numThreads);

	// �J�[�l�����s
	pbdsphDensityFluctuation<<< numBlocks, numThreads >>>(dErr, dDens, rest_dens, nprts);

	RX_CUERROR("pbdsphDensityFluctuation kernel execution failed");	// �J�[�l�����s�G���[�`�F�b�N
	RX_CUCHECK(cudaThreadSynchronize());		// �S�ẴX���b�h���I���̂�҂�

	// �e�p�[�e�B�N���̖��x�ϓ���Scan
	CuScanf(dErrScan, dErr, nprts);

	// Exclusive scan (�Ō�̗v�f��0�Ԗڂ���n-2�Ԗڂ܂ł̍��v�ɂȂ��Ă���)�Ȃ̂ŁC
	// Scan�O�z��̍Ō�(n-1�Ԗ�)�̗v�f�ƍ��v���邱�ƂŖ��x�ϓ��̍��v���v�Z
	float lval, lsval;
	RX_CUCHECK(cudaMemcpy((void*)&lval, (void*)(dErr+nprts-1), sizeof(float), cudaMemcpyDeviceToHost));
	RX_CUCHECK(cudaMemcpy((void*)&lsval, (void*)(dErrScan+nprts-1), sizeof(float), cudaMemcpyDeviceToHost));
	float dens_var = lval+lsval;

	return dens_var/(float)nprts;
}

/*!
 * �ʒu�C���ʂ̌v�Z
 * @param[in] dPos �p�[�e�B�N�����S���W
 * @param[in] dScl �X�P�[�����O�t�@�N�^
 * @param[out] dDp �ʒu�C����
 * @param[in] cell �p�[�e�B�N���O���b�h�f�[�^
 */
void CuPbdSphPositionCorrection(float* dPos, float* dScl, float* dDp, rxParticleCell cell)
{
#if USE_TEX
	RX_CUCHECK(cudaBindTexture(0, dSortedPosTex, cell.dSortedPos, cell.uNumParticles*sizeof(float4)));
	RX_CUCHECK(cudaBindTexture(0, dCellStartTex, cell.dCellStart, cell.uNumCells*sizeof(uint)));
	RX_CUCHECK(cudaBindTexture(0, dCellEndTex, cell.dCellEnd, cell.uNumCells*sizeof(uint)));	
#endif
	//RX_CUCHECK(cudaMemset((void*)dNewDens, 0, sizeof(float2)*cell.uNumParticles));

	// 1�X���b�h/�p�[�e�B�N��
	uint numThreads, numBlocks;
	computeGridSize(cell.uNumParticles, THREAD_NUM, numBlocks, numThreads);

	// �J�[�l�����s
	pbdsphPositionCorrection<<< numBlocks, numThreads >>>((float4*)dPos, dScl, (float4*)dDp, cell);

	RX_CUERROR("pbdsphPositionCorrection kernel execution failed");	// �J�[�l�����s�G���[�`�F�b�N
	RX_CUCHECK(cudaThreadSynchronize());		// �S�ẴX���b�h���I���̂�҂�

#if USE_TEX
	RX_CUCHECK(cudaUnbindTexture(dSortedPosTex));
	RX_CUCHECK(cudaUnbindTexture(dCellStartTex));
	RX_CUCHECK(cudaUnbindTexture(dCellEndTex));
#endif
}

/*!
 * �p�[�e�B�N���ʒu���X�V
 * @param[inout] dPos �p�[�e�B�N���ʒu
 * @param[in] dDp �ʒu�C����
 * @param[in] nprts �p�[�e�B�N����
 */
void CuPbdSphCorrectPosition(float* dPos, float* dDp, uint nprts)
{
	// 1�X���b�h/�p�[�e�B�N��
	uint numThreads, numBlocks;
	computeGridSize(nprts, THREAD_NUM, numBlocks, numThreads);

	// �J�[�l�����s
	pbdsphCorrectPosition<<< numBlocks, numThreads >>>((float4*)dPos, (float4*)dDp, nprts);
	
	RX_CUERROR("pbdsphCorrectPosition kernel execution failed");	// �J�[�l�����s�G���[�`�F�b�N
	RX_CUCHECK(cudaThreadSynchronize());		// �S�ẴX���b�h���I���̂�҂�
}



/*!
 * ���E�p�[�e�B�N�����x���]���̃p�[�e�B�N�����x�ɉ�����
 * @param[inout] dDens ���̃p�[�e�B�N�����x
 * @param[in] dPos  ���̃p�[�e�B�N������
 * @param[in] dVolB ���E�p�[�e�B�N���̐�
 * @param[in] bcell ���E�p�[�e�B�N���O���b�h�f�[�^
 */
void CuPbdSphBoundaryDensity(float* dDens, float* dPos, float* dVolB, rxParticleCell bcell, uint pnum)
{
	// 1�X���b�h/�p�[�e�B�N��
	uint numThreads, numBlocks;
	computeGridSize(pnum, THREAD_NUM, numBlocks, numThreads);

	// �J�[�l�����s
	pbdsphCalBoundaryDensity<<< numBlocks, numThreads >>>(dDens, (float4*)dPos, dVolB, bcell, pnum);

	RX_CUERROR("kernel execution failed : sphCalBoundaryDensity");	// �J�[�l�����s�G���[�`�F�b�N
	RX_CUCHECK(cudaThreadSynchronize());		// �S�ẴX���b�h���I���̂�҂�
}


/*!
 * �X�P�[�����O�t�@�N�^�̌v�Z(���E�p�[�e�B�N���܂�)
 * @param[in] dPos ���̃p�[�e�B�N�����S���W
 * @param[out] dDens ���̃p�[�e�B�N�����x
 * @param[out] dScl ���̃p�[�e�B�N���̃X�P�[�����O�t�@�N�^
 * @param[in] eps �ɘa�W��
 * @param[in] cell ���̃p�[�e�B�N���O���b�h�f�[�^
 * @param[in] dVolB ���E�p�[�e�B�N���̐�
 * @param[out] dSclB ���E�p�[�e�B�N���̃X�P�[�����O�t�@�N�^
 * @param[in] bcell ���E�p�[�e�B�N���O���b�h�f�[�^
 */
void CuPbdSphScalingFactorWithBoundary(float* dPos, float* dDens, float* dScl, float eps, rxParticleCell cell, 
									   float* dVolB, float* dSclB, rxParticleCell bcell)
{
#if USE_TEX
	RX_CUCHECK(cudaBindTexture(0, dSortedPosTex, cell.dSortedPos, cell.uNumParticles*sizeof(float4)));
	RX_CUCHECK(cudaBindTexture(0, dCellStartTex, cell.dCellStart, cell.uNumCells*sizeof(uint)));
	RX_CUCHECK(cudaBindTexture(0, dCellEndTex, cell.dCellEnd, cell.uNumCells*sizeof(uint)));	
#endif
	//RX_CUCHECK(cudaMemset((void*)dNewDens, 0, sizeof(float2)*cell.uNumParticles));

	// ���̃p�[�e�B�N���̐������X���b�h�𗧂Ă�
	uint numThreads, numBlocks;
	computeGridSize(cell.uNumParticles, THREAD_NUM, numBlocks, numThreads);

	// �J�[�l�����s
	pbdsphCalScalingFactorWithBoundary<<< numBlocks, numThreads >>>((float4*)dPos, dDens, dScl, eps, cell, dVolB, bcell);

	RX_CUERROR("kernel execution failed : pbdsphCalScalingFactorWithBoundary");	// �J�[�l�����s�G���[�`�F�b�N
	RX_CUCHECK(cudaThreadSynchronize());		// �S�ẴX���b�h���I���̂�҂�

#if USE_TEX
	RX_CUCHECK(cudaUnbindTexture(dSortedPosTex));
	RX_CUCHECK(cudaUnbindTexture(dCellStartTex));
	RX_CUCHECK(cudaUnbindTexture(dCellEndTex));
#endif

	// ���E�p�[�e�B�N���̃X�P�[�����O�t�@�N�^�̌v�Z
	// ���E�p�[�e�B�N���̐������X���b�h�𗧂Ă�
	computeGridSize(bcell.uNumParticles, THREAD_NUM, numBlocks, numThreads);

	// �J�[�l�����s
	pbdsphCalBoundaryScalingFactor<<< numBlocks, numThreads >>>((float4*)dPos, dDens, eps, cell, dVolB, dSclB, bcell);

	RX_CUERROR("kernel execution failed : pbdsphCalScalingFactorWithBoundary");	// �J�[�l�����s�G���[�`�F�b�N
	RX_CUCHECK(cudaThreadSynchronize());		// �S�ẴX���b�h���I���̂�҂�

}

/*!
 * �ʒu�C���ʂ̌v�Z(���E�p�[�e�B�N���܂�)
 * @param[in] dPos ���̃p�[�e�B�N�����S���W
 * @param[in] dScl ���̃p�[�e�B�N���̃X�P�[�����O�t�@�N�^
 * @param[out] dDens ���̃p�[�e�B�N���ʒu�C����
 * @param[in] cell ���̃p�[�e�B�N���O���b�h�f�[�^
 * @param[in] dVolB ���E�p�[�e�B�N���̐�
 * @param[in] dSclB ���E�p�[�e�B�N���̃X�P�[�����O�t�@�N�^
 * @param[in] bcell ���E�p�[�e�B�N���O���b�h�f�[�^
 */
void CuPbdSphPositionCorrectionWithBoundary(float* dPos, float* dScl, float* dDp, rxParticleCell cell, 
											float* dVolB, float* dSclB, rxParticleCell bcell)
{
#if USE_TEX
	RX_CUCHECK(cudaBindTexture(0, dSortedPosTex, cell.dSortedPos, cell.uNumParticles*sizeof(float4)));
	RX_CUCHECK(cudaBindTexture(0, dCellStartTex, cell.dCellStart, cell.uNumCells*sizeof(uint)));
	RX_CUCHECK(cudaBindTexture(0, dCellEndTex, cell.dCellEnd, cell.uNumCells*sizeof(uint)));	
#endif
	//RX_CUCHECK(cudaMemset((void*)dNewDens, 0, sizeof(float2)*cell.uNumParticles));

	// 1�X���b�h/�p�[�e�B�N��
	uint numThreads, numBlocks;
	computeGridSize(cell.uNumParticles, THREAD_NUM, numBlocks, numThreads);

	// �J�[�l�����s
	pbdsphPositionCorrectionWithBoundary<<< numBlocks, numThreads >>>((float4*)dPos, dScl, (float4*)dDp, cell, 
																	  dVolB, dSclB, bcell);

	RX_CUERROR("kernel execution failed");	// �J�[�l�����s�G���[�`�F�b�N
	RX_CUCHECK(cudaThreadSynchronize());		// �S�ẴX���b�h���I���̂�҂�

#if USE_TEX
	RX_CUCHECK(cudaUnbindTexture(dSortedPosTex));
	RX_CUCHECK(cudaUnbindTexture(dCellStartTex));
	RX_CUCHECK(cudaUnbindTexture(dCellEndTex));
#endif
}



/*!
 * �p�[�e�B�N���ʒu�C���x�̍X�V
 * @param[inout] pos �p�[�e�B�N���ʒu
 * @param[inout] vel �p�[�e�B�N�����x
 * @param[inout] velOld �O�X�e�b�v�̃p�[�e�B�N�����x
 * @param[in] frc �p�[�e�B�N���ɂ������
 * @param[in] dens �p�[�e�B�N�����x
 * @param[in] dt ���ԃX�e�b�v��
 * @param[in] nprts �p�[�e�B�N����
 */
void CuPbdSphIntegrate(float* pos, float* vel, float* acc, int* attr, 
					   float* new_pos, float* new_vel, float dt, uint nprts)
{
	// 1�X���b�h/�p�[�e�B�N��
	uint numThreads, numBlocks;
	computeGridSize(nprts, THREAD_NUM, numBlocks, numThreads);

	// �J�[�l�����s
	pbdsphIntegrate<<< numBlocks, numThreads >>>((float4*)pos, (float4*)vel, (float4*)acc, attr, 
												 (float4*)new_pos, (float4*)new_vel, dt, nprts);
	
	RX_CUERROR("pbdsphIntegrate kernel execution failed");	// �J�[�l�����s�G���[�`�F�b�N
	RX_CUCHECK(cudaThreadSynchronize());		// �S�ẴX���b�h���I���̂�҂�
}


/*!
 * �p�[�e�B�N���ʒu�C���x�̍X�V
 * @param[inout] pos �p�[�e�B�N���ʒu
 * @param[inout] vel �p�[�e�B�N�����x
 * @param[inout] velOld �O�X�e�b�v�̃p�[�e�B�N�����x
 * @param[in] frc �p�[�e�B�N���ɂ������
 * @param[in] dens �p�[�e�B�N�����x
 * @param[in] dt ���ԃX�e�b�v��
 * @param[in] nprts �p�[�e�B�N����
 */
void CuPbdSphIntegrateWithPolygon(float* pos, float* vel, float* acc, int* attr, 
								  float* new_pos, float* new_vel, 
								  float* vrts, int* tris, int tri_num, float dt, rxParticleCell cell)
{
	// 1�X���b�h/�p�[�e�B�N��
	uint numThreads, numBlocks;
	computeGridSize(cell.uNumParticles, THREAD_NUM, numBlocks, numThreads);

	// �J�[�l�����s
	pbdsphIntegrateWithPolygon<<< numBlocks, numThreads >>>((float4*)pos, (float4*)vel, (float4*)acc, attr, 
															(float4*)new_pos, (float4*)new_vel, (float3*)vrts, (int3*)tris, tri_num, dt, cell);
	
	RX_CUERROR("pbdsphIntegrateWithPolygon kernel execution failed");	// �J�[�l�����s�G���[�`�F�b�N
	RX_CUCHECK(cudaThreadSynchronize());		// �S�ẴX���b�h���I���̂�҂�
}


/*!
 * �p�[�e�B�N���ʒu�C���x�̍X�V
 * @param[in] pos �X�V���ꂽ�p�[�e�B�N���ʒu
 * @param[inout] new_pos �X�e�b�v�ŏ��̃p�[�e�B�N���ʒu/�V�����p�[�e�B�N�����x
 * @param[out] new_vel �V�����p�[�e�B�N�����x
 * @param[in] dt ���ԃX�e�b�v��
 * @param[in] nprts �p�[�e�B�N����
 */
void CuPbdSphUpdatePosition(float* pos, float* new_pos, float* new_vel, float dt, uint nprts)
{
	// 1�X���b�h/�p�[�e�B�N��
	uint numThreads, numBlocks;
	computeGridSize(nprts, THREAD_NUM, numBlocks, numThreads);

	// �J�[�l�����s
	pbdsphUpdatePosition<<< numBlocks, numThreads >>>((float4*)pos, (float4*)new_pos, (float4*)new_vel, dt, nprts);
	
	RX_CUERROR("CuPbdSphUpdatePosition kernel execution failed");	// �J�[�l�����s�G���[�`�F�b�N
	RX_CUCHECK(cudaThreadSynchronize());		// �S�ẴX���b�h���I���̂�҂�
}

/*!
 * �p�[�e�B�N���ʒu�C���x�̍X�V
 * @param[in] pos �X�V���ꂽ�p�[�e�B�N���ʒu
 * @param[inout] new_pos �X�e�b�v�ŏ��̃p�[�e�B�N���ʒu/�V�����p�[�e�B�N�����x
 * @param[out] new_vel �V�����p�[�e�B�N�����x
 * @param[in] dt ���ԃX�e�b�v��
 * @param[in] nprts �p�[�e�B�N����
 */
void CuPbdSphUpdateVelocity(float* pos, float* new_pos, float* new_vel, float dt, uint nprts)
{
	// 1�X���b�h/�p�[�e�B�N��
	uint numThreads, numBlocks;
	computeGridSize(nprts, THREAD_NUM, numBlocks, numThreads);

	// �J�[�l�����s
	pbdsphUpdateVelocity<<< numBlocks, numThreads >>>((float4*)pos, (float4*)new_pos, (float4*)new_vel, dt, nprts);
	
	RX_CUERROR("pbdsphUpdateVelocity kernel execution failed");	// �J�[�l�����s�G���[�`�F�b�N
	RX_CUCHECK(cudaThreadSynchronize());		// �S�ẴX���b�h���I���̂�҂�
}

/*!
 * XSPH�ɂ��S���v�Z
 * @param[in] dPos �p�[�e�B�N�����S���W
 * @param[in] dVel �p�[�e�B�N�����x
 * @param[out] dNewVel �X�V���ꂽ�p�[�e�B�N�����x
 * @param[in] c �S���v�Z�p�p�����[�^
 * @param[in] cell �p�[�e�B�N���O���b�h�f�[�^
 */
void CuXSphViscosity(float* dPos, float* dVel, float* dNewVel, float* dDens, float c, rxParticleCell cell)
{
#if USE_TEX
	RX_CUCHECK(cudaBindTexture(0, dSortedPosTex, cell.dSortedPos, cell.uNumParticles*sizeof(float4)));
	RX_CUCHECK(cudaBindTexture(0, dCellStartTex, cell.dCellStart, cell.uNumCells*sizeof(uint)));
	RX_CUCHECK(cudaBindTexture(0, dCellEndTex, cell.dCellEnd, cell.uNumCells*sizeof(uint)));	
#endif

	// 1�X���b�h/�p�[�e�B�N��
	uint numThreads, numBlocks;
	computeGridSize(cell.uNumParticles, THREAD_NUM, numBlocks, numThreads);

	// �J�[�l�����s
	xsphVisocosity<<< numBlocks, numThreads >>>((float4*)dPos, (float4*)dVel, (float4*)dNewVel, dDens, c, cell);

	RX_CUERROR("pbdsphPositionCorrection kernel execution failed");	// �J�[�l�����s�G���[�`�F�b�N
	RX_CUCHECK(cudaThreadSynchronize());		// �S�ẴX���b�h���I���̂�҂�

#if USE_TEX
	RX_CUCHECK(cudaUnbindTexture(dSortedPosTex));
	RX_CUCHECK(cudaUnbindTexture(dCellStartTex));
	RX_CUCHECK(cudaUnbindTexture(dCellEndTex));
#endif
}



/*!
 * �p�[�e�B�N���ʒu���`�F�b�N���č폜�̈���Ȃ�Α�����-1�ɂ��āC�͈͊O�ɔ�΂�
 * @param[inout] pos �p�[�e�B�N���ʒu
 * @param[inout] vel �p�[�e�B�N�����x
 * @param[in] dens �p�[�e�B�N�����x
 * @param[in] dt ���ԃX�e�b�v��
 * @param[in] nprts �p�[�e�B�N����
 */
void CuPbdSphCheckDelete(float* pos, float* vel, int* attr, float minp[3], float maxp[3], float farpoint[3], uint nprts)
{
	// 1�X���b�h/�p�[�e�B�N��
	uint numThreads, numBlocks;
	computeGridSize(nprts, THREAD_NUM, numBlocks, numThreads);

	float3 minp3, maxp3, farpoint3;
	minp3.x = minp[0];
	minp3.y = minp[1];
	minp3.z = minp[2];
	maxp3.x = maxp[0];
	maxp3.y = maxp[1];
	maxp3.z = maxp[2];
	farpoint3.x = farpoint[0];
	farpoint3.y = farpoint[1];
	farpoint3.z = farpoint[2];


	// �J�[�l�����s
	checkDeletePB<<< numBlocks, numThreads >>>((float4*)pos, (float4*)vel, attr, minp3, maxp3, farpoint3, nprts);
	
	RX_CUERROR("sphIntegrate kernel execution failed");	// �J�[�l�����s�G���[�`�F�b�N
	RX_CUCHECK(cudaThreadSynchronize());		// �S�ẴX���b�h���I���̂�҂�
}

/*!
 * �p�[�e�B�N����x���W�l���`�F�b�N���Ă���l�ȏ�Ȃ�Α�����-1�ɂ��āC�͈͊O�ɔ�΂�
 * @param[inout] pos �p�[�e�B�N���ʒu
 * @param[inout] vel �p�[�e�B�N�����x
 * @param[in] dens �p�[�e�B�N�����x
 * @param[in] dt ���ԃX�e�b�v��
 * @param[in] nprts �p�[�e�B�N����
 */
void CuPbdSphCheckDeleteX(float* pos, float* vel, int* attr, float xmax, float farpoint[3], uint nprts)
{
	// 1�X���b�h/�p�[�e�B�N��
	uint numThreads, numBlocks;
	computeGridSize(nprts, THREAD_NUM, numBlocks, numThreads);

	float3 farpoint3;
	farpoint3.x = farpoint[0];
	farpoint3.y = farpoint[1];
	farpoint3.z = farpoint[2];

	// �J�[�l�����s
	checkDeleteXPB<<< numBlocks, numThreads >>>((float4*)pos, (float4*)vel, attr, xmax, farpoint3, nprts);
	
	RX_CUERROR("sphIntegrate kernel execution failed");	// �J�[�l�����s�G���[�`�F�b�N
	RX_CUCHECK(cudaThreadSynchronize());		// �S�ẴX���b�h���I���̂�҂�
}

/*!
 * �O���b�h��̖��x���Z�o
 * @param[out] dGridD �O���b�h��̖��x�l
 * @param[in] cell �p�[�e�B�N���O���b�h�f�[�^
 * @param[in] nx,ny �O���b�h��
 * @param[in] x0,y0 �O���b�h�ŏ����W
 * @param[in] dx,dy �O���b�h��
 */
void CuPbdSphGridDensity(float *dGridD, rxParticleCell cell, 
					  int nx, int ny, int nz, float x0, float y0, float z0, float dx, float dy, float dz)
{
#if USE_TEX
	RX_CUCHECK(cudaBindTexture(0, dSortedPosTex, cell.dSortedPos, cell.uNumParticles*sizeof(float4)));
	RX_CUCHECK(cudaBindTexture(0, dCellStartTex, cell.dCellStart, cell.uNumCells*sizeof(uint)));
	RX_CUCHECK(cudaBindTexture(0, dCellEndTex, cell.dCellEnd, cell.uNumCells*sizeof(uint)));	
#endif

	uint3  gnum = make_uint3(nx, ny, nz);
	float3 gmin = make_float3(x0, y0, z0);
	float3 glen = make_float3(dx, dy, dz);

	int numcell = gnum.x*gnum.y*gnum.z;

	int threads = 128;
	dim3 grid((numcell+threads-1)/threads, 1, 1);
	if(grid.x > 65535){
		grid.y = (grid.x+32768-1)/32768;
		grid.x = 32768;
	}

	// �J�[�l�����s
	pbdsphCalDensityInGrid<<<grid, threads>>>(dGridD, cell, gnum, gmin, glen);

	RX_CUERROR("Kernel execution failed");	// �J�[�l�����s�G���[�`�F�b�N
	RX_CUCHECK(cudaThreadSynchronize());		// �S�ẴX���b�h���I���̂�҂�

#if USE_TEX
	RX_CUCHECK(cudaUnbindTexture(dSortedPosTex));
	RX_CUCHECK(cudaUnbindTexture(dCellStartTex));
	RX_CUCHECK(cudaUnbindTexture(dCellEndTex));
#endif
}

/*!
 * �p�[�e�B�N���@���̌v�Z
 * @param[out] dNewDens �p�[�e�B�N�����x
 * @param[out] dNewPres �p�[�e�B�N������
 * @param[in]  cell �p�[�e�B�N���O���b�h�f�[�^
 */
void CuPbdSphNormal(float* dNrms, float* dDens, rxParticleCell cell)
{
	// MRK:CuSphNormal

	// 1�X���b�h/�p�[�e�B�N��
	uint numThreads, numBlocks;
	computeGridSize(cell.uNumParticles, THREAD_NUM, numBlocks, numThreads);

	// �J�[�l�����s
	pbdsphCalNormal<<< numBlocks, numThreads >>>((float4*)dNrms, dDens, cell);

	RX_CUERROR("sphCalNormal kernel execution failed");	// �J�[�l�����s�G���[�`�F�b�N
	RX_CUCHECK(cudaThreadSynchronize());		// �S�ẴX���b�h���I���̂�҂�

}





//-----------------------------------------------------------------------------
// Anisotropic Kernel
//-----------------------------------------------------------------------------
/*!
 * �J�[�l�����S�ʒu�̍X�V�Əd�ݕt�����ς̌v�Z(�J�[�l���֐�)
 * @param[out] dUpPos �X�V�J�[�l�����S
 * @param[out] dPosW �d�ݕt�����σp�[�e�B�N�����W 
 * @param[in]  lambda �������̂��߂̒萔
 * @param[in]  h �T�����a
 * @param[in]  cell �p�[�e�B�N���O���b�h�f�[�^
 */
void CuSphCalUpdatedPosition(float* dUpPos, float* dPosW, float lambda, float h, rxParticleCell cell)
{
	// MRK:CuSphCalUpdatedPosition
#if USE_TEX
	RX_CUCHECK(cudaBindTexture(0, dSortedPosTex, cell.dSortedPos, cell.uNumParticles*sizeof(float4)));
	RX_CUCHECK(cudaBindTexture(0, dCellStartTex, cell.dCellStart, cell.uNumCells*sizeof(uint)));
	RX_CUCHECK(cudaBindTexture(0, dCellEndTex, cell.dCellEnd, cell.uNumCells*sizeof(uint)));	
#endif

	// 1�X���b�h/�p�[�e�B�N��
	uint numThreads, numBlocks;
	computeGridSize(cell.uNumParticles, THREAD_NUM, numBlocks, numThreads);

	// �J�[�l�����s
	sphCalUpdatedPosition<<< numBlocks, numThreads >>>((float4*)dUpPos, (float4*)dPosW, lambda, h, cell);

	RX_CUERROR("sphCalUpdatedPosition kernel execution failed");	// �J�[�l�����s�G���[�`�F�b�N
	RX_CUCHECK(cudaThreadSynchronize());		// �S�ẴX���b�h���I���̂�҂�

#if USE_TEX
	RX_CUCHECK(cudaUnbindTexture(dSortedPosTex));
	RX_CUCHECK(cudaUnbindTexture(dCellStartTex));
	RX_CUCHECK(cudaUnbindTexture(dCellEndTex));
#endif
}

/*!
 * �������ʒu�ł̏d�ݕt�����ψʒu�̍Čv�Z��covariance matrix�̌v�Z
 * @param[out] dPosW �d�ݕt�����σp�[�e�B�N�����W 
 * @param[out] dCMat Covariance Matrix
 * @param[in]  h �T�����a
 * @param[in]  cell �p�[�e�B�N���O���b�h�f�[�^
 */
void CuSphCalCovarianceMatrix(float* dPosW, float* dCMat, float h, rxParticleCell cell)
{
	// MRK:CuSphCalCovarianceMatrix
#if USE_TEX
	RX_CUCHECK(cudaBindTexture(0, dSortedPosTex, cell.dSortedPos, cell.uNumParticles*sizeof(float4)));
	RX_CUCHECK(cudaBindTexture(0, dCellStartTex, cell.dCellStart, cell.uNumCells*sizeof(uint)));
	RX_CUCHECK(cudaBindTexture(0, dCellEndTex, cell.dCellEnd, cell.uNumCells*sizeof(uint)));	
#endif

	// 1�X���b�h/�p�[�e�B�N��
	uint numThreads, numBlocks;
	computeGridSize(cell.uNumParticles, THREAD_NUM, numBlocks, numThreads);

	// �J�[�l�����s
	sphCalCovarianceMatrix<<< numBlocks, numThreads >>>((float4*)dPosW, (matrix3x3*)dCMat, h, cell);

	RX_CUERROR("sphCalCovarianceMatrix kernel execution failed");	// �J�[�l�����s�G���[�`�F�b�N
	RX_CUCHECK(cudaThreadSynchronize());		// �S�ẴX���b�h���I���̂�҂�

#if USE_TEX
	RX_CUCHECK(cudaUnbindTexture(dSortedPosTex));
	RX_CUCHECK(cudaUnbindTexture(dCellStartTex));
	RX_CUCHECK(cudaUnbindTexture(dCellEndTex));
#endif
}

/*!
 * ���ْl�����ɂ��ŗL�l���v�Z
 * @param[in]  dC Covariance Matrix
 * @param[in]  dPosW �d�ݕt�����ψʒu
 * @param[out] dEigen �ŗL�l
 * @param[out] dR �ŗL�x�N�g��(��]�s��)
 * @param[in]  numParticles �p�[�e�B�N����
 */
void CuSphSVDecomposition(float* dC, float* dPosW, float* dEigen, float* dR, uint numParticles)
{
	// MRK:CuSphCalTransformMatrix

	// 1�X���b�h/�p�[�e�B�N��
	uint numThreads, numBlocks;
	computeGridSize(numParticles, THREAD_NUM, numBlocks, numThreads);

	// �J�[�l�����s
	sphSVDecomposition<<< numBlocks, numThreads >>>((matrix3x3*)dC, (float4*)dPosW, (float3*)dEigen, (matrix3x3*)dR, numParticles);

	RX_CUERROR("sphCalTransformMatrix kernel execution failed");	// �J�[�l�����s�G���[�`�F�b�N
	RX_CUCHECK(cudaThreadSynchronize());		// �S�ẴX���b�h���I���̂�҂�
}

/*!
 * �ŗL�l�C�ŗL�x�N�g��(��]�s��)����ό`�s����v�Z
 * @param[in]  dEigen �ŗL�l
 * @param[in]  dR �ŗL�x�N�g��(��]�s��)
 * @param[out] dG �ό`�s��
 * @param[in]  numParticles �p�[�e�B�N����
 */
void CuSphCalTransformMatrix(float* dEigen, float* dR, float *dG, uint numParticles)
{
	// MRK:CuSphCalTransformMatrix

	// 1�X���b�h/�p�[�e�B�N��
	uint numThreads, numBlocks;
	computeGridSize(numParticles, THREAD_NUM, numBlocks, numThreads);

	// �J�[�l�����s
	sphCalTransformMatrix<<< numBlocks, numThreads >>>((float3*)dEigen, (matrix3x3*)dR, (matrix3x3*)dG, numParticles);

	RX_CUERROR("sphCalTransformMatrix kernel execution failed");	// �J�[�l�����s�G���[�`�F�b�N
	RX_CUCHECK(cudaThreadSynchronize());		// �S�ẴX���b�h���I���̂�҂�
}


/*!
 * �O���b�h��̖��x���Z�o
 * @param[out] dGridD �O���b�h��̖��x�l
 * @param[in] cell �p�[�e�B�N���O���b�h�f�[�^
 * @param[in] nx,ny �O���b�h��
 * @param[in] x0,y0 �O���b�h�ŏ����W
 * @param[in] dx,dy �O���b�h��
 */
void CuSphGridDensityAniso(float *dGridD, float *dG, float Emax, rxParticleCell cell, 
						   int nx, int ny, int nz, float x0, float y0, float z0, float dx, float dy, float dz)
{
#if USE_TEX
	RX_CUCHECK(cudaBindTexture(0, dSortedPosTex, cell.dSortedPos, cell.uNumParticles*sizeof(float4)));
	RX_CUCHECK(cudaBindTexture(0, dCellStartTex, cell.dCellStart, cell.uNumCells*sizeof(uint)));
	RX_CUCHECK(cudaBindTexture(0, dCellEndTex, cell.dCellEnd, cell.uNumCells*sizeof(uint)));	
#endif

	uint3  gnum = make_uint3(nx, ny, nz);
	float3 gmin = make_float3(x0, y0, z0);
	float3 glen = make_float3(dx, dy, dz);

	int numcell = gnum.x*gnum.y*gnum.z;

	int threads = THREAD_NUM;
	dim3 grid((numcell+threads-1)/threads, 1, 1);
	if(grid.x > 65535){
		grid.y = (grid.x+32768-1)/32768;
		grid.x = 32768;
	}

	// �J�[�l�����s
	sphCalDensityAnisoInGrid<<<grid, threads>>>(dGridD, (matrix3x3*)dG, Emax, cell, gnum, gmin, glen);

	RX_CUERROR("Kernel execution failed");	// �J�[�l�����s�G���[�`�F�b�N
	RX_CUCHECK(cudaThreadSynchronize());		// �S�ẴX���b�h���I���̂�҂�

#if USE_TEX
	RX_CUCHECK(cudaUnbindTexture(dSortedPosTex));
	RX_CUCHECK(cudaUnbindTexture(dCellStartTex));
	RX_CUCHECK(cudaUnbindTexture(dCellEndTex));
#endif
}




//-----------------------------------------------------------------------------
// MARK:�E�F�[�u���b�g����
//-----------------------------------------------------------------------------

/*!
 * �p�[�e�B�N�����x�̃G�l���M�[�X�y�N�g�������v�Z
 * @param[out] dEt
 */
void CuSphES(float* dEt, float scale, float coef_et, float max_et, rxParticleCell cell, float3 cell_size)
{
	// MRK:CuSphES
#if USE_TEX
	RX_CUCHECK(cudaBindTexture(0, dSortedPosTex, cell.dSortedPos, cell.uNumParticles*sizeof(float4)));
	RX_CUCHECK(cudaBindTexture(0, dSortedVelTex, cell.dSortedVel, cell.uNumParticles*sizeof(float4)));
	RX_CUCHECK(cudaBindTexture(0, dCellStartTex, cell.dCellStart, cell.uNumCells*sizeof(uint)));
	RX_CUCHECK(cudaBindTexture(0, dCellEndTex, cell.dCellEnd, cell.uNumCells*sizeof(uint)));	
#endif

	float h = scale;//*MEXICAN_HAT_R;

	int na = cell.uNumArdGrid;
	na = (int)(h/cell_size.x)+1;
	//printf("CuSphES : na = %d, h = %f, cell = %f\n", na, h, cell_size.x);

	scale = h/MEXICAN_HAT_R;

	// 1�X���b�h/�p�[�e�B�N��
	uint numThreads, numBlocks;
	computeGridSize(cell.uNumParticles, THREAD_NUM, numBlocks, numThreads);

	// �J�[�l�����s
	calParticleES<<< numBlocks, numThreads >>>((float*)dEt, h, scale, coef_et, max_et, cell, na);

	// �J�[�l�����s�G���[�`�F�b�N
	RX_CUERROR("Kernel execution failed");
	RX_CUCHECK(cudaThreadSynchronize());		// �S�ẴX���b�h���I���̂�҂�

#if USE_TEX
	RX_CUCHECK(cudaUnbindTexture(dSortedPosTex));
	RX_CUCHECK(cudaUnbindTexture(dSortedVelTex));
	RX_CUCHECK(cudaUnbindTexture(dCellStartTex));
	RX_CUCHECK(cudaUnbindTexture(dCellEndTex));
#endif
}

//! �E�F�[�u���b�g�m�C�Y�^�C���̃Z�b�g
void CuSetWaveletTile3D(float *tile, int n, int d)
{
	if(!g_dNoiseTile[d] && g_udNoiseTileSize != n){
		// �f�o�C�X�������̊m�ۂƃz�X�g����̓]��
		int size = n*n*n*sizeof(float);
		RX_CUCHECK(cudaMalloc((void**)&g_dNoiseTile[d], size));
		RX_CUCHECK(cudaMemcpy((void*)g_dNoiseTile[d], (void*)tile, size, cudaMemcpyHostToDevice));

		g_udNoiseTileSize = n;
	}
}


//! �E�F�[�u���b�g�m�C�Y�^�C���̃Z�b�g
void CuSetWaveletTile3DB(float *tile, int nx, int ny, int nz, int d)
{
	if(g_dNoiseTile[d]) RX_CUCHECK(cudaFree(g_dNoiseTile[d]));

	// �f�o�C�X�������̊m�ۂƃz�X�g����̓]��
	int size = nx*ny*nz*sizeof(float);
	RX_CUCHECK(cudaMalloc((void**)&g_dNoiseTile[d], size));
	RX_CUCHECK(cudaMemcpy((void*)g_dNoiseTile[d], (void*)tile, size, cudaMemcpyHostToDevice));

	g_uNoiseTileNum[3*d+0] = nx;
	g_uNoiseTileNum[3*d+1] = ny;
	g_uNoiseTileNum[3*d+2] = nz;
}

//! �E�F�[�u���b�g�m�C�Y�^�C���̃Z�b�g
void CuSetWaveletTile3DBT(float *tile, int nx, int ny, int nz, int d)
{
	if(g_dNoiseTile[d]) RX_CUCHECK(cudaFree(g_dNoiseTile[d]));

	// �f�o�C�X�������̊m�ۂƃz�X�g����̓]��
	int size = nx*ny*nz*sizeof(float);
	RX_CUCHECK(cudaMalloc((void**)&g_dNoiseTile[d], size));
	RX_CUCHECK(cudaMemcpy((void*)g_dNoiseTile[d], (void*)tile, size, cudaMemcpyHostToDevice));

	cudaChannelFormatDesc cdesc = cudaCreateChannelDesc(32, 0, 0, 0, cudaChannelFormatKindFloat);
	RX_CUCHECK(cudaBindTexture(0, g_TexNoiseTile3D, g_dNoiseTile[d], cdesc) );

	g_uNoiseTileNum[3*d+0] = nx;
	g_uNoiseTileNum[3*d+1] = ny;
	g_uNoiseTileNum[3*d+2] = nz;
}

/*!
 * Wavelet turbulence���ԑ��x��ɒǉ�
 * @param[out] dVturb �������x��
 * @param[out] dFrc   �����ɂ���
 * @param[in] dEt   �Q��̃G�l���M�[�X�y�N�g����
 * @param[in] dDens �p�[�e�B�N�����x
 * @param[in] first,nbands ���d�ш�m�C�Y�̍ŏ��̑ш�Ƒш敝
 * @param[in] dPos  �p�[�e�B�N�����W
 * @param[in] pdim  �p�[�e�B�N���v�Z��Ԃ̑傫��
 * @param[in] nprts �p�[�e�B�N����
 * @param[in] dt    ���ԃX�e�b�v��
 */
void CuAddWTurb3D(float *dVturb, float *dFrc, float *dEt, float *dDens, int d, 
				  int first, int nbands, float* dPos, float pmin[3], float pdim[3], uint nprts, float dt)
{
	if(!g_dNoiseTile[0]) return;

	int3 tile_n;
	tile_n.x = g_uNoiseTileNum[3*d+0];
	tile_n.y = g_uNoiseTileNum[3*d+1];
	tile_n.z = g_uNoiseTileNum[3*d+2];
	
	int tile_size = tile_n.x*tile_n.y*tile_n.z;

	float3 pmin3, pdim3;
	pmin3.x = pmin[0];
	pmin3.y = pmin[1];
	pmin3.z = pmin[2];
	pdim3.x = pdim[0];
	pdim3.y = pdim[1];
	pdim3.z = pdim[2];

	// 1�X���b�h/�p�[�e�B�N��
	uint numThreads, numBlocks;
	computeGridSize(nprts, THREAD_NUM, numBlocks, numThreads);

	// �J�[�l�����s
	calWaveletTurbulence3D<<< numBlocks, numThreads >>>((float4*)dVturb, (float4*)dFrc, (float*)dEt, (float*)dDens, 
		first, nbands, g_dNoiseTile[0], tile_n, tile_size, (float4*)dPos, pmin3, pdim3, nprts, dt);

	// �J�[�l�����s�G���[�̃`�F�b�N
	RX_CUERROR("Kernel execution failed");
	RX_CUCHECK(cudaThreadSynchronize());
}




//-----------------------------------------------------------------------------
// MARK:Sub Particle
//-----------------------------------------------------------------------------
/*!
 * �T�u���q�̃G�l���M�[�X�y�N�g���̍X�V
 * @param[inout]dSubEt		�T�u���q�̃G�l���M�[�X�y�N�g���̔z��
 * @param[in]	dEt			���̃G�l���M�[�X�y�N�g��
 * @param[in]	scale		���̃X�P�[��
 * @param[in]	radius		���x��0�ł̔��a
 * @param[in]	maxSubLevel	�T�u���q�̍ő僌�x��
 * @param[in]	cell		�p�[�e�B�N���O���b�h�f�[�^
 
void CuSubUpdateEt(float* dSubEt,float* dEt,
			  float	scale,	float radius,
			  rxParticleCell cell)
{
	// 1�X���b�h/�p�[�e�B�N��
	uint numThreads, numBlocks;
	computeGridSize(cell.uNumParticles, 64, numBlocks, numThreads);

	// �J�[�l�����s
	updateSubEt<<< numBlocks, numThreads >>>(dSubEt, dEt, scale, radius, cell.uNumParticles);

	RX_CUERROR("updateSubEt kernel execution failed");	// �J�[�l�����s�G���[�`�F�b�N
	RX_CUCHECK(cudaThreadSynchronize());		// �S�ẴX���b�h���I���̂�҂�
}*/

/*!
 * �T�u���q��Et�̌v�Z
 * @param[inout]dSubEt		�T�u���q�̃G�l���M�[�X�y�N�g���ւ̃A�h���X
 * @param[in]�@	dEt			���q(���x��0)�̃G�l���M�[�X�y�N�g���ւ̃A�h���X
 * @param[in]	scale		���̃X�P�[��
 * @param[in]	radius		���x��0�ł̔��a
 * @param[in]	maxSubLevel	�T�u���q�̍ő僌�x��
 * @param[in]	cell		�p�[�e�B�N���O���b�h�f�[�^
 */
void CuSubCalEt(	float*	dSubEt,
					float*	dEt,
					float   et_coef, 
					float	scale, 
					float	radius,
					uint	maxSubLevel,
					rxParticleCell cell)
{
	// 1�X���b�h/�p�[�e�B�N��
	uint numThreads, numBlocks;
	computeGridSize(cell.uNumParticles, THREAD_NUM, numBlocks, numThreads);

	// �J�[�l�����s
	setSubEt<<< numBlocks, numThreads >>>(dSubEt, dEt, et_coef, scale, radius, cell.uNumParticles);
	RX_CUERROR("setSubPos kernel execution failed");	// �J�[�l�����s�G���[�`�F�b�N
	RX_CUCHECK(cudaThreadSynchronize());		// �S�ẴX���b�h���I���̂�҂�

	// �J�[�l�����s
	for(uint level = 0; level < maxSubLevel; level++){

		uint subBlockDim = calUintPow(2,level);

		for(uint subBlockIdx = 0; subBlockIdx < subBlockDim; subBlockIdx++){

			uint	subIndexP		= cell.uNumParticles * (calUintPow(2,level) - 1 + subBlockIdx);
			uint	subIndexC		= cell.uNumParticles * (calUintPow(2,level+1) - 1 + subBlockIdx*2);

			float* dSubBlockEtC	= &dSubEt[subIndexC];
			float* dSubBlockEtP	= &dSubEt[subIndexP];

			updateSubEt<<< numBlocks, numThreads >>>(dSubBlockEtC, dSubBlockEtP, cell.uNumParticles);
		}
		RX_CUERROR("updateSubEt kernel execution failed");	// �J�[�l�����s�G���[�`�F�b�N
		RX_CUCHECK(cudaThreadSynchronize());		// �S�ẴX���b�h���I���̂�҂�
		
	}

}

/*!
 * �T�u���q�̐�΍��W�̌v�Z
 * @param[inout]dSubPos		�T�u���q�̐�΍��W�ւ̃A�h���X
 * @param[in]	dPos		���q(���x��0)�̐�΍��W�ւ̃A�h���X
 * @param[in]	dSubChild	�T�u���q�̎q�ւ̒P�ʃx�N�g���ւ̃A�h���X
 * @param[in]	radius		���x��0�ł̔��a
 * @param[in]	maxSubLevel	�T�u���q�̍ő僌�x��
 * @param[in]	numParticles �p�[�e�B�N���O���b�h�f�[�^
 */
void CuSubCalPos(	float *dSubPos,
					float *dPos,
					float *dSubChild,
					float	radius,
					uint	maxSubLevel,
					uint	nprts)
{
	// 1�X���b�h/�p�[�e�B�N��
	uint numThreads, numBlocks;
	computeGridSize(nprts, THREAD_NUM, numBlocks, numThreads);

	// �J�[�l�����s
	setSubPos<<< numBlocks, numThreads >>>((float4*)dSubPos, (float4*)dPos, nprts);

	RX_CUERROR("setSubPos kernel execution failed");	// �J�[�l�����s�G���[�`�F�b�N
	RX_CUCHECK(cudaThreadSynchronize());		// �S�ẴX���b�h���I���̂�҂�

	
	for(uint level = 0; level < maxSubLevel; level++){

		uint subBlockDim = calUintPow(2,level);
		float radius_level = radius * pow(2.0f,-(float)level/3.0f);

		for(uint subBlockIdx = 0; subBlockIdx < subBlockDim; subBlockIdx++){

			uint	subIndexP		= nprts * (calUintPow(2,level)   - 1 + subBlockIdx);
			uint	subIndexC		= nprts * (calUintPow(2,level+1) - 1 + subBlockIdx*2);

			float4* dSubBlockPosC	= (float4*)&dSubPos[4*subIndexC];
			float4* dSubBlockPosP	= (float4*)&dSubPos[4*subIndexP];
			float4* dSubBlockChild	= (float4*)&dSubChild[4*subIndexP];
			// �J�[�l�����s
			calSubPos<<< numBlocks, numThreads >>>(dSubBlockPosC, dSubBlockPosP, dSubBlockChild,
				0.33f*radius_level, nprts);

		}
		RX_CUERROR("calSubPos kernel execution failed");	// �J�[�l�����s�G���[�`�F�b�N
		RX_CUCHECK(cudaThreadSynchronize());		// �S�ẴX���b�h���I���̂�҂�		
	}

}

/*!
 * �T�u���q�̎q1�ւ̒P�ʃx�N�g���Ǝ��̍X�V
 * @param[out]	dSubChild	�T�u���q�̎q�ւ̒P�ʃx�N�g���ւ̃A�h���X
 * @param[in]	dSubAxis	�T�u���q�̎��ւ̃A�h���X
 * @param[in]	dSubEt		�T�u���q�̃G�l���M�[�X�y�N�g���ւ̃A�h���X
 * @param[in]	radius		���x��0�ł̔��a
 * @param[in]	maxSubLevel	�T�u���q�̍ő僌�x��
 * @param[in]	cell		�p�[�e�B�N���O���b�h�f�[�^
 * @param[in]	dt			���ԃX�e�b�v��
 */
void CuSubUpdateChildAxis(	float	*dSubChild,
							float	*dSubAxis,
							float	*dSubEt,
							uint	*dSubRand,
							float	radius,
							float	ratioAxisAngle,
							uint	maxSubLevel,
							rxParticleCell cell,
							float	dt)
{
	// 1�X���b�h/�p�[�e�B�N��
	uint numThreads, numBlocks;
	computeGridSize(cell.uNumParticles, THREAD_NUM, numBlocks, numThreads);

	// �J�[�l�����s
	for(uint level = 0; level < maxSubLevel; level++){

		uint subBlockDim = calUintPow(2,level);
		float radius_level = radius * pow(2.0f,-(float)level/3.0f);

		for(uint subBlockIdx = 0; subBlockIdx < subBlockDim; subBlockIdx++){

			uint	subIndex		= cell.uNumParticles * (calUintPow(2,level) - 1 + subBlockIdx);

			float4* dSubBlockChild	= (float4*)&dSubChild[4*subIndex];
			float4* dSubBlockAxis	= (float4*)&dSubAxis[4*subIndex];
			float*	dSubBlockEt		= &dSubEt[subIndex];

			updateSubChildAxis<<< numBlocks, numThreads >>>(dSubBlockChild, dSubBlockAxis, dSubBlockEt, dSubRand, 
				radius_level, cell.uNumParticles, ratioAxisAngle, dt);
			
		}		
	}
	RX_CUERROR("updateSubChildAxis kernel execution failed");	// �J�[�l�����s�G���[�`�F�b�N
	RX_CUCHECK(cudaThreadSynchronize());		// �S�ẴX���b�h���I���̂�҂�

}

/*!
 * �T�u���q�̏�����
 * @param[inout]dSubPos		�T�u���q�̐�΍��W�ւ̃A�h���X
 * @param[in]	dPos		���q(���x��0)�̐�΍��W�ւ̃A�h���X
 * @param[inout]dSubAxis	�T�u���q�̎��ւ̃A�h���X
 * @param[inout]dSubChild	�T�u���q�̎q�ւ̒P�ʃx�N�g���ւ̃A�h���X
 * @param[in]	radius		���x��0�ł̔��a
 * @param[in]	maxSubLevel	�T�u���q�̍ő僌�x��
 * @param[in]	nprts		�p�[�e�B�N����
 */
void CuSubInit(	float *dSubPos,
				float *dPos,
				float *dSubAxis,
				float *dSubChild,
				uint*	dSubRand,
				float	radius,
				uint	maxSubLevel,
				uint nprts)
{
	// 1�X���b�h/�p�[�e�B�N��
	uint numThreads, numBlocks;
	computeGridSize(nprts, THREAD_NUM, numBlocks, numThreads);

	initSubRand<<< numBlocks, numThreads >>>(dSubRand, nprts);

	// �J�[�l�����s
	for(uint level = 0; level < maxSubLevel; level++){

		uint subBlockDim = calUintPow(2,level);

		for(uint subBlockIdx = 0; subBlockIdx < subBlockDim; subBlockIdx++){

			uint	subIndex		= nprts * (calUintPow(2,level) - 1 + subBlockIdx);
			float4* dSubBlockChild	= (float4*)(&dSubChild[4*subIndex]);
			float4* dSubBlockAxis	= (float4*)(&dSubAxis[4*subIndex]);

			// �J�[�l�����s
			initSub<<< numBlocks, numThreads >>>(dSubBlockChild, dSubBlockAxis, dSubRand, nprts);		
		}
	}
	RX_CUERROR("initSub kernel execution failed");	// �J�[�l�����s�G���[�`�F�b�N
	RX_CUCHECK(cudaThreadSynchronize());		// �S�ẴX���b�h���I���̂�҂�

	//��΍��W���v�Z
	CuSubCalPos(dSubPos, dPos, dSubChild, radius, maxSubLevel, nprts);

}

//-----------------------------------------------------------------------------

/*!
 * �T�u���q�̏�����
 * @param[inout]dSubPos		�T�u���q�̐�΍��W�ւ̃A�h���X
 * @param[in]	dPos		���q(���x��0)�̐�΍��W�ւ̃A�h���X
 * @param[inout]dSubAxis	�T�u���q�̎��ւ̃A�h���X
 * @param[inout]dSubChild	�T�u���q�̎q�ւ̒P�ʃx�N�g���ւ̃A�h���X
 * @param[in]	radius		���x��0�ł̔��a
 * @param[in]	maxSubLevel	�T�u���q�̍ő僌�x��
 * @param[in]	nprts		�p�[�e�B�N����
 */
void CuSubInit2(float	*dSubPos,
				float	*dPos,
				float	*dSubAxis,
				float	*dSubChild,
				uint	*dSubRand,
				float	radius,
				uint	maxSubLevel,
				uint	nprts, 
				uint	uMaxParticles)
{
	// 1�X���b�h/�p�[�e�B�N��
	uint numThreads, numBlocks;
	computeGridSize(nprts, THREAD_NUM, numBlocks, numThreads);

	// ���x��0�̃p�[�e�B�N���ɗ����̊��ݒ�
	initSubRand<<< numBlocks, numThreads >>>(dSubRand, nprts);
	RX_CUCHECK(cudaThreadSynchronize());

	// ���x��0�̃p�[�e�B�N���ʒu��ݒ�
	setSubPos<<< numBlocks, numThreads >>>((float4*)dSubPos, (float4*)dPos, nprts);
	RX_CUCHECK(cudaThreadSynchronize());		// �S�ẴX���b�h���I���̂�҂�

	float4 *dSubPosP, *dSubPosC1, *dSubPosC2, *dSubChildP, *dSubAxisP;

	// �J�[�l�����s
	for(uint level = 0; level < maxSubLevel; level++){
		uint subBlockDim = calUintPow(2, level);

		for(uint subBlockIdx = 0; subBlockIdx < subBlockDim; subBlockIdx++){

			uint block_idx = subBlockDim-1+subBlockIdx;
			dSubPosP	= (float4*)(&dSubPos[4 * uMaxParticles * block_idx]);		// �T�u�p�[�e�B�N���Q(L,j)
			dSubPosC1	= (float4*)(&dSubPos[4 * uMaxParticles * (2*block_idx+1)]);	// �q(L+1,0)
			dSubPosC2	= (float4*)(&dSubPos[4 * uMaxParticles * (2*block_idx+2)]);	// �q(L+1,1)
			dSubChildP	= (float4*)(&dSubChild[4 * uMaxParticles * block_idx]);		// �q(j=0)�ւ̒P�ʃx�N�g��
			dSubAxisP	= (float4*)(&dSubAxis[4 * uMaxParticles * block_idx]);		// ��]��

			// �J�[�l�����s
			initSubParticle<<< numBlocks, numThreads >>>(dSubPosP, dSubPosC1, dSubPosC2, dSubChildP, dSubAxisP,
				dSubRand, radius, nprts, uMaxParticles);

		}
		RX_CUERROR("initSub2 kernel execution failed");	// �J�[�l�����s�G���[�`�F�b�N
		RX_CUCHECK(cudaThreadSynchronize());		// �S�ẴX���b�h���I���̂�҂�
	}

}
/*!
 * �T�u���q�̏�����
 * @param[inout]dSubPos		�T�u���q�̐�΍��W�ւ̃A�h���X
 * @param[in]	dPos		���q(���x��0)�̐�΍��W�ւ̃A�h���X
 * @param[inout]dSubAxis	�T�u���q�̎��ւ̃A�h���X
 * @param[inout]dSubChild	�T�u���q�̎q�ւ̒P�ʃx�N�g���ւ̃A�h���X
 * @param[in]	radius		���x��0�ł̔��a
 * @param[in]	maxSubLevel	�T�u���q�̍ő僌�x��
 * @param[in]	nprts		�p�[�e�B�N����
 */
void CuSubAdd(float	*dSubPos,
			  float	*dPos,
			  float	*dSubAxis,
			  float	*dSubChild,
			  uint	*dSubRand,
			  float	radius,
			  uint	maxSubLevel,
			  uint	nprts, 
			  int	uMaxParticles, 
			  uint	uStart, 
			  uint	uCount)
{
	// 1�X���b�h/�p�[�e�B�N��
	uint numThreads, numBlocks;
	// FIX:�u���b�N�������Ȃ��Ƃ��܂������Ȃ��H
	computeGridSize(uMaxParticles, THREAD_NUM, numBlocks, numThreads);
	//printf("block %d, thread %d\n", numBlocks, numThreads);

	// ���x��0�̃p�[�e�B�N���ɗ����̊��ݒ�
	addSubRand<<< numBlocks, numThreads >>>(dSubRand, nprts, uStart, uCount);
	RX_CUCHECK(cudaThreadSynchronize());

	// ���x��0�̃p�[�e�B�N���ʒu��ݒ�
	addSubPos<<< numBlocks, numThreads >>>((float4*)dSubPos, (float4*)dPos, nprts, uStart, uCount);
	RX_CUCHECK(cudaThreadSynchronize());		// �S�ẴX���b�h���I���̂�҂�

	float4 *dSubPosP, *dSubPosC1, *dSubPosC2, *dSubChildP, *dSubAxisP;

	// �J�[�l�����s
	for(uint level = 0; level < maxSubLevel; level++){
		uint subBlockDim = calUintPow(2, level);

		for(uint subBlockIdx = 0; subBlockIdx < subBlockDim; subBlockIdx++){

			uint block_idx = subBlockDim-1+subBlockIdx;
			dSubPosP	= (float4*)(&dSubPos[4 * uMaxParticles * block_idx]);		// �T�u�p�[�e�B�N���Q(L,j)
			dSubPosC1	= (float4*)(&dSubPos[4 * uMaxParticles * (2*block_idx+1)]);	// �q(L+1,0)
			dSubPosC2	= (float4*)(&dSubPos[4 * uMaxParticles * (2*block_idx+2)]);	// �q(L+1,1)
			dSubChildP	= (float4*)(&dSubChild[4 * uMaxParticles * block_idx]);		// �q(j=0)�ւ̒P�ʃx�N�g��
			dSubAxisP	= (float4*)(&dSubAxis[4 * uMaxParticles * block_idx]);		// ��]��

			// �J�[�l�����s
			addSubParticle<<< numBlocks, numThreads >>>(dSubPosP, dSubPosC1, dSubPosC2, dSubChildP, dSubAxisP, dSubRand, 
														radius, nprts, uMaxParticles, uStart, uCount);

			RX_CUERROR("addSubParticle kernel execution failed");	// �J�[�l�����s�G���[�`�F�b�N
			RX_CUCHECK(cudaThreadSynchronize());		// �S�ẴX���b�h���I���̂�҂�
		}
	}
}

/*!
 * �T�u���q�̎q1�ւ̒P�ʃx�N�g���Ǝ��̍X�V
 * @param[out]	dSubChild	�T�u���q�̎q�ւ̒P�ʃx�N�g���ւ̃A�h���X
 * @param[in]	dSubAxis	�T�u���q�̎��ւ̃A�h���X
 * @param[in]	dSubEt		�T�u���q�̃G�l���M�[�X�y�N�g���ւ̃A�h���X
 * @param[in]	radius		���x��0�ł̔��a
 * @param[in]	maxSubLevel	�T�u���q�̍ő僌�x��
 * @param[in]	cell		�p�[�e�B�N���O���b�h�f�[�^
 * @param[in]	dt			���ԃX�e�b�v��
 */
void CuSubUpdate(	float	*dPos,
					float	*dEt,
					float	*dSubPos,
					float	*dSubChild,
					float	*dSubAxis,
					float	*dSubEt,
					uint	*dSubRand,
					float   et_coef, 
					float	radius,
					float	ratioAxisAngle,
					uint	maxSubLevel,
					uint	nprts,
					uint	uMaxParticles, 
					float	scale,
					float	dt)
{
	// 1�X���b�h/�p�[�e�B�N��
	uint grid, block;	// �O���b�h���u���b�N��,�u���b�N���X���b�h��
	block = THREAD_NUM;
	grid = DivCeil(nprts, block);

	// ���x��0�̃T�u�p�[�e�B�N���̈ʒu�Ɨ����G�l���M�[�l���X�V
	setSubEt<<<grid, block>>>(dSubEt, dEt, et_coef, scale, radius, nprts);
	setSubPos<<<grid, block>>>((float4*)dSubPos, (float4*)dPos, nprts);

	RX_CUERROR("setSubPos kernel execution failed");	// �J�[�l�����s�G���[�`�F�b�N
	RX_CUCHECK(cudaThreadSynchronize());		// �S�ẴX���b�h���I���̂�҂�

	float4 *dSubPosP, *dSubPosC1, *dSubPosC2, *dSubChildP, *dSubAxisP;

	// �J�[�l�����s
	for(uint level = 0; level < maxSubLevel; level++){
		uint subBlockDim = calUintPow(2,level);

		float radius_level = radius * pow(2.0f,-(float)level/3.0f);
		float ratio_et	= pow(2.0f, -(float)level*5.0f/9.0f);

		for(uint subBlockIdx = 0; subBlockIdx < subBlockDim; subBlockIdx++){

			//printf("level = %d / subBlockIdx = %d \n", level, subBlockIdx);
			uint blockIdx = subBlockDim-1+subBlockIdx;
			dSubPosP	= (float4*)(&dSubPos[4 * uMaxParticles * blockIdx]);		// �T�u�p�[�e�B�N���Q(L,j)
			dSubPosC1	= (float4*)(&dSubPos[4 * uMaxParticles * (2*blockIdx+1)]);	// �q(L+1,0)
			dSubPosC2	= (float4*)(&dSubPos[4 * uMaxParticles * (2*blockIdx+2)]);	// �q(L+1,1)
			dSubChildP	= (float4*)(&dSubChild[4 * uMaxParticles * blockIdx]);		// �q(j=0)�ւ̒P�ʃx�N�g��
			dSubAxisP	= (float4*)(&dSubAxis[ 4 * uMaxParticles * blockIdx]);		// ��]��

			// (level, j)�̃T�u�p�[�e�B�N���̈ʒu���X�V
			updateSubParticle<<<grid, block>>>(dSubPosP, dSubPosC1, dSubPosC2, dSubChildP, dSubAxisP,
				dSubEt, dSubRand, ratio_et, radius_level, nprts, uMaxParticles, dt);

			
		}
		RX_CUERROR("CuSubUpdate kernel execution failed");
		RX_CUCHECK(cudaThreadSynchronize());	// �S�ẴX���b�h���I���̂�҂�
		
	}
}


/*!
 * �����_�����O�Ɏg�p����T�u�p�[�e�B�N�����x���C�e���W�����v�Z
 * @param[in] dEt ���x��0�̃p�[�e�B�N���̃G�l���M�[�l
 * @param[in] dSubPos �T�u�p�[�e�B�N���̈ʒu
 * @param[out] subCell �T�u�p�[�e�B�N�����
 * @param[in] et_coef ���x��0�̃p�[�e�B�N���̃G�l���M�[�l�Ɋ|����W��
 */
void CuSubSetUnsortArray(float* dEt, float* dSubPos, rxSubParticleCell &subCell, float et_coef)
{
	// 1�X���b�h/�e�p�[�e�B�N��
	uint grid, block;	// �O���b�h���u���b�N��,�u���b�N���X���b�h��
	block = THREAD_NUM;
	grid = DivCeil(subCell.uNumParticles, block);

	RX_CUCHECK(cudaMemset(subCell.dSubOcc, 0, subCell.uSubNumAllParticles*sizeof(uint)));

	setSubUnsortArray<<<grid, block>>>(subCell, dEt, (float4*)dSubPos, et_coef);

	RX_CUERROR("CuSubSetUnsortArray kernel execution failed");
	RX_CUCHECK(cudaThreadSynchronize());	// �S�ẴX���b�h���I���̂�҂�

	subCell.uSubNumValidParticles = subCell.uSubNumAllParticles;


	int size = subCell.uSubNumAllParticles;

	// �T�u�p�[�e�B�N���L��/������Scan
	thrust::exclusive_scan(thrust::device_ptr<unsigned int>(subCell.dSubOcc), 
						   thrust::device_ptr<unsigned int>(subCell.dSubOcc+size),
						   thrust::device_ptr<unsigned int>(subCell.dSubOccScan));

	// Exclusive scan (�Ō�̗v�f��0�Ԗڂ���n-2�Ԗڂ܂ł̍��v�ɂȂ��Ă���)�Ȃ̂ŁC
	// Scan�O�z��̍Ō�(n-1�Ԗ�)�̗v�f�ƍ��v���邱�ƂŗL���T�u�p�[�e�B�N�������v�Z
	uint lval, lsval;
	RX_CUCHECK(cudaMemcpy((void*)&lval, (void*)(subCell.dSubOcc+size-1), sizeof(uint), cudaMemcpyDeviceToHost));
	RX_CUCHECK(cudaMemcpy((void*)&lsval, (void*)(subCell.dSubOccScan+size-1), sizeof(uint), cudaMemcpyDeviceToHost));
	uint num_valid_particles = lval+lsval;

	//printf("num of valid sub-particles = %d / %d\n", num_valid_particles, subCell.uSubNumAllParticles);

	// 1�X���b�h/�S�T�u�p�[�e�B�N��
	block = THREAD_NUM;
	grid = DivCeil(subCell.uSubNumAllParticles, block);

	// �L���ȃp�[�e�B�N�������l�߂�
	compactSubParticles<<<grid, block>>>(subCell, size);
	
	RX_CUERROR("compactSubParticles kernel execution failed");	// �J�[�l�����s�G���[�`�F�b�N
	RX_CUCHECK(cudaThreadSynchronize());		// �S�ẴX���b�h���I���̂�҂�

	RX_CUCHECK(cudaMemcpy(subCell.dSubUnsortPos, subCell.dSubSortedPos, num_valid_particles*sizeof(float4), cudaMemcpyDeviceToDevice));
	RX_CUCHECK(cudaMemcpy(subCell.dSubUnsortRad, subCell.dSubSortedRad, num_valid_particles*sizeof(float),  cudaMemcpyDeviceToDevice));
	RX_CUCHECK(cudaMemcpy(subCell.dSubUnsortRat, subCell.dSubSortedRat, num_valid_particles*sizeof(float),  cudaMemcpyDeviceToDevice));
	subCell.uSubNumValidParticles = num_valid_particles;
}


void CuSubCalcHash(rxSubParticleCell subCell)
{
	uint numThreads, numBlocks;
	computeGridSize(subCell.uSubNumValidParticles, THREAD_NUM, numBlocks, numThreads);

	// �J�[�l�����s
	calcHashD3<<< numBlocks, numThreads >>>(subCell.dSubGridParticleHash,
										   subCell.dSubSortedIndex,
										   subCell.dSubUnsortPos,
										   subCell.uSubNumValidParticles);
	
	RX_CUERROR("Kernel execution failed");	// �J�[�l���G���[�`�F�b�N
}

void CuSubCheckRatio(rxSubParticleCell subCell)
{
	uint numThreads, numBlocks;
	computeGridSize(subCell.uSubNumValidParticles, THREAD_NUM, numBlocks, numThreads);

	checkSubRatio<<< numBlocks, numThreads >>>(subCell.dSubGridParticleHash,
				subCell.dSubUnsortRat, 
				subCell.uSubNumCells,
				subCell.uSubNumValidParticles);

	RX_CUERROR("Kernel execution failed: CuSubCheckRatio");
	RX_CUCHECK(cudaThreadSynchronize());

}


/*!
 * �p�[�e�B�N���z����\�[�g���ꂽ���Ԃɕ��ёւ��C
 * �e�Z���̎n�܂�ƏI���̃C���f�b�N�X������
 * @param[in] cell �p�[�e�B�N���O���b�h�f�[�^
 * @param[in] oldPos �p�[�e�B�N���ʒu
 * @param[in] oldVel �p�[�e�B�N�����x
 */
void CuSubReorderDataAndFindCellStart(rxSubParticleCell subCell)
{
	uint numThreads, numBlocks;
	computeGridSize(subCell.uSubNumValidParticles, THREAD_NUM, numBlocks, numThreads);

	RX_CUCHECK(cudaMemset(subCell.dSubCellStart, 0xffffffff, subCell.uSubNumCells*sizeof(uint)));

#if USE_TEX//�e�N�X�`���̖��O�ɒ���
	RX_CUCHECK(cudaBindTexture(0, dSubUnsortPosTex, subCell.dSubUnsortPos, subCell.uSubNumValidParticles*sizeof(float4)));
	RX_CUCHECK(cudaBindTexture(0, dSubUnsortRadTex, subCell.dSubUnsortRad, subCell.uSubNumValidParticles*sizeof(float)));
	RX_CUCHECK(cudaBindTexture(0, dSubUnsortRatTex, subCell.dSubUnsortRat, subCell.uSubNumValidParticles*sizeof(float)));
#endif

	uint smemSize = sizeof(uint)*(numThreads+1);

	
	// �J�[�l�����s
	reorderDataAndFindCellStartF4F1F1<<< numBlocks, numThreads, smemSize>>>(subCell);

	RX_CUERROR("Kernel execution failed: CuSubReorderDataAndFindCellStart");
	RX_CUCHECK(cudaThreadSynchronize());		// �S�ẴX���b�h���I���̂�҂�


#if USE_TEX
	RX_CUCHECK(cudaUnbindTexture(dSubUnsortPosTex));
	RX_CUCHECK(cudaUnbindTexture(dSubUnsortRadTex));
	RX_CUCHECK(cudaUnbindTexture(dSubUnsortRatTex));
#endif
}

/*!
 * �O���b�h��̖��x���Z�o
 * @param[out] dGridD �O���b�h��̖��x�l
 * @param[in] cell �p�[�e�B�N���O���b�h�f�[�^
 * @param[in] nx,ny �O���b�h��
 * @param[in] x0,y0 �O���b�h�ŏ����W
 * @param[in] dx,dy �O���b�h��
 */
void CuSubSphGridDensity(float *dGridD, rxSubParticleCell subCell, 
					  int nx, int ny, int nz, float x0, float y0, float z0, float dx, float dy, float dz)
{
#if USE_TEX
	RX_CUCHECK(cudaBindTexture(0, dSubSortedPosTex,	subCell.dSubSortedPos,	subCell.uSubNumValidParticles*sizeof(float4)));
	RX_CUCHECK(cudaBindTexture(0, dSubSortedRadTex,	subCell.dSubSortedRad,	subCell.uSubNumValidParticles*sizeof(float)));
	RX_CUCHECK(cudaBindTexture(0, dSubSortedRatTex,	subCell.dSubSortedRat,	subCell.uSubNumValidParticles*sizeof(float)));
	RX_CUCHECK(cudaBindTexture(0, dSubCellStartTex,	subCell.dSubCellStart,	subCell.uSubNumCells*sizeof(uint)));
	RX_CUCHECK(cudaBindTexture(0, dSubCellEndTex,	subCell.dSubCellEnd,	subCell.uSubNumCells*sizeof(uint)));	
#endif

	uint3  gnum = make_uint3(nx, ny, nz);
	float3 gmin = make_float3(x0, y0, z0);
	float3 glen = make_float3(dx, dy, dz);

	int numcell = gnum.x * gnum.y * gnum.z;

	int threads = 128;
	dim3 grid(numcell/threads, 1, 1);
	if(grid.x > 65535){
		grid.y = grid.x/32768;
		grid.x = 32768;
	}

	// �J�[�l�����s
	calSubGridDensity<<<grid, threads>>>(dGridD, subCell, gnum, gmin, glen);

	RX_CUERROR("Kernel execution failed");	// �J�[�l�����s�G���[�`�F�b�N
	RX_CUCHECK(cudaThreadSynchronize());		// �S�ẴX���b�h���I���̂�҂�

#if USE_TEX
	RX_CUCHECK(cudaUnbindTexture(dSubSortedPosTex));
	RX_CUCHECK(cudaUnbindTexture(dSubSortedRadTex));
	RX_CUCHECK(cudaUnbindTexture(dSubSortedRatTex));
	RX_CUCHECK(cudaUnbindTexture(dSubCellStartTex));
	RX_CUCHECK(cudaUnbindTexture(dSubCellEndTex));
#endif
	
}

uint CuCheckNumValidHashData(rxSubParticleCell subCell)
{
	uint numThreads, numBlocks;
	computeGridSize( subCell.uSubNumAllParticles, THREAD_NUM, numBlocks, numThreads);

	uint *dNum;
	uint hNum;

	hNum = subCell.uSubNumAllParticles;

	RX_CUCHECK(cudaMalloc((void**)&dNum, sizeof(uint)));

	RX_CUCHECK(cudaMemcpy(dNum, &hNum, sizeof(uint), cudaMemcpyHostToDevice));

	uint smemSize = sizeof(uint)*(numThreads+1);
	checkNumUintData<<< numBlocks, numThreads, smemSize>>>(subCell.dSubGridParticleHash, dNum, subCell.uSubNumCells, subCell.uSubNumAllParticles);

	RX_CUERROR("Kernel execution failed");	// �J�[�l�����s�G���[�`�F�b�N
	RX_CUCHECK(cudaThreadSynchronize());		// �S�ẴX���b�h���I���̂�҂�

	RX_CUCHECK(cudaMemcpy(&hNum, dNum, sizeof(uint), cudaMemcpyDeviceToHost));

	
	//free(hNum);
	RX_CUCHECK(cudaFree(dNum));
	//printf("[CuFindNumValidHashData]Number of Valid Hash Data : %d\n",subCell.uSubNumMCParticles);
	return hNum;

}





}   // extern "C"
