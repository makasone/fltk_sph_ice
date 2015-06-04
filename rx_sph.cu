/*! 
  @file rx_sph.cu
	
  @brief CUDAによるSPH

  @author Makoto Fujisawa
  @date 2009-08, 2011-06
*/
// FILE --rx_sph.cu--


//-----------------------------------------------------------------------------
// インクルードファイル
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
// MARK:グローバル変数
//-----------------------------------------------------------------------------
cudaArray *g_caNoiseTile = 0;
float *g_dNoiseTile[3] = {0, 0, 0};
uint g_udNoiseTileSize = 0;
uint g_uNoiseTileNum[3*3] = {0, 0, 0,  0, 0, 0,  0, 0, 0};


//-----------------------------------------------------------------------------
// CUDA関数
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
 * thrust::exclusive_scanの呼び出し
 * @param[out] dScanData scan後のデータ
 * @param[in] dData 元データ
 * @param[in] num データ数
 */
void CuScanf(float* dScanData, float* dData, unsigned int num)
{
	thrust::exclusive_scan(thrust::device_ptr<float>(dData), 
						   thrust::device_ptr<float>(dData+num),
						   thrust::device_ptr<float>(dScanData));
}


/*!
 * 分割セルのハッシュを計算
 * @param[in] 
 * @return 
 */
void CuCalcHash(uint* dGridParticleHash, uint* dSortedIndex, float* dPos, int *attr, int nprts)
{
	uint numThreads, numBlocks;
	computeGridSize(nprts, THREAD_NUM, numBlocks, numThreads);

	// カーネル実行
	calcHashD2<<< numBlocks, numThreads >>>(dGridParticleHash,
										   dSortedIndex,
										   (float4*)dPos,
										   attr, 
										   nprts);
	//calcHashD<<< numBlocks, numThreads >>>(dGridParticleHash,
	//									   dSortedIndex,
	//									   (float4*)dPos,
	//									   nprts);
	
	RX_CUERROR("Kernel execution failed");	// カーネルエラーチェック
}


/*!
 * パーティクル配列をソートされた順番に並び替え，
 * 各セルの始まりと終わりのインデックスを検索
 * @param[in] cell パーティクルグリッドデータ
 * @param[in] oldPos パーティクル位置
 * @param[in] oldVel パーティクル速度
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

	// カーネル実行
	reorderDataAndFindCellStartD<<< numBlocks, numThreads, smemSize>>>(cell, (float4*)oldPos, (float4*)oldVel);

	RX_CUERROR("Kernel execution failed: CuReorderDataAndFindCellStartD");
	RX_CUCHECK(cudaThreadSynchronize());		// 全てのスレッドが終わるのを待つ

#if USE_TEX
	RX_CUCHECK(cudaUnbindTexture(dSortedPosTex));
	RX_CUCHECK(cudaUnbindTexture(dSortedVelTex));
#endif
}


/*!
 * 分割セルのハッシュを計算
 * @param[in] 
 * @return 
 */
void CuCalcHashB(uint* dGridParticleHash, uint* dSortedIndex, float* dPos, 
				 float3 world_origin, float3 cell_width, uint3 grid_size, int nprts)
{
	uint numThreads, numBlocks;
	computeGridSize(nprts, THREAD_NUM, numBlocks, numThreads);

	// カーネル実行
	calcHashB<<< numBlocks, numThreads >>>(dGridParticleHash,
										   dSortedIndex,
										   (float4*)dPos,
										   world_origin, 
										   cell_width, 
										   grid_size, 
										   nprts);
	
	RX_CUERROR("Kernel execution failed : calcHashB");	// カーネルエラーチェック
}

/*!
 * パーティクル配列をソートされた順番に並び替え，
 * 各セルの始まりと終わりのインデックスを検索
 * @param[in] cell パーティクルグリッドデータ
 * @param[in] oldPos パーティクル位置
 */
void CuReorderDataAndFindCellStartB(rxParticleCell cell, float* oldPos)
{
	uint numThreads, numBlocks;
	computeGridSize(cell.uNumParticles, THREAD_NUM, numBlocks, numThreads);

	RX_CUCHECK(cudaMemset(cell.dCellStart, 0xffffffff, cell.uNumCells*sizeof(uint)));

	uint smemSize = sizeof(uint)*(numThreads+1);

	// カーネル実行
	reorderDataAndFindCellStartB<<< numBlocks, numThreads, smemSize>>>(cell, (float4*)oldPos);

	RX_CUERROR("Kernel execution failed: CuReorderDataAndFindCellStartB");
	RX_CUCHECK(cudaThreadSynchronize());		// 全てのスレッドが終わるのを待つ
}


//-----------------------------------------------------------------------------
// MARK:3D SPH
//-----------------------------------------------------------------------------



/*!
 * パーティクル密度の計算(カーネル呼び出し)
 * @param[out] dDens パーティクル密度
 * @param[out] dPres パーティクル圧力
 * @param[in]  cell パーティクルグリッドデータ
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

	// 1スレッド/パーティクル
	uint numThreads, numBlocks;
	computeGridSize(cell.uNumParticles, THREAD_NUM, numBlocks, numThreads);

	// カーネル実行
	sphCalDensity<<< numBlocks, numThreads >>>(dDens, dPres, cell);

	RX_CUERROR("sphCalDensity kernel execution failed");	// カーネル実行エラーチェック
	RX_CUCHECK(cudaThreadSynchronize());		// 全てのスレッドが終わるのを待つ

#if USE_TEX
	RX_CUCHECK(cudaUnbindTexture(dSortedPosTex));
	RX_CUCHECK(cudaUnbindTexture(dCellStartTex));
	RX_CUCHECK(cudaUnbindTexture(dCellEndTex));
#endif
}

/*!
 * パーティクル法線の計算
 * @param[out] dNewDens パーティクル密度
 * @param[out] dNewPres パーティクル圧力
 * @param[in]  cell パーティクルグリッドデータ
 */
void CuSphNormal(float* dNrms, float* dDens, rxParticleCell cell)
{
	// MRK:CuSphNormal

	// 1スレッド/パーティクル
	uint numThreads, numBlocks;
	computeGridSize(cell.uNumParticles, THREAD_NUM, numBlocks, numThreads);

	// カーネル実行
	sphCalNormal<<< numBlocks, numThreads >>>((float4*)dNrms, dDens, cell);

	RX_CUERROR("sphCalNormal kernel execution failed");	// カーネル実行エラーチェック
	RX_CUCHECK(cudaThreadSynchronize());		// 全てのスレッドが終わるのを待つ

}

/*!
 * パーティクルにかかる力の計算(カーネル呼び出し)
 * @param[in] dDens パーティクル密度
 * @param[in] dPres パーティクル圧力
 * @param[out] dFrc パーティクルにかかる力
 * @param[in] cell パーティクルグリッドデータ
 * @param[in] dt 時間ステップ幅
 */
void CuSphForces(float* dDens, float* dPres, float* dFrc, rxParticleCell cell, float dt)
{
#if USE_TEX
	RX_CUCHECK(cudaBindTexture(0, dSortedPosTex, cell.dSortedPos, cell.uNumParticles*sizeof(float4)));
	RX_CUCHECK(cudaBindTexture(0, dSortedVelTex, cell.dSortedVel, cell.uNumParticles*sizeof(float4)));
	RX_CUCHECK(cudaBindTexture(0, dCellStartTex, cell.dCellStart, cell.uNumCells*sizeof(uint)));
	RX_CUCHECK(cudaBindTexture(0, dCellEndTex, cell.dCellEnd, cell.uNumCells*sizeof(uint)));	
#endif

	// 1スレッド/パーティクル
	uint numThreads, numBlocks;
	computeGridSize(cell.uNumParticles, THREAD_NUM, numBlocks, numThreads);

	// カーネル実行
	sphCalForces<<< numBlocks, numThreads >>>(dDens, dPres, (float4*)dFrc, cell);

	RX_CUERROR("calForcesSPH kernel execution failed");	// カーネル実行エラーチェック
	RX_CUCHECK(cudaThreadSynchronize());		// 全てのスレッドが終わるのを待つ

#if USE_TEX
	RX_CUCHECK(cudaUnbindTexture(dSortedPosTex));
	RX_CUCHECK(cudaUnbindTexture(dSortedVelTex));
	RX_CUCHECK(cudaUnbindTexture(dCellStartTex));
	RX_CUCHECK(cudaUnbindTexture(dCellEndTex));
#endif
}

/*!
 * パーティクル位置，速度の更新
 * @param[inout] pos パーティクル位置
 * @param[inout] vel パーティクル速度
 * @param[inout] velOld 前ステップのパーティクル速度
 * @param[in] frc パーティクルにかかる力
 * @param[in] dens パーティクル密度
 * @param[in] dt 時間ステップ幅
 * @param[in] nprts パーティクル数
 */
void CuSphIntegrate(float* pos, float* vel, float* frc, float* dens, int* attr, 
					float dt, uint nprts)
{
	// 1スレッド/パーティクル
	uint numThreads, numBlocks;
	computeGridSize(nprts, THREAD_NUM, numBlocks, numThreads);

	// カーネル実行
	sphIntegrate<<< numBlocks, numThreads >>>((float4*)pos, (float4*)vel, (float4*)frc, dens, attr, 
											  dt, nprts);
	
	RX_CUERROR("sphIntegrate kernel execution failed");	// カーネル実行エラーチェック
	RX_CUCHECK(cudaThreadSynchronize());		// 全てのスレッドが終わるのを待つ
}


/*!
 * パーティクル位置，速度の更新
 * @param[inout] pos パーティクル位置
 * @param[inout] vel パーティクル速度
 * @param[inout] velOld 前ステップのパーティクル速度
 * @param[in] frc パーティクルにかかる力
 * @param[in] dens パーティクル密度
 * @param[in] dt 時間ステップ幅
 * @param[in] nprts パーティクル数
 */
void CuSphIntegrateWithPolygon(float* pos, float* vel, float* frc, float* dens, int* attr, 
							   float* vrts, int* tris, int tri_num, float dt, rxParticleCell cell)
{
	// 1スレッド/パーティクル
	uint numThreads, numBlocks;
	computeGridSize(cell.uNumParticles, THREAD_NUM, numBlocks, numThreads);

	// カーネル実行
	sphIntegrateWithPolygon<<< numBlocks, numThreads >>>((float4*)pos, (float4*)vel, (float4*)frc, dens, attr, 
											   (float3*)vrts, (int3*)tris, tri_num, dt, cell);
	
	RX_CUERROR("sphIntegrate kernel execution failed");	// カーネル実行エラーチェック
	RX_CUCHECK(cudaThreadSynchronize());		// 全てのスレッドが終わるのを待つ
}


/*!
 * パーティクル位置をチェックして削除領域内ならば属性を-1にして，範囲外に飛ばす
 * @param[inout] pos パーティクル位置
 * @param[inout] vel パーティクル速度
 * @param[in] dens パーティクル密度
 * @param[in] dt 時間ステップ幅
 * @param[in] nprts パーティクル数
 */
void CuSphCheckDelete(float* pos, float* vel, int* attr, float minp[3], float maxp[3], float farpoint[3], uint nprts)
{
	// 1スレッド/パーティクル
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


	// カーネル実行
	checkDelete<<< numBlocks, numThreads >>>((float4*)pos, (float4*)vel, attr, minp3, maxp3, farpoint3, nprts);
	
	RX_CUERROR("sphIntegrate kernel execution failed");	// カーネル実行エラーチェック
	RX_CUCHECK(cudaThreadSynchronize());		// 全てのスレッドが終わるのを待つ
}

/*!
 * パーティクルのx座標値をチェックしてある値以上ならば属性を-1にして，範囲外に飛ばす
 * @param[inout] pos パーティクル位置
 * @param[inout] vel パーティクル速度
 * @param[in] dens パーティクル密度
 * @param[in] dt 時間ステップ幅
 * @param[in] nprts パーティクル数
 */
void CuSphCheckDeleteX(float* pos, float* vel, int* attr, float xmax, float farpoint[3], uint nprts)
{
	// 1スレッド/パーティクル
	uint numThreads, numBlocks;
	computeGridSize(nprts, THREAD_NUM, numBlocks, numThreads);

	float3 farpoint3;
	farpoint3.x = farpoint[0];
	farpoint3.y = farpoint[1];
	farpoint3.z = farpoint[2];

	// カーネル実行
	checkDeleteX<<< numBlocks, numThreads >>>((float4*)pos, (float4*)vel, attr, xmax, farpoint3, nprts);
	
	RX_CUERROR("sphIntegrate kernel execution failed");	// カーネル実行エラーチェック
	RX_CUCHECK(cudaThreadSynchronize());		// 全てのスレッドが終わるのを待つ
}

/*!
 * グリッド上の密度を算出
 * @param[out] dGridD グリッド上の密度値
 * @param[in] cell パーティクルグリッドデータ
 * @param[in] nx,ny グリッド数
 * @param[in] x0,y0 グリッド最小座標
 * @param[in] dx,dy グリッド幅
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

	// カーネル実行
	sphCalDensityInGrid<<<grid, threads>>>(dGridD, cell, gnum, gmin, glen);

	RX_CUERROR("Kernel execution failed");	// カーネル実行エラーチェック
	RX_CUCHECK(cudaThreadSynchronize());		// 全てのスレッドが終わるのを待つ

#if USE_TEX
	RX_CUCHECK(cudaUnbindTexture(dSortedPosTex));
	RX_CUCHECK(cudaUnbindTexture(dCellStartTex));
	RX_CUCHECK(cudaUnbindTexture(dCellEndTex));
#endif
}

/*!追加　氷用
 * グリッド上の密度を算出
 * @param[out] dGridD グリッド上の密度値
 * @param[in] cell パーティクルグリッドデータ
 * @param[in] nx,ny グリッド数
 * @param[in] x0,y0 グリッド最小座標
 * @param[in] dx,dy グリッド幅
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

	// カーネル実行
	sphCalDensityInGridIceMesh<<<grid, threads>>>(dGridD, cell, gnum, gmin, glen, bIceFlag);

	RX_CUERROR("Kernel execution failed");	// カーネル実行エラーチェック
	RX_CUCHECK(cudaThreadSynchronize());		// 全てのスレッドが終わるのを待つ

#if USE_TEX
	RX_CUCHECK(cudaUnbindTexture(dSortedPosTex));
	RX_CUCHECK(cudaUnbindTexture(dCellStartTex));
	RX_CUCHECK(cudaUnbindTexture(dCellEndTex));
#endif
}


//追加：表面パーティクルの検出
#define NEIGHT_MAX 30	//近傍最大粒子数
void CuSphDetectSurfaceParticles(int* neights, int* surface, rxParticleCell celldata, float* pos, float radius)
{
	uint numThreads, numBlocks;
	computeGridSize(celldata.uNumParticles, THREAD_NUM, numBlocks, numThreads);

	detectSurfaceParticles<<<numBlocks, numThreads>>>(neights, surface, celldata, NEIGHT_MAX, (float4*)pos, radius);

	RX_CUERROR("Kernel execution failed");	// カーネル実行エラーチェック
	cudaThreadSynchronize();
}

//追加：近傍パーティクルの検出
void CuSphDetectNeighborParticles(int* neights, rxParticleCell celldata, float radius, int prtNum)
{
	uint numThreads, numBlocks;
	computeGridSize(celldata.uNumParticles, THREAD_NUM, numBlocks, numThreads);

	detectNeighborParticles<<<numBlocks, numThreads>>>(neights, celldata, radius, NEIGHT_MAX);

	RX_CUERROR("Kernel execution failed");	// カーネル実行エラーチェック
	cudaThreadSynchronize();
}

/*!
 * 境界パーティクルの体積を計算
 *  - "Versatile Rigid-Fluid Coupling for Incompressible SPH", 2.2 式(3)の上
 * @param[out] dVolB 境界パーティクルの体積
 * @param[in]  mass パーティクル質量
 * @param[in]  cell パーティクルグリッドデータ
 */
void CuSphBoundaryVolume(float* dVolB, float mass, rxParticleCell cell)
{
	// 1スレッド/パーティクル
	uint numThreads, numBlocks;
	computeGridSize(cell.uNumParticles, THREAD_NUM, numBlocks, numThreads);

	// カーネル実行
	sphCalBoundaryVolume<<< numBlocks, numThreads >>>(dVolB, cell);

	RX_CUERROR("kernel execution failed : sphCalBoundaryVolume");	// カーネル実行エラーチェック
	RX_CUCHECK(cudaThreadSynchronize());		// 全てのスレッドが終わるのを待つ
}

/*!
 * パーティクル密度の計算(カーネル呼び出し)
 * @param[out] dDens パーティクル密度
 * @param[out] dPres パーティクル圧力
 * @param[in]  cell パーティクルグリッドデータ
 */
void CuSphBoundaryDensity(float* dDens, float* dPres, float* dPos, float* dVolB, rxParticleCell bcell, uint pnum)
{
	// 1スレッド/パーティクル
	uint numThreads, numBlocks;
	computeGridSize(pnum, THREAD_NUM, numBlocks, numThreads);

	// カーネル実行
	sphCalBoundaryDensity<<< numBlocks, numThreads >>>(dDens, dPres, (float4*)dPos, dVolB, bcell, pnum);

	RX_CUERROR("kernel execution failed : sphCalBoundaryDensity");	// カーネル実行エラーチェック
	RX_CUCHECK(cudaThreadSynchronize());		// 全てのスレッドが終わるのを待つ
}

/*!
 * パーティクルにかかる力の計算(カーネル呼び出し)
 * @param[in] dDens パーティクル密度
 * @param[in] dPres パーティクル圧力
 * @param[out] dFrc パーティクルにかかる力
 * @param[in] cell パーティクルグリッドデータ
 * @param[in] dt 時間ステップ幅
 */
void CuSphBoundaryForces(float* dDens, float* dPres, float* dPos, float* dVolB, float* dFrc, rxParticleCell bcell, uint pnum)
{
	// 1スレッド/パーティクル
	uint numThreads, numBlocks;
	computeGridSize(pnum, THREAD_NUM, numBlocks, numThreads);

	// カーネル実行
	sphCalBoundaryForce<<< numBlocks, numThreads >>>(dDens, dPres, (float4*)dPos, dVolB, (float4*)dFrc, bcell, pnum);

	RX_CUERROR("kernel execution failed : sphCalBoundaryForce");	// カーネル実行エラーチェック
	RX_CUCHECK(cudaThreadSynchronize());		// 全てのスレッドが終わるのを待つ
}

//-----------------------------------------------------------------------------
// PBDSPH
//-----------------------------------------------------------------------------

/*!
 * パーティクル密度の計算(カーネル呼び出し)
 * @param[out] dDens パーティクル密度
 * @param[out] dPres パーティクル圧力
 * @param[in]  cell パーティクルグリッドデータ
 */
void CuPbdSphDensity(float* dDens, rxParticleCell cell)
{
#if USE_TEX
	RX_CUCHECK(cudaBindTexture(0, dSortedPosTex, cell.dSortedPos, cell.uNumParticles*sizeof(float4)));
	RX_CUCHECK(cudaBindTexture(0, dCellStartTex, cell.dCellStart, cell.uNumCells*sizeof(uint)));
	RX_CUCHECK(cudaBindTexture(0, dCellEndTex, cell.dCellEnd, cell.uNumCells*sizeof(uint)));	
#endif
	//RX_CUCHECK(cudaMemset((void*)dNewDens, 0, sizeof(float2)*cell.uNumParticles));

	// 1スレッド/パーティクル
	uint numThreads, numBlocks;
	computeGridSize(cell.uNumParticles, THREAD_NUM, numBlocks, numThreads);

	// カーネル実行
	pbdsphCalDensity<<< numBlocks, numThreads >>>(dDens, cell);

	RX_CUERROR("pbdsphCalDensity kernel execution failed");	// カーネル実行エラーチェック
	RX_CUCHECK(cudaThreadSynchronize());		// 全てのスレッドが終わるのを待つ

#if USE_TEX
	RX_CUCHECK(cudaUnbindTexture(dSortedPosTex));
	RX_CUCHECK(cudaUnbindTexture(dCellStartTex));
	RX_CUCHECK(cudaUnbindTexture(dCellEndTex));
#endif
}

/*!
 * パーティクルにかかる力の計算(カーネル呼び出し)
 * @param[in] dDens パーティクル密度
 * @param[out] dFrc パーティクルにかかる力
 * @param[in] cell パーティクルグリッドデータ
 * @param[in] dt 時間ステップ幅
 */
void CuPbdSphExternalForces(float* dDens, float* dFrc, rxParticleCell cell, float dt)
{
#if USE_TEX
	RX_CUCHECK(cudaBindTexture(0, dSortedPosTex, cell.dSortedPos, cell.uNumParticles*sizeof(float4)));
	RX_CUCHECK(cudaBindTexture(0, dSortedVelTex, cell.dSortedVel, cell.uNumParticles*sizeof(float4)));
	RX_CUCHECK(cudaBindTexture(0, dCellStartTex, cell.dCellStart, cell.uNumCells*sizeof(uint)));
	RX_CUCHECK(cudaBindTexture(0, dCellEndTex, cell.dCellEnd, cell.uNumCells*sizeof(uint)));	
#endif

	// 1スレッド/パーティクル
	uint numThreads, numBlocks;
	computeGridSize(cell.uNumParticles, THREAD_NUM, numBlocks, numThreads);

	// カーネル実行
	pbdsphCalExternalForces<<< numBlocks, numThreads >>>(dDens, (float4*)dFrc, cell);

	RX_CUERROR("calForcesSPH kernel execution failed");	// カーネル実行エラーチェック
	RX_CUCHECK(cudaThreadSynchronize());		// 全てのスレッドが終わるのを待つ

#if USE_TEX
	RX_CUCHECK(cudaUnbindTexture(dSortedPosTex));
	RX_CUCHECK(cudaUnbindTexture(dSortedVelTex));
	RX_CUCHECK(cudaUnbindTexture(dCellStartTex));
	RX_CUCHECK(cudaUnbindTexture(dCellEndTex));
#endif
}

/*!
 * スケーリングファクタの計算
 * @param[in] dPos パーティクル中心座標
 * @param[out] dDens パーティクル密度
 * @param[out] dScl スケーリングファクタ
 * @param[in] eps 緩和係数
 * @param[in] cell パーティクルグリッドデータ
 */
void CuPbdSphScalingFactor(float* dPos, float* dDens, float* dScl, float eps, rxParticleCell cell)
{
#if USE_TEX
	RX_CUCHECK(cudaBindTexture(0, dSortedPosTex, cell.dSortedPos, cell.uNumParticles*sizeof(float4)));
	RX_CUCHECK(cudaBindTexture(0, dCellStartTex, cell.dCellStart, cell.uNumCells*sizeof(uint)));
	RX_CUCHECK(cudaBindTexture(0, dCellEndTex, cell.dCellEnd, cell.uNumCells*sizeof(uint)));	
#endif
	//RX_CUCHECK(cudaMemset((void*)dNewDens, 0, sizeof(float2)*cell.uNumParticles));

	// 1スレッド/パーティクル
	uint numThreads, numBlocks;
	computeGridSize(cell.uNumParticles, THREAD_NUM, numBlocks, numThreads);

	// カーネル実行
	pbdsphCalScalingFactor<<< numBlocks, numThreads >>>((float4*)dPos, dDens, dScl, eps, cell);

	RX_CUERROR("pbdsphCalScalingFactor kernel execution failed");	// カーネル実行エラーチェック
	RX_CUCHECK(cudaThreadSynchronize());		// 全てのスレッドが終わるのを待つ

#if USE_TEX
	RX_CUCHECK(cudaUnbindTexture(dSortedPosTex));
	RX_CUCHECK(cudaUnbindTexture(dCellStartTex));
	RX_CUCHECK(cudaUnbindTexture(dCellEndTex));
#endif
}

/*!
 * 平均密度変動の計算
 *  - すべてのパーティクル密度の初期密度との差をカーネルで計算し，Prefix Sum (Scan)でその合計を求める
 * @param[out] dErrScan 変動値のScan結果を格納する配列
 * @param[out] dErr パーティクル密度変動値
 * @param[in] dDens パーティクル密度
 * @param[in] rest_dens 初期密度
 * @param[in] nprts パーティクル数
 * @return 平均密度変動
 */
float CuPbdSphCalDensityFluctuation(float* dErrScan, float* dErr, float* dDens, float rest_dens, uint nprts)
{
	// 1スレッド/パーティクル
	uint numThreads, numBlocks;
	computeGridSize(nprts, THREAD_NUM, numBlocks, numThreads);

	// カーネル実行
	pbdsphDensityFluctuation<<< numBlocks, numThreads >>>(dErr, dDens, rest_dens, nprts);

	RX_CUERROR("pbdsphDensityFluctuation kernel execution failed");	// カーネル実行エラーチェック
	RX_CUCHECK(cudaThreadSynchronize());		// 全てのスレッドが終わるのを待つ

	// 各パーティクルの密度変動をScan
	CuScanf(dErrScan, dErr, nprts);

	// Exclusive scan (最後の要素が0番目からn-2番目までの合計になっている)なので，
	// Scan前配列の最後(n-1番目)の要素と合計することで密度変動の合計を計算
	float lval, lsval;
	RX_CUCHECK(cudaMemcpy((void*)&lval, (void*)(dErr+nprts-1), sizeof(float), cudaMemcpyDeviceToHost));
	RX_CUCHECK(cudaMemcpy((void*)&lsval, (void*)(dErrScan+nprts-1), sizeof(float), cudaMemcpyDeviceToHost));
	float dens_var = lval+lsval;

	return dens_var/(float)nprts;
}

/*!
 * 位置修正量の計算
 * @param[in] dPos パーティクル中心座標
 * @param[in] dScl スケーリングファクタ
 * @param[out] dDp 位置修正量
 * @param[in] cell パーティクルグリッドデータ
 */
void CuPbdSphPositionCorrection(float* dPos, float* dScl, float* dDp, rxParticleCell cell)
{
#if USE_TEX
	RX_CUCHECK(cudaBindTexture(0, dSortedPosTex, cell.dSortedPos, cell.uNumParticles*sizeof(float4)));
	RX_CUCHECK(cudaBindTexture(0, dCellStartTex, cell.dCellStart, cell.uNumCells*sizeof(uint)));
	RX_CUCHECK(cudaBindTexture(0, dCellEndTex, cell.dCellEnd, cell.uNumCells*sizeof(uint)));	
#endif
	//RX_CUCHECK(cudaMemset((void*)dNewDens, 0, sizeof(float2)*cell.uNumParticles));

	// 1スレッド/パーティクル
	uint numThreads, numBlocks;
	computeGridSize(cell.uNumParticles, THREAD_NUM, numBlocks, numThreads);

	// カーネル実行
	pbdsphPositionCorrection<<< numBlocks, numThreads >>>((float4*)dPos, dScl, (float4*)dDp, cell);

	RX_CUERROR("pbdsphPositionCorrection kernel execution failed");	// カーネル実行エラーチェック
	RX_CUCHECK(cudaThreadSynchronize());		// 全てのスレッドが終わるのを待つ

#if USE_TEX
	RX_CUCHECK(cudaUnbindTexture(dSortedPosTex));
	RX_CUCHECK(cudaUnbindTexture(dCellStartTex));
	RX_CUCHECK(cudaUnbindTexture(dCellEndTex));
#endif
}

/*!
 * パーティクル位置を更新
 * @param[inout] dPos パーティクル位置
 * @param[in] dDp 位置修正量
 * @param[in] nprts パーティクル数
 */
void CuPbdSphCorrectPosition(float* dPos, float* dDp, uint nprts)
{
	// 1スレッド/パーティクル
	uint numThreads, numBlocks;
	computeGridSize(nprts, THREAD_NUM, numBlocks, numThreads);

	// カーネル実行
	pbdsphCorrectPosition<<< numBlocks, numThreads >>>((float4*)dPos, (float4*)dDp, nprts);
	
	RX_CUERROR("pbdsphCorrectPosition kernel execution failed");	// カーネル実行エラーチェック
	RX_CUCHECK(cudaThreadSynchronize());		// 全てのスレッドが終わるのを待つ
}



/*!
 * 境界パーティクル密度を従来のパーティクル密度に加える
 * @param[inout] dDens 流体パーティクル密度
 * @param[in] dPos  流体パーティクル圧力
 * @param[in] dVolB 境界パーティクル体積
 * @param[in] bcell 境界パーティクルグリッドデータ
 */
void CuPbdSphBoundaryDensity(float* dDens, float* dPos, float* dVolB, rxParticleCell bcell, uint pnum)
{
	// 1スレッド/パーティクル
	uint numThreads, numBlocks;
	computeGridSize(pnum, THREAD_NUM, numBlocks, numThreads);

	// カーネル実行
	pbdsphCalBoundaryDensity<<< numBlocks, numThreads >>>(dDens, (float4*)dPos, dVolB, bcell, pnum);

	RX_CUERROR("kernel execution failed : sphCalBoundaryDensity");	// カーネル実行エラーチェック
	RX_CUCHECK(cudaThreadSynchronize());		// 全てのスレッドが終わるのを待つ
}


/*!
 * スケーリングファクタの計算(境界パーティクル含む)
 * @param[in] dPos 流体パーティクル中心座標
 * @param[out] dDens 流体パーティクル密度
 * @param[out] dScl 流体パーティクルのスケーリングファクタ
 * @param[in] eps 緩和係数
 * @param[in] cell 流体パーティクルグリッドデータ
 * @param[in] dVolB 境界パーティクル体積
 * @param[out] dSclB 境界パーティクルのスケーリングファクタ
 * @param[in] bcell 境界パーティクルグリッドデータ
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

	// 流体パーティクルの数だけスレッドを立てる
	uint numThreads, numBlocks;
	computeGridSize(cell.uNumParticles, THREAD_NUM, numBlocks, numThreads);

	// カーネル実行
	pbdsphCalScalingFactorWithBoundary<<< numBlocks, numThreads >>>((float4*)dPos, dDens, dScl, eps, cell, dVolB, bcell);

	RX_CUERROR("kernel execution failed : pbdsphCalScalingFactorWithBoundary");	// カーネル実行エラーチェック
	RX_CUCHECK(cudaThreadSynchronize());		// 全てのスレッドが終わるのを待つ

#if USE_TEX
	RX_CUCHECK(cudaUnbindTexture(dSortedPosTex));
	RX_CUCHECK(cudaUnbindTexture(dCellStartTex));
	RX_CUCHECK(cudaUnbindTexture(dCellEndTex));
#endif

	// 境界パーティクルのスケーリングファクタの計算
	// 境界パーティクルの数だけスレッドを立てる
	computeGridSize(bcell.uNumParticles, THREAD_NUM, numBlocks, numThreads);

	// カーネル実行
	pbdsphCalBoundaryScalingFactor<<< numBlocks, numThreads >>>((float4*)dPos, dDens, eps, cell, dVolB, dSclB, bcell);

	RX_CUERROR("kernel execution failed : pbdsphCalScalingFactorWithBoundary");	// カーネル実行エラーチェック
	RX_CUCHECK(cudaThreadSynchronize());		// 全てのスレッドが終わるのを待つ

}

/*!
 * 位置修正量の計算(境界パーティクル含む)
 * @param[in] dPos 流体パーティクル中心座標
 * @param[in] dScl 流体パーティクルのスケーリングファクタ
 * @param[out] dDens 流体パーティクル位置修正量
 * @param[in] cell 流体パーティクルグリッドデータ
 * @param[in] dVolB 境界パーティクル体積
 * @param[in] dSclB 境界パーティクルのスケーリングファクタ
 * @param[in] bcell 境界パーティクルグリッドデータ
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

	// 1スレッド/パーティクル
	uint numThreads, numBlocks;
	computeGridSize(cell.uNumParticles, THREAD_NUM, numBlocks, numThreads);

	// カーネル実行
	pbdsphPositionCorrectionWithBoundary<<< numBlocks, numThreads >>>((float4*)dPos, dScl, (float4*)dDp, cell, 
																	  dVolB, dSclB, bcell);

	RX_CUERROR("kernel execution failed");	// カーネル実行エラーチェック
	RX_CUCHECK(cudaThreadSynchronize());		// 全てのスレッドが終わるのを待つ

#if USE_TEX
	RX_CUCHECK(cudaUnbindTexture(dSortedPosTex));
	RX_CUCHECK(cudaUnbindTexture(dCellStartTex));
	RX_CUCHECK(cudaUnbindTexture(dCellEndTex));
#endif
}



/*!
 * パーティクル位置，速度の更新
 * @param[inout] pos パーティクル位置
 * @param[inout] vel パーティクル速度
 * @param[inout] velOld 前ステップのパーティクル速度
 * @param[in] frc パーティクルにかかる力
 * @param[in] dens パーティクル密度
 * @param[in] dt 時間ステップ幅
 * @param[in] nprts パーティクル数
 */
void CuPbdSphIntegrate(float* pos, float* vel, float* acc, int* attr, 
					   float* new_pos, float* new_vel, float dt, uint nprts)
{
	// 1スレッド/パーティクル
	uint numThreads, numBlocks;
	computeGridSize(nprts, THREAD_NUM, numBlocks, numThreads);

	// カーネル実行
	pbdsphIntegrate<<< numBlocks, numThreads >>>((float4*)pos, (float4*)vel, (float4*)acc, attr, 
												 (float4*)new_pos, (float4*)new_vel, dt, nprts);
	
	RX_CUERROR("pbdsphIntegrate kernel execution failed");	// カーネル実行エラーチェック
	RX_CUCHECK(cudaThreadSynchronize());		// 全てのスレッドが終わるのを待つ
}


/*!
 * パーティクル位置，速度の更新
 * @param[inout] pos パーティクル位置
 * @param[inout] vel パーティクル速度
 * @param[inout] velOld 前ステップのパーティクル速度
 * @param[in] frc パーティクルにかかる力
 * @param[in] dens パーティクル密度
 * @param[in] dt 時間ステップ幅
 * @param[in] nprts パーティクル数
 */
void CuPbdSphIntegrateWithPolygon(float* pos, float* vel, float* acc, int* attr, 
								  float* new_pos, float* new_vel, 
								  float* vrts, int* tris, int tri_num, float dt, rxParticleCell cell)
{
	// 1スレッド/パーティクル
	uint numThreads, numBlocks;
	computeGridSize(cell.uNumParticles, THREAD_NUM, numBlocks, numThreads);

	// カーネル実行
	pbdsphIntegrateWithPolygon<<< numBlocks, numThreads >>>((float4*)pos, (float4*)vel, (float4*)acc, attr, 
															(float4*)new_pos, (float4*)new_vel, (float3*)vrts, (int3*)tris, tri_num, dt, cell);
	
	RX_CUERROR("pbdsphIntegrateWithPolygon kernel execution failed");	// カーネル実行エラーチェック
	RX_CUCHECK(cudaThreadSynchronize());		// 全てのスレッドが終わるのを待つ
}


/*!
 * パーティクル位置，速度の更新
 * @param[in] pos 更新されたパーティクル位置
 * @param[inout] new_pos ステップ最初のパーティクル位置/新しいパーティクル速度
 * @param[out] new_vel 新しいパーティクル速度
 * @param[in] dt 時間ステップ幅
 * @param[in] nprts パーティクル数
 */
void CuPbdSphUpdatePosition(float* pos, float* new_pos, float* new_vel, float dt, uint nprts)
{
	// 1スレッド/パーティクル
	uint numThreads, numBlocks;
	computeGridSize(nprts, THREAD_NUM, numBlocks, numThreads);

	// カーネル実行
	pbdsphUpdatePosition<<< numBlocks, numThreads >>>((float4*)pos, (float4*)new_pos, (float4*)new_vel, dt, nprts);
	
	RX_CUERROR("CuPbdSphUpdatePosition kernel execution failed");	// カーネル実行エラーチェック
	RX_CUCHECK(cudaThreadSynchronize());		// 全てのスレッドが終わるのを待つ
}

/*!
 * パーティクル位置，速度の更新
 * @param[in] pos 更新されたパーティクル位置
 * @param[inout] new_pos ステップ最初のパーティクル位置/新しいパーティクル速度
 * @param[out] new_vel 新しいパーティクル速度
 * @param[in] dt 時間ステップ幅
 * @param[in] nprts パーティクル数
 */
void CuPbdSphUpdateVelocity(float* pos, float* new_pos, float* new_vel, float dt, uint nprts)
{
	// 1スレッド/パーティクル
	uint numThreads, numBlocks;
	computeGridSize(nprts, THREAD_NUM, numBlocks, numThreads);

	// カーネル実行
	pbdsphUpdateVelocity<<< numBlocks, numThreads >>>((float4*)pos, (float4*)new_pos, (float4*)new_vel, dt, nprts);
	
	RX_CUERROR("pbdsphUpdateVelocity kernel execution failed");	// カーネル実行エラーチェック
	RX_CUCHECK(cudaThreadSynchronize());		// 全てのスレッドが終わるのを待つ
}

/*!
 * XSPHによる粘性計算
 * @param[in] dPos パーティクル中心座標
 * @param[in] dVel パーティクル速度
 * @param[out] dNewVel 更新されたパーティクル速度
 * @param[in] c 粘性計算用パラメータ
 * @param[in] cell パーティクルグリッドデータ
 */
void CuXSphViscosity(float* dPos, float* dVel, float* dNewVel, float* dDens, float c, rxParticleCell cell)
{
#if USE_TEX
	RX_CUCHECK(cudaBindTexture(0, dSortedPosTex, cell.dSortedPos, cell.uNumParticles*sizeof(float4)));
	RX_CUCHECK(cudaBindTexture(0, dCellStartTex, cell.dCellStart, cell.uNumCells*sizeof(uint)));
	RX_CUCHECK(cudaBindTexture(0, dCellEndTex, cell.dCellEnd, cell.uNumCells*sizeof(uint)));	
#endif

	// 1スレッド/パーティクル
	uint numThreads, numBlocks;
	computeGridSize(cell.uNumParticles, THREAD_NUM, numBlocks, numThreads);

	// カーネル実行
	xsphVisocosity<<< numBlocks, numThreads >>>((float4*)dPos, (float4*)dVel, (float4*)dNewVel, dDens, c, cell);

	RX_CUERROR("pbdsphPositionCorrection kernel execution failed");	// カーネル実行エラーチェック
	RX_CUCHECK(cudaThreadSynchronize());		// 全てのスレッドが終わるのを待つ

#if USE_TEX
	RX_CUCHECK(cudaUnbindTexture(dSortedPosTex));
	RX_CUCHECK(cudaUnbindTexture(dCellStartTex));
	RX_CUCHECK(cudaUnbindTexture(dCellEndTex));
#endif
}



/*!
 * パーティクル位置をチェックして削除領域内ならば属性を-1にして，範囲外に飛ばす
 * @param[inout] pos パーティクル位置
 * @param[inout] vel パーティクル速度
 * @param[in] dens パーティクル密度
 * @param[in] dt 時間ステップ幅
 * @param[in] nprts パーティクル数
 */
void CuPbdSphCheckDelete(float* pos, float* vel, int* attr, float minp[3], float maxp[3], float farpoint[3], uint nprts)
{
	// 1スレッド/パーティクル
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


	// カーネル実行
	checkDeletePB<<< numBlocks, numThreads >>>((float4*)pos, (float4*)vel, attr, minp3, maxp3, farpoint3, nprts);
	
	RX_CUERROR("sphIntegrate kernel execution failed");	// カーネル実行エラーチェック
	RX_CUCHECK(cudaThreadSynchronize());		// 全てのスレッドが終わるのを待つ
}

/*!
 * パーティクルのx座標値をチェックしてある値以上ならば属性を-1にして，範囲外に飛ばす
 * @param[inout] pos パーティクル位置
 * @param[inout] vel パーティクル速度
 * @param[in] dens パーティクル密度
 * @param[in] dt 時間ステップ幅
 * @param[in] nprts パーティクル数
 */
void CuPbdSphCheckDeleteX(float* pos, float* vel, int* attr, float xmax, float farpoint[3], uint nprts)
{
	// 1スレッド/パーティクル
	uint numThreads, numBlocks;
	computeGridSize(nprts, THREAD_NUM, numBlocks, numThreads);

	float3 farpoint3;
	farpoint3.x = farpoint[0];
	farpoint3.y = farpoint[1];
	farpoint3.z = farpoint[2];

	// カーネル実行
	checkDeleteXPB<<< numBlocks, numThreads >>>((float4*)pos, (float4*)vel, attr, xmax, farpoint3, nprts);
	
	RX_CUERROR("sphIntegrate kernel execution failed");	// カーネル実行エラーチェック
	RX_CUCHECK(cudaThreadSynchronize());		// 全てのスレッドが終わるのを待つ
}

/*!
 * グリッド上の密度を算出
 * @param[out] dGridD グリッド上の密度値
 * @param[in] cell パーティクルグリッドデータ
 * @param[in] nx,ny グリッド数
 * @param[in] x0,y0 グリッド最小座標
 * @param[in] dx,dy グリッド幅
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

	// カーネル実行
	pbdsphCalDensityInGrid<<<grid, threads>>>(dGridD, cell, gnum, gmin, glen);

	RX_CUERROR("Kernel execution failed");	// カーネル実行エラーチェック
	RX_CUCHECK(cudaThreadSynchronize());		// 全てのスレッドが終わるのを待つ

#if USE_TEX
	RX_CUCHECK(cudaUnbindTexture(dSortedPosTex));
	RX_CUCHECK(cudaUnbindTexture(dCellStartTex));
	RX_CUCHECK(cudaUnbindTexture(dCellEndTex));
#endif
}

/*!
 * パーティクル法線の計算
 * @param[out] dNewDens パーティクル密度
 * @param[out] dNewPres パーティクル圧力
 * @param[in]  cell パーティクルグリッドデータ
 */
void CuPbdSphNormal(float* dNrms, float* dDens, rxParticleCell cell)
{
	// MRK:CuSphNormal

	// 1スレッド/パーティクル
	uint numThreads, numBlocks;
	computeGridSize(cell.uNumParticles, THREAD_NUM, numBlocks, numThreads);

	// カーネル実行
	pbdsphCalNormal<<< numBlocks, numThreads >>>((float4*)dNrms, dDens, cell);

	RX_CUERROR("sphCalNormal kernel execution failed");	// カーネル実行エラーチェック
	RX_CUCHECK(cudaThreadSynchronize());		// 全てのスレッドが終わるのを待つ

}





//-----------------------------------------------------------------------------
// Anisotropic Kernel
//-----------------------------------------------------------------------------
/*!
 * カーネル中心位置の更新と重み付き平均の計算(カーネル関数)
 * @param[out] dUpPos 更新カーネル中心
 * @param[out] dPosW 重み付き平均パーティクル座標 
 * @param[in]  lambda 平滑化のための定数
 * @param[in]  h 探索半径
 * @param[in]  cell パーティクルグリッドデータ
 */
void CuSphCalUpdatedPosition(float* dUpPos, float* dPosW, float lambda, float h, rxParticleCell cell)
{
	// MRK:CuSphCalUpdatedPosition
#if USE_TEX
	RX_CUCHECK(cudaBindTexture(0, dSortedPosTex, cell.dSortedPos, cell.uNumParticles*sizeof(float4)));
	RX_CUCHECK(cudaBindTexture(0, dCellStartTex, cell.dCellStart, cell.uNumCells*sizeof(uint)));
	RX_CUCHECK(cudaBindTexture(0, dCellEndTex, cell.dCellEnd, cell.uNumCells*sizeof(uint)));	
#endif

	// 1スレッド/パーティクル
	uint numThreads, numBlocks;
	computeGridSize(cell.uNumParticles, THREAD_NUM, numBlocks, numThreads);

	// カーネル実行
	sphCalUpdatedPosition<<< numBlocks, numThreads >>>((float4*)dUpPos, (float4*)dPosW, lambda, h, cell);

	RX_CUERROR("sphCalUpdatedPosition kernel execution failed");	// カーネル実行エラーチェック
	RX_CUCHECK(cudaThreadSynchronize());		// 全てのスレッドが終わるのを待つ

#if USE_TEX
	RX_CUCHECK(cudaUnbindTexture(dSortedPosTex));
	RX_CUCHECK(cudaUnbindTexture(dCellStartTex));
	RX_CUCHECK(cudaUnbindTexture(dCellEndTex));
#endif
}

/*!
 * 平滑化位置での重み付き平均位置の再計算とcovariance matrixの計算
 * @param[out] dPosW 重み付き平均パーティクル座標 
 * @param[out] dCMat Covariance Matrix
 * @param[in]  h 探索半径
 * @param[in]  cell パーティクルグリッドデータ
 */
void CuSphCalCovarianceMatrix(float* dPosW, float* dCMat, float h, rxParticleCell cell)
{
	// MRK:CuSphCalCovarianceMatrix
#if USE_TEX
	RX_CUCHECK(cudaBindTexture(0, dSortedPosTex, cell.dSortedPos, cell.uNumParticles*sizeof(float4)));
	RX_CUCHECK(cudaBindTexture(0, dCellStartTex, cell.dCellStart, cell.uNumCells*sizeof(uint)));
	RX_CUCHECK(cudaBindTexture(0, dCellEndTex, cell.dCellEnd, cell.uNumCells*sizeof(uint)));	
#endif

	// 1スレッド/パーティクル
	uint numThreads, numBlocks;
	computeGridSize(cell.uNumParticles, THREAD_NUM, numBlocks, numThreads);

	// カーネル実行
	sphCalCovarianceMatrix<<< numBlocks, numThreads >>>((float4*)dPosW, (matrix3x3*)dCMat, h, cell);

	RX_CUERROR("sphCalCovarianceMatrix kernel execution failed");	// カーネル実行エラーチェック
	RX_CUCHECK(cudaThreadSynchronize());		// 全てのスレッドが終わるのを待つ

#if USE_TEX
	RX_CUCHECK(cudaUnbindTexture(dSortedPosTex));
	RX_CUCHECK(cudaUnbindTexture(dCellStartTex));
	RX_CUCHECK(cudaUnbindTexture(dCellEndTex));
#endif
}

/*!
 * 特異値分解により固有値を計算
 * @param[in]  dC Covariance Matrix
 * @param[in]  dPosW 重み付き平均位置
 * @param[out] dEigen 固有値
 * @param[out] dR 固有ベクトル(回転行列)
 * @param[in]  numParticles パーティクル数
 */
void CuSphSVDecomposition(float* dC, float* dPosW, float* dEigen, float* dR, uint numParticles)
{
	// MRK:CuSphCalTransformMatrix

	// 1スレッド/パーティクル
	uint numThreads, numBlocks;
	computeGridSize(numParticles, THREAD_NUM, numBlocks, numThreads);

	// カーネル実行
	sphSVDecomposition<<< numBlocks, numThreads >>>((matrix3x3*)dC, (float4*)dPosW, (float3*)dEigen, (matrix3x3*)dR, numParticles);

	RX_CUERROR("sphCalTransformMatrix kernel execution failed");	// カーネル実行エラーチェック
	RX_CUCHECK(cudaThreadSynchronize());		// 全てのスレッドが終わるのを待つ
}

/*!
 * 固有値，固有ベクトル(回転行列)から変形行列を計算
 * @param[in]  dEigen 固有値
 * @param[in]  dR 固有ベクトル(回転行列)
 * @param[out] dG 変形行列
 * @param[in]  numParticles パーティクル数
 */
void CuSphCalTransformMatrix(float* dEigen, float* dR, float *dG, uint numParticles)
{
	// MRK:CuSphCalTransformMatrix

	// 1スレッド/パーティクル
	uint numThreads, numBlocks;
	computeGridSize(numParticles, THREAD_NUM, numBlocks, numThreads);

	// カーネル実行
	sphCalTransformMatrix<<< numBlocks, numThreads >>>((float3*)dEigen, (matrix3x3*)dR, (matrix3x3*)dG, numParticles);

	RX_CUERROR("sphCalTransformMatrix kernel execution failed");	// カーネル実行エラーチェック
	RX_CUCHECK(cudaThreadSynchronize());		// 全てのスレッドが終わるのを待つ
}


/*!
 * グリッド上の密度を算出
 * @param[out] dGridD グリッド上の密度値
 * @param[in] cell パーティクルグリッドデータ
 * @param[in] nx,ny グリッド数
 * @param[in] x0,y0 グリッド最小座標
 * @param[in] dx,dy グリッド幅
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

	// カーネル実行
	sphCalDensityAnisoInGrid<<<grid, threads>>>(dGridD, (matrix3x3*)dG, Emax, cell, gnum, gmin, glen);

	RX_CUERROR("Kernel execution failed");	// カーネル実行エラーチェック
	RX_CUCHECK(cudaThreadSynchronize());		// 全てのスレッドが終わるのを待つ

#if USE_TEX
	RX_CUCHECK(cudaUnbindTexture(dSortedPosTex));
	RX_CUCHECK(cudaUnbindTexture(dCellStartTex));
	RX_CUCHECK(cudaUnbindTexture(dCellEndTex));
#endif
}




//-----------------------------------------------------------------------------
// MARK:ウェーブレット乱流
//-----------------------------------------------------------------------------

/*!
 * パーティクル速度のエネルギースペクトラムを計算
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

	// 1スレッド/パーティクル
	uint numThreads, numBlocks;
	computeGridSize(cell.uNumParticles, THREAD_NUM, numBlocks, numThreads);

	// カーネル実行
	calParticleES<<< numBlocks, numThreads >>>((float*)dEt, h, scale, coef_et, max_et, cell, na);

	// カーネル実行エラーチェック
	RX_CUERROR("Kernel execution failed");
	RX_CUCHECK(cudaThreadSynchronize());		// 全てのスレッドが終わるのを待つ

#if USE_TEX
	RX_CUCHECK(cudaUnbindTexture(dSortedPosTex));
	RX_CUCHECK(cudaUnbindTexture(dSortedVelTex));
	RX_CUCHECK(cudaUnbindTexture(dCellStartTex));
	RX_CUCHECK(cudaUnbindTexture(dCellEndTex));
#endif
}

//! ウェーブレットノイズタイルのセット
void CuSetWaveletTile3D(float *tile, int n, int d)
{
	if(!g_dNoiseTile[d] && g_udNoiseTileSize != n){
		// デバイスメモリの確保とホストからの転送
		int size = n*n*n*sizeof(float);
		RX_CUCHECK(cudaMalloc((void**)&g_dNoiseTile[d], size));
		RX_CUCHECK(cudaMemcpy((void*)g_dNoiseTile[d], (void*)tile, size, cudaMemcpyHostToDevice));

		g_udNoiseTileSize = n;
	}
}


//! ウェーブレットノイズタイルのセット
void CuSetWaveletTile3DB(float *tile, int nx, int ny, int nz, int d)
{
	if(g_dNoiseTile[d]) RX_CUCHECK(cudaFree(g_dNoiseTile[d]));

	// デバイスメモリの確保とホストからの転送
	int size = nx*ny*nz*sizeof(float);
	RX_CUCHECK(cudaMalloc((void**)&g_dNoiseTile[d], size));
	RX_CUCHECK(cudaMemcpy((void*)g_dNoiseTile[d], (void*)tile, size, cudaMemcpyHostToDevice));

	g_uNoiseTileNum[3*d+0] = nx;
	g_uNoiseTileNum[3*d+1] = ny;
	g_uNoiseTileNum[3*d+2] = nz;
}

//! ウェーブレットノイズタイルのセット
void CuSetWaveletTile3DBT(float *tile, int nx, int ny, int nz, int d)
{
	if(g_dNoiseTile[d]) RX_CUCHECK(cudaFree(g_dNoiseTile[d]));

	// デバイスメモリの確保とホストからの転送
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
 * Wavelet turbulenceを補間速度場に追加
 * @param[out] dVturb 乱流速度場
 * @param[out] dFrc   乱流による力
 * @param[in] dEt   渦場のエネルギースペクトラム
 * @param[in] dDens パーティクル密度
 * @param[in] first,nbands 多重帯域ノイズの最初の帯域と帯域幅
 * @param[in] dPos  パーティクル座標
 * @param[in] pdim  パーティクル計算空間の大きさ
 * @param[in] nprts パーティクル数
 * @param[in] dt    時間ステップ幅
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

	// 1スレッド/パーティクル
	uint numThreads, numBlocks;
	computeGridSize(nprts, THREAD_NUM, numBlocks, numThreads);

	// カーネル実行
	calWaveletTurbulence3D<<< numBlocks, numThreads >>>((float4*)dVturb, (float4*)dFrc, (float*)dEt, (float*)dDens, 
		first, nbands, g_dNoiseTile[0], tile_n, tile_size, (float4*)dPos, pmin3, pdim3, nprts, dt);

	// カーネル実行エラーのチェック
	RX_CUERROR("Kernel execution failed");
	RX_CUCHECK(cudaThreadSynchronize());
}




//-----------------------------------------------------------------------------
// MARK:Sub Particle
//-----------------------------------------------------------------------------
/*!
 * サブ粒子のエネルギースペクトルの更新
 * @param[inout]dSubEt		サブ粒子のエネルギースペクトルの配列
 * @param[in]	dEt			元のエネルギースペクトル
 * @param[in]	scale		そのスケール
 * @param[in]	radius		レベル0での半径
 * @param[in]	maxSubLevel	サブ粒子の最大レベル
 * @param[in]	cell		パーティクルグリッドデータ
 
void CuSubUpdateEt(float* dSubEt,float* dEt,
			  float	scale,	float radius,
			  rxParticleCell cell)
{
	// 1スレッド/パーティクル
	uint numThreads, numBlocks;
	computeGridSize(cell.uNumParticles, 64, numBlocks, numThreads);

	// カーネル実行
	updateSubEt<<< numBlocks, numThreads >>>(dSubEt, dEt, scale, radius, cell.uNumParticles);

	RX_CUERROR("updateSubEt kernel execution failed");	// カーネル実行エラーチェック
	RX_CUCHECK(cudaThreadSynchronize());		// 全てのスレッドが終わるのを待つ
}*/

/*!
 * サブ粒子のEtの計算
 * @param[inout]dSubEt		サブ粒子のエネルギースペクトルへのアドレス
 * @param[in]　	dEt			粒子(レベル0)のエネルギースペクトルへのアドレス
 * @param[in]	scale		そのスケール
 * @param[in]	radius		レベル0での半径
 * @param[in]	maxSubLevel	サブ粒子の最大レベル
 * @param[in]	cell		パーティクルグリッドデータ
 */
void CuSubCalEt(	float*	dSubEt,
					float*	dEt,
					float   et_coef, 
					float	scale, 
					float	radius,
					uint	maxSubLevel,
					rxParticleCell cell)
{
	// 1スレッド/パーティクル
	uint numThreads, numBlocks;
	computeGridSize(cell.uNumParticles, THREAD_NUM, numBlocks, numThreads);

	// カーネル実行
	setSubEt<<< numBlocks, numThreads >>>(dSubEt, dEt, et_coef, scale, radius, cell.uNumParticles);
	RX_CUERROR("setSubPos kernel execution failed");	// カーネル実行エラーチェック
	RX_CUCHECK(cudaThreadSynchronize());		// 全てのスレッドが終わるのを待つ

	// カーネル実行
	for(uint level = 0; level < maxSubLevel; level++){

		uint subBlockDim = calUintPow(2,level);

		for(uint subBlockIdx = 0; subBlockIdx < subBlockDim; subBlockIdx++){

			uint	subIndexP		= cell.uNumParticles * (calUintPow(2,level) - 1 + subBlockIdx);
			uint	subIndexC		= cell.uNumParticles * (calUintPow(2,level+1) - 1 + subBlockIdx*2);

			float* dSubBlockEtC	= &dSubEt[subIndexC];
			float* dSubBlockEtP	= &dSubEt[subIndexP];

			updateSubEt<<< numBlocks, numThreads >>>(dSubBlockEtC, dSubBlockEtP, cell.uNumParticles);
		}
		RX_CUERROR("updateSubEt kernel execution failed");	// カーネル実行エラーチェック
		RX_CUCHECK(cudaThreadSynchronize());		// 全てのスレッドが終わるのを待つ
		
	}

}

/*!
 * サブ粒子の絶対座標の計算
 * @param[inout]dSubPos		サブ粒子の絶対座標へのアドレス
 * @param[in]	dPos		粒子(レベル0)の絶対座標へのアドレス
 * @param[in]	dSubChild	サブ粒子の子への単位ベクトルへのアドレス
 * @param[in]	radius		レベル0での半径
 * @param[in]	maxSubLevel	サブ粒子の最大レベル
 * @param[in]	numParticles パーティクルグリッドデータ
 */
void CuSubCalPos(	float *dSubPos,
					float *dPos,
					float *dSubChild,
					float	radius,
					uint	maxSubLevel,
					uint	nprts)
{
	// 1スレッド/パーティクル
	uint numThreads, numBlocks;
	computeGridSize(nprts, THREAD_NUM, numBlocks, numThreads);

	// カーネル実行
	setSubPos<<< numBlocks, numThreads >>>((float4*)dSubPos, (float4*)dPos, nprts);

	RX_CUERROR("setSubPos kernel execution failed");	// カーネル実行エラーチェック
	RX_CUCHECK(cudaThreadSynchronize());		// 全てのスレッドが終わるのを待つ

	
	for(uint level = 0; level < maxSubLevel; level++){

		uint subBlockDim = calUintPow(2,level);
		float radius_level = radius * pow(2.0f,-(float)level/3.0f);

		for(uint subBlockIdx = 0; subBlockIdx < subBlockDim; subBlockIdx++){

			uint	subIndexP		= nprts * (calUintPow(2,level)   - 1 + subBlockIdx);
			uint	subIndexC		= nprts * (calUintPow(2,level+1) - 1 + subBlockIdx*2);

			float4* dSubBlockPosC	= (float4*)&dSubPos[4*subIndexC];
			float4* dSubBlockPosP	= (float4*)&dSubPos[4*subIndexP];
			float4* dSubBlockChild	= (float4*)&dSubChild[4*subIndexP];
			// カーネル実行
			calSubPos<<< numBlocks, numThreads >>>(dSubBlockPosC, dSubBlockPosP, dSubBlockChild,
				0.33f*radius_level, nprts);

		}
		RX_CUERROR("calSubPos kernel execution failed");	// カーネル実行エラーチェック
		RX_CUCHECK(cudaThreadSynchronize());		// 全てのスレッドが終わるのを待つ		
	}

}

/*!
 * サブ粒子の子1への単位ベクトルと軸の更新
 * @param[out]	dSubChild	サブ粒子の子への単位ベクトルへのアドレス
 * @param[in]	dSubAxis	サブ粒子の軸へのアドレス
 * @param[in]	dSubEt		サブ粒子のエネルギースペクトルへのアドレス
 * @param[in]	radius		レベル0での半径
 * @param[in]	maxSubLevel	サブ粒子の最大レベル
 * @param[in]	cell		パーティクルグリッドデータ
 * @param[in]	dt			時間ステップ幅
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
	// 1スレッド/パーティクル
	uint numThreads, numBlocks;
	computeGridSize(cell.uNumParticles, THREAD_NUM, numBlocks, numThreads);

	// カーネル実行
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
	RX_CUERROR("updateSubChildAxis kernel execution failed");	// カーネル実行エラーチェック
	RX_CUCHECK(cudaThreadSynchronize());		// 全てのスレッドが終わるのを待つ

}

/*!
 * サブ粒子の初期化
 * @param[inout]dSubPos		サブ粒子の絶対座標へのアドレス
 * @param[in]	dPos		粒子(レベル0)の絶対座標へのアドレス
 * @param[inout]dSubAxis	サブ粒子の軸へのアドレス
 * @param[inout]dSubChild	サブ粒子の子への単位ベクトルへのアドレス
 * @param[in]	radius		レベル0での半径
 * @param[in]	maxSubLevel	サブ粒子の最大レベル
 * @param[in]	nprts		パーティクル数
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
	// 1スレッド/パーティクル
	uint numThreads, numBlocks;
	computeGridSize(nprts, THREAD_NUM, numBlocks, numThreads);

	initSubRand<<< numBlocks, numThreads >>>(dSubRand, nprts);

	// カーネル実行
	for(uint level = 0; level < maxSubLevel; level++){

		uint subBlockDim = calUintPow(2,level);

		for(uint subBlockIdx = 0; subBlockIdx < subBlockDim; subBlockIdx++){

			uint	subIndex		= nprts * (calUintPow(2,level) - 1 + subBlockIdx);
			float4* dSubBlockChild	= (float4*)(&dSubChild[4*subIndex]);
			float4* dSubBlockAxis	= (float4*)(&dSubAxis[4*subIndex]);

			// カーネル実行
			initSub<<< numBlocks, numThreads >>>(dSubBlockChild, dSubBlockAxis, dSubRand, nprts);		
		}
	}
	RX_CUERROR("initSub kernel execution failed");	// カーネル実行エラーチェック
	RX_CUCHECK(cudaThreadSynchronize());		// 全てのスレッドが終わるのを待つ

	//絶対座標を計算
	CuSubCalPos(dSubPos, dPos, dSubChild, radius, maxSubLevel, nprts);

}

//-----------------------------------------------------------------------------

/*!
 * サブ粒子の初期化
 * @param[inout]dSubPos		サブ粒子の絶対座標へのアドレス
 * @param[in]	dPos		粒子(レベル0)の絶対座標へのアドレス
 * @param[inout]dSubAxis	サブ粒子の軸へのアドレス
 * @param[inout]dSubChild	サブ粒子の子への単位ベクトルへのアドレス
 * @param[in]	radius		レベル0での半径
 * @param[in]	maxSubLevel	サブ粒子の最大レベル
 * @param[in]	nprts		パーティクル数
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
	// 1スレッド/パーティクル
	uint numThreads, numBlocks;
	computeGridSize(nprts, THREAD_NUM, numBlocks, numThreads);

	// レベル0のパーティクルに乱数の基数を設定
	initSubRand<<< numBlocks, numThreads >>>(dSubRand, nprts);
	RX_CUCHECK(cudaThreadSynchronize());

	// レベル0のパーティクル位置を設定
	setSubPos<<< numBlocks, numThreads >>>((float4*)dSubPos, (float4*)dPos, nprts);
	RX_CUCHECK(cudaThreadSynchronize());		// 全てのスレッドが終わるのを待つ

	float4 *dSubPosP, *dSubPosC1, *dSubPosC2, *dSubChildP, *dSubAxisP;

	// カーネル実行
	for(uint level = 0; level < maxSubLevel; level++){
		uint subBlockDim = calUintPow(2, level);

		for(uint subBlockIdx = 0; subBlockIdx < subBlockDim; subBlockIdx++){

			uint block_idx = subBlockDim-1+subBlockIdx;
			dSubPosP	= (float4*)(&dSubPos[4 * uMaxParticles * block_idx]);		// サブパーティクル群(L,j)
			dSubPosC1	= (float4*)(&dSubPos[4 * uMaxParticles * (2*block_idx+1)]);	// 子(L+1,0)
			dSubPosC2	= (float4*)(&dSubPos[4 * uMaxParticles * (2*block_idx+2)]);	// 子(L+1,1)
			dSubChildP	= (float4*)(&dSubChild[4 * uMaxParticles * block_idx]);		// 子(j=0)への単位ベクトル
			dSubAxisP	= (float4*)(&dSubAxis[4 * uMaxParticles * block_idx]);		// 回転軸

			// カーネル実行
			initSubParticle<<< numBlocks, numThreads >>>(dSubPosP, dSubPosC1, dSubPosC2, dSubChildP, dSubAxisP,
				dSubRand, radius, nprts, uMaxParticles);

		}
		RX_CUERROR("initSub2 kernel execution failed");	// カーネル実行エラーチェック
		RX_CUCHECK(cudaThreadSynchronize());		// 全てのスレッドが終わるのを待つ
	}

}
/*!
 * サブ粒子の初期化
 * @param[inout]dSubPos		サブ粒子の絶対座標へのアドレス
 * @param[in]	dPos		粒子(レベル0)の絶対座標へのアドレス
 * @param[inout]dSubAxis	サブ粒子の軸へのアドレス
 * @param[inout]dSubChild	サブ粒子の子への単位ベクトルへのアドレス
 * @param[in]	radius		レベル0での半径
 * @param[in]	maxSubLevel	サブ粒子の最大レベル
 * @param[in]	nprts		パーティクル数
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
	// 1スレッド/パーティクル
	uint numThreads, numBlocks;
	// FIX:ブロック数が少ないとうまく動かない？
	computeGridSize(uMaxParticles, THREAD_NUM, numBlocks, numThreads);
	//printf("block %d, thread %d\n", numBlocks, numThreads);

	// レベル0のパーティクルに乱数の基数を設定
	addSubRand<<< numBlocks, numThreads >>>(dSubRand, nprts, uStart, uCount);
	RX_CUCHECK(cudaThreadSynchronize());

	// レベル0のパーティクル位置を設定
	addSubPos<<< numBlocks, numThreads >>>((float4*)dSubPos, (float4*)dPos, nprts, uStart, uCount);
	RX_CUCHECK(cudaThreadSynchronize());		// 全てのスレッドが終わるのを待つ

	float4 *dSubPosP, *dSubPosC1, *dSubPosC2, *dSubChildP, *dSubAxisP;

	// カーネル実行
	for(uint level = 0; level < maxSubLevel; level++){
		uint subBlockDim = calUintPow(2, level);

		for(uint subBlockIdx = 0; subBlockIdx < subBlockDim; subBlockIdx++){

			uint block_idx = subBlockDim-1+subBlockIdx;
			dSubPosP	= (float4*)(&dSubPos[4 * uMaxParticles * block_idx]);		// サブパーティクル群(L,j)
			dSubPosC1	= (float4*)(&dSubPos[4 * uMaxParticles * (2*block_idx+1)]);	// 子(L+1,0)
			dSubPosC2	= (float4*)(&dSubPos[4 * uMaxParticles * (2*block_idx+2)]);	// 子(L+1,1)
			dSubChildP	= (float4*)(&dSubChild[4 * uMaxParticles * block_idx]);		// 子(j=0)への単位ベクトル
			dSubAxisP	= (float4*)(&dSubAxis[4 * uMaxParticles * block_idx]);		// 回転軸

			// カーネル実行
			addSubParticle<<< numBlocks, numThreads >>>(dSubPosP, dSubPosC1, dSubPosC2, dSubChildP, dSubAxisP, dSubRand, 
														radius, nprts, uMaxParticles, uStart, uCount);

			RX_CUERROR("addSubParticle kernel execution failed");	// カーネル実行エラーチェック
			RX_CUCHECK(cudaThreadSynchronize());		// 全てのスレッドが終わるのを待つ
		}
	}
}

/*!
 * サブ粒子の子1への単位ベクトルと軸の更新
 * @param[out]	dSubChild	サブ粒子の子への単位ベクトルへのアドレス
 * @param[in]	dSubAxis	サブ粒子の軸へのアドレス
 * @param[in]	dSubEt		サブ粒子のエネルギースペクトルへのアドレス
 * @param[in]	radius		レベル0での半径
 * @param[in]	maxSubLevel	サブ粒子の最大レベル
 * @param[in]	cell		パーティクルグリッドデータ
 * @param[in]	dt			時間ステップ幅
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
	// 1スレッド/パーティクル
	uint grid, block;	// グリッド内ブロック数,ブロック内スレッド数
	block = THREAD_NUM;
	grid = DivCeil(nprts, block);

	// レベル0のサブパーティクルの位置と乱流エネルギー値を更新
	setSubEt<<<grid, block>>>(dSubEt, dEt, et_coef, scale, radius, nprts);
	setSubPos<<<grid, block>>>((float4*)dSubPos, (float4*)dPos, nprts);

	RX_CUERROR("setSubPos kernel execution failed");	// カーネル実行エラーチェック
	RX_CUCHECK(cudaThreadSynchronize());		// 全てのスレッドが終わるのを待つ

	float4 *dSubPosP, *dSubPosC1, *dSubPosC2, *dSubChildP, *dSubAxisP;

	// カーネル実行
	for(uint level = 0; level < maxSubLevel; level++){
		uint subBlockDim = calUintPow(2,level);

		float radius_level = radius * pow(2.0f,-(float)level/3.0f);
		float ratio_et	= pow(2.0f, -(float)level*5.0f/9.0f);

		for(uint subBlockIdx = 0; subBlockIdx < subBlockDim; subBlockIdx++){

			//printf("level = %d / subBlockIdx = %d \n", level, subBlockIdx);
			uint blockIdx = subBlockDim-1+subBlockIdx;
			dSubPosP	= (float4*)(&dSubPos[4 * uMaxParticles * blockIdx]);		// サブパーティクル群(L,j)
			dSubPosC1	= (float4*)(&dSubPos[4 * uMaxParticles * (2*blockIdx+1)]);	// 子(L+1,0)
			dSubPosC2	= (float4*)(&dSubPos[4 * uMaxParticles * (2*blockIdx+2)]);	// 子(L+1,1)
			dSubChildP	= (float4*)(&dSubChild[4 * uMaxParticles * blockIdx]);		// 子(j=0)への単位ベクトル
			dSubAxisP	= (float4*)(&dSubAxis[ 4 * uMaxParticles * blockIdx]);		// 回転軸

			// (level, j)のサブパーティクルの位置を更新
			updateSubParticle<<<grid, block>>>(dSubPosP, dSubPosC1, dSubPosC2, dSubChildP, dSubAxisP,
				dSubEt, dSubRand, ratio_et, radius_level, nprts, uMaxParticles, dt);

			
		}
		RX_CUERROR("CuSubUpdate kernel execution failed");
		RX_CUCHECK(cudaThreadSynchronize());	// 全てのスレッドが終わるのを待つ
		
	}
}


/*!
 * レンダリングに使用するサブパーティクルレベル，影響係数を計算
 * @param[in] dEt レベル0のパーティクルのエネルギー値
 * @param[in] dSubPos サブパーティクルの位置
 * @param[out] subCell サブパーティクル情報
 * @param[in] et_coef レベル0のパーティクルのエネルギー値に掛ける係数
 */
void CuSubSetUnsortArray(float* dEt, float* dSubPos, rxSubParticleCell &subCell, float et_coef)
{
	// 1スレッド/親パーティクル
	uint grid, block;	// グリッド内ブロック数,ブロック内スレッド数
	block = THREAD_NUM;
	grid = DivCeil(subCell.uNumParticles, block);

	RX_CUCHECK(cudaMemset(subCell.dSubOcc, 0, subCell.uSubNumAllParticles*sizeof(uint)));

	setSubUnsortArray<<<grid, block>>>(subCell, dEt, (float4*)dSubPos, et_coef);

	RX_CUERROR("CuSubSetUnsortArray kernel execution failed");
	RX_CUCHECK(cudaThreadSynchronize());	// 全てのスレッドが終わるのを待つ

	subCell.uSubNumValidParticles = subCell.uSubNumAllParticles;


	int size = subCell.uSubNumAllParticles;

	// サブパーティクル有効/無効のScan
	thrust::exclusive_scan(thrust::device_ptr<unsigned int>(subCell.dSubOcc), 
						   thrust::device_ptr<unsigned int>(subCell.dSubOcc+size),
						   thrust::device_ptr<unsigned int>(subCell.dSubOccScan));

	// Exclusive scan (最後の要素が0番目からn-2番目までの合計になっている)なので，
	// Scan前配列の最後(n-1番目)の要素と合計することで有効サブパーティクル数を計算
	uint lval, lsval;
	RX_CUCHECK(cudaMemcpy((void*)&lval, (void*)(subCell.dSubOcc+size-1), sizeof(uint), cudaMemcpyDeviceToHost));
	RX_CUCHECK(cudaMemcpy((void*)&lsval, (void*)(subCell.dSubOccScan+size-1), sizeof(uint), cudaMemcpyDeviceToHost));
	uint num_valid_particles = lval+lsval;

	//printf("num of valid sub-particles = %d / %d\n", num_valid_particles, subCell.uSubNumAllParticles);

	// 1スレッド/全サブパーティクル
	block = THREAD_NUM;
	grid = DivCeil(subCell.uSubNumAllParticles, block);

	// 有効なパーティクルだけ詰める
	compactSubParticles<<<grid, block>>>(subCell, size);
	
	RX_CUERROR("compactSubParticles kernel execution failed");	// カーネル実行エラーチェック
	RX_CUCHECK(cudaThreadSynchronize());		// 全てのスレッドが終わるのを待つ

	RX_CUCHECK(cudaMemcpy(subCell.dSubUnsortPos, subCell.dSubSortedPos, num_valid_particles*sizeof(float4), cudaMemcpyDeviceToDevice));
	RX_CUCHECK(cudaMemcpy(subCell.dSubUnsortRad, subCell.dSubSortedRad, num_valid_particles*sizeof(float),  cudaMemcpyDeviceToDevice));
	RX_CUCHECK(cudaMemcpy(subCell.dSubUnsortRat, subCell.dSubSortedRat, num_valid_particles*sizeof(float),  cudaMemcpyDeviceToDevice));
	subCell.uSubNumValidParticles = num_valid_particles;
}


void CuSubCalcHash(rxSubParticleCell subCell)
{
	uint numThreads, numBlocks;
	computeGridSize(subCell.uSubNumValidParticles, THREAD_NUM, numBlocks, numThreads);

	// カーネル実行
	calcHashD3<<< numBlocks, numThreads >>>(subCell.dSubGridParticleHash,
										   subCell.dSubSortedIndex,
										   subCell.dSubUnsortPos,
										   subCell.uSubNumValidParticles);
	
	RX_CUERROR("Kernel execution failed");	// カーネルエラーチェック
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
 * パーティクル配列をソートされた順番に並び替え，
 * 各セルの始まりと終わりのインデックスを検索
 * @param[in] cell パーティクルグリッドデータ
 * @param[in] oldPos パーティクル位置
 * @param[in] oldVel パーティクル速度
 */
void CuSubReorderDataAndFindCellStart(rxSubParticleCell subCell)
{
	uint numThreads, numBlocks;
	computeGridSize(subCell.uSubNumValidParticles, THREAD_NUM, numBlocks, numThreads);

	RX_CUCHECK(cudaMemset(subCell.dSubCellStart, 0xffffffff, subCell.uSubNumCells*sizeof(uint)));

#if USE_TEX//テクスチャの名前に注意
	RX_CUCHECK(cudaBindTexture(0, dSubUnsortPosTex, subCell.dSubUnsortPos, subCell.uSubNumValidParticles*sizeof(float4)));
	RX_CUCHECK(cudaBindTexture(0, dSubUnsortRadTex, subCell.dSubUnsortRad, subCell.uSubNumValidParticles*sizeof(float)));
	RX_CUCHECK(cudaBindTexture(0, dSubUnsortRatTex, subCell.dSubUnsortRat, subCell.uSubNumValidParticles*sizeof(float)));
#endif

	uint smemSize = sizeof(uint)*(numThreads+1);

	
	// カーネル実行
	reorderDataAndFindCellStartF4F1F1<<< numBlocks, numThreads, smemSize>>>(subCell);

	RX_CUERROR("Kernel execution failed: CuSubReorderDataAndFindCellStart");
	RX_CUCHECK(cudaThreadSynchronize());		// 全てのスレッドが終わるのを待つ


#if USE_TEX
	RX_CUCHECK(cudaUnbindTexture(dSubUnsortPosTex));
	RX_CUCHECK(cudaUnbindTexture(dSubUnsortRadTex));
	RX_CUCHECK(cudaUnbindTexture(dSubUnsortRatTex));
#endif
}

/*!
 * グリッド上の密度を算出
 * @param[out] dGridD グリッド上の密度値
 * @param[in] cell パーティクルグリッドデータ
 * @param[in] nx,ny グリッド数
 * @param[in] x0,y0 グリッド最小座標
 * @param[in] dx,dy グリッド幅
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

	// カーネル実行
	calSubGridDensity<<<grid, threads>>>(dGridD, subCell, gnum, gmin, glen);

	RX_CUERROR("Kernel execution failed");	// カーネル実行エラーチェック
	RX_CUCHECK(cudaThreadSynchronize());		// 全てのスレッドが終わるのを待つ

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

	RX_CUERROR("Kernel execution failed");	// カーネル実行エラーチェック
	RX_CUCHECK(cudaThreadSynchronize());		// 全てのスレッドが終わるのを待つ

	RX_CUCHECK(cudaMemcpy(&hNum, dNum, sizeof(uint), cudaMemcpyDeviceToHost));

	
	//free(hNum);
	RX_CUCHECK(cudaFree(dNum));
	//printf("[CuFindNumValidHashData]Number of Valid Hash Data : %d\n",subCell.uSubNumMCParticles);
	return hNum;

}





}   // extern "C"
