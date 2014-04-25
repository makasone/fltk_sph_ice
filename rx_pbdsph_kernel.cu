/*! 
  @file rx_sph_kernel.cu
	
  @brief CUDAによるSPH
 
  @author Makoto Fujisawa
  @date 2009-08, 2011-06
*/
// FILE --rx_sph_kernel.cu--

#ifndef _RX_PBDSPH_KERNEL_CU_
#define _RX_PBDSPH_KERNEL_CU_


//-----------------------------------------------------------------------------
// インクルードファイル
//-----------------------------------------------------------------------------
#include "rx_cu_common.cu"


//-----------------------------------------------------------------------------
// PBDSPH
//-----------------------------------------------------------------------------
/*!
 * 与えられたセル内のパーティクルとの距離から密度を計算
 * @param[in] gridPos グリッド位置
 * @param[in] index パーティクルインデックス
 * @param[in] pos 計算座標
 * @param[in] cell パーティクルグリッドデータ
 * @return セル内のパーティクルから計算した密度値
 */
__device__
float calDensityCellPB(int3 gridPos, uint i, float3 pos0, rxParticleCell cell)
{
	uint gridHash = calcGridHash(gridPos);

	// セル内のパーティクルのスタートインデックス
	uint startIndex = FETCHC(dCellStart, gridHash);

	float h = params.EffectiveRadius;
	float dens = 0.0f;
	if(startIndex != 0xffffffff){	// セルが空でないかのチェック
		// セル内のパーティクルで反復
		uint endIndex = FETCHC(dCellEnd, gridHash);
		for(uint j = startIndex; j < endIndex; ++j){
			//if(j == i) continue;

			float3 pos1 = make_float3(FETCHC(dSortedPos, j));

			float3 rij = pos0-pos1;
			float r = length(rij);

			if(r <= h){
				float q = h*h-r*r;
				dens += params.Mass*params.Wpoly6*q*q*q;
			}
		}
	}

	return dens;
}



/*!
 * パーティクル密度計算(カーネル関数)
 * @param[out] newDens パーティクル密度
 * @param[out] newPres パーティクル圧力
 * @param[in]  cell パーティクルグリッドデータ
 */
__global__
void pbdsphCalDensity(float* newDens, rxParticleCell cell)
{
	// パーティクルインデックス
	uint index = __mul24(blockIdx.x,blockDim.x)+threadIdx.x;
	if(index >= cell.uNumParticles) return;	
	
	float3 pos = make_float3(FETCHC(dSortedPos, index));	// パーティクル位置
	//int3 grid_pos = calcGridPos(pos);	// パーティクルが属するグリッド位置
	float h = params.EffectiveRadius;

	// パーティクル周囲のグリッド
	int3 grid_pos0, grid_pos1;
	grid_pos0 = calcGridPos(pos-make_float3(h));
	grid_pos1 = calcGridPos(pos+make_float3(h));

	// 周囲のグリッドも含めて近傍探索，密度計算
	float dens = 0.0f;
	for(int z = grid_pos0.z; z <= grid_pos1.z; ++z){
		for(int y = grid_pos0.y; y <= grid_pos1.y; ++y){
			for(int x = grid_pos0.x; x <= grid_pos1.x; ++x){
				int3 n_grid_pos = make_int3(x, y, z);
				dens += calDensityCellPB(n_grid_pos, index, pos, cell);
			}
		}
	}

	// 密度と圧力値を結果に書き込み
	uint oIdx = cell.dSortedIndex[index];
	newDens[oIdx] = dens;
}

/*!
 * 与えられたセル内のパーティクルとの距離から力場を計算
 * @param[in] gridPos グリッド位置
 * @param[in] i パーティクルインデックス
 * @param[in] pos0 計算座標
 * @param[in] vel0 計算座標の速度
 * @param[in] dens0 計算座標の密度
 * @param[in] dens パーティクル密度
 * @param[in] cell パーティクルグリッドデータ
 * @return セル内のパーティクルから計算した力場
 */
__device__
float3 calExtForceCell(int3 gridPos, uint i, float3 pos0, float3 vel0, float dens0, float* dens, rxParticleCell cell)
{
	uint gridHash = calcGridHash(gridPos);

	// セル内のパーティクルのスタートインデックス
	uint startIndex = FETCHC(dCellStart, gridHash);

	float h = params.EffectiveRadius;

	float3 frc = make_float3(0.0f);
	if(startIndex != 0xffffffff){	// セルが空でないかのチェック
		// セル内のパーティクルで反復
		uint endIndex = FETCHC(dCellEnd, gridHash);
		for(uint j = startIndex; j < endIndex; ++j){
			if(j != i){
				// 近傍パーティクルのパラメータ
				float3 pos1 = make_float3(FETCHC(dSortedPos, j));
				float3 vel1 = make_float3(FETCHC(dSortedVel, j));

				float3 rij = pos0-pos1;
				float r = length(rij);

				if(r <= h && r > 0.0001){
					float dens1 = dens[cell.dSortedIndex[j]];

					float3 vij = vel1-vel0;

					float q = h-r;

					// 粘性項
					frc += params.Viscosity*params.Mass*(vij/dens1)*params.LWvisc*q;
				}
			}
		}
	}

	return frc;
}

/*!
 * パーティクルにかかる外力の計算(カーネル関数)
 * @param[in] dens パーティクル密度
 * @param[out] outFrc パーティクルにかかる力
 * @param[in] cell パーティクルグリッドデータ
 */
__global__
void pbdsphCalExternalForces(float* dens, float4* outFrc, rxParticleCell cell)
{
	uint index = __mul24(blockIdx.x,blockDim.x)+threadIdx.x;
	if(index >= cell.uNumParticles) return;	
	
	// ソート済み配列からパーティクルデータを取得
	float3 pos0 = make_float3(FETCHC(dSortedPos, index));
	float3 vel0 = make_float3(FETCHC(dSortedVel, index));
	float h = params.EffectiveRadius;

	// パーティクルのソートなし配列上でのインデックス
	uint oIdx = cell.dSortedIndex[index];

	float3 frc = make_float3(0.0f);
	float dens0 = dens[oIdx];

	int3 grid_pos0, grid_pos1;
	grid_pos0 = calcGridPos(pos0-make_float3(h));
	grid_pos1 = calcGridPos(pos0+make_float3(h));

	// 周囲のグリッドも含めて近傍探索，圧力項，粘性項を計算
	for(int z = grid_pos0.z; z <= grid_pos1.z; ++z){
		for(int y = grid_pos0.y; y <= grid_pos1.y; ++y){
			for(int x = grid_pos0.x; x <= grid_pos1.x; ++x){
				int3 n_grid_pos = make_int3(x, y, z);

				frc += calExtForceCell(n_grid_pos, index, pos0, vel0, dens0, dens, cell);
			}
		}
	}

	// 外力(重力)
	frc += params.Gravity;

	outFrc[oIdx] = make_float4(frc, 0.0f);
}


/*!
 * 与えられたセル内のパーティクルとの距離からスケーリングファクタの分母項計算
 * @param[in] gridPos グリッド位置
 * @param[in] index パーティクルインデックス
 * @param[in] pos 計算座標
 * @param[in] cell パーティクルグリッドデータ
 * @return セル内のパーティクルから計算した密度値
 */
__device__
float calScalingFactorCell(int3 gridPos, uint i, float3 pos0, rxParticleCell cell)
{
	uint gridHash = calcGridHash(gridPos);

	// セル内のパーティクルのスタートインデックス
	uint startIndex = FETCHC(dCellStart, gridHash);

	float h = params.EffectiveRadius;
	float r0 = params.Density;
	float sd = 0.0f;
	if(startIndex != 0xffffffff){	// セルが空でないかのチェック
		// セル内のパーティクルで反復
		uint endIndex = FETCHC(dCellEnd, gridHash);
		for(uint j = startIndex; j < endIndex; ++j){
			if(j == i) continue;

			float3 pos1 = make_float3(FETCHC(dSortedPos, j));

			float3 rij = pos0-pos1;
			float r = length(rij);

			if(r <= h && r > 0.0){
				float q = h-r;

				// Spikyカーネルで位置変動を計算
				float3 dp = (params.GWspiky*q*q*rij/r)/r0;

				sd += dot(dp, dp);
			}

		}
	}

	return sd;
}

/*!
 * スケーリングファクタの計算
 * @param[in] ppos パーティクル中心座標
 * @param[out] pdens パーティクル密度
 * @param[out] pscl スケーリングファクタ
 * @param[in] eps 緩和係数
 * @param[in]  cell パーティクルグリッドデータ
 */
__global__
void pbdsphCalScalingFactor(float4* ppos, float* pdens, float* pscl, float eps, rxParticleCell cell)
{
	// パーティクルインデックス
	uint index = __mul24(blockIdx.x,blockDim.x)+threadIdx.x;
	if(index >= cell.uNumParticles) return;	
	
	float3 pos = make_float3(FETCHC(dSortedPos, index));	// パーティクル位置
	//int3 grid_pos = calcGridPos(pos);	// パーティクルが属するグリッド位置

	float h = params.EffectiveRadius;
	float r0 = params.Density;

	// パーティクル周囲のグリッド
	int3 grid_pos0, grid_pos1;
	grid_pos0 = calcGridPos(pos-make_float3(h));
	grid_pos1 = calcGridPos(pos+make_float3(h));

	// 周囲のグリッドも含めて近傍探索，密度計算
	float dens = 0.0f;
	for(int z = grid_pos0.z; z <= grid_pos1.z; ++z){
		for(int y = grid_pos0.y; y <= grid_pos1.y; ++y){
			for(int x = grid_pos0.x; x <= grid_pos1.x; ++x){
				int3 n_grid_pos = make_int3(x, y, z);
				dens += calDensityCellPB(n_grid_pos, index, pos, cell);
			}
		}
	}

	// 密度拘束条件(式(1))
	float C = dens/r0-1.0;

	// 周囲のグリッドも含めて近傍探索，スケーリングファクタの分母項計算
	float sd = 0.0f;
	for(int z = grid_pos0.z; z <= grid_pos1.z; ++z){
		for(int y = grid_pos0.y; y <= grid_pos1.y; ++y){
			for(int x = grid_pos0.x; x <= grid_pos1.x; ++x){
				int3 n_grid_pos = make_int3(x, y, z);
				sd += calScalingFactorCell(n_grid_pos, index, pos, cell);
			}
		}
	}

	// パーティクルのソートなし配列上でのインデックス
	uint oIdx = cell.dSortedIndex[index];

	// スケーリングファクタの計算(式(11))
	pscl[oIdx] = -C/(sd+eps);

	// 更新された密度
	pdens[oIdx] = dens;
}


/*!
 * 与えられたセル内のパーティクルとの距離からスケーリングファクタの分母項計算
 * @param[in] gridPos グリッド位置
 * @param[in] index パーティクルインデックス
 * @param[in] pos 計算座標
 * @param[in] cell パーティクルグリッドデータ
 * @return セル内のパーティクルから計算した密度値
 */
__device__
float3 calPositionCorrectionCell(int3 gridPos, uint i, float3 pos0, float* pscl, rxParticleCell cell)
{
	uint gridHash = calcGridHash(gridPos);

	// セル内のパーティクルのスタートインデックス
	uint startIndex = FETCHC(dCellStart, gridHash);

	float k = params.AP_K;
	float n = params.AP_N;
	float wq = params.AP_WQ;

	float h = params.EffectiveRadius;
	float r0 = params.Density;
	float3 dp = make_float3(0.0);

	float dt = params.Dt;

	float si = pscl[cell.dSortedIndex[i]];

	if(startIndex != 0xffffffff){	// セルが空でないかのチェック
		// セル内のパーティクルで反復
		uint endIndex = FETCHC(dCellEnd, gridHash);
		for(uint j = startIndex; j < endIndex; ++j){
			if(j == i) continue;

			float3 pos1 = make_float3(FETCHC(dSortedPos, j));

			float3 rij = pos0-pos1;
			float r = length(rij);

			if(r <= h && r > 0.0){
				float scorr = 0.0f;

				if(params.AP){
					float q1 = h*h-r*r;
					float ww = params.Wpoly6*q1*q1*q1/wq;
					scorr = -k*pow(ww, n)*dt*dt;
				}
				float q = h-r;
				float sj = pscl[cell.dSortedIndex[j]];

				// Spikyカーネルで位置修正量を計算
				dp += (si+sj+scorr)*(params.GWspiky*q*q*rij/r)/r0;
			}

		}
	}

	return dp;
}

/*!
 * スケーリングファクタの計算
 * @param[in] ppos パーティクル中心座標
 * @param[out] pdens パーティクル密度
 * @param[out] pscl スケーリングファクタ
 * @param[in] eps 緩和係数
 * @param[in]  cell パーティクルグリッドデータ
 */
__global__
void pbdsphPositionCorrection(float4* ppos, float* pscl, float4* pdp, rxParticleCell cell)
{
	// パーティクルインデックス
	uint index = __mul24(blockIdx.x,blockDim.x)+threadIdx.x;
	if(index >= cell.uNumParticles) return;	
	
	float3 pos = make_float3(FETCHC(dSortedPos, index));	// パーティクル位置
	//int3 grid_pos = calcGridPos(pos);	// パーティクルが属するグリッド位置

	float h = params.EffectiveRadius;

	// パーティクル周囲のグリッド
	int3 grid_pos0, grid_pos1;
	grid_pos0 = calcGridPos(pos-make_float3(h));
	grid_pos1 = calcGridPos(pos+make_float3(h));

	// 周囲のグリッドも含めて近傍探索，位置修正量を計算
	float3 dpij = make_float3(0.0f);
	for(int z = grid_pos0.z; z <= grid_pos1.z; ++z){
		for(int y = grid_pos0.y; y <= grid_pos1.y; ++y){
			for(int x = grid_pos0.x; x <= grid_pos1.x; ++x){
				int3 n_grid_pos = make_int3(x, y, z);
				dpij += calPositionCorrectionCell(n_grid_pos, index, pos, pscl, cell);
			}
		}
	}

	// パーティクルのソートなし配列上でのインデックス
	uint oIdx = cell.dSortedIndex[index];

	// 位置修正量
	pdp[oIdx] = make_float4(dpij, 0.0);
}

/*!
 * パーティクル位置修正
 * @param[inout] pos パーティクル位置
 * @param[in] pdp 位置修正量
 * @param[in] nprts パーティクル数
 */
__global__
void pbdsphCorrectPosition(float4* ppos, float4* pdp, uint nprts)
{
	uint index = __umul24(blockIdx.x,blockDim.x)+threadIdx.x;
	if(index >= nprts) return;

	// 位置修正
	ppos[index] += pdp[index];
}

/*!
 * 密度変動の計算
 * @param[inout] pos パーティクル位置
 * @param[in] pdp 位置修正量
 * @param[in] nprts パーティクル数
 */
__global__
void pbdsphDensityFluctuation(float* perr, float* pdens, float rest_dens, uint nprts)
{
	uint index = __umul24(blockIdx.x,blockDim.x)+threadIdx.x;
	if(index >= nprts) return;

	// 密度変動
	//perr[index] = fabs(pdens[index]-rest_dens)/rest_dens;
	float err = pdens[index]-rest_dens;
	perr[index] = (err >= 0.0f ? err : 0.0f)/rest_dens;
}




/*!
 * 与えられたセル内のパーティクルとの距離から密度を計算
 * @param[in] gridPos グリッド位置
 * @param[in] index パーティクルインデックス
 * @param[in] pos 計算座標
 * @param[in] cell パーティクルグリッドデータ
 * @return セル内のパーティクルから計算した密度値
 */
__device__
float calBoundaryDensityCellPB(int3 gridPos, uint i, float3 pos0, float* dVolB, rxParticleCell bcell)
{
	uint gridHash = calcGridHashB(gridPos, params.GridSizeB);

	// セル内のパーティクルのスタートインデックス
	uint startIndex = bcell.dCellStart[gridHash];

	float h = params.EffectiveRadius;
	float dens = 0.0f;
	if(startIndex != 0xffffffff){	// セルが空でないかのチェック
		// セル内のパーティクルで反復
		uint endIndex = bcell.dCellEnd[gridHash];
		for(uint j = startIndex; j < endIndex; ++j){
			//if(j == i) continue;

			float3 pos1 = make_float3(bcell.dSortedPos[j]);
			uint jdx = bcell.dSortedIndex[j];

			float3 rij = pos0-pos1;
			float r = length(rij);

			if(r <= h){
				float q = h*h-r*r;
				dens += params.Density*dVolB[jdx]*params.Wpoly6*q*q*q;
			}
		}
	}

	return dens;
}

/*!
 * パーティクル密度計算(カーネル関数)
 * @param[out] newDens パーティクル密度
 * @param[out] newPres パーティクル圧力
 * @param[in]  cell パーティクルグリッドデータ
 */
__global__
void pbdsphCalBoundaryDensity(float* newDens, float4* dPos, float* dVolB, rxParticleCell bcell, uint pnum)
{
	// パーティクルインデックス
	uint index = __mul24(blockIdx.x,blockDim.x)+threadIdx.x;
	if(index >= pnum) return;	
	
	float3 pos = make_float3(dPos[index]);	// パーティクル位置
	float h = params.EffectiveRadius;

	// パーティクル周囲のグリッド
	int3 grid_pos0, grid_pos1;
	grid_pos0 = calcGridPosB(pos-make_float3(h), params.WorldOriginB, params.CellWidthB, params.GridSizeB);
	grid_pos1 = calcGridPosB(pos+make_float3(h), params.WorldOriginB, params.CellWidthB, params.GridSizeB);

	// 周囲のグリッドも含めて近傍探索，密度計算
	float dens = 0.0f;
	for(int z = grid_pos0.z; z <= grid_pos1.z; ++z){
		for(int y = grid_pos0.y; y <= grid_pos1.y; ++y){
			for(int x = grid_pos0.x; x <= grid_pos1.x; ++x){
				int3 n_grid_pos = make_int3(x, y, z);
				dens += calBoundaryDensityCellPB(n_grid_pos, index, pos, dVolB, bcell);
			}
		}
	}

	// 密度を結果に書き込み
	newDens[index] += dens;
}




/*!
 * 与えられたセル内のパーティクルとの距離からスケーリングファクタの分母項計算
 * @param[in] gridPos グリッド位置
 * @param[in] index パーティクルインデックス
 * @param[in] pos 計算座標
 * @param[in] cell パーティクルグリッドデータ
 * @return セル内のパーティクルから計算した密度値
 */
__device__
float calBoundaryScalingFactorCell(int3 gridPos, uint i, float3 pos0, float* dVolB, rxParticleCell bcell)
{
	uint gridHash = calcGridHashB(gridPos, params.GridSizeB);

	// セル内のパーティクルのスタートインデックス
	uint startIndex = bcell.dCellStart[gridHash];

	float h = params.EffectiveRadius;
	float r0 = params.Density;
	float sd = 0.0f;
	if(startIndex != 0xffffffff){	// セルが空でないかのチェック
		// セル内のパーティクルで反復
		uint endIndex = bcell.dCellEnd[gridHash];
		for(uint j = startIndex; j < endIndex; ++j){
			float3 pos1 = make_float3(bcell.dSortedPos[j]);
			uint jdx = bcell.dSortedIndex[j];

			float3 rij = pos0-pos1;
			float r = length(rij);

			if(r <= h && r > 0.0){
				float q = h-r;

				// Spikyカーネルで位置変動を計算
				float3 dp = (params.Density*dVolB[jdx]/params.Mass)*(params.GWspiky*q*q*rij/r)/r0;

				sd += dot(dp, dp);
			}

		}
	}

	return sd;
}

/*!
 * スケーリングファクタの計算(境界パーティクル含む)
 * @param[in] ppos パーティクル中心座標
 * @param[out] pdens パーティクル密度
 * @param[out] pscl スケーリングファクタ
 * @param[in] eps 緩和係数
 * @param[in] cell パーティクルグリッドデータ
 */
__global__
void pbdsphCalScalingFactorWithBoundary(float4* ppos, float* pdens, float* pscl, float eps, rxParticleCell cell, 
										float* bvol, rxParticleCell bcell)
{
	// パーティクルインデックス
	uint index = __mul24(blockIdx.x,blockDim.x)+threadIdx.x;
	if(index >= cell.uNumParticles) return;	
	
	float3 pos = make_float3(FETCHC(dSortedPos, index));	// パーティクル位置
	//int3 grid_pos = calcGridPos(pos);	// パーティクルが属するグリッド位置

	float h = params.EffectiveRadius;
	float r0 = params.Density;

	// パーティクル周囲のグリッド
	int3 grid_pos0, grid_pos1;
	grid_pos0 = calcGridPos(pos-make_float3(h));
	grid_pos1 = calcGridPos(pos+make_float3(h));

	// 流体パーティクルによる密度
	float dens = 0.0f;
	for(int z = grid_pos0.z; z <= grid_pos1.z; ++z){
		for(int y = grid_pos0.y; y <= grid_pos1.y; ++y){
			for(int x = grid_pos0.x; x <= grid_pos1.x; ++x){
				int3 n_grid_pos = make_int3(x, y, z);
				dens += calDensityCellPB(n_grid_pos, index, pos, cell);
			}
		}
	}

	// パーティクル周囲のグリッド
	int3 grid_pos2, grid_pos3;
	grid_pos2 = calcGridPosB(pos-make_float3(h), params.WorldOriginB, params.CellWidthB, params.GridSizeB);
	grid_pos3 = calcGridPosB(pos+make_float3(h), params.WorldOriginB, params.CellWidthB, params.GridSizeB);

	// 境界パーティクルによる密度
	for(int z = grid_pos2.z; z <= grid_pos3.z; ++z){
		for(int y = grid_pos2.y; y <= grid_pos3.y; ++y){
			for(int x = grid_pos2.x; x <= grid_pos3.x; ++x){
				int3 n_grid_pos = make_int3(x, y, z);
				dens += calBoundaryDensityCellPB(n_grid_pos, index, pos, bvol, bcell);
			}
		}
	}

	// 密度拘束条件(式(1))
	float C = dens/r0-1.0;

	// 流体パーティクルによるスケーリングファクタの分母項計算
	float sd = 0.0f;
	for(int z = grid_pos0.z; z <= grid_pos1.z; ++z){
		for(int y = grid_pos0.y; y <= grid_pos1.y; ++y){
			for(int x = grid_pos0.x; x <= grid_pos1.x; ++x){
				int3 n_grid_pos = make_int3(x, y, z);
				sd += calScalingFactorCell(n_grid_pos, index, pos, cell);
			}
		}
	}

	// 境界パーティクルによるスケーリングファクタの分母項計算
	for(int z = grid_pos2.z; z <= grid_pos3.z; ++z){
		for(int y = grid_pos2.y; y <= grid_pos3.y; ++y){
			for(int x = grid_pos2.x; x <= grid_pos3.x; ++x){
				int3 n_grid_pos = make_int3(x, y, z);
				sd += calBoundaryScalingFactorCell(n_grid_pos, index, pos, bvol, bcell);
			}
		}
	}

	// パーティクルのソートなし配列上でのインデックス
	uint oIdx = cell.dSortedIndex[index];

	// スケーリングファクタの計算(式(11))
	pscl[oIdx] = -C/(sd+eps);

	// 更新された密度
	pdens[oIdx] = dens;
}



/*!
 * スケーリングファクタの計算(境界パーティクル含む)
 * @param[in] ppos パーティクル中心座標
 * @param[out] pdens パーティクル密度
 * @param[out] pscl スケーリングファクタ
 * @param[in] eps 緩和係数
 * @param[in] cell パーティクルグリッドデータ
 */
__global__
void pbdsphCalBoundaryScalingFactor(float4* ppos, float* pdens, float eps, rxParticleCell cell, 
									float* bvol, float* bscl, rxParticleCell bcell)
{
	// パーティクルインデックス
	uint index = __mul24(blockIdx.x,blockDim.x)+threadIdx.x;
	if(index >= bcell.uNumParticles) return;	
	
	float3 pos = make_float3(bcell.dSortedPos[index]);	// パーティクル位置

	float h = params.EffectiveRadius;
	float r0 = params.Density;

	// パーティクル周囲のグリッド(流体パーティクル用)
	int3 grid_pos0, grid_pos1;
	grid_pos0 = calcGridPos(pos-make_float3(h));
	grid_pos1 = calcGridPos(pos+make_float3(h));

	// 流体パーティクルによる密度
	float dens = 0.0f;
	for(int z = grid_pos0.z; z <= grid_pos1.z; ++z){
		for(int y = grid_pos0.y; y <= grid_pos1.y; ++y){
			for(int x = grid_pos0.x; x <= grid_pos1.x; ++x){
				int3 n_grid_pos = make_int3(x, y, z);
				dens += calDensityCellPB(n_grid_pos, index, pos, cell);
			}
		}
	}

	// パーティクル周囲のグリッド(境界パーティクル用)
	int3 grid_pos2, grid_pos3;
	grid_pos2 = calcGridPosB(pos-make_float3(h), params.WorldOriginB, params.CellWidthB, params.GridSizeB);
	grid_pos3 = calcGridPosB(pos+make_float3(h), params.WorldOriginB, params.CellWidthB, params.GridSizeB);

	// 境界パーティクルによる密度
	for(int z = grid_pos2.z; z <= grid_pos3.z; ++z){
		for(int y = grid_pos2.y; y <= grid_pos3.y; ++y){
			for(int x = grid_pos2.x; x <= grid_pos3.x; ++x){
				int3 n_grid_pos = make_int3(x, y, z);
				dens += calBoundaryDensityCellPB(n_grid_pos, index, pos, bvol, bcell);
			}
		}
	}

	// 密度拘束条件(式(1))
	float C = dens/r0-1.0;

	// 流体パーティクルによるスケーリングファクタの分母項計算
	float sd = 0.0f;
	for(int z = grid_pos0.z; z <= grid_pos1.z; ++z){
		for(int y = grid_pos0.y; y <= grid_pos1.y; ++y){
			for(int x = grid_pos0.x; x <= grid_pos1.x; ++x){
				int3 n_grid_pos = make_int3(x, y, z);
				sd += calScalingFactorCell(n_grid_pos, index, pos, cell);
			}
		}
	}

	// 境界パーティクルによるスケーリングファクタの分母項計算
	for(int z = grid_pos2.z; z <= grid_pos3.z; ++z){
		for(int y = grid_pos2.y; y <= grid_pos3.y; ++y){
			for(int x = grid_pos2.x; x <= grid_pos3.x; ++x){
				int3 n_grid_pos = make_int3(x, y, z);
				sd += calBoundaryScalingFactorCell(n_grid_pos, index, pos, bvol, bcell);
			}
		}
	}

	// パーティクルのソートなし配列上でのインデックス
	uint oIdx = bcell.dSortedIndex[index];

	// スケーリングファクタの計算(式(11))
	bscl[oIdx] = -C/(sd+eps);
}



/*!
 * 与えられたセル内のパーティクルとの距離からスケーリングファクタの分母項計算
 * @param[in] gridPos グリッド位置
 * @param[in] index パーティクルインデックス
 * @param[in] pos 計算座標
 * @param[in] cell パーティクルグリッドデータ
 * @return セル内のパーティクルから計算した密度値
 */
__device__
float3 calBoundaryPositionCorrectionCell(int3 gridPos, uint i, float3 pos0, float si, float* bscl, float* bvol, rxParticleCell bcell)
{
	uint gridHash = calcGridHashB(gridPos, params.GridSizeB);

	// セル内のパーティクルのスタートインデックス
	uint startIndex = bcell.dCellStart[gridHash];

	float k = params.AP_K;
	float n = params.AP_N;
	float wq = params.AP_WQ;

	float h = params.EffectiveRadius;
	float r0 = params.Density;
	float3 dp = make_float3(0.0);

	float dt = params.Dt;

	if(startIndex != 0xffffffff){	// セルが空でないかのチェック
		// セル内のパーティクルで反復
		uint endIndex = bcell.dCellEnd[gridHash];
		for(uint j = startIndex; j < endIndex; ++j){
			float3 pos1 = make_float3(bcell.dSortedPos[j]);
			uint jdx = bcell.dSortedIndex[j];

			float3 rij = pos0-pos1;
			float r = length(rij);

			if(r <= h && r > 0.0){
				float scorr = 0.0f;

				if(params.AP){
					float q1 = h*h-r*r;
					float ww = (params.Density*bvol[jdx]/params.Mass)*params.Wpoly6*q1*q1*q1/wq;
					scorr = -k*pow(ww, n)*dt*dt;
				}
				float q = h-r;
				float sj = bscl[jdx];

				// Spikyカーネルで位置修正量を計算
				dp += (si+sj+scorr)*(params.GWspiky*q*q*rij/r)/r0;
			}

		}
	}

	return dp;
}

/*!
 * スケーリングファクタの計算
 * @param[in] ppos パーティクル中心座標
 * @param[out] pdens パーティクル密度
 * @param[out] pscl スケーリングファクタ
 * @param[in] eps 緩和係数
 * @param[in]  cell パーティクルグリッドデータ
 */
__global__
void pbdsphPositionCorrectionWithBoundary(float4* ppos, float* pscl, float4* pdp, rxParticleCell cell, 
										  float* bvol, float* bscl, rxParticleCell bcell)
{
	// パーティクルインデックス
	uint index = __mul24(blockIdx.x,blockDim.x)+threadIdx.x;
	if(index >= cell.uNumParticles) return;	
	
	float3 pos = make_float3(FETCHC(dSortedPos, index));	// パーティクル位置
	//int3 grid_pos = calcGridPos(pos);	// パーティクルが属するグリッド位置

	float h = params.EffectiveRadius;

	float si = pscl[cell.dSortedIndex[index]];


	// パーティクル周囲のグリッド(流体パーティクル用)
	int3 grid_pos0, grid_pos1;
	grid_pos0 = calcGridPos(pos-make_float3(h));
	grid_pos1 = calcGridPos(pos+make_float3(h));

	// 流体パーティクルによる位置修正量を計算
	float3 dpij = make_float3(0.0f);
	for(int z = grid_pos0.z; z <= grid_pos1.z; ++z){
		for(int y = grid_pos0.y; y <= grid_pos1.y; ++y){
			for(int x = grid_pos0.x; x <= grid_pos1.x; ++x){
				int3 n_grid_pos = make_int3(x, y, z);
				dpij += calPositionCorrectionCell(n_grid_pos, index, pos, pscl, cell);
			}
		}
	}

	// パーティクル周囲のグリッド(境界パーティクル用)
	int3 grid_pos2, grid_pos3;
	grid_pos2 = calcGridPosB(pos-make_float3(h), params.WorldOriginB, params.CellWidthB, params.GridSizeB);
	grid_pos3 = calcGridPosB(pos+make_float3(h), params.WorldOriginB, params.CellWidthB, params.GridSizeB);

	// 境界パーティクルによる位置修正量を計算
	for(int z = grid_pos2.z; z <= grid_pos3.z; ++z){
		for(int y = grid_pos2.y; y <= grid_pos3.y; ++y){
			for(int x = grid_pos2.x; x <= grid_pos3.x; ++x){
				int3 n_grid_pos = make_int3(x, y, z);
				dpij += calBoundaryPositionCorrectionCell(n_grid_pos, index, pos, si, bscl, bvol, bcell);
			}
		}
	}

	// パーティクルのソートなし配列上でのインデックス
	uint oIdx = cell.dSortedIndex[index];

	// 位置修正量
	pdp[oIdx] = make_float4(dpij, 0.0);
}


__device__
void calCollisionSolidPB(float3 &pos, float3 &vel, float dt)
{
	float d;
	float3 n;
	float3 cp;

	// ボックス形状のオブジェクトとの衝突
#if MAX_BOX_NUM
	for(int i = 0; i < params.BoxNum; ++i){
		if(params.BoxFlg[i] == 0) continue;
		
		collisionPointBox(pos, params.BoxCen[i], params.BoxExt[i], params.BoxRot[i], params.BoxInvRot[i], cp, d, n);

		if(d < 0.0){
			float res = params.Restitution;
			res = (res > 0) ? (res*fabs(d)/(dt*length(vel))) : 0.0f;
			vel -= (1+res)*n*dot(n, vel);
			pos = cp;
		}
	}
#endif

	// 球形状のオブジェクトとの衝突
#if MAX_SPHERE_NUM
	for(int i = 0; i < params.SphereNum; ++i){
		if(params.SphereFlg[i] == 0) continue;

		collisionPointSphere(pos, params.SphereCen[i], params.SphereRad[i], cp, d, n);

		if(d < 0.0){
			float res = params.Restitution;
			res = (res > 0) ? (res*fabs(d)/(dt*length(vel))) : 0.0f;
			vel -= (1+res)*n*dot(n, vel);
			pos = cp;
		}
	}
#endif

	// 周囲の境界との衝突判定
	float3 l0 = params.BoundaryMin;
	float3 l1 = params.BoundaryMax;
	collisionPointAABB(pos, 0.5*(l1+l0), 0.5*(l1-l0), cp, d, n);

	if(d < 0.0){
		float res = params.Restitution;
		res = (res > 0) ? (res*fabs(d)/(dt*length(vel))) : 0.0f;
		vel -= (1+res)*n*dot(n, vel);
		pos = cp;
	}
}

__device__
inline bool calCollisionPolygonPB(float3 &pos0, float3 &pos1, float3 &vel, float3 v0, float3 v1, float3 v2, float dt)
{
	float3 cp, n;
	if(intersectSegmentTriangle(pos0, pos1, v0, v1, v2, cp, n, params.ParticleRadius) == 1){
		float d = length(pos1-cp);
		n = normalize(n);

		float res = params.Restitution;
		res = (res > 0) ? (res*fabs(d)/(dt*length(vel))) : 0.0f;
		float3 vr = -(1+res)*n*dot(n, vel);

		float l = length(pos1-pos0);
		pos1 = cp+vr*(dt*d/l);
		vel += vr;//+params.PolyVel[0];
		//vel.x = 1.0;

		return true;
	}
	return false;
}



/*!
 * パーティクル位置，速度の更新
 * @param[inout] ppos パーティクル位置
 * @param[inout] pvel パーティクル速度
 * @param[in] pfrc パーティクルにかかる力
 * @param[in] dens パーティクル密度
 * @param[in] dt 時間ステップ幅
 * @param[in] nprts パーティクル数
 */
__global__
void pbdsphIntegrate(float4* ppos, float4* pvel, float4* pacc, int* attr, 
					 float4* new_ppos, float4* new_pvel, float dt, uint nprts)
{
	uint index = __umul24(blockIdx.x,blockDim.x)+threadIdx.x;
	if(index >= nprts) return;
	//if(attr[index] == -1) return;

	float3 x = make_float3(ppos[index]);
	float3 v = make_float3(pvel[index]);
	float3 a = make_float3(pacc[index]);
	//float3 v_old = v;

	// 更新位置，速度の更新
	v += dt*a;
	x += dt*v;

	// 固体・境界との衝突
	calCollisionSolidPB(x, v, dt);

	// 位置と速度の更新
	new_ppos[index] = make_float4(x);
	new_pvel[index] = make_float4(v);
}



/*!
 * パーティクル位置，速度の更新(Leap-Frog)
 * @param[inout] ppos パーティクル位置
 * @param[inout] pvel パーティクル速度
 * @param[in] pfrc パーティクルにかかる力
 * @param[in] dens パーティクル密度
 * @param[in] vrts
 * @param[in] dt 時間ステップ幅
 * @param[in] nprts パーティクル数
 */
__global__
void pbdsphIntegrateWithPolygon(float4* ppos, float4* pvel, float4* pacc, int* attr, 
								float4* new_ppos, float4* new_pvel, 
								float3* vrts, int3* tris, int tri_num, float dt, rxParticleCell cell)
{
	uint index = __umul24(blockIdx.x,blockDim.x)+threadIdx.x;
	if(index >= cell.uNumParticles) return;
	//if(attr[index] == -1) return;

	float3 x = make_float3(ppos[index]);
	float3 v = make_float3(pvel[index]);
	float3 a = make_float3(pacc[index]);
	//float3 v_old = v;
	float3 x_old = x;

	// 更新位置，速度の更新
	v += dt*a;
	x += dt*v;

	// ポリゴンオブジェクトとの衝突
	int3 gridPos[2];
	gridPos[0] = calcGridPos(x_old);	// 位置更新前のパーティクルが属するグリッド
	gridPos[1] = calcGridPos(x);		// 位置更新後のパーティクルが属するグリッド
	for(int i = 0; i < 2; ++i){
		uint grid_hash = calcGridHash(gridPos[i]);
		uint start_index = cell.dPolyCellStart[grid_hash];
		if(start_index != 0xffffffff){	// セルが空でないかのチェック

			uint end_index = cell.dPolyCellEnd[grid_hash];
			for(uint j = start_index; j < end_index; ++j){
				uint pidx = cell.dSortedPolyIdx[j];

				int3 idx = tris[pidx];
				if(calCollisionPolygonPB(x_old, x, v, vrts[idx.x], vrts[idx.y], vrts[idx.z], dt)){
				}
			}
		}
	}

	// 固体・境界との衝突
	calCollisionSolidPB(x, v, dt);

	// 位置と速度の更新
	new_ppos[index] = make_float4(x);
	new_pvel[index] = make_float4(v);
}




/*!
 * パーティクル位置，速度の更新
 * @param[in] ppos 更新されたパーティクル位置
 * @param[inout] new_ppos ステップ最初のパーティクル位置/新しいパーティクル速度
 * @param[out] new_pvel 新しいパーティクル速度
 * @param[in] dt 時間ステップ幅
 * @param[in] nprts パーティクル数
 */
__global__
void pbdsphUpdatePosition(float4* ppos, float4* new_ppos, float4* new_pvel, float dt, uint nprts)
{
	uint index = __umul24(blockIdx.x,blockDim.x)+threadIdx.x;
	if(index >= nprts) return;

	float3 x0 = make_float3(new_ppos[index]);
	float3 x1 = make_float3(ppos[index]);
	float3 v = (x1-x0)/dt;

	// 位置と速度の更新
	new_pvel[index] = make_float4(v);
	new_ppos[index] = make_float4(x1);
}

/*!
 * パーティクル速度の更新
 * @param[in] ppos 更新されたパーティクル位置
 * @param[in] new_ppos ステップ最初のパーティクル位置/新しいパーティクル速度
 * @param[out] new_pvel 新しいパーティクル速度
 * @param[in] dt 時間ステップ幅
 * @param[in] nprts パーティクル数
 */
__global__
void pbdsphUpdateVelocity(float4* ppos, float4* new_ppos, float4* new_pvel, float dt, uint nprts)
{
	uint index = __umul24(blockIdx.x,blockDim.x)+threadIdx.x;
	if(index >= nprts) return;

	float3 x0 = make_float3(new_ppos[index]);
	float3 x1 = make_float3(ppos[index]);
	float3 v = (x1-x0)/dt;

	// 位置と速度の更新
	new_pvel[index] = make_float4(v);
}



/*!
 * 与えられたセル内のパーティクルとの距離から密度を計算
 * @param[in] gridPos グリッド位置
 * @param[in] index パーティクルインデックス
 * @param[in] pos 計算座標
 * @param[in] cell パーティクルグリッドデータ
 * @return セル内のパーティクルから計算した密度値
 */
__device__
float3 calXsphViscosityCell(int3 gridPos, uint i, float3 pos0, float3 vel0, float4* pvel, float* dens, rxParticleCell cell)
{
	uint gridHash = calcGridHash(gridPos);

	// セル内のパーティクルのスタートインデックス
	uint startIndex = FETCHC(dCellStart, gridHash);

	float h = params.EffectiveRadius;
	float3 v = make_float3(0.0);
	if(startIndex != 0xffffffff){	// セルが空でないかのチェック
		// セル内のパーティクルで反復
		uint endIndex = FETCHC(dCellEnd, gridHash);
		for(uint j = startIndex; j < endIndex; ++j){
			//if(j == i) continue;

			float3 pos1 = make_float3(FETCHC(dSortedPos, j));

			float3 rij = pos0-pos1;
			float r = length(rij);

			if(r <= h){
				float3 vel1 = make_float3(pvel[cell.dSortedIndex[j]]);
				float3 rho1 = make_float3(dens[cell.dSortedIndex[j]]);

				float q = h*h-r*r;
				v += (params.Mass/rho1)*(vel1-vel0)*params.Wpoly6*q*q*q;
			}
		}
	}

	return v;
}

/*!
 * パーティクル密度計算(カーネル関数)
 * @param[out] newDens パーティクル密度
 * @param[out] newPres パーティクル圧力
 * @param[in]  cell パーティクルグリッドデータ
 */
__global__
void xsphVisocosity(float4* ppos, float4* pvel, float4* new_pvel, float* dens, float c, rxParticleCell cell)
{
	// パーティクルインデックス
	uint index = __mul24(blockIdx.x,blockDim.x)+threadIdx.x;
	if(index >= cell.uNumParticles) return;	
	
	float3 pos0 = make_float3(FETCHC(dSortedPos, index));	// パーティクル位置
	float3 vel0 = make_float3(pvel[cell.dSortedIndex[index]]);	// パーティクル速度
	//int3 grid_pos = calcGridPos(pos0);	// パーティクルが属するグリッド位置
	float h = params.EffectiveRadius;

	// パーティクル周囲のグリッド
	int3 grid_pos0, grid_pos1;
	grid_pos0 = calcGridPos(pos0-make_float3(h));
	grid_pos1 = calcGridPos(pos0+make_float3(h));

	// 周囲のグリッドも含めて近傍探索，密度計算
	float3 v = make_float3(0.0);
	for(int z = grid_pos0.z; z <= grid_pos1.z; ++z){
		for(int y = grid_pos0.y; y <= grid_pos1.y; ++y){
			for(int x = grid_pos0.x; x <= grid_pos1.x; ++x){
				int3 n_grid_pos = make_int3(x, y, z);
				v += calXsphViscosityCell(n_grid_pos, index, pos0, vel0, pvel, dens, cell);
			}
		}
	}

	// 密度と圧力値を結果に書き込み
	uint oIdx = cell.dSortedIndex[index];
	new_pvel[oIdx] = make_float4(vel0+c*v);
	//new_pvel[oIdx] = make_float4(vel0);
}


/*!
 * パーティクル位置をチェックして削除領域内ならば属性を-1にして，範囲外に飛ばす
 * @param[inout] ppos パーティクル位置
 * @param[in] dt 時間ステップ幅
 * @param[in] nprts パーティクル数
 */
__global__
void checkDeletePB(float4* ppos, float4* pvel, int* attr, float3 minp, float3 maxp, float3 far_pos, uint nprts)
{
	uint index = __umul24(blockIdx.x,blockDim.x)+threadIdx.x;
	if(index >= nprts) return;
	if(attr[index] == -1) return;	// すでに削除されていた場合は飛ばす

	float3 x = make_float3(ppos[index]);

	if((x.x > minp.x && x.x < maxp.x) && (x.y > minp.y && x.y < maxp.y) && (x.z > minp.z && x.z < maxp.z)){
		// 属性を-1にして，far_pointに飛ばす
		ppos[index] = make_float4(far_pos);
		pvel[index] = make_float4(0.0);
		attr[index] = -1;
	}
}
/*!
 * パーティクル位置をチェックして削除領域内ならば属性を-1にして，範囲外に飛ばす
 * @param[inout] ppos パーティクル位置
 * @param[in] dt 時間ステップ幅
 * @param[in] nprts パーティクル数
 */
__global__
void checkDeleteXPB(float4* ppos, float4* pvel, int* attr, float xmax, float3 far_pos, uint nprts)
{
	uint index = __umul24(blockIdx.x,blockDim.x)+threadIdx.x;
	if(index >= nprts) return;
	if(attr[index] == -1) return;	// すでに削除されていた場合は飛ばす

	float3 x = make_float3(ppos[index]);

	if(x.x > xmax){
		// 属性を-1にして，far_pointに飛ばす
		ppos[index] = make_float4(far_pos);
		pvel[index] = make_float4(0.0);
		attr[index] = -1;
	}
}




/*!
 * 与えられたセル内のパーティクルとの距離から密度を計算
 * @param[in] gridPos グリッド位置
 * @param[in] pos 計算座標
 * @param[in] cell パーティクルグリッドデータ
 * @return セル内のパーティクルから計算した密度値
 */
__device__
float calDensityCellGPB(int3 gridPos, float3 pos0, rxParticleCell cell)
{
	uint gridHash = calcGridHash(gridPos);

	// セル内のパーティクルのスタートインデックス
	uint startIndex = FETCHC(dCellStart, gridHash);

	float h = params.EffectiveRadius;
	float d = 0.0f;
	if(startIndex != 0xffffffff){	// セルが空でないかのチェック
		// セル内のパーティクルで反復
		uint endIndex = FETCHC(dCellEnd, gridHash);

		for(uint j = startIndex; j < endIndex; ++j){
			//if(j != index){
				float3 pos1 = make_float3(FETCHC(dSortedPos, j));

				float3 rij = pos0-pos1;
				float r = length(rij);

				if(r <= h){
					float q = h*h-r*r;

					d += params.Mass*params.Wpoly6*q*q*q;
				}

			//}
		}
	}

	return d;
}

/*!
 * グリッド上での密度を計算
 * @param[out] GridD グリッド密度
 * @param[in] cell パーティクルグリッドデータ
 * @param[in] gnum グリッド数
 * @param[in] gmin グリッド最小座標
 * @param[in] glen グリッド幅
 */
__global__
void pbdsphCalDensityInGrid(float* GridD, rxParticleCell cell, 
					uint3 gnum, float3 gmin, float3 glen)
{
	uint blockId = __mul24(blockIdx.y, gridDim.x)+blockIdx.x;
	uint i = __mul24(blockId, blockDim.x)+threadIdx.x;

	uint3 gridPos = calcGridPosU(i, gnum);

	if(gridPos.x < gnum.x && gridPos.y < gnum.y && gridPos.z < gnum.z){
		float3 gpos;
		gpos.x = gmin.x+(gridPos.x)*glen.x;
		gpos.y = gmin.y+(gridPos.y)*glen.y;
		gpos.z = gmin.z+(gridPos.z)*glen.z;

		float d = 0.0f;

		int3 pgpos = calcGridPos(gpos);

		float h = params.EffectiveRadius;
		int3 grid_pos0, grid_pos1;
		grid_pos0 = calcGridPos(gpos-make_float3(h));
		grid_pos1 = calcGridPos(gpos+make_float3(h));

		for(int z = grid_pos0.z; z <= grid_pos1.z; ++z){
			for(int y = grid_pos0.y; y <= grid_pos1.y; ++y){
				for(int x = grid_pos0.x; x <= grid_pos1.x; ++x){
					int3 neighbourPos = make_int3(x, y, z);

					d += calDensityCellGPB(neighbourPos, gpos, cell);
				}
			}
		}

		GridD[gridPos.x+gridPos.y*gnum.x+gridPos.z*gnum.x*gnum.y] = d;
	}

}

/*!
 * 与えられたセル内のパーティクルとの距離から法線を計算
 * @param[in] gridPos グリッド位置
 * @param[in] i パーティクルインデックス
 * @param[in] pos 計算座標
 * @param[in] cell パーティクルグリッドデータ
 * @return セル内のパーティクルから計算した密度値
 */
__device__
float3 calNormalCellPB(int3 gridPos, uint i, float3 pos0, float* dens, rxParticleCell cell)
{
	uint gridHash = calcGridHash(gridPos);

	// セル内のパーティクルのスタートインデックス
	uint startIndex = FETCHC(dCellStart, gridHash);

	float h = params.EffectiveRadius;
	float3 nrm = make_float3(0.0f);
	if(startIndex != 0xffffffff){	// セルが空でないかのチェック
		// セル内のパーティクルで反復
		uint endIndex = FETCHC(dCellEnd, gridHash);

		for(uint j = startIndex; j < endIndex; ++j){
			if(j != i){
				float3 pos1 = make_float3(FETCHC(dSortedPos, j));

				float3 rij = pos0-pos1;
				float r = length(rij);

				if(r <= h && r > 0.0001){
					float d1 = dens[cell.dSortedIndex[j]];
					float q = h*h-r*r;

					nrm += (params.Mass/d1)*params.GWpoly6*q*q*rij;
				}

			}
		}
	}

	return nrm;
}


/*!
 * パーティクル法線計算(カーネル関数)
 * @param[out] newNrms パーティクル法線
 * @param[in] dens パーティクル密度
 * @param[in] cell パーティクルグリッドデータ
 */
__global__
void pbdsphCalNormal(float4* newNrms, float* dens, rxParticleCell cell)
{
	uint index = __mul24(blockIdx.x,blockDim.x)+threadIdx.x;
	if(index >= cell.uNumParticles) return;	
	
	float3 pos = make_float3(FETCHC(dSortedPos, index));	// パーティクル位置
	float h = params.EffectiveRadius;
	//int3 grid_pos = calcGridPos(pos);	// パーティクルが属するグリッド位置

	// パーティクル周囲のグリッド
	int3 grid_pos0, grid_pos1;
	grid_pos0 = calcGridPos(pos-make_float3(h));
	grid_pos1 = calcGridPos(pos+make_float3(h));

	// 周囲のグリッドも含めて近傍探索，密度計算
	float3 nrm = make_float3(0.0f);
	for(int z = grid_pos0.z; z <= grid_pos1.z; ++z){
		for(int y = grid_pos0.y; y <= grid_pos1.y; ++y){
			for(int x = grid_pos0.x; x <= grid_pos1.x; ++x){
				int3 n_grid_pos = make_int3(x, y, z);
				nrm += calNormalCellPB(n_grid_pos, index, pos, dens, cell);
			}
		}
	}

	float l = length(nrm);
	if(l > 0){
		nrm /= l;
	}

	uint oIdx = cell.dSortedIndex[index];
	newNrms[oIdx] = make_float4(nrm, 0.0f);
}





#endif // #ifndef _RX_PBDSPH_KERNEL_CU_



