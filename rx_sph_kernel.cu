/*! 
  @file rx_sph_kernel.cu
	
  @brief CUDAによるSPH
 
  @author Makoto Fujisawa
  @date 2009-08, 2011-06
*/
// FILE --rx_sph_kernel.cu--

#ifndef _RX_CUSPH_KERNEL_CU_
#define _RX_CUSPH_KERNEL_CU_


//-----------------------------------------------------------------------------
// インクルードファイル
//-----------------------------------------------------------------------------
#include "rx_cu_common.cu"


//-----------------------------------------------------------------------------
// ハッシュ
//-----------------------------------------------------------------------------
/*!
 * 各パーティクルのグリッドハッシュ値
 * @param[out] gridParticleHash
 * @param[out] dSortedIndex
 * @param[in] pos パーティクル位置を格納した配列
 * @param[in] nprts パーティクル数
 */
__global__
void calcHashD(uint*   dGridParticleHash, 
			   uint*   dSortedIndex, 
			   float4* dPos, 
			   uint	   nprts)
{
	uint index = __umul24(blockIdx.x, blockDim.x)+threadIdx.x;
	if(index >= nprts) return;
	
	volatile float4 p = dPos[index];
	int3 gridPos = calcGridPos(make_float3(p.x, p.y, p.z));
	uint hash = calcGridHash(gridPos);

	dGridParticleHash[index] = hash;
	dSortedIndex[index] = index;
}
/*!
 * 各パーティクルのグリッドハッシュ値
 * @param[out] gridParticleHash
 * @param[out] dSortedIndex
 * @param[in] pos パーティクル位置を格納した配列
 * @param[in] nprts パーティクル数
 */
__global__
void calcHashD2(uint*   dGridParticleHash, 
				uint*   dSortedIndex, 
				float4* dPos, 
				int*    dAttr, 
				uint	   nprts)
{
	uint index = __umul24(blockIdx.x, blockDim.x)+threadIdx.x;
	if(index >= nprts) return;
	
	uint hash = 2147483647;
	if(dAttr[index] != -1){
		volatile float4 p = dPos[index];
		int3 gridPos = calcGridPos(make_float3(p.x, p.y, p.z));
		hash = calcGridHash(gridPos);
	}

	dGridParticleHash[index] = hash;
	dSortedIndex[index] = index;
}
/*!
 * 各パーティクルのグリッドハッシュ値
 * @param[out] gridParticleHash
 * @param[out] dSortedIndex
 * @param[in] pos パーティクル位置を格納した配列
 * @param[in] nprts パーティクル数
 */
__global__
void calcHashB(uint*   dGridParticleHash, 
			   uint*   dSortedIndex, 
			   float4*  dPos, 
			   float3  world_origin, 
			   float3  cell_width, 
			   uint3   grid_size, 
			   uint	   nprts)
{
	uint index = __umul24(blockIdx.x, blockDim.x)+threadIdx.x;
	if(index >= nprts) return;
	
	float3 p = make_float3(dPos[index]);
	int3 gridPos = calcGridPosB(make_float3(p.x, p.y, p.z), world_origin, cell_width, grid_size);
	uint hash = calcGridHashB(gridPos, grid_size);

	dGridParticleHash[index] = hash;
	dSortedIndex[index] = index;
}


/*!
 * パーティクルデータをソートして，ハッシュ内の各セルの最初のアドレスを検索
 * @param[in] cell パーティクルグリッドデータ
 * @param[in] oldPos パーティクル位置
 * @param[in] oldVel パーティクル速度
 */
__global__
void reorderDataAndFindCellStartD(rxParticleCell cell, float4* dSortedPos, float4* dSortedVel)
{
	extern __shared__ uint sharedHash[];	// サイズ : blockSize+1
	uint index = __umul24(blockIdx.x,blockDim.x)+threadIdx.x;
	
	uint hash;
	if(index < cell.uNumParticles){
		hash = cell.dGridParticleHash[index];	// ハッシュ値

		sharedHash[threadIdx.x+1] = hash;	// ハッシュ値をシェアードメモリに格納

		if(index > 0 && threadIdx.x == 0){
			// 各シェアードメモリの最初は隣のグリッドのパーティクルのハッシュ値を格納
			sharedHash[0] = cell.dGridParticleHash[index-1];
		}
	}

	__syncthreads();
	
	if(index < cell.uNumParticles){
		// インデックス0である，もしくは，一つ前のパーティクルのグリッドハッシュ値が異なる場合，
		// パーティクルは分割領域の最初
		if(index == 0 || hash != sharedHash[threadIdx.x]){
			cell.dCellStart[hash] = index;
			if(index > 0){
				// 一つ前のパーティクルは，一つ前の分割領域の最後
				cell.dCellEnd[sharedHash[threadIdx.x]] = index;
			}
		}

		// インデックスが最後ならば，分割領域の最後
		if(index == cell.uNumParticles-1){
			cell.dCellEnd[hash] = index+1;
		}

		// 位置と速度のデータを並び替え
		// ソートしたインデックスで参照も可能だが探索時のグローバルメモリアクセスを極力抑えるためにデータそのものを並び替える
		uint sortedIndex = cell.dSortedIndex[index];
		float4 pos = FETCH(dSortedPos, sortedIndex);
		float4 vel = FETCH(dSortedVel, sortedIndex);

		cell.dSortedPos[index] = pos;
		cell.dSortedVel[index] = vel;
	}
}

/*!
 * パーティクルデータをソートして，ハッシュ内の各セルの最初のアドレスを検索
 * @param[in] cell パーティクルグリッドデータ
 * @param[in] oldPos パーティクル位置
 * @param[in] oldVel パーティクル速度
 */
__global__
void reorderDataAndFindCellStartB(rxParticleCell cell, float4* dPos)
{
	extern __shared__ uint sharedHash[];	// サイズ : blockSize+1
	uint index = __umul24(blockIdx.x,blockDim.x)+threadIdx.x;
	
	uint hash;
	if(index < cell.uNumParticles){
		hash = cell.dGridParticleHash[index];	// ハッシュ値

		sharedHash[threadIdx.x+1] = hash;	// ハッシュ値をシェアードメモリに格納

		if(index > 0 && threadIdx.x == 0){
			// 各シェアードメモリの最初は隣のグリッドのパーティクルのハッシュ値を格納
			sharedHash[0] = cell.dGridParticleHash[index-1];
		}
	}

	__syncthreads();
	
	if(index < cell.uNumParticles){
		// インデックス0である，もしくは，一つ前のパーティクルのグリッドハッシュ値が異なる場合，
		// パーティクルは分割領域の最初
		if(index == 0 || hash != sharedHash[threadIdx.x]){
			cell.dCellStart[hash] = index;
			if(index > 0){
				// 一つ前のパーティクルは，一つ前の分割領域の最後
				cell.dCellEnd[sharedHash[threadIdx.x]] = index;
			}
		}

		// インデックスが最後ならば，分割領域の最後
		if(index == cell.uNumParticles-1){
			cell.dCellEnd[hash] = index+1;
		}

		// 位置と速度のデータを並び替え
		// ソートしたインデックスで参照も可能だが探索時のグローバルメモリアクセスを極力抑えるためにデータそのものを並び替える
		uint sortedIndex = cell.dSortedIndex[index];
		float4 pos = dPos[sortedIndex];
		cell.dSortedPos[index] = pos;
	}
}

//-----------------------------------------------------------------------------
// MARK:SPH
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
float calDensityCell(int3 gridPos, uint i, float3 pos0, rxParticleCell cell)
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
void sphCalDensity(float* newDens, float* newPres, rxParticleCell cell)
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
				dens += calDensityCell(n_grid_pos, index, pos, cell);
			}
		}
	}

	// ガス定数を使った圧力算出
	float pres;
	pres = params.GasStiffness*(dens-params.Density);

	// 密度と圧力値を結果に書き込み
	uint oIdx = cell.dSortedIndex[index];
	newDens[oIdx] = dens;
	newPres[oIdx] = pres;
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
float3 calNormalCell(int3 gridPos, uint i, float3 pos0, float* dens, rxParticleCell cell)
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
void sphCalNormal(float4* newNrms, float* dens, rxParticleCell cell)
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
				nrm += calNormalCell(n_grid_pos, index, pos, dens, cell);
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




/*!
 * 与えられたセル内のパーティクルとの距離から力場を計算
 * @param[in] gridPos グリッド位置
 * @param[in] i パーティクルインデックス
 * @param[in] pos0 計算座標
 * @param[in] vel0 計算座標の速度
 * @param[in] dens0 計算座標の密度
 * @param[in] pres0 計算座標の圧力
 * @param[in] dens パーティクル密度
 * @param[in] pres パーティクル圧力
 * @param[in] cell パーティクルグリッドデータ
 * @return セル内のパーティクルから計算した力場
 */
__device__
float3 calForceCell(int3 gridPos, uint i, float3 pos0, float3 vel0, float dens0, float pres0, 
					float* dens, float* pres, rxParticleCell cell)
{
	uint gridHash = calcGridHash(gridPos);

	// セル内のパーティクルのスタートインデックス
	uint startIndex = FETCHC(dCellStart, gridHash);

	float h = params.EffectiveRadius;

	float3 frc = make_float3(0.0f);
	if(startIndex != 0xffffffff){	// セルが空でないかのチェック
		// セル内のパーティクルで反復
		uint endIndex = FETCHC(dCellEnd, gridHash);
		float prsi = pres0/(dens0*dens0);
		for(uint j = startIndex; j < endIndex; ++j){
			if(j != i){
				// 近傍パーティクルのパラメータ
				float3 pos1 = make_float3(FETCHC(dSortedPos, j));
				float3 vel1 = make_float3(FETCHC(dSortedVel, j));

				float3 rij = pos0-pos1;
				float r = length(rij);

				if(r <= h && r > 0.0001){
					//float3 vel1 = make_float3(vel[cell.dSortedIndex[j]]);
					float dens1 = dens[cell.dSortedIndex[j]];
					float pres1 = pres[cell.dSortedIndex[j]];

					float3 vji = vel1-vel0;

					float prsj = pres1/(dens1*dens1);
					float q = h-r;

					// 圧力項
					frc += -params.Mass*(prsi+prsj)*params.GWspiky*q*q*rij/r;

					// 粘性項
					frc += params.Viscosity*params.Mass*(vji/dens1)*params.LWvisc*q;
				}
			}
		}
	}

	return frc;
}

/*!
 * パーティクルにかかる力の計算(カーネル関数)
 * @param[in] dens パーティクル密度
 * @param[in] pres パーティクル圧力
 * @param[out] outFrc パーティクルにかかる力
 * @param[in] cell パーティクルグリッドデータ
 */
__global__
void sphCalForces(float* dens, float* pres, float4* outFrc, rxParticleCell cell)
{
	uint index = __mul24(blockIdx.x,blockDim.x)+threadIdx.x;
	if(index >= cell.uNumParticles) return;	
	
	// ソート済み配列からパーティクルデータを取得
	float3 pos0 = make_float3(FETCHC(dSortedPos, index));
	float3 vel0 = make_float3(FETCHC(dSortedVel, index));

	int3 gridPos0 = calcGridPos(pos0);

	// パーティクルのソートなし配列上でのインデックス
	uint oIdx = cell.dSortedIndex[index];

	float dens0 = dens[oIdx];
	float pres0 = pres[oIdx];

	float h = params.EffectiveRadius;
	int3 grid_pos0, grid_pos1;
	grid_pos0 = calcGridPos(pos0-make_float3(h));
	grid_pos1 = calcGridPos(pos0+make_float3(h));

	// 周囲のグリッドも含めて近傍探索，圧力項，粘性項を計算
	float3 frc = make_float3(0.0f);
	for(int z = grid_pos0.z; z <= grid_pos1.z; ++z){
		for(int y = grid_pos0.y; y <= grid_pos1.y; ++y){
			for(int x = grid_pos0.x; x <= grid_pos1.x; ++x){
				int3 n_grid_pos = make_int3(x, y, z);

				frc += calForceCell(n_grid_pos, index, pos0, vel0, dens0, pres0, dens, pres, cell);
			}
		}
	}

	// 外力(重力)
	frc += params.Gravity;

	outFrc[oIdx] = make_float4(frc, 0.0f);
}


/*!
 * 与えられたセル内のパーティクルとの距離から境界パーティクルの体積を計算
 * @param[in] gridPos グリッド位置
 * @param[in] index パーティクルインデックス
 * @param[in] pos 計算座標
 * @param[in] cell パーティクルグリッドデータ
 * @return セル内のパーティクルから計算した密度値
 */
__device__
float calBoundaryVolumeCell(int3 gridPos, uint i, float3 pos0, rxParticleCell cell)
{
	uint gridHash = calcGridHashB(gridPos, params.GridSizeB);

	// セル内のパーティクルのスタートインデックス
	uint startIndex = cell.dCellStart[gridHash];

	float h = params.EffectiveRadius;
	float mw = 0.0f;
	if(startIndex != 0xffffffff){	// セルが空でないかのチェック
		// セル内のパーティクルで反復
		uint endIndex = cell.dCellEnd[gridHash];
		for(uint j = startIndex; j < endIndex; ++j){
			float3 pos1 = make_float3(cell.dSortedPos[j]);

			float3 rij = pos0-pos1;
			float r = length(rij);

			if(r <= h){
				float q = h*h-r*r;
				mw += params.Mass*params.Wpoly6*q*q*q;
			}
		}
	}

	return mw;
}

/*!
 * 境界パーティクルの体積計算(カーネル関数)
 * @param[out] newVolB パーティクル体積
 * @param[in]  cell 境界パーティクルグリッドデータ
 */
__global__
void sphCalBoundaryVolume(float* newVolB, rxParticleCell cell)
{
	// パーティクルインデックス
	uint index = __mul24(blockIdx.x,blockDim.x)+threadIdx.x;
	if(index >= cell.uNumParticles) return;	
	
	float3 pos = make_float3(cell.dSortedPos[index]);	// パーティクル位置
	//int3 grid_pos = calcGridPos(pos);	// パーティクルが属するグリッド位置
	float h = params.EffectiveRadius;

	// パーティクル周囲のグリッド
	int3 grid_pos0, grid_pos1;
	grid_pos0 = calcGridPosB(pos-make_float3(h), params.WorldOriginB, params.CellWidthB, params.GridSizeB);
	grid_pos1 = calcGridPosB(pos+make_float3(h), params.WorldOriginB, params.CellWidthB, params.GridSizeB);


	// 周囲のグリッドも含めて近傍探索
	float mw = 0.0f;
	for(int z = grid_pos0.z; z <= grid_pos1.z; ++z){
		for(int y = grid_pos0.y; y <= grid_pos1.y; ++y){
			for(int x = grid_pos0.x; x <= grid_pos1.x; ++x){
				int3 n_grid_pos = make_int3(x, y, z);
				mw += calBoundaryVolumeCell(n_grid_pos, index, pos, cell);
			}
		}
	}

	// 体積を結果に書き込み
	uint oIdx = cell.dSortedIndex[index];
	newVolB[oIdx] = params.Mass/mw;
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
float calBoundaryDensityCell(int3 gridPos, uint i, float3 pos0, float* dVolB, rxParticleCell bcell)
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
void sphCalBoundaryDensity(float* newDens, float* newPres, float4* dPos, float* dVolB, rxParticleCell bcell, uint pnum)
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
				dens += calBoundaryDensityCell(n_grid_pos, index, pos, dVolB, bcell);
			}
		}
	}

	dens += newDens[index];

	// ガス定数を使った圧力算出
	float pres;
	pres = params.GasStiffness*(dens-params.Density);

	// 密度と圧力値を結果に書き込み
	newDens[index] = dens;
	newPres[index] = pres;
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
float3 calBoundaryForceCell(int3 gridPos, uint i, float3 pos0, float* dVolB, float dens0, float pres0, rxParticleCell bcell)
{
	uint gridHash = calcGridHashB(gridPos, params.GridSizeB);

	// セル内のパーティクルのスタートインデックス
	uint startIndex = bcell.dCellStart[gridHash];

	float h = params.EffectiveRadius;
	float3 bp = make_float3(0.0f);
	if(startIndex != 0xffffffff){	// セルが空でないかのチェック
		// セル内のパーティクルで反復
		uint endIndex = bcell.dCellEnd[gridHash];
		float prsi = pres0/(dens0*dens0);
		for(uint j = startIndex; j < endIndex; ++j){
			float3 pos1 = make_float3(bcell.dSortedPos[j]);
			uint jdx = bcell.dSortedIndex[j];

			float3 rij = pos0-pos1;
			float r = length(rij);

			if(r <= h && r > 0.0001){
				float q = h-r;
				bp += -params.Density*dVolB[jdx]*prsi*params.GWspiky*q*q*rij/r;
			}
		}
	}

	return bp;
}

/*!
 * 境界パーティクルによる力の計算(カーネル関数)
 * @param[out] newDens パーティクル密度
 * @param[out] newPres パーティクル圧力
 * @param[in]  cell パーティクルグリッドデータ
 */
__global__
void sphCalBoundaryForce(float* dDens, float* dPres, float4* dPos, float* dVolB, float4* outFrc, rxParticleCell bcell, uint pnum)
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

	// 密度と圧力
	float dens0 = dDens[index];
	float pres0 = dPres[index];

	// 周囲のグリッドも含めて近傍探索，密度計算
	float3 frc = make_float3(0.0f);
	for(int z = grid_pos0.z; z <= grid_pos1.z; ++z){
		for(int y = grid_pos0.y; y <= grid_pos1.y; ++y){
			for(int x = grid_pos0.x; x <= grid_pos1.x; ++x){
				int3 n_grid_pos = make_int3(x, y, z);
				frc += calBoundaryForceCell(n_grid_pos, index, pos, dVolB, dens0, pres0, bcell);
			}
		}
	}

	// 密度と圧力値を結果に書き込み
	outFrc[index] += make_float4(frc, 0.0f);
}

__device__
void calCollisionSolid(float3 &pos, float3 &vel, float dt)
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
inline bool calCollisionPolygon(float3 &pos0, float3 &pos1, float3 &vel, float3 v0, float3 v1, float3 v2, float dt)
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
void sphIntegrate(float4* ppos,	float4* pvel, 
				  float4* pacc, float* dens, int* attr, float dt, uint nprts)
{
	uint index = __umul24(blockIdx.x,blockDim.x)+threadIdx.x;
	if(index >= nprts) return;
	if(attr[index] == -1) return;

	float3 x = make_float3(ppos[index]);
	float3 v = make_float3(pvel[index]);
	float3 a = make_float3(pacc[index]);
	//float3 v_old = v;

	// 更新位置，速度の更新
	v += dt*a;
	x += dt*v;

	// 固体・境界との衝突
	calCollisionSolid(x, v, dt);

	// 位置と速度の更新
	ppos[index] = make_float4(x);
	pvel[index] = make_float4(v);
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
void sphIntegrateWithPolygon(float4* ppos, float4* pvel, float4* pacc, float* dens, int* attr, 
							 float3* vrts, int3* tris, int tri_num, float dt, rxParticleCell cell)
{
	uint index = __umul24(blockIdx.x,blockDim.x)+threadIdx.x;
	if(index >= cell.uNumParticles) return;
	if(attr[index] == -1) return;

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
				if(calCollisionPolygon(x_old, x, v, vrts[idx.x], vrts[idx.y], vrts[idx.z], dt)){
				}
			}
		}
	}

	// 固体・境界との衝突
	calCollisionSolid(x, v, dt);

	// 位置と速度の更新
	ppos[index] = make_float4(x);
	pvel[index] = make_float4(v);
}

/*!
 * パーティクル位置をチェックして削除領域内ならば属性を-1にして，範囲外に飛ばす
 * @param[inout] ppos パーティクル位置
 * @param[in] dt 時間ステップ幅
 * @param[in] nprts パーティクル数
 */
__global__
void checkDelete(float4* ppos, float4* pvel, int* attr, float3 minp, float3 maxp, float3 far_pos, uint nprts)
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
void checkDeleteX(float4* ppos, float4* pvel, int* attr, float xmax, float3 far_pos, uint nprts)
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
float calDensityCellG(int3 gridPos, float3 pos0, rxParticleCell cell)
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
				//cell.dSortedIndex[j];

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

/*! 追加：氷用
 * 与えられたセル内のパーティクルとの距離から密度を計算
 * @param[in] gridPos グリッド位置
 * @param[in] pos 計算座標
 * @param[in] cell パーティクルグリッドデータ
 * @return セル内のパーティクルから計算した密度値
 */
__device__
float calDensityCellGIceMesh(int3 gridPos, float3 pos0, rxParticleCell cell, float* bIceCheck)
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
				//cell.dSortedIndex[j];

				float3 rij = pos0-pos1;
				float r = length(rij);

				if(r <= h){
					uint pIndx = cell.dSortedIndex[j];
					if( bIceCheck[pIndx] < 0.0f ) continue;		//氷でないならメッシュを作らない
//					if( bIceCheck[j] < 0.0f ) continue;		//氷でないならメッシュを作らない
//					printf("pIndx = %d\n", pIndx);
//					printf("j = %d\n", j);
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
void sphCalDensityInGrid(float* GridD, rxParticleCell cell, 
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

					d += calDensityCellG(neighbourPos, gpos, cell);
				}
			}
		}

		GridD[gridPos.x+gridPos.y*gnum.x+gridPos.z*gnum.x*gnum.y] = d;
	}

}


/*!追加：氷用
 * グリッド上での密度を計算
 * @param[out] GridD グリッド密度
 * @param[in] cell パーティクルグリッドデータ
 * @param[in] gnum グリッド数
 * @param[in] gmin グリッド最小座標
 * @param[in] glen グリッド幅
 */
__global__
void sphCalDensityInGridIceMesh(float* GridD, rxParticleCell cell, 
					uint3 gnum, float3 gmin, float3 glen, float* bIceFlag)
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

					d += calDensityCellGIceMesh(neighbourPos, gpos, cell, bIceFlag);
				}
			}
		}

		GridD[gridPos.x+gridPos.y*gnum.x+gridPos.z*gnum.x*gnum.y] = d;
	}

}

#endif // #ifndef _RX_CUSPH_KERNEL_CU_



