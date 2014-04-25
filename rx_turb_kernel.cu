/*! 
  @file rx_turb_kernel.cu
	
  @brief CUDAによるSPH乱流計算
 
  @author Makoto Fujisawa
  @date 2009-08, 2011-06
*/

#ifndef _RX_TURBULENCE_KERNEL_CU_
#define _RX_TURBULENCE_KERNEL_CU_


//-----------------------------------------------------------------------------
// インクルードファイル
//-----------------------------------------------------------------------------
#include "rx_cu_common.cu"

//-----------------------------------------------------------------------------
// 変数
//-----------------------------------------------------------------------------
#if USE_TEX
texture<float4, cudaTextureType1D, cudaReadModeElementType> dSortedPosTex;
texture<float4, cudaTextureType1D, cudaReadModeElementType> dSortedVelTex;
texture<uint,   cudaTextureType1D, cudaReadModeElementType> dCellStartTex;
texture<uint,   cudaTextureType1D, cudaReadModeElementType> dCellEndTex;

//サブパーティクル
texture<float4, cudaTextureType1D, cudaReadModeElementType> dSubUnsortPosTex;
texture<float,  cudaTextureType1D, cudaReadModeElementType> dSubUnsortRadTex;
texture<float,  cudaTextureType1D, cudaReadModeElementType> dSubUnsortRatTex;
texture<float4, cudaTextureType1D, cudaReadModeElementType> dSubSortedPosTex;
texture<float,  cudaTextureType1D, cudaReadModeElementType> dSubSortedRadTex;
texture<float,  cudaTextureType1D, cudaReadModeElementType> dSubSortedRatTex;
texture<uint,   cudaTextureType1D, cudaReadModeElementType> dSubCellStartTex;
texture<uint,   cudaTextureType1D, cudaReadModeElementType> dSubCellEndTex;
#endif

texture<float, cudaTextureType1D, cudaReadModeElementType> g_TexNoiseTile3D;


__constant__ int RXNA[] = {0, -1, 1, -2, 2, -3, 3, -4, 4, -5, 5, -6, 6, -7, 7, -8, 8};

/*!
 * 各パーティクルのグリッドハッシュ値
 * @param[out] gridParticleHash
 * @param[out] dSortedIndex
 * @param[in] pos パーティクル位置を格納した配列
 * @param[in] nprts パーティクル数
 */
__global__
void calcHashD3(uint*   dGridParticleHash, 
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


//-----------------------------------------------------------------------------
// MARK:ウェーブレット乱流
//-----------------------------------------------------------------------------

/*!
 * 与えられたセル内のパーティクルから連続ウェーブレット変換を計算
 * @param[in] gridPos グリッド位置
 * @param[in] index パーティクルインデックス
 * @param[in] pos0 計算座標
 * @param[in] vel0 計算座標での速度
 * @param[in] cell パーティクルグリッドデータ
 */
__device__
float4 calParticleWt(int3 gridPos, uint i, int &c, float3 pos0, float3 vel0, rxParticleCell cell, float h, float scale)
{
	uint gridHash = calcGridHash(gridPos);

	// セル内のパーティクルのスタートインデックス
	uint startIndex = FETCHC(dCellStart, gridHash);

	//float h = scale*MEXICAN_HAT_R;
	float4 wt = make_float4(0.0f);
	if(startIndex != 0xffffffff){	// セルが空でないかのチェック
		// セル内のパーティクルで反復
		uint endIndex = FETCHC(dCellEnd, gridHash);

		for(uint j = startIndex; j < endIndex; ++j){
			if(j != i){
				float3 pos1 = make_float3(FETCHC(dSortedPos, j));
				float3 vel1 = make_float3(FETCHC(dSortedVel, j));

				float3 rij = pos0-pos1;
				float r = length(rij);

				if(r <= h){
					float Tx = rij.x/scale;
					float Ty = rij.y/scale;
					float Tz = rij.z/scale;

					float w = MexicanHat3D(Tx, Ty, Tz);

					wt.x += vel1.x*w;
					wt.y += vel1.y*w;
					wt.z += vel1.z*w;

					wt.w += fabs(w);

					c++;
				}

			}
		}
	}

	return wt;
}

/*!
 * パーティクル速度のエネルギースペクトラムを計算
 * @param[out] Et エネルギースペクトラム
 * @param[in] scale ウェーブレットスケール
 * @param[in] coef_et エネルギースペクトラム値にかける係数
 * @param[in] max_et  エネルギースペクトラム値の最大値(max_et以上ならクランプ)
 * @param[in] cell パーティクルグリッドデータ
 */
__global__
void calParticleES(float* Et, float h, float scale, float coef_et, float max_et, rxParticleCell cell, int na)
{
	uint index = __mul24(blockIdx.x,blockDim.x)+threadIdx.x;
	if(index >= cell.uNumParticles) return;	
	
	// read particle data from sorted arrays
	float3 pos = make_float3(FETCHC(dSortedPos, index));
	float3 vel = make_float3(FETCHC(dSortedVel, index));

	//int3 gridPos = calcGridPos(pos);
	//float h = params.EffectiveRadius;
	int3 grid_pos0, grid_pos1;
	grid_pos0 = calcGridPos(pos-make_float3(h));
	grid_pos1 = calcGridPos(pos+make_float3(h));

	int c = 0;
	float4 wt = make_float4(0.0f);
	for(int z = grid_pos0.z; z <= grid_pos1.z; ++z){
		for(int y = grid_pos0.y; y <= grid_pos1.y; ++y){
			for(int x = grid_pos0.x; x <= grid_pos1.x; ++x){
				int3 n_grid_pos = make_int3(x, y, z);
				wt += calParticleWt(n_grid_pos, index, c, pos, vel, cell, h, scale);
			}
		}
	}

	if(wt.w > 1.0e-6){
		float s = (sqrtf(scale)*wt.w)/(c/50.0);
		wt.x /= s;
		wt.y /= s;
		wt.z /= s;
	}

	float ev = 0.5*dot(vel, vel);
	float et = coef_et*ev*0.5*(wt.x*wt.x+wt.y*wt.y+wt.z*wt.z);
	if(et > max_et) et = max_et;

	uint oIdx = cell.dSortedIndex[index];	// ソートしていないときのインデックス
	Et[oIdx] = et;
}

__device__
float wdnoise3d(float *tile, int n[3], int offset, int tile_size, float p[3], int d)
{
	int f[3], c[3];		// filter, noise coef. indices
	int mid[3];

	float w[3][3], t, result = 0;

	// 2次のB-スプライン(quadratic B-spline)基底関数を計算
	//  [t^2/2, 1/2+t-t^2, (1-t)^2/2]
	for(int k = 0; k < 3; ++k){
		mid[k] = ceil(p[k]-0.5);
		t = mid[k]-(p[k]-0.5);

		if(k == d){
			w[k][0] = -t; 
			w[k][1] = 2*t-1;
			w[k][2] = 1-t;
		}
		else{
			w[k][0] = t*t/2; 
			w[k][2] = (1-t)*(1-t)/2; 
			w[k][1] = 1-w[k][0]-w[k][2];
		}
	}

	// ノイズタイルを基底関数の値で重み付け補間する
	for(f[2] = -1; f[2] <= 1; ++f[2]){
		for(f[1] = -1; f[1] <= 1; ++f[1]){
			for(f[0] = -1; f[0] <= 1; ++f[0]){
				float weight = 1;
				for(int k = 0; k < 3; ++k){
					c[k] = Mod(mid[k]+f[k], n[k]);
					weight *= w[k][f[k]+1];
				}

				int tl = c[2]*n[0]*n[1]+c[1]*n[0]+c[0]+offset;
				if(tl >= tile_size) tl -= tile_size;

				//result += weight*tile[tl];
				result += weight*tex1Dfetch(g_TexNoiseTile3D, tl);
			}
		}
	}


	return result;
}

__global__
void calWaveletTurbulence3D(float4 *turb, float4 *frc, float *Et, float *dens, 
							int first, int nbands, float *wtile, int3 tile_n, int tile_size, 
							float4* posArray, float3 pmin, float3 pdim, uint nprts, float dt)
{
	uint index = __umul24(blockIdx.x,blockDim.x)+threadIdx.x;
	if(index >= nprts) return;

	float e = Et[index];
	float3 p = make_float3(posArray[index]);


	float3 y;
	y.x = 0.0;
	y.y = 0.0;
	y.z = 0.0;

	int3 offset;
	offset.x = 0;
	offset.y = 0.457*tile_n.x;
	offset.z = 0.896*tile_n.x;

	int n[3];
	n[0] = tile_n.x;
	n[1] = tile_n.y;
	n[2] = tile_n.z;

	float q[3];
	for(int b = 0; b < nbands; ++b){
		float scl = powf(2.0f, (float)(first+b));

		q[0] = ((p.x-pmin.x)/pdim.x)*scl;
		q[1] = ((p.y-pmin.y)/pdim.y)*scl;
		q[2] = ((p.z-pmin.z)/pdim.z)*scl;

		float3 w;
		w.x = wdnoise3d(wtile, n, offset.x, tile_size, q, 1)-wdnoise3d(wtile, n, offset.y, tile_size, q, 2);
		w.y = wdnoise3d(wtile, n, offset.z, tile_size, q, 2)-wdnoise3d(wtile, n, offset.x, tile_size, q, 0);
		w.z = wdnoise3d(wtile, n, offset.y, tile_size, q, 0)-wdnoise3d(wtile, n, offset.z, tile_size, q, 1);
		//w.x = pmin.x;//wnoise3d_(wtile, tile_n, tile_size, q);
		//w.y = pdim.x;//wnoise3d_(wtile, tile_n, tile_size, q);
		//w.z = scl;//wnoise3d_(wtile, tile_n, tile_size, q);

		float wl = sqrtf(w.x*w.x+w.y*w.y+w.z*w.z);

		if(wl > 1e-6){
			w.x /= wl;
			w.y /= wl;
			w.z /= wl;
		}
		else{
			w = make_float3(0.0);
		}

		float c = powf(2.0f, (float)(-5.0/6.0*b));
		y.x += c*w.x;
		y.y += c*w.y;
		y.z += c*w.z;

		//y += c*w;
	}

	float dens0 = dens[index];

	float3 tr;
	tr.x = e*y.x;
	tr.y = e*y.y;
	tr.z = e*y.z;

	//turb[index].x = tr.x*dens0;
	//turb[index].y = tr.y*dens0;
	//turb[index].z = tr.z*dens0;

	frc[index].x += tr.x*dens0/dt;
	frc[index].y += tr.y*dens0/dt;
	frc[index].z += tr.z*dens0/dt;
}




//-----------------------------------------------------------------------------
// MARK:SPS乱流
//-----------------------------------------------------------------------------
/*!
 * ある位置とグリッドの最短距離を計算
 * @param[in] p		グリッド座標
 * @param[in] gridPos グリッド座標
 * @return 距離
 */
__device__ 
float calcGridPosDisMin(float3 p, int3 gridPos)
{
	float3 dis2,gpos1,gpos2;

	gpos1.x = params.WorldOrigin.x + params.CellWidth.x * (float)gridPos.x;
	gpos1.y = params.WorldOrigin.y + params.CellWidth.y * (float)gridPos.y; 
	gpos1.z = params.WorldOrigin.z + params.CellWidth.z * (float)gridPos.z;
	gpos2.x = gpos1.x + params.CellWidth.x;
	gpos2.y = gpos1.y + params.CellWidth.y;
	gpos2.z = gpos1.z + params.CellWidth.z;
	dis2.x = min( (gpos1.x-p.x)*(gpos1.x-p.x) , (gpos2.x-p.x)*(gpos2.x-p.x));
	dis2.y = min( (gpos1.y-p.y)*(gpos1.y-p.y) , (gpos2.y-p.y)*(gpos2.y-p.y));
	dis2.z = min( (gpos1.z-p.z)*(gpos1.z-p.z) , (gpos2.z-p.z)*(gpos2.z-p.z));
	
	if(p.x>= gpos1.x && p.x<=gpos2.x) dis2.x=0.0;
	if(p.y>= gpos1.y && p.y<=gpos2.y) dis2.y=0.0;
	if(p.z>= gpos1.z && p.z<=gpos2.z) dis2.z=0.0;

	return pow((dis2.x+dis2.y+dis2.z), 0.50f);
}

__device__
float3 getWorldMaxPos(void)
{
	float3 maxPos;
	maxPos.x = params.WorldOrigin.x + params.CellWidth.x * (float)params.GridSize.x;
	maxPos.y = params.WorldOrigin.y + params.CellWidth.y * (float)params.GridSize.y;
	maxPos.z = params.WorldOrigin.z + params.CellWidth.z * (float)params.GridSize.z;
	return maxPos;
}

__global__
void checkSubRatio(uint* dGridParticleHash,
				   float* dSubRat,
				   uint uHashMax,
					uint uNum)
{
	uint index = __umul24(blockIdx.x,blockDim.x)+threadIdx.x;
	if(index >= uNum) return;

	float rat = dSubRat[index];
	if(rat <= 0.0f) dGridParticleHash[index] = uHashMax;
}

/*!
 * Quaternionの掛け算
 */
__device__
float4 calMultiQuaternion(float4 a, float4 b){
	float4 ab;

	ab.x = a.x*b.x - a.y*b.y - a.z*b.z - a.w*b.w;
	ab.y = a.x*b.y + a.y*b.x + a.z*b.w - a.w*b.z;
	ab.z = a.x*b.z - a.y*b.w + a.z*b.x + a.w*b.y;
	ab.w = a.x*b.w + a.y*b.z - a.z*b.y + a.w*b.x;

	return ab;
}

__device__
float4 calRotateQuaternion(float4 q, float3 a, float theta){

	float inv_a,cos_half_theta,sin_half_theta,tmp;
	float4 dq,new_q;

	inv_a = sqrt(a.x*a.x + a.y*a.y + a.z*a.z);

	if(inv_a <= 0.0){
		return q;
	}
	else{
		inv_a = 1.0/inv_a;

		cos_half_theta = cos(0.50*theta);
		sin_half_theta = sin(0.50*theta);

		tmp = inv_a * sin_half_theta;
		dq = make_float4(cos_half_theta, tmp*a.x, tmp*a.y, tmp*a.z);

		new_q = calMultiQuaternion(dq,q);

		dq = make_float4(cos_half_theta, -tmp*a.x, -tmp*a.y, -tmp*a.z);

		new_q = calMultiQuaternion(new_q,dq);
	}

	return new_q;

}

__device__
float4 calRotateQuaternionPos(float4 pos,float3 a,const float theta){

	float4 posq = make_float4(0.0, pos.x, pos.y, pos.z);

	posq = calRotateQuaternion(posq,a,theta);

	float tmp = sqrt( posq.y*posq.y + posq.z*posq.z + posq.w*posq.w);
	
	float4 new_pos;
	if(tmp > 0.0){
		tmp = 1.0/tmp;
		new_pos = make_float4(tmp * posq.y, tmp * posq.z, tmp * posq.w, 0.0);
	}
	else  new_pos = pos;

	return new_pos;
}

/*!
 * サブ粒子(レベル0)のエネルギースペクトル更新
 * @param[out]	dSubEt	サブ粒子のエネルギースペクトルへのアドレス
 * @param[in]	dEt		粒子のエネルギースペクトルへのアドレス
 * @param[in]	scale		そのスケール
 * @param[in]	radius		レベル0での半径
 * @param[in]	nprts パーティクル数
 */
__global__
void setSubEt(	float* dSubEt,
				float* dEt,
				float et_coef, 
				float	scale, 
				float	radius,
				uint nprts)
{
	uint index = __umul24(blockIdx.x,blockDim.x)+threadIdx.x;
	if(index >= nprts) return;

	float et = dEt[index];

	dSubEt[index] = et_coef*et*pow(scale/(2.0f*radius), -5.0f/3.0f);//* POW_2_M5D9
} 

/*!
 * エネルギースペクトルの更新
 * @param[out]	dSubBlockEtC	サブ粒子(子)のエネルギースペクトルへのアドレス
 * @param[in]	dSubBlockEtP	サブ粒子(親)のエネルギースペクトルへのアドレス
 * @param[in]	nprts パーティクル数
 */
__global__
void updateSubEt(	float* dSubBlockEtC,
					float* dSubBlockEtP,
					uint nprts)
{
	uint index = __umul24(blockIdx.x,blockDim.x)+threadIdx.x;
	if(index >= nprts) return;

	uint	indexc1	= 2*index;
	uint	indexc2	= 2*index+1;
	
	float etp	= dSubBlockEtP[index];

	float etc	= POW_2_M5D9 * etp;

	dSubBlockEtC[indexc1] = etc;
	dSubBlockEtC[indexc2] = etc;
} 

/*!
 * サブ粒子(レベル0)の位置更新
 * @param[in]	dSubPos	サブ粒子の絶対座標へのアドレス
 * @param[in]	dPos	粒子の絶対座標へのアドレス
 * @param[in]	nprts パーティクル数
 */
__global__
void setSubPos(	float4* dSubPos,
				float4* dPos,
				uint nprts)
{
	uint index = __umul24(blockIdx.x,blockDim.x)+threadIdx.x;
	if(index >= nprts) return;

	float4 posData = dPos[index];

	//float4 pos = make_float4(posData.x, posData.y, posData.z, 0.0);

	dSubPos[index] = posData;
} 


/*!
 * サブ粒子の位置更新
 * @param[in]	dSubBlockPosC	サブ粒子(子)の絶対座標へのアドレス
 * @param[in]	dSubBlockPosP	サブ粒子(親)の絶対座標へのアドレス
 * @param[in]	dSubBlockChild	サブ粒子(親)の子1への単位ベクトルへのアドレス
 * @param[in]	radius_level	サブ粒子(親)の半径
 * @param[in]	nprts パーティクル数
 */
__global__
void calSubPos(	float4* dSubBlockPosC,
				float4* dSubBlockPosP,
				float4* dSubBlockChild,
				float radius_level,
				uint nprts)
{
	uint index = __umul24(blockIdx.x,blockDim.x)+threadIdx.x;
	if(index >= nprts) return;

	uint	indexc1	= 2*index;
	uint	indexc2	= 2*index+1;
	
	float4 pos		= dSubBlockPosP[index];
	float4 child	= dSubBlockChild[index];

	//float4 pos_c1	= pos + radius_level * child;
	//float4 pos_c2	= pos - radius_level * child;

	float4 pos_c1	= make_float4(pos.x + radius_level * child.x, pos.y + radius_level * child.y, pos.z + radius_level * child.z, 0.0);
	float4 pos_c2	= make_float4(pos.x - radius_level * child.x, pos.y - radius_level * child.y, pos.z - radius_level * child.z, 0.0);

	dSubBlockPosC[indexc1] = pos_c1;
	dSubBlockPosC[indexc2] = pos_c2;
}


/*!
 * サブ粒子の回転(単位ベクトル更新)
 * @param[inout]dSubBlockChild	サブ粒子の子1への単位ベクトルの行列
 * @param[in]	dSubBlockAxis	サブ粒子のエネルギースペクトルの配列
 * @param[in]	dSubBlockEt		サブ粒子のエネルギースペクトルの配列
 * @param[in]	radius_level	このサブ粒子のの半径
 * @param[in]	nprts パーティクル数
 * @param[in]	dt				微少時間
 */
__global__
void updateSubChildAxis(float4* dSubBlockChild,
						float4* dSubBlockAxis,
						float*	dSubBlockEt,
						uint*	dSubRand,
						float	radius_level,
						uint	nprts,
						float	ratioAxisAngle,
						float	dt)
{
	uint index = __umul24(blockIdx.x,blockDim.x)+threadIdx.x;
	if(index >= nprts) return;

	float4 child	= dSubBlockChild[index];
	float4 axis		= dSubBlockAxis[index];
	float  et		= dSubBlockEt[index];
	uint   pre_rand	= dSubRand[index];

	uint rand	= Rand2(pre_rand);
	float theta = sqrt(2.0f*et* POW_2_M5D9) * dt / radius_level;
	float rand_theta = (-1.0+2.0*((float)rand)/RAND2_MAX) * theta * ratioAxisAngle;
	//float rand_theta = XorFrand(-1.0,1.0) * theta * ratioAxisAngle;
	
	float3 axisF3 = make_float3(axis.x,axis.y,axis.z);
	float4 newChild = calRotateQuaternionPos(child,axisF3,theta);

	float3 newChildF3 = make_float3(newChild.x, newChild.y, newChild.z);
	float4 newAxis = calRotateQuaternionPos(axis,newChildF3,rand_theta);

	dSubBlockChild[index]	= newChild;
	dSubBlockAxis[index]	= newAxis;
	dSubRand[index]			= rand;
}


/*!
 * サブ粒子の回転(単位ベクトル更新)
 * @param[inout]dSubBlockChild	サブ粒子の子1への単位ベクトルの行列
 * @param[in]	dSubBlockAxis	サブ粒子のエネルギースペクトルの配列
 * @param[in]	dSubRand		uint乱数
 * @param[in]	nprts パーティクル数
 */
__global__
void initSub(	float4* dSubBlockChild,
				float4* dSubBlockAxis,
				uint*	dSubRand,
				uint nprts)
{
	uint index = __umul24(blockIdx.x,blockDim.x)+threadIdx.x;
	if(index >= nprts) return;
	
	float tmp = 0.0f;

	uint rand = dSubRand[index];

	float3 randf3 = make_float3(0.0,0.0,0.0);

	while(randf3.x == 0.0){
		rand = Rand2(rand);
		randf3.x = -1.0+2.0*(float)rand/RAND2_MAX;
	}
	while(randf3.y == 0.0){
		rand = Rand2(rand);
		randf3.y = -1.0+2.0*(float)rand/RAND2_MAX;
	}
	while(randf3.z == 0.0){
		rand = Rand2(rand);
		randf3.z = -1.0+2.0*(float)rand/RAND2_MAX;
	}

	tmp = 1.0f/sqrt(randf3.x*randf3.x + randf3.y*randf3.y + randf3.z*randf3.z);
	randf3 = make_float3(tmp*randf3.x, tmp*randf3.y, tmp*randf3.z);


	float3 newChildF3 = randf3;
	float4 newChild = make_float4(newChildF3.x, newChildF3.y, newChildF3.z, 0.0);

	tmp		= 1.0 / sqrt(newChild.x*newChild.x + newChild.y*newChild.y);
	
	float4 randAxis = make_float4(tmp * newChild.y, -tmp * newChild.x, 0.0, 0.0);
	//float4 randAxis = make_float4(0.0, tmp * newChild.z, -tmp * newChild.y, 0.0);

	rand = Rand2(rand);
	float rand_theta = (-1.0+2.0*(float)rand/RAND2_MAX) * M_PI;//-PI~PI

	float4 newAxis = calRotateQuaternionPos(randAxis, newChildF3, rand_theta);

	dSubBlockChild[index]	= newChild;
	dSubBlockAxis[index]	= newAxis;
	dSubRand[index]			= rand;

}

//---------------------------------------------------------------------------

/*!
 * サブ粒子の更新
 * @param[in]	dSubPosP		サブ粒子の親の位置
 * @param[in]	dSubPosC1		サブ粒子の子1の位置
 * @param[in]	dSubPosC2		サブ粒子の子2の位置
 * @param[in]	dSubChild		サブ粒子の子1への単位ベクトル
 * @param[in]	dSubAxis		サブ粒子の軸
 * @param[in]	dSubEt			サブ粒子のエネルギースペクトル(スケールは直径)
 * @param[in]	dSubRand		uint乱数
 * @param[in]	ratio_et		エネルギースペクトルに掛ける係数
 * @param[in]	radius_sub		このサブ粒子の半径
 * @param[in]	nprts	パーティクル数
 * @param[in]	dt				微少時間
 */
__global__
void updateSubParticle(	float4* dSubPosP,
						float4* dSubPosC1,
						float4* dSubPosC2,
						float4* dSubChild,
						float4* dSubAxis,
						float*	dSubEt,
						uint*	dSubRand,
						float	ratio_et,
						float	radius_sub,
						uint	nprts, 
						uint	uMaxParticles, 
						float	dt)
{
	uint index = __umul24(blockIdx.x,blockDim.x)+threadIdx.x;
	if(index >= nprts) return;

	float4 posP		= dSubPosP[index];
	float4 child	= dSubChild[index];
	float4 axis		= dSubAxis[index];
	float  et		= dSubEt[index];
	uint   pre_rand	= dSubRand[index];


	float radiusC = radius_sub * POW_2_M1D3;
	float etP	= et * ratio_et;
	float etC	= etP * POW_2_M5D9;

	uint rand	= Rand2(pre_rand);
	float randF = -1.0 + 2.0 * (float)rand /RAND2_MAX;//65535.0;

	float child_theta = sqrt(2.0 * etC) * dt / radiusC;
	float axis_theta =  RATIO_AXIS_ANGLE * randF * child_theta ;
	
	float3 axisF3		= make_float3(axis.x,axis.y,axis.z);
	float4 newChild		= calRotateQuaternionPos(child,axisF3,child_theta);
	float3 newChildF3	= make_float3(newChild.x, newChild.y, newChild.z);
	float4 newAxis		= calRotateQuaternionPos(axis,newChildF3,axis_theta);

	float4 posC1	= posP + radiusC * newChild;
	float4 posC2	= posP - radiusC * newChild;


	dSubChild[index]	= newChild;
	dSubAxis[index]		= newAxis;
	dSubRand[index]		= rand;
	dSubPosC1[index]	= posC1;
	dSubPosC2[index]	= posC2;

}

/*!
 * 
 * @param[in]	dSubRand		uint乱数
 * @param[in]	nprts	パーティクル数
 */
__global__
void initSubRand(	uint*	dSubRand,
					uint nprts)		
{
	uint index = __umul24(blockIdx.x, blockDim.x)+threadIdx.x;
	if(index >= nprts) return;

	uint rand = Rand2(index);

	dSubRand[index] = rand;
}

/*!
 * 
 * @param[in]	dSubRand		uint乱数
 * @param[in]	nprts	パーティクル数
 */
__global__
void addSubRand(uint* dSubRand,
				uint  nprts, 
				uint  uStart, 
				uint  uCount)
{
	uint index = __umul24(blockIdx.x, blockDim.x)+threadIdx.x;
	//if(index >= uCount || uStart+index >= nprts) return;
	if(index < uStart || index >= uStart+uCount) return;

	uint rand = Rand2(uStart+index);

	dSubRand[uStart+index] = rand;
}


/*!
 * サブ粒子(レベル0)の位置更新
 * @param[in]	dSubPos	サブ粒子の絶対座標へのアドレス
 * @param[in]	dPos	粒子の絶対座標へのアドレス
 * @param[in]	nprts パーティクル数
 */
__global__
void addSubPos(float4* dSubPos,
			   float4* dPos,
			   uint nprts, 
			   uint uStart, 
			   uint uCount)
{
	uint index = __umul24(blockIdx.x,blockDim.x)+threadIdx.x;
	//if(index >= uCount || uStart+index >= nprts) return;
	if(index < uStart || index >= uStart+uCount) return;

	float4 posData = dPos[uStart+index];

	//float4 pos = make_float4(posData.x, posData.y, posData.z, 0.0);

	dSubPos[uStart+index] = posData;
} 


/*!
 * サブ粒子の追加
 * @param[in] dSubPosP		サブ粒子の親の位置
 * @param[in] dSubPosC1		サブ粒子の子1の位置
 * @param[in] dSubPosC2		サブ粒子の子2の位置
 * @param[in] dSubChild		サブ粒子の子1への単位ベクトル
 * @param[in] dSubAxis		サブ粒子の軸
 * @param[in] dSubRand		uint乱数
 * @param[in] radius_sub	このサブ粒子の半径
 * @param[in] nprts	パーティクル数
 * @param[in] uMaxParticles	最大パーティクル数
 * @param[in] uStart		追加開始インデックス
 * @param[in] uCount		追加数
 */
__global__
void addSubParticle(float4* dSubPosP,
					float4* dSubPosC1,
					float4* dSubPosC2,
					float4* dSubChild,
					float4* dSubAxis,
					uint*	dSubRand,
					float radius_sub,
					uint nprts, 
					uint uMaxParticles, 
					uint uStart, 
					uint uCount)		
{
	uint index = __umul24(blockIdx.x,blockDim.x)+threadIdx.x;
	//if(index >= uCount || uStart+index >= nprts) return;
	if(index < uStart || index >= uStart+uCount) return;

	//index += uStart;

	float4 posP		= dSubPosP[index];
	uint rand = dSubRand[index];

	float3 randf3 = make_float3(0.0,0.0,0.0);
	float tmp = 0.0;
	float radiusC = radius_sub * POW_2_M1D3;


	while(randf3.x == 0.0){
		rand = Rand2(rand);
		randf3.x = -1.0+2.0*(float)rand/RAND2_MAX;
	}
	while(randf3.y == 0.0){
		rand = Rand2(rand);
		randf3.y = -1.0+2.0*(float)rand/RAND2_MAX;
	}
	while(randf3.z == 0.0){
		rand = Rand2(rand);
		randf3.z = -1.0+2.0*(float)rand/RAND2_MAX;
	}

	tmp = 1.0f/sqrt(randf3.x*randf3.x + randf3.y*randf3.y + randf3.z*randf3.z);
	randf3 = make_float3(tmp*randf3.x, tmp*randf3.y, tmp*randf3.z);


	float3 newChildF3 = randf3;
	//float3 newChildF3 = make_float3(0, 1, 0);
	float4 newChild = make_float4(newChildF3.x, newChildF3.y, newChildF3.z, 0.0);

	tmp		= 1.0 / sqrt(newChild.x*newChild.x + newChild.y*newChild.y);
	
	float4 randAxis = make_float4(tmp * newChild.y, -tmp * newChild.x, 0.0, 0.0);

	rand = Rand2(rand);
	float rand_theta = (-1.0+2.0*(float)rand/RAND2_MAX) * M_PI;//-PI~PI

	float4 newAxis = calRotateQuaternionPos(randAxis, newChildF3, rand_theta);


	float4 posC1	= posP + radiusC * newChild;
	float4 posC2	= posP - radiusC * newChild;

	dSubChild[index]	= newChild;
	dSubAxis[index]		= newAxis;
	dSubPosC1[index]	= posC1;
	dSubPosC2[index]	= posC2;
	dSubRand[index]		= rand;

}

/*!
 * サブ粒子の更新
 * @param[in]	dSubPosP		サブ粒子の親の位置
 * @param[in]	dSubPosC1		サブ粒子の子1の位置
 * @param[in]	dSubPosC2		サブ粒子の子2の位置
 * @param[in]	dSubChild		サブ粒子の子1への単位ベクトル
 * @param[in]	dSubAxis		サブ粒子の軸
 * @param[in]	dSubRand		uint乱数
 * @param[in]	radius_sub		このサブ粒子の半径
 * @param[in]	nprts	パーティクル数
 * @param[in]	uMaxParticles	最大パーティクル数
 */
__global__
void initSubParticle(	float4* dSubPosP,
						float4* dSubPosC1,
						float4* dSubPosC2,
						float4* dSubChild,
						float4* dSubAxis,
						uint*	dSubRand,
						float radius_sub,
						uint nprts, 
						uint uMaxParticles)		
{
	uint index = __umul24(blockIdx.x,blockDim.x)+threadIdx.x;
	if(index >= nprts) return;


	float4 posP		= dSubPosP[index];
	uint rand = dSubRand[index];

	float3 randf3 = make_float3(0.0,0.0,0.0);
	float tmp = 0.0;
	float radiusC = radius_sub * POW_2_M1D3;


	while(randf3.x == 0.0){
		rand = Rand2(rand);
		randf3.x = -1.0+2.0*(float)rand/RAND2_MAX;
	}
	while(randf3.y == 0.0){
		rand = Rand2(rand);
		randf3.y = -1.0+2.0*(float)rand/RAND2_MAX;
	}
	while(randf3.z == 0.0){
		rand = Rand2(rand);
		randf3.z = -1.0+2.0*(float)rand/RAND2_MAX;
	}

	tmp = 1.0f/sqrt(randf3.x*randf3.x + randf3.y*randf3.y + randf3.z*randf3.z);
	randf3 = make_float3(tmp*randf3.x, tmp*randf3.y, tmp*randf3.z);


	float3 newChildF3 = randf3;
	float4 newChild = make_float4(newChildF3.x, newChildF3.y, newChildF3.z, 0.0);


	tmp		= 1.0 / sqrt(newChild.x*newChild.x + newChild.y*newChild.y);
	
	float4 randAxis = make_float4(tmp * newChild.y, -tmp * newChild.x, 0.0, 0.0);

	rand = Rand2(rand);
	float rand_theta = (-1.0+2.0*(float)rand/RAND2_MAX) * M_PI;//-PI~PI

	float4 newAxis = calRotateQuaternionPos(randAxis, newChildF3, rand_theta);


	float4 posC1	= posP + radiusC * newChild;
	float4 posC2	= posP - radiusC * newChild;

	dSubChild[index]	= newChild;
	dSubAxis[index]		= newAxis;
	dSubPosC1[index]	= posC1;
	dSubPosC2[index]	= posC2;
	dSubRand[index]		= rand;

}


/*!
 * エネルギーに基づき使用する子パーティクルを選定
 * @param[out] subCell 子パーティクル情報
 * @param[in] dEt レベル0のパーティクルのエネルギー値
 * @param[in] dSubPos 子パーティクルの位置
 * @param[in] ratio_e0 レベル0のパーティクルのエネルギー値に掛ける係数
 */
__global__
void setSubUnsortArray(rxSubParticleCell subCell, float* dEt, float4* dSubPos, float ratio_et0)
{
	uint index = __umul24(blockIdx.x,blockDim.x)+threadIdx.x;
	if(index >= subCell.uNumParticles) return;

	uint sub_i, subIndex, subIndex1;
	float ratioL[MAX_SUB_LEVEL+1];	// 各レベルの影響係数
	float ratio[MAX_SUB_NUM];
	float radius[MAX_SUB_NUM];
	//float4 posData[MAX_SUB_NUM];

	float ets = dEt[index];			// パーティクル乱流エネルギー
	float et_cri = subCell.fEtcri;	// エネルギー基準値

	// 各レベルの影響係数を初期化
	for(uint level = 0; level <= MAX_SUB_LEVEL; level++){
		ratioL[level] = 0.0;
	}

	// レベル0のエネルギー
	float et0 = ratio_et0*ets;

	//比率決定
	float setLevel = log(et_cri/et0)*INV_LOG_2_M5D9;

	// 各レベルの影響係数を計算
	if(setLevel <= 0.0f){	// レベル0のパーティクルのみ
		ratioL[0] = 1.0f;
	}
	else if(setLevel >= (float)MAX_SUB_LEVEL){	// 最下層の子パーティクルのみ
		ratioL[MAX_SUB_LEVEL] = powf(2.0f, -(float)MAX_SUB_LEVEL);
	}
	else {	// 子パーティクル間(leveldとleveld+1の間)
		uint leveld = (uint)setLevel;
		ratioL[leveld]   = (setLevel-(float)leveld)*powf(2.0f, -(float)leveld) ;
		ratioL[leveld+1] = (1.0-ratioL[leveld])*powf(2.0f, -(float)(1+leveld)) ;
	}

	// 影響係数と半径をローカルメモリに格納
	sub_i = 0;
	for(uint level = 0; level <= MAX_SUB_LEVEL; ++level){
		for(uint lIndex = 0; lIndex < subCell.uSubNumEach[level]; ++lIndex){
			ratio[sub_i] = ratioL[level];
			radius[sub_i] = subCell.fSubRad[level];
			sub_i++;
		}
	}

	// サブパーティクル情報を格納
	subIndex  = index;
	subIndex1 = index;
	for(uint i = 0; i < MAX_SUB_NUM; ++i){
		subCell.dSubUnsortRad[subIndex]	= radius[i];
		subCell.dSubUnsortRat[subIndex]	= ratio[i];
		subCell.dSubUnsortPos[subIndex] = dSubPos[subIndex1];

		if(ratio[i] > 1e-10){
			subCell.dSubOcc[subIndex] = 1;
		}

		subIndex  += subCell.uNumParticles;
		subIndex1 += subCell.uMaxParticles;
	}
	
}

/*!
 * サブパーティクルを詰める 
 */
__global__ 
void compactSubParticles(rxSubParticleCell subCell, uint num)
{
	// インデックス
	uint blockId = __mul24(blockIdx.y, gridDim.x)+blockIdx.x;
	uint i = __mul24(blockId, blockDim.x)+threadIdx.x;

	if(subCell.dSubOcc[i] && (i < num)) {
		uint idx = subCell.dSubOccScan[i];
		subCell.dSubSortedPos[idx] = subCell.dSubUnsortPos[i];
		subCell.dSubSortedRad[idx] = subCell.dSubUnsortRad[i];
		subCell.dSubSortedRat[idx] = subCell.dSubUnsortRat[i];
	}
}


/*!
 * レベルごとにサブパーティクルデータをソートして，
 * ハッシュ内の各セルの最初のアドレスを検索
 * @param[in] cell パーティクルグリッドデータ
 * @param[in] oldPos パーティクル位置
 * @param[in] oldVel パーティクル速度
 */
__global__
void reorderDataAndFindCellStartF4F1F1(rxSubParticleCell cell)
{
	extern __shared__ uint sharedHash[];	// blockSize+1 elements
	uint index = __umul24(blockIdx.x,blockDim.x)+threadIdx.x;
	
	uint hash;
	if(index < cell.uSubNumValidParticles){
		hash = cell.dSubGridParticleHash[index];	// ハッシュ値

		sharedHash[threadIdx.x+1] = hash;	// ハッシュ値をシェアードメモリに格納

		if(index > 0 && threadIdx.x == 0){
			// 各シェアードメモリの最初は隣のパーティクルのハッシュ値を格納
			sharedHash[0] = cell.dSubGridParticleHash[index-1];
		}
	}

	__syncthreads();
	
	if(index < cell.uSubNumValidParticles && hash < cell.uSubNumCells){
		if(index == 0 || hash != sharedHash[threadIdx.x]){
			cell.dSubCellStart[hash] = index;
			if(index > 0){
				cell.dSubCellEnd[sharedHash[threadIdx.x]] = index;
			}
		}

		if(index == cell.uSubNumValidParticles-1){
			cell.dSubCellEnd[hash] = index+1;
		}

		// 位置と速度のデータを並び替えた物を格納
		uint sortedIndex = cell.dSubSortedIndex[index];
		cell.dSubSortedPos[index] = FETCHC(dSubUnsortPos, sortedIndex);
		cell.dSubSortedRad[index] = FETCHC(dSubUnsortRad, sortedIndex);
		cell.dSubSortedRat[index] = FETCHC(dSubUnsortRat, sortedIndex);
	}
}

/*!
 * 与えられたセル内のパーティクルとの距離から密度を計算
 * @param[in] gridPos グリッド位置
 * @param[in] pos0 計算座標
 * @param[in] cell パーティクルグリッドデータ
 * @return セル内のパーティクルから計算した密度値
 */
__device__
float calSubDensityCellG(int3 gridPos,
						 float3 pos0, 
						 rxSubParticleCell cell)
{
	uint gridHash = calcGridHash(gridPos);

	// セル内のパーティクルのスタートインデックス
	uint startIndex = FETCHC(dSubCellStart, gridHash);

	float d = 0.0f;
	if(startIndex != 0xffffffff){	// セルが空でないかのチェック
		// セル内のパーティクルで反復
		uint endIndex = FETCHC(dSubCellEnd, gridHash);

		for(uint j = startIndex; j < endIndex; ++j){
			//if(j != index){
				float3 pos1 = make_float3(FETCHC(dSubSortedPos, j));
				float  rad  = FETCHC(dSubSortedRad, j);

				float3 rij = pos0-pos1;
				float h = cell.fSubEffectiveFactor * rad;
				float r = length(rij);

				if(r <= h){
					float rat = FETCHC(dSubSortedRat, j);
					float q = h*h-r*r;

					//d += rat*params.Mass*params.Wpoly6*q*q*q;
					d += rat*params.Mass*315.0f/64.0f/3.14f/powf(h,9.0f)*q*q*q;
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
void calSubGridDensity(	float* GridD, 
						rxSubParticleCell subCell, 
						uint3 gnum, 
						float3 gmin, 
						float3 glen)
{
	uint blockId = __mul24(blockIdx.y, gridDim.x)+blockIdx.x;
	uint i = __mul24(blockId, blockDim.x)+threadIdx.x;

	uint3 gridPos = calcGridPosU(i, gnum);
	//uint3 gridPos = calcGridPosU3(i, gnum.x, gnum.y, gnum.z);

	if(gridPos.x < gnum.x && gridPos.y < gnum.y && gridPos.z < gnum.z){
		float3 gpos;
		gpos.x = gmin.x+(gridPos.x)*glen.x;
		gpos.y = gmin.y+(gridPos.y)*glen.y;
		gpos.z = gmin.z+(gridPos.z)*glen.z;

		float d = 0.0f;

		int3 pgpos = calcGridPos(gpos);
		
		float3 gpos_max,gpos_min;
		int3 pgpos_max,pgpos_min;

		gpos_max.x = gmin.x+(gridPos.x)*glen.x+params.EffectiveRadius;
		gpos_max.y = gmin.y+(gridPos.y)*glen.y+params.EffectiveRadius;
		gpos_max.z = gmin.z+(gridPos.z)*glen.z+params.EffectiveRadius;

		gpos_min.x = gmin.x+(gridPos.x)*glen.x-params.EffectiveRadius;
		gpos_min.y = gmin.y+(gridPos.y)*glen.y-params.EffectiveRadius;
		gpos_min.z = gmin.z+(gridPos.z)*glen.z-params.EffectiveRadius;

		pgpos_max = calcGridPos(gpos_max);
		pgpos_min = calcGridPos(gpos_min);

		for(int z = pgpos_min.z; z <= pgpos_max.z; ++z){
			for(int y = pgpos_min.y; y <= pgpos_max.y; ++y){
				for(int x = pgpos_min.x; x <= pgpos_max.x; ++x){
					int3 neighbourPos = make_int3(x, y, z);

					float minDis = calcGridPosDisMin(gpos, neighbourPos);
					if(minDis < params.EffectiveRadius){
						d += calSubDensityCellG(neighbourPos, gpos, subCell);
					}
				}
			}
		}
		
		GridD[gridPos.x+gridPos.y*gnum.x+gridPos.z*gnum.x*gnum.y] = d;
	}

}


//
__global__
void checkNumUintData(uint*	dUintData,
						uint*	dNumValidUintData,
						uint	uUintMax,
						uint	uNumUintData)
{
	extern __shared__ uint sharedUintData[];	// blockSize+1 elements

	uint index = __umul24(blockIdx.x,blockDim.x)+threadIdx.x;
	
	uint uintIndex;
	if(index < uNumUintData){
		uintIndex = dUintData[index];	// ハッシュ値

		sharedUintData[threadIdx.x+1] = uintIndex;	// ハッシュ値をシェアードメモリに格納

		if(index > 0 && threadIdx.x == 0){
			// 各シェアードメモリの最初は隣のパーティクルのハッシュ値を格納
			sharedUintData[0] = dUintData[index-1];
		}
	}

	__syncthreads();
	
	if(index < uNumUintData){
		if(index == 0){
			if(sharedUintData[1] >= uUintMax){
				dNumValidUintData[0] = 0;
			}
		}	
		//else if(index == uNumF4Data-1){
		//}
		else if( sharedUintData[threadIdx.x] < uUintMax && sharedUintData[threadIdx.x+1] >= uUintMax ){
			uint validNum = dNumValidUintData[0];
			if(validNum > index) dNumValidUintData[0] = index;
		}
	}

}






#endif // #ifndef _RX_TURBULENCE_KERNEL_CU_



