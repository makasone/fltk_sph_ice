/*! 
  @file rx_ssm_kernel.cu
	
  @brief SSM法によるメッシュ生成
		- M. Muller et al., Screen space meshes, SCA2007, 2007. 

  @author Makoto Fujisawa
  @date 2011-05
*/
// FILE --rx_ssm_kernel.cu--

#ifndef _RX_SSM_KERNEL_CU_
#define _RX_SSM_KERNEL_CU_


//-----------------------------------------------------------------------------
// インクルードファイル
//-----------------------------------------------------------------------------
#include "rx_cu_common.cu"


//-----------------------------------------------------------------------------
// 変数
//-----------------------------------------------------------------------------
// テクスチャに格納したテーブル
texture<uint, 1, cudaReadModeElementType> g_texSSMeshTable;
texture<uint, 1, cudaReadModeElementType> g_texSSEdgeTable;
texture<uint, 1, cudaReadModeElementType> g_texSSNodeTable;
texture<uint, 1, cudaReadModeElementType> g_texSSVRotTable;

// シミュレーションパラメータ
__constant__ rxSsmParams g_cSSMParams;

// Binomialフィルタ係数(n_filter = RX_MAX_FILTER_SIZEまで対応)
__constant__ float g_cBinomials[RX_BINOMIALS_SIZE];

// フィルタ用ビット列
__constant__ int BITS[] = {1, 2, 4, 8,  16, 32, 64, 128};


//-----------------------------------------------------------------------------
// MARK:SSM
//-----------------------------------------------------------------------------
/*!
 * float型の配列を初期化
 * @param[in] farray 値を代入したいfloat型の配列
 * @param[in] val 代入する値
 * @param[in] n 配列のサイズ
 */
__global__
void initFloatArray(float* farray, float val, int n)
{
	uint index = __mul24(blockIdx.x,blockDim.x)+threadIdx.x;
	if(index >= n) return;

	farray[index] = val;
}

/*!
 * パーティクル半径の配列を初期化
 * @param[in] rad 出力半径配列(サイズ=n*3)
 * @param[in] val 入力半径配列(サイズ=n)
 * @param[in] n 配列のサイズ
 */
__global__
void initRadArray(float3* rad, float* val, int n)
{
	uint index = __mul24(blockIdx.x,blockDim.x)+threadIdx.x;
	if(index >= n) return;

	float r = val[index];
	rad[index] = make_float3(r, r, r);
}


/*!
 * デプスマップを∞で初期化
 */
__global__
void initDepthMap(float* dmap, int n)
{
	uint index = __mul24(blockIdx.x,blockDim.x)+threadIdx.x;
	if(index >= n) return;

	dmap[index] = 1.0e10;
}


/*!
 * デプスマップの計算
 */
__global__
void calDepthMap(float* dmap, float3* prts_pos, float3* prts_rad, int pnum, float3 tr, matrix4x4 PMV, 
				 float W, float H, float dw, float dh, int nx, int ny)
{
	uint index = __mul24(blockIdx.x,blockDim.x)+threadIdx.x;
	if(index >= pnum) return;	

	float4 x = make_float4(prts_pos[index], 1.0);
	float3 r = prts_rad[index];

	// 投影変換
	float4 xd;
	xd.x = dot(PMV.e[0], x);
	xd.y = dot(PMV.e[1], x);
	xd.z = dot(PMV.e[2], x);
	xd.w = dot(PMV.e[3], x);

	// wで割ることで[-1, 1]の正規化座標系に変換
	float3 xp;
	xp.x = W*(0.5+0.5*xd.x/xd.w);
	xp.y = H*(0.5+0.5*xd.y/xd.w);
	xp.z = xd.z;

	prts_pos[index] = xp;

	// 正規化座標系での半径値
	float3 rp;
	rp.x = r.x*tr.x/xd.w/2;
	rp.y = r.y*tr.y/xd.w/2;
	rp.z = r.z*tr.z;

	prts_rad[index] = rp;

	float rrp = rp.x*rp.x;

	if(xp.x < 0 || xp.x >= W || xp.y < 0 || xp.y >= H){
		return;
	}

	// デプスマップ上でのパーティクルの範囲
	int cen[2];	// パーティクル中心
	cen[0] = int(xp.x/dw)+1;
	cen[1] = int(xp.y/dh)+1;

	int minp[2], maxp[2];
	minp[0] = cen[0]-(rp.x/dw+1);
	minp[1] = cen[1]-(rp.y/dh+1);
	maxp[0] = cen[0]+(rp.x/dw+1);
	maxp[1] = cen[1]+(rp.y/dh+1);

	// 範囲がマップ外にならないようにクランプ
	minp[0] = (minp[0] < 0 ? 0 : (minp[0] > nx ? nx : minp[0]));
	minp[1] = (minp[1] < 0 ? 0 : (minp[1] > ny ? ny : minp[1]));
	maxp[0] = (maxp[0] < 0 ? 0 : (maxp[0] > nx ? nx : maxp[0]));
	maxp[1] = (maxp[1] < 0 ? 0 : (maxp[1] > ny ? ny : maxp[1]));

	// パーティクルデプス値更新
	for(int j = minp[1]; j <= maxp[1]; ++j){
		for(int i = minp[0]; i <= maxp[0]; ++i){
			float rr = (i*dw-xp.x)*(i*dw-xp.x)+(j*dh-xp.y)*(j*dh-xp.y);
			if(rr <= rrp){
				float hij = sqrt(1.0-rr/rrp);
				float z = xp.z-rp.z*hij;

				if(z >= 0){
					atomicFloatMin(&dmap[i+j*(nx+1)], z);
				}

				//float zij = dmap[i+j*(nx+1)];
				//if(z >= 0 && z < zij){
				//	dmap[i+j*(nx+1)] = z;
				//}
			}
		}
	}
}


/*!
 * デプスマップの計算
 */
__global__
void calDepthMap(float* dmap, float4* prts_pos, float3* prts_rad, int pnum, float3 tr, matrix4x4 PMV, 
				 float W, float H, float dw, float dh, int nx, int ny)
{
	uint index = __mul24(blockIdx.x,blockDim.x)+threadIdx.x;
	if(index >= pnum) return;	

	float4 x = prts_pos[index];
	x.w = 1.0;
	float3 r = prts_rad[index];

	// 投影変換
	float4 xd;
	xd.x = dot(PMV.e[0], x);
	xd.y = dot(PMV.e[1], x);
	xd.z = dot(PMV.e[2], x);
	xd.w = dot(PMV.e[3], x);

	// wで割ることで[-1, 1]の正規化座標系に変換
	float4 xp;
	xp.x = W*(0.5+0.5*xd.x/xd.w);
	xp.y = H*(0.5+0.5*xd.y/xd.w);
	xp.z = xd.z;
	xp.w = xd.w;

	prts_pos[index] = xp;

	// 正規化座標系での半径値
	float3 rp;
	rp.x = r.x*tr.x/xd.w/2;
	rp.y = r.y*tr.y/xd.w/2;
	rp.z = r.z*tr.z;

	prts_rad[index] = rp;

	float rrp = rp.x*rp.x;

	if(xp.x < 0 || xp.x >= W || xp.y < 0 || xp.y >= H){
		return;
	}

	// デプスマップ上でのパーティクルの範囲
	int cen[2];	// パーティクル中心
	cen[0] = int(xp.x/dw)+1;
	cen[1] = int(xp.y/dh)+1;

	int minp[2], maxp[2];
	minp[0] = cen[0]-(rp.x/dw+1);
	minp[1] = cen[1]-(rp.y/dh+1);
	maxp[0] = cen[0]+(rp.x/dw+1);
	maxp[1] = cen[1]+(rp.y/dh+1);

	// 範囲がマップ外にならないようにクランプ
	minp[0] = (minp[0] < 0 ? 0 : (minp[0] > nx ? nx : minp[0]));
	minp[1] = (minp[1] < 0 ? 0 : (minp[1] > ny ? ny : minp[1]));
	maxp[0] = (maxp[0] < 0 ? 0 : (maxp[0] > nx ? nx : maxp[0]));
	maxp[1] = (maxp[1] < 0 ? 0 : (maxp[1] > ny ? ny : maxp[1]));

	// パーティクルデプス値更新
	for(int j = minp[1]; j <= maxp[1]; ++j){
		for(int i = minp[0]; i <= maxp[0]; ++i){
			float rr = (i*dw-xp.x)*(i*dw-xp.x)+(j*dh-xp.y)*(j*dh-xp.y);
			if(rr <= rrp){
				float hij = sqrt(1.0-rr/rrp);
				float z = xp.z-rp.z*hij;

				if(z >= 0){
					atomicFloatMin(&dmap[i+j*(nx+1)], z);
				}

				//float zij = dmap[i+j*(nx+1)];
				//if(z >= 0 && z < zij){
				//	dmap[i+j*(nx+1)] = z;
				//}
			}
		}
	}

}

/*!
 * デプスマップからデプス値を取得する(境界処理あり)
 * @param[in] data デプスマップ
 * @param[in] x,y  ピクセル位置
 * @param[in] w,h  デプスマップの大きさ
 * @return デプス値
 */
__device__ 
float GetDepth(float *data, int x, int y, int w, int h)
{
	x = CuClamp(x, 0, w-1);
	y = CuClamp(y, 0, h-1);
	return data[y*w+x];
}


/*!
 * デプスマップを平滑化(x方向)
 *  - シェアードメモリを使用
 *     _____________
 *  r |   :     :   |
 *    |_ _:_____:_ _|
 *    |   |     |   |
 * bh |   |     |   |
 *    |_ _|_____|_ _|
 *  r |   :     :   |
 *    |___:_____:___|
 *      r    bw   r
 *    <----tilew---->
 */
__global__
void smoothDepthMapX(float* dmap, int nx, int ny, int r, int tilew, float zmax)
{
	extern __shared__ float d[];

	// Shared memory上のインデックス取得
	#define SMEM(X, Y) d[(Y)*tilew+(X)]

	// ブロック内のスレッド位置
	int tx = threadIdx.x;
	int ty = threadIdx.y;

	// ブロックサイズ
	int bw = blockDim.x;
	int bh = blockDim.y;

	// グリッド内のインデックス(画像上の位置)
	int x = blockIdx.x*bw+tx;
	int y = blockIdx.y*bh+ty;

	if(x >= nx || y >= ny) return;

	// フィルター範囲の値をシェアードメモリに転送
	// 中心領域
	SMEM(r+tx, r+ty) = GetDepth(dmap, x, y, nx, ny);

	// エッジ領域
	if(threadIdx.x < r){
		SMEM(tx, r+ty)      = GetDepth(dmap, x-r,  y, nx, ny);	// 右
		SMEM(r+bw+tx, r+ty) = GetDepth(dmap, x+bw, y, nx, ny);	// 左
	}
	if(threadIdx.y < r){
		SMEM(r+tx, ty)      = GetDepth(dmap, x, y-r,  nx, ny);
		SMEM(r+tx, r+bh+ty) = GetDepth(dmap, x, y+bh, nx, ny);
	}

	// コーナー領域
	if((threadIdx.x < r) && (threadIdx.y < r)){
		SMEM(tx, ty)           = GetDepth(dmap, x-r,  y-r,  nx, ny);
		SMEM(tx, r+bh+ty)      = GetDepth(dmap, x-r,  y+bh, nx, ny);
		SMEM(r+bw+tx, ty)      = GetDepth(dmap, x+bh, y-r,  nx, ny);
		SMEM(r+bw+tx, r+bh+ty) = GetDepth(dmap, x+bw, y+bh, nx, ny);
	}

	// 全スレッドの転送が終わるまで待つ
	__syncthreads();


	// シェアードメモリ内でのグリッドインデックス
	int i = r+tx;
	int j = r+ty;

	// x方向フィルタ
	float d0 = SMEM(i, j);
	if(d0 < 0.99e10){	// != ∞のピクセルのみに適用
		// 周囲rグリッドで同一レイヤーに属するグリッド数を調べる
		int n = 0;
		for(int k = 1; k <= r; ++k){
			if(fabs(d0-SMEM(i-k, j)) > zmax || fabs(d0-SMEM(i+k, j)) > zmax){
				break;
			}
			n++;
		}

		// Binomial係数を掛けて積算
		int offset = n*n;
		float new_depth = g_cBinomials[offset+n]*d0;
		for(int k = 1; k <= n; ++k){
			new_depth += g_cBinomials[offset+(n-k)]*SMEM(i-k, j);
			new_depth += g_cBinomials[offset+(n+k)]*SMEM(i+k, j);
		}

		dmap[x+y*nx] = new_depth;
	}
}

/*!
 * デプスマップを平滑化(y方向)
 *  - シェアードメモリを使用
 *     _____________
 *  r |   :     :   |
 *    |_ _:_____:_ _|
 *    |   |     |   |
 * bh |   |     |   |
 *    |_ _|_____|_ _|
 *  r |   :     :   |
 *    |___:_____:___|
 *      r    bw   r
 *    <----tilew---->
 */
__global__
void smoothDepthMapY(float* dmap, int nx, int ny, int r, int tilew, float zmax)
{
	extern __shared__ float d[];

	// Shared memory上のインデックス取得
	#define SMEM(X, Y) d[(Y)*tilew+(X)]

	// ブロック内のスレッド位置
	int tx = threadIdx.x;
	int ty = threadIdx.y;

	// ブロックサイズ
	int bw = blockDim.x;
	int bh = blockDim.y;

	// グリッド内のインデックス(画像上の位置)
	int x = blockIdx.x*bw+tx;
	int y = blockIdx.y*bh+ty;

	if(x >= nx || y >= ny) return;

	// 更新したデプス値を再度シェアードメモリに転送
	// 中心領域
	SMEM(r+tx, r+ty) = GetDepth(dmap, x, y, nx, ny);

	// エッジ領域
	if(threadIdx.x < r){
		SMEM(tx, r+ty)      = GetDepth(dmap, x-r,  y, nx, ny);	// 右
		SMEM(r+bw+tx, r+ty) = GetDepth(dmap, x+bw, y, nx, ny);	// 左
	}
	if(threadIdx.y < r){
		SMEM(r+tx, ty)      = GetDepth(dmap, x, y-r,  nx, ny);
		SMEM(r+tx, r+bh+ty) = GetDepth(dmap, x, y+bh, nx, ny);
	}

	// コーナー領域
	if((threadIdx.x < r) && (threadIdx.y < r)){
		SMEM(tx, ty)           = GetDepth(dmap, x-r,  y-r,  nx, ny);
		SMEM(tx, r+bh+ty)      = GetDepth(dmap, x-r,  y+bh, nx, ny);
		SMEM(r+bw+tx, ty)      = GetDepth(dmap, x+bh, y-r,  nx, ny);
		SMEM(r+bw+tx, r+bh+ty) = GetDepth(dmap, x+bw, y+bh, nx, ny);
	}

	// 全スレッドの転送が終わるまで待つ
	__syncthreads();


	// シェアードメモリ内でのグリッドインデックス
	int i = r+tx;
	int j = r+ty;

	// y方向フィルタ
	float d0 = SMEM(i, j);
	if(d0 < 0.99e10){	// != ∞のピクセルのみに適用
		// 周囲rグリッドで同一レイヤーに属するグリッド数を調べる
		int n = 0;
		for(int k = 1; k <= r; ++k){
			if(fabs(d0-SMEM(i, j-k)) > zmax || fabs(d0-SMEM(i, j+k)) > zmax){
				break;
			}
			n++;
		}

		// Binomial係数を掛けて積算
		int offset = n*n;
		float new_depth = g_cBinomials[offset+n]*d0;
		for(int k = 1; k <= n; ++k){
			new_depth += g_cBinomials[offset+(n-k)]*SMEM(i, j-k);
			new_depth += g_cBinomials[offset+(n+k)]*SMEM(i, j+k);
		}

		dmap[x+y*nx] = new_depth;
	}

}

/*!
 * x方向輪郭エッジを抽出
 */
__global__
void detectSilhouetteEdgeX(float* dmap, int nx, int ny, float dw, float dh, float zmax, rxSSEdgeG *dedge, uint *dedgesil)
{
	// グリッド内のインデックス
	int x = blockIdx.x*blockDim.x+threadIdx.x;
	int y = blockIdx.y*blockDim.y+threadIdx.y;

	if(x >= nx || y >= ny) return;

	rxSSEdgeG e;
	e.x0 = make_float3(dw*x, dh*y, dmap[x+y*(nx+1)]);
	e.x1 = make_float3(dw*(x+1), dh*y, dmap[(x+1)+y*(nx+1)]);
	e.depth = 0.5*(e.x0.z+e.x1.z);
	e.front_vertex = -1;
	e.dx = -1.0;

	int silhouette = 0;
	if(fabs(e.x0.z - e.x1.z) > zmax){
		silhouette = 1;
	}
	e.silhouette = silhouette;

	dedge[x+y*nx] = e;
	dedgesil[x+y*nx] = silhouette;
}

/*!
 * y方向輪郭エッジを抽出
 */
__global__
void detectSilhouetteEdgeY(float* dmap, int nx, int ny, float dw, float dh, float zmax, rxSSEdgeG *dedge, uint *dedgesil)
{
	// グリッド内のインデックス
	int x = blockIdx.x*blockDim.x+threadIdx.x;
	int y = blockIdx.y*blockDim.y+threadIdx.y;

	if(x >= nx || y >= ny) return;

	rxSSEdgeG e;
	e.x0 = make_float3(dw*x, dh*y, dmap[x+y*nx]);
	e.x1 = make_float3(dw*x, dh*(y+1), dmap[x+(y+1)*nx]);
	e.depth = 0.5*(e.x0.z+e.x1.z);
	e.front_vertex = -1;
	e.dx = -1.0;

	int silhouette = 0;
	if(fabs(e.x0.z - e.x1.z) > zmax){
		silhouette = 1;
	}
	e.silhouette = silhouette;

	dedge[x+y*nx] = e;
	dedgesil[x+y*nx] = silhouette;
}

/*!
 * 輪郭エッジ情報を詰める 
 * @param[out] compacted_edges 輪郭エッジを詰めた配列
 * @param[in] silhouette 輪郭エッジ情報(輪郭エッジ:1, それ以外:0)
 * @param[in] silhouette_Scan 輪郭エッジ情報から作成したPrefix Sum(Scan)
 * @param[in] edges 全エッジ
 * @param[in] num 総エッジ数
 */
__global__ 
void compactSilhouetteEdges(rxSSEdgeG *compacted_edges, uint *silhouette, uint *silhouette_scan, rxSSEdgeG *edges, uint num)
{
	// インデックス
	uint blockId = __mul24(blockIdx.y, gridDim.x)+blockIdx.x;
	uint i = __mul24(blockId, blockDim.x)+threadIdx.x;

	if(silhouette[i] && (i < num)) {
		compacted_edges[silhouette_scan[i]] = edges[i];
	}
}



/*!
 * 輪郭エッジのfront vertexを算出
 */
__global__
void calFrontEdgeVertex(float3* prts_pos, float3* prts_rad, int pnum, 
						rxSSEdgeG *edges, uint *silhouette, int yoffset, 
						int nx, int ny, float dw, float dh, float W, float H, 
						float3* edge_vertices, uint *edge_vertices_occ)
{
	uint index = __mul24(blockIdx.x,blockDim.x)+threadIdx.x;
	if(index >= pnum) return;	

	float3 xp = prts_pos[index];
	float3 rp = prts_rad[index];

	if(xp.x < 0 || xp.x >= W || xp.y < 0 || xp.y >= H){
		return;
	}

	// デプスマップ上でのパーティクルの範囲
	int cen[2];	// パーティクル中心
	cen[0] = int(xp.x/dw)+1;
	cen[1] = int(xp.y/dh)+1;

	int minp[2], maxp[2];
	minp[0] = cen[0]-(rp.x/dw+1);
	minp[1] = cen[1]-(rp.y/dh+1);
	maxp[0] = cen[0]+(rp.x/dw+1);
	maxp[1] = cen[1]+(rp.y/dh+1);

	// 範囲がマップ外にならないようにクランプ
	minp[0] = CuClamp(minp[0], 0, nx-1);
	minp[1] = CuClamp(minp[1], 0, ny-1);
	maxp[0] = CuClamp(maxp[0], 0, nx-1);
	maxp[1] = CuClamp(maxp[1], 0, ny-1);

	// 範囲内のエッジを調査(x方向)
	for(int j = minp[1]; j <= maxp[1]+1; ++j){
		for(int i = minp[0]; i <= maxp[0]; ++i){
			rxSSEdgeG &e = edges[i+j*(nx)];
			int sil = silhouette[i+j*(nx)];
			if(!sil) continue;
			if(e.depth < xp.z) continue;

			// 円とエッジの交点
			float2 A, B, C, P[2];
			float r, t[2];
			if(e.x0.z <= e.x1.z){
				A = make_float2(e.x0.x, e.x0.y);
				B = make_float2(e.x1.x, e.x1.y);
			}
			else{
				A = make_float2(e.x1.x, e.x1.y);
				B = make_float2(e.x0.x, e.x0.y);
			}
			C = make_float2(xp.x, xp.y);
			r = rp.x;
			int inter = CuLineCircleIntersection(A, B, C, r, P, t);
			if(inter == 1){
				// HACK:要改良
				atomicFloatMax(&e.dx, t[0]);
				if(e.dx == t[0]){
					edge_vertices_occ[i+j*(nx)] = 1;
					edge_vertices[i+j*(nx)] = make_float3(P[0].x, P[0].y, xp.z);
					e.front_vertex = 1;
				}
			}
		}
	}

	// 範囲内のエッジを調査(y方向)
	for(int i = minp[0]; i <= maxp[0]+1; ++i){
		for(int j = minp[1]; j <= maxp[1]; ++j){
			rxSSEdgeG &e = edges[yoffset+i+j*(nx+1)];
			int sil = silhouette[yoffset+i+j*(nx+1)];
			if(!sil) continue;
			if(e.depth < xp.z) continue;

			// 円とエッジの交点
			float2 A, B, C, P[2];
			float r, t[2];
			if(e.x0.z <= e.x1.z){
				A = make_float2(e.x0.x, e.x0.y);
				B = make_float2(e.x1.x, e.x1.y);
			}
			else{
				A = make_float2(e.x1.x, e.x1.y);
				B = make_float2(e.x0.x, e.x0.y);
			}
			C = make_float2(xp.x, xp.y);
			r = rp.x;
			int inter = CuLineCircleIntersection(A, B, C, r, P, t);
			if(inter == 1){
				atomicFloatMax(&e.dx, t[0]);
				if(e.dx == t[0]){
					edge_vertices_occ[yoffset+i+j*(nx+1)] = 1;
					edge_vertices[yoffset+i+j*(nx+1)] = make_float3(P[0].x, P[0].y, xp.z);
					e.front_vertex = 1;
				}
			}
		}
	}

}


						
/*!
 * 輪郭エッジのfront vertexを算出
 */
__global__
void calFrontEdgeVertex(float4* prts_pos, float3* prts_rad, int pnum, 
						rxSSEdgeG *edges, uint *silhouette, int yoffset, 
						int nx, int ny, float dw, float dh, float W, float H, 
						float3* edge_vertices, uint *edge_vertices_occ)
{
	uint index = __mul24(blockIdx.x,blockDim.x)+threadIdx.x;
	if(index >= pnum) return;	

	float3 xp = make_float3(prts_pos[index]);
	float3 rp = prts_rad[index];

	if(xp.x < 0 || xp.x >= W || xp.y < 0 || xp.y >= H){
		return;
	}

	// デプスマップ上でのパーティクルの範囲
	int cen[2];	// パーティクル中心
	cen[0] = int(xp.x/dw)+1;
	cen[1] = int(xp.y/dh)+1;

	int minp[2], maxp[2];
	minp[0] = cen[0]-(rp.x/dw+1);
	minp[1] = cen[1]-(rp.y/dh+1);
	maxp[0] = cen[0]+(rp.x/dw+1);
	maxp[1] = cen[1]+(rp.y/dh+1);

	// 範囲がマップ外にならないようにクランプ
	minp[0] = CuClamp(minp[0], 0, nx-1);
	minp[1] = CuClamp(minp[1], 0, ny-1);
	maxp[0] = CuClamp(maxp[0], 0, nx-1);
	maxp[1] = CuClamp(maxp[1], 0, ny-1);

	// 範囲内のエッジを調査(x方向)
	for(int j = minp[1]; j <= maxp[1]+1; ++j){
		for(int i = minp[0]; i <= maxp[0]; ++i){
			rxSSEdgeG &e = edges[i+j*(nx)];
			int sil = silhouette[i+j*(nx)];
			if(!sil) continue;
			if(e.depth < xp.z) continue;

			// 円とエッジの交点
			float2 A, B, C, P[2];
			float r, t[2];
			if(e.x0.z <= e.x1.z){
				A = make_float2(e.x0.x, e.x0.y);
				B = make_float2(e.x1.x, e.x1.y);
			}
			else{
				A = make_float2(e.x1.x, e.x1.y);
				B = make_float2(e.x0.x, e.x0.y);
			}
			C = make_float2(xp.x, xp.y);
			r = rp.x;
			int inter = CuLineCircleIntersection(A, B, C, r, P, t);
			if(inter == 1){
				// HACK:要改良
				atomicFloatMax(&e.dx, t[0]);
				if(e.dx == t[0]){
					edge_vertices_occ[i+j*(nx)] = 1;
					edge_vertices[i+j*(nx)] = make_float3(P[0].x, P[0].y, xp.z);
					e.front_vertex = 1;
				}
			}
		}
	}

	// 範囲内のエッジを調査(y方向)
	for(int i = minp[0]; i <= maxp[0]+1; ++i){
		for(int j = minp[1]; j <= maxp[1]; ++j){
			rxSSEdgeG &e = edges[yoffset+i+j*(nx+1)];
			int sil = silhouette[yoffset+i+j*(nx+1)];
			if(!sil) continue;
			if(e.depth < xp.z) continue;

			// 円とエッジの交点
			float2 A, B, C, P[2];
			float r, t[2];
			if(e.x0.z <= e.x1.z){
				A = make_float2(e.x0.x, e.x0.y);
				B = make_float2(e.x1.x, e.x1.y);
			}
			else{
				A = make_float2(e.x1.x, e.x1.y);
				B = make_float2(e.x0.x, e.x0.y);
			}
			C = make_float2(xp.x, xp.y);
			r = rp.x;
			int inter = CuLineCircleIntersection(A, B, C, r, P, t);
			if(inter == 1){
				atomicFloatMax(&e.dx, t[0]);
				if(e.dx == t[0]){
					edge_vertices_occ[yoffset+i+j*(nx+1)] = 1;
					edge_vertices[yoffset+i+j*(nx+1)] = make_float3(P[0].x, P[0].y, xp.z);
					e.front_vertex = 1;
				}
			}
		}
	}

}

/*!
 * エッジ頂点を詰める 
 * @param[out] compacted_ev エッジ頂点を詰めた配列
 * @param[in] ev_occ エッジ頂点有無情報(エッジ頂点有り:1, 無し:0)
 * @param[in] ev_occ_scan エッジ頂点有無情報から作成したPrefix Sum(Scan)
 * @param[in] ev エッジ頂点
 * @param[in] edges エッジ
 * @param[in] num 総エッジ数
 */
__global__ 
void compactEdgeVertices(float3 *compacted_ev, uint *ev_occ, uint *ev_occ_scan, float3 *ev, rxSSEdgeG *edges, uint num)
{
	// インデックス
	uint blockId = __mul24(blockIdx.y, gridDim.x)+blockIdx.x;
	uint i = __mul24(blockId, blockDim.x)+threadIdx.x;

	if(ev_occ[i] && (i < num)) {
		compacted_ev[ev_occ_scan[i]] = ev[i];
		edges[i].front_vertex = ev_occ_scan[i];
	}
}



/*!
 * ノード頂点の生成
 */
__global__
void calNodeVertex(float* dmap, int nx, int ny, float dw, float dh, float zmax, float3 *node, uint *node_occ)
{
	// グリッド内のインデックス
	int x = blockIdx.x*blockDim.x+threadIdx.x;
	int y = blockIdx.y*blockDim.y+threadIdx.y;

	if(x >= nx || y >= ny) return;

	int idx = x+y*nx;
	float d = dmap[idx];
	if(d < 0.99e10){	// != ∞のノードに頂点生成
		node[idx] = make_float3(dw*x, dh*y, d);
		node_occ[idx] = 1;
	}
}

/*!
 * ノード頂点を詰める 
 */
__global__ 
void compactNodeVertex(float3 *compacted_nv, uint *nv_occ, uint *nv_occ_scan, float3 *nv, uint num)
{
	// インデックス
	uint blockId = __mul24(blockIdx.y, gridDim.x)+blockIdx.x;
	uint i = __mul24(blockId, blockDim.x)+threadIdx.x;

	if(nv_occ[i] && (i < num)) {
		compacted_nv[nv_occ_scan[i]] = nv[i];
	}
}


/*!
 * ノード頂点インデックスをメッシュグリッドに格納
 */
__global__
void storeNodeVertex(rxSSGridG* mgrid, int nx, int ny, float3 *node_vertices, uint *node_occ, uint *node_occ_scan)
{
	// グリッドインデックス
	int x = blockIdx.x*blockDim.x+threadIdx.x;
	int y = blockIdx.y*blockDim.y+threadIdx.y;

	if(x >= nx || y >= ny) return;

	int idx = x+y*nx;
	rxSSGridG &g = mgrid[idx];

	g.i = x;
	g.j = y;
	g.num_nv = 0;
	g.vrot = 0;
	g.table_index0 = 0;
	g.table_index1 = 0;
	g.mesh_num = 0;

	// 左下ノード
	if(node_occ[x+y*(nx+1)]){
		int vidx = node_occ_scan[x+y*(nx+1)];
		g.node_vrts[0] = vidx;
		g.num_nv++;
		g.node_depth[0] = node_vertices[vidx].z;
	}
	else{
		g.node_vrts[0] = -1;
		g.node_depth[0] = 1e10;
	}

	// 右下ノード
	if(node_occ[(x+1)+y*(nx+1)]){
		int vidx = node_occ_scan[(x+1)+y*(nx+1)];
		g.node_vrts[1] = vidx;
		g.num_nv++;
		g.node_depth[1] = node_vertices[vidx].z;
	}
	else{
		g.node_vrts[1] = -1;
		g.node_depth[1] = 1e10;
	}

	// 左上ノード
	if(node_occ[(x+1)+(y+1)*(nx+1)]){
		int vidx = node_occ_scan[(x+1)+(y+1)*(nx+1)];
		g.node_vrts[2] = vidx;
		g.num_nv++;
		g.node_depth[2] = node_vertices[vidx].z;
	}
	else{
		g.node_vrts[2] = -1;
		g.node_depth[2] = 1e10;
	}

	// 右上ノード
	if(node_occ[x+(y+1)*(nx+1)]){
		int vidx = node_occ_scan[x+(y+1)*(nx+1)];
		g.node_vrts[3] = vidx;
		g.num_nv++;
		g.node_depth[3] = node_vertices[vidx].z;
	}
	else{
		g.node_vrts[3] = -1;
		g.node_depth[3] = 1e10;
	}
}




/*!
 * 輪郭エッジのback vertexを算出
 */
__global__
void calBackEdgeVertexX(float* dmap, rxSSEdgeG *edges, uint *silhouette, 
						int nx, int ny, float dw, float dh, float W, float H, 
						float3* front_vertices, float3* back_vertices, uint *back_vertices_occ)
{
	int x = blockIdx.x*blockDim.x+threadIdx.x;
	int y = blockIdx.y*blockDim.y+threadIdx.y;

	if(x >= nx || y >= ny) return;
	if(!silhouette[x+y*nx]) return;

	rxSSEdgeG e = edges[x+y*nx];


	// back vertexがあるかどうかをチェック
	if(e.x0.z < 0.99e10 && e.x1.z < 0.99e10){	// エッジ端点が両方とも != ∞ ならば back vertex が存在
		int back_node;		// back layerに属するエッジ端点
		//int nn_back_node;	// 隣接ノード

		// デプス値が大きい方が back layer に属する
		if(e.x0.z > e.x1.z){
			back_node = x+y*(nx+1);
			//nn_back_node = (x == 0 ? x : x-1)+y*(nx+1);
		}
		else{
			back_node = (x+1)+y*(nx+1);
			//nn_back_node = ((x+1) == nx ? x+1 : x+2)+y*(nx+1);
		}

		// back layerに属するエッジ端点とその隣接ノードのデプス値
		float back_node_depth = dmap[back_node];
		//float nn_back_node_depth = dmap[nn_back_node];

		// デプス値を外挿により近似
		float back_depth = back_node_depth;
		//float back_depth = back_node_depth*((2*l-dx)/l)-nn_back_node_depth*((l-dx)/l);

		float3 vrt_pos = front_vertices[x+y*nx];

		// back vertexを設定
		float3 back_vrt;
		back_vrt.x = vrt_pos.x;
		back_vrt.y = vrt_pos.y;
		back_vrt.z = back_depth;

		back_vertices_occ[x+y*nx] = 1;
		back_vertices[x+y*nx] = back_vrt;
	}
}

/*!
 * 輪郭エッジのback vertexを算出
 */
__global__
void calBackEdgeVertexY(float* dmap, rxSSEdgeG *edges, uint *silhouette, 
						int nx, int ny, float dw, float dh, float W, float H, 
						float3* front_vertices, float3* back_vertices, uint *back_vertices_occ, int offset)
{
	int x = blockIdx.x*blockDim.x+threadIdx.x;
	int y = blockIdx.y*blockDim.y+threadIdx.y;

	if(x >= nx || y >= ny) return;
	if(!silhouette[x+y*nx+offset]) return;

	rxSSEdgeG e = edges[x+y*nx+offset];

	// back vertexがあるかどうかをチェック
	if(e.x0.z < 0.99e10 && e.x1.z < 0.99e10){	// エッジ端点が両方とも != ∞ ならば back vertex が存在
		int back_node;		// back layerに属するエッジ端点
		//int nn_back_node;	// 隣接ノード

		// デプス値が大きい方が back layer に属する
		if(e.x0.z > e.x1.z){
			back_node = x+y*nx;
			//nn_back_node = (x == 0 ? x : x-1)+y*nx;
		}
		else{
			back_node = x+(y+1)*nx;
			//nn_back_node = ((x+1) == nx ? x+1 : x+2)+y*nx;
		}

		// back layerに属するエッジ端点とその隣接ノードのデプス値
		float back_node_depth = dmap[back_node];
		//float nn_back_node_depth = dmap[nn_back_node];

		// デプス値を外挿により近似
		float back_depth = back_node_depth;
		//float back_depth = back_node_depth*((2*l-dx)/l)-nn_back_node_depth*((l-dx)/l);

		float3 vrt_pos = front_vertices[x+y*nx+offset];

		// back vertexを設定
		float3 back_vrt;
		back_vrt.x = vrt_pos.x;
		back_vrt.y = vrt_pos.y;
		back_vrt.z = back_depth;

		back_vertices_occ[x+y*nx+offset] = 1;
		back_vertices[x+y*nx+offset] = back_vrt;
	}
}

/*!
 * エッジ頂点を詰める 
 * @param[out] compacted_ev エッジ頂点を詰めた配列
 * @param[in] ev_occ エッジ頂点有無情報(エッジ頂点有り:1, 無し:0)
 * @param[in] ev_occ_scan エッジ頂点有無情報から作成したPrefix Sum(Scan)
 * @param[in] ev エッジ頂点
 * @param[in] num 総エッジ数
 */
__global__ 
void compactBackEdgeVertices(float3 *compacted_ev, uint *ev_occ, uint *ev_occ_scan, float3 *ev, uint num)
{
	// インデックス
	uint blockId = __mul24(blockIdx.y, gridDim.x)+blockIdx.x;
	uint i = __mul24(blockId, blockDim.x)+threadIdx.x;

	if(ev_occ[i] && (i < num)) {
		compacted_ev[ev_occ_scan[i]] = ev[i];
	}
}


/*!
 * エッジ頂点インデックスをメッシュグリッドに格納
 */
__global__
void storeEdgeVertex(rxSSGridG* mgrid, int nx, int ny, int yoffset, int num_nv, int num_nev, 
					 float3 *edge_vrts, uint *edge_occ, uint *edge_occ_scan, 
					 float3 *back_vrts, uint *back_occ, uint *back_occ_scan)
{
	// グリッドインデックス
	int x = blockIdx.x*blockDim.x+threadIdx.x;
	int y = blockIdx.y*blockDim.y+threadIdx.y;

	if(x >= nx || y >= ny) return;

	int idx = x+y*nx;
	rxSSGridG &g = mgrid[idx];

	g.num_ev = 0;
	g.num_bv = 0;
	g.back2 = -1;

	//
	// 前面エッジ頂点
	//
	// 下エッジ
	if(edge_occ[x+y*nx]){
		int eidx = edge_occ_scan[x+y*nx]+num_nv;
		g.edge_vrts[0] = eidx;
		g.num_ev++;
	}
	else{
		g.edge_vrts[0] = -1;
	}

	// 上エッジ
	if(edge_occ[x+(y+1)*nx]){
		int eidx = edge_occ_scan[x+(y+1)*nx]+num_nv;
		g.edge_vrts[2] = eidx;
		g.num_ev++;
	}
	else{
		g.edge_vrts[2] = -1;
	}

	// 左エッジ
	if(edge_occ[yoffset+x+y*(nx+1)]){
		int eidx = edge_occ_scan[yoffset+x+y*(nx+1)]+num_nv;
		g.edge_vrts[3] = eidx;
		g.num_ev++;
	}
	else{
		g.edge_vrts[3] = -1;
	}

	// 右エッジ
	if(edge_occ[yoffset+(x+1)+y*(nx+1)]){
		int eidx = edge_occ_scan[yoffset+(x+1)+y*(nx+1)]+num_nv;
		g.edge_vrts[1] = eidx;
		g.num_ev++;
	}
	else{
		g.edge_vrts[1] = -1;
	}

	//
	// 背面エッジ頂点
	//
	// 下エッジ
	if(back_occ[x+y*nx]){
		int eidx = back_occ_scan[x+y*nx]+num_nev;
		g.back_vrts[0] = eidx;
		g.num_bv++;
	}
	else{
		g.back_vrts[0] = -1;
	}

	// 上エッジ
	if(back_occ[x+(y+1)*nx]){
		int eidx = back_occ_scan[x+(y+1)*nx]+num_nev;
		g.back_vrts[2] = eidx;
		g.num_bv++;
	}
	else{
		g.back_vrts[2] = -1;
	}

	// 左エッジ
	if(back_occ[yoffset+x+y*(nx+1)]){
		int eidx = back_occ_scan[yoffset+x+y*(nx+1)]+num_nev;
		g.back_vrts[3] = eidx;
		g.num_bv++;
	}
	else{
		g.back_vrts[3] = -1;
	}

	// 右エッジ
	if(back_occ[yoffset+(x+1)+y*(nx+1)]){
		int eidx = back_occ_scan[yoffset+(x+1)+y*(nx+1)]+num_nev;
		g.back_vrts[1] = eidx;
		g.num_bv++;
	}
	else{
		g.back_vrts[1] = -1;
	}
	g.back_vrts[4] = -1;
	g.back_vrts[5] = -1;
}




/*!
 * パターン0 内部メッシュ用のテーブルインデックス更新
 * @param[inout] table_index 元のテーブルインデックス(上位4ビットがエッジ頂点,下位4ビットがノード頂点)
 * @param[inout] vrot 頂点ローテーション数
 * @param[in] v グリッド内の頂点インデックス
 * @param[inout] g メッシュを生成するグリッド
 */
__device__ 
void CuUpdateTableIndexE0N4(int &table_index, int &vrot, int v[], rxSSGridG *g)
{
	// グリッド番号により90度回転させる
	table_index -= ((g->i+g->j) & 0x01) ? 1 : 0;
}

/*!
 * パターン1,2,4,8 内部輪郭(輪郭の始点)用のテーブルインデックス更新
 * @param[inout] table_index 元のテーブルインデックス(上位4ビットがエッジ頂点,下位4ビットがノード頂点)
 * @param[inout] vrot 頂点ローテーション数
 * @param[in] v グリッド内の頂点インデックス
 * @param[inout] g メッシュを生成するグリッド
 */
__device__ 
void CuUpdateTableIndexE1(int &table_index, int &vrot, int v[], rxSSGridG *g)
{
	// ノード頂点で back layer に属する物を探す
	// 頂点があるエッジの端点でデプス値が大きい方が back layer
	int kb = 0;
	for(int k = 0; k < 4; ++k){
		if(v[k+4] != -1){
			int k1 = (k == 3 ? 0 : k+1);
			kb = ((g->node_depth[k] > g->node_depth[k1]) ? BITS[k] : BITS[k1]);
			break;
		}
	}

	kb = (~kb & 0x0F);
	

	// バックレイヤーにあるノード頂点ビットを0にする
	table_index = (table_index & 0xF0)+kb;
}

/*!
 * パターン3,5,6,9,10,12 内部輪郭用のテーブルインデックス更新
 * @param[inout] table_index 元のテーブルインデックス(上位4ビットがエッジ頂点,下位4ビットがノード頂点)
 * @param[inout] vrot 頂点ローテーション数
 * @param[in] v グリッド内の頂点インデックス
 * @param[inout] g メッシュを生成するグリッド
 */
__device__ 
void CuUpdateTableIndexE2N4(int &table_index, int &vrot, int v[], rxSSGridG *g)
{
	// back vertexのノード頂点ビットを0としたビット列を生成
	int btable = 0;
	int k0 = 0, k1;
	for(int k = 0; k <= 4; ++k){
		k1 = (k0 == 3 ? 0 : k0+1);
		if(v[k0+4] != -1){	// ノード間にエッジ頂点あり
			btable |= (g->node_depth[k0] > g->node_depth[k1]) ? BITS[k0] : BITS[k1];
		}
		else{
			btable |= ((btable & BITS[k0]) ? BITS[k1] : 0);
		}
		k0++;
		if(k0 == 4) k0 = 0;
	}

	// 外部輪郭と区別するために，内部輪郭の場合ビット列に+2する
	btable = (btable+2 & 0x0F);

	table_index = (table_index & 0xF0)+btable;
}

/*!
 * 7,11,13,14 外部/内部輪郭混在用のテーブルインデックス更新
 * @param[inout] table_index 元のテーブルインデックス(上位4ビットがエッジ頂点,下位4ビットがノード頂点)
 * @param[inout] vrot 頂点ローテーション数
 * @param[in] v グリッド内の頂点インデックス
 * @param[inout] g メッシュを生成するグリッド
 */
__device__ 
void CuUpdateTableIndexE3N23(int &table_index, int &vrot, int v[], rxSSGridG *g, float3 *vrts, float zmax, rxVrtAdd &add)
{
	int btable = 0;	// ノード頂点数2,3の場合，1番目のビットは0(0xxx)
	int pattern = (table_index >> 4);
	int node = tex1Dfetch(g_texSSEdgeTable, pattern*4)-4;
	int R[4];	// 頂点がないエッジから反時計回りにノード頂点番号を並べたリスト
	for(int k = 0; k < 4; ++k){
		R[k] = node++;
		if(node == 4) node = 0;
	}

	int ntable = (table_index & 0x0F);	// ノード頂点ビット列

	// 2番目のビット
	btable |= (((ntable >> R[1]) & 1) ? 4 : 0);

	// 3番目のビット
	btable |= (((ntable >> R[2]) & 1) ? 2 : 0);

	// Rに従いbtableを並び替え
	int btable0 = 0;
	btable0 |= (((ntable >> R[0]) & 1) ? 8 : 0);
	btable0 |= (((ntable >> R[1]) & 1) ? 4 : 0);
	btable0 |= (((ntable >> R[2]) & 1) ? 2 : 0);
	btable0 |= (((ntable >> R[3]) & 1) ? 1 : 0);

	// 並び替えたbtableの下位2ビット-1
	int n0 = (btable0 & 3)-1;
	int n1 = n0+1;

	add.num = 0;
	if(n0 != 1){
		// 追加のback vertexをg_EdgeTable[pattern][1]の位置に追加
		int add_edge = tex1Dfetch(g_texSSEdgeTable, pattern*4+1)-4;	// 追加するエッジの位置
		float3 add_vrt = vrts[g->edge_vrts[add_edge]];		// 追加頂点

		// 追加頂点デプス値
		int ref_edge = tex1Dfetch(g_texSSEdgeTable, pattern*4)-4;		// 追加頂点のデプス値を参照するエッジ頂点
		int ref_node = (btable & 4) ? R[1] : R[2];
		if(fabs(vrts[g->edge_vrts[ref_edge]].z-g->node_depth[ref_node]) > zmax){
			// edge vertexを使用
			add_vrt.z = vrts[g->edge_vrts[ref_edge]].z;	// 追加頂点のデプス値
		}
		else{
			// back vertexを使用
			if(g->back_vrts[ref_edge] == -1){
				ref_edge = tex1Dfetch(g_texSSEdgeTable, pattern*4+2)-4;
			}
			add_vrt.z = vrts[g->back_vrts[ref_edge]].z;	// 追加頂点のデプス値
		}

		add.edge[0] = add_edge;
		add.vrts[0] = add_vrt;
		add.num = 1;

		g->back_vrts[4] = add_edge;
		
		// グリッドの頂点リストに追加
		if(add_vrt.z < vrts[v[add_edge+4]].z){
			add.layer = 0;
		}
		else{
			add.layer = 1;
		}
			
	}

	// 4番目のビット : ノードのデプス値を比較
	btable |= ((g->node_depth[R[n0]] >= g->node_depth[R[n1]]) ? 0 : 1);

	// テーブルインデックスのエッジ頂点部分を更新
	table_index &= 0xF0;
	table_index |= btable;
}

/*!
 * パターン7,11,13,14 内部輪郭用のテーブルインデックス更新
 *  - さらなるback vertexが必要
 * @param[inout] table_index 元のテーブルインデックス(上位4ビットがエッジ頂点,下位4ビットがノード頂点)
 * @param[inout] vrot 頂点ローテーション数
 * @param[in] v グリッド内の頂点インデックス
 * @param[inout] g メッシュを生成するグリッド
 */
__device__ 
void CuUpdateTableIndexE3N4(int &table_index, int &vrot, int v[], rxSSGridG *g, float3 *vrts, float zmax, rxVrtAdd &add)
{
	int btable = 8;	// ノード頂点数4の場合，1番目のビットは1(1xxx)
	int pattern = (table_index >> 4);

	// ノードのデプス値を比べる順番の決定
	int node = tex1Dfetch(g_texSSEdgeTable, pattern*4)-4;
	int R[4];	// 頂点がないエッジから反時計回りにノード頂点番号を並べたリスト
	for(int k = 0; k < 4; ++k){
		R[k] = node++;
		if(node == 4) node = 0;
	}

	// ノードのデプス値の大小でビット列を変更(2-4番目のビット)
	for(int k = 0; k < 3; ++k){
		// R[k] > R[k+1]ならば対応するビットを1にする
		if(g->node_depth[R[k]] > g->node_depth[R[k+1]]){
			btable |= BITS[2-k];
		}
	}

	// テーブルインデックスのエッジ頂点部分を更新
	table_index &= 0xF0;
	table_index |= btable;

	//
	// g_EdgeTable[pattern][1]の位置に頂点を追加
	//
	int add_edge = tex1Dfetch(g_texSSEdgeTable, pattern*4+1)-4;	// 追加するエッジの位置
	float3 add_vrt = vrts[g->edge_vrts[add_edge]];		// 追加頂点

	// メッシュCの頂点から補間でデプス値を求める
	float ref_depths[4];
	// ／＼
	// 2ー3
	// |＼|
	// 0ー1

	ref_depths[0] = g->node_depth[R[3]];
	ref_depths[1] = g->node_depth[R[0]];

	int e2 = tex1Dfetch(g_texSSEdgeTable, pattern*4+2)-4;
	int e3 = tex1Dfetch(g_texSSEdgeTable, pattern*4+0)-4;

	if(fabs(vrts[g->edge_vrts[e2]].z-ref_depths[0]) < zmax){
		ref_depths[2] = vrts[g->edge_vrts[e2]].z;
	}
	else{
		ref_depths[2] = vrts[g->back_vrts[e2]].z;
	}
	if(fabs(vrts[g->edge_vrts[e3]].z-ref_depths[1]) < zmax){
		ref_depths[3] = vrts[g->edge_vrts[e3]].z;
	}
	else{
		ref_depths[3] = vrts[g->back_vrts[e3]].z;
	}

	// 追加頂点のデプス値
	add_vrt.z = 0.5*(ref_depths[2]+ref_depths[3]);

	add.edge[0] = add_edge;
	add.vrts[0] = add_vrt;
	add.num = 1;

	g->back_vrts[4] = add_edge;

	// グリッドの頂点リストに追加
	if(add_vrt.z < vrts[v[add_edge+8]].z){
		if(add_vrt.z < vrts[v[add_edge+4]].z){
			// front vertexとして挿入
			add.layer = 0;
		}
		else{
			// back vertexとして挿入
			add.layer = 1;
		}
	}
	else{
		// back-2 vertexとして挿入
		add.layer = 2;
	}	
}

/*!
 * パターン15 外部輪郭用のテーブルインデックス更新
 * @param[inout] table_index 元のテーブルインデックス(上位4ビットがエッジ頂点,下位4ビットがノード頂点)
 * @param[inout] vrot 頂点ローテーション数
 * @param[in] v グリッド内の頂点インデックス
 * @param[inout] g メッシュを生成するグリッド
 */
__device__ 
void CuUpdateTableIndexE4N2(int &table_index, int &vrot, int v[], rxSSGridG *g)
{
	int ntable = (table_index & 0x0F);	// ノード頂点ビット列
	ntable = (ntable == 5 ? 0 : 15);

	// テーブルインデックスのノード頂点部分を更新
	table_index &= 0xF0;
	table_index |= ntable;
}

/*!
 * パターン15 外部/内部輪郭混在用のテーブルインデックス更新
 * @param[inout] table_index 元のテーブルインデックス(上位4ビットがエッジ頂点,下位4ビットがノード頂点)
 * @param[inout] vrot 頂点ローテーション数
 * @param[in] v グリッド内の頂点インデックス
 * @param[inout] g メッシュを生成するグリッド
 */
__device__ 
void CuUpdateTableIndexE4N3(int &table_index, int &vrot, int v[], rxSSGridG *g)
{
	int btable = 0;	// ノード頂点数3の場合，1番目のビットは0(0xxx)
	int ntable = (table_index & 0x0F);	// ノード頂点ビット列

	// 頂点がないノード
	int zero_node = log((double)(~ntable & 0x0F))/log(2.0);

	// ノードのデプス値の大小でビット列を変更(2-4番目のビット)
	for(int k = 0; k < 3; ++k){
		int k0 = tex1Dfetch(g_texSSNodeTable, zero_node*8+k);
		int k1 = tex1Dfetch(g_texSSNodeTable, zero_node*8+k+1);

		// k0 > k1ならば対応するビットを1にする
		if(g->node_depth[k0] > g->node_depth[k1]){
			btable |= BITS[2-k];
		}
	}

	// 頂点ローテーション
	vrot = zero_node;

	// テーブルインデックスのノード頂点部分を更新
	table_index &= 0xF0;
	table_index |= btable;
}


/*!
 * パターン15 内部輪郭用のテーブルインデックス更新
 *  - 追加のback vertexが2つ必要
 * @param[inout] table_index 元のテーブルインデックス(上位4ビットがエッジ頂点,下位4ビットがノード頂点)
 * @param[inout] vrot 頂点ローテーション数
 * @param[in] v グリッド内の頂点インデックス
 * @param[inout] g メッシュを生成するグリッド
 */
__device__ 
void CuUpdateTableIndexE4N4(int &table_index, int &vrot, int v[], rxSSGridG *g, float3 *vrts, float zmax, rxVrtAdd &add)
{
	int btable = 8;	// ノード頂点数4の場合，1番目のビットは1(1xxx)
	//int ntable = (table_index & 0x0F);	// ノード頂点ビット列

	// デプス値がもっとも大きいノード
	int zero_node = 0;
	double max_depth = 0.0;
	for(int k = 1; k < 4; ++k){
		if(g->node_depth[k] > max_depth){
			max_depth = g->node_depth[k];
			zero_node = k;
		}
	}

	// ノードのデプス値の大小でビット列を変更(2-4番目のビット)
	for(int k = 0; k < 3; ++k){
		int k0 = tex1Dfetch(g_texSSNodeTable, zero_node*8+k);
		int k1 = tex1Dfetch(g_texSSNodeTable, zero_node*8+k+1);

		// k0 > k1ならば対応するビットを1にする
		if(g->node_depth[k0] > g->node_depth[k1]){
			btable |= BITS[2-k];
		}
	}

	// 頂点ローテーション
	vrot = zero_node;

	// テーブルインデックスのエッジ頂点部分を更新
	table_index &= 0xF0;
	table_index |= btable;

	// 
	// g_NodeTable[zero_node][5,6]の位置に頂点(back-2 vertex)を追加
	//
	int add_edge1 = tex1Dfetch(g_texSSNodeTable, zero_node*8+5)-4;	// 追加するエッジの位置
	int add_edge2 = tex1Dfetch(g_texSSNodeTable, zero_node*8+6)-4;	// 追加するエッジの位置
	float3 add_vrt1 = vrts[g->edge_vrts[add_edge1]];		// 追加頂点
	float3 add_vrt2 = vrts[g->edge_vrts[add_edge2]];		// 追加頂点

	// g_NodeTable[zero_node][4,7]の位置のback vertexのデプスを設定
	int ref_edge1 = tex1Dfetch(g_texSSNodeTable, zero_node*8+4)-4;	// 追加するエッジの位置
	int ref_edge2 = tex1Dfetch(g_texSSNodeTable, zero_node*8+7)-4;	// 追加するエッジの位置
	add_vrt1.z = vrts[g->back_vrts[ref_edge1]].z;
	add_vrt2.z = vrts[g->back_vrts[ref_edge2]].z;

	add.num = 2;
	add.edge[0] = add_edge1;
	add.vrts[0] = add_vrt1;
	add.edge[1] = add_edge2;
	add.vrts[1] = add_vrt2;

	g->back_vrts[4] = add_edge1;
	g->back_vrts[5] = add_edge2;

	add.layer = 2;
}


/*!
 * グリッドごとにメッシュ生成
 */
__global__
void calGridMesh(rxSSGridG* mgrid, int nx, int ny, float zmax, float3 *vrts, int num_vrts, uint *tri_num, 
				 float3 *back2_vrts, uint *back2_occ, int yoffset)
{
	// グリッドインデックス
	int x = blockIdx.x*blockDim.x+threadIdx.x;
	int y = blockIdx.y*blockDim.y+threadIdx.y;

	if(x >= nx || y >= ny) return;

	rxSSGridG &g = mgrid[x+y*nx];

	int v[14];
	// ノード頂点
	// 3 - 2
	// |   |
	// 0 - 1
		 
	// エッジ頂点
	// - 6 -
	// 7   5
	// - 4 -
		 
	// エッジ頂点(back vertex)
	//  - 10 -
	// 11     9
	//  -  8 -
	// back-2 vertex : 12,13

	int table_index = 0;
	for(int k = 0; k < 4; ++k){
		v[k]   = g.node_vrts[k];
		v[k+4] = g.edge_vrts[k];
		v[k+8] = g.back_vrts[k];

		table_index |= ((v[k] != -1) ? BITS[k] : 0);		// ノード頂点下位4ビット
		table_index |= ((v[k+4] != -1) ? BITS[k]*16 : 0);	// エッジ頂点上位4ビット
	}
	v[12] = -1;
	v[13] = -1;

	int rotation = 0;

	g.table_index0 = table_index;	// デバッグ用

	int fidx = g.num_ev*5+g.num_nv;

	rxVrtAdd add;
	add.num = 0;
	if(fidx == 4){
		CuUpdateTableIndexE0N4(table_index, rotation, v, &g);
	}
	else if(fidx >= 5 &&  fidx <= 9){
		CuUpdateTableIndexE1(table_index, rotation, v, &g);
	}
	else if(fidx == 14){
		CuUpdateTableIndexE2N4(table_index, rotation, v, &g);
	}
	else if(fidx == 17 || fidx == 18){
		// 頂点追加
		CuUpdateTableIndexE3N23(table_index, rotation, v, &g, vrts, zmax, add);
	}
	else if(fidx == 19){
		// 頂点追加
		CuUpdateTableIndexE3N4(table_index, rotation, v, &g, vrts, zmax, add);
	}
	else if(fidx == 22){
		CuUpdateTableIndexE4N2(table_index, rotation, v, &g);
	}
	else if(fidx == 23){
		CuUpdateTableIndexE4N3(table_index, rotation, v, &g);
	}
	else if(fidx == 24){
		// 頂点追加
		CuUpdateTableIndexE4N4(table_index, rotation, v, &g, vrts, zmax, add);
	}

	g.back2 = add.num;
	if(add.num){
		if(add.edge[0]%2){
			// y方向エッジ
			int idx = (x+(add.edge[0]-1)/2)+y*(nx+1)+yoffset;
			back2_occ[idx] = 1;
			back2_vrts[idx] = add.vrts[0];
		}
		else{
			// x方向エッジ
			int idx = x+(y+add.edge[0]/2)*nx;
			back2_occ[idx] = 1;
			back2_vrts[idx] = add.vrts[0];
		}
		g.back2 += (add.layer << 2);

		if(add.num == 2){
			if(add.edge[1]%2){
				// y方向エッジ
				int idx = (x+(add.edge[1]-1)/2)+y*(nx+1)+yoffset;
				back2_occ[idx] = 1;
				back2_vrts[idx] = add.vrts[1];
			}
			else{
				// x方向エッジ
				int idx = x+(y+add.edge[1]/2)*nx;
				back2_occ[idx] = 1;
				back2_vrts[idx] = add.vrts[1];
			}
			g.back2 += (add.layer << 4);
		}
	}

	g.table_index1 = table_index;	// デバッグ用

	int num_tri = tex1Dfetch(g_texSSMeshTable, table_index*19);
	if(num_tri > 0){	// グリッド内のメッシュ数が0より大きかったらメッシュ生成
		g.mesh_num = num_tri;
		g.vrot = rotation;
		tri_num[x+y*nx] = num_tri;

		for(int k = 0; k < 14; ++k){
			g.v[k] = v[k];
		}
	}
}

/*!
 * グリッドごとにメッシュ生成
 */
__global__
void genGridMesh(rxSSGridG* mgrid, int nx, int ny, float zmax, float3 *vrts, int num_vrts, uint* tri_array, uint* tri_num_scan, 
				 float3 *back2_vrts, uint *back2_occ, uint *back2_occ_scan, int yoffset, int voffset)
{
	// グリッドインデックス
	int x = blockIdx.x*blockDim.x+threadIdx.x;
	int y = blockIdx.y*blockDim.y+threadIdx.y;

	if(x >= nx || y >= ny) return;

	rxSSGridG &g = mgrid[x+y*nx];
	if(g.mesh_num == 0) return;

	int v[14];
	for(int k = 0; k < 14; ++k){
		v[k] = g.v[k];
	}

	int back2 = g.back2;
	int add_num = back2 & 0x03;
	if(add_num){
		for(int k = 1; k <= add_num; ++k){
			int layer = (back2 >> (2*k)) & 0x03;
			int edge = g.back_vrts[3+k];
			int odd = edge%2;

			// エッジインデックス
			int idx;
			if(odd){	// y方向エッジ
				idx = (x+(edge-1)/2)+y*(nx+1)+yoffset;
			}
			else{		// x方向エッジ
				idx = x+(y+edge/2)*nx;
			}

			// 追加頂点インデックス
			int bidx = back2_occ_scan[idx]+voffset;

			if(layer == 2){
				// 最背面エッジ頂点として追加
				v[11+k] = bidx;
				g.back_vrts[3+k] = bidx;
			}
			else if(layer == 1){
				// 背面エッジ頂点として追加
				if(v[edge+8] == -1){
					v[edge+8] = bidx;
					g.back_vrts[edge] = bidx;
					g.num_bv++;

					v[11+k] = -1;
					g.back_vrts[3+k] = -1;
				}
				else{
					v[11+k] = v[edge+8];
					g.back_vrts[3+k] = v[edge+8];

					v[edge+8] = bidx;
					g.back_vrts[edge] = bidx;
				}
			}
			else{
				// 前面エッジ頂点として追加
				if(v[edge+8] == -1){
					v[edge+8] = v[edge+4];
					g.back_vrts[edge] = v[edge+4];
					g.num_bv++;

					v[edge+4] = bidx;
					g.edge_vrts[edge] = bidx;

					v[11+k] = -1;
					g.back_vrts[3+k] = -1;
				}
				else{
					v[11+k] = v[edge+8];
					g.back_vrts[3+k] = v[edge+8];

					v[edge+8] = v[edge+4];
					g.back_vrts[edge] = v[edge+4];

					v[edge+4] = bidx;
					g.edge_vrts[edge] = bidx;
				}
			}
		}
	}
	
	// デバッグ用
	for(int k = 0; k < 14; ++k){
		g.v[k] = v[k];
	}

	int m = g.mesh_num;
	if(m > 0){	// グリッド内のメッシュ数が0より大きかったらメッシュ生成
		uint tri[3];
		uint midx = tri_num_scan[x+y*nx];

		for(int k = 0; k < m; ++k){
			for(int l = 0; l < 3; ++l){
				int tidx0 = tex1Dfetch(g_texSSMeshTable, g.table_index1*19+(k*3+l+1));
				int tidx = tex1Dfetch(g_texSSVRotTable, g.vrot*14+tidx0);
				if(v[tidx] == -1) v[tidx] = 0;

				tri[l] = v[tidx];
			}

			g.mesh[k] = midx+k;

			tri_array[3*(midx+k)+0] = tri[0];
			tri_array[3*(midx+k)+1] = tri[1];
			tri_array[3*(midx+k)+2] = tri[2];
		}
	}
}


/*!
 * 輪郭の平滑化 : 平均位置の計算
 */
__global__
void smoothSilhouette(float3 *vrts, int num_nvrts, uint* tri_array, int num_tris, float4 *avg_vrts)
{
	// インデックス
	uint blockId = __mul24(blockIdx.y, gridDim.x)+blockIdx.x;
	uint i = __mul24(blockId, blockDim.x)+threadIdx.x;
	
	if(i >= num_tris) return;

	for(int j = 0; j < 3; ++j){
		int idx = tri_array[3*i+j];
		if(idx >= num_nvrts){
			float4 avg = avg_vrts[idx-num_nvrts];

			// 頂点自身
			if(avg.w == 0){
				avg += make_float4(vrts[idx], 1.0);
			}

			// 隣接頂点
			int jn = j;
			for(int k = 0; k < 2; ++k){
				jn++;
				if(jn == 3) jn = 0;

				int nidx = tri_array[3*i+jn];
				if(nidx >= num_nvrts){
					avg += make_float4(vrts[nidx], 1.0);
				}
				else{
					avg += make_float4(0.5*vrts[nidx], 0.5);
				}
			}

			avg_vrts[idx-num_nvrts] = avg;
		}
	}
}

/*!
 * 輪郭の平滑化 : エッジ頂点を平均位置に移動
 */
__global__
void smoothSilhouette2(float3 *vrts, int num_vrts, int num_nvrts, float4 *avg_vrts)
{
	// インデックス
	uint blockId = __mul24(blockIdx.y, gridDim.x)+blockIdx.x;
	uint i = __mul24(blockId, blockDim.x)+threadIdx.x;
	
	if(i >= num_vrts-num_nvrts) return;

	float4 avg = avg_vrts[i];
	if(avg.w > 0.01){
		avg.x /= avg.w;
		avg.y /= avg.w;
		avg.z /= avg.w;
		vrts[num_nvrts+i] = make_float3(avg.x, avg.y, avg.z);
	}
}


/*!
 * スクリーンスペースから元の3D空間に戻す
 */
__global__
void transfromBack(float3 *ssvrts, float3 *vrts, int num_vrts, matrix4x4 IMVQ, float4 Q, float W, float H)
{
	uint idx = __mul24(blockIdx.x,blockDim.x)+threadIdx.x;
	if(idx >= num_vrts) return;	

	float3 xp = ssvrts[idx];
	float4 xd;
	xd.x = -1.0+2.0*xp.x/W;
	xd.y = -1.0+2.0*xp.y/H;
	xd.w = (1.0-Q.z*xp.z)/(Q.x*xd.x+Q.y*xd.y+Q.w);

	xd.x *= xd.w;
	xd.y *= xd.w;
	xd.z = xp.z;

	// 逆投影変換
	float4 x;
	x.x = dot(IMVQ.e[0], xd);
	x.y = dot(IMVQ.e[1], xd);
	x.z = dot(IMVQ.e[2], xd);
	x.w = dot(IMVQ.e[3], xd);

	vrts[idx] = make_float3(x.x, x.y, x.z);
}


/*!
 * 頂点法線の計算 : 面法線の蓄積
 */
__global__
void sumFaceNormal(float3 *vrts, int num_nvrts, uint* tri_array, int num_tris, float *nrms)
{
	// インデックス
	uint blockId = __mul24(blockIdx.y, gridDim.x)+blockIdx.x;
	uint i = __mul24(blockId, blockDim.x)+threadIdx.x;
	
	if(i >= num_tris) return;

	uint id0, id1, id2;
	id0 = tri_array[3*i+0];
	id1 = tri_array[3*i+1];
	id2 = tri_array[3*i+2];

	float3 vec1, vec2, normal;
	vec1 = vrts[id1]-vrts[id0];
	vec2 = vrts[id2]-vrts[id0];
	normal = cross(vec1, vec2);

	atomicFloatAdd(&nrms[3*id0],   normal.x);
	atomicFloatAdd(&nrms[3*id0+1], normal.y);
	atomicFloatAdd(&nrms[3*id0+2], normal.z);

	atomicFloatAdd(&nrms[3*id1],   normal.x);
	atomicFloatAdd(&nrms[3*id1+1], normal.y);
	atomicFloatAdd(&nrms[3*id1+2], normal.z);

	atomicFloatAdd(&nrms[3*id2],   normal.x);
	atomicFloatAdd(&nrms[3*id2+1], normal.y);
	atomicFloatAdd(&nrms[3*id2+2], normal.z);
}

/*!
 * 頂点法線の計算 : 頂点法線の正規化
 */
__global__
void normalizeNormal(float3 *nrms, int num_nrms)
{
	// インデックス
	uint blockId = __mul24(blockIdx.y, gridDim.x)+blockIdx.x;
	uint i = __mul24(blockId, blockDim.x)+threadIdx.x;
	
	if(i >= num_nrms) return;

	float3 normal = nrms[i];
	nrms[i] = normalize(normal);
}






#endif // #ifndef _RX_SSM_KERNEL_CU_



