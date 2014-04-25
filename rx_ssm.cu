/*! 
  @file rx_ssm.cu
	
  @brief SSM法によるメッシュ生成
		- M. Muller et al., Screen space meshes, SCA2007, 2007. 

  @author Makoto Fujisawa
  @date 2011-05
*/
// FILE --rx_ssm.cu--



//-----------------------------------------------------------------------------
// インクルードファイル
//-----------------------------------------------------------------------------
#include <cstdio>
#include <GL/glew.h>
#include <GL/glut.h>

#include "rx_ssm_kernel.cu"
#include "rx_ssm_tables.h"

#include "rx_cu_funcs.cuh"


//-----------------------------------------------------------------------------
// MARK:グローバル変数
//-----------------------------------------------------------------------------
// MS法のテーブル
uint* g_puSSMeshTable = 0;
uint* g_puSSEdgeTable = 0;
uint* g_puSSNodeTable = 0;
uint* g_puSSVRotTable = 0;


//-----------------------------------------------------------------------------
// CUDA関数
//-----------------------------------------------------------------------------
extern "C"
{


/*!
 * 頂点情報の初期化
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
 * 頂点情報の初期化
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
 * 頂点情報の消去
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
 * 頂点情報の消去
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
 * MS法のテーブル
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
 * MS法のテーブルの破棄
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
 * パラメータをコンスタントメモリに転送
 */
void CuSetSSMParameters(rxSsmParams *hostParams)
{
	RX_CUCHECK( cudaMemcpyToSymbol(g_cSSMParams, hostParams, sizeof(rxSsmParams)) );
}
/*!
 * float型の配列を初期化
 * @param[in] dArray 値を代入したいfloat型の配列
 * @param[in] val 代入する値
 * @param[in] n 配列のサイズ
 */
void CuInitFloatArray(float* dArray, float val, int n)
{
	// 1スレッド/要素
	uint grid, block;	// グリッド内ブロック数，ブロック内スレッド数
	block = THREAD_NUM;
	grid = DivCeil(n, block);

	// カーネル実行
	initFloatArray<<<grid, block>>>(dArray, val, n);

	RX_CUERROR("initDepthMap kernel execution failed");	// カーネル実行エラーチェック
	RX_CUCHECK(cudaThreadSynchronize());		// 全てのスレッドが終わるのを待つ
}

/*!
 * パーティクル半径の配列を初期化
 * @param[in] drad 出力半径配列(サイズ=n*3)
 * @param[in] dval 入力半径配列(サイズ=n)
 * @param[in] n 配列のサイズ
 */
void CuInitRadArray(float* drad, float* dval, int n)
{
	// 1スレッド/要素
	uint grid, block;	// グリッド内ブロック数，ブロック内スレッド数
	block = THREAD_NUM;
	grid = DivCeil(n, block);

	// カーネル実行
	initRadArray<<<grid, block>>>((float3*)drad, dval, n);

	RX_CUERROR("initDepthMap kernel execution failed");	// カーネル実行エラーチェック
	RX_CUCHECK(cudaThreadSynchronize());		// 全てのスレッドが終わるのを待つ
}

/*!
 * デプスマップを∞で初期化
 * @param[in] dDMap デプスマップ((nx+1)x(ny+1))
 * @param[in] nx,ny グリッド解像度
 */
void CuInitDepthMap(float* dDMap, int nx, int ny)
{
	int n = (nx+1)*(ny+1);

	// 1スレッド/ピクセル
	uint grid, block;	// グリッド内ブロック数，ブロック内スレッド数
	block = THREAD_NUM;
	grid = DivCeil(n, block);

	// カーネル実行
	initDepthMap<<<grid, block>>>(dDMap, n);

	RX_CUERROR("initDepthMap kernel execution failed");	// カーネル実行エラーチェック
	RX_CUCHECK(cudaThreadSynchronize());		// 全てのスレッドが終わるのを待つ
}

/*!
 * パーティクルからデプスマップの計算
 * @param[in] dDMap デプスマップ((nx+1)x(ny+1)xfloat)
 * @param[inout] dPrtPos パーティクル座標(pnumxfloat3)，スクリーンスペースでの座標を返す
 * @param[out] dPrtRad スクリーンスペースでのパーティクル半径(pnumxfloat3)
 * @param[in] pnum パーティクル数
 * @param[in] pdim 1パーティクルのメモリサイズ(3 or 4)
 * @param[in] tr パーティクル半径計算用係数(要素数3)
 * @param[in] pmv 透視投影とモデルビュー行列を掛けた行列(要素数16=4x4)
 * @param[in] W,H スクリーンの解像度
 * @param[in] dw,dh グリッド幅
 * @param[in] nx,ny グリッド解像度
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

	// 1スレッド/パーティクル
	uint grid, block;	// グリッド内ブロック数，ブロック内スレッド数
	block = THREAD_NUM;
	grid = DivCeil(pnum, block);

	// カーネル実行
	if(pdim == 3){
		calDepthMap<<<grid, block>>>(dDMap, (float3*)dPrtPos, (float3*)dPrtRad, pnum, tr3, mat_pmv, 
									 (float)W, (float)H, dw, dh, nx, ny);
	}
	else{
		calDepthMap<<<grid, block>>>(dDMap, (float4*)dPrtPos, (float3*)dPrtRad, pnum, tr3, mat_pmv, 
									 (float)W, (float)H, dw, dh, nx, ny);
	}

	RX_CUERROR("calDepthMap kernel execution failed");	// カーネル実行エラーチェック
	RX_CUCHECK(cudaThreadSynchronize());		// 全てのスレッドが終わるのを待つ
}

/*!
 * デプスマップに平滑化を施す
 * @param[in] dDMap デプスマップ((nx)x(ny))
 * @param[in] nx,ny マップ解像度
 * @param[in] n_filter フィルタ幅
 * @param[in] binomials 二項演算係数(n_filter = 0, 1, 2,..,RX_MAX_FILTER_SIZE の係数が順番に格納されている)
 * @param[in] zmax 輪郭となるデプス差の閾値
 */
void CuDepthSmoothing(float* dDMap, int nx, int ny, int n_filter, float *binomials, float zmax)
{
	int b = (n_filter+1)*(n_filter+1);
	cudaMemcpyToSymbol(g_cBinomials, binomials, sizeof(float)*b);

	// 1スレッド/ピクセル(各スレッドがマップのピクセルに対応するように2次元的に配置)
	dim3 block(BLOCK_SIZE, BLOCK_SIZE);
	dim3 grid((nx+block.x-1)/block.x, (ny+block.y-1)/block.y);

	// ブロックごとのシェアードメモリサイズ(最大16KB)
	int sbytes = (block.x+2*n_filter)*(block.y+2*n_filter)*sizeof(float);

	// カーネル実行(x方向平滑化)
	smoothDepthMapX<<< grid, block, sbytes >>>(dDMap, nx, ny, n_filter, block.x+2*n_filter, zmax);

	RX_CUERROR("smoothDepthMapX kernel execution failed");	// カーネル実行エラーチェック
	RX_CUCHECK(cudaThreadSynchronize());		// 全てのスレッドが終わるのを待つ

	// カーネル実行(y方向平滑化)
	smoothDepthMapY<<< grid, block, sbytes >>>(dDMap, nx, ny, n_filter, block.x+2*n_filter, zmax);

	RX_CUERROR("smoothDepthMapY kernel execution failed");	// カーネル実行エラーチェック
	RX_CUCHECK(cudaThreadSynchronize());		// 全てのスレッドが終わるのを待つ
}


/*!
 * 輪郭エッジの検出とfront edge vertexの計算
 * @param[in] dDMap デプスマップ(サイズ = (nx+1)x(ny+1))
 * @param[in] nx,ny マップ解像度
 * @param[in] dw,dh グリッド幅
 * @param[in] zmax 輪郭となるデプス差の閾値
 * @param[out] dNodeVrts ノード頂点(サイズ = (nx+1)x(ny+1))
 * @param[out] dMGrid メッシュ生成用グリッド(サイズ = (nx)x(ny))
 * @param[out] num_node_vertex ノード頂点数
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
	// ノード頂点
	//
	// 1スレッド/ノード
	block2 = dim3(BLOCK_SIZE, BLOCK_SIZE);
	grid2 = dim3(((nx+1)+block2.x-1)/block2.x, ((ny+1)+block2.y-1)/block2.y);
	
	// カーネル実行
	calNodeVertex<<< grid2, block2 >>>(dDMap, nx+1, ny+1, dw, dh, zmax, (float3*)dNodeVrts.dPos, dNodeVrts.dOcc);

	RX_CUERROR("calNodeVertex kernel execution failed");	// カーネル実行エラーチェック
	RX_CUCHECK(cudaThreadSynchronize());		// 全てのスレッドが終わるのを待つ


	//
	// ノード頂点のパッキング
	//
	// ノード頂点があるなら1,そうでないなら0が格納された配列をScan
	CuScan(dNodeVrts.dOccScan, dNodeVrts.dOcc, size);

	// Exclusive scan (最後の要素が0番目からn-2番目までの合計になっている)なので，
	// Scan前配列の最後(n-1番目)の要素と合計することで輪郭エッジ数を計算
	RX_CUCHECK(cudaMemcpy((void*)&last_val, (void*)(dNodeVrts.dOcc+size-1), sizeof(uint), cudaMemcpyDeviceToHost));
	RX_CUCHECK(cudaMemcpy((void *) &last_scan_val, (void*)(dNodeVrts.dOccScan+size-1), sizeof(uint), cudaMemcpyDeviceToHost));
	num_node_vertex = last_val+last_scan_val;

	// 詰めた輪郭エッジ情報を格納する領域の確保
	if(dNodeVrts.dCompactedPos) RX_CUCHECK(cudaFree(dNodeVrts.dCompactedPos));
	RX_CUCHECK(cudaMalloc((void**)&dNodeVrts.dCompactedPos, num_node_vertex*3*sizeof(float)));

	// 輪郭エッジを詰める
	block1 = THREAD_NUM;
	grid1 = DivCeil(size, block1);

	// カーネル実行
	compactNodeVertex<<<grid1, block1>>>((float3*)(dNodeVrts.dCompactedPos), dNodeVrts.dOcc, dNodeVrts.dOccScan, (float3*)(dNodeVrts.dPos), size);
	
	RX_CUERROR("compactNodeVertex kernel execution failed");	// カーネル実行エラーチェック
	RX_CUCHECK(cudaThreadSynchronize());		// 全てのスレッドが終わるのを待つ


	//
	// ノード頂点インデックスをメッシュグリッドに格納
	//
	// 1スレッド/グリッド
	block2 = dim3(BLOCK_SIZE, BLOCK_SIZE);
	grid2 = dim3(((nx)+block2.x-1)/block2.x, ((ny)+block2.y-1)/block2.y);
		
	// カーネル実行
	storeNodeVertex<<< grid2, block2 >>>(dMGrid, nx, ny, (float3*)(dNodeVrts.dCompactedPos), dNodeVrts.dOcc, dNodeVrts.dOccScan);

	RX_CUERROR("storeNodeVertex kernel execution failed");	// カーネル実行エラーチェック
	RX_CUCHECK(cudaThreadSynchronize());		// 全てのスレッドが終わるのを待つ
}



/*!
 * 輪郭エッジの検出とfront edge vertexの計算
 * @param[in] dDMap デプスマップ((nx+1)x(ny+1))
 * @param[in] nx,ny マップ解像度
 * @param[in] dw,dh グリッド幅
 * @param[in] zmax 輪郭となるデプス差の閾値
 * @param[in] dPrtPos パーティクル座標(pnumxfloat3)，スクリーンスペースでの座標を返す
 * @param[in] dPrtRad スクリーンスペースでのパーティクル半径(pnumxfloat3)
 * @param[in] pnum パーティクル数
 * @param[in] pdim 1パーティクルのメモリサイズ(3 or 4)
 * @param[in] W,H スクリーンの解像度
 * @param[out] dEdge エッジ情報(nx*(ny+1)+(nx+1)*ny)
 * @param[out] dCompactedEdge 輪郭エッジ情報(エッジ情報から輪郭エッジのみを詰めた物)
 * @param[out] dEdgeSil 輪郭エッジならば1,そうでなければ0を格納する配列(nx*(ny+1)+(nx+1)*ny)
 * @param[out] dEdgeSilScan dEdgeSilのPrefix sum (nx*(ny+1)+(nx+1)*ny)
 * @param[out] dEdgeVrts エッジ頂点(サイズ = (nx*(ny+1)+(nx+1)*ny))
 * @param[out] dBackVrts 背面エッジ頂点(サイズ = (nx*(ny+1)+(nx+1)*ny))
 * @param[out] dMGrid メッシュ生成用グリッド(サイズ = (nx)x(ny))
 * @param[in] num_nv ノード頂点数
 * @param[out] num_edge 輪郭エッジ数
 * @param[out] num_front_edge エッジ頂点数(前面)
 * @param[out] num_back_edge エッジ頂点数(背面)
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
	// 輪郭エッジ検出
	//
	RX_CUCHECK(cudaMemset((void*)dEdge.dOcc, 0, size*sizeof(uint)));

	// x方向輪郭エッジ検出
	// 1スレッド/エッジ
	block2 = dim3(BLOCK_SIZE, BLOCK_SIZE);
	grid2 = dim3(((nx)+block2.x-1)/block2.x, ((ny+1)+block2.y-1)/block2.y);
	
	// カーネル実行
	detectSilhouetteEdgeX<<< grid2, block2 >>>(dDMap, nx, ny+1, dw, dh, zmax, dEdge.dPos, dEdge.dOcc);

	RX_CUERROR("detectSilhouetteEdgeX kernel execution failed");	// カーネル実行エラーチェック
	RX_CUCHECK(cudaThreadSynchronize());		// 全てのスレッドが終わるのを待つ
	
	// y方向輪郭エッジ検出
	// 1スレッド/エッジ
	block2 = dim3(BLOCK_SIZE, BLOCK_SIZE);
	grid2 = dim3(((nx+1)+block2.x-1)/block2.x, ((ny)+block2.y-1)/block2.y);
		
	// カーネル実行
	detectSilhouetteEdgeY<<< grid2, block2 >>>(dDMap, nx+1, ny, dw, dh, zmax, dEdge.dPos+offset, dEdge.dOcc+offset);

	RX_CUERROR("detectSilhouetteEdgeY kernel execution failed");	// カーネル実行エラーチェック
	RX_CUCHECK(cudaThreadSynchronize());		// 全てのスレッドが終わるのを待つ
	
	//
	// 輪郭エッジの可否でPrefix Sum作成して，エッジ情報を詰める
	//
	// 輪郭エッジなら1,そうでないなら0が格納された配列をScan
	CuScan(dEdge.dOccScan, dEdge.dOcc, size);

	// Exclusive scan (最後の要素が0番目からn-2番目までの合計になっている)なので，
	// Scan前配列の最後(n-1番目)の要素と合計することで輪郭エッジ数を計算
	RX_CUCHECK(cudaMemcpy((void*)&last_val, (void*)(dEdge.dOcc+size-1), sizeof(uint), cudaMemcpyDeviceToHost));
	RX_CUCHECK(cudaMemcpy((void*)&last_scan_val, (void*)(dEdge.dOccScan+size-1), sizeof(uint), cudaMemcpyDeviceToHost));
	num_edge = last_val+last_scan_val;

	// 詰めた輪郭エッジ情報を格納する領域の確保
	if(dEdge.dCompactedPos) RX_CUCHECK(cudaFree(dEdge.dCompactedPos));
	RX_CUCHECK(cudaMalloc((void**)&dEdge.dCompactedPos, num_edge*sizeof(rxSSEdgeG)));

	// 輪郭エッジを詰める
	block1 = THREAD_NUM;
	grid1 = DivCeil(size, block1);

	// カーネル実行
	compactSilhouetteEdges<<<grid1, block1>>>(dEdge.dCompactedPos, dEdge.dOcc, dEdge.dOccScan, dEdge.dPos, size);
	
	RX_CUERROR("CompactSilhouetteEdges kernel execution failed");	// カーネル実行エラーチェック
	RX_CUCHECK(cudaThreadSynchronize());		// 全てのスレッドが終わるのを待つ



	//
	// 輪郭エッジのfront vertexを算出
	//
	// エッジ頂点有無情報配列の初期化
	RX_CUCHECK(cudaMemset((void*)dEdgeVrts.dOcc, 0, size*sizeof(uint)));

	// 1スレッド/パーティクル
	block1 = THREAD_NUM;
	grid1 = DivCeil(pnum, block1);

	// カーネル実行
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

	RX_CUERROR("calFrontEdgeVertex kernel execution failed");	// カーネル実行エラーチェック
	RX_CUCHECK(cudaThreadSynchronize());		// 全てのスレッドが終わるのを待つ
	
	//
	// エッジ頂点有無情報のPrefix Sumを作成して，エッジ頂点を詰める
	//
	// 輪郭エッジなら1,そうでないなら0が格納された配列をScan
	CuScan(dEdgeVrts.dOccScan, dEdgeVrts.dOcc, size);

	// 輪郭エッジ数を計算
	RX_CUCHECK(cudaMemcpy((void*)&last_val, (void*)(dEdgeVrts.dOcc+size-1), sizeof(uint), cudaMemcpyDeviceToHost));
	RX_CUCHECK(cudaMemcpy((void*)&last_scan_val, (void*)(dEdgeVrts.dOccScan+size-1), sizeof(uint), cudaMemcpyDeviceToHost));
	num_front_edge = last_val+last_scan_val;

	// 詰めた輪郭エッジ情報を格納する領域の確保
	if(dEdgeVrts.dCompactedPos) RX_CUCHECK(cudaFree(dEdgeVrts.dCompactedPos));
	RX_CUCHECK(cudaMalloc((void**)&dEdgeVrts.dCompactedPos, num_front_edge*sizeof(float)*3));

	// 輪郭エッジを詰める
	block1 = THREAD_NUM;
	grid1 = DivCeil(size, block1);

	// カーネル実行
	compactEdgeVertices<<<grid1, block1>>>((float3*)dEdgeVrts.dCompactedPos, dEdgeVrts.dOcc, dEdgeVrts.dOccScan, (float3*)dEdgeVrts.dPos, dEdge.dPos, size);
	
	RX_CUERROR("compactEdgeVertices kernel execution failed");	// カーネル実行エラーチェック
	RX_CUCHECK(cudaThreadSynchronize());		// 全てのスレッドが終わるのを待つ



	//
	// 輪郭エッジのback vertexを算出
	//
	// エッジ頂点有無情報配列の初期化
	RX_CUCHECK(cudaMemset((void*)dBackVrts.dOcc, 0, size*sizeof(uint)));

	// x方向エッジ
	// 1スレッド/エッジ
	block2 = dim3(BLOCK_SIZE, BLOCK_SIZE);
	grid2 = dim3(((nx)+block2.x-1)/block2.x, ((ny+1)+block2.y-1)/block2.y);
	
	// カーネル実行
	calBackEdgeVertexX<<<grid2, block2>>>(dDMap, dEdge.dPos, dEdge.dOcc, nx, ny+1, dw, dh, W, H, 
										  (float3*)dEdgeVrts.dPos, (float3*)dBackVrts.dPos, dBackVrts.dOcc);

	RX_CUERROR("calBackEdgeVertexX kernel execution failed");	// カーネル実行エラーチェック
	RX_CUCHECK(cudaThreadSynchronize());		// 全てのスレッドが終わるのを待つ

	offset = (nx)*(ny+1);

	// y方向エッジ
	// 1スレッド/エッジ
	block2 = dim3(BLOCK_SIZE, BLOCK_SIZE);
	grid2 = dim3(((nx+1)+block2.x-1)/block2.x, ((ny)+block2.y-1)/block2.y);
	
	// カーネル実行
	calBackEdgeVertexY<<<grid2, block2>>>(dDMap, dEdge.dPos, dEdge.dOcc, nx+1, ny, dw, dh, W, H, 
										  (float3*)dEdgeVrts.dPos, (float3*)dBackVrts.dPos, dBackVrts.dOcc, offset);

	RX_CUERROR("calBackEdgeVertexY kernel execution failed");	// カーネル実行エラーチェック
	RX_CUCHECK(cudaThreadSynchronize());		// 全てのスレッドが終わるのを待つ

	//
	// エッジ頂点有無情報のPrefix Sumを作成して，エッジ頂点を詰める
	//
	// 輪郭エッジなら1,そうでないなら0が格納された配列をScan
	CuScan(dBackVrts.dOccScan, dBackVrts.dOcc, size);

	// 輪郭エッジ数を計算
	RX_CUCHECK(cudaMemcpy((void*)&last_val, (void*)(dBackVrts.dOcc+size-1), sizeof(uint), cudaMemcpyDeviceToHost));
	RX_CUCHECK(cudaMemcpy((void*)&last_scan_val, (void*)(dBackVrts.dOccScan+size-1), sizeof(uint), cudaMemcpyDeviceToHost));
	num_back_edge = last_val+last_scan_val;

	if(num_back_edge){
		// 詰めた輪郭エッジ情報を格納する領域の確保
		if(dBackVrts.dCompactedPos) RX_CUCHECK(cudaFree(dBackVrts.dCompactedPos));
		RX_CUCHECK(cudaMalloc((void**)&dBackVrts.dCompactedPos, num_back_edge*sizeof(float)*3));

		// 輪郭エッジを詰める
		block1 = THREAD_NUM;
		grid1 = DivCeil(size, block1);

		// カーネル実行
		compactBackEdgeVertices<<<grid1, block1>>>((float3*)dBackVrts.dCompactedPos, dBackVrts.dOcc, dBackVrts.dOccScan, (float3*)dBackVrts.dPos, size);
	
		RX_CUERROR("compactEdgeVertices kernel execution failed");	// カーネル実行エラーチェック
		RX_CUCHECK(cudaThreadSynchronize());		// 全てのスレッドが終わるのを待つ
	}


	//
	// エッジ頂点インデックスをメッシュグリッドに格納
	//
	// 1スレッド/グリッド
	block2 = dim3(BLOCK_SIZE, BLOCK_SIZE);
	grid2 = dim3(((nx)+block2.x-1)/block2.x, ((ny)+block2.y-1)/block2.y);
		
	// カーネル実行
	storeEdgeVertex<<< grid2, block2 >>>(dMGrid, nx, ny, offset, num_nv, num_nv+num_front_edge, 
										 (float3*)(dEdgeVrts.dCompactedPos), dEdgeVrts.dOcc, dEdgeVrts.dOccScan, 
										 (float3*)(dBackVrts.dCompactedPos), dBackVrts.dOcc, dBackVrts.dOccScan);

	RX_CUERROR("storeEdgeVertex kernel execution failed");	// カーネル実行エラーチェック
	RX_CUCHECK(cudaThreadSynchronize());		// 全てのスレッドが終わるのを待つ
}



/*!
 * メッシュ生成
 * @param[in] dDMap デプスマップ(サイズ = (nx+1)x(ny+1))
 * @param[in] nx,ny マップ解像度
 * @param[in] dw,dh グリッド幅
 * @param[in] zmax 輪郭となるデプス差の閾値
 * @param[out] dTriNum 各グリッドのメッシュ数
 * @param[out] dTriNumScan 各グリッドのメッシュ数のScan
 * @param[out] dBack2Vrts 最背面エッジ頂点(サイズ = (nx*(ny+1)+(nx+1)*ny))
 * @param[in] dVrts メッシュ頂点列
 * @param[in] num_vrts メッシュ頂点数
 * @param[out] num_back2_vrts 最背面エッジ頂点数
 * @param[out] dTriArray メッシュ
 * @param[out] num_mesh メッシュ数
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


	// エッジ頂点有無情報配列の初期化
	RX_CUCHECK(cudaMemset((void*)dBack2Vrts.dOcc, 0, size_e*sizeof(uint)));

	// メッシュ数配列の初期化
	RX_CUCHECK(cudaMemset((void*)dTriNum, 0, size*sizeof(uint)));

	//
	// 各グリッドでメッシュ生成
	//
	// 1スレッド/グリッド
	block2 = dim3(BLOCK_SIZE, BLOCK_SIZE);
	grid2 = dim3(((nx)+block2.x-1)/block2.x, ((ny)+block2.y-1)/block2.y);

	// カーネル実行
	calGridMesh<<< grid2, block2 >>>(dMGrid, nx, ny, zmax, (float3*)dVrts, num_vrts, dTriNum, 
									 (float3*)dBack2Vrts.dPos, dBack2Vrts.dOcc, offset);

	RX_CUERROR("calGridMesh kernel execution failed");	// カーネル実行エラーチェック
	RX_CUCHECK(cudaThreadSynchronize());		// 全てのスレッドが終わるのを待つ
	
	//
	// エッジ頂点有無情報のPrefix Sumを作成して，エッジ頂点を詰める
	//
	// エッジ頂点があるなら1,そうでないなら0が格納された配列をScan
	CuScan(dBack2Vrts.dOccScan, dBack2Vrts.dOcc, size_e);

	// 輪郭エッジ数を計算
	RX_CUCHECK(cudaMemcpy((void*)&last_val, (void*)(dBack2Vrts.dOcc+size_e-1), sizeof(uint), cudaMemcpyDeviceToHost));
	RX_CUCHECK(cudaMemcpy((void*)&last_scan_val, (void*)(dBack2Vrts.dOccScan+size_e-1), sizeof(uint), cudaMemcpyDeviceToHost));
	num_back2_vrts = last_val+last_scan_val;

	if(num_back2_vrts){
		// 詰めた輪郭エッジ情報を格納する領域の確保
		if(dBack2Vrts.dCompactedPos) RX_CUCHECK(cudaFree(dBack2Vrts.dCompactedPos));
		RX_CUCHECK(cudaMalloc((void**)&dBack2Vrts.dCompactedPos, num_back2_vrts*sizeof(float)*3));

		// 輪郭エッジを詰める
		block1 = THREAD_NUM;
		grid1 = DivCeil(size_e, block1);

		// カーネル実行
		compactBackEdgeVertices<<<grid1, block1>>>((float3*)dBack2Vrts.dCompactedPos, dBack2Vrts.dOcc, dBack2Vrts.dOccScan, (float3*)dBack2Vrts.dPos, size_e);
	
		RX_CUERROR("compactBackEdgeVertices kernel execution failed");	// カーネル実行エラーチェック
		RX_CUCHECK(cudaThreadSynchronize());		// 全てのスレッドが終わるのを待つ
	}

	// メッシュ数配列をScan
	CuScan(dTriNumScan, dTriNum, size);
	
	// 総メッシュ数を計算
	RX_CUCHECK(cudaMemcpy((void*)&last_val, (void*)(dTriNum+size-1), sizeof(uint), cudaMemcpyDeviceToHost));
	RX_CUCHECK(cudaMemcpy((void*)&last_scan_val, (void*)(dTriNumScan+size-1), sizeof(uint), cudaMemcpyDeviceToHost));
	num_mesh = last_val+last_scan_val;

	// メッシュを格納する領域の確保
	if(dTriArray) RX_CUCHECK(cudaFree(dTriArray));
	RX_CUCHECK(cudaMalloc((void**)&dTriArray, num_mesh*sizeof(uint)*3));

	// 1スレッド/グリッド
	block2 = dim3(BLOCK_SIZE, BLOCK_SIZE);
	grid2 = dim3(((nx)+block2.x-1)/block2.x, ((ny)+block2.y-1)/block2.y);

	//printf("num_vrts = %d\n", num_vrts);

	// カーネル実行
	genGridMesh<<< grid2, block2 >>>(dMGrid, nx, ny, zmax, (float3*)dVrts, num_vrts, dTriArray, dTriNumScan, 
									 (float3*)dBack2Vrts.dPos, dBack2Vrts.dOcc, dBack2Vrts.dOccScan, offset, num_vrts);

	RX_CUERROR("genGridMesh kernel execution failed");	// カーネル実行エラーチェック
	RX_CUCHECK(cudaThreadSynchronize());		// 全てのスレッドが終わるのを待つ
}

/*!
 * 輪郭平滑化
 * @param[inout] dVrts メッシュ頂点列
 * @param[in] num_vrts メッシュ頂点総数
 * @param[in] num_node_vrts ノード頂点数
 * @param[in] dTriArray メッシュ
 * @param[in] num_mesh メッシュ数
 * @param[in] n_iter 平滑化反復数
 */
void CuSilhouetteSmoothing(float* dVrts, int num_vrts, int num_node_vrts, uint* dTriArray, int num_mesh, int n_iter)
{
	uint block1, grid1;

	float4 *dAvgPos;
	RX_CUCHECK(cudaMalloc(&dAvgPos, (num_vrts-num_node_vrts)*sizeof(float4)));


	for(int l = 0; l < n_iter; ++l){
		RX_CUCHECK(cudaMemset(dAvgPos, 0, (num_vrts-num_node_vrts)*sizeof(float4)));

		// 1スレッド/メッシュ
		block1 = THREAD_NUM;
		grid1 = DivCeil(num_mesh, block1);

		// カーネル実行(平均位置の算出)
		smoothSilhouette<<< grid1, block1 >>>((float3*)dVrts, num_node_vrts, dTriArray, num_mesh, dAvgPos);

		RX_CUERROR("smoothSilhouette kernel execution failed");	// カーネル実行エラーチェック
		RX_CUCHECK(cudaThreadSynchronize());		// 全てのスレッドが終わるのを待つ


		// 1スレッド/頂点
		block1 = THREAD_NUM;
		grid1 = DivCeil(num_vrts-num_node_vrts, block1);

		// カーネル実行(平均位置の算出)
		smoothSilhouette2<<< grid1, block1 >>>((float3*)dVrts, num_vrts, num_node_vrts, dAvgPos);

		RX_CUERROR("smoothSilhouette2 kernel execution failed");	// カーネル実行エラーチェック
		RX_CUCHECK(cudaThreadSynchronize());		// 全てのスレッドが終わるのを待つ
	}
}

/*!
 * スクリーンスペースから頂点列を元の3D空間に戻す
 * @param[in] dSSVrts スクリーンスペースメッシュ頂点列
 * @param[in] num_vrts メッシュ頂点総数
 * @param[in] mvq 透視投影とモデルビュー逆行列を掛けた行列(要素数16=4x4)
 * @param[in] W,H スクリーンの解像度
 * @param[out] dVrts 3Dメッシュ頂点列
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

	// 1スレッド/頂点
	block1 = THREAD_NUM;
	grid1 = DivCeil(num_vrts, block1);

	// カーネル実行(平均位置の算出)
	transfromBack<<< grid1, block1 >>>((float3*)dSSVrts, (float3*)dVrts, num_vrts, mat_mvq, vec_q, W, H);

	RX_CUERROR("transfromBack kernel execution failed");	// カーネル実行エラーチェック
	RX_CUCHECK(cudaThreadSynchronize());		// 全てのスレッドが終わるのを待つ

}

/*!
 * 頂点法線計算
 * @param[in] dVrts メッシュ頂点列
 * @param[in] num_vrts メッシュ頂点総数
 * @param[in] dTriArray メッシュ
 * @param[in] num_mesh メッシュ数
 * @param[out] dNrms 頂点法線
 */
void CuCalVertexNormal(float* dVrts, int num_vrts, uint* dTriArray, int num_mesh, float* dNrms)
{
	uint block1, grid1;

	// 法線配列の初期化
	RX_CUCHECK(cudaMemset(dNrms, 0, num_vrts*sizeof(float3)));

	// 1スレッド/メッシュ
	block1 = THREAD_NUM;
	grid1 = DivCeil(num_mesh, block1);

	// カーネル実行(平均位置の算出)
	sumFaceNormal<<< grid1, block1 >>>((float3*)dVrts, num_vrts, dTriArray, num_mesh, dNrms);

	RX_CUERROR("sumFaceNormal kernel execution failed");	// カーネル実行エラーチェック
	RX_CUCHECK(cudaThreadSynchronize());		// 全てのスレッドが終わるのを待つ


	// 1スレッド/頂点
	block1 = THREAD_NUM;
	grid1 = DivCeil(num_vrts, block1);

	// カーネル実行(平均位置の算出)
	normalizeNormal<<< grid1, block1 >>>((float3*)dNrms, num_vrts);

	RX_CUERROR("normalizeNormal kernel execution failed");	// カーネル実行エラーチェック
	RX_CUCHECK(cudaThreadSynchronize());		// 全てのスレッドが終わるのを待つ
}


}   // extern "C"
