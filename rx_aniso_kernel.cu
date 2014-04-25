/*! 
  @file rx_aniso_kernel.cu
	
  @brief Anisotropic Kernelの計算
  @ref J. Yu and G. Turk, Reconstructing Surfaces of Particle-Based Fluids Using Anisotropic Kernels, SCA2010. 
 
  @author Makoto Fujisawa
  @date 2009-08, 2011-06
*/

#ifndef _RX_ANISOTROPIC_KERNEL_CU_
#define _RX_ANISOTROPIC_KERNEL_CU_


//-----------------------------------------------------------------------------
// インクルードファイル
//-----------------------------------------------------------------------------
#include "rx_cu_common.cu"


//-----------------------------------------------------------------------------
// MARK:Anisotropic Kernel
//-----------------------------------------------------------------------------
/*!
 * 与えられたセル内のパーティクルとの距離から重み付き平均を計算
 * @param[in] gridPos グリッド位置
 * @param[in] i パーティクルインデックス
 * @param[in] pos0 計算座標
 * @param[in] cell パーティクルグリッドデータ
 * @param[in] h 探索半径
 * @return セル内のパーティクルから計算した重み付き平均値
 */
__device__
float4 calWeighedAvgPositionCell(int3 gridPos, uint i, float3 pos0, float h, rxParticleCell cell)
{
	uint gridHash = calcGridHash(gridPos);

	// セル内のパーティクルのスタートインデックス
	uint startIndex = FETCHC(dCellStart, gridHash);

	float4 posw = make_float4(0.0f);	
	if(startIndex != 0xffffffff){	// セルが空でないかのチェック
		// セル内のパーティクルで反復
		uint endIndex = FETCHC(dCellEnd, gridHash);
		for(uint j = startIndex; j < endIndex; ++j){
			//if(j == i) continue;

			float3 pos1 = make_float3(FETCHC(dSortedPos, j));

			float3 rij = pos0-pos1;
			float r = length(rij);

			if(r < h){
				float q = (1.0f-r/h);
				float wij = q*q*q;
				posw += make_float4(pos1*wij, wij);
			}
		}
	}

	return posw;
}

/*!
 * カーネル中心位置の平滑化と重み付き平均の計算(カーネル関数)
 * @param[out] upPos 平滑化カーネル中心
 * @param[out] wPos 重み付き平均パーティクル座標 
 * @param[in]  lambda 平滑化のための定数
 * @param[in]  h 探索半径
 * @param[in]  cell パーティクルグリッドデータ
 */
__global__
void sphCalUpdatedPosition(float4* upPos, float4* wPos, float lambda, float h, rxParticleCell cell)
{
	// パーティクルインデックス
	uint index = __mul24(blockIdx.x,blockDim.x)+threadIdx.x;
	if(index >= cell.uNumParticles) return;	
	
	float3 pos = make_float3(FETCHC(dSortedPos, index));	// パーティクル位置

	// パーティクル周囲のグリッド
	int3 grid_pos0, grid_pos1;
	grid_pos0 = calcGridPos(pos-make_float3(h));
	grid_pos1 = calcGridPos(pos+make_float3(h));

	// 周囲のグリッドも含めて近傍探索，重み付き平均位置を計算
	float4 posw  = make_float4(0.0f);
	for(int z = grid_pos0.z; z <= grid_pos1.z; ++z){
		for(int y = grid_pos0.y; y <= grid_pos1.y; ++y){
			for(int x = grid_pos0.x; x <= grid_pos1.x; ++x){
				int3 n_grid_pos = make_int3(x, y, z);
				posw += calWeighedAvgPositionCell(n_grid_pos, index, pos, h, cell);
			}
		}
	}

	// 結果を書き込み
	uint oIdx = cell.dSortedIndex[index];
	posw.x /= posw.w;
	posw.y /= posw.w;
	posw.z /= posw.w;
	posw.w = 0.0f;

	wPos[oIdx] = posw;
	upPos[oIdx]  = make_float4((1.0-lambda)*pos+lambda*make_float3(posw));
}


/*!
 * 与えられたセル内のパーティクルとの距離からCovariance Matrixを計算
 * @param[in] gridPos グリッド位置
 * @param[in] i パーティクルインデックス
 * @param[in] pos0 計算座標
 * @param[in] cell パーティクルグリッドデータ
 * @param[in] h 探索半径
 * @return セル内のパーティクルから計算した重み
 */
__device__
float2 calCovarianceMatrixCell(int3 gridPos, uint i, float3 pos0, float h, matrix3x3 &c, float4 xiw, rxParticleCell cell)
{
	uint gridHash = calcGridHash(gridPos);

	// セル内のパーティクルのスタートインデックス
	uint startIndex = FETCHC(dCellStart, gridHash);

	float2 wn = make_float2(0.0f);
	if(startIndex != 0xffffffff){	// セルが空でないかのチェック
		// セル内のパーティクルで反復
		uint endIndex = FETCHC(dCellEnd, gridHash);
		for(uint j = startIndex; j < endIndex; ++j){
			//if(j == i) continue;

			float3 pos1 = make_float3(FETCHC(dSortedPos, j));

			float3 rij = pos0-pos1;
			float r = length(rij);

			if(r < h){
				float q = (1.0f-r/h);
				float wij = q*q*q;

				float3 dxj = pos1-make_float3(xiw);

				c.e[0].x += wij*dxj.x*dxj.x;
				c.e[0].y += wij*dxj.x*dxj.y;
				c.e[0].z += wij*dxj.x*dxj.z;
				c.e[1].x += wij*dxj.y*dxj.x;
				c.e[1].y += wij*dxj.y*dxj.y;
				c.e[1].z += wij*dxj.y*dxj.z;
				c.e[2].x += wij*dxj.z*dxj.x;
				c.e[2].y += wij*dxj.z*dxj.y;
				c.e[2].z += wij*dxj.z*dxj.z;
				wn.x += wij;
				wn.y += 1.0;
			}
		}
	}

	return wn;
}


/*!
 * Covariance Matrixの計算
 * @param[out] PosW 重み付き平均パーティクル座標 
 * @param[out] CMat Covariance Matrix
 * @param[in]  lambda 平滑化のための定数
 * @param[in]  h 探索半径
 * @param[in]  cell パーティクルグリッドデータ
 */
__global__
void sphCalCovarianceMatrix(float4* PosW, matrix3x3 *CMat, float h, rxParticleCell cell)
{
	// パーティクルインデックス
	uint index = __mul24(blockIdx.x,blockDim.x)+threadIdx.x;
	if(index >= cell.uNumParticles) return;	
	
	float3 pos = make_float3(FETCHC(dSortedPos, index));	// パーティクル位置

	// パーティクル周囲のグリッド
	int3 grid_pos0, grid_pos1;
	grid_pos0 = calcGridPos(pos-make_float3(h));
	grid_pos1 = calcGridPos(pos+make_float3(h));

	// 周囲のグリッドも含めて近傍探索，重み付き平均位置を計算
	float4 xiw  = make_float4(0.0f);
	for(int z = grid_pos0.z; z <= grid_pos1.z; ++z){
		for(int y = grid_pos0.y; y <= grid_pos1.y; ++y){
			for(int x = grid_pos0.x; x <= grid_pos1.x; ++x){
				int3 n_grid_pos = make_int3(x, y, z);
				xiw += calWeighedAvgPositionCell(n_grid_pos, index, pos, h, cell);
			}
		}
	}

	xiw.x /= xiw.w;
	xiw.y /= xiw.w;
	xiw.z /= xiw.w;

	//__syncthreads();

	// 行列の要素の初期化
	matrix3x3 c;
	for(int k = 0; k < 3; ++k){
		c.e[k].x = 0.0;
		c.e[k].y = 0.0;
		c.e[k].z = 0.0;
	}

	// 周囲のグリッドも含めて近傍探索，Covariance Matrixを計算
	float2 wn = make_float2(0.0f);
	for(int z = grid_pos0.z; z <= grid_pos1.z; ++z){
		for(int y = grid_pos0.y; y <= grid_pos1.y; ++y){
			for(int x = grid_pos0.x; x <= grid_pos1.x; ++x){
				int3 n_grid_pos = make_int3(x, y, z);
				wn += calCovarianceMatrixCell(n_grid_pos, index, pos, h, c, xiw, cell);
			}
		}
	}

	for(int k = 0; k < 3; ++k){
		c.e[k].x /= wn.x;
		c.e[k].y /= wn.x;
		c.e[k].z /= wn.x;
	}

	xiw.w = wn.y;

	// 結果を書き込み
	uint oIdx = cell.dSortedIndex[index];
	PosW[oIdx] = xiw;
	CMat[oIdx] = c;
}

__device__
inline float RxPythag(const float a, const float b)
{
	float absa = abs(a), absb = abs(b);
	return (absa > absb ? absa*(float)sqrt((double)(1.0+(absb/absa)*(absb/absa))) :
		   (absb == 0.0 ? 0.0 : absb*(float)sqrt((double)(1.0+(absa/absb)*(absa/absb)))));
}

//! 最小値判定(2値)
__device__
inline float RXD_MIN(const float &a, const float &b){ return ((a < b) ? a : b); }

//! 最大値判定(2値)
__device__
inline float RXD_MAX(const float &a, const float &b){ return ((a > b) ? a : b); }

//! aの符号をbの符号にあわせる
__device__
inline float RXD_SIGN2(const float &a, const float &b){ return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a); }

__device__
int svdecomp3(float w[3], float u[9], float v[9], float eps)
{
	bool flag;
	int i, its, j, jj, k, l, nm;
	float anorm, c, f, g, h, s, scale, x, y, z;
	float rv1[3];
	g = scale = anorm = 0.0;
	for(i = 0; i < 3; ++i){
		l = i+2;
		rv1[i] = scale*g;
		g = s = scale = 0.0;
		for(k = i; k < 3; ++k) scale += abs(u[k*3+i]);
		if(scale != 0.0){
			for(k = i; k < 3; ++k){
				u[k*3+i] /= scale;
				s += u[k*3+i]*u[k*3+i];
			}
			f = u[i*3+i];
			g = -RXD_SIGN2(sqrt(s), f);
			h = f*g-s;
			u[i*3+i] = f-g;
			for(j = l-1; j < 3; ++j){
				for(s = 0.0, k = i; k < 3; ++k) s += u[k*3+i]*u[k*3+j];
				f = s/h;
				for(k = i; k < 3; ++k) u[k*3+j] += f*u[k*3+i];
			}
			for(k = i; k < 3; ++k) u[k*3+i] *= scale;
		}

		w[i] = scale*g;
		g = s = scale = 0.0;
		if(i+1 <= 3 && i+1 != 3){
			for(k = l-1; k < 3; ++k) scale += abs(u[i*3+k]);
			if(scale != 0.0){
				for(k = l-1; k < 3; ++k){
					u[i*3+k] /= scale;
					s += u[i*3+k]*u[i*3+k];
				}
				f = u[i*3+l-1];
				g = -RXD_SIGN2(sqrt(s), f);
				h = f*g-s;
				u[i*3+l-1] = f-g;
				for(k = l-1; k < 3; ++k) rv1[k] = u[i*3+k]/h;
				for(j = l-1; j < 3; ++j){
					for(s = 0.0,k = l-1; k < 3; ++k) s += u[j*3+k]*u[i*3+k];
					for(k = l-1; k < 3; ++k) u[j*3+k] += s*rv1[k];
				}
				for(k = l-1; k < 3; ++k) u[i*3+k] *= scale;
			}
		}
		anorm = RXD_MAX(anorm, (abs(w[i])+abs(rv1[i])));
	}
	for(i = 2; i >= 0; --i){
		if(i < 2){
			if(g != 0.0){
				for(j = l; j < 3; ++j){
					v[j*3+i] = (u[i*3+j]/u[i*3+l])/g;
				}
				for(j = l; j < 3; ++j){
					for(s = 0.0, k = l; k < 3; ++k) s += u[i*3+k]*v[k*3+j];
					for(k = l; k < 3; ++k) v[k*3+j] += s*v[k*3+i];
				}
			}
			for(j = l; j < 3; ++j) v[i*3+j] = v[j*3+i] = 0.0;
		}
		v[i*3+i] = 1.0;
		g = rv1[i];
		l = i;
	}
	for(i = 2; i >= 0; --i){
		l = i+1;
		g = w[i];
		for(j = l; j < 3; ++j) u[i*3+j] = 0.0;
		if(g != 0.0){
			g = 1.0/g;
			for(j = l; j < 3; ++j){
				for(s = 0.0, k = l; k < 3; ++k) s += u[k*3+i]*u[k*3+j];
				f = (s/u[i*3+i])*g;
				for(k = i; k < 3; ++k) u[k*3+j] += f*u[k*3+i];
			}
			for(j = i; j < 3; ++j) u[j*3+i] *= g;
		}
		else{
			for(j = i; j < 3; ++j) u[j*3+i] = 0.0;
		}
		++u[i*3+i];
	}
	for(k = 2; k >= 0; --k){
		for(its = 0; its < 30; ++its){
			flag = true;
			for(l = k; l >= 0; --l){
				nm = l-1;
				if(l == 0 || abs(rv1[l]) <= eps*anorm){
					flag = false;
					break;
				}
				if(abs(w[nm]) <= eps*anorm) break;
			}
			if(flag){
				c = 0.0;
				s = 1.0;
				for(i = l; i < k+1; ++i){
					f = s*rv1[i];
					rv1[i] = c*rv1[i];
					if(abs(f) <= eps*anorm) break;
					g = w[i];
					h = RxPythag(f, g);
					w[i] = h;
					h = 1.0/h;
					c = g*h;
					s = -f*h;
					for(j = 0; j < 3; ++j){
						y = u[j*3+nm];
						z = u[j*3+i];
						u[j*3+nm] = y*c+z*s;
						u[j*3+i] = z*c-y*s;
					}
				}
			}
			z = w[k];
			if(l == k){
				if(z < 0.0){
					w[k] = -z;
					for(j = 0; j < 3; ++j) v[j*3+k] = -v[j*3+k];
				}
				break;
			}
			if(its == 29){
				//printf("no convergence in 30 svdcmp iterations");
				return 0;
			}
			x = w[l];
			nm = k-1;
			y = w[nm];
			g = rv1[nm];
			h = rv1[k];
			f = ((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
			g = RxPythag(f, 1.0f);
			f = ((x-z)*(x+z)+h*((y/(f+RXD_SIGN2(g, f)))-h))/x;
			c = s = 1.0;
			for(j = l; j <= nm; ++j){
				i = j+1;
				g = rv1[i];
				y = w[i];
				h = s*g;
				g = c*g;
				z = RxPythag(f, h);
				rv1[j] = z;
				c = f/z;
				s = h/z;
				f = x*c+g*s;
				g = g*c-x*s;
				h = y*s;
				y *= c;
				for(jj = 0; jj < 3; ++jj){
					x = v[jj*3+j];
					z = v[jj*3+i];
					v[jj*3+j] = x*c+z*s;
					v[jj*3+i] = z*c-x*s;
				}
				z = RxPythag(f, h);
				w[j] = z;
				if(z){
					z = 1.0/z;
					c = f*z;
					s = h*z;
				}
				f = c*g+s*y;
				x = c*y-s*g;
				for(jj = 0; jj < 3; ++jj){
					y = u[jj*3+j];
					z = u[jj*3+i];
					u[jj*3+j] = y*c+z*s;
					u[jj*3+i] = z*c-y*s;
				}
			}
			rv1[l] = 0.0;
			rv1[k] = f;
			w[k] = x;
		}
	}

	// reorder
	int inc = 1;
	float sw;
	float su[3], sv[3];

	do{
		inc *= 3;
		inc++; 
	}while(inc <= 3);

	do{
		inc /= 3;
		for(i = inc; i < 3; ++i){
			sw = w[i];
			for(k = 0; k < 3; ++k) su[k] = u[k*3+i];
			for(k = 0; k < 3; ++k) sv[k] = v[k*3+i];
			j = i;
			while (w[j-inc] < sw){
				w[j] = w[j-inc];
				for(k = 0; k < 3; ++k) u[k*3+j] = u[k*3+j-inc];
				for(k = 0; k < 3; ++k) v[k*3+j] = v[k*3+j-inc];
				j -= inc;
				if (j < inc) break;
			}
			w[j] = sw;
			for(k = 0; k < 3; ++k) u[k*3+j] = su[k];
			for(k = 0; k < 3; ++k) v[k*3+j] = sv[k];

		}
	}while(inc > 1);

	for(k = 0; k < 3; ++k){
		s = 0;
		for(i = 0; i < 3; ++i) if(u[i*3+k] < 0.) s++;
		for(j = 0; j < 3; ++j) if(v[j*3+k] < 0.) s++;
		if(s > 3){
			for(i = 0; i < 3; ++i) u[i*3+k] = -u[i*3+k];
			for(j = 0; j < 3; ++j) v[j*3+k] = -v[j*3+k];
		}
	}

	return 1;
}

/*!
 * 特異値分解により固有値を計算
 * @param[in]  CMat Covariance Matrix
 * @param[in]  posw 重み付き平均位置
 * @param[out] eigen 固有値
 * @param[out] RMat 固有ベクトル(回転行列)
 * @param[in]  numParticles パーティクル数
 */
__global__
void sphSVDecomposition(matrix3x3 *CMat, float4* posw, float3* eigen, matrix3x3 *RMat, uint numParticles)
{
	// パーティクルインデックス
	uint index = __mul24(blockIdx.x,blockDim.x)+threadIdx.x;
	if(index >= numParticles) return;

	int n = (int)posw[index].w;
	matrix3x3 tmp = CMat[index];
	float u[9], v[9], w[3];
	for(int i = 0; i < 3; ++i){
		u[i*3+0] = tmp.e[i].x;
		u[i*3+1] = tmp.e[i].y;
		u[i*3+2] = tmp.e[i].z;
	}

	// 特異値分解
	svdecomp3(w, u, v, 1.0e-10);
	
	float3 sigma;
	sigma.x = w[0];
	sigma.y = w[1];
	sigma.z = w[2];
	for(int i = 0; i < 3; ++i){
		tmp.e[i].x = u[i*3+0];
		tmp.e[i].y = u[i*3+1];
		tmp.e[i].z = u[i*3+2];
	}
	
	int ne = 10;
	float ks = 1400;
	float kn = 0.5;
	float kr = 4.0;
	if(n > ne){
		float s0 = sigma.x/kr;
		sigma.y = (sigma.y >= s0 ? sigma.y : s0);
		sigma.z = (sigma.z >= s0 ? sigma.z : s0);
		sigma *= ks;
	}
	else{
		sigma = make_float3(kn*1.0f);
	}

	eigen[index] = sigma;
	RMat[index]  = tmp;
}

__device__
inline float calDeterminant(matrix3x3 &mat)
{
	float d = mat.e[0].x*mat.e[1].y*mat.e[2].z- 
			  mat.e[0].x*mat.e[2].y*mat.e[1].z+ 
			  mat.e[1].x*mat.e[2].y*mat.e[0].z- 
			  mat.e[1].x*mat.e[0].y*mat.e[2].z+ 
			  mat.e[2].x*mat.e[0].y*mat.e[1].z- 
			  mat.e[2].x*mat.e[1].y*mat.e[0].z;
	return d;
}

__device__
inline matrix3x3 calInverse(matrix3x3 &mat)
{
	matrix3x3 inv_mat;

	float d = mat.e[0].x*mat.e[1].y*mat.e[2].z- 
			  mat.e[0].x*mat.e[2].y*mat.e[1].z+ 
			  mat.e[1].x*mat.e[2].y*mat.e[0].z- 
			  mat.e[1].x*mat.e[0].y*mat.e[2].z+ 
			  mat.e[2].x*mat.e[0].y*mat.e[1].z- 
			  mat.e[2].x*mat.e[1].y*mat.e[0].z;

	if(d == 0) d = 1;

	inv_mat.e[0].x =  (mat.e[1].y*mat.e[2].z-mat.e[1].z*mat.e[2].y)/d;
	inv_mat.e[1].x = -(mat.e[0].y*mat.e[2].z-mat.e[0].z*mat.e[2].y)/d;
	inv_mat.e[2].x =  (mat.e[0].y*mat.e[1].z-mat.e[0].z*mat.e[1].y)/d;
	inv_mat.e[0].y = -(mat.e[1].x*mat.e[2].z-mat.e[1].z*mat.e[2].x)/d;
	inv_mat.e[1].y =  (mat.e[0].x*mat.e[2].z-mat.e[0].z*mat.e[2].x)/d;
	inv_mat.e[2].y = -(mat.e[0].x*mat.e[1].z-mat.e[0].z*mat.e[1].x)/d;
	inv_mat.e[0].z =  (mat.e[1].x*mat.e[2].y-mat.e[1].y*mat.e[2].x)/d;
	inv_mat.e[1].z = -(mat.e[0].x*mat.e[2].y-mat.e[0].y*mat.e[2].x)/d;
	inv_mat.e[2].z =  (mat.e[0].x*mat.e[1].y-mat.e[0].y*mat.e[1].x)/d;

	return inv_mat;
}

/*!
 * 固有値，固有ベクトル(回転行列)から変形行列を計算
 * @param[in]  eigen 固有値
 * @param[in]  RMat 固有ベクトル(回転行列)
 * @param[out] GMat 変形行列
 * @param[in]  numParticles パーティクル数
 */
__global__
void sphCalTransformMatrix(float3* eigen, matrix3x3 *RMat, matrix3x3 *GMat, uint numParticles)
{
	// パーティクルインデックス
	uint index = __mul24(blockIdx.x,blockDim.x)+threadIdx.x;
	if(index >= numParticles) return;

	float3 sigma = eigen[index];
	matrix3x3 R = RMat[index];
	matrix3x3 G;

	for(int j = 0; j < 3; ++j){
		G.e[j].x = R.e[j].x*R.e[0].x/sigma.x+R.e[j].y*R.e[0].y/sigma.y+R.e[j].z*R.e[0].z/sigma.z;
		G.e[j].y = R.e[j].x*R.e[1].x/sigma.x+R.e[j].y*R.e[1].y/sigma.y+R.e[j].z*R.e[1].z/sigma.z;
		G.e[j].z = R.e[j].x*R.e[2].x/sigma.x+R.e[j].y*R.e[2].y/sigma.y+R.e[j].z*R.e[2].z/sigma.z;
	}

	float max_diag = -100000000.0f;
	if(G.e[0].x > max_diag) max_diag = G.e[0].x;
	if(G.e[1].y > max_diag) max_diag = G.e[1].y;
	if(G.e[2].z > max_diag) max_diag = G.e[2].z;

	for(int j = 0; j < 3; ++j){
		G.e[j].x /= max_diag;
		G.e[j].y /= max_diag;
		G.e[j].z /= max_diag;
	}

	//G = calInverse(G);

	GMat[index] = G;
}



/*!
 * 与えられたセル内のパーティクルとの距離から密度を計算
 * @param[in] gridPos グリッド位置
 * @param[in] pos 計算座標
 * @param[in] cell パーティクルグリッドデータ
 * @return セル内のパーティクルから計算した密度値
 */
__device__
float calDensityAnisoCellG(int3 gridPos, float3 pos0, float h, matrix3x3 *GMat, rxParticleCell cell)
{
	uint gridHash = calcGridHash(gridPos);

	// セル内のパーティクルのスタートインデックス
	uint startIndex = FETCHC(dCellStart, gridHash);

	float d = 0.0f;
	if(startIndex != 0xffffffff){	// セルが空でないかのチェック
		// セル内のパーティクルで反復
		uint endIndex = FETCHC(dCellEnd, gridHash);

		for(uint j = startIndex; j < endIndex; ++j){
			//if(j != index){
				float3 pos1 = make_float3(FETCHC(dSortedPos, j));
				uint jdx = cell.dSortedIndex[j];

				float3 rij = pos0-pos1;
				matrix3x3 Gj = GMat[jdx];
				float3 rg;
				rg.x = Gj.e[0].x*rij.x+Gj.e[0].y*rij.y+Gj.e[0].z*rij.z;
				rg.y = Gj.e[1].x*rij.x+Gj.e[1].y*rij.y+Gj.e[1].z*rij.z;
				rg.z = Gj.e[2].x*rij.x+Gj.e[2].y*rij.y+Gj.e[2].z*rij.z;

				float r = length(rg);

				if(r <= h){
					float q = h*h-r*r;
					float detG = calDeterminant(Gj);

					d += detG*params.Mass*params.Wpoly6*q*q*q;
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
void sphCalDensityAnisoInGrid(float* GridD, matrix3x3 *GMat, float Emax, 
							  rxParticleCell cell, uint3 gnum, float3 gmin, float3 glen)
{
	uint blockId = __mul24(blockIdx.y, gridDim.x)+blockIdx.x;
	uint i = __mul24(blockId, blockDim.x)+threadIdx.x;

	uint3 gridPos = calcGridPosU(i, gnum);

	if(gridPos.x < gnum.x && gridPos.y < gnum.y && gridPos.z < gnum.z){
		float3 gpos;
		gpos.x = gmin.x+(gridPos.x)*glen.x;
		gpos.y = gmin.y+(gridPos.y)*glen.y;
		gpos.z = gmin.z+(gridPos.z)*glen.z;

		matrix3x3 G = GMat[i];
		float detG = calDeterminant(G);
		int3 pgpos = calcGridPos(gpos);

		float h = params.EffectiveRadius;
		int3 grid_pos0, grid_pos1;
		grid_pos0 = calcGridPos(gpos-make_float3(h*3));
		grid_pos1 = calcGridPos(gpos+make_float3(h*3));

		float P = 0.0f;
		for(int z = grid_pos0.z; z <= grid_pos1.z; ++z){
			for(int y = grid_pos0.y; y <= grid_pos1.y; ++y){
				for(int x = grid_pos0.x; x <= grid_pos1.x; ++x){
					int3 neighbourPos = make_int3(x, y, z);

					P += calDensityAnisoCellG(neighbourPos, gpos, h, GMat, cell);
				}
			}
		}

		GridD[gridPos.x+gridPos.y*gnum.x+gridPos.z*gnum.x*gnum.y] = P;
	}

}


#endif // #ifndef _RX_ANISOTROPIC_KERNEL_CU_



