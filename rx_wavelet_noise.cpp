/*!
 @file rx_wavelet_noise.cpp

 @brief Wavelet noise 生成

		Robert L. Cook and Tony DeRose
		Wavelet noise, 
		SIGGRAPH 2005, 2005. 
  
 @author Makoto Fujisawa
 @date 2010-01
*/

#pragma warning (disable: 4305)

//-----------------------------------------------------------------------------
// インクルードファイル
//-----------------------------------------------------------------------------
#include "rx_wavelet_noise.h"

#include "mt19937ar.h"


//-----------------------------------------------------------------------------
// 定義
//-----------------------------------------------------------------------------
const int ARAD = 16;


//-----------------------------------------------------------------------------
// MARK:グローバル変数
//-----------------------------------------------------------------------------
RXREAL *g_pNoiseTileData = NULL;	//!< ノイズタイルデータ
int g_iNoiseTileSize = 0;


//-----------------------------------------------------------------------------
// 関数
//-----------------------------------------------------------------------------
static RXREAL GaussianNoise()
{
	RXREAL x1, x2;
	RXREAL ret;
	
	RXREAL r2;

	do {

		x1 = 2.0 * genrand_real2() - 1.0;	/* [-1, 1) */
		x2 = 2.0 * genrand_real2() - 1.0;

		r2 = x1*x1 + x2*x2;

	} while ((r2 == 0) || (r2 > 1.0));

	ret = x1 * sqrt((-2.0 * log(r2))/r2);
	ret *= 0.25;		// Possibility of ( N(0, 1) < 4.0 ) = 100%

	if (ret < -1.0) ret = -1.0; /* Account for loss of precision. */
	if (ret >  1.0) ret = 1.0;

	return ret;
}

inline int ModW(int x, int n)
{
	int m = x%n; 
	return ((m < 0) ? m+n : m);
}


void Downsample(const RXREAL *from, RXREAL *to, int n, int stride)
{
	RXREAL *a;
	RXREAL aCoeffs[2*ARAD+1] = { 0.000334,-0.001528, 0.000410, 0.003545,-0.000938,-0.008233, 0.002172, 0.019120, \
							  -0.005040,-0.044412, 0.011655, 0.103311,-0.025936,-0.243780, 0.033979, 0.655340, \
							   0.655340, 0.033979,-0.243780,-0.025936, 0.103311, 0.011655,-0.044412,-0.005040, \
							   0.019120, 0.002172,-0.008233,-0.000938, 0.003546, 0.000410,-0.001528, 0.000334, 0.0};

	a = &aCoeffs[ARAD];
	for(int i = 0; i < n/2; ++i){
		to[i*stride] = 0;
		for(int k = 2*i-ARAD; k <= 2*i+ARAD; ++k){
			to[i*stride] += a[k-2*i]*from[ModW(k, n)*stride];
		}
	}
}

void Upsample(const RXREAL *from, RXREAL *to, int n, int stride)
{
	RXREAL *p, pCoeffs[4] = { 0.25, 0.75, 0.75, 0.25 };
	p = &pCoeffs[2];
	for(int i = 0; i < n; ++i){
		to[i*stride] = 0;
		for(int k = i/2; k <= i/2+1; ++k){
			to[i*stride] += p[i-2*k]*from[ModW(k, n/2)*stride];
		}
	}
}


/*!
 * 2D wavelet noise
 * @param[in] nx,ny 解像度
 */
void rxWaveletNoise::GenerateNoiseTile2(int n)
{
	// MARK:GenerateNoiseTile2
	if(g_pNoiseTileData != NULL && n == g_iNoiseTileSize) return;

	if(n%2) n++; // tile size must be even

	int sz = n*n;

	RXREAL *temp1 = new RXREAL[sz];
	RXREAL *temp2 = new RXREAL[sz];
	RXREAL *noise = new RXREAL[sz];

	// Step 1. Fill the tile with random numbers in the range -1 to 1.
	for(int i = 0; i < n*n; ++i){
		noise[i] = GaussianNoise();
		temp1[i] = temp2[i] = 0;
	}

	// Steps 2 and 3. Downsample and Upsample the tile
	for(int iy = 0; iy < n; ++iy){ // each x row
		int idx = iy*n;
		Downsample(&noise[idx], &temp1[idx], n, 1);
		Upsample(  &temp1[idx], &temp2[idx], n, 1);
	}

	for(int ix = 0; ix < n; ++ix){ // each y row
		int idx = ix; 
		Downsample(&temp2[idx], &temp1[idx], n, n);
		Upsample(  &temp1[idx], &temp2[idx], n, n);
	}

	// Step 4. Subtract out the coarse-scale contribution
	for(int i = 0; i < n*n; ++i){
		noise[i] -= temp2[i];
	}

	// Avoid even/odd variance difference by adding odd-offset version of noise to itself.
	int offset = n/2;
	if(offset%2 == 0) offset++;

	for(int i = 0, ix = 0; ix < n; ++ix){
		for(int iy = 0; iy < n; ++iy){
			temp1[i++] = noise[ModW(ix+offset, n)+ModW(iy+offset, n)*n];
		}
	}

	for(int i = 0; i < n*n; ++i){
		noise[i] += temp1[i];
	}

	g_pNoiseTileData = noise;
	g_iNoiseTileSize = n;

	delete [] temp1;
	delete [] temp2;
}

/*!
 * Non-projected 2D noise
 * @param[in] p[2] 2D座標
 * @return Waveletノイズ値
 */
RXREAL rxWaveletNoise::Noise2(RXREAL p[2])
{
	int f[3], c[3];		// filter, noise coef. indices
	int mid[3], n = g_iNoiseTileSize;

	RXREAL w[2][3], t, result = 0;

 	// 2次のB-スプライン(quadratic B-spline)基底関数を計算
	//  [t^2/2, 1/2+t-t^2, (1-t)^2/2]
	for(int i = 0; i < 2; ++i){
		mid[i] = ceil(p[i]-0.5);
		t = mid[i]-(p[i]-0.5);

		w[i][0] = t*t/2; 
		w[i][2] = (1-t)*(1-t)/2; 
		w[i][1] = 1-w[i][0]-w[i][2];
	}

	// ノイズタイルを基底関数の値で重み付け補間する
	for(f[1] = -1; f[1] <= 1; ++f[1]){
		for(f[0] = -1; f[0] <= 1; ++f[0]){
			RXREAL weight = 1;
			for(int i = 0; i < 2; ++i){
				c[i] = ModW(mid[i]+f[i], n);
				weight *= w[i][f[i]+1];
			}
			result += weight*g_pNoiseTileData[c[1]*n+c[0]];
		}
	}

	return result;
}



/*!
 * Multiband Noise 2D
 * @param[in] p[2] 2D座標
 * @param[in] s
 * @param[in] firstBand 最初のバンド
 * @param[in] nbands バンド数
 * @param[in] w 重み
 * @return 
 */
RXREAL rxWaveletNoise::MultibandNoise2(RXREAL p[2], int first, int nbands, RXREAL *w0)
{
	// HACK:WMultibandNoise2
	RXREAL result = 0;
	RXREAL w = pow((RXREAL)2.0, (RXREAL)first);
	RXREAL q[2];

	for(int b = 0; b < nbands; ++b){
		q[0] = p[0]*w;
		q[1] = p[1]*w;
		result += w0[b]*Noise2(q);
		w *= 2.0;
	}

	RXREAL sigma_m = 0;
	for(int b = 0; b < nbands; ++b){
		sigma_m += w0[b]*w0[b];
	}

	// B-Splineの平均分散σNは2Dで0.265，3Dで0.210，2D上へ投影した3Dノイズで0.296
	RXREAL sigma_n = 0.265;
	sigma_m = sqrt(sigma_n*sigma_m);

	// 分散が1になるようにノイズを調整
	if(sigma_m) result /= sigma_m;

	return result;

}


/*!
 * derivative of 2D wavelet noise in x direction
 * @param[in] p[2] 2D座標
 * @return Waveletノイズ値
 */
RXREAL rxWaveletNoise::DNoise2(RXREAL p[2], int d)
{
	int f[3], c[3];		// filter, noise coef. indices
	int mid[3], n = g_iNoiseTileSize;

	RXREAL w[2][3], t, result = 0;

 	// 2次のB-スプライン(quadratic B-spline)基底関数を計算
	//  [t^2/2, 1/2+t-t^2, (1-t)^2/2]
	for(int i = 0; i < 2; ++i){
		mid[i] = ceil(p[i]-0.5);
		t = mid[i]-(p[i]-0.5);

		if(i == d){
			w[i][0] = -t; 
			w[i][1] = 2*t-1;
			w[i][2] = 1-t;
		}
		else{
			w[i][0] = t*t/2; 
			w[i][2] = (1-t)*(1-t)/2; 
			w[i][1] = 1-w[i][0]-w[i][2];
		}
	}

	// ノイズタイルを基底関数の値で重み付け補間する
	for(f[1] = -1; f[1] <= 1; ++f[1]){
		for(f[0] = -1; f[0] <= 1; ++f[0]){
			RXREAL weight = 1;
			for(int i = 0; i < 2; ++i){
				c[i] = ModW(mid[i]+f[i], n);
				weight *= w[i][f[i]+1];
			}
			result += weight*g_pNoiseTileData[c[1]*n+c[0]];
		}
	}

	return result;
}



/*!
 * 
 * @param[in] 
 * @return 
 */
void rxWaveletNoise::GenerateNoiseTile3(int n)
{
	if(g_pNoiseTileData != NULL && n == g_iNoiseTileSize) return;

	if(n%2) n++; // tile size must be even

	int sz = n*n*n;
	RXREAL *temp1 = new RXREAL[sz];
	RXREAL *temp2 = new RXREAL[sz];
	RXREAL *noise = new RXREAL[sz];

	// Step 1. Fill the tile with random numbers in the range -1 to 1.
	for(int i = 0; i < sz; ++i){
		noise[i] = GaussianNoise();
	}

	// Steps 2 and 3. Downsample and Upsample the tile
	for(int iy = 0; iy < n; ++iy){
		for(int iz = 0; iz < n; ++iz){
			// each x row
			int idx = iy*n+iz*n*n;
			Downsample(&noise[idx], &temp1[idx], n, 1);
			Upsample(  &temp1[idx], &temp2[idx], n, 1);
		}
	}

	for(int ix = 0; ix < n; ++ix){
		for(int iz = 0; iz < n; ++iz){
			// each y row
			int idx = ix+iz*n*n;
			Downsample(&temp2[idx], &temp1[idx], n, n);
			Upsample(  &temp1[idx], &temp2[idx], n, n);
		}
	}

	for(int ix = 0; ix < n; ++ix){
		for(int iy = 0; iy < n; ++iy){
			// each z row
			int idx = ix+iy*n;
			Downsample(&temp2[idx], &temp1[idx], n, n*n);
			Upsample(  &temp1[idx], &temp2[idx], n, n*n);
		}
	}

	// Step 4. Subtract out the coarse-scale contribution
	for(int i = 0; i < sz; ++i){
		noise[i] -= temp2[i];
	}

	// Avoid even/odd variance difference by adding odd-offset version of noise to itself.
	int offset = n/2;
	if(offset%2 == 0) offset++;

	for(int i = 0, ix = 0; ix < n; ++ix){
		for(int iy = 0; iy < n; ++iy){
			for(int iz = 0; iz < n; ++iz){
				temp1[i++] = noise[ModW(ix+offset, n)+ModW(iy+offset, n)*n+ModW(iz+offset, n)*n*n];
			}
		}
	}

	for(int i = 0; i < sz; ++i){
		noise[i] += temp1[i];
	}

	g_pNoiseTileData = noise;
	g_iNoiseTileSize = n;

	delete [] temp1;
	delete [] temp2;
}

/*!
 * 
 * @param[in] 
 * @return 
 */
RXREAL* rxWaveletNoise::GenerateNoiseTile3r(int &nx, int &ny, int &nz)
{
	if(nx%2) nx++; // tile size must be even
	if(ny%2) ny++; // tile size must be even
	if(nz%2) nz++; // tile size must be even

	int sz = nx*ny*nz;
	RXREAL *temp1 = new RXREAL[sz];
	RXREAL *temp2 = new RXREAL[sz];
	RXREAL *noise = new RXREAL[sz];

	init_genrand((unsigned)time(NULL));
	//init_genrand(1234);

	// Step 1. Fill the tile with random numbers in the range -1 to 1.
	for(int i = 0; i < sz; ++i){
		noise[i] = GaussianNoise();
	}

	// Steps 2 and 3. Downsample and Upsample the tile
	for(int iy = 0; iy < ny; ++iy){
		for(int iz = 0; iz < nz; ++iz){
			// each x row
			int idx = iy*nx+iz*nx*ny;
			Downsample(&noise[idx], &temp1[idx], nx, 1);
			Upsample(  &temp1[idx], &temp2[idx], nx, 1);
		}
	}

	for(int ix = 0; ix < nx; ++ix){
		for(int iz = 0; iz < nz; ++iz){
			// each y row
			int idx = ix+iz*nx*ny;
			Downsample(&temp2[idx], &temp1[idx], ny, nx);
			Upsample(  &temp1[idx], &temp2[idx], ny, nx);
		}
	}

	for(int ix = 0; ix < nx; ++ix){
		for(int iy = 0; iy < ny; ++iy){
			// each z row
			int idx = ix+iy*nx;
			Downsample(&temp2[idx], &temp1[idx], nz, nx*ny);
			Upsample(  &temp1[idx], &temp2[idx], nz, nx*ny);
		}
	}

	// Step 4. Subtract out the coarse-scale contribution
	for(int i = 0; i < sz; ++i){
		noise[i] -= temp2[i];
	}

	// Avoid even/odd variance difference by adding odd-offset version of noise to itself.
	int offset = nx/2;
	if(offset%2 == 0) offset++;

	for(int i = 0, ix = 0; ix < nx; ++ix){
		for(int iy = 0; iy < ny; ++iy){
			for(int iz = 0; iz < nz; ++iz){
				temp1[i++] = noise[ModW(ix+offset, nx)+ModW(iy+offset, ny)*nx+ModW(iz+offset, nz)*nx*ny];
			}
		}
	}

	for(int i = 0; i < sz; ++i){
		noise[i] += temp1[i];
	}

	delete [] temp1;
	delete [] temp2;

	return noise;
}

/*!
 * Non-projected 3D noise
 * @param[in] 
 * @return 
 */
RXREAL rxWaveletNoise::Noise3(RXREAL p[3])
{
	int f[3], c[3];	// filter, noise coef. indices
	int mid[3], n = g_iNoiseTileSize;
	RXREAL w[3][3], t, result = 0;

	// Evaluate quadratic B-spline basis functions
	for(int i = 0; i < 3; ++i){
		mid[i] = ceil(p[i]-0.5);
		t = mid[i]-(p[i]-0.5);

		w[i][0] = t*t/2;
		w[i][2] = (1-t)*(1-t)/2;
		w[i][1] = 1-w[i][0]-w[i][2];
	}

	// Evaluate noise by weighting noise coefficients by basis function values
	for(f[2] = -1; f[2] <= 1; ++f[2]){
		for(f[1] = -1; f[1] <= 1; ++f[1]){
			for(f[0] = -1; f[0] <= 1; ++f[0]){
				RXREAL weight = 1;
				for(int i = 0; i < 3; ++i){
					c[i] = ModW(mid[i]+f[i], n);
					weight *= w[i][f[i]+1];
				}
				result += weight*g_pNoiseTileData[c[2]*n*n+c[1]*n+c[0]];
			}
		}
	}
	
	return result;
}


/*!
 * derivative of 2D wavelet noise
 * @param[in] p[2] 2D座標
 * @return Waveletノイズ値
 */
RXREAL rxWaveletNoise::DNoise3(RXREAL p[3], int d)
{
	int f[3], c[3];		// filter, noise coef. indices
	int mid[3], n = g_iNoiseTileSize;

	RXREAL w[3][3], t, result = 0;

 	// 2次のB-スプライン(quadratic B-spline)基底関数を計算
	//  [t^2/2, 1/2+t-t^2, (1-t)^2/2]
	for(int i = 0; i < 3; ++i){
		mid[i] = ceil(p[i]-0.5);
		t = mid[i]-(p[i]-0.5);

		if(i == d){
			w[i][0] = -t; 
			w[i][1] = 2*t-1;
			w[i][2] = 1-t;
		}
		else{
			w[i][0] = t*t/2; 
			w[i][2] = (1-t)*(1-t)/2; 
			w[i][1] = 1-w[i][0]-w[i][2];
		}
	}

	// ノイズタイルを基底関数の値で重み付け補間する
	for(f[2] = -1; f[2] <= 1; ++f[2]){
		for(f[1] = -1; f[1] <= 1; ++f[1]){
			for(f[0] = -1; f[0] <= 1; ++f[0]){
				RXREAL weight = 1;
				for(int i = 0; i < 3; ++i){
					c[i] = ModW(mid[i]+f[i], n);
					weight *= w[i][f[i]+1];
				}
				result += weight*g_pNoiseTileData[c[2]*n*n+c[1]*n+c[0]];
			}
		}
	}

	return result;
}

/*!
* 3D noise projected onto 2D
* @param[in] 
* @return 
*/
RXREAL rxWaveletNoise::ProjectedNoise(RXREAL p[3], RXREAL normal[3])
{
	int c[3];	// noise coef. location
	int _min[3], _max[3];
	int n = g_iNoiseTileSize;

	RXREAL support, result = 0;

	// Bound the support of the basis functions for this projection direction
	for(int i = 0; i < 3; ++i){
		support = 3*abs(normal[i]) + 3*sqrt((1-normal[i]*normal[i])/2);
		_min[i] = ceil( p[i] - (3*abs(normal[i]) + 3*sqrt((1-normal[i]*normal[i])/2)) );
		_max[i] = floor( p[i] + (3*abs(normal[i]) + 3*sqrt((1-normal[i]*normal[i])/2)) );
	}


	// Loop over the noise coefficients within the bound.
	for(c[2] = _min[2]; c[2] <= _max[2]; ++c[2]){
		for(c[1] = _min[1]; c[1] <= _max[1]; ++c[1]){
			for(c[0] = _min[0]; c[0] <= _max[0]; ++c[0]){
				RXREAL t, t1, t2, t3, dot = 0, weight = 1;

				// Dot the normal with the vector from c to p */
				for(int i = 0; i < 3; ++i){
					dot += normal[i]*(p[i]-c[i]);
				}

				// Evaluate the basis function at c moved halfway to p along the normal. 
				for(int i = 0; i < 3; ++i){
					t = (c[i]+normal[i]*dot/2)-(p[i]-1.5);
					t1 = t-1; t2 = 2-t; t3 = 3-t;
					
					weight *= (t <= 0 || t >= 3) ? 0 : ((t < 1) ? t*t/2 : ((t < 2) ? 1-(t1*t1+t2*t2)/2 : t3*t3/2));
				}

				// Evaluate noise by weighting noise coefficients by basis function values.
				result += weight*g_pNoiseTileData[ModW(c[2], n)*n*n+ModW(c[1], n)*n+ModW(c[0], n)];
			}
		}
	}

	return result;
}


/*!
 * Multiband Noise 3D
 * @param[in] 
 * @return 
 */
RXREAL rxWaveletNoise::MultibandNoise3(RXREAL p[3], RXREAL s, RXREAL *normal, int firstBand, int nbands, RXREAL *w)
{
	RXREAL q[3], result = 0, variance = 0;
	for(int b = 0; b < nbands && s+firstBand+b < 0; ++b){
		for(int i = 0; i < 3; ++i){
			q[i] = 2*p[i]*pow((RXREAL)2.0, (RXREAL)(firstBand+b));
		}
		
		result += (normal) ? w[b]*ProjectedNoise(q, normal) : w[b]*Noise3(q);
	}
	for(int b = 0; b < nbands; ++b){
		variance += w[b]*w[b];
	}

	// Adjust the noise so it has a variance of 1.
	if(variance) result /= sqrt(variance*((normal) ? 0.296 : 0.210));

	return result;
}



/*!
 * 3Dノイズタイル生成
 * @param[in] n タイルグリッド数
 */
void rxWaveletNoise::GenerateNoiseTile4(int n, int nt)
{
	//if(g_pNoiseTileData != NULL && n == g_iNoiseTileSize) return;

	if(g_pNoiseTileData != NULL) delete [] g_pNoiseTileData;

	if(n%2) n++; // tile size must be even

	long sz = n*n*n*nt;
	RXREAL *temp1 = new RXREAL[sz];
	RXREAL *temp2 = new RXREAL[sz];
	RXREAL *noise = new RXREAL[sz];

	// Step 1. Fill the tile with random numbers in the range -1 to 1.
	for(int i = 0; i < sz; ++i){
		noise[i] = GaussianNoise();
	}

	// Steps 2 and 3. Downsample and Upsample the tile
	for(int iy = 0; iy < n; ++iy){
		for(int iz = 0; iz < n; ++iz){
			for(int it = 0; it < nt; ++it){
				// each x row
				long idx = iy*n+iz*n*n+it*n*n*n;
				Downsample(&noise[idx], &temp1[idx], n, 1);
				Upsample(  &temp1[idx], &temp2[idx], n, 1);
			}
		}
	}

	for(int ix = 0; ix < n; ++ix){
		for(int iz = 0; iz < n; ++iz){
			for(int it = 0; it < nt; ++it){
				// each y row
				long idx = ix+iz*n*n+it*n*n*n;
				Downsample(&temp2[idx], &temp1[idx], n, n);
				Upsample(  &temp1[idx], &temp2[idx], n, n);
			}
		}
	}

	for(int ix = 0; ix < n; ++ix){
		for(int iy = 0; iy < n; ++iy){
			for(int it = 0; it < nt; ++it){
				// each z row
				long idx = ix+iy*n+it*n*n*n;
				Downsample(&temp2[idx], &temp1[idx], n, n*n);
				Upsample(  &temp1[idx], &temp2[idx], n, n*n);
			}
		}
	}

	for(int ix = 0; ix < n; ++ix){
		for(int iy = 0; iy < n; ++iy){
			for(int iz = 0; iz < n; ++iz){
				// each t row
				long idx = ix+iy*n+iz*n*n;
				Downsample(&temp2[idx], &temp1[idx], nt, n*n*n);
				Upsample(  &temp1[idx], &temp2[idx], nt, n*n*n);
			}
		}
	}

	// Step 4. Subtract out the coarse-scale contribution
	for(int i = 0; i < sz; ++i){
		noise[i] -= temp2[i];
	}

	// Avoid even/odd variance difference by adding odd-offset version of noise to itself.
	int offset = n/2;
	if(offset%2 == 0) offset++;

	for(int i = 0, ix = 0; ix < n; ++ix){
		for(int iy = 0; iy < n; ++iy){
			for(int iz = 0; iz < n; ++iz){
				for(int it = 0; it < nt; ++it){
					temp1[i++] = noise[ModW(ix+offset, n)+ModW(iy+offset, n)*n+ModW(iz+offset, n)*n*n+ModW(it+offset, nt)*n*n*n];
				}
			}
		}
	}

	for(int i = 0; i < sz; ++i){
		noise[i] += temp1[i];
	}

	g_pNoiseTileData = noise;
	g_iNoiseTileSize = n;

	delete [] temp1;
	delete [] temp2;
}


/*!
 * 3Dノイズタイル生成
 * @param[in] n タイルグリッド数
 */
RXREAL* rxWaveletNoise::GenerateNoiseTile4r(int &n, int &nt)
{
	if(n%2) n++; // tile size must be even

	long sz = n*n*n*nt;
	RXREAL *temp1 = new RXREAL[sz];
	RXREAL *temp2 = new RXREAL[sz];
	RXREAL *noise = new RXREAL[sz];

	// Step 1. Fill the tile with random numbers in the range -1 to 1.
	for(int i = 0; i < sz; ++i){
		noise[i] = GaussianNoise();
	}

	// Steps 2 and 3. Downsample and Upsample the tile
	for(int iy = 0; iy < n; ++iy){
		for(int iz = 0; iz < n; ++iz){
			for(int it = 0; it < nt; ++it){
				// each x row
				long idx = iy*n+iz*n*n+it*n*n*n;
				Downsample(&noise[idx], &temp1[idx], n, 1);
				Upsample(  &temp1[idx], &temp2[idx], n, 1);
			}
		}
	}

	for(int ix = 0; ix < n; ++ix){
		for(int iz = 0; iz < n; ++iz){
			for(int it = 0; it < nt; ++it){
				// each y row
				long idx = ix+iz*n*n+it*n*n*n;
				Downsample(&temp2[idx], &temp1[idx], n, n);
				Upsample(  &temp1[idx], &temp2[idx], n, n);
			}
		}
	}

	for(int ix = 0; ix < n; ++ix){
		for(int iy = 0; iy < n; ++iy){
			for(int it = 0; it < nt; ++it){
				// each z row
				long idx = ix+iy*n+it*n*n*n;
				Downsample(&temp2[idx], &temp1[idx], n, n*n);
				Upsample(  &temp1[idx], &temp2[idx], n, n*n);
			}
		}
	}

	for(int ix = 0; ix < n; ++ix){
		for(int iy = 0; iy < n; ++iy){
			for(int iz = 0; iz < n; ++iz){
				// each t row
				long idx = ix+iy*n+iz*n*n;
				Downsample(&temp2[idx], &temp1[idx], nt, n*n*n);
				Upsample(  &temp1[idx], &temp2[idx], nt, n*n*n);
			}
		}
	}

	// Step 4. Subtract out the coarse-scale contribution
	for(int i = 0; i < sz; ++i){
		noise[i] -= temp2[i];
	}

	// Avoid even/odd variance difference by adding odd-offset version of noise to itself.
	int offset = n/2;
	if(offset%2 == 0) offset++;

	for(int i = 0, ix = 0; ix < n; ++ix){
		for(int iy = 0; iy < n; ++iy){
			for(int iz = 0; iz < n; ++iz){
				for(int it = 0; it < nt; ++it){
					temp1[i++] = noise[ModW(ix+offset, n)+ModW(iy+offset, n)*n+ModW(iz+offset, n)*n*n+ModW(it+offset, nt)*n*n*n];
				}
			}
		}
	}

	for(int i = 0; i < sz; ++i){
		noise[i] += temp1[i];
	}

	delete [] temp1;
	delete [] temp2;

	return noise;
}

/*!
 * Non-projected 3D noise
 * @param[in] 
 * @return 
 */
RXREAL rxWaveletNoise::Noise4(RXREAL p[4])
{
	int f[4], c[4];	// filter, noise coef. indices
	int mid[4], n = g_iNoiseTileSize;

	RXREAL w[4][3], t, result = 0;

	// Evaluate quadratic B-spline basis functions
	for(int i = 0; i < 4; ++i){
		mid[i] = ceil(p[i]-0.5);
		t = mid[i]-(p[i]-0.5);

		w[i][0] = t*t/2;
		w[i][2] = (1-t)*(1-t)/2;
		w[i][1] = 1-w[i][0]-w[i][2];
	}

	// Evaluate noise by weighting noise coefficients by basis function values
	for(f[3] = -1; f[3] <= 1; ++f[3]){
		for(f[2] = -1; f[2] <= 1; ++f[2]){
			for(f[1] = -1; f[1] <= 1; ++f[1]){
				for(f[0] = -1; f[0] <= 1; ++f[0]){
					RXREAL weight = 1;
					for(int i = 0; i < 4; ++i){
						c[i] = ModW(mid[i]+f[i], n);
						weight *= w[i][f[i]+1];
					}
					result += weight*g_pNoiseTileData[c[3]*n*n*n+c[2]*n*n+c[1]*n+c[0]];
				}
			}
		}
	}
	
	return result;
}

/*!
 * derivative of 4D wavelet noise
 * @param[in] p[4] 4D座標
 * @return Waveletノイズ値
 */
RXREAL rxWaveletNoise::DNoise4(RXREAL p[4], int d)
{
	int f[4], c[4];	// filter, noise coef. indices
	int mid[4], n = g_iNoiseTileSize;

	RXREAL w[4][3], t, result = 0;

	// Evaluate quadratic B-spline basis functions
	for(int i = 0; i < 4; ++i){
		mid[i] = ceil(p[i]-0.5);
		t = mid[i]-(p[i]-0.5);

		if(i == d){
			w[i][0] = -t; 
			w[i][1] = 2*t-1;
			w[i][2] = 1-t;
		}
		else{
			w[i][0] = t*t/2; 
			w[i][2] = (1-t)*(1-t)/2; 
			w[i][1] = 1-w[i][0]-w[i][2];
		}
	}

	// Evaluate noise by weighting noise coefficients by basis function values
	for(f[3] = -1; f[3] <= 1; ++f[3]){
		for(f[2] = -1; f[2] <= 1; ++f[2]){
			for(f[1] = -1; f[1] <= 1; ++f[1]){
				for(f[0] = -1; f[0] <= 1; ++f[0]){
					RXREAL weight = 1;
					for(int i = 0; i < 4; ++i){
						c[i] = ModW(mid[i]+f[i], n);
						weight *= w[i][f[i]+1];
					}
					result += weight*g_pNoiseTileData[c[3]*n*n*n+c[2]*n*n+c[1]*n+c[0]];
				}
			}
		}
	}
	
	return result;
}

/*!
 * Multiband Noise 2D
 * @param[in] p[2] 2D座標
 * @param[in] s
 * @param[in] firstBand 最初のバンド
 * @param[in] nbands バンド数
 * @param[in] w 重み
 * @return 
 */
RXREAL rxWaveletNoise::MultibandNoise4(RXREAL p[4], int first, int nbands, RXREAL *w0)
{
	// HACK:WMultibandNoise2
	RXREAL result = 0;
	RXREAL w = pow((RXREAL)2.0, (RXREAL)first);
	RXREAL q[4];

	for(int b = 0; b < nbands; ++b){
		q[0] = p[0]*w;
		q[1] = p[1]*w;
		q[2] = p[2]*w;
		q[3] = p[3]*w;
		result += w0[b]*Noise4(q);
		w *= 2.0;
	}

	RXREAL sigma_m = 0;
	for(int b = 0; b < nbands; ++b){
		sigma_m += w0[b]*w0[b];
	}

	// B-Splineの平均分散σNは2Dで0.265，3Dで0.210，2D上へ投影した3Dノイズで0.296
	RXREAL sigma_n = 0.210;
	sigma_m = sqrt(sigma_n*sigma_m);

	// 分散が1になるようにノイズを調整
	if(sigma_m) result /= sigma_m;

	return result;

}

