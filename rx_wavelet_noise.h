/*! 
 @file rx_wavelet_noise.h

 @brief Wavelet noise 生成

		Robert L. Cook and Tony DeRose
		Wavelet noise, 
		SIGGRAPH 2005, 2005. 
  
 @author Makoto Fujisawa
 @date 2010-01
*/


/*! @file
	
	@brief 
	@author Makoto Fujisawa
	@date 2010-01
*/



#ifndef _RX_WAVELET_NOISE_H_
#define _RX_WAVELET_NOISE_H_


//-----------------------------------------------------------------------------
// インクルードファイル
//-----------------------------------------------------------------------------
#include "rx_sph_commons.h"

//-----------------------------------------------------------------------------
// 名前空間
//-----------------------------------------------------------------------------
using namespace std;

//-----------------------------------------------------------------------------
// 変数
//-----------------------------------------------------------------------------
extern RXREAL *g_pNoiseTileData;
extern int g_iNoiseTileSize;

//-----------------------------------------------------------------------------
// MARK:rxWaveletNoise
//-----------------------------------------------------------------------------
namespace rxWaveletNoise
{
	void GenerateNoiseTile2(int n);
	void GenerateNoiseTile3(int n);
	void GenerateNoiseTile4(int n, int nt);

	//RXREAL* GenerateNoiseTile2r(int n);
	RXREAL* GenerateNoiseTile3r(int &n, int &ny, int &nz);
	RXREAL* GenerateNoiseTile4r(int &n, int &nt);

	RXREAL Noise3(RXREAL p[3]);
	RXREAL Noise2(RXREAL p[2]);
	RXREAL Noise4(RXREAL p[4]);

	RXREAL DNoise2(RXREAL p[2], int d);
	RXREAL DNoise3(RXREAL p[3], int d);
	RXREAL DNoise4(RXREAL p[4], int d);

	RXREAL MultibandNoise2(RXREAL p[2], int firstBand, int nbands, RXREAL *w);
	RXREAL MultibandNoise3(RXREAL p[3], RXREAL s, RXREAL *normal, int firstBand, int nbands, RXREAL *w);
	RXREAL MultibandNoise4(RXREAL p[4], int firstBand, int nbands, RXREAL *w0);

	RXREAL ProjectedNoise(RXREAL p[3], RXREAL normal[3]);
};





#endif // #ifndef _RX_WAVELET_NOISE_H_