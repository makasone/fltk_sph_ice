/*! 
  @file rx_cu_common.cu
	
  @brief CUDA���ʃf�o�C�X�֐�
 
  @author Makoto Fujisawa
  @date 2009-08, 2011-06
*/
// FILE --rx_cu_common.cu--

#ifndef _RX_CU_COMMON_CU_
#define _RX_CU_COMMON_CU_


//-----------------------------------------------------------------------------
// �C���N���[�h�t�@�C��
//-----------------------------------------------------------------------------
#include <stdio.h>
#include <math.h>

#include "helper_math.h"
#include <math_constants.h>

#include "rx_cuda_utils.h"

#include "rx_cu_common.cuh"


// �V�~�����[�V�����p�����[�^(�R���X�^���g������)
__constant__ rxSimParams params;



//-----------------------------------------------------------------------------
// �֐�
//-----------------------------------------------------------------------------
__device__ __host__
inline uint calUintPow(uint x, uint y)
{
	uint x_y = 1;
	for(uint i=0; i < y;i++) x_y *= x;
	return x_y;
}

/*!
 * a/b�̌v�Z���ʂ�؂�グ
 * @param[in] a,b a/b
 * @return �؂�グ�����Z����
 */
__device__ __host__
inline uint DivCeil(uint a, uint b)
{
	return (a % b != 0) ? (a / b + 1) : (a / b);
}


/*!
 * [a,b]�ɃN�����v
 * @param[in] x �N�����v���������l
 * @param[in] a,b �N�����v���E
 * @return �N�����v���ꂽ���l
 */
__device__
inline float CuClamp(float x, float a, float b)
{
	return max(a, min(b, x));
}
__device__
inline int CuClamp(int x, int a, int b)
{
	return max(a, min(b, x));
}

/*!
 * �[������ for float3
 * @param[in] v �l
 */
__device__
inline int CuIsZero(float3 v)
{
	if(fabsf(v.x) < 1.0e-10 && fabsf(v.y) < 1.0e-10 && fabsf(v.z) < 1.0e-10){
		return 1;
	}
	else{
		return 0;
	}
}

/*!
 * �s��ƃx�N�g���̐�
 * @param[in] m 3x3�s��
 * @param[in] v 3D�x�N�g��
 * @return �ς̌���
 */
__device__
inline float3 CuMulMV(matrix3x3 m, float3 v)
{
	return make_float3(dot(m.e[0], v), dot(m.e[1], v), dot(m.e[2], v));
}



// �O���b�h���u���b�N���C�u���b�N���X���b�h���̌v�Z
__device__ __host__
inline void computeGridSize(uint n, uint thread_per_block, uint &numBlocks, uint &numThreads)
{
	numThreads = min(thread_per_block, n);
	numBlocks = DivCeil(n, numThreads);
}


//-----------------------------------------------------------------------------
// �O���b�h
//-----------------------------------------------------------------------------
/*!
 * �O���b�h�ʒu�v�Z
 * @param[in] p ���W
 * @return �O���b�h���W
 */
__device__ 
inline int3 calcGridPos(float3 p)
{
	int3 gridPos;
	gridPos.x = floor((p.x-params.WorldOrigin.x)/params.CellWidth.x);
	gridPos.y = floor((p.y-params.WorldOrigin.y)/params.CellWidth.y);
	gridPos.z = floor((p.z-params.WorldOrigin.z)/params.CellWidth.z);

	gridPos.x = min(max(gridPos.x, 0), params.GridSize.x-1);
	gridPos.y = min(max(gridPos.y, 0), params.GridSize.y-1);
	gridPos.z = min(max(gridPos.z, 0), params.GridSize.z-1);

	return gridPos;
}

/*!
 * �O���b�h���W����1�����z�񒆂ł̈ʒu���v�Z
 * @param[in] gridPos �O���b�h���W
 * @return �A�h���X
 */
__device__ 
inline uint calcGridHash(int3 gridPos)
{
	return __umul24(__umul24(gridPos.z, params.GridSize.y), params.GridSize.x)+__umul24(gridPos.y, params.GridSize.x)+gridPos.x;
}

/*!
 * �O���b�h�ʒu�v�Z
 * @param[in] p ���W
 * @param[in] origin �O���b�h�̍ŏ����W
 * @param[in] cell_width 1�O���b�h�Z���̕�
 * @param[in] grid_size �O���b�h��
 * @return �O���b�h���W
 */
__device__ 
inline int3 calcGridPosB(float3 p, float3 origin, float3 cell_width, uint3 grid_size)
{
	int3 gridPos;
	gridPos.x = floor((p.x-origin.x)/cell_width.x);
	gridPos.y = floor((p.y-origin.y)/cell_width.y);
	gridPos.z = floor((p.z-origin.z)/cell_width.z);

	gridPos.x = min(max(gridPos.x, 0), grid_size.x-1);
	gridPos.y = min(max(gridPos.y, 0), grid_size.y-1);
	gridPos.z = min(max(gridPos.z, 0), grid_size.z-1);

	return gridPos;
}

/*!
 * �O���b�h���W����1�����z�񒆂ł̈ʒu���v�Z
 * @param[in] gridPos �O���b�h���W
 * @return �A�h���X
 */
__device__ 
inline uint calcGridHashB(int3 gridPos, uint3 grid_size)
{
	return __umul24(__umul24(gridPos.z, grid_size.y), grid_size.x)+__umul24(gridPos.y, grid_size.x)+gridPos.x;
}





//-----------------------------------------------------------------------------
// �A�g�~�b�N�֐�
//-----------------------------------------------------------------------------
#ifdef RX_USE_ATOMIC_FUNC

/*!
 * float��atomicAdd
 */
__device__ 
inline void atomicFloatAdd(float *address, float val)
{
	int i_val = __float_as_int(val);
	int tmp0 = 0;
	int tmp1;
 
	while( (tmp1 = atomicCAS((int *)address, tmp0, i_val)) != tmp0)
	{
		tmp0 = tmp1;
		i_val = __float_as_int(val + __int_as_float(tmp1));
	}
}
/*!
 * double��atomicAdd
 */
__device__ 
inline double atomicDoubleAdd(double *address, double val)
{
	unsigned long long int *address_as_ull = (unsigned long long int*)address;
	unsigned long long int old = *address_as_ull, assumed;
	do{
		assumed = old;
		old = atomicCAS(address_as_ull, assumed, __double_as_longlong(val+__longlong_as_double(assumed)));
	}while(assumed != old);
	return __longlong_as_double(old);
}
/*!
 * float��atomicMin
 */
__device__ 
inline float atomicFloatMin(float *address, float val)
{
	int *address_as_int = (int*)address;
	int old = atomicMin(address_as_int, __float_as_int(val));
	return __int_as_float(old);
}

/*!
 * float��atomicMax
 */
__device__ 
inline float atomicFloatMax(float *address, float val)
{
	int *address_as_int = (int*)address;
	int old = atomicMax(address_as_int, __float_as_int(val));
	return __int_as_float(old);
}

#endif // #ifdef RX_USE_ATOMIC_FUNC


//-----------------------------------------------------------------------------
// �O���b�h
//-----------------------------------------------------------------------------
/*!
 * 1D�C���f�b�N�X����3D�C���f�b�N�X�ւ̕ϊ�(�O���b�h���͔C��)
 * @param[in] i 1D�C���f�b�N�X
 * @param[in] gridSize �O���b�h��
 * @return 3D�C���f�b�N�X
 */
__device__
inline uint3 calcGridPosU(uint i, uint3 ngrid)
{
	uint3 gridPos;
	uint w = i%(ngrid.x*ngrid.y);
	gridPos.x = w%ngrid.x;
	gridPos.y = w/ngrid.x;
	gridPos.z = i/(ngrid.x*ngrid.y);
	return gridPos;
}
/*!
 * 3D�C���f�b�N�X����1D�C���f�b�N�X�ւ̕ϊ�(�O���b�h���͔C��)
 * @param[in] p 3D�C���f�b�N�X
 * @param[in] gridSize �O���b�h��
 * @return 1D�C���f�b�N�X
 */
__device__
inline uint calcGridPos3(uint3 p, uint3 ngrid)
{
	p.x = min(p.x, ngrid.x-1);
	p.y = min(p.y, ngrid.y-1);
	p.z = min(p.z, ngrid.z-1);
	return (p.z*ngrid.x*ngrid.y)+(p.y*ngrid.x)+p.x;
}



//-----------------------------------------------------------------------------
// CWT�f�o�C�X�֐�
//-----------------------------------------------------------------------------
/*!
 * ���L�V�J���n�b�g
 * @param[in] t ���W
 * @return �E�F�[�u���b�g��֐��l
 */
__device__
inline float MexicanHat(float t)
{
	t = t*t;
	return MEXICAN_HAT_C*(1.0-t)*exp(-t/2.0);
}
__device__
inline float MexicanHatIm(float t)
{
	return 0.0f;
}

/*!
 * ���L�V�J���n�b�g(�g�����)
 * @param[in] w �g��
 * @return �E�F�[�u���b�g��֐��l
 */
__device__
inline float MexicanHatWave(float w)
{
	w = w*w;
	return MEXICAN_HAT_C*M_SQRT2PI*w*exp(-w/2.0);
}
inline float MexicanHatWaveIm(float w)
{
	return 0.0f;
}

/*!
 * ���L�V�J���n�b�g(2D)
 * @param[in] x,y ���W
 * @return �E�F�[�u���b�g��֐��l
 */
__device__
inline float MexicanHat2D(float x, float y)
{
	x = x*x;
	y = y*y;
	return MEXICAN_HAT_C*(x+y-2)*exp(-(x+y)/2.0);
}
__device__
inline float MexicanHat2DIm(float x, float y)
{
	return 0.0f;
}

/*!
 * ���L�V�J���n�b�g(3D)
 * @param[in] x,y ���W
 * @return �E�F�[�u���b�g��֐��l
 */
__device__ __host__
inline float MexicanHat3D(float x, float y, float z)
{
	x = x*x;
	y = y*y;
	z = z*z;
	return MEXICAN_HAT_C*(x+y+z-3.0f)*exp(-(x+y+z)/2.0f);
}
__device__ __host__
inline float MexicanHat3DIm(float x, float y)
{
	return 0.0f;
}

__device__
inline int Mod(int x, int n)
{
	int m = (int)fmodf((float)x, (float)n); 
	return ((m < 0) ? m+n : m);
}


//-----------------------------------------------------------------------------
// ����
//-----------------------------------------------------------------------------
__device__ static mt_struct_stripped ds_MT[MT_RNG_COUNT];
static mt_struct_stripped h_MT[MT_RNG_COUNT];

/*!
 * Mersenne Twister �ɂ�闐������ (CUDA�T���v�����)
 * @param[out] d_Random ������������
 * @param[in] NPerRng ������
 */
__global__
static void RandomGPU(float *d_Random, int NPerRng)
{
	const int	  tid = blockDim.x * blockIdx.x + threadIdx.x;
	const int THREAD_N = blockDim.x * gridDim.x;

	int iState, iState1, iStateM, iOut;
	unsigned int mti, mti1, mtiM, x;
	unsigned int mt[MT_NN];

	for(int iRng = tid; iRng < MT_RNG_COUNT; iRng += THREAD_N){
		//Load bit-vector Mersenne Twister parameters
		mt_struct_stripped config = ds_MT[iRng];

		//Initialize current state
		mt[0] = config.seed;
		for(iState = 1; iState < MT_NN; iState++)
			mt[iState] = (1812433253U * (mt[iState - 1] ^ (mt[iState - 1] >> 30)) + iState) & MT_WMASK;

		iState = 0;
		mti1 = mt[0];
		for(iOut = 0; iOut < NPerRng; iOut++){
			//iState1 = (iState +	 1) % MT_NN
			//iStateM = (iState + MT_MM) % MT_NN
			iState1 = iState + 1;
			iStateM = iState + MT_MM;
			if(iState1 >= MT_NN) iState1 -= MT_NN;
			if(iStateM >= MT_NN) iStateM -= MT_NN;
			mti  = mti1;
			mti1 = mt[iState1];
			mtiM = mt[iStateM];

			x	= (mti & MT_UMASK) | (mti1 & MT_LMASK);
			x	=  mtiM ^ (x >> 1) ^ ((x & 1) ? config.matrix_a : 0);
			mt[iState] = x;
			iState = iState1;

			//Tempering transformation
			x ^= (x >> MT_SHIFT0);
			x ^= (x << MT_SHIFTB) & config.mask_b;
			x ^= (x << MT_SHIFTC) & config.mask_c;
			x ^= (x >> MT_SHIFT1);

			//Convert to (0, 1] float and write to global memory
			d_Random[iRng + iOut * MT_RNG_COUNT] = ((float)x + 1.0f) / 4294967296.0f;
		}
	}
}


// ���^�����@�ɂ�闐������(C����ȂǂƓ���)
__device__ static unsigned int randx = 1;

__device__
inline void Srand(unsigned int s)
{
	randx = s;
}

__device__
inline unsigned int Rand()
{
	randx = randx*1103515245+12345;
	return randx&2147483647;
}

__device__
inline unsigned int Rand2(unsigned int x)
{
	x = x*1103515245+12345;
	return x&2147483647;
}

#define RAND2_MAX (2147483647)



// XORShift�ɂ�闐��
__device__ static unsigned long xors_x = 123456789;
__device__ static unsigned long xors_y = 362436069;
__device__ static unsigned long xors_z = 521288629;
__device__ static unsigned long xors_w = 88675123;

/*!
 * G. Marsaglia, "Xorshift RNGs", Journal of Statistical Software, Vol. 8(14), pp.1-6, 2003. 
 *  - http://www.jstatsoft.org/v08/i14/
 * @param[in] 
 * @return 
 */
__device__
inline unsigned long Xorshift128()
{ 
	unsigned long t; 
	t = (xors_x^(xors_x<<11));
	xors_x = xors_y; xors_y = xors_z; xors_z = xors_w; 
	return ( xors_w = (xors_w^(xors_w>>19))^(t^(t>>8)) ); 
}
__device__
inline long Xorshift128(long l, long h)
{ 
	unsigned long t; 
	t = (xors_x^(xors_x<<11));
	xors_x = xors_y; xors_y = xors_z; xors_z = xors_w; 
	xors_w = (xors_w^(xors_w>>19))^(t^(t>>8));
	return l+(xors_w%(h-l));
}


__device__
inline float XorFrand(float l, float h)
{
	return l+(h-l)*(Xorshift128(0, 1000000)/1000000.0f);
}

__device__
inline void Random(float2 &x, float a, float b)
{
	x.x = XorFrand(a, b);
	x.y = XorFrand(a, b);
}

__device__
inline void Random(float3 &x, float a, float b)
{
	x.x = XorFrand(a, b);
	x.y = XorFrand(a, b);
	x.z = XorFrand(a, b);
}

// �K�E�X�m�C�Y
__device__
inline float GaussianNoise(void)
{
	float x1, x2;
	float ret;
	float r2;

	do {
		x1 = 2.0 * XorFrand(0.0, 1.0-(1e-10)) - 1.0;	/* [-1, 1) */
		x2 = 2.0 * XorFrand(0.0, 1.0-(1e-10)) - 1.0;

		r2 = x1*x1 + x2*x2;

	} while ((r2 == 0) || (r2 > 1.0));

	ret = x1 * sqrtf((-2.0 * logf(r2))/r2);
	ret *= 0.25;		// Possibility of ( N(0, 1) < 4.0 ) = 100%

	if (ret < -1.0) ret = -1.0; /* Account for loss of precision. */
	if (ret >  1.0) ret = 1.0;

	return ret;
}




//-----------------------------------------------------------------------------
// ��������
//-----------------------------------------------------------------------------

/*!
 * �����Ɖ~�̌�������(2D, A��)
 * @param[in] A,B �����̗��[�_���W
 * @param[in] C �~�̒��S
 * @param[in] r �~�̔��a
 * @param[out] P ��_���W
 * @return ��_��
 */
__device__ 
static int CuLineCircleIntersection(float2 A, float2 B, float2 C, float r, float2 P[2], float t[2])
{
	float rr = r*r;
	float2 AC = C-A;
	float2 BC = C-B;

	float2 v = B-A;
	float l = length(v);
	v /= l;

	float td = dot(v, AC);
	float2 D = A+td*v;
	float dd = dot(D-C, D-C);

	if(dd < rr){
		float dt = sqrtf(rr-dd);

		float da = rr-dot(AC, AC);
		float db = rr-dot(BC, BC);

		int inter = 0;
		float t1 = td-dt;
		float t2 = td+dt;
		if(t1 >= 0 && t1 <= l){
			P[inter] = A+t1*v;
			t[inter] = t1;
			inter++;
		}
		if(t2 >= 0 && t2 <= l){
			P[inter] = A+t2*v;
			t[inter] = t2;
			inter++;
		}

		return inter;
	}
	else{
		return 0;
	}
}


/*!
 * AABB�Ƌ��̋���
 * @param[in] spos �����S
 * @param[in] r �����a
 * @param[in] sgn
 * @param[in] box_min,box_max AABB�ŏ��C�ő���W�l
 * @param[out] cp AABB�\�ʂ̍ŋߖT�_
 * @param[out] d ���s��AABB�̋���
 * @param[out] n ��_�ɂ�����P�ʖ@���x�N�g��
 */
__device__
inline int collisionSphereAABB(float3 spos, float r, int sgn, float3 box_min, float3 box_max, float3 &cp, float &d, float3 &n)
{
	float3 dist_min;	// box_min�Ƃ̋���
	float3 dist_max;	// box_max�Ƃ̋���
	float d0 = 0.0f;
	float3 n0 = make_float3(0.0f, 0.0f, 0.0f);
	int bout = 0;
	int count = 0;

	// �e�����Ƃɍŏ��ƍő勫�E�O�ɂȂ��Ă��Ȃ������ׂ�
	if((dist_min.x = (spos.x-r)-box_min.x) < 0.0){ bout |= 0x0001; count++; d0 = dist_min.x; n0 = make_float3( 1.0,  0.0,  0.0);}
	if((dist_min.y = (spos.y-r)-box_min.y) < 0.0){ bout |= 0x0002; count++; d0 = dist_min.y; n0 = make_float3( 0.0,  1.0,  0.0);}
	if((dist_min.z = (spos.z-r)-box_min.z) < 0.0){ bout |= 0x0004; count++; d0 = dist_min.z; n0 = make_float3( 0.0,  0.0,  1.0);}
	if((dist_max.x = box_max.x-(spos.x+r)) < 0.0){ bout |= 0x0008; count++; d0 = dist_max.x; n0 = make_float3(-1.0,  0.0,  0.0);}
	if((dist_max.y = box_max.y-(spos.y+r)) < 0.0){ bout |= 0x0010; count++; d0 = dist_max.y; n0 = make_float3( 0.0, -1.0,  0.0);}
	if((dist_max.z = box_max.z-(spos.z+r)) < 0.0){ bout |= 0x0020; count++; d0 = dist_max.z; n0 = make_float3( 0.0,  0.0, -1.0);}

	// �����̓�(�S���ŋ��E��)
	if(bout == 0){
		float min_d = 1e10;
		if(dist_min.x < min_d){ min_d = dist_min.x; n = make_float3( 1.0,  0.0,  0.0); }
		if(dist_min.y < min_d){ min_d = dist_min.y; n = make_float3( 0.0,  1.0,  0.0); }
		if(dist_min.z < min_d){ min_d = dist_min.z; n = make_float3( 0.0,  0.0,  1.0); }

		if(dist_max.x < min_d){ min_d = dist_max.x; n = make_float3(-1.0,  0.0,  0.0); }
		if(dist_max.y < min_d){ min_d = dist_max.y; n = make_float3( 0.0, -1.0,  0.0); }
		if(dist_max.z < min_d){ min_d = dist_max.z; n = make_float3( 0.0,  0.0, -1.0); }

		d = (float)sgn*min_d;
		n *= (float)sgn;
		cp = spos+n*fabs(d);
		return 1;
	}

	// �����̊O
	// sgn = 1:���C-1:�I�u�W�F�N�g
	if(count == 1){
		// ���ʋߖT
		d = (float)sgn*d0;
		n = (float)sgn*n0;
		cp = spos+n*fabs(d);
	}
	else{
		// �G�b�W/�R�[�i�[�ߖT
		float3 x = make_float3(0.0f, 0.0f, 0.0f);
		if(bout & 0x0001) x.x =  dist_min.x;
		if(bout & 0x0002) x.y =  dist_min.y;
		if(bout & 0x0004) x.z =  dist_min.z;
		if(bout & 0x0008) x.x = -dist_max.x;
		if(bout & 0x0010) x.y = -dist_max.y;
		if(bout & 0x0020) x.z = -dist_max.z;

		d = length(x);
		n = normalize(x);

		d *= -(float)sgn;
		n *= -(float)sgn;

		cp = spos+n*fabs(d);

		float3 disp = make_float3(0.00001);
		//Random(disp, 0, 0.00001);
		disp = disp*n;
		cp += disp;
	}

	return 0;
}


/*!
 * AABB�Ɠ_�̋���
 * @param[in] p �_���W
 * @param[in] box_cen AABB�̒��S
 * @param[in] box_ext AABB�̊e�ӂ̒�����1/2
 * @param[out] cp AABB�\�ʂ̍ŋߖT�_
 * @param[out] d ���s��AABB�̋���
 * @param[out] n ��_�ɂ�����P�ʖ@���x�N�g��
 */
__device__
inline int collisionPointAABB(float3 p, float3 box_cen, float3 box_ext, float3 &cp, float &d, float3 &n)
{
	cp = p-box_cen;

	float3 tmp = fabs(cp)-box_ext;
	float res = ((tmp.x > tmp.y && tmp.x > tmp.z) ? tmp.x : (tmp.y > tmp.z ? tmp.y : tmp.z));

	float sgn = (res > 0.0) ? -1.0 : 1.0;

	int coli = 0;
	n = make_float3(0.0f);

	if(cp.x > box_ext.x){
		cp.x = box_ext.x;
		n.x -= 1.0;
		coli++;
	}
	else if(cp.x < -box_ext.x){
		cp.x = -box_ext.x;
		n.x += 1.0;
		coli++;
	}

	if(cp.y > box_ext.y){
		cp.y = box_ext.y;
		n.y -= 1.0;
		coli++;
	}
	else if(cp.y < -box_ext.y){
		cp.y = -box_ext.y;
		n.y += 1.0;
		coli++;
	}

	if(cp.z > box_ext.z){
		cp.z = box_ext.z;
		n.z -= 1.0;
		coli++;
	}
	else if(cp.z < -box_ext.z){
		cp.z = -box_ext.z;
		n.z += 1.0;
		coli++;
	}

	n = normalize(n);

	//if(coli > 1){
	//	float3 disp;
	//	Random(disp, 0, 0.00001);
	//	disp = disp*n;
	//	cp += disp;
	//}

	cp += box_cen;
	d = sgn*length(cp-p);

	return 0;
}


/*!
 * �_��BOX�̋���
 * @param[in] p �_���W
 * @param[in] box_cen BOX�̒��S
 * @param[in] box_ext BOX�̊e�ӂ̒�����1/2
 * @param[in] box_rot BOX�̕����s��(3x3��]�s��)
 * @param[in] box_inv_rot BOX�̕����s��̋t�s��(3x3)
 * @param[out] cp BOX�\�ʂ̍ŋߖT�_
 * @param[out] d �_��BOX�̋���
 * @param[out] n ��_�ɂ�����P�ʖ@���x�N�g��
 */
__device__
inline int collisionPointBox(float3 p, float3 box_cen, float3 box_ext, matrix3x3 box_rot, matrix3x3 box_inv_rot, float3 &cp, float &d, float3 &n)
{
	cp = p-box_cen;
	cp = CuMulMV(box_rot, cp);

	float3 tmp = fabs(cp)-box_ext;

	int coli = 0;
	n = make_float3(0.0f);

	if(tmp.x < 0.0 && tmp.y < 0.0 && tmp.z < 0.0){
		tmp = fabs(tmp);

		if(tmp.x <= tmp.y && tmp.x <= tmp.z){	// x���ʂɋ߂�
			if(cp.x > 0){
				cp.x = box_ext.x;
				n.x += 1.0;
			}
			else{
				cp.x = -box_ext.x;
				n.x -= 1.0;
			}
		}
		else if(tmp.y <= tmp.x && tmp.y <= tmp.z){ // y���ʂɋ߂�
			if(cp.y > 0){
				cp.y = box_ext.y;
				n.y += 1.0;
			}
			else{
				cp.y = -box_ext.y;
				n.y -= 1.0;
			}
		}
		else{ // z���ʂɋ߂�
			if(cp.z > 0){
				cp.z = box_ext.z;
				n.z += 1.0;
			}
			else{
				cp.z = -box_ext.z;
				n.z -= 1.0;
			}
		}

		coli++;
	}

	cp = CuMulMV(box_inv_rot, cp);
	n  = CuMulMV(box_inv_rot, n);

	n = normalize(n);
	cp += box_cen;

	float sgn = (coli) ? -1.0 : 1.0;
	d = sgn*(length(cp-p));

	return 0;
}

/*!
 * �_�Ƌ��̋���
 * @param[in] p �_���W
 * @param[in] sphere_cen ���̒��S
 * @param[in] sphere_rad ���̔��a
 * @param[out] cp �_�Ƌ����S�����Ԑ����Ƌ��̌�_
 * @param[out] d �_�Ƌ��\�ʂ̋���
 * @param[out] n �����S����_�ւ̒P�ʃx�N�g��
 */
__device__
inline int collisionPointSphere(float3 p, float3 sphere_cen, float sphere_rad, float3 &cp, float &d, float3 &n)
{
	n = make_float3(0.0f);

	float3 l = p-sphere_cen;
	float ll = length(l);

	d = ll-sphere_rad;
	if(d < 0.0){
		n = normalize(p-sphere_cen);
		cp = sphere_cen+n*sphere_rad;
	}

	return 0;
}

/*!
 * �_�ƕ��ʂ̋���
 * @param[in] v  �_�̍��W
 * @param[in] px ���ʏ�̓_
 * @param[in] pn ���ʂ̖@��
 * @return ����
 */
__device__ 
inline float distPointPlane(float3 v, float3 px, float3 pn)
{
	return dot((v-px), pn)/length(pn);
}

/*!
 * �O�p�`�Ɠ_�̋����ƍŋߖT�_
 * @param[in] v0,v1,v2	�O�p�`�̒��_
 * @param[in] n			�O�p�`�̖@��
 * @param[in] p			�_
 * @return 
 */
__device__ 
inline int distPointTriangle(float3 v0, float3 v1, float3 v2, float3 n, float3 p, float &dist, float3 &p0)
{
	// �|���S�����܂ޕ��ʂƓ_�̋���
	float l = distPointPlane(p, v0, n);
	
	// ���ʂƂ̍ŋߖT�_���W
	float3 np = p-l*n;

	// �ߖT�_���O�p�`�����ǂ����̔���
	float3 n1 = cross((v0-p), (v1-p));
	float3 n2 = cross((v1-p), (v2-p));
	float3 n3 = cross((v2-p), (v0-p));

	if(dot(n1, n2) > 0 && dot(n2, n3) > 0){
		// �O�p�`��
		dist = l;
		p0 = np;
		return 1;
	}
	else{
		// �O�p�`�O
		return 0;
	}
}


/*!
 * ���C/�����ƎO�p�`�̌���
 * @param[in] P0,P1 ���C/�����̒[�_or���C��̓_
 * @param[in] V0,V1,V2 �O�p�`�̒��_���W
 * @param[out] I ��_���W
 * @retval 1 ��_I�Ō��� 
 * @retval 0 ��_�Ȃ�
 * @retval 2 �O�p�`�̕��ʓ�
 * @retval -1 �O�p�`��"degenerate"�ł���(�ʐς�0�C�܂�C�������_�ɂȂ��Ă���)
 */
inline __device__ 
int intersectSegmentTriangle(float3 P0, float3 P1, 
							 float3 V0, float3 V1, float3 V2, 
							 float3 &I, float3 &n, float rp = 0.01)
{
	// �O�p�`�̃G�b�W�x�N�g���Ɩ@��
	float3 u = V1-V0;		
	float3 v = V2-V0;			
	n = normalize(cross(u, v));
	if(CuIsZero(n)){
		return -1;	// �O�p�`��"degenerate"�ł���(�ʐς�0)
	}

	// ����
	float3 dir = P1-P0;
	float a = dot(n, P0-V0);
	float b = dot(n, dir);
	if(fabs(b) < 1e-10){	// �����ƎO�p�`���ʂ����s
		if(a == 0){
			return 2;	// ���������ʏ�
		}
		else{
			return 0;	// ��_�Ȃ�
		}
	}


	// ��_�v�Z

	// 2�[�_�����ꂼ��قȂ�ʂɂ��邩�ǂ����𔻒�
	float r = -a/b;
	if(r < 0.0 || fabs(a) > fabs(b) || b > 0){
		return 0;
	}
	//if(r < 0.0){
	//	return 0;
	//}
	//else{
	//	if(fabs(a) > fabs(b)){
	//		return 0;
	//	}
	//	else{
	//		if(b > 0){
	//			return 0;
	//		}
	//	}
	//}

	// �����ƕ��ʂ̌�_
	I = P0+r*dir;

	// ��_���O�p�`���ɂ��邩�ǂ����̔���
	float uu, uv, vv, wu, wv, D;
	uu = dot(u, u);
	uv = dot(u, v);
	vv = dot(v, v);
	float3 w = I-V0;
	wu = dot(w, u);
	wv = dot(w, v);
	D = uv*uv-uu*vv;

	float s, t;
	s = (uv*wv-vv*wu)/D;
	if(s < 0.0 || s > 1.0){
		return 0;
	}
	
	t = (uv*wu-uu*wv)/D;
	if(t < 0.0 || (s+t) > 1.0){
		return 0;
	}

	return 1;
}



#endif // #ifndef _RX_CU_COMMON_CU_



