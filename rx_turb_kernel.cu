/*! 
  @file rx_turb_kernel.cu
	
  @brief CUDA�ɂ��SPH�����v�Z
 
  @author Makoto Fujisawa
  @date 2009-08, 2011-06
*/

#ifndef _RX_TURBULENCE_KERNEL_CU_
#define _RX_TURBULENCE_KERNEL_CU_


//-----------------------------------------------------------------------------
// �C���N���[�h�t�@�C��
//-----------------------------------------------------------------------------
#include "rx_cu_common.cu"

//-----------------------------------------------------------------------------
// �ϐ�
//-----------------------------------------------------------------------------
#if USE_TEX
texture<float4, cudaTextureType1D, cudaReadModeElementType> dSortedPosTex;
texture<float4, cudaTextureType1D, cudaReadModeElementType> dSortedVelTex;
texture<uint,   cudaTextureType1D, cudaReadModeElementType> dCellStartTex;
texture<uint,   cudaTextureType1D, cudaReadModeElementType> dCellEndTex;

//�T�u�p�[�e�B�N��
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
 * �e�p�[�e�B�N���̃O���b�h�n�b�V���l
 * @param[out] gridParticleHash
 * @param[out] dSortedIndex
 * @param[in] pos �p�[�e�B�N���ʒu���i�[�����z��
 * @param[in] nprts �p�[�e�B�N����
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
// MARK:�E�F�[�u���b�g����
//-----------------------------------------------------------------------------

/*!
 * �^����ꂽ�Z�����̃p�[�e�B�N������A���E�F�[�u���b�g�ϊ����v�Z
 * @param[in] gridPos �O���b�h�ʒu
 * @param[in] index �p�[�e�B�N���C���f�b�N�X
 * @param[in] pos0 �v�Z���W
 * @param[in] vel0 �v�Z���W�ł̑��x
 * @param[in] cell �p�[�e�B�N���O���b�h�f�[�^
 */
__device__
float4 calParticleWt(int3 gridPos, uint i, int &c, float3 pos0, float3 vel0, rxParticleCell cell, float h, float scale)
{
	uint gridHash = calcGridHash(gridPos);

	// �Z�����̃p�[�e�B�N���̃X�^�[�g�C���f�b�N�X
	uint startIndex = FETCHC(dCellStart, gridHash);

	//float h = scale*MEXICAN_HAT_R;
	float4 wt = make_float4(0.0f);
	if(startIndex != 0xffffffff){	// �Z������łȂ����̃`�F�b�N
		// �Z�����̃p�[�e�B�N���Ŕ���
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
 * �p�[�e�B�N�����x�̃G�l���M�[�X�y�N�g�������v�Z
 * @param[out] Et �G�l���M�[�X�y�N�g����
 * @param[in] scale �E�F�[�u���b�g�X�P�[��
 * @param[in] coef_et �G�l���M�[�X�y�N�g�����l�ɂ�����W��
 * @param[in] max_et  �G�l���M�[�X�y�N�g�����l�̍ő�l(max_et�ȏ�Ȃ�N�����v)
 * @param[in] cell �p�[�e�B�N���O���b�h�f�[�^
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

	uint oIdx = cell.dSortedIndex[index];	// �\�[�g���Ă��Ȃ��Ƃ��̃C���f�b�N�X
	Et[oIdx] = et;
}

__device__
float wdnoise3d(float *tile, int n[3], int offset, int tile_size, float p[3], int d)
{
	int f[3], c[3];		// filter, noise coef. indices
	int mid[3];

	float w[3][3], t, result = 0;

	// 2����B-�X�v���C��(quadratic B-spline)���֐����v�Z
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

	// �m�C�Y�^�C�������֐��̒l�ŏd�ݕt����Ԃ���
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
// MARK:SPS����
//-----------------------------------------------------------------------------
/*!
 * ����ʒu�ƃO���b�h�̍ŒZ�������v�Z
 * @param[in] p		�O���b�h���W
 * @param[in] gridPos �O���b�h���W
 * @return ����
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
 * Quaternion�̊|���Z
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
 * �T�u���q(���x��0)�̃G�l���M�[�X�y�N�g���X�V
 * @param[out]	dSubEt	�T�u���q�̃G�l���M�[�X�y�N�g���ւ̃A�h���X
 * @param[in]	dEt		���q�̃G�l���M�[�X�y�N�g���ւ̃A�h���X
 * @param[in]	scale		���̃X�P�[��
 * @param[in]	radius		���x��0�ł̔��a
 * @param[in]	nprts �p�[�e�B�N����
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
 * �G�l���M�[�X�y�N�g���̍X�V
 * @param[out]	dSubBlockEtC	�T�u���q(�q)�̃G�l���M�[�X�y�N�g���ւ̃A�h���X
 * @param[in]	dSubBlockEtP	�T�u���q(�e)�̃G�l���M�[�X�y�N�g���ւ̃A�h���X
 * @param[in]	nprts �p�[�e�B�N����
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
 * �T�u���q(���x��0)�̈ʒu�X�V
 * @param[in]	dSubPos	�T�u���q�̐�΍��W�ւ̃A�h���X
 * @param[in]	dPos	���q�̐�΍��W�ւ̃A�h���X
 * @param[in]	nprts �p�[�e�B�N����
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
 * �T�u���q�̈ʒu�X�V
 * @param[in]	dSubBlockPosC	�T�u���q(�q)�̐�΍��W�ւ̃A�h���X
 * @param[in]	dSubBlockPosP	�T�u���q(�e)�̐�΍��W�ւ̃A�h���X
 * @param[in]	dSubBlockChild	�T�u���q(�e)�̎q1�ւ̒P�ʃx�N�g���ւ̃A�h���X
 * @param[in]	radius_level	�T�u���q(�e)�̔��a
 * @param[in]	nprts �p�[�e�B�N����
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
 * �T�u���q�̉�](�P�ʃx�N�g���X�V)
 * @param[inout]dSubBlockChild	�T�u���q�̎q1�ւ̒P�ʃx�N�g���̍s��
 * @param[in]	dSubBlockAxis	�T�u���q�̃G�l���M�[�X�y�N�g���̔z��
 * @param[in]	dSubBlockEt		�T�u���q�̃G�l���M�[�X�y�N�g���̔z��
 * @param[in]	radius_level	���̃T�u���q�̂̔��a
 * @param[in]	nprts �p�[�e�B�N����
 * @param[in]	dt				��������
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
 * �T�u���q�̉�](�P�ʃx�N�g���X�V)
 * @param[inout]dSubBlockChild	�T�u���q�̎q1�ւ̒P�ʃx�N�g���̍s��
 * @param[in]	dSubBlockAxis	�T�u���q�̃G�l���M�[�X�y�N�g���̔z��
 * @param[in]	dSubRand		uint����
 * @param[in]	nprts �p�[�e�B�N����
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
 * �T�u���q�̍X�V
 * @param[in]	dSubPosP		�T�u���q�̐e�̈ʒu
 * @param[in]	dSubPosC1		�T�u���q�̎q1�̈ʒu
 * @param[in]	dSubPosC2		�T�u���q�̎q2�̈ʒu
 * @param[in]	dSubChild		�T�u���q�̎q1�ւ̒P�ʃx�N�g��
 * @param[in]	dSubAxis		�T�u���q�̎�
 * @param[in]	dSubEt			�T�u���q�̃G�l���M�[�X�y�N�g��(�X�P�[���͒��a)
 * @param[in]	dSubRand		uint����
 * @param[in]	ratio_et		�G�l���M�[�X�y�N�g���Ɋ|����W��
 * @param[in]	radius_sub		���̃T�u���q�̔��a
 * @param[in]	nprts	�p�[�e�B�N����
 * @param[in]	dt				��������
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
 * @param[in]	dSubRand		uint����
 * @param[in]	nprts	�p�[�e�B�N����
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
 * @param[in]	dSubRand		uint����
 * @param[in]	nprts	�p�[�e�B�N����
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
 * �T�u���q(���x��0)�̈ʒu�X�V
 * @param[in]	dSubPos	�T�u���q�̐�΍��W�ւ̃A�h���X
 * @param[in]	dPos	���q�̐�΍��W�ւ̃A�h���X
 * @param[in]	nprts �p�[�e�B�N����
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
 * �T�u���q�̒ǉ�
 * @param[in] dSubPosP		�T�u���q�̐e�̈ʒu
 * @param[in] dSubPosC1		�T�u���q�̎q1�̈ʒu
 * @param[in] dSubPosC2		�T�u���q�̎q2�̈ʒu
 * @param[in] dSubChild		�T�u���q�̎q1�ւ̒P�ʃx�N�g��
 * @param[in] dSubAxis		�T�u���q�̎�
 * @param[in] dSubRand		uint����
 * @param[in] radius_sub	���̃T�u���q�̔��a
 * @param[in] nprts	�p�[�e�B�N����
 * @param[in] uMaxParticles	�ő�p�[�e�B�N����
 * @param[in] uStart		�ǉ��J�n�C���f�b�N�X
 * @param[in] uCount		�ǉ���
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
 * �T�u���q�̍X�V
 * @param[in]	dSubPosP		�T�u���q�̐e�̈ʒu
 * @param[in]	dSubPosC1		�T�u���q�̎q1�̈ʒu
 * @param[in]	dSubPosC2		�T�u���q�̎q2�̈ʒu
 * @param[in]	dSubChild		�T�u���q�̎q1�ւ̒P�ʃx�N�g��
 * @param[in]	dSubAxis		�T�u���q�̎�
 * @param[in]	dSubRand		uint����
 * @param[in]	radius_sub		���̃T�u���q�̔��a
 * @param[in]	nprts	�p�[�e�B�N����
 * @param[in]	uMaxParticles	�ő�p�[�e�B�N����
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
 * �G�l���M�[�Ɋ�Â��g�p����q�p�[�e�B�N����I��
 * @param[out] subCell �q�p�[�e�B�N�����
 * @param[in] dEt ���x��0�̃p�[�e�B�N���̃G�l���M�[�l
 * @param[in] dSubPos �q�p�[�e�B�N���̈ʒu
 * @param[in] ratio_e0 ���x��0�̃p�[�e�B�N���̃G�l���M�[�l�Ɋ|����W��
 */
__global__
void setSubUnsortArray(rxSubParticleCell subCell, float* dEt, float4* dSubPos, float ratio_et0)
{
	uint index = __umul24(blockIdx.x,blockDim.x)+threadIdx.x;
	if(index >= subCell.uNumParticles) return;

	uint sub_i, subIndex, subIndex1;
	float ratioL[MAX_SUB_LEVEL+1];	// �e���x���̉e���W��
	float ratio[MAX_SUB_NUM];
	float radius[MAX_SUB_NUM];
	//float4 posData[MAX_SUB_NUM];

	float ets = dEt[index];			// �p�[�e�B�N�������G�l���M�[
	float et_cri = subCell.fEtcri;	// �G�l���M�[��l

	// �e���x���̉e���W����������
	for(uint level = 0; level <= MAX_SUB_LEVEL; level++){
		ratioL[level] = 0.0;
	}

	// ���x��0�̃G�l���M�[
	float et0 = ratio_et0*ets;

	//�䗦����
	float setLevel = log(et_cri/et0)*INV_LOG_2_M5D9;

	// �e���x���̉e���W�����v�Z
	if(setLevel <= 0.0f){	// ���x��0�̃p�[�e�B�N���̂�
		ratioL[0] = 1.0f;
	}
	else if(setLevel >= (float)MAX_SUB_LEVEL){	// �ŉ��w�̎q�p�[�e�B�N���̂�
		ratioL[MAX_SUB_LEVEL] = powf(2.0f, -(float)MAX_SUB_LEVEL);
	}
	else {	// �q�p�[�e�B�N����(leveld��leveld+1�̊�)
		uint leveld = (uint)setLevel;
		ratioL[leveld]   = (setLevel-(float)leveld)*powf(2.0f, -(float)leveld) ;
		ratioL[leveld+1] = (1.0-ratioL[leveld])*powf(2.0f, -(float)(1+leveld)) ;
	}

	// �e���W���Ɣ��a�����[�J���������Ɋi�[
	sub_i = 0;
	for(uint level = 0; level <= MAX_SUB_LEVEL; ++level){
		for(uint lIndex = 0; lIndex < subCell.uSubNumEach[level]; ++lIndex){
			ratio[sub_i] = ratioL[level];
			radius[sub_i] = subCell.fSubRad[level];
			sub_i++;
		}
	}

	// �T�u�p�[�e�B�N�������i�[
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
 * �T�u�p�[�e�B�N�����l�߂� 
 */
__global__ 
void compactSubParticles(rxSubParticleCell subCell, uint num)
{
	// �C���f�b�N�X
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
 * ���x�����ƂɃT�u�p�[�e�B�N���f�[�^���\�[�g���āC
 * �n�b�V�����̊e�Z���̍ŏ��̃A�h���X������
 * @param[in] cell �p�[�e�B�N���O���b�h�f�[�^
 * @param[in] oldPos �p�[�e�B�N���ʒu
 * @param[in] oldVel �p�[�e�B�N�����x
 */
__global__
void reorderDataAndFindCellStartF4F1F1(rxSubParticleCell cell)
{
	extern __shared__ uint sharedHash[];	// blockSize+1 elements
	uint index = __umul24(blockIdx.x,blockDim.x)+threadIdx.x;
	
	uint hash;
	if(index < cell.uSubNumValidParticles){
		hash = cell.dSubGridParticleHash[index];	// �n�b�V���l

		sharedHash[threadIdx.x+1] = hash;	// �n�b�V���l���V�F�A�[�h�������Ɋi�[

		if(index > 0 && threadIdx.x == 0){
			// �e�V�F�A�[�h�������̍ŏ��ׂ͗̃p�[�e�B�N���̃n�b�V���l���i�[
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

		// �ʒu�Ƒ��x�̃f�[�^����ёւ��������i�[
		uint sortedIndex = cell.dSubSortedIndex[index];
		cell.dSubSortedPos[index] = FETCHC(dSubUnsortPos, sortedIndex);
		cell.dSubSortedRad[index] = FETCHC(dSubUnsortRad, sortedIndex);
		cell.dSubSortedRat[index] = FETCHC(dSubUnsortRat, sortedIndex);
	}
}

/*!
 * �^����ꂽ�Z�����̃p�[�e�B�N���Ƃ̋������疧�x���v�Z
 * @param[in] gridPos �O���b�h�ʒu
 * @param[in] pos0 �v�Z���W
 * @param[in] cell �p�[�e�B�N���O���b�h�f�[�^
 * @return �Z�����̃p�[�e�B�N������v�Z�������x�l
 */
__device__
float calSubDensityCellG(int3 gridPos,
						 float3 pos0, 
						 rxSubParticleCell cell)
{
	uint gridHash = calcGridHash(gridPos);

	// �Z�����̃p�[�e�B�N���̃X�^�[�g�C���f�b�N�X
	uint startIndex = FETCHC(dSubCellStart, gridHash);

	float d = 0.0f;
	if(startIndex != 0xffffffff){	// �Z������łȂ����̃`�F�b�N
		// �Z�����̃p�[�e�B�N���Ŕ���
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
 * �O���b�h��ł̖��x���v�Z
 * @param[out] GridD �O���b�h���x
 * @param[in] cell �p�[�e�B�N���O���b�h�f�[�^
 * @param[in] gnum �O���b�h��
 * @param[in] gmin �O���b�h�ŏ����W
 * @param[in] glen �O���b�h��
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
		uintIndex = dUintData[index];	// �n�b�V���l

		sharedUintData[threadIdx.x+1] = uintIndex;	// �n�b�V���l���V�F�A�[�h�������Ɋi�[

		if(index > 0 && threadIdx.x == 0){
			// �e�V�F�A�[�h�������̍ŏ��ׂ͗̃p�[�e�B�N���̃n�b�V���l���i�[
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



