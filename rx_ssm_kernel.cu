/*! 
  @file rx_ssm_kernel.cu
	
  @brief SSM�@�ɂ�郁�b�V������
		- M. Muller et al., Screen space meshes, SCA2007, 2007. 

  @author Makoto Fujisawa
  @date 2011-05
*/
// FILE --rx_ssm_kernel.cu--

#ifndef _RX_SSM_KERNEL_CU_
#define _RX_SSM_KERNEL_CU_


//-----------------------------------------------------------------------------
// �C���N���[�h�t�@�C��
//-----------------------------------------------------------------------------
#include "rx_cu_common.cu"


//-----------------------------------------------------------------------------
// �ϐ�
//-----------------------------------------------------------------------------
// �e�N�X�`���Ɋi�[�����e�[�u��
texture<uint, 1, cudaReadModeElementType> g_texSSMeshTable;
texture<uint, 1, cudaReadModeElementType> g_texSSEdgeTable;
texture<uint, 1, cudaReadModeElementType> g_texSSNodeTable;
texture<uint, 1, cudaReadModeElementType> g_texSSVRotTable;

// �V�~�����[�V�����p�����[�^
__constant__ rxSsmParams g_cSSMParams;

// Binomial�t�B���^�W��(n_filter = RX_MAX_FILTER_SIZE�܂őΉ�)
__constant__ float g_cBinomials[RX_BINOMIALS_SIZE];

// �t�B���^�p�r�b�g��
__constant__ int BITS[] = {1, 2, 4, 8,  16, 32, 64, 128};


//-----------------------------------------------------------------------------
// MARK:SSM
//-----------------------------------------------------------------------------
/*!
 * float�^�̔z���������
 * @param[in] farray �l����������float�^�̔z��
 * @param[in] val �������l
 * @param[in] n �z��̃T�C�Y
 */
__global__
void initFloatArray(float* farray, float val, int n)
{
	uint index = __mul24(blockIdx.x,blockDim.x)+threadIdx.x;
	if(index >= n) return;

	farray[index] = val;
}

/*!
 * �p�[�e�B�N�����a�̔z���������
 * @param[in] rad �o�͔��a�z��(�T�C�Y=n*3)
 * @param[in] val ���͔��a�z��(�T�C�Y=n)
 * @param[in] n �z��̃T�C�Y
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
 * �f�v�X�}�b�v�����ŏ�����
 */
__global__
void initDepthMap(float* dmap, int n)
{
	uint index = __mul24(blockIdx.x,blockDim.x)+threadIdx.x;
	if(index >= n) return;

	dmap[index] = 1.0e10;
}


/*!
 * �f�v�X�}�b�v�̌v�Z
 */
__global__
void calDepthMap(float* dmap, float3* prts_pos, float3* prts_rad, int pnum, float3 tr, matrix4x4 PMV, 
				 float W, float H, float dw, float dh, int nx, int ny)
{
	uint index = __mul24(blockIdx.x,blockDim.x)+threadIdx.x;
	if(index >= pnum) return;	

	float4 x = make_float4(prts_pos[index], 1.0);
	float3 r = prts_rad[index];

	// ���e�ϊ�
	float4 xd;
	xd.x = dot(PMV.e[0], x);
	xd.y = dot(PMV.e[1], x);
	xd.z = dot(PMV.e[2], x);
	xd.w = dot(PMV.e[3], x);

	// w�Ŋ��邱�Ƃ�[-1, 1]�̐��K�����W�n�ɕϊ�
	float3 xp;
	xp.x = W*(0.5+0.5*xd.x/xd.w);
	xp.y = H*(0.5+0.5*xd.y/xd.w);
	xp.z = xd.z;

	prts_pos[index] = xp;

	// ���K�����W�n�ł̔��a�l
	float3 rp;
	rp.x = r.x*tr.x/xd.w/2;
	rp.y = r.y*tr.y/xd.w/2;
	rp.z = r.z*tr.z;

	prts_rad[index] = rp;

	float rrp = rp.x*rp.x;

	if(xp.x < 0 || xp.x >= W || xp.y < 0 || xp.y >= H){
		return;
	}

	// �f�v�X�}�b�v��ł̃p�[�e�B�N���͈̔�
	int cen[2];	// �p�[�e�B�N�����S
	cen[0] = int(xp.x/dw)+1;
	cen[1] = int(xp.y/dh)+1;

	int minp[2], maxp[2];
	minp[0] = cen[0]-(rp.x/dw+1);
	minp[1] = cen[1]-(rp.y/dh+1);
	maxp[0] = cen[0]+(rp.x/dw+1);
	maxp[1] = cen[1]+(rp.y/dh+1);

	// �͈͂��}�b�v�O�ɂȂ�Ȃ��悤�ɃN�����v
	minp[0] = (minp[0] < 0 ? 0 : (minp[0] > nx ? nx : minp[0]));
	minp[1] = (minp[1] < 0 ? 0 : (minp[1] > ny ? ny : minp[1]));
	maxp[0] = (maxp[0] < 0 ? 0 : (maxp[0] > nx ? nx : maxp[0]));
	maxp[1] = (maxp[1] < 0 ? 0 : (maxp[1] > ny ? ny : maxp[1]));

	// �p�[�e�B�N���f�v�X�l�X�V
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
 * �f�v�X�}�b�v�̌v�Z
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

	// ���e�ϊ�
	float4 xd;
	xd.x = dot(PMV.e[0], x);
	xd.y = dot(PMV.e[1], x);
	xd.z = dot(PMV.e[2], x);
	xd.w = dot(PMV.e[3], x);

	// w�Ŋ��邱�Ƃ�[-1, 1]�̐��K�����W�n�ɕϊ�
	float4 xp;
	xp.x = W*(0.5+0.5*xd.x/xd.w);
	xp.y = H*(0.5+0.5*xd.y/xd.w);
	xp.z = xd.z;
	xp.w = xd.w;

	prts_pos[index] = xp;

	// ���K�����W�n�ł̔��a�l
	float3 rp;
	rp.x = r.x*tr.x/xd.w/2;
	rp.y = r.y*tr.y/xd.w/2;
	rp.z = r.z*tr.z;

	prts_rad[index] = rp;

	float rrp = rp.x*rp.x;

	if(xp.x < 0 || xp.x >= W || xp.y < 0 || xp.y >= H){
		return;
	}

	// �f�v�X�}�b�v��ł̃p�[�e�B�N���͈̔�
	int cen[2];	// �p�[�e�B�N�����S
	cen[0] = int(xp.x/dw)+1;
	cen[1] = int(xp.y/dh)+1;

	int minp[2], maxp[2];
	minp[0] = cen[0]-(rp.x/dw+1);
	minp[1] = cen[1]-(rp.y/dh+1);
	maxp[0] = cen[0]+(rp.x/dw+1);
	maxp[1] = cen[1]+(rp.y/dh+1);

	// �͈͂��}�b�v�O�ɂȂ�Ȃ��悤�ɃN�����v
	minp[0] = (minp[0] < 0 ? 0 : (minp[0] > nx ? nx : minp[0]));
	minp[1] = (minp[1] < 0 ? 0 : (minp[1] > ny ? ny : minp[1]));
	maxp[0] = (maxp[0] < 0 ? 0 : (maxp[0] > nx ? nx : maxp[0]));
	maxp[1] = (maxp[1] < 0 ? 0 : (maxp[1] > ny ? ny : maxp[1]));

	// �p�[�e�B�N���f�v�X�l�X�V
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
 * �f�v�X�}�b�v����f�v�X�l���擾����(���E��������)
 * @param[in] data �f�v�X�}�b�v
 * @param[in] x,y  �s�N�Z���ʒu
 * @param[in] w,h  �f�v�X�}�b�v�̑傫��
 * @return �f�v�X�l
 */
__device__ 
float GetDepth(float *data, int x, int y, int w, int h)
{
	x = CuClamp(x, 0, w-1);
	y = CuClamp(y, 0, h-1);
	return data[y*w+x];
}


/*!
 * �f�v�X�}�b�v�𕽊���(x����)
 *  - �V�F�A�[�h���������g�p
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

	// Shared memory��̃C���f�b�N�X�擾
	#define SMEM(X, Y) d[(Y)*tilew+(X)]

	// �u���b�N���̃X���b�h�ʒu
	int tx = threadIdx.x;
	int ty = threadIdx.y;

	// �u���b�N�T�C�Y
	int bw = blockDim.x;
	int bh = blockDim.y;

	// �O���b�h���̃C���f�b�N�X(�摜��̈ʒu)
	int x = blockIdx.x*bw+tx;
	int y = blockIdx.y*bh+ty;

	if(x >= nx || y >= ny) return;

	// �t�B���^�[�͈͂̒l���V�F�A�[�h�������ɓ]��
	// ���S�̈�
	SMEM(r+tx, r+ty) = GetDepth(dmap, x, y, nx, ny);

	// �G�b�W�̈�
	if(threadIdx.x < r){
		SMEM(tx, r+ty)      = GetDepth(dmap, x-r,  y, nx, ny);	// �E
		SMEM(r+bw+tx, r+ty) = GetDepth(dmap, x+bw, y, nx, ny);	// ��
	}
	if(threadIdx.y < r){
		SMEM(r+tx, ty)      = GetDepth(dmap, x, y-r,  nx, ny);
		SMEM(r+tx, r+bh+ty) = GetDepth(dmap, x, y+bh, nx, ny);
	}

	// �R�[�i�[�̈�
	if((threadIdx.x < r) && (threadIdx.y < r)){
		SMEM(tx, ty)           = GetDepth(dmap, x-r,  y-r,  nx, ny);
		SMEM(tx, r+bh+ty)      = GetDepth(dmap, x-r,  y+bh, nx, ny);
		SMEM(r+bw+tx, ty)      = GetDepth(dmap, x+bh, y-r,  nx, ny);
		SMEM(r+bw+tx, r+bh+ty) = GetDepth(dmap, x+bw, y+bh, nx, ny);
	}

	// �S�X���b�h�̓]�����I���܂ő҂�
	__syncthreads();


	// �V�F�A�[�h���������ł̃O���b�h�C���f�b�N�X
	int i = r+tx;
	int j = r+ty;

	// x�����t�B���^
	float d0 = SMEM(i, j);
	if(d0 < 0.99e10){	// != ���̃s�N�Z���݂̂ɓK�p
		// ����r�O���b�h�œ��ꃌ�C���[�ɑ�����O���b�h���𒲂ׂ�
		int n = 0;
		for(int k = 1; k <= r; ++k){
			if(fabs(d0-SMEM(i-k, j)) > zmax || fabs(d0-SMEM(i+k, j)) > zmax){
				break;
			}
			n++;
		}

		// Binomial�W�����|���ĐώZ
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
 * �f�v�X�}�b�v�𕽊���(y����)
 *  - �V�F�A�[�h���������g�p
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

	// Shared memory��̃C���f�b�N�X�擾
	#define SMEM(X, Y) d[(Y)*tilew+(X)]

	// �u���b�N���̃X���b�h�ʒu
	int tx = threadIdx.x;
	int ty = threadIdx.y;

	// �u���b�N�T�C�Y
	int bw = blockDim.x;
	int bh = blockDim.y;

	// �O���b�h���̃C���f�b�N�X(�摜��̈ʒu)
	int x = blockIdx.x*bw+tx;
	int y = blockIdx.y*bh+ty;

	if(x >= nx || y >= ny) return;

	// �X�V�����f�v�X�l���ēx�V�F�A�[�h�������ɓ]��
	// ���S�̈�
	SMEM(r+tx, r+ty) = GetDepth(dmap, x, y, nx, ny);

	// �G�b�W�̈�
	if(threadIdx.x < r){
		SMEM(tx, r+ty)      = GetDepth(dmap, x-r,  y, nx, ny);	// �E
		SMEM(r+bw+tx, r+ty) = GetDepth(dmap, x+bw, y, nx, ny);	// ��
	}
	if(threadIdx.y < r){
		SMEM(r+tx, ty)      = GetDepth(dmap, x, y-r,  nx, ny);
		SMEM(r+tx, r+bh+ty) = GetDepth(dmap, x, y+bh, nx, ny);
	}

	// �R�[�i�[�̈�
	if((threadIdx.x < r) && (threadIdx.y < r)){
		SMEM(tx, ty)           = GetDepth(dmap, x-r,  y-r,  nx, ny);
		SMEM(tx, r+bh+ty)      = GetDepth(dmap, x-r,  y+bh, nx, ny);
		SMEM(r+bw+tx, ty)      = GetDepth(dmap, x+bh, y-r,  nx, ny);
		SMEM(r+bw+tx, r+bh+ty) = GetDepth(dmap, x+bw, y+bh, nx, ny);
	}

	// �S�X���b�h�̓]�����I���܂ő҂�
	__syncthreads();


	// �V�F�A�[�h���������ł̃O���b�h�C���f�b�N�X
	int i = r+tx;
	int j = r+ty;

	// y�����t�B���^
	float d0 = SMEM(i, j);
	if(d0 < 0.99e10){	// != ���̃s�N�Z���݂̂ɓK�p
		// ����r�O���b�h�œ��ꃌ�C���[�ɑ�����O���b�h���𒲂ׂ�
		int n = 0;
		for(int k = 1; k <= r; ++k){
			if(fabs(d0-SMEM(i, j-k)) > zmax || fabs(d0-SMEM(i, j+k)) > zmax){
				break;
			}
			n++;
		}

		// Binomial�W�����|���ĐώZ
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
 * x�����֊s�G�b�W�𒊏o
 */
__global__
void detectSilhouetteEdgeX(float* dmap, int nx, int ny, float dw, float dh, float zmax, rxSSEdgeG *dedge, uint *dedgesil)
{
	// �O���b�h���̃C���f�b�N�X
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
 * y�����֊s�G�b�W�𒊏o
 */
__global__
void detectSilhouetteEdgeY(float* dmap, int nx, int ny, float dw, float dh, float zmax, rxSSEdgeG *dedge, uint *dedgesil)
{
	// �O���b�h���̃C���f�b�N�X
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
 * �֊s�G�b�W�����l�߂� 
 * @param[out] compacted_edges �֊s�G�b�W���l�߂��z��
 * @param[in] silhouette �֊s�G�b�W���(�֊s�G�b�W:1, ����ȊO:0)
 * @param[in] silhouette_Scan �֊s�G�b�W��񂩂�쐬����Prefix Sum(Scan)
 * @param[in] edges �S�G�b�W
 * @param[in] num ���G�b�W��
 */
__global__ 
void compactSilhouetteEdges(rxSSEdgeG *compacted_edges, uint *silhouette, uint *silhouette_scan, rxSSEdgeG *edges, uint num)
{
	// �C���f�b�N�X
	uint blockId = __mul24(blockIdx.y, gridDim.x)+blockIdx.x;
	uint i = __mul24(blockId, blockDim.x)+threadIdx.x;

	if(silhouette[i] && (i < num)) {
		compacted_edges[silhouette_scan[i]] = edges[i];
	}
}



/*!
 * �֊s�G�b�W��front vertex���Z�o
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

	// �f�v�X�}�b�v��ł̃p�[�e�B�N���͈̔�
	int cen[2];	// �p�[�e�B�N�����S
	cen[0] = int(xp.x/dw)+1;
	cen[1] = int(xp.y/dh)+1;

	int minp[2], maxp[2];
	minp[0] = cen[0]-(rp.x/dw+1);
	minp[1] = cen[1]-(rp.y/dh+1);
	maxp[0] = cen[0]+(rp.x/dw+1);
	maxp[1] = cen[1]+(rp.y/dh+1);

	// �͈͂��}�b�v�O�ɂȂ�Ȃ��悤�ɃN�����v
	minp[0] = CuClamp(minp[0], 0, nx-1);
	minp[1] = CuClamp(minp[1], 0, ny-1);
	maxp[0] = CuClamp(maxp[0], 0, nx-1);
	maxp[1] = CuClamp(maxp[1], 0, ny-1);

	// �͈͓��̃G�b�W�𒲍�(x����)
	for(int j = minp[1]; j <= maxp[1]+1; ++j){
		for(int i = minp[0]; i <= maxp[0]; ++i){
			rxSSEdgeG &e = edges[i+j*(nx)];
			int sil = silhouette[i+j*(nx)];
			if(!sil) continue;
			if(e.depth < xp.z) continue;

			// �~�ƃG�b�W�̌�_
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
				// HACK:�v����
				atomicFloatMax(&e.dx, t[0]);
				if(e.dx == t[0]){
					edge_vertices_occ[i+j*(nx)] = 1;
					edge_vertices[i+j*(nx)] = make_float3(P[0].x, P[0].y, xp.z);
					e.front_vertex = 1;
				}
			}
		}
	}

	// �͈͓��̃G�b�W�𒲍�(y����)
	for(int i = minp[0]; i <= maxp[0]+1; ++i){
		for(int j = minp[1]; j <= maxp[1]; ++j){
			rxSSEdgeG &e = edges[yoffset+i+j*(nx+1)];
			int sil = silhouette[yoffset+i+j*(nx+1)];
			if(!sil) continue;
			if(e.depth < xp.z) continue;

			// �~�ƃG�b�W�̌�_
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
 * �֊s�G�b�W��front vertex���Z�o
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

	// �f�v�X�}�b�v��ł̃p�[�e�B�N���͈̔�
	int cen[2];	// �p�[�e�B�N�����S
	cen[0] = int(xp.x/dw)+1;
	cen[1] = int(xp.y/dh)+1;

	int minp[2], maxp[2];
	minp[0] = cen[0]-(rp.x/dw+1);
	minp[1] = cen[1]-(rp.y/dh+1);
	maxp[0] = cen[0]+(rp.x/dw+1);
	maxp[1] = cen[1]+(rp.y/dh+1);

	// �͈͂��}�b�v�O�ɂȂ�Ȃ��悤�ɃN�����v
	minp[0] = CuClamp(minp[0], 0, nx-1);
	minp[1] = CuClamp(minp[1], 0, ny-1);
	maxp[0] = CuClamp(maxp[0], 0, nx-1);
	maxp[1] = CuClamp(maxp[1], 0, ny-1);

	// �͈͓��̃G�b�W�𒲍�(x����)
	for(int j = minp[1]; j <= maxp[1]+1; ++j){
		for(int i = minp[0]; i <= maxp[0]; ++i){
			rxSSEdgeG &e = edges[i+j*(nx)];
			int sil = silhouette[i+j*(nx)];
			if(!sil) continue;
			if(e.depth < xp.z) continue;

			// �~�ƃG�b�W�̌�_
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
				// HACK:�v����
				atomicFloatMax(&e.dx, t[0]);
				if(e.dx == t[0]){
					edge_vertices_occ[i+j*(nx)] = 1;
					edge_vertices[i+j*(nx)] = make_float3(P[0].x, P[0].y, xp.z);
					e.front_vertex = 1;
				}
			}
		}
	}

	// �͈͓��̃G�b�W�𒲍�(y����)
	for(int i = minp[0]; i <= maxp[0]+1; ++i){
		for(int j = minp[1]; j <= maxp[1]; ++j){
			rxSSEdgeG &e = edges[yoffset+i+j*(nx+1)];
			int sil = silhouette[yoffset+i+j*(nx+1)];
			if(!sil) continue;
			if(e.depth < xp.z) continue;

			// �~�ƃG�b�W�̌�_
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
 * �G�b�W���_���l�߂� 
 * @param[out] compacted_ev �G�b�W���_���l�߂��z��
 * @param[in] ev_occ �G�b�W���_�L�����(�G�b�W���_�L��:1, ����:0)
 * @param[in] ev_occ_scan �G�b�W���_�L����񂩂�쐬����Prefix Sum(Scan)
 * @param[in] ev �G�b�W���_
 * @param[in] edges �G�b�W
 * @param[in] num ���G�b�W��
 */
__global__ 
void compactEdgeVertices(float3 *compacted_ev, uint *ev_occ, uint *ev_occ_scan, float3 *ev, rxSSEdgeG *edges, uint num)
{
	// �C���f�b�N�X
	uint blockId = __mul24(blockIdx.y, gridDim.x)+blockIdx.x;
	uint i = __mul24(blockId, blockDim.x)+threadIdx.x;

	if(ev_occ[i] && (i < num)) {
		compacted_ev[ev_occ_scan[i]] = ev[i];
		edges[i].front_vertex = ev_occ_scan[i];
	}
}



/*!
 * �m�[�h���_�̐���
 */
__global__
void calNodeVertex(float* dmap, int nx, int ny, float dw, float dh, float zmax, float3 *node, uint *node_occ)
{
	// �O���b�h���̃C���f�b�N�X
	int x = blockIdx.x*blockDim.x+threadIdx.x;
	int y = blockIdx.y*blockDim.y+threadIdx.y;

	if(x >= nx || y >= ny) return;

	int idx = x+y*nx;
	float d = dmap[idx];
	if(d < 0.99e10){	// != ���̃m�[�h�ɒ��_����
		node[idx] = make_float3(dw*x, dh*y, d);
		node_occ[idx] = 1;
	}
}

/*!
 * �m�[�h���_���l�߂� 
 */
__global__ 
void compactNodeVertex(float3 *compacted_nv, uint *nv_occ, uint *nv_occ_scan, float3 *nv, uint num)
{
	// �C���f�b�N�X
	uint blockId = __mul24(blockIdx.y, gridDim.x)+blockIdx.x;
	uint i = __mul24(blockId, blockDim.x)+threadIdx.x;

	if(nv_occ[i] && (i < num)) {
		compacted_nv[nv_occ_scan[i]] = nv[i];
	}
}


/*!
 * �m�[�h���_�C���f�b�N�X�����b�V���O���b�h�Ɋi�[
 */
__global__
void storeNodeVertex(rxSSGridG* mgrid, int nx, int ny, float3 *node_vertices, uint *node_occ, uint *node_occ_scan)
{
	// �O���b�h�C���f�b�N�X
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

	// �����m�[�h
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

	// �E���m�[�h
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

	// ����m�[�h
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

	// �E��m�[�h
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
 * �֊s�G�b�W��back vertex���Z�o
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


	// back vertex�����邩�ǂ������`�F�b�N
	if(e.x0.z < 0.99e10 && e.x1.z < 0.99e10){	// �G�b�W�[�_�������Ƃ� != �� �Ȃ�� back vertex ������
		int back_node;		// back layer�ɑ�����G�b�W�[�_
		//int nn_back_node;	// �אڃm�[�h

		// �f�v�X�l���傫������ back layer �ɑ�����
		if(e.x0.z > e.x1.z){
			back_node = x+y*(nx+1);
			//nn_back_node = (x == 0 ? x : x-1)+y*(nx+1);
		}
		else{
			back_node = (x+1)+y*(nx+1);
			//nn_back_node = ((x+1) == nx ? x+1 : x+2)+y*(nx+1);
		}

		// back layer�ɑ�����G�b�W�[�_�Ƃ��̗אڃm�[�h�̃f�v�X�l
		float back_node_depth = dmap[back_node];
		//float nn_back_node_depth = dmap[nn_back_node];

		// �f�v�X�l���O�}�ɂ��ߎ�
		float back_depth = back_node_depth;
		//float back_depth = back_node_depth*((2*l-dx)/l)-nn_back_node_depth*((l-dx)/l);

		float3 vrt_pos = front_vertices[x+y*nx];

		// back vertex��ݒ�
		float3 back_vrt;
		back_vrt.x = vrt_pos.x;
		back_vrt.y = vrt_pos.y;
		back_vrt.z = back_depth;

		back_vertices_occ[x+y*nx] = 1;
		back_vertices[x+y*nx] = back_vrt;
	}
}

/*!
 * �֊s�G�b�W��back vertex���Z�o
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

	// back vertex�����邩�ǂ������`�F�b�N
	if(e.x0.z < 0.99e10 && e.x1.z < 0.99e10){	// �G�b�W�[�_�������Ƃ� != �� �Ȃ�� back vertex ������
		int back_node;		// back layer�ɑ�����G�b�W�[�_
		//int nn_back_node;	// �אڃm�[�h

		// �f�v�X�l���傫������ back layer �ɑ�����
		if(e.x0.z > e.x1.z){
			back_node = x+y*nx;
			//nn_back_node = (x == 0 ? x : x-1)+y*nx;
		}
		else{
			back_node = x+(y+1)*nx;
			//nn_back_node = ((x+1) == nx ? x+1 : x+2)+y*nx;
		}

		// back layer�ɑ�����G�b�W�[�_�Ƃ��̗אڃm�[�h�̃f�v�X�l
		float back_node_depth = dmap[back_node];
		//float nn_back_node_depth = dmap[nn_back_node];

		// �f�v�X�l���O�}�ɂ��ߎ�
		float back_depth = back_node_depth;
		//float back_depth = back_node_depth*((2*l-dx)/l)-nn_back_node_depth*((l-dx)/l);

		float3 vrt_pos = front_vertices[x+y*nx+offset];

		// back vertex��ݒ�
		float3 back_vrt;
		back_vrt.x = vrt_pos.x;
		back_vrt.y = vrt_pos.y;
		back_vrt.z = back_depth;

		back_vertices_occ[x+y*nx+offset] = 1;
		back_vertices[x+y*nx+offset] = back_vrt;
	}
}

/*!
 * �G�b�W���_���l�߂� 
 * @param[out] compacted_ev �G�b�W���_���l�߂��z��
 * @param[in] ev_occ �G�b�W���_�L�����(�G�b�W���_�L��:1, ����:0)
 * @param[in] ev_occ_scan �G�b�W���_�L����񂩂�쐬����Prefix Sum(Scan)
 * @param[in] ev �G�b�W���_
 * @param[in] num ���G�b�W��
 */
__global__ 
void compactBackEdgeVertices(float3 *compacted_ev, uint *ev_occ, uint *ev_occ_scan, float3 *ev, uint num)
{
	// �C���f�b�N�X
	uint blockId = __mul24(blockIdx.y, gridDim.x)+blockIdx.x;
	uint i = __mul24(blockId, blockDim.x)+threadIdx.x;

	if(ev_occ[i] && (i < num)) {
		compacted_ev[ev_occ_scan[i]] = ev[i];
	}
}


/*!
 * �G�b�W���_�C���f�b�N�X�����b�V���O���b�h�Ɋi�[
 */
__global__
void storeEdgeVertex(rxSSGridG* mgrid, int nx, int ny, int yoffset, int num_nv, int num_nev, 
					 float3 *edge_vrts, uint *edge_occ, uint *edge_occ_scan, 
					 float3 *back_vrts, uint *back_occ, uint *back_occ_scan)
{
	// �O���b�h�C���f�b�N�X
	int x = blockIdx.x*blockDim.x+threadIdx.x;
	int y = blockIdx.y*blockDim.y+threadIdx.y;

	if(x >= nx || y >= ny) return;

	int idx = x+y*nx;
	rxSSGridG &g = mgrid[idx];

	g.num_ev = 0;
	g.num_bv = 0;
	g.back2 = -1;

	//
	// �O�ʃG�b�W���_
	//
	// ���G�b�W
	if(edge_occ[x+y*nx]){
		int eidx = edge_occ_scan[x+y*nx]+num_nv;
		g.edge_vrts[0] = eidx;
		g.num_ev++;
	}
	else{
		g.edge_vrts[0] = -1;
	}

	// ��G�b�W
	if(edge_occ[x+(y+1)*nx]){
		int eidx = edge_occ_scan[x+(y+1)*nx]+num_nv;
		g.edge_vrts[2] = eidx;
		g.num_ev++;
	}
	else{
		g.edge_vrts[2] = -1;
	}

	// ���G�b�W
	if(edge_occ[yoffset+x+y*(nx+1)]){
		int eidx = edge_occ_scan[yoffset+x+y*(nx+1)]+num_nv;
		g.edge_vrts[3] = eidx;
		g.num_ev++;
	}
	else{
		g.edge_vrts[3] = -1;
	}

	// �E�G�b�W
	if(edge_occ[yoffset+(x+1)+y*(nx+1)]){
		int eidx = edge_occ_scan[yoffset+(x+1)+y*(nx+1)]+num_nv;
		g.edge_vrts[1] = eidx;
		g.num_ev++;
	}
	else{
		g.edge_vrts[1] = -1;
	}

	//
	// �w�ʃG�b�W���_
	//
	// ���G�b�W
	if(back_occ[x+y*nx]){
		int eidx = back_occ_scan[x+y*nx]+num_nev;
		g.back_vrts[0] = eidx;
		g.num_bv++;
	}
	else{
		g.back_vrts[0] = -1;
	}

	// ��G�b�W
	if(back_occ[x+(y+1)*nx]){
		int eidx = back_occ_scan[x+(y+1)*nx]+num_nev;
		g.back_vrts[2] = eidx;
		g.num_bv++;
	}
	else{
		g.back_vrts[2] = -1;
	}

	// ���G�b�W
	if(back_occ[yoffset+x+y*(nx+1)]){
		int eidx = back_occ_scan[yoffset+x+y*(nx+1)]+num_nev;
		g.back_vrts[3] = eidx;
		g.num_bv++;
	}
	else{
		g.back_vrts[3] = -1;
	}

	// �E�G�b�W
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
 * �p�^�[��0 �������b�V���p�̃e�[�u���C���f�b�N�X�X�V
 * @param[inout] table_index ���̃e�[�u���C���f�b�N�X(���4�r�b�g���G�b�W���_,����4�r�b�g���m�[�h���_)
 * @param[inout] vrot ���_���[�e�[�V������
 * @param[in] v �O���b�h���̒��_�C���f�b�N�X
 * @param[inout] g ���b�V���𐶐�����O���b�h
 */
__device__ 
void CuUpdateTableIndexE0N4(int &table_index, int &vrot, int v[], rxSSGridG *g)
{
	// �O���b�h�ԍ��ɂ��90�x��]������
	table_index -= ((g->i+g->j) & 0x01) ? 1 : 0;
}

/*!
 * �p�^�[��1,2,4,8 �����֊s(�֊s�̎n�_)�p�̃e�[�u���C���f�b�N�X�X�V
 * @param[inout] table_index ���̃e�[�u���C���f�b�N�X(���4�r�b�g���G�b�W���_,����4�r�b�g���m�[�h���_)
 * @param[inout] vrot ���_���[�e�[�V������
 * @param[in] v �O���b�h���̒��_�C���f�b�N�X
 * @param[inout] g ���b�V���𐶐�����O���b�h
 */
__device__ 
void CuUpdateTableIndexE1(int &table_index, int &vrot, int v[], rxSSGridG *g)
{
	// �m�[�h���_�� back layer �ɑ����镨��T��
	// ���_������G�b�W�̒[�_�Ńf�v�X�l���傫������ back layer
	int kb = 0;
	for(int k = 0; k < 4; ++k){
		if(v[k+4] != -1){
			int k1 = (k == 3 ? 0 : k+1);
			kb = ((g->node_depth[k] > g->node_depth[k1]) ? BITS[k] : BITS[k1]);
			break;
		}
	}

	kb = (~kb & 0x0F);
	

	// �o�b�N���C���[�ɂ���m�[�h���_�r�b�g��0�ɂ���
	table_index = (table_index & 0xF0)+kb;
}

/*!
 * �p�^�[��3,5,6,9,10,12 �����֊s�p�̃e�[�u���C���f�b�N�X�X�V
 * @param[inout] table_index ���̃e�[�u���C���f�b�N�X(���4�r�b�g���G�b�W���_,����4�r�b�g���m�[�h���_)
 * @param[inout] vrot ���_���[�e�[�V������
 * @param[in] v �O���b�h���̒��_�C���f�b�N�X
 * @param[inout] g ���b�V���𐶐�����O���b�h
 */
__device__ 
void CuUpdateTableIndexE2N4(int &table_index, int &vrot, int v[], rxSSGridG *g)
{
	// back vertex�̃m�[�h���_�r�b�g��0�Ƃ����r�b�g��𐶐�
	int btable = 0;
	int k0 = 0, k1;
	for(int k = 0; k <= 4; ++k){
		k1 = (k0 == 3 ? 0 : k0+1);
		if(v[k0+4] != -1){	// �m�[�h�ԂɃG�b�W���_����
			btable |= (g->node_depth[k0] > g->node_depth[k1]) ? BITS[k0] : BITS[k1];
		}
		else{
			btable |= ((btable & BITS[k0]) ? BITS[k1] : 0);
		}
		k0++;
		if(k0 == 4) k0 = 0;
	}

	// �O���֊s�Ƌ�ʂ��邽�߂ɁC�����֊s�̏ꍇ�r�b�g���+2����
	btable = (btable+2 & 0x0F);

	table_index = (table_index & 0xF0)+btable;
}

/*!
 * 7,11,13,14 �O��/�����֊s���ݗp�̃e�[�u���C���f�b�N�X�X�V
 * @param[inout] table_index ���̃e�[�u���C���f�b�N�X(���4�r�b�g���G�b�W���_,����4�r�b�g���m�[�h���_)
 * @param[inout] vrot ���_���[�e�[�V������
 * @param[in] v �O���b�h���̒��_�C���f�b�N�X
 * @param[inout] g ���b�V���𐶐�����O���b�h
 */
__device__ 
void CuUpdateTableIndexE3N23(int &table_index, int &vrot, int v[], rxSSGridG *g, float3 *vrts, float zmax, rxVrtAdd &add)
{
	int btable = 0;	// �m�[�h���_��2,3�̏ꍇ�C1�Ԗڂ̃r�b�g��0(0xxx)
	int pattern = (table_index >> 4);
	int node = tex1Dfetch(g_texSSEdgeTable, pattern*4)-4;
	int R[4];	// ���_���Ȃ��G�b�W���甽���v���Ƀm�[�h���_�ԍ�����ׂ����X�g
	for(int k = 0; k < 4; ++k){
		R[k] = node++;
		if(node == 4) node = 0;
	}

	int ntable = (table_index & 0x0F);	// �m�[�h���_�r�b�g��

	// 2�Ԗڂ̃r�b�g
	btable |= (((ntable >> R[1]) & 1) ? 4 : 0);

	// 3�Ԗڂ̃r�b�g
	btable |= (((ntable >> R[2]) & 1) ? 2 : 0);

	// R�ɏ]��btable����ёւ�
	int btable0 = 0;
	btable0 |= (((ntable >> R[0]) & 1) ? 8 : 0);
	btable0 |= (((ntable >> R[1]) & 1) ? 4 : 0);
	btable0 |= (((ntable >> R[2]) & 1) ? 2 : 0);
	btable0 |= (((ntable >> R[3]) & 1) ? 1 : 0);

	// ���ёւ���btable�̉���2�r�b�g-1
	int n0 = (btable0 & 3)-1;
	int n1 = n0+1;

	add.num = 0;
	if(n0 != 1){
		// �ǉ���back vertex��g_EdgeTable[pattern][1]�̈ʒu�ɒǉ�
		int add_edge = tex1Dfetch(g_texSSEdgeTable, pattern*4+1)-4;	// �ǉ�����G�b�W�̈ʒu
		float3 add_vrt = vrts[g->edge_vrts[add_edge]];		// �ǉ����_

		// �ǉ����_�f�v�X�l
		int ref_edge = tex1Dfetch(g_texSSEdgeTable, pattern*4)-4;		// �ǉ����_�̃f�v�X�l���Q�Ƃ���G�b�W���_
		int ref_node = (btable & 4) ? R[1] : R[2];
		if(fabs(vrts[g->edge_vrts[ref_edge]].z-g->node_depth[ref_node]) > zmax){
			// edge vertex���g�p
			add_vrt.z = vrts[g->edge_vrts[ref_edge]].z;	// �ǉ����_�̃f�v�X�l
		}
		else{
			// back vertex���g�p
			if(g->back_vrts[ref_edge] == -1){
				ref_edge = tex1Dfetch(g_texSSEdgeTable, pattern*4+2)-4;
			}
			add_vrt.z = vrts[g->back_vrts[ref_edge]].z;	// �ǉ����_�̃f�v�X�l
		}

		add.edge[0] = add_edge;
		add.vrts[0] = add_vrt;
		add.num = 1;

		g->back_vrts[4] = add_edge;
		
		// �O���b�h�̒��_���X�g�ɒǉ�
		if(add_vrt.z < vrts[v[add_edge+4]].z){
			add.layer = 0;
		}
		else{
			add.layer = 1;
		}
			
	}

	// 4�Ԗڂ̃r�b�g : �m�[�h�̃f�v�X�l���r
	btable |= ((g->node_depth[R[n0]] >= g->node_depth[R[n1]]) ? 0 : 1);

	// �e�[�u���C���f�b�N�X�̃G�b�W���_�������X�V
	table_index &= 0xF0;
	table_index |= btable;
}

/*!
 * �p�^�[��7,11,13,14 �����֊s�p�̃e�[�u���C���f�b�N�X�X�V
 *  - ����Ȃ�back vertex���K�v
 * @param[inout] table_index ���̃e�[�u���C���f�b�N�X(���4�r�b�g���G�b�W���_,����4�r�b�g���m�[�h���_)
 * @param[inout] vrot ���_���[�e�[�V������
 * @param[in] v �O���b�h���̒��_�C���f�b�N�X
 * @param[inout] g ���b�V���𐶐�����O���b�h
 */
__device__ 
void CuUpdateTableIndexE3N4(int &table_index, int &vrot, int v[], rxSSGridG *g, float3 *vrts, float zmax, rxVrtAdd &add)
{
	int btable = 8;	// �m�[�h���_��4�̏ꍇ�C1�Ԗڂ̃r�b�g��1(1xxx)
	int pattern = (table_index >> 4);

	// �m�[�h�̃f�v�X�l���ׂ鏇�Ԃ̌���
	int node = tex1Dfetch(g_texSSEdgeTable, pattern*4)-4;
	int R[4];	// ���_���Ȃ��G�b�W���甽���v���Ƀm�[�h���_�ԍ�����ׂ����X�g
	for(int k = 0; k < 4; ++k){
		R[k] = node++;
		if(node == 4) node = 0;
	}

	// �m�[�h�̃f�v�X�l�̑召�Ńr�b�g���ύX(2-4�Ԗڂ̃r�b�g)
	for(int k = 0; k < 3; ++k){
		// R[k] > R[k+1]�Ȃ�ΑΉ�����r�b�g��1�ɂ���
		if(g->node_depth[R[k]] > g->node_depth[R[k+1]]){
			btable |= BITS[2-k];
		}
	}

	// �e�[�u���C���f�b�N�X�̃G�b�W���_�������X�V
	table_index &= 0xF0;
	table_index |= btable;

	//
	// g_EdgeTable[pattern][1]�̈ʒu�ɒ��_��ǉ�
	//
	int add_edge = tex1Dfetch(g_texSSEdgeTable, pattern*4+1)-4;	// �ǉ�����G�b�W�̈ʒu
	float3 add_vrt = vrts[g->edge_vrts[add_edge]];		// �ǉ����_

	// ���b�V��C�̒��_�����ԂŃf�v�X�l�����߂�
	float ref_depths[4];
	// �^�_
	// 2�[3
	// |�_|
	// 0�[1

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

	// �ǉ����_�̃f�v�X�l
	add_vrt.z = 0.5*(ref_depths[2]+ref_depths[3]);

	add.edge[0] = add_edge;
	add.vrts[0] = add_vrt;
	add.num = 1;

	g->back_vrts[4] = add_edge;

	// �O���b�h�̒��_���X�g�ɒǉ�
	if(add_vrt.z < vrts[v[add_edge+8]].z){
		if(add_vrt.z < vrts[v[add_edge+4]].z){
			// front vertex�Ƃ��đ}��
			add.layer = 0;
		}
		else{
			// back vertex�Ƃ��đ}��
			add.layer = 1;
		}
	}
	else{
		// back-2 vertex�Ƃ��đ}��
		add.layer = 2;
	}	
}

/*!
 * �p�^�[��15 �O���֊s�p�̃e�[�u���C���f�b�N�X�X�V
 * @param[inout] table_index ���̃e�[�u���C���f�b�N�X(���4�r�b�g���G�b�W���_,����4�r�b�g���m�[�h���_)
 * @param[inout] vrot ���_���[�e�[�V������
 * @param[in] v �O���b�h���̒��_�C���f�b�N�X
 * @param[inout] g ���b�V���𐶐�����O���b�h
 */
__device__ 
void CuUpdateTableIndexE4N2(int &table_index, int &vrot, int v[], rxSSGridG *g)
{
	int ntable = (table_index & 0x0F);	// �m�[�h���_�r�b�g��
	ntable = (ntable == 5 ? 0 : 15);

	// �e�[�u���C���f�b�N�X�̃m�[�h���_�������X�V
	table_index &= 0xF0;
	table_index |= ntable;
}

/*!
 * �p�^�[��15 �O��/�����֊s���ݗp�̃e�[�u���C���f�b�N�X�X�V
 * @param[inout] table_index ���̃e�[�u���C���f�b�N�X(���4�r�b�g���G�b�W���_,����4�r�b�g���m�[�h���_)
 * @param[inout] vrot ���_���[�e�[�V������
 * @param[in] v �O���b�h���̒��_�C���f�b�N�X
 * @param[inout] g ���b�V���𐶐�����O���b�h
 */
__device__ 
void CuUpdateTableIndexE4N3(int &table_index, int &vrot, int v[], rxSSGridG *g)
{
	int btable = 0;	// �m�[�h���_��3�̏ꍇ�C1�Ԗڂ̃r�b�g��0(0xxx)
	int ntable = (table_index & 0x0F);	// �m�[�h���_�r�b�g��

	// ���_���Ȃ��m�[�h
	int zero_node = log((double)(~ntable & 0x0F))/log(2.0);

	// �m�[�h�̃f�v�X�l�̑召�Ńr�b�g���ύX(2-4�Ԗڂ̃r�b�g)
	for(int k = 0; k < 3; ++k){
		int k0 = tex1Dfetch(g_texSSNodeTable, zero_node*8+k);
		int k1 = tex1Dfetch(g_texSSNodeTable, zero_node*8+k+1);

		// k0 > k1�Ȃ�ΑΉ�����r�b�g��1�ɂ���
		if(g->node_depth[k0] > g->node_depth[k1]){
			btable |= BITS[2-k];
		}
	}

	// ���_���[�e�[�V����
	vrot = zero_node;

	// �e�[�u���C���f�b�N�X�̃m�[�h���_�������X�V
	table_index &= 0xF0;
	table_index |= btable;
}


/*!
 * �p�^�[��15 �����֊s�p�̃e�[�u���C���f�b�N�X�X�V
 *  - �ǉ���back vertex��2�K�v
 * @param[inout] table_index ���̃e�[�u���C���f�b�N�X(���4�r�b�g���G�b�W���_,����4�r�b�g���m�[�h���_)
 * @param[inout] vrot ���_���[�e�[�V������
 * @param[in] v �O���b�h���̒��_�C���f�b�N�X
 * @param[inout] g ���b�V���𐶐�����O���b�h
 */
__device__ 
void CuUpdateTableIndexE4N4(int &table_index, int &vrot, int v[], rxSSGridG *g, float3 *vrts, float zmax, rxVrtAdd &add)
{
	int btable = 8;	// �m�[�h���_��4�̏ꍇ�C1�Ԗڂ̃r�b�g��1(1xxx)
	//int ntable = (table_index & 0x0F);	// �m�[�h���_�r�b�g��

	// �f�v�X�l�������Ƃ��傫���m�[�h
	int zero_node = 0;
	double max_depth = 0.0;
	for(int k = 1; k < 4; ++k){
		if(g->node_depth[k] > max_depth){
			max_depth = g->node_depth[k];
			zero_node = k;
		}
	}

	// �m�[�h�̃f�v�X�l�̑召�Ńr�b�g���ύX(2-4�Ԗڂ̃r�b�g)
	for(int k = 0; k < 3; ++k){
		int k0 = tex1Dfetch(g_texSSNodeTable, zero_node*8+k);
		int k1 = tex1Dfetch(g_texSSNodeTable, zero_node*8+k+1);

		// k0 > k1�Ȃ�ΑΉ�����r�b�g��1�ɂ���
		if(g->node_depth[k0] > g->node_depth[k1]){
			btable |= BITS[2-k];
		}
	}

	// ���_���[�e�[�V����
	vrot = zero_node;

	// �e�[�u���C���f�b�N�X�̃G�b�W���_�������X�V
	table_index &= 0xF0;
	table_index |= btable;

	// 
	// g_NodeTable[zero_node][5,6]�̈ʒu�ɒ��_(back-2 vertex)��ǉ�
	//
	int add_edge1 = tex1Dfetch(g_texSSNodeTable, zero_node*8+5)-4;	// �ǉ�����G�b�W�̈ʒu
	int add_edge2 = tex1Dfetch(g_texSSNodeTable, zero_node*8+6)-4;	// �ǉ�����G�b�W�̈ʒu
	float3 add_vrt1 = vrts[g->edge_vrts[add_edge1]];		// �ǉ����_
	float3 add_vrt2 = vrts[g->edge_vrts[add_edge2]];		// �ǉ����_

	// g_NodeTable[zero_node][4,7]�̈ʒu��back vertex�̃f�v�X��ݒ�
	int ref_edge1 = tex1Dfetch(g_texSSNodeTable, zero_node*8+4)-4;	// �ǉ�����G�b�W�̈ʒu
	int ref_edge2 = tex1Dfetch(g_texSSNodeTable, zero_node*8+7)-4;	// �ǉ�����G�b�W�̈ʒu
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
 * �O���b�h���ƂɃ��b�V������
 */
__global__
void calGridMesh(rxSSGridG* mgrid, int nx, int ny, float zmax, float3 *vrts, int num_vrts, uint *tri_num, 
				 float3 *back2_vrts, uint *back2_occ, int yoffset)
{
	// �O���b�h�C���f�b�N�X
	int x = blockIdx.x*blockDim.x+threadIdx.x;
	int y = blockIdx.y*blockDim.y+threadIdx.y;

	if(x >= nx || y >= ny) return;

	rxSSGridG &g = mgrid[x+y*nx];

	int v[14];
	// �m�[�h���_
	// 3 - 2
	// |   |
	// 0 - 1
		 
	// �G�b�W���_
	// - 6 -
	// 7   5
	// - 4 -
		 
	// �G�b�W���_(back vertex)
	//  - 10 -
	// 11     9
	//  -  8 -
	// back-2 vertex : 12,13

	int table_index = 0;
	for(int k = 0; k < 4; ++k){
		v[k]   = g.node_vrts[k];
		v[k+4] = g.edge_vrts[k];
		v[k+8] = g.back_vrts[k];

		table_index |= ((v[k] != -1) ? BITS[k] : 0);		// �m�[�h���_����4�r�b�g
		table_index |= ((v[k+4] != -1) ? BITS[k]*16 : 0);	// �G�b�W���_���4�r�b�g
	}
	v[12] = -1;
	v[13] = -1;

	int rotation = 0;

	g.table_index0 = table_index;	// �f�o�b�O�p

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
		// ���_�ǉ�
		CuUpdateTableIndexE3N23(table_index, rotation, v, &g, vrts, zmax, add);
	}
	else if(fidx == 19){
		// ���_�ǉ�
		CuUpdateTableIndexE3N4(table_index, rotation, v, &g, vrts, zmax, add);
	}
	else if(fidx == 22){
		CuUpdateTableIndexE4N2(table_index, rotation, v, &g);
	}
	else if(fidx == 23){
		CuUpdateTableIndexE4N3(table_index, rotation, v, &g);
	}
	else if(fidx == 24){
		// ���_�ǉ�
		CuUpdateTableIndexE4N4(table_index, rotation, v, &g, vrts, zmax, add);
	}

	g.back2 = add.num;
	if(add.num){
		if(add.edge[0]%2){
			// y�����G�b�W
			int idx = (x+(add.edge[0]-1)/2)+y*(nx+1)+yoffset;
			back2_occ[idx] = 1;
			back2_vrts[idx] = add.vrts[0];
		}
		else{
			// x�����G�b�W
			int idx = x+(y+add.edge[0]/2)*nx;
			back2_occ[idx] = 1;
			back2_vrts[idx] = add.vrts[0];
		}
		g.back2 += (add.layer << 2);

		if(add.num == 2){
			if(add.edge[1]%2){
				// y�����G�b�W
				int idx = (x+(add.edge[1]-1)/2)+y*(nx+1)+yoffset;
				back2_occ[idx] = 1;
				back2_vrts[idx] = add.vrts[1];
			}
			else{
				// x�����G�b�W
				int idx = x+(y+add.edge[1]/2)*nx;
				back2_occ[idx] = 1;
				back2_vrts[idx] = add.vrts[1];
			}
			g.back2 += (add.layer << 4);
		}
	}

	g.table_index1 = table_index;	// �f�o�b�O�p

	int num_tri = tex1Dfetch(g_texSSMeshTable, table_index*19);
	if(num_tri > 0){	// �O���b�h���̃��b�V������0���傫�������烁�b�V������
		g.mesh_num = num_tri;
		g.vrot = rotation;
		tri_num[x+y*nx] = num_tri;

		for(int k = 0; k < 14; ++k){
			g.v[k] = v[k];
		}
	}
}

/*!
 * �O���b�h���ƂɃ��b�V������
 */
__global__
void genGridMesh(rxSSGridG* mgrid, int nx, int ny, float zmax, float3 *vrts, int num_vrts, uint* tri_array, uint* tri_num_scan, 
				 float3 *back2_vrts, uint *back2_occ, uint *back2_occ_scan, int yoffset, int voffset)
{
	// �O���b�h�C���f�b�N�X
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

			// �G�b�W�C���f�b�N�X
			int idx;
			if(odd){	// y�����G�b�W
				idx = (x+(edge-1)/2)+y*(nx+1)+yoffset;
			}
			else{		// x�����G�b�W
				idx = x+(y+edge/2)*nx;
			}

			// �ǉ����_�C���f�b�N�X
			int bidx = back2_occ_scan[idx]+voffset;

			if(layer == 2){
				// �Ŕw�ʃG�b�W���_�Ƃ��Ēǉ�
				v[11+k] = bidx;
				g.back_vrts[3+k] = bidx;
			}
			else if(layer == 1){
				// �w�ʃG�b�W���_�Ƃ��Ēǉ�
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
				// �O�ʃG�b�W���_�Ƃ��Ēǉ�
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
	
	// �f�o�b�O�p
	for(int k = 0; k < 14; ++k){
		g.v[k] = v[k];
	}

	int m = g.mesh_num;
	if(m > 0){	// �O���b�h���̃��b�V������0���傫�������烁�b�V������
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
 * �֊s�̕����� : ���ψʒu�̌v�Z
 */
__global__
void smoothSilhouette(float3 *vrts, int num_nvrts, uint* tri_array, int num_tris, float4 *avg_vrts)
{
	// �C���f�b�N�X
	uint blockId = __mul24(blockIdx.y, gridDim.x)+blockIdx.x;
	uint i = __mul24(blockId, blockDim.x)+threadIdx.x;
	
	if(i >= num_tris) return;

	for(int j = 0; j < 3; ++j){
		int idx = tri_array[3*i+j];
		if(idx >= num_nvrts){
			float4 avg = avg_vrts[idx-num_nvrts];

			// ���_���g
			if(avg.w == 0){
				avg += make_float4(vrts[idx], 1.0);
			}

			// �אڒ��_
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
 * �֊s�̕����� : �G�b�W���_�𕽋ψʒu�Ɉړ�
 */
__global__
void smoothSilhouette2(float3 *vrts, int num_vrts, int num_nvrts, float4 *avg_vrts)
{
	// �C���f�b�N�X
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
 * �X�N���[���X�y�[�X���猳��3D��Ԃɖ߂�
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

	// �t���e�ϊ�
	float4 x;
	x.x = dot(IMVQ.e[0], xd);
	x.y = dot(IMVQ.e[1], xd);
	x.z = dot(IMVQ.e[2], xd);
	x.w = dot(IMVQ.e[3], xd);

	vrts[idx] = make_float3(x.x, x.y, x.z);
}


/*!
 * ���_�@���̌v�Z : �ʖ@���̒~��
 */
__global__
void sumFaceNormal(float3 *vrts, int num_nvrts, uint* tri_array, int num_tris, float *nrms)
{
	// �C���f�b�N�X
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
 * ���_�@���̌v�Z : ���_�@���̐��K��
 */
__global__
void normalizeNormal(float3 *nrms, int num_nrms)
{
	// �C���f�b�N�X
	uint blockId = __mul24(blockIdx.y, gridDim.x)+blockIdx.x;
	uint i = __mul24(blockId, blockDim.x)+threadIdx.x;
	
	if(i >= num_nrms) return;

	float3 normal = nrms[i];
	nrms[i] = normalize(normal);
}






#endif // #ifndef _RX_SSM_KERNEL_CU_



