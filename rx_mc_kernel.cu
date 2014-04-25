/*! 
  @file rx_mc_kernel.cu
	
  @brief CUDA�ɂ�郁�b�V������(MC�@)
 
  @author Makoto Fujisawa
  @date 2009-08, 2011-06
*/
// FILE --rx_mc_kernel.cu--

#ifndef _RX_CUMC_KERNEL_CU_
#define _RX_CUMC_KERNEL_CU_


//-----------------------------------------------------------------------------
// �C���N���[�h�t�@�C��
//-----------------------------------------------------------------------------
#include "rx_cu_common.cu"




//-----------------------------------------------------------------------------
// �ϐ�
//-----------------------------------------------------------------------------
// �e�N�X�`���Ɋi�[�����e�[�u��
texture<uint, cudaTextureType1D, cudaReadModeElementType> edgeTex;
texture<uint, cudaTextureType1D, cudaReadModeElementType> triTex;
texture<uint, cudaTextureType1D, cudaReadModeElementType> numVertsTex;

// �e�N�X�`���Ɋi�[�����T���v���{�����[���f�[�^
texture<uchar, cudaTextureType1D, cudaReadModeNormalizedFloat> volumeTex;



//-----------------------------------------------------------------------------
// MARK:MC�@
//-----------------------------------------------------------------------------
/*!
 * �t�B�[���h�֐�(��)
 * @param[in] x,y,z ���W[-1,1]
 * @return �֐��l
 */
__device__
float tangle(float x, float y, float z)
{
	x *= 3.0f;
	y *= 3.0f;
	z *= 3.0f;
	return (x*x*x*x-5.0f*x*x+y*y*y*y-5.0f*y*y+z*z*z*z-5.0f*z*z+11.8f)*0.2f+0.5f;
}

/*!
 * �t�B�[���h�֐�(��)
 * @param[in] p ���W[-1,1]
 * @return �֐��l
 */
__device__
float fieldFunc(float3 p)
{
	return tangle(p.x, p.y, p.z);
}

/*!
 * �t�B�[���h�֐��l�Ɩ@��(��)
 * @param[in] p ���W[-1,1]
 * @return �֐��l
 */
__device__
float4 fieldFunc4(float3 p)
{
	float v = tangle(p.x, p.y, p.z);
	const float d = 0.001f;
	float dx = tangle(p.x+d, p.y, p.z)-v;
	float dy = tangle(p.x, p.y+d, p.z)-v;
	float dz = tangle(p.x, p.y, p.z+d)-v;
	return make_float4(dx, dy, dz, v);
}

/*!
 * �{�����[���f�[�^�̎Q��
 * @param[in] data �{�����[���f�[�^
 * @param[in] p �O���b�h�C���f�b�N�X
 * @param[in] gridSize �O���b�h��
 * @return �{�����[���l
 */
__device__
float sampleVolume(float *data, uint3 p, uint3 gridSize)
{
	p.x = min(p.x, gridSize.x-1);
	p.y = min(p.y, gridSize.y-1);
	p.z = min(p.z, gridSize.z-1);
	uint i = (p.z*gridSize.x*gridSize.y)+(p.y*gridSize.x)+p.x;
//	return tex1Dfetch(volumeTex, i);
	return data[i];
}

/*!
 * �{�����[���f�[�^�̎Q��
 * @param[in] data �{�����[���f�[�^
 * @param[in] p �O���b�h�C���f�b�N�X
 * @param[in] gridSize �O���b�h��
 * @return �{�����[���l
 */
__device__
float sampleVolume2(float *data, uint3 p, uint3 gridSize)
{
	p.x = min(p.x, gridSize.x-1);
	p.y = min(p.y, gridSize.y-1);
	p.z = min(p.z, gridSize.z-1);
	uint i = (p.z*gridSize.x*gridSize.y)+(p.y*gridSize.x)+p.x;
	return data[i];
}

/*!
 * 1D�C���f�b�N�X����3D�C���f�b�N�X�ւ̕ϊ�(�O���b�h����2^n�̂Ƃ��̂ݗL��)
 * @param[in] i 1D�C���f�b�N�X
 * @param[in] gridSizeShift (0,n,n*n)
 * @param[in] gridSizeMask  �}�X�N
 * @return 3D�C���f�b�N�X
__device__
uint3 calcGridIdx(uint i, uint3 gridSizeShift, uint3 gridSizeMask)
{
	uint3 gridPos;
	gridPos.x = i & gridSizeMask.x;
	gridPos.y = (i >> gridSizeShift.y) & gridSizeMask.y;
	gridPos.z = (i >> gridSizeShift.z) & gridSizeMask.z;
	return gridPos;
}
 */
/*!
 * 1D�C���f�b�N�X����3D�C���f�b�N�X�ւ̕ϊ�(�O���b�h���͔C��)
 * @param[in] i 1D�C���f�b�N�X
 * @param[in] gridSize �O���b�h��
 * @return 3D�C���f�b�N�X
 */
__device__
uint3 calcGridIdxU(uint i, uint3 ngrid)
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
uint calcGridIdx3(uint3 p, uint3 ngrid)
{
	p.x = min(p.x, ngrid.x-1);
	p.y = min(p.y, ngrid.y-1);
	p.z = min(p.z, ngrid.z-1);
	return (p.z*ngrid.x*ngrid.y)+(p.y*ngrid.x)+p.x;
}



/*!
 * �T���v���{�����[���̐��^��Ԃɂ��C�ӈʒu�ł̒l�擾
 * @param[in] data �{�����[���f�[�^
 * @param[in] pos �l�擾���W
 * @param[in] gridSize �{�����[���O���b�h��
 * @param[in] voxelSize �{�����[���O���b�h��
 * @return �l
 */
__device__
float linearInterp(float *data, float3 pos, uint3 gridSize, float3 voxelSize)
{
	pos.x = (pos.x+1.0f)/voxelSize.x;
	pos.y = (pos.y+1.0f)/voxelSize.y;
	pos.z = (pos.z+1.0f)/voxelSize.z;

	uint3 grid;
	grid.x = (uint)(pos.x);
	grid.y = (uint)(pos.y);
	grid.z = (uint)(pos.z);

	float f[8];
	f[0] = sampleVolume(data, grid, gridSize);
	f[1] = sampleVolume(data, grid+make_uint3(1, 0, 0), gridSize);
	f[2] = sampleVolume(data, grid+make_uint3(1, 1, 0), gridSize);
	f[3] = sampleVolume(data, grid+make_uint3(0, 1, 0), gridSize);
	f[4] = sampleVolume(data, grid+make_uint3(0, 0, 1), gridSize);
	f[5] = sampleVolume(data, grid+make_uint3(1, 0, 1), gridSize);
	f[6] = sampleVolume(data, grid+make_uint3(1, 1, 1), gridSize);
	f[7] = sampleVolume(data, grid+make_uint3(0, 1, 1), gridSize);

	float s0 = pos.x-(float)grid.x;
	float t0 = pos.y-(float)grid.y;
	float u0 = pos.z-(float)grid.z;

	float s1 = 1.0-s0;
	float t1 = 1.0-t0;
	float u1 = 1.0-u0;

	float fz0, fz1;

	fz0 = s1*(t1*f[0]+t0*f[3])+s0*(t1*f[1]+t0*f[2]);
	fz1 = s1*(t1*f[4]+t0*f[7])+s0*(t1*f[5]+t0*f[6]);

	return u1*fz0+u0*fz1;
}

/*!
 * ���`���
 * @param[in] isolavel 臒l
 * @param[in] p0,p1 ���[�_�l
 * @param[in] f0,f1 ���[�_�̉A�֐��l
 * @return ��Ԓl
 */
__device__
float linearInterp1(float isolevel, float p0, float p1, float f0, float f1)
{
	float t = (isolevel-f0)/(f1-f0);
	return lerp(p0, p1, t);
} 

/*!
 * �G�b�W�ɉ��������`��Ԃɂ��A�֐��l0�ƂȂ钸�_�ʒu�̌v�Z
 * @param[in] isolavel 臒l
 * @param[in] p0,p1 ���[�_���W
 * @param[in] f0,f1 ���[�_�̉A�֐��l
 * @return ���_���W
 */
__device__
float3 vertexInterp(float isolevel, float3 p0, float3 p1, float f0, float f1)
{
	float t = (isolevel-f0)/(f1-f0);
	return lerp(p0, p1, t);
} 

/*!
 * �G�b�W�ɉ��������`��Ԃɂ��A�֐��l0�ƂȂ钸�_�ʒu�Ɩ@���̌v�Z
 * @param[in] isolavel 臒l
 * @param[in] p0,p1 ���[�_���W
 * @param[in] f0,f1 ���[�_�̉A�֐��l
 * @param[out] p ���_���W
 * @param[out] n ���_�@��
 */
__device__
void vertexInterp2(float isolevel, float3 p0, float3 p1, float4 f0, float4 f1, float3 &p, float3 &n)
{
	float t = (isolevel-f0.w)/(f1.w-f0.w);
	p = lerp(p0, p1, t);
	n.x = lerp(f0.x, f1.x, t);
	n.y = lerp(f0.y, f1.y, t);
	n.z = lerp(f0.z, f1.z, t);
//	n = normalize(n);
} 

/*!
 * �O�p�`���b�V���̖@���v�Z(�O�ώg�p)
 * @param[in] v0,v1,v2 ���b�V����3���_
 * @return �@��(�P�ʖ@���ł͂Ȃ�)
 */
__device__
float3 calcNormal(float3 *v0, float3 *v1, float3 *v2)
{
	float3 edge0 = *v1-*v0;
	float3 edge1 = *v2-*v0;
	return cross(edge0, edge1);
}

/*!
 * �{�����[���f�[�^����̖@���v�Z
 * @param[in] v0,v1,v2 ���b�V����3���_
 * @return �@��(�P�ʖ@���ł͂Ȃ�)
 */
__device__
float3 calcNormalF(float *data, float3 pos, uint3 gridSize, float3 voxelSize)
{
	float3 nrm;
	float3 p, m;

	float3 d = 0.2*voxelSize;

	p.x = linearInterp(data, pos+make_float3(d.x, 0, 0), gridSize, voxelSize);
	p.y = linearInterp(data, pos+make_float3(0, d.y, 0), gridSize, voxelSize);
	p.z = linearInterp(data, pos+make_float3(0, 0, d.z), gridSize, voxelSize);
	m.x = linearInterp(data, pos-make_float3(d.x, 0, 0), gridSize, voxelSize);
	m.y = linearInterp(data, pos-make_float3(0, d.y, 0), gridSize, voxelSize);
	m.z = linearInterp(data, pos-make_float3(0, 0, d.z), gridSize, voxelSize);

	nrm.x = -(p.x-m.x)/(2.0*d.x);
	nrm.y = -(p.y-m.y)/(2.0*d.y);
	nrm.z = -(p.z-m.z)/(2.0*d.z);

	float l = length(nrm);
	if(l > 0){
		nrm /= l;
	}

	return nrm;
}



//-----------------------------------------------------------------------------
// MC�@�J�[�l��
//-----------------------------------------------------------------------------
/*!
 * �{�N�Z�����Ƃɓ��l���_�ʒu(voxelCubeIdx),���l���_��(voxelVerts),���b�V����(voxelTris)���v�Z
 * @param[out] voxelCubeIdx �{�N�Z���̓��l���_�ʒu(8bit�}�X�N)
 * @param[out] voxelVerts �{�N�Z���̓��l���_��
 * @param[out] voxelTris �{�N�Z���̃��b�V����
 * @param[in] volume �{�����[���f�[�^
 * @param[in] gridSize �O���b�h��
 * @param[in] gridSizeShift,gridSizeMask �V�t�g�ƃ}�X�N
 * @param[in] numVoxels ���{�N�Z����
 * @param[in] voxelSize �{�N�Z����
 * @param[in] isoValue 臒l
 */
__global__ 
void ClassifyVoxel2(uint* voxelCubeIdx, uint *voxelVerts, uint *voxelTris, uint *voxelOccupied, float *volume,
					uint3 ngrids, uint nvoxels, float3 voxel_h, float threshold)
{
	uint blockId = __mul24(blockIdx.y, gridDim.x)+blockIdx.x;
	uint i = __mul24(blockId, blockDim.x)+threadIdx.x;

	if(i < nvoxels){
		uint3 gridPos = calcGridIdxU(i, ngrids);

		// �O���b�h8���_�̉A�֐��l���Q�Ƃ���
#if SAMPLE_VOLUME
		float field[8];
		field[0] = sampleVolume2(volume, gridPos, ngrids);
		field[1] = sampleVolume2(volume, gridPos+make_uint3(1, 0, 0), ngrids);
		field[2] = sampleVolume2(volume, gridPos+make_uint3(1, 1, 0), ngrids);
		field[3] = sampleVolume2(volume, gridPos+make_uint3(0, 1, 0), ngrids);
		field[4] = sampleVolume2(volume, gridPos+make_uint3(0, 0, 1), ngrids);
		field[5] = sampleVolume2(volume, gridPos+make_uint3(1, 0, 1), ngrids);
		field[6] = sampleVolume2(volume, gridPos+make_uint3(1, 1, 1), ngrids);
		field[7] = sampleVolume2(volume, gridPos+make_uint3(0, 1, 1), ngrids);
#else
		float3 p;
		p.x = -1.0f+(gridPos.x*voxel_h.x);
		p.y = -1.0f+(gridPos.y*voxel_h.y);
		p.z = -1.0f+(gridPos.z*voxel_h.z);

		float field[8];
		field[0] = fieldFunc(p);
		field[1] = fieldFunc(p+make_float3(voxel_h.x, 0, 0));
		field[2] = fieldFunc(p+make_float3(voxel_h.x, voxel_h.y, 0));
		field[3] = fieldFunc(p+make_float3(0, voxel_h.y, 0));
		field[4] = fieldFunc(p+make_float3(0, 0, voxel_h.z));
		field[5] = fieldFunc(p+make_float3(voxel_h.x, 0, voxel_h.z));
		field[6] = fieldFunc(p+make_float3(voxel_h.x, voxel_h.y, voxel_h.z));
		field[7] = fieldFunc(p+make_float3(0, voxel_h.y, voxel_h.z));
#endif

		// �O���b�h���̒��_���e�[�u���p�̃C���f�b�N�X�쐬
		uint cubeindex;
		cubeindex =  uint(field[0] < threshold); 
		cubeindex += uint(field[1] < threshold)*2; 
		cubeindex += uint(field[2] < threshold)*4; 
		cubeindex += uint(field[3] < threshold)*8; 
		cubeindex += uint(field[4] < threshold)*16; 
		cubeindex += uint(field[5] < threshold)*32; 
		cubeindex += uint(field[6] < threshold)*64; 
		cubeindex += uint(field[7] < threshold)*128;


		uint numVerts = tex1Dfetch(numVertsTex, cubeindex);

		voxelCubeIdx[i] = cubeindex;	// ��̌v�Z�̂��߂�cubeindex���L�^���Ă���
		voxelVerts[i] = numVerts;		// �O���b�h���̒��_��
		voxelTris[i] = numVerts/3;		// �O���b�h���̎O�p�`�|���S����

#if SKIP_EMPTY_VOXELS
		voxelOccupied[i] = (numVerts > 0);	// �O���b�h���Ƀ��b�V�����܂ނ��ǂ���
#endif
	}

}

/*!
 * �G�b�W���Ƃ�0���l�ʒ��_��T��
 * @param[out] edgeVrts ���_��
 * @param[out] edgeOccupied �G�b�W�̒��_��L�z��
 * @param[in] volume �{�����[���f�[�^
 * @param[in] dir �T������
 * @param[in] gridSize1 �G�b�W��(�O���b�h��+1)
 * @param[in] gridSize �O���b�h��
 * @param[in] gridSizeShift,gridSizeMask �V�t�g�ƃ}�X�N
 * @param[in] numVoxels ���{�N�Z����
 * @param[in] voxelSize �{�N�Z����
 * @param[in] isoValue 臒l
 */
__global__ 
void CalVertexEdge(float4* edgeVrts, uint *edgeOccupied, float *volume, uint3 dir, 
				   uint3 edgeSize, uint3 ngrids, uint nvoxels, uint nedge, 
				   float3 voxel_h, float3 voxel_min, float threshold)
{
	uint blockId = __mul24(blockIdx.y, gridDim.x)+blockIdx.x;
	uint i = __mul24(blockId, blockDim.x)+threadIdx.x;
	if(i >= nedge) return;

	uint3 gridPos = calcGridIdxU(i, edgeSize);

	float3 p;
	p.x = voxel_min.x+(gridPos.x*voxel_h.x);
	p.y = voxel_min.y+(gridPos.y*voxel_h.y);
	p.z = voxel_min.z+(gridPos.z*voxel_h.z);

	// calculate cell vertex positions
	float3 v[2];
	v[0] = p;
	v[1] = p+make_float3(dir.x*voxel_h.x, dir.y*voxel_h.y, dir.z*voxel_h.z);

	// read field values at neighbouring grid vertices
#if SAMPLE_VOLUME
	float field[2];
	field[0] = sampleVolume2(volume, gridPos,	  ngrids);
	field[1] = sampleVolume2(volume, gridPos+dir, ngrids);

	uint cubeindex;
	cubeindex =  uint(field[0] < threshold); 
	cubeindex += uint(field[1] < threshold)*2; 
#else
	// evaluate field values
	float4 field[2];
	field[0] = fieldFunc4(p);
	field[1] = fieldFunc4(p+make_float3(dir.x*voxelSize.x, dir.y*voxelSize.y, dir.z*voxelSize.z));

	uint cubeindex;
	cubeindex =  uint(field[0].w < isoValue); 
	cubeindex += uint(field[1].w < isoValue)*2; 
#endif


	if(cubeindex == 1 || cubeindex == 2){
		float3 vertex, normal = make_float3(0.0f, 5.0f, 0.0f);

#if SAMPLE_VOLUME
		vertex = vertexInterp(threshold, v[0], v[1], field[0], field[1]);
#else
		vertexInterp2(isoValue, v[0], v[1], field[0], field[1], vertex, normal);
#endif

		edgeVrts[i] = make_float4(vertex, 1.0f);
		edgeOccupied[i] = 1;
	}
	else{
		//edgeVrts[i] = make_float4(0.0f, 5.0f, 0.0f, 1.0f);
		edgeOccupied[i] = 0;
	}
}

/*!
 * �{�N�Z�����Ƃ�Scan���ʂɊ�Â��ăG�b�W��̒��_�����l�߂� 
 * @param[out] compactedVrts ��G�b�W���l�߂����_���W�z��
 * @param[in] occupied �G�b�W��L���(���_�����݂���:1, ���Ȃ�:0)
 * @param[in] occupiedScan �G�b�W��L��񂩂�쐬����Prefix Sum(Scan)
 * @param[in] vrts �e�G�b�W�̒��_���W(���_�����݂��Ȃ��ꍇ�͂��ׂĂ̗v�f��0)
 * @param[in] num ���G�b�W��(3x(�{�N�Z����+1))
 */
__global__ 
void CompactVoxels(uint *compacted, uint *occupied, uint *occupiedScan, uint num)
{
	// �C���f�b�N�X
	uint blockId = __mul24(blockIdx.y, gridDim.x)+blockIdx.x;
	uint i = __mul24(blockId, blockDim.x)+threadIdx.x;

	if(occupied[i] && (i < num)) {
		compacted[occupiedScan[i]] = i;
	}
}

/*!
 * �{�N�Z�����Ƃ�Scan���ʂɊ�Â��ăG�b�W��̒��_�����l�߂� 
 * @param[out] compactedVrts ��G�b�W���l�߂����_���W�z��
 * @param[in] occupied �G�b�W��L���(���_�����݂���:1, ���Ȃ�:0)
 * @param[in] occupiedScan �G�b�W��L��񂩂�쐬����Prefix Sum(Scan)
 * @param[in] vrts �e�G�b�W�̒��_���W(���_�����݂��Ȃ��ꍇ�͂��ׂĂ̗v�f��0)
 * @param[in] num ���G�b�W��(3x(�{�N�Z����+1))
 */
__global__ 
void CompactEdges(float4 *compactedVrts, uint *occupied, uint *occupiedScan, float4 *vrts, uint num)
{
	// �C���f�b�N�X
	uint blockId = __mul24(blockIdx.y, gridDim.x)+blockIdx.x;
	uint i = __mul24(blockId, blockDim.x)+threadIdx.x;

	if(occupied[i] && (i < num)) {
		compactedVrts[occupiedScan[i]] = vrts[i];
	}
}



/*!
 * �{�N�Z�����ƂɃ��b�V�����_�C���f�b�N�X���쐬
 * @param[out] voxelIdx ���b�V�����_�C���f�b�N�X
 * @param[in] voxelTrisScan �{�N�Z�����Ƃ̃��b�V����Scan
 * @param[in] edgeOccupiedScan �G�b�W���Ƃ̒��_����Scan
 * @param[in] voxelCubeIdx �{�N�Z�����Ƃ̓��l�ʒ��_���݃}�X�N
 * @param[in] gridSize1 �G�b�W��(�O���b�h��+1)
 * @param[in] gridSize �O���b�h��
 * @param[in] gridSizeShift,gridSizeMask �V�t�g�ƃ}�X�N
 * @param[in] numVoxels ���{�N�Z����
 * @param[in] voxelSize �{�N�Z����
 * @param[in] isoValue 臒l
 */
__global__ 
void GenerateTriangles3(uint3 *vertIdx, uint *voxelTrisScan, uint *edgeOccupiedScan, uint *voxelCubeIdx, 
						uint3 edgeSizeX, uint3 edgeSizeY, uint3 edgeSizeZ, uint3 edgeNum, 
						uint *compactedVoxelArray, uint3 ngrids, uint activeVoxels,
						float3 voxelSize, float isoValue, uint maxVerts, uint numMesh)
{
	uint blockId = __mul24(blockIdx.y, gridDim.x)+blockIdx.x;
	uint idx = __mul24(blockId, blockDim.x)+threadIdx.x;

	if(idx > activeVoxels-1){
		return;
		//idx = activeVoxels-1;
	}

#if SKIP_EMPTY_VOXELS
	uint voxel = compactedVoxelArray[idx];
#else
	uint voxel = idx;
#endif

	uint3 gpos = calcGridIdxU(voxel, ngrids);

	uint cubeindex = voxelCubeIdx[voxel];

#if USE_SHARED
	__shared__ uint vertlist[12*NTHREADS];
	vertlist[12*threadIdx.x+0]  = edgeOccupiedScan[calcGridIdx3(make_uint3(gpos.x,   gpos.y,   gpos.z),   edgeSizeX)];
	vertlist[12*threadIdx.x+1]  = edgeOccupiedScan[calcGridIdx3(make_uint3(gpos.x+1, gpos.y,   gpos.z),   edgeSizeY)+edgeNum.x];
	vertlist[12*threadIdx.x+2]  = edgeOccupiedScan[calcGridIdx3(make_uint3(gpos.x,   gpos.y+1, gpos.z),   edgeSizeX)];
	vertlist[12*threadIdx.x+3]  = edgeOccupiedScan[calcGridIdx3(make_uint3(gpos.x,   gpos.y,   gpos.z),   edgeSizeY)+edgeNum.x];
	vertlist[12*threadIdx.x+4]  = edgeOccupiedScan[calcGridIdx3(make_uint3(gpos.x,   gpos.y,   gpos.z+1), edgeSizeX)];
	vertlist[12*threadIdx.x+5]  = edgeOccupiedScan[calcGridIdx3(make_uint3(gpos.x+1, gpos.y,   gpos.z+1), edgeSizeY)+edgeNum.x];
	vertlist[12*threadIdx.x+6]  = edgeOccupiedScan[calcGridIdx3(make_uint3(gpos.x,   gpos.y+1, gpos.z+1), edgeSizeX)];
	vertlist[12*threadIdx.x+7]  = edgeOccupiedScan[calcGridIdx3(make_uint3(gpos.x,   gpos.y,   gpos.z+1), edgeSizeY)+edgeNum.x];
	vertlist[12*threadIdx.x+8]  = edgeOccupiedScan[calcGridIdx3(make_uint3(gpos.x,   gpos.y,   gpos.z), edgeSizeZ)+edgeNum.x+edgeNum.y];
	vertlist[12*threadIdx.x+9]  = edgeOccupiedScan[calcGridIdx3(make_uint3(gpos.x+1, gpos.y,   gpos.z), edgeSizeZ)+edgeNum.x+edgeNum.y];
	vertlist[12*threadIdx.x+10] = edgeOccupiedScan[calcGridIdx3(make_uint3(gpos.x+1, gpos.y+1, gpos.z), edgeSizeZ)+edgeNum.x+edgeNum.y];
	vertlist[12*threadIdx.x+11] = edgeOccupiedScan[calcGridIdx3(make_uint3(gpos.x,   gpos.y+1, gpos.z), edgeSizeZ)+edgeNum.x+edgeNum.y];
	__syncthreads();
#else
	uint vertlist[12];
	vertlist[0]  = edgeOccupiedScan[calcGridIdx3(make_uint3(gpos.x,   gpos.y,   gpos.z),   edgeSizeX)];
	vertlist[2]  = edgeOccupiedScan[calcGridIdx3(make_uint3(gpos.x,   gpos.y+1, gpos.z),   edgeSizeX)];
	vertlist[4]  = edgeOccupiedScan[calcGridIdx3(make_uint3(gpos.x,   gpos.y,   gpos.z+1), edgeSizeX)];
	vertlist[6]  = edgeOccupiedScan[calcGridIdx3(make_uint3(gpos.x,   gpos.y+1, gpos.z+1), edgeSizeX)];

	vertlist[1]  = edgeOccupiedScan[calcGridIdx3(make_uint3(gpos.x+1, gpos.y,   gpos.z),   edgeSizeY)+edgeNum.x];
	vertlist[3]  = edgeOccupiedScan[calcGridIdx3(make_uint3(gpos.x,   gpos.y,   gpos.z),   edgeSizeY)+edgeNum.x];
	vertlist[5]  = edgeOccupiedScan[calcGridIdx3(make_uint3(gpos.x+1, gpos.y,   gpos.z+1), edgeSizeY)+edgeNum.x];
	vertlist[7]  = edgeOccupiedScan[calcGridIdx3(make_uint3(gpos.x,   gpos.y,   gpos.z+1), edgeSizeY)+edgeNum.x];

	vertlist[8]  = edgeOccupiedScan[calcGridIdx3(make_uint3(gpos.x,   gpos.y,   gpos.z), edgeSizeZ)+edgeNum.x+edgeNum.y];
	vertlist[9]  = edgeOccupiedScan[calcGridIdx3(make_uint3(gpos.x+1, gpos.y,   gpos.z), edgeSizeZ)+edgeNum.x+edgeNum.y];
	vertlist[10] = edgeOccupiedScan[calcGridIdx3(make_uint3(gpos.x+1, gpos.y+1, gpos.z), edgeSizeZ)+edgeNum.x+edgeNum.y];
	vertlist[11] = edgeOccupiedScan[calcGridIdx3(make_uint3(gpos.x,   gpos.y+1, gpos.z), edgeSizeZ)+edgeNum.x+edgeNum.y];
#endif

	// output triangle
	uint numTri = tex1Dfetch(numVertsTex, cubeindex)/3;

	for(int i = 0; i < numTri; ++i){
		uint index = voxelTrisScan[voxel]+i;

		uint vidx[3];
		uint edge[3];
		edge[0] = tex1Dfetch(triTex, (cubeindex*16)+3*i);
		edge[1] = tex1Dfetch(triTex, (cubeindex*16)+3*i+1);
		edge[2] = tex1Dfetch(triTex, (cubeindex*16)+3*i+2);

#if USE_SHARED
		vidx[0] = min(vertlist[12*threadIdx.x+edge[0]], maxVerts-1);
		vidx[2] = min(vertlist[12*threadIdx.x+edge[1]], maxVerts-1);
		vidx[1] = min(vertlist[12*threadIdx.x+edge[2]], maxVerts-1);
#else
		vidx[0] = min(vertlist[edge[0]], maxVerts-1);
		vidx[2] = min(vertlist[edge[1]], maxVerts-1);
		vidx[1] = min(vertlist[edge[2]], maxVerts-1);
#endif
		if(index < numMesh){
			vertIdx[index] = make_uint3(vidx[2], vidx[1], vidx[0]);
		}
	}
}

__global__ 
void CalVertexNormalKernel(float4 *vrts, uint3 *tris, float3 *nrms, int nvrts, int ntris)
{
	uint blockId = __mul24(blockIdx.y, gridDim.x)+blockIdx.x;
	uint idx = __mul24(blockId, blockDim.x)+threadIdx.x;
	if(idx >= ntris) return;

	uint3 id = tris[idx];

	float3 normal;
	normal = cross(make_float3(vrts[id.y])-make_float3(vrts[id.x]), 
				   make_float3(vrts[id.z])-make_float3(vrts[id.x]));

	//normal = normalize(normal);

#ifdef RX_USE_ATOMIC_FUNC
	atomicFloatAdd(&nrms[id.x].x, normal.x);
	atomicFloatAdd(&nrms[id.x].y, normal.y);
	atomicFloatAdd(&nrms[id.x].z, normal.z);
	atomicFloatAdd(&nrms[id.y].x, normal.x);
	atomicFloatAdd(&nrms[id.y].y, normal.y);
	atomicFloatAdd(&nrms[id.y].z, normal.z);
	atomicFloatAdd(&nrms[id.z].x, normal.x);
	atomicFloatAdd(&nrms[id.z].y, normal.y);
	atomicFloatAdd(&nrms[id.z].z, normal.z);
#else
	nrms[id.x] += normal;
	__threadfence();
	nrms[id.y] += normal;
	__threadfence();
	nrms[id.z] += normal;
#endif
}

__global__ 
void NormalizeKernel(float3 *v, int nv)
{
	uint blockId = __mul24(blockIdx.y, gridDim.x)+blockIdx.x;
	uint idx = __mul24(blockId, blockDim.x)+threadIdx.x;
	if(idx >= nv) return;

	v[idx] = normalize(v[idx]);
}






//-----------------------------------------------------------------------------
// ���_�E�@���\���Ń��b�V���쐬
//-----------------------------------------------------------------------------

// classify voxel based on number of vertices it will generate
// one thread per voxel
__global__ void
classifyVoxel(uint* voxelVerts, uint *voxelOccupied, float *volume,
			  uint3 ngrids, uint nvoxels, float3 voxel_h, float threshold)
{
	uint blockId = __mul24(blockIdx.y, gridDim.x) + blockIdx.x;
	uint i = __mul24(blockId, blockDim.x) + threadIdx.x;

	uint3 gridPos = calcGridIdxU(i, ngrids);

#if SAMPLE_VOLUME
	float field[8];
	field[0] = sampleVolume(volume, gridPos, ngrids);
	field[1] = sampleVolume(volume, gridPos + make_uint3(1, 0, 0), ngrids);
	field[2] = sampleVolume(volume, gridPos + make_uint3(1, 1, 0), ngrids);
	field[3] = sampleVolume(volume, gridPos + make_uint3(0, 1, 0), ngrids);
	field[4] = sampleVolume(volume, gridPos + make_uint3(0, 0, 1), ngrids);
	field[5] = sampleVolume(volume, gridPos + make_uint3(1, 0, 1), ngrids);
	field[6] = sampleVolume(volume, gridPos + make_uint3(1, 1, 1), ngrids);
	field[7] = sampleVolume(volume, gridPos + make_uint3(0, 1, 1), ngrids);
#else
	float3 p;
	p.x = -1.0f + (gridPos.x * voxel_h.x);
	p.y = -1.0f + (gridPos.y * voxel_h.y);
	p.z = -1.0f + (gridPos.z * voxel_h.z);

	float field[8];
	field[0] = fieldFunc(p);
	field[1] = fieldFunc(p + make_float3(voxel_h.x, 0, 0));
	field[2] = fieldFunc(p + make_float3(voxel_h.x, voxel_h.y, 0));
	field[3] = fieldFunc(p + make_float3(0, voxel_h.y, 0));
	field[4] = fieldFunc(p + make_float3(0, 0, voxel_h.z));
	field[5] = fieldFunc(p + make_float3(voxel_h.x, 0, voxel_h.z));
	field[6] = fieldFunc(p + make_float3(voxel_h.x, voxel_h.y, voxel_h.z));
	field[7] = fieldFunc(p + make_float3(0, voxel_h.y, voxel_h.z));
#endif

	// calculate flag indicating if each vertex is inside or outside isosurface
	uint cubeindex;
	cubeindex =  uint(field[0] < threshold); 
	cubeindex += uint(field[1] < threshold)*2; 
	cubeindex += uint(field[2] < threshold)*4; 
	cubeindex += uint(field[3] < threshold)*8; 
	cubeindex += uint(field[4] < threshold)*16; 
	cubeindex += uint(field[5] < threshold)*32; 
	cubeindex += uint(field[6] < threshold)*64; 
	cubeindex += uint(field[7] < threshold)*128;

	// read number of vertices from texture
	uint numVerts = tex1Dfetch(numVertsTex, cubeindex);

	if (i < nvoxels) {
		voxelVerts[i] = numVerts;
		voxelOccupied[i] = (numVerts > 0);
	}
}


// compact voxel array
__global__ void
compactVoxels(uint *compactedVoxelArray, uint *voxelOccupied, uint *voxelOccupiedScan, uint numVoxels)
{
	uint blockId = __mul24(blockIdx.y, gridDim.x) + blockIdx.x;
	uint i = __mul24(blockId, blockDim.x) + threadIdx.x;

	if (voxelOccupied[i] && (i < numVoxels)) {
		compactedVoxelArray[ voxelOccupiedScan[i] ] = i;
	}
}


// generate triangles for each voxel using marching cubes
// interpolates normals from field function
__global__ void
generateTriangles(float4 *pos, float4 *norm, uint *compactedVoxelArray, uint *numVertsScanned,
				  uint3 ngrids, float3 voxel_h, float3 voxel_min, float isoValue, uint activeVoxels, uint maxVerts)
{
	uint blockId = __mul24(blockIdx.y, gridDim.x) + blockIdx.x;
	uint i = __mul24(blockId, blockDim.x) + threadIdx.x;

	if (i > activeVoxels - 1) {
		// can't return here because of syncthreads()
		i = activeVoxels - 1;
	}

#if SKIP_EMPTY_VOXELS
	uint voxel = compactedVoxelArray[i];
#else
	uint voxel = i;
#endif

	// compute position in 3d grid
	uint3 gridPos = calcGridIdxU(voxel, ngrids);

	float3 p;
	p.x = voxel_min.x+(gridPos.x*voxel_h.x);
	p.y = voxel_min.y+(gridPos.y*voxel_h.y);
	p.z = voxel_min.z+(gridPos.z*voxel_h.z);
	//p.x = -1.0f + (gridPos.x * voxelSize.x);
	//p.y = -1.0f + (gridPos.y * voxelSize.y);
	//p.z = -1.0f + (gridPos.z * voxelSize.z);

	// calculate cell vertex positions
	float3 v[8];
	v[0] = p;
	v[1] = p + make_float3(voxel_h.x,		 0,		 0);
	v[2] = p + make_float3(voxel_h.x, voxel_h.y,		 0);
	v[3] = p + make_float3(		0, voxel_h.y,		 0);
	v[4] = p + make_float3(		0,		 0, voxel_h.z);
	v[5] = p + make_float3(voxel_h.x,		 0, voxel_h.z);
	v[6] = p + make_float3(voxel_h.x, voxel_h.y, voxel_h.z);
	v[7] = p + make_float3(		0, voxel_h.y, voxel_h.z);

	// evaluate field values
	float4 field[8];
	field[0] = fieldFunc4(v[0]);
	field[1] = fieldFunc4(v[1]);
	field[2] = fieldFunc4(v[2]);
	field[3] = fieldFunc4(v[3]);
	field[4] = fieldFunc4(v[4]);
	field[5] = fieldFunc4(v[5]);
	field[6] = fieldFunc4(v[6]);
	field[7] = fieldFunc4(v[7]);

	// recalculate flag
	// (this is faster than storing it in global memory)
	uint cubeindex;
	cubeindex =  uint(field[0].w < isoValue); 
	cubeindex += uint(field[1].w < isoValue)*2; 
	cubeindex += uint(field[2].w < isoValue)*4; 
	cubeindex += uint(field[3].w < isoValue)*8; 
	cubeindex += uint(field[4].w < isoValue)*16; 
	cubeindex += uint(field[5].w < isoValue)*32; 
	cubeindex += uint(field[6].w < isoValue)*64; 
	cubeindex += uint(field[7].w < isoValue)*128;

	// find the vertices where the surface intersects the cube 

#if USE_SHARED
	// use partioned shared memory to avoid using local memory
	__shared__ float3 vertlist[12*NTHREADS];
	__shared__ float3 normlist[12*NTHREADS];

	//vertexInterp2(isoValue, v[0], v[1], field[0], field[1], vertlist[threadIdx.x], normlist[threadIdx.x]);
	//vertexInterp2(isoValue, v[1], v[2], field[1], field[2], vertlist[threadIdx.x+NTHREADS], normlist[threadIdx.x+NTHREADS]);
	//vertexInterp2(isoValue, v[2], v[3], field[2], field[3], vertlist[threadIdx.x+(NTHREADS*2)], normlist[threadIdx.x+(NTHREADS*2)]);
	//vertexInterp2(isoValue, v[3], v[0], field[3], field[0], vertlist[threadIdx.x+(NTHREADS*3)], normlist[threadIdx.x+(NTHREADS*3)]);
	//vertexInterp2(isoValue, v[4], v[5], field[4], field[5], vertlist[threadIdx.x+(NTHREADS*4)], normlist[threadIdx.x+(NTHREADS*4)]);
	//vertexInterp2(isoValue, v[5], v[6], field[5], field[6], vertlist[threadIdx.x+(NTHREADS*5)], normlist[threadIdx.x+(NTHREADS*5)]);
	//vertexInterp2(isoValue, v[6], v[7], field[6], field[7], vertlist[threadIdx.x+(NTHREADS*6)], normlist[threadIdx.x+(NTHREADS*6)]);
	//vertexInterp2(isoValue, v[7], v[4], field[7], field[4], vertlist[threadIdx.x+(NTHREADS*7)], normlist[threadIdx.x+(NTHREADS*7)]);
	//vertexInterp2(isoValue, v[0], v[4], field[0], field[4], vertlist[threadIdx.x+(NTHREADS*8)], normlist[threadIdx.x+(NTHREADS*8)]);
	//vertexInterp2(isoValue, v[1], v[5], field[1], field[5], vertlist[threadIdx.x+(NTHREADS*9)], normlist[threadIdx.x+(NTHREADS*9)]);
	//vertexInterp2(isoValue, v[2], v[6], field[2], field[6], vertlist[threadIdx.x+(NTHREADS*10)], normlist[threadIdx.x+(NTHREADS*10)]);
	//vertexInterp2(isoValue, v[3], v[7], field[3], field[7], vertlist[threadIdx.x+(NTHREADS*11)], normlist[threadIdx.x+(NTHREADS*11)]);
	vertexInterp2(isoValue, v[0], v[1], field[0], field[1], vertlist[12*threadIdx.x+0], normlist[12*threadIdx.x+0]);
	vertexInterp2(isoValue, v[1], v[2], field[1], field[2], vertlist[12*threadIdx.x+1], normlist[12*threadIdx.x+1]);
	vertexInterp2(isoValue, v[2], v[3], field[2], field[3], vertlist[12*threadIdx.x+2], normlist[12*threadIdx.x+2]);
	vertexInterp2(isoValue, v[3], v[0], field[3], field[0], vertlist[12*threadIdx.x+3], normlist[12*threadIdx.x+3]);
	vertexInterp2(isoValue, v[4], v[5], field[4], field[5], vertlist[12*threadIdx.x+4], normlist[12*threadIdx.x+4]);
	vertexInterp2(isoValue, v[5], v[6], field[5], field[6], vertlist[12*threadIdx.x+5], normlist[12*threadIdx.x+5]);
	vertexInterp2(isoValue, v[6], v[7], field[6], field[7], vertlist[12*threadIdx.x+6], normlist[12*threadIdx.x+6]);
	vertexInterp2(isoValue, v[7], v[4], field[7], field[4], vertlist[12*threadIdx.x+7], normlist[12*threadIdx.x+7]);
	vertexInterp2(isoValue, v[0], v[4], field[0], field[4], vertlist[12*threadIdx.x+8], normlist[12*threadIdx.x+8]);
	vertexInterp2(isoValue, v[1], v[5], field[1], field[5], vertlist[12*threadIdx.x+9], normlist[12*threadIdx.x+9]);
	vertexInterp2(isoValue, v[2], v[6], field[2], field[6], vertlist[12*threadIdx.x+10], normlist[12*threadIdx.x+10]);
	vertexInterp2(isoValue, v[3], v[7], field[3], field[7], vertlist[12*threadIdx.x+11], normlist[12*threadIdx.x+11]);
	__syncthreads();

#else
	float3 vertlist[12];
	float3 normlist[12];

	vertexInterp2(isoValue, v[0], v[1], field[0], field[1], vertlist[0], normlist[0]);
	vertexInterp2(isoValue, v[1], v[2], field[1], field[2], vertlist[1], normlist[1]);	
	vertexInterp2(isoValue, v[2], v[3], field[2], field[3], vertlist[2], normlist[2]);
	vertexInterp2(isoValue, v[3], v[0], field[3], field[0], vertlist[3], normlist[3]);

	vertexInterp2(isoValue, v[4], v[5], field[4], field[5], vertlist[4], normlist[4]); 
	vertexInterp2(isoValue, v[5], v[6], field[5], field[6], vertlist[5], normlist[5]);
	vertexInterp2(isoValue, v[6], v[7], field[6], field[7], vertlist[6], normlist[6]);
	vertexInterp2(isoValue, v[7], v[4], field[7], field[4], vertlist[7], normlist[7]);

	vertexInterp2(isoValue, v[0], v[4], field[0], field[4], vertlist[8], normlist[8]); 
	vertexInterp2(isoValue, v[1], v[5], field[1], field[5], vertlist[9], normlist[9]);
	vertexInterp2(isoValue, v[2], v[6], field[2], field[6], vertlist[10], normlist[10]);
	vertexInterp2(isoValue, v[3], v[7], field[3], field[7], vertlist[11], normlist[11]);
#endif

	// output triangle vertices
	uint numVerts = tex1Dfetch(numVertsTex, cubeindex);
	for(int i=0; i<numVerts; i++) {
		uint edge = tex1Dfetch(triTex, cubeindex*16 + i);

		uint index = numVertsScanned[voxel] + i;
		if (index < maxVerts) {
#if USE_SHARED
			//pos[index]  = make_float4(vertlist[(edge*NTHREADS)+threadIdx.x], 1.0f);
			//norm[index] = make_float4(normlist[(edge*NTHREADS)+threadIdx.x], 0.0f);
			pos[index]  = make_float4(vertlist[12*threadIdx.x+edge], 1.0f);
			norm[index] = make_float4(normlist[12*threadIdx.x+edge], 0.0f);
#else
			pos[index] = make_float4(vertlist[edge], 1.0f);
			norm[index] = make_float4(normlist[edge], 0.0f);
#endif
		}
	}
}

					



// version that calculates flat surface normal for each triangle
__global__ void
generateTriangles2(float4 *pos, float4 *norm, uint *compactedVoxelArray, uint *numVertsScanned, float *volume,
				   uint3 ngrids, float3 voxel_h, float3 voxel_min, float isoValue, uint activeVoxels, uint maxVerts)
{
	uint blockId = __mul24(blockIdx.y, gridDim.x) + blockIdx.x;
	uint i = __mul24(blockId, blockDim.x) + threadIdx.x;

	if(i > activeVoxels-1){
		i = activeVoxels-1;
	}

#if SKIP_EMPTY_VOXELS
	uint voxel = compactedVoxelArray[i];
#else
	uint voxel = i;
#endif

	uint3 gridPos = calcGridIdxU(voxel, ngrids);

	float3 p;
	p.x = voxel_min.x+(gridPos.x*voxel_h.x);
	p.y = voxel_min.y+(gridPos.y*voxel_h.y);
	p.z = voxel_min.z+(gridPos.z*voxel_h.z);
	//p.x = -1.0f + (gridPos.x * voxelSize.x);
	//p.y = -1.0f + (gridPos.y * voxelSize.y);
	//p.z = -1.0f + (gridPos.z * voxelSize.z);

	// calculate cell vertex positions
	float3 v[8];
	v[0] = p;
	v[1] = p + make_float3(voxel_h.x,		 0,		 0);
	v[2] = p + make_float3(voxel_h.x, voxel_h.y,		 0);
	v[3] = p + make_float3(		0, voxel_h.y,		 0);
	v[4] = p + make_float3(		0,		 0, voxel_h.z);
	v[5] = p + make_float3(voxel_h.x,		 0, voxel_h.z);
	v[6] = p + make_float3(voxel_h.x, voxel_h.y, voxel_h.z);
	v[7] = p + make_float3(		0, voxel_h.y, voxel_h.z);

#if SAMPLE_VOLUME
	float field[8];
	field[0] = sampleVolume(volume, gridPos, ngrids);
	field[1] = sampleVolume(volume, gridPos + make_uint3(1, 0, 0), ngrids);
	field[2] = sampleVolume(volume, gridPos + make_uint3(1, 1, 0), ngrids);
	field[3] = sampleVolume(volume, gridPos + make_uint3(0, 1, 0), ngrids);
	field[4] = sampleVolume(volume, gridPos + make_uint3(0, 0, 1), ngrids);
	field[5] = sampleVolume(volume, gridPos + make_uint3(1, 0, 1), ngrids);
	field[6] = sampleVolume(volume, gridPos + make_uint3(1, 1, 1), ngrids);
	field[7] = sampleVolume(volume, gridPos + make_uint3(0, 1, 1), ngrids);
#else
	// evaluate field values
	float field[8];
	field[0] = fieldFunc(v[0]);
	field[1] = fieldFunc(v[1]);
	field[2] = fieldFunc(v[2]);
	field[3] = fieldFunc(v[3]);
	field[4] = fieldFunc(v[4]);
	field[5] = fieldFunc(v[5]);
	field[6] = fieldFunc(v[6]);
	field[7] = fieldFunc(v[7]);
#endif

	// recalculate flag
	uint cubeindex;
	cubeindex =  uint(field[0] < isoValue); 
	cubeindex += uint(field[1] < isoValue)*2; 
	cubeindex += uint(field[2] < isoValue)*4; 
	cubeindex += uint(field[3] < isoValue)*8; 
	cubeindex += uint(field[4] < isoValue)*16; 
	cubeindex += uint(field[5] < isoValue)*32; 
	cubeindex += uint(field[6] < isoValue)*64; 
	cubeindex += uint(field[7] < isoValue)*128;

	// find the vertices where the surface intersects the cube 

#if USE_SHARED
	// use shared memory to avoid using local
	__shared__ float3 vertlist[12*NTHREADS];

	//vertlist[threadIdx.x] = vertexInterp(isoValue, v[0], v[1], field[0], field[1]);
	//vertlist[NTHREADS+threadIdx.x] = vertexInterp(isoValue, v[1], v[2], field[1], field[2]);
	//vertlist[(NTHREADS*2)+threadIdx.x] = vertexInterp(isoValue, v[2], v[3], field[2], field[3]);
	//vertlist[(NTHREADS*3)+threadIdx.x] = vertexInterp(isoValue, v[3], v[0], field[3], field[0]);
	//vertlist[(NTHREADS*4)+threadIdx.x] = vertexInterp(isoValue, v[4], v[5], field[4], field[5]);
	//vertlist[(NTHREADS*5)+threadIdx.x] = vertexInterp(isoValue, v[5], v[6], field[5], field[6]);
	//vertlist[(NTHREADS*6)+threadIdx.x] = vertexInterp(isoValue, v[6], v[7], field[6], field[7]);
	//vertlist[(NTHREADS*7)+threadIdx.x] = vertexInterp(isoValue, v[7], v[4], field[7], field[4]);
	//vertlist[(NTHREADS*8)+threadIdx.x] = vertexInterp(isoValue, v[0], v[4], field[0], field[4]);
	//vertlist[(NTHREADS*9)+threadIdx.x] = vertexInterp(isoValue, v[1], v[5], field[1], field[5]);
	//vertlist[(NTHREADS*10)+threadIdx.x] = vertexInterp(isoValue, v[2], v[6], field[2], field[6]);
	//vertlist[(NTHREADS*11)+threadIdx.x] = vertexInterp(isoValue, v[3], v[7], field[3], field[7]);
	vertlist[12*threadIdx.x+0]  = vertexInterp(isoValue, v[0], v[1], field[0], field[1]);
	vertlist[12*threadIdx.x+1]  = vertexInterp(isoValue, v[1], v[2], field[1], field[2]);
	vertlist[12*threadIdx.x+2]  = vertexInterp(isoValue, v[2], v[3], field[2], field[3]);
	vertlist[12*threadIdx.x+3]  = vertexInterp(isoValue, v[3], v[0], field[3], field[0]);
	vertlist[12*threadIdx.x+4]  = vertexInterp(isoValue, v[4], v[5], field[4], field[5]);
	vertlist[12*threadIdx.x+5]  = vertexInterp(isoValue, v[5], v[6], field[5], field[6]);
	vertlist[12*threadIdx.x+6]  = vertexInterp(isoValue, v[6], v[7], field[6], field[7]);
	vertlist[12*threadIdx.x+7]  = vertexInterp(isoValue, v[7], v[4], field[7], field[4]);
	vertlist[12*threadIdx.x+8]  = vertexInterp(isoValue, v[0], v[4], field[0], field[4]);
	vertlist[12*threadIdx.x+9]  = vertexInterp(isoValue, v[1], v[5], field[1], field[5]);
	vertlist[12*threadIdx.x+10] = vertexInterp(isoValue, v[2], v[6], field[2], field[6]);
	vertlist[12*threadIdx.x+11] = vertexInterp(isoValue, v[3], v[7], field[3], field[7]);
   __syncthreads();
#else

	float3 vertlist[12];

	vertlist[0] = vertexInterp(isoValue, v[0], v[1], field[0], field[1]);
	vertlist[1] = vertexInterp(isoValue, v[1], v[2], field[1], field[2]);
	vertlist[2] = vertexInterp(isoValue, v[2], v[3], field[2], field[3]);
	vertlist[3] = vertexInterp(isoValue, v[3], v[0], field[3], field[0]);

	vertlist[4] = vertexInterp(isoValue, v[4], v[5], field[4], field[5]);
	vertlist[5] = vertexInterp(isoValue, v[5], v[6], field[5], field[6]);
	vertlist[6] = vertexInterp(isoValue, v[6], v[7], field[6], field[7]);
	vertlist[7] = vertexInterp(isoValue, v[7], v[4], field[7], field[4]);

	vertlist[8] = vertexInterp(isoValue, v[0], v[4], field[0], field[4]);
	vertlist[9] = vertexInterp(isoValue, v[1], v[5], field[1], field[5]);
	vertlist[10] = vertexInterp(isoValue, v[2], v[6], field[2], field[6]);
	vertlist[11] = vertexInterp(isoValue, v[3], v[7], field[3], field[7]);
#endif

	// output triangle vertices
	uint numVerts = tex1Dfetch(numVertsTex, cubeindex);
	for(int i=0; i<numVerts; i+=3) {
		uint index = numVertsScanned[voxel] + i;


/*
		uint edge;
		edge = tex1Dfetch(triTex, (cubeindex*16) + i);
		float3 *v[3];
#if USE_SHARED
		//v[0] = &vertlist[(edge*NTHREADS)+threadIdx.x];
		v[0] = &vertlist[12*threadIdx.x+edge];
#else
		v[0] = &vertlist[edge];
#endif

		edge = tex1Dfetch(triTex, (cubeindex*16) + i + 1);
#if USE_SHARED
		//v[1] = &vertlist[(edge*NTHREADS)+threadIdx.x];
		v[1] = &vertlist[12*threadIdx.x+edge];
#else
		v[1] = &vertlist[edge];
#endif

		edge = tex1Dfetch(triTex, (cubeindex*16) + i + 2);
#if USE_SHARED
		//v[2] = &vertlist[(edge*NTHREADS)+threadIdx.x];
		v[2] = &vertlist[12*threadIdx.x+edge];
#else
		v[2] = &vertlist[edge];
#endif
*/
		// calculate triangle surface normal
		float3 n;

		if(index < (maxVerts-3)){
			//n = calcNormal(v[0], v[1], v[2]);
/*
			n = calcNormalF(volume, *v[0], ngrids, voxel_h);
			pos[index] = make_float4(*v[0], 1.0f);
			norm[index] = make_float4(n, 0.0f);

			n = calcNormalF(volume, *v[1], ngrids, voxel_h);
			pos[index+1] = make_float4(*v[1], 1.0f);
			norm[index+1] = make_float4(n, 0.0f);

			n = calcNormalF(volume, *v[2], ngrids, voxel_h);
			pos[index+2] = make_float4(*v[2], 1.0f);
			norm[index+2] = make_float4(n, 0.0f);
*/

			uint edge;
			edge = tex1Dfetch(triTex, (cubeindex*16) + i);
#if USE_SHARED
			n = calcNormalF(volume, vertlist[12*threadIdx.x+edge], ngrids, voxel_h);
			pos[index] = make_float4(vertlist[12*threadIdx.x+edge], 1.0f);
#else
			n = calcNormalF(volume, vertlist[edge], ngrids, voxel_h);
			pos[index] = make_float4(vertlist[edge], 1.0f);
#endif
			norm[index] = make_float4(n, 0.0f);

			edge = tex1Dfetch(triTex, (cubeindex*16) + i + 1);
#if USE_SHARED
			n = calcNormalF(volume, vertlist[12*threadIdx.x+edge], ngrids, voxel_h);
			pos[index+1] = make_float4(vertlist[12*threadIdx.x+edge], 1.0f);
#else
			n = calcNormalF(volume, vertlist[edge], ngrids, voxel_h);
			pos[index+1] = make_float4(vertlist[edge], 1.0f);
#endif
			norm[index+1] = make_float4(n, 0.0f);

			edge = tex1Dfetch(triTex, (cubeindex*16) + i + 2);
#if USE_SHARED
			n = calcNormalF(volume, vertlist[12*threadIdx.x+edge], ngrids, voxel_h);
			pos[index+2] = make_float4(vertlist[12*threadIdx.x+edge], 1.0f);
#else
			n = calcNormalF(volume, vertlist[edge], ngrids, voxel_h);
			pos[index+2] = make_float4(vertlist[edge], 1.0f);
#endif
			norm[index+2] = make_float4(n, 0.0f);
		}
	}
}






#endif // #ifndef _RX_CUMC_KERNEL_CU_



