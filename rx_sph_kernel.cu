/*! 
  @file rx_sph_kernel.cu
	
  @brief CUDA�ɂ��SPH
 
  @author Makoto Fujisawa
  @date 2009-08, 2011-06
*/
// FILE --rx_sph_kernel.cu--

#ifndef _RX_CUSPH_KERNEL_CU_
#define _RX_CUSPH_KERNEL_CU_


//-----------------------------------------------------------------------------
// �C���N���[�h�t�@�C��
//-----------------------------------------------------------------------------
#include "rx_cu_common.cu"


//-----------------------------------------------------------------------------
// �n�b�V��
//-----------------------------------------------------------------------------
/*!
 * �e�p�[�e�B�N���̃O���b�h�n�b�V���l
 * @param[out] gridParticleHash
 * @param[out] dSortedIndex
 * @param[in] pos �p�[�e�B�N���ʒu���i�[�����z��
 * @param[in] nprts �p�[�e�B�N����
 */
__global__
void calcHashD(uint*   dGridParticleHash, 
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
/*!
 * �e�p�[�e�B�N���̃O���b�h�n�b�V���l
 * @param[out] gridParticleHash
 * @param[out] dSortedIndex
 * @param[in] pos �p�[�e�B�N���ʒu���i�[�����z��
 * @param[in] nprts �p�[�e�B�N����
 */
__global__
void calcHashD2(uint*   dGridParticleHash, 
				uint*   dSortedIndex, 
				float4* dPos, 
				int*    dAttr, 
				uint	   nprts)
{
	uint index = __umul24(blockIdx.x, blockDim.x)+threadIdx.x;
	if(index >= nprts) return;
	
	uint hash = 2147483647;
	if(dAttr[index] != -1){
		volatile float4 p = dPos[index];
		int3 gridPos = calcGridPos(make_float3(p.x, p.y, p.z));
		hash = calcGridHash(gridPos);
	}

	dGridParticleHash[index] = hash;
	dSortedIndex[index] = index;
}
/*!
 * �e�p�[�e�B�N���̃O���b�h�n�b�V���l
 * @param[out] gridParticleHash
 * @param[out] dSortedIndex
 * @param[in] pos �p�[�e�B�N���ʒu���i�[�����z��
 * @param[in] nprts �p�[�e�B�N����
 */
__global__
void calcHashB(uint*   dGridParticleHash, 
			   uint*   dSortedIndex, 
			   float4*  dPos, 
			   float3  world_origin, 
			   float3  cell_width, 
			   uint3   grid_size, 
			   uint	   nprts)
{
	uint index = __umul24(blockIdx.x, blockDim.x)+threadIdx.x;
	if(index >= nprts) return;
	
	float3 p = make_float3(dPos[index]);
	int3 gridPos = calcGridPosB(make_float3(p.x, p.y, p.z), world_origin, cell_width, grid_size);
	uint hash = calcGridHashB(gridPos, grid_size);

	dGridParticleHash[index] = hash;
	dSortedIndex[index] = index;
}


/*!
 * �p�[�e�B�N���f�[�^���\�[�g���āC�n�b�V�����̊e�Z���̍ŏ��̃A�h���X������
 * @param[in] cell �p�[�e�B�N���O���b�h�f�[�^
 * @param[in] oldPos �p�[�e�B�N���ʒu
 * @param[in] oldVel �p�[�e�B�N�����x
 */
__global__
void reorderDataAndFindCellStartD(rxParticleCell cell, float4* dSortedPos, float4* dSortedVel)
{
	extern __shared__ uint sharedHash[];	// �T�C�Y : blockSize+1
	uint index = __umul24(blockIdx.x,blockDim.x)+threadIdx.x;
	
	uint hash;
	if(index < cell.uNumParticles){
		hash = cell.dGridParticleHash[index];	// �n�b�V���l

		sharedHash[threadIdx.x+1] = hash;	// �n�b�V���l���V�F�A�[�h�������Ɋi�[

		if(index > 0 && threadIdx.x == 0){
			// �e�V�F�A�[�h�������̍ŏ��ׂ͗̃O���b�h�̃p�[�e�B�N���̃n�b�V���l���i�[
			sharedHash[0] = cell.dGridParticleHash[index-1];
		}
	}

	__syncthreads();
	
	if(index < cell.uNumParticles){
		// �C���f�b�N�X0�ł���C�������́C��O�̃p�[�e�B�N���̃O���b�h�n�b�V���l���قȂ�ꍇ�C
		// �p�[�e�B�N���͕����̈�̍ŏ�
		if(index == 0 || hash != sharedHash[threadIdx.x]){
			cell.dCellStart[hash] = index;
			if(index > 0){
				// ��O�̃p�[�e�B�N���́C��O�̕����̈�̍Ō�
				cell.dCellEnd[sharedHash[threadIdx.x]] = index;
			}
		}

		// �C���f�b�N�X���Ō�Ȃ�΁C�����̈�̍Ō�
		if(index == cell.uNumParticles-1){
			cell.dCellEnd[hash] = index+1;
		}

		// �ʒu�Ƒ��x�̃f�[�^����ёւ�
		// �\�[�g�����C���f�b�N�X�ŎQ�Ƃ��\�����T�����̃O���[�o���������A�N�Z�X���ɗ͗}���邽�߂Ƀf�[�^���̂��̂���ёւ���
		uint sortedIndex = cell.dSortedIndex[index];
		float4 pos = FETCH(dSortedPos, sortedIndex);
		float4 vel = FETCH(dSortedVel, sortedIndex);

		cell.dSortedPos[index] = pos;
		cell.dSortedVel[index] = vel;
	}
}

/*!
 * �p�[�e�B�N���f�[�^���\�[�g���āC�n�b�V�����̊e�Z���̍ŏ��̃A�h���X������
 * @param[in] cell �p�[�e�B�N���O���b�h�f�[�^
 * @param[in] oldPos �p�[�e�B�N���ʒu
 * @param[in] oldVel �p�[�e�B�N�����x
 */
__global__
void reorderDataAndFindCellStartB(rxParticleCell cell, float4* dPos)
{
	extern __shared__ uint sharedHash[];	// �T�C�Y : blockSize+1
	uint index = __umul24(blockIdx.x,blockDim.x)+threadIdx.x;
	
	uint hash;
	if(index < cell.uNumParticles){
		hash = cell.dGridParticleHash[index];	// �n�b�V���l

		sharedHash[threadIdx.x+1] = hash;	// �n�b�V���l���V�F�A�[�h�������Ɋi�[

		if(index > 0 && threadIdx.x == 0){
			// �e�V�F�A�[�h�������̍ŏ��ׂ͗̃O���b�h�̃p�[�e�B�N���̃n�b�V���l���i�[
			sharedHash[0] = cell.dGridParticleHash[index-1];
		}
	}

	__syncthreads();
	
	if(index < cell.uNumParticles){
		// �C���f�b�N�X0�ł���C�������́C��O�̃p�[�e�B�N���̃O���b�h�n�b�V���l���قȂ�ꍇ�C
		// �p�[�e�B�N���͕����̈�̍ŏ�
		if(index == 0 || hash != sharedHash[threadIdx.x]){
			cell.dCellStart[hash] = index;
			if(index > 0){
				// ��O�̃p�[�e�B�N���́C��O�̕����̈�̍Ō�
				cell.dCellEnd[sharedHash[threadIdx.x]] = index;
			}
		}

		// �C���f�b�N�X���Ō�Ȃ�΁C�����̈�̍Ō�
		if(index == cell.uNumParticles-1){
			cell.dCellEnd[hash] = index+1;
		}

		// �ʒu�Ƒ��x�̃f�[�^����ёւ�
		// �\�[�g�����C���f�b�N�X�ŎQ�Ƃ��\�����T�����̃O���[�o���������A�N�Z�X���ɗ͗}���邽�߂Ƀf�[�^���̂��̂���ёւ���
		uint sortedIndex = cell.dSortedIndex[index];
		float4 pos = dPos[sortedIndex];
		cell.dSortedPos[index] = pos;
	}
}

//-----------------------------------------------------------------------------
// MARK:SPH
//-----------------------------------------------------------------------------
/*!
 * �^����ꂽ�Z�����̃p�[�e�B�N���Ƃ̋������疧�x���v�Z
 * @param[in] gridPos �O���b�h�ʒu
 * @param[in] index �p�[�e�B�N���C���f�b�N�X
 * @param[in] pos �v�Z���W
 * @param[in] cell �p�[�e�B�N���O���b�h�f�[�^
 * @return �Z�����̃p�[�e�B�N������v�Z�������x�l
 */
__device__
float calDensityCell(int3 gridPos, uint i, float3 pos0, rxParticleCell cell)
{
	uint gridHash = calcGridHash(gridPos);

	// �Z�����̃p�[�e�B�N���̃X�^�[�g�C���f�b�N�X
	uint startIndex = FETCHC(dCellStart, gridHash);

	float h = params.EffectiveRadius;
	float dens = 0.0f;
	if(startIndex != 0xffffffff){	// �Z������łȂ����̃`�F�b�N
		// �Z�����̃p�[�e�B�N���Ŕ���
		uint endIndex = FETCHC(dCellEnd, gridHash);
		for(uint j = startIndex; j < endIndex; ++j){
			//if(j == i) continue;

			float3 pos1 = make_float3(FETCHC(dSortedPos, j));

			float3 rij = pos0-pos1;
			float r = length(rij);

			if(r <= h){
				float q = h*h-r*r;
				dens += params.Mass*params.Wpoly6*q*q*q;
			}
		}
	}

	return dens;
}

/*!
 * �p�[�e�B�N�����x�v�Z(�J�[�l���֐�)
 * @param[out] newDens �p�[�e�B�N�����x
 * @param[out] newPres �p�[�e�B�N������
 * @param[in]  cell �p�[�e�B�N���O���b�h�f�[�^
 */
__global__
void sphCalDensity(float* newDens, float* newPres, rxParticleCell cell)
{
	// �p�[�e�B�N���C���f�b�N�X
	uint index = __mul24(blockIdx.x,blockDim.x)+threadIdx.x;
	if(index >= cell.uNumParticles) return;	
	
	float3 pos = make_float3(FETCHC(dSortedPos, index));	// �p�[�e�B�N���ʒu
	//int3 grid_pos = calcGridPos(pos);	// �p�[�e�B�N����������O���b�h�ʒu
	float h = params.EffectiveRadius;

	// �p�[�e�B�N�����͂̃O���b�h
	int3 grid_pos0, grid_pos1;
	grid_pos0 = calcGridPos(pos-make_float3(h));
	grid_pos1 = calcGridPos(pos+make_float3(h));

	// ���͂̃O���b�h���܂߂ċߖT�T���C���x�v�Z
	float dens = 0.0f;
	for(int z = grid_pos0.z; z <= grid_pos1.z; ++z){
		for(int y = grid_pos0.y; y <= grid_pos1.y; ++y){
			for(int x = grid_pos0.x; x <= grid_pos1.x; ++x){
				int3 n_grid_pos = make_int3(x, y, z);
				dens += calDensityCell(n_grid_pos, index, pos, cell);
			}
		}
	}

	// �K�X�萔���g�������͎Z�o
	float pres;
	pres = params.GasStiffness*(dens-params.Density);

	// ���x�ƈ��͒l�����ʂɏ�������
	uint oIdx = cell.dSortedIndex[index];
	newDens[oIdx] = dens;
	newPres[oIdx] = pres;
}




/*!
 * �^����ꂽ�Z�����̃p�[�e�B�N���Ƃ̋�������@�����v�Z
 * @param[in] gridPos �O���b�h�ʒu
 * @param[in] i �p�[�e�B�N���C���f�b�N�X
 * @param[in] pos �v�Z���W
 * @param[in] cell �p�[�e�B�N���O���b�h�f�[�^
 * @return �Z�����̃p�[�e�B�N������v�Z�������x�l
 */
__device__
float3 calNormalCell(int3 gridPos, uint i, float3 pos0, float* dens, rxParticleCell cell)
{
	uint gridHash = calcGridHash(gridPos);

	// �Z�����̃p�[�e�B�N���̃X�^�[�g�C���f�b�N�X
	uint startIndex = FETCHC(dCellStart, gridHash);

	float h = params.EffectiveRadius;
	float3 nrm = make_float3(0.0f);
	if(startIndex != 0xffffffff){	// �Z������łȂ����̃`�F�b�N
		// �Z�����̃p�[�e�B�N���Ŕ���
		uint endIndex = FETCHC(dCellEnd, gridHash);

		for(uint j = startIndex; j < endIndex; ++j){
			if(j != i){
				float3 pos1 = make_float3(FETCHC(dSortedPos, j));

				float3 rij = pos0-pos1;
				float r = length(rij);

				if(r <= h && r > 0.0001){
					float d1 = dens[cell.dSortedIndex[j]];
					float q = h*h-r*r;

					nrm += (params.Mass/d1)*params.GWpoly6*q*q*rij;
				}

			}
		}
	}

	return nrm;
}


/*!
 * �p�[�e�B�N���@���v�Z(�J�[�l���֐�)
 * @param[out] newNrms �p�[�e�B�N���@��
 * @param[in] dens �p�[�e�B�N�����x
 * @param[in] cell �p�[�e�B�N���O���b�h�f�[�^
 */
__global__
void sphCalNormal(float4* newNrms, float* dens, rxParticleCell cell)
{
	uint index = __mul24(blockIdx.x,blockDim.x)+threadIdx.x;
	if(index >= cell.uNumParticles) return;	
	
	float3 pos = make_float3(FETCHC(dSortedPos, index));	// �p�[�e�B�N���ʒu
	float h = params.EffectiveRadius;
	//int3 grid_pos = calcGridPos(pos);	// �p�[�e�B�N����������O���b�h�ʒu

	// �p�[�e�B�N�����͂̃O���b�h
	int3 grid_pos0, grid_pos1;
	grid_pos0 = calcGridPos(pos-make_float3(h));
	grid_pos1 = calcGridPos(pos+make_float3(h));

	// ���͂̃O���b�h���܂߂ċߖT�T���C���x�v�Z
	float3 nrm = make_float3(0.0f);
	for(int z = grid_pos0.z; z <= grid_pos1.z; ++z){
		for(int y = grid_pos0.y; y <= grid_pos1.y; ++y){
			for(int x = grid_pos0.x; x <= grid_pos1.x; ++x){
				int3 n_grid_pos = make_int3(x, y, z);
				nrm += calNormalCell(n_grid_pos, index, pos, dens, cell);
			}
		}
	}

	float l = length(nrm);
	if(l > 0){
		nrm /= l;
	}

	uint oIdx = cell.dSortedIndex[index];
	newNrms[oIdx] = make_float4(nrm, 0.0f);
}




/*!
 * �^����ꂽ�Z�����̃p�[�e�B�N���Ƃ̋�������͏���v�Z
 * @param[in] gridPos �O���b�h�ʒu
 * @param[in] i �p�[�e�B�N���C���f�b�N�X
 * @param[in] pos0 �v�Z���W
 * @param[in] vel0 �v�Z���W�̑��x
 * @param[in] dens0 �v�Z���W�̖��x
 * @param[in] pres0 �v�Z���W�̈���
 * @param[in] dens �p�[�e�B�N�����x
 * @param[in] pres �p�[�e�B�N������
 * @param[in] cell �p�[�e�B�N���O���b�h�f�[�^
 * @return �Z�����̃p�[�e�B�N������v�Z�����͏�
 */
__device__
float3 calForceCell(int3 gridPos, uint i, float3 pos0, float3 vel0, float dens0, float pres0, 
					float* dens, float* pres, rxParticleCell cell)
{
	uint gridHash = calcGridHash(gridPos);

	// �Z�����̃p�[�e�B�N���̃X�^�[�g�C���f�b�N�X
	uint startIndex = FETCHC(dCellStart, gridHash);

	float h = params.EffectiveRadius;

	float3 frc = make_float3(0.0f);
	if(startIndex != 0xffffffff){	// �Z������łȂ����̃`�F�b�N
		// �Z�����̃p�[�e�B�N���Ŕ���
		uint endIndex = FETCHC(dCellEnd, gridHash);
		float prsi = pres0/(dens0*dens0);
		for(uint j = startIndex; j < endIndex; ++j){
			if(j != i){
				// �ߖT�p�[�e�B�N���̃p�����[�^
				float3 pos1 = make_float3(FETCHC(dSortedPos, j));
				float3 vel1 = make_float3(FETCHC(dSortedVel, j));

				float3 rij = pos0-pos1;
				float r = length(rij);

				if(r <= h && r > 0.0001){
					//float3 vel1 = make_float3(vel[cell.dSortedIndex[j]]);
					float dens1 = dens[cell.dSortedIndex[j]];
					float pres1 = pres[cell.dSortedIndex[j]];

					float3 vji = vel1-vel0;

					float prsj = pres1/(dens1*dens1);
					float q = h-r;

					// ���͍�
					frc += -params.Mass*(prsi+prsj)*params.GWspiky*q*q*rij/r;

					// �S����
					frc += params.Viscosity*params.Mass*(vji/dens1)*params.LWvisc*q;
				}
			}
		}
	}

	return frc;
}

/*!
 * �p�[�e�B�N���ɂ�����͂̌v�Z(�J�[�l���֐�)
 * @param[in] dens �p�[�e�B�N�����x
 * @param[in] pres �p�[�e�B�N������
 * @param[out] outFrc �p�[�e�B�N���ɂ������
 * @param[in] cell �p�[�e�B�N���O���b�h�f�[�^
 */
__global__
void sphCalForces(float* dens, float* pres, float4* outFrc, rxParticleCell cell)
{
	uint index = __mul24(blockIdx.x,blockDim.x)+threadIdx.x;
	if(index >= cell.uNumParticles) return;	
	
	// �\�[�g�ςݔz�񂩂�p�[�e�B�N���f�[�^���擾
	float3 pos0 = make_float3(FETCHC(dSortedPos, index));
	float3 vel0 = make_float3(FETCHC(dSortedVel, index));

	int3 gridPos0 = calcGridPos(pos0);

	// �p�[�e�B�N���̃\�[�g�Ȃ��z���ł̃C���f�b�N�X
	uint oIdx = cell.dSortedIndex[index];

	float dens0 = dens[oIdx];
	float pres0 = pres[oIdx];

	float h = params.EffectiveRadius;
	int3 grid_pos0, grid_pos1;
	grid_pos0 = calcGridPos(pos0-make_float3(h));
	grid_pos1 = calcGridPos(pos0+make_float3(h));

	// ���͂̃O���b�h���܂߂ċߖT�T���C���͍��C�S�������v�Z
	float3 frc = make_float3(0.0f);
	for(int z = grid_pos0.z; z <= grid_pos1.z; ++z){
		for(int y = grid_pos0.y; y <= grid_pos1.y; ++y){
			for(int x = grid_pos0.x; x <= grid_pos1.x; ++x){
				int3 n_grid_pos = make_int3(x, y, z);

				frc += calForceCell(n_grid_pos, index, pos0, vel0, dens0, pres0, dens, pres, cell);
			}
		}
	}

	// �O��(�d��)
	frc += params.Gravity;

	outFrc[oIdx] = make_float4(frc, 0.0f);
}


/*!
 * �^����ꂽ�Z�����̃p�[�e�B�N���Ƃ̋������狫�E�p�[�e�B�N���̑̐ς��v�Z
 * @param[in] gridPos �O���b�h�ʒu
 * @param[in] index �p�[�e�B�N���C���f�b�N�X
 * @param[in] pos �v�Z���W
 * @param[in] cell �p�[�e�B�N���O���b�h�f�[�^
 * @return �Z�����̃p�[�e�B�N������v�Z�������x�l
 */
__device__
float calBoundaryVolumeCell(int3 gridPos, uint i, float3 pos0, rxParticleCell cell)
{
	uint gridHash = calcGridHashB(gridPos, params.GridSizeB);

	// �Z�����̃p�[�e�B�N���̃X�^�[�g�C���f�b�N�X
	uint startIndex = cell.dCellStart[gridHash];

	float h = params.EffectiveRadius;
	float mw = 0.0f;
	if(startIndex != 0xffffffff){	// �Z������łȂ����̃`�F�b�N
		// �Z�����̃p�[�e�B�N���Ŕ���
		uint endIndex = cell.dCellEnd[gridHash];
		for(uint j = startIndex; j < endIndex; ++j){
			float3 pos1 = make_float3(cell.dSortedPos[j]);

			float3 rij = pos0-pos1;
			float r = length(rij);

			if(r <= h){
				float q = h*h-r*r;
				mw += params.Mass*params.Wpoly6*q*q*q;
			}
		}
	}

	return mw;
}

/*!
 * ���E�p�[�e�B�N���̑̐όv�Z(�J�[�l���֐�)
 * @param[out] newVolB �p�[�e�B�N���̐�
 * @param[in]  cell ���E�p�[�e�B�N���O���b�h�f�[�^
 */
__global__
void sphCalBoundaryVolume(float* newVolB, rxParticleCell cell)
{
	// �p�[�e�B�N���C���f�b�N�X
	uint index = __mul24(blockIdx.x,blockDim.x)+threadIdx.x;
	if(index >= cell.uNumParticles) return;	
	
	float3 pos = make_float3(cell.dSortedPos[index]);	// �p�[�e�B�N���ʒu
	//int3 grid_pos = calcGridPos(pos);	// �p�[�e�B�N����������O���b�h�ʒu
	float h = params.EffectiveRadius;

	// �p�[�e�B�N�����͂̃O���b�h
	int3 grid_pos0, grid_pos1;
	grid_pos0 = calcGridPosB(pos-make_float3(h), params.WorldOriginB, params.CellWidthB, params.GridSizeB);
	grid_pos1 = calcGridPosB(pos+make_float3(h), params.WorldOriginB, params.CellWidthB, params.GridSizeB);


	// ���͂̃O���b�h���܂߂ċߖT�T��
	float mw = 0.0f;
	for(int z = grid_pos0.z; z <= grid_pos1.z; ++z){
		for(int y = grid_pos0.y; y <= grid_pos1.y; ++y){
			for(int x = grid_pos0.x; x <= grid_pos1.x; ++x){
				int3 n_grid_pos = make_int3(x, y, z);
				mw += calBoundaryVolumeCell(n_grid_pos, index, pos, cell);
			}
		}
	}

	// �̐ς����ʂɏ�������
	uint oIdx = cell.dSortedIndex[index];
	newVolB[oIdx] = params.Mass/mw;
}

/*!
 * �^����ꂽ�Z�����̃p�[�e�B�N���Ƃ̋������疧�x���v�Z
 * @param[in] gridPos �O���b�h�ʒu
 * @param[in] index �p�[�e�B�N���C���f�b�N�X
 * @param[in] pos �v�Z���W
 * @param[in] cell �p�[�e�B�N���O���b�h�f�[�^
 * @return �Z�����̃p�[�e�B�N������v�Z�������x�l
 */
__device__
float calBoundaryDensityCell(int3 gridPos, uint i, float3 pos0, float* dVolB, rxParticleCell bcell)
{
	uint gridHash = calcGridHashB(gridPos, params.GridSizeB);

	// �Z�����̃p�[�e�B�N���̃X�^�[�g�C���f�b�N�X
	uint startIndex = bcell.dCellStart[gridHash];

	float h = params.EffectiveRadius;
	float dens = 0.0f;
	if(startIndex != 0xffffffff){	// �Z������łȂ����̃`�F�b�N
		// �Z�����̃p�[�e�B�N���Ŕ���
		uint endIndex = bcell.dCellEnd[gridHash];
		for(uint j = startIndex; j < endIndex; ++j){
			//if(j == i) continue;

			float3 pos1 = make_float3(bcell.dSortedPos[j]);
			uint jdx = bcell.dSortedIndex[j];

			float3 rij = pos0-pos1;
			float r = length(rij);

			if(r <= h){
				float q = h*h-r*r;
				dens += params.Density*dVolB[jdx]*params.Wpoly6*q*q*q;
			}
		}
	}

	return dens;
}

/*!
 * �p�[�e�B�N�����x�v�Z(�J�[�l���֐�)
 * @param[out] newDens �p�[�e�B�N�����x
 * @param[out] newPres �p�[�e�B�N������
 * @param[in]  cell �p�[�e�B�N���O���b�h�f�[�^
 */
__global__
void sphCalBoundaryDensity(float* newDens, float* newPres, float4* dPos, float* dVolB, rxParticleCell bcell, uint pnum)
{
	// �p�[�e�B�N���C���f�b�N�X
	uint index = __mul24(blockIdx.x,blockDim.x)+threadIdx.x;
	if(index >= pnum) return;	
	
	float3 pos = make_float3(dPos[index]);	// �p�[�e�B�N���ʒu
	float h = params.EffectiveRadius;

	// �p�[�e�B�N�����͂̃O���b�h
	int3 grid_pos0, grid_pos1;
	grid_pos0 = calcGridPosB(pos-make_float3(h), params.WorldOriginB, params.CellWidthB, params.GridSizeB);
	grid_pos1 = calcGridPosB(pos+make_float3(h), params.WorldOriginB, params.CellWidthB, params.GridSizeB);

	// ���͂̃O���b�h���܂߂ċߖT�T���C���x�v�Z
	float dens = 0.0f;
	for(int z = grid_pos0.z; z <= grid_pos1.z; ++z){
		for(int y = grid_pos0.y; y <= grid_pos1.y; ++y){
			for(int x = grid_pos0.x; x <= grid_pos1.x; ++x){
				int3 n_grid_pos = make_int3(x, y, z);
				dens += calBoundaryDensityCell(n_grid_pos, index, pos, dVolB, bcell);
			}
		}
	}

	dens += newDens[index];

	// �K�X�萔���g�������͎Z�o
	float pres;
	pres = params.GasStiffness*(dens-params.Density);

	// ���x�ƈ��͒l�����ʂɏ�������
	newDens[index] = dens;
	newPres[index] = pres;
}


/*!
 * �^����ꂽ�Z�����̃p�[�e�B�N���Ƃ̋������疧�x���v�Z
 * @param[in] gridPos �O���b�h�ʒu
 * @param[in] index �p�[�e�B�N���C���f�b�N�X
 * @param[in] pos �v�Z���W
 * @param[in] cell �p�[�e�B�N���O���b�h�f�[�^
 * @return �Z�����̃p�[�e�B�N������v�Z�������x�l
 */
__device__
float3 calBoundaryForceCell(int3 gridPos, uint i, float3 pos0, float* dVolB, float dens0, float pres0, rxParticleCell bcell)
{
	uint gridHash = calcGridHashB(gridPos, params.GridSizeB);

	// �Z�����̃p�[�e�B�N���̃X�^�[�g�C���f�b�N�X
	uint startIndex = bcell.dCellStart[gridHash];

	float h = params.EffectiveRadius;
	float3 bp = make_float3(0.0f);
	if(startIndex != 0xffffffff){	// �Z������łȂ����̃`�F�b�N
		// �Z�����̃p�[�e�B�N���Ŕ���
		uint endIndex = bcell.dCellEnd[gridHash];
		float prsi = pres0/(dens0*dens0);
		for(uint j = startIndex; j < endIndex; ++j){
			float3 pos1 = make_float3(bcell.dSortedPos[j]);
			uint jdx = bcell.dSortedIndex[j];

			float3 rij = pos0-pos1;
			float r = length(rij);

			if(r <= h && r > 0.0001){
				float q = h-r;
				bp += -params.Density*dVolB[jdx]*prsi*params.GWspiky*q*q*rij/r;
			}
		}
	}

	return bp;
}

/*!
 * ���E�p�[�e�B�N���ɂ��͂̌v�Z(�J�[�l���֐�)
 * @param[out] newDens �p�[�e�B�N�����x
 * @param[out] newPres �p�[�e�B�N������
 * @param[in]  cell �p�[�e�B�N���O���b�h�f�[�^
 */
__global__
void sphCalBoundaryForce(float* dDens, float* dPres, float4* dPos, float* dVolB, float4* outFrc, rxParticleCell bcell, uint pnum)
{
	// �p�[�e�B�N���C���f�b�N�X
	uint index = __mul24(blockIdx.x,blockDim.x)+threadIdx.x;
	if(index >= pnum) return;	
	
	float3 pos = make_float3(dPos[index]);	// �p�[�e�B�N���ʒu
	float h = params.EffectiveRadius;

	// �p�[�e�B�N�����͂̃O���b�h
	int3 grid_pos0, grid_pos1;
	grid_pos0 = calcGridPosB(pos-make_float3(h), params.WorldOriginB, params.CellWidthB, params.GridSizeB);
	grid_pos1 = calcGridPosB(pos+make_float3(h), params.WorldOriginB, params.CellWidthB, params.GridSizeB);

	// ���x�ƈ���
	float dens0 = dDens[index];
	float pres0 = dPres[index];

	// ���͂̃O���b�h���܂߂ċߖT�T���C���x�v�Z
	float3 frc = make_float3(0.0f);
	for(int z = grid_pos0.z; z <= grid_pos1.z; ++z){
		for(int y = grid_pos0.y; y <= grid_pos1.y; ++y){
			for(int x = grid_pos0.x; x <= grid_pos1.x; ++x){
				int3 n_grid_pos = make_int3(x, y, z);
				frc += calBoundaryForceCell(n_grid_pos, index, pos, dVolB, dens0, pres0, bcell);
			}
		}
	}

	// ���x�ƈ��͒l�����ʂɏ�������
	outFrc[index] += make_float4(frc, 0.0f);
}

__device__
void calCollisionSolid(float3 &pos, float3 &vel, float dt)
{
	float d;
	float3 n;
	float3 cp;

	// �{�b�N�X�`��̃I�u�W�F�N�g�Ƃ̏Փ�
#if MAX_BOX_NUM
	for(int i = 0; i < params.BoxNum; ++i){
		if(params.BoxFlg[i] == 0) continue;
		
		collisionPointBox(pos, params.BoxCen[i], params.BoxExt[i], params.BoxRot[i], params.BoxInvRot[i], cp, d, n);

		if(d < 0.0){
			float res = params.Restitution;
			res = (res > 0) ? (res*fabs(d)/(dt*length(vel))) : 0.0f;
			vel -= (1+res)*n*dot(n, vel);
			pos = cp;
		}
	}
#endif

	// ���`��̃I�u�W�F�N�g�Ƃ̏Փ�
#if MAX_SPHERE_NUM
	for(int i = 0; i < params.SphereNum; ++i){
		if(params.SphereFlg[i] == 0) continue;

		collisionPointSphere(pos, params.SphereCen[i], params.SphereRad[i], cp, d, n);

		if(d < 0.0){
			float res = params.Restitution;
			res = (res > 0) ? (res*fabs(d)/(dt*length(vel))) : 0.0f;
			vel -= (1+res)*n*dot(n, vel);
			pos = cp;
		}
	}
#endif

	// ���͂̋��E�Ƃ̏Փ˔���
	float3 l0 = params.BoundaryMin;
	float3 l1 = params.BoundaryMax;
	collisionPointAABB(pos, 0.5*(l1+l0), 0.5*(l1-l0), cp, d, n);

	if(d < 0.0){
		float res = params.Restitution;
		res = (res > 0) ? (res*fabs(d)/(dt*length(vel))) : 0.0f;
		vel -= (1+res)*n*dot(n, vel);
		pos = cp;
	}
}

__device__
inline bool calCollisionPolygon(float3 &pos0, float3 &pos1, float3 &vel, float3 v0, float3 v1, float3 v2, float dt)
{
	float3 cp, n;
	if(intersectSegmentTriangle(pos0, pos1, v0, v1, v2, cp, n, params.ParticleRadius) == 1){
		float d = length(pos1-cp);
		n = normalize(n);

		float res = params.Restitution;
		res = (res > 0) ? (res*fabs(d)/(dt*length(vel))) : 0.0f;
		float3 vr = -(1+res)*n*dot(n, vel);

		float l = length(pos1-pos0);
		pos1 = cp+vr*(dt*d/l);
		vel += vr;//+params.PolyVel[0];
		//vel.x = 1.0;

		return true;
	}
	return false;
}



/*!
 * �p�[�e�B�N���ʒu�C���x�̍X�V
 * @param[inout] ppos �p�[�e�B�N���ʒu
 * @param[inout] pvel �p�[�e�B�N�����x
 * @param[in] pfrc �p�[�e�B�N���ɂ������
 * @param[in] dens �p�[�e�B�N�����x
 * @param[in] dt ���ԃX�e�b�v��
 * @param[in] nprts �p�[�e�B�N����
 */
__global__
void sphIntegrate(float4* ppos,	float4* pvel, 
				  float4* pacc, float* dens, int* attr, float dt, uint nprts)
{
	uint index = __umul24(blockIdx.x,blockDim.x)+threadIdx.x;
	if(index >= nprts) return;
	if(attr[index] == -1) return;

	float3 x = make_float3(ppos[index]);
	float3 v = make_float3(pvel[index]);
	float3 a = make_float3(pacc[index]);
	//float3 v_old = v;

	// �X�V�ʒu�C���x�̍X�V
	v += dt*a;
	x += dt*v;

	// �ő́E���E�Ƃ̏Փ�
	calCollisionSolid(x, v, dt);

	// �ʒu�Ƒ��x�̍X�V
	ppos[index] = make_float4(x);
	pvel[index] = make_float4(v);
}



/*!
 * �p�[�e�B�N���ʒu�C���x�̍X�V(Leap-Frog)
 * @param[inout] ppos �p�[�e�B�N���ʒu
 * @param[inout] pvel �p�[�e�B�N�����x
 * @param[in] pfrc �p�[�e�B�N���ɂ������
 * @param[in] dens �p�[�e�B�N�����x
 * @param[in] vrts
 * @param[in] dt ���ԃX�e�b�v��
 * @param[in] nprts �p�[�e�B�N����
 */
__global__
void sphIntegrateWithPolygon(float4* ppos, float4* pvel, float4* pacc, float* dens, int* attr, 
							 float3* vrts, int3* tris, int tri_num, float dt, rxParticleCell cell)
{
	uint index = __umul24(blockIdx.x,blockDim.x)+threadIdx.x;
	if(index >= cell.uNumParticles) return;
	if(attr[index] == -1) return;

	float3 x = make_float3(ppos[index]);
	float3 v = make_float3(pvel[index]);
	float3 a = make_float3(pacc[index]);
	//float3 v_old = v;
	float3 x_old = x;

	// �X�V�ʒu�C���x�̍X�V
	v += dt*a;
	x += dt*v;

	// �|���S���I�u�W�F�N�g�Ƃ̏Փ�
	int3 gridPos[2];
	gridPos[0] = calcGridPos(x_old);	// �ʒu�X�V�O�̃p�[�e�B�N����������O���b�h
	gridPos[1] = calcGridPos(x);		// �ʒu�X�V��̃p�[�e�B�N����������O���b�h
	for(int i = 0; i < 2; ++i){
		uint grid_hash = calcGridHash(gridPos[i]);
		uint start_index = cell.dPolyCellStart[grid_hash];
		if(start_index != 0xffffffff){	// �Z������łȂ����̃`�F�b�N

			uint end_index = cell.dPolyCellEnd[grid_hash];
			for(uint j = start_index; j < end_index; ++j){
				uint pidx = cell.dSortedPolyIdx[j];

				int3 idx = tris[pidx];
				if(calCollisionPolygon(x_old, x, v, vrts[idx.x], vrts[idx.y], vrts[idx.z], dt)){
				}
			}
		}
	}

	// �ő́E���E�Ƃ̏Փ�
	calCollisionSolid(x, v, dt);

	// �ʒu�Ƒ��x�̍X�V
	ppos[index] = make_float4(x);
	pvel[index] = make_float4(v);
}

/*!
 * �p�[�e�B�N���ʒu���`�F�b�N���č폜�̈���Ȃ�Α�����-1�ɂ��āC�͈͊O�ɔ�΂�
 * @param[inout] ppos �p�[�e�B�N���ʒu
 * @param[in] dt ���ԃX�e�b�v��
 * @param[in] nprts �p�[�e�B�N����
 */
__global__
void checkDelete(float4* ppos, float4* pvel, int* attr, float3 minp, float3 maxp, float3 far_pos, uint nprts)
{
	uint index = __umul24(blockIdx.x,blockDim.x)+threadIdx.x;
	if(index >= nprts) return;
	if(attr[index] == -1) return;	// ���łɍ폜����Ă����ꍇ�͔�΂�

	float3 x = make_float3(ppos[index]);

	if((x.x > minp.x && x.x < maxp.x) && (x.y > minp.y && x.y < maxp.y) && (x.z > minp.z && x.z < maxp.z)){
		// ������-1�ɂ��āCfar_point�ɔ�΂�
		ppos[index] = make_float4(far_pos);
		pvel[index] = make_float4(0.0);
		attr[index] = -1;
	}
}
/*!
 * �p�[�e�B�N���ʒu���`�F�b�N���č폜�̈���Ȃ�Α�����-1�ɂ��āC�͈͊O�ɔ�΂�
 * @param[inout] ppos �p�[�e�B�N���ʒu
 * @param[in] dt ���ԃX�e�b�v��
 * @param[in] nprts �p�[�e�B�N����
 */
__global__
void checkDeleteX(float4* ppos, float4* pvel, int* attr, float xmax, float3 far_pos, uint nprts)
{
	uint index = __umul24(blockIdx.x,blockDim.x)+threadIdx.x;
	if(index >= nprts) return;
	if(attr[index] == -1) return;	// ���łɍ폜����Ă����ꍇ�͔�΂�

	float3 x = make_float3(ppos[index]);

	if(x.x > xmax){
		// ������-1�ɂ��āCfar_point�ɔ�΂�
		ppos[index] = make_float4(far_pos);
		pvel[index] = make_float4(0.0);
		attr[index] = -1;
	}
}




/*!
 * �^����ꂽ�Z�����̃p�[�e�B�N���Ƃ̋������疧�x���v�Z
 * @param[in] gridPos �O���b�h�ʒu
 * @param[in] pos �v�Z���W
 * @param[in] cell �p�[�e�B�N���O���b�h�f�[�^
 * @return �Z�����̃p�[�e�B�N������v�Z�������x�l
 */
__device__
float calDensityCellG(int3 gridPos, float3 pos0, rxParticleCell cell)
{
	uint gridHash = calcGridHash(gridPos);

	// �Z�����̃p�[�e�B�N���̃X�^�[�g�C���f�b�N�X
	uint startIndex = FETCHC(dCellStart, gridHash);

	float h = params.EffectiveRadius;
	float d = 0.0f;
	if(startIndex != 0xffffffff){	// �Z������łȂ����̃`�F�b�N
		// �Z�����̃p�[�e�B�N���Ŕ���
		uint endIndex = FETCHC(dCellEnd, gridHash);

		for(uint j = startIndex; j < endIndex; ++j){
			//if(j != index){
				float3 pos1 = make_float3(FETCHC(dSortedPos, j));
				//cell.dSortedIndex[j];

				float3 rij = pos0-pos1;
				float r = length(rij);

				if(r <= h){
					float q = h*h-r*r;

					d += params.Mass*params.Wpoly6*q*q*q;
				}

			//}
		}
	}

	return d;
}

/*! �ǉ��F�X�p
 * �^����ꂽ�Z�����̃p�[�e�B�N���Ƃ̋������疧�x���v�Z
 * @param[in] gridPos �O���b�h�ʒu
 * @param[in] pos �v�Z���W
 * @param[in] cell �p�[�e�B�N���O���b�h�f�[�^
 * @return �Z�����̃p�[�e�B�N������v�Z�������x�l
 */
__device__
float calDensityCellGIceMesh(int3 gridPos, float3 pos0, rxParticleCell cell, float* bIceCheck)
{
	uint gridHash = calcGridHash(gridPos);

	// �Z�����̃p�[�e�B�N���̃X�^�[�g�C���f�b�N�X
	uint startIndex = FETCHC(dCellStart, gridHash);

	float h = params.EffectiveRadius;
	float d = 0.0f;
	if(startIndex != 0xffffffff){	// �Z������łȂ����̃`�F�b�N
		// �Z�����̃p�[�e�B�N���Ŕ���
		uint endIndex = FETCHC(dCellEnd, gridHash);

		for(uint j = startIndex; j < endIndex; ++j){
			//if(j != index){
				float3 pos1 = make_float3(FETCHC(dSortedPos, j));
				//cell.dSortedIndex[j];

				float3 rij = pos0-pos1;
				float r = length(rij);

				if(r <= h){
					uint pIndx = cell.dSortedIndex[j];
					if( bIceCheck[pIndx] < 0.0f ) continue;		//�X�łȂ��Ȃ烁�b�V�������Ȃ�
//					if( bIceCheck[j] < 0.0f ) continue;		//�X�łȂ��Ȃ烁�b�V�������Ȃ�
//					printf("pIndx = %d\n", pIndx);
//					printf("j = %d\n", j);
					float q = h*h-r*r;
					d += params.Mass*params.Wpoly6*q*q*q;
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
void sphCalDensityInGrid(float* GridD, rxParticleCell cell, 
					uint3 gnum, float3 gmin, float3 glen)
{
	uint blockId = __mul24(blockIdx.y, gridDim.x)+blockIdx.x;
	uint i = __mul24(blockId, blockDim.x)+threadIdx.x;

	uint3 gridPos = calcGridPosU(i, gnum);

	if(gridPos.x < gnum.x && gridPos.y < gnum.y && gridPos.z < gnum.z){
		float3 gpos;
		gpos.x = gmin.x+(gridPos.x)*glen.x;
		gpos.y = gmin.y+(gridPos.y)*glen.y;
		gpos.z = gmin.z+(gridPos.z)*glen.z;

		float d = 0.0f;

		int3 pgpos = calcGridPos(gpos);

		float h = params.EffectiveRadius;
		int3 grid_pos0, grid_pos1;
		grid_pos0 = calcGridPos(gpos-make_float3(h));
		grid_pos1 = calcGridPos(gpos+make_float3(h));

		for(int z = grid_pos0.z; z <= grid_pos1.z; ++z){
			for(int y = grid_pos0.y; y <= grid_pos1.y; ++y){
				for(int x = grid_pos0.x; x <= grid_pos1.x; ++x){
					int3 neighbourPos = make_int3(x, y, z);

					d += calDensityCellG(neighbourPos, gpos, cell);
				}
			}
		}

		GridD[gridPos.x+gridPos.y*gnum.x+gridPos.z*gnum.x*gnum.y] = d;
	}

}


/*!�ǉ��F�X�p
 * �O���b�h��ł̖��x���v�Z
 * @param[out] GridD �O���b�h���x
 * @param[in] cell �p�[�e�B�N���O���b�h�f�[�^
 * @param[in] gnum �O���b�h��
 * @param[in] gmin �O���b�h�ŏ����W
 * @param[in] glen �O���b�h��
 */
__global__
void sphCalDensityInGridIceMesh(float* GridD, rxParticleCell cell, 
					uint3 gnum, float3 gmin, float3 glen, float* bIceFlag)
{
	uint blockId = __mul24(blockIdx.y, gridDim.x)+blockIdx.x;
	uint i = __mul24(blockId, blockDim.x)+threadIdx.x;

	uint3 gridPos = calcGridPosU(i, gnum);

	if(gridPos.x < gnum.x && gridPos.y < gnum.y && gridPos.z < gnum.z){
		float3 gpos;
		gpos.x = gmin.x+(gridPos.x)*glen.x;
		gpos.y = gmin.y+(gridPos.y)*glen.y;
		gpos.z = gmin.z+(gridPos.z)*glen.z;

		float d = 0.0f;

		int3 pgpos = calcGridPos(gpos);

		float h = params.EffectiveRadius;
		int3 grid_pos0, grid_pos1;
		grid_pos0 = calcGridPos(gpos-make_float3(h));
		grid_pos1 = calcGridPos(gpos+make_float3(h));

		for(int z = grid_pos0.z; z <= grid_pos1.z; ++z){
			for(int y = grid_pos0.y; y <= grid_pos1.y; ++y){
				for(int x = grid_pos0.x; x <= grid_pos1.x; ++x){
					int3 neighbourPos = make_int3(x, y, z);

					d += calDensityCellGIceMesh(neighbourPos, gpos, cell, bIceFlag);
				}
			}
		}

		GridD[gridPos.x+gridPos.y*gnum.x+gridPos.z*gnum.x*gnum.y] = d;
	}

}

#endif // #ifndef _RX_CUSPH_KERNEL_CU_



