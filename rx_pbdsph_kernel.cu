/*! 
  @file rx_sph_kernel.cu
	
  @brief CUDA�ɂ��SPH
 
  @author Makoto Fujisawa
  @date 2009-08, 2011-06
*/
// FILE --rx_sph_kernel.cu--

#ifndef _RX_PBDSPH_KERNEL_CU_
#define _RX_PBDSPH_KERNEL_CU_


//-----------------------------------------------------------------------------
// �C���N���[�h�t�@�C��
//-----------------------------------------------------------------------------
#include "rx_cu_common.cu"


//-----------------------------------------------------------------------------
// PBDSPH
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
float calDensityCellPB(int3 gridPos, uint i, float3 pos0, rxParticleCell cell)
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
void pbdsphCalDensity(float* newDens, rxParticleCell cell)
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
				dens += calDensityCellPB(n_grid_pos, index, pos, cell);
			}
		}
	}

	// ���x�ƈ��͒l�����ʂɏ�������
	uint oIdx = cell.dSortedIndex[index];
	newDens[oIdx] = dens;
}

/*!
 * �^����ꂽ�Z�����̃p�[�e�B�N���Ƃ̋�������͏���v�Z
 * @param[in] gridPos �O���b�h�ʒu
 * @param[in] i �p�[�e�B�N���C���f�b�N�X
 * @param[in] pos0 �v�Z���W
 * @param[in] vel0 �v�Z���W�̑��x
 * @param[in] dens0 �v�Z���W�̖��x
 * @param[in] dens �p�[�e�B�N�����x
 * @param[in] cell �p�[�e�B�N���O���b�h�f�[�^
 * @return �Z�����̃p�[�e�B�N������v�Z�����͏�
 */
__device__
float3 calExtForceCell(int3 gridPos, uint i, float3 pos0, float3 vel0, float dens0, float* dens, rxParticleCell cell)
{
	uint gridHash = calcGridHash(gridPos);

	// �Z�����̃p�[�e�B�N���̃X�^�[�g�C���f�b�N�X
	uint startIndex = FETCHC(dCellStart, gridHash);

	float h = params.EffectiveRadius;

	float3 frc = make_float3(0.0f);
	if(startIndex != 0xffffffff){	// �Z������łȂ����̃`�F�b�N
		// �Z�����̃p�[�e�B�N���Ŕ���
		uint endIndex = FETCHC(dCellEnd, gridHash);
		for(uint j = startIndex; j < endIndex; ++j){
			if(j != i){
				// �ߖT�p�[�e�B�N���̃p�����[�^
				float3 pos1 = make_float3(FETCHC(dSortedPos, j));
				float3 vel1 = make_float3(FETCHC(dSortedVel, j));

				float3 rij = pos0-pos1;
				float r = length(rij);

				if(r <= h && r > 0.0001){
					float dens1 = dens[cell.dSortedIndex[j]];

					float3 vij = vel1-vel0;

					float q = h-r;

					// �S����
					frc += params.Viscosity*params.Mass*(vij/dens1)*params.LWvisc*q;
				}
			}
		}
	}

	return frc;
}

/*!
 * �p�[�e�B�N���ɂ�����O�͂̌v�Z(�J�[�l���֐�)
 * @param[in] dens �p�[�e�B�N�����x
 * @param[out] outFrc �p�[�e�B�N���ɂ������
 * @param[in] cell �p�[�e�B�N���O���b�h�f�[�^
 */
__global__
void pbdsphCalExternalForces(float* dens, float4* outFrc, rxParticleCell cell)
{
	uint index = __mul24(blockIdx.x,blockDim.x)+threadIdx.x;
	if(index >= cell.uNumParticles) return;	
	
	// �\�[�g�ςݔz�񂩂�p�[�e�B�N���f�[�^���擾
	float3 pos0 = make_float3(FETCHC(dSortedPos, index));
	float3 vel0 = make_float3(FETCHC(dSortedVel, index));
	float h = params.EffectiveRadius;

	// �p�[�e�B�N���̃\�[�g�Ȃ��z���ł̃C���f�b�N�X
	uint oIdx = cell.dSortedIndex[index];

	float3 frc = make_float3(0.0f);
	float dens0 = dens[oIdx];

	int3 grid_pos0, grid_pos1;
	grid_pos0 = calcGridPos(pos0-make_float3(h));
	grid_pos1 = calcGridPos(pos0+make_float3(h));

	// ���͂̃O���b�h���܂߂ċߖT�T���C���͍��C�S�������v�Z
	for(int z = grid_pos0.z; z <= grid_pos1.z; ++z){
		for(int y = grid_pos0.y; y <= grid_pos1.y; ++y){
			for(int x = grid_pos0.x; x <= grid_pos1.x; ++x){
				int3 n_grid_pos = make_int3(x, y, z);

				frc += calExtForceCell(n_grid_pos, index, pos0, vel0, dens0, dens, cell);
			}
		}
	}

	// �O��(�d��)
	frc += params.Gravity;

	outFrc[oIdx] = make_float4(frc, 0.0f);
}


/*!
 * �^����ꂽ�Z�����̃p�[�e�B�N���Ƃ̋�������X�P�[�����O�t�@�N�^�̕��ꍀ�v�Z
 * @param[in] gridPos �O���b�h�ʒu
 * @param[in] index �p�[�e�B�N���C���f�b�N�X
 * @param[in] pos �v�Z���W
 * @param[in] cell �p�[�e�B�N���O���b�h�f�[�^
 * @return �Z�����̃p�[�e�B�N������v�Z�������x�l
 */
__device__
float calScalingFactorCell(int3 gridPos, uint i, float3 pos0, rxParticleCell cell)
{
	uint gridHash = calcGridHash(gridPos);

	// �Z�����̃p�[�e�B�N���̃X�^�[�g�C���f�b�N�X
	uint startIndex = FETCHC(dCellStart, gridHash);

	float h = params.EffectiveRadius;
	float r0 = params.Density;
	float sd = 0.0f;
	if(startIndex != 0xffffffff){	// �Z������łȂ����̃`�F�b�N
		// �Z�����̃p�[�e�B�N���Ŕ���
		uint endIndex = FETCHC(dCellEnd, gridHash);
		for(uint j = startIndex; j < endIndex; ++j){
			if(j == i) continue;

			float3 pos1 = make_float3(FETCHC(dSortedPos, j));

			float3 rij = pos0-pos1;
			float r = length(rij);

			if(r <= h && r > 0.0){
				float q = h-r;

				// Spiky�J�[�l���ňʒu�ϓ����v�Z
				float3 dp = (params.GWspiky*q*q*rij/r)/r0;

				sd += dot(dp, dp);
			}

		}
	}

	return sd;
}

/*!
 * �X�P�[�����O�t�@�N�^�̌v�Z
 * @param[in] ppos �p�[�e�B�N�����S���W
 * @param[out] pdens �p�[�e�B�N�����x
 * @param[out] pscl �X�P�[�����O�t�@�N�^
 * @param[in] eps �ɘa�W��
 * @param[in]  cell �p�[�e�B�N���O���b�h�f�[�^
 */
__global__
void pbdsphCalScalingFactor(float4* ppos, float* pdens, float* pscl, float eps, rxParticleCell cell)
{
	// �p�[�e�B�N���C���f�b�N�X
	uint index = __mul24(blockIdx.x,blockDim.x)+threadIdx.x;
	if(index >= cell.uNumParticles) return;	
	
	float3 pos = make_float3(FETCHC(dSortedPos, index));	// �p�[�e�B�N���ʒu
	//int3 grid_pos = calcGridPos(pos);	// �p�[�e�B�N����������O���b�h�ʒu

	float h = params.EffectiveRadius;
	float r0 = params.Density;

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
				dens += calDensityCellPB(n_grid_pos, index, pos, cell);
			}
		}
	}

	// ���x�S������(��(1))
	float C = dens/r0-1.0;

	// ���͂̃O���b�h���܂߂ċߖT�T���C�X�P�[�����O�t�@�N�^�̕��ꍀ�v�Z
	float sd = 0.0f;
	for(int z = grid_pos0.z; z <= grid_pos1.z; ++z){
		for(int y = grid_pos0.y; y <= grid_pos1.y; ++y){
			for(int x = grid_pos0.x; x <= grid_pos1.x; ++x){
				int3 n_grid_pos = make_int3(x, y, z);
				sd += calScalingFactorCell(n_grid_pos, index, pos, cell);
			}
		}
	}

	// �p�[�e�B�N���̃\�[�g�Ȃ��z���ł̃C���f�b�N�X
	uint oIdx = cell.dSortedIndex[index];

	// �X�P�[�����O�t�@�N�^�̌v�Z(��(11))
	pscl[oIdx] = -C/(sd+eps);

	// �X�V���ꂽ���x
	pdens[oIdx] = dens;
}


/*!
 * �^����ꂽ�Z�����̃p�[�e�B�N���Ƃ̋�������X�P�[�����O�t�@�N�^�̕��ꍀ�v�Z
 * @param[in] gridPos �O���b�h�ʒu
 * @param[in] index �p�[�e�B�N���C���f�b�N�X
 * @param[in] pos �v�Z���W
 * @param[in] cell �p�[�e�B�N���O���b�h�f�[�^
 * @return �Z�����̃p�[�e�B�N������v�Z�������x�l
 */
__device__
float3 calPositionCorrectionCell(int3 gridPos, uint i, float3 pos0, float* pscl, rxParticleCell cell)
{
	uint gridHash = calcGridHash(gridPos);

	// �Z�����̃p�[�e�B�N���̃X�^�[�g�C���f�b�N�X
	uint startIndex = FETCHC(dCellStart, gridHash);

	float k = params.AP_K;
	float n = params.AP_N;
	float wq = params.AP_WQ;

	float h = params.EffectiveRadius;
	float r0 = params.Density;
	float3 dp = make_float3(0.0);

	float dt = params.Dt;

	float si = pscl[cell.dSortedIndex[i]];

	if(startIndex != 0xffffffff){	// �Z������łȂ����̃`�F�b�N
		// �Z�����̃p�[�e�B�N���Ŕ���
		uint endIndex = FETCHC(dCellEnd, gridHash);
		for(uint j = startIndex; j < endIndex; ++j){
			if(j == i) continue;

			float3 pos1 = make_float3(FETCHC(dSortedPos, j));

			float3 rij = pos0-pos1;
			float r = length(rij);

			if(r <= h && r > 0.0){
				float scorr = 0.0f;

				if(params.AP){
					float q1 = h*h-r*r;
					float ww = params.Wpoly6*q1*q1*q1/wq;
					scorr = -k*pow(ww, n)*dt*dt;
				}
				float q = h-r;
				float sj = pscl[cell.dSortedIndex[j]];

				// Spiky�J�[�l���ňʒu�C���ʂ��v�Z
				dp += (si+sj+scorr)*(params.GWspiky*q*q*rij/r)/r0;
			}

		}
	}

	return dp;
}

/*!
 * �X�P�[�����O�t�@�N�^�̌v�Z
 * @param[in] ppos �p�[�e�B�N�����S���W
 * @param[out] pdens �p�[�e�B�N�����x
 * @param[out] pscl �X�P�[�����O�t�@�N�^
 * @param[in] eps �ɘa�W��
 * @param[in]  cell �p�[�e�B�N���O���b�h�f�[�^
 */
__global__
void pbdsphPositionCorrection(float4* ppos, float* pscl, float4* pdp, rxParticleCell cell)
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

	// ���͂̃O���b�h���܂߂ċߖT�T���C�ʒu�C���ʂ��v�Z
	float3 dpij = make_float3(0.0f);
	for(int z = grid_pos0.z; z <= grid_pos1.z; ++z){
		for(int y = grid_pos0.y; y <= grid_pos1.y; ++y){
			for(int x = grid_pos0.x; x <= grid_pos1.x; ++x){
				int3 n_grid_pos = make_int3(x, y, z);
				dpij += calPositionCorrectionCell(n_grid_pos, index, pos, pscl, cell);
			}
		}
	}

	// �p�[�e�B�N���̃\�[�g�Ȃ��z���ł̃C���f�b�N�X
	uint oIdx = cell.dSortedIndex[index];

	// �ʒu�C����
	pdp[oIdx] = make_float4(dpij, 0.0);
}

/*!
 * �p�[�e�B�N���ʒu�C��
 * @param[inout] pos �p�[�e�B�N���ʒu
 * @param[in] pdp �ʒu�C����
 * @param[in] nprts �p�[�e�B�N����
 */
__global__
void pbdsphCorrectPosition(float4* ppos, float4* pdp, uint nprts)
{
	uint index = __umul24(blockIdx.x,blockDim.x)+threadIdx.x;
	if(index >= nprts) return;

	// �ʒu�C��
	ppos[index] += pdp[index];
}

/*!
 * ���x�ϓ��̌v�Z
 * @param[inout] pos �p�[�e�B�N���ʒu
 * @param[in] pdp �ʒu�C����
 * @param[in] nprts �p�[�e�B�N����
 */
__global__
void pbdsphDensityFluctuation(float* perr, float* pdens, float rest_dens, uint nprts)
{
	uint index = __umul24(blockIdx.x,blockDim.x)+threadIdx.x;
	if(index >= nprts) return;

	// ���x�ϓ�
	//perr[index] = fabs(pdens[index]-rest_dens)/rest_dens;
	float err = pdens[index]-rest_dens;
	perr[index] = (err >= 0.0f ? err : 0.0f)/rest_dens;
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
float calBoundaryDensityCellPB(int3 gridPos, uint i, float3 pos0, float* dVolB, rxParticleCell bcell)
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
void pbdsphCalBoundaryDensity(float* newDens, float4* dPos, float* dVolB, rxParticleCell bcell, uint pnum)
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
				dens += calBoundaryDensityCellPB(n_grid_pos, index, pos, dVolB, bcell);
			}
		}
	}

	// ���x�����ʂɏ�������
	newDens[index] += dens;
}




/*!
 * �^����ꂽ�Z�����̃p�[�e�B�N���Ƃ̋�������X�P�[�����O�t�@�N�^�̕��ꍀ�v�Z
 * @param[in] gridPos �O���b�h�ʒu
 * @param[in] index �p�[�e�B�N���C���f�b�N�X
 * @param[in] pos �v�Z���W
 * @param[in] cell �p�[�e�B�N���O���b�h�f�[�^
 * @return �Z�����̃p�[�e�B�N������v�Z�������x�l
 */
__device__
float calBoundaryScalingFactorCell(int3 gridPos, uint i, float3 pos0, float* dVolB, rxParticleCell bcell)
{
	uint gridHash = calcGridHashB(gridPos, params.GridSizeB);

	// �Z�����̃p�[�e�B�N���̃X�^�[�g�C���f�b�N�X
	uint startIndex = bcell.dCellStart[gridHash];

	float h = params.EffectiveRadius;
	float r0 = params.Density;
	float sd = 0.0f;
	if(startIndex != 0xffffffff){	// �Z������łȂ����̃`�F�b�N
		// �Z�����̃p�[�e�B�N���Ŕ���
		uint endIndex = bcell.dCellEnd[gridHash];
		for(uint j = startIndex; j < endIndex; ++j){
			float3 pos1 = make_float3(bcell.dSortedPos[j]);
			uint jdx = bcell.dSortedIndex[j];

			float3 rij = pos0-pos1;
			float r = length(rij);

			if(r <= h && r > 0.0){
				float q = h-r;

				// Spiky�J�[�l���ňʒu�ϓ����v�Z
				float3 dp = (params.Density*dVolB[jdx]/params.Mass)*(params.GWspiky*q*q*rij/r)/r0;

				sd += dot(dp, dp);
			}

		}
	}

	return sd;
}

/*!
 * �X�P�[�����O�t�@�N�^�̌v�Z(���E�p�[�e�B�N���܂�)
 * @param[in] ppos �p�[�e�B�N�����S���W
 * @param[out] pdens �p�[�e�B�N�����x
 * @param[out] pscl �X�P�[�����O�t�@�N�^
 * @param[in] eps �ɘa�W��
 * @param[in] cell �p�[�e�B�N���O���b�h�f�[�^
 */
__global__
void pbdsphCalScalingFactorWithBoundary(float4* ppos, float* pdens, float* pscl, float eps, rxParticleCell cell, 
										float* bvol, rxParticleCell bcell)
{
	// �p�[�e�B�N���C���f�b�N�X
	uint index = __mul24(blockIdx.x,blockDim.x)+threadIdx.x;
	if(index >= cell.uNumParticles) return;	
	
	float3 pos = make_float3(FETCHC(dSortedPos, index));	// �p�[�e�B�N���ʒu
	//int3 grid_pos = calcGridPos(pos);	// �p�[�e�B�N����������O���b�h�ʒu

	float h = params.EffectiveRadius;
	float r0 = params.Density;

	// �p�[�e�B�N�����͂̃O���b�h
	int3 grid_pos0, grid_pos1;
	grid_pos0 = calcGridPos(pos-make_float3(h));
	grid_pos1 = calcGridPos(pos+make_float3(h));

	// ���̃p�[�e�B�N���ɂ�閧�x
	float dens = 0.0f;
	for(int z = grid_pos0.z; z <= grid_pos1.z; ++z){
		for(int y = grid_pos0.y; y <= grid_pos1.y; ++y){
			for(int x = grid_pos0.x; x <= grid_pos1.x; ++x){
				int3 n_grid_pos = make_int3(x, y, z);
				dens += calDensityCellPB(n_grid_pos, index, pos, cell);
			}
		}
	}

	// �p�[�e�B�N�����͂̃O���b�h
	int3 grid_pos2, grid_pos3;
	grid_pos2 = calcGridPosB(pos-make_float3(h), params.WorldOriginB, params.CellWidthB, params.GridSizeB);
	grid_pos3 = calcGridPosB(pos+make_float3(h), params.WorldOriginB, params.CellWidthB, params.GridSizeB);

	// ���E�p�[�e�B�N���ɂ�閧�x
	for(int z = grid_pos2.z; z <= grid_pos3.z; ++z){
		for(int y = grid_pos2.y; y <= grid_pos3.y; ++y){
			for(int x = grid_pos2.x; x <= grid_pos3.x; ++x){
				int3 n_grid_pos = make_int3(x, y, z);
				dens += calBoundaryDensityCellPB(n_grid_pos, index, pos, bvol, bcell);
			}
		}
	}

	// ���x�S������(��(1))
	float C = dens/r0-1.0;

	// ���̃p�[�e�B�N���ɂ��X�P�[�����O�t�@�N�^�̕��ꍀ�v�Z
	float sd = 0.0f;
	for(int z = grid_pos0.z; z <= grid_pos1.z; ++z){
		for(int y = grid_pos0.y; y <= grid_pos1.y; ++y){
			for(int x = grid_pos0.x; x <= grid_pos1.x; ++x){
				int3 n_grid_pos = make_int3(x, y, z);
				sd += calScalingFactorCell(n_grid_pos, index, pos, cell);
			}
		}
	}

	// ���E�p�[�e�B�N���ɂ��X�P�[�����O�t�@�N�^�̕��ꍀ�v�Z
	for(int z = grid_pos2.z; z <= grid_pos3.z; ++z){
		for(int y = grid_pos2.y; y <= grid_pos3.y; ++y){
			for(int x = grid_pos2.x; x <= grid_pos3.x; ++x){
				int3 n_grid_pos = make_int3(x, y, z);
				sd += calBoundaryScalingFactorCell(n_grid_pos, index, pos, bvol, bcell);
			}
		}
	}

	// �p�[�e�B�N���̃\�[�g�Ȃ��z���ł̃C���f�b�N�X
	uint oIdx = cell.dSortedIndex[index];

	// �X�P�[�����O�t�@�N�^�̌v�Z(��(11))
	pscl[oIdx] = -C/(sd+eps);

	// �X�V���ꂽ���x
	pdens[oIdx] = dens;
}



/*!
 * �X�P�[�����O�t�@�N�^�̌v�Z(���E�p�[�e�B�N���܂�)
 * @param[in] ppos �p�[�e�B�N�����S���W
 * @param[out] pdens �p�[�e�B�N�����x
 * @param[out] pscl �X�P�[�����O�t�@�N�^
 * @param[in] eps �ɘa�W��
 * @param[in] cell �p�[�e�B�N���O���b�h�f�[�^
 */
__global__
void pbdsphCalBoundaryScalingFactor(float4* ppos, float* pdens, float eps, rxParticleCell cell, 
									float* bvol, float* bscl, rxParticleCell bcell)
{
	// �p�[�e�B�N���C���f�b�N�X
	uint index = __mul24(blockIdx.x,blockDim.x)+threadIdx.x;
	if(index >= bcell.uNumParticles) return;	
	
	float3 pos = make_float3(bcell.dSortedPos[index]);	// �p�[�e�B�N���ʒu

	float h = params.EffectiveRadius;
	float r0 = params.Density;

	// �p�[�e�B�N�����͂̃O���b�h(���̃p�[�e�B�N���p)
	int3 grid_pos0, grid_pos1;
	grid_pos0 = calcGridPos(pos-make_float3(h));
	grid_pos1 = calcGridPos(pos+make_float3(h));

	// ���̃p�[�e�B�N���ɂ�閧�x
	float dens = 0.0f;
	for(int z = grid_pos0.z; z <= grid_pos1.z; ++z){
		for(int y = grid_pos0.y; y <= grid_pos1.y; ++y){
			for(int x = grid_pos0.x; x <= grid_pos1.x; ++x){
				int3 n_grid_pos = make_int3(x, y, z);
				dens += calDensityCellPB(n_grid_pos, index, pos, cell);
			}
		}
	}

	// �p�[�e�B�N�����͂̃O���b�h(���E�p�[�e�B�N���p)
	int3 grid_pos2, grid_pos3;
	grid_pos2 = calcGridPosB(pos-make_float3(h), params.WorldOriginB, params.CellWidthB, params.GridSizeB);
	grid_pos3 = calcGridPosB(pos+make_float3(h), params.WorldOriginB, params.CellWidthB, params.GridSizeB);

	// ���E�p�[�e�B�N���ɂ�閧�x
	for(int z = grid_pos2.z; z <= grid_pos3.z; ++z){
		for(int y = grid_pos2.y; y <= grid_pos3.y; ++y){
			for(int x = grid_pos2.x; x <= grid_pos3.x; ++x){
				int3 n_grid_pos = make_int3(x, y, z);
				dens += calBoundaryDensityCellPB(n_grid_pos, index, pos, bvol, bcell);
			}
		}
	}

	// ���x�S������(��(1))
	float C = dens/r0-1.0;

	// ���̃p�[�e�B�N���ɂ��X�P�[�����O�t�@�N�^�̕��ꍀ�v�Z
	float sd = 0.0f;
	for(int z = grid_pos0.z; z <= grid_pos1.z; ++z){
		for(int y = grid_pos0.y; y <= grid_pos1.y; ++y){
			for(int x = grid_pos0.x; x <= grid_pos1.x; ++x){
				int3 n_grid_pos = make_int3(x, y, z);
				sd += calScalingFactorCell(n_grid_pos, index, pos, cell);
			}
		}
	}

	// ���E�p�[�e�B�N���ɂ��X�P�[�����O�t�@�N�^�̕��ꍀ�v�Z
	for(int z = grid_pos2.z; z <= grid_pos3.z; ++z){
		for(int y = grid_pos2.y; y <= grid_pos3.y; ++y){
			for(int x = grid_pos2.x; x <= grid_pos3.x; ++x){
				int3 n_grid_pos = make_int3(x, y, z);
				sd += calBoundaryScalingFactorCell(n_grid_pos, index, pos, bvol, bcell);
			}
		}
	}

	// �p�[�e�B�N���̃\�[�g�Ȃ��z���ł̃C���f�b�N�X
	uint oIdx = bcell.dSortedIndex[index];

	// �X�P�[�����O�t�@�N�^�̌v�Z(��(11))
	bscl[oIdx] = -C/(sd+eps);
}



/*!
 * �^����ꂽ�Z�����̃p�[�e�B�N���Ƃ̋�������X�P�[�����O�t�@�N�^�̕��ꍀ�v�Z
 * @param[in] gridPos �O���b�h�ʒu
 * @param[in] index �p�[�e�B�N���C���f�b�N�X
 * @param[in] pos �v�Z���W
 * @param[in] cell �p�[�e�B�N���O���b�h�f�[�^
 * @return �Z�����̃p�[�e�B�N������v�Z�������x�l
 */
__device__
float3 calBoundaryPositionCorrectionCell(int3 gridPos, uint i, float3 pos0, float si, float* bscl, float* bvol, rxParticleCell bcell)
{
	uint gridHash = calcGridHashB(gridPos, params.GridSizeB);

	// �Z�����̃p�[�e�B�N���̃X�^�[�g�C���f�b�N�X
	uint startIndex = bcell.dCellStart[gridHash];

	float k = params.AP_K;
	float n = params.AP_N;
	float wq = params.AP_WQ;

	float h = params.EffectiveRadius;
	float r0 = params.Density;
	float3 dp = make_float3(0.0);

	float dt = params.Dt;

	if(startIndex != 0xffffffff){	// �Z������łȂ����̃`�F�b�N
		// �Z�����̃p�[�e�B�N���Ŕ���
		uint endIndex = bcell.dCellEnd[gridHash];
		for(uint j = startIndex; j < endIndex; ++j){
			float3 pos1 = make_float3(bcell.dSortedPos[j]);
			uint jdx = bcell.dSortedIndex[j];

			float3 rij = pos0-pos1;
			float r = length(rij);

			if(r <= h && r > 0.0){
				float scorr = 0.0f;

				if(params.AP){
					float q1 = h*h-r*r;
					float ww = (params.Density*bvol[jdx]/params.Mass)*params.Wpoly6*q1*q1*q1/wq;
					scorr = -k*pow(ww, n)*dt*dt;
				}
				float q = h-r;
				float sj = bscl[jdx];

				// Spiky�J�[�l���ňʒu�C���ʂ��v�Z
				dp += (si+sj+scorr)*(params.GWspiky*q*q*rij/r)/r0;
			}

		}
	}

	return dp;
}

/*!
 * �X�P�[�����O�t�@�N�^�̌v�Z
 * @param[in] ppos �p�[�e�B�N�����S���W
 * @param[out] pdens �p�[�e�B�N�����x
 * @param[out] pscl �X�P�[�����O�t�@�N�^
 * @param[in] eps �ɘa�W��
 * @param[in]  cell �p�[�e�B�N���O���b�h�f�[�^
 */
__global__
void pbdsphPositionCorrectionWithBoundary(float4* ppos, float* pscl, float4* pdp, rxParticleCell cell, 
										  float* bvol, float* bscl, rxParticleCell bcell)
{
	// �p�[�e�B�N���C���f�b�N�X
	uint index = __mul24(blockIdx.x,blockDim.x)+threadIdx.x;
	if(index >= cell.uNumParticles) return;	
	
	float3 pos = make_float3(FETCHC(dSortedPos, index));	// �p�[�e�B�N���ʒu
	//int3 grid_pos = calcGridPos(pos);	// �p�[�e�B�N����������O���b�h�ʒu

	float h = params.EffectiveRadius;

	float si = pscl[cell.dSortedIndex[index]];


	// �p�[�e�B�N�����͂̃O���b�h(���̃p�[�e�B�N���p)
	int3 grid_pos0, grid_pos1;
	grid_pos0 = calcGridPos(pos-make_float3(h));
	grid_pos1 = calcGridPos(pos+make_float3(h));

	// ���̃p�[�e�B�N���ɂ��ʒu�C���ʂ��v�Z
	float3 dpij = make_float3(0.0f);
	for(int z = grid_pos0.z; z <= grid_pos1.z; ++z){
		for(int y = grid_pos0.y; y <= grid_pos1.y; ++y){
			for(int x = grid_pos0.x; x <= grid_pos1.x; ++x){
				int3 n_grid_pos = make_int3(x, y, z);
				dpij += calPositionCorrectionCell(n_grid_pos, index, pos, pscl, cell);
			}
		}
	}

	// �p�[�e�B�N�����͂̃O���b�h(���E�p�[�e�B�N���p)
	int3 grid_pos2, grid_pos3;
	grid_pos2 = calcGridPosB(pos-make_float3(h), params.WorldOriginB, params.CellWidthB, params.GridSizeB);
	grid_pos3 = calcGridPosB(pos+make_float3(h), params.WorldOriginB, params.CellWidthB, params.GridSizeB);

	// ���E�p�[�e�B�N���ɂ��ʒu�C���ʂ��v�Z
	for(int z = grid_pos2.z; z <= grid_pos3.z; ++z){
		for(int y = grid_pos2.y; y <= grid_pos3.y; ++y){
			for(int x = grid_pos2.x; x <= grid_pos3.x; ++x){
				int3 n_grid_pos = make_int3(x, y, z);
				dpij += calBoundaryPositionCorrectionCell(n_grid_pos, index, pos, si, bscl, bvol, bcell);
			}
		}
	}

	// �p�[�e�B�N���̃\�[�g�Ȃ��z���ł̃C���f�b�N�X
	uint oIdx = cell.dSortedIndex[index];

	// �ʒu�C����
	pdp[oIdx] = make_float4(dpij, 0.0);
}


__device__
void calCollisionSolidPB(float3 &pos, float3 &vel, float dt)
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
inline bool calCollisionPolygonPB(float3 &pos0, float3 &pos1, float3 &vel, float3 v0, float3 v1, float3 v2, float dt)
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
void pbdsphIntegrate(float4* ppos, float4* pvel, float4* pacc, int* attr, 
					 float4* new_ppos, float4* new_pvel, float dt, uint nprts)
{
	uint index = __umul24(blockIdx.x,blockDim.x)+threadIdx.x;
	if(index >= nprts) return;
	//if(attr[index] == -1) return;

	float3 x = make_float3(ppos[index]);
	float3 v = make_float3(pvel[index]);
	float3 a = make_float3(pacc[index]);
	//float3 v_old = v;

	// �X�V�ʒu�C���x�̍X�V
	v += dt*a;
	x += dt*v;

	// �ő́E���E�Ƃ̏Փ�
	calCollisionSolidPB(x, v, dt);

	// �ʒu�Ƒ��x�̍X�V
	new_ppos[index] = make_float4(x);
	new_pvel[index] = make_float4(v);
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
void pbdsphIntegrateWithPolygon(float4* ppos, float4* pvel, float4* pacc, int* attr, 
								float4* new_ppos, float4* new_pvel, 
								float3* vrts, int3* tris, int tri_num, float dt, rxParticleCell cell)
{
	uint index = __umul24(blockIdx.x,blockDim.x)+threadIdx.x;
	if(index >= cell.uNumParticles) return;
	//if(attr[index] == -1) return;

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
				if(calCollisionPolygonPB(x_old, x, v, vrts[idx.x], vrts[idx.y], vrts[idx.z], dt)){
				}
			}
		}
	}

	// �ő́E���E�Ƃ̏Փ�
	calCollisionSolidPB(x, v, dt);

	// �ʒu�Ƒ��x�̍X�V
	new_ppos[index] = make_float4(x);
	new_pvel[index] = make_float4(v);
}




/*!
 * �p�[�e�B�N���ʒu�C���x�̍X�V
 * @param[in] ppos �X�V���ꂽ�p�[�e�B�N���ʒu
 * @param[inout] new_ppos �X�e�b�v�ŏ��̃p�[�e�B�N���ʒu/�V�����p�[�e�B�N�����x
 * @param[out] new_pvel �V�����p�[�e�B�N�����x
 * @param[in] dt ���ԃX�e�b�v��
 * @param[in] nprts �p�[�e�B�N����
 */
__global__
void pbdsphUpdatePosition(float4* ppos, float4* new_ppos, float4* new_pvel, float dt, uint nprts)
{
	uint index = __umul24(blockIdx.x,blockDim.x)+threadIdx.x;
	if(index >= nprts) return;

	float3 x0 = make_float3(new_ppos[index]);
	float3 x1 = make_float3(ppos[index]);
	float3 v = (x1-x0)/dt;

	// �ʒu�Ƒ��x�̍X�V
	new_pvel[index] = make_float4(v);
	new_ppos[index] = make_float4(x1);
}

/*!
 * �p�[�e�B�N�����x�̍X�V
 * @param[in] ppos �X�V���ꂽ�p�[�e�B�N���ʒu
 * @param[in] new_ppos �X�e�b�v�ŏ��̃p�[�e�B�N���ʒu/�V�����p�[�e�B�N�����x
 * @param[out] new_pvel �V�����p�[�e�B�N�����x
 * @param[in] dt ���ԃX�e�b�v��
 * @param[in] nprts �p�[�e�B�N����
 */
__global__
void pbdsphUpdateVelocity(float4* ppos, float4* new_ppos, float4* new_pvel, float dt, uint nprts)
{
	uint index = __umul24(blockIdx.x,blockDim.x)+threadIdx.x;
	if(index >= nprts) return;

	float3 x0 = make_float3(new_ppos[index]);
	float3 x1 = make_float3(ppos[index]);
	float3 v = (x1-x0)/dt;

	// �ʒu�Ƒ��x�̍X�V
	new_pvel[index] = make_float4(v);
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
float3 calXsphViscosityCell(int3 gridPos, uint i, float3 pos0, float3 vel0, float4* pvel, float* dens, rxParticleCell cell)
{
	uint gridHash = calcGridHash(gridPos);

	// �Z�����̃p�[�e�B�N���̃X�^�[�g�C���f�b�N�X
	uint startIndex = FETCHC(dCellStart, gridHash);

	float h = params.EffectiveRadius;
	float3 v = make_float3(0.0);
	if(startIndex != 0xffffffff){	// �Z������łȂ����̃`�F�b�N
		// �Z�����̃p�[�e�B�N���Ŕ���
		uint endIndex = FETCHC(dCellEnd, gridHash);
		for(uint j = startIndex; j < endIndex; ++j){
			//if(j == i) continue;

			float3 pos1 = make_float3(FETCHC(dSortedPos, j));

			float3 rij = pos0-pos1;
			float r = length(rij);

			if(r <= h){
				float3 vel1 = make_float3(pvel[cell.dSortedIndex[j]]);
				float3 rho1 = make_float3(dens[cell.dSortedIndex[j]]);

				float q = h*h-r*r;
				v += (params.Mass/rho1)*(vel1-vel0)*params.Wpoly6*q*q*q;
			}
		}
	}

	return v;
}

/*!
 * �p�[�e�B�N�����x�v�Z(�J�[�l���֐�)
 * @param[out] newDens �p�[�e�B�N�����x
 * @param[out] newPres �p�[�e�B�N������
 * @param[in]  cell �p�[�e�B�N���O���b�h�f�[�^
 */
__global__
void xsphVisocosity(float4* ppos, float4* pvel, float4* new_pvel, float* dens, float c, rxParticleCell cell)
{
	// �p�[�e�B�N���C���f�b�N�X
	uint index = __mul24(blockIdx.x,blockDim.x)+threadIdx.x;
	if(index >= cell.uNumParticles) return;	
	
	float3 pos0 = make_float3(FETCHC(dSortedPos, index));	// �p�[�e�B�N���ʒu
	float3 vel0 = make_float3(pvel[cell.dSortedIndex[index]]);	// �p�[�e�B�N�����x
	//int3 grid_pos = calcGridPos(pos0);	// �p�[�e�B�N����������O���b�h�ʒu
	float h = params.EffectiveRadius;

	// �p�[�e�B�N�����͂̃O���b�h
	int3 grid_pos0, grid_pos1;
	grid_pos0 = calcGridPos(pos0-make_float3(h));
	grid_pos1 = calcGridPos(pos0+make_float3(h));

	// ���͂̃O���b�h���܂߂ċߖT�T���C���x�v�Z
	float3 v = make_float3(0.0);
	for(int z = grid_pos0.z; z <= grid_pos1.z; ++z){
		for(int y = grid_pos0.y; y <= grid_pos1.y; ++y){
			for(int x = grid_pos0.x; x <= grid_pos1.x; ++x){
				int3 n_grid_pos = make_int3(x, y, z);
				v += calXsphViscosityCell(n_grid_pos, index, pos0, vel0, pvel, dens, cell);
			}
		}
	}

	// ���x�ƈ��͒l�����ʂɏ�������
	uint oIdx = cell.dSortedIndex[index];
	new_pvel[oIdx] = make_float4(vel0+c*v);
	//new_pvel[oIdx] = make_float4(vel0);
}


/*!
 * �p�[�e�B�N���ʒu���`�F�b�N���č폜�̈���Ȃ�Α�����-1�ɂ��āC�͈͊O�ɔ�΂�
 * @param[inout] ppos �p�[�e�B�N���ʒu
 * @param[in] dt ���ԃX�e�b�v��
 * @param[in] nprts �p�[�e�B�N����
 */
__global__
void checkDeletePB(float4* ppos, float4* pvel, int* attr, float3 minp, float3 maxp, float3 far_pos, uint nprts)
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
void checkDeleteXPB(float4* ppos, float4* pvel, int* attr, float xmax, float3 far_pos, uint nprts)
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
float calDensityCellGPB(int3 gridPos, float3 pos0, rxParticleCell cell)
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

/*!
 * �O���b�h��ł̖��x���v�Z
 * @param[out] GridD �O���b�h���x
 * @param[in] cell �p�[�e�B�N���O���b�h�f�[�^
 * @param[in] gnum �O���b�h��
 * @param[in] gmin �O���b�h�ŏ����W
 * @param[in] glen �O���b�h��
 */
__global__
void pbdsphCalDensityInGrid(float* GridD, rxParticleCell cell, 
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

					d += calDensityCellGPB(neighbourPos, gpos, cell);
				}
			}
		}

		GridD[gridPos.x+gridPos.y*gnum.x+gridPos.z*gnum.x*gnum.y] = d;
	}

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
float3 calNormalCellPB(int3 gridPos, uint i, float3 pos0, float* dens, rxParticleCell cell)
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
void pbdsphCalNormal(float4* newNrms, float* dens, rxParticleCell cell)
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
				nrm += calNormalCellPB(n_grid_pos, index, pos, dens, cell);
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





#endif // #ifndef _RX_PBDSPH_KERNEL_CU_



