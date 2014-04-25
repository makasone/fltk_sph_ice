/*! 
  @file rx_cu_funcs.cuh
	
  @brief CUDA関数の宣言
 
  @author Makoto Fujisawa
  @date 2009-08, 2011-06
*/
// FILE --rx_cu_funcs.cuh--

#ifndef _RX_CU_FUNCS_CUH_
#define _RX_CU_FUNCS_CUH_


//-----------------------------------------------------------------------------
// インクルードファイル
//-----------------------------------------------------------------------------
#include "rx_cu_common.cuh"


//-----------------------------------------------------------------------------
// CUDA関数
//-----------------------------------------------------------------------------
extern "C"
{
void CuInit(int argc, char **argv);
void CuSetDevice(int id);

void CuDeviceProp(void);

void CuAllocateArray(void **devPtr, int size);
void CuSetArrayValue(void *devPtr, int val, size_t size);
void CuFreeArray(void *devPtr);

void CuCopyArrayD2D(void *dDst, void *dSrc, int size);
void CuThreadSync(void);

void CuCopyArrayFromDevice(void* host, const void* device, cudaGraphicsResource **resource, int size);
void CuCopyArrayToDevice(void* device, const void* host, int offset, int size);
void CuRegisterGLBufferObject(unsigned int vbo, cudaGraphicsResource **resource);
void CuUnregisterGLBufferObject(cudaGraphicsResource *resource);
void *CuMapGLBufferObject(cudaGraphicsResource **resource);
void CuUnmapGLBufferObject(cudaGraphicsResource *resource);

void CuScan(unsigned int* output, unsigned int* input, unsigned int numElements);
void CuSort(unsigned int *dHash, uint *dIndex, uint numParticles);

void CuScanf(float* output, float* input, unsigned int numElements);


//-----------------------------------------------------------------------------
// 3D SPH
//-----------------------------------------------------------------------------
void CuSPHInit(int max_particles);
void CuSPHClean(void);

void CuSetParameters(rxSimParams *hostParams);
void CuClearData(void);

// 近傍パーティクル検索用グリッド
void CuCalcHash(uint*  gridParticleHash, uint*  gridParticleIndex, float* pos, int *attr, int numParticles);
void CuReorderDataAndFindCellStart(rxParticleCell cell, float* oldPos, float* oldVel);

void CuCalcHashB(uint* dGridParticleHash, uint* dSortedIndex, float* dPos, 
				 float3 world_origin, float3 cell_width, uint3 grid_size, int nprts);
void CuReorderDataAndFindCellStartB(rxParticleCell cell, float* oldPos);


// SPH計算
void CuSphDensity(float* dDens, float* dPres, rxParticleCell cell);
void CuSphForces(float* dDens, float* dPres, float* dFrc, rxParticleCell cell, float dt);
void CuSphNormal(float* dNrms, float* dDens, rxParticleCell cell);
void CuSphIntegrate(float* pos, float* vel, float* frc, float* dens, int* attr, float dt, uint numParticles);
void CuSphIntegrateWithPolygon(float* pos, float* vel, float* frc, float* dens, int* attr, 
							   float* vrts, int* tris, int tri_num, float dt, rxParticleCell cell);

void CuSphCheckDelete(float* pos, float* vel, int* attr, float minp[3], float maxp[3], float far_point[3], uint nprts);
void CuSphCheckDeleteX(float* pos, float* vel, int* attr, float xmax, float farpoint[3], uint nprts);

// 境界パーティクル
void CuSphBoundaryVolume(float* dVolB, float mass, rxParticleCell cell);
void CuSphBoundaryDensity(float* dDens, float* dPres, float* dPos, float* dVolB, rxParticleCell bcell, uint pnum);
void CuSphBoundaryForces(float* dDens, float* dPres, float* dPos, float* dVolB, float* dFrc, rxParticleCell bcell, uint pnum);


// グリッド
void CuSphGridDensity(float *dGridD, rxParticleCell cell, 
					  int nx, int ny, int nz, float x0, float y0, float z0, float dx, float dy, float dz);

// グリッド　追加：氷用
void CuIceMeshMake(float *dGridD, rxParticleCell cell, 
					  int nx, int ny, int nz, float x0, float y0, float z0, float dx, float dy, float dz, float *bIceFlag);

//-----------------------------------------------------------------------------
// PBDSPH
//-----------------------------------------------------------------------------
void CuPbdSphDensity(float* dDens, rxParticleCell cell);
void CuPbdSphExternalForces(float* dDens, float* dFrc, rxParticleCell cell, float dt);

void CuPbdSphScalingFactor(float* dPos, float* dDens, float* dScl, float eps, rxParticleCell cell);
void CuPbdSphPositionCorrection(float* dPos, float* dScl, float* dDp, rxParticleCell cell);

void CuPbdSphCorrectPosition(float* dPos, float* dDp, uint nprts);
void CuPbdSphUpdatePosition(float* pos, float* new_pos, float* new_vel, float dt, uint nprts);
void CuPbdSphUpdateVelocity(float* pos, float* new_pos, float* new_vel, float dt, uint nprts);

void CuXSphViscosity(float* dPos, float* dVel, float* dNewVel, float* dDens, float c, rxParticleCell cell);

float CuPbdSphCalDensityFluctuation(float* dErrScan, float* dErr, float* dDens, float rest_dens, uint nprts);

void CuPbdSphIntegrate(float* pos, float* vel, float* acc, int* attr, 
					   float* new_pos, float* new_vel, float dt, uint numParticles);
void CuPbdSphIntegrateWithPolygon(float* pos, float* vel, float* acc, int* attr, 
								  float* new_pos, float* new_vel, 
								  float* vrts, int* tris, int tri_num, float dt, rxParticleCell cell);

void CuPbdSphCheckDelete(float* pos, float* vel, int* attr, float minp[3], float maxp[3], float far_point[3], uint nprts);
void CuPbdSphCheckDeleteX(float* pos, float* vel, int* attr, float xmax, float farpoint[3], uint nprts);

void CuPbdSphGridDensity(float *dGridD, rxParticleCell cell, 
						 int nx, int ny, int nz, float x0, float y0, float z0, float dx, float dy, float dz);
void CuPbdSphNormal(float* dNrms, float* dDens, rxParticleCell cell);

// 境界パーティクル
void CuPbdSphBoundaryDensity(float* dDens, float* dPos, float* dVolB, rxParticleCell bcell, uint pnum);
void CuPbdSphScalingFactorWithBoundary(float* dPos, float* dDens, float* dScl, float eps, rxParticleCell cell, 
									   float* dVolB, float* dSclB, rxParticleCell bcell);
void CuPbdSphPositionCorrectionWithBoundary(float* dPos, float* dScl, float* dDp, rxParticleCell cell, 
											float* dVolB, float* dSclB, rxParticleCell bcell);


//-----------------------------------------------------------------------------
// Anisotropic Kernel
//-----------------------------------------------------------------------------
void CuSphCalUpdatedPosition(float* dUpPos, float* dPosW, float lambda, float h, rxParticleCell cell);
void CuSphCalCovarianceMatrix(float* dPosW, float* dCMat, float h, rxParticleCell cell);
void CuSphSVDecomposition(float* dC, float* dPosW, float* dEigen, float* dR, uint numParticles);
void CuSphCalTransformMatrix(float* dEigen, float* dR, float *dG, uint numParticles);

void CuSphGridDensityAniso(float *dGridD, float *dG, float Emax, rxParticleCell cell, 
						   int nx, int ny, int nz, float x0, float y0, float z0, float dx, float dy, float dz);


//-----------------------------------------------------------------------------
// ウェーブレット乱流
//-----------------------------------------------------------------------------
// SPH速度場からEnergy Spectrumを計算
void CuSphES(float* dEt, float scale, float coef_et, float max_et, rxParticleCell cell, float3 cell_size);
void CuAddWTurb3D(float *dVturb, float *dFrc, float *dEt, float *dDens, int d, 
				  int first, int nbands, float* dPos, float pmin[3], float pdim[3], uint numParticles, float dt);

void CuSetWaveletTile3D(float *tile, int n, int d);
void CuSetWaveletTile3DB(float *tile, int nx, int ny, int nz, int d);
void CuSetWaveletTile3DBT(float *tile, int nx, int ny, int nz, int d);


//-----------------------------------------------------------------------------
// SPS乱流
//-----------------------------------------------------------------------------
void CuSubInit(float *dSubPos,float *dPos,float *dSubAxis,float *dSubChild,uint* dSubRand, float radius, uint maxSubLevel, uint nprts);
void CuSubCalEt(float *dSubEt,float *dEt, float et_coef, float scale , float radius, uint maxSubLevel, rxParticleCell cell);
void CuSubCalPos(float *dSubPos,float *dPos, float *dSubChild, float radius, uint maxSubLevel, uint	nprts);
void CuSubUpdateChildAxis(float *dSubChild,float *dSubAxis,float* dSubEt, uint	*dSubRand,float radius, float ratioAxisAngle, uint maxSubLevel, rxParticleCell cell, float dt);

void CuSubInit2(float *dSubPos,float *dPos,float *dSubAxis,float *dSubChild,uint *dSubRand,float radius,uint maxSubLevel,uint nprts, uint uMaxParticles);

void CuSubAdd(float *dSubPos, float *dPos, float *dSubAxis, float *dSubChild, uint *dSubRand, 
			  float radius, uint maxSubLevel, uint nprts, uint uMaxParticles, uint uStart, uint uCount);

void CuSubUpdate(float *dPos, float	*dEt, float *dSubPos, float	*dSubChild, float *dSubAxis, float *dSubEt, uint *dSubRand,
				 float et_coef, float radius, float	ratioAxisAngle,uint maxSubLevel,uint nprts, uint uMaxParticles, float scale,float dt);

//グリッド
void CuSubSetUnsortArray(float *dEt, float *dSubPos, rxSubParticleCell &subCell, float et_coef);
void CuSubCalcHash(rxSubParticleCell subCell);
void CuSubReorderDataAndFindCellStart(rxSubParticleCell subCell);
void CuSubSphGridDensity(float *dGridDl, rxSubParticleCell subCell, 
					  int nx, int ny, int nz, float x0, float y0, float z0, float dx, float dy, float dz);
void CuSubCheckRatio(rxSubParticleCell subCell);

uint CuCheckNumValidHashData(rxSubParticleCell subCell);


//-----------------------------------------------------------------------------
// MC法によるメッシュ化
//-----------------------------------------------------------------------------
#ifdef RX_CUMC_USE_GEOMETRY
void CuMCCreateMesh(float threshold, unsigned int &nvrts, unsigned int &ntris);
#else
void CuMCCreateMesh(GLuint pvbo, GLuint nvbo, float threshold, unsigned int &nvrts, unsigned int &ntris);
#endif


bool CuMCInit(void);

void CuInitMCTable(void);
void CuCleanMCTable(void);

void CuMCCalTriNum(float *dVolume, uint *dVoxBit, uint *dVoxVNum, uint *dVoxVNumScan, 
				   uint *dVoxTNum, uint *dVoxTNumScan, uint *dVoxOcc, uint *dVoxOccScan, uint *dCompactedVox, 
				   uint3 grid_size, uint num_voxels, float3 grid_width, float threshold, 
				   uint &num_active_voxels, uint &nvrts, uint &ntris);
void CuMCCalEdgeVrts(float *dVolume, float *dEdgeVrts, float *dCompactedEdgeVrts, 
					 uint *dEdgeOcc, uint *dEdgeOccScan, uint3 edge_size[3], uint num_edge[4], 
					 uint3 grid_size, uint num_voxels, float3 grid_width, float3 grid_min, float threshold, 
					 uint &nvrts);
void CuMCCalTri(uint *dTris, uint *dVoxBit, uint *dVoxTNumScan, uint *dCompactedVox, 
				uint *dEdgeOccScan, uint3 edge_size[3], uint num_edge[4], 
				uint3 grid_size, uint num_voxels, float3 grid_width, float threshold, 
				uint num_active_voxels, uint nvrts, uint ntris);
void CuMCCalNrm(float *dNrms, uint *dTris, float *dCompactedEdgeVrts, uint nvrts, uint ntris);


//-----------------------------------------------------------------------------
// Screen Space Mesh
//-----------------------------------------------------------------------------
void CuInitFloatArray(float* dArray, float val, int n);
void CuInitRadArray(float* drad, float* dval, int n);

void CuInitDepthMap(float* dDMap, int nx, int ny);
void CuCreateDepthMap(float* dDMap, float* dPrtPos, float* dPrtRad, int pnum, int pdim, float* tr, float* pmv, 
					  int W, int H, float dw, float dh, int nx, int ny);
void CuDepthSmoothing(float* dDMap, int nx, int ny, int n_filter, float *binomials, float zmax);
void CuDetectSilhouetteEdgeVertex(float* dDMap, int nx, int ny, float dw, float dh, float zmax, 
								  float* dPrtPos, float* dPrtRad, int pnum, int pdim, int W, int H, 
								  rxVPacke &dEdge, rxVPackf &dEdgeVrts, rxVPackf &dBackEdgeVrts, rxSSGridG* dMGrid, 
								  int num_nv, int &num_edge, int &num_front_edge, int &num_back_edge);
void CuCalNodeVertex(float* dDMap, int nx, int ny, float dw, float dh, float zmax, 
					 rxVPackf &dNodeVrts, rxSSGridG* dMGrid, int &num_node_vertex);
void CuCalMesh(rxSSGridG* dMGrid, int nx, int ny, float dw, float dh, float zmax, uint* dTriNum, uint* dTriNumScan, 
			   rxVPackf &dBack2Vrts, float* dVrts, int num_vrts, int &num_back2_vrts, uint* &dTriArray, int &num_mesh);
void CuSilhouetteSmoothing(float* dVrts, int num_vrts, int num_node_vrts, uint* dTriArray, int num_mesh, int n_iter);
void CuTransformBack(float* dSSVrts, int num_vrts, float* mvq, float* q, int W, int H, float* dVrts);
void CuCalVertexNormal(float* dVrts, int num_vrts, uint* dTriArray, int num_mesh, float* dNrms);

void CuSetSSMParameters(rxSsmParams *hostParams);

void CuInitVPackf(rxVPackf &pack, int size);
void CuCleanVPackf(rxVPackf &pack);
void CuInitVPacke(rxVPacke &pack, int size);
void CuCleanVPacke(rxVPacke &pack);

void CuInitTable(void);
void CuCleanTable(void);


} // extern "C"


#endif // #ifdef _RX_CU_FUNCS_CUH_