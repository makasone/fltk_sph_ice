/*! 
  @file rx_cu_funcs.cu
	
  @brief CUDA�֐� - �������֌W�Ȃ�

  @author Makoto Fujisawa
  @date 2009-08, 2011-06
*/
// FILE --rx_cu_funcs.cu--


//-----------------------------------------------------------------------------
// �C���N���[�h�t�@�C��
//-----------------------------------------------------------------------------
#include <cstdio>
#include <GL/glew.h>
#include <GL/glut.h>

#include <thrust/device_vector.h>
#include <thrust/scan.h>
#include <thrust/sort.h>

#include <GL/freeglut.h>

#include "rx_cu_common.cu"



//-----------------------------------------------------------------------------
// CUDA�֐�
//-----------------------------------------------------------------------------
extern "C"
{
/*!
 * CUDA�f�o�C�X�̐ݒ�
 *  - �R�}���h���C�������Ɋ�Â�CUDA�f�o�C�X��ݒ�((��)-device 0)
 * @param[in] argc �R�}���h���C�������̐�
 * @param[in] argv �R�}���h���C���������X�g(argv[0]�͎��s�t�@�C����)
 */
void CuInit(int argc, char **argv)
{   
	if(checkCmdLineFlag(argc, (const char**)argv, "device")){
		int id = getCmdLineArgumentInt(argc, (const char**)argv, "device=");
		if(id < 0){
			id = gpuGetMaxGflopsDeviceId();
			cudaSetDevice(id);
		}
		else{
			cudaSetDevice(id);
		}
	}
	else{
		cudaSetDevice( gpuGetMaxGflopsDeviceId() );
	}
}

/*!
 * CUDA�f�o�C�X�̐ݒ� - id�𒼐ڎw��
 * @param[in] id �f�o�C�XID
 */
void CuSetDevice(int id)
{ 
	int device_count = 0;
	cudaGetDeviceCount(&device_count);
	if(id < 0 || id >= device_count){
		id = gpuGetMaxGflopsDeviceId();
		cudaSetDevice(0);
	}
	else{
		cudaSetDevice(id);
	}
}

/*!
 * �f�o�C�X�������̊m��
 * @param[out] dPtr �f�o�C�X�������ւ̃|�C���^
 * @param[in] size �m�ۃT�C�Y(��������̃T�C�Y)
 */
void CuAllocateArray(void **dPtr, size_t size)
{
	RX_CUCHECK(cudaMalloc(dPtr, size));
}

/*!
 * �f�o�C�X�������̉��
 * @param[in] devPtr �f�o�C�X�������ւ̃|�C���^
 */
void CuFreeArray(void *dPtr)
{
	RX_CUCHECK(cudaFree(dPtr));
}

/*!
 * �f�o�C�X�������̈�̏�����
 * @param[in] dPtr �f�o�C�X�������ւ̃|�C���^
 * @param[in] val �����l
 * @param[in] size ����������̈�̃T�C�Y(��������̃T�C�Y)
 */
void CuSetArrayValue(void *dPtr, int val, size_t size)
{
	RX_CUCHECK(cudaMemset(dPtr, val, size));
}

/*!
 * �f�o�C�X�������ԃR�s�[
 * @param[in] dDst �R�s�[��
 * @param[in] dSrc �R�s�[��
 * @param[in] size �R�s�[�T�C�Y(��������̃T�C�Y)
 */
void CuCopyArrayD2D(void *dDst, void *dSrc, int size)
{
	RX_CUCHECK(cudaMemcpy(dDst, dSrc, size, cudaMemcpyDeviceToDevice));
}


/*!
 * VBO���}�b�s���O
 * @param[in] vbo VBO,PBO��
 */
void *CuMapGLBufferObject(cudaGraphicsResource **resource)
{
	void *ptr;
	RX_CUCHECK(cudaGraphicsMapResources(1, resource, 0));
	size_t num_bytes;
	RX_CUCHECK(cudaGraphicsResourceGetMappedPointer((void**)&ptr, &num_bytes, *resource));
	return ptr;
}

/*!
 * VBO���A���}�b�v
 * @param[in] vbo VBO,PBO��
 */
void CuUnmapGLBufferObject(cudaGraphicsResource *resource)
{
	RX_CUCHECK(cudaGraphicsUnmapResources(1, &resource, 0));
}

/*!
 * PBO,VBO�o�b�t�@��CUDA�ɓo�^
 * @param[in] vbo VBO,PBO��
 */
void CuRegisterGLBufferObject(uint vbo, cudaGraphicsResource **resource)
{
	RX_CUCHECK(cudaGraphicsGLRegisterBuffer(resource, vbo, cudaGraphicsMapFlagsNone));
}

/*!
 * PBO,VBO�o�b�t�@��CUDA����폜
 * @param[in] vbo VBO,PBO��
 */
void CuUnregisterGLBufferObject(cudaGraphicsResource *resource)
{
	RX_CUCHECK(cudaGraphicsUnregisterResource(resource));
}

/*!
 * �f�o�C�X����z�X�g�������ւ̃R�s�[
 * @param[in] hDst �R�s�[��z�X�g������(�Œ�size���m�ۂ���Ă��邱��)
 * @param[in] dSrc �R�s�[���f�o�C�X������
 * @param[in] vbo dSrc��VBO�̏ꍇ�CVBO��ID�D�����łȂ��ꍇ��0���w��
 * @param[in] size �R�s�[�T�C�Y(��������̃T�C�Y)
 */
void CuCopyArrayFromDevice(void* hDst, const void* dSrc, cudaGraphicsResource **resource, int size)
{   
	if(resource) dSrc = CuMapGLBufferObject(resource);

	RX_CUCHECK(cudaMemcpy(hDst, dSrc, size, cudaMemcpyDeviceToHost));
	
	if(resource) CuUnmapGLBufferObject(*resource);
}

/*!
 * �z�X�g����f�o�C�X�������ւ̃R�s�[
 * @param[in] dDst �R�s�[��f�o�C�X������(�Œ�size���m�ۂ���Ă��邱��)
 * @param[in] hSrc �R�s�[���z�X�g������
 * @param[in] offset �R�s�[��I�t�Z�b�g
 * @param[in] size �R�s�[�T�C�Y(��������̃T�C�Y)
 */
void CuCopyArrayToDevice(void* dDst, const void* hSrc, int offset, int size)
{
	RX_CUCHECK(cudaMemcpy((char*)dDst+offset, hSrc, size, cudaMemcpyHostToDevice));
}

/*!
 * �X���b�h����
 */
void CuThreadSync(void)
{
	RX_CUCHECK(cudaThreadSynchronize());
}

/*!
 * �f�o�C�X�v���p�e�B�̕\��
 */
void CuDeviceProp(void)
{
	int n;	//�f�o�C�X��
	RX_CUCHECK(cudaGetDeviceCount(&n));

	for(int i = 0; i < n; ++i){
		cudaDeviceProp dev;

		// �f�o�C�X�v���p�e�B�擾
		RX_CUCHECK(cudaGetDeviceProperties(&dev, i));

		printf("device %d\n", i);
		printf(" device name : %s\n", dev.name);
		printf(" total global memory : %d (MB)\n", dev.totalGlobalMem/1024/1024);
		printf(" shared memory / block : %d (KB)\n", dev.sharedMemPerBlock/1024);
		printf(" register / block : %d\n", dev.regsPerBlock);
		printf(" warp size : %d\n", dev.warpSize);
		printf(" max pitch : %d (B)\n", dev.memPitch);
		printf(" max threads / block : %d\n", dev.maxThreadsPerBlock);
		printf(" max size of each dim. of block : (%d, %d, %d)\n", dev.maxThreadsDim[0], dev.maxThreadsDim[1], dev.maxThreadsDim[2]);
		printf(" max size of each dim. of grid  : (%d, %d, %d)\n", dev.maxGridSize[0], dev.maxGridSize[1], dev.maxGridSize[2]);
		printf(" clock rate : %d (MHz)\n", dev.clockRate/1000);
		printf(" total constant memory : %d (KB)\n", dev.totalConstMem/1024);
		printf(" compute capability : %d.%d\n", dev.major, dev.minor);
		printf(" alignment requirement for texture : %d\n", dev.textureAlignment);
		printf(" device overlap : %s\n", (dev.deviceOverlap ? "ok" : "not"));
		printf(" num. of multiprocessors : %d\n", dev.multiProcessorCount);
		printf(" kernel execution timeout : %s\n", (dev.kernelExecTimeoutEnabled ? "on" : "off"));
		printf(" integrated : %s\n", (dev.integrated ? "on" : "off"));
		printf(" host memory mapping : %s\n", (dev.canMapHostMemory ? "on" : "off"));

		printf(" compute mode : ");
		if(dev.computeMode == cudaComputeModeDefault) printf("default mode (multiple threads can use) \n");
		else if(dev.computeMode == cudaComputeModeExclusive) printf("exclusive mode (only one thread will be able to use)\n");
		else if(dev.computeMode == cudaComputeModeProhibited) printf("prohibited mode (no threads can use)\n");
		
	}

	printf("Device with Maximum GFLOPS : %d\n", gpuGetMaxGflopsDeviceId());
}

/*!
 * thrust::exclusive_scan�̌Ăяo��
 * @param[out] dScanData scan��̃f�[�^
 * @param[in] dData ���f�[�^
 * @param[in] num �f�[�^��
 */
void CuScan(unsigned int* dScanData, unsigned int* dData, unsigned int num)
{
	thrust::exclusive_scan(thrust::device_ptr<unsigned int>(dData), 
						   thrust::device_ptr<unsigned int>(dData+num),
						   thrust::device_ptr<unsigned int>(dScanData));
}

/*!
 * thrust::sort_by_key�ɂ��n�b�V���l�Ɋ�Â��\�[�g
 * @param[in] dHash �n�b�V���l
 * @param[in] dIndex �C���f�b�N�X(�p�[�e�B�N���C�|���S���Ȃ�)
 * @param[in] num �f�[�^��
 */
void CuSort(unsigned int *dHash, uint *dIndex, uint num)
{
	thrust::sort_by_key(thrust::device_ptr<unsigned int>(dHash),
						thrust::device_ptr<unsigned int>(dHash+num),
						thrust::device_ptr<unsigned int>(dIndex));
}




}   // extern "C"
