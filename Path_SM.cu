#ifndef _GPU_PATH_SM_
#define _GPU_PATH_SM_

#include "Path_SM.cuh"

//TODO: timeStepを引数として渡せるように
void LaunchUpdatePrefixSumGPU
(
	int prtNum,
	int PosSizeX,
	int PosSizeY,
	int ApqSizeX,
	int ApqSizeY,
	int* md_2DiPTHtoPRT,
	int* md_2DiPRTtoPTH,
	float* md_2Df3PrfxPos,
	float* md_2Df9PrfxApq,
	int* md_3DiPTHandPrfxSet,
	float* md_f3OrgPos,
	float* md_f3OrgCm,
	unsigned* dgroupPos,
	unsigned* dgroupApq,
	const float* md_fPos,
	const float* md_fVel
)
{	
	//cout << __FUNCTION__ << endl;

	int side = pow(prtNum, 1.0/3.0) + 0.5;	//立方体の１辺の頂点数

	dim3 grid(side, side);
	dim3 block(side, 1, 1);

	//prefixSumの更新
	//posとvelからデータ作成
	UpdatePrefixSumGPU<<<grid ,block>>>
	(
		prtNum,
		PosSizeX, PosSizeY,
		ApqSizeX, ApqSizeY,
		md_2DiPTHtoPRT, md_2DiPRTtoPTH,
		md_2Df3PrfxPos,
		md_2Df9PrfxApq,
		md_3DiPTHandPrfxSet,
		md_f3OrgPos,	md_f3OrgCm,
		md_fPos,		md_fVel
	);

	//ThrustでprefixSum計算
	thrust::inclusive_scan_by_key(
		thrust::device_ptr<unsigned>(dgroupPos),
		thrust::device_ptr<unsigned>(dgroupPos + prtNum*3),
		thrust::device_ptr<float>(md_2Df3PrfxPos),
		thrust::device_ptr<float>(md_2Df3PrfxPos)
	);

	//Apq
	thrust::inclusive_scan_by_key(
		thrust::device_ptr<unsigned>(dgroupApq),
		thrust::device_ptr<unsigned>(dgroupApq + prtNum*9),
		thrust::device_ptr<float>(md_2Df9PrfxApq),
		thrust::device_ptr<float>(md_2Df9PrfxApq)
	);

	cudaThreadSynchronize();
}

//Trustのテスト
void ThrustTest()
{
	cout << __FUNCTION__ << endl;

	float* dPosOutData;
	float* dApqOutData;
	cudaMalloc((void**)&dPosOutData,	sizeof(float) * 3*3);
	cudaMalloc((void**)&dApqOutData,	sizeof(float) * 3*3);

	float* posData = new float[3 * 3];
	for(int i = 0; i < 3 * 3; i++)
	{
		posData[i] = (float)i;
	}

	for(int i = 0; i < 3 * 3; i++)
	{
		cout << "in:" << posData[i] << endl;
	}

	cudaMemcpy(dPosOutData, posData, sizeof(float) * 3*3, cudaMemcpyHostToDevice);	

	for(int indx = 0; indx < 3; indx++)
	{
		int startIndx = 3 * indx;
		int endIndx = 3 * (indx+1);
		thrust::inclusive_scan(
			thrust::device_ptr<float>(dPosOutData + startIndx),
			thrust::device_ptr<float>(dPosOutData + endIndx), 
			thrust::device_ptr<float>(dApqOutData + startIndx)); // in-place scan
	}

	cudaMemcpy(posData, dApqOutData, sizeof(float) * 3*3, cudaMemcpyDeviceToHost);
	
	for(int i = 0; i < 3*3; i++)
	{
		cout << "out:" << posData[i] << endl;
	}

	delete[] posData;
	cudaFree(dPosOutData);
	cudaFree(dApqOutData);
}

void ThrustTest2()
{
	int data[10] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
	int keys[10] = {0, 0, 0, 1, 1, 2, 3, 3, 3, 3};
	int out[10];

	int* group = new int[10];
	for(int i = 0; i < 10; i++)
	{
		if( i < 5)
		{
			group[i] = 0;
		}
		else if( i >= 5 && i < 8)
		{
			group[i] = 1;
		}
		else
		{
			group[i] = 2;
		}

		cout << "group:" << group[i] << endl;
	}

	thrust::inclusive_scan_by_key(group, group + 10, data, out); // in-place scan

	for(int i = 0; i < 10; i++)
	{
		cout << "out:" << out[i] << endl;
	}

	delete[] group;
}

__global__
void UpdatePrefixSumGPU(
	int prtNum,
	int PosSizeX,
	int PosSizeY,
	int ApqSizeX,
	int ApqSizeY,
	int* md_2DiPTHtoPRT,
	int* md_2DiPRTtoPTH,
	float* prfxPos,
	float* prfxApq,
	int* md_3DiPTHandPrfxSet,
	float* orgPos,
	float* md_f3OrgCm,
	const float* pos,
	const float* vel
	)
{
	//UpdatePrefixSumPosGPU(md_2DiPTHtoPRT, prfxPos, pos, vel, prtNum);
	//UpdatePrefixSumApqGPU(md_2DiPTHtoPRT, prfxApq, orgPos, pos, vel, prtNum);

	//TODO: パスは一本であるという前提の式になっている
	int thid = blockIdx.x * (blockDim.x * blockDim.y)
				+ blockIdx.y * (blockDim.x * blockDim.y) * gridDim.x
				+ threadIdx.x;

	int pvIndx = thid*4;

	float3 p;
	float3 q;

	p.x = pos[pvIndx+0] + vel[pvIndx+0] * 0.01f;
	p.y = pos[pvIndx+1] + vel[pvIndx+1] * 0.01f;
	p.z = pos[pvIndx+2] + vel[pvIndx+2] * 0.01f;

	q.x = orgPos[thid*3+0];
	q.y = orgPos[thid*3+1];
	q.z = orgPos[thid*3+2];

//Pos
	//posとvelのデータをprfixPosに入れるだけ
	//prfixPosのデータの並びは，x[0][1][2]...y[0][1][2]...z[0][1][2]...とする
	//x[0]y[0]z[0] x[1]y[1]z[1]...ではないのに注意
	prfxPos[prtNum*0+thid] = p.x;	//x
	prfxPos[prtNum*1+thid] = p.y;	//y
	prfxPos[prtNum*2+thid] = p.z;	//z

//Apq
	prfxApq[prtNum*0+thid] = p.x * q.x;
	prfxApq[prtNum*1+thid] = p.x * q.y;
	prfxApq[prtNum*2+thid] = p.x * q.z;
	
	prfxApq[prtNum*3+thid] = p.y * q.x;
	prfxApq[prtNum*4+thid] = p.y * q.y;
	prfxApq[prtNum*5+thid] = p.y * q.z;
	
	prfxApq[prtNum*6+thid] = p.z * q.x;
	prfxApq[prtNum*7+thid] = p.z * q.y;
	prfxApq[prtNum*8+thid] = p.z * q.z;
}

__device__
void UpdatePrefixSumPosGPU
(
	int* PTHtoPRT,
	float* prfixPos,
	const float* pos, 
	const float* vel,
	int num
)
{
	//TODO: パスは一本であるという前提の式になっている
	int thid = blockIdx.x * (blockDim.x * blockDim.y)
				+ blockIdx.y * (blockDim.x * blockDim.y) * gridDim.x
				+ threadIdx.x;

	//posとvelのデータをprfixPosに入れるだけ
	//prfixPosのデータの並びは，x[0][1][2]...y[0][1][2]...z[0][1][2]...とする
	//x[0]y[0]z[0] x[1]y[1]z[1]...ではないのに注意
	int pvIndx = thid*4;
	prfixPos[num*0+thid] = pos[pvIndx+0] + vel[pvIndx+0] * 0.01f;	//x
	prfixPos[num*1+thid] = pos[pvIndx+1] + vel[pvIndx+1] * 0.01f;	//y
	prfixPos[num*2+thid] = pos[pvIndx+2] + vel[pvIndx+2] * 0.01f;	//z
}

__device__
void UpdatePrefixSumApqGPU
(
	int* PTHtoPRT,
	float* prfixApq,
	float* orgPos,
	const float* pos, 
	const float* vel,
	int num
)
{
	int idxA = blockIdx.x * (blockDim.x * blockDim.y);
	int idxB = blockIdx.y * (blockDim.x * blockDim.y) * gridDim.x;
	int idxC = threadIdx.x;
	int thid = idxA + idxB + idxC;

	float3 p;
	float3 q;

	p.x = pos[thid*4+0] + vel[thid*4+0] * 0.01f;
	p.y = pos[thid*4+1] + vel[thid*4+1] * 0.01f;
	p.z = pos[thid*4+2] + vel[thid*4+2] * 0.01f;

	q.x = orgPos[thid*3+0];
	q.y = orgPos[thid*3+1];
	q.z = orgPos[thid*3+2];

	prfixApq[num*0+thid] = p.x * q.x;
	prfixApq[num*1+thid] = p.x * q.y;
	prfixApq[num*2+thid] = p.x * q.z;
	
	prfixApq[num*3+thid] = p.y * q.x;
	prfixApq[num*4+thid] = p.y * q.y;
	prfixApq[num*5+thid] = p.y * q.z;
	
	prfixApq[num*6+thid] = p.z * q.x;
	prfixApq[num*7+thid] = p.z * q.y;
	prfixApq[num*8+thid] = p.z * q.z;
}

//--------------------------------------------------------------------------------------------
__device__
void scan_A
(
	int* g_idata,
	int* g_odata,
	int n
)
{
	__shared__ float temp[TESTDATANUM*2];

	int thid = threadIdx.x;
	int pout = 0, pin = 1;

	// load input into shared memory.
	// Exclusive scan: shift right by one and set first element to 0
	temp[thid] = g_idata[thid];
	temp[thid+TESTDATANUM] = g_idata[thid];

	__syncthreads();

	for( int offset = 1; offset < n; offset <<= 1 )
	{
		pout = 1 - pout; // swap double buffer indices
		pin = 1 - pout;
		if (thid >= offset)
		temp[pout*n+thid] = temp[pin*n+thid] + temp[pin*n+thid - offset];
		else
		temp[pout*n+thid] = temp[pin*n+thid];
		__syncthreads();
	}

	g_odata[thid] = temp[pout*n+thid]; // write output
}

__device__
	void scan_B
(
	int* g_idata,
	int* g_odata,
	int n
)
{
	//printf("data:%d\n", inData[threadIdx.x]);

	__shared__ float temp[TESTDATANUM*2];
	int thid = threadIdx.x;
	int offset = 1;

	temp[2*thid] = g_idata[2*thid];
	temp[2*thid+1] = g_idata[2*thid+1];

	//右シフトして値を小さくしていく
	for(int d = n>>1; d > 0; d >>= 1)
	{
		__syncthreads();
		
		if(thid < d)
		{
			int ai = offset*(2*thid+1)-1;
			int bi = offset*(2*thid+2)-1;

			temp[bi] += temp[ai];
		}

		offset *= 2;
	}
	
	int last = -1;	//最後の要素を保存しておく
	if(thid == 0){ last = temp[n-1]; temp[n-1] = 0; }

	for (int d = 1; d < n; d *= 2) // traverse down tree & build scan
	{
		offset >>= 1;
		__syncthreads();
	
		if (thid < d)
		{
			int ai = offset*(2*thid+1)-1;
			int bi = offset*(2*thid+2)-1;
	
			float t = temp[ai];
			temp[ai] = temp[bi];
			temp[bi] += t;
		}
	}
	
	__syncthreads();
	
	if(thid >= n){	return;	}

	if(thid > 0)
	{
		g_odata[2*thid   -1] = temp[2*thid]; // write results to device memory
		g_odata[2*thid+1 -1] = temp[2*thid+1];
	}
	else if(thid == 0)
	{
		g_odata[0] = temp[2*thid+1];	//最初の要素
		g_odata[n-1] = last;			//最後の要素
	}
}

__device__
void scan_C
(
	int* g_idata,
	int* g_odata,
	int n
)
{
	//int thid = threadIdx.x;

	int idxA = blockIdx.x * (blockDim.x * blockDim.y);
	int idxB = blockIdx.y * (blockDim.x * blockDim.y) * gridDim.x;
	int idxC = threadIdx.x;
	int thid = idxA + idxB + idxC;

	int pout = 0, pin = 1;

	g_odata[thid] = g_idata[thid];
	g_odata[thid+n] = g_idata[thid];

	__syncthreads();
	//__threadfence();

	for( int offset = 1; offset < n; offset <<= 1 )
	{
		pout = 1 - pout; // swap double buffer indices
		pin = 1 - pout;
		if (thid >= offset)
			g_odata[pout*n+thid] = g_odata[pin*n+thid] + g_odata[pin*n+thid - offset];
		else
			g_odata[pout*n+thid] = g_odata[pin*n+thid];

		__syncthreads();
		//__threadfence();
		//__threadfence_block();
	}

	g_idata[thid] = g_odata[pout*n+thid]; // write output
}

__device__
	void scan_D
(
	int* g_idata,
	int* g_odata,
	int n
)
{
	int thid = threadIdx.x;
	int offset = 1;

	g_odata[2*thid] = g_idata[2*thid];
	g_odata[2*thid+1] = g_idata[2*thid+1];

	//右シフトして値を小さくしていく
	for(int d = n>>1; d > 0; d >>= 1)
	{
		__syncthreads();
		
		if(thid < d)
		{
			int ai = offset*(2*thid+1)-1;
			int bi = offset*(2*thid+2)-1;

			g_odata[bi] += g_odata[ai];
		}

		offset *= 2;
	}
	
	int last = -1;	//最後の要素を保存しておく
	if(thid == 0){ last = g_odata[n-1]; g_odata[n-1] = 0; }

	for (int d = 1; d < n; d *= 2) // traverse down tree & build scan
	{
		offset >>= 1;
		__syncthreads();
	
		if (thid < d)
		{
			int ai = offset*(2*thid+1)-1;
			int bi = offset*(2*thid+2)-1;
	
			float t = g_odata[ai];
			g_odata[ai] = g_odata[bi];
			g_odata[bi] += t;
		}
	}
	
	__syncthreads();
	
	if(thid >= n){	return;	}

	if(thid > 0)
	{
		g_idata[2*thid   -1] = g_odata[2*thid]; // write results to device memory
		g_idata[2*thid+1 -1] = g_odata[2*thid+1];
	}
	else if(thid == 0)
	{
		g_idata[0] = g_odata[1];		//最初の要素
		g_idata[n-1] = last;			//最後の要素
	}
}

#endif