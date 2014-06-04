#include <stdio.h>
#include <math.h>


void LaunchShapeMathcingGPU();
__global__ void d_Integrate();


//GPU処理
void LaunchShapeMatchingGPU(
	/*unsigned int num_particles,
	float (*pos)[2],
	float time,
	float dt*/)
{
	dim3 grid(1, 1);
	dim3 block(1, 1, 1);

	//printf("LaunchGPUKernel");

	d_Integrate <<< grid , block >>> (/*num_particles, pos, time, dt*/);
}


//GPUの位置・速度更新
__global__ void d_Integrate(
	/*int num_particles,	
	float (*pos)[2],
	float time, 
	float dt*/)
{
	//printf("d_Integrate\n");	//めちゃくちゃ出るので注意

	//unsigned int index;
	//float xn, yn, p1, q1, p2, q2, p3, q3, p4, q4;
	//float x, y, t;

 //   // 処理対象の粒子の決定．
 //   index = blockDim.x * blockIdx.x + threadIdx.x;
 //   if (index >= num_particles)
	//	return;
	//xn = pos[index][0];
 //   yn = pos[index][1];

 //   // 1段目．
	//p1 = d_U(xn, yn, time);
	//q1 = d_V(xn, yn, time);

	//// 2段目．
	//x = xn + 0.5f * p1 * dt;
	//y = yn + 0.5f * q1 * dt;
	//t = time + 0.5f * dt;
	//p2 = d_U(x, y, t);
	//q2 = d_V(x, y, t);

	//// 3段目．
	//x = xn + 0.5f * p2 * dt;
	//y = yn + 0.5f * q2 * dt;
	//t = time + 0.5f * dt;
	//p3 = d_U(x, y, t);
	//q3 = d_V(x, y, t);

	//// 4段目．
	//x = xn + p3 * dt;
	//y = yn + q3 * dt;
	//t = time + dt;
	//p4 = d_U(x, y, t);
	//q4 = d_V(x, y, t);
 //  
 //   // 粒子位置の更新．
 //   pos[index][0] = xn + (p1 + p2 + p3 + p4) / 6.0f * dt;
 //   pos[index][1] = yn + (q1 + q2 + q3 + q4) / 6.0f * dt;
}