#include <stdio.h>
#include <math.h>


void LaunchShapeMathcingGPU();
__global__ void d_Integrate();


//GPU����
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


//GPU�̈ʒu�E���x�X�V
__global__ void d_Integrate(
	/*int num_particles,	
	float (*pos)[2],
	float time, 
	float dt*/)
{
	//printf("d_Integrate\n");	//�߂��Ⴍ����o��̂Œ���

	//unsigned int index;
	//float xn, yn, p1, q1, p2, q2, p3, q3, p4, q4;
	//float x, y, t;

 //   // �����Ώۂ̗��q�̌���D
 //   index = blockDim.x * blockIdx.x + threadIdx.x;
 //   if (index >= num_particles)
	//	return;
	//xn = pos[index][0];
 //   yn = pos[index][1];

 //   // 1�i�ځD
	//p1 = d_U(xn, yn, time);
	//q1 = d_V(xn, yn, time);

	//// 2�i�ځD
	//x = xn + 0.5f * p1 * dt;
	//y = yn + 0.5f * q1 * dt;
	//t = time + 0.5f * dt;
	//p2 = d_U(x, y, t);
	//q2 = d_V(x, y, t);

	//// 3�i�ځD
	//x = xn + 0.5f * p2 * dt;
	//y = yn + 0.5f * q2 * dt;
	//t = time + 0.5f * dt;
	//p3 = d_U(x, y, t);
	//q3 = d_V(x, y, t);

	//// 4�i�ځD
	//x = xn + p3 * dt;
	//y = yn + q3 * dt;
	//t = time + dt;
	//p4 = d_U(x, y, t);
	//q4 = d_V(x, y, t);
 //  
 //   // ���q�ʒu�̍X�V�D
 //   pos[index][0] = xn + (p1 + p2 + p3 + p4) / 6.0f * dt;
 //   pos[index][1] = yn + (q1 + q2 + q3 + q4) / 6.0f * dt;
}