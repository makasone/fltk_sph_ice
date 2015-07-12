//mk_Quaternion�FEigen�������t�@�C���Ŗ]�񂾂悤�ɃC���N���[�h�ł��Ȃ��������߂ɍ�����N���X�@�o����Ύg�������Ȃ�
//mk_ExpMap�F�ΐ��N�H�[�^�j�I��

#ifndef _MK_QUATERNION_
#define _MK_QUATERNION_

#include <vector>
#include <math.h>

using namespace std;

class mk_Quaternion;
class mk_ExpMap;

mk_ExpMap QuaternionToExpMap(const mk_Quaternion& q);
mk_Quaternion ExpMapToQuaterinon(const mk_ExpMap& eq);

//�ꎞ���̂��N�H�[�^�j�I��
class mk_Quaternion{
public:
	float	w, x, y, z;

public:
	mk_Quaternion(){
		mk_Quaternion(1.0f, 0.0f, 0.0f, 0.0f);
	}

	mk_Quaternion(float xx, float yy, float zz, float ww){
		x = xx; y = yy; z = zz; w = ww;
	}
};

//ExpMap�@�N�H�[�^�j�I���̑ΐ�
//�ʂɃN���X�𕪂���K�v�͂Ȃ������̂ł́c�H
class mk_ExpMap{
public:
	float x, y, z;

public:
	//���܂��ĂׂȂ��݂����H�@�Ȃ�������������Ȃ�
	mk_ExpMap(){
		mk_ExpMap(0.0f, 0.0f, 0.0f);
	}

	mk_ExpMap(float vx, float vy, float vz){
		x = vx; y = vy; z = vz;
	}

	mk_ExpMap(const mk_Quaternion& q){
		*this = QuaternionToExpMap(q);
	}

//���Z�q
	inline mk_ExpMap operator+(const mk_ExpMap& a)
	{
		return mk_ExpMap(x + a.x, y + a.y, z + a.z);
	}

	inline mk_ExpMap operator*(const float a)
	{
		return mk_ExpMap(x * a, y * a, z * a);
	}

	inline mk_ExpMap& operator=(const mk_ExpMap& a)
	{
		x = a.x;
		y = a.y;
		z = a.z;
		return *this;
	}

	inline mk_ExpMap& operator+=(const mk_ExpMap& a)
	{
		x += a.x;
		y += a.y;
		z += a.z;
		return *this;
	}

	//ExpMap�̕��
	inline mk_ExpMap ExpLinerInterpolation(mk_ExpMap exp, float weight){
		return (*this) * weight + exp * (1.0f-weight);
	}

	//ExpMap�̕��
	inline mk_ExpMap ExpLinerInterpolation(const vector<mk_ExpMap>& expes, const vector<float>& weights, unsigned prtNum){
		
		mk_ExpMap exp = mk_ExpMap(0.0f, 0.0f, 0.0f);

		for(unsigned i = 0; i < prtNum; i++){
			mk_ExpMap q = expes.at(i);
			exp += q * weights.at(i);
		}

		return exp;
	}
};

	//�N�H�[�^�j�I����ExpMap
	inline static mk_ExpMap QuaternionToExpMap(const mk_Quaternion& q){

		//�P�ʃN�H�[�^�j�I����ExpMap��Zero��Ή��t��
		if( fabs(q.w) >= 1.0f - 0.00001f ){
			return mk_ExpMap(0.0f, 0.0f, 0.0f);
		}
		
		const float theta = acos(q.w);
		const float thetaInV = theta / sin(theta);

		return mk_ExpMap(q.x * thetaInV, q.y * thetaInV, q.z * thetaInV);
	}

	//ExpMap���N�H�[�^�j�I��
	inline static mk_Quaternion ExpMapToQuaterinon(const mk_ExpMap& eq){

		const float theta = sqrtf(eq.x * eq.x + eq.y * eq.y + eq.z * eq.z);

		//�[�����Z�������
		if(theta < 0.001f){
			return mk_Quaternion(0.0f, 0.0f, 0.0f, 1.0f);
		}

		const float sintheta = sinf(theta) / theta;

		return mk_Quaternion(eq.x * sintheta, eq.y * sintheta, eq.z * sintheta, cosf(theta));
	}

#endif