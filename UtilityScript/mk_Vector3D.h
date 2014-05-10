//�R�����Ɍ�����P����vector
//�v�f�����܂���Ƀ\�[�g����ȂǁC��������ɂ͑Ή����Ă��Ȃ�
//(i,j)�ŃC�e���[�^�ɐG���ƕ֗��I
//���S�ɔj���������Ȃ�C�������������ǃ|�C���^�ō����delete���邵���Ȃ��D

#ifndef _MK_VECTOR3D_
#define _MK_VECTOR3D_

#include <vector>
using namespace std;

template<class Type> class mk_Vector3D
{
public:
	mk_Vector3D():mVector(0)
	{
	}

	mk_Vector3D(int sizeX, int sizeY, int sizeZ) : mVector(0), mSizeX(sizeX), mSizeY(sizeY), mSizeZ(sizeZ)
	{
		SetSize(sizeX, sizeY, sizeZ);
	}

	~mk_Vector3D()
	{
		vector<Type>().swap(mVector);		//STL�̃X���b�v�Z�@�ŉ������
	}

	//�z��m��
	void SetSize(int sizeX, int sizeY, int sizeZ)
	{
		mSizeX = sizeX;
		mSizeY = sizeY;
		mSizeZ = sizeZ;
		mVector.resize(sizeX*sizeY*sizeZ);
	}

	//vector�̃����o�֐����g�������ꍇ�́C������o�R���Ďg���D
	vector<Type>& Get()
	{
		return mVector;
	}

	int GetSizeX()
	{
		return mSizeX;
	}

	int GetSizeY()
	{
		return mSizeY;
	}

	int GetSizeZ()
	{
		return mSizeZ;
	}

	//���Z�q
	//Y�ŘA�������Ă���.
	Type& operator()(int X, int Y, int Z)
	{
		return mVector[mSizeY*mSizeZ*X + mSizeY*Z + Y];
	}

	const Type& operator()(int X, int Y, int Z) const
	{
		return mVector[mSizeY*mSizeZ*X + mSizeY*Z + Y];
	}

	Type& operator[](int i)
	{
		return mVector[i];
	}

	const Type& operator()(int i) const
	{
		return mVector[i];
	}

private:
	vector<Type> mVector;
public:
	int mSizeX;
	int mSizeY;
	int mSizeZ;
};

#endif