//３次元に見える１次元vector
//要素をうまい具合にソートするなど，難しい処理には対応していない
//(i,j)でイテレータに触れると便利！
//完全に破棄したいなら，かっこ悪いけどポインタで作ってdeleteするしかない．

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
		vector<Type>().swap(mVector);		//STLのスワップ技法で解放する
	}

	//配列確保
	void SetSize(int sizeX, int sizeY, int sizeZ)
	{
		mSizeX = sizeX;
		mSizeY = sizeY;
		mSizeZ = sizeZ;
		mVector.resize(sizeX*sizeY*sizeZ);
	}

	//vectorのメンバ関数を使いたい場合は，これを経由して使う．
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

	//演算子
	//Yで連続させている.
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