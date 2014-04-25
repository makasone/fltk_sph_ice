//２次元に見える１次元vector
//要素をうまい具合にソートするなど，難しい処理には対応していない
//(i,j)でイテレータに触れると便利！
//完全に破棄したいなら，かっこ悪いけどポインタで作ってdeleteするしかない．

#include <vector>
using namespace std;

template<class Type> class mk_Vector2D
{
public:
	mk_Vector2D():mVector(0)
	{
	}

	mk_Vector2D(int sizeX, int sizeY):mVector(0), mSizeX(sizeX), mSizeY(sizeY)
	{
		SetSize(sizeX, sizeY);
	}

	~mk_Vector2D()
	{
		vector<Type>().swap(mVector);
	}

	//配列確保
	void SetSize(int sizeX, int sizeY)
	{
		mSizeX = sizeX;
		mSizeY = sizeY;
		mVector.resize(sizeX*sizeY);
	}

	//vectorのメンバ関数を使いたい場合は，これを経由して使う．
	vector<Type>& Get()
	{
		return mVector;
	}

	//演算子
	Type& operator()(int indx0, int indx1)
	{
		return mVector[indx1*mSizeX+indx0];
	}

	const Type& operator()(int indx0, int indx1) const
	{
		return mVector[indx1*mSizeX+indx0];
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
	int mSizeX;
	int mSizeY;
};