//�Q�����Ɍ�����P����vector
//�v�f�����܂���Ƀ\�[�g����ȂǁC��������ɂ͑Ή����Ă��Ȃ�
//(i,j)�ŃC�e���[�^�ɐG���ƕ֗��I
//���S�ɔj���������Ȃ�C�������������ǃ|�C���^�ō����delete���邵���Ȃ��D

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

	//�z��m��
	void SetSize(int sizeX, int sizeY)
	{
		mSizeX = sizeX;
		mSizeY = sizeY;
		mVector.resize(sizeX*sizeY);
	}

	//vector�̃����o�֐����g�������ꍇ�́C������o�R���Ďg���D
	vector<Type>& Get()
	{
		return mVector;
	}

	//���Z�q
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