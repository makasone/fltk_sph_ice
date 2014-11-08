//-----------------------------------------------------------------------------
// �C���N���[�h�t�@�C��
//-----------------------------------------------------------------------------
#include "HeatTransfar.h"
#include <stdio.h>
#include <iostream>

using namespace std;

#define RX_PI          (double)(3.1415926535897932384626433832795)   

HeatTransfar::HeatTransfar(int num)
{	
	mNumVertices = num;

	mAirTemp = 0.0f;

	mTemps = new float[mNumVertices];				//���x
	mTempsDelta = new float[mNumVertices];			//�ω����x
	mHeats = new float[mNumVertices];				//�M��

	mPhase = new int[mNumVertices];					//���݂̏��
	mPhaseChange = new int[mNumVertices];			//���ω����s�����̃t���O

	mSurfaceParticleNums = new int[mNumVertices];

	mNeighborhoodsId.clear();
	mNeighborhoodsDis.clear();

	//���q�̉��x�C�M�ʁC��Ԃ̏�����
	for( int i = 0; i < mNumVertices; i++)
	{
		mTemps[i] = 0.0f;
		mHeats[i] = 0.0f;
		if( mTemps[i] < 250 )
		{
			mPhase[i]	 = -2;	//�X
			mPhaseChange[i] = 0;
		}
		else
		{
			mPhase[i]	 = 2;	//��
			mPhaseChange[i] = 0;
		}
	}
}

HeatTransfar::~HeatTransfar(void)
{
}


void HeatTransfar::initState()
{

}

void HeatTransfar::MeltParticle(int pIndx)
{
	//���M�ω��I��
	setTemps(pIndx, 1000);
	setHeats(pIndx, 1000);
	calcTempAndHeat(pIndx);						//�M�ʂ̉��x�ϊ��C���x�̔M�ʕϊ�

	//���M�ω��I��
	setTemps(pIndx, 1000);
	setHeats(pIndx, 1000);
	calcTempAndHeat(pIndx);						//�M�ʂ̉��x�ϊ��C���x�̔M�ʕϊ�
}

void HeatTransfar::WarmParticle(int pIndx, float temp, float heat)
{
	//���M�ω��I��
	float newTemp = getTemps()[pIndx] + temp;
	float newHeat = getHeats()[pIndx] + heat;

	setTemps(pIndx, newTemp);
	setHeats(pIndx, newHeat);
	calcTempAndHeat(pIndx);						//�M�ʂ̉��x�ϊ��C���x�̔M�ʕϊ�

	//���M�ω��I��
	calcTempAndHeat(pIndx);
}

void HeatTransfar::AddParticle(int nowVerticesNum)
{
	//���q�̉��x�C�M�ʁC��Ԃ̏�����
	for( int i = mNumVertices; i < nowVerticesNum; i++)
	{
		mTemps[i] = 1000.0f;
		mHeats[i] = 1000.0f;

		if( mTemps[i] < 250 )
		{
			mPhase[i]	 = -2;	//�X
			mPhaseChange[i] = 0;
		}
		else
		{
			mPhase[i]	 = 2;	//��
			mPhaseChange[i] = 0;
		}
	}

	mNumVertices = nowVerticesNum;
}

//����܂�����w�b�_�Ɉړ�
void HeatTransfar::resetNeighborhoodsId()
{
	mNeighborhoodsId.clear();
}

void HeatTransfar::resetNeighborhoodsDist()
{
	mNeighborhoodsDis.clear();
}

void HeatTransfar::AddNeighborhoodsId(std::vector<int> ids)
{
	mNeighborhoodsId.push_back( ids );
}

void HeatTransfar::AddNeighborhoodsDist(std::vector<float> dists)
{
	mNeighborhoodsDis.push_back( dists );
}

//�e�����a����J�[�l���֐��̒萔��ݒ肷��
void HeatTransfar::setCarnelConstant(float radius)
{
	// �J�[�l���֐��̒萔
	m_fWpoly6	=   4.0/(RX_PI*pow((double)radius, (double)8.0));
	m_fGWpoly6	= -24.0/(RX_PI*pow((double)radius, (double)8.0));
	m_fLWpoly6	= -24.0/(RX_PI*pow((double)radius, (double)8.0));

	m_fWspiky	=  10.0/(RX_PI*pow((double)radius, (double)5.0));
	m_fGWspiky	= -30.0/(RX_PI*pow((double)radius, (double)5.0));
	m_fLWspiky	= -60.0/(RX_PI*pow((double)radius, (double)5.0));

	m_fWvisc	= 10.0/(3.0*RX_PI*pow((double)radius, (double)2.0));
	m_fGWvisc	= 10.0/(3.0*RX_PI*pow((double)radius, (double)2.0));
	m_fLWvisc	= 20.0/(3.0*RX_PI*pow((double)radius, (double)5.0));

	m_fWspline  = 10.0/(7.0*RX_PI*pow((double)radius, (double)2.0));
	m_fGWspline = 15.0/(14.0*RX_PI*pow((double)radius, (double)4.0));
	m_fLWspline = 30.0/(7.0*RX_PI*pow((double)radius, (double)4.0));
}

/*!
 * Spline�J�[�l���֐��l�̌v�Z
 * @param[in] r ����
 * @param[in] h �L�����a
 * @return �֐��l
 */
inline double HeatTransfar::KernelSpline(const double &r, const double &h)
{
    double q = r/h;
    if(q >= 0.0 && q < 1.0)
	{
        return m_fWspline*(1.0-1.5*q*q+0.75*q*q*q);
    }
    else if(q >= 1.0 && q < 2.0)
	{
        return m_fWspline*0.25*(2.0-q)*(2.0-q)*(2.0-q);
    }
    else{
        return 0.0;
    }
}

/*!
 * Spline�J�[�l���֐����z�l�̌v�Z
 * @param[in] r ����
 * @param[in] rij ���Έʒu�x�N�g��
 * @param[in] h �L�����a
 * @return ���z�l
 */
/*
inline Vec2 HeatTransfar::KernelSplineG(const double &r, const Vec2 &rij, const double &h)
{
    double q = r/h;
    if(q >= 0.0 && q < 1.0){
        return  m_fGWspline*(0.75*q-1.0)*rij;
    }
    else if(q >= 1.0 && q < 2.0){
        return -m_fGWspline*(2.0-q)*(2.0-q)*rij/q;
    }
    else{
        return Vec2(0.0);
    }
}
 */

/*!
 * Spline�J�[�l���֐����v���V�A���̌v�Z
 * @param[in] r ����
 * @param[in] h �L�����a
 * @return ���v���V�A���̒l
 */
inline double HeatTransfar::KernelSplineL(const double &r, const double &h)
{
    double q = r/h;
    if(q >= 0.0 && q < 1.0)
	{
        return m_fLWspline*(1.5*q-1.0);
    }
    else if(q >= 1.0 && q < 2.0)
	{
        return m_fLWspline*0.5*(2.0-q);
    }
    else{
        return 0.0;
    }
}

//���x�ƔM�ʂ̏����@���M�E���M�̌v�Z�C���ω�����
void HeatTransfar::calcTempAndHeat(int i)
{
	//���ԏ�Ԃւ̕ω����o
	if( mPhase[i] == -2 && mTemps[i] > 250.0f )					//�X�̏ꍇ
	{
		mPhase[i] = -1;											//�X���ԏ��
		mTemps[i] = 250.0f;
		mHeats[i] = 0;
	}
	else if( mPhase[i] == 2 && mTemps[i] < 250.0f )				//���̏ꍇ
	{
		mPhase[i] = 1;											//�����ԏ��
		mTemps[i] = 250.0f;
		mHeats[i] = mLatentHeat;								//�Z����M
	}

	//���M�E���M�̌v�Z
	if( mPhase[i] == -1 || mPhase[i] == 1 )						//�����Ԃ��X���Ԃ̏ꍇ
	{
		//���M�v�Z
		//�ω����x��M�ʂɕϊ����Đ��M�v�Z�@�i���x�͕ω����Ȃ��j
		float spcfHt = 2.1f + (2.1f * mHeats[i] / mLatentHeat);	//��M�@�ܗL�M�ʂŐ��ƕX�̔�M����
		mHeats[i]	+= mTempsDelta[i] * spcfHt * 1.0f;			//���x�ω���M�ʂɊ��Z�@���ʂ͂P�ŌŒ�

//		cout << "mHeats[" << i << "] = " << mHeats[i] << "mPhase[i] = " << mPhase[i] << endl;1

		//���M�ω����猰�M�ω��֖߂锻��
		if( mPhase[i] == -1 && mHeats[i] < -mLatentHeat/300.0f )//���P�@<0�@�ɂ���ƁC�n���Ȃ��Ȃ�
		{
			//�X���ԏ�ԁ��X
			mPhase[i] = -2;										//�X�ւ̑��ω�
			mTemps[i] = 249.0f;
			mHeats[i] = 0;										//�M�i�|�e���V�����G�l���M�[�j����o
//			cout << "������֖߂���  i= " << i << endl;
		}
		else if( mPhase[i] == 1 && mHeats[i] > mLatentHeat )
		{
			//�����ԏ�ԁ���
			mPhase[i] = 2;										//���ւ̑��ω�
			mTemps[i] = 251.0f;
			mHeats[i] = 0;										//�M�i�|�e���V�����G�l���M�[�j���z��
//			cout << "�݂��֖߂��� i = " << i << endl;
		}

		//���ω�����
		if( mPhase[i] == -1 && mHeats[i] > mLatentHeat )		//�ܗL�M�ʂ��Z����M������
		{
			mPhase[i] = 2;										//���ւ̑��ω�
			mTemps[i] = 251.0f;
			mHeats[i] = 0;										//�M�i�|�e���V�����G�l���M�[�j���z��
			mPhaseChange[i] = 1;
//			cout << "�݂��� i = " << i << endl;
		}
		else if( mPhase[i] == 1 && mHeats[i] < 0.0f )			//�ܗL�M�ʂ��ÌŐ��M���g���؂�
		{
			mPhase[i] = -2;										//�X�ւ̑��ω�
			mTemps[i] = 249.0f;
			mHeats[i] = 0;										//�M�i�|�e���V�����G�l���M�[�j����o
			mPhaseChange[i] = 1;
//			cout << "�������  i= " << i << endl;
		}
	}
	else
	{	
		//���M�v�Z
		//�ω��M�ʂ̓K�p
		float spcfHt = (mTemps[i] > 250.0f)? 4.2f : 2.1f;		//��M�@���ƕX�Ŕ�M��ϓ��@��4.2�@�X2.1
		mTemps[i] =	mTemps[i] + mHeats[i] / (spcfHt * 1.0f);	//�M�ʂ����x�Ɋ��Z�@���ʂ͌Œ肵�Ă���
		mHeats[i] = 0.0f;										//�������@�t���[�����ƂɔM�ʁi���M�j�͒~�ς���Ȃ�

		//�ω����x��K�p
		mTemps[i] += mTempsDelta[i];
		if( mTemps[i] > mTempMax) mTemps[i] = mTempMax;
		if( mTemps[i] < mTempMin) mTemps[i] = mTempMin;
	}
}

//���x�ƔM�ʂ̏����@���M�E���M�̌v�Z�C���ω�����
void HeatTransfar::calcTempAndHeat()
{
	for( int i = 0; i < mNumVertices; i++)
	{
		//���ԏ�Ԃւ̕ω����o
		if( mPhase[i] == -2 && mTemps[i] > 250.0f )					//�X�̏ꍇ
		{
			mPhase[i] = -1;											//�X���ԏ��
			mTemps[i] = 250.0f;
			mHeats[i] = 0;
		}
		else if( mPhase[i] == 2 && mTemps[i] < 250.0f )				//���̏ꍇ
		{
			mPhase[i] = 1;											//�����ԏ��
			mTemps[i] = 250.0f;
			mHeats[i] = mLatentHeat;								//�Z����M
		}

		//���M�E���M�̌v�Z
		if( mPhase[i] == -1 || mPhase[i] == 1 )						//�����Ԃ��X���Ԃ̏ꍇ
		{
			//���M�v�Z
			//�ω����x��M�ʂɕϊ����Đ��M�v�Z�@�i���x�͕ω����Ȃ��j
			float spcfHt = 2.1f + (2.1f * mHeats[i] / mLatentHeat);	//��M�@�ܗL�M�ʂŐ��ƕX�̔�M����
			mHeats[i]	+= mTempsDelta[i] * spcfHt * 1.0f;			//���x�ω���M�ʂɊ��Z�@���ʂ͂P�ŌŒ�

//			cout << "mHeats[" << i << "] = " << mHeats[i] << "mPhase[i] = " << mPhase[i] << endl;1

			//���M�ω����猰�M�ω��֖߂锻��
			if( mPhase[i] == -1 && mHeats[i] < -mLatentHeat/300.0f )//���P�@<0�@�ɂ���ƁC�n���Ȃ��Ȃ�
			{
				//�X���ԏ�ԁ��X
				mPhase[i] = -2;										//�X�ւ̑��ω�
				mTemps[i] = 249.0f;
				mHeats[i] = 0;										//�M�i�|�e���V�����G�l���M�[�j����o
//				cout << "������֖߂���  i= " << i << endl;
			}
			else if( mPhase[i] == 1 && mHeats[i] > mLatentHeat )
			{
				//�����ԏ�ԁ���
				mPhase[i] = 2;										//���ւ̑��ω�
				mTemps[i] = 251.0f;
				mHeats[i] = 0;										//�M�i�|�e���V�����G�l���M�[�j���z��
//				cout << "�݂��֖߂��� i = " << i << endl;
			}

			//���ω�����
			if( mPhase[i] == -1 && mHeats[i] > mLatentHeat )		//�ܗL�M�ʂ��Z����M������
			{
				mPhase[i] = 2;										//���ւ̑��ω�
				mTemps[i] = 251.0f;
				mHeats[i] = 0;										//�M�i�|�e���V�����G�l���M�[�j���z��
				mPhaseChange[i] = 1;
//				cout << "�݂��� i = " << i << endl;
			}
			else if( mPhase[i] == 1 && mHeats[i] < 0.0f )			//�ܗL�M�ʂ��ÌŐ��M���g���؂�
			{
				mPhase[i] = -2;										//�X�ւ̑��ω�
				mTemps[i] = 249.0f;
				mHeats[i] = 0;										//�M�i�|�e���V�����G�l���M�[�j����o
				mPhaseChange[i] = 1;
//				cout << "�������  i= " << i << endl;
			}
		}
		else
		{	
			//���M�v�Z
			//�ω��M�ʂ̓K�p
			float spcfHt = (mTemps[i] > 250.0f)? 4.2f : 2.1f;		//��M�@���ƕX�Ŕ�M��ϓ��@��4.2�@�X2.1
			mTemps[i] =	mTemps[i] + mHeats[i] / (spcfHt * 1.0f);	//�M�ʂ����x�Ɋ��Z�@���ʂ͌Œ肵�Ă���
			mHeats[i] = 0.0f;										//�������@�t���[�����ƂɔM�ʁi���M�j�͒~�ς���Ȃ�

			//�ω����x��K�p
			mTemps[i] += mTempsDelta[i];
			if( mTemps[i] > mTempMax) mTemps[i] = mTempMax;
			if( mTemps[i] < mTempMin) mTemps[i] = mTempMin;
		}
	}
}

//��C�Ɨ��q�̔M����
void HeatTransfar::heatAirAndParticle()
{
	int *surfaceId = getSurfaceParticleNums();

	for( int i = 0; i < mNumVertices; i++ )
	{
		if( surfaceId[i] == -1 ) continue;
//		if( i >  2200 ) continue;
//		cout << "surcaceid = " << surfaceId[i] << endl;

		double airNum = 20.0-(double)surfaceId[i];
		if(airNum < 0) airNum = 0.0;

		double surfaceArea = airNum/20.0;						//��C�ƐG��Ă���\�ʐρ@0�`1.0 15�͓K��

		if( surfaceArea < 0.0) surfaceArea = 0.0;
		if( surfaceArea > 1.0) surfaceArea = 1.0;

		double qHeat = mHT * ( mAirTemp - mTemps[i])*surfaceArea;		//�j���[�g���̗�p�@���̎�����M�ʂ��v�Z
		mHeats[i] += qHeat;												//�M�ʂ����Z
//		cout << "i = " << i << "qHeat=" << qHeat << " mHT=" << mHT 
//			<< " mAirTemp=" << mAirTemp << " mTemps[" << i << "]=" << mTemps[i] << " surfaceArea=" << surfaceArea << endl;
	}
}

//���q���m�̔M����
void HeatTransfar::heatParticleAndParticle(const float* d, double h)	//d:���x�z��@h:�e�����a
{
	double tmp = 0.0;

	for( unsigned i = 0; i < mNeighborhoodsId.size(); i++)
	{
		tmp = 0.0;

		//�ߖT���q�̐������܂킷
		for( unsigned j = 0; j < mNeighborhoodsId[i].size(); j++)
		{
			int id = mNeighborhoodsId[i][j];
			double dis = mNeighborhoodsDis[i][j];
			if( i == id )	continue;							//�ߖT���q�Ɏ������܂܂�邱�Ƃ����邽��
			float densty = d[id];
			if( densty < 0.05f ) densty = 0.05f;				//���x������������or�O�̏ꍇ������̂Œ���
			tmp += timeStep * 1.0 * (mTemps[id] - mTemps[i]) / densty * KernelSpline(dis, h);	//�_�����Q�l�Ɂ@�Ƃ肠�������ʂP�ň��
		}

		mTempsDelta[i] = tmp * mTD;		//�M�g�U�W��
	}

}