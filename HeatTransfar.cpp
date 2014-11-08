//-----------------------------------------------------------------------------
// インクルードファイル
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

	mTemps = new float[mNumVertices];				//温度
	mTempsDelta = new float[mNumVertices];			//変化温度
	mHeats = new float[mNumVertices];				//熱量

	mPhase = new int[mNumVertices];					//現在の状態
	mPhaseChange = new int[mNumVertices];			//相変化を行うかのフラグ

	mSurfaceParticleNums = new int[mNumVertices];

	mNeighborhoodsId.clear();
	mNeighborhoodsDis.clear();

	//粒子の温度，熱量，状態の初期化
	for( int i = 0; i < mNumVertices; i++)
	{
		mTemps[i] = 0.0f;
		mHeats[i] = 0.0f;
		if( mTemps[i] < 250 )
		{
			mPhase[i]	 = -2;	//氷
			mPhaseChange[i] = 0;
		}
		else
		{
			mPhase[i]	 = 2;	//水
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
	//顕熱変化終了
	setTemps(pIndx, 1000);
	setHeats(pIndx, 1000);
	calcTempAndHeat(pIndx);						//熱量の温度変換，温度の熱量変換

	//潜熱変化終了
	setTemps(pIndx, 1000);
	setHeats(pIndx, 1000);
	calcTempAndHeat(pIndx);						//熱量の温度変換，温度の熱量変換
}

void HeatTransfar::WarmParticle(int pIndx, float temp, float heat)
{
	//顕熱変化終了
	float newTemp = getTemps()[pIndx] + temp;
	float newHeat = getHeats()[pIndx] + heat;

	setTemps(pIndx, newTemp);
	setHeats(pIndx, newHeat);
	calcTempAndHeat(pIndx);						//熱量の温度変換，温度の熱量変換

	//潜熱変化終了
	calcTempAndHeat(pIndx);
}

void HeatTransfar::AddParticle(int nowVerticesNum)
{
	//粒子の温度，熱量，状態の初期化
	for( int i = mNumVertices; i < nowVerticesNum; i++)
	{
		mTemps[i] = 1000.0f;
		mHeats[i] = 1000.0f;

		if( mTemps[i] < 250 )
		{
			mPhase[i]	 = -2;	//氷
			mPhaseChange[i] = 0;
		}
		else
		{
			mPhase[i]	 = 2;	//水
			mPhaseChange[i] = 0;
		}
	}

	mNumVertices = nowVerticesNum;
}

//じゃまだからヘッダに移動
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

//影響半径からカーネル関数の定数を設定する
void HeatTransfar::setCarnelConstant(float radius)
{
	// カーネル関数の定数
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
 * Splineカーネル関数値の計算
 * @param[in] r 距離
 * @param[in] h 有効半径
 * @return 関数値
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
 * Splineカーネル関数勾配値の計算
 * @param[in] r 距離
 * @param[in] rij 相対位置ベクトル
 * @param[in] h 有効半径
 * @return 勾配値
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
 * Splineカーネル関数ラプラシアンの計算
 * @param[in] r 距離
 * @param[in] h 有効半径
 * @return ラプラシアンの値
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

//温度と熱量の処理　顕熱・潜熱の計算，相変化判定
void HeatTransfar::calcTempAndHeat(int i)
{
	//中間状態への変化検出
	if( mPhase[i] == -2 && mTemps[i] > 250.0f )					//氷の場合
	{
		mPhase[i] = -1;											//氷中間状態
		mTemps[i] = 250.0f;
		mHeats[i] = 0;
	}
	else if( mPhase[i] == 2 && mTemps[i] < 250.0f )				//水の場合
	{
		mPhase[i] = 1;											//水中間状態
		mTemps[i] = 250.0f;
		mHeats[i] = mLatentHeat;								//融解潜熱
	}

	//顕熱・潜熱の計算
	if( mPhase[i] == -1 || mPhase[i] == 1 )						//水中間か氷中間の場合
	{
		//潜熱計算
		//変化温度を熱量に変換して潜熱計算　（温度は変化しない）
		float spcfHt = 2.1f + (2.1f * mHeats[i] / mLatentHeat);	//比熱　含有熱量で水と氷の比熱を補間
		mHeats[i]	+= mTempsDelta[i] * spcfHt * 1.0f;			//温度変化を熱量に換算　質量は１で固定

//		cout << "mHeats[" << i << "] = " << mHeats[i] << "mPhase[i] = " << mPhase[i] << endl;1

		//潜熱変化から顕熱変化へ戻る判定
		if( mPhase[i] == -1 && mHeats[i] < -mLatentHeat/300.0f )//改善　<0　にすると，溶けなくなる
		{
			//氷中間状態→氷
			mPhase[i] = -2;										//氷への相変化
			mTemps[i] = 249.0f;
			mHeats[i] = 0;										//熱（ポテンシャルエネルギー）を放出
//			cout << "こおりへ戻った  i= " << i << endl;
		}
		else if( mPhase[i] == 1 && mHeats[i] > mLatentHeat )
		{
			//水中間状態→水
			mPhase[i] = 2;										//水への相変化
			mTemps[i] = 251.0f;
			mHeats[i] = 0;										//熱（ポテンシャルエネルギー）を吸収
//			cout << "みずへ戻った i = " << i << endl;
		}

		//相変化判定
		if( mPhase[i] == -1 && mHeats[i] > mLatentHeat )		//含有熱量が融解潜熱を上回る
		{
			mPhase[i] = 2;										//水への相変化
			mTemps[i] = 251.0f;
			mHeats[i] = 0;										//熱（ポテンシャルエネルギー）を吸収
			mPhaseChange[i] = 1;
//			cout << "みずへ i = " << i << endl;
		}
		else if( mPhase[i] == 1 && mHeats[i] < 0.0f )			//含有熱量が凝固潜熱を使い切る
		{
			mPhase[i] = -2;										//氷への相変化
			mTemps[i] = 249.0f;
			mHeats[i] = 0;										//熱（ポテンシャルエネルギー）を放出
			mPhaseChange[i] = 1;
//			cout << "こおりへ  i= " << i << endl;
		}
	}
	else
	{	
		//顕熱計算
		//変化熱量の適用
		float spcfHt = (mTemps[i] > 250.0f)? 4.2f : 2.1f;		//比熱　水と氷で比熱を変動　水4.2　氷2.1
		mTemps[i] =	mTemps[i] + mHeats[i] / (spcfHt * 1.0f);	//熱量を温度に換算　質量は固定している
		mHeats[i] = 0.0f;										//初期化　フレームごとに熱量（顕熱）は蓄積されない

		//変化温度を適用
		mTemps[i] += mTempsDelta[i];
		if( mTemps[i] > mTempMax) mTemps[i] = mTempMax;
		if( mTemps[i] < mTempMin) mTemps[i] = mTempMin;
	}
}

//温度と熱量の処理　顕熱・潜熱の計算，相変化判定
void HeatTransfar::calcTempAndHeat()
{
	for( int i = 0; i < mNumVertices; i++)
	{
		//中間状態への変化検出
		if( mPhase[i] == -2 && mTemps[i] > 250.0f )					//氷の場合
		{
			mPhase[i] = -1;											//氷中間状態
			mTemps[i] = 250.0f;
			mHeats[i] = 0;
		}
		else if( mPhase[i] == 2 && mTemps[i] < 250.0f )				//水の場合
		{
			mPhase[i] = 1;											//水中間状態
			mTemps[i] = 250.0f;
			mHeats[i] = mLatentHeat;								//融解潜熱
		}

		//顕熱・潜熱の計算
		if( mPhase[i] == -1 || mPhase[i] == 1 )						//水中間か氷中間の場合
		{
			//潜熱計算
			//変化温度を熱量に変換して潜熱計算　（温度は変化しない）
			float spcfHt = 2.1f + (2.1f * mHeats[i] / mLatentHeat);	//比熱　含有熱量で水と氷の比熱を補間
			mHeats[i]	+= mTempsDelta[i] * spcfHt * 1.0f;			//温度変化を熱量に換算　質量は１で固定

//			cout << "mHeats[" << i << "] = " << mHeats[i] << "mPhase[i] = " << mPhase[i] << endl;1

			//潜熱変化から顕熱変化へ戻る判定
			if( mPhase[i] == -1 && mHeats[i] < -mLatentHeat/300.0f )//改善　<0　にすると，溶けなくなる
			{
				//氷中間状態→氷
				mPhase[i] = -2;										//氷への相変化
				mTemps[i] = 249.0f;
				mHeats[i] = 0;										//熱（ポテンシャルエネルギー）を放出
//				cout << "こおりへ戻った  i= " << i << endl;
			}
			else if( mPhase[i] == 1 && mHeats[i] > mLatentHeat )
			{
				//水中間状態→水
				mPhase[i] = 2;										//水への相変化
				mTemps[i] = 251.0f;
				mHeats[i] = 0;										//熱（ポテンシャルエネルギー）を吸収
//				cout << "みずへ戻った i = " << i << endl;
			}

			//相変化判定
			if( mPhase[i] == -1 && mHeats[i] > mLatentHeat )		//含有熱量が融解潜熱を上回る
			{
				mPhase[i] = 2;										//水への相変化
				mTemps[i] = 251.0f;
				mHeats[i] = 0;										//熱（ポテンシャルエネルギー）を吸収
				mPhaseChange[i] = 1;
//				cout << "みずへ i = " << i << endl;
			}
			else if( mPhase[i] == 1 && mHeats[i] < 0.0f )			//含有熱量が凝固潜熱を使い切る
			{
				mPhase[i] = -2;										//氷への相変化
				mTemps[i] = 249.0f;
				mHeats[i] = 0;										//熱（ポテンシャルエネルギー）を放出
				mPhaseChange[i] = 1;
//				cout << "こおりへ  i= " << i << endl;
			}
		}
		else
		{	
			//顕熱計算
			//変化熱量の適用
			float spcfHt = (mTemps[i] > 250.0f)? 4.2f : 2.1f;		//比熱　水と氷で比熱を変動　水4.2　氷2.1
			mTemps[i] =	mTemps[i] + mHeats[i] / (spcfHt * 1.0f);	//熱量を温度に換算　質量は固定している
			mHeats[i] = 0.0f;										//初期化　フレームごとに熱量（顕熱）は蓄積されない

			//変化温度を適用
			mTemps[i] += mTempsDelta[i];
			if( mTemps[i] > mTempMax) mTemps[i] = mTempMax;
			if( mTemps[i] < mTempMin) mTemps[i] = mTempMin;
		}
	}
}

//空気と粒子の熱処理
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

		double surfaceArea = airNum/20.0;						//空気と触れている表面積　0～1.0 15は適当

		if( surfaceArea < 0.0) surfaceArea = 0.0;
		if( surfaceArea > 1.0) surfaceArea = 1.0;

		double qHeat = mHT * ( mAirTemp - mTemps[i])*surfaceArea;		//ニュートンの冷却法則の式から熱量を計算
		mHeats[i] += qHeat;												//熱量を加算
//		cout << "i = " << i << "qHeat=" << qHeat << " mHT=" << mHT 
//			<< " mAirTemp=" << mAirTemp << " mTemps[" << i << "]=" << mTemps[i] << " surfaceArea=" << surfaceArea << endl;
	}
}

//粒子同士の熱処理
void HeatTransfar::heatParticleAndParticle(const float* d, double h)	//d:密度配列　h:影響半径
{
	double tmp = 0.0;

	for( unsigned i = 0; i < mNeighborhoodsId.size(); i++)
	{
		tmp = 0.0;

		//近傍粒子の数だけまわす
		for( unsigned j = 0; j < mNeighborhoodsId[i].size(); j++)
		{
			int id = mNeighborhoodsId[i][j];
			double dis = mNeighborhoodsDis[i][j];
			if( i == id )	continue;							//近傍粒子に自分が含まれることがあるため
			float densty = d[id];
			if( densty < 0.05f ) densty = 0.05f;				//密度が小さすぎるor０の場合があるので調整
			tmp += timeStep * 1.0 * (mTemps[id] - mTemps[i]) / densty * KernelSpline(dis, h);	//論文を参考に　とりあえず質量１で一定
		}

		mTempsDelta[i] = tmp * mTD;		//熱拡散係数
	}

}