//---------------------------------------------------------------------------
//このファイルのクラス・関数が二回以上コンパイルされるのを防ぐための処理（インクルードガード）
#ifndef HeatTransfarH
#define HeatTransfarH
//---------------------------------------------------------------------------

#include <cuda_runtime.h>

#include <vector>
#include "Math2d\math2d.h"

//GPU処理
extern void LaunchHeatTransferGPU
(
	float* heats,
	float* temps,
	float* dtemps,
	const int* surfParticles,
	float airTemp,
	float cffCnHt,
	int prtNum
	//const vector<vector<rxNeigh>>& neights,
	//const vector<int>& objNeight,
	//float floor,
	//float effRadius,
	//const float* pos,
	//const float* dens
);

class HeatTransfar
{
public:
	HeatTransfar(int num);
	~HeatTransfar(void);

	void InitGPU();

	void MeltParticle(int pIndx);
	void WarmParticle(int pIndx, float temp, float heat);

//----------------------------------------GPU----------------------------------------------
	float* getHostHeats(){	return sd_Heats;	}
	float* getHostTemps(){	return sd_Temps;	}
	float* getHostDTemps(){	return sd_DTemps;	}

	int* getHostPhase(){	return sd_Phase;	}
	int* getHostPhaseChangeFlag(){	return sd_PhaseChangeFlag;	}
//----------------------------------------GPU----------------------------------------------

	float* getTemps(){	return mTemps;	}								//各粒子の温度配列を取得
	float* getHeats(){	return mHeats;	}								//各粒子の熱量配列を取得
	int* getSurfaceParticleNums(){	return mSurfaceParticleNums;	}	//各表面粒子の番号を取得

	int getAirTemp() {	return mAirTemp;	}							//空気温度を取得

	float getTempMax(){	return mTempMax;	}							//上界（最大）温度を取得
	float getTempMin(){ return mTempMin;	}							//下界（最小）温度を取得

	float getCffCntHt(){ return mHT;	}
	float getCffCntTd(){ return mTD;	}

	float getLatentHeat(){ return mLatentHeat; }						//融解潜熱を取得

	int getPhase(int i){		return mPhase[i];	}					//現在の状態を返す（水・氷）
	int getPhaseChange(int i){	return mPhaseChange[i]; }

	void setNumVertices(int num){	mNumVertices = num;	};				//粒子の個数を設定
	
	void setTimeStep( float time ){ timeStep = time; };
	
	void setTempMax(float max){	mTempMax = max; }
	void setTempMin(float min){ mTempMin = min; }
	
	void setLatentHeat(float heat){	mLatentHeat = heat; }

	void setTemps(int nr, double tmp){	mTemps[nr] = tmp;	};			//ある添え字の温度を設定
	void setHeats(int nr, double heat){	mHeats[nr] = heat;	};									//ある添え字の熱量を設定
	void setAirTemp(float temp){	mAirTemp = temp;	};				//空気温度を設定

	void setCffCntHt(float h){ mHT = h; };
	void setCffCntTd(float h){ mTD = h;	};

	void setPhase(int i, int phase){ mPhase[i] = phase; }
	void setPhaseChange(int i, int phase){ mPhaseChange[i] = phase; }

	void setSurfaceParticleNums(int nr, int num){	mSurfaceParticleNums[nr] = num;	}	//ある表面粒子の添え字を格納

	void resetNeighborhoodsId();										//各近傍粒子の添え字を初期化
	void resetNeighborhoodsDist();										//各近傍粒子の距離を初期化

	void AddNeighborhoodsId(std::vector<int> ids);						//各近傍粒子の添え字を設定
	void AddNeighborhoodsDist(std::vector<float> dists);				//各近傍粒子の距離を設定

	void AddParticle(int nowVerticesNum);								//sph法で粒子が追加された際の処理

	void heatAirAndParticle();											//空気とパーティクルの熱処理
	void heatObjAndParticle(const std::vector<int>& neight);					//オブジェクトとパーティクルの熱処理
	void heatParticleAndParticle(const float* d, double h);				//パーティクル間の熱処理

	void calcTempAndHeat();												//熱量から温度を求める
	void calcTempAndHeat(int pIndx);									//1個の粒子だけ熱量から温度を求める

	// splineカーネル　カーネル関数
	inline double KernelSpline(const double &r, const double &h);
//	inline Vec2   KernelSplineG(const double &r, const Vec2 &rij, const double &h);
	inline double KernelSplineL(const double &r, const double &h);

	void setCarnelConstant(float radius);

private:
	void initState();						//初期化

//----------------------------------------GPU----------------------------------------------
	static float* sd_Heats;
	static float* sd_Temps;
	static float* sd_DTemps;

	static int* sd_Phase;
	static int* sd_PhaseChangeFlag;
//----------------------------------------GPU----------------------------------------------

	float timeStep;							//タイムステップ

	float mAirHeat;							//空気の熱量	未使用
	float mAirTemp;							//空気の温度

	float mTempMax;							//温度の上限値
	float mTempMin;							//温度の下限値

	float mLatentHeat;						//融解潜熱	readOnlyにする？

	double mHT;								//熱伝達係数	coefficient of heat transfer
	double mTD;								//熱拡散係数	coefficient of thermal diffusivity


	int mNumVertices;						//頂点数

	float *mTemps;							//温度
	float *mTempsDelta;						//変化温度
	float *mHeats;							//熱量

	int *mPhase;							//状態　-2＝氷　-1＝氷中間　1＝水中間　2＝水
	int *mPhaseChange;						//前フレームの状態　０＝相変化なし　１＝相変化あり

	int *mSurfaceParticleNums;				//表面粒子の近傍粒子数

	std::vector < std::vector<int> >	mNeighborhoodsId;				//各粒子の近傍となる粒子のID
	std::vector < std::vector<float> >	mNeighborhoodsDis;				//各粒子の近傍となる粒子の距離

	// カーネル関数の計算の際に用いられる定数係数
	double m_fWpoly6;				//!< Pory6カーネルの定数係数
	double m_fGWpoly6;				//!< Pory6カーネルの勾配の定数係数
	double m_fLWpoly6;				//!< Pory6カーネルのラプラシアンの定数係数
	double m_fWspiky;				//!< Spikyカーネルの定数係数
	double m_fGWspiky;				//!< Spikyカーネルの勾配の定数係数
	double m_fLWspiky;				//!< Spikyカーネルのラプラシアンの定数係数
	double m_fWvisc;				//!< Viscosityカーネルの定数係数
	double m_fGWvisc;				//!< Viscosityカーネルの勾配の定数係数
	double m_fLWvisc;				//!< Viscosityカーネルのラプラシアンの定数係数

	double m_fWspline;				//!< Splineカーネルの定数係数
	double m_fGWspline;				//!< Splineカーネルの勾配の定数係数
	double m_fLWspline;				//!< Splineカーネルのラプラシアンの定数係数

};
#endif