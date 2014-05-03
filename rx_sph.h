/*!
  @file rx_sph.h
	
  @brief SPH法
 
  @author Makoto Fujisawa
  @date 2008-10,2011-06
*/
// FILE --rx_sph.h--

#ifndef _RX_SPH_H_
#define _RX_SPH_H_


//-----------------------------------------------------------------------------
// インクルードファイル
//-----------------------------------------------------------------------------
#include "rx_sph_commons.h"

#include "rx_ps.h"			// パーティクルシステム基底クラス
#include "rx_nnsearch.h"	// グリッド分割による近傍探索

#include "rx_sph_solid.h"
#include "rx_wavelet_noise.h"

#include "rx_kernel.h"

#include "rx_cu_common.cuh"



//-----------------------------------------------------------------------------
// 定義
//-----------------------------------------------------------------------------
//#define RX_USE_BOUNDARY_PARTICLE	// 境界パーティクルの有無

//#define GL_REAL GL_DOUBLE
#define GL_REAL GL_FLOAT

// 時間計測
class rxTimerAvg;
extern rxTimerAvg g_Time;

#define RX_USE_TIMER

#ifdef RX_USE_TIMER
#define RXTIMER_CLEAR g_Time.ClearTime()
#define RXTIMER_RESET g_Time.ResetTime()
#define RXTIMER(x) g_Time.Split(x)
#define RXTIMER_PRINT g_Time.Print()
#define RXTIMER_STRING(x) g_Time.PrintToString(x)
#else
#define RXTIMER_CLEAR
#define RXTIMER_RESET
#define RXTIMER(x) 
#define RXTIMER_PRINT
#define RXTIMER_STRING(x)
#endif

const RXREAL KOL2 = (RXREAL)0.561231024;	// 2^(-5/6)

// グローバル変数の宣言
extern double g_fCoefEt;
extern double g_fMaxEt;
extern double g_fWaveletScale;
extern double g_fMaxEnergySpectrum;
extern double g_fMaxWaveletTurb;
extern double g_fESScale;
extern double g_fVCCoef;
extern double g_fCoefTurb;
extern double g_fCoefTurbForMesh;
extern double g_fEtCri;

extern double g_fNoiseScl;
extern double g_fNoiseEthr;
extern double g_fNoiseMag;


inline RXREAL3 MAKE_RXREAL3(RXREAL x, RXREAL y, RXREAL z)
{
	return make_float3(x, y, z);
}
inline RXREAL2 MAKE_RXREAL2(RXREAL x, RXREAL y)
{
	return make_float2(x, y);
}
inline RXREAL3 MAKE_FLOAT3V(Vec3 x)
{
	return make_float3((FLOAT)x[0], (FLOAT)x[1], (FLOAT)x[2]);
}



//! SPHシーンのパラメータ
struct rxSPHEnviroment
{
	#define MAX_DELETE_REGIONS 64

	int max_particles;			//!< 最大パーティクル数
	Vec3 boundary_cen;			//!< 境界の中心
	Vec3 boundary_ext;			//!< 境界の大きさ(各辺の長さの1/2)
	RXREAL dens;				//!< 初期密度
	RXREAL mass;				//!< パーティクルの質量
	RXREAL kernel_particles;	//!< 有効半径h以内のパーティクル数
	RXREAL dt;					//!< 時間ステップ幅
	RXREAL viscosity;			//!< 動粘性係数
	RXREAL gas_k;				//!< ガス定数

	int use_inlet;				//!< 流入境界条件の有無
	RXREAL et_cri;				//!< 乱流形成用の係数

	RXREAL epsilon;				//!< CFMの緩和係数
	RXREAL eta;					//!< 密度変動許容量
	int min_iter;				//!< ヤコビ反復最小数
	int max_iter;				//!< ヤコビ反復最大数

	int use_ap;					//!< 人工圧力ON/OFF (0 or 1)
	RXREAL ap_k;				//!< 人工圧力のための係数k (倍数)
	RXREAL ap_n;				//!< 人工圧力のための係数n (n乗)
	RXREAL ap_q;				//!< 人工圧力計算時の基準カーネル値計算用係数(有効半径hに対する係数, [0,1])

	int use_delete_region;		//!< パーティクル削除領域の有無(数)
	Vec3 delete_region[MAX_DELETE_REGIONS][2];	//!< 削除領域の範囲(最小，最大座標)

	// 表面メッシュ
	Vec3 mesh_boundary_cen;		//!< メッシュ生成境界の中心
	Vec3 mesh_boundary_ext;		//!< メッシュ生成境界の大きさ(各辺の長さの1/2)
	int mesh_vertex_store;		//!< 頂点数からポリゴン数を予測するときの係数
	int mesh_max_n;				//!< MC法用グリッドの最大分割数

	//追加：：ファイルから読み込んだパラメータ
	//熱処理
	float htTimeStep;
	float tempMax;
	float tempMin;

	float latentHeat;

	double cffcntHt;
	double cffcntTd;

	//SM法
	float smTimeStep;

	//氷構造
	int layer;

	rxSPHEnviroment()
	{
		max_particles = 50000;
		boundary_cen = Vec3(0.0);
		boundary_ext = Vec3(2.0, 0.8, 0.8);
		dens = (RXREAL)998.29;
		mass = (RXREAL)0.04;
		kernel_particles = (RXREAL)20.0;
		dt = 0.01;
		viscosity = 1.0e-3;
		gas_k = 3.0;

		mesh_vertex_store = 10;
		use_inlet = 0;
		mesh_max_n = 128;
		et_cri = 1.0;

		epsilon = 0.001;
		eta = 0.05;
		min_iter = 2;
		max_iter = 10;

		use_ap = true;
		ap_k = 0.1;
		ap_n = 4.0;
		ap_q = 0.2;

		use_delete_region = 0;
		for(int i = 0; i < MAX_DELETE_REGIONS; ++i){
			delete_region[i][0] = Vec3(0.0);
			delete_region[i][1] = Vec3(0.0);
		}
	}
};



//! 表面パーティクル
struct rxSurfaceParticle
{
	Vec3 pos;					//!< 中心座標
	Vec3 nrm;					//!< 法線
	Vec3 vel;					//!< 速度
	RXREAL d;					//!< 探索中心からの距離
	int idx;					//!< パーティクルインデックス
};

extern double g_fSurfThr[2];


//-----------------------------------------------------------------------------
// MARK:rxSPHクラスの宣言
//-----------------------------------------------------------------------------
class rxSPH : public rxParticleSystemBase
{
private:
	// パーティクル
	RXREAL *m_hNrm;					//!< パーティクル法線
	RXREAL *m_hFrc;					//!< パーティクルにかかる力
	RXREAL *m_hDens;				//!< パーティクル密度
	RXREAL *m_hPres;				//!< パーティクル圧力

	// 表面生成用(Anisotropic kernel)
	RXREAL *m_hUpPos;				//!< 平滑化パーティクル位置
	RXREAL *m_hPosW;				//!< 重み付き平均座標
	RXREAL *m_hCMatrix;				//!< 共分散行列
	RXREAL *m_hEigen;				//!< 共分散行列の特異値
	RXREAL *m_hRMatrix;				//!< 回転行列(共分散行列の特異ベクトル)
	RXREAL *m_hG;					//!< 変形行列

	uint *m_hSurf;					//!< 表面パーティクル

	// 境界・固体
	rxSolid *m_pBoundary;			//!< シミュレーション空間の境界
	vector<rxSolid*> m_vSolids;		//!< 固体物体
	RXREAL *m_hVrts;				//!< 固体ポリゴンの頂点
	int m_iNumVrts;					//!< 固体ポリゴンの頂点数
	int *m_hTris;					//!< 固体ポリゴン
	int m_iNumTris;					//!< 固体ポリゴンの数
	RXREAL *m_hSVels;				//!< 固体ポリゴン速度

	// 空間分割格子関連
	rxNNGrid *m_pNNGrid;			//!< 分割グリッドによる近傍探索
	vector< vector<rxNeigh> > m_vNeighs;	//!< 近傍パーティクル

	rxNNGrid *m_pNNGridB;			//!< 境界パーティクル用分割グリッド

	// 粒子パラメータ
	uint m_iKernelParticles;		//!< カーネル内のパーティクル数
	RXREAL m_fRestDens, m_fMass;	//!< 密度，質量
	RXREAL m_fEffectiveRadius;		//!< 有効半径
	RXREAL m_fKernelRadius;			//!< カーネルの影響範囲

	// シミュレーションパラメータ
	RXREAL m_fGasStiffness;			//!< ガス定数
	RXREAL m_fViscosity;			//!< 粘性係数
	RXREAL m_fBuoyancy;				//!< 浮力

	// カーネル関数の計算の際に用いられる定数係数
	double m_fAw;
	double m_fAg;
	double m_fAl;

	// カーネル関数
	double (*m_fpW)(double, double, double);
	Vec3 (*m_fpGW)(double, double, double, Vec3);
	double (*m_fpLW)(double, double, double, double);

protected:
	rxSPH(){}

public:
	//! コンストラクタ
	rxSPH(bool use_opengl);

	//! デストラクタ
	virtual ~rxSPH();

	// パーティクル半径
	float GetEffectiveRadius(){ return m_fEffectiveRadius; }

	// 近傍パーティクル
	uint* GetNeighborList(const int &i, int &n);

public:
	//
	// 仮想関数
	//
	// パーティクルデータ
	virtual RXREAL* GetParticle(void);
	virtual RXREAL* GetParticleDevice(void);
	virtual void UnmapParticle(void){}

	// シミュレーションステップ
	virtual bool Update(RXREAL dt, int step = 0);

	// シーンの設定
	virtual void SetPolygonObstacle(const vector<Vec3> &vrts, const vector<Vec3> &nrms, const vector< vector<int> > &tris, Vec3 vel);
	virtual void SetBoxObstacle(Vec3 cen, Vec3 ext, Vec3 ang, Vec3 vel, int flg);
	virtual void SetSphereObstacle(Vec3 cen, double rad, Vec3 vel, int flg);
	virtual void MoveSphereObstacle(int b, Vec3 disp);
	virtual Vec3 GetSphereObstaclePos(int b = -1);

	// ホスト<->VBO間転送
	virtual RXREAL* GetArrayVBO(rxParticleArray type, bool d2h = true, int num = -1);
	virtual void SetArrayVBO(rxParticleArray type, const RXREAL* data, int start, int count);
	virtual void SetArrayVBO(rxParticleArray type, const int* data, int start, int count);
	virtual void SetColorVBO(int type);

	// 陰関数値計算
	virtual double GetImplicit(double x, double y, double z);
	virtual void CalImplicitField(int n[3], Vec3 minp, Vec3 d, RXREAL *hF);
	virtual void CalImplicitFieldDevice(int n[3], Vec3 minp, Vec3 d, RXREAL *dF);
	//追加
	double GetImplicitSolid(double x, double y, double z);
	void CalImplicitFieldSolid(int n[3], Vec3 minp, Vec3 d, RXREAL *hF);
	void CalImplicitFieldDeviceSolid(int n[3], Vec3 minp, Vec3 d, RXREAL *dF, float *m_fIntrps);

	// SPH情報出力
	virtual void OutputSetting(string fn);

	// 描画関数
	virtual void DrawCell(int i, int j, int k);
	virtual void DrawCells(Vec3 col, Vec3 col2, int sel = 0);
	virtual void DrawObstacles(void);


public:
	// SPH初期化
	void Initialize(const rxSPHEnviroment &env);
	void Allocate(int max_particles);
	void Finalize(void);

	// 近傍取得
	void GetNearestNeighbors(Vec3 pos, vector<rxNeigh> &neighs, RXREAL h = -1.0);
	void GetNearestNeighbors(int idx, RXREAL *p, vector<rxNeigh> &neighs, RXREAL h = -1.0);

	void GetNearestNeighborsB(Vec3 pos, vector<rxNeigh> &neighs, RXREAL h = -1.0);

	// 分割セルにパーティクルを格納
	virtual void SetParticlesToCell(void);
	virtual void SetParticlesToCell(RXREAL *prts, int n, RXREAL h);

	// メタボールによる陰関数値
	double CalColorField(double x, double y, double z);
	//追加
	double CalColorFieldSolid(double x, double y, double z);

	// 分割セルにポリゴンを格納
	virtual void SetPolygonsToCell(void);

	int  GetPolygonsInCell(int gi, int gj, int gk, vector<int> &polys);
	bool IsPolygonsInCell(int gi, int gj, int gk);

	
	// 表面パーティクル検出
	void DetectSurfaceParticles(void);					// 表面パーティクルの検出
	double CalDistToNormalizedMassCenter(const int i);	// 近傍パーティクルの重心までの距離計算
	uint* GetArraySurf(void);							// 表面パーティクル情報の取得

	// 表面パーティクル情報の取得
	int GetSurfaceParticles(const Vec3 pos, RXREAL h, vector<rxSurfaceParticle> &sp);

	// 法線計算
	void CalNormalFromDensity(void);
	void CalNormal(void);

	//追加：：近傍粒子の取得
	vector<vector<rxNeigh>>& GetNeights(void){ return m_vNeighs; }

	//
	// Anisotropic kernel
	//
	virtual void CalAnisotropicKernel(void);

	//
	// SPS乱流
	//
	int	GetMaxSubLevel() const { return 0; }

	// サブパーティクルデータ(デバイスメモリ)
	RXREAL* GetSubParticleDev(void){ return 0; }
	void    UnmapSubParticle(void){}
	RXREAL* GetAllSubParticleDev(void){ return 0; }
	void    UnmapAllSubParticle(void){}
	RXREAL* GetSubParticlePosDev(void){ return 0; }
	RXREAL* GetSubParticleRadDev(void){ return 0; }
	RXREAL* GetSubParticleRatDev(void){ return 0; }

	int GetNumAllSubParticles(void){ return 0; }
	int GetNumValidSubParticles(void){ return 0; }

	// サブパーティクルデータ(ホストメモリ)
	RXREAL* GetSubParticlePos(void){ return 0; }
	RXREAL* GetSubParticleRad(void){ return 0; }
	RXREAL* GetSubParticleRat(void){ return 0; }
	unsigned int* GetSubParticleOcc(void){ return 0; }
	unsigned int* GetSubParticleOccScan(void){ return 0; }

	RXREAL GetSubRadius(int level){ return 0.0f; }

	//サブ粒子
	void SetEtcri(const double &et_cri){}
	void InitSubParticles(const double &et_cri){}
	void AddSubParticles(int start, int count){}
	void UpdateSubParticle(RXREAL *dPos, RXREAL *dSubPos, RXREAL scale, RXREAL dt){}
	void CalRadiusAndRatio(bool force = false){}

	void CalSubImplicitFieldDevice(int n[3], Vec3 minp, Vec3 d, RXREAL *dF, int turb){}
	void CalSubImplicitField(int n[3], Vec3 minp, Vec3 d, RXREAL *hF, const double &et_cri){}

	void AddSphereField(int n[3], Vec3 minp, Vec3 d, RXREAL *hF, 
						const Vec3 &pos, const double &eRad, const double& ratio){}



protected:
	// CPUによるSPH計算
	void calDensity(const RXREAL *pos, RXREAL *dens, RXREAL h);
	void calNormal(void);
	void calForce(const RXREAL *ppos, const RXREAL *pvel, const RXREAL *pdens, RXREAL *ppres, RXREAL *pfrc, RXREAL h);

	// rest densityの計算
	RXREAL calRestDensity(RXREAL h);

	// 個々の境界パーティクルの体積を計算
	void calBoundaryVolumes(const RXREAL *bpos, RXREAL *bvol, RXREAL mass, uint n, RXREAL h);

	// 時間ステップ幅の修正
	RXREAL calTimeStep(RXREAL &dt, RXREAL eta_avg, const RXREAL *pfrc, const RXREAL *pvel, const RXREAL *pdens);

	// 位置と速度の更新
	void integrate(const RXREAL *pos, const RXREAL *vel, const RXREAL *dens, const RXREAL *acc, 
				   RXREAL *pos_new, RXREAL *vel_new, RXREAL dt);

	// 衝突判定
	int calCollisionPolygon(uint grid_hash, Vec3 &pos0, Vec3 &pos1, Vec3 &vel, RXREAL dt);
	int calCollisionSolid(Vec3 &pos0, Vec3 &pos1, Vec3 &vel, RXREAL dt);
};


//-----------------------------------------------------------------------------
// MARK:rxSPH_GPUクラスの宣言
//-----------------------------------------------------------------------------
class rxSPH_GPU : public rxParticleSystemBase
{
private:
	//
	// メンバ変数(GPU変数)
	//
	RXREAL *m_dPos;			//!< パーティクル位置
	RXREAL *m_dVel;			//!< パーティクル速度
	RXREAL *m_dNrm;			//!< パーティクル法線
	RXREAL *m_dFrc;			//!< パーティクルにかかる力
	RXREAL *m_dDens;		//!< パーティクル密度
	RXREAL *m_dPres;		//!< パーティクル圧力

	RXREAL *m_dPosB;		//!< 境界パーティクル
	RXREAL *m_dVolB;		//!< 境界パーティクルの体積

	int *m_dAttr;			//!< パーティクル属性

	cudaGraphicsResource *m_pPosResource;	//!< OpenGL(のVBO)-CUDA間のデータ転送を扱うためのハンドル

	// ウェーブレット乱流
	RXREAL *m_dEt;			//!< 速度エネルギースペクトル
	RXREAL *m_dTurb;		//!< ウェーブレット乱流による速度場

	//サブ粒子 GPU data
	RXREAL *m_dSubPos;		//!< サブ粒子の絶対座標 float4(0.0,x,y,z)
	RXREAL *m_dSubChild;	//!< サブ粒子の子1への単位ベクトル
	RXREAL *m_dSubAxis;		//!< サブ粒子の回転軸(単位ベクトル)
	RXREAL *m_dSubEt;		//!< サブ粒子のエネルギースペクトル
	uint   *m_dSubRand;		//!< 乱数テーブル

	uint m_subposVBO;		//!< サブ粒子座標VBO
	cudaGraphicsResource *m_pSubPosResource;	//!< OpenGL(のVBO)-CUDA間のデータ転送を扱うためのハンドル

	// SPS乱流
	RXREAL *m_dUpPos;		//!< 平滑化パーティクル位置
	RXREAL *m_dPosW;		//!< 重み付き平均座標
	RXREAL *m_dCMatrix;		//!< 共分散行列
	RXREAL *m_dEigen;		//!< 共分散行列の特異値
	RXREAL *m_dRMatrix;		//!< 回転行列(共分散行列の特異ベクトル)
	RXREAL *m_dG;			//!< 変形行列

	// 表面メッシュ
	RXREAL *m_dVrts;		//!< 固体メッシュ頂点
	int    *m_dTris;		//!< 固体メッシュ

	// シミュレーションパラメータ
	rxSimParams m_params;	//!< シミュレーションパラメータ(GPUへのデータ渡し用)
	uint3 m_gridSize;		//!< 近傍探索グリッドの各軸の分割数
	uint m_numGridCells;	//!< 近傍探索グリッド分割総数
	
	// 空間分割(GPU)
	rxParticleCell m_dCellData;	//!< 近傍探索グリッド
	rxSubParticleCell m_dSubCellData;	//!< サブ粒子用近傍探索グリッド
	uint m_gridSortBits;		//!< ハッシュ値による基数ソート時の基数桁数

	rxParticleCell m_dCellDataB;	//!< 境界パーティクル用近傍探索グリッド
	uint3 m_gridSizeB;		//!< 近傍探索グリッドの各軸の分割数

	//
	// メンバ変数(CPU変数)
	//
	RXREAL *m_hNrm;			//!< パーティクル法線
	RXREAL *m_hFrc;			//!< パーティクルにかかる力
	RXREAL *m_hDens;		//!< パーティクル密度
	RXREAL *m_hPres;		//!< パーティクル圧力

	uint *m_hSurf;			//!< 表面パーティクル

	// 表面生成用(Anisotropic kernel)
	RXREAL *m_hUpPos;		//!< 平滑化パーティクル位置
	RXREAL *m_hPosW;		//!< 重み付き平均座標
	RXREAL *m_hCMatrix;		//!< 共分散行列
	RXREAL *m_hEigen;		//!< 共分散行列の特異値
	RXREAL *m_hRMatrix;		//!< 回転行列(共分散行列の特異ベクトル)
	RXREAL *m_hG;			//!< 変形行列
	RXREAL  m_fEigenMax;	//!< 特異値の最大値(探索半径拡張に用いる)

	// ウェーブレット乱流
	RXREAL *m_hVwt;			//!< パーティクル速度のウェーブレット変換
	RXREAL *m_hEt;			//!< 速度のエネルギースペクトラム
	RXREAL *m_hTurb;		//!< ウェーブレット乱流による速度場
	RXREAL *m_hNoiseTile[3];//!< ウェーブレットノイズタイル
	RXREAL m_fNoiseTileWidth;	//!< ノイズタイルの幅
	uint   m_iNoiseTileN[3];	//!< ノイズタイルの解像度

	// SPS乱流
	RXREAL *m_hSubPos;		//!< サブ粒子の絶対座標
	RXREAL *m_hSubChild;	//!< サブ粒子の子1への単位ベクトル
	RXREAL *m_hSubAxis;		//!< サブ粒子の回転軸(単位ベクトル)
	RXREAL *m_hSubEt;		//!< サブ粒子のエネルギースペクトル

	uint m_uNumSubParticles;//!< サブ粒子数
	uint m_uMaxSubParticles;//!< サブ粒子最大数
	uint m_uMaxSubLevel;	//!< サブ粒子最大レベル

	RXREAL *m_hSubPosPack;
	RXREAL *m_hSubRadPack;
	RXREAL *m_hSubRatPack;
	unsigned int *m_hSubOcc;
	unsigned int *m_hSubOccScan;
	bool m_bSubPacked;		//!< サブパーティクルをPackしたらtrue

	// 表面メッシュ
	vector<RXREAL> m_vVrts;	//!< 固体メッシュ頂点
	int m_iNumVrts;			//!< 固体メッシュ頂点数
	vector<int> m_vTris;	//!< 固体メッシュ
	int m_iNumTris;			//!< 固体メッシュ数

	bool m_bCalNormal;		//!< 法線計算フラグ

	//追加：近傍粒子のGPUからのコピーデータ
	uint* m_hSortedIndex;		//!< ハッシュ値でソートしたパーティクルインデックス
	uint* m_hGridParticleHash;	//!< 各パーティクルのグリッドハッシュ値
	uint* m_hCellStart;			//!< ソートリスト内の各セルのスタートインデックス
	uint* m_hCellEnd;			//!< ソートリスト内の各セルのエンドインデックス
	uint  m_uNumCells;			//!< 総セル数
	vector< vector<rxNeigh> > m_vNeighs;	//!< 近傍パーティクル

	// カーネル関数の計算の際に用いられる定数係数
	double m_fAw;
	double m_fAg;
	double m_fAl;

	// カーネル関数
	double (*m_fpW)(double, double, double);
	Vec3 (*m_fpGW)(double, double, double, Vec3);
	double (*m_fpLW)(double, double, double, double);

protected:
	rxSPH_GPU(){}

public:
	//! コンストラクタ
	rxSPH_GPU(bool use_opengl);

	//! デストラクタ
	~rxSPH_GPU();

	// パーティクル半径
	//追加：
	void SetEffectiveRadius(float r){	m_params.EffectiveRadius = r;	}
	float GetEffectiveRadius(){ return m_params.EffectiveRadius; }
	
	// 近傍パーティクル
	uint* GetNeighborList(const int &i, int &n);

	// シミュレーションパラメータ
	rxSimParams GetParams(void){ return m_params; }
	void UpdateParams(void);

	// フラグ切替
	void ToggleNormalCalc(int t = -1){ RX_TOGGLE(m_bCalNormal, t); }			//!< パーティクル法線の計算
	bool IsNormalCalc(void) const { return m_bCalNormal; }
#if MAX_BOX_NUM
	void ToggleSolidFlg(int t = -1);
#endif

public:
	//
	// 仮想関数
	//
	// パーティクルデータ
	virtual RXREAL* GetParticle(void);
	virtual RXREAL* GetParticleDevice(void);
	virtual void UnmapParticle(void);

	// シミュレーションステップ
	virtual bool Update(RXREAL dt, int step = 0);

	// シーンの設定
	virtual void SetPolygonObstacle(const vector<Vec3> &vrts, const vector<Vec3> &nrms, const vector< vector<int> > &tris, Vec3 vel);
	virtual void SetBoxObstacle(Vec3 cen, Vec3 ext, Vec3 ang, Vec3 vel, int flg);
	virtual void SetSphereObstacle(Vec3 cen, double rad, Vec3 vel, int flg);
	virtual void MoveSphereObstacle(int b, Vec3 disp);
	virtual Vec3 GetSphereObstaclePos(int b = -1);

	// ホスト<->VBO間転送
	virtual RXREAL* GetArrayVBO(rxParticleArray type, bool d2h = true, int num = -1);
	virtual void SetArrayVBO(rxParticleArray type, const RXREAL* data, int start, int count);
	virtual void SetArrayVBO(rxParticleArray type, const int* data, int start, int count);
	virtual void SetColorVBO(int type);


	// 陰関数値計算
	virtual double GetImplicit(double x, double y, double z);
	virtual void CalImplicitField(int n[3], Vec3 minp, Vec3 d, RXREAL *hF);
	virtual void CalImplicitFieldDevice(int n[3], Vec3 minp, Vec3 d, RXREAL *dF);

	//追加
	double GetImplicitSolid(double x, double y, double z, float* fIceCheck);
	void CalImplicitFieldSolid(int n[3], Vec3 minp, Vec3 d, RXREAL *hF, float* fIceCheck);
	void CalImplicitFieldDeviceSolid(int n[3], Vec3 minp, Vec3 d, RXREAL *dF, float* fIceCheck);

	// SPH情報出力
	virtual void OutputSetting(string fn);

	// 描画関数
	virtual void DrawCell(int i, int j, int k);
	virtual void DrawCells(Vec3 col, Vec3 col2, int sel = 0);
	virtual void DrawObstacles(void);

protected:
	void setObjectToCell(RXREAL *p);


public:
	// SPH初期化
	void Initialize(const rxSPHEnviroment &env);
	void Allocate(int max_particles);
	void Finalize(void);

	// 分割セルにパーティクルを格納
	virtual void SetParticlesToCell(void);
	virtual void SetParticlesToCell(RXREAL *prts, int n, RXREAL h);
	virtual void SetPolygonsToCell(void);

	int  GetPolygonsInCell(int gi, int gj, int gk, vector<int> &polys){ return 0; }
	bool IsPolygonsInCell(int gi, int gj, int gk){ return false; }
	
	//追加：：カラーフィールド値計算？
	double CalColorFieldSolid(double x, double y, double z, float* fIceCheck);
	
	void CalMaxDensity(int k);

	// 表面パーティクル検出
	void DetectSurfaceParticles(void);					// 表面パーティクルの検出
	double CalDistToNormalizedMassCenter(const int i);	// 近傍パーティクルの重心までの距離計算
	uint* GetArraySurf(void);							// 表面パーティクル情報の取得

	// 表面パーティクル情報の取得
	int GetSurfaceParticles(const Vec3 pos, RXREAL h, vector<rxSurfaceParticle> &sp);

	// 法線計算
	void CalNormalFromDensity(void);
	void CalNormal(void);

	// 境界パーティクルの初期化
	void InitBoundary(void);

	//追加：：近傍粒子取得
	vector<vector<rxNeigh>>& GetNeights(void){ return m_vNeighs; }

	//
	// Anisotropic kernel
	//
	virtual void CalAnisotropicKernel(void);

	//
	// ウェーブレット乱流
	//
	void CalParticleCWT(RXREAL scale);
	void CalWaveletTurbulence(RXREAL scale, RXREAL *dPos, RXREAL dt);

	//
	// SPS乱流
	//
	int	GetMaxSubLevel() const { return m_uMaxSubLevel; }

	// サブパーティクルデータ(デバイスメモリ)
	RXREAL* GetSubParticleDev(void);
	void    UnmapSubParticle(void);
	RXREAL* GetAllSubParticleDev(void);
	void    UnmapAllSubParticle(void);
	RXREAL* GetSubParticlePosDev(void);
	RXREAL* GetSubParticleRadDev(void);
	RXREAL* GetSubParticleRatDev(void);

	int GetNumAllSubParticles(void);
	int GetNumValidSubParticles(void);

	// サブパーティクルデータ(ホストメモリ)
	RXREAL* GetSubParticlePos(void);
	RXREAL* GetSubParticleRad(void);
	RXREAL* GetSubParticleRat(void);
	unsigned int* GetSubParticleOcc(void);
	unsigned int* GetSubParticleOccScan(void);

	RXREAL GetSubRadius(int level){ return m_dSubCellData.fSubRad[level]; }

	//サブ粒子
	void SetEtcri(const double &et_cri){ m_dSubCellData.fEtcri = et_cri;}
	void InitSubParticles(const double &et_cri);
	void AddSubParticles(int start, int count);
	void UpdateSubParticle(RXREAL *dPos, RXREAL *dSubPos, RXREAL scale, RXREAL dt);
	void CalRadiusAndRatio(bool force = false);

	void CalSubImplicitFieldDevice(int n[3], Vec3 minp, Vec3 d, RXREAL *dF, int turb);
	void CalSubImplicitField(int n[3], Vec3 minp, Vec3 d, RXREAL *hF, const double &et_cri);

	void AddSphereField(int n[3], Vec3 minp, Vec3 d, RXREAL *hF, 
						const Vec3 &pos, const double &eRad, const double& ratio);

protected:
	// rest densityの計算
	RXREAL calRestDensity(RXREAL h);

	// 分割セルの初期設定
	void setupCells(rxParticleCell &cell, uint3 &gridsize, double &cell_width, Vec3 vMin, Vec3 vMax, double h);

	// グリッドハッシュの計算
	uint calGridHash(uint x, uint y, uint z);
	uint calGridHash(Vec3 pos);

	// ポリゴンを分割セルに格納
	void setPolysToCell(RXREAL *vrts, int nv, int* tris, int nt);

	// CPU用の近傍粒子探索
	void searchNeighbors(void);
	void getNeighborsInCell(Vec3 pos, RXREAL *p, int gi, int gj, int gk, vector<rxNeigh> &neighs, RXREAL h);
};




//-----------------------------------------------------------------------------
// MARK:rxDDSPHクラスの宣言
//  - Double Density Relaxation によるSPH
//  - S.Clavet et al., "Particle-based Viscoelastic Fluid Simulation", SCA2005, 2005. 
//  - http://www.iro.umontreal.ca/labs/infographie/papers/Clavet-2005-PVFS/index.html
//-----------------------------------------------------------------------------
class rxDDSPH : public rxParticleSystemBase
{
	//! 粒子間バネ
	class rxSPHSpring
	{
	public:
		bool enable;
		double L;
		double r;
		int pi, pj;

		rxSPHSpring()
		{
			enable = false;
			L = r = 0.0;
		}

		~rxSPHSpring()
		{
		}
	};

protected:
	// パーティクル
	RXREAL *m_hNrm;					//!< パーティクル法線
	RXREAL *m_hFrc;					//!< パーティクルにかかる力
	RXREAL *m_hDens;				//!< パーティクル密度
	RXREAL *m_hPres;				//!< パーティクル圧力
	RXREAL *m_hVelOld;				//!< 前ステップ速度

	RXREAL *m_hPredictPos;			//!< 予測位置
	RXREAL *m_hDist;				//!< パーティクル間距離

	uint *m_hSurf;					//!< 表面パーティクル


	// 境界・固体
	rxSolid *m_pBoundary;			//!< シミュレーション空間の境界
	vector<rxSolid*> m_vSolids;		//!< 固体物体
	RXREAL *m_hVrts;				//!< 固体ポリゴンの頂点
	int m_iNumVrts;					//!< 固体ポリゴンの頂点数
	int *m_hTris;					//!< 固体ポリゴン
	int m_iNumTris;					//!< 固体ポリゴンの数

	// 近傍粒子探索
	rxNNGrid *m_pNNGrid;			//!< 分割グリッドによる近傍探索
	vector< vector<rxNeigh> > m_vNeighs;	//!< 近傍パーティクル

	// 粒子パラメータ
	uint m_iKernelParticles;		//!< カーネル内のパーティクル数
	RXREAL m_fMass;					//!< 質量
	RXREAL m_fEffectiveRadius;		//!< 有効半径

	// シミュレーションパラメータ
	RXREAL m_fBuoyancy;				//!< 浮力

	// Double Density用パラメータ
	RXREAL m_fK;					//!< ガス定数
	RXREAL m_fKnear;				//!< near density用ガス定数
	RXREAL m_fViscC;
	RXREAL m_fViscBeta;
	RXREAL m_fElasKspring;
	RXREAL m_fPlasAlpha;
	RXREAL m_fPlasGamma;

	// 粘弾性力用粒子間スプリング
	map<int, rxSPHSpring> m_mapSpring;

	// 密度
	double m_fInitDens;				//!< 初期平均密度
	double m_fInitDensNear;			//!< 初期平均密度(near)
	double m_fInitMaxDens;			//!< 初期最大密度

	double m_fMaxDens;				//!< 最大密度
	double m_fSurfDens;				//!< 表面粒子とする密度の比率

	bool m_bCalDens;				//!< rest density 計算フラグ

	// カーネル関数の計算の際に用いられる定数係数
	double m_fWpoly6;				//!< Pory6カーネルの定数係数
	double m_fGWpoly6;				//!< Pory6カーネルの勾配の定数係数
	double m_fLWpoly6;				//!< Pory6カーネルのラプラシアンの定数係数

protected:
	rxDDSPH(){}

public:
	//! コンストラクタ
	rxDDSPH(bool use_opengl);

	//! デストラクタ
	virtual ~rxDDSPH();

	// パーティクル半径
	float GetEffectiveRadius(){ return m_fEffectiveRadius; }

	// 近傍パーティクル
	uint* GetNeighborList(const int &i, int &n);

public:
	//
	// 仮想関数
	//
	// パーティクルデータ
	virtual RXREAL* GetParticle(void);
	virtual RXREAL* GetParticleDevice(void);
	virtual void UnmapParticle(void){}

	// シミュレーションステップ
	virtual bool Update(RXREAL dt, int step = 0);

	// シーンの設定
	virtual void SetPolygonObstacle(const vector<Vec3> &vrts, const vector<Vec3> &nrms, const vector< vector<int> > &tris, Vec3 vel);
	virtual void SetBoxObstacle(Vec3 cen, Vec3 ext, Vec3 ang, Vec3 vel, int flg);
	virtual void SetSphereObstacle(Vec3 cen, double rad, Vec3 vel, int flg);
	virtual void MoveSphereObstacle(int b, Vec3 disp);
	virtual Vec3 GetSphereObstaclePos(int b = -1);

	// ホスト<->VBO間転送
	virtual RXREAL* GetArrayVBO(rxParticleArray type, bool d2h = true, int num = -1);
	virtual void SetArrayVBO(rxParticleArray type, const RXREAL* data, int start, int count);
	virtual void SetArrayVBO(rxParticleArray type, const int* data, int start, int count){}
	virtual void SetColorVBO(int type);
	
	// 陰関数値計算
	virtual double GetImplicit(double x, double y, double z);
	virtual void CalImplicitField(int n[3], Vec3 minp, Vec3 d, RXREAL *hF);
	virtual void CalImplicitFieldDevice(int n[3], Vec3 minp, Vec3 d, RXREAL *dF);

	// SPH情報出力
	virtual void OutputSetting(string fn);

	// 描画関数
	virtual void DrawCell(int i, int j, int k);
	virtual void DrawCells(Vec3 col, Vec3 col2, int sel = 0);
	virtual void DrawObstacles(void);

	
public:
	// SPH初期化
	void Initialize(const rxSPHEnviroment &env);
	void Allocate(int max_particles);
	void Finalize(void);

	// 近傍取得
	void GetNearestNeighbors(Vec3 pos, vector<rxNeigh> &neighs, RXREAL h = -1.0);
	void GetNearestNeighbors(int idx, RXREAL *p, vector<rxNeigh> &neighs, RXREAL h = -1.0);

	// 分割セルにパーティクルを格納
	virtual void SetParticlesToCell(void);
	virtual void SetParticlesToCell(RXREAL *prts, int n, RXREAL h);

	// メタボールによる陰関数値
	double CalColorField(double x, double y, double z);

	// 分割セルにポリゴンを格納
	virtual void SetPolygonsToCell(void);

	int  GetPolygonsInCell(int gi, int gj, int gk, vector<int> &polys);
	bool IsPolygonsInCell(int gi, int gj, int gk);

	// 表面パーティクル検出
	void DetectSurfaceParticles(void);					// 表面パーティクルの検出
	double CalDistToNormalizedMassCenter(const int i);	// 近傍パーティクルの重心までの距離計算
	uint* GetArraySurf(void);							// 表面パーティクル情報の取得

	// 表面パーティクル情報の取得
	int GetSurfaceParticles(const Vec3 pos, RXREAL h, vector<rxSurfaceParticle> &sp);

	// 法線計算
	void CalNormalFromDensity(void);
	void CalNormal(void);

	//追加：：近傍粒子の取得
	vector<vector<rxNeigh>>& GetNeights(void){ return m_vNeighs; }

protected:
	// Double density relaxation method
	void calExternalForceToVel(RXREAL dt);
	void calDoubleDensity(RXREAL *ppos);
	void calDoubleDensity(RXREAL *ppos, double &avg_dens, double &avg_dens_near, double &max_dens);
	void calDoubleDensityRelaxation(RXREAL *ppos, RXREAL dt);
	void adjustSprings(RXREAL dt);
	void applySpringDisplacements(RXREAL dt);
	
	// 衝突判定
	int calCollisionPolygon(uint grid_hash, Vec3 &pos0, Vec3 &pos1, Vec3 &vel, RXREAL dt);
	int calCollisionSolid(Vec3 &pos0, Vec3 &pos1, Vec3 &vel, RXREAL dt);

public:
	//
	// Anisotropic kernel - DD用は未実装
	//
	virtual void CalAnisotropicKernel(void){}

	//
	// SPS乱流 - CPU用は未実装
	//
	int	GetMaxSubLevel() const { return 0; }

	// サブパーティクルデータ(デバイスメモリ)
	RXREAL* GetSubParticleDev(void){ return 0; }
	void    UnmapSubParticle(void){}
	RXREAL* GetAllSubParticleDev(void){ return 0; }
	void    UnmapAllSubParticle(void){}
	RXREAL* GetSubParticlePosDev(void){ return 0; }
	RXREAL* GetSubParticleRadDev(void){ return 0; }
	RXREAL* GetSubParticleRatDev(void){ return 0; }

	int GetNumAllSubParticles(void){ return 0; }
	int GetNumValidSubParticles(void){ return 0; }

	// サブパーティクルデータ(ホストメモリ)
	RXREAL* GetSubParticlePos(void){ return 0; }
	RXREAL* GetSubParticleRad(void){ return 0; }
	RXREAL* GetSubParticleRat(void){ return 0; }
	unsigned int* GetSubParticleOcc(void){ return 0; }
	unsigned int* GetSubParticleOccScan(void){ return 0; }

	RXREAL GetSubRadius(int level){ return 0.0f; }

	//サブ粒子
	void SetEtcri(const double &et_cri){}
	void InitSubParticles(const double &et_cri){}
	void AddSubParticles(int start, int count){}
	void UpdateSubParticle(RXREAL *dPos, RXREAL *dSubPos, RXREAL scale, RXREAL dt){}
	void CalRadiusAndRatio(bool force = false){}

	void CalSubImplicitFieldDevice(int n[3], Vec3 minp, Vec3 d, RXREAL *dF, int turb){}
	void CalSubImplicitField(int n[3], Vec3 minp, Vec3 d, RXREAL *hF, const double &et_cri){}

	void AddSphereField(int n[3], Vec3 minp, Vec3 d, RXREAL *hF, 
						const Vec3 &pos, const double &eRad, const double& ratio){}
};




//-----------------------------------------------------------------------------
// MARK:rxPBDSPHクラスの宣言
//  - Miles Macklin and Matthias Muller, "Position Based Fluids", Proc. SIGGRAPH 2013, 2013. 
//  - http://blog.mmacklin.com/publications/
//-----------------------------------------------------------------------------
class rxPBDSPH : public rxParticleSystemBase
{
private:
	// パーティクル
	RXREAL *m_hNrm;					//!< パーティクル法線
	RXREAL *m_hFrc;					//!< パーティクルにかかる力
	RXREAL *m_hDens;				//!< パーティクル密度
	RXREAL *m_hPres;				//!< パーティクル圧力

	RXREAL *m_hS;					//!< Scaling factor for CFM
	RXREAL *m_hDp;					//!< 位置修正量

	RXREAL *m_hPredictPos;			//!< 予測位置
	RXREAL *m_hPredictVel;			//!< 予測速度

	RXREAL *m_hSb;					//!< 境界パーティクルのScaling factor

	// 表面生成用(Anisotropic kernel)
	RXREAL *m_hUpPos;				//!< 平滑化パーティクル位置
	RXREAL *m_hPosW;				//!< 重み付き平均座標
	RXREAL *m_hCMatrix;				//!< 共分散行列
	RXREAL *m_hEigen;				//!< 共分散行列の特異値
	RXREAL *m_hRMatrix;				//!< 回転行列(共分散行列の特異ベクトル)
	RXREAL *m_hG;					//!< 変形行列

	uint *m_hSurf;					//!< 表面パーティクル

	// 境界・固体
	rxSolid *m_pBoundary;			//!< シミュレーション空間の境界
	vector<rxSolid*> m_vSolids;		//!< 固体物体
	RXREAL *m_hVrts;				//!< 固体ポリゴンの頂点
	int m_iNumVrts;					//!< 固体ポリゴンの頂点数
	int *m_hTris;					//!< 固体ポリゴン
	int m_iNumTris;					//!< 固体ポリゴンの数
	RXREAL *m_hSVels;				//!< 固体ポリゴン速度

	// 空間分割格子関連
	rxNNGrid *m_pNNGrid;			//!< 分割グリッドによる近傍探索
	vector< vector<rxNeigh> > m_vNeighs;	//!< 近傍パーティクル

	rxNNGrid *m_pNNGridB;			//!< 境界パーティクル用分割グリッド


	// 粒子パラメータ
	uint m_iKernelParticles;		//!< カーネル内のパーティクル数
	RXREAL m_fRestDens, m_fMass;	//!< 密度，質量
	RXREAL m_fEffectiveRadius;		//!< 有効半径
	RXREAL m_fKernelRadius;			//!< カーネルの影響範囲

	// シミュレーションパラメータ
	RXREAL m_fGasStiffness;			//!< ガス定数
	RXREAL m_fViscosity;			//!< 粘性係数
	RXREAL m_fBuoyancy;				//!< 浮力

	RXREAL m_fEpsilon;				//!< CFMの緩和係数
	RXREAL m_fEta;					//!< 密度変動率
	int m_iMinIterations;			//!< ヤコビ反復最小反復回数
	int m_iMaxIterations;			//!< ヤコビ反復最大反復回数

	bool m_bArtificialPressure;		//!< クラスタリングを防ぐためのArtificial Pressure項を追加するフラグ
	RXREAL m_fApK;					//!< 人工圧力のための係数k
	RXREAL m_fApN;					//!< 人工圧力のための係数n
	RXREAL m_fApQ;					//!< 人工圧力計算時の基準カーネル値計算用係数(有効半径hに対する係数, [0,1])


	// カーネル関数の計算の際に用いられる定数係数
	double m_fAw;
	double m_fAg;
	double m_fAl;

	// カーネル関数
	double (*m_fpW)(double, double, double);
	Vec3 (*m_fpGW)(double, double, double, Vec3);
	double (*m_fpLW)(double, double, double, double);

protected:
	rxPBDSPH(){}

public:
	//! コンストラクタ
	rxPBDSPH(bool use_opengl);

	//! デストラクタ
	virtual ~rxPBDSPH();

	// パーティクル半径
	float GetEffectiveRadius(){ return m_fEffectiveRadius; }

	// 近傍パーティクル
	uint* GetNeighborList(const int &i, int &n);

public:
	//
	// 仮想関数
	//
	// パーティクルデータ
	virtual RXREAL* GetParticle(void);
	virtual RXREAL* GetParticleDevice(void);
	virtual void UnmapParticle(void){}

	// シミュレーションステップ
	virtual bool Update(RXREAL dt, int step = 0);

	// シーンの設定
	virtual void SetPolygonObstacle(const vector<Vec3> &vrts, const vector<Vec3> &nrms, const vector< vector<int> > &tris, Vec3 vel);
	virtual void SetBoxObstacle(Vec3 cen, Vec3 ext, Vec3 ang, Vec3 vel, int flg);
	virtual void SetSphereObstacle(Vec3 cen, double rad, Vec3 vel, int flg);
	virtual void MoveSphereObstacle(int b, Vec3 disp);
	virtual Vec3 GetSphereObstaclePos(int b = -1);

	// ホスト<->VBO間転送
	virtual RXREAL* GetArrayVBO(rxParticleArray type, bool d2h = true, int num = -1);
	virtual void SetArrayVBO(rxParticleArray type, const RXREAL* data, int start, int count);
	virtual void SetArrayVBO(rxParticleArray type, const int* data, int start, int count);
	virtual void SetColorVBO(int type);

	// 陰関数値計算
	virtual double GetImplicit(double x, double y, double z);
	virtual void CalImplicitField(int n[3], Vec3 minp, Vec3 d, RXREAL *hF);
	virtual void CalImplicitFieldDevice(int n[3], Vec3 minp, Vec3 d, RXREAL *dF);

	// SPH情報出力
	virtual void OutputSetting(string fn);

	// 描画関数
	virtual void DrawCell(int i, int j, int k);
	virtual void DrawCells(Vec3 col, Vec3 col2, int sel = 0);
	virtual void DrawObstacles(void);


public:
	// SPH初期化
	void Initialize(const rxSPHEnviroment &env);
	void Allocate(int max_particles);
	void Finalize(void);

	// 近傍取得
	void GetNearestNeighbors(Vec3 pos, vector<rxNeigh> &neighs, RXREAL h = -1.0);
	void GetNearestNeighbors(int idx, RXREAL *p, vector<rxNeigh> &neighs, RXREAL h = -1.0);

	void GetNearestNeighborsB(Vec3 pos, vector<rxNeigh> &neighs, RXREAL h = -1.0);

	// 分割セルにパーティクルを格納
	virtual void SetParticlesToCell(void);
	virtual void SetParticlesToCell(RXREAL *prts, int n, RXREAL h);

	// メタボールによる陰関数値
	double CalColorField(double x, double y, double z);

	// 分割セルにポリゴンを格納
	virtual void SetPolygonsToCell(void);

	int  GetPolygonsInCell(int gi, int gj, int gk, vector<int> &polys);
	bool IsPolygonsInCell(int gi, int gj, int gk);

	
	// 表面パーティクル検出
	void DetectSurfaceParticles(void);					// 表面パーティクルの検出
	double CalDistToNormalizedMassCenter(const int i);	// 近傍パーティクルの重心までの距離計算
	uint* GetArraySurf(void);							// 表面パーティクル情報の取得

	// 表面パーティクル情報の取得
	int GetSurfaceParticles(const Vec3 pos, RXREAL h, vector<rxSurfaceParticle> &sp);

	// 法線計算
	void CalNormalFromDensity(void);
	void CalNormal(void);

	//追加：：近傍粒子の取得
	vector<vector<rxNeigh>>& GetNeights(void){ return m_vNeighs; }

	// 人口圧力項
	bool& GetArtificialPressure(void){ return m_bArtificialPressure; }

	//
	// Anisotropic kernel
	//
	virtual void CalAnisotropicKernel(void);

	//
	// SPS乱流
	//
	int	GetMaxSubLevel() const { return 0; }

	// サブパーティクルデータ(デバイスメモリ)
	RXREAL* GetSubParticleDev(void){ return 0; }
	void    UnmapSubParticle(void){}
	RXREAL* GetAllSubParticleDev(void){ return 0; }
	void    UnmapAllSubParticle(void){}
	RXREAL* GetSubParticlePosDev(void){ return 0; }
	RXREAL* GetSubParticleRadDev(void){ return 0; }
	RXREAL* GetSubParticleRatDev(void){ return 0; }

	int GetNumAllSubParticles(void){ return 0; }
	int GetNumValidSubParticles(void){ return 0; }

	// サブパーティクルデータ(ホストメモリ)
	RXREAL* GetSubParticlePos(void){ return 0; }
	RXREAL* GetSubParticleRad(void){ return 0; }
	RXREAL* GetSubParticleRat(void){ return 0; }
	unsigned int* GetSubParticleOcc(void){ return 0; }
	unsigned int* GetSubParticleOccScan(void){ return 0; }

	RXREAL GetSubRadius(int level){ return 0.0f; }

	//サブ粒子
	void SetEtcri(const double &et_cri){}
	void InitSubParticles(const double &et_cri){}
	void AddSubParticles(int start, int count){}
	void UpdateSubParticle(RXREAL *dPos, RXREAL *dSubPos, RXREAL scale, RXREAL dt){}
	void CalRadiusAndRatio(bool force = false){}

	void CalSubImplicitFieldDevice(int n[3], Vec3 minp, Vec3 d, RXREAL *dF, int turb){}
	void CalSubImplicitField(int n[3], Vec3 minp, Vec3 d, RXREAL *hF, const double &et_cri){}

	void AddSphereField(int n[3], Vec3 minp, Vec3 d, RXREAL *hF, 
						const Vec3 &pos, const double &eRad, const double& ratio){}



protected:
	// CPUによるSPH計算
	void calDensity(const RXREAL *pos, RXREAL *dens, RXREAL h);
	void calNormal(void);
	void calForceExtAndVisc(const RXREAL *pos, const RXREAL *vel, const RXREAL *dens, RXREAL *frc, RXREAL h);

	// Scaling factorの計算
	void calScalingFactor(const RXREAL *ppos, RXREAL *pdens, RXREAL *pscl, RXREAL h, RXREAL dt);

	// Scaling factorの計算
	void calPositionCorrection(const RXREAL *ppos, const RXREAL *pscl, RXREAL *pdp, RXREAL h, RXREAL dt);

	// rest densityの計算
	RXREAL calRestDensity(RXREAL h);

	// 個々の境界パーティクルの体積を計算
	void calBoundaryVolumes(const RXREAL *bpos, RXREAL *bvol, RXREAL mass, uint n, RXREAL h);

	// 時間ステップ幅の修正
	RXREAL calTimeStep(RXREAL &dt, RXREAL eta_avg, const RXREAL *pfrc, const RXREAL *pvel, const RXREAL *pdens);
	// 位置と速度の更新(Leap-Frog)
	void integrate(const RXREAL *pos, const RXREAL *vel, const RXREAL *dens, const RXREAL *acc, 
				   RXREAL *pos_new, RXREAL *vel_new, RXREAL dt);

	// 衝突判定
	int calCollisionPolygon(uint grid_hash, Vec3 &pos0, Vec3 &pos1, Vec3 &vel, RXREAL dt);
	int calCollisionSolid(Vec3 &pos0, Vec3 &pos1, Vec3 &vel, RXREAL dt);
};




//-----------------------------------------------------------------------------
// MARK:rxPBDSPH_GPUクラスの宣言
//  - Miles Macklin and Matthias Muller, "Position Based Fluids", Proc. SIGGRAPH 2013, 2013. 
//  - http://blog.mmacklin.com/publications/
//-----------------------------------------------------------------------------
class rxPBDSPH_GPU : public rxParticleSystemBase
{
private:
	//
	// メンバ変数(GPU変数)
	//
	RXREAL *m_dPos;			//!< パーティクル位置
	RXREAL *m_dVel;			//!< パーティクル速度
	RXREAL *m_dNrm;			//!< パーティクル法線
	RXREAL *m_dFrc;			//!< パーティクルにかかる力
	RXREAL *m_dDens;		//!< パーティクル密度
	RXREAL *m_dPres;		//!< パーティクル圧力

	RXREAL *m_dPosB;		//!< 境界パーティクル
	RXREAL *m_dVolB;		//!< 境界パーティクルの体積

	RXREAL *m_dS;			//!< Scaling factor for CFM
	RXREAL *m_dDp;			//!< 位置修正量
	RXREAL *m_dPredictPos;	//!< 予測位置
	RXREAL *m_dPredictVel;	//!< 予測速度

	RXREAL *m_dSb;			//!< 境界パーティクルのScaling factor

	RXREAL *m_dErr;			//!< 密度変動値	
	RXREAL *m_dErrScan;		//!< 密度変動値のScan結果

	int *m_dAttr;			//!< パーティクル属性

	cudaGraphicsResource *m_pPosResource;	//!< OpenGL(のVBO)-CUDA間のデータ転送を扱うためのハンドル

	// ウェーブレット乱流
	RXREAL *m_dEt;			//!< 速度エネルギースペクトル
	RXREAL *m_dTurb;		//!< ウェーブレット乱流による速度場

	//サブ粒子 GPU data
	RXREAL *m_dSubPos;		//!< サブ粒子の絶対座標 float4(0.0,x,y,z)
	RXREAL *m_dSubChild;	//!< サブ粒子の子1への単位ベクトル
	RXREAL *m_dSubAxis;		//!< サブ粒子の回転軸(単位ベクトル)
	RXREAL *m_dSubEt;		//!< サブ粒子のエネルギースペクトル
	uint   *m_dSubRand;		//!< 乱数テーブル

	uint m_subposVBO;		//!< サブ粒子座標VBO
	cudaGraphicsResource *m_pSubPosResource;	//!< OpenGL(のVBO)-CUDA間のデータ転送を扱うためのハンドル

	// SPS乱流
	RXREAL *m_dUpPos;		//!< 平滑化パーティクル位置
	RXREAL *m_dPosW;		//!< 重み付き平均座標
	RXREAL *m_dCMatrix;		//!< 共分散行列
	RXREAL *m_dEigen;		//!< 共分散行列の特異値
	RXREAL *m_dRMatrix;		//!< 回転行列(共分散行列の特異ベクトル)
	RXREAL *m_dG;			//!< 変形行列

	// 表面メッシュ
	RXREAL *m_dVrts;		//!< 固体メッシュ頂点
	int    *m_dTris;		//!< 固体メッシュ

	// シミュレーションパラメータ
	rxSimParams m_params;	//!< シミュレーションパラメータ(GPUへのデータ渡し用)
	uint3 m_gridSize;		//!< 近傍探索グリッドの各軸の分割数
	uint m_numGridCells;	//!< 近傍探索グリッド分割総数
	
	// 空間分割(GPU)
	rxParticleCell m_dCellData;			//!< 近傍探索グリッド
	rxSubParticleCell m_dSubCellData;	//!< サブ粒子用近傍探索グリッド
	uint m_gridSortBits;				//!< ハッシュ値による基数ソート時の基数桁数

	rxParticleCell m_dCellDataB;		//!< 境界パーティクル用近傍探索グリッド
	uint3 m_gridSizeB;					//!< 境界パーティクル用近傍探索グリッドの各軸の分割数

	//
	// メンバ変数(CPU変数)
	//
	RXREAL *m_hNrm;			//!< パーティクル法線
	RXREAL *m_hFrc;			//!< パーティクルにかかる力
	RXREAL *m_hDens;		//!< パーティクル密度
	RXREAL *m_hPres;		//!< パーティクル圧力

	RXREAL *m_hS;			//!< Scaling factor for CFM
	RXREAL *m_hDp;			//!< 位置修正量
	RXREAL *m_hPredictPos;	//!< 予測位置
	RXREAL *m_hPredictVel;	//!< 予測速度

	RXREAL *m_hSb;			//!< 境界パーティクルのScaling factor

	uint *m_hSurf;			//!< 表面パーティクル

	// 表面生成用(Anisotropic kernel)
	RXREAL *m_hUpPos;		//!< 平滑化パーティクル位置
	RXREAL *m_hPosW;		//!< 重み付き平均座標
	RXREAL *m_hCMatrix;		//!< 共分散行列
	RXREAL *m_hEigen;		//!< 共分散行列の特異値
	RXREAL *m_hRMatrix;		//!< 回転行列(共分散行列の特異ベクトル)
	RXREAL *m_hG;			//!< 変形行列
	RXREAL  m_fEigenMax;	//!< 特異値の最大値(探索半径拡張に用いる)

	// ウェーブレット乱流
	RXREAL *m_hVwt;			//!< パーティクル速度のウェーブレット変換
	RXREAL *m_hEt;			//!< 速度のエネルギースペクトラム
	RXREAL *m_hTurb;		//!< ウェーブレット乱流による速度場
	RXREAL *m_hNoiseTile[3];//!< ウェーブレットノイズタイル
	RXREAL m_fNoiseTileWidth;	//!< ノイズタイルの幅
	uint   m_iNoiseTileN[3];	//!< ノイズタイルの解像度

	// SPS乱流
	RXREAL *m_hSubPos;		//!< サブ粒子の絶対座標
	RXREAL *m_hSubChild;	//!< サブ粒子の子1への単位ベクトル
	RXREAL *m_hSubAxis;		//!< サブ粒子の回転軸(単位ベクトル)
	RXREAL *m_hSubEt;		//!< サブ粒子のエネルギースペクトル

	uint m_uNumSubParticles;//!< サブ粒子数
	uint m_uMaxSubParticles;//!< サブ粒子最大数
	uint m_uMaxSubLevel;	//!< サブ粒子最大レベル

	RXREAL *m_hSubPosPack;
	RXREAL *m_hSubRadPack;
	RXREAL *m_hSubRatPack;
	unsigned int *m_hSubOcc;
	unsigned int *m_hSubOccScan;
	bool m_bSubPacked;		//!< サブパーティクルをPackしたらtrue

	// 表面メッシュ
	vector<RXREAL> m_vVrts;	//!< 固体メッシュ頂点
	int m_iNumVrts;			//!< 固体メッシュ頂点数
	vector<int> m_vTris;	//!< 固体メッシュ
	int m_iNumTris;			//!< 固体メッシュ数

	bool m_bCalNormal;		//!< 法線計算フラグ

	RXREAL m_fEpsilon;				//!< CFMの緩和係数
	RXREAL m_fEta;					//!< 密度変動率
	int m_iMinIterations;			//!< ヤコビ反復最小反復回数
	int m_iMaxIterations;			//!< ヤコビ反復最大反復回数

	bool m_bArtificialPressure;		//!< クラスタリングを防ぐためのArtificial Pressure項を追加するフラグ

	//追加：近傍粒子のGPUからのコピーデータ
	uint* m_hSortedIndex;		//!< ハッシュ値でソートしたパーティクルインデックス
	uint* m_hGridParticleHash;	//!< 各パーティクルのグリッドハッシュ値
	uint* m_hCellStart;			//!< ソートリスト内の各セルのスタートインデックス
	uint* m_hCellEnd;			//!< ソートリスト内の各セルのエンドインデックス
	uint  m_uNumCells;			//!< 総セル数
	vector< vector<rxNeigh> > m_vNeighs;	//!< 近傍パーティクル

protected:
	rxPBDSPH_GPU(){}

public:
	//! コンストラクタ
	rxPBDSPH_GPU(bool use_opengl);

	//! デストラクタ
	~rxPBDSPH_GPU();

	// パーティクル半径
	float GetEffectiveRadius(){ return m_params.EffectiveRadius; }

	// 近傍パーティクル
	uint* GetNeighborList(const int &i, int &n);

	// シミュレーションパラメータ
	rxSimParams GetParams(void){ return m_params; }
	void UpdateParams(void);

	// フラグ切替
	void ToggleNormalCalc(int t = -1){ RX_TOGGLE(m_bCalNormal, t); }			//!< パーティクル法線の計算
	bool IsNormalCalc(void) const { return m_bCalNormal; }
#if MAX_BOX_NUM
	void ToggleSolidFlg(int t = -1);
#endif

public:
	//
	// 仮想関数
	//
	// パーティクルデータ
	virtual RXREAL* GetParticle(void);
	virtual RXREAL* GetParticleDevice(void);
	virtual void UnmapParticle(void);

	// シミュレーションステップ
	virtual bool Update(RXREAL dt, int step = 0);

	// シーンの設定
	virtual void SetPolygonObstacle(const vector<Vec3> &vrts, const vector<Vec3> &nrms, const vector< vector<int> > &tris, Vec3 vel);
	virtual void SetBoxObstacle(Vec3 cen, Vec3 ext, Vec3 ang, Vec3 vel, int flg);
	virtual void SetSphereObstacle(Vec3 cen, double rad, Vec3 vel, int flg);
	virtual void MoveSphereObstacle(int b, Vec3 disp);
	virtual Vec3 GetSphereObstaclePos(int b = -1);

	// ホスト<->VBO間転送
	virtual RXREAL* GetArrayVBO(rxParticleArray type, bool d2h = true, int num = -1);
	virtual void SetArrayVBO(rxParticleArray type, const RXREAL* data, int start, int count);
	virtual void SetArrayVBO(rxParticleArray type, const int* data, int start, int count);
	virtual void SetColorVBO(int type);


	// 陰関数値計算
	virtual double GetImplicit(double x, double y, double z);
	virtual void CalImplicitField(int n[3], Vec3 minp, Vec3 d, RXREAL *hF);
	virtual void CalImplicitFieldDevice(int n[3], Vec3 minp, Vec3 d, RXREAL *dF);

	// SPH情報出力
	virtual void OutputSetting(string fn);

	// 描画関数
	virtual void DrawCell(int i, int j, int k);
	virtual void DrawCells(Vec3 col, Vec3 col2, int sel = 0);
	virtual void DrawObstacles(void);

protected:
	void setObjectToCell(RXREAL *p, RXREAL *v);


public:
	// SPH初期化
	void Initialize(const rxSPHEnviroment &env);
	void Allocate(int max_particles);
	void Finalize(void);

	// 分割セルにパーティクルを格納
	virtual void SetParticlesToCell(void);
	virtual void SetParticlesToCell(RXREAL *prts, int n, RXREAL h);
	virtual void SetPolygonsToCell(void);

	int  GetPolygonsInCell(int gi, int gj, int gk, vector<int> &polys){ return 0; }
	bool IsPolygonsInCell(int gi, int gj, int gk){ return false; }

	void CalMaxDensity(int k);

	// 表面パーティクル検出
	void DetectSurfaceParticles(void);					// 表面パーティクルの検出
	double CalDistToNormalizedMassCenter(const int i);	// 近傍パーティクルの重心までの距離計算
	uint* GetArraySurf(void);							// 表面パーティクル情報の取得

	// 表面パーティクル情報の取得
	int GetSurfaceParticles(const Vec3 pos, RXREAL h, vector<rxSurfaceParticle> &sp);

	// 法線計算
	void CalNormalFromDensity(void);
	void CalNormal(void);

	// 人口圧力項
	bool& GetArtificialPressure(void){ return m_bArtificialPressure; }

	// 境界パーティクルの初期化
	void InitBoundary(void);

	//追加：：近傍粒子取得
	vector<vector<rxNeigh>>& GetNeights(void){ return m_vNeighs; }

	//
	// Anisotropic kernel
	//
	virtual void CalAnisotropicKernel(void);

	//
	// ウェーブレット乱流
	//
	void CalParticleCWT(RXREAL scale);
	void CalWaveletTurbulence(RXREAL scale, RXREAL *dPos, RXREAL dt);

	//
	// SPS乱流
	//
	int	GetMaxSubLevel() const { return m_uMaxSubLevel; }

	// サブパーティクルデータ(デバイスメモリ)
	RXREAL* GetSubParticleDev(void);
	void    UnmapSubParticle(void);
	RXREAL* GetAllSubParticleDev(void);
	void    UnmapAllSubParticle(void);
	RXREAL* GetSubParticlePosDev(void);
	RXREAL* GetSubParticleRadDev(void);
	RXREAL* GetSubParticleRatDev(void);

	int GetNumAllSubParticles(void);
	int GetNumValidSubParticles(void);

	// サブパーティクルデータ(ホストメモリ)
	RXREAL* GetSubParticlePos(void);
	RXREAL* GetSubParticleRad(void);
	RXREAL* GetSubParticleRat(void);
	unsigned int* GetSubParticleOcc(void);
	unsigned int* GetSubParticleOccScan(void);

	RXREAL GetSubRadius(int level){ return m_dSubCellData.fSubRad[level]; }

	//サブ粒子
	void SetEtcri(const double &et_cri){ m_dSubCellData.fEtcri = et_cri;}
	void InitSubParticles(const double &et_cri);
	void AddSubParticles(int start, int count);
	void UpdateSubParticle(RXREAL *dPos, RXREAL *dSubPos, RXREAL scale, RXREAL dt);
	void CalRadiusAndRatio(bool force = false);

	void CalSubImplicitFieldDevice(int n[3], Vec3 minp, Vec3 d, RXREAL *dF, int turb);
	void CalSubImplicitField(int n[3], Vec3 minp, Vec3 d, RXREAL *hF, const double &et_cri);

	void AddSphereField(int n[3], Vec3 minp, Vec3 d, RXREAL *hF, 
						const Vec3 &pos, const double &eRad, const double& ratio);

protected:
	// rest densityの計算
	RXREAL calRestDensity(RXREAL h);

	// 分割セルの初期設定
	void setupCells(rxParticleCell &cell, uint3 &gridsize, double &cell_width, Vec3 vMin, Vec3 vMax, double h);

	// グリッドハッシュの計算
	uint calGridHash(uint x, uint y, uint z);
	uint calGridHash(Vec3 pos);

	// ポリゴンを分割セルに格納
	void setPolysToCell(RXREAL *vrts, int nv, int* tris, int nt);

	// CPU用の近傍粒子探索　追加：：GPU→CPUのデータ転送
	void searchNeighbors(void);
	void getNeighborsInCell(Vec3 pos, RXREAL *p, int gi, int gj, int gk, vector<rxNeigh> &neighs, RXREAL h);
};




#endif	// _SPH_H_

