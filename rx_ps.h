/*!
  @file rx_ps.h
	
  @brief パーティクルを扱うシミュレーションの基底クラス
 
  @author Makoto Fujisawa
  @date 2011-06
*/
// FILE --rx_ps.h--

#ifndef _RX_PS_H_
#define _RX_PS_H_


//-----------------------------------------------------------------------------
// MARK:インクルードファイル
//-----------------------------------------------------------------------------
#include "rx_sph_commons.h"

#include "rx_cu_common.cuh"

#include "rx_sph_solid.h"

//#include <helper_functions.h>



//-----------------------------------------------------------------------------
// 定義
//-----------------------------------------------------------------------------
#ifndef DIM
	#define DIM 4
#endif

const int RX_MAX_STEPS = 100000;


//-----------------------------------------------------------------------------
// パーティクル流入ライン
//-----------------------------------------------------------------------------
struct rxInletLine
{
	Vec3 pos1, pos2;	//!< ラインの端点
	Vec3 vel;			//!< 追加するパーティクルの速度
	Vec3 up;			//!< パーティクル堆積方向
	int accum;			//!< パーティクル堆積数
	int span;			//!< 時間的なスパン
	double spacing;		//!< 空間的なスパン
};


//-----------------------------------------------------------------------------
// パーティクルを扱うシミュレーションの基底クラス
//-----------------------------------------------------------------------------
class rxParticleSystemBase
{
public:

	enum rxParticleConfig
	{
		RX_CONFIG_RANDOM,
		RX_CONFIG_GRID,
		RX_CONFIG_BOX,
		RX_CONFIG_NONE, 
		_NUM_CONFIGS
	};

	enum rxParticleArray
	{
		RX_POSITION = 0,
		RX_VELOCITY,
		RX_NORMAL, 
		RX_FORCE, 
		RX_DENSITY, 
		RX_PRESSURE, 

		RX_SURFACE, 
		RX_ATTRIBUTE, 

		RX_PREDICT_POS, 
		RX_PREDICT_VEL, 
		RX_SCALING_FACTOR, 
		RX_CORRECTION, 

		RX_BOUNDARY_PARTICLE, 
		RX_BOUNDARY_PARTICLE_VOL, 

		RX_TEST, 
		RX_DEBUG_VECTOR, 

		RX_TURB_VELOCITY, 
		RX_ENERGY_SPECTRUM, 
		RX_VORTICITY, 

		RX_UPDATED_POSITION, 
		RX_EIGEN_VALUE, 
		RX_ROTATION_MATRIX, 
		RX_TRANSFORMATION, 
		RX_SUBPOSITION,

		RX_CONSTANT, 
		RX_RAMP, 
		RX_NONE, 

		//追加
		RX_TEMP,			//熱処理のために新しく追加　温度
		RX_SHPMTCHNG,		//クラスタの計算のために新しく追加　クラスタ計算
		RX_ICE_CONNECT,		//接続クラスタの判別のために追加　氷塊
		RX_ICE_CALC,		//計算クラスタの判別のために追加　氷塊
		RX_EDGE,			//表面氷粒子取得のために追加　氷塊
		RX_ICE_FAST_PATH,	//高速計算用のパスのために追加
		RX_ICE_HIGH_CLUSTER,//

		RX_PSDATA_END, 
	};

protected:
	bool m_bInitialized;
	bool m_bUseOpenGL;

	uint m_uNumParticles;	//!< 現在のパーティクル数
	uint m_uMaxParticles;	//!< 最大パーティクル数

	uint m_uNumBParticles;	//!< 境界パーティクルの数

	uint m_uNumArdGrid;

	uint m_solverIterations;

	RXREAL m_fParticleRadius;
	Vec3   m_v3Gravity;
	RXREAL m_fDamping;
	RXREAL m_fRestitution;

	RXREAL *m_hPos;		//!< パーティクル位置
	RXREAL *m_hVel;		//!< パーティクル速度

	int *m_hAttr;		//!< パーティクル属性

	RXREAL *m_hPosB;	//!< 境界パーティクル
	RXREAL *m_hVolB;	//!< 境界パーティクルの体積

	RXREAL *m_hSb;		//!< 境界パーティクルのScaling factor

	uint m_posVBO;		//!< パーティクル座標VBO
	uint m_colorVBO;	//!< カラーVBO

	RXREAL *m_hTmp;		//!< 一時的な値の格納場所
	RXREAL m_fTmpMax;
	RXREAL *m_hDebugVec;	//!< テストデータの格納(ベクトルデータ)


	Vec3 m_v3EnvMin;	//!< 環境のAABB最小座標
	Vec3 m_v3EnvMax;	//!< 環境のAABB最大座標
	
	int m_iColorType;

	RXREAL m_fTime;

	bool m_bUseVorticity;	//!< 渦度強制使用フラグ
	bool m_bUseWaveletTurb;	//!< ウェーブレット乱流使用フラグ
	bool m_bGridVelocity;	//!< グリッドへの速度場投影フラグ
	bool m_bCalNormal;		//!< 法線計算フラグ
	bool m_bUpsampling;		//!< パーティクルの再サンプリングフラグ
	bool m_bSubParticle;	//!< サブパーティクル使用フラグ

	int m_iDeleteRegion;	//!< 削除領域使用フラグ
	vector<Vec3> m_vDeleteRegion;	//!< パーティクル削除領域

	vector<rxInletLine> m_vInletLines;	//!< 流入ライン
	int m_iInletStart;		//!< パーティクル追加開始インデックス

public:	
	vector<RXREAL> m_vFuncB;

protected:
	rxParticleSystemBase(){}

public:
	//! コンストラクタ
	rxParticleSystemBase(bool bUseOpenGL) : 
		m_bInitialized(false),
		m_bUseOpenGL(bUseOpenGL), 
		m_hPos(0),
		m_hVel(0), 
		m_hAttr(0), 
		m_hTmp(0), 
		m_hDebugVec(0), 
		m_hPosB(0), 
		m_hVolB(0)
	{
		m_v3Gravity = Vec3(0.0, -9.82, 0.0);
		m_fDamping = 0.0;
		m_fRestitution = 0.0;
		m_fParticleRadius = 0.1;
		m_fTime = 0.0;
		m_bCalAnisotropic = false;
		m_iDeleteRegion = 0;
		m_iInletStart = -1;
		m_fTmpMax = 1.0;
	}

	//! デストラクタ
	virtual ~rxParticleSystemBase(){}


	// シミュレーション空間
	Vec3 GetMax(void) const { return m_v3EnvMax; }
	Vec3 GetMin(void) const { return m_v3EnvMin; }
	Vec3 GetDim(void) const { return m_v3EnvMax-m_v3EnvMin; }
	Vec3 GetCen(void) const { return 0.5*(m_v3EnvMax+m_v3EnvMin); }

	// パーティクル数
	int	GetNumParticles() const { return m_uNumParticles; }
	int	GetMaxParticles() const { return m_uMaxParticles; }
	int GetNumBoundaryParticles() const { return m_uNumBParticles; }

	// シミュレーション反復回数
	void SetIterations(int i) { m_solverIterations = i; }
		
	// パーティクル半径
	float GetParticleRadius(){ return m_fParticleRadius; }

	// シミュレーション設定
	void SetDamping(RXREAL x){ m_fDamping = x; }	//!< 固体境界での反発
	void SetGravity(RXREAL x){ m_v3Gravity = Vec3(0.0, x, 0.0); }	//!< 重力

	// パーティクルVBO
	unsigned int GetCurrentReadBuffer() const { return m_posVBO; }
	unsigned int GetColorBuffer()	   const { return m_colorVBO; }

	// 描画用カラー設定
	void SetColorType(int type){ m_iColorType = type; }
	int  GetColorType(void) const { return m_iColorType; }

	// フラグ切替
	void ToggleUseVorticity(int t = -1){ RX_TOGGLE(m_bUseVorticity, t); }	//!< 渦度強制
	bool IsUseVorticity(void) const { return m_bUseVorticity; }
	void ToggleWaveletTurb(int t = -1){ RX_TOGGLE(m_bUseWaveletTurb, t); }	//!< パーティクルDWTによる乱流
	bool IsWaveletTurb(void) const { return m_bUseWaveletTurb; }
	void ToggleGridVelocity(int t = -1){ RX_TOGGLE(m_bGridVelocity, t); }	//!< グリッドへの速度の投影
	bool IsGridVelocity(void) const { return m_bGridVelocity; }
	void ToggleNormalCalc(int t = -1){ RX_TOGGLE(m_bCalNormal, t); }		//!< パーティクル法線の計算
	bool IsNormalCalc(void) const { return m_bCalNormal; }
	void ToggleUpsampling(int t = -1){ RX_TOGGLE(m_bUpsampling, t); }		//!< 表面パーティクルのアップサンプリング
	bool IsUpsampling(void) const { return m_bUpsampling; }
	void ToggleSubParticle(int t = -1){ RX_TOGGLE(m_bSubParticle, t); }		//!< 乱流サブパーティクル
	bool IsSubParticle(void) const { return m_bSubParticle; }

public:
	// 純粋仮想関数
	virtual bool Update(RXREAL dt, int step = 0) = 0;

	virtual RXREAL* GetArrayVBO(rxParticleArray type, bool d2h = true, int num = -1) = 0;
	virtual void SetArrayVBO(rxParticleArray type, const RXREAL* data, int start, int count) = 0;
	virtual void SetArrayVBO(rxParticleArray type, const int* data, int start, int count) = 0;
	virtual void SetColorVBO(int type) = 0;

	virtual RXREAL* GetParticle(void) = 0;
	virtual RXREAL* GetParticleDevice(void) = 0;

public:
	// 仮想関数
	virtual void UnmapParticle(void){}

	virtual void SetPolygonObstacle(const vector<Vec3> &vrts, const vector<Vec3> &nrms, const vector< vector<int> > &tris, Vec3 vel){}
	virtual void SetBoxObstacle(Vec3 cen, Vec3 ext, Vec3 ang, Vec3 vel, int flg){}
	virtual void SetSphereObstacle(Vec3 cen, double rad, Vec3 vel, int flg){}
	virtual void MoveSphereObstacle(int b, Vec3 disp){}
	virtual Vec3 GetSphereObstaclePos(int b = -1){ return Vec3(0.0); }

	virtual void SetParticlesToCell(void) = 0;
	virtual void SetParticlesToCell(RXREAL *prts, int n, RXREAL h) = 0;

	virtual void SetPolygonsToCell(void){}

	// 陰関数値計算
	virtual double GetImplicit(double x, double y, double z){ return 0.0; }
	virtual void CalImplicitField(int n[3], Vec3 minp, Vec3 d, RXREAL *hF){}
	virtual void CalImplicitFieldDevice(int n[3], Vec3 minp, Vec3 d, RXREAL *dF){}
	//追加
	virtual void CalImplicitFieldDeviceSolid(int n[3], Vec3 minp, Vec3 d, RXREAL *dF, float *fIceCheck){}

	// 描画関数
	virtual void DrawCell(int i, int j, int k){}
	virtual void DrawCells(Vec3 col, Vec3 col2, int sel = 0){}

	virtual void DrawObstacles(void){}

	// シミュレーション設定の出力
	virtual void OutputSetting(string fn){}

	// Anisotropic Kernel
	virtual void CalAnisotropicKernel(void){}
	bool m_bCalAnisotropic;

public:
	void Reset(rxParticleConfig config);
	bool Set(const vector<Vec3> &ppos, const vector<Vec3> &pvel);

	void AddSphere(int start, RXREAL *pos, RXREAL *vel, int r, RXREAL spacing, int attr = 0);
	void AddBox(int start, Vec3 cen, Vec3 dim, Vec3 vel, RXREAL spacing, int attr = 0);

	//追加：粒子を表面に配置
	void AddBoxSurface(int start, Vec3 cen, Vec3 dim, Vec3 vel, RXREAL spacing, int attr = 0);

	//追加：粒子配置
	bool SetParticle(uint& indx, uint& count, RXREAL dx[3], Vec3& cen, Vec3& vel, RXREAL jitter, int attr);

	int  AddLine(rxInletLine line);

	RXREAL SetColorVBOFromArray(RXREAL *hVal, int d, bool use_max = true, RXREAL vmax = 1.0);
	void SetColorVBO(void){ SetColorVBO(m_iColorType); }

	int OutputParticles(string fn);
	int InputParticles(string fn);

protected:
	int  addParticles(int &start, rxInletLine line, int attr = 0);

	uint createVBO(uint size)
	{
		GLuint vbo;
		glGenBuffers(1, &vbo);
		glBindBuffer(GL_ARRAY_BUFFER, vbo);
		glBufferData(GL_ARRAY_BUFFER, size, 0, GL_DYNAMIC_DRAW);
		glBindBuffer(GL_ARRAY_BUFFER, 0);
		return vbo;
	}

};



#endif	// _PS_H_

