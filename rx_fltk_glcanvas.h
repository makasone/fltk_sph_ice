/*!
  @file rx_fltk_glcanvas.h
	
  @brief FLTKによるOpenGLウィンドウクラス
 
  @author Makoto Fujisawa 
  @date   2011-09
*/

#ifndef _RX_FLTK_GLCANVAS_H_
#define _RX_FLTK_GLCANVAS_H_

//-----------------------------------------------------------------------------
// インクルードライブラリ
//-----------------------------------------------------------------------------
#pragma comment(lib, "libtet.lib")

//-----------------------------------------------------------------------------
// インクルードファイル
//-----------------------------------------------------------------------------
#include <iostream>

// STL
#include <vector>
#include <string>

// FLTK
#include <FL/Fl.H>
#include <FL/Fl_Gl_Window.H>
#include <FL/Fl_Menu_Bar.H>
#include <FL/Fl_Spinner.H>
#include <FL/Fl_Value_Slider.H>
#include <FL/Fl_Round_Button.H>
#include <FL/Fl_Check_Button.H>
#include <FL/Fl_Text_Editor.H>

#include "rx_sph_commons.h"
#include "rx_sph_config.h"
#include "rx_fltk_widgets.h"

#include "rx_trackball.h"
#include "rx_model.h"
#include "rx_pov.h"

#include "rx_texture.h"

//追加：：
#include "HeatTransfar.h"
#include "Ice_SM.h"
#include "IceStructure.h"
#include "IceObject.h"
#include "tetgen.h"
#include <UtilityScript\mk_Vector2D.h>
#include "QueryCounter.h"
#include "mk_ArrayScript.h"
#include <time.h>

#include <omp.h>
#include <fstream>
#include <algorithm>

#include <boost/bind.hpp>
#include <boost/function.hpp>

#include "test.h"

using namespace std;

//-----------------------------------------------------------------------------
// 定義/定数
//-----------------------------------------------------------------------------
class rxFlWindow;
class rxParticleSystemBase;
struct rxCell;

class rxMCMeshCPU;
class rxMCMeshGPU;
class rxSSMeshCPU;
class rxSSMeshGPU;

//　ここの定義でSPH法のGPU，CPUなどを切り替える
#define RX_USE_GPU

//#define RX_USE_DD
//#define RX_USE_PBD

#if defined(RX_USE_GPU)
	#if defined(RX_USE_PBD)
		#define RXSPH rxPBDSPH_GPU
	#else
		#define RXSPH rxSPH_GPU
	#endif
#else
	#if defined(RX_USE_PBD)
		#define RXSPH rxPBDSPH
	#elif defined(RX_USE_DD)
		#define RXSPH rxDDSPH
	#else
		#define RXSPH rxSPH
	#endif
#endif

//クラスタを作成する際の粒子数，四面体数を指定
//もう使ってない
//#define ICENUM 27
//#define ICENUM 125
//#define ICENUM	729
//#define ICENUM	1331	//11_11_11 or バニーモデル
//#define TETRANUM 10900	//バニーモデルの場合の四面体数
//#define ICENUM	2197	//13_13_13 or バニーモデル
//#define SIDE	13
//#define ICENUM	3463	//バニーモデル
//#define ICENUM	4913	//17_17_17
//#define TETRANUM	26000
//#define ICENUM	6859	//19_19_19
//#define ICENUM	9261	//21_21_21
//#define ICENUM	12167
//#define ICENUM	15625	//25_25_25
//#define ICENUM	19683	//27_27_27
//#define TETRANUM	110000
//#define ICENUM	24389	//29_29_29

//(1000, 1000, 1000);
//(3000, 3000, 12000);
//(5000, 5000, 26000);			//粒子数4913個の場合のパラメータ
//(7000, 7000, 37000);
//(13000, 13000, 50000);
//(10000, 10000, 1);
//(16000, 16000, 1);

//(16000, 16000, 90000);
//(20000, 20000, 110000);
//(24400, 24400, 137000);		//動かなかった．

//表面のみの場合
//#define ICENUM	6			//1_1_1
//#define ICENUM	27
//#define ICENUM	54			//3_3_3
//#define ICENUM	1014
//#define ICENUM	2646	//21_21_21
//#define ICENUM	5046	//29_29_29

#define HIGHNUM 8

// 描画フラグ
enum
{
	RXD_PARTICLE		= 0x0001,	//!< パーティクル
	RXD_VELOCITY		= 0x0002,	//!< 速度場
	RXD_NORMAL			= 0x0004,	//!< 法線
	RXD_BBOX			= 0x0008,	//!< AABB(シミュレーション空間)
	RXD_CELLS			= 0x0010,	//!< 近傍探索用分割セル
	RXD_MESH			= 0x0020,	//!< メッシュ
	RXD_SOLID			= 0x0040,	//!< 固体
	RXD_REFRAC			= 0x0080,	//!< 屈折描画

	RXD_TURB_VELOC		= 0x0100,	//!< 乱流速度場
	RXD_FOAM			= 0x0200,	//!< 泡

	RXD_PARAMS			= 0x0400,	//!< パラメータ画面描画
	RXD_ANISOTROPICS	= 0x0800,	//!< 異方性カーネル
	RXD_UPDATED_PRTS	= 0x1000,	//!< 
	RXD_AXIS			= 0x2000,   //!< 軸

	RXD_BPARTICLE		= 0x4000,	//!< 境界パーティクル
	RXD_DEBUG_VECTOR	= 0x8000,	//!< デバッグ用のベクトル場
};
const string RX_DRAW_STR[] = {
	"Particle",					"p", 
	"Velocity",					"v", 
	"Normal",					"n",
	"AABB (Simulation Space)", 	"b",
	"Cells", 					"d",
	"Mesh", 					"m",
	"Solid", 					"o",
	"Refrac", 					"r",
	"Turbulence Velocity", 		"t",
	"Foam", 					"F",
	"Params", 					"P",
	"Anisotoropics", 			"A",
	"Updated Particles", 		"u",
	"Axis",				 		"",
	"Boundary Particle",		"B",
	"Debug Vectors",				"",
	"-1"
};

// パーティクル描画方法
enum
{
	RXP_POINTSPRITE = 0, 
	RXP_POINT, 
	RXP_POINT_UPSAMPLE, 
	RXP_POINT_NONE, 

	IDP_END, 
};
const string RX_PARTICLE_DRAW[] = {
	"Point Sprite", "^1", 
	"GL_POINT", 	"^2", 
	"Upsampling", 	"^3", 
	"None", 		"^4", 
	"-1"
};


// 三角形メッシュ生成法
enum
{
	RXM_MC_CPU = 0, 
	RXM_MC_GPU, 
	RXM_SSM_CPU, 
	RXM_SSM_GPU, 
};
const string RX_TRIANGULATION_METHOD[] = {
	"Marching Cubes (CPU)",    "", 
	"Marching Cubes (GPU)",    "", 
	"Screen Space Mesh (CPU)", "", 
	"Screen Space Mesh (GPU)", "", 
	"-1"
};

// 固体描画
enum
{
	RXS_VERTEX		= 0x0001, 
	RXS_EDGE		= 0x0002, 
	RXS_FACE		= 0x0004, 
	RXS_NORMAL		= 0x0008, 
	RXS_MOVE		= 0x0010,

	RXS_END
};
const string RXS_STR[] = {
	"Vertex", "", 
	"Edge",   "", 
	"Face",   "", 
	"Normal", "", 
	"Move",   "", 
	"-1"
};


//! 描画領域サイズ候補
const string RX_CANVAS_SIZE_STR[] = {
	"1920x1080",	"",
	"1280x720",		"", 
	"1024x768",		"", 
	"800x800",		"", 
	"800x600",		"", 
	"640x480",		"", 
	"-1", 
};


//! SPH設定
enum SettingMode
{
	ID_SPH_MESH = 0,	// メッシュ生成
	ID_SPH_INPUT,		// パーティクルデータ出力
	ID_SPH_OUTPUT,		// パーティクルデータ入力
	ID_SPH_MESH_OUTPUT, // メッシュ出力
	ID_SPH_INLET,		// 流入境界
	ID_SPH_VC, 
	ID_SPH_PS_TURB, 
	ID_SPH_SPS_TURB, 

	ID_HEAT,			//追加　熱処理
	ID_SM,				//追加　SM法
	ID_ICE,				//追加　氷構造

	ID_SPH_ANISOTROPIC, // 異方性カーネル

	ID_SPH_END, 
};

//-----------------------------------------------------------------------------
//! rxFlGLWindowクラス - fltkによるOpenGLウィンドウ
//-----------------------------------------------------------------------------
class rxFlGLWindow : public Fl_Gl_Window
{
protected:
	int m_iWinW;					//!< 描画ウィンドウの幅
	int m_iWinH;					//!< 描画ウィンドウの高さ
	int m_iMouseButton;				//!< マウスボタンの状態
	int m_iKeyMod;					//!< 修飾キーの状態
	rxTrackball m_tbView;			//!< トラックボール

	double m_fBGColor[3];			//!< 背景色
	bool m_bAnimation;				//!< アニメーションON/OFF
	bool m_bFullscreen;				//!< フルスクリーンON/OFF

	rxFlWindow *m_pParent;			//!< 親クラス

	vector<rxPolygons> m_vPolys;	//!< ポリゴンオブジェクト

	// FTGL
	unsigned long m_ulFontSize;		//!< フォントサイズ

	//
	// 粒子法関連変数
	//
	rxParticleSystemBase *m_pPS;	//!< SPH
	double m_fDt;					//!< タイムステップ幅
	double m_fGravity;				//!< 重力加速度

	int m_iCurrentStep;				//!< 現在のステップ数
	bool m_bPause;					//!< シミュレーションのポーズフラグ

public:
	//
	// 追加　熱処理関連変数
	//
	vector<float> m_fIntrps;		//SPH法とSM法の線形補間のパラメータ配列　各パーティクルごとに適用してやる

	Vec2 m_ht_vStartPoint;			//矩形内の粒子温度を上げるための始点
	Vec2 m_ht_vEndPoint;			//終点

	bool m_ht_bRectFlag;			//矩形内温度変化機能を利用するかのフラグ

	vector<int> m_ht_vSelectedVertices;

	float m_fMakeParticleTemp;

	//
	//追加　相変化オブジェクト
	//
	IceObject* m_iceObj;

	//
	// 追加　SM法関連変数
	//
	vector<Ice_SM*> m_sm_connects;	//関係情報クラスタ
	vector<Ice_SM*> m_sm_calcs;		//計算処理クラスタ

	int m_iIcePrtNum;				//初期の固体粒子数
	int m_iIceTtrNum;				//初期の固体粒子の四面体数
	int m_iClusteresNum;			//クラスタの数　使用中のクラスタで最大の値
	int m_iTetraNum;				//未使用　四面体の数
	int m_iTetraNumNum;				//デバッグ用

	int m_iMaxClusterNum;			//最大クラスタ数，最大粒子数
	int m_iIceItr;					//固体運動計算の反復回数

	//
	// 追加　氷構造関連変数
	//
	IceStructure *m_ice;

	int meltPIndx;
	int debugIndx;

	float *m_fIceFlag;						//氷のフラグ　デバイスメモリ

	bool *testFlag;							//実行速度を確かめるための実験

	vector<vector<int>> m_vviTetraList;		//四面体の組み合わせリスト

	int m_iLayer;							//計算処理クラスタ作成のために取得する近傍クラスタの，レイヤー数
	int m_iShowClusterIndx;					//GUIで表示するクラスタ番号　０～クラスタ数　値はctr+shift+qで変化させる
	int m_iShowHighClstrIndx;
	int m_iShowTetraIndx;					//GUIで表示する四面体番号　０～四面体数　値はctr+shift+aで変化させる
	
	//objファイル
	rxPolygons m_poly;

	//
	//追加　処理切り替え
	//
	int m_bMode;							//現在の処理モード
	enum
	{
		MODE_SPH = 0,
		MODE_HEAT,
		MODE_SM,
		MODE_ICE,
		MODE_NUM,
	};

	//追加　レンダリングパラメータ
	double m_etaRatio;
	double m_fresnelBias;
	double m_fresnelPower;
	double m_fresnelScale;

	//落下開始フラグ
	bool m_bFall;

	//テストクラス
	//mk_CGAL test;

protected:
	// シーン
	//string m_strCurrentScene;		//!< 現在のシーンの名前
	//vector<string> m_vSceneFiles;	//!< シーンファイルリスト
	//int m_iSceneFileNum;			//!< シーンファイルの数
	rxSPHConfig m_Scene;
	int m_iSimuSetting;				//!< ミュレーション設定保存用

	// パーティクル情報出力
	string m_strSphOutputName0 ;
	string m_strSphOutputHeader;

	// パーティクル情報入力
	string m_strSphInputName0 ;
	string m_strSphInputHeader;

	// 固体移動フラグ
	bool m_bSolidMove;

	// 固体の動き
	Vec3 m_vMovePos[2];
	double m_fMoveMaxVel;
	bool m_bMoveSolid;
	int m_iMoveStart;


	//
	// メッシュ
	//
	uint m_iNumVrts, m_iNumTris;	//!< 生成されたメッシュの頂点数とメッシュ数
	int m_iVertexStore;				//!< サンプリングボリューム数に対する予想される頂点数(nx*ny*store)

	int m_iMeshMaxN;				//!< メッシュ化グリッド数(境界がもっとも長い軸方向の分割数)
	int m_iMeshN[3];				//!< メッシュ化グリッド数

	Vec3 m_vMeshBoundaryExt;		//!< メッシュ境界ボックスの各辺の長さの1/2
	Vec3 m_vMeshBoundaryCen;		//!< メッシュ境界ボックスの中心座標

	rxPolygons m_Poly;				//!< メッシュ
	//vector<rxPolygons*> m_vSolidPoly;//!< 固体メッシュ
	rxMaterialOBJ m_matPoly;

	GLuint m_uVrtVBO;				//!< メッシュ頂点(VBO)
	GLuint m_uTriVBO;				//!< メッシュポリゴン(VBO)
	GLuint m_uNrmVBO;				//!< メッシュ頂点法線(VBO)

	//追加：固体用
	uint m_iNumVrts_solid, m_iNumTris_solid;	//!< 生成されたメッシュの頂点数とメッシュ数
	GLuint m_uVrtVBO_solid;			//!< メッシュ頂点(VBO)
	GLuint m_uTriVBO_solid;			//!< メッシュポリゴン(VBO)
	GLuint m_uNrmVBO_solid;			//!< メッシュ頂点法線(VBO)

	int m_iDimVBO;

	// メッシュ出力
	int m_iSaveMeshSpacing;

	// 背景画像
	bool m_bUseCubeMap;				//!< キューブマップ使用フラグ
	rxCubeMapData m_CubeMap;		//!< キューブマップ

	// メッシュ生成
	rxMCMeshCPU *m_pMCMeshCPU;
	rxMCMeshGPU *m_pMCMeshGPU;
	rxSSMeshCPU *m_pSSMCPU;			//!< Screen Space Mesh
	rxSSMeshGPU *m_pSSMGPU;			//!< Screen Space Mesh by CUDA
	int m_iDepthFiltering;			//!< 平滑化(デプスマップ0x01, 輪郭0x02)
	int m_iSSDebugOutput;

	double m_fSpacing;				//!< デプスマップのサンプリング間隔
	double m_fPrtRad;				//!< パーティクルの半径
	double m_fZmax;					//!< 輪郭となるデプス差の閾値
	int m_iNfilter;					//!< デプス値平滑化のフィルタサイズ
	int m_iNiters;					//!< 輪郭平滑化の反復回数

	double m_fProjectionMatrix[16];	//!< 透視投影行列
	double m_fModelviewMatrix[16];	//!< モデルビュー行列

	int m_iPickedParticle;			//!< マウスピックされたパーティクル


public:
	// 描画フラグ
	int m_iDraw;					//!< 描画フラグ
	int m_iDrawPS;					//!< パーティクル描画方法
	int m_iColorType;				//!< パーティクル描画時の色
	int m_iTriangulationMethod;		//!< 三角形メッシュ生成法
	int m_iSolidDraw;				//!< 固体描画

	// シミュレーション設定
	bitset<32> m_bsSimuSetting;		//!< ミュレーション設定フラグ
	double m_fVScale;				//!< ベクトル場描画時のスケール
	double m_fMeshThr;				//!< 陰関数メッシュ化時の閾値

	// シーンリスト
	vector<string> m_vSceneTitles;	//!< シーンファイルリスト
	int m_iCurrentSceneIdx;			//!< 現在のシーンファイル

	// 画像出力
	int m_iSaveImageSpacing;		//!< 画像保存間隔(=-1なら保存しない)

public:
	//! コンストラクタ
	rxFlGLWindow(int x, int y, int w, int h, const char* l, void *parent);

	//! デストラクタ
	~rxFlGLWindow();

public:
	// OpenGL初期化
	void InitGL(void);

	// OpenGL描画
	void Projection(void);
	vector<string> SetDrawString(void);
	void ReDisplay(void);

	// GUIコールバック
	void Display(void);
	void Resize(int w, int h);
	void Mouse(int button, int state, int x, int y);
	void Motion(int x, int y);
	void PassiveMotion(int x, int y);
	void Idle(void);
	void Timer(void);
	void Keyboard(int key, int x, int y);
	void SpecialKey(int key, int x, int y);

	// マウスピック用
	static void Projection_s(void* x);
	static void DisplayForPick_s(void* x);
	void DisplayForPick(void);
	void PickParticle(int x, int y);

	// 視点
	void InitView(void);

	// アニメーション
	static void OnTimer_s(void* x);
	static void OnIdle_s(void* x);
	
	// アニメーション切り替え
	bool SwitchIdle(int on);

	// フルスクリーン切り替え
	void SwitchFullScreen(int win = 1);
	int  IsFullScreen(void);


	// ファイル入出力
	void OpenFile(const string &fn);
	void SaveFile(const string &fn);
	void SaveDisplay(const string &fn);
	void SaveDisplay(const int &stp);
	
	void SaveMesh(const string fn, rxPolygons &polys);
	void SaveMesh(const int &stp, rxPolygons &polys);

	// FTGLフォント設定
	int SetupFonts(const char* file);


public:
	// メッシュ生成
	bool CalMeshSPH(int nmax, double thr = 300.0);
	bool ResetMesh(void);
	void SetMeshCreation(void);
	RXREAL GetImplicitSPH(double x, double y, double z);

protected:
	bool calMeshSPH_CPU(int nmax, double thr = 300.0);
	bool calMeshSPH_GPU(int nmax, double thr = 300.0);
	bool calMeshSPH_SSM(int nmax, double thr = 300.0);
	bool calMeshSPH_SSM_GPU(int nmax, double thr = 300.0);

	// SPH
	void InitSPH(rxSPHConfig &sph_scene);
	void InitSPH(void){ InitSPH(m_Scene); }
	void AddSphere(void);

	void StepPS(double dt);
	void ComputeFPS(void);

	void DivideString(const string &org, vector<string> &div);

	void DrawParticleVector(RXREAL *prts, RXREAL *vels, int n, int d, double *c0, double *c1, double len = 0.1);
	void DrawParticlePoints(unsigned int vbo, int n, unsigned int color_vbo = 0, RXREAL *data = 0);
	void DrawSubParticles(void);

	void CreateVBO(GLuint* vbo, unsigned int size);
	void DeleteVBO(GLuint* vbo);

	void DrawLiquidSurface(void);

	//追加：氷用
	void DrawSolidSurface(void);

	//追加：高速化パス
	void DrawFastPath(RXREAL *prts, RXREAL *vels, int n, int d, double *c0, double *c1, double len = 0.1);

	void SetParticleColorType(int type, int change = 0);

	void RenderSphScene(void);

	//追加：：氷
	void InitIceObj(void);
	void StepTimeEvent(void);

	void CountTetraHedra(int tIndx, vector<int>& pList);
	void MakeTetraInfo(int tIndx, int* PtoTNum);
	void MakeTetraInfo(int tIndx, vector<int> pList);

	//追加：：粒子ベース：：クラスタ
	void MakeCluster(int pIndx);
	void MakeClusterFromNeight();

	void StepIceObj();
	void StepIceStructure();
	void StepHeatTransfer();

	void StepClusterCPU(double dt);
	void StepClusterGPU(double dt);

	//追加：：粒子ベース：：固体情報（クラスタ情報）
	void MakeClusterInfo(int cIndx);

	void SearchReconstructCluster_Freeze(const vector<int>& pList, vector<int>& cList, vector<int>& lList);

	void SearchReconstructTetra_Freeze(const vector<int>& pList, vector<int>& tList, vector<int>& lList);

	void UpdateInfo_Delete_TandP(const vector<int>& tList, const vector<int>& deleteList);

	void SetClusterInfo(const vector<int>& pList, const vector<int>& cList, const vector<int>& lList);
	void SetTetraInfo(const vector<int>& pList, const vector<int>& cList, const vector<int>& lList);
	void SetFreezeTetraInfo(vector<int>& pList);
	void SetFreezeClusterInfo(const vector<int>& pList);

	void CheckDeleteCluster(void);
	void CheckDeleteTetra(vector<int>& tList, vector<int>& lList);
	void CheckSameTetra(int tIndx, const vector<int>& searchList, vector<int>& deleteList);
	void CheckIncludeTetra(int tIndx, const vector<int>& searchList, vector<int>& deleteList);

	void RemoveReconstructTetra(vector<int>& tList, vector<int>& lList, vector<int>& deleteTList);

	//追加：：四面体ベース：：SM法
	void InitSM(rxSPHConfig &sph_scene);
	void InitSM_Layer(rxSPHConfig &sph_scene);
	void StepSM(double dt);

	void AddVertexToCluster(int pIndx, int cIndx);

	//追加：：四面体ベース：：氷構造
	void InitICE(rxSPHConfig &sph_scene);
	void InitICE_Layer(rxSPHConfig &sph_scene);

	void StepICE_Melt(void);

	void StepICE_Freeze(void);
	bool StepICE_Freeze_MakeCluster(vector<int>& pList);
	bool StepICE_Freeze_AddParticle(vector<int>& pList);
	bool StepICE_Freeze_PhaseChangeParticles(vector<int>& pList);

	void StepICE_CalcParametor(double dt);
	void StepICE_Interpolation(double dt);
	void StepICE_ClusterConstruct(void);

	void ReConstructCalcCluster(void);
	void ReConstructCalcCluster_Melt(vector<int>& pIndxList);
	void ReConstructCalcCluster_Freeze(vector<int>& pIndxList);

	vector<int> GetNeighborDistanceIceParticles(int pIndx, int layer, int num);

	vector<int> CheckSameCluster_Connect(int cIndx);
	void CheckZeroCluster_Connect(int cIndx);

	void DeleteIceInfoParticle(int pIndx);
	void DeleteIceInfoCluster(int cIndx);

	//追加：：その他処理
	void Collision(Vec3 &p, Vec3 &np, Vec3 &v, int obj);	//衝突判定クラス
	void ClearPick(void);
	void ClearRect(void);
	void StepParticleColor(void);
	void UpdateInfo();

	//追加：：四面体作成のための処理
	void MakeFreezeTetrahedra(vector<int>& pList, vector<int>& tList);
	void MakeFreezeTetrahedra_OnlyFreezeParticle(const vector<int>& pList, vector<int>& tList);
	void AddFreezeTetrahedra(const vector<int>& pList, vector<int>& tList);

	void DumpParticleData();

private:
	//! 描画コールバック関数のオーバーライド
	void draw(void)
	{
		if(!context_valid()) InitGL();
		if(!valid()) Resize(w(), h());
		Display();    // OpenGL描画
	}

	//! リサイズコールバック関数のオーバーライド
	void resize(int x_, int y_, int w_, int h_)
	{
		Fl_Gl_Window::resize(x_, y_, w_, h_);
		//Resize(w_, h_);
	}


public:
	// イベントハンドラ
	int handle(int e);	// handle関数のオーバーライド

	// メニューイベントハンドラ
	static void OnMenuDraw_s(Fl_Widget *widget, void* x);
	inline void OnMenuDraw(double val, string label);

	static void OnMenuSimulation_s(Fl_Widget *widget, void* x);
	inline void OnMenuSimulation(double val, string label);

	static void OnMenuParticle_s(Fl_Widget *widget, void* x);
	inline void OnMenuParticle(double val, string label);
	inline void OnMenuParticleColor(double val, string label);
	inline void OnMenuParticleDraw(double val, string label);

	static void OnMenuSolid_s(Fl_Widget *widget, void* x);
	inline void OnMenuSolid(double val, string label);

	static void OnMenuTriangulation_s(Fl_Widget *widget, void* x);
	inline void OnMenuTriangulation(double val, string label);

	static void OnMenuScene_s(Fl_Widget *widget, void* x);
	inline void OnMenuScene(double val, string label);
};




#endif // #ifdef _RX_FLTK_GLCANVAS_H_
