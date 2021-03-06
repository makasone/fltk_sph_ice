/*!
  @file rx_fltk_glcanvas.cpp
	
  @brief FLTKによるOpenGLウィンドウクラス
 
  @author Makoto Fujisawa 
  @date   2011-09
*/
// FILE --rx_fltk_glcanvas.cpp--

#pragma warning (disable: 4996)
#pragma warning (disable: 4819)


//-----------------------------------------------------------------------------
// インクルードファイル
//-----------------------------------------------------------------------------
#include "rx_fltk_glcanvas.h"
#include "rx_fltk_window.h"

// テクスチャ・画像
#include "rx_texture.h"
//#include "rx_jpeg.h"
#include "rx_png.h"
#include "rx_bitmap.h"

// シミュレーション
#include "rx_sph.h"

// 設定ファイル
#include "rx_atom_ini.h"

// OpenGL描画関連
//#include "rx_gltexture.h"	// GLテクスチャ & PNGによる画像保存
#include "rx_trackball.h"	// 視点変更用トラックボールクラス
#include "rx_glgui.h"		// GLを使ったGUI
#include "rx_gldraw.h"		// GL描画関数群
#include "rx_shaders.h"		// GLSL関数
#include "rx_shadow.h"		// シャドウマップによる影付け

#include "rx_pick.h"

// CUDA
#include "rx_cu_funcs.cuh"

//#include <helper_cuda.h>
#include <helper_timer.h>

// メッシュ化
#include "rx_model.h"
#include "rx_mesh.h"		// メッシュ構造体，クラス定義
#include "rx_mc.h"			// Marching Cubes
#include "rx_ssm.h"			// Screen Space Mesh

// OpenCV
//#include <opencv2/opencv.hpp>

// FTGL
#include <FTGL/ftgl.h>


//-----------------------------------------------------------------------------
// 定数・変数
//-----------------------------------------------------------------------------
const uint NUM_PARTICLES = 50000;

const GLfloat RX_LIGHT0_POS[4] = {  0.0f, 1.0f, 0.0f, 1.0f };
const GLfloat RX_LIGHT1_POS[4] = { -1.0f, -10.0f, -1.0f, 0.0f };

const GLfloat RX_LIGHT_AMBI[4] = { 0.2f, 0.2f, 0.2f, 1.0f };
const GLfloat RX_LIGHT_DIFF[4] = { 1.0f, 1.0f, 1.0f, 1.0f };
const GLfloat RX_LIGHT_SPEC[4] = { 0.4f, 0.4f, 0.4f, 1.0f };

const GLfloat RX_GREEN[] = { 0.1f, 1.0f, 0.1f, 1.0f };
const GLfloat RX_RED[]   = { 1.0f, 0.1f, 0.1f, 1.0f };
const GLfloat RX_BLUE[]  = { 0.1f, 0.1f, 1.0f, 1.0f };
const GLfloat RX_WHITE[] = { 1.0f, 1.0f, 1.0f, 1.0f };

const GLfloat RX_FOV = 45.0f;
const double RX_SPACE_INC = 0.002;

const int RX_SAVE_IMAGE_SPACING = 1;
const int RX_DISPLAY_STEP = 2000;

// 計算結果出力のデフォルトフォルダ
const string RX_DEFAULT_RESULT_DIR = "result/";
const string RX_DEFAULT_IMAGE_DIR  = RX_DEFAULT_RESULT_DIR+"images/";
const string RX_DEFAULT_MESH_DIR   = RX_DEFAULT_RESULT_DIR+"mesh/";
const string RX_DEFAULT_DATA_DIR   = RX_DEFAULT_RESULT_DIR+"data/";

// 設定ファイル
extern rxINI *g_pINI;

// 設定ファイルへの保存用
double g_fTBTran[3] = {0, 0, -5};	//!< 視点移動用トラックボールの平行移動量
double g_fTBQuat[4] = {1, 0, 0, 0};	//!< 視点移動用トラックボールの回転量

// FTGL
#define FONT_FILE "Inconsolata.ttf"
static FTPixmapFont* g_pFont;
	
int g_iIterations = 1;				//!< 修正反復回数
double g_fEta = 0.0;				//!< 密度変動率

// 描画
rxGLSL g_glslPointSprite;			//!< GLSLを使った描画
rxGLSL g_glslFresnel;				//!< GLSLを使った描画

// 時間計測
rxTimer g_TimerFPS;					//!< FPS測定用タイマー

// 平均計算時間計測
rxTimerAvg g_Time;

double g_fTotalTime;
double g_fAvgTime;
int g_iTimeCount;

double g_fCoefEt = 0.00019;
double g_fMaxEt = 0.02;
double g_fWaveletScale = 3.0;
double g_fMaxEnergySpectrum = 0.0;
double g_fMaxWaveletTurb = 0.0;
double g_fESScale = 0.01;
double g_fVCCoef = 0.0005;
double g_fCoefTurb = 100.0;
double g_fCoefTurbForMesh = 10.0;
double g_fEtCri = 1.0;

double g_fNoiseScl = 300.0;
double g_fNoiseEthr = 0.01;
double g_fNoiseMag = 100.0;

//! 追加::マウスピック
vector<int> g_vSelectedVertices;
int g_iPickedObj = -1;
double g_fPickDist = 1.0;	//!< ピックされた点までの距離


//-----------------------------------------------------------------------------
// 関数プロトタイプ宣言
//-----------------------------------------------------------------------------
int DrawAxis(double len, double line_width = 5.0);
void DrawStrings(vector<string> &static_str, int w, int h, int offsetx, int offsety);
void DrawStringsBottom(vector<string> &static_str, int w, int h, int offsetx, int offsety);

void CalVertexNormalsFromVBO(GLuint vrts_vbo, GLuint tris_vbo, GLuint nrms_vbo, uint nvrts, uint ntris);
static int LoadGLTextureCV(const string &fn, GLuint &tex_name, bool mipmap, bool compress);
//void ReadOBJ(const string filename, rxPolygons &polys, Vec3 cen, Vec3 ext, Vec3 ang);

RXREAL CalDeterminant3x3(const RXREAL *m);
void CalInverse3x3(const RXREAL *m, RXREAL *invm);

unsigned char* ReadImageFile(const std::string &fn, int &w, int &h, int &c);
int WriteImageFile(const std::string &fn, unsigned char *img, int w, int h, int c, int quality = 200);
bool SaveFrameBuffer(const string &fn, int w, int h);
bool SaveTexture(const string &fn, rxTexObj2D &tex_obj, int jpeg_quality = 200);
bool ReadTexture(const string &fn, rxTexObj2D &tex_obj);
//int LoadGLTexture(const string &fn, GLuint &tex_name, bool mipmap, bool compress);
bool LoadTPNG(const string &fn, GLuint &tex_name, int &w, int &h, int &c);
bool LoadTPNG(const string &fn, GLuint &tex_name);
bool LoadCubeMap(rxCubeMapData &cube_map, string base, string ext);




//-----------------------------------------------------------------------------
// rxFlWindowクラスの実装
//-----------------------------------------------------------------------------

//! コンストラクタ
rxFlGLWindow::rxFlGLWindow(int x_, int y_, int w_, int h_, const char* l, void *parent)
	: Fl_Gl_Window(x_, y_, w_, h_, l), m_iWinW(w_), m_iWinH(h_)
{
	m_pParent = (rxFlWindow*)parent;
	m_bAnimation = false;
	resizable(this);
	end();


	// 描画フラグ
	m_iDraw = 0;

	// フォントサイズ
	m_ulFontSize = 14;


	//
	// シミュレーション用変数の初期化
	//
	m_pPS = 0;
	m_fDt = 0.005;
	m_fGravity = 9.80665;

	// アニメーションステップ
	m_iCurrentStep = 0;
	m_bPause = false;

	// シーンファイルの読み取り
	m_Scene.ReadSceneFiles();

	// シーンタイトルリスト
	m_vSceneTitles = m_Scene.GetSceneTitles();

	// 現在のシーン
	m_iCurrentSceneIdx = m_Scene.GetCurrentSceneIdx();

	m_iSimuSetting = 0;

	m_iSaveImageSpacing = -1;

	m_iColorType = rxParticleSystemBase::RX_RAMP;
	m_fVScale = 0.02;

	// パーティクル情報出力
	m_strSphOutputName0  = RX_DEFAULT_DATA_DIR+"sph_setting.dat";
	m_strSphOutputHeader = RX_DEFAULT_DATA_DIR+"sph_particles_";

	// パーティクル情報入力
	m_strSphInputName0  = RX_DEFAULT_DATA_DIR+"sph_setting.dat";
	m_strSphInputHeader = RX_DEFAULT_DATA_DIR+"sph_particles_";

	// 固体移動フラグ
	m_bSolidMove = false;

	// 固体の動き
	m_vMovePos[0] = Vec3(0.0);
	m_vMovePos[1] = Vec3(0.0);
	m_fMoveMaxVel = 0.0;
	m_bMoveSolid = false;
	m_iMoveStart = -1;

	m_fIceFlag = 0;

	//
	// メッシュ
	//
	m_iNumVrts = 0; m_iNumTris = 0;		// 生成されたメッシュの頂点数とメッシュ数
	m_iVertexStore = 5;					// サンプリングボリューム数に対する予想される頂点数(nx*ny*store)

	m_fMeshThr = 400.0;					// 陰関数メッシュ化時の閾値
	m_iMeshMaxN = 128;					// メッシュ化グリッド数(境界がもっとも長い軸方向の分割数)
	m_iMeshN[0] = 64;					// メッシュ化グリッド数
	m_iMeshN[1] = 64;
	m_iMeshN[2] = 64;

	m_vMeshBoundaryExt = Vec3(1.0);		// メッシュ境界ボックスの各辺の長さの1/2
	m_vMeshBoundaryCen = Vec3(0.0);		// メッシュ境界ボックスの中心座標

	m_iSolidDraw = 2;

	m_uVrtVBO = 0;
	m_uTriVBO = 0;
	m_uNrmVBO = 0;
	m_iDimVBO = 4;

	// メッシュ出力
	m_iSaveMeshSpacing = RX_SAVE_IMAGE_SPACING;

	// 背景画像
	m_bUseCubeMap = false;		// キューブマップ使用フラグ

	g_fTotalTime = 0;
	g_fAvgTime = 0;
	g_iTimeCount = 0;

	m_pMCMeshCPU = 0;
	m_pMCMeshGPU = 0;
	m_pSSMCPU = 0;				// Screen Space Mesh
	m_pSSMGPU = 0;				// Screen Space Mesh by CUDA
	m_iDepthFiltering = 1;		// 平滑化(デプスマップ0x01, 輪郭0x02)
	m_iSSDebugOutput = 0;
	m_iTriangulationMethod = RXM_MC_GPU;

	m_fSpacing = 2;				// デプスマップのサンプリング間隔
	m_fPrtRad = 0.01;			// パーティクルの半径
	m_fZmax = 1.8*m_fPrtRad;	// 輪郭となるデプス差の閾値
	m_iNfilter = 1;				// デプス値平滑化のフィルタサイズ
	m_iNiters = 3;				// 輪郭平滑化の反復回数


	m_iPickedParticle = -1;
	

	// 設定ファイル
	if(g_pINI){
		//g_pINI->Set("gl", "draw", &m_iDraw, m_iDraw);
		g_pINI->Set("gl", "font", &m_ulFontSize, m_ulFontSize);

		g_pINI->Set("trackball", "tx",  &g_fTBTran[0],  0.0);
		g_pINI->Set("trackball", "ty",  &g_fTBTran[1],  0.0);
		g_pINI->Set("trackball", "tz",  &g_fTBTran[2], -4.0);
		g_pINI->Set("trackball", "q0",  &g_fTBQuat[0],  1.0);
		g_pINI->Set("trackball", "q1",  &g_fTBQuat[1],  0.0);
		g_pINI->Set("trackball", "q2",  &g_fTBQuat[2],  0.0);
		g_pINI->Set("trackball", "q3",  &g_fTBQuat[3],  0.0);

		g_pINI->Set("sph", "scene",  &m_iCurrentSceneIdx, -1);
		g_pINI->Set("sph", "mesh_threshold", &m_fMeshThr, m_fMeshThr);
		g_pINI->Set("sph", "setting", &m_iSimuSetting, 0);

		g_pINI->Set("sph", "et_coef", &g_fCoefEt, g_fCoefEt);
		g_pINI->Set("sph", "et_max", &g_fMaxEt, g_fMaxEt);
		g_pINI->Set("sph", "et_scale", &g_fWaveletScale, g_fWaveletScale);

		g_pINI->Set("sph", "sp_et_coef", &g_fCoefTurb, g_fCoefTurb);
		g_pINI->Set("sph", "sp_et_coef_mesh", &g_fCoefTurbForMesh, g_fCoefTurbForMesh);
		g_pINI->Set("sph", "sp_et_cri", &g_fEtCri, g_fEtCri);

		g_pINI->Set("sph", "noise_scale", &g_fNoiseScl, g_fNoiseScl);
		g_pINI->Set("sph", "noise_ethr", &g_fNoiseEthr, g_fNoiseEthr);
		g_pINI->Set("sph", "noise_emag", &g_fNoiseMag, g_fNoiseMag);

		g_pINI->Set("sph", "solid_mesh", &m_iSolidDraw, m_iSolidDraw);

		g_pINI->Set("sph", "setting", &m_iSimuSetting, 0);
	}
}

//! デストラクタ
rxFlGLWindow::~rxFlGLWindow()
{
	if(m_pMCMeshGPU) delete m_pMCMeshGPU;
	if(m_pMCMeshCPU) delete m_pMCMeshCPU;
	if(m_pSSMCPU) delete m_pSSMCPU;
	if(m_pSSMGPU) delete m_pSSMGPU;

	if(m_uVrtVBO) glDeleteBuffers(1, &m_uVrtVBO);
	if(m_uTriVBO) glDeleteBuffers(1, &m_uTriVBO);
	if(m_uNrmVBO) glDeleteBuffers(1, &m_uNrmVBO);

	if(m_fIceFlag) CuFreeArray(m_fIceFlag);

	if(m_pPS) delete m_pPS;

	if(m_ice) delete m_ice;
}


/*! 
 * GLの初期化関数
 */
void rxFlGLWindow::InitGL(void)
{
	// MARK:InitGL
	make_current();
	static int init = false;
	if(init) return;

	init = true;

	CuDeviceProp();
	CuSetDevice(0);

	CheckAndMakeDir(RX_DEFAULT_IMAGE_DIR);
	CheckAndMakeDir(RX_DEFAULT_DATA_DIR);
	CheckAndMakeDir(RX_DEFAULT_MESH_DIR);

	RXCOUT << "OpenGL Ver. " << glGetString(GL_VERSION) << endl;

	GLenum err = glewInit();
	if(err == GLEW_OK){
		RXCOUT << "GLEW OK : Glew Ver. " << glewGetString(GLEW_VERSION) << endl;
	}
	else{
		RXCOUT << "GLEW Error : " << glewGetErrorString(err) << endl;
	}

	// OpenGL拡張のバージョンチェック
	if(!glewIsSupported("GL_VERSION_2_0 " 
						"GL_ARB_pixel_buffer_object "
						"GL_EXT_framebuffer_object "
						"GL_ARB_multitexture "
						"GL_ARB_vertex_buffer_object "
						)){
		RXCOUT << "ERROR: Support for necessary OpenGL extensions missing." << endl;
		exit(1);
	}

	GLint buf, sbuf;
	glGetIntegerv(GL_SAMPLE_BUFFERS, &buf);
	RXCOUT << "number of sample buffers is " << buf << endl;
	glGetIntegerv(GL_SAMPLES, &sbuf);
	RXCOUT << "number of samples is " << sbuf << endl;

	m_fBGColor[0] = 1.0; m_fBGColor[1] = 1.0; m_fBGColor[2] = 1.0;
	m_iMouseButton = 0;

	glClearColor((GLfloat)m_fBGColor[0], (GLfloat)m_fBGColor[1], (GLfloat)m_fBGColor[2], 1.0f);
	glClearDepth(1.0f);

	glEnable(GL_DEPTH_TEST);
	glEnable(GL_CULL_FACE);

	glEnable(GL_AUTO_NORMAL);
	glEnable(GL_NORMALIZE);

	// 光源設定
	glEnable(GL_LIGHT0);
	glEnable(GL_LIGHTING);
	glLightfv(GL_LIGHT0, GL_POSITION, RX_LIGHT0_POS);
	glLightfv(GL_LIGHT0, GL_DIFFUSE,  RX_LIGHT_DIFF);
	glLightfv(GL_LIGHT0, GL_SPECULAR, RX_LIGHT_SPEC);
	glLightfv(GL_LIGHT0, GL_AMBIENT,  RX_LIGHT_AMBI);

	glShadeModel(GL_SMOOTH);

	glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
	glEnable(GL_COLOR_MATERIAL);

	// 視点初期化
	//InitView();

	// デフォルト材質
	m_matPoly.name      = "";
	m_matPoly.diffuse   = Vec4(0.2, 0.2, 1.0, 1.0);
	m_matPoly.specular  = Vec4(0.3, 0.3, 0.3, 1.0);
	m_matPoly.ambient   = Vec4(0.1, 0.1, 0.1, 1.0);
	m_matPoly.color     = Vec4(0.0);
	m_matPoly.emission  = Vec4(0.0);
	m_matPoly.shininess = 10.0;
	m_matPoly.illum = 2;
	m_matPoly.tex_file = "";
	m_matPoly.tex_name = 0;

	// キューブマップの読み込み
//	if(LoadCubeMap(m_CubeMap, "texture/terragen0_", ".png")){
//	if(LoadCubeMap(m_CubeMap, "texture/cubemap_", ".png")){
//	if(LoadCubeMap(m_CubeMap, "texture/nvlobby_new_", ".png")){
	if(LoadCubeMap(m_CubeMap, "texture/nvlobby_new_gauss_", ".png")){
//	if(LoadCubeMap(m_CubeMap, "texture/IceLake_", ".jpg")){	//jpgだとうまくいかない？
		m_bUseCubeMap = true;
	}
	else{
		RXCOUT << "error : can't load the cube map" << endl;
	}

	// フォント読み込み
	SetupFonts(FONT_FILE);
	//RXCOUT << "char width : " << glutBitmapWidth(GLUT_BITMAP_HELVETICA_12, 'W') << endl;

	//g_InletLine.span = -1;
	m_bsSimuSetting.set(ID_SPH_ANISOTROPIC, false);

	// トラックボール初期姿勢
	m_tbView.SetScaling(-5.0);
	//m_tbView.SetTranslation(0.0, -2.0);

	m_tbView.SetTranslation(g_fTBTran[0], g_fTBTran[1]);
	m_tbView.SetScaling(g_fTBTran[2]);
	m_tbView.SetQuaternion(g_fTBQuat);

	int settings = m_iSimuSetting;
	for(int i = 0; i < ID_SPH_END; ++i){
		int flag = settings & 1;

		m_bsSimuSetting.set(i, (flag == 1) ? true : false);

		settings >>= 1;
	}
	m_bsSimuSetting.set(ID_SPH_SPS_TURB, false);

	m_bsSimuSetting.set(ID_HEAT, false);	//追加
	m_bsSimuSetting.set(ID_SM, false);		//追加
	m_bsSimuSetting.set(ID_ICE, false);		//追加

	// 描画フラグ初期化
	m_iDraw = 0;
	m_iDraw |= RXD_BBOX;
	m_iDraw |= RXD_PARTICLE;
	m_iDraw |= RXD_SOLID;
	m_iDraw |= RXD_PARAMS;
	if(m_bsSimuSetting.at(ID_SPH_MESH)) m_iDraw |= RXD_MESH;

	m_iDrawPS = RXP_POINTSPRITE;


	// カレントのシーン設定
	m_Scene.SetCurrentScene(m_iCurrentSceneIdx);

	// SPH初期化
	InitSPH(m_Scene);

	//相変化オブジェクトの初期化
	InitIceObj();					

	// GLSLのコンパイル
	g_glslPointSprite = CreateGLSL(ps_vs, ps_fs, "point sprite");
	g_glslFresnel     = CreateGLSL(fresnel2_vs, fresnel2_fs, "fresnel");
	//g_glslFresnel     = CreateGLSL(fresnel_vs, fresnel_fs, "fresnel");


	m_pParent->UpdateMenuState();
}

/*!
 * フォントの設定
 * @param[in] file フォントファイルパス
 */
int rxFlGLWindow::SetupFonts(const char* file)
{
	g_pFont = new FTPixmapFont(file);
	if(g_pFont->Error()){
		RXCOUT << "Failed to open font " << file << endl;
		delete g_pFont;
		g_pFont = 0;
		return 1;
	}

	g_pFont->FaceSize(m_ulFontSize);

	return 0;
}
/*! 
 * リサイズイベント処理関数
 * @param[in] w キャンバス幅(ピクセル数)
 * @param[in] h キャンバス高さ(ピクセル数)
 */
void rxFlGLWindow::Resize(int w, int h)
{
	m_iWinW = w;
	m_iWinH = h;

	cout << "m_iWinW = " << m_iWinW << " m_iWinH = " << m_iWinH << endl;

	glViewport(0, 0, w, h);
	m_tbView.SetRegion(w, h);

	// 透視変換行列の設定
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	Projection();

	// モデルビュー変換行列の設定
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
}

/*!
 * 透視投影変換
 */
void rxFlGLWindow::Projection(void)
{
	gluPerspective(RX_FOV, (float)m_iWinW/(float)m_iWinH, 0.2f, 1000.0f);
	//glOrtho(-1, 1, -1, 1, -1, 1);
}

/*!
 * 再描画命令
 */
void rxFlGLWindow::ReDisplay(void)
{
	redraw();
}

/*!
 * 視点の初期化
 */
void rxFlGLWindow::InitView(void)
{
	double q[4] = {1, 0, 0, 0};
	m_tbView.SetQuaternion(q);
	m_tbView.SetScaling(-6.0);
	m_tbView.SetTranslation(0.0, 0.0);
}

/*!
 * 再描画イベント処理関数
 */
void rxFlGLWindow::Display(void)
{
	make_current();

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glViewport(0, 0, m_iWinW, m_iWinH);

	// 透視変換行列の設定
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(RX_FOV, (float)m_iWinW/(float)m_iWinH, 0.001f, 1000.0f);

	// モデルビュー行列初期化
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glPushMatrix();

	//glScalef(10.0, 10.0, 10.0);

	// マウスによる回転・平行移動の適用
	m_tbView.Apply();

	if(m_bUseCubeMap && (m_iDraw & RXD_REFRAC)){
		glPushMatrix();

		glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
		glColor3d(0.0, 0.0, 0.0);
	
		glDisable(GL_CULL_FACE);
		//TrackballApplyTranslation();
		//TrackballApplyRotation();
		DrawCubeMap(m_CubeMap, 100.0);
		glPopMatrix();
	}


	glPushMatrix();
		
	// OpenGL投影変換行列取得
	glGetDoublev(GL_PROJECTION_MATRIX, m_fProjectionMatrix);

	// OpenGLモデルビュー変換行列取得
	glGetDoublev(GL_MODELVIEW_MATRIX, m_fModelviewMatrix);

	RenderSphScene();

	glPopMatrix();

	glPopMatrix();

	//矩形描画
	if( m_ht_vStartPoint[0] != 0.0 && m_ht_vStartPoint[1] != 0.0
	 && m_ht_vEndPoint[0] != 0.0   && m_ht_vEndPoint[1] != 0.0 )
	{
		DrawRubber(1, m_ht_vStartPoint, m_ht_vEndPoint, m_iWinW, m_iWinH);
	}

	// 画面文字列描画
	if(m_iDraw & RXD_PARAMS){
		glColor3d(0.0, 0.0, 0.5);
		DrawStringsBottom(SetDrawString(), m_iWinW, m_iWinH, 0, 0);
	}
}

void rxFlGLWindow::DisplayForPick_s(void* x)
{
	((rxFlGLWindow*)x)->DisplayForPick();
}
void rxFlGLWindow::DisplayForPick(void)
{
	glPushMatrix();

	// マウスによる回転・平行移動の適用
	m_tbView.Apply();

	glDisable(GL_LIGHTING);
	glColor4d(0.0, 0.0, 1.0, 1.0);
	glEnable(GL_POINT_SMOOTH);
	glPointSize(4.0);

	int pnum = m_pPS->GetNumParticles();
	RXREAL *data = 0;
	data = m_pPS->GetArrayVBO(rxParticleSystemBase::RX_POSITION);
	int k = 0;
	for(int i = 0; i < pnum; ++i){
		glLoadName(i);
		glBegin(GL_POINTS);
		glVertex3d(data[k], data[k+1], data[k+2]);
		glEnd();
		k += DIM;
	}

	glPopMatrix();
}

void rxFlGLWindow::Projection_s(void *x)
{
	((rxFlGLWindow*)x)->Projection();
}

//粒子選択
void rxFlGLWindow::PickParticle(int x, int y)
{
	rxGLPick pick;
	pick.Set(rxFlGLWindow::DisplayForPick_s, this, rxFlGLWindow::Projection_s, this);
	
	int hit = pick.Pick(x, y);
	if(hit >= 1){
		m_iPickedParticle = hit-1;
	}
	else
	{
		m_iPickedParticle = -1;
	}
	cout << "picked partcile : " << m_iPickedParticle << endl;
}


/*!
 * マウスイベント処理関数
 * @param[in] button マウスボタン(FL_LEFT_MOUSE,FL_MIDDLE_MOUSE,FL_RIGHT_MOUSE)
 * @param[in] state マウスボタンの状態(1:down, 0:up)
 * @param[in] x,y マウス座標(スクリーン座標系)
 */
void rxFlGLWindow::Mouse(int button, int state, int x, int y)
{
	make_current();
	m_iKeyMod = (Fl::event_state(FL_SHIFT) ? 2 : (Fl::event_state(FL_CTRL) ? 3 : 1));
	int mask = 0x01 << (button-1);
	m_iMouseButton = (state ? (m_iMouseButton | mask) : (m_iMouseButton & ~mask));

	if(button == FL_LEFT_MOUSE){
		if(!m_bSolidMove){
			if(state)
			{	// ボタンダウン				
				PickParticle(x, y);			//粒子取得

				if( m_iPickedParticle != -1 )
				{
					cout << "vertex " << m_iPickedParticle << endl;

					vector<int> vrts;
					vrts.push_back(m_iPickedParticle);

					g_vSelectedVertices.resize(1);

					// 視点からピック点までの距離計算
					Vec3 ray_from;
					Vec3 init_pos = Vec3(0.0);
					m_tbView.CalLocalPos(ray_from, init_pos);
					RXREAL *p  = m_pPS->GetArrayVBO(rxParticleSystemBase::RX_POSITION);
					Vec3 pos = Vec3(	p[DIM*m_iPickedParticle+0], 
										p[DIM*m_iPickedParticle+1],
										p[DIM*m_iPickedParticle+2]
									);

					g_fPickDist = length(pos-ray_from);
				}
				else
				{
					m_tbView.Start(x, y, m_iKeyMod);
					ClearPick();
					ClearRect();
				}
			}
			else
			{	// ボタンアップ
				m_tbView.Stop(x, y);

				m_tbView.GetTranslation(g_fTBTran);
				m_tbView.GetScaling(g_fTBTran[2]);
				m_tbView.GetQuaternion(g_fTBQuat);
				//m_bsSimuSetting.set(ID_SPH_INLET, false);
				m_iSimuSetting = (int)m_bsSimuSetting.to_ulong();

				ClearPick();
			}
			ClearRect();
		}
	}
	else if(button == FL_MIDDLE_MOUSE){

	}
	else if(button == FL_RIGHT_MOUSE){

	}
}

/*!
 * 氷粒子選択のリセット
 */
void rxFlGLWindow::ClearPick(void)
{//	cout<< "clearPick" << endl;
	g_vSelectedVertices.clear();

	if( m_iPickedParticle != -1 )
	{
		if(m_iceObj == 0)	return;
		if(m_iceObj->GetParticleNum() <= m_iPickedParticle){	return;	}	//融解のみの実験のときに必要になる．

		//クラスタの情報を更新
		for( int i = 0; i < m_iceObj->GetPtoCIndx(m_iPickedParticle); i++ )
		{
			int cIndx = m_iceObj->GetPtoC(m_iPickedParticle, i, 0);
			int oIndx = m_iceObj->GetPtoC(m_iPickedParticle, i, 1);

			if(cIndx == -1 || oIndx == -1)
			{
				continue;
			}

			m_iceObj->GetMoveObj(cIndx)->UnFixVertex(oIndx);	//粒子を非選択に
			m_iceObj->GetOrientedObj(cIndx)->UnFixVertex(oIndx);
		}

		m_iPickedParticle = -1;
	}
}

/*!
 * 矩形選択領域のリセット
 */
void rxFlGLWindow::ClearRect(void)
{
	m_ht_vStartPoint = Vec2(0.0, 0.0);
	m_ht_vEndPoint = Vec2(0.0, 0.0);
	m_ht_vSelectedVertices.clear();
	m_ht_bRectFlag = false;
}

/*!
 * モーションイベント処理関数(マウスボタンを押したままドラッグ)
 * @param[in] x,y マウス座標(スクリーン座標系)
 */
void rxFlGLWindow::Motion(int x, int y)
{
	if(x < 0 || y < 0) return;
	make_current();

	static int x0 = -1, y0 = -1;
	if(x0 == -1 || y0 == -1){
		x0 = x;
		y0 = y;
	}

	if(m_bSolidMove)
	{
		//if(((m_iMouseButton >> 3) == GLUT_LEFT_BUTTON) && ((m_iMouseButton & 1) == GLUT_DOWN)){
		//	double dm = 0.4/(double)m_iWinW;

		//	m_pPS->MoveSphereObstacle(0, Vec3(0.0, dm*(y0-y), dm*(x-x0)));
		//}

		Vec3 ray_from, ray_to;
		Vec3 init_pos = Vec3(0.0);
		m_tbView.CalLocalPos(ray_from, init_pos);
		m_tbView.GetRayTo(x, y, RX_FOV, ray_to);

		Vec3 spherePos = m_pPS->GetSphereObstaclePos();
		float dist = length(spherePos-ray_from);
		Vec3 dir = Unit(ray_to-ray_from);	// 視点からマウス位置へのベクトル
		Vec3 new_pos = ray_from+dir*dist;

		m_pPS->MoveSphereObstacle(0, new_pos-spherePos);
	}
	else if( m_ht_bRectFlag )
	{
		m_ht_vEndPoint = Vec2(x, y);
	}
	else
	{
		if(g_vSelectedVertices.empty())
		{
			m_tbView.Motion(x, y);
		}
		//クラスタに属する粒子の移動
		else if( m_iPickedParticle != -1 )
		{
			if(m_iceObj == 0)	return;
			if(m_iceObj->GetParticleNum() <= m_iPickedParticle){	return;	}	//融解のみの実験のときに必要になる．

			Vec3 ray_from, ray_to;
			Vec3 init_pos = Vec3(0.0);
			m_tbView.CalLocalPos(ray_from, init_pos);
			m_tbView.GetRayTo(x, y, RX_FOV, ray_to);

			Vec3 dir = Unit(ray_to-ray_from);	// 視点からマウス位置へのベクトル
			Vec3 new_pos = ray_from+dir*g_fPickDist;
			int v = g_vSelectedVertices[0];

			//粒子が属するクラスタの数値を更新
			for( int i = 0; i < m_iceObj->GetPtoCIndx(m_iPickedParticle); i++ )
			{
				int cIndx = m_iceObj->GetPtoC(m_iPickedParticle, i, 0);
				int oIndx = m_iceObj->GetPtoC(m_iPickedParticle, i, 1);
				if(cIndx == -1 || oIndx == -1) continue;
	
				m_iceObj->GetMoveObj(cIndx)->FixVertex(oIndx, new_pos);
				m_iceObj->GetMoveObj(cIndx)->SetCurrentPos(oIndx, new_pos);

				//楕円体の実験のために追加
				m_iceObj->GetOrientedObj(cIndx)->FixVertex(oIndx, new_pos);
				m_iceObj->GetOrientedObj(cIndx)->SetCurrentPos(oIndx, new_pos);
			}

			//粒子の移動
			RXREAL *p  = m_pPS->GetArrayVBO(rxParticleSystemBase::RX_POSITION);
			p[DIM*m_iPickedParticle+0] = new_pos[0];
			p[DIM*m_iPickedParticle+1] = new_pos[1];
			p[DIM*m_iPickedParticle+2] = new_pos[2];

			m_pPS->SetArrayVBO(rxParticleSystemBase::RX_POSITION, p, 0, m_pPS->GetNumParticles());
		}
	}
	//RXCOUT << "        " << (m_iMouseButton >> 3) << ", " << (m_iMouseButton & 1) << endl;

	x0 = x;
	y0 = y;

	redraw();
}

/*!
 * モーションイベント処理関数(マウスボタンを押さない移動)
 * @param[in] x,y マウス座標(スクリーン座標系)
 */
void rxFlGLWindow::PassiveMotion(int x, int y)
{
	if(x < 0 || y < 0) return;
	//make_current();
	//redraw();
}

/*!
 * キーボードイベント処理関数
 * @param[in] key キーの種類
 * @param[in] x,y キーが押されたときのマウス座標(スクリーン座標系)
 */
void rxFlGLWindow::Keyboard(int key, int x, int y)
{
	make_current();
	m_iKeyMod = (Fl::event_state(FL_SHIFT) ? 2 : (Fl::event_state(FL_CTRL) ? 3 : 1));

	switch(key){
	case 'i':
		InitView();
		break;

	case 'S':	// 画面の画像保存
		SaveDisplay( ((m_iCurrentStep >= 0) ? m_iCurrentStep : 0) );
		break;

	case 'g':
		SwitchFullScreen();
		break;

	// 
	// シミュレーション設定
	// 
	case '*':
		AddSphere();
		break;

	case FL_Left:
	case FL_Right:
	case FL_Up:
	case FL_Down:
	case FL_Home:
	case FL_End:
	case FL_Page_Up:
	case FL_Page_Down:
	//case FL_Tab:
	case FL_BackSpace:
	case FL_Scroll_Lock:
	case FL_Delete:
		SpecialKey(key, x, y);
		return;

	default:
		break;
	}


	m_pParent->UpdateMenuState();
	ReDisplay();
}

/*!
 * 特殊キーボードイベント処理関数
 * @param[in] key キーの種類
 * @param[in] x,y キーが押されたときのマウス座標(スクリーン座標系)
 */
void rxFlGLWindow::SpecialKey(int key, int x, int y)
{
	if(m_iKeyMod == GLUT_ACTIVE_CTRL){
		if(key == FL_Left){
			RXCOUT << "clockwise" << endl;
			m_tbView.AddRotation(90, 0, 1, 0);
		}
		else if(key == FL_Right){
			RXCOUT << "counterclockwise" << endl;
			m_tbView.AddRotation(90, 0, -1, 0);
		}
		else if(key == FL_Up){
			RXCOUT << "up" << endl;
			m_tbView.AddRotation(90, 1, 0, 0);
		}
		else if(key == FL_Down){
			RXCOUT << "down" << endl;
			m_tbView.AddRotation(90, -1, 0, 0);
		}

		ReDisplay();
		return;
	}

	if(m_bSolidMove){
		if(key == FL_Left){
			m_pPS->MoveSphereObstacle(0, Vec3( 0.01, 0.0, 0.0));
		}
		else if(key == FL_Right){
			m_pPS->MoveSphereObstacle(0, Vec3(-0.01, 0.0, 0.0));
		}
		else if(key == FL_Up){
			m_pPS->MoveSphereObstacle(0, Vec3(0.0, 0.0, -0.01));
		}
		else if(key == FL_Down){
			m_pPS->MoveSphereObstacle(0, Vec3(0.0, 0.0,  0.01));
		}
	}
	else{
		if(key == FL_Left){
			m_fVScale -= (m_fVScale <= 0 ? 0.0 : 0.005);
			RXCOUT << "vscale : " << m_fVScale << endl;
		}
		else if(key == FL_Right){
			m_fVScale += (m_fVScale >= 1 ? 0.0 : 0.005);
			RXCOUT << "vscale : " << m_fVScale << endl;
		}
		else if(key == FL_Up){
			//g_fCoefEt += (g_fCoefEt >= 100 ? 0.0 : 0.1);
			g_fWaveletScale += (g_fWaveletScale >= 1 ? 0.0 : 0.005);
		}
		else if(key == FL_Down){
			//g_fCoefEt -= (g_fCoefEt <= 0 ? 0.0 : 0.1);
			g_fWaveletScale -= (g_fWaveletScale <= 0 ? 0.0 : 0.005);
		}
	}

	m_pParent->UpdateMenuState();
	ReDisplay();
}

/*!
 * アイドルコールバック関数
 * @param[in] x ユーザ定義変数
 */
void rxFlGLWindow::OnIdle_s(void* x)
{
	((rxFlGLWindow*)x)->Idle();
}
void rxFlGLWindow::Idle(void)
{
	make_current();

	if(m_iCurrentStep%RX_DISPLAY_STEP == 0 && m_iCurrentStep != 0){
		RXCOUT << "step " << m_iCurrentStep << endl;
	}
	// 出力用ステップ数
	int stp = m_iCurrentStep-2;
	stp = (m_iCurrentStep == 0) ? 0 : stp;
	stp = (stp < 0) ? 0 : stp;

	//
	// 前のフレームの描画を画像ファイルに保存
	//
	if(m_iSaveImageSpacing > 0){
		if(stp%m_iSaveImageSpacing == 0){
			SaveDisplay(stp);
		}
	}

	//
	// 前のフレームのデータをファイル保存
	//
	if(m_bsSimuSetting.at(ID_SPH_OUTPUT)){
		m_pPS->OutputParticles(CreateFileName(m_strSphOutputHeader, "dat", stp, 5));
	}

	g_TimerFPS.Reset();
	g_TimerFPS.Start();
	RXTIMER_RESET;

	if(m_iCurrentStep >= m_iMoveStart) m_bMoveSolid = true;

	//
	// 固体を動かす
	//
	if(m_bMoveSolid){
		Vec3 spos = m_pPS->GetSphereObstaclePos(0);

		static bool solid_init = false;
		static Vec3 sdir = Vec3(0.0);
		static double svel = 0.0;
		double acc = 0.1;
		if(!solid_init){
			if(norm2(m_vMovePos[0]-spos) > norm2(m_vMovePos[1]-spos)){
				sdir = Unit(m_vMovePos[0]-spos);
			}
			else{
				sdir = Unit(m_vMovePos[1]-spos);
			}
			solid_init = true;
		}

		if(dot(m_vMovePos[0]-spos, m_vMovePos[1]-spos) > 0){
			// 反転
			sdir *= -1;
		}
		
		svel = m_fMoveMaxVel;

		m_pPS->MoveSphereObstacle(0, sdir*svel*m_fDt);

	}

	//
	// シミュレーションタイムステップを進める
	//
	if(m_bsSimuSetting.at(ID_SPH_INPUT)){
		// パーティクル位置の入力
		if(!m_pPS->InputParticles(CreateFileName(m_strSphInputHeader, "dat", m_iCurrentStep, 5))){
			SwitchIdle(0);
			m_bsSimuSetting.set(ID_SPH_INPUT, false);
		}
	}
	else if( m_bMode == MODE_SPH )
	{
		RXTIMER_RESET;
		StepPS(m_fDt);
	}

	//
	//追加：：シミュレーションタイムステップを進める　熱処理
	//
	if( m_bMode == MODE_HEAT )
	{
		RXTIMER_RESET;
		
		StepHeatTransfer();
		StepPS(m_fDt);
//		StepSM(m_fDt);
	}
	//
	//追加：：シミュレーションタイムステップを進める　SM法
	//
	else if( m_bMode == MODE_SM )
	{
		RXTIMER_RESET;

//		StepPS(m_fDt);
//		StepSM(m_fDt);						//SM法　粒子位置のフィードバック，計算，線形補間
	}
	//
	//追加：：シミュレーションタイムステップを進める　氷処理
	//
	else if(m_bMode == MODE_ICE)
	{
		StepTimeEvent();			RXTIMER("StepTimeEvent");		//タイムイベント
		
		StepPS(m_fDt);				RXTIMER("StepPS");				//液体運動
		StepIceObj();												//固体運動
		StepHeatTransfer();			RXTIMER("StepHeatTransfer");	//熱処理
		StepIceStructure();			RXTIMER("StepIceStructure");	//相変化
	}

	//
	//追加：：粒子の色設定
	//
	StepParticleColor();
	RXTIMER("StepParticleColor");

	if(m_bsSimuSetting.at(ID_SPH_ANISOTROPIC)){
		// 異方性カーネル
		RXTIMER_RESET;
		m_pPS->CalAnisotropicKernel();
		RXTIMER("anisotropic");
	}

	//
	// 流体表面メッシュ生成
	//
	if(m_bsSimuSetting.at(ID_SPH_MESH)){
		CalMeshSPH(m_iMeshMaxN, m_fMeshThr);
	}
	RXTIMER("mesh");

	//
	// FPSの計算
	//
	g_TimerFPS.Stop();
	if(m_iCurrentStep) ComputeFPS();

	if(m_iCurrentStep%RX_DISPLAY_STEP == 0 && m_iCurrentStep != 0){
		string tstr;
		RXTIMER_STRING(tstr);
		RXCOUT << "time : " << tstr << endl;
	}

	//
	// メッシュ保存
	//
	if(m_bsSimuSetting.at(ID_SPH_MESH)){
		if(m_bsSimuSetting.at(ID_SPH_MESH_OUTPUT) && m_iCurrentStep%m_iSaveMeshSpacing){
			SaveMesh(m_iCurrentStep, m_Poly);
		}
	}

	// 最大ステップ数を超えたらアイドルを止める
	if(m_iCurrentStep > RX_MAX_STEPS) SwitchIdle(0);

	m_iCurrentStep++;		// 現在のステップ数

	redraw();
}

/*!
 * タイマーコールバック関数
 * @param[in] x ユーザ定義変数
 */
void rxFlGLWindow::OnTimer_s(void* x)
{
	((rxFlGLWindow*)x)->Timer();
}
void rxFlGLWindow::Timer(void)
{
	Idle();
	if(m_bAnimation){
		Fl::repeat_timeout(0.033, OnTimer_s, this);
	}
}


/*!
 * アイドル関数のON/OFF
 * @param[in] on trueでON, falseでOFF
 */
bool rxFlGLWindow::SwitchIdle(int on)
{
	m_bAnimation = (on == -1) ? !m_bAnimation : (on ? true : false);
	if(m_bAnimation){
		//Fl::add_timeout(0.033, CB_Timer_s, this);
		if(!Fl::has_idle(rxFlGLWindow::OnIdle_s, this)){
			Fl::add_idle(rxFlGLWindow::OnIdle_s, this);
		}
	}
	else{
		if(Fl::has_idle(rxFlGLWindow::OnIdle_s, this)){
			Fl::remove_idle(rxFlGLWindow::OnIdle_s, this);
		}
	}
	return m_bAnimation;
}

/*!
 * フルスクリーン/ウィンドウ表示の切り替え
 * @param[in] win 1でGLキャンパスをウィンドウ内でフル化, 0でウィンドウ含めてフル化
 */
void rxFlGLWindow::SwitchFullScreen(int win)
{
	static int pos0[2] = { 0, 0 };
	static int win0[2] = { 500, 500 };
	if(win){
		// ウィンドウ内でフル化
		if(m_bFullscreen){
			fullscreen_off(pos0[0], pos0[1], win0[0], win0[1]);
			m_bFullscreen = false;
		}
		else if(m_pParent->IsFullScreen()){
			pos0[0] = x();
			pos0[1] = y();
			win0[0] = w();
			win0[1] = h();
			fullscreen();
			m_bFullscreen = true;
		}
	}
	else{
		// ウィンドウ含めてフル化
		if(m_bFullscreen){
			fullscreen_off(pos0[0], pos0[1], win0[0], win0[1]);
			m_bFullscreen = false;
			if(m_pParent->IsFullScreen()) m_pParent->SwitchFullScreen();
		}
		else{
			if(!m_pParent->IsFullScreen()) m_pParent->SwitchFullScreen();
			pos0[0] = x();
			pos0[1] = y();
			win0[0] = w();
			win0[1] = h();
			fullscreen();
			m_bFullscreen = true;
		}
	}
}

/*!
 * フルスクリーン/ウィンドウ表示の状態取得
 */
int rxFlGLWindow::IsFullScreen(void)
{
	return (m_bFullscreen ? 1 : 0);
}


/*!
 * ファイル読み込み
 * @param[in] fn ファイルパス
 */
void rxFlGLWindow::OpenFile(const string &fn)
{	cout << __FUNCTION__ << endl;

	// 3Dファイル読み込み
	RxModel::Read(fn, m_poly);

	int vertex_count = 0;
	int index_count = 0; 

	vertex_count = (int)m_poly.vertices.size(); // 総頂点数
	index_count = (int)m_poly.faces.size(); // 総ポリゴン数

	cout << "vertex = " << vertex_count << " index_count = " << index_count << endl;
	
	redraw();
}

/*!
 * ファイル書き込み
 * @param[in] fn ファイルパス
 */
void rxFlGLWindow::SaveFile(const string &fn)
{
}



/*!
 * 現在の画面描画を画像ファイルとして保存
 * @param[in] fn ファイルパス
 */
void rxFlGLWindow::SaveDisplay(const string &fn)
{
	static int count = 0;
	//string fn = CreateFileName("view_", ".png", count, 5);

	int w_ = w();
	int h_ = h();
	int c_ = 4;

	make_current();
	unsigned char* data = new unsigned char[w_*h_*c_];

	glReadPixels(0, 0, w_, h_, GL_RGBA, GL_UNSIGNED_BYTE, data);

	// 上下反転
	int stride = w_*c_;
	for(int j = 0; j < h_/2; ++j){
		for(int i = 0; i < stride; ++i){
			unsigned char tmp = data[j*stride+i];
			data[j*stride+i] = data[(h_-j-1)*stride+i];
			data[(h_-j-1)*stride+i] = tmp;
		}
	}

	string ext = GetExtension(fn);
	if(ext == "bmp"){
		WriteBitmapFile(fn, data, w_, h_, c_, RX_BMP_WINDOWS_V3);
		RXCOUT << "saved the screen image to " << fn << endl;
		count++;
	}
	else if(ext == "png"){
		WritePngFile(fn, data, w_, h_, c_);
		RXCOUT << "saved the screen image to " << fn << endl;
		count++;
	}

	delete [] data;
}


/*!
 * 現在の画面描画を画像ファイルとして保存
 * @param[in] stp 現在のステップ数(ファイル名として使用)
 */
void rxFlGLWindow::SaveDisplay(const int &stp)
{
	// MRK:SaveDisplay
	string image_name = "sph";
	string head = RX_DEFAULT_IMAGE_DIR+image_name+"_";
	string fn = CreateFileName(head, ".png", stp, 5);

	SaveDisplay(fn);
}


/*!
 * メッシュをファイルとして保存
 * @param[in] fn ファイル名
 * @param[in] polys ポリゴンオブジェクト
 */
void rxFlGLWindow::SaveMesh(const string fn, rxPolygons &polys)
{
	rxPOV pov;
	if(pov.SaveListData(fn, polys)){
		RXCOUT << "saved the mesh to " << fn << endl;
	}
}

/*!
 * メッシュをファイルとして保存
 * @param[in] stp 現在のステップ数(ファイル名として使用)
 * @param[in] polys ポリゴンオブジェクト
 */
void rxFlGLWindow::SaveMesh(const int &stp, rxPolygons &polys)
{
	SaveMesh(CreateFileName(RX_DEFAULT_MESH_DIR+"sph_", "inc", stp, 5), polys);
}


/*!
 * イベントハンドラ
 * @param[in] ev イベントID
 */
int rxFlGLWindow::handle(int e)
{
	switch(e){
	case FL_DND_ENTER:
	case FL_DND_RELEASE:
	case FL_DND_LEAVE:
	case FL_DND_DRAG:
	case FL_PASTE:
		// DnDBoxにペースト処理を渡すために，
		// これらのイベントが来たときはFl_Gl_Window::handleに処理を渡さない．
		return 1;

	default:
		break;
	}

	return Fl_Gl_Window::handle(e);
}


/*!
 * Fl_Menu_Barのコールバック関数 - Draw
 * @param[in] widget ウィジットの親クラスオブジェクト
 * @param[in] x ユーザ定義変数
 */
void rxFlGLWindow::OnMenuDraw_s(Fl_Widget *widget, void* x)
{
	Fl_Menu_Bar *menubar = (Fl_Menu_Bar*)widget;
	char picked[80];
	menubar->item_pathname(picked, sizeof(picked)-1);	// メニュー名

	string label = picked;
	string menu_name = "Draw/";
	label = label.substr(menu_name.size(), string::npos);

	((rxFlGLWindow*)x)->OnMenuDraw(1.0, label);
}
void rxFlGLWindow::OnMenuDraw(double val, string label)
{
	int idx = 0;
	while(label.find(RX_DRAW_STR[2*idx]) != 0) idx++;

	int flag = (0x01 << idx);
	m_iDraw ^= flag;

	switch(flag){
	case RXD_NORMAL:	// 法線
		m_pPS->ToggleNormalCalc((m_iDraw & RXD_NORMAL));
		break;

	case RXD_CELLS:	// 分割セル
		if(m_iDraw & RXD_CELLS){
			m_pPS->SetParticlesToCell();
			m_pPS->SetPolygonsToCell();
		}
		break;

	case RXD_REFRAC:	// 屈折描画
		// if(m_iDraw & RXD_REFRAC) m_iDraw |= RXD_FOAM;
		// m_bSphMesh = (m_iDraw & RXD_REFRAC);
		break;

	case RXD_MESH:	// メッシュ生成と描画
		SetMeshCreation();
		break;

	case RXD_UPDATED_PRTS:
		if(m_iDraw & RXD_UPDATED_PRTS){
			m_pPS->CalAnisotropicKernel();
		}
		break;

	case RXD_ANISOTROPICS:
		m_bsSimuSetting.flip(ID_SPH_ANISOTROPIC);
		if(m_bsSimuSetting.at(ID_SPH_ANISOTROPIC)){
			RXCOUT << "Anisotropic Kernel" << endl;
			m_pPS->CalAnisotropicKernel();
			RXTIMER_CLEAR;
		}
		else{
			RXCOUT << "Normal Kernel" << endl;
			m_pPS->m_bCalAnisotropic = false;
		}
		if(m_bsSimuSetting.at(ID_SPH_MESH)){
			CalMeshSPH(m_iMeshMaxN, m_fMeshThr);
		}
		break;

	default:
		break;
	}
	
	RXCOUT << "draw : " << label << " - " << GetBitArray(m_iDraw, 16) << endl;
	m_pParent->UpdateMenuState();
	redraw();
}


/*!
 * Fl_Menu_Barのコールバック関数 - Simulation
 * @param[in] widget ウィジットの親クラスオブジェクト
 * @param[in] x ユーザ定義変数
 */
void rxFlGLWindow::OnMenuSimulation_s(Fl_Widget *widget, void* x)
{
	Fl_Menu_Bar *menubar = (Fl_Menu_Bar*)widget;
	char picked[80];
	menubar->item_pathname(picked, sizeof(picked)-1);	// メニュー名

	string label = picked;
	string menu_name = "Simulation/";
	label = label.substr(menu_name.size(), string::npos);

	((rxFlGLWindow*)x)->OnMenuSimulation(1.0, label);
}
void rxFlGLWindow::OnMenuSimulation(double val, string label)
{
	if(label.find("Reset") != string::npos){
		// シーンリセット
		if(m_pPS) delete m_pPS;
		m_pPS = 0;
		InitSPH(m_Scene);
	}

	if(label.find("Wavelet Turbulence") != string::npos){	
		// ウェーブレット乱流
		m_pPS->ToggleWaveletTurb();
		m_bsSimuSetting.set(ID_SPH_PS_TURB, m_pPS->IsWaveletTurb());
		RXCOUT << "Turb : " << (m_pPS->IsWaveletTurb() ? "on" : "off") << endl;
	}
	else if(label.find("SPS Turbulence") != string::npos){
		// ウェーブレット乱流(Sub-Particle-Scale)
		m_bsSimuSetting.flip(ID_SPH_SPS_TURB);
		RXCOUT << "SPS : " << (m_bsSimuSetting[ID_SPH_SPS_TURB] ? "on" : "off") << endl;
		if(m_bsSimuSetting.at(ID_SPH_MESH) && m_iCurrentStep){
			CalMeshSPH(m_iMeshMaxN, m_fMeshThr);
		}
		m_pPS->ToggleSubParticle();
		if(m_bsSimuSetting.at(ID_SPH_SPS_TURB)){
			((RXSPH*)m_pPS)->InitSubParticles(g_fEtCri);
		}
	}
	else if(label.find("Vorticity Confinement") != string::npos){
		// Vorticity Confinement
		m_pPS->ToggleUseVorticity();
		m_bsSimuSetting.set(ID_SPH_VC, m_pPS->IsUseVorticity());
	}
	//追加
	else if( label.find("SPH only") != string::npos )
	{
		m_bMode = MODE_SPH;
	}
	else if( label.find("HeatTransfar") != string::npos )
	{
		//熱処理のみ
		cout << "Change HeatTransfar" << endl;
//		RXCOUT << "SPS : " << (m_bsSimuSetting[ID_SPH_SPS_TURB] ? "on" : "off") << endl;
		m_bsSimuSetting.set(ID_HEAT, !m_bsSimuSetting.at(ID_HEAT));
		m_bsSimuSetting.set(ID_SM, false);
		m_bsSimuSetting.set(ID_ICE, false);
		m_bMode = MODE_HEAT;
	}
	//追加
	else if( label.find("ShapeMatching") != string::npos )
	{
		//SM法のみ
		cout << "Change ShapeMatching" << endl;
		m_bsSimuSetting.set(ID_HEAT, false);
		m_bsSimuSetting.set(ID_SM, !m_bsSimuSetting.at(ID_SM));
		m_bsSimuSetting.set(ID_ICE, false);
		m_bMode = MODE_SM;
	}
	//追加
	else if( label.find("IceStructure") != string::npos )
	{
		//氷構造
		cout << "Change IceStructure" << endl;
		m_bsSimuSetting.set(ID_HEAT, false);
		m_bsSimuSetting.set(ID_SM, false);
		m_bsSimuSetting.set(ID_ICE, !m_bsSimuSetting.at(ID_ICE));
		m_bMode = MODE_ICE;
	}

#ifdef RX_USE_PBD
	else if(label.find("Artificial Pressure") != string::npos){
		// クラスタリングを防ぐための人工圧力のON/OFF (Tensile Instability)
		static_cast<RXSPH*>(m_pPS)->GetArtificialPressure() ^= 1;
	}
#endif
	else if(label.find("Particle Data Input") != string::npos){
		// パーティクルデータのインプット
		m_bsSimuSetting.flip(ID_SPH_INPUT);
	}
	else if(label.find("Particle Data Output") != string::npos){
		// パーティクルデータのアウトプット
		m_bsSimuSetting.flip(ID_SPH_OUTPUT);
		if(m_bsSimuSetting.at(ID_SPH_OUTPUT)) m_pPS->OutputSetting(m_strSphOutputName0);
	}
	else if(label.find("Mesh Saving") != string::npos){
		// メッシュ保存	
		m_bsSimuSetting.flip(ID_SPH_MESH_OUTPUT);
	}
	else if(label.find("Image Saving") != string::npos){
		// 画面の定期画像保存
		if(m_iSaveImageSpacing != -1){
			m_iSaveImageSpacing = -1;
		}
		else{
			m_iSaveImageSpacing = RX_SAVE_IMAGE_SPACING;
			RXCOUT << "save frames per " << m_iSaveImageSpacing << " steps" << endl;
		}
	}
	
	RXCOUT << "simulation : " << label << endl;
	m_pParent->UpdateMenuState();
	redraw();
}


/*!
 * Fl_Menu_Barのコールバック関数 - Particle/Color/
 * @param[in] widget ウィジットの親クラスオブジェクト
 * @param[in] x ユーザ定義変数
 */
void rxFlGLWindow::OnMenuParticle_s(Fl_Widget *widget, void* x)
{
	Fl_Menu_Bar *menubar = (Fl_Menu_Bar*)widget;
	char picked[80];
	menubar->item_pathname(picked, sizeof(picked)-1);	// メニュー名

	string label = picked;
	string menu_name = "Particle/";
	label = label.substr(menu_name.size(), string::npos);

	if(label.find("Color/") != string::npos){
		menu_name = "Color/";
		label = label.substr(menu_name.size(), string::npos);
		((rxFlGLWindow*)x)->OnMenuParticleColor(1.0, label);
	}
	else if(label.find("Draw/") != string::npos){
		menu_name = "Draw/";
		label = label.substr(menu_name.size(), string::npos);
		((rxFlGLWindow*)x)->OnMenuParticleDraw(1.0, label);
	}
	else{
		((rxFlGLWindow*)x)->OnMenuParticle(1.0, label);
	}
}
/*!
 * Particle/メニューのイベントハンドラ
 */
void rxFlGLWindow::OnMenuParticle(double val, string label)
{
	if(label.find("Anisotoropic Kernel") != string::npos){
		m_bsSimuSetting.flip(ID_SPH_ANISOTROPIC);
		if(m_bsSimuSetting.at(ID_SPH_ANISOTROPIC)){
			RXCOUT << "Anisotropic Kernel" << endl;
			m_pPS->CalAnisotropicKernel();
			RXTIMER_CLEAR;
		}
		else{
			RXCOUT << "Normal Kernel" << endl;
			m_pPS->m_bCalAnisotropic = false;
		}
		if(m_bsSimuSetting.at(ID_SPH_MESH)){
			CalMeshSPH(m_iMeshMaxN, m_fMeshThr);
		}
	}

	RXCOUT << "particle : " << label << endl;
	m_pParent->UpdateMenuState();
	redraw();
}
/*!
 * Particle/Color/メニューのイベントハンドラ
 * 粒子の色を切り替え
 */
void rxFlGLWindow::OnMenuParticleColor(double val, string label)
{	cout << __FUNCTION__ << endl;

	if(label.find("Ramp") != string::npos){
		SetParticleColorType(rxParticleSystemBase::RX_RAMP);
	}
	else if(label.find("Constant") != string::npos){
		SetParticleColorType(rxParticleSystemBase::RX_CONSTANT);
	}
	else if(label.find("Density") != string::npos){
		SetParticleColorType(rxParticleSystemBase::RX_DENSITY);
	}
	else if(label.find("Energy Spectrum") != string::npos){
		SetParticleColorType(rxParticleSystemBase::RX_ENERGY_SPECTRUM);
	}
	else if(label.find("Pressure") != string::npos){
		SetParticleColorType(rxParticleSystemBase::RX_PRESSURE);
	}
	//追加　温度
	else if(label.find("Temperature") != string::npos)
	{
		((RXSPH*)m_pPS)->DetectSurfaceParticles();
		m_pPS->SetColorType(rxParticleSystemBase::RX_TEMP);
		m_pPS->SetColorVBOFromArray(m_iceObj->GetTemps(), 1, false, 1.5f * 500.0f);
		m_iColorType = rxParticleSystemBase::RX_TEMP;
	}
	//追加　氷
	else if(label.find("Ice_Cnct") != string::npos)
	{
		((RXSPH*)m_pPS)->DetectSurfaceParticles();
		SetParticleColorType(rxParticleSystemBase::RX_ICE_CONNECT);
		m_iColorType = rxParticleSystemBase::RX_ICE_CONNECT;

		m_iShowTetraIndx++;
		
		int tetraNum = m_iceObj->GetTetraNum();

		//情報を表示する四面体番号の範囲
		if( m_iShowTetraIndx > tetraNum || m_iShowTetraIndx < 0 )
		{
			m_iShowTetraIndx = 0;
		}

		cout << "Tetrahedra :: m_iShowTetraIndx = " << m_iShowTetraIndx << "/" << tetraNum << endl;
		cout << "	Indxes :: " << endl;
		for( int i = 0; i < m_iceObj->GetTtoPIndx(m_iShowTetraIndx); i++ )
		{
			cout << "                  i = " << i << " pIndx = " << m_iceObj->GetTtoP(m_iShowTetraIndx, i) << endl;
		}
	}
	//クラスタ
	else if(label.find("Ice_Calc") != string::npos)
	{
		((RXSPH*)m_pPS)->DetectSurfaceParticles();
		SetParticleColorType(rxParticleSystemBase::RX_ICE_CALC);
		m_iColorType = rxParticleSystemBase::RX_ICE_CALC;

		m_iShowClusterIndx++;

		if(m_iShowClusterIndx > m_iClusteresNum || m_iShowClusterIndx < 0){	m_iShowClusterIndx = 0;	}

		cout << "ClusterCountUp :: m_iShowClusterIndx = " << m_iShowClusterIndx << "/" << m_iClusteresNum << endl;

		if(m_iShowClusterIndx != m_iClusteresNum && m_iShowClusterIndx >= 0)
		{
			cout << "Cluster :: Indxes :: " << endl;

			//クラスタに関する情報
			m_iceObj->GetMoveObj(m_iShowClusterIndx)->DebugClusterInfo();

			//近傍特徴クラスタ情報
			 m_iceObj->GetMoveObj(m_iShowClusterIndx)->DebugNeighborFeatureClusterInfo();
		}

		StepParticleColor();
	}

	//追加　表面氷粒子　未実装
	else if(label.find("Edge") != string::npos)
	{
		((RXSPH*)m_pPS)->DetectSurfaceParticles();
		SetParticleColorType(rxParticleSystemBase::RX_RAMP);
//		m_pPS->SetColorType(rxParticleSystemBase::RX_EDGE);
		m_iColorType = rxParticleSystemBase::RX_EDGE;
	}
	//追加　高速化用パス
	else if(label.find("FAST_PATH") != string::npos)
	{
		((RXSPH*)m_pPS)->DetectSurfaceParticles();
		SetParticleColorType(rxParticleSystemBase::RX_ICE_FAST_PATH);
		m_iColorType = rxParticleSystemBase::RX_ICE_FAST_PATH;
		StepParticleColor();
	}
	else if(label.find("Surface") != string::npos)
	{
		((RXSPH*)m_pPS)->DetectSurfaceParticles();
		SetParticleColorType(rxParticleSystemBase::RX_SURFACE);
	}
	//追加　選択的運動計算
	else if(label.find("SELECTED") != string::npos)
	{
		((RXSPH*)m_pPS)->DetectSurfaceParticles();
		SetParticleColorType(rxParticleSystemBase::RX_ICE_SELECTED);
		m_iColorType = rxParticleSystemBase::RX_ICE_SELECTED;

		int selectNum = 0;

		cout << "selectClusters:" << endl;
		for(int cIndx = 0; cIndx < m_pPS->GetNumParticles(); cIndx++)
		{
			if(m_iceObj->GetMotionCalcClsuter(cIndx) == 0){	continue;	}
			//cout << "	cIndx = " << cIndx << " , " << m_iceObj->GetMotionCalcClsuter(cIndx) << endl;
			selectNum++;
		}

		cout << "selectNum = " << selectNum << endl;

		StepParticleColor();
	}
	else if(label.find("DEFORMATION") != string::npos)
	{
		((RXSPH*)m_pPS)->DetectSurfaceParticles();
		SetParticleColorType(rxParticleSystemBase::RX_ICE_DEFORMATION);
		m_iColorType = rxParticleSystemBase::RX_ICE_DEFORMATION;
		StepParticleColor();
	}
	else if(label.find("None") != string::npos){
		SetParticleColorType(rxParticleSystemBase::RX_NONE);
	}
	
	RXCOUT << "color of particle : " << label << endl;
	m_pParent->UpdateMenuState();
	redraw();
}
/*!
 * Particle/Draw/メニューのイベントハンドラ
 */
void rxFlGLWindow::OnMenuParticleDraw(double val, string label)
{
	int idx = 0;
	while(label.find(RX_PARTICLE_DRAW[2*idx]) == string::npos) idx++;

	m_iDrawPS = idx;

#if RX_USE_UPSAMPLING
	m_pPS->ToggleUpsampling((m_iDrawPS == RXP_POINT_UPSAMPLE));
	if(m_pPS->IsUpsampling()){
#if RX_USE_UPSAMPLING_GPU
		m_pPS->CalUpsamplePointsGPU();
#else
		m_pPS->DetectSurfaceParticlesGPU();
		m_pPS->CalUpsamplePoints();
#endif
	}
#endif	// #if RX_USE_UPSAMPLING

	RXCOUT << "drawing method for particles : " << label << endl;
	m_pParent->UpdateMenuState();
	redraw();
}


/*!
 * Fl_Menu_Barのコールバック関数 - Solid
 * @param[in] widget ウィジットの親クラスオブジェクト
 * @param[in] x ユーザ定義変数
 */
void rxFlGLWindow::OnMenuSolid_s(Fl_Widget *widget, void* x)
{
	Fl_Menu_Bar *menubar = (Fl_Menu_Bar*)widget;
	char picked[80];
	menubar->item_pathname(picked, sizeof(picked)-1);	// メニュー名

	string label = picked;
	string menu_name = "Solid/";
	label = label.substr(menu_name.size(), string::npos);

	((rxFlGLWindow*)x)->OnMenuSolid(1.0, label);
}
void rxFlGLWindow::OnMenuSolid(double val, string label)
{
	int idx = 0;
	while(label.find(RXS_STR[2*idx]) == string::npos) idx++;

	int flag = (0x01 << idx);
	m_iSolidDraw ^= flag;

	switch(flag){
	case RXS_MOVE:	// 固体移動
		//m_bSolidMove = !m_bSolidMove;
		RXCOUT << "moving : " << (m_bSolidMove ? "on" : "off") << endl;
		if(m_fMoveMaxVel > RX_FEQ_EPS){
			m_bMoveSolid = !m_bMoveSolid;
			//RXCOUT << "moving : " << (m_bMoveSolid ? "on" : "off") << endl;
		}
		break;

	default:
		break;
	}
	
	RXCOUT << "solid : " << label << endl;
	m_pParent->UpdateMenuState();
	redraw();
}


/*!
 * Fl_Menu_Barのコールバック関数 - Mesh
 * @param[in] widget ウィジットの親クラスオブジェクト
 * @param[in] x ユーザ定義変数
 */
void rxFlGLWindow::OnMenuTriangulation_s(Fl_Widget *widget, void* x)
{
	Fl_Menu_Bar *menubar = (Fl_Menu_Bar*)widget;
	char picked[80];
	menubar->item_pathname(picked, sizeof(picked)-1);	// メニュー名

	string label = picked;
	string menu_name = "Mesh/";
	label = label.substr(menu_name.size(), string::npos);

	((rxFlGLWindow*)x)->OnMenuTriangulation(1.0, label);
}
void rxFlGLWindow::OnMenuTriangulation(double val, string label)
{
	int idx = 0;
	while(label.find(RX_TRIANGULATION_METHOD[2*idx]) == string::npos) idx++;

	m_iTriangulationMethod = idx;
	if(m_bsSimuSetting.at(ID_SPH_MESH))	CalMeshSPH(m_iMeshMaxN, m_fMeshThr);

	RXCOUT << "triangulation method : " << label << endl;
	m_pParent->UpdateMenuState();
	redraw();
}


/*!
 * Fl_Menu_Barのコールバック関数 - Scene
 * @param[in] widget ウィジットの親クラスオブジェクト
 * @param[in] x ユーザ定義変数
 */
void rxFlGLWindow::OnMenuScene_s(Fl_Widget *widget, void* x)
{
	Fl_Menu_Bar *menubar = (Fl_Menu_Bar*)widget;
	char picked[80];
	menubar->item_pathname(picked, sizeof(picked)-1);	// メニュー名

	string label = picked;
	string menu_name = "Scene/";
	label = label.substr(menu_name.size(), string::npos);

	((rxFlGLWindow*)x)->OnMenuScene(1.0, label);
}
void rxFlGLWindow::OnMenuScene(double val, string label)
{
	if(m_Scene.SetCurrentSceneFromTitle(label)){
		m_iCurrentSceneIdx = m_Scene.GetCurrentSceneIdx();
		InitSPH(m_Scene);
	}
	
	m_pParent->UpdateMenuState();
	redraw();
}



//-----------------------------------------------------------------------------
// MARK:メッシュ作成
//-----------------------------------------------------------------------------
/*!
 * 三角形メッシュ生成
 * @param[in] nmax メッシュ化グリッド解像度(最大)
 * @param[in] thr メッシュ化閾値
 */
bool rxFlGLWindow::CalMeshSPH(int nmax, double thr)
{
	switch(m_iTriangulationMethod){
	case RXM_MC_CPU:
		return calMeshSPH_CPU(nmax, thr);

	case RXM_MC_GPU:
		return calMeshSPH_GPU(nmax, thr);

	case RXM_SSM_CPU:
		return calMeshSPH_SSM(nmax, thr);

	case RXM_SSM_GPU:
		return calMeshSPH_SSM_GPU(nmax, thr);

	default:
		return calMeshSPH_GPU(nmax, thr);
	};

	return false;
}

/*!
 * メッシュ情報の初期化
 */
bool rxFlGLWindow::ResetMesh(void)
{
	// ポリゴン初期化
	if(!m_Poly.vertices.empty()){
		m_Poly.vertices.clear();
		m_Poly.normals.clear();
		m_Poly.faces.clear();
		m_Poly.materials.clear();
	}
	m_iNumVrts = 0;
	m_iNumTris = 0;

	// 頂点VBO
	if(!m_uVrtVBO) glDeleteBuffers(1, &m_uVrtVBO);
	m_uVrtVBO = 0;

	// 法線VBO
	if(!m_uNrmVBO) glDeleteBuffers(1, &m_uNrmVBO);
	m_uNrmVBO = 0;

	// メッシュVBO
	if(m_uTriVBO) glDeleteBuffers(1, &m_uTriVBO);
	m_uTriVBO = 0;

	if(m_pMCMeshCPU) delete m_pMCMeshCPU;
	if(m_pMCMeshGPU) delete m_pMCMeshGPU;
	if(m_pSSMCPU) delete m_pSSMCPU;
	if(m_pSSMGPU) delete m_pSSMGPU;


	m_pMCMeshCPU = 0;
	m_pMCMeshGPU = 0;
	m_pSSMCPU = 0;
	m_pSSMGPU = 0;

	return true;
}


/*!
 * メッシュ生成のセット
 */
void rxFlGLWindow::SetMeshCreation(void)
{
	m_bsSimuSetting.set(ID_SPH_MESH, ((m_iDraw & RXD_MESH) ? true : false));
	if(m_bsSimuSetting.at(ID_SPH_MESH) && m_iCurrentStep){
		CalMeshSPH(m_iMeshMaxN, m_fMeshThr);
	}
}

/*!
 * 三角形メッシュの生成(MC法,CPU)
 * @param[in] nmax メッシュ化グリッド解像度(最大)
 * @param[in] thr メッシュ化閾値
 */
bool rxFlGLWindow::calMeshSPH_CPU(int nmax, double thr)
{
	m_pPS->SetParticlesToCell();

	Vec3 minp = m_pPS->GetMin();
	Vec3 maxp = m_pPS->GetMax();

	double h;
	int n[3];
	CalMeshDiv(minp, maxp, nmax, h, n, 0.05);
	for(int i = 0; i < 3; ++i) m_iMeshN[i] = n[i];
	//cout << "mc : " << n[0] << " x " << n[1] << " x " << n[2] << endl;

	// ポリゴン初期化
	if(!m_Poly.vertices.empty()){
		m_Poly.vertices.clear();
		m_Poly.normals.clear();
		m_Poly.faces.clear();
		m_Poly.materials.clear();
	}

	if(!m_pMCMeshCPU){
		m_pMCMeshCPU = new rxMCMeshCPU;
	}

	// HACK:calMeshSPH_CPU
	//m_pMCMeshCPU->CreateMesh(GetImplicitSPH, minp, h, n, thr, m_Poly.vertices, m_Poly.normals, m_Poly.faces);

	m_iNumVrts = (int)m_Poly.vertices.size();
	m_iNumTris = (int)m_Poly.faces.size();
	
	if(m_Poly.normals.empty()){
		CalVertexNormals(m_Poly);
	}
	m_Poly.materials[m_matPoly.name] = m_matPoly;

	m_iDimVBO = 3;

	AssignArrayBuffers(m_iNumVrts, 3, m_uVrtVBO, m_uNrmVBO, m_uTriVBO);
	SetFBOFromArray(m_uVrtVBO, m_uNrmVBO, m_uTriVBO, m_Poly.vertices, m_Poly.normals, m_Poly.faces);

	return true;
}

/*!
 * 三角形メッシュの生成(MC法,GPU)
 * @param[in] nmax メッシュ化グリッド解像度(最大)
 * @param[in] thr メッシュ化閾値
 */
bool rxFlGLWindow::calMeshSPH_GPU(int nmax, double thr)
{
	if(m_pPS == NULL) return false;

	Vec3 minp = m_vMeshBoundaryCen-m_vMeshBoundaryExt;
	Vec3 maxp = m_vMeshBoundaryCen+m_vMeshBoundaryExt;

	double h;	//セル幅
	int n[3];	//各軸のセル数
	CalMeshDiv(minp, maxp, nmax, h, n, 0.05);
	for(int i = 0; i < 3; ++i) m_iMeshN[i] = n[i];
	//cout << "mc : " << n[0] << " x " << n[1] << " x " << n[2] << endl;

	m_iDimVBO = 4;

	if(!m_pMCMeshGPU){
		m_pMCMeshGPU = new rxMCMeshGPU;
		m_pMCMeshGPU->Set(minp, Vec3(h), n, m_iVertexStore);

		// メッシュデータ格納用のVBOの確保
		AssignArrayBuffers(m_pMCMeshGPU->GetMaxVrts(), 4, m_uVrtVBO, m_uNrmVBO, m_uTriVBO);

		////追加：氷のフラグ用データをGPUに確保
		//CuAllocateArray((void**)&m_fIceFlag, m_pPS->GetMaxParticles()*sizeof(RXREAL));
		//CuSetArrayValue((void*)m_fIceFlag, 0, m_pPS->GetMaxParticles()*sizeof(RXREAL));

		////追加：固体用データの確保
		//AssignArrayBuffers(m_pMCMeshGPU->GetMaxVrts(), 4, m_uVrtVBO_solid, m_uNrmVBO_solid, m_uTriVBO_solid);
	}
	
	if(!m_Poly.vertices.empty()){
		m_Poly.vertices.clear();
		m_Poly.normals.clear();
		m_Poly.faces.clear();
		m_Poly.materials.clear();
	}
	
	//cout << "check start" << endl;
	// サンプリングボリュームの計算(陰関数値を格納したグリッドデータ)
	m_pPS->CalImplicitFieldDevice(n, minp, Vec3(h, h, h), m_pMCMeshGPU->GetSampleVolumeDevice());
	//cout << "check end" << endl;
	unsigned int nvrts = 0;
	unsigned int ntris = 0;
	m_pMCMeshGPU->CreateMeshV(minp, h, n, thr, nvrts, ntris);
	
	m_iNumVrts = nvrts;
	m_iNumTris = ntris;
	if(m_iCurrentStep%RX_DISPLAY_STEP == 0 && m_iCurrentStep != 0){
		RXCOUT << nvrts << " verts and " << ntris << " tri were created." << endl;
	}

	// FBOにデータをコピー
	m_pMCMeshGPU->SetDataToFBO(m_uVrtVBO, m_uNrmVBO, m_uTriVBO);
	
/*
	//追加：氷用処理
	float* fIceCheck = new float[m_pPS->GetNumParticles()];
	for(int i=0; i<m_pPS->GetNumParticles(); i++)
	{//	cout << "i = " << i << " m_ice->GetCtoPNum_Connect = " << m_ice->GetCtoPNum_Connect(i) << endl;
		if(m_ice->GetCtoPNum_Connect(i) != 0)
		{
			fIceCheck[i] = 1.0f;
		}
		else
		{
			fIceCheck[i] = -1.0f;
		}
	}

	CuCopyArrayToDevice(m_fIceFlag, fIceCheck, 0, m_pPS->GetNumParticles()*sizeof(RXREAL));		//氷のフラグ情報の更新
	m_pPS->CalImplicitFieldDeviceSolid(n, minp, Vec3(h, h, h), m_pMCMeshGPU->GetSampleVolumeDevice(), m_fIceFlag);
	delete[] fIceCheck;

	nvrts = 0;
	ntris = 0;
	//ここのthrで閾値を変えられる
	m_pMCMeshGPU->CreateMeshV(minp, h, n, thr*1.2, nvrts, ntris);

	m_iNumVrts_solid = nvrts;
	m_iNumTris_solid = ntris;
	if(m_iCurrentStep%RX_DISPLAY_STEP == 0 && m_iCurrentStep != 0){
		RXCOUT << "iceMesh " << nvrts << " verts and " << ntris << " tri were created (solid)." << endl;
	}

	//// FBOにデータをコピー
	m_pMCMeshGPU->SetDataToFBO(m_uVrtVBO_solid, m_uNrmVBO_solid, m_uTriVBO_solid);
	//終了：氷用処理 DrawSolidSurface()で描画
*/
	// ファイル保存のためにメッシュデータをホスト側配列にコピー
	if(m_bsSimuSetting.at(ID_SPH_MESH_OUTPUT)){
		m_pMCMeshGPU->SetDataToArray(m_Poly.vertices, m_Poly.normals, m_Poly.faces);
	}

	return true;
}


/*!
 * 三角形メッシュの生成(Screen Space Mesh, CPU)
 * @param[in] nmax メッシュ化グリッド解像度(最大)
 * @param[in] thr メッシュ化閾値
 */
bool rxFlGLWindow::calMeshSPH_SSM(int nmax, double thr)
{
	// MRK:calMeshSPH_SSM
	if(!m_pPS) return false;
	if(!m_pPS->GetNumParticles()) return false;

	int pnum = m_pPS->GetNumParticles();
	m_fPrtRad = m_pPS->GetParticleRadius()*1.5;
	m_fZmax = 10.0*m_fPrtRad;

	// SSM初期化
	if(!m_pSSMCPU){
		m_pSSMCPU = new rxSSMeshCPU(m_fZmax, m_fSpacing, m_fPrtRad, m_iNfilter, m_iNiters);
	}
	else{
		m_pSSMCPU->SetSpacing(m_fSpacing);
		m_pSSMCPU->SetRadius(m_fPrtRad);
		m_pSSMCPU->SetZMax(m_fZmax);
		m_pSSMCPU->SetFilterRadius(m_iNfilter);
		m_pSSMCPU->SetSmoothIter(m_iNiters);
	}

	// ポリゴン初期化
	if(!m_Poly.vertices.empty()){
		m_Poly.vertices.clear();
		m_Poly.normals.clear();
		m_Poly.faces.clear();
		m_Poly.materials.clear();
	}

	// パーティクルデータ生成
	RXREAL *data = 0;
	data = m_pPS->GetArrayVBO(rxParticleSystemBase::RX_POSITION);

	vector<Vec3> prts;
	prts.resize(pnum);
	int k = 0;
	for(int i = 0; i < pnum; ++i){
		prts[i] = Vec3(data[k], data[k+1], data[k+2]);
		k += 4;
	}

	// スクリーンスペースメッシュ生成
	vector< vector<int> > tris;
	m_pSSMCPU->CreateMesh(m_fProjectionMatrix, m_fModelviewMatrix, m_iWinW, m_iWinH, 
					   prts, pnum, m_Poly.vertices, m_Poly.normals, tris, 
					   m_iDepthFiltering, m_iSSDebugOutput);

	// ポリゴン生成
	int n = (int)tris.size();
	m_Poly.faces.resize(n);
	for(int i = 0; i < n; ++i){
		m_Poly.faces[i].vert_idx.resize(3);

		for(int j = 0; j < 3; ++j){
			m_Poly.faces[i][j] = tris[i][j];
		}
	}
	if(m_Poly.normals.empty()){
		CalVertexNormals(m_Poly);
	}

	// ポリゴン材質
	m_Poly.materials[m_matPoly.name] = m_matPoly;

	// VBOへの移動
	int nv = (int)m_Poly.vertices.size();
	int nm = (int)tris.size();
	int nn = nv;

	m_iNumVrts = nv;
	m_iNumTris = nm;

	m_iDimVBO = 3;
	AssignArrayBuffers(nv, 3, m_uVrtVBO, m_uNrmVBO, m_uTriVBO);
	SetFBOFromArray(m_uVrtVBO, m_uNrmVBO, m_uTriVBO, m_Poly.vertices, m_Poly.normals, m_Poly.faces);

	if(m_iCurrentStep%RX_DISPLAY_STEP == 0 && m_iCurrentStep != 0){
		RXCOUT << m_iNumVrts << " verts and " << m_iNumTris << " tri were created." << endl;
	}

	//int nx, ny;
	//RXREAL *dm;

	//m_pSSMCPU->GetDepthMapSize(nx, ny);
	//dm = m_pSSMCPU->GetDepthMap();

	//// デプス値テクスチャ生成
	//g_texDepth.SetSize(nx-1, ny-1, 4);
	//for(int j = 0; j < ny-1; ++j){
	//	for(int i = 0; i < nx-1; ++i){
	//		//rxSSDepth depth = 0.25*(dm[i+j*nx]+dm[(i+1)+j*nx]+dm[i+(j+1)*nx]+dm[(i+1)+(j+1)*nx]);
	//		//unsigned int ud = RX_DEPTH2COLOR(depth);
	//		unsigned int ud = RX_DEPTH2COLOR(dm[i+j*nx]);

	//		g_texDepth.SetColor(i, j, ud, ud, ud, 255);
	//	}
	//}
	//BindTexture(g_texDepth);


	// デプス値画像出力
	//SaveTexture("depth.png", g_texDepth);

	return true;
}

/*!
 * 三角形メッシュの生成(Screen Space Mesh, GPU)
 * @param[in] nmax メッシュ化グリッド解像度(最大)
 * @param[in] thr メッシュ化閾値
 */
bool rxFlGLWindow::calMeshSPH_SSM_GPU(int nmax, double thr)
{
	// MRK:calMeshSPH_SSM_GPU
	if(!m_pPS) return false;
	if(!m_pPS->GetNumParticles()) return false;

	int pnum = m_pPS->GetNumParticles();
	m_fPrtRad = m_pPS->GetParticleRadius();
	m_fZmax = 10.0*m_fPrtRad;

	m_iDimVBO = 3;

	// SSM初期化
	if(!m_pSSMGPU){
		m_pSSMGPU = new rxSSMeshGPU(m_fZmax, m_fSpacing, m_fPrtRad, m_iNfilter, m_iNiters);

		//if(m_uVrtVBO) glDeleteBuffers(1, &m_uVrtVBO);
		//if(m_uTriVBO) glDeleteBuffers(1, &m_uTriVBO);
		//if(m_uNrmVBO) glDeleteBuffers(1, &m_uNrmVBO);

		//int max_vrts = (m_pSSMGPU->GetNx()+1)*(m_pSSMGPU->GetNy()+1)*1.5;
		//AssignArrayBuffers(max_vrts, 3, m_uVrtVBO, m_uNrmVBO, m_uTriVBO);
	}
	else{
		m_pSSMGPU->SetSpacing(m_fSpacing);
		m_pSSMGPU->SetRadius(m_fPrtRad);
		m_pSSMGPU->SetZMax(m_fZmax);
		m_pSSMGPU->SetFilterRadius(m_iNfilter);
		m_pSSMGPU->SetSmoothIter(m_iNiters);
	}

	// ポリゴン初期化
	if(!m_Poly.vertices.empty()){
		m_Poly.vertices.clear();
		m_Poly.normals.clear();
		m_Poly.faces.clear();
		m_Poly.materials.clear();
	}

	// パーティクルデータ生成
	RXREAL *dprts = 0;
	RXREAL *drads = 0;

	if(m_pPS->IsSubParticle()){
		// 有効なサブパーティクルリストの作成
		((RXSPH*)m_pPS)->CalRadiusAndRatio();

		// パーティクルデータの取得
		dprts = ((RXSPH*)m_pPS)->GetSubParticlePosDev();		// 座標
		drads = ((RXSPH*)m_pPS)->GetSubParticleRadDev();		// 半径
		pnum  = ((RXSPH*)m_pPS)->GetNumValidSubParticles();	// パーティクル数
	}
	else{
		dprts = m_pPS->GetParticleDevice();
	}

	// スクリーンスペースメッシュ生成
	m_pSSMGPU->CreateMeshVBO(m_fProjectionMatrix, m_fModelviewMatrix, m_iWinW, m_iWinH, dprts, drads, pnum, 4, 
								m_uVrtVBO, m_uNrmVBO, m_uTriVBO, m_iDepthFiltering, m_iSSDebugOutput);

	if(!m_pPS->IsSubParticle()){
		m_pPS->UnmapParticle();
	}

	if(m_iCurrentStep%RX_DISPLAY_STEP == 0 && m_iCurrentStep != 0){
		RXCOUT << m_iNumVrts << " verts and " << m_iNumTris << " tri were created." << endl;
	}

	m_iNumVrts = m_pSSMGPU->GetVertexNum();
	m_iNumTris = m_pSSMGPU->GetMeshNum();

	if(m_bsSimuSetting.at(ID_SPH_MESH_OUTPUT)){
		SetArrayFromFBO(m_uVrtVBO, m_uNrmVBO, m_uTriVBO, m_Poly.vertices, m_Poly.normals, m_Poly.faces, m_iNumVrts, m_iNumTris);
	}


	return true;
}

//-----------------------------------------------------------------------------
// SPH関数
//-----------------------------------------------------------------------------
/*!
 * SPHの初期化
 * @param[in] fn_scene シーン記述ファイル名
 */
void rxFlGLWindow::InitSPH(rxSPHConfig &sph_scene)
{
	// SPHクラスの初期化とシーンクラスへの設定
	if(m_pPS) delete m_pPS;
	m_pPS = new RXSPH(true);
	sph_scene.Clear();
	sph_scene.SetPS(m_pPS);
		
	// シーン全体情報の読み込み
	if(m_bsSimuSetting.at(ID_SPH_INPUT) && ExistFile(m_strSphInputName0)){
		sph_scene.LoadSpaceFromFile(m_strSphInputName0);
	}
	else{
		if(sph_scene.LoadSpaceFromFile()){
			rxSPHEnviroment sph_env = sph_scene.GetSphEnv();
			((RXSPH*)m_pPS)->Initialize(sph_env);

			m_fDt = sph_env.dt;
			m_iVertexStore = sph_env.mesh_vertex_store;
			m_iMeshMaxN = sph_env.mesh_max_n;
			g_fEtCri = sph_env.et_cri;
			m_bsSimuSetting.set(ID_SPH_INLET, (sph_env.use_inlet ? true : false));
		}
	}
	m_iCurrentStep = 0;

	// メッシュ生成境界の設定
	m_vMeshBoundaryCen = sph_scene.GetSphEnv().mesh_boundary_cen;
	m_vMeshBoundaryExt = sph_scene.GetSphEnv().mesh_boundary_ext;

	// 乱流に関するフラグの設定
	m_pPS->ToggleWaveletTurb(m_bsSimuSetting.at(ID_SPH_PS_TURB));
	m_pPS->ToggleUseVorticity(m_bsSimuSetting.at(ID_SPH_VC));

	// パーティクル描画のためのカラータイプの設定
	m_pPS->SetColorType(m_iColorType);
	m_pPS->SetColorVBO();

	// シーンの個別設定(パーティクル，固体)の読み込み
	if(m_bsSimuSetting.at(ID_SPH_INPUT)){
		// パーティクル情報をファイルから読み込み
		m_pPS->InputParticles(CreateFileName(m_strSphInputHeader, "dat", 0, 5));
		CalMeshSPH(m_iMeshMaxN, m_fMeshThr);
	}
	else{
		m_pPS->Reset(rxParticleSystemBase::RX_CONFIG_NONE);
		sph_scene.LoadSceneFromFile();
	}

	// サブパーティクル設定
	((RXSPH*)m_pPS)->SetEtcri(g_fEtCri);
	((RXSPH*)m_pPS)->InitSubParticles(g_fEtCri);

	// 異方性カーネル設定
	if(m_bsSimuSetting.at(ID_SPH_ANISOTROPIC)){
		m_pPS->CalAnisotropicKernel();
	}

	// 表面メッシュのリセット
	ResetMesh();


	// 時間計測用変数の初期化
	g_fTotalTime = 0.0f;
	g_fAvgTime = 0.0f;
	g_iTimeCount = 0;
}

/*!
 * シーンに水滴を追加
 */
void rxFlGLWindow::AddSphere(void)
{
	int ballr = 10;

	RXREAL pr = m_pPS->GetParticleRadius();
	RXREAL tr = pr+(pr*2.0f)*ballr;
	RXREAL pos[4], vel[4];
	pos[0] = -1.0+tr+RX_FRAND()*(2.0-tr*2.0);
	pos[1] =  1.0-tr;
	pos[2] = -1.0+tr+RX_FRAND()*(2.0-tr*2.0);
	pos[3] =  0.0;
	vel[0] = vel[1] = vel[2] = vel[3] = 0.0;
	m_pPS->AddSphere(0, pos, vel, ballr, pr*2.0f);
}


/*!
 * SPHのタイムステップを進める
 * @param[in] dt タイムステップ幅
 */
void rxFlGLWindow::StepPS(double dt)
{
	if(!m_bPause)
	{
		float damping = 1.0f;
		int iterations = 1;

		// シミュレーションパラメータ
		m_pPS->SetIterations(iterations);
		m_pPS->SetDamping(damping);
		m_pPS->SetGravity((RXREAL)(-m_fGravity));

		// シミュレーションステップを進める
		m_pPS->Update((RXREAL)dt, m_iCurrentStep);
		//m_pPS->UpdateWithoutPosAndVel((RXREAL)dt, m_iCurrentStep);

		//追加：粒子数をチェックして，各情報を更新
		UpdateInfo();
	}
}

/*!
 * 追加された粒子に対する処理
 * TODO::IceObject.cppを追加しても正常に動くように修正
 */
void rxFlGLWindow::UpdateInfo()
{
	int beforeSize = m_fIntrps.size();
	int nowSize = m_pPS->GetNumParticles();

	if( beforeSize == 0 || beforeSize >= nowSize ) return;		//粒子数に変化がないなら終了

	//cout << __FUNCTION__ << " AddInfo beforeSize = " << beforeSize << " nowSize = " << nowSize << endl;

	//熱処理
	//m_ht->AddParticle( nowSize );

	//線形補間
	//for( int i = beforeSize; i < nowSize; i++ )
	//{
	//	m_fIntrps.push_back( 1.0f-(m_ht->getTemps()[i] / m_ht->getTempMax()) );
	//}

	//融解の場合は，一切凝固が行われないとして追加しない
	return;

	//RXREAL *p = m_pPS->GetArrayVBO(rxParticleSystemBase::RX_POSITION);
	//RXREAL *v = m_pPS->GetArrayVBO(rxParticleSystemBase::RX_VELOCITY);

	////クラスタ
	//rxSPHEnviroment sph_env = m_Scene.GetSphEnv();				// Shape Matchingの設定　パラメータ読み込み
	//for(int i = beforeSize; i < nowSize; i++)
	//{
	//	m_sm_cluster.push_back(new Ice_SM(i));
	//	m_sm_cluster[i]->SetSimulationSpace(-sph_env.boundary_ext, sph_env.boundary_ext);
	//	m_sm_cluster[i]->SetTimeStep(sph_env.smTimeStep);
	//	m_sm_cluster[i]->SetCollisionFunc(0);
	//	m_sm_cluster[i]->SetStiffness(1.0, 0.0);
	//	m_iClusteresNum++;
	//}

	////固体構造
	//m_ice->SetParticleNum(nowSize);
	//m_ice->SetClusterNum(nowSize);

	////デバイスメモリの更新

}

//-----------------------------------------------------------------------------
// ShapeMatching関数
//-----------------------------------------------------------------------------
/*!
 * 衝突処理関数
 * @param[in] p 現在の座標
 * @param[out] np 衝突後の座標
 * @param[in] v 速度
 * @param[in] obj オブジェクト番号
 */
void rxFlGLWindow::Collision(Vec3 &p, Vec3 &np, Vec3 &v, int obj)
{
	//// 球体との衝突判定
	//Vec3 rpos = p-g_v3Cen;
	//double d = norm(rpos)-g_fRad;
	//if(d < 0.0){
	//	Vec3 n = Unit(rpos);
	//	np = g_v3Cen+n*g_fRad;
	//}

	//// 変形メッシュ同士の衝突
	//uint grid_hash = m_pNNGrid->CalGridHash(np);
	//vector<int> polys_in_cell;
	//if(m_pNNGrid->GetPolygonsInCell(grid_hash, polys_in_cell)){
	//	for(int j = 0; j < (int)polys_in_cell.size(); ++j){
	//		int pidx = polys_in_cell[j];

	//		int vidx[3];
	//		vidx[0] = g_vTris[3*pidx+0];
	//		vidx[1] = g_vTris[3*pidx+1];
	//		vidx[2] = g_vTris[3*pidx+2];

	//		if(obj >= 0 && (vidx[0] >= g_vObj[obj].vstart && vidx[0] <= g_vObj[obj].vend)){
	//			continue;
	//		}

	//		Vec3 vrts[3];
	//		vrts[0] = Vec3(g_vVrts[3*vidx[0]], g_vVrts[3*vidx[0]+1], g_vVrts[3*vidx[0]+2]);
	//		vrts[1] = Vec3(g_vVrts[3*vidx[1]], g_vVrts[3*vidx[1]+1], g_vVrts[3*vidx[1]+2]);
	//		vrts[2] = Vec3(g_vVrts[3*vidx[2]], g_vVrts[3*vidx[2]+1], g_vVrts[3*vidx[2]+2]);

	//		Vec3 cp, n;
	//		if(IntersectSegmentTriangle(p, np, vrts[0], vrts[1], vrts[2], cp, n, 0.0) == 1){
	//			double d = length(np-cp);
	//			n = Unit(n);

	//			RXREAL res = 1.0;
	//			res = (res > 0) ? (res*fabs(d)/(g_fDt*length(v))) : 0.0f;
	//			Vec3 vr = -(1.0+res)*n*dot(n, v);

	//			double l = norm(np-p);
	//			np = cp+vr*g_fDt;//*d/l);
	//			v += vr;
	//		}
	//	}
	//}
}


//-------------------------------------------------------------------------------------------------------------
//　物体情報
//-------------------------------------------------------------------------------------------------------------
/*!
 * 物体情報の初期化
 */
void rxFlGLWindow::InitIceObj(void)
{	cout << __FUNCTION__ << endl;

#ifdef _OPENMP
    cout << "OpenMP : On, threads =" << omp_get_max_threads() << endl;	//OpenMPによる簡易並列処理
#endif

	////シミュレーションテストのための準備
	//m_stSimuFiles.push_back("test1");
	//m_stSimuFiles.push_back("test2");
	//m_stSimuFiles.push_back("test3");

	//m_uSimulationNum = m_stSimuFiles.size();

	//パラメータ設定
	//TODO::ある粒子は液体，ある粒子は固体，とファイルから指定する　球は固体，立方体は液体，のように
	m_iIcePrtNum = m_pPS->GetNumParticles();
	m_iIceTtrNum = m_iIcePrtNum * 5.4 + 4000;					//四面体数　Excelでてきとうに線形近似して得られた式
	m_iLayer = m_Scene.GetSphEnv().layer;
	m_iIceItr = m_Scene.GetSphEnv().smItr;

	m_iMaxClusterNum = m_pPS->GetNumParticles();

	//ファイルから読み込んだパラメータ
	rxSPHEnviroment sph_env = m_Scene.GetSphEnv();

	//ポインタの初期化
	m_ice = 0;
	m_iceObj = 0;

	RXREAL *hp = m_pPS->GetArrayVBO(rxParticleSystemBase::RX_POSITION);
	RXREAL *hv = m_pPS->GetArrayVBO(rxParticleSystemBase::RX_VELOCITY);
	RXREAL *dp = ((RXSPH*)m_pPS)->GetDevicePointer_Pos();
	RXREAL *dv = ((RXSPH*)m_pPS)->GetDevicePointer_Vel();

	//ソリッドモデル
	m_iceObj = new IceObject(
		m_iIcePrtNum+1000,
		m_iIcePrtNum+1000,
		m_iIceTtrNum,
		m_iIcePrtNum,
		hp, hv, dp, dv,
		m_iLayer,
		m_iMaxClusterNum
	);

	ApplyNowData();

	//表面粒子
	((RXSPH*)m_pPS)->DetectSurfaceParticles();								//表面粒子検出
	int* surfaceParticles = (int *)( ((RXSPH*)m_pPS)->GetArraySurf() );		//表面粒子

	//近傍粒子
	float radius = ((RXSPH*)m_pPS)->GetEffectiveRadius();
	((RXSPH*)m_pPS)->SetEffectiveRadius(radius * 1.2f);
	StepPS(m_fDt*0.0);														//一度タイムステップを勧めないと，近傍粒子が取得されないみたい
	((RXSPH*)m_pPS)->searchNeighbors();
	vector<vector<rxNeigh>>& neights = ((RXSPH*)m_pPS)->GetNeights();
	((RXSPH*)m_pPS)->SetEffectiveRadius(radius);

	//sm法の初期化
	m_iceObj->InitCluster(
		sph_env.boundary_ext,
		-sph_env.boundary_ext,
		sph_env.smTimeStep,
		m_iIceItr,
		neights,
		surfaceParticles
	);

	m_iceObj->InitStrct();										//粒子とクラスタの関係情報を初期化

	//熱処理初期化
	m_iceObj->InitHeatTransfer(
		((RXSPH*)m_pPS)->GetEffectiveRadius(),
		sph_env.htTimeStep,
		sph_env.tempMax,
		sph_env.tempMin,
		sph_env.latentHeat,
		sph_env.cffcntHt,
		sph_env.cffcntTd
	);

	//運動計算部分の初期化
	m_iceObj->InitMoveMethod();

	m_iClusteresNum = m_iceObj->GetClusterNum();
	m_iShowClusterIndx = m_iClusteresNum;

	m_iceObj->SaveInitPos();	//初期位置の保存

//デバッグ
	//DumpParticleData();										//初期粒子情報をダンプ
}

/*
 *	あるタイムステップを迎えると行われるイベント
 *  データ取得や動画作成に使う
 */
void rxFlGLWindow::StepTimeEvent(void)
{
	////粒子が中心の列から融解
	//if(200 <= g_iTimeCount && g_iTimeCount <= 300)
	//{		
	//	//特定領域の粒子を融解　真ん中から真っ二つ
	//	//上の列の粒子を先に温める
	//	/*for(int i = CUBE_X*CUBE_Y*(CUBE_Z/2+1); i < CUBE_X*CUBE_Y*(CUBE_Z/2+2); i++)*/
	//	//for(int i = SIDE*SIDE*(SIDE/2); i < SIDE*SIDE*(SIDE/2+1); i++)	//一列
	//	for(int i = SIDE*SIDE*(SIDE/2 -1); i < SIDE*SIDE*(SIDE/2+2); i++)	//三列
	//	{
	//		int pIndx = i;
	//		WarmParticle(pIndx, 20.0f, 20.0f);
	//	}
	//}

	////特定の番号の粒子が融解
	//if(g_iTimeCount == 50)
	//{
	//	for(int i = 0; i < 1; i++)
	//	{
	//		int pIndx = i;
	//		m_iceObj->MeltParticle(pIndx);
	//	}
	//}

	////粒子が中心の列から融解
	//if(200 <= g_iTimeCount && g_iTimeCount <= 300)
	//{		
	//	//特定領域の粒子を融解　真ん中から真っ二つ
	//	int SIDE = pow(m_iIcePrtNum, 1.0/3.0) + 0.5;	//立方体の１辺の頂点数
	//	int rine = 1;

	//	int first = SIDE*SIDE*(SIDE/2 - (rine-1) );
	//	int end = SIDE*SIDE*(SIDE/2 + rine);

	//	for(int i = first; i < end; i++)	//rine列
	//	{
	//		m_iceObj->MeltParticle(i);
	//	}
	//}

	////あるフレーム以上になるとデータを作成して終了
	//if(250 < g_iTimeCount)
	//{
	//	m_iceObj->DebugMakeGraph();
	//	exit(1);
	//}

	////シミュレーションが終わると，次のファイル名を渡してシミュレーションさせる
	//if(300 < g_iTimeCount)
	//{
	//	ResetSimulation();
	//	g_iTimeCount = 0;
	//	m_iCurrentStep = 0;
	//	m_uSimulationNum--;
	//	if(m_uSimulationNum <= 0){	exit(1);	} //全てのシミュレーションが終了

	//	//次のシミュレーション開始
	//	string nextSimulation = m_stSimuFiles[m_uSimulationNum];
	//	m_iceObj->TestSimulationFromFile(nextSimulation);
	//}
}

/*!
 * 四面体に含まれている粒子数のカウント，粒子が所属する四面体数のカウント
 * @param[in] tNum 四面体番号
 */
void rxFlGLWindow::CountTetraHedra(int tIndx, vector<int>& pList)
{
	for(unsigned i = 0; i < pList.size(); i++)
	{
		int pIndx = pList[i];		

		m_iceObj->CountPtoT(pIndx);
		m_iceObj->CountTtoP(tIndx);

		//Indxの更新
		if(m_iceObj->GetPtoTNum(pIndx) >= m_iceObj->GetPtoTIndx(pIndx))
		{
			m_iceObj->SetPtoTIndx(pIndx);
		}
	}
	
	//Indxの更新
	if(m_iceObj->GetTtoPNum(tIndx) >= m_iceObj->GetTtoPIndx(tIndx))
	{
		m_iceObj->SetTtoPIndx(tIndx);
	}
}

/*!
 * 四面体に含まれている情報の登録　初期化専用関数
 * @param[in] tNum 四面体番号
 * @param[in] PtoTNum 四面体が属する粒子数
 */
void rxFlGLWindow::MakeTetraInfo(int tIndx, int* PtoTNum)
{
	IceTetrahedra &tetra = IceTetrahedra::GetInstance();

	//粒子が属している四面体の番号を登録するための準備
	//pCountListには，tIndx番目の四面体に含まれる各粒子が，それぞれいくつの四面体に属するかを求めて保存する
	int* pCountList = new int[tetra.GetTetraList(tIndx).size()];

	for(unsigned j = 0; j < tetra.GetTetraList(tIndx).size(); j++)
	{
		int pIndx = tetra.GetTetraList(tIndx)[j];
		pCountList[j] = m_iceObj->GetPtoTNum(pIndx)-PtoTNum[pIndx];
		PtoTNum[pIndx]--;
	}

	//粒子と四面体の情報登録
	vector<int>& pIndxList = tetra.GetTetraList(tIndx);

	for(int i = 0; i < m_iceObj->GetTtoPNum(tIndx); i++)
	{
		m_iceObj->SetPtoT(pIndxList[i], pCountList[i], tIndx, i);
	}

	m_iceObj->SetTtoP(tIndx, pIndxList);

	delete[] pCountList;
}

/*!
 * 四面体に含まれている情報の登録
 * @param[in] tNum 四面体番号
 */
void rxFlGLWindow::MakeTetraInfo(int tIndx, vector<int> pList)
{//	cout << __FUNCTION__ << endl;
	
	CountTetraHedra(tIndx, pList);									//各四面体に含まれる粒子数のカウント，添え字の更新

	//粒子と四面体の情報登録
	for(unsigned i = 0; i < pList.size(); i++)
	{
		int freeIndx = m_ice->GetPtoTFreeIndx(pList[i]);
		m_ice->SetPtoT(pList[i], freeIndx, tIndx, i);				//粒子が所属している四面体を登録
	}

	m_ice->SetTtoP(tIndx, pList);									//四面体が含んでいる粒子を登録
	
	m_iTetraNum++;
	m_ice->SetTetraNum(m_iTetraNum);

//デバッグ
	//cout << "Debug After" << endl;
	//for(unsigned j = 0; j < pList.size(); j++)
	//{
	//	int jpIndx = pList[j];
	//	m_ice->DebugPtoT(jpIndx);
	//}
	//
	//m_ice->DebugTtoP(m_iTetraNum-1);
}

//-------------------------------------------------------------------------------------------------------------
//　クラスタ情報
//-------------------------------------------------------------------------------------------------------------
/*! 
 * クラスタ作成
 * @param[in] pIndx 粒子番号
 */
void rxFlGLWindow::MakeCluster(int pIndx)
{
	RXREAL *p = m_pPS->GetArrayVBO(rxParticleSystemBase::RX_POSITION);

	//((RXSPH*)m_pPS)->DetectSurfaceParticles();								//表面粒子検出
	//int *surfaceParticles = (int *)( ((RXSPH*)m_pPS)->GetArraySurf() );		//表面粒子
	//float weightedMass = 300.0f;
	//int layer = 1;

	//ある粒子が含まれる四面体を近傍クラスタとし，各四面体に含まれる粒子でクラスタを作成
	for(int i = 0; i < m_iceObj->GetPtoTIndx(pIndx); i++)
	{
		int itIndx = m_iceObj->GetPtoT(pIndx, i, 0);
		int ioIndx = m_iceObj->GetPtoT(pIndx, i, 1);

		if(itIndx == -1 || ioIndx == -1){ continue;	}

		//四面体の全ての粒子をクラスタに登録
		for(int j = 0; j < m_iceObj->GetTtoPIndx(itIndx); j++)
		{
			int jpIndx = m_iceObj->GetTtoP(itIndx, j);

			if(jpIndx == -1){	continue;	}
			if(m_iceObj->GetMoveObj(pIndx)->CheckIndx(jpIndx)){	continue;	}

			int pNum = m_iceObj->GetMoveObj(pIndx)->GetNumVertices();
			float mass = 1.0f;
			
			m_iceObj->GetMoveObj(pIndx)->AddVertex( 
					Vec3(p[jpIndx*4+0], p[jpIndx*4+1], p[jpIndx*4+2] ), mass, jpIndx
				);
			m_iceObj->GetMoveObj(pIndx)->SetAlphas(pNum, 1.0);
			m_iceObj->GetMoveObj(pIndx)->SetBetas (pNum, 0.0);
			m_iceObj->GetMoveObj(pIndx)->SetLayer (pNum, 0);
		}
	}

	//近傍四面体のlayerたどり，粒子を追加していく
	//TODO::不安定になるなら，layerが遠いほどbetaを下げる
	for(int i = 0; i < m_iceObj->GetPtoTIndx(pIndx); i++)
	{
		int itIndx = m_iceObj->GetPtoT(pIndx, i, 0);
		int ioIndx = m_iceObj->GetPtoT(pIndx, i, 1);

		if(itIndx == -1 || ioIndx == -1){	continue;	}

		for(int j = 0; j < m_iceObj->GetNTNum(itIndx); j++)
		{
			int jtIndx = m_iceObj->GetNeighborTetra(itIndx, j, 0);
			int jlIndx = m_iceObj->GetNeighborTetra(itIndx, j, 1);

			if(jtIndx == -1 || jlIndx == -1){	continue;	}

			//四面体の全ての粒子をクラスタに登録
			for(int k = 0; k < m_iceObj->GetTtoPIndx(jtIndx); k++)
			{
				int kpIndx = m_iceObj->GetTtoP(jtIndx, k);

				if(kpIndx == -1){	continue;	}
				if(m_iceObj->GetMoveObj(pIndx)->CheckIndx(kpIndx)){	/*cout << "contineue kpIndx = " << kpIndx << endl;*/ continue;	}

				int pNum = m_iceObj->GetMoveObj(pIndx)->GetNumVertices();
				float mass = 1.0f;

				m_iceObj->GetMoveObj(pIndx)->AddVertex(
						Vec3(p[kpIndx*4+0], p[kpIndx*4+1], p[kpIndx*4+2] ), mass, kpIndx
					);
				m_iceObj->GetMoveObj(pIndx)->SetAlphas(pNum, 1.0);
				m_iceObj->GetMoveObj(pIndx)->SetBetas (pNum, 0.0);
				m_iceObj->GetMoveObj(pIndx)->SetLayer (pNum, jlIndx);
			}
		}
	}
}

/*!
 * 近傍情報のみによるクラスタ作成
 * TODO::IceObject.cpp追加に合わせえて修正
 * @param[in]
 */
void rxFlGLWindow::MakeClusterFromNeight()
{
	////初期化のために影響半径を広くしてみる
	//float radius = ((RXSPH*)m_pPS)->GetEffectiveRadius();
	//((RXSPH*)m_pPS)->SetEffectiveRadius(radius * 3.0f);							//3.0で近傍粒子が500を超えるのに注意
	//StepPS(m_fDt);																//一度タイムステップを勧めないと，近傍粒子が取得されないみたい
	//((RXSPH*)m_pPS)->SetEffectiveRadius(radius);
	//
	//vector<vector<rxNeigh>>& neights = ((RXSPH*)m_pPS)->GetNeights();			//近傍粒子を取得
	//RXREAL *p = m_pPS->GetArrayVBO(rxParticleSystemBase::RX_POSITION);

	//rxSPHEnviroment sph_env = m_Scene.GetSphEnv();								// Shape Matchingの設定　パラメータ読み込み

	//for(int i = 0; i < m_iIcePrtNum; i++)
	//{
	//	//クラスタ初期化
	//	m_sm_cluster.push_back(new Ice_SM(m_iClusteresNum));
	//	m_sm_cluster[m_iClusteresNum]->SetSimulationSpace(-sph_env.boundary_ext, sph_env.boundary_ext);
	//	m_sm_cluster[m_iClusteresNum]->SetTimeStep(sph_env.smTimeStep);
	//	m_sm_cluster[m_iClusteresNum]->SetCollisionFunc(0);
	//	m_sm_cluster[m_iClusteresNum]->SetStiffness(1.0, 0.0);

	//	//近傍粒子をクラスタに追加
	//	for(int id_np = 0; id_np < neights[i].size(); id_np++)
	//	{
	//		int np_pIndx = neights[i][id_np].Idx;
	//		int pNum = m_sm_cluster[m_iClusteresNum]->GetNumVertices();

	//		m_sm_cluster[m_iClusteresNum]->AddVertex( Vec3(p[np_pIndx*4+0], p[np_pIndx*4+1], p[np_pIndx*4+2] ), 1.0, np_pIndx);
	//		m_sm_cluster[m_iClusteresNum]->SetAlphas(pNum, 1.0);
	//		m_sm_cluster[m_iClusteresNum]->SetBetas (pNum, 0.0);
	//		m_sm_cluster[m_iClusteresNum]->SetLayer (pNum, 0);
	//	}

	//	m_iClusteresNum++;
	//}
}

//相変化オブジェクトの運動計算
void rxFlGLWindow::StepIceObj()
{
#if defined(MK_USE_GPU)		//GPU使用
	StepClusterGPU(m_fDt);				RXTIMER("StepClusterGPU");
#else						//CPUのみ
	StepClusterCPU(m_fDt);				RXTIMER("StepClusterCPU");			//クラスタの運動計算
#endif
}

//熱処理
void rxFlGLWindow::StepHeatTransfer()
{
	//RXREAL *d  = m_pPS->GetArrayVBO(rxParticleSystemBase::RX_DENSITY);		//各粒子の密度
	//RXREAL *p  = m_pPS->GetArrayVBO(rxParticleSystemBase::RX_POSITION);		//各粒子の位置

	//((RXSPH*)m_pPS)->searchNeighbors();
	//((RXSPH*)m_pPS)->DetectSurfaceParticles();								//表面粒子検出

	//int* surfaceParticles = (int *)( ((RXSPH*)m_pPS)->GetArraySurf() );		//表面粒子
	//vector<vector<rxNeigh>>& neights = ((RXSPH*)m_pPS)->GetNeights();

	////ofstream ofs( "result/data/neightsDistance.txt");
	////ofs << "distance" << endl;
	////for(int i = 0; i < m_iIcePrtNum; i++){
	////	ofs << i << ":";
	////	for(int j = 0; j < neights[i].size(); j++){
	////		rxNeigh n = neights[i][j];
	////		int indx = n.Idx;
	////		double distance = n.Dist;
	////		ofs << " " << indx << ":" << distance << endl;
	////	}
	////}

	//float floor = -m_Scene.GetSphEnv().boundary_ext[2];
	//float radius = ((RXSPH*)m_pPS)->GetEffectiveRadius();

	//vector<int> nhs;
	//RXREAL h = radius * 3.0f;	//TODO: パラメータ化
	//Vec3 pos = m_pPS->GetSphereObstaclePos();

	//Vec3 min = ((RXSPH*)m_pPS)->GetMin();
	//rxSimParams param = ((RXSPH*)m_pPS)->GetParams();

	//uint* cellStart = ((RXSPH*)m_pPS)->GetCellStartInfo();
	//uint* cellEnd = ((RXSPH*)m_pPS)->GetCellEndInfo();
	//uint* sortedIndx = ((RXSPH*)m_pPS)->GetSortedIndx();

	//// 分割セルインデックスの算出
	//int x = (pos[0]-min[0])/param.CellWidth.x;
	//int y = (pos[1]-min[1])/param.CellWidth.y;
	//int z = (pos[2]-min[2])/param.CellWidth.z;

	//int numArdGrid = (int)(h/param.CellWidth.x+1);
	//for(int k = -numArdGrid; k <= numArdGrid; ++k){
	//	for(int j = -numArdGrid; j <= numArdGrid; ++j){
	//		for(int i = -numArdGrid; i <= numArdGrid; ++i){
	//			int i1 = x+i;
	//			int j1 = y+j;
	//			int k1 = z+k;
	//			if(i1 < 0 || (unsigned)i1 >= param.GridSize.x || j1 < 0 || (unsigned)j1 >= param.GridSize.y || k1 < 0 || (unsigned)k1 >= param.GridSize.z){
	//				continue;
	//			}

	//			RXREAL h2 = h*h;

	//			uint grid_hash = k1*param.GridSize.y*param.GridSize.x+j1*param.GridSize.x+i1;

	//			uint start_index = cellStart[grid_hash];
	//			if(start_index != 0xffffffff){	// セルが空でないかのチェック
	//				uint end_index = cellEnd[grid_hash];
	//				for(uint j = start_index; j < end_index; ++j){
	//					uint idx = sortedIndx[j];

	//					Vec3 xij;
	//					xij[0] = pos[0]-p[DIM*idx+0];
	//					xij[1] = pos[1]-p[DIM*idx+1];
	//					xij[2] = pos[2]-p[DIM*idx+2];

	//					rxNeigh neigh;
	//					neigh.Dist2 = norm2(xij);

	//					if(neigh.Dist2 <= h2){
	//						neigh.Idx = idx;
	//						nhs.push_back(neigh.Idx);
	//					}
	//				}
	//			}


	//		}
	//	}
	//}

	//m_iceObj->StepHeatTransfer(surfaceParticles, neights, nhs, floor, radius, p, d);

//GPU
	((RXSPH*)m_pPS)->DetectNeighborParticlesGPU();			//近傍粒子探索
	((RXSPH*)m_pPS)->DetectSurfaceParticlesGPU();			//表面粒子検出

	int* neight = ((RXSPH*)m_pPS)->GetNeightsGPU();
	int* surfPrt = ((RXSPH*)m_pPS)->GetArraySurfGPU();
	float* densty = ((RXSPH*)m_pPS)->GetDenstyGPU();
	float radius = ((RXSPH*)m_pPS)->GetEffectiveRadius();

	//オブジェクトの近傍粒子
	RXREAL *p  = m_pPS->GetArrayVBO(rxParticleSystemBase::RX_POSITION);		//各粒子の位置
	((RXSPH*)m_pPS)->copyCellInfo();

	int* nhs = new int[1000];
	int nhsIndx = 1;
	for(int i = 0; i < 1000; i++){	nhs[i] = -1;	}
	nhs[0] = 0;

	RXREAL h = radius * 3.0f;	//TODO: パラメータ化
	Vec3 pos = m_pPS->GetSphereObstaclePos();

	Vec3 min = ((RXSPH*)m_pPS)->GetMin();
	rxSimParams param = ((RXSPH*)m_pPS)->GetParams();

	uint* cellStart = ((RXSPH*)m_pPS)->GetCellStartInfo();
	uint* cellEnd = ((RXSPH*)m_pPS)->GetCellEndInfo();
	uint* sortedIndx = ((RXSPH*)m_pPS)->GetSortedIndx();

	// 分割セルインデックスの算出
	int x = (pos[0]-min[0])/param.CellWidth.x;
	int y = (pos[1]-min[1])/param.CellWidth.y;
	int z = (pos[2]-min[2])/param.CellWidth.z;

	int numArdGrid = (int)(h/param.CellWidth.x+1);
	for(int k = -numArdGrid; k <= numArdGrid; ++k){
		for(int j = -numArdGrid; j <= numArdGrid; ++j){
			for(int i = -numArdGrid; i <= numArdGrid; ++i){
				int i1 = x+i;
				int j1 = y+j;
				int k1 = z+k;
				if(i1 < 0 || (unsigned)i1 >= param.GridSize.x || j1 < 0 || (unsigned)j1 >= param.GridSize.y || k1 < 0 || (unsigned)k1 >= param.GridSize.z){
					continue;
				}

				RXREAL h2 = h*h;

				uint grid_hash = k1*param.GridSize.y*param.GridSize.x+j1*param.GridSize.x+i1;

				uint start_index = cellStart[grid_hash];
				if(start_index != 0xffffffff){	// セルが空でないかのチェック
					uint end_index = cellEnd[grid_hash];
					for(uint j = start_index; j < end_index; ++j){
						uint idx = sortedIndx[j];

						Vec3 xij;
						xij[0] = pos[0]-p[DIM*idx+0];
						xij[1] = pos[1]-p[DIM*idx+1];
						xij[2] = pos[2]-p[DIM*idx+2];

						float dist = norm2(xij);

						if(dist <= h2){
							if(nhsIndx >= 1000){ break; }
							nhs[nhsIndx++] = idx;
						}
					}
				}


			}
		}
	}
	
	nhs[0] = nhsIndx-1;	//近傍粒子数

	m_iceObj->StepHeatTransferGPU(surfPrt, neight, densty, nhs, radius);

	delete[] nhs;
}

//相変化処理
void rxFlGLWindow::StepIceStructure()
{
	m_iceObj->StepIceStructure();
}

/*!
 * ShapeMatching法のタイムステップを進める　CPU
 * @param[in] dt タイムステップ幅
 */
void rxFlGLWindow::StepClusterCPU(double dt)
{//	cout << __FUNCTION__ << endl;

	//デバイスのデータをホストへコピー
	RXREAL *p = m_pPS->GetArrayVBO(rxParticleSystemBase::RX_POSITION);
	RXREAL *v = m_pPS->GetArrayVBO(rxParticleSystemBase::RX_VELOCITY);
	
	m_iceObj->SetSPHHostPointer(p, v);

	//int* surfaceParticles = (int *)( ((RXSPH*)m_pPS)->GetArraySurf() );		//表面粒子

	//固体の運動計算
	m_iceObj->StepObjMoveCPU();				

	//テスト
	//m_iceObj->TestFixUpperPos();	//粒子を上部で固定
	//m_iceObj->TestFixSidePos();	//粒子を側面で固定

	//SPHのデータの更新　位置・速度
	m_pPS->SetArrayVBO(rxParticleSystemBase::RX_POSITION, p, 0, m_pPS->GetNumParticles());
	m_pPS->SetArrayVBO(rxParticleSystemBase::RX_VELOCITY, v, 0, m_pPS->GetNumParticles());
}

/*!
 * ShapeMatching法のタイムステップを進める　GPUで並列処理
 * @param[in] dt タイムステップ幅
 */
void rxFlGLWindow::StepClusterGPU(double dt)
{
	//位置情報はVBOを使っているので，マップしないといけない
	//頂点バッファオブジェクトをマップ
	float* dvbo_pos = ((RXSPH*)m_pPS)->GetDevicePointer_Pos();
	cudaGraphicsResource* d_resource_pos = ((RXSPH*)m_pPS)->GetDevicePointer_PosVBO();

	cudaGraphicsMapResources(1, &d_resource_pos, 0);
	cudaGraphicsResourceGetMappedPointer((void**)&dvbo_pos, NULL, d_resource_pos);

	//マップしたポインタを渡す
	m_iceObj->SetSPHDevicePointer(dvbo_pos, ((RXSPH*)m_pPS)->GetDevicePointer_Vel());

	//運動計算
	m_iceObj->StepObjMoveGPU();

	//位置情報をアンマップ
	cudaGraphicsUnmapResources(1, &d_resource_pos, 0);
}

void rxFlGLWindow::ResetSimulation()
{
	m_iceObj->ResetPosAndVel();
	ApplyNowData();

	g_iTimeCount = 0;
	m_iCurrentStep = 0;
	g_fTotalTime = 0;
	g_fAvgTime = 0;

	redraw();
}

void rxFlGLWindow::ApplyNowData()
{
	RXREAL *p = m_iceObj->GetSPHHostPosPointer();
	RXREAL *v = m_iceObj->GetSPHHostVelPointer();

	//SPHのデータの更新　位置・速度
	m_pPS->SetArrayVBO(rxParticleSystemBase::RX_POSITION, p, 0, m_pPS->GetNumParticles());
	m_pPS->SetArrayVBO(rxParticleSystemBase::RX_VELOCITY, v, 0, m_pPS->GetNumParticles());
}

//------------------------------------------------------------------------------------------------------------
// 固体構造関数
//------------------------------------------------------------------------------------------------------------
/*!
 * クラスタに含まれている情報の登録
 * @param[in] tNum クラスタ番号
 */
//バグだらけのコードだった
void rxFlGLWindow::MakeClusterInfo(int cIndx)
{	//cout << __FUNCTION__ << " start cIndx = " << cIndx << endl;
	//カウント
	//CountSolid(cIndx);

	//粒子とクラスタの情報登録
	vector<int> pIndxList;
	int pNum = m_iceObj->GetMoveObj(cIndx)->GetNumVertices();
	int* pLayerList = new int[pNum];									//粒子の所属レイヤー

	//layerのためのコピー
	for(int i = 0; i < pNum; i++)
	{
		pLayerList[i] = m_iceObj->GetMoveObj(cIndx)->GetLayer(i);		//粒子が何層目の近傍なのかを取得
	}

	for(int i = 0; i < m_ice->GetCtoPNum(cIndx); i++)
	{
		int pIndx = m_iceObj->GetMoveObj(cIndx)->GetParticleIndx(i);
		int freeIndx = m_ice->GetPtoCFreeIndx(pIndx);

		m_ice->SetPtoC(pIndx, freeIndx, cIndx, i, pLayerList[i]);		//粒子が所属している四面体を登録
		pIndxList.push_back(pIndx);
	}

	m_ice->SetCtoP(cIndx, pIndxList, pLayerList);						//四面体が含んでいる粒子を登録

	delete[] pLayerList;
}

/*!
 * 四面体削除時の粒子・四面体情報の更新
 * @param tList 四面体番号参照リスト
 */
void rxFlGLWindow::UpdateInfo_Delete_TandP(const vector<int>& tList, const vector<int>& deleteList)
{//	cout << __FUNCTION__ << endl;
	if(tList.size() == 0){	return; }
	
	for(unsigned i = 0; i < deleteList.size(); i++)
	{
		int itIndx = tList[deleteList[i]];

		for(int j = 0; j < m_ice->GetTtoPIndx(itIndx); j++)
		{
			int jpIndx = m_ice->GetTtoP(itIndx, j);
			if(jpIndx == -1){	continue;	}

			for(int k = 0; k < m_ice->GetPtoTIndx(jpIndx); k++)
			{
				int ktIndx = m_ice->GetPtoT(jpIndx, k, 0);
				int koIndx = m_ice->GetPtoT(jpIndx, k, 1);

				if(ktIndx == -1 || koIndx == -1 || koIndx != itIndx){	continue;	}

				m_ice->DeletePtoT(jpIndx, k);
				break;
			}
		}

		m_ice->ClearTtoP(itIndx);
		m_iTetraNumNum--;
	}
}

/*!
 * 同一，包含関係にある四面体を削除		//TODO::この関数はIceStructureに持たせるべき
 */
void rxFlGLWindow::CheckDeleteTetra(vector<int>& tList, vector<int>& lList)
{
	//cout << "before tList.size =  " << tList.size() << " m_iTetraNumNum = " << m_iTetraNumNum << endl;

	//四面体で内容が同一，または包含関係にあるものを削除．
	//TODO::並列化をそれぞれで行う
	for(vector<int>::iterator indx = tList.begin(); indx !=tList.end(); indx++)
	{
		if(*indx == -1){	continue;	}

		vector<int> deleteTList;									//削除する四面体の添え字

		CheckSameTetra(*indx, tList, deleteTList);					//同一内容の四面体の添え字取得
		CheckIncludeTetra(*indx, tList, deleteTList);				//内包している四面体の添え字取得
		UpdateInfo_Delete_TandP(tList, deleteTList);				//四面体に関する情報の削除
		RemoveReconstructTetra(tList, lList, deleteTList);			//削除した四面体が再定義する四面体に含まれるなら取り除く
	}

	//-1である四面体データを削除
	tList.erase(std::remove(tList.begin(), tList.end(), -1), tList.end());
	lList.erase(std::remove(lList.begin(), lList.end(), -1), lList.end());

	//cout << "after tList.size =  " << tList.size() << " m_iTetraNumNum = " << m_iTetraNumNum << endl;
}

/*!
 * 同一の四面体が存在するかの判定
 * @param[in] tIndx 判定する四面体番号
 * @param[in] searchTList 探索対象となる四面体リスト
 * @param[in] deleteList 結果を追加する配列
 */
void rxFlGLWindow::CheckSameTetra(int tIndx, const vector<int>& searchTList, vector<int>& deleteList)
{
	int tNum = m_ice->GetTtoPNum(tIndx);

	if(tNum == 4 || tNum == 1){	return;	}								//所属粒子数４なら同一は存在しない（はず）

	//同一内容を探索
	for(unsigned i = 0; i < searchTList.size(); i++)
	{
		int itIndx = searchTList[i];
		if(itIndx < 0)						{	continue;	}
		if(tIndx == itIndx)					{	continue;	}
		if(m_ice->GetTtoPNum(itIndx) == 0)	{	continue;	}
		if(m_ice->GetTtoPNum(itIndx) == 4)	{	continue;	}
		if(tNum != m_ice->GetTtoPNum(itIndx)){	continue;	}	//所属する粒子数のチェック

		//tIndxとi番目の四面体の内容が，完全に一致しているかのチェック
		bool check = false;
		for(int j = 0; j < m_ice->GetTtoPIndx(tIndx); j++)
		{
			check = false;

			int jtIndx = m_ice->GetTtoP(tIndx, j);
			if(jtIndx == -1){	continue;	}

			//j番目の粒子が共通しているかをチェック
			for(int k = 0; k < m_ice->GetTtoPIndx(itIndx); k++)
			{
				int ktIndx = m_ice->GetTtoP(itIndx, k);
				if(ktIndx == -1){	continue;	}

				if(jtIndx == ktIndx){	check = true; break;}
				else				{	check = false;		}
			}

			if(!check){	break;	}				//1つでも含んでいないなら不一致 trueなら次の粒子判定へ
		}

		//最後までtrueなら同一内容．
		if(check)
		{	
			deleteList.push_back(i);
		}
	}

	//デバッグ
	//if(tList.size() == 0){	return;	}
	//cout << "Debug::" << __FUNCTION__ << endl;
	//m_ice->DebugTtoP(tIndx);
	//for(int i = 0; i < tList.size(); i++)
	//{
	//	m_ice->DebugTtoP(tList[i]);
	//}
}

/*!
 * 内包する四面体が存在するかの判定
 * @param[in] tIndx 判定する四面体番号
 * @param[in] searchTList　探索する四面体リスト
 * @param[in] deleteList 結果を追加する配列
 */
void rxFlGLWindow::CheckIncludeTetra(int tIndx, const vector<int>& searchTList, vector<int>& deleteList)
{
	int tNum = m_ice->GetTtoPNum(tIndx);

	if(tNum == 1){	return;	}								//所属粒子数1なら内包は存在しない（はず）

	//内包する四面体を再定義する四面体から探索
	for(unsigned i = 0; i < searchTList.size(); i++)
	{
		int itIndx = searchTList[i];
		if(itIndx < 0)						{	continue;	}
		if(tIndx == itIndx)					{	continue;	}
		if(m_ice->GetTtoPNum(itIndx) == 0)	{	continue;	}
		if(m_ice->GetTtoPNum(itIndx) == 4)	{	continue;	}
		if(tNum <= m_ice->GetTtoPNum(itIndx)){	continue;	}	//所属する粒子数のチェック

		//i番目の四面体の内容がtIndx番目の四面体に全て含まれているかのチェック
		bool check = false;
		for(int j = 0; j < m_ice->GetTtoPIndx(itIndx); j++)
		{
			check = false;

			int jtIndx = m_ice->GetTtoP(itIndx, j);
			if(jtIndx == -1){	continue;	}

			//j番目の粒子が共通しているかをチェック
			for(int k = 0; k < m_ice->GetTtoPIndx(tIndx); k++)
			{
				int ktIndx = m_ice->GetTtoP(tIndx, k);
				if(ktIndx == -1){	continue;	}

				if(jtIndx == ktIndx){	check = true; break;}
				else				{	check = false;		}
			}

			if(!check){	break;	}				//1つでも含んでいないなら不一致 trueなら次の粒子判定へ
		}

		//最後までtrueなら含んでいる．
		if(check)
		{
			deleteList.push_back(i);
		}		
	}

	//デバッグ
	//if(tList.size() == 0){	return;	}
	//cout << "Debug::" << __FUNCTION__ << endl;
	//m_ice->DebugTtoP(tIndx);
	//for(int i = 0; i < tList.size(); i++)
	//{
	//	m_ice->DebugTtoP(tList[i]);
	//}
}

/*
 * 再定義する四面体から不要とみなされた四面体を-1に書き換える
 */
void rxFlGLWindow::RemoveReconstructTetra(vector<int>& tList, vector<int>& lList, vector<int>& deleteTList)
{
	for(vector<int>::iterator indx = deleteTList.begin(); indx != deleteTList.end(); indx++)
	{
		tList[*indx] = -1;
		lList[*indx] = -1;
	}
}

/*!
 * 粒子・クラスタ情報の再定義
 * @param pList 粒子リスト
 * @param cList 再定義クラスタリスト
 * @param lList レイヤーリスト
 */
void rxFlGLWindow::SetClusterInfo(const vector<int>& pList, const vector<int>& cList, const vector<int>& lList)
{	//cout << __FUNCTION__ << endl;
	if(pList.size() == 0 || cList.size() == 0){	return;	}
	
	vector<int> checkTList;
	int j = 0, k = 0, l = 0;
	int icIndx = 0;
	int jtIndx = 0, joIndx = 0;
	int kpIndx = 0, ktIndx = 0, klIndx = 0;
	int lpIndx = 0;
	int pNum = 0;

	RXREAL *p = m_pPS->GetArrayVBO(rxParticleSystemBase::RX_POSITION);

	//TODO::まずは，lListを用いないでやってみる．おそらく計算時間はそれほどかからないはずだが…

	//SM法クラスタの初期化
	#pragma omp parallel
	{
	#pragma omp for
		for(int i = 0; i < (int)cList.size(); ++i)
		{
			if(std::find(pList.begin(), pList.end(), cList[i]) != pList.end()){	continue;	}
			m_iceObj->GetMoveObj(cList[i])->Clear();
		}
	}//end #pragma omp parallel

	//SM法クラスタの再定義
	#pragma omp parallel
	{
	#pragma omp for private(checkTList, j, k, l, icIndx, jtIndx, joIndx, kpIndx, ktIndx, klIndx, lpIndx, pNum)
		for(int i = 0; i < (int)cList.size(); ++i)
		{
			checkTList.clear();
			icIndx = cList[i];
			if(std::find(pList.begin(), pList.end(), icIndx) != pList.end()){	continue;	}

			//クラスタを再定義する際，基本となる粒子が属する四面体から初期粒子を得る．
			//クラスタ番号＝＝粒子番号なのに注意
			//以下を関数にすると，エラーが出てうまくいかない
			for(j = 0; j < m_ice->GetPtoTIndx(icIndx); j++)
			{
				jtIndx = m_ice->GetPtoT(icIndx, j, 0);
				joIndx = m_ice->GetPtoT(icIndx, j, 1);

				if(jtIndx == -1 || joIndx == -1){ continue;	}
				if(std::find(checkTList.begin(), checkTList.end(), jtIndx) != checkTList.end())
				{	continue;	}
				else
				{	checkTList.push_back(jtIndx);	}

				//四面体の全ての粒子をクラスタに登録
				for(k = 0; k < m_ice->GetTtoPIndx(jtIndx); k++)
				{
					kpIndx = m_ice->GetTtoP(jtIndx, k);

					if(kpIndx == -1){	continue;	}
					if(m_iceObj->GetMoveObj(icIndx)->CheckIndx(kpIndx)){	continue;	}

					pNum = m_iceObj->GetMoveObj(icIndx)->GetNumVertices();

					m_iceObj->GetMoveObj(icIndx)->AddVertex(Vec3(p[kpIndx*4+0], p[kpIndx*4+1], p[kpIndx*4+2]), 1.0, kpIndx);
					m_iceObj->GetMoveObj(icIndx)->SetAlphas(pNum, 1.0);
					m_iceObj->GetMoveObj(icIndx)->SetBetas (pNum, 0.0);
					m_iceObj->GetMoveObj(icIndx)->SetLayer (pNum, 0);
				}
			}

			//近傍四面体のlayerたどり，粒子を追加していく
			//TODO::不安定になるなら，layerが遠いほどbetaを下げる
			for(j = 0; j < m_ice->GetPtoTIndx(icIndx); ++j)
			{
				jtIndx = m_ice->GetPtoT(icIndx, j, 0);
				joIndx = m_ice->GetPtoT(icIndx, j, 1);

				if(jtIndx == -1 || joIndx == -1){ continue;	}

				for(k = 0; k < m_ice->GetNTNum(jtIndx); k++)
				{
					ktIndx = m_ice->GetNeighborTetra(jtIndx, k, 0);
					klIndx = m_ice->GetNeighborTetra(jtIndx, k, 1);

					if(ktIndx == -1 || klIndx == -1){	continue;	}
					if(std::find(checkTList.begin(), checkTList.end(), ktIndx) != checkTList.end())
					{	continue;	}
					else
					{	checkTList.push_back(ktIndx);	}

					//四面体の全ての粒子をクラスタに登録
					for(l = 0; l < m_ice->GetTtoPIndx(ktIndx); l++)
					{
						lpIndx = m_ice->GetTtoP(ktIndx, l);

						if(lpIndx == -1){	continue;	}
						if(m_iceObj->GetMoveObj(icIndx)->CheckIndx(lpIndx)){	/*cout << "contineue kpIndx = " << kpIndx << endl;*/ continue;	}

						pNum = m_iceObj->GetMoveObj(icIndx)->GetNumVertices();

						m_iceObj->GetMoveObj(icIndx)->AddVertex( Vec3(p[lpIndx*4+0], p[lpIndx*4+1], p[lpIndx*4+2] ), 1.0, lpIndx);
						m_iceObj->GetMoveObj(icIndx)->SetAlphas(pNum, 1.0);
						m_iceObj->GetMoveObj(icIndx)->SetBetas (pNum, 0.0);
						m_iceObj->GetMoveObj(icIndx)->SetLayer (pNum, klIndx);
					}
				}
			}
		}
	}//end #pragma omp parallel

	//固体情報の初期化
	//#pragma omp parallel
	//{
	//#pragma omp for
		for(int i = 0; i < (int)cList.size(); ++i)
		{
			if(std::find(pList.begin(), pList.end(), cList[i]) != pList.end())
			{
				continue;
			}
			m_ice->ClearCtoP(cList[i]);
		}
	//}//end #pragma omp parallel

	//固体情報の登録
	for(unsigned i = 0; i < cList.size(); ++i)
	{
		if(std::find(pList.begin(), pList.end(), cList[i]) != pList.end())
		{
			continue;	
		}
		MakeClusterInfo(cList[i]);
	}

}

/*!
 * 粒子融解時の粒子・四面体情報の再定義
 * @param pList 融解粒子リスト
 * @param tList 再定義四面体リスト
 * @param lList レイヤーリスト
 */
void rxFlGLWindow::SetTetraInfo(const vector<int>& pList, const vector<int>& tList, const vector<int>& lList)
{	
	if(pList.size() == 0 || tList.size() == 0){	return;	}

	//cout << __FUNCTION__ << endl;

	int itIndx = 0;
	int ilayer = 0;

	#pragma omp parallel
	{
	#pragma omp for private(itIndx,ilayer)
		//近傍四面体の再定義
		for(int i = 0; i < (int)tList.size(); i++)
		{
			itIndx = tList[i];
			ilayer = lList[i];

			m_ice->ClearNeighborTetraFromLayer(itIndx, ilayer);
			m_ice->SetNeighborTetraFromLayer(itIndx, m_iLayer, ilayer);	//ここが非常に重い
		}
	}//end #pragma omp parallel
}

/*!
 * 凝固粒子に関する情報の作成　クラスタ
 */
void rxFlGLWindow::SetFreezeClusterInfo(const vector<int>& pList)
{	//cout << __FUNCTION__ << " start pList.size() = " << pList.size() << endl;

	for(unsigned i = 0; i < pList.size(); i++)
	{
		int pIndx = pList[i];
		if(m_ice->GetPtoTNum(pIndx) == 0){	continue;	}	//四面体に属しない場合は戻る

		//pIndx == cIndxa なのに注意
		MakeCluster(pIndx);
		MakeClusterInfo(pIndx);
	}
}

/*!
 * 凝固後，再定義するクラスタの探索
 * @param[in] pList   凝固粒子リスト
 * @param[in] cList   再定義するクラスタの参照リスト
 * @param[in] lList   再定義するクラスタのレイヤー参照リスト
 */
void rxFlGLWindow::SearchReconstructCluster_Freeze(const vector<int>& pList, vector<int>& cList, vector<int>& lList)
{
	if(pList.size() == 0){	return; }
	//cout << __FUNCTION__ << endl;

	//凝固粒子が所属するクラスタに含まれている各粒子のクラスタが，再定義クラスタ．
	for(unsigned i = 0; i < pList.size(); i++)
	{
		int ipIndx = pList[i];
		if(m_ice->GetPtoCNum(ipIndx) == 0){		continue;	}

		//凝固粒子が含まれているクラスタを取得
		for(int j = 0; j < m_ice->GetPtoCIndx(ipIndx); j++)
		{
			int jcIndx = m_ice->GetPtoC(ipIndx, j, 0);
			int joIndx = m_ice->GetPtoC(ipIndx, j, 1);
			int jlIndx = m_ice->GetPtoC(ipIndx, j, 2);

			if(jcIndx == -1 || joIndx == -1){	continue;	}
			if(jcIndx == ipIndx)			{	continue;	}
			if(std::find(cList.begin(), cList.end(), jcIndx) != cList.end()){	continue;	}
			//if(std::find(pList.begin(), pList.end(), jcIndx) != pList.end()){	continue;	}

			cList.push_back(jcIndx);
			lList.push_back(1);

			//クラスタに含まれている粒子＝クラスタを再定義クラスタとする
			//for(int k = 0; k < m_ice->GetCtoPNum(jcIndx); k++)
			//{
			//	int kcIndx = m_ice->GetCtoP(jcIndx, k, 0);
			//	int koIndx = m_ice->GetCtoP(jcIndx, k, 1);
			//	int klIndx = m_ice->GetCtoP(jcIndx, k, 2);

			//	if(kcIndx == -1 || koIndx == -1){			continue;	}
			//	if(kcIndx == ipIndx || kcIndx == jcIndx){	continue;	}
			//	if(std::find(cList.begin(), cList.end(), kcIndx) != cList.end()){	continue;	}

			//	cList.push_back(kcIndx);
			//	lList.push_back(1);
			//}
		}
	}
}

/*!
 * 凝固後，近傍を再定義する四面体の探索
 * @param[in] pList   凝固粒子リスト
 * @param[in] tList   再定義する四面体のリスト
 * @param[in] lList   再定義する四面体のレイヤーリスト
 */
void rxFlGLWindow::SearchReconstructTetra_Freeze(const vector<int>& pList, vector<int>& tList, vector<int>& lList)
{	
	if(pList.size() == 0){	return; }

	//今回新しく作られた四面体群Aの近傍四面体群をA'とすると，A'の近傍四面体へAを追加する必要がある．
	//しかし追加する際のlayerが昇順でないという問題があるので，計算時間はかかるかもしれないが単純に再定義する．
	
	//凝固粒子が所属する四面体の，近傍四面体全てを再定義する．

	m_ice->ResetTFlag(m_iTetraNum);							//四面体探索フラグの初期化

	//１，凝固粒子が含まれていた四面体
	for(unsigned i = 0; i < pList.size(); i++)
	{
		int ipIndx = pList[i];

		for(int j = 0; j < m_ice->GetPtoTIndx(ipIndx); j++)
		{
			int jtIndx = m_ice->GetPtoT(ipIndx, j, 0);
			int jcIndx = m_ice->GetPtoT(ipIndx, j, 1);

			if(jtIndx == -1 || jcIndx == -1){	continue;	}

			//if(std::find(tList.begin(), tList.end(), jtlSet[0]) != tList.end()){	continue;	}
			if(m_ice->GetTFlag(jtIndx))	{	continue;	}
			else						{	m_ice->SetTFlag(jtIndx, true);	}

			tList.push_back(jtIndx);
			lList.push_back(1);								//0か1かの判断はできないので1に合わせる．メモを参考に．
		}
	}

	//２，１の四面体の近傍四面体
	int tetraNum = tList.size();
	for(int i = 0; i < tetraNum; i++)
	{
		int itIndx = tList[i];

		for(int j = 0; j < m_ice->GetNTNum(itIndx); j++)
		{
			int jtIndx = m_ice->GetNeighborTetra(itIndx, j, 0);
			int jlIndx = m_ice->GetNeighborTetra(itIndx, j, 1);

			//既に含まれているのなら，layerを比べて小さいほうを優先する
			if(m_ice->GetTFlag(jtIndx))
			{
				vector<int>::iterator check = std::find(tList.begin(), tList.end(), jtIndx);

				int layerIndx = check-tList.begin();
				if(lList[layerIndx] > jlIndx)
				{
					lList[layerIndx] = jlIndx;
				}
				continue;
			}
			else
			{
				m_ice->SetTFlag(jtIndx, true);	
			}

			tList.push_back(jtIndx);
			lList.push_back(jlIndx);
		}
	}
}

/*!
 * 最近傍最短粒子の取得　layer = 0なら最大27程度
 * @param[in] pIndx　粒子番号
 * @param[in] layer　近傍層数
 * @param[in] num　　最大取得粒子数
 * @param[out] vector<int> 最近傍最短粒子の集合　距離の近い順にソートされている
 */
vector<int> rxFlGLWindow::GetNeighborDistanceIceParticles(int pIndx, int layer, int num)
{//	cout << __FUNCTION__ << " start" << endl;
	vector<int> pList;
	vector<float> distanceList;
	RXREAL *p = m_pPS->GetArrayVBO(rxParticleSystemBase::RX_POSITION);

	//四面体ベース
	//ある粒子が含まれている接続クラスタ全てを探索する．
	//for( int i = 0; i < m_ice->GetPtoCIndx_Connect(pIndx); i++ )
	//{
	//	int* coList = m_ice->GetPtoC_Connect( pIndx, i );
	//	int cIndx = coList[0];
	//	int oIndx = coList[1];
	//	if(cIndx == -1 || oIndx == -1) continue;

	//	for( int j = 0; j < m_sm_connects[cIndx]->GetNumVertices(); j++ )
	//	{
	//		int jpIndx = m_sm_connects[cIndx]->GetParticleIndx(j);

	//		if( pIndx == jpIndx ) continue;

	//		vector<int>::iterator check = std::find( pList.begin(),
	//												 pList.end()  ,
	//												 jpIndx			);
	//		if( check != pList.end() ){	continue;	}	//すでに含んでいたら戻る

	//		float distance = dot(Vec3(p[pIndx*DIM+0],  p[pIndx*DIM+1],  p[pIndx*DIM+2]),
	//							 Vec3(p[jpIndx*DIM+0], p[jpIndx*DIM+1], p[jpIndx*DIM+2]));
	//		distanceList.push_back(distance);
	//		pList.push_back(jpIndx);

	//		if( pList.size() >= (unsigned)num ){ break;}
	//	}
	//	if( pList.size() >= (unsigned)num ){ break;}
	//}

	//粒子ベース
	//ある粒子が含まれている四面体全てを探索する
	for(int i = 0; i < m_ice->GetPtoTIndx(pIndx); i++)
	{
		int tIndx = m_ice->GetPtoT(pIndx, i, 0);
		int oIndx = m_ice->GetPtoT(pIndx, i, 1);

		if(tIndx == -1 || oIndx == -1) continue;

		for(int j = 0; j < m_ice->GetTtoPIndx(tIndx); j++)
		{
			int jpIndx = m_ice->GetTtoP(tIndx, j);

			if(jpIndx == -1)		continue;
			if(pIndx == jpIndx)		continue;

			vector<int>::iterator check = std::find( pList.begin(),
													 pList.end()  ,
													 jpIndx			);
			if( check != pList.end() ){	continue;	}	//すでに含んでいたら戻る

			float distance = dot(Vec3(p[pIndx*DIM+0],  p[pIndx*DIM+1],  p[pIndx*DIM+2]),
								 Vec3(p[jpIndx*DIM+0], p[jpIndx*DIM+1], p[jpIndx*DIM+2]));
			distanceList.push_back(distance);
			pList.push_back(jpIndx);

			if( pList.size() >= (unsigned)num ){ break;}
		}
		if( pList.size() >= (unsigned)num ){ break;}
	}

	//ある粒子が含まれているクラスタの，近傍クラスタを何層か探索する．layer回探索．
	for( int i = 0; i < layer; i++ )
	{
	}

	//距離の近い順に粒子番号をソート
	for( unsigned i = 0; i < distanceList.size(); i++ )
	{
		for( unsigned j = i+1; j < distanceList.size(); j++ )
		{
			if( distanceList[j] < distanceList[i] )
			{
				float dist = distanceList[i];
				distanceList[i] = distanceList[j];
				distanceList[j] = dist;

				int pIndx = pList[i];
				pList[i] = pList[j];
				pList[j] = pIndx;
			}
		}
	}

	return pList;
}

/*
 * 粒子の色設定　温度，熱量，四面体，クラスタを切り替える．
 * 
 */
void rxFlGLWindow::StepParticleColor()
{
	//温度
	if( m_iColorType == rxParticleSystemBase::RX_TEMP )
	{
		//真っ黒では見えにくかったので，ちょっと青いのがデフォルトとした．

		////CPU
		//float* tempColor = new float[m_pPS->GetNumParticles()];
		//for( int i = 0; i < m_pPS->GetNumParticles(); i++ ){
		//	tempColor[i] = m_iceObj->GetTemps()[i] + 100.0f;
		//}

		//m_pPS->SetColorVBOFromArray( tempColor, 1, false, 1000.0f );
		//delete[] tempColor;

		//GPU
		float* tempColor = new float[m_pPS->GetNumParticles()];
		m_iceObj->CopyTempDeviceToHost(tempColor);
		for( int i = 0; i < m_pPS->GetNumParticles(); i++ ){
			tempColor[i] += 100.0f;
		}

		m_pPS->SetColorVBOFromArray( tempColor, 1, false, 1000.0f );
		delete[] tempColor;

	}
	//表面粒子
	else if(m_iColorType == rxParticleSystemBase::RX_SURFACE)
	{
		////CPU
		//((RXSPH*)m_pPS)->DetectSurfaceParticles();								//表面粒子検出

		//int* surfaceParticles = (int *)( ((RXSPH*)m_pPS)->GetArraySurf() );		//表面粒子
		//vector<vector<rxNeigh>>& neights = ((RXSPH*)m_pPS)->GetNeights();

		//float* tempColor = new float[m_pPS->GetNumParticles()];

		//for( int i = 0; i < m_pPS->GetNumParticles(); i++ ){
		//	if(surfaceParticles[i] == 1){
		//		tempColor[i] = /*(float)neights[i].size()*/ 90.0f;
		//	}
		//	else{
		//		tempColor[i] = 0.0f;
		//	}
		//}

		//m_pPS->SetColorVBOFromArray( tempColor, 1, false, 100.0f );
		//delete[] tempColor;

		//GPU
		//((RXSPH*)m_pPS)->DetectSurfaceParticlesGPU();		//表面粒子検出
		int* dSurfPrt = ((RXSPH*)m_pPS)->GetArraySurfGPU();
		int* hSurfPrt = new int[m_pPS->GetNumParticles()];
		
		cudaMemcpy(hSurfPrt, dSurfPrt, sizeof(int) * m_pPS->GetNumParticles(), cudaMemcpyDeviceToHost);

		float* tempColor = new float[m_pPS->GetNumParticles()];

		for(int i = 0; i < m_pPS->GetNumParticles(); i++){
			if(hSurfPrt[i] == 1){
				tempColor[i] = 90.0f;
			}
			else{
				tempColor[i] = 0.0f;
			}
		}

		m_pPS->SetColorVBOFromArray(tempColor, 1, false, 100.0f);

		delete[] tempColor;
		delete[] hSurfPrt;
	}
	//四面体情報
	//TODO::いずれGUIから入力して操作できるようにする
	else if( m_iColorType == rxParticleSystemBase::RX_ICE_CONNECT )
	{
		//クラスタにより接続関係にある粒子に色をつける
		float* tempColor = new float[m_pPS->GetNumParticles()];
		for( int i = 0; i < m_pPS->GetNumParticles(); i++ ){	tempColor[i] = 100.0f;		}

		//描画モードで分岐
		if(m_iShowTetraIndx == m_iceObj->GetTetraNum())
		{	
			//四面体全体を表示
			for(int i = 0; i < m_iceObj->GetTetraNum(); i++)
			{	
				for(int j = 0; j < m_iceObj->GetTtoPIndx(i); j++)
				{
					int jpIndx = m_iceObj->GetTtoP(i, j);

					if(jpIndx == -1){	continue;	}
					tempColor[jpIndx] = 600.0f;
				}
			}
		}
		else
		{	
			//１つの四面体のみを表示
			for( int j = 0; j < m_iceObj->GetTtoPIndx(m_iShowTetraIndx); j++ )
			{	
				int jpIndx = m_iceObj->GetTtoP(m_iShowTetraIndx, j);

				if(jpIndx == -1){	continue;	}
				tempColor[jpIndx] = 600.0f;
			}
		}

		m_pPS->SetColorVBOFromArray( tempColor, 1, false, 1000.0f );
		delete[] tempColor;
	}
	//クラスタ
	else if( m_iColorType == rxParticleSystemBase::RX_ICE_CALC )
	{
		//クラスタにより接続関係にある粒子に色をつける
		float* tempColor = new float[m_pPS->GetNumParticles()];
		for( int i = 0; i < m_pPS->GetNumParticles(); i++ )
		{
			if(m_iceObj->GetPtoCNum(i) == 0)
			{
				tempColor[i] = 0.0f;
			}
			else
			{
				tempColor[i] = 100.0f;
			}
		}

		//粒子ベース版
		if( m_iShowClusterIndx == m_iClusteresNum )
		{	
			//クラスタ全体を表示
			for(int i = 0; i < m_iClusteresNum; i++)
			{
				for(unsigned j = 0; j < m_iceObj->GetMoveObj(i)->GetIndxNum(); j++ )
				{
					int jpIndx = m_iceObj->GetMoveObj(i)->GetParticleIndx(j);

					if(jpIndx == MAXINT){	continue;	}

					tempColor[jpIndx] = 700.0f;
				}
			}
		}
		else
		{	
			//１つのクラスタのみを表示
			for(unsigned j = 0; j < m_iceObj->GetMoveObj(m_iShowClusterIndx)->GetIndxNum(); j++ )
			{
				int jpIndx = m_iceObj->GetMoveObj(m_iShowClusterIndx)->GetParticleIndx(j);

				if(jpIndx == MAXINT){	continue;	}


				//クラスタと対になる粒子の色は変える
				if(m_iShowClusterIndx == jpIndx)
				{
					tempColor[jpIndx] = 666.0f;
				}
				else
				{
					tempColor[jpIndx] = 700.0f;
				}
			}
		}

		m_pPS->SetColorVBOFromArray( tempColor, 1, false, 1000.0f );
		delete[] tempColor;
	}		
	//高速化用パス
	else if( m_iColorType == rxParticleSystemBase::RX_ICE_FAST_PATH )
	{
	}
	//選択的運動計算
	else if( m_iColorType == rxParticleSystemBase::RX_ICE_SELECTED )
	{
		//運動計算する粒子を赤，それ以外を黒に
		float* tempColor = new float[m_pPS->GetNumParticles()];

		for( int i = 0; i < m_pPS->GetNumParticles(); i++ )
		{
			//適当に，奇数番だけを運動計算してみる
			if(m_iceObj->GetMotionCalcClsuter(i) != 0)
			{
				tempColor[i] = 1000.0f;
			}
			else
			{
				tempColor[i] = 0.0f;
			}
		}

		m_pPS->SetColorVBOFromArray( tempColor, 1, false, 1000.0f );

		delete[] tempColor;
	}
	//変形クラスタの検出
	else if( m_iColorType == rxParticleSystemBase::RX_ICE_DEFORMATION )
	{
		float* tempColor = new float[m_pPS->GetNumParticles()];

		float threshold = 0.1f;
		int defNum = 0;

		////相対的な変形量
		////ある程度正確に検出できるが，耳の付け根をガッツリ検出する，といったことができない．
		////平均や標準偏差の定義上，近傍に含まれるうちの幾つかしか検出することができない．
		//for(int i = 0; i < m_pPS->GetNumParticles(); i++){
		//	//平均値
		//	float defAmountSum = 0.0f;
		//	for(int j = 0; j < m_iceObj->GetOrientedObj(i)->GetIndxNum(); j++){
		//		int pIndx = m_iceObj->GetOrientedObj(i)->GetParticleIndx(j);

		//		float defAmount = m_iceObj->GetOrientedObj(pIndx)->FinalDefAmount();
		//		if(defAmount < 0.0f) defAmount = 0.0f;

		//		defAmountSum += defAmount;
		//	}

		//	float defAmountAve = defAmountSum / m_iceObj->GetOrientedObj(i)->GetNumVertices();

		//	//分散
		//	float variance = 0.0f;
		//	for(int j = 0; j < m_iceObj->GetOrientedObj(i)->GetIndxNum(); j++){
		//		int pIndx = m_iceObj->GetOrientedObj(i)->GetParticleIndx(j);

		//		float diff = m_iceObj->GetOrientedObj(pIndx)->FinalDefAmount() - defAmountAve;
		//		variance += diff * diff;
		//	}

		//	variance /= m_iceObj->GetOrientedObj(i)->GetNumVertices();

		//	//標準偏差
		//	float deviation = sqrt(variance);

		//	//平均より大きく上ならば変形部位
		//	//float defAmount = m_iceObj->GetOrientedObj(i)->FinalDefAmount();
		//	//if(defAmount < 0.0f) defAmount = 0.0f;
		//
		//	//if(defAmount > defAmountAve + deviation /* || pNumAve > 20*/){
		//	//	tempColor[i] = 0.83f;	//赤
		//	//	defNum++;
		//	//}
		//	//else{
		//	//	tempColor[i] = 0.15f;	//青
		//	//}

		//	//標準偏差より下なら非変形部位
		//	float defAmount = m_iceObj->GetOrientedObj(i)->FinalDefAmount();
		//	if(defAmount < 0.0f) defAmount = 0.0f;
		//
		//	if(defAmount < defAmountAve + deviation * 1.0f /* || pNumAve > 20*/){
		//		tempColor[i] = 0.15f;	//青
		//	}
		//	else{
		//		tempColor[i] = 0.83f;	//赤
		//		defNum++;
		//	}
		//}

		//cout << defNum << endl;

		////非変形部位の検出
		//for(int i = 0; i < m_pPS->GetNumParticles(); i++){

		//	if(m_iceObj->GetOrientedObj(i)->IsThereOutlier()){
		//		tempColor[i] = 0.15f;	//青
		//	}
		//	else{
		//		tempColor[i] = 0.95f;	//赤
		//		defNum++;
		//	}
		//}

		////外れ値の検出
		//for(int i = 0; i < m_pPS->GetNumParticles(); i++){

		//	if(m_iceObj->GetOrientedObj(i)->IsThereOutlier()){
		//		tempColor[i] = 0.95f;	//赤
		//		defNum++;
		//	}
		//	else{
		//		tempColor[i] = 0.15f;	//青
		//	}
		//}

		////体積変化率
		////妙な色の分布になる　内部に色がついて，表面粒子は黒　大きな変形には反応する
		//for(int i = 0; i < m_pPS->GetNumParticles(); i++){
		//	float defAmount = abs(m_iceObj->GetOrientedObj(i)->VolumeChange()*2000.0f);
		//	if(defAmount < 0.0f) defAmount = 0.0f;
		//	if(defAmount > 1.0f) defAmount = 1.0f;

		//	tempColor[i] = defAmount;

		//	//defAmount= defAmount;
		//	//if(defAmount < 0.6f){
		//	//	tempColor[i] = 0.15f;	//青
		//	//}
		//	//else{
		//	//	tempColor[i] = 0.95f;	//赤
		//	//	defNum++;
		//	//}
		//}

		//粒子数と変形量の組み合わせ　ギリギリ及第点
		for(int i = 0; i < m_pPS->GetNumParticles(); i++){
			int pNumAve = 0;

			//クラスタに含まれる粒子数
			//for(int j = 0; j < m_iceObj->GetOrientedObj(i)->GetNumVertices(); j++){
			//	int pIndx = m_iceObj->GetOrientedObj(i)->GetParticleIndx(j);
			//	int particleNum = m_iceObj->GetOrientedObj(pIndx)->GetNumVertices();
			//	pNumAve += particleNum;
			//}

			//pNumAve /= m_iceObj->GetOrientedObj(i)->GetNumVertices();

			//閾値の代わりに変形量の標準偏差を用いて判定
			float defAmount = m_iceObj->GetOrientedObj(i)->FinalDefAmount() * 1000.0f;
			if(defAmount < 0.0f) defAmount = 0.0f;
			if(defAmount > 1.0f) defAmount = 1.0f;
		
			defAmount= defAmount;

			float def2 = m_iceObj->GetOrientedObj(i)->DefAmount();

			if(defAmount < threshold/* || pNumAve > 20*/){
			//if(defAmount < def2){
				tempColor[i] = 0.15f;	//青
			}
			else{
				tempColor[i] = 0.83f;	//赤
				defNum++;
			}
		}

		cout << defNum << endl;

		//最終位置と初期位置の差分の平均
		//for(int i = 0; i < m_pPS->GetNumParticles(); i++){
		//	float defAmount = m_iceObj->GetOrientedObj(i)->FinalDefAmount() * 100.0f;
		//	if(defAmount < 0.0f) defAmount = 0.0f;
		//	if(defAmount > 1.0f) defAmount = 1.0f;

		//	//tempColor[i] = defAmount;

		//	defAmount= defAmount;
		//	if(defAmount < threshold || m_iceObj->GetMoveObj(i)->GetNumVertices() > 17){
		//		tempColor[i] = 0.15f;	//青
		//	}
		//	else{
		//		tempColor[i] = 0.95f;	//赤
		//		defNum++;
		//	}
		//}

		//cout << defNum << endl;

		//角速度の分布
		//for(int i = 0; i < m_pPS->GetNumParticles(); i++){
		//	float vel = norm(m_iceObj->GetOrientedObj(i)->Particle(0)->AngularVel());
		//	if(vel < 0.0f) vel = 0.0f;
		//	tempColor[i] = vel/5.0f;
		//}

		//速度の分布
		//for(int i = 0; i < m_pPS->GetNumParticles(); i++){
		//	float vel = norm(m_iceObj->GetOrientedObj(i)->GetVertexVel(0));
		//	if(vel < 0.0f) vel = 0.0f;
		//	tempColor[i] = vel/5.0f;
		//}

		////変形量がどんな分布になっているか
		//for(int i = 0; i < m_pPS->GetNumParticles(); i++){
		//	float defAmount = m_iceObj->GetOrientedObj(i)->DefAmount();
		//	if(defAmount < 0.0f) defAmount = 0.0f;
		//	tempColor[i] = defAmount;
		//}

		//方向ベクトルの内積総和　quaternionの場合と全然変わらない
		//for(int i = 0; i < m_pPS->GetNumParticles(); i++){
		//	float dotSum = m_iceObj->GetOrientedObj(i)->QuatInnerDotSum();
		//	if(dotSum < 0.0f) dotSum = 0.0f;
		//	if(dotSum > 1.0f) dotSum = 1.0f;
		//	tempColor[i] = dotSum;
		//}

		//quaternionの内積総和がどんな分布になっているか
		//全然色変わらない
		//for(int i = 0; i < m_pPS->GetNumParticles(); i++){
		//	float quatDef = m_iceObj->GetOrientedObj(i)->QuatInnerDotSum();
		//	if(quatDef < 0.0f) quatDef = 0.0f;
		//	if(quatDef > 1.0f) quatDef = 1.0f;
		//	tempColor[i] = quatDef * quatDef * quatDef * quatDef * quatDef * quatDef;
		//}

		//quaternionの内積総和を閾値とした場合
		//for(int i = 0; i < m_pPS->GetNumParticles(); i++){
		//	float quatDef = m_iceObj->GetOrientedObj(i)->QuatInnerDotSum();
		//	if(quatDef < threshold){
		//		tempColor[i] = 0.95f;	//赤
		//		cout << quatDef << endl;
		//	}
		//	else{
		//		tempColor[i] = 0.15f;	//青
		//	}
		//}

		////変形量を閾値とした場合
		//for( int i = 0; i < m_pPS->GetNumParticles(); i++ ){
		//	float defAmount = m_iceObj->GetOrientedObj(i)->DefAmount();
		//	if(defAmount > threshold){				
		//		tempColor[i] = 0.95f;	//赤
		//	}
		//	else{
		//		tempColor[i] = 0.25f;	//青
		//	}
		//}

		m_pPS->SetColorVBOFromArray( tempColor, 1, false, 1.0f );

		delete[] tempColor;
	}
}

//現在の粒子の位置と速度をファイルにダンプ
void rxFlGLWindow::DumpParticleData()
{	cout << __FUNCTION__ << endl;

	ofstream ofs( "DampParticleData.txt" );
	ofs << "DampParticleData" << endl;

	RXREAL *p = m_pPS->GetArrayVBO(rxParticleSystemBase::RX_POSITION);
	RXREAL *v = m_pPS->GetArrayVBO(rxParticleSystemBase::RX_VELOCITY);

	//位置
	ofs << "Pos::" << endl;

	for(int i = 0; i < m_iIcePrtNum; i++)
	{
		ofs << i << ", (" << p[i*4+0] << ", " << p[i*4+1] << ", " << p[i*4+2] << ")" << endl;
	}

	//速度
	ofs << "Vel::" << endl;

	for(int i = 0; i < m_iIcePrtNum; i++)
	{
		ofs << i << ", (" << v[i*4+0] << ", " << v[i*4+1] << ", " << v[i*4+2] << ")" << endl;
	}

}

//
/*!
 * 粒子情報から四面体を作成
 * @param[in] pList 四面体を作成する粒子リスト
 * @param[in] tList 作成した四面体の番号リスト
 */
void rxFlGLWindow::MakeFreezeTetrahedra(vector<int>& pList, vector<int>& tList)
{
	RXREAL *p = m_pPS->GetArrayVBO(rxParticleSystemBase::RX_POSITION);
	vector<vector<rxNeigh>>& neights = ((RXSPH*)m_pPS)->GetNeights();		//近傍粒子を取得

	//凝固粒子毎に近傍固体粒子を探索，四面体を作成
	for(unsigned i = 0; i < pList.size(); i++)
	{	
		int ipIndx = pList[i];
		vector<int> npList;													//近傍固体粒子

		//cout << __FUNCTION__ << " check1 i = " << i << " ipIndx = " << ipIndx << endl;

		//凝固粒子の近傍固体粒子を探索
		for(unsigned j = 0; j < neights[ipIndx].size(); j++)
		{
			int jpIndx = neights[ipIndx][j].Idx;							//近傍粒子の番号
			float distance	= neights[ipIndx][j].Dist2;
			float radius	= pow(((RXSPH*)m_pPS)->GetEffectiveRadius()*0.95f, 2); 

			//cout << "distance = " << distance << " radius = " << radius << endl;

			if(ipIndx == jpIndx )								continue;	//自分自身は除く
			//if(m_ht->getTemps()[jpIndx] > 250)					continue;	//粒子の温度が一定以下
			if(m_ice->GetParticleNum() <= jpIndx)			{	continue;	}	//融解のみの実験のときに必要になる．
			if(m_ice->GetPtoCNum(jpIndx) <= 0)					continue;	//クラスタに所属していないなら戻る
			if(m_ice->GetPtoTNum(jpIndx) <= 0)					continue;	//四面体に属していないなら戻る
			if(distance > radius)								continue;	//凝固粒子との距離が影響半径の半分

			npList.push_back(jpIndx);										//粒子番号を登録
		}

		//３　近傍氷塊粒子数で分岐　終了か最探索
		int iParticleSize = npList.size();
		if(iParticleSize == 0)
		{
			//cout << __FUNCTION__ << " check2a iParticleSize = " << iParticleSize << endl;
			continue;
		}
		else if(1 <= iParticleSize && iParticleSize <= 2)
		{
			//最近傍氷塊粒子が３個になるまで探索　近傍クラスタ0層目まで探索
			//cout << __FUNCTION__ << " check2b iParticleSize = " << iParticleSize << endl;

			for(int j = 0; j < iParticleSize; j++)
			{
				int jpIndx = npList[j];
				vector<int> jpList = GetNeighborDistanceIceParticles(jpIndx, 0, 20);			//最大20個取得

				for(unsigned k = 0; k < jpList.size(); k++)
				{
					int kpIndx = jpList[k];
					if(std::find(npList.begin(), npList.end(), kpIndx) != npList.end()){	continue;	}
					npList.push_back(kpIndx);

					if(npList.size() == 3){	break;	}
				}
				if(npList.size() == 3){	break;	}
			}
		}
		else if(iParticleSize >= 3)
		{
			//cout << __FUNCTION__ << " check2c iParticleSize = " << iParticleSize << endl;
		}

		if(npList.size() == 0) continue;

		//四面体作成　相転移粒子＋氷塊粒子
		npList.push_back(ipIndx);
		iParticleSize = npList.size();

		if(iParticleSize >= 4)
		{
			//cout << __FUNCTION__ << " check2a iParticleSize = " << iParticleSize << endl;
			//試作　複数の四面体を考えない
			while(npList.size() > 4)
			{	npList.erase(npList.begin());	}

			tList.push_back(m_iTetraNum);

			MakeTetraInfo(m_iTetraNum, npList);								//四面体作成，各種情報登録

			//デバッグ
			//cout << __FUNCTION__ << " Debug ipIndx = " << ipIndx << " npList.size() = " << npList.size() << endl;
			//for( unsigned j = 0; j < npList.size(); j++ )
			//{
			//	cout << " " << npList[j];
			//}
			//cout << endl;

			m_fIntrps[ipIndx] = 1.0f;										//線形補間もしない
			//m_ht->setPhaseChange(ipIndx, 0);								//相転移し終わったことを伝える
		}
		else
		{
			//cout << __FUNCTION__ << " check2b iParticleSize = " << iParticleSize << endl;
		}
	}//for(unsigned i = 0; i < pList.size(); i++)

}

/*!
 * 凝固粒子のみで四面体を作成
 * @param[in] pList 四面体を作成する粒子リスト
 * @param[in] tList 追加された四面体の番号リスト
 */
void rxFlGLWindow::MakeFreezeTetrahedra_OnlyFreezeParticle(const vector<int>& pList, vector<int>& tList)
{

}

/*!
 * 凝固粒子を四面体に追加
 * @param[in] pList 四面体を作成する粒子リスト
 * @param[in] tList 追加された四面体の番号リスト
 */
 void rxFlGLWindow::AddFreezeTetrahedra(const vector<int>& pList, vector<int>& tList)
{

}


/*!
 * FPSを計算してウィンドウタイトルに表示
 */
void rxFlGLWindow::ComputeFPS(void)
{
	static unsigned int frame_count = 0;
	static bool verify = false;
	static int fps_count = 0;
	static int fps_limit = 1;

	frame_count++;
	fps_count++;
	if(fps_count == fps_limit-1){
		verify = true;
	}
	if(fps_count == fps_limit){
		char fps[256];
	
		g_iTimeCount++;

		//int startFrame = 299;

		////300フレームから計測
		////if(g_iTimeCount > startFrame)
		{
			g_fTotalTime += g_TimerFPS.GetTime(0);
		}

		g_fAvgTime = g_fTotalTime/(float)g_iTimeCount;

		////平均計算時間を表示
		//if(g_iTimeCount%100 == 0)
		//{
		//	cout << "count = " << g_iTimeCount << " AvgTime :: " << g_fAvgTime << endl;
		//}

		int num_particles = m_pPS->GetNumParticles();
		double ifps = 1.0/g_fAvgTime;
		sprintf(fps, "CUDA Particles (%d particles): %d step, %3.1f fps", num_particles, m_iCurrentStep, ifps);  

		//glutSetWindowTitle(fps);
		m_pParent->SetStatusLabel(fps);
		fps_count = 0; 
	
		g_TimerFPS.Reset();
	}
}

/*!
 * "\n"が含まれるstringを複数のstringに分解する
 * @param[in] org 元の文字列
 * @param[in] div 分解後の文字列配列
 */
void rxFlGLWindow::DivideString(const string &org, vector<string> &div)
{
	size_t pos0 = 0, pos1 = 0;
	while(pos1 != string::npos){
		pos1 = org.find("\n", pos0);

		div.push_back(org.substr(pos0, pos1-pos0));

		pos0 = pos1+1;
	}
}

/*!
 * 画面出力用の文字列の生成
 * @return 文字列
 */
vector<string> rxFlGLWindow::SetDrawString(void)
{
	vector<string> str;

	str.push_back("Color : ");
	int type = m_pPS->GetColorType();
	if(type == rxParticleSystemBase::RX_DENSITY) str.back() << STR(RX_DENSITY);
	else if(type == rxParticleSystemBase::RX_PRESSURE) str.back() << STR(RX_PRESSURE);
	else if(type == rxParticleSystemBase::RX_ENERGY_SPECTRUM) str.back() << STR(RX_ENERGY_SPECTRUM);
	else if(type == rxParticleSystemBase::RX_RAMP) str.back() << STR(RX_RAMP);
	else str.back() << STR(RX_NONE);

	str.back() << ", Drawing : ";
	if((m_iDraw & RXD_REFRAC)) str.back() << "refraction mesh";
	else if(m_iDrawPS == RXP_POINT) str.back() << STR(RXP_POINT);
	else if(m_iDrawPS == RXP_POINTSPRITE) str.back() << STR(RXP_POINTSPRITE);
	else if(m_iDrawPS == RXP_POINT_UPSAMPLE) str.back() << STR(RXP_POINT_UPSAMPLE);
	else if(m_iDrawPS == RXP_POINT_NONE) str.back() << STR(RXP_POINT_NONE);

	str.push_back("");
#if defined(RX_USE_GPU)
	#if defined(RX_USE_PBD)
		str.back() << "PBDSPH(GPU)";
	#else
		str.back() << "SPH(GPU)";
	#endif
#else
	#if defined(RX_USE_PBD)
		str.back() << "PBDSPH";
	#elif defined(RX_USE_DD)
		str.back() << "SPH(DD)";
	#else
		str.back() << "SPH(CPU)";
	#endif
#endif

	if(m_bsSimuSetting.at(ID_SPH_MESH)){
		str.back() << " - ";
		if(m_iTriangulationMethod == RXM_SSM_CPU) str.back() << "SSM";
		else if(m_iTriangulationMethod == RXM_SSM_GPU) str.back() << "SSM(GPU)";
		else if(m_iTriangulationMethod == RXM_MC_GPU) str.back() << "MC(GPU)";
		else str.back() << "MC";
	}
	str.back() << ", WT " << (m_pPS->IsWaveletTurb() ? "on" : "off");
	str.back() << ", SPS " << (m_pPS->IsSubParticle() ? "on" : "off");
	str.back() << ", VC " << (m_pPS->IsUseVorticity() ? "on" : "off");

	if((m_iSaveImageSpacing > 0) || m_bsSimuSetting.at(ID_SPH_OUTPUT) || m_bsSimuSetting.at(ID_SPH_MESH_OUTPUT)){
		str.push_back("");
		str.back() << "Output : " << (m_iSaveImageSpacing > 0 ? "img " : "");
		str.back() << (m_bsSimuSetting.at(ID_SPH_OUTPUT) ? "prts " : "");
		str.back() << (m_bsSimuSetting.at(ID_SPH_MESH_OUTPUT) ? "mesh " : "");
		if(m_bsSimuSetting.at(ID_SPH_INPUT)){
			str.back() << ", Input : " << (m_bsSimuSetting.at(ID_SPH_INPUT) ? "prts " : "");
		}
	}

	str.push_back("num particles : ");
	str.back() << m_pPS->GetNumParticles();
	if(m_pPS->GetNumBoundaryParticles()){
		str.back() << " (boundary : " << m_pPS->GetNumBoundaryParticles() << ")";
	}

	//str.push_back("");
	//str.back() << "wavelet scale : ";
	//str.back() << RXFRMT(g_fWaveletScale);
	//str.back() << ", coef. et : ";
	//str.back() << RXFRMTE(g_fCoefEt);
#ifdef RX_USE_PBD
	// 反復
	str.push_back("");
	str.back() << "Iterations : " << g_iIterations;
	str.push_back("");
	str.back() << "Eta : " << g_fEta;
	str.push_back("Artificial pressure [t] : ");
	str.back() << (static_cast<RXSPH*>(m_pPS)->GetArtificialPressure() ? "on" : "off");
#endif

	if(m_bsSimuSetting[ID_SPH_SPS_TURB]){
		str.push_back("valid sub-particles : ");
		str.back() << ((RXSPH*)m_pPS)->GetNumValidSubParticles() << " / " << ((RXSPH*)m_pPS)->GetNumAllSubParticles();
	}

	// メッシュ
	if(m_bsSimuSetting.at(ID_SPH_MESH)){
		if(m_iTriangulationMethod == RXM_SSM_CPU || m_iTriangulationMethod == RXM_SSM_GPU){
			str.push_back("SSM : ");
			str.back() << "grid (" << m_iWinW/m_fSpacing << "x" << m_iWinH/m_fSpacing << "), ";
		}
		else{
			str.push_back("MC : ");
			str.back() << "grid (" << m_iMeshN[0] << "x" << m_iMeshN[1] << "x" << m_iMeshN[2] << "), ";
		}
		str.back() << "vertex " << m_iNumVrts << ", face " << m_iNumTris;
		str.push_back("MC(solid) : ");
		str.back() << "vertex " << m_iNumVrts_solid << ", face " << m_iNumTris_solid;
	}

	// 計測された計算時間
	str.push_back("time : ");
	string tstr;
	RXTIMER_STRING(tstr);
	DivideString(tstr, str);

	return str;
}


/*!
 * パーティクル速度をGL_LINESで描画
 * @param[in] prts パーティクル位置
 * @param[in] vels パーティクル速度
 * @param[in] n パーティクル数
 * @param[in] d 配列のステップ
 * @param[in] len 線の長さ
 */
void rxFlGLWindow::DrawParticleVector(RXREAL *prts, RXREAL *vels, int n, int d, double *c0, double *c1, double len)
{
	//RXCOUT << "GL_POINTS" << endl;
	glBegin(GL_LINES);
	int k = 0;
	for(int i = 0; i < n; ++i){
		glColor3dv(c0);
		glVertex3d(prts[k], prts[k+1], prts[k+2]);
		glColor3dv(c1);
		glVertex3d(prts[k]+len*vels[k], prts[k+1]+len*vels[k+1], prts[k+2]+len*vels[k+2]);
		k += d;
	}
	glEnd();
}


//追加
/*!
 * パーティクル速度をGL_LINESで描画
 * @param[in] prts パーティクル位置
 * @param[in] vels パーティクル速度
 * @param[in] n パーティクル数
 * @param[in] d 配列のステップ
 * @param[in] len 線の長さ
 */
void rxFlGLWindow::DrawFastPath(RXREAL *prts, RXREAL *vels, int n, int d, double *c0, double *c1, double len)
{//	RXCOUT << __FUNCTION__ << endl;
	glBegin(GL_LINE_STRIP);
		int k = 0;
		glColor3dv(c0);

		for(int i = 0; i < n; ++i)
		{
			glVertex3d(prts[k], prts[k+1], prts[k+2]);
			k += d;
		}
	glEnd();
}


/*!
 * パーティクルをGL_POINTSで描画
 *  - VBOがあれば用いる
 * @param[in] vbo パーティクル位置を格納したVBO
 * @param[in] n パーティクル数
 * @param[in] color_vbo パーティクルの色を格納したVBO
 * @param[in] data パーティクル座標(ホスト側メモリ，vboが0の時に使用)
 */
void rxFlGLWindow::DrawParticlePoints(unsigned int vbo, int n, unsigned int color_vbo, RXREAL *data)
{
	if(!vbo){
		glBegin(GL_POINTS);
		int k = 0;
		for(int i = 0; i < n; ++i){
			glVertex3d(data[k], data[k+1], data[k+2]);
			k += 4;
		}
		glEnd();
	}
	else{
		glBindBuffer(GL_ARRAY_BUFFER, vbo);
		glEnableClientState(GL_VERTEX_ARRAY);
		glVertexPointer(4, GL_REAL, 0, 0);

		if(color_vbo){
			glBindBufferARB(GL_ARRAY_BUFFER_ARB, color_vbo);
			glColorPointer(4, GL_REAL, 0, 0);
			glEnableClientState(GL_COLOR_ARRAY);
		}

		glDrawArrays(GL_POINTS, 0, n);

		glDisableClientState(GL_VERTEX_ARRAY); 
		glDisableClientState(GL_COLOR_ARRAY); 
		glBindBuffer(GL_ARRAY_BUFFER, 0);
	}
}


/*!
 * サブパーティクルの描画
 * @param[in] pVBO 位置のVBO
 * @param[in] uColorVBO 色のVBO
 * @param[in] data 位置(CPU)
 */
void rxFlGLWindow::DrawSubParticles(void)
{
	// Etに基づいて各サブパーティクルの混合比率と半径を計算して詰める
	((RXSPH*)m_pPS)->CalRadiusAndRatio();

	// サブパーティクル数
	int num_all_sp = ((RXSPH*)m_pPS)->GetNumAllSubParticles();
	int num_sp = ((RXSPH*)m_pPS)->GetNumValidSubParticles();
	int num_p = m_pPS->GetNumParticles();

	if(!num_sp) return;

	// サブパーティクルの位置と半径(有効な物だけ詰めた配列)
	RXREAL *sppos = ((RXSPH*)m_pPS)->GetSubParticlePos();
	RXREAL *sprad = ((RXSPH*)m_pPS)->GetSubParticleRad();

	// パーティクルの大きさ
	RXREAL prad_level[MAX_SUB_LEVEL+1];
	for(int l= 0; l <= MAX_SUB_LEVEL; ++l){
		prad_level[l] = ((RXSPH*)m_pPS)->GetSubRadius(l);
	}

	//cout << num_sp << "/" << num_all_sp << " - " << num_p << endl;
	RXREAL rad = prad_level[0];
	glUniform1f( glGetUniformLocation(g_glslPointSprite.Prog, "pointRadius"), rad);
	int c = 0;
	int level = 0;
	for(int i = 0; i < num_sp; ++i){
		while(fabs(rad-sprad[i]) > 0.01*rad){
			level++;
			rad = prad_level[level];
		}

		glUniform1f( glGetUniformLocation(g_glslPointSprite.Prog, "pointRadius"), sprad[i]);
		Vec3 col(0.0);
		col[(level <= 2 ? level : 2)] = 1.0;
		glColor3dv(col.data);

		glBegin(GL_POINTS);
		glVertex3f(sppos[4*i+0], sppos[4*i+1], sppos[4*i+2]);
		glEnd();
	}
}


/*!
 * Vertex Buffer Object(VBO)の作成
 * @param[inout] vbo バッファID
 * @param[in] size バッファサイズ
 */
void rxFlGLWindow::CreateVBO(GLuint* vbo, unsigned int size)
{
	glGenBuffers(1, vbo);
	glBindBuffer(GL_ARRAY_BUFFER, *vbo);

	// バッファオブジェクトの初期化
	glBufferData(GL_ARRAY_BUFFER, size, 0, GL_DYNAMIC_COPY);
	glBindBuffer(GL_ARRAY_BUFFER, 0);

	//CuRegisterGLBufferObject(*vbo);
}

/*!
 * Vertex Buffer Object(VBO)の削除
 * @param[inout] vbo バッファID
 */
void rxFlGLWindow::DeleteVBO(GLuint* vbo)
{
	//CuUnregisterGLBufferObject(*vbo);

	glBindBuffer(1, *vbo);
	glDeleteBuffers(1, vbo);

	*vbo = 0;
}


/*!
 * VBOを用いた等値面ポリゴン描画
 */
void rxFlGLWindow::DrawLiquidSurface(void)
{
	if(m_uVrtVBO){
		glEnableClientState(GL_VERTEX_ARRAY);
		glEnableClientState(GL_NORMAL_ARRAY);

		glBindBuffer(GL_ARRAY_BUFFER, m_uVrtVBO);
		if(m_iTriangulationMethod == RXM_SSM_CPU || m_iTriangulationMethod == RXM_SSM_GPU){
			glVertexPointer(3, GL_FLOAT, 0, 0);
		}
		else{
			glVertexPointer(m_iDimVBO, GL_FLOAT, 0, 0);
		}
		glBindBuffer(GL_ARRAY_BUFFER, m_uNrmVBO);
		glNormalPointer(GL_FLOAT, 0, 0);
	 
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_uTriVBO);
		glDrawElements(GL_TRIANGLES, m_iNumTris*3, GL_UNSIGNED_INT, 0);
		 
		glBindBuffer(GL_ARRAY_BUFFER, 0);
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

		glDisableClientState(GL_VERTEX_ARRAY); 
		glDisableClientState(GL_NORMAL_ARRAY); 
	}
	else{
		if(m_Poly.vertices.empty()) return;
		int ntris = (int)m_Poly.faces.size();
		for(int i = 0; i < ntris; ++i){
			// 面
			glColor4d(0.0, 0.0, 1.0, 1.0);
			glBegin(GL_POLYGON);
			for(int j = 0; j < 3; ++j){
				Vec3 pos = m_Poly.vertices[m_Poly.faces[i][j]];
				Vec3 nrm = m_Poly.normals[m_Poly.faces[i][j]];
				glNormal3dv(nrm.data);
				glVertex3dv(pos.data);
			}
			glEnd();
		}
	}
}

//追加：氷用
/*
 * VBOを用いた等値面ポリゴン描画(固体)
 */
void rxFlGLWindow::DrawSolidSurface(void)
{
	if(m_uVrtVBO_solid){
		glEnableClientState(GL_VERTEX_ARRAY);
		glEnableClientState(GL_NORMAL_ARRAY);

		glBindBuffer(GL_ARRAY_BUFFER, m_uVrtVBO_solid);
		if(m_iTriangulationMethod == RXM_SSM_CPU || m_iTriangulationMethod == RXM_SSM_GPU){
			glVertexPointer(3, GL_FLOAT, 0, 0);
		}
		else{
			glVertexPointer(m_iDimVBO, GL_FLOAT, 0, 0);
		}
		glBindBuffer(GL_ARRAY_BUFFER, m_uNrmVBO_solid);
		glNormalPointer(GL_FLOAT, 0, 0);
	 
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_uTriVBO_solid);
		glDrawElements(GL_TRIANGLES, m_iNumTris_solid*3, GL_UNSIGNED_INT, 0);
		 
		glBindBuffer(GL_ARRAY_BUFFER, 0);
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

		glDisableClientState(GL_VERTEX_ARRAY); 
		glDisableClientState(GL_NORMAL_ARRAY); 
	}
}


void rxFlGLWindow::SetParticleColorType(int type, int change)
{
	m_pPS->SetColorType(type);
	m_pPS->SetColorVBO();
	m_iColorType = type;
}

/*!
 * シーン描画
 */
void rxFlGLWindow::RenderSphScene(void)
{
	// MARK:RenderScene
	glEnable(GL_LIGHTING);
	glEnable(GL_DEPTH_TEST);
	glDisable(GL_CULL_FACE);
	glDepthFunc(GL_LESS);

	int pnum = m_pPS->GetNumParticles();
	Vec3 cen = m_pPS->GetCen();
	Vec3 dim = m_pPS->GetDim();

	//
	// 周囲境界
	//
	if((m_iDraw & RXD_BBOX)){
		glDisable(GL_LIGHTING);
		glColor4d(0.0, 0.0, 0.0, 1.0);
//		glColor4d(1.0, 1.0, 1.0, 1.0);
		glLineWidth(1.0);
		glPushMatrix();
		glTranslatef(cen[0], cen[1], cen[2]);
		glScalef(dim[0], dim[1], dim[2]);
		glutWireCube(1.0);
		glPopMatrix();
	}

	// 軸
	if((m_iDraw & RXD_AXIS)){
		Vec3 len = 0.6*dim;
		glLineWidth((GLfloat)3.0);

		// x軸
		glColor3f(1.0, 0.0, 0.0);
		glBegin(GL_LINES);
		glVertex3d(0.0, 0.0, 0.0);
		glVertex3d(len[0], 0.0, 0.0);
		glEnd();

		// y軸
		glColor3f(0.0, 1.0, 0.0);
		glBegin(GL_LINES);
		glVertex3d(0.0, 0.0, 0.0);
		glVertex3d(0.0, len[1], 0.0);
		glEnd();

		// z軸
		glColor3f(0.0, 0.0, 1.0);
		glBegin(GL_LINES);
		glVertex3d(0.0, 0.0, 0.0);
		glVertex3d(0.0, 0.0, len[2]);
		glEnd();
	}


	//
	// 障害物
	//
	if((m_iDraw & RXD_SOLID)){
		glDisable(GL_LIGHTING);
		glColor4d(0.0, 1.0, 0.0, 1.0);
		glLineWidth(1.0);

		m_pPS->DrawObstacles();

		vector<rxPolygons>& solid_poly = m_Scene.GetSolidPolys();
		
		vector<rxPolygons>::iterator itr = solid_poly.begin();
		for(; itr != solid_poly.end(); ++itr){
			if(itr->open){
				// 描画フラグ : 下位ビットから頂点,エッジ,面,法線 - 1,2,4,8)
				itr->Draw(m_iSolidDraw);
			}
		}
	}

	// 
	// パーティクル
	// 
	if(!m_bsSimuSetting.at(ID_SPH_ANISOTROPIC) && (m_iDraw & RXD_PARTICLE) && !(m_iDraw & RXD_REFRAC)){
		glDisable(GL_LIGHTING);
		glColor4d(0.0, 0.0, 1.0, 1.0);
		glEnable(GL_POINT_SMOOTH);
		glPointSize(1.0);

		unsigned int pvbo = m_pPS->GetCurrentReadBuffer();
		unsigned int cvbo = m_pPS->GetColorBuffer();

		if(m_iDrawPS == RXP_POINT){				// GL_POINTSで描画
			RXREAL *data = 0;
			data = m_pPS->GetArrayVBO(rxParticleSystemBase::RX_POSITION);
			DrawParticlePoints(0, pnum, cvbo, data);
		}
		else if(m_iDrawPS == RXP_POINTSPRITE){	// ポイントスプライトで描画
			float prad = (float)m_pPS->GetParticleRadius();			// パーティクル半径
			float pscale = m_iWinH/tanf(RX_FOV*0.5f*(float)M_PI/180.0f);

			glEnable(GL_POINT_SPRITE_ARB);
			glTexEnvi(GL_POINT_SPRITE_ARB, GL_COORD_REPLACE_ARB, GL_TRUE);
			glEnable(GL_VERTEX_PROGRAM_POINT_SIZE_NV);
			glDepthMask(GL_TRUE);
			glEnable(GL_DEPTH_TEST);

			glUseProgram(g_glslPointSprite.Prog);
			glUniform1f( glGetUniformLocation(g_glslPointSprite.Prog, "pointScale"), pscale );
			glUniform1f( glGetUniformLocation(g_glslPointSprite.Prog, "pointRadius"), prad );

			if(m_pPS->IsSubParticle()){
				DrawSubParticles();
			}
			else{
				DrawParticlePoints(pvbo, pnum, cvbo, 0);
			}

			glUseProgram(0);
			glDisable(GL_POINT_SPRITE_ARB);
		}

		if(m_iPickedParticle != -1){
			RXREAL prad = m_pPS->GetParticleRadius();
			int k = m_iPickedParticle;
			RXREAL *data = m_pPS->GetArrayVBO(rxParticleSystemBase::RX_POSITION);
			Vec3 pos(data[DIM*k+0], data[DIM*k+1], data[DIM*k+2]);
			glDisable(GL_LIGHTING);
			glColor3d(1.0, 1.0, 0.0);
			glPushMatrix();
			glTranslatef(pos[0], pos[1], pos[2]);
			glutSolidSphere(1.5*prad, 8, 4);
			glPopMatrix();


		}
	}

	//
	// 境界パーティクル
	//
	if(m_iDraw & RXD_BPARTICLE){
		glDisable(GL_LIGHTING);
		glColor4d(0.0, 1.0, 0.0, 1.0);
		glEnable(GL_POINT_SMOOTH);
		glPointSize(1.0);

		int bnum = m_pPS->GetNumBoundaryParticles();

		RXREAL *data = m_pPS->GetArrayVBO(rxParticleSystemBase::RX_BOUNDARY_PARTICLE, true, bnum);

		if(m_pPS->m_vFuncB.empty()){
			DrawParticlePoints(0, bnum, 0, data);
		}
		else{
			glBegin(GL_POINTS);
			int k = 0;
			for(int i = 0; i < bnum; ++i){
				Vec3 col, col0(0, 1, 0), col1(0, 0, 0);
				RXFunc::Gradation(col.data, col0.data, col1.data, m_pPS->m_vFuncB[i]);
				glColor3dv(col.data);
				glVertex3d(data[k], data[k+1], data[k+2]);
				k += DIM;
			}
			glEnd();
		}
	}


	//
	// Anisotropic Kernel
	//
	if(m_bsSimuSetting.at(ID_SPH_ANISOTROPIC) && (m_iDraw & RXD_PARTICLE) && !(m_iDraw & RXD_REFRAC)){
		glDisable(GL_LIGHTING);

		RXREAL *ppos = m_pPS->GetArrayVBO(rxParticleSystemBase::RX_POSITION);

		glEnable(GL_LIGHTING);

		GLfloat mat_diff[] = { 0.1f, 0.1f, 1.0f, 1.0f };
		GLfloat mat_spec[] = { 0.2f, 0.2f, 0.2f, 1.0f };
		GLfloat mat_ambi[] = { 0.1f, 0.1f, 0.1f, 1.0f };

		glColor3d(0.0, 0.0, 1.0);
		glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE,  mat_diff);
		glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, mat_spec);
		glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT,  mat_ambi);
		glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 100);

		RXREAL *G = m_pPS->GetArrayVBO(rxParticleSystemBase::RX_TRANSFORMATION);
		RXREAL trans[9];
		GLfloat gltrans[16];
		for(int i = 0; i < pnum; ++i){
#ifdef RX_USE_GPU
			CalInverse3x3(&G[9*i], trans);
#else
			for(int j = 0; j < 9; ++j) trans[j] = G[9*i+j];
#endif
			GetGLMatrix(trans, gltrans);

			glPushMatrix();
			glTranslatef(ppos[DIM*i+0], ppos[DIM*i+1], ppos[DIM*i+2]);
			glMultMatrixf(gltrans);
			glutSolidSphere(0.01, 8, 4);
			glPopMatrix();
		}
	}

	//
	// メッシュ描画
	//
	if((m_iDraw & RXD_MESH)){
		if((m_iDraw & RXD_REFRAC) && m_bUseCubeMap){
			glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
			glEnable(GL_LIGHTING);

			glEnable(GL_COLOR_MATERIAL);

			Vec3 eye_pos(0.0);		// 視点位置
			m_tbView.CalLocalPos(eye_pos.data, Vec3(0.0).data);
			
			//ここをコメントアウトすれば，環境マッピングなどは適用されない
			glUseProgram(g_glslFresnel.Prog);

			// パラメータ設定
			// バーテックスシェーダ用パラメータ
			//glUniform1f(glGetUniformLocation(g_glslFresnel.Prog, "etaRatio"), 0.53);
			//glUniform1f(glGetUniformLocation(g_glslFresnel.Prog, "fresnelBias"), 0.005);
			//glUniform1f(glGetUniformLocation(g_glslFresnel.Prog, "fresnelPower"), 0.98);
			//glUniform1f(glGetUniformLocation(g_glslFresnel.Prog, "fresnelScale"), 1.0);
			//glUniform3f(glGetUniformLocation(g_glslFresnel.Prog, "eyePosition"), eye_pos[0], eye_pos[1], eye_pos[2]);

			//glUniform1f(glGetUniformLocation(g_glslFresnel.Prog, "etaRatio"), 0.53);
			//glUniform1f(glGetUniformLocation(g_glslFresnel.Prog, "fresnelBias"), 0.85);
			//glUniform1f(glGetUniformLocation(g_glslFresnel.Prog, "fresnelPower"), 0.11);
			//glUniform1f(glGetUniformLocation(g_glslFresnel.Prog, "fresnelScale"), 0.5);
			//glUniform3f(glGetUniformLocation(g_glslFresnel.Prog, "eyePosition"), eye_pos[0], eye_pos[1], eye_pos[2]);

			glUniform1f(glGetUniformLocation(g_glslFresnel.Prog, "etaRatio"),		m_etaRatio);
			glUniform1f(glGetUniformLocation(g_glslFresnel.Prog, "fresnelBias"),	m_fresnelBias);
			glUniform1f(glGetUniformLocation(g_glslFresnel.Prog, "fresnelPower"),	m_fresnelPower);
			glUniform1f(glGetUniformLocation(g_glslFresnel.Prog, "fresnelScale"),	m_fresnelScale);
			glUniform3f(glGetUniformLocation(g_glslFresnel.Prog, "eyePosition"), eye_pos[0], eye_pos[1], eye_pos[2]);

			// フラグメントシェーダ用パラメータ
			glUniform1i(glGetUniformLocation(g_glslFresnel.Prog, "envmap"), 0);

			//氷と水で色を変える
			//glColor3d(0.0, 0.6, 1.0);
			DrawLiquidSurface();

//			glColor3d(0.6, 0.6, 1.0);
//			DrawSolidSurface();

			glUseProgram(0);
		}
		else{
			glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
			glDisable(GL_LIGHTING);
			glEnable(GL_CULL_FACE);
//			glColor4d(0.3, 0.3, 1.0, 1.0);
			//DrawLiquidSurface();
//			glColor4d(0.3, 1.0, 0.3, 1.0);
//			DrawSolidSurface();
		}
	}

	//
	// パーティクル速度
	//
	if((m_iDraw & RXD_VELOCITY)){
		RXREAL *p = m_pPS->GetArrayVBO(rxParticleSystemBase::RX_POSITION);
		RXREAL *v = m_pPS->GetArrayVBO(rxParticleSystemBase::RX_VELOCITY);
		if(p && v){
			glDisable(GL_LIGHTING);
			glLineWidth(1.0);
			glColor3d(0.0, 1.0, 1.0);
			DrawParticleVector(p, v, pnum, DIM, Vec3(0.8, 0.8, 1.0).data, Vec3(0.0, 0.0, 1.0).data, norm(dim)*m_fVScale);
		}
	}

	//
	// パーティクル乱流速度
	//
	if((m_iDraw & RXD_TURB_VELOC)){
		RXREAL *p = m_pPS->GetArrayVBO(rxParticleSystemBase::RX_POSITION);
		RXREAL *v = m_pPS->GetArrayVBO(rxParticleSystemBase::RX_TURB_VELOCITY, true);
		if(p && v){
			glDisable(GL_LIGHTING);
			glLineWidth(1.0);
			glColor3d(0.0, 1.0, 1.0);
			DrawParticleVector(p, v, pnum, DIM, Vec3(1.0, 0.8, 0.8).data, Vec3(1.0, 0.0, 0.0).data, norm(dim)*m_fVScale*0.02);
		}
	}

	// 
	// 分割セル
	//
	if((m_iDraw & RXD_CELLS)){
		glDisable(GL_LIGHTING);
		glLineWidth(1.0);
		m_pPS->DrawCells(Vec3(0.0, 1.0, 0.0), Vec3(1.0, 0.0, 0.0));
	}

	//
	// デバッグ用ベクトル場
	//
	if(m_iDraw & RXD_DEBUG_VECTOR){
		RXREAL *p = m_pPS->GetArrayVBO(rxParticleSystemBase::RX_POSITION);
		RXREAL *v = m_pPS->GetArrayVBO(rxParticleSystemBase::RX_DEBUG_VECTOR, true);
		if(p && v){
			glDisable(GL_LIGHTING);
			glLineWidth(2.0);
			glColor3d(1.0, 0.0, 0.0);
			DrawParticleVector(p, v, pnum, DIM, Vec3(1.0, 0.8, 0.8).data, Vec3(1.0, 0.0, 0.0).data, norm(dim)*m_fVScale*0.02);
		}
	}

	//
	//高速化用パス
	//
	if(m_iColorType == rxParticleSystemBase::RX_ICE_FAST_PATH)
	{
		RXREAL *p = m_pPS->GetArrayVBO(rxParticleSystemBase::RX_POSITION);
		RXREAL *v = m_pPS->GetArrayVBO(rxParticleSystemBase::RX_VELOCITY);

		if(p && v)
		{
			glDisable(GL_LIGHTING);
			glLineWidth(3.0);
			glColor3d(0.0, 1.0, 1.0);
			DrawFastPath(p, v, pnum, DIM, Vec3(0.8, 0.8, 1.0).data, Vec3(0.0, 0.0, 1.0).data, norm(dim)*m_fVScale);
		}
	}

	//iceObjectのデバッグ情報の描画
	m_iceObj->Display();
}



//-----------------------------------------------------------------------------
// 各種関数
//-----------------------------------------------------------------------------
/*!
 * xyz軸描画(x軸:赤,y軸:緑,z軸:青)
 * @param[in] len 軸の長さ
 */
int DrawAxis(double len, double line_width)
{
	glLineWidth((GLfloat)line_width);

	// x軸
	glColor3f(1.0, 0.0, 0.0);
	glBegin(GL_LINES);
	glVertex3d(0.0, 0.0, 0.0);
	glVertex3d(len, 0.0, 0.0);
	glEnd();

	// y軸
	glColor3f(0.0, 1.0, 0.0);
	glBegin(GL_LINES);
	glVertex3d(0.0, 0.0, 0.0);
	glVertex3d(0.0, len, 0.0);
	glEnd();

	// z軸
	glColor3f(0.0, 0.0, 1.0);
	glBegin(GL_LINES);
	glVertex3d(0.0, 0.0, 0.0);
	glVertex3d(0.0, 0.0, len);
	glEnd();

	return 1;
}

/*!
 * 任意ベクトル周りの回転
 * @param[in] pos 元の座標値
 * @param[in] axis 回転軸
 * @param[in] ang 回転角度(rad)
 * @return 回転した座標値
 */
inline Vec3 Rotate(const Vec3 &pos, const Vec3 &axis, const double &ang)
{
	Vec3 pos1;    // 回転後の座標値

	double c = cos(ang);
	double s = sin(ang);
	double x, y, z;
	x = axis[0]; y = axis[1]; z = axis[2];
	
	// | xx(1-c)+c	xy(1-c)-zs	xz(1-c)+ys	0 |
	// | yx(1-c)-zs	yy(1-c)+c	yz(1-c)-xs	0 |
	// | xz(1-c)-ys	yz(1-c)+xs	zz(1-c)+c	0 |
	// | 0			0			0			1 |
	pos1[0] =   (x*x*(1.0-c)+c)*pos[0] +(x*y*(1.0-c)-z*s)*pos[1] +(x*z*(1.0-c)+y*s)*pos[2];
	pos1[1] = (y*x*(1.0-c)-z*s)*pos[0]   +(y*y*(1.0-c)+c)*pos[1] +(y*z*(1.0-c)-x*s)*pos[2];
	pos1[2] = (x*z*(1.0-c)-y*s)*pos[0] +(y*z*(1.0-c)+x*s)*pos[1]   +(z*z*(1.0-c)+c)*pos[2];
 
	return pos1;
}


/*!
 * 文字列描画
 * @param[in] cir_str 文字列循環バッファ
 * @param[in] static_str 静的な文字列バッファ
 * @param[in] w,h ウィンドウサイズ
 */
void DrawStrings(vector<string> &static_str, int w, int h, int offsetx, int offsety)
{
	glDisable(GL_LIGHTING);
	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();
	gluOrtho2D(0, w, 0, h);
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadIdentity();

	float x0 = offsetx;
	float y0 = h-offsety;

	glRasterPos2f(x0, y0);

	if(g_pFont){
		// FTGLで文字列を描画
		for(int j = 0; j < (int)static_str.size(); ++j){
			glRasterPos2f(x0, y0);
			g_pFont->Render(static_str[j].c_str());
			y0 -= g_pFont->LineHeight();
		}
	}
	else{
		// glutBitmapStringで文字列描画
		for(int j = 0; j < (int)static_str.size(); ++j){
			glRasterPos2f(x0, y0);
			glutBitmapString(GLUT_BITMAP_HELVETICA_12, (unsigned char*)static_str[j].c_str());
			y0 -= 20;
		}
	}

	glPopMatrix();
	glMatrixMode(GL_PROJECTION);
	glPopMatrix();
	glMatrixMode(GL_MODELVIEW);
}


/*!
 * 文字列描画
 * @param[in] cir_str 文字列循環バッファ
 * @param[in] static_str 静的な文字列バッファ
 * @param[in] w,h ウィンドウサイズ
 */
void DrawStringsBottom(vector<string> &static_str, int w, int h, int offsetx, int offsety)
{
	//w *= 0.5;
	//h *= 0.5;
	
	glDisable(GL_LIGHTING);
	//glColor3f(0.0, 0.0, 0.0);
	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();
	gluOrtho2D(0, w, 0, h);
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadIdentity();

	float x0 = 20.0f;

	if(g_pFont){
		float y0 = static_str.size()*(g_pFont->LineHeight())+5;
		glRasterPos2f(20.0f+offsetx, y0);

		// FTGLで文字列描画
		for(int j = 0; j < (int)static_str.size(); ++j){
			glRasterPos2i(20+offsetx, y0);
			g_pFont->Render(static_str[j].c_str());

			y0 -= g_pFont->LineHeight();
		}
	}
	else{
		float y0 = static_str.size()*20.0f;
		glRasterPos2f(20.0f+offsetx, y0);

		// glutBitmapStringで文字列描画
		for(int j = 0; j < (int)static_str.size(); ++j){
			glRasterPos2i(20+offsetx, y0);

			glutBitmapString(GLUT_BITMAP_HELVETICA_12, (unsigned char*)static_str[j].c_str());

			y0 -= 20;
		}
	}

	glPopMatrix();
	glMatrixMode(GL_PROJECTION);
	glPopMatrix();
	glMatrixMode(GL_MODELVIEW);
}



RXREAL CalDeterminant3x3(const RXREAL *m)
{
	return m[0]*m[4]*m[8]+m[3]*m[7]*m[2]+m[6]*m[1]*m[5]-m[0]*m[7]*m[5]-m[6]*m[4]*m[2]-m[3]*m[1]*m[8];
}

void CalInverse3x3(const RXREAL *m, RXREAL *invm)
{
	RXREAL d = m[0]*m[4]*m[8]+m[3]*m[7]*m[2]+m[6]*m[1]*m[5]-m[0]*m[7]*m[5]-m[6]*m[4]*m[2]-m[3]*m[1]*m[8];

	if(d == 0) d = 1;

	RXREAL inv_det = 1.0/d;
	invm[0] = inv_det*(m[4]*m[8]-m[5]*m[7]);
	invm[1] = inv_det*(m[2]*m[7]-m[1]*m[8]);
	invm[2] = inv_det*(m[1]*m[5]-m[2]*m[4]);
	
	invm[3] = inv_det*(m[5]*m[6]-m[3]*m[8]);
	invm[4] = inv_det*(m[0]*m[8]-m[2]*m[6]);
	invm[5] = inv_det*(m[2]*m[3]-m[0]*m[5]);
	
	invm[6] = inv_det*(m[3]*m[7]-m[4]*m[6]);
	invm[7] = inv_det*(m[1]*m[6]-m[0]*m[7]);
	invm[8] = inv_det*(m[0]*m[4]-m[1]*m[3]);
}






//-----------------------------------------------------------------------------
// 画像の読込
//-----------------------------------------------------------------------------
static int m_iBitmapType = RX_BMP_WINDOWS_V3;
inline void SetBitmapType(int type){ m_iBitmapType = type; }


unsigned char* ReadImageFile(const std::string &fn, int &w, int &h, int &c)
{
	// 拡張子抽出
	string ext;
	size_t pos1 = fn.rfind('.');
	if(pos1 != string::npos){
		ext = fn.substr(pos1+1, fn.size()-pos1);
		string::iterator itr = ext.begin();
		while(itr != ext.end()){
			*itr = tolower(*itr);
			itr++;
		}
		itr = ext.end()-1;
		while(itr != ext.begin()){	// パスの最後に\0やスペースがあったときの対策
			if(*itr == 0 || *itr == 32){
				ext.erase(itr--);
			}
			else{
				itr--;
			}
		}
	}

	// 画像読み込み
	unsigned char* pimg;
	if(ext == "bmp"){
		pimg = ReadBitmapFile(fn, w, h, c);
	}
	else if(ext == "png"){
		pimg = ReadPngFile(fn, w, h, c);
	}
	//else if(ext == "jpg" || ext == "jpeg"){
	//	pimg = ReadJpegFile(fn, w, h, c);
	//}
	else{
		return 0;
	}

	return pimg;
}


int WriteImageFile(const std::string &fn, unsigned char *img, int w, int h, int c, int quality)
{
	// 拡張子抽出
	string ext;
	size_t pos1 = fn.rfind('.');
	if(pos1 != string::npos){
		ext = fn.substr(pos1+1, fn.size()-pos1);
		string::iterator itr = ext.begin();
		while(itr != ext.end()){
			*itr = tolower(*itr);
			itr++;
		}
		itr = ext.end()-1;
		while(itr != ext.begin()){	// パスの最後に\0やスペースがあったときの対策
			if(*itr == 0 || *itr == 32){
				ext.erase(itr--);
			}
			else{
				itr--;
			}
		}
	}

	// 画像読み込み
	if(ext == "bmp"){
		WriteBitmapFile(fn, img, w, h, c, m_iBitmapType);
	}
	else if(ext == "png"){
		WritePngFile(fn, img, w, h, c);
	}
	//else if(ext == "jpg" || ext == "jpeg"){
	//	WriteJpegFile(fn, img, w, h, c, quality);
	//}
	else{
		return 0;
	}

	return 1;
}


/*! 
 * フレームバッファのRGB情報を一時的なバッファに確保
 * @retval true 保存成功
 * @retval false 保存失敗
 */
bool SaveFrameBuffer(const string &fn, int w, int h)
{
	vector<unsigned char> imm_buf(w*h*3);

	glReadPixels(0, 0, w, h, GL_RGB, GL_UNSIGNED_BYTE, &imm_buf[0]);

	// 上下反転
	int stride = w*3;
	for(int j = 0; j < h/2; ++j){
		for(int i = 0; i < stride; ++i){
			unsigned char tmp = imm_buf[j*stride+i];
			imm_buf[j*stride+i] = imm_buf[(h-j-1)*stride+i];
			imm_buf[(h-j-1)*stride+i] = tmp;
		}
	}

	WriteImageFile(fn, &imm_buf[0], w, h, 3);

	return true;
}


/*! テクスチャをビットマップで保存
	@param[in] file_name 保存ファイル名
	@param[in] tex_obj 保存したいテクスチャオブジェクト
	@param[in] jpeg_quality JPEGで保存する場合の保存品質(0-255)
	@retval true 保存成功
	@retval false 保存失敗
*/
bool SaveTexture(const string &fn, rxTexObj2D &tex_obj, int jpeg_quality)
{
	int w = tex_obj.m_iW;
	int h = tex_obj.m_iH;
	int c = 3;
	vector<unsigned char> imm_buf(w*h*c);

	int ic, jc, idx;
	for(jc = 0; jc < h; ++jc){
		for(ic = 0; ic < w; ++ic){
			idx = 3*(jc*w+ic);

			imm_buf[idx+0] = tex_obj.m_pImage[idx];
			imm_buf[idx+1] = tex_obj.m_pImage[idx+1];
			imm_buf[idx+2] = tex_obj.m_pImage[idx+2];
		}
	}

	WriteImageFile(fn, &imm_buf[0], w, h, c, jpeg_quality);

	return true;
}

/*! 
 * ReadTexture テクスチャの読み込み
 * @param[in] path テクスチャ画像のパス
 * @param[out] tex_obj テクスチャを格納する
 * @return テクスチャが読み込めたかどうか
 */
bool ReadTexture(const string &fn, rxTexObj2D &tex_obj)
{
	int w0, h0, c0;
	unsigned char* pimg;
	pimg = ReadImageFile(fn, w0, h0, c0);
	if(pimg == 0){
		return false;
	}

	tex_obj.SetSize(w0, h0, c0);

	int ic, jc;
	for(jc = 0; jc < tex_obj.m_iH; ++jc){
		for(ic = 0; ic < tex_obj.m_iW; ++ic){
			int idx = 3*(jc*w0+ic);
			tex_obj.SetColor(ic, h0-jc-1, pimg[idx+0], pimg[idx+1], pimg[idx+2]);
		}
	}

	// テクスチャ登録
	tex_obj.Bind();

	// テクスチャパラメータの設定
	tex_obj.SetParameter(GL_TEXTURE_WRAP_S, GL_REPEAT);
	tex_obj.SetParameter(GL_TEXTURE_WRAP_T, GL_REPEAT);
	tex_obj.SetParameter(GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	tex_obj.SetParameter(GL_TEXTURE_MIN_FILTER, GL_NEAREST);

	tex_obj.SetTexture();

	return true;
}


/*!
 * OpenCVで画像読込み -> OpenGLテクスチャ登録
 * @param[in] fn ファイル名
 * @param[inout] tex_name テクスチャ名(0なら新たに生成)
 * @param[in] mipmap ミップマップ使用フラグ
 * @param[in] compress テクスチャ圧縮使用フラグ
int LoadGLTexture(const string &fn, GLuint &tex_name, bool mipmap, bool compress)
{
	// 画像読み込み
	int w, h, c;
	unsigned char* pimg;
	pimg = ReadImageFile(fn, w, h, c);
	if(!pimg){
		return 0;
	}

	cout << "image : " << w << " x " << h << " x " << c << endl;
	GLuint iformat, format;

	// 画像フォーマット
	format = GL_RGBA;
	if(c == 1){
		format = GL_LUMINANCE;
	}
	else if(c == 3){
		format = GL_RGB;
	}
 
	// OpenGL内部の格納フォーマット
	if(compress){
		iformat = GL_COMPRESSED_RGBA_S3TC_DXT1_EXT;
		if(c == 1){
			iformat = GL_COMPRESSED_LUMINANCE_ARB;
		}
		else if(c == 3){
			iformat = GL_COMPRESSED_RGB_S3TC_DXT1_EXT ;
		}
	}
	else{
		iformat = GL_RGBA;
		if(c == 1){
			iformat = GL_LUMINANCE;
		}
		else if(c == 3){
			iformat = GL_RGB;
		}
	}
 
	//glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
 
	// テクスチャ作成
	if(tex_name == 0){
		glGenTextures(1, &tex_name);
 
		glBindTexture(GL_TEXTURE_2D, tex_name);
		
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, (mipmap ? GL_LINEAR_MIPMAP_LINEAR : GL_LINEAR));
		glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
		glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
 
		if(mipmap){
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_BASE_LEVEL, 0);
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAX_LEVEL, 6);
		}
 
		glTexImage2D(GL_TEXTURE_2D, 0, iformat, w, h, 0, format, GL_UNSIGNED_BYTE, pimg);
 
		if(mipmap){
			glGenerateMipmapEXT(GL_TEXTURE_2D);
		}
	}
	else{
		glBindTexture(GL_TEXTURE_2D, tex_name);
		//glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, w, h, format, GL_UNSIGNED_BYTE, pimg);
		glTexImage2D(GL_TEXTURE_2D, 0, iformat, w, h, 0, format, GL_UNSIGNED_BYTE, pimg);

		if(mipmap){
			glGenerateMipmapEXT(GL_TEXTURE_2D);
		}
	}
 
	glBindTexture(GL_TEXTURE_2D, 0);

	delete [] pimg;

	return 1;
}
 */

/*!
 * 透過PNG画像読み込み
 * @param[in] fn ファイル名
 * @param[inout] tex_name テクスチャ名(0なら新たに生成)
 * @param[out] w,h 読み込まれた画像のサイズ
 */
bool LoadTPNG(const string &fn, GLuint &tex_name, int &w, int &h, int &c)
{
	int w0, h0, c0;
	unsigned char* pimg;
	pimg = ReadPngFile(fn, w0, h0, c0);
	if(pimg == NULL){
		return false;
	}

	w = w0;
	h = h0;
	c = c0;

	//RXCOUT << "read(libpng) : " << w << "x" << h << "x" << c << endl;

	//cvFlip(img, NULL, 0);

	GLuint format;
	format = GL_RGBA;
	if(c == 1){
		format = GL_LUMINANCE;
	}
	else if(c == 3){
		format = GL_RGB;
	}

	GLuint iformat;
	iformat = GL_RGBA;
	if(c == 1){
		iformat = GL_LUMINANCE;
	}
	else if(c == 3){
		iformat = GL_RGB;
	}

	// テクスチャ作成
	if(tex_name == 0){
		glGenTextures(1, &tex_name);
	}

	glBindTexture(GL_TEXTURE_2D, tex_name);

	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);

	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	//glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);

	//glTexImage2D(GL_TEXTURE_2D, 0, format, w, h, 0, iformat, GL_UNSIGNED_BYTE, pimg);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_BASE_LEVEL, 0);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAX_LEVEL, 6);
	gluBuild2DMipmaps(GL_TEXTURE_2D, format, w, h, iformat, GL_UNSIGNED_BYTE, pimg); 

	glBindTexture(GL_TEXTURE_2D, 0);	

	free(pimg);

	return true;
}
bool LoadTPNG(const string &fn, GLuint &tex_name)
{
	int w, h, c;
	return LoadTPNG(fn, tex_name, w, h, c);
}



/*! 
 * 環境マップ用のキューブマップテクスチャの読み込み
 * @param[in] fn[6] テクスチャ画像(6枚)のパス(x+,x-,y+,y-,z+,z-)(右,左,上,下,後,前)
 * @param[out] cube_map rxCubeMapData型
 * @retval true  キューブマップ用画像の読み込み成功
 * @retval false キューブマップ用画像の読み込み失敗
 */
bool LoadCubeMapTexture(const string fn[6], rxCubeMapData &cube_map)
{
	GLuint tex_name;
	glGenTextures(1, &tex_name);
	glBindTexture(GL_TEXTURE_CUBE_MAP, tex_name);

	// キューブマップテクスチャパラメータの設定
	glTexParameterf(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);		// 画像境界の扱いの指定
	glTexParameterf(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
	glTexParameterf(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);
	glTexParameterf(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MAG_FILTER, GL_NEAREST);	// 画像フィルタの指定
	glTexParameterf(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	
	glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_BASE_LEVEL, 0);
	glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MAX_LEVEL, 6);

	glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
	
	GLenum target[6] = { GL_TEXTURE_CUBE_MAP_POSITIVE_X, GL_TEXTURE_CUBE_MAP_NEGATIVE_X, 
						 GL_TEXTURE_CUBE_MAP_POSITIVE_Y, GL_TEXTURE_CUBE_MAP_NEGATIVE_Y, 
						 GL_TEXTURE_CUBE_MAP_POSITIVE_Z, GL_TEXTURE_CUBE_MAP_NEGATIVE_Z };

	for(int i = 0; i < 6; ++i){
		int w, h, c;
		unsigned char* pimg;
		pimg = ReadImageFile(fn[i], w, h, c);
		if(!pimg){
			return false;
		}

		GLuint format;
		format = GL_RGBA;
		if(c == 1){
			format = GL_LUMINANCE;
		}
		else if(c == 3){
			format = GL_RGB;
		}

		GLuint iformat;
		iformat = GL_RGBA;
		if(c == 1){
			iformat = GL_LUMINANCE;
		}
		else if(c == 3){
			iformat = GL_RGB;
		}

		gluBuild2DMipmaps(target[i], format, w, h, iformat, GL_UNSIGNED_BYTE, pimg); 


		free(pimg);	
	}

	glBindTexture(GL_TEXTURE_2D, 0);	

	cube_map.id = tex_name;

	return true;
}

/*! 
 * 環境マップ用のキューブマップテクスチャの読み込み
 * @param cube_map キューブマップデータ
 * @param base キューブマップ用画像のファイル名のベース部分
 * @param ext キューブマップ用画像のファイルの拡張子
 * @retval true  キューブマップ用画像の読み込み成功
 * @retval false キューブマップ用画像の読み込み失敗
 */
bool LoadCubeMap(rxCubeMapData &cube_map, string base, string ext)
{
	// キューブマップ用画像の読み込み(x+,x-,y+,y-,z+,z-)(右,左,上,下,後,前)
	string fn[6];
	fn[0] = base+"posx"+ext;
	fn[1] = base+"negx"+ext;
	fn[2] = base+"posy"+ext;
	fn[3] = base+"negy"+ext;
	fn[4] = base+"posz"+ext;
	fn[5] = base+"negz"+ext;
	
	if(!LoadCubeMapTexture(fn, cube_map)){
		return false;
	}

	return true;
}
