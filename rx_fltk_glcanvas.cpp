/*!
  @file rx_fltk_glcanvas.cpp
	
  @brief FLTK�ɂ��OpenGL�E�B���h�E�N���X
 
  @author Makoto Fujisawa 
  @date   2011-09
*/
// FILE --rx_fltk_glcanvas.cpp--

#pragma warning (disable: 4996)
#pragma warning (disable: 4819)


//-----------------------------------------------------------------------------
// �C���N���[�h�t�@�C��
//-----------------------------------------------------------------------------
#include "rx_fltk_glcanvas.h"
#include "rx_fltk_window.h"

// �e�N�X�`���E�摜
#include "rx_texture.h"
//#include "rx_jpeg.h"
#include "rx_png.h"
#include "rx_bitmap.h"

// �V�~�����[�V����
#include "rx_sph.h"

// �ݒ�t�@�C��
#include "rx_atom_ini.h"

// OpenGL�`��֘A
//#include "rx_gltexture.h"	// GL�e�N�X�`�� & PNG�ɂ��摜�ۑ�
#include "rx_trackball.h"	// ���_�ύX�p�g���b�N�{�[���N���X
#include "rx_glgui.h"		// GL���g����GUI
#include "rx_gldraw.h"		// GL�`��֐��Q
#include "rx_shaders.h"		// GLSL�֐�
#include "rx_shadow.h"		// �V���h�E�}�b�v�ɂ��e�t��

#include "rx_pick.h"

// CUDA
#include "rx_cu_funcs.cuh"

//#include <helper_cuda.h>
#include <helper_timer.h>

// ���b�V����
#include "rx_mesh.h"		// ���b�V���\���́C�N���X��`
#include "rx_mc.h"			// Marching Cubes
#include "rx_ssm.h"			// Screen Space Mesh

// OpenCV
//#include <opencv2/opencv.hpp>

// FTGL
#include <FTGL/ftgl.h>


//-----------------------------------------------------------------------------
// �萔�E�ϐ�
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
const int RX_DISPLAY_STEP = 200;

// �v�Z���ʏo�͂̃f�t�H���g�t�H���_
const string RX_DEFAULT_RESULT_DIR = "result/";
const string RX_DEFAULT_IMAGE_DIR  = RX_DEFAULT_RESULT_DIR+"images/";
const string RX_DEFAULT_MESH_DIR   = RX_DEFAULT_RESULT_DIR+"mesh/";
const string RX_DEFAULT_DATA_DIR   = RX_DEFAULT_RESULT_DIR+"data/";

// �ݒ�t�@�C��
extern rxINI *g_pINI;

// �ݒ�t�@�C���ւ̕ۑ��p
double g_fTBTran[3] = {0, 0, -5};	//!< ���_�ړ��p�g���b�N�{�[���̕��s�ړ���
double g_fTBQuat[4] = {1, 0, 0, 0};	//!< ���_�ړ��p�g���b�N�{�[���̉�]��

// FTGL
#define FONT_FILE "Inconsolata.ttf"
static FTPixmapFont* g_pFont;
	
int g_iIterations = 1;				//!< �C��������
double g_fEta = 0.0;				//!< ���x�ϓ���

// �`��
rxGLSL g_glslPointSprite;			//!< GLSL���g�����`��
rxGLSL g_glslFresnel;				//!< GLSL���g�����`��

// ���Ԍv��
rxTimer g_TimerFPS;					//!< FPS����p�^�C�}�[

// ���όv�Z���Ԍv��
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

//! �ǉ�::�}�E�X�s�b�N
vector<int> g_vSelectedVertices;
int g_iPickedObj = -1;
double g_fPickDist = 1.0;	//!< �s�b�N���ꂽ�_�܂ł̋���


//-----------------------------------------------------------------------------
// �֐��v���g�^�C�v�錾
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
// rxFlWindow�N���X�̎���
//-----------------------------------------------------------------------------

//! �R���X�g���N�^
rxFlGLWindow::rxFlGLWindow(int x_, int y_, int w_, int h_, const char* l, void *parent)
	: Fl_Gl_Window(x_, y_, w_, h_, l), m_iWinW(w_), m_iWinH(h_)
{
	m_pParent = (rxFlWindow*)parent;
	m_bAnimation = false;
	resizable(this);
	end();


	// �`��t���O
	m_iDraw = 0;

	// �t�H���g�T�C�Y
	m_ulFontSize = 14;


	//
	// �V�~�����[�V�����p�ϐ��̏�����
	//
	m_pPS = 0;
	m_fDt = 0.005;
	m_fGravity = 9.80665;

	// �A�j���[�V�����X�e�b�v
	m_iCurrentStep = 0;
	m_bPause = false;

	// �V�[���t�@�C���̓ǂݎ��
	m_Scene.ReadSceneFiles();

	// �V�[���^�C�g�����X�g
	m_vSceneTitles = m_Scene.GetSceneTitles();

	// ���݂̃V�[��
	m_iCurrentSceneIdx = m_Scene.GetCurrentSceneIdx();

	m_iSimuSetting = 0;

	m_iSaveImageSpacing = -1;

	m_iColorType = rxParticleSystemBase::RX_RAMP;
	m_fVScale = 0.02;

	// �p�[�e�B�N�����o��
	m_strSphOutputName0  = RX_DEFAULT_DATA_DIR+"sph_setting.dat";
	m_strSphOutputHeader = RX_DEFAULT_DATA_DIR+"sph_particles_";

	// �p�[�e�B�N��������
	m_strSphInputName0  = RX_DEFAULT_DATA_DIR+"sph_setting.dat";
	m_strSphInputHeader = RX_DEFAULT_DATA_DIR+"sph_particles_";

	// �ő̈ړ��t���O
	m_bSolidMove = false;

	// �ő̂̓���
	m_vMovePos[0] = Vec3(0.0);
	m_vMovePos[1] = Vec3(0.0);
	m_fMoveMaxVel = 0.0;
	m_bMoveSolid = false;
	m_iMoveStart = -1;

	m_fIceFlag = 0;

	//
	// ���b�V��
	//
	m_iNumVrts = 0; m_iNumTris = 0;		// �������ꂽ���b�V���̒��_���ƃ��b�V����
	m_iVertexStore = 5;					// �T���v�����O�{�����[�����ɑ΂���\�z����钸�_��(nx*ny*store)

	m_fMeshThr = 400.0;					// �A�֐����b�V��������臒l
	m_iMeshMaxN = 128;					// ���b�V�����O���b�h��(���E�������Ƃ������������̕�����)
	m_iMeshN[0] = 64;					// ���b�V�����O���b�h��
	m_iMeshN[1] = 64;
	m_iMeshN[2] = 64;

	m_vMeshBoundaryExt = Vec3(1.0);		// ���b�V�����E�{�b�N�X�̊e�ӂ̒�����1/2
	m_vMeshBoundaryCen = Vec3(0.0);		// ���b�V�����E�{�b�N�X�̒��S���W

	m_iSolidDraw = 2;

	m_uVrtVBO = 0;
	m_uTriVBO = 0;
	m_uNrmVBO = 0;
	m_iDimVBO = 4;

	// ���b�V���o��
	m_iSaveMeshSpacing = RX_SAVE_IMAGE_SPACING;

	// �w�i�摜
	m_bUseCubeMap = false;		// �L���[�u�}�b�v�g�p�t���O

	g_fTotalTime = 0;
	g_fAvgTime = 0;
	g_iTimeCount = 0;

	m_pMCMeshCPU = 0;
	m_pMCMeshGPU = 0;
	m_pSSMCPU = 0;				// Screen Space Mesh
	m_pSSMGPU = 0;				// Screen Space Mesh by CUDA
	m_iDepthFiltering = 1;		// ������(�f�v�X�}�b�v0x01, �֊s0x02)
	m_iSSDebugOutput = 0;
	m_iTriangulationMethod = RXM_MC_GPU;

	m_fSpacing = 2;				// �f�v�X�}�b�v�̃T���v�����O�Ԋu
	m_fPrtRad = 0.01;			// �p�[�e�B�N���̔��a
	m_fZmax = 1.8*m_fPrtRad;	// �֊s�ƂȂ�f�v�X����臒l
	m_iNfilter = 1;				// �f�v�X�l�������̃t�B���^�T�C�Y
	m_iNiters = 3;				// �֊s�������̔�����


	m_iPickedParticle = -1;
	

	// �ݒ�t�@�C��
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

//! �f�X�g���N�^
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
	if(m_ht)  delete m_ht;
}


/*! 
 * GL�̏������֐�
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

	// OpenGL�g���̃o�[�W�����`�F�b�N
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

	// �����ݒ�
	glEnable(GL_LIGHT0);
	glEnable(GL_LIGHTING);
	glLightfv(GL_LIGHT0, GL_POSITION, RX_LIGHT0_POS);
	glLightfv(GL_LIGHT0, GL_DIFFUSE,  RX_LIGHT_DIFF);
	glLightfv(GL_LIGHT0, GL_SPECULAR, RX_LIGHT_SPEC);
	glLightfv(GL_LIGHT0, GL_AMBIENT,  RX_LIGHT_AMBI);

	glShadeModel(GL_SMOOTH);

	glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
	glEnable(GL_COLOR_MATERIAL);

	// ���_������
	//InitView();

	// �f�t�H���g�ގ�
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

	// �L���[�u�}�b�v�̓ǂݍ���
//	if(LoadCubeMap(m_CubeMap, "texture/terragen0_", ".png")){
//	if(LoadCubeMap(m_CubeMap, "texture/cubemap_", ".png")){
	if(LoadCubeMap(m_CubeMap, "texture/nvlobby_new_", ".png")){
		m_bUseCubeMap = true;
	}
	else{
		RXCOUT << "error : can't load the cube map" << endl;
	}

	// �t�H���g�ǂݍ���
	SetupFonts(FONT_FILE);
	//RXCOUT << "char width : " << glutBitmapWidth(GLUT_BITMAP_HELVETICA_12, 'W') << endl;

	//g_InletLine.span = -1;
	m_bsSimuSetting.set(ID_SPH_ANISOTROPIC, false);

	// �g���b�N�{�[�������p��
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

	m_bsSimuSetting.set(ID_HEAT, false);	//�ǉ�
	m_bsSimuSetting.set(ID_SM, false);		//�ǉ�
	m_bsSimuSetting.set(ID_ICE, false);		//�ǉ�

	// �`��t���O������
	m_iDraw = 0;
	m_iDraw |= RXD_BBOX;
	m_iDraw |= RXD_PARTICLE;
	m_iDraw |= RXD_SOLID;
	m_iDraw |= RXD_PARAMS;
	if(m_bsSimuSetting.at(ID_SPH_MESH)) m_iDraw |= RXD_MESH;

	m_iDrawPS = RXP_POINTSPRITE;


	// �J�����g�̃V�[���ݒ�
	m_Scene.SetCurrentScene(m_iCurrentSceneIdx);

	// SPH������
	InitSPH(m_Scene);
	
	//OpenMP�ɂ��ȈՕ��񏈗�
#ifdef _OPENMP
    cout << "OpenMP : On, threads =" << omp_get_max_threads() << endl;
#endif
	
	m_ht = 0;
	m_ice = 0;

	//�ǉ�	����������
	InitHT(m_Scene);		//�M����������
	InitICE();				//�X������
//	InitTetra();			//�l�ʑ̏�����
	InitCluster();			//�N���X�^������
	InitICE_Cluster();		//���q�ƃN���X�^�̊֌W����������

	// GLSL�̃R���p�C��
	g_glslPointSprite = CreateGLSL(ps_vs, ps_fs, "point sprite");
	g_glslFresnel     = CreateGLSL(fresnel2_vs, fresnel2_fs, "fresnel");


	m_pParent->UpdateMenuState();
}

/*!
 * �t�H���g�̐ݒ�
 * @param[in] file �t�H���g�t�@�C���p�X
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
 * ���T�C�Y�C�x���g�����֐�
 * @param[in] w �L�����o�X��(�s�N�Z����)
 * @param[in] h �L�����o�X����(�s�N�Z����)
 */
void rxFlGLWindow::Resize(int w, int h)
{
	m_iWinW = w;
	m_iWinH = h;

	cout << "m_iWinW = " << m_iWinW << " m_iWinH = " << m_iWinH << endl;

	glViewport(0, 0, w, h);
	m_tbView.SetRegion(w, h);

	// �����ϊ��s��̐ݒ�
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	Projection();

	// ���f���r���[�ϊ��s��̐ݒ�
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
}

/*!
 * �������e�ϊ�
 */
void rxFlGLWindow::Projection(void)
{
	gluPerspective(RX_FOV, (float)m_iWinW/(float)m_iWinH, 0.2f, 1000.0f);
	//glOrtho(-1, 1, -1, 1, -1, 1);
}

/*!
 * �ĕ`�施��
 */
void rxFlGLWindow::ReDisplay(void)
{
	redraw();
}

/*!
 * ���_�̏�����
 */
void rxFlGLWindow::InitView(void)
{
	double q[4] = {1, 0, 0, 0};
	m_tbView.SetQuaternion(q);
	m_tbView.SetScaling(-6.0);
	m_tbView.SetTranslation(0.0, 0.0);
}

/*!
 * �ĕ`��C�x���g�����֐�
 */
void rxFlGLWindow::Display(void)
{
	make_current();

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glViewport(0, 0, m_iWinW, m_iWinH);

	// �����ϊ��s��̐ݒ�
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(RX_FOV, (float)m_iWinW/(float)m_iWinH, 0.001f, 1000.0f);

	// ���f���r���[�s�񏉊���
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glPushMatrix();

	//glScalef(10.0, 10.0, 10.0);

	// �}�E�X�ɂ���]�E���s�ړ��̓K�p
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
		
	// OpenGL���e�ϊ��s��擾
	glGetDoublev(GL_PROJECTION_MATRIX, m_fProjectionMatrix);

	// OpenGL���f���r���[�ϊ��s��擾
	glGetDoublev(GL_MODELVIEW_MATRIX, m_fModelviewMatrix);

	RenderSphScene();

	glPopMatrix();

	glPopMatrix();

	//��`�`��
	if( m_ht_vStartPoint[0] != 0.0 && m_ht_vStartPoint[1] != 0.0
	 && m_ht_vEndPoint[0] != 0.0   && m_ht_vEndPoint[1] != 0.0 )
	{
		DrawRubber(1, m_ht_vStartPoint, m_ht_vEndPoint, m_iWinW, m_iWinH);
	}

	// ��ʕ�����`��
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

	// �}�E�X�ɂ���]�E���s�ړ��̓K�p
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

//���q�I��
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
 * �}�E�X�C�x���g�����֐�
 * @param[in] button �}�E�X�{�^��(FL_LEFT_MOUSE,FL_MIDDLE_MOUSE,FL_RIGHT_MOUSE)
 * @param[in] state �}�E�X�{�^���̏��(1:down, 0:up)
 * @param[in] x,y �}�E�X���W(�X�N���[�����W�n)
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
			{	// �{�^���_�E��				
				PickParticle(x, y);			//���q�擾

				if( m_iPickedParticle != -1 )
				{
					vector<int> vrts;
					vrts.push_back(m_iPickedParticle);
					cout << "vertex " << m_iPickedParticle << endl;
					g_vSelectedVertices.resize(1);

//					// ���_����s�b�N�_�܂ł̋����v�Z
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
			{	// �{�^���A�b�v
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
		if(state)
		{

		}
		else
		{
			//�f�o�b�O
			if(debugIndx == 0)
			{
				for(int i = 0; i < m_pPS->GetNumParticles(); i++)
				{
					//if(m_ice->GetPtoCNum_Connect(i)==0){ continue;}
					//m_ice->DebugPtoC_Connect(i);				//���q���ڑ��N���X�^
				}
				debugIndx++;
			}
			else if(debugIndx == 1)
			{
				for(int i = 0; i < m_iClusteresNum; i++)
				{
					//m_ice->DebugCtoP_Connect(i);				//�ڑ��N���X�^�����q
				}
				debugIndx++;
			}
			else if(debugIndx == 2)
			{
				for(int i = 0; i < m_iClusteresNum; i++)
				{
					//m_ice->DebugNeighborCluster(i);				//�e�N���X�^�̋ߖT�N���X�^
				}
				debugIndx++;
			}
			else if(debugIndx == 3)
			{
				for(int i = 0; i < m_pPS->GetNumParticles(); i++)
				{
					//if(m_ice->GetPtoCNum_Calc(i)==0){ continue;}
					//m_ice->DebugPtoC_Calc(i);							//���q���v�Z�N���X�^
					if(m_ice->GetParticleNum() <= i){	continue;	}	//�Z���݂̂̎����̂Ƃ��ɕK�v�ɂȂ�D
					if(m_ice->GetPtoCNum(i) == 0)	{	continue;	}
					m_ice->DebugPtoC(i);
				}
				//debugIndx++;
				debugIndx = 4;
			}
			else if(debugIndx == 4)
			{
				for(int i = 0; i < m_iClusteresNum; i++)
				{
					//m_ice->DebugCtoP_Calc(i);					//�v�Z�N���X�^�����q
					if(m_ice->GetParticleNum() <= i){	continue;	}	//�Z���݂̂̎����̂Ƃ��ɕK�v�ɂȂ�D
					if(m_ice->GetCtoPNum(i) == 0){	continue;	}
					m_ice->DebugCtoP(i);
				}
				//debugIndx++;
				debugIndx = 3;
			}
			else if(debugIndx == 5)
			{
				//for(int i = 0; i < m_iClusteresNum; i++)
				//{
				//	m_ice->DebugCalcToConnect(i);				//�ڑ��N���X�^�̗��q���v�Z�N���X�^�̗��q
				//}
				debugIndx++;
			}
			else if(debugIndx == 6)
			{
				//for(int i = 0; i < m_iClusteresNum; i++)
				//{
				//	m_sm_connects[i]->DebugLayer();				//�ڑ��N���X�^�̗��q�̃��C���[
				//}
				debugIndx++;
			}
			else if(debugIndx == 7)
			{
				for(int i = 0; i < m_iClusteresNum; i++)
				//{
				//	m_sm_calcs[i]->DebugLayer();				//�v�Z�N���X�^�̗��q�̃��C���[
				//}
				debugIndx = 0;
			}
		}

	}
	else if(button == FL_RIGHT_MOUSE){
		//��`�͈͂��쐬���C�͈͓��̗��q��I���C���x�E�M�ʂ��㏸������
		//if(state)
		//{
		//	m_ht_vStartPoint = Vec2(x, y);
		//	m_ht_bRectFlag = true;
		//}
		//else
		//{
		//	m_ht_vEndPoint = Vec2(x, y);

		//	vector<int> vrts;
		//	vector<rxPickInfo> hits;
		//	rxGLPick pick;
		//	pick.Set(rxFlGLWindow::DisplayForPick_s, this, rxFlGLWindow::Projection_s, this);
		//	hits = pick.Pick(
		//						(int)m_ht_vStartPoint[0],
		//						(int)m_ht_vStartPoint[1],
		//						(int)m_ht_vEndPoint[0],
		//						(int)m_ht_vEndPoint[1]
		//	);
		//	for(vector<rxPickInfo>::iterator i = hits.begin(); i != hits.end(); ++i){
		//		vrts.push_back(i->name-1);
		//		cout << "pIndx = " << vrts[vrts.size()-1] << endl;
		//	}

		//	cout << vrts.size() << " vertices are selected." << endl;
		//	m_ht_vSelectedVertices = vrts;

		//	//���x�E�M�ʏ㏸
		//	for( unsigned i = 0; i < vrts.size(); i++ )
		//	{
		//		m_ht->setTemps(vrts[i], 1000);
		//		m_ht->setHeats(vrts[i], 1000);
		//	}
		//	ClearRect();
		//}

		//ClearPick();

		//�f�o�b�O�@�w�肵�����q��Z��
		//if(state)
		//{
		//	m_ht->setTemps(meltPIndx, 1000);
		//	m_ht->setHeats(meltPIndx, 1000);
		//	m_ht->calcTempAndHeat();								//�M�ʂ̉��x�ϊ��C���x�̔M�ʕϊ�
		//}
		//else
		//{
		//	m_ht->setTemps(meltPIndx, 1000);
		//	m_ht->setHeats(meltPIndx, 1000);
		//	m_ht->calcTempAndHeat();								//�M�ʂ̉��x�ϊ��C���x�̔M�ʕϊ�

		//	StepSolid_Melt(m_fDt);

		//	if(meltPIndx == 0)
		//	{
		//		meltPIndx = 3;
		//	}
		//	else if(meltPIndx == 1)
		//	{
		//		meltPIndx = 0;
		//	}
		//	else if(meltPIndx == 3)
		//	{
		//		meltPIndx = 9;
		//	}
		//	else if(meltPIndx == 9)
		//	{
		//		meltPIndx = 1;
		//	}
		//}

	redraw();
	}
}

/*!
 * �X���q�I���̃��Z�b�g
 */
void rxFlGLWindow::ClearPick(void)
{//	cout<< "clearPick" << endl;
	g_vSelectedVertices.clear();

	if( m_iPickedParticle != -1 )
	{
		if(m_ice == 0)	return;
		if(m_ht == 0)	return;
		if(m_ice->GetParticleNum() <= m_iPickedParticle){	return;	}	//�Z���݂̂̎����̂Ƃ��ɕK�v�ɂȂ�D

		//���q�x�[�X�N���X�^
		for( int i = 0; i < m_ice->GetPtoCIndx(m_iPickedParticle); i++ )
		{
			int cIndx = m_ice->GetPtoC(m_iPickedParticle, i, 0);
			int oIndx = m_ice->GetPtoC(m_iPickedParticle, i, 1);

			if(cIndx == -1 || oIndx == -1)
			{
				continue;
			}

			m_sm_cluster[cIndx]->UnFixVertex(oIndx);
		}

		m_iPickedParticle = -1;
	}
}

/*!
 * ��`�I��̈�̃��Z�b�g
 */
void rxFlGLWindow::ClearRect(void)
{
	m_ht_vStartPoint = Vec2(0.0, 0.0);
	m_ht_vEndPoint = Vec2(0.0, 0.0);
	m_ht_vSelectedVertices.clear();
	m_ht_bRectFlag = false;
}

/*!
 * ���[�V�����C�x���g�����֐�(�}�E�X�{�^�����������܂܃h���b�O)
 * @param[in] x,y �}�E�X���W(�X�N���[�����W�n)
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
		if(((m_iMouseButton >> 3) == GLUT_LEFT_BUTTON) && ((m_iMouseButton & 1) == GLUT_DOWN)){
			double dm = 0.4/(double)m_iWinW;

			m_pPS->MoveSphereObstacle(0, Vec3(0.0, dm*(y0-y), dm*(x-x0)));
		}
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
		//�N���X�^�ɑ����闱�q�̈ړ�
		else if( m_iPickedParticle != -1 )
		{
			if(m_ice == 0)	return;
			if(m_ht == 0)	return;
			if(m_ice->GetParticleNum() <= m_iPickedParticle){	return;	}	//�Z���݂̂̎����̂Ƃ��ɕK�v�ɂȂ�D

			Vec3 ray_from, ray_to;
			Vec3 init_pos = Vec3(0.0);
			m_tbView.CalLocalPos(ray_from, init_pos);
			m_tbView.GetRayTo(x, y, RX_FOV, ray_to);

			Vec3 dir = Unit(ray_to-ray_from);	// ���_����}�E�X�ʒu�ւ̃x�N�g��
			Vec3 new_pos = ray_from+dir*g_fPickDist;
			//int v = g_vSelectedVertices[0];

			//���q��������N���X�^�̐��l���X�V
			for( int i = 0; i < m_ice->GetPtoCIndx(m_iPickedParticle); i++ )
			{
				int cIndx = m_ice->GetPtoC(m_iPickedParticle, i, 0);
				int oIndx = m_ice->GetPtoC(m_iPickedParticle, i, 1);
				if(cIndx == -1 || oIndx == -1) continue;
	
				m_sm_cluster[cIndx]->FixVertex(oIndx, new_pos);
			}

			//���q�̉��x�㏸
			m_ht->setTemps(m_iPickedParticle, m_ht->getAirTemp());
			m_ht->setHeats(m_iPickedParticle, m_ht->getAirTemp());

			//���q�̈ړ�
			RXREAL *p  = m_pPS->GetArrayVBO(rxParticleSystemBase::RX_POSITION);
			p[DIM*m_iPickedParticle+0] = new_pos[0];
			p[DIM*m_iPickedParticle+1] = new_pos[1];
			p[DIM*m_iPickedParticle+2] = new_pos[2];
		}
	}
	//RXCOUT << "        " << (m_iMouseButton >> 3) << ", " << (m_iMouseButton & 1) << endl;

	x0 = x;
	y0 = y;

	redraw();
}

/*!
 * ���[�V�����C�x���g�����֐�(�}�E�X�{�^���������Ȃ��ړ�)
 * @param[in] x,y �}�E�X���W(�X�N���[�����W�n)
 */
void rxFlGLWindow::PassiveMotion(int x, int y)
{
	if(x < 0 || y < 0) return;
	//make_current();
	//redraw();
}

/*!
 * �L�[�{�[�h�C�x���g�����֐�
 * @param[in] key �L�[�̎��
 * @param[in] x,y �L�[�������ꂽ�Ƃ��̃}�E�X���W(�X�N���[�����W�n)
 */
void rxFlGLWindow::Keyboard(int key, int x, int y)
{
	make_current();
	m_iKeyMod = (Fl::event_state(FL_SHIFT) ? 2 : (Fl::event_state(FL_CTRL) ? 3 : 1));

	switch(key){
	case 'i':
		InitView();
		break;

	case 'S':	// ��ʂ̉摜�ۑ�
		SaveDisplay( ((m_iCurrentStep >= 0) ? m_iCurrentStep : 0) );
		break;

	case 'g':
		SwitchFullScreen();
		break;

	// 
	// �V�~�����[�V�����ݒ�
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
 * ����L�[�{�[�h�C�x���g�����֐�
 * @param[in] key �L�[�̎��
 * @param[in] x,y �L�[�������ꂽ�Ƃ��̃}�E�X���W(�X�N���[�����W�n)
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
 * �A�C�h���R�[���o�b�N�֐�
 * @param[in] x ���[�U��`�ϐ�
 */
void rxFlGLWindow::OnIdle_s(void* x)
{
	((rxFlGLWindow*)x)->Idle();
}
void rxFlGLWindow::Idle(void)
{
	make_current();

	if(m_iCurrentStep%RX_DISPLAY_STEP == 0){
		RXCOUT << "step " << m_iCurrentStep << endl;
	}
	// �o�͗p�X�e�b�v��
	int stp = m_iCurrentStep-2;
	stp = (m_iCurrentStep == 0) ? 0 : stp;
	stp = (stp < 0) ? 0 : stp;

	//
	// �O�̃t���[���̕`����摜�t�@�C���ɕۑ�
	//
	if(m_iSaveImageSpacing > 0){
		if(stp%m_iSaveImageSpacing == 0){
			SaveDisplay(stp);
		}
	}

	//
	// �O�̃t���[���̃f�[�^���t�@�C���ۑ�
	//
	if(m_bsSimuSetting.at(ID_SPH_OUTPUT)){
		m_pPS->OutputParticles(CreateFileName(m_strSphOutputHeader, "dat", stp, 5));
	}

	g_TimerFPS.Reset();
	g_TimerFPS.Start();
	RXTIMER_RESET;

	if(m_iCurrentStep >= m_iMoveStart) m_bMoveSolid = true;

	//
	// �ő̂𓮂���
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
			// ���]
			sdir *= -1;
		}
		
		svel = m_fMoveMaxVel;

		m_pPS->MoveSphereObstacle(0, sdir*svel*m_fDt);

	}

	//
	// �V�~�����[�V�����^�C���X�e�b�v��i�߂�
	//
	if(m_bsSimuSetting.at(ID_SPH_INPUT)){
		// �p�[�e�B�N���ʒu�̓���
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
	//�ǉ��F�F�V�~�����[�V�����^�C���X�e�b�v��i�߂�@�M����
	//
	if( m_bMode == MODE_HEAT )
	{
		RXTIMER_RESET;

		StepHT(m_fDt);						//�M����
		StepPS(m_fDt);
//		StepSM(m_fDt);
	}
	//
	//�ǉ��F�F�V�~�����[�V�����^�C���X�e�b�v��i�߂�@SM�@
	//
	else if( m_bMode == MODE_SM )
	{
		RXTIMER_RESET;

//		StepPS(m_fDt);
//		StepSM(m_fDt);						//SM�@�@���q�ʒu�̃t�B�[�h�o�b�N�C�v�Z�C���`���
	}
	//
	//�ǉ��F�F�V�~�����[�V�����^�C���X�e�b�v��i�߂�@�X����
	//
	else if( m_bMode == MODE_ICE )
	{
		StepPS(m_fDt);						//���q�@�̉^��
		RXTIMER_RESET;

		StepHT(m_fDt);						//�M����

		//���q�x�[�X
		StepSolid_Melt(m_fDt);				//�Z������
		StepSolid_Freeze(m_fDt);			//�Ìŏ���

		//StepCalcParam(m_fDt);				//���x�ɂ����`��ԌW������@���ԏ�Ԃ���
		StepCluster(m_fDt);					//�N���X�^�̌v�Z

		StepInterpolation(m_fDt);			//�t�̂ƌő̂̉^������`���
	}

	//
	//�ǉ��F�F���q�̐F�ݒ�
	//
	StepParticleColor();

	if(m_bsSimuSetting.at(ID_SPH_ANISOTROPIC)){
		// �ٕ����J�[�l��
		RXTIMER_RESET;
		m_pPS->CalAnisotropicKernel();
	}

//	RXTIMER("anisotropic");

	//
	// ���̕\�ʃ��b�V������
	//
	if(m_bsSimuSetting.at(ID_SPH_MESH)){
		CalMeshSPH(m_iMeshMaxN, m_fMeshThr);
	}
	RXTIMER("mesh mc");

	//
	// FPS�̌v�Z
	//
	g_TimerFPS.Stop();
	if(m_iCurrentStep) ComputeFPS();

	if(m_iCurrentStep%RX_DISPLAY_STEP == 0){
		string tstr;
		RXTIMER_STRING(tstr);
		RXCOUT << "time : " << tstr << endl;
	}

	//
	// ���b�V���ۑ�
	//
	if(m_bsSimuSetting.at(ID_SPH_MESH)){
		if(m_bsSimuSetting.at(ID_SPH_MESH_OUTPUT) && m_iCurrentStep%m_iSaveMeshSpacing){
			SaveMesh(m_iCurrentStep, m_Poly);
		}
	}

	// �ő�X�e�b�v���𒴂�����A�C�h�����~�߂�
	if(m_iCurrentStep > RX_MAX_STEPS) SwitchIdle(0);

	m_iCurrentStep++;		// ���݂̃X�e�b�v��

	redraw();
}

/*!
 * �^�C�}�[�R�[���o�b�N�֐�
 * @param[in] x ���[�U��`�ϐ�
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
 * �A�C�h���֐���ON/OFF
 * @param[in] on true��ON, false��OFF
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
 * �t���X�N���[��/�E�B���h�E�\���̐؂�ւ�
 * @param[in] win 1��GL�L�����p�X���E�B���h�E���Ńt����, 0�ŃE�B���h�E�܂߂ăt����
 */
void rxFlGLWindow::SwitchFullScreen(int win)
{
	static int pos0[2] = { 0, 0 };
	static int win0[2] = { 500, 500 };
	if(win){
		// �E�B���h�E���Ńt����
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
		// �E�B���h�E�܂߂ăt����
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
 * �t���X�N���[��/�E�B���h�E�\���̏�Ԏ擾
 */
int rxFlGLWindow::IsFullScreen(void)
{
	return (m_bFullscreen ? 1 : 0);
}


/*!
 * �t�@�C���ǂݍ���
 * @param[in] fn �t�@�C���p�X
 */
void rxFlGLWindow::OpenFile(const string &fn)
{
	redraw();
}

/*!
 * �t�@�C����������
 * @param[in] fn �t�@�C���p�X
 */
void rxFlGLWindow::SaveFile(const string &fn)
{
}



/*!
 * ���݂̉�ʕ`����摜�t�@�C���Ƃ��ĕۑ�
 * @param[in] fn �t�@�C���p�X
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

	// �㉺���]
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
 * ���݂̉�ʕ`����摜�t�@�C���Ƃ��ĕۑ�
 * @param[in] stp ���݂̃X�e�b�v��(�t�@�C�����Ƃ��Ďg�p)
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
 * ���b�V�����t�@�C���Ƃ��ĕۑ�
 * @param[in] fn �t�@�C����
 * @param[in] polys �|���S���I�u�W�F�N�g
 */
void rxFlGLWindow::SaveMesh(const string fn, rxPolygons &polys)
{
	rxPOV pov;
	if(pov.SaveListData(fn, polys)){
		RXCOUT << "saved the mesh to " << fn << endl;
	}
}

/*!
 * ���b�V�����t�@�C���Ƃ��ĕۑ�
 * @param[in] stp ���݂̃X�e�b�v��(�t�@�C�����Ƃ��Ďg�p)
 * @param[in] polys �|���S���I�u�W�F�N�g
 */
void rxFlGLWindow::SaveMesh(const int &stp, rxPolygons &polys)
{
	SaveMesh(CreateFileName(RX_DEFAULT_MESH_DIR+"sph_", "inc", stp, 5), polys);
}


/*!
 * �C�x���g�n���h��
 * @param[in] ev �C�x���gID
 */
int rxFlGLWindow::handle(int e)
{
	switch(e){
	case FL_DND_ENTER:
	case FL_DND_RELEASE:
	case FL_DND_LEAVE:
	case FL_DND_DRAG:
	case FL_PASTE:
		// DnDBox�Ƀy�[�X�g������n�����߂ɁC
		// �����̃C�x���g�������Ƃ���Fl_Gl_Window::handle�ɏ�����n���Ȃ��D
		return 1;

	default:
		break;
	}

	return Fl_Gl_Window::handle(e);
}


/*!
 * Fl_Menu_Bar�̃R�[���o�b�N�֐� - Draw
 * @param[in] widget �E�B�W�b�g�̐e�N���X�I�u�W�F�N�g
 * @param[in] x ���[�U��`�ϐ�
 */
void rxFlGLWindow::OnMenuDraw_s(Fl_Widget *widget, void* x)
{
	Fl_Menu_Bar *menubar = (Fl_Menu_Bar*)widget;
	char picked[80];
	menubar->item_pathname(picked, sizeof(picked)-1);	// ���j���[��

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
	case RXD_NORMAL:	// �@��
		m_pPS->ToggleNormalCalc((m_iDraw & RXD_NORMAL));
		break;

	case RXD_CELLS:	// �����Z��
		if(m_iDraw & RXD_CELLS){
			m_pPS->SetParticlesToCell();
			m_pPS->SetPolygonsToCell();
		}
		break;

	case RXD_REFRAC:	// ���ܕ`��
		// if(m_iDraw & RXD_REFRAC) m_iDraw |= RXD_FOAM;
		// m_bSphMesh = (m_iDraw & RXD_REFRAC);
		break;

	case RXD_MESH:	// ���b�V�������ƕ`��
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
 * Fl_Menu_Bar�̃R�[���o�b�N�֐� - Simulation
 * @param[in] widget �E�B�W�b�g�̐e�N���X�I�u�W�F�N�g
 * @param[in] x ���[�U��`�ϐ�
 */
void rxFlGLWindow::OnMenuSimulation_s(Fl_Widget *widget, void* x)
{
	Fl_Menu_Bar *menubar = (Fl_Menu_Bar*)widget;
	char picked[80];
	menubar->item_pathname(picked, sizeof(picked)-1);	// ���j���[��

	string label = picked;
	string menu_name = "Simulation/";
	label = label.substr(menu_name.size(), string::npos);

	((rxFlGLWindow*)x)->OnMenuSimulation(1.0, label);
}
void rxFlGLWindow::OnMenuSimulation(double val, string label)
{
	if(label.find("Reset") != string::npos){
		// �V�[�����Z�b�g
		if(m_pPS) delete m_pPS;
		m_pPS = 0;
		InitSPH(m_Scene);
		InitHT(m_Scene);
	}

	if(label.find("Wavelet Turbulence") != string::npos){	
		// �E�F�[�u���b�g����
		m_pPS->ToggleWaveletTurb();
		m_bsSimuSetting.set(ID_SPH_PS_TURB, m_pPS->IsWaveletTurb());
		RXCOUT << "Turb : " << (m_pPS->IsWaveletTurb() ? "on" : "off") << endl;
	}
	else if(label.find("SPS Turbulence") != string::npos){
		// �E�F�[�u���b�g����(Sub-Particle-Scale)
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
	//�ǉ�
	else if( label.find("SPH only") != string::npos )
	{
		m_bMode = MODE_SPH;
	}
	else if( label.find("HeatTransfar") != string::npos )
	{
		//�M�����̂�
		cout << "Change HeatTransfar" << endl;
//		RXCOUT << "SPS : " << (m_bsSimuSetting[ID_SPH_SPS_TURB] ? "on" : "off") << endl;
		m_bsSimuSetting.set(ID_HEAT, !m_bsSimuSetting.at(ID_HEAT));
		m_bsSimuSetting.set(ID_SM, false);
		m_bsSimuSetting.set(ID_ICE, false);
		m_bMode = MODE_HEAT;
	}
	//�ǉ�
	else if( label.find("ShapeMatching") != string::npos )
	{
		//SM�@�̂�
		cout << "Change ShapeMatching" << endl;
		m_bsSimuSetting.set(ID_HEAT, false);
		m_bsSimuSetting.set(ID_SM, !m_bsSimuSetting.at(ID_SM));
		m_bsSimuSetting.set(ID_ICE, false);
		m_bMode = MODE_SM;
	}
	//�ǉ�
	else if( label.find("IceStructure") != string::npos )
	{
		//�X�\��
		cout << "Change IceStructure" << endl;
		m_bsSimuSetting.set(ID_HEAT, false);
		m_bsSimuSetting.set(ID_SM, false);
		m_bsSimuSetting.set(ID_ICE, !m_bsSimuSetting.at(ID_ICE));
		m_bMode = MODE_ICE;
	}

#ifdef RX_USE_PBD
	else if(label.find("Artificial Pressure") != string::npos){
		// �N���X�^�����O��h�����߂̐l�H���͂�ON/OFF (Tensile Instability)
		static_cast<RXSPH*>(m_pPS)->GetArtificialPressure() ^= 1;
	}
#endif
	else if(label.find("Particle Data Input") != string::npos){
		// �p�[�e�B�N���f�[�^�̃C���v�b�g
		m_bsSimuSetting.flip(ID_SPH_INPUT);
	}
	else if(label.find("Particle Data Output") != string::npos){
		// �p�[�e�B�N���f�[�^�̃A�E�g�v�b�g
		m_bsSimuSetting.flip(ID_SPH_OUTPUT);
		if(m_bsSimuSetting.at(ID_SPH_OUTPUT)) m_pPS->OutputSetting(m_strSphOutputName0);
	}
	else if(label.find("Mesh Saving") != string::npos){
		// ���b�V���ۑ�	
		m_bsSimuSetting.flip(ID_SPH_MESH_OUTPUT);
	}
	else if(label.find("Image Saving") != string::npos){
		// ��ʂ̒���摜�ۑ�
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
 * Fl_Menu_Bar�̃R�[���o�b�N�֐� - Particle/Color/
 * @param[in] widget �E�B�W�b�g�̐e�N���X�I�u�W�F�N�g
 * @param[in] x ���[�U��`�ϐ�
 */
void rxFlGLWindow::OnMenuParticle_s(Fl_Widget *widget, void* x)
{
	Fl_Menu_Bar *menubar = (Fl_Menu_Bar*)widget;
	char picked[80];
	menubar->item_pathname(picked, sizeof(picked)-1);	// ���j���[��

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
 * Particle/���j���[�̃C�x���g�n���h��
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
 * Particle/Color/���j���[�̃C�x���g�n���h��
 * ���q�̐F��؂�ւ�
 */
void rxFlGLWindow::OnMenuParticleColor(double val, string label)
{
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
	//�ǉ��@���x
	else if(label.find("Temperature") != string::npos)
	{
		((RXSPH*)m_pPS)->DetectSurfaceParticles();
		m_pPS->SetColorType(rxParticleSystemBase::RX_TEMP);
		m_pPS->SetColorVBOFromArray(m_ht->getTemps(), 1, false, 1.5f * m_ht->getTempMax());
		m_iColorType = rxParticleSystemBase::RX_TEMP;
	}
	//�ǉ��@�X
	else if(label.find("Ice_Cnct") != string::npos)
	{
		((RXSPH*)m_pPS)->DetectSurfaceParticles();
		SetParticleColorType(rxParticleSystemBase::RX_ICE_CONNECT);
//		m_pPS->SetColorVBOFromArray(m_ht->getTemps(), 1, false, 1.5f * m_ht->getTempMax());
//		m_pPS->SetColorType(rxParticleSystemBase::RX_ICE_CONNECT);
		m_iColorType = rxParticleSystemBase::RX_ICE_CONNECT;

		//
		//�@�ǉ��F�N���X�^�̕\���؂�ւ�
		//�l�ʑ̃x�[�X��
		//m_iShowClusterIndx++;

		//if( (unsigned)m_iShowClusterIndx > m_sm_connects.size() )
		//{
		//	m_iShowClusterIndx = 0;
		//}
		//else if( m_iShowClusterIndx < 0 )
		//{
		//	m_iShowClusterIndx = 0;
		//}

		//cout << "ClusterCountUp :: m_iShowClusterIndx = " << m_iShowClusterIndx << endl;

		//if( m_iShowClusterIndx != m_sm_connects.size() && m_iShowClusterIndx >= 0)
		//{
		//	cout << "Cnct :: Indxes :: " << endl;
		//	for( int i = 0; i < m_sm_connects[m_iShowClusterIndx]->GetNumVertices(); i++ )
		//	{
		//		cout << "                  i = " << i << " pIndx = " << m_sm_connects[m_iShowClusterIndx]->GetParticleIndx(i) << endl;
		//	}
		//}

		//���q�x�[�X��
		m_iShowTetraIndx++;
		
		if( m_iShowTetraIndx > m_iTetraNum )
		{
			m_iShowTetraIndx = 0;
		}
		else if( m_iShowTetraIndx < 0 )
		{
			m_iShowTetraIndx = 0;
		}

		cout << "Tetrahedra :: m_iShowTetraIndx = " << m_iShowTetraIndx << endl;

		if( m_iShowTetraIndx != m_iTetraNum && m_iShowTetraIndx >= 0)
		{
			cout << "Tetrahedra :: Indxes :: " << endl;
			for( int i = 0; i < m_ice->GetTtoPIndx(m_iShowTetraIndx); i++ )
			{
				cout << "                  i = " << i << " pIndx = " << m_ice->GetTtoP(m_iShowTetraIndx, i) << endl;
			}
		}

		StepParticleColor();
	}
	else if(label.find("Ice_Calc") != string::npos)
	{
		((RXSPH*)m_pPS)->DetectSurfaceParticles();
		SetParticleColorType(rxParticleSystemBase::RX_ICE_CALC);
//		m_pPS->SetColorVBOFromArray(m_ht->getTemps(), 1, false, 1.5f * m_ht->getTempMax());
//		m_pPS->SetColorType(rxParticleSystemBase::RX_ICE_CALC);
		m_iColorType = rxParticleSystemBase::RX_ICE_CALC;
		//
		//�@�ǉ��F�N���X�^�̕\���؂�ւ�
		//�l�ʑ̃x�[�X��
		m_iShowClusterIndx++;

		if( m_iShowClusterIndx > m_iClusteresNum )
		{
			m_iShowClusterIndx = 0;
		}
		else if( m_iShowClusterIndx < 0 )
		{
			m_iShowClusterIndx = 0;
		}

		cout << "ClusterCountUp :: m_iShowClusterIndx = " << m_iShowClusterIndx << endl;

		if( m_iShowClusterIndx != m_iClusteresNum && m_iShowClusterIndx >= 0)
		{
			cout << "Cluster :: Indxes :: " << endl;
			for( int i = 0; i < m_sm_cluster[m_iShowClusterIndx]->GetNumVertices(); i++ )
			{
				cout << "                  i = " << i << " pIndx = " << m_sm_cluster[m_iShowClusterIndx]->GetParticleIndx(i) << endl;
			}
		}

		StepParticleColor();
	}

	//�ǉ��@�\�ʕX���q�@������
	else if(label.find("Edge") != string::npos){
		((RXSPH*)m_pPS)->DetectSurfaceParticles();
		SetParticleColorType(rxParticleSystemBase::RX_RAMP);
//		m_pPS->SetColorType(rxParticleSystemBase::RX_EDGE);
		m_iColorType = rxParticleSystemBase::RX_EDGE;
	}
	//�ǉ��@�������p�p�X
	else if(label.find("FAST_PATH") != string::npos){
		((RXSPH*)m_pPS)->DetectSurfaceParticles();
		SetParticleColorType(rxParticleSystemBase::RX_ICE_FAST_PATH);
		m_iColorType = rxParticleSystemBase::RX_ICE_FAST_PATH;
		StepParticleColor();
	}
	else if(label.find("Surface") != string::npos){
		((RXSPH*)m_pPS)->DetectSurfaceParticles();
		SetParticleColorType(rxParticleSystemBase::RX_SURFACE);
	}
	else if(label.find("None") != string::npos){
		SetParticleColorType(rxParticleSystemBase::RX_NONE);
	}
	
	RXCOUT << "color of particle : " << label << endl;
	m_pParent->UpdateMenuState();
	redraw();
}
/*!
 * Particle/Draw/���j���[�̃C�x���g�n���h��
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
 * Fl_Menu_Bar�̃R�[���o�b�N�֐� - Solid
 * @param[in] widget �E�B�W�b�g�̐e�N���X�I�u�W�F�N�g
 * @param[in] x ���[�U��`�ϐ�
 */
void rxFlGLWindow::OnMenuSolid_s(Fl_Widget *widget, void* x)
{
	Fl_Menu_Bar *menubar = (Fl_Menu_Bar*)widget;
	char picked[80];
	menubar->item_pathname(picked, sizeof(picked)-1);	// ���j���[��

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
	case RXS_MOVE:	// �ő̈ړ�
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
 * Fl_Menu_Bar�̃R�[���o�b�N�֐� - Mesh
 * @param[in] widget �E�B�W�b�g�̐e�N���X�I�u�W�F�N�g
 * @param[in] x ���[�U��`�ϐ�
 */
void rxFlGLWindow::OnMenuTriangulation_s(Fl_Widget *widget, void* x)
{
	Fl_Menu_Bar *menubar = (Fl_Menu_Bar*)widget;
	char picked[80];
	menubar->item_pathname(picked, sizeof(picked)-1);	// ���j���[��

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
 * Fl_Menu_Bar�̃R�[���o�b�N�֐� - Scene
 * @param[in] widget �E�B�W�b�g�̐e�N���X�I�u�W�F�N�g
 * @param[in] x ���[�U��`�ϐ�
 */
void rxFlGLWindow::OnMenuScene_s(Fl_Widget *widget, void* x)
{
	Fl_Menu_Bar *menubar = (Fl_Menu_Bar*)widget;
	char picked[80];
	menubar->item_pathname(picked, sizeof(picked)-1);	// ���j���[��

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
		InitHT(m_Scene);
	}
	
	m_pParent->UpdateMenuState();
	redraw();
}



//-----------------------------------------------------------------------------
// MARK:���b�V���쐬
//-----------------------------------------------------------------------------
/*!
 * �O�p�`���b�V������
 * @param[in] nmax ���b�V�����O���b�h�𑜓x(�ő�)
 * @param[in] thr ���b�V����臒l
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
 * ���b�V�����̏�����
 */
bool rxFlGLWindow::ResetMesh(void)
{
	// �|���S��������
	if(!m_Poly.vertices.empty()){
		m_Poly.vertices.clear();
		m_Poly.normals.clear();
		m_Poly.faces.clear();
		m_Poly.materials.clear();
	}
	m_iNumVrts = 0;
	m_iNumTris = 0;

	// ���_VBO
	if(!m_uVrtVBO) glDeleteBuffers(1, &m_uVrtVBO);
	m_uVrtVBO = 0;

	// �@��VBO
	if(!m_uNrmVBO) glDeleteBuffers(1, &m_uNrmVBO);
	m_uNrmVBO = 0;

	// ���b�V��VBO
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
 * ���b�V�������̃Z�b�g
 */
void rxFlGLWindow::SetMeshCreation(void)
{
	m_bsSimuSetting.set(ID_SPH_MESH, ((m_iDraw & RXD_MESH) ? true : false));
	if(m_bsSimuSetting.at(ID_SPH_MESH) && m_iCurrentStep){
		CalMeshSPH(m_iMeshMaxN, m_fMeshThr);
	}
}

/*!
 * �O�p�`���b�V���̐���(MC�@,CPU)
 * @param[in] nmax ���b�V�����O���b�h�𑜓x(�ő�)
 * @param[in] thr ���b�V����臒l
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

	// �|���S��������
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
 * �O�p�`���b�V���̐���(MC�@,GPU)
 * @param[in] nmax ���b�V�����O���b�h�𑜓x(�ő�)
 * @param[in] thr ���b�V����臒l
 */
bool rxFlGLWindow::calMeshSPH_GPU(int nmax, double thr)
{
	if(m_pPS == NULL) return false;

	Vec3 minp = m_vMeshBoundaryCen-m_vMeshBoundaryExt;
	Vec3 maxp = m_vMeshBoundaryCen+m_vMeshBoundaryExt;

	double h;	//�Z����
	int n[3];	//�e���̃Z����
	CalMeshDiv(minp, maxp, nmax, h, n, 0.05);
	for(int i = 0; i < 3; ++i) m_iMeshN[i] = n[i];
	//cout << "mc : " << n[0] << " x " << n[1] << " x " << n[2] << endl;

	m_iDimVBO = 4;

	if(!m_pMCMeshGPU){
		m_pMCMeshGPU = new rxMCMeshGPU;
		m_pMCMeshGPU->Set(minp, Vec3(h), n, m_iVertexStore);

		// ���b�V���f�[�^�i�[�p��VBO�̊m��
		AssignArrayBuffers(m_pMCMeshGPU->GetMaxVrts(), 4, m_uVrtVBO, m_uNrmVBO, m_uTriVBO);

		////�ǉ��F�X�̃t���O�p�f�[�^��GPU�Ɋm��
		//CuAllocateArray((void**)&m_fIceFlag, m_pPS->GetMaxParticles()*sizeof(RXREAL));
		//CuSetArrayValue((void*)m_fIceFlag, 0, m_pPS->GetMaxParticles()*sizeof(RXREAL));

		////�ǉ��F�ő̗p�f�[�^�̊m��
		//AssignArrayBuffers(m_pMCMeshGPU->GetMaxVrts(), 4, m_uVrtVBO_solid, m_uNrmVBO_solid, m_uTriVBO_solid);
	}
	
	if(!m_Poly.vertices.empty()){
		m_Poly.vertices.clear();
		m_Poly.normals.clear();
		m_Poly.faces.clear();
		m_Poly.materials.clear();
	}
	
	//cout << "check start" << endl;
	// �T���v�����O�{�����[���̌v�Z(�A�֐��l���i�[�����O���b�h�f�[�^)
	m_pPS->CalImplicitFieldDevice(n, minp, Vec3(h, h, h), m_pMCMeshGPU->GetSampleVolumeDevice());
	//cout << "check end" << endl;
	unsigned int nvrts = 0;
	unsigned int ntris = 0;
	m_pMCMeshGPU->CreateMeshV(minp, h, n, thr, nvrts, ntris);
	
	m_iNumVrts = nvrts;
	m_iNumTris = ntris;
	if(m_iCurrentStep%RX_DISPLAY_STEP == 0){
		RXCOUT << nvrts << " verts and " << ntris << " tri were created." << endl;
	}

	// FBO�Ƀf�[�^���R�s�[
	m_pMCMeshGPU->SetDataToFBO(m_uVrtVBO, m_uNrmVBO, m_uTriVBO);
	
/*
	//�ǉ��F�X�p����
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

	CuCopyArrayToDevice(m_fIceFlag, fIceCheck, 0, m_pPS->GetNumParticles()*sizeof(RXREAL));		//�X�̃t���O���̍X�V
	m_pPS->CalImplicitFieldDeviceSolid(n, minp, Vec3(h, h, h), m_pMCMeshGPU->GetSampleVolumeDevice(), m_fIceFlag);
	delete[] fIceCheck;

	nvrts = 0;
	ntris = 0;
	//������thr��臒l��ς�����
	m_pMCMeshGPU->CreateMeshV(minp, h, n, thr*1.2, nvrts, ntris);

	m_iNumVrts_solid = nvrts;
	m_iNumTris_solid = ntris;
	if(m_iCurrentStep%RX_DISPLAY_STEP == 0){
		RXCOUT << "iceMesh " << nvrts << " verts and " << ntris << " tri were created (solid)." << endl;
	}

	//// FBO�Ƀf�[�^���R�s�[
	m_pMCMeshGPU->SetDataToFBO(m_uVrtVBO_solid, m_uNrmVBO_solid, m_uTriVBO_solid);
	//�I���F�X�p���� DrawSolidSurface()�ŕ`��
*/
	// �t�@�C���ۑ��̂��߂Ƀ��b�V���f�[�^���z�X�g���z��ɃR�s�[
	if(m_bsSimuSetting.at(ID_SPH_MESH_OUTPUT)){
		m_pMCMeshGPU->SetDataToArray(m_Poly.vertices, m_Poly.normals, m_Poly.faces);
	}

	return true;
}


/*!
 * �O�p�`���b�V���̐���(Screen Space Mesh, CPU)
 * @param[in] nmax ���b�V�����O���b�h�𑜓x(�ő�)
 * @param[in] thr ���b�V����臒l
 */
bool rxFlGLWindow::calMeshSPH_SSM(int nmax, double thr)
{
	// MRK:calMeshSPH_SSM
	if(!m_pPS) return false;
	if(!m_pPS->GetNumParticles()) return false;

	int pnum = m_pPS->GetNumParticles();
	m_fPrtRad = m_pPS->GetParticleRadius()*1.5;
	m_fZmax = 10.0*m_fPrtRad;

	// SSM������
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

	// �|���S��������
	if(!m_Poly.vertices.empty()){
		m_Poly.vertices.clear();
		m_Poly.normals.clear();
		m_Poly.faces.clear();
		m_Poly.materials.clear();
	}

	// �p�[�e�B�N���f�[�^����
	RXREAL *data = 0;
	data = m_pPS->GetArrayVBO(rxParticleSystemBase::RX_POSITION);

	vector<Vec3> prts;
	prts.resize(pnum);
	int k = 0;
	for(int i = 0; i < pnum; ++i){
		prts[i] = Vec3(data[k], data[k+1], data[k+2]);
		k += 4;
	}

	// �X�N���[���X�y�[�X���b�V������
	vector< vector<int> > tris;
	m_pSSMCPU->CreateMesh(m_fProjectionMatrix, m_fModelviewMatrix, m_iWinW, m_iWinH, 
					   prts, pnum, m_Poly.vertices, m_Poly.normals, tris, 
					   m_iDepthFiltering, m_iSSDebugOutput);

	// �|���S������
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

	// �|���S���ގ�
	m_Poly.materials[m_matPoly.name] = m_matPoly;

	// VBO�ւ̈ړ�
	int nv = (int)m_Poly.vertices.size();
	int nm = (int)tris.size();
	int nn = nv;

	m_iNumVrts = nv;
	m_iNumTris = nm;

	m_iDimVBO = 3;
	AssignArrayBuffers(nv, 3, m_uVrtVBO, m_uNrmVBO, m_uTriVBO);
	SetFBOFromArray(m_uVrtVBO, m_uNrmVBO, m_uTriVBO, m_Poly.vertices, m_Poly.normals, m_Poly.faces);

	if(m_iCurrentStep%RX_DISPLAY_STEP == 0){
		RXCOUT << m_iNumVrts << " verts and " << m_iNumTris << " tri were created." << endl;
	}

	//int nx, ny;
	//RXREAL *dm;

	//m_pSSMCPU->GetDepthMapSize(nx, ny);
	//dm = m_pSSMCPU->GetDepthMap();

	//// �f�v�X�l�e�N�X�`������
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


	// �f�v�X�l�摜�o��
	//SaveTexture("depth.png", g_texDepth);

	return true;
}

/*!
 * �O�p�`���b�V���̐���(Screen Space Mesh, GPU)
 * @param[in] nmax ���b�V�����O���b�h�𑜓x(�ő�)
 * @param[in] thr ���b�V����臒l
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

	// SSM������
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

	// �|���S��������
	if(!m_Poly.vertices.empty()){
		m_Poly.vertices.clear();
		m_Poly.normals.clear();
		m_Poly.faces.clear();
		m_Poly.materials.clear();
	}

	// �p�[�e�B�N���f�[�^����
	RXREAL *dprts = 0;
	RXREAL *drads = 0;

	if(m_pPS->IsSubParticle()){
		// �L���ȃT�u�p�[�e�B�N�����X�g�̍쐬
		((RXSPH*)m_pPS)->CalRadiusAndRatio();

		// �p�[�e�B�N���f�[�^�̎擾
		dprts = ((RXSPH*)m_pPS)->GetSubParticlePosDev();		// ���W
		drads = ((RXSPH*)m_pPS)->GetSubParticleRadDev();		// ���a
		pnum  = ((RXSPH*)m_pPS)->GetNumValidSubParticles();	// �p�[�e�B�N����
	}
	else{
		dprts = m_pPS->GetParticleDevice();
	}

	// �X�N���[���X�y�[�X���b�V������
	m_pSSMGPU->CreateMeshVBO(m_fProjectionMatrix, m_fModelviewMatrix, m_iWinW, m_iWinH, dprts, drads, pnum, 4, 
								m_uVrtVBO, m_uNrmVBO, m_uTriVBO, m_iDepthFiltering, m_iSSDebugOutput);

	if(!m_pPS->IsSubParticle()){
		m_pPS->UnmapParticle();
	}

	if(m_iCurrentStep%RX_DISPLAY_STEP == 0){
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
// SPH�֐�
//-----------------------------------------------------------------------------
/*!
 * SPH�̏�����
 * @param[in] fn_scene �V�[���L�q�t�@�C����
 */
void rxFlGLWindow::InitSPH(rxSPHConfig &sph_scene)
{
	// SPH�N���X�̏������ƃV�[���N���X�ւ̐ݒ�
	if(m_pPS) delete m_pPS;
	m_pPS = new RXSPH(true);
	sph_scene.Clear();
	sph_scene.SetPS(m_pPS);
		
	// �V�[���S�̏��̓ǂݍ���
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

	// ���b�V���������E�̐ݒ�
	m_vMeshBoundaryCen = sph_scene.GetSphEnv().mesh_boundary_cen;
	m_vMeshBoundaryExt = sph_scene.GetSphEnv().mesh_boundary_ext;

	// �����Ɋւ���t���O�̐ݒ�
	m_pPS->ToggleWaveletTurb(m_bsSimuSetting.at(ID_SPH_PS_TURB));
	m_pPS->ToggleUseVorticity(m_bsSimuSetting.at(ID_SPH_VC));

	// �p�[�e�B�N���`��̂��߂̃J���[�^�C�v�̐ݒ�
	m_pPS->SetColorType(m_iColorType);
	m_pPS->SetColorVBO();

	// �V�[���̌ʐݒ�(�p�[�e�B�N���C�ő�)�̓ǂݍ���
	if(m_bsSimuSetting.at(ID_SPH_INPUT)){
		// �p�[�e�B�N�������t�@�C������ǂݍ���
		m_pPS->InputParticles(CreateFileName(m_strSphInputHeader, "dat", 0, 5));
		CalMeshSPH(m_iMeshMaxN, m_fMeshThr);
	}
	else{
		m_pPS->Reset(rxParticleSystemBase::RX_CONFIG_NONE);
		sph_scene.LoadSceneFromFile();
	}

	// �T�u�p�[�e�B�N���ݒ�
	((RXSPH*)m_pPS)->SetEtcri(g_fEtCri);
	((RXSPH*)m_pPS)->InitSubParticles(g_fEtCri);

	// �ٕ����J�[�l���ݒ�
	if(m_bsSimuSetting.at(ID_SPH_ANISOTROPIC)){
		m_pPS->CalAnisotropicKernel();
	}

	// �\�ʃ��b�V���̃��Z�b�g
	ResetMesh();


	// ���Ԍv���p�ϐ��̏�����
	g_fTotalTime = 0.0f;
	g_fAvgTime = 0.0f;
	g_iTimeCount = 0;
}

/*!
 * �V�[���ɐ��H��ǉ�
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
 * SPH�̃^�C���X�e�b�v��i�߂�
 * @param[in] dt �^�C���X�e�b�v��
 */
void rxFlGLWindow::StepPS(double dt)
{
	if(!m_bPause)
	{
		float damping = 1.0f;
		int iterations = 1;

		// �V�~�����[�V�����p�����[�^
		m_pPS->SetIterations(iterations);
		m_pPS->SetDamping(damping);
		m_pPS->SetGravity((RXREAL)(-m_fGravity));

		// �V�~�����[�V�����X�e�b�v��i�߂�
		m_pPS->Update((RXREAL)dt, m_iCurrentStep);

		//�ǉ��F���q�����`�F�b�N���āC�e�����X�V
		UpdateInfo();
	}
}

/*!
 * �ǉ����ꂽ���q�ɑ΂��鏈��
 * 
 */
void rxFlGLWindow::UpdateInfo()
{
	int beforeSize = m_fIntrps.size();
	int nowSize = m_pPS->GetNumParticles();

	if( beforeSize == 0 || beforeSize >= nowSize ) return;		//���q���ɕω����Ȃ��Ȃ�I��

	cout << __FUNCTION__ << " AddInfo beforeSize = " << beforeSize << " nowSize = " << nowSize << endl;

	//�M����
	m_ht->AddParticle( nowSize );

	//���`���
	for( int i = beforeSize; i < nowSize; i++ )
	{
		m_fIntrps.push_back( 1.0f-(m_ht->getTemps()[i] / m_ht->getTempMax()) );
	}

	//�Z���̏ꍇ�́C��؋Ìł��s���Ȃ��Ƃ��Ēǉ����Ȃ�
	//return;

	//�N���X�^
	rxSPHEnviroment sph_env = m_Scene.GetSphEnv();				// Shape Matching�̐ݒ�@�p�����[�^�ǂݍ���
	for(int i = beforeSize; i < nowSize; i++)
	{
		m_sm_cluster.push_back(new Ice_SM(i));
		m_sm_cluster[i]->SetSimulationSpace(-sph_env.boundary_ext, sph_env.boundary_ext);
		m_sm_cluster[i]->SetTimeStep(sph_env.smTimeStep);
		m_sm_cluster[i]->SetCollisionFunc(0);
		m_sm_cluster[i]->SetStiffness(1.0, 0.0);
		m_iClusteresNum++;
	}

	//�ő̍\��
	m_ice->SetParticleNum(nowSize);
	m_ice->SetClusterNum(nowSize);

	//�f�o�C�X�������̍X�V

}

//-----------------------------------------------------------------------------
// �M�����֐�
//-----------------------------------------------------------------------------
/*!
 * �M�����̏�����
 * @param[in] fn_scene �V�[���L�q�t�@�C����
 */
void rxFlGLWindow::InitHT(rxSPHConfig &sph_scene)
{	cout << __FUNCTION__ <<	endl;
	//	if( m_ht ) delete m_ht;
	m_ht = new HeatTransfar( m_Scene.GetSphEnv().max_particles );		//�ŏ��ɍő吔���m�ۂ��Ă����āC�g���͍̂쐬���ꂽ�p�[�e�B�N���܂łƂ���
	m_ht->setCarnelConstant( ((RXSPH*)m_pPS)->GetEffectiveRadius() );	//�J�[�l���֐��̒萔�̂��߂̏���
	m_ht->setNumVertices( ICENUM );					//�p�[�e�B�N���̐����擾

	//�t�@�C������p�����[�^�̓ǂݍ���
	rxSPHEnviroment sph_env = sph_scene.GetSphEnv();

	m_ht->setTimeStep(sph_env.htTimeStep);
	m_ht->setTempMax(sph_env.tempMax);
	m_ht->setTempMin(sph_env.tempMin);
	m_ht->setLatentHeat(sph_env.latentHeat);
	m_ht->setCffCntHt(sph_env.cffcntHt);
	m_ht->setCffCntTd(sph_env.cffcntTd);
	meltPIndx = 0;
	debugIndx = 0;
	ClearRect();
}
/*!
 * �M�����̃^�C���X�e�b�v��i�߂�
 * @param[in] dt �^�C���X�e�b�v��
 */
void rxFlGLWindow::StepHT(double dt)
{//	cout << "StepHT" << endl;

	//�M����
	RXREAL *d  = m_pPS->GetArrayVBO(rxParticleSystemBase::RX_DENSITY);		//�e���q�̖��x
	RXREAL *p  = m_pPS->GetArrayVBO(rxParticleSystemBase::RX_POSITION);		//�e���q�̈ʒu

	((RXSPH*)m_pPS)->DetectSurfaceParticles();								//�\�ʗ��q���o

	int *surfaceParticles = (int *)( ((RXSPH*)m_pPS)->GetArraySurf() );		//�\�ʗ��q
	vector<int> ids;														//�\�ʗ��q�̓Y����
	vector<float> dists;													//�\�ʗ��q�̋���

	//������
	m_ht->resetNeighborhoodsId();
	m_ht->resetNeighborhoodsDist();
	
	//�ߖT���q���擾
	vector<vector<rxNeigh>>& neights = ((RXSPH*)m_pPS)->GetNeights();

//	cout << __FUNCTION__ << "Step1" << endl;
//	RXTIMER("ht1");

	//�ߖT���q�̓Y�����Ƌ����̐ݒ�
	for( int i = 0; i < m_pPS->GetNumParticles(); i++ )
	{
//		if( m_fIntrps.size() <= (unsigned)i ) continue;			//�l�ʑ̃x�[�X�ł̂�

		//�\�ʗ��q������@�P�Ȃ�Ε\�ʗ��q�C�O�Ȃ�������q�@���̍ۂ̐��l�͋ߖT���q������\��
		//�ߖT���q�����ŕ\�ʐς��ߎ�
		if( surfaceParticles[i] == 1 )
		{
			//���ɋ߂��C���͂̍������q�́C�\�ʂł͂Ȃ���ʂƂ���
			double floor = -m_Scene.GetSphEnv().boundary_ext[1];		//���̍����@���̒l
//			cout << "d[" << i << "] = " << d[i] << endl;

			//���q���ɂ���ăp�����[�^��ς��Ȃ��Ƃ����Ȃ��D
			//�\�ʗ��q����ɕs�������̂ŏC�����K�v�D
			//0.75�ŉ��Q�i�Ƃ��
			if(p[i*4+1] < floor+((RXSPH*)m_pPS)->GetEffectiveRadius()*0.2)				//1331�@���P�i
			{
				if(d[i] < 950.0)
				{
					m_ht->setSurfaceParticleNums(i, (int)( neights[i].size() ));		//�\�ʈ���
				}
				else
				{
					m_ht->setSurfaceParticleNums(i, -1);								//��ʈ���
				}
			}
			else
			{
				m_ht->setSurfaceParticleNums(i, (int)( neights[i].size() ));
			}
		}
		else
		{
			m_ht->setSurfaceParticleNums(i, -1);
		}
		
		//������
		ids.clear();
		dists.clear();

		for( unsigned j = 0; j < neights[i].size(); j++)
		{
			if( i == (int)( neights[i][j].Idx ) ) continue;							//�������g���Ȃ�
			ids.push_back( (int)( neights[i][j].Idx ) );
			dists.push_back( (float)( neights[i][j].Dist ) );
		}

		m_ht->AddNeighborhoodsId( ids );
		m_ht->AddNeighborhoodsDist( dists );
	}
//	cout << __FUNCTION__ << "Step2" << endl;
//	RXTIMER("ht2");

	//�M�����v�Z
	m_ht->heatAirAndParticle(); 		 										//�M�����@��C�Ɨ��q
	m_ht->heatParticleAndParticle(d, ((RXSPH*)m_pPS)->GetEffectiveRadius());		//�M�����@���q��
	m_ht->calcTempAndHeat();														//�M�ʂ̉��x�ϊ��C���x�̔M�ʕϊ�

//	RXTIMER("ht3");
//	cout << __FUNCTION__ << "Step3" << endl;
}

//-----------------------------------------------------------------------------
// ShapeMatching�֐�
//-----------------------------------------------------------------------------
/*!
 * �Փˏ����֐�
 * @param[in] p ���݂̍��W
 * @param[out] np �Փˌ�̍��W
 * @param[in] v ���x
 * @param[in] obj �I�u�W�F�N�g�ԍ�
 */
void rxFlGLWindow::Collision(Vec3 &p, Vec3 &np, Vec3 &v, int obj)
{
	//// ���̂Ƃ̏Փ˔���
	//Vec3 rpos = p-g_v3Cen;
	//double d = norm(rpos)-g_fRad;
	//if(d < 0.0){
	//	Vec3 n = Unit(rpos);
	//	np = g_v3Cen+n*g_fRad;
	//}

	//// �ό`���b�V�����m�̏Փ�
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
//�@���̏��
//-------------------------------------------------------------------------------------------------------------
/*!
 * ���̏��̏�����
 */
void rxFlGLWindow::InitICE(void)
{	cout << __FUNCTION__ << endl;

	m_iLayer = m_Scene.GetSphEnv().layer;

	//m_ice = new IceStructure(5000, 5000, 26000);				//���q��4913�̏ꍇ�̃p�����[�^
	//m_ice = new IceStructure(2500, 2500, 12000);				//�ő嗱�q���@�ő�N���X�^���@�ő�l�ʑ̐�
	m_ice = new IceStructure(6000, 6000, 1);					//�\�ʗ��q�݂̂̏ꍇ
	m_ice->SetParticleNum(ICENUM);								//���q���̓o�^
}


//-------------------------------------------------------------------------------------------------------------
//�@�l�ʑ̏��
//-------------------------------------------------------------------------------------------------------------
/*!
 * �l�ʑ̏��̏�����
 */
void rxFlGLWindow::InitTetra()
{	cout << __FUNCTION__ << endl;
	MakeTetrahedra();											//tetgen�𗘗p�����l�ʑ̍쐬
//	MakeTetrahedraOnlySurface();								//�\�ʗ��q�ɂ��e�X�g��

//	Load_ELE_File("test.1.ele");								//ele�t�@�C����ǂݍ��݁C���X�g���쐬
		
	m_ice->SetTetraNum(m_vviTetraList.size());					//���l�ʑ̐���o�^

	//�e�l�ʑ̂Ɋ܂܂�闱�q���̃J�E���g
	for(unsigned i = 0; i < m_vviTetraList.size(); i++)
	{
		CountTetraHedra(i, m_vviTetraList[i]);
	}

	//�������m��
	m_ice->InitTetraInfo();

	//���q���������Ă���N���X�^���̔z����R�s�[
	int *PtoTNum = new int[ICENUM];

	for(int i = 0; i < ICENUM; i++)
	{
		PtoTNum[i] = m_ice->GetPtoTNum(i);
	}

	//�l�ʑ̃f�[�^�o�^
	for(unsigned i = 0; i < m_vviTetraList.size(); i++)
	{
		MakeTetraInfo(i, PtoTNum);
	}
	delete[] PtoTNum;

	//�ߖT�l�ʑ̃f�[�^�o�^
	for(unsigned i = 0; i < m_vviTetraList.size(); i++)
	{
		m_ice->SetNeighborTetra(i, m_iLayer);
	}

	m_iTetraNum = m_vviTetraList.size();
	m_iTetraNumNum = m_iTetraNum;		//�f�o�b�O�p

	//�f�o�b�O
	//DebugTetra();
}

/*!
 * �l�ʑ̏��̃f�o�b�O
 * @param[in]
 */
void rxFlGLWindow::DebugTetra()
{
	//�l�ʑ́����q
	for(unsigned i = 0; i < m_vviTetraList.size(); i++ )
	{
		m_ice->DebugTtoP(i);
	}

	//���q���l�ʑ�
	for(int i = 0; i < m_pPS->GetNumParticles(); i++)
	{
		m_ice->DebugPtoT(i);
	}

	//�ߖT�l�ʑ�
	for(unsigned i = 0; i < m_vviTetraList.size(); i++ )
	{
		m_ice->DebugNeighborTetra(i);
	}
}

/*!
 * �l�ʑ̂Ɋ܂܂�Ă��闱�q���̃J�E���g�C���q����������l�ʑ̐��̃J�E���g
 * @param[in] tNum �l�ʑ̔ԍ�
 */
void rxFlGLWindow::CountTetraHedra(int tIndx, vector<int>& pList)
{
	for(unsigned i = 0; i < pList.size(); i++)
	{
		int pIndx = pList[i];		
		m_ice->CountPtoT(pIndx);
		m_ice->CountTtoP(tIndx);

		//Indx�̍X�V
		if(m_ice->GetPtoTNum(pIndx) >= m_ice->GetPtoTIndx(pIndx))
		{
			m_ice->SetPtoTIndx(pIndx, m_ice->GetPtoTNum(pIndx));
		}
	}
	
	//Indx�̍X�V
	if(m_ice->GetTtoPNum(tIndx) >= m_ice->GetTtoPIndx(tIndx))
	{
		m_ice->SetTtoPIndx(tIndx, m_ice->GetTtoPNum(tIndx));
	}
}

/*!
 * �l�ʑ̂Ɋ܂܂�Ă�����̓o�^�@��������p�֐�
 * @param[in] tNum �l�ʑ̔ԍ�
 * @param[in] PtoTNum �l�ʑ̂������闱�q��
 */
void rxFlGLWindow::MakeTetraInfo(int tIndx, int* PtoTNum)
{
	//���q�������Ă���l�ʑ̂̔ԍ���o�^���邽�߂̏���
	//pCountList�ɂ́CtIndx�Ԗڂ̎l�ʑ̂Ɋ܂܂��e���q���C���ꂼ�ꂢ���̎l�ʑ̂ɑ����邩�����߂ĕۑ�����
	int* pCountList = new int[m_vviTetraList[tIndx].size()];

	for(int j = 0; j < m_vviTetraList[tIndx].size(); j++)
	{
		int pIndx = m_vviTetraList[tIndx][j];
		pCountList[j] = m_ice->GetPtoTNum(pIndx)-PtoTNum[pIndx];	//���q���������鉽�Ԗڂ̃N���X�^�Ȃ̂������߂�
		PtoTNum[pIndx]--;
	}

	//���q�Ǝl�ʑ̂̏��o�^
	vector<int>& pIndxList = m_vviTetraList[tIndx];

	for(int i = 0; i < m_ice->GetTtoPNum(tIndx); i++)
	{
		m_ice->SetPtoT(pIndxList[i], pCountList[i], tIndx, i);	//���q���������Ă���l�ʑ̂�o�^
	}

	m_ice->SetTtoP(tIndx, pIndxList);							//�l�ʑ̂��܂�ł��闱�q��o�^

	delete[] pCountList;
}

/*!
 * �l�ʑ̂Ɋ܂܂�Ă�����̓o�^
 * @param[in] tNum �l�ʑ̔ԍ�
 */
void rxFlGLWindow::MakeTetraInfo(int tIndx, vector<int> pList)
{//	cout << __FUNCTION__ << endl;
	
	CountTetraHedra(tIndx, pList);									//�e�l�ʑ̂Ɋ܂܂�闱�q���̃J�E���g�C�Y�����̍X�V

	//���q�Ǝl�ʑ̂̏��o�^
	for(unsigned i = 0; i < pList.size(); i++)
	{
		int freeIndx = m_ice->GetPtoTFreeIndx(pList[i]);
		m_ice->SetPtoT(pList[i], freeIndx, tIndx, i);				//���q���������Ă���l�ʑ̂�o�^
	}

	m_ice->SetTtoP(tIndx, pList);									//�l�ʑ̂��܂�ł��闱�q��o�^
	
	m_iTetraNum++;
	m_ice->SetTetraNum(m_iTetraNum);

	//�f�o�b�O
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
//�@�N���X�^���
//-------------------------------------------------------------------------------------------------------------
/*! 
 * �N���X�^�\���̏�����
 */
void rxFlGLWindow::InitCluster()
{	cout << __FUNCTION__ << endl;

	//�ϐ��̏�����
	for( vector<Ice_SM*>::iterator it = m_sm_cluster.begin(); it != m_sm_cluster.end(); ++it )
	{
		if(*it) delete *it;
	}

	//���q���ɒ�`������̂�p��
	for( int i = 0; i < ICENUM; i++)
	{	
		m_fIntrps.push_back( 1.0f-(m_ht->getTemps()[i] / m_ht->getTempMax()) );	//���`��ԌW��
	}
	
	//�N���X�^�쐬
	m_iClusteresNum = 0;														//�N���X�^��

	////�p�^�[���P�F�l�ʑ̃��X�g�����ɁC���q���ɃN���X�^�쐬
	//for(int i = 0; i < ICENUM; i++)
	//{
	//	// Shape Matching�̐ݒ�@�p�����[�^�ǂݍ���
	//	rxSPHEnviroment sph_env = m_Scene.GetSphEnv();

	//	//�N���X�^������
	//	m_sm_cluster.push_back(new Ice_SM(m_iClusteresNum));
	//	m_sm_cluster[m_iClusteresNum]->SetSimulationSpace(-sph_env.boundary_ext, sph_env.boundary_ext);
	//	m_sm_cluster[m_iClusteresNum]->SetTimeStep(sph_env.smTimeStep);
	//	m_sm_cluster[m_iClusteresNum]->SetCollisionFunc(0);
	//	m_sm_cluster[m_iClusteresNum]->SetStiffness(1.0, 0.0);

	//	//�l�ʑ̃��X�g�����ɁC���q���ɃN���X�^�쐬
	//	MakeCluster(i);

	//	m_iClusteresNum++;
	//}

	//�p�^�[���Q�F�ߖT���݂̂ŃN���X�^�쐬
	MakeClusterFromNeight();

	//TODO::���q���ʂ�������@���͂𐶂ނ���

	m_ice->SetClusterNum(m_iClusteresNum);				//�N���X�^���̓o�^
	m_iShowClusterIndx = m_iClusteresNum;

	//�f�o�b�O
	//�N���X�^�Ɋ܂܂�闱�q
	//for(int i = 0; i < ICENUM; i++)
	//{
	//	cout << "pIndx = " << endl;

	//	for(int j = 0; j < m_sm_cluster[i]->GetNumVertices(); j++)
	//	{
	//		cout << " j = " << m_sm_cluster[i]->GetParticleIndx(j);
	//	}
	//	cout << endl;
	//}
}

/*! 
 * �N���X�^�쐬
 * @param[in] pIndx ���q�ԍ�
 */
void rxFlGLWindow::MakeCluster(int pIndx)
{//	cout << __FUNCTION__ << " pIndx = " << pIndx << endl;
	RXREAL *p = m_pPS->GetArrayVBO(rxParticleSystemBase::RX_POSITION);

	//���闱�q���܂܂��l�ʑ̂��ߖT�N���X�^�Ƃ��C�e�l�ʑ̂Ɋ܂܂�闱�q�ŃN���X�^���쐬
	for(int i = 0; i < m_ice->GetPtoTIndx(pIndx); i++)
	{
		int itIndx = m_ice->GetPtoT(pIndx, i, 0);
		int ioIndx = m_ice->GetPtoT(pIndx, i, 1);

		if(itIndx == -1 || ioIndx == -1){ continue;	}

		//�l�ʑ̂̑S�Ă̗��q���N���X�^�ɓo�^
		for(int j = 0; j < m_ice->GetTtoPIndx(itIndx); j++)
		{
			int jpIndx = m_ice->GetTtoP(itIndx, j);

			if(jpIndx == -1){	continue;	}
			if(m_sm_cluster[pIndx]->CheckIndx(jpIndx)){	continue;	}

			int pNum = m_sm_cluster[pIndx]->GetNumVertices();

			m_sm_cluster[pIndx]->AddVertex( Vec3(p[jpIndx*4+0], p[jpIndx*4+1], p[jpIndx*4+2] ), 1.0, jpIndx);
			m_sm_cluster[pIndx]->SetAlphas(pNum, 1.0);
			m_sm_cluster[pIndx]->SetBetas (pNum, 0.0);
			m_sm_cluster[pIndx]->SetLayer (pNum, 0);
		}
	}

	//�ߖT�l�ʑ̂�layer���ǂ�C���q��ǉ����Ă���
	//TODO::�s����ɂȂ�Ȃ�Clayer�������ق�beta��������
	for(int i = 0; i < m_ice->GetPtoTIndx(pIndx); i++)
	{
		int itIndx = m_ice->GetPtoT(pIndx, i, 0);
		int ioIndx = m_ice->GetPtoT(pIndx, i, 1);

		if(itIndx == -1 || ioIndx == -1){	continue;	}

		for(int j = 0; j < m_ice->GetNTNum(itIndx); j++)
		{
			int jtIndx = m_ice->GetNeighborTetra(itIndx, j, 0);
			int jlIndx = m_ice->GetNeighborTetra(itIndx, j, 1);

			if(jtIndx == -1 || jlIndx == -1){	continue;	}

			//�l�ʑ̂̑S�Ă̗��q���N���X�^�ɓo�^
			for(int k = 0; k < m_ice->GetTtoPIndx(jtIndx); k++)
			{
				int kpIndx = m_ice->GetTtoP(jtIndx, k);

				if(kpIndx == -1){	continue;	}
				if(m_sm_cluster[pIndx]->CheckIndx(kpIndx)){	/*cout << "contineue kpIndx = " << kpIndx << endl;*/ continue;	}

				int pNum = m_sm_cluster[pIndx]->GetNumVertices();

				m_sm_cluster[pIndx]->AddVertex( Vec3(p[kpIndx*4+0], p[kpIndx*4+1], p[kpIndx*4+2] ), 1.0, kpIndx);
				m_sm_cluster[pIndx]->SetAlphas(pNum, 1.0);
				m_sm_cluster[pIndx]->SetBetas (pNum, 0.0);
				m_sm_cluster[pIndx]->SetLayer (pNum, jlIndx);
			}
		}
	}
}

/*!
 * �ߖT���݂̂ɂ��N���X�^�쐬
 * @param[in]
 */
void rxFlGLWindow::MakeClusterFromNeight()
{
	//�������̂��߂ɉe�����a���L�����Ă݂�
	float radius = ((RXSPH*)m_pPS)->GetEffectiveRadius();
	((RXSPH*)m_pPS)->SetEffectiveRadius(radius * 2.0f);
	StepPS(m_fDt);																//��x�^�C���X�e�b�v�����߂Ȃ��ƁC�ߖT���q���擾����Ȃ��݂���
	((RXSPH*)m_pPS)->SetEffectiveRadius(radius);
	
	vector<vector<rxNeigh>>& neights = ((RXSPH*)m_pPS)->GetNeights();			//�ߖT���q���擾
	RXREAL *p = m_pPS->GetArrayVBO(rxParticleSystemBase::RX_POSITION);

	rxSPHEnviroment sph_env = m_Scene.GetSphEnv();								// Shape Matching�̐ݒ�@�p�����[�^�ǂݍ���

	for(int i = 0; i < ICENUM; i++)
	{
		//�N���X�^������
		m_sm_cluster.push_back(new Ice_SM(m_iClusteresNum));
		m_sm_cluster[m_iClusteresNum]->SetSimulationSpace(-sph_env.boundary_ext, sph_env.boundary_ext);
		m_sm_cluster[m_iClusteresNum]->SetTimeStep(sph_env.smTimeStep);
		m_sm_cluster[m_iClusteresNum]->SetCollisionFunc(0);
		m_sm_cluster[m_iClusteresNum]->SetStiffness(1.0, 0.0);

		//�ߖT���q���N���X�^�ɒǉ�
		for(int id_np = 0; id_np < neights[i].size(); id_np++)
		{
			int np_pIndx = neights[i][id_np].Idx;
			int pNum = m_sm_cluster[m_iClusteresNum]->GetNumVertices();

			m_sm_cluster[m_iClusteresNum]->AddVertex( Vec3(p[np_pIndx*4+0], p[np_pIndx*4+1], p[np_pIndx*4+2] ), 1.0, np_pIndx);
			m_sm_cluster[m_iClusteresNum]->SetAlphas(pNum, 1.0);
			m_sm_cluster[m_iClusteresNum]->SetBetas (pNum, 0.0);
			m_sm_cluster[m_iClusteresNum]->SetLayer (pNum, 0);
		}

		m_iClusteresNum++;
	}
}

/*!
 * ShapeMatching�@�̃^�C���X�e�b�v��i�߂�
 * @param[in] dt �^�C���X�e�b�v��
 */
void rxFlGLWindow::StepCalcParam(double dt)
{
	if(m_bPause) return;

	//�M�����ɂ��ϓ��������x��,�e�p�����[�^�ɓK�p
//	for( int i = 0; i < m_pPS->GetNumParticles(); i++ )
//	{
//		if( m_fIntrps.size() <= (unsigned)i )							continue;			//���݂��Ȃ��ꍇ�͔�΂�
//		if( m_ht->getPhase(i) == 2 || m_ht->getPhase(i) == -2 )			continue;			//���ԏ�ԏo�Ȃ��ꍇ�͖߂�
//
//		for( int j = 0; j < m_ice->GetPtoCIndx_Connect(i); j++ )
//		{
//			int* coSet = m_ice->GetPtoC_Connect(i, j);
//			int cIndx = coSet[0];
//			int oIndx = coSet[1];
//			if(cIndx == -1 || oIndx == -1) continue;
//			
//			//���x�͂O�`�T�O�O�x�@�Q�T�O�ŕX
////			m_sm_connects[cIndx]-SetAlphas( oIndx, 1.0 - (m_ht->getTemps()[i] / m_ht->getTempMax()) );		//alpha�i�����j�@���x��
//			m_sm_connects[cIndx]->SetAlphas( oIndx, 1.0 - (m_ht->getHeats()[i]/m_ht->getLatentHeat()) );	//alpha�i�����j�@�M�ʔ�
//
//			//�l��̐���
//			if( m_sm_connects[cIndx]->GetAlphas(oIndx) > 1.0 )
//				m_sm_connects[cIndx]->SetAlphas( oIndx, 1.0 );
//
//			if( m_sm_connects[cIndx]->GetAlphas(oIndx) < 0.0 ) 
//				m_sm_connects[cIndx]->SetAlphas( oIndx, 0.0 );
//		}
//	}

	//���`��ԃp�����[�^
	for( int i = 0; i < m_pPS->GetNumParticles(); i++ )
	{
		//if( m_fIntrps.size() <= (unsigned)i )					continue;				//���݂��Ȃ��ꍇ�͖߂�@�l�ʑ̃x�[�X�̂�
		if( m_ice->GetPtoCNum(i) == 0)							continue;
		if( m_ht->getPhase(i) == 2 || m_ht->getPhase(i) == -2 ) continue;				//���ԏ�ԏo�Ȃ��ꍇ�͖߂�

		m_fIntrps[i] = 1.0f - (m_ht->getTemps()[i] / m_ht->getTempMax());				//���x�Ō���@�������̂ق������R
//		m_fIntrps[i] = 1.0f - (g_ht->getTemps()[i] / g_ht->getLatentHeat());			//�M�ʂŌ���

		//�l��̐���
		if( m_fIntrps[i] > 1.0f ) m_fIntrps[i] = 1.0f;
		if( m_fIntrps[i] < 0.0f ) m_fIntrps[i] = 0.0f;
	}
//	cout << __FUNCTION__ << "Step3" << endl;
}

/*!
 * ShapeMatching�@�̃^�C���X�e�b�v��i�߂�
 * @param[in] dt �^�C���X�e�b�v��
 */
void rxFlGLWindow::StepCluster(double dt)
{//	cout << __FUNCTION__ << endl;

	if(m_bPause){	return;}

	int j = 0;
	int jpIndx = 0, jlIndx = 0, jcIndx = 0, joIndx = 0;
	double shapeNum = 0.0;
	int* coSet = 0;
	Vec3 pos,vel,veltemp;
	
	RXREAL *p = m_pPS->GetArrayVBO(rxParticleSystemBase::RX_POSITION);
	RXREAL *v = m_pPS->GetArrayVBO(rxParticleSystemBase::RX_VELOCITY);

	//�N���X�^�̐��l�X�V�@�ʒu�E���x
	#pragma omp parallel
	{
	#pragma omp for private(j, jpIndx, jlIndx)
		for(int i = 0; i < m_iClusteresNum; i++)
		{
			if(m_ice->GetPtoCNum(i) == 0){	continue;	}

			for(j = 0; j < m_ice->GetCtoPIndx(i); j++)
			{
				jpIndx = m_ice->GetCtoP(i, j, 0);								//�ǂ̌ő̂���̗��q
				jlIndx = m_ice->GetCtoP(i, j, 1);								//���w�ڂ̗��q��
				
				if(jpIndx == -1 || jlIndx == -1){	continue;	}
	
				m_sm_cluster[i]->SetCurrentPos	( j, Vec3(p[jpIndx*4+0], p[jpIndx*4+1], p[jpIndx*4+2]) );	//������Ȃ����ƁC���q�̊֘A��ۂ��ʒu���ς��D
				m_sm_cluster[i]->SetVelocity	( j, Vec3(v[jpIndx*4+0], v[jpIndx*4+1], v[jpIndx*4+2]) );	//���x�͖��X�V	������X�V����ƁC�ςɂȂ�

	//			m_sm_cluster[i]->SetOriginalPos	( j, m_sm_connects[cIndx]->GetOriginalPos(oIndx) );			//������Ȃ����ƒǉ������܂�������
	//			m_sm_cluster[i]->SetGoalPos		( j, Vec3(p[jpIndx*4+0], p[jpIndx*4+1], p[jpIndx*4+2]) );
	//			m_sm_cluster[i]->SetNewPos		( j, Vec3(p[jpIndx*4+0], p[jpIndx*4+1], p[jpIndx*4+2]) );
	//			m_sm_cluster[i]->parames[j].alpha = m_sm_connects[cIndx]->parames[oIndx].alpha;
			}
		}
	}//#pragma omp parallel

	//�N���X�^�̉^������
	#pragma omp parallel
	{
	#pragma omp for
		for(int i = 0; i < m_iClusteresNum; i++)
		{	
			if(m_ice->GetPtoCNum(i) == 0){	continue;	}
			m_sm_cluster[i]->Update();											//�^���v�Z
		}
	}//#pragma omp parallel
}

/* 
 * �t�̉^���ƌő̉^���̐��`��ԂƔ��f
 */
void rxFlGLWindow::StepInterpolation(double dt)
{
	if(m_bPause){	return;}

	int j = 0;
	int jpIndx = 0, jlIndx = 0, jcIndx = 0, joIndx = 0;
	double shapeNum = 0.0;
	int* coSet = 0;
	Vec3 pos,vel,veltemp;
	
	RXREAL *p = m_pPS->GetArrayVBO(rxParticleSystemBase::RX_POSITION);
	RXREAL *v = m_pPS->GetArrayVBO(rxParticleSystemBase::RX_VELOCITY);

	//���`���
	#pragma omp parallel	//�e����
	{
	#pragma omp for private(j, coSet, jcIndx, joIndx, pos, vel, shapeNum, veltemp)
		for(int i = 0; i < m_pPS->GetNumParticles(); i++)
		{
			//���q�̑��x����
			veltemp = Vec3(v[i*4+0], v[i*4+1], v[i*4+2]);
			while(norm2(veltemp)>=50)
			{
				//cout << "vel = " << veltemp << endl;
				v[i*4+0] = v[i*4+0] * 0.75;
				v[i*4+1] = v[i*4+1] * 0.75;
				v[i*4+2] = v[i*4+2] * 0.75;
				veltemp = Vec3(v[i*4+0], v[i*4+1], v[i*4+2]);
			}

			if(m_ice->GetParticleNum() <= i){		continue;	}	//�Z���݂̂̎����̂Ƃ��ɕK�v�ɂȂ�D
			if(m_ice->GetPtoCNum(i) <= 0)	{		continue;	}

			//���ꂼ��̃x�N�g�������������ς��Ƃ�
			pos = Vec3(0.0, 0.0, 0.0);
			vel = Vec3(0.0, 0.0, 0.0);
			shapeNum = 0.0;											//�N���X�^�̐�

			//�l�̎擾�C����
			for(j = 0; j < m_ice->GetPtoCIndx(i); j++)
			{
				jcIndx = m_ice->GetPtoC(i, j, 0);
				joIndx = m_ice->GetPtoC(i, j, 1);

				if(jcIndx == -1 || joIndx == -1){	continue;	}

				pos += m_sm_cluster[jcIndx]->GetVertexPos(joIndx);
				vel += m_sm_cluster[jcIndx]->GetVertexVel(joIndx);

				shapeNum += 1.0;
			}

			//�N���X�^�̐��Ŋ���
			//�ǂ̃N���X�^�ɂ��܂܂�Ă��Ȃ��ꍇ�C�^����SPH�@�ɏ]��
			if(shapeNum != 0.0)
			{
				pos /= shapeNum;
				vel /= shapeNum;
			}		
			else
			{
				pos = Vec3(p[i*4+0], p[i*4+1], p[i*4+2]);
				vel = Vec3(v[i*4+0], v[i*4+1], v[i*4+2]);
			}

			//SPH�@��SM�@�ŋ��߂����x�ƈʒu����
			v[i*4+0] = vel[0] * m_fIntrps[i] + v[i*4+0] * (1-m_fIntrps[i]);
			v[i*4+1] = vel[1] * m_fIntrps[i] + v[i*4+1] * (1-m_fIntrps[i]);
			v[i*4+2] = vel[2] * m_fIntrps[i] + v[i*4+2] * (1-m_fIntrps[i]);

			p[i*4+0] = pos[0] * m_fIntrps[i] + p[i*4+0] * (1-m_fIntrps[i]);
			p[i*4+1] = pos[1] * m_fIntrps[i] + p[i*4+1] * (1-m_fIntrps[i]);
			p[i*4+2] = pos[2] * m_fIntrps[i] + p[i*4+2] * (1-m_fIntrps[i]);

			//���肵�Ȃ�����
			//p[i*4+0] += v[i*4*0] * dt;
			//p[i*4+1] += v[i*4*1] * dt;
			//p[i*4+2] += v[i*4*2] * dt;
		}
	}//end #pragma omp parallel

	//SPH�̃f�[�^�̍X�V�@�ʒu�E���x
	m_pPS->SetArrayVBO(rxParticleSystemBase::RX_POSITION, p, 0, m_pPS->GetNumParticles());
	m_pPS->SetArrayVBO(rxParticleSystemBase::RX_VELOCITY, v, 0, m_pPS->GetNumParticles());
	
	p = m_pPS->GetArrayVBO(rxParticleSystemBase::RX_POSITION);
	v = m_pPS->GetArrayVBO(rxParticleSystemBase::RX_VELOCITY);
}

//------------------------------------------------------------------------------------------------------------
// �ő̍\���֐�
//------------------------------------------------------------------------------------------------------------
/*! 
 * �ő̍\���̏�����
 */
void rxFlGLWindow::InitICE_Cluster()
{//	cout << __FUNCTION__ << endl;

	//�J�E���g
	for(int i = 0; i < ICENUM; i++)
	{
		CountSolid(i);
	}

	//�������m��
	m_ice->InitClusterInfo();

	//���q���������Ă���N���X�^���̔z����R�s�[
	int *PtoCNum = new int[ICENUM];

	for(int i = 0; i < ICENUM; i++)
	{
		PtoCNum[i] = m_ice->GetPtoCNum(i);
	}

	//�N���X�^�Ɨ��q�̊֘A���̓o�^
	for(int i = 0; i < ICENUM; i++)
	{
		MakeClusterInfo(i, PtoCNum);	//�J�E���g��O�ōs���Ă��邽�߁C��������g��
	}

	delete[] PtoCNum;

	//TODO::�N���X�^�Ǝl�ʑ̂̊֘A���̓o�^

	RXREAL *p = m_pPS->GetArrayVBO(rxParticleSystemBase::RX_POSITION);

	m_ice->InitPath(p, m_sm_cluster, ICENUM);			//�������̂��߂̃p�X�쐬



	//�f�o�b�O
//	for(int i = 0; i < m_iClusteresNum; i++)
//	{
//		m_ice->DebugCtoP(i);
//	}
//
//	for(int i = 0; i < m_pPS->GetNumParticles(); i++)
//	{
//		m_ice->DebugPtoC(i);
//	}
}

/*!
 * �N���X�^�Ɋ܂܂�Ă��闱�q���̃J�E���g�C���q����������N���X�^���̃J�E���g
 * @param[in] cIndx �N���X�^�ԍ�
 */
void rxFlGLWindow::CountSolid(int cIndx)
{//	cout << __FUNCTION__ << " start cIndx = " << cIndx << endl;
	for(int j = 0; j < m_sm_cluster[cIndx]->GetNumVertices(); j++)
	{
		int jpIndx = m_sm_cluster[cIndx]->GetParticleIndx(j);
		m_ice->CountPtoC(jpIndx);										//���q���ڑ��N���X�^�ɏ���������̃J�E���g
		m_ice->CountCtoP(cIndx);										//�ڑ��N���X�^�����q���܂ތ��̃J�E���g
	
		//Indx�̍X�V
		if(m_ice->GetPtoCNum(jpIndx) >= m_ice->GetPtoCIndx(jpIndx))
		{
			m_ice->SetPtoCIndx(jpIndx, m_ice->GetPtoCNum(jpIndx));
		}
	}

	//Indx�̍X�V
	if(m_ice->GetCtoPNum(cIndx) >= m_ice->GetCtoPIndx(cIndx))
	{
		m_ice->SetCtoPIndx(cIndx, m_ice->GetCtoPNum(cIndx));
	}
//	cout << __FUNCTION__ << " end" << endl;
}

/*!
 * �N���X�^�Ɋ܂܂�Ă�����̓o�^�@��������p�֐�
 * @param[in] tNum �N���X�^�ԍ�
 */
void rxFlGLWindow::MakeClusterInfo(int cIndx, int* PtoCNum)
{
	//���q�������Ă���l�ʑ̂̔ԍ���o�^���邽�߂̏���
	//pCountList�ɂ́CcIndx�Ԗڂ̃N���X�^�Ɋ܂܂��e���q���C���ꂼ�ꂢ���̃N���X�^�ɑ����邩�����߂ĕۑ�����
	int* pCountList = new int[m_sm_cluster[cIndx]->GetNumVertices()];
	int* pLayerList = new int[m_sm_cluster[cIndx]->GetNumVertices()];			//���q�̏������C���[
	
	for(int i = 0; i < m_sm_cluster[cIndx]->GetNumVertices(); i++)
	{
		int pIndx = m_sm_cluster[cIndx]->GetParticleIndx(i);

		//��������z�肵�Ă��Ȃ����߁C����ď㏑�����Ă��܂��Ă���
		//-1��T�����ď㏑������̂ɐ؂�ւ���D
		for(int j = 0; j < m_ice->GetPtoCMax(); j++)
		{
			if(m_ice->GetPtoC(pIndx, j, 0) != -1 || m_ice->GetPtoC(pIndx, j, 1) != -1){	continue;	}

			if(j >= m_ice->GetPtoCIndx(pIndx))
			{
				m_ice->SetPtoCIndx(pIndx, j+1);								//���݂�Indx���傫���Ȃ�X�V
			}

			pCountList[i] = j;
			break;
		}

		pLayerList[i] = m_sm_cluster[cIndx]->GetLayer(i);					//���q�����w�ڂ̋ߖT�Ȃ̂����擾
	}

	//���q�ƃN���X�^�̏��o�^
	const vector<int>& pIndxList = m_sm_cluster[cIndx]->GetVertexIndxList();

	for(int i = 0; i < m_ice->GetCtoPNum(cIndx); i++)
	{
		m_ice->SetPtoC(pIndxList[i], pCountList[i], cIndx, i, pLayerList[i]);	//���q���������Ă���l�ʑ̂�o�^
	}

	m_ice->SetCtoP(cIndx, pIndxList, pLayerList);								//�l�ʑ̂��܂�ł��闱�q��o�^

	delete[] pCountList;
	delete[] pLayerList;
}

/*!
 * �N���X�^�Ɋ܂܂�Ă�����̓o�^
 * @param[in] tNum �N���X�^�ԍ�
 */
void rxFlGLWindow::MakeClusterInfo(int cIndx)
{	//cout << __FUNCTION__ << " start cIndx = " << cIndx << endl;
	//�J�E���g
	CountSolid(cIndx);

	//���q�ƃN���X�^�̏��o�^
	const vector<int>& pIndxList = m_sm_cluster[cIndx]->GetVertexIndxList();
	int* pLayerList = new int[m_sm_cluster[cIndx]->GetNumVertices()];			//���q�̏������C���[

	//layer�̂��߂̃R�s�[
	for(int i = 0; i < m_sm_cluster[cIndx]->GetNumVertices(); i++)
	{
		pLayerList[i] = m_sm_cluster[cIndx]->GetLayer(i);					//���q�����w�ڂ̋ߖT�Ȃ̂����擾
	}

	for(int i = 0; i < m_ice->GetCtoPNum(cIndx); i++)
	{
		int freeIndx = m_ice->GetPtoCFreeIndx(pIndxList[i]);
		m_ice->SetPtoC(pIndxList[i], freeIndx, cIndx, i, pLayerList[i]);		//���q���������Ă���l�ʑ̂�o�^
	}
	
	m_ice->SetCtoP(cIndx, pIndxList, pLayerList);								//�l�ʑ̂��܂�ł��闱�q��o�^

	delete[] pLayerList;
}

/*!
 * �ő̗̂Z������
 * @param dt �^�C���X�e�b�v
 */
void rxFlGLWindow::StepSolid_Melt(double dt)
{//	cout << __FUNCTION__ << endl;

	vector<int> viParticleList;												//�Z���������q�W��
	vector<int> viClusterList;												//�Ē�`����N���X�^�̏W��
	vector<int> viCLayerList;												//�Ē�`����N���X�^�̃��C���[
	vector<int> viTetraList;												//�Ē�`����l�ʑ̂̏W��
	vector<int> viTLayerList;												//�Ē�`����l�ʑ̂̃��C���[

	SearchMeltParticle(viParticleList);										//�Z�𗱎q�̒T��
	//RXTIMER("SearchReconstructTetra_Melt start");
	SearchReconstructTetra_Melt(viParticleList, viTetraList, viTLayerList);	//�Ē�`�l�ʑ̂̒T��
	//RXTIMER("SearchReconstructTetra_Melt end");
	SearchReconstructCluster_Melt(viParticleList, viClusterList, viCLayerList);	//�Ē�`�N���X�^�̒T��

	UpdateInfo_Melt_PandT(viParticleList);									//���q�E�l�ʑ̏��̍X�V
	//RXTIMER("UpdateInfo_Melt_PandC start");
	UpdateInfo_Melt_PandC(viParticleList, viClusterList);					//���q�E�N���X�^���̍X�V
	//RXTIMER("UpdateInfo_Melt_PandC end");
	
	//CheckDeleteCluster();													//����C��܊֌W�ɂ���N���X�^���폜
	//RXTIMER("CheckDeleteTetra start");
	//CheckDeleteTetra(viTetraList, viTLayerList);							//����C��܊֌W�ɂ���l�ʑ̂��폜
	//RXTIMER("CheckDeleteTetra end");

	//RXTIMER("SetTetraInfo start");
	SetTetraInfo(viParticleList, viTetraList, viTLayerList);				//���q�E�ߖT�l�ʑ̏��̍Ē�`
	//RXTIMER("SetTetraInfo end");
	//RXTIMER("SetClusterInfo start");
	SetClusterInfo(viParticleList, viClusterList, viCLayerList);			//���q�E�N���X�^���̍Ē�`
	//RXTIMER("SetClusterInfo end");

	//�f�o�b�O
	//if(viParticleList.size() == 0){	return;	}
	//cout << "Debug" << __FUNCTION__ << endl;
	//cout << "viParticleList.size = " << viParticleList.size() << " ";
	//for(unsigned i = 0; i < viParticleList.size(); i++)
	//{
	//	cout << " " << viParticleList[i];
	//}
	//cout << endl;

	//cout << "viClusterList.size =  " << viClusterList.size() << " ";
	//for(unsigned i = 0; i < viClusterList.size(); i++)
	//{
	//	cout << " " << viClusterList[i];
	//}
	//cout << endl;

	//cout << "viCLayerList:: ";
	//for(unsigned i = 0; i < viCLayerList.size(); i++)
	//{
	//	cout << " " << viCLayerList[i];
	//}
	//cout << endl;

	//cout << "viTetraList.size = " << viTetraList.size() << " ";
	//for(unsigned i = 0; i < viTetraList.size(); i++)
	//{
	//	cout << " " << viTetraList[i];
	//}
	//cout << endl;

	//cout << "viTLayerList:: ";
	//for(unsigned i = 0; i < viTLayerList.size(); i++)
	//{
	//	cout << " " << viTLayerList[i];
	//}
	//cout << endl;

	////�N���X�^�����q
	//for(int i = 0; i < m_iClusteresNum; i++){	m_ice->DebugCtoP(i);	}

	//���q���N���X�^
	//for(int i = 0; i < m_pPS->GetNumParticles(); i++){	m_ice->DebugPtoC(i);	}

	//SM�N���X�^�Ɋ܂܂�闱�q�͋@�\�Ŋm�F�ł���
	//for(int i = 0; i < ICENUM; i++)
	//{
	//	cout << "pIndx = " << endl;

	//	for(int j = 0; j < m_sm_cluster[i]->GetNumVertices(); j++)
	//	{
	//		cout << " j = " << m_sm_cluster[i]->GetParticleIndx(j);
	//	}
	//	cout << endl;
	//}

	//�l�ʑ́����q�͋@�\�Ŋm�F�ł���

	//���q���l�ʑ�
	//for(int i = 0; i < m_pPS->GetNumParticles(); i++){	m_ice->DebugPtoT(i);	}

	//�ߖT�l�ʑ�
	//for(unsigned i = 0; i < m_vviTetraList.size(); i++ ){	m_ice->DebugNeighborTetra(i);	}
}

/*!
 * �Z�𗱎q�̒T��
 * @param pList ���q�ԍ����X�g
 */
void rxFlGLWindow::SearchMeltParticle(vector<int>& pList)
{//	cout << __FUNCTION__ << endl;
	for(int i = 0; i < m_pPS->GetNumParticles(); i++)
	{	
		if( m_ht->getPhaseChange(i) != 1 )				continue;	//���]�ڂ̏����𖞂����Ă��Ȃ��ꍇ�͖߂�
		if( m_ht->getPhase(i) != 2 )					continue;	//���ւƑ��]�ڂ��Ă��Ȃ��ꍇ�͖߂�
		if( m_ice->GetParticleNum() <= i)				continue;	//�Z���݂̂̎����̂Ƃ��ɕK�v�ɂȂ�D
		if( m_ice->GetPtoCNum(i) == 0 )					continue;	//�N���X�^�Ɋ܂܂�Ă���
		if( m_ice->GetPtoTNum(i) == 0 )					continue;	//�N���X�^�Ɋ܂܂�Ă���

//		if(pList.size() > 1){	break;	}							//�Z�𗱎q���̐���

		m_fIntrps[i] = 0.0f;										//���`��Ԃ����Ȃ�
		m_ht->setPhaseChange(i, 0);									//���]�ڂ��I��������Ƃ�`����
		pList.push_back(i);											//�Z�𗱎q�̋L�^
	}
}

/*!
 * �Ē�`����N���X�^�̒T��
 * @param[in] pList   �Z�𗱎q���X�g
 * @param[in] cList   �Ē�`����N���X�^�̎Q�ƃ��X�g
 * @param[in] lList   �Ē�`����N���X�^�̃��C���[�Q�ƃ��X�g
 */
void rxFlGLWindow::SearchReconstructCluster_Melt(const vector<int>& pList, vector<int>& cList, vector<int>& lList)
{
	if(pList.size() == 0){	return; }

	//�Z�𗱎q���������Ă����N���X�^���Ē�`�N���X�^�D
	for(unsigned i = 0; i < pList.size(); i++)
	{
		int ipIndx = pList[i];

		for(int j = 0; j < m_ice->GetPtoCIndx(ipIndx); j++)
		{
			int jcIndx = m_ice->GetPtoC(ipIndx, j, 0);
			int joIndx = m_ice->GetPtoC(ipIndx, j, 1);
			int jlIndx = m_ice->GetPtoC(ipIndx, j, 2);

			if( jcIndx == -1 || joIndx == -1)
			{
				continue;
			}

			vector<int>::iterator check = std::find(cList.begin(), cList.end(), jcIndx);

			//���Ɋ܂܂�Ă���̂Ȃ�Clayer���ׂď������ق���D�悷��
			if(check != cList.end())
			{
				int layerIndx = check-cList.begin();
				if(lList[layerIndx] > jlIndx)
				{
					lList[layerIndx] = jlIndx;
				}
				continue;
			}

			cList.push_back(jcIndx);
			lList.push_back(jlIndx);
		}
	}
}

/*!
 * �Ē�`����l�ʑ̂̒T��
 * @param[in] pList   �Z�𗱎q���X�g
 * @param[in] tList   �Ē�`����l�ʑ̂̃��X�g
 * @param[in] lList   �Ē�`����l�ʑ̂̃��C���[���X�g
 */
void rxFlGLWindow::SearchReconstructTetra_Melt(const vector<int>& pList, vector<int>& tList, vector<int>& lList)
{//	cout << __FUNCTION__ << endl;
	if(pList.size() == 0){	return; }

	//�P�C���q���܂܂�Ă����C����l�ʑ�
	//�Q�C�P�̎l�ʑ̂̋ߖT�l�ʑ�
	//�܂�̓N���X�^���\�������l�ʑ́@TODO::�o��������Ȃ̂ŁC�����̌v�Z�R�X�g��������ΏC���\
	
	m_ice->ResetTFlag(m_iTetraNum);							//�l�ʑ̒T���t���O�̏�����

	//�P�C�Z�𗱎q���܂܂�Ă����l�ʑ�
	for(unsigned i = 0; i < pList.size(); i++)
	{
		int ipIndx = pList[i];

		for(int j = 0; j < m_ice->GetPtoTIndx(ipIndx); j++)
		{
			if(m_ice->GetPtoT(ipIndx, j, 0) == -1
			|| m_ice->GetPtoT(ipIndx, j, 1) == -1)
			{
				continue;
			}

			//if(std::find(tList.begin(), tList.end(), jtlSet[0]) != tList.end()){ continue;	}
			if( m_ice->GetTFlag(m_ice->GetPtoT(ipIndx, j, 0)) )	{	continue;	}
			else												{	m_ice->SetTFlag(m_ice->GetPtoT(ipIndx, j, 0), true);	}

			tList.push_back(m_ice->GetPtoT(ipIndx, j, 0));
			lList.push_back(1);								//0��1���̔��f�͂ł��Ȃ��̂�1�ɍ��킹��D
		}
	}

	//�Q�C�P�̎l�ʑ̂̋ߖT�l�ʑ�
	int tetraNum = tList.size();
	for(int i = 0; i < tetraNum; i++)
	{
		int itIndx = tList[i];

		for(int j = 0; j < m_ice->GetNTNum(itIndx); j++)
		{
			int jtIndx = m_ice->GetNeighborTetra(itIndx, j, 0);
			int jlIndx = m_ice->GetNeighborTetra(itIndx, j, 1);

			//���Ɋ܂܂�Ă���̂Ȃ�Clayer���ׂď������ق���D�悷��
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
 * ���q�Z�����̗��q�E�l�ʑ̏��̍X�V
 * @param pList ���q�ԍ��Q�ƃ��X�g
 */
void rxFlGLWindow::UpdateInfo_Melt_PandT(const vector<int>& pList)
{//	cout << __FUNCTION__ << endl;
	if(pList.size() == 0){	return; }

	for(unsigned i = 0; i < pList.size(); i++)
	{
		int ipIndx = pList[i];

		for(int j = 0; j < m_ice->GetPtoTIndx(ipIndx); j++)
		{
			int tIndx = m_ice->GetPtoT(ipIndx, j, 0);
			int oIndx = m_ice->GetPtoT(ipIndx, j, 1);										//���̏ꍇ�͂��̂܂ܓY����
			
			if(tIndx == -1 || oIndx == -1){ continue;	}

			m_ice->DeleteTtoP(tIndx, oIndx);
		}

		m_ice->ClearPtoT(ipIndx);
	}
}

/*!
 * �l�ʑ̍폜���̗��q�E�l�ʑ̏��̍X�V
 * @param tList �l�ʑ̔ԍ��Q�ƃ��X�g
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
 * ���q�Z�����̗��q�E�N���X�^���̍X�V
 * @param pList �Z�𗱎q�z��
 * @param cList �Ē�`�N���X�^
 */
void rxFlGLWindow::UpdateInfo_Melt_PandC(const vector<int>& pList, const vector<int>& cList)
{//	cout << __FUNCTION__ << endl;
	if(pList.size() == 0){	return; }

	//���񏈗��ŗp����ϐ����܂Ƃ߂Ē�`
	int j= 0, k = 0;
	int icIndx = 0;
	int jpIndx = 0;
	int* coSet;

	//�Z���������q�̃N���X�^�̏����C���̗��q�����菜��
	#pragma omp parallel
	{
	#pragma omp for private(j,k,icIndx,jpIndx,coSet)
		for(int i = 0; i < pList.size(); i++)
		{
			icIndx = pList[i];
	
			for(j = 0; j < m_ice->GetCtoPIndx(icIndx); j++)
			{
				jpIndx = m_ice->GetCtoP(icIndx, j, 0);
				
				if(jpIndx == -1){	continue;	}

				for(k = 0; k < m_ice->GetPtoCIndx(jpIndx); k++)
				{					
					if(m_ice->GetPtoC(jpIndx, k, 0) == -1
					|| m_ice->GetPtoC(jpIndx, k, 1) == -1
					|| m_ice->GetPtoC(jpIndx, k, 0) != icIndx)
					{
						continue;
					}
	
	#pragma omp critical (DeletePtoC)	//TODO�F�F��ɃJ�E���g�����ق������񉻂ł��Ă悢
	{
					m_ice->DeletePtoC(jpIndx, k);
	}
					break;			//�����N���X�^�ɕ����������邱�Ƃ͖����̂ŁCbreak
				}
			}
		}
	}//end #pragma omp parallel

	//�Z���������q���N���X�^�����폜
	#pragma omp parallel
	{
	#pragma omp for
		for(int i = 0; i < pList.size(); i++)
		{
			m_ice->ClearPtoC(pList[i]);
		}
	}//end #pragma omp parallel

	//�Z�������N���X�^�����q�����폜�@���q�ԍ����N���X�^�ԍ�
	#pragma omp parallel
	{
	#pragma omp for
		for(int i = 0; i < pList.size(); i++)
		{
			m_ice->ClearCtoP(pList[i]);
			m_sm_cluster[pList[i]]->Clear();
		}
	}//end #pragma omp parallel

	//�Ē�`����N���X�^�Ɋ܂܂�闱�q�̏����N���X�^�����X�V
	//���ꂩ��Ē�`����N���X�^�̏����폜
	#pragma omp parallel
	{
	#pragma omp for private(j,k,icIndx,jpIndx,coSet)
		for(int i = 0; i < cList.size(); i++)
		{
			icIndx = cList[i];
	
			for(j = 0; j < m_ice->GetCtoPIndx(icIndx); j++)
			{
				jpIndx = m_ice->GetCtoP(icIndx, j, 0);
				
				if(jpIndx == -1){	continue;	}

				for(k = 0; k < m_ice->GetPtoCIndx(jpIndx); k++)
				{
					if(m_ice->GetPtoC(jpIndx, k, 0) == -1
					|| m_ice->GetPtoC(jpIndx, k, 1) == -1
					|| m_ice->GetPtoC(jpIndx, k, 0) != icIndx)
					{
						continue;
					}


	#pragma omp critical (DeletePtoC)	//TODO�F�F��ɃJ�E���g�����ق������񉻂ł��Ă悢
	{
					m_ice->DeletePtoC(jpIndx, k);
	}
					break;			//�����N���X�^�ɕ����������邱�Ƃ͖����̂ŁCbreak
				}
			}
		}
	}//end #pragma omp parallel
}

/*!
 * ����C��܊֌W�ɂ���l�ʑ̂��폜		//TODO::���̊֐���IceStructure�Ɏ�������ׂ�
 */
void rxFlGLWindow::CheckDeleteTetra(vector<int>& tList, vector<int>& lList)
{
	//cout << "before tList.size =  " << tList.size() << " m_iTetraNumNum = " << m_iTetraNumNum << endl;

	//�l�ʑ̂œ��e������C�܂��͕�܊֌W�ɂ�����̂��폜�D
	//TODO::���񉻂����ꂼ��ōs��
	for(vector<int>::iterator indx = tList.begin(); indx !=tList.end(); indx++)
	{
		if(*indx == -1){	continue;	}

		vector<int> deleteTList;									//�폜����l�ʑ̂̓Y����

		CheckSameTetra(*indx, tList, deleteTList);					//������e�̎l�ʑ̂̓Y�����擾
		CheckIncludeTetra(*indx, tList, deleteTList);				//����Ă���l�ʑ̂̓Y�����擾
		UpdateInfo_Delete_TandP(tList, deleteTList);				//�l�ʑ̂Ɋւ�����̍폜
		RemoveReconstructTetra(tList, lList, deleteTList);			//�폜�����l�ʑ̂��Ē�`����l�ʑ̂Ɋ܂܂��Ȃ��菜��
	}

	//-1�ł���l�ʑ̃f�[�^���폜
	tList.erase(std::remove(tList.begin(), tList.end(), -1), tList.end());
	lList.erase(std::remove(lList.begin(), lList.end(), -1), lList.end());

	//cout << "after tList.size =  " << tList.size() << " m_iTetraNumNum = " << m_iTetraNumNum << endl;
}

/*!
 * ����̎l�ʑ̂����݂��邩�̔���
 * @param[in] tIndx ���肷��l�ʑ̔ԍ�
 * @param[in] searchTList �T���ΏۂƂȂ�l�ʑ̃��X�g
 * @param[in] deleteList ���ʂ�ǉ�����z��
 */
void rxFlGLWindow::CheckSameTetra(int tIndx, const vector<int>& searchTList, vector<int>& deleteList)
{
	int tNum = m_ice->GetTtoPNum(tIndx);

	if(tNum == 4 || tNum == 1){	return;	}								//�������q���S�Ȃ瓯��͑��݂��Ȃ��i�͂��j

	//������e��T��
	for(unsigned i = 0; i < searchTList.size(); i++)
	{
		int itIndx = searchTList[i];
		if(itIndx < 0)						{	continue;	}
		if(tIndx == itIndx)					{	continue;	}
		if(m_ice->GetTtoPNum(itIndx) == 0)	{	continue;	}
		if(m_ice->GetTtoPNum(itIndx) == 4)	{	continue;	}
		if(tNum != m_ice->GetTtoPNum(itIndx)){	continue;	}	//�������闱�q���̃`�F�b�N

		//tIndx��i�Ԗڂ̎l�ʑ̂̓��e���C���S�Ɉ�v���Ă��邩�̃`�F�b�N
		bool check = false;
		for(int j = 0; j < m_ice->GetTtoPIndx(tIndx); j++)
		{
			check = false;

			int jtIndx = m_ice->GetTtoP(tIndx, j);
			if(jtIndx == -1){	continue;	}

			//j�Ԗڂ̗��q�����ʂ��Ă��邩���`�F�b�N
			for(int k = 0; k < m_ice->GetTtoPIndx(itIndx); k++)
			{
				int ktIndx = m_ice->GetTtoP(itIndx, k);
				if(ktIndx == -1){	continue;	}

				if(jtIndx == ktIndx){	check = true; break;}
				else				{	check = false;		}
			}

			if(!check){	break;	}				//1�ł��܂�ł��Ȃ��Ȃ�s��v true�Ȃ玟�̗��q�����
		}

		//�Ō�܂�true�Ȃ瓯����e�D
		if(check)
		{	
			deleteList.push_back(i);
		}
	}

	//�f�o�b�O
	//if(tList.size() == 0){	return;	}
	//cout << "Debug::" << __FUNCTION__ << endl;
	//m_ice->DebugTtoP(tIndx);
	//for(int i = 0; i < tList.size(); i++)
	//{
	//	m_ice->DebugTtoP(tList[i]);
	//}
}

/*!
 * �����l�ʑ̂����݂��邩�̔���
 * @param[in] tIndx ���肷��l�ʑ̔ԍ�
 * @param[in] searchTList�@�T������l�ʑ̃��X�g
 * @param[in] deleteList ���ʂ�ǉ�����z��
 */
void rxFlGLWindow::CheckIncludeTetra(int tIndx, const vector<int>& searchTList, vector<int>& deleteList)
{
	int tNum = m_ice->GetTtoPNum(tIndx);

	if(tNum == 1){	return;	}								//�������q��1�Ȃ����͑��݂��Ȃ��i�͂��j

	//�����l�ʑ̂��Ē�`����l�ʑ̂���T��
	for(unsigned i = 0; i < searchTList.size(); i++)
	{
		int itIndx = searchTList[i];
		if(itIndx < 0)						{	continue;	}
		if(tIndx == itIndx)					{	continue;	}
		if(m_ice->GetTtoPNum(itIndx) == 0)	{	continue;	}
		if(m_ice->GetTtoPNum(itIndx) == 4)	{	continue;	}
		if(tNum <= m_ice->GetTtoPNum(itIndx)){	continue;	}	//�������闱�q���̃`�F�b�N

		//i�Ԗڂ̎l�ʑ̂̓��e��tIndx�Ԗڂ̎l�ʑ̂ɑS�Ċ܂܂�Ă��邩�̃`�F�b�N
		bool check = false;
		for(int j = 0; j < m_ice->GetTtoPIndx(itIndx); j++)
		{
			check = false;

			int jtIndx = m_ice->GetTtoP(itIndx, j);
			if(jtIndx == -1){	continue;	}

			//j�Ԗڂ̗��q�����ʂ��Ă��邩���`�F�b�N
			for(int k = 0; k < m_ice->GetTtoPIndx(tIndx); k++)
			{
				int ktIndx = m_ice->GetTtoP(tIndx, k);
				if(ktIndx == -1){	continue;	}

				if(jtIndx == ktIndx){	check = true; break;}
				else				{	check = false;		}
			}

			if(!check){	break;	}				//1�ł��܂�ł��Ȃ��Ȃ�s��v true�Ȃ玟�̗��q�����
		}

		//�Ō�܂�true�Ȃ�܂�ł���D
		if(check)
		{
			deleteList.push_back(i);
		}		
	}

	//�f�o�b�O
	//if(tList.size() == 0){	return;	}
	//cout << "Debug::" << __FUNCTION__ << endl;
	//m_ice->DebugTtoP(tIndx);
	//for(int i = 0; i < tList.size(); i++)
	//{
	//	m_ice->DebugTtoP(tList[i]);
	//}
}

/*
 * �Ē�`����l�ʑ̂���s�v�Ƃ݂Ȃ��ꂽ�l�ʑ̂�-1�ɏ���������
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
 * ���q�E�N���X�^���̍Ē�`
 * @param pList ���q���X�g
 * @param cList �Ē�`�N���X�^���X�g
 * @param lList ���C���[���X�g
 */
void rxFlGLWindow::SetClusterInfo(const vector<int>& pList, const vector<int>& cList, const vector<int>& lList)
{//	cout << __FUNCTION__ << endl;
	if(pList.size() == 0 || cList.size() == 0){	return;	}
	
	vector<int> checkTList;
	bool check = false;
	int j = 0, k = 0, l = 0;
	int icIndx = 0;
	int jtIndx = 0, joIndx = 0;
	int kpIndx = 0, ktIndx = 0, klIndx = 0;
	int lpIndx = 0;
	int pNum = 0;

	RXREAL *p = m_pPS->GetArrayVBO(rxParticleSystemBase::RX_POSITION);

	//TODO::�܂��́ClList��p���Ȃ��ł���Ă݂�D�����炭�v�Z���Ԃ͂���قǂ�����Ȃ��͂������c

	//SM�@�N���X�^�̏�����
	#pragma omp parallel
	{
	#pragma omp for
		for(int i = 0; i < cList.size(); i++)
		{
			if(std::find(pList.begin(), pList.end(), cList[i]) != pList.end()){	continue;	}
			m_sm_cluster[cList[i]]->Clear();
		}
	}//end #pragma omp parallel

	//SM�@�N���X�^�̍Ē�`
	#pragma omp parallel
	{
	#pragma omp for private(checkTList, check, j, k, l, icIndx, jtIndx, joIndx, kpIndx, ktIndx, klIndx, lpIndx, pNum)
		for(int i = 0; i < cList.size(); i++)
		{
			checkTList.clear();
			icIndx = cList[i];
			if(std::find(pList.begin(), pList.end(), icIndx) != pList.end()){	continue;	}

			//�N���X�^���Ē�`����ہC��{�ƂȂ闱�q��������l�ʑ̂��珉�����q�𓾂�D
			//�N���X�^�ԍ��������q�ԍ��Ȃ̂ɒ���
			//�ȉ����֐��ɂ���ƁC�G���[���o�Ă��܂������Ȃ�
			for(j = 0; j < m_ice->GetPtoTIndx(icIndx); j++)
			{
				jtIndx = m_ice->GetPtoT(icIndx, j, 0);
				joIndx = m_ice->GetPtoT(icIndx, j, 1);

				if(jtIndx == -1 || joIndx == -1){ continue;	}
				if(std::find(checkTList.begin(), checkTList.end(), jtIndx) != checkTList.end())
				{	continue;	}
				else
				{	checkTList.push_back(jtIndx);	}

				//�l�ʑ̂̑S�Ă̗��q���N���X�^�ɓo�^
				for(k = 0; k < m_ice->GetTtoPIndx(jtIndx); k++)
				{
					kpIndx = m_ice->GetTtoP(jtIndx, k);

					if(kpIndx == -1){	continue;	}
					if(m_sm_cluster[icIndx]->CheckIndx(kpIndx)){	continue;	}

					pNum = m_sm_cluster[icIndx]->GetNumVertices();

					m_sm_cluster[icIndx]->AddVertex(Vec3(p[kpIndx*4+0], p[kpIndx*4+1], p[kpIndx*4+2]), 1.0, kpIndx);
					m_sm_cluster[icIndx]->SetAlphas(pNum, 1.0);
					m_sm_cluster[icIndx]->SetBetas (pNum, 0.0);
					m_sm_cluster[icIndx]->SetLayer (pNum, 0);
				}
			}

			//�ߖT�l�ʑ̂�layer���ǂ�C���q��ǉ����Ă���
			//TODO::�s����ɂȂ�Ȃ�Clayer�������ق�beta��������
			for(j = 0; j < m_ice->GetPtoTIndx(icIndx); j++)
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

					//�l�ʑ̂̑S�Ă̗��q���N���X�^�ɓo�^
					for(l = 0; l < m_ice->GetTtoPIndx(ktIndx); l++)
					{
						lpIndx = m_ice->GetTtoP(ktIndx, l);

						if(lpIndx == -1){	continue;	}
						if(m_sm_cluster[icIndx]->CheckIndx(lpIndx)){	/*cout << "contineue kpIndx = " << kpIndx << endl;*/ continue;	}

						pNum = m_sm_cluster[icIndx]->GetNumVertices();

						m_sm_cluster[icIndx]->AddVertex( Vec3(p[lpIndx*4+0], p[lpIndx*4+1], p[lpIndx*4+2] ), 1.0, lpIndx);
						m_sm_cluster[icIndx]->SetAlphas(pNum, 1.0);
						m_sm_cluster[icIndx]->SetBetas (pNum, 0.0);
						m_sm_cluster[icIndx]->SetLayer (pNum, klIndx);
					}
				}
			}
		}
	}//end #pragma omp parallel

	//�ő̏��̏�����
	//#pragma omp parallel
	//{
	//#pragma omp for
		for(int i = 0; i < cList.size(); i++)
		{
			if(std::find(pList.begin(), pList.end(), cList[i]) != pList.end())
			{
				continue;
			}
			m_ice->ClearCtoP(cList[i]);
		}
	//}//end #pragma omp parallel

	//�ő̏��̓o�^
	for(unsigned i = 0; i < cList.size(); i++)
	{
		if(std::find(pList.begin(), pList.end(), cList[i]) != pList.end())
		{
			continue;	
		}
		MakeClusterInfo(cList[i]);
	}
}

/*!
 * ���q�Z�����̗��q�E�l�ʑ̏��̍Ē�`
 * @param pList �Z�𗱎q���X�g
 * @param tList �Ē�`�l�ʑ̃��X�g
 * @param lList ���C���[���X�g
 */
void rxFlGLWindow::SetTetraInfo(const vector<int>& pList, const vector<int>& tList, const vector<int>& lList)
{//	cout << __FUNCTION__ << endl;
	if(pList.size() == 0 || tList.size() == 0){	return;	}

	int itIndx = 0;
	int ilayer = 0;

	#pragma omp parallel
	{
	#pragma omp for private(itIndx,ilayer)
		//�ߖT�l�ʑ̂̍Ē�`
		for(int i = 0; i < tList.size(); i++)
		{
			itIndx = tList[i];
			ilayer = lList[i];

			m_ice->ClearNeighborTetraFromLayer(itIndx, ilayer);
			m_ice->SetNeighborTetraFromLayer(itIndx, m_iLayer, ilayer);	//���������ɏd��
		}
	}//end #pragma omp parallel
}

/*!
 * ����C��܊֌W�ɂ���N���X�^���폜
 */
void rxFlGLWindow::CheckDeleteCluster()
{
}

/*!
 * �Ìŏ���
 * @param dt �^�C���X�e�b�v
 */
void rxFlGLWindow::StepSolid_Freeze(double dt)
{
	vector<int> viParticleList;														//�Ìł������q�W��
	vector<int> viClusterList;														//�Ē�`����N���X�^�̏W��
	vector<int> viCLayerList;														//�Ē�`����N���X�^�̃��C���[
	vector<int> viTetraList;														//�Ē�`����l�ʑ̂̏W��
	vector<int> viTLayerList;														//�Ē�`����l�ʑ̂̃��C���[
	
	SearchFreezeParticle(viParticleList);											//�Ìŗ��q�̒T��
	SetFreezeTetraInfo(viParticleList);												//�Ìŗ��q�Ɋւ���l�ʑ̂̍쐬
	SetFreezeClusterInfo(viParticleList);											//�Ìŗ��q�Ɋւ���N���X�^�̍쐬

	SearchReconstructTetra_Freeze(viParticleList, viTetraList, viTLayerList);		//�Ē�`�l�ʑ̂̒T��
	SearchReconstructCluster_Freeze(viParticleList, viClusterList, viCLayerList);	//�Ē�`�N���X�^�̒T��

	//CheckDeleteCluster();															//����C��܊֌W�ɂ���N���X�^���폜
	//CheckDeleteTetra(viTetraList, viTLayerList);									//����C��܊֌W�ɂ���l�ʑ̂��폜

	SetTetraInfo(viParticleList, viTetraList, viTLayerList);						//���q�E�ߖT�l�ʑ̏��̍Ē�`
	SetClusterInfo(viParticleList, viClusterList, viCLayerList);					//���q�E�N���X�^���̍Ē�`

	//�f�o�b�O
	if(viParticleList.size() == 0 || viClusterList.size() == 0){	return;	}
	cout << "Debug " << __FUNCTION__ << "viParticleList.size = " << viParticleList.size() << " " << endl;
	//for(unsigned i = 0; i < viParticleList.size(); i++)
	//{
	//	cout << " " << viParticleList[i];
	//}
	//cout << endl;

	//cout << "viClusterList.size =  " << viClusterList.size() << " ";
	//for(unsigned i = 0; i < viClusterList.size(); i++)
	//{
	//	cout << " " << viClusterList[i];
	//}
	//cout << endl;

	//cout << "viCLayerList:: ";
	//for(unsigned i = 0; i < viCLayerList.size(); i++)
	//{
	//	cout << " " << viCLayerList[i];
	//}
	//cout << endl;

	//cout << "viTetraList.size = " << viTetraList.size() << " ";
	//for(unsigned i = 0; i < viTetraList.size(); i++)
	//{
	//	cout << " " << viTetraList[i];
	//}
	//cout << endl;

	//cout << "viTLayerList:: "3
	//for(unsigned i = 0; i < viTLayerList.size(); i++)
	//{
	//	cout << " " << viTLayerList[i];
	//}
	//cout << endl;

	////�N���X�^�����q
	//for(int i = 0; i < m_iClusteresNum; i++){	m_ice->DebugCtoP(i);	}

	//���q���N���X�^
	//for(int i = 0; i < m_pPS->GetNumParticles(); i++){	m_ice->DebugPtoC(i);	}

	//SM�N���X�^�Ɋ܂܂�闱�q�͋@�\�Ŋm�F�ł���
	//for(int i = 0; i < ICENUM; i++)
	//{
	//	cout << "pIndx = " << endl;

	//	for(int j = 0; j < m_sm_cluster[i]->GetNumVertices(); j++)
	//	{
	//		cout << " j = " << m_sm_cluster[i]->GetParticleIndx(j);
	//	}
	//	cout << endl;
	//}

	//�l�ʑ́����q�͋@�\�Ŋm�F�ł���

	//���q���l�ʑ�
	//for(int i = 0; i < m_pPS->GetNumParticles(); i++){	m_ice->DebugPtoT(i);	}

	//�ߖT�l�ʑ�
	//for(unsigned i = 0; i < m_vviTetraList.size(); i++ ){	m_ice->DebugNeighborTetra(i);	}
}

/*!
 * �Ìŗ��q�̒T��
 * @param pList ���q�ԍ����X�g
 */
void rxFlGLWindow::SearchFreezeParticle(vector<int>& pList)
{//	cout << __FUNCTION__ << endl;
	for(int i = 0; i < m_pPS->GetNumParticles(); i++)
	{	
		if(m_ht->getPhaseChange(i) != 1)				continue;	//���]�ڂ̏����𖞂����Ă��Ȃ��ꍇ�͖߂�
		if(m_ht->getPhase(i) != -2)						continue;	//�X�ւƑ��]�ڂ��Ă��Ȃ��ꍇ�͖߂�
		if(m_ice->GetParticleNum() <= i)				continue;	//�Z���݂̂̎����̂Ƃ��ɕK�v�ɂȂ�D
		if(m_ice->GetPtoCNum(i) != 0)					continue;	//�N���X�^�Ɋ܂܂�Ă���
		if(m_ice->GetPtoTNum(i) != 0)					continue;	//�N���X�^�Ɋ܂܂�Ă���
//		if(pList.size() > 1){	break;	}							//�Ìŗ��q���̐���
		
		pList.push_back(i);											//�Ìŗ��q�̋L�^
	}
}

/*!
 * �Ìŗ��q�Ɋւ�����̍쐬�@�l�ʑ�
 */
void rxFlGLWindow::SetFreezeTetraInfo(const vector<int>& pList)
{	//cout << __FUNCTION__ << " start" << endl;
	vector<int> tList;

	MakeFreezeTetrahedra(pList, tList);								//�l�ʑ̂̍쐬�@�ő̗��q�Ƃ̋Ì�
	//TODO::�ȉ��͂Q�͖�����
	AddFreezeTetrahedra(pList, tList);								//�l�ʑ̂֒ǉ��@�������q�����R�ȉ��̎l�ʑ̂ւ̒ǉ�
	MakeFreezeTetrahedra_OnlyFreezeParticle(pList, tList);			//�l�ʑ̂̍쐬�@�Ìŗ��q�݂̂ł̋Ì�

	for(unsigned i = 0; i < tList.size(); i++)
	{
		m_ice->SetNeighborTetra(tList[i], m_iLayer);				//�ߖT�l�ʑ̏��̐ݒ�
		//�f�o�b�O
		//m_ice->DebugNeighborTetra(tList[i]);
	}

	//cout << __FUNCTION__ << " end" << endl;
}

/*!
 * �Ìŗ��q�Ɋւ�����̍쐬�@�N���X�^
 */
void rxFlGLWindow::SetFreezeClusterInfo(const vector<int>& pList)
{	//cout << __FUNCTION__ << " start pList.size() = " << pList.size() << endl;

	for(unsigned i = 0; i < pList.size(); i++)
	{
		int pIndx = pList[i];
		if(m_ice->GetPtoTNum(pIndx) == 0){	continue;	}	//�l�ʑ̂ɑ����Ȃ��ꍇ�͖߂�

		//pIndx == cIndxa �Ȃ̂ɒ���
		MakeCluster(pIndx);
		MakeClusterInfo(pIndx);
	}
}

/*!
 * �ÌŌ�C�Ē�`����N���X�^�̒T��
 * @param[in] pList   �Ìŗ��q���X�g
 * @param[in] cList   �Ē�`����N���X�^�̎Q�ƃ��X�g
 * @param[in] lList   �Ē�`����N���X�^�̃��C���[�Q�ƃ��X�g
 */
void rxFlGLWindow::SearchReconstructCluster_Freeze(const vector<int>& pList, vector<int>& cList, vector<int>& lList)
{
	if(pList.size() == 0){	return; }
	//cout << __FUNCTION__ << endl;

	//�Ìŗ��q����������N���X�^�Ɋ܂܂�Ă���e���q�̃N���X�^���C�Ē�`�N���X�^�D
	for(unsigned i = 0; i < pList.size(); i++)
	{
		int ipIndx = pList[i];
		if(m_ice->GetPtoCNum(ipIndx) == 0){		continue;	}

		//�Ìŗ��q���܂܂�Ă���N���X�^���擾
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
		}
	}
}

/*!
 * �ÌŌ�C�ߖT���Ē�`����l�ʑ̂̒T��
 * @param[in] pList   �Ìŗ��q���X�g
 * @param[in] tList   �Ē�`����l�ʑ̂̃��X�g
 * @param[in] lList   �Ē�`����l�ʑ̂̃��C���[���X�g
 */
void rxFlGLWindow::SearchReconstructTetra_Freeze(const vector<int>& pList, vector<int>& tList, vector<int>& lList)
{	
	if(pList.size() == 0){	return; }

	//����V�������ꂽ�l�ʑ̌QA�̋ߖT�l�ʑ̌Q��A'�Ƃ���ƁCA'�̋ߖT�l�ʑ̂�A��ǉ�����K�v������D
	//�������ǉ�����ۂ�layer�������łȂ��Ƃ�����肪����̂ŁC�v�Z���Ԃ͂����邩������Ȃ����P���ɍĒ�`����D
	
	//�Ìŗ��q����������l�ʑ̂́C�ߖT�l�ʑ̑S�Ă��Ē�`����D

	m_ice->ResetTFlag(m_iTetraNum);							//�l�ʑ̒T���t���O�̏�����

	//�P�C�Ìŗ��q���܂܂�Ă����l�ʑ�
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
			else							{	m_ice->SetTFlag(jtIndx, true);	}

			tList.push_back(jtIndx);
			lList.push_back(1);								//0��1���̔��f�͂ł��Ȃ��̂�1�ɍ��킹��D�������Q�l�ɁD
		}
	}

	//�Q�C�P�̎l�ʑ̂̋ߖT�l�ʑ�
	int tetraNum = tList.size();
	for(int i = 0; i < tetraNum; i++)
	{
		int itIndx = tList[i];

		for(int j = 0; j < m_ice->GetNTNum(itIndx); j++)
		{
			int jtIndx = m_ice->GetNeighborTetra(itIndx, j, 0);
			int jlIndx = m_ice->GetNeighborTetra(itIndx, j, 1);

			//���Ɋ܂܂�Ă���̂Ȃ�Clayer���ׂď������ق���D�悷��
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
 * �ŋߖT�ŒZ���q�̎擾�@layer = 0�Ȃ�ő�27���x
 * @param[in] pIndx�@���q�ԍ�
 * @param[in] layer�@�ߖT�w��
 * @param[in] num�@�@�ő�擾���q��
 * @param[out] vector<int> �ŋߖT�ŒZ���q�̏W���@�����̋߂����Ƀ\�[�g����Ă���
 */
vector<int> rxFlGLWindow::GetNeighborDistanceIceParticles(int pIndx, int layer, int num)
{//	cout << __FUNCTION__ << " start" << endl;
	vector<int> pList;
	vector<float> distanceList;
	RXREAL *p = m_pPS->GetArrayVBO(rxParticleSystemBase::RX_POSITION);

	//�l�ʑ̃x�[�X
	//���闱�q���܂܂�Ă���ڑ��N���X�^�S�Ă�T������D
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
	//		if( check != pList.end() ){	continue;	}	//���łɊ܂�ł�����߂�

	//		float distance = dot(Vec3(p[pIndx*DIM+0],  p[pIndx*DIM+1],  p[pIndx*DIM+2]),
	//							 Vec3(p[jpIndx*DIM+0], p[jpIndx*DIM+1], p[jpIndx*DIM+2]));
	//		distanceList.push_back(distance);
	//		pList.push_back(jpIndx);

	//		if( pList.size() >= (unsigned)num ){ break;}
	//	}
	//	if( pList.size() >= (unsigned)num ){ break;}
	//}

	//���q�x�[�X
	//���闱�q���܂܂�Ă���l�ʑ̑S�Ă�T������
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
			if( check != pList.end() ){	continue;	}	//���łɊ܂�ł�����߂�

			float distance = dot(Vec3(p[pIndx*DIM+0],  p[pIndx*DIM+1],  p[pIndx*DIM+2]),
								 Vec3(p[jpIndx*DIM+0], p[jpIndx*DIM+1], p[jpIndx*DIM+2]));
			distanceList.push_back(distance);
			pList.push_back(jpIndx);

			if( pList.size() >= (unsigned)num ){ break;}
		}
		if( pList.size() >= (unsigned)num ){ break;}
	}

	//���闱�q���܂܂�Ă���N���X�^�́C�ߖT�N���X�^�����w���T������Dlayer��T���D
	for( int i = 0; i < layer; i++ )
	{
	}

	//�����̋߂����ɗ��q�ԍ����\�[�g
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
 * ���q�̐F�ݒ�@���x�C�M�ʁC�ڑ��N���X�^�C�v�Z�N���X�^��؂�ւ���D
 * 
 */
void rxFlGLWindow::StepParticleColor()
{
	//���x
	if( m_iColorType == rxParticleSystemBase::RX_TEMP )
	{
		//�^�����ł͌����ɂ��������̂ŁC������Ɛ��̂��f�t�H���g�Ƃ����D
		float* tempColor = new float[m_pPS->GetNumParticles()];
		for( int i = 0; i < m_pPS->GetNumParticles(); i++ )
		{
			//if( m_fIntrps.size() <= (unsigned)i )
			//if(m_ice->GetPtoCNum(i) == 0 || m_ice->GetPtoTNum(i) == 0)
			if(false)
			{
				tempColor[i] = 100.0f;
			}
			else
			{
				tempColor[i] = m_ht->getTemps()[i] + 100.0f;
			}
		}
		m_pPS->SetColorVBOFromArray( tempColor, 1, false, 1.5f * m_ht->getTempMax() );
		delete[] tempColor;
	}
	//SM�@
	else if( m_iColorType == rxParticleSystemBase::RX_SHPMTCHNG )
	{

	}
	//�ڑ��N���X�^
	else if( m_iColorType == rxParticleSystemBase::RX_ICE_CONNECT )
	{
		//�N���X�^�ɂ��ڑ��֌W�ɂ��闱�q�ɐF������
		//�z��̑S�v�f�����������Ȃ��ƁC�`�悪���������Ȃ�̂ɒ��ӁD
		float* tempColor = new float[m_pPS->GetNumParticles()];
		for( int i = 0; i < m_pPS->GetNumParticles(); i++ )
		{
//			if( m_fIntrps.size() <= (unsigned)i )
			//if(m_ice->GetPtoCNum(i) == 0 || m_ice->GetPtoTNum(i) == 0)
			if(false)
			{
				tempColor[i] = 0.0f;
			}
			else
			{
				tempColor[i] = 100.0f;
			}
		}

		//�l�ʑ̃x�[�X��
/*		if( m_iShowClusterIndx == m_sm_connects.size() )
		{	//�N���X�^�S�̂�\��
			for( unsigned i = 0; i < m_sm_connects.size(); i++ )
			{	
				for( int j = 0; j < m_sm_connects[i]->GetNumVertices(); j++ )
				{
					int jpIndx = m_sm_connects[i]->GetParticleIndx(j);
					tempColor[jpIndx] = 600.0f;
				}
			}
		}
		else
		{	//�P�̃N���X�^�݂̂�\��
			for( int j = 0; j < m_sm_connects[m_iShowClusterIndx]->GetNumVertices(); j++ )
			{	
				int jpIndx = m_sm_connects[m_iShowClusterIndx]->GetParticleIndx(j);
				tempColor[jpIndx] = 600.0f;
			}
		}
*/

		//���q�x�[�X��
		if(m_iShowTetraIndx == m_iTetraNum)
		{	//�l�ʑ̑S�̂�\��
			for(int i = 0; i < m_iTetraNum; i++)
			{	
				for(int j = 0; j < m_ice->GetTtoPIndx(i); j++)
				{
					int jpIndx = m_ice->GetTtoP(i, j);

					if(jpIndx == -1){	continue;	}
					tempColor[jpIndx] = 600.0f;
				}
			}
		}
		else
		{	//�P�̎l�ʑ݂̂̂�\��
			for( int j = 0; j < m_ice->GetTtoPIndx(m_iShowTetraIndx); j++ )
			{	
				int jpIndx = m_ice->GetTtoP(m_iShowTetraIndx, j);
				if(jpIndx == -1){	continue;	}
				tempColor[jpIndx] = 600.0f;
			}
		}

		m_pPS->SetColorVBOFromArray( tempColor, 1, false, 1.5f * m_ht->getTempMax() );
		delete[] tempColor;
	}
	//�v�Z�N���X�^
	else if( m_iColorType == rxParticleSystemBase::RX_ICE_CALC )
	{
		//�N���X�^�ɂ��ڑ��֌W�ɂ��闱�q�ɐF������
		//�z��̑S�v�f�����������Ȃ��ƁC�`�悪���������Ȃ�̂ɒ��ӁD
		float* tempColor = new float[m_pPS->GetNumParticles()];
		for( int i = 0; i < m_pPS->GetNumParticles(); i++ )
		{
			//if( m_fIntrps.size() <= (unsigned)i )
			if(m_ice->GetPtoCNum(i) == 0)
			{
				tempColor[i] = 0.0f;
			}
			else
			{
				tempColor[i] = 100.0f;
			}
		}

/*		//�l�ʑ̃x�[�X��
		if( m_iShowClusterIndx == m_iClusteresNum )
		{	//�N���X�^�S�̂�\��
			for( unsigned i = 0; i < m_sm_calcs.size(); i++ )
			{
				for( int j = 0; j < m_sm_calcs[i]->GetNumVertices(); j++ )
				{
					int jpIndx = m_sm_calcs[i]->GetParticleIndx(j);
					tempColor[jpIndx] = 700.0f;
				}
			}
		}
		else
		{	//�P�̃N���X�^�݂̂�\��
			for( int j = 0; j < m_sm_calcs[m_iShowClusterIndx]->GetNumVertices(); j++ )
			{
				int jpIndx = m_sm_calcs[m_iShowClusterIndx]->GetParticleIndx(j);
				tempColor[jpIndx] = 700.0f;
			}
		}
*/

		//���q�x�[�X��
		if( m_iShowClusterIndx == m_iClusteresNum )
		{	//�N���X�^�S�̂�\��
			for( unsigned i = 0; i < m_sm_cluster.size(); i++ )
			{
				for( int j = 0; j < m_sm_cluster[i]->GetNumVertices(); j++ )
				{
					int jpIndx = m_sm_cluster[i]->GetParticleIndx(j);
					tempColor[jpIndx] = 700.0f;
				}
			}
		}
		else
		{	//�P�̃N���X�^�݂̂�\��
			for( int j = 0; j < m_sm_cluster[m_iShowClusterIndx]->GetNumVertices(); j++ )
			{
				int jpIndx = m_sm_cluster[m_iShowClusterIndx]->GetParticleIndx(j);

				//�N���X�^�Ƒ΂ɂȂ闱�q�̐F�͕ς���
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

		m_pPS->SetColorVBOFromArray( tempColor, 1, false, 1.5f * m_ht->getTempMax() );
		delete[] tempColor;
	}		
	//�������p�p�X
	else if( m_iColorType == rxParticleSystemBase::RX_ICE_FAST_PATH )
	{
	}
}


/*!
 * ������Ԃ�tetgen�p�ɕϊ������t�@�C�����쐬�D���_���̂݁D
 * �t�@�C���� src/fltk_sph_turb/bin�@�ɍ쐬�����D
 */
void rxFlGLWindow::Save_NODE_File()
{
	cout << "Save_NODE_File" << endl;
	RXREAL *p = m_pPS->GetArrayVBO(rxParticleSystemBase::RX_POSITION);
	//�t�@�C�����쐬���C�t�H�[�}�b�g�ɂ��������Ĉʒu�����������ށ@#include <fstream>�̈ʒu�ɒ���
	std::ofstream ofs( "sph_pos_before.node" );

	//���_�S�̂̏��@���_���C�������i�R�ŌŒ�j�Cattribute�Cboundarymark�D
	ofs << "       " << m_pPS->GetNumParticles() << " ";
	ofs << "3 ";
	ofs << "0 ";
	ofs << "0 ";
	ofs << endl;

	//���_���ꂼ��̏��
	for(int i = 0; i < m_pPS->GetNumParticles(); i++ )
	{
		ofs << "       " << i+1 << " " << p[DIM*i+0] << " " << p[DIM*i+1] << " " << p[DIM*i+2] << endl;
	}
}

/*!
 * ������Ԃ�tetgen�p�ɕϊ������t�@�C�����쐬�D���_���C�ʏ��D
 * �t�@�C���� src/fltk_sph_turb/bin�@�ɍ쐬�����D
 */
void rxFlGLWindow::Save_POLY_File()
{
	cout << "Save_POLY_File" << endl;
	RXREAL *p = m_pPS->GetArrayVBO(rxParticleSystemBase::RX_POSITION);

	//�t�@�C�����쐬���C�t�H�[�}�b�g�ɂ��������Ĉʒu���Ɩʏ����������ށ@#include <fstream>�̈ʒu�ɒ���
	ofstream ofs( "sph_pos_5_5.poly" );

	//Poly�t�@�C���p�e�X�g
//�P�@���_�̈ʒu���
	//���_�S�̂̏��@���_���C�������i�R�ŌŒ�j�Cattribute�Cboundarymark�D
	ofs << "# Part 1 - node list" << endl;
	ofs << "       " << m_pPS->GetNumParticles() << " ";
	ofs << "3 ";
	ofs << "0 ";
	ofs << "0 ";
	ofs << endl;

	//���_���ꂼ��̏��
	for(int i = 0; i < m_pPS->GetNumParticles(); i++ )
	{
		ofs << "       " << i+1 << " " << p[DIM*i+0] << " " << p[DIM*i+1] << " " << p[DIM*i+2] << endl;
	}
//�Q�@�ʂ̐ڑ����
	//���_�ō쐬�����ʂ̏��
	ofs << "# Part 2 - facet list" << endl;

	int n = pow( m_pPS->GetNumParticles(), 1.0/3.0 ) + 0.5;	//�����̂̂P�ӂ̒��_��
	if( n == 1 )
	{
		cout << "error::n == 1" << endl;
	}

	////�\�ʗ��q�̂݃o�[�W����
//	ofs << "\t" << (n-1) * (n-1) * 6 << "\t";		//������Ԃ͕K�������́C�Ƃ����ꍇ�̖ʂ̐�
	ofs << "\t" << (n-1) * (n-1) * n * 3 << "\t";		//������Ԃ͕K�������́C�Ƃ����ꍇ�̖ʂ̐�
	ofs << "\t" << 1 << "\t";
	ofs << endl;

	//�㉺��
	for( int k = 0; k < n; k++ )
	{
		for( int i = 1; i < n; i++ )
		{
			for( int j = 1; j < n; j++ )
			{
				ofs << "\t" << 1 << "\t" << 0 << "\t" << 1 << endl;		//�p�����[�^
				ofs << "\t" << 4					<< "\t"				//�ʂ̐�
							<< k*n*n+(i-1)*n+j		<< "\t" 
							<< k*n*n+(i-1)*n+j+1	<< "\t"
							<< k*n*n+(i-1)*n+j+1+n	<< "\t"
							<< k*n*n+(i-1)*n+j+n
				<< endl;
			}
		}
	}

	//���E��
	for( int k = 0; k < n; k++ )
	{
		for( int i = 1; i < n; i++ )
		{
			for( int j = 1; j < n; j++ )
			{
				ofs << "\t" << 1 << "\t" << 0 << "\t" << 1 << endl;		//�p�����[�^
				ofs << "\t" << 4						<< "\t"			//�ʂ̐�
							<< k+(i-1)*n*n+(j-1)*n+1	<< "\t" 
							<< k+(i-1)*n*n+(j-1)*n+1+n	<< "\t"
							<< k+i*n*n+(j-1)*n+1+n		<< "\t"
							<< k+i*n*n+(j-1)*n+1
				<< endl;
			}
		}
	}

	//�O���
	for( int k = 0; k < n; k++ )
	{
		for( int i = 1; i < n; i++ )
		{
			for( int j = 1; j < n; j++ )
			{
				ofs << "\t" << 1 << "\t" << 0 << "\t" << 1 << endl;		//�p�����[�^
				ofs << "\t" << 4						<< "\t"			//�ʂ̐�
							<< k*n+(i-1)*n*n+j			<< "\t" 
							<< k*n+(i-1)*n*n+j+1		<< "\t"
							<< k*n+(i-1)*n*n+j+1+n*n	<< "\t"
							<< k*n+(i-1)*n*n+j+n*n
				<< endl;
			}
		}
	}

//�R�@
	ofs << "# Part 3 - hole list" << endl;
	ofs << "0";
	ofs << endl;

//�S
	ofs << "# Part 4 - region list" << endl;
	ofs << "0";
	ofs << endl;
}

/*!
 * tetgen�œ���ꂽ�t�@�C����ǂݍ��݁C������Ԃ��쐬�D�l�ʑ̏��D
 * �t�@�C���� src/fltk_sph_turb/bin�@����ǂݍ��ށD
 */
void rxFlGLWindow::Load_ELE_File(string name)
{
	cout << "Load_ELE_File" << endl;
	//�t�@�C����ǂݍ��݁C�l�ʑ̂ƂȂ�_�̑g�ݍ��킹��List�ɓ����D
	ifstream ifs( name );
	string str;

	//�t�@�C���̑��݊m�F
	if(ifs.fail()) 
	{
		cerr << "File do not exist.\n";
		exit(0);
	}

	//�ϐ��̗p�ӁC������
	int a=0, b=0, c=0, d=0, e=0, f=0;
	bool line_1 = false;

	m_vviTetraList.clear();

	//������̓ǂݍ���
	while( getline(ifs, str) )
	{
		a=0; b=0; c=0; d=0; e=0;
		//������肾���ǂƂ肠����
		if( !line_1 )
		{
			line_1 = true;
			sscanf(str.data(), "%d %d %d", &a, &b, &c);
			
			//cout << "a = " << a << "\t";
			//cout << "b = " << b << "\t";
			//cout << "c = " << c << endl;
		}
		else
		{
			if( str[0] == '#' )
			{
				continue;
			}

			sscanf(str.data(), "%d %d %d %d %d", &a, &b, &c, &d, &e);
			//cout << "a = " << a << "\t";
			//cout << "b = " << b << "\t";
			//cout << "c = " << c << "\t";
			//cout << "d = " << d << "\t";
			//cout << "e = " << e << endl;

			vector<int> list;
			list.push_back(b);
			list.push_back(c);
			list.push_back(d);
			list.push_back(e);

			m_vviTetraList.push_back( list );
		}
	}
}

/*!
 * tetgen��p���ē_�Q�̈ʒu���C���b�V����񂩂珉���̎l�ʑ̍\�����쐬.
 * �t�@�C���� src/fltk_sph_turb/bin�@�ɍ����D
 */
void rxFlGLWindow::MakeTetrahedra()
{	cout << __FUNCTION__ << endl;
	tetgenio in, out;	// ���̓��b�V���Əo�͎l�ʑ̃��b�V��
	RXREAL *p = m_pPS->GetArrayVBO(rxParticleSystemBase::RX_POSITION);

	// �|���S�����_�C���f�b�N�X�̃X�^�[�g(0�X�^�[�g��1�X�^�[�g)
	in.firstnumber = 1;

	// ���b�V�����_�̐ݒ�
//	in.numberofpoints = m_pPS->GetNumParticles();
	in.numberofpoints = ICENUM;
	in.pointlist = new double[in.numberofpoints*3];
	for(int i = 0; i < in.numberofpoints; ++i){
		for(int j = 0; j < 3; ++j){
			//���q�ʒu���̓o�^
			in.pointlist[3*i+j] = (double)p[DIM*i+j];
		}
	}

	// �|���S���̐ݒ�
	int n = pow( in.numberofpoints, 1.0/3.0 ) + 0.5;	//�����̂̂P�ӂ̒��_��
	if( n == 1 ){	cout << "error::n == 1" << endl;	}

	//�e�|���S���̒��_�ԍ����v�Z
	//�i�q��ɔz�u���ꂽ���q�ō����ʂ�����Ă���D
	//�i�q�̓����ɑ��݂���ʂ��l������K�v������̂ł�₱�����Ȃ��Ă���D
	vector<vector<int>> poligonList;
	int poligonNum = 0;

	//�㉺��
	for( int k = 0; k < n; k++ ){
		for( int i = 1; i < n; i++ ){
			for( int j = 1; j < n; j++ ){
				vector<int> list;
				list.push_back(k*n*n+(i-1)*n+j);
				list.push_back(k*n*n+(i-1)*n+j+1);
				list.push_back(k*n*n+(i-1)*n+j+1+n);
				list.push_back(k*n*n+(i-1)*n+j+n);
				poligonList.push_back(list);
	}}}

	//���E��
	for( int k = 0; k < n; k++ ){
		for( int i = 1; i < n; i++ ){
			for( int j = 1; j < n; j++ ){
				vector<int> list;
				list.push_back(k+(i-1)*n*n+(j-1)*n+1);
				list.push_back(k+(i-1)*n*n+(j-1)*n+1+n);
				list.push_back(k+i*n*n+(j-1)*n+1+n);
				list.push_back(k+i*n*n+(j-1)*n+1);
				poligonList.push_back(list);
	}}}

	//�O���
	for( int k = 0; k < n; k++ ){
		for( int i = 1; i < n; i++ ){
			for( int j = 1; j < n; j++ ){
				vector<int> list;
				list.push_back(k*n+(i-1)*n*n+j);
				list.push_back(k*n+(i-1)*n*n+j+1);
				list.push_back(k*n+(i-1)*n*n+j+1+n*n);
				list.push_back(k*n+(i-1)*n*n+j+n*n);
				poligonList.push_back(list);
	}}}

	in.numberoffacets = (n-1) * (n-1) * n * 3;				//�ʂ̑��� (��ӂ̒��_��-1)��2��*(��ӂ̒��_��)*�i���̐��j
	in.facetlist = new tetgenio::facet[in.numberoffacets];
	in.facetmarkerlist = new int[in.numberoffacets];

	for(int i = 0; i < in.numberoffacets; ++i){
		tetgenio::facet *f = &in.facetlist[i];

		// �|���S�����X�g
		f->numberofpolygons = 1;
		f->polygonlist = new tetgenio::polygon[f->numberofpolygons];

		// hole(��)���X�g
		f->numberofholes = 0;
		f->holelist = NULL;

		// �|���S�����_�C���f�b�N�X
		tetgenio::polygon *p = &f->polygonlist[0];

		p->numberofvertices = 4;
		p->vertexlist = new int[p->numberofvertices];
		for(int j = 0; j < p->numberofvertices; ++j)
		{
			p->vertexlist[j] = poligonList[i][j];
		}

		in.facetmarkerlist[i] =  1;	//??
	}

	// ���̓��b�V�������t�@�C���Ƀ_���v
	//in.save_poly("test_poly");	// test.poly
	
	// �������́C"p":PLC�ǂݍ��݁C"q":quality mesh generation(q�̌��quality bound�𐔒l�Ŏw��)�C
	// "a":�ő�̐ϐ���(a�̌�ɑ̐ς𐔒l�Ŏw��)
	tetgenbehavior b;
//	if(!b.parse_commandline("-pq1.414a0.1n")){
	if(!b.parse_commandline("-pn")){
		terminatetetgen(0);
	}

	tetrahedralize(&b, &in, &out);				// �l�ʑ̃��b�V������

	//std::cout << " the number of node   : "    << out.numberofpoints << endl;
	//std::cout << " the number of element   : " << out.numberoftetrahedra << endl;
	//std::cout << " the number of face  : "     << out.numberoftrifaces << endl;

	// �o�̓��b�V�������t�@�C���Ƀ_���v
//	out.save_elements("test_out");	// .ele

	// �o�͎l�ʑ̒��_�ԍ���z��Ɋi�[
	int nelem = out.numberoftetrahedra;
	int nc = out.numberofcorners;
	for(int i = 0; i < nelem; ++i)
	{
		vector<int> list;
		for(int j = 0; j < nc; ++j)
		{
			int pIndx = out.tetrahedronlist[nc*i+j]-out.firstnumber;
			list.push_back(pIndx);
		}
		m_vviTetraList.push_back( list );
	}

//	delete[] in.pointlist;	//�K�v�H
}

/*!
 * tetgen��p���ē_�Q�̈ʒu���C���b�V����񂩂珉���̎l�ʑ̍\�����쐬.
 * �t�@�C���� src/fltk_sph_turb/bin�@�ɍ����D
 */
void rxFlGLWindow::MakeTetrahedraOnlySurface()
{	cout << __FUNCTION__ << endl;
	tetgenio in, out;	// ���̓��b�V���Əo�͎l�ʑ̃��b�V��
	RXREAL *p = m_pPS->GetArrayVBO(rxParticleSystemBase::RX_POSITION);

	// �|���S�����_�C���f�b�N�X�̃X�^�[�g(0�X�^�[�g��1�X�^�[�g)
	in.firstnumber = 1;

	// ���b�V�����_�̐ݒ�
//	in.numberofpoints = m_pPS->GetNumParticles();
	in.numberofpoints = ICENUM;
	in.pointlist = new double[in.numberofpoints*3];
	for(int i = 0; i < in.numberofpoints; ++i){
		for(int j = 0; j < 3; ++j){
			in.pointlist[3*i+j] = (double)p[DIM*i+j];
		}
	}

	// �|���S���̐ݒ�
	int n = pow( in.numberofpoints, 1.0/3.0 ) + 0.5;	//�����̂̂P�ӂ̒��_��
	if( n == 1 ){	cout << "error::n == 1" << endl;	}

	//�e�|���S���̒��_�ԍ����v�Z
	vector<vector<int>> poligonList;
	int poligonNum = 0;

	//�㉺��
	for( int k = 0; k < n; k++ ){
		for( int i = 1; i < n; i++ ){
			for( int j = 1; j < n; j++ ){
				vector<int> list;
				list.push_back(k*n*n+(i-1)*n+j);
				list.push_back(k*n*n+(i-1)*n+j+1);
				list.push_back(k*n*n+(i-1)*n+j+1+n);
				list.push_back(k*n*n+(i-1)*n+j+n);
				poligonList.push_back(list);
	}}}

	//���E��
	for( int k = 0; k < n; k++ ){
		for( int i = 1; i < n; i++ ){
			for( int j = 1; j < n; j++ ){
				vector<int> list;
				list.push_back(k+(i-1)*n*n+(j-1)*n+1);
				list.push_back(k+(i-1)*n*n+(j-1)*n+1+n);
				list.push_back(k+i*n*n+(j-1)*n+1+n);
				list.push_back(k+i*n*n+(j-1)*n+1);
				poligonList.push_back(list);
	}}}

	//�O���
	for( int k = 0; k < n; k++ ){
		for( int i = 1; i < n; i++ ){
			for( int j = 1; j < n; j++ ){
				vector<int> list;
				list.push_back(k*n+(i-1)*n*n+j);
				list.push_back(k*n+(i-1)*n*n+j+1);
				list.push_back(k*n+(i-1)*n*n+j+1+n*n);
				list.push_back(k*n+(i-1)*n*n+j+n*n);
				poligonList.push_back(list);
	}}}

	in.numberoffacets = (n-1) * (n-1) * n * 3;
	in.facetlist = new tetgenio::facet[in.numberoffacets];
	in.facetmarkerlist = new int[in.numberoffacets];

	for(int i = 0; i < in.numberoffacets; ++i){
		tetgenio::facet *f = &in.facetlist[i];

		// �|���S�����X�g
		f->numberofpolygons = 1;
		f->polygonlist = new tetgenio::polygon[f->numberofpolygons];

		// hole(��)���X�g
		f->numberofholes = 0;
		f->holelist = NULL;

		// �|���S�����_�C���f�b�N�X
		tetgenio::polygon *p = &f->polygonlist[0];

		p->numberofvertices = 4;
		p->vertexlist = new int[p->numberofvertices];
		for(int j = 0; j < p->numberofvertices; ++j)
		{
			p->vertexlist[j] = poligonList[i][j];
		}

		in.facetmarkerlist[i] =  1;	//??
	}

	// ���̓��b�V�������t�@�C���Ƀ_���v
	//in.save_poly("test_poly");	// test.poly
	
	// �������́C"p":PLC�ǂݍ��݁C"q":quality mesh generation(q�̌��quality bound�𐔒l�Ŏw��)�C
	// "a":�ő�̐ϐ���(a�̌�ɑ̐ς𐔒l�Ŏw��)
	tetgenbehavior b;
//	if(!b.parse_commandline("-pq1.414a0.1n")){
	if(!b.parse_commandline("-pn")){
		terminatetetgen(0);
	}

	tetrahedralize(&b, &in, &out);				// �l�ʑ̃��b�V������

	//std::cout << " the number of node   : "    << out.numberofpoints << endl;
	//std::cout << " the number of element   : " << out.numberoftetrahedra << endl;
	//std::cout << " the number of face  : "     << out.numberoftrifaces << endl;

	// �o�̓��b�V�������t�@�C���Ƀ_���v
//	out.save_elements("test_out");	// .ele

	// �o�͎l�ʑ̒��_�ԍ���z��Ɋi�[
	int nelem = out.numberoftetrahedra;
	int nc = out.numberofcorners;
	for(int i = 0; i < nelem; ++i)
	{
		vector<int> list;
		for(int j = 0; j < nc; ++j)
		{
			int pIndx = out.tetrahedronlist[nc*i+j]-out.firstnumber;
			list.push_back(pIndx);
		}
		m_vviTetraList.push_back( list );
	}
}


/*!
 * ���q��񂩂�l�ʑ̂��쐬
 * @param[in] pList �l�ʑ̂��쐬���闱�q���X�g
 * @param[in] tList �쐬�����l�ʑ̂̔ԍ����X�g
 *
 * tetgen��p���ē_�Q�̈ʒu���C���b�V����񂩂珉���̎l�ʑ̍\�����쐬.
 * �t�@�C���� src/fltk_sph_turb/bin�@�ɍ����D
 */
void rxFlGLWindow::MakeFreezeTetrahedra(const vector<int>& pList, vector<int>& tList)
{
	RXREAL *p = m_pPS->GetArrayVBO(rxParticleSystemBase::RX_POSITION);
	vector<vector<rxNeigh>>& neights = ((RXSPH*)m_pPS)->GetNeights();		//�ߖT���q���擾

	//�Ìŗ��q���ɋߖT�ő̗��q��T���C�l�ʑ̂��쐬
	for(unsigned i = 0; i < pList.size(); i++)
	{	
		int ipIndx = pList[i];
		vector<int> npList;													//�ߖT�ő̗��q

		//cout << __FUNCTION__ << " check1 i = " << i << " ipIndx = " << ipIndx << endl;

		//�Ìŗ��q�̋ߖT�ő̗��q��T��
		for(unsigned j = 0; j < neights[ipIndx].size(); j++)
		{
			int jpIndx = neights[ipIndx][j].Idx;							//�ߖT���q�̔ԍ�
			float distance	= neights[ipIndx][j].Dist2;
			float radius	= pow(((RXSPH*)m_pPS)->GetEffectiveRadius()*0.95f, 2); 

			//cout << "distance = " << distance << " radius = " << radius << endl;

			if(ipIndx == jpIndx )								continue;	//�������g�͏���
			if(m_ht->getTemps()[jpIndx] > 250)					continue;	//���q�̉��x�����ȉ�
			if(m_ice->GetParticleNum() <= jpIndx)			{	continue;	}	//�Z���݂̂̎����̂Ƃ��ɕK�v�ɂȂ�D
			if(m_ice->GetPtoCNum(jpIndx) <= 0)					continue;	//�N���X�^�ɏ������Ă��Ȃ��Ȃ�߂�
			if(m_ice->GetPtoTNum(jpIndx) <= 0)					continue;	//�l�ʑ̂ɑ����Ă��Ȃ��Ȃ�߂�
			if(distance > radius)								continue;	//�Ìŗ��q�Ƃ̋������e�����a�̔���

			npList.push_back(jpIndx);										//���q�ԍ���o�^
		}

		//�R�@�ߖT�X�򗱎q���ŕ���@�I�����ŒT��
		int iParticleSize = npList.size();
		if(iParticleSize == 0)
		{
			//cout << __FUNCTION__ << " check2a iParticleSize = " << iParticleSize << endl;
			continue;
		}
		else if(1 <= iParticleSize && iParticleSize <= 2)
		{
			//�ŋߖT�X�򗱎q���R�ɂȂ�܂ŒT���@�ߖT�N���X�^0�w�ڂ܂ŒT��
			//cout << __FUNCTION__ << " check2b iParticleSize = " << iParticleSize << endl;

			for(int j = 0; j < iParticleSize; j++)
			{
				int jpIndx = npList[j];
				vector<int> jpList = GetNeighborDistanceIceParticles(jpIndx, 0, 20);			//�ő�20�擾

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

		//�l�ʑ̍쐬�@���]�ڗ��q�{�X�򗱎q
		npList.push_back(ipIndx);
		iParticleSize = npList.size();

		if(iParticleSize >= 4)
		{
			//cout << __FUNCTION__ << " check2a iParticleSize = " << iParticleSize << endl;
			//����@�����̎l�ʑ̂��l���Ȃ�
			while(npList.size() > 4)
			{	npList.erase(npList.begin());	}

			tList.push_back(m_iTetraNum);

			MakeTetraInfo(m_iTetraNum, npList);								//�l�ʑ̍쐬�C�e����o�^
			
			//�f�o�b�O
			//cout << __FUNCTION__ << " Debug ipIndx = " << ipIndx << " npList.size() = " << npList.size() << endl;
			//for( unsigned j = 0; j < npList.size(); j++ )
			//{
			//	cout << " " << npList[j];
			//}
			//cout << endl;

			m_fIntrps[ipIndx] = 1.0f;										//���`��Ԃ����Ȃ�
			m_ht->setPhaseChange(ipIndx, 0);								//���]�ڂ��I��������Ƃ�`����
		}
		else
		{
			//cout << __FUNCTION__ << " check2b iParticleSize = " << iParticleSize << endl;
		}
	}//for(unsigned i = 0; i < pList.size(); i++)
}

/*!
 * �Ìŗ��q�݂̂Ŏl�ʑ̂��쐬
 * @param[in] pList �l�ʑ̂��쐬���闱�q���X�g
 * @param[in] tList �ǉ����ꂽ�l�ʑ̂̔ԍ����X�g
 */
void rxFlGLWindow::MakeFreezeTetrahedra_OnlyFreezeParticle(const vector<int>& pList, vector<int>& tList)
{

}

/*!
 * �Ìŗ��q���l�ʑ̂ɒǉ�
 * @param[in] pList �l�ʑ̂��쐬���闱�q���X�g
 * @param[in] tList �ǉ����ꂽ�l�ʑ̂̔ԍ����X�g
 */
 void rxFlGLWindow::AddFreezeTetrahedra(const vector<int>& pList, vector<int>& tList)
{

}


/*!
 * FPS���v�Z���ăE�B���h�E�^�C�g���ɕ\��
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
		
		g_fTotalTime += g_TimerFPS.GetTime(0);
		g_iTimeCount++;
		g_fAvgTime = g_fTotalTime/(float)g_iTimeCount;

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
 * "\n"���܂܂��string�𕡐���string�ɕ�������
 * @param[in] org ���̕�����
 * @param[in] div ������̕�����z��
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
 * ��ʏo�͗p�̕�����̐���
 * @return ������
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
	// ����
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

	// ���b�V��
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

	// �v�����ꂽ�v�Z����
	str.push_back("time : ");
	string tstr;
	RXTIMER_STRING(tstr);
	DivideString(tstr, str);

	return str;
}


/*!
 * �p�[�e�B�N�����x��GL_LINES�ŕ`��
 * @param[in] prts �p�[�e�B�N���ʒu
 * @param[in] vels �p�[�e�B�N�����x
 * @param[in] n �p�[�e�B�N����
 * @param[in] d �z��̃X�e�b�v
 * @param[in] len ���̒���
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


//�ǉ�
/*!
 * �p�[�e�B�N�����x��GL_LINES�ŕ`��
 * @param[in] prts �p�[�e�B�N���ʒu
 * @param[in] vels �p�[�e�B�N�����x
 * @param[in] n �p�[�e�B�N����
 * @param[in] d �z��̃X�e�b�v
 * @param[in] len ���̒���
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
 * �p�[�e�B�N����GL_POINTS�ŕ`��
 *  - VBO������Ηp����
 * @param[in] vbo �p�[�e�B�N���ʒu���i�[����VBO
 * @param[in] n �p�[�e�B�N����
 * @param[in] color_vbo �p�[�e�B�N���̐F���i�[����VBO
 * @param[in] data �p�[�e�B�N�����W(�z�X�g���������Cvbo��0�̎��Ɏg�p)
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
 * �T�u�p�[�e�B�N���̕`��
 * @param[in] pVBO �ʒu��VBO
 * @param[in] uColorVBO �F��VBO
 * @param[in] data �ʒu(CPU)
 */
void rxFlGLWindow::DrawSubParticles(void)
{
	// Et�Ɋ�Â��Ċe�T�u�p�[�e�B�N���̍����䗦�Ɣ��a���v�Z���ċl�߂�
	((RXSPH*)m_pPS)->CalRadiusAndRatio();

	// �T�u�p�[�e�B�N����
	int num_all_sp = ((RXSPH*)m_pPS)->GetNumAllSubParticles();
	int num_sp = ((RXSPH*)m_pPS)->GetNumValidSubParticles();
	int num_p = m_pPS->GetNumParticles();

	if(!num_sp) return;

	// �T�u�p�[�e�B�N���̈ʒu�Ɣ��a(�L���ȕ������l�߂��z��)
	RXREAL *sppos = ((RXSPH*)m_pPS)->GetSubParticlePos();
	RXREAL *sprad = ((RXSPH*)m_pPS)->GetSubParticleRad();

	// �p�[�e�B�N���̑傫��
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
 * Vertex Buffer Object(VBO)�̍쐬
 * @param[inout] vbo �o�b�t�@ID
 * @param[in] size �o�b�t�@�T�C�Y
 */
void rxFlGLWindow::CreateVBO(GLuint* vbo, unsigned int size)
{
	glGenBuffers(1, vbo);
	glBindBuffer(GL_ARRAY_BUFFER, *vbo);

	// �o�b�t�@�I�u�W�F�N�g�̏�����
	glBufferData(GL_ARRAY_BUFFER, size, 0, GL_DYNAMIC_COPY);
	glBindBuffer(GL_ARRAY_BUFFER, 0);

	//CuRegisterGLBufferObject(*vbo);
}

/*!
 * Vertex Buffer Object(VBO)�̍폜
 * @param[inout] vbo �o�b�t�@ID
 */
void rxFlGLWindow::DeleteVBO(GLuint* vbo)
{
	//CuUnregisterGLBufferObject(*vbo);

	glBindBuffer(1, *vbo);
	glDeleteBuffers(1, vbo);

	*vbo = 0;
}


/*!
 * VBO��p�������l�ʃ|���S���`��
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
			// ��
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

//�ǉ��F�X�p
/*
 * VBO��p�������l�ʃ|���S���`��(�ő�)
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
 * �V�[���`��
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
	// ���͋��E
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

	// ��
	if((m_iDraw & RXD_AXIS)){
		Vec3 len = 0.6*dim;
		glLineWidth((GLfloat)3.0);

		// x��
		glColor3f(1.0, 0.0, 0.0);
		glBegin(GL_LINES);
		glVertex3d(0.0, 0.0, 0.0);
		glVertex3d(len[0], 0.0, 0.0);
		glEnd();

		// y��
		glColor3f(0.0, 1.0, 0.0);
		glBegin(GL_LINES);
		glVertex3d(0.0, 0.0, 0.0);
		glVertex3d(0.0, len[1], 0.0);
		glEnd();

		// z��
		glColor3f(0.0, 0.0, 1.0);
		glBegin(GL_LINES);
		glVertex3d(0.0, 0.0, 0.0);
		glVertex3d(0.0, 0.0, len[2]);
		glEnd();
	}


	//
	// ��Q��
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
				// �`��t���O : ���ʃr�b�g���璸�_,�G�b�W,��,�@�� - 1,2,4,8)
				itr->Draw(m_iSolidDraw);
			}
		}
	}
	
	// 
	// �p�[�e�B�N��
	// 
	if(!m_bsSimuSetting.at(ID_SPH_ANISOTROPIC) && (m_iDraw & RXD_PARTICLE) && !(m_iDraw & RXD_REFRAC)){
		glDisable(GL_LIGHTING);
		glColor4d(0.0, 0.0, 1.0, 1.0);
		glEnable(GL_POINT_SMOOTH);
		glPointSize(1.0);

		unsigned int pvbo = m_pPS->GetCurrentReadBuffer();
		unsigned int cvbo = m_pPS->GetColorBuffer();

		if(m_iDrawPS == RXP_POINT){				// GL_POINTS�ŕ`��
			RXREAL *data = 0;
			data = m_pPS->GetArrayVBO(rxParticleSystemBase::RX_POSITION);
			DrawParticlePoints(0, pnum, cvbo, data);
		}
		else if(m_iDrawPS == RXP_POINTSPRITE){	// �|�C���g�X�v���C�g�ŕ`��
			float prad = (float)m_pPS->GetParticleRadius();			// �p�[�e�B�N�����a
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
	// ���E�p�[�e�B�N��
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
	// ���b�V���`��
	//
	if((m_iDraw & RXD_MESH)){
		if((m_iDraw & RXD_REFRAC) && m_bUseCubeMap){
			glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
			glEnable(GL_LIGHTING);

			glEnable(GL_COLOR_MATERIAL);

			Vec3 eye_pos(0.0);		// ���_�ʒu
			m_tbView.CalLocalPos(eye_pos.data, Vec3(0.0).data);
			
			glUseProgram(g_glslFresnel.Prog);

			// �p�����[�^�ݒ�
			// �o�[�e�b�N�X�V�F�[�_�p�p�����[�^
			glUniform1f(glGetUniformLocation(g_glslFresnel.Prog, "etaRatio"), 0.93);
			glUniform1f(glGetUniformLocation(g_glslFresnel.Prog, "fresnelBias"), 0.005);
			glUniform1f(glGetUniformLocation(g_glslFresnel.Prog, "fresnelPower"), 0.98);
			glUniform1f(glGetUniformLocation(g_glslFresnel.Prog, "fresnelScale"), 1.0);
			glUniform3f(glGetUniformLocation(g_glslFresnel.Prog, "eyePosition"), eye_pos[0], eye_pos[1], eye_pos[2]);

			// �t���O�����g�V�F�[�_�p�p�����[�^
			glUniform1i(glGetUniformLocation(g_glslFresnel.Prog, "envmap"), 0);

			//�X�Ɛ��ŐF��ς���
//			glColor3d(0.3, 0.3, 1.0);
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
			DrawLiquidSurface();
//			glColor4d(0.3, 1.0, 0.3, 1.0);
//			DrawSolidSurface();
		}
	}

	//
	// �p�[�e�B�N�����x
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
	// �p�[�e�B�N���������x
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
	// �����Z��
	//
	if((m_iDraw & RXD_CELLS)){
		glDisable(GL_LIGHTING);
		glLineWidth(1.0);
		m_pPS->DrawCells(Vec3(0.0, 1.0, 0.0), Vec3(1.0, 0.0, 0.0));
	}

	//
	// �f�o�b�O�p�x�N�g����
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
	//�������p�p�X
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

}



//-----------------------------------------------------------------------------
// �e��֐�
//-----------------------------------------------------------------------------
/*!
 * xyz���`��(x��:��,y��:��,z��:��)
 * @param[in] len ���̒���
 */
int DrawAxis(double len, double line_width)
{
	glLineWidth((GLfloat)line_width);

	// x��
	glColor3f(1.0, 0.0, 0.0);
	glBegin(GL_LINES);
	glVertex3d(0.0, 0.0, 0.0);
	glVertex3d(len, 0.0, 0.0);
	glEnd();

	// y��
	glColor3f(0.0, 1.0, 0.0);
	glBegin(GL_LINES);
	glVertex3d(0.0, 0.0, 0.0);
	glVertex3d(0.0, len, 0.0);
	glEnd();

	// z��
	glColor3f(0.0, 0.0, 1.0);
	glBegin(GL_LINES);
	glVertex3d(0.0, 0.0, 0.0);
	glVertex3d(0.0, 0.0, len);
	glEnd();

	return 1;
}

/*!
 * �C�Ӄx�N�g������̉�]
 * @param[in] pos ���̍��W�l
 * @param[in] axis ��]��
 * @param[in] ang ��]�p�x(rad)
 * @return ��]�������W�l
 */
inline Vec3 Rotate(const Vec3 &pos, const Vec3 &axis, const double &ang)
{
	Vec3 pos1;    // ��]��̍��W�l

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
 * ������`��
 * @param[in] cir_str ������z�o�b�t�@
 * @param[in] static_str �ÓI�ȕ�����o�b�t�@
 * @param[in] w,h �E�B���h�E�T�C�Y
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
		// FTGL�ŕ������`��
		for(int j = 0; j < (int)static_str.size(); ++j){
			glRasterPos2f(x0, y0);
			g_pFont->Render(static_str[j].c_str());
			y0 -= g_pFont->LineHeight();
		}
	}
	else{
		// glutBitmapString�ŕ�����`��
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
 * ������`��
 * @param[in] cir_str ������z�o�b�t�@
 * @param[in] static_str �ÓI�ȕ�����o�b�t�@
 * @param[in] w,h �E�B���h�E�T�C�Y
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

		// FTGL�ŕ�����`��
		for(int j = 0; j < (int)static_str.size(); ++j){
			glRasterPos2i(20+offsetx, y0);
			g_pFont->Render(static_str[j].c_str());

			y0 -= g_pFont->LineHeight();
		}
	}
	else{
		float y0 = static_str.size()*20.0f;
		glRasterPos2f(20.0f+offsetx, y0);

		// glutBitmapString�ŕ�����`��
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
// �摜�̓Ǎ�
//-----------------------------------------------------------------------------
static int m_iBitmapType = RX_BMP_WINDOWS_V3;
inline void SetBitmapType(int type){ m_iBitmapType = type; }


unsigned char* ReadImageFile(const std::string &fn, int &w, int &h, int &c)
{
	// �g���q���o
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
		while(itr != ext.begin()){	// �p�X�̍Ō��\0��X�y�[�X���������Ƃ��̑΍�
			if(*itr == 0 || *itr == 32){
				ext.erase(itr--);
			}
			else{
				itr--;
			}
		}
	}

	// �摜�ǂݍ���
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
	// �g���q���o
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
		while(itr != ext.begin()){	// �p�X�̍Ō��\0��X�y�[�X���������Ƃ��̑΍�
			if(*itr == 0 || *itr == 32){
				ext.erase(itr--);
			}
			else{
				itr--;
			}
		}
	}

	// �摜�ǂݍ���
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
 * �t���[���o�b�t�@��RGB�����ꎞ�I�ȃo�b�t�@�Ɋm��
 * @retval true �ۑ�����
 * @retval false �ۑ����s
 */
bool SaveFrameBuffer(const string &fn, int w, int h)
{
	vector<unsigned char> imm_buf(w*h*3);

	glReadPixels(0, 0, w, h, GL_RGB, GL_UNSIGNED_BYTE, &imm_buf[0]);

	// �㉺���]
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


/*! �e�N�X�`�����r�b�g�}�b�v�ŕۑ�
	@param[in] file_name �ۑ��t�@�C����
	@param[in] tex_obj �ۑ��������e�N�X�`���I�u�W�F�N�g
	@param[in] jpeg_quality JPEG�ŕۑ�����ꍇ�̕ۑ��i��(0-255)
	@retval true �ۑ�����
	@retval false �ۑ����s
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
 * ReadTexture �e�N�X�`���̓ǂݍ���
 * @param[in] path �e�N�X�`���摜�̃p�X
 * @param[out] tex_obj �e�N�X�`�����i�[����
 * @return �e�N�X�`�����ǂݍ��߂����ǂ���
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

	// �e�N�X�`���o�^
	tex_obj.Bind();

	// �e�N�X�`���p�����[�^�̐ݒ�
	tex_obj.SetParameter(GL_TEXTURE_WRAP_S, GL_REPEAT);
	tex_obj.SetParameter(GL_TEXTURE_WRAP_T, GL_REPEAT);
	tex_obj.SetParameter(GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	tex_obj.SetParameter(GL_TEXTURE_MIN_FILTER, GL_NEAREST);

	tex_obj.SetTexture();

	return true;
}


/*!
 * OpenCV�ŉ摜�Ǎ��� -> OpenGL�e�N�X�`���o�^
 * @param[in] fn �t�@�C����
 * @param[inout] tex_name �e�N�X�`����(0�Ȃ�V���ɐ���)
 * @param[in] mipmap �~�b�v�}�b�v�g�p�t���O
 * @param[in] compress �e�N�X�`�����k�g�p�t���O
int LoadGLTexture(const string &fn, GLuint &tex_name, bool mipmap, bool compress)
{
	// �摜�ǂݍ���
	int w, h, c;
	unsigned char* pimg;
	pimg = ReadImageFile(fn, w, h, c);
	if(!pimg){
		return 0;
	}

	cout << "image : " << w << " x " << h << " x " << c << endl;
	GLuint iformat, format;

	// �摜�t�H�[�}�b�g
	format = GL_RGBA;
	if(c == 1){
		format = GL_LUMINANCE;
	}
	else if(c == 3){
		format = GL_RGB;
	}
 
	// OpenGL�����̊i�[�t�H�[�}�b�g
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
 
	// �e�N�X�`���쐬
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
 * ����PNG�摜�ǂݍ���
 * @param[in] fn �t�@�C����
 * @param[inout] tex_name �e�N�X�`����(0�Ȃ�V���ɐ���)
 * @param[out] w,h �ǂݍ��܂ꂽ�摜�̃T�C�Y
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

	// �e�N�X�`���쐬
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
 * ���}�b�v�p�̃L���[�u�}�b�v�e�N�X�`���̓ǂݍ���
 * @param[in] fn[6] �e�N�X�`���摜(6��)�̃p�X(x+,x-,y+,y-,z+,z-)(�E,��,��,��,��,�O)
 * @param[out] cube_map rxCubeMapData�^
 * @retval true  �L���[�u�}�b�v�p�摜�̓ǂݍ��ݐ���
 * @retval false �L���[�u�}�b�v�p�摜�̓ǂݍ��ݎ��s
 */
bool LoadCubeMapTexture(const string fn[6], rxCubeMapData &cube_map)
{
	GLuint tex_name;
	glGenTextures(1, &tex_name);
	glBindTexture(GL_TEXTURE_CUBE_MAP, tex_name);

	// �L���[�u�}�b�v�e�N�X�`���p�����[�^�̐ݒ�
	glTexParameterf(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);		// �摜���E�̈����̎w��
	glTexParameterf(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
	glTexParameterf(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);
	glTexParameterf(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MAG_FILTER, GL_NEAREST);	// �摜�t�B���^�̎w��
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
 * ���}�b�v�p�̃L���[�u�}�b�v�e�N�X�`���̓ǂݍ���
 * @param cube_map �L���[�u�}�b�v�f�[�^
 * @param base �L���[�u�}�b�v�p�摜�̃t�@�C�����̃x�[�X����
 * @param ext �L���[�u�}�b�v�p�摜�̃t�@�C���̊g���q
 * @retval true  �L���[�u�}�b�v�p�摜�̓ǂݍ��ݐ���
 * @retval false �L���[�u�}�b�v�p�摜�̓ǂݍ��ݎ��s
 */
bool LoadCubeMap(rxCubeMapData &cube_map, string base, string ext)
{
	// �L���[�u�}�b�v�p�摜�̓ǂݍ���(x+,x-,y+,y-,z+,z-)(�E,��,��,��,��,�O)
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
