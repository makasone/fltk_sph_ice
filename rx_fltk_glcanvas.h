/*!
  @file rx_fltk_glcanvas.h
	
  @brief FLTK�ɂ��OpenGL�E�B���h�E�N���X
 
  @author Makoto Fujisawa 
  @date   2011-09
*/

#ifndef _RX_FLTK_GLCANVAS_H_
#define _RX_FLTK_GLCANVAS_H_

//-----------------------------------------------------------------------------
// �C���N���[�h���C�u����
//-----------------------------------------------------------------------------
#pragma comment(lib, "libtet.lib")

//-----------------------------------------------------------------------------
// �C���N���[�h�t�@�C��
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

//�ǉ��F�F
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
// ��`/�萔
//-----------------------------------------------------------------------------
class rxFlWindow;
class rxParticleSystemBase;
struct rxCell;

class rxMCMeshCPU;
class rxMCMeshGPU;
class rxSSMeshCPU;
class rxSSMeshGPU;

//�@�����̒�`��SPH�@��GPU�CCPU�Ȃǂ�؂�ւ���
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

//�N���X�^���쐬����ۂ̗��q���C�l�ʑ̐����w��
//�����g���ĂȂ�
//#define ICENUM 27
//#define ICENUM 125
//#define ICENUM	729
//#define ICENUM	1331	//11_11_11 or �o�j�[���f��
//#define TETRANUM 10900	//�o�j�[���f���̏ꍇ�̎l�ʑ̐�
//#define ICENUM	2197	//13_13_13 or �o�j�[���f��
//#define SIDE	13
//#define ICENUM	3463	//�o�j�[���f��
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
//(5000, 5000, 26000);			//���q��4913�̏ꍇ�̃p�����[�^
//(7000, 7000, 37000);
//(13000, 13000, 50000);
//(10000, 10000, 1);
//(16000, 16000, 1);

//(16000, 16000, 90000);
//(20000, 20000, 110000);
//(24400, 24400, 137000);		//�����Ȃ������D

//�\�ʂ݂̂̏ꍇ
//#define ICENUM	6			//1_1_1
//#define ICENUM	27
//#define ICENUM	54			//3_3_3
//#define ICENUM	1014
//#define ICENUM	2646	//21_21_21
//#define ICENUM	5046	//29_29_29

#define HIGHNUM 8

// �`��t���O
enum
{
	RXD_PARTICLE		= 0x0001,	//!< �p�[�e�B�N��
	RXD_VELOCITY		= 0x0002,	//!< ���x��
	RXD_NORMAL			= 0x0004,	//!< �@��
	RXD_BBOX			= 0x0008,	//!< AABB(�V�~�����[�V�������)
	RXD_CELLS			= 0x0010,	//!< �ߖT�T���p�����Z��
	RXD_MESH			= 0x0020,	//!< ���b�V��
	RXD_SOLID			= 0x0040,	//!< �ő�
	RXD_REFRAC			= 0x0080,	//!< ���ܕ`��

	RXD_TURB_VELOC		= 0x0100,	//!< �������x��
	RXD_FOAM			= 0x0200,	//!< �A

	RXD_PARAMS			= 0x0400,	//!< �p�����[�^��ʕ`��
	RXD_ANISOTROPICS	= 0x0800,	//!< �ٕ����J�[�l��
	RXD_UPDATED_PRTS	= 0x1000,	//!< 
	RXD_AXIS			= 0x2000,   //!< ��

	RXD_BPARTICLE		= 0x4000,	//!< ���E�p�[�e�B�N��
	RXD_DEBUG_VECTOR	= 0x8000,	//!< �f�o�b�O�p�̃x�N�g����
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

// �p�[�e�B�N���`����@
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


// �O�p�`���b�V�������@
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

// �ő̕`��
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


//! �`��̈�T�C�Y���
const string RX_CANVAS_SIZE_STR[] = {
	"1920x1080",	"",
	"1280x720",		"", 
	"1024x768",		"", 
	"800x800",		"", 
	"800x600",		"", 
	"640x480",		"", 
	"-1", 
};


//! SPH�ݒ�
enum SettingMode
{
	ID_SPH_MESH = 0,	// ���b�V������
	ID_SPH_INPUT,		// �p�[�e�B�N���f�[�^�o��
	ID_SPH_OUTPUT,		// �p�[�e�B�N���f�[�^����
	ID_SPH_MESH_OUTPUT, // ���b�V���o��
	ID_SPH_INLET,		// �������E
	ID_SPH_VC, 
	ID_SPH_PS_TURB, 
	ID_SPH_SPS_TURB, 

	ID_HEAT,			//�ǉ��@�M����
	ID_SM,				//�ǉ��@SM�@
	ID_ICE,				//�ǉ��@�X�\��

	ID_SPH_ANISOTROPIC, // �ٕ����J�[�l��

	ID_SPH_END, 
};

//-----------------------------------------------------------------------------
//! rxFlGLWindow�N���X - fltk�ɂ��OpenGL�E�B���h�E
//-----------------------------------------------------------------------------
class rxFlGLWindow : public Fl_Gl_Window
{
protected:
	int m_iWinW;					//!< �`��E�B���h�E�̕�
	int m_iWinH;					//!< �`��E�B���h�E�̍���
	int m_iMouseButton;				//!< �}�E�X�{�^���̏��
	int m_iKeyMod;					//!< �C���L�[�̏��
	rxTrackball m_tbView;			//!< �g���b�N�{�[��

	double m_fBGColor[3];			//!< �w�i�F
	bool m_bAnimation;				//!< �A�j���[�V����ON/OFF
	bool m_bFullscreen;				//!< �t���X�N���[��ON/OFF

	rxFlWindow *m_pParent;			//!< �e�N���X

	vector<rxPolygons> m_vPolys;	//!< �|���S���I�u�W�F�N�g

	// FTGL
	unsigned long m_ulFontSize;		//!< �t�H���g�T�C�Y

	//
	// ���q�@�֘A�ϐ�
	//
	rxParticleSystemBase *m_pPS;	//!< SPH
	double m_fDt;					//!< �^�C���X�e�b�v��
	double m_fGravity;				//!< �d�͉����x

	int m_iCurrentStep;				//!< ���݂̃X�e�b�v��
	bool m_bPause;					//!< �V�~�����[�V�����̃|�[�Y�t���O

public:
	//
	// �ǉ��@�M�����֘A�ϐ�
	//
	vector<float> m_fIntrps;		//SPH�@��SM�@�̐��`��Ԃ̃p�����[�^�z��@�e�p�[�e�B�N�����ƂɓK�p���Ă��

	Vec2 m_ht_vStartPoint;			//��`���̗��q���x���グ�邽�߂̎n�_
	Vec2 m_ht_vEndPoint;			//�I�_

	bool m_ht_bRectFlag;			//��`�����x�ω��@�\�𗘗p���邩�̃t���O

	vector<int> m_ht_vSelectedVertices;

	float m_fMakeParticleTemp;

	//
	//�ǉ��@���ω��I�u�W�F�N�g
	//
	IceObject* m_iceObj;

	//
	// �ǉ��@SM�@�֘A�ϐ�
	//
	vector<Ice_SM*> m_sm_connects;	//�֌W���N���X�^
	vector<Ice_SM*> m_sm_calcs;		//�v�Z�����N���X�^

	int m_iIcePrtNum;				//�����̌ő̗��q��
	int m_iIceTtrNum;				//�����̌ő̗��q�̎l�ʑ̐�
	int m_iClusteresNum;			//�N���X�^�̐��@�g�p���̃N���X�^�ōő�̒l
	int m_iTetraNum;				//���g�p�@�l�ʑ̂̐�
	int m_iTetraNumNum;				//�f�o�b�O�p

	int m_iMaxClusterNum;			//�ő�N���X�^���C�ő嗱�q��
	int m_iIceItr;					//�ő̉^���v�Z�̔�����

	//
	// �ǉ��@�X�\���֘A�ϐ�
	//
	IceStructure *m_ice;

	int meltPIndx;
	int debugIndx;

	float *m_fIceFlag;						//�X�̃t���O�@�f�o�C�X������

	bool *testFlag;							//���s���x���m���߂邽�߂̎���

	vector<vector<int>> m_vviTetraList;		//�l�ʑ̂̑g�ݍ��킹���X�g

	int m_iLayer;							//�v�Z�����N���X�^�쐬�̂��߂Ɏ擾����ߖT�N���X�^�́C���C���[��
	int m_iShowClusterIndx;					//GUI�ŕ\������N���X�^�ԍ��@�O�`�N���X�^���@�l��ctr+shift+q�ŕω�������
	int m_iShowHighClstrIndx;
	int m_iShowTetraIndx;					//GUI�ŕ\������l�ʑ̔ԍ��@�O�`�l�ʑ̐��@�l��ctr+shift+a�ŕω�������
	
	//obj�t�@�C��
	rxPolygons m_poly;

	//
	//�ǉ��@�����؂�ւ�
	//
	int m_bMode;							//���݂̏������[�h
	enum
	{
		MODE_SPH = 0,
		MODE_HEAT,
		MODE_SM,
		MODE_ICE,
		MODE_NUM,
	};

	//�ǉ��@�����_�����O�p�����[�^
	double m_etaRatio;
	double m_fresnelBias;
	double m_fresnelPower;
	double m_fresnelScale;

	//�����J�n�t���O
	bool m_bFall;

	//�e�X�g�N���X
	//mk_CGAL test;

protected:
	// �V�[��
	//string m_strCurrentScene;		//!< ���݂̃V�[���̖��O
	//vector<string> m_vSceneFiles;	//!< �V�[���t�@�C�����X�g
	//int m_iSceneFileNum;			//!< �V�[���t�@�C���̐�
	rxSPHConfig m_Scene;
	int m_iSimuSetting;				//!< �~�����[�V�����ݒ�ۑ��p

	// �p�[�e�B�N�����o��
	string m_strSphOutputName0 ;
	string m_strSphOutputHeader;

	// �p�[�e�B�N��������
	string m_strSphInputName0 ;
	string m_strSphInputHeader;

	// �ő̈ړ��t���O
	bool m_bSolidMove;

	// �ő̂̓���
	Vec3 m_vMovePos[2];
	double m_fMoveMaxVel;
	bool m_bMoveSolid;
	int m_iMoveStart;


	//
	// ���b�V��
	//
	uint m_iNumVrts, m_iNumTris;	//!< �������ꂽ���b�V���̒��_���ƃ��b�V����
	int m_iVertexStore;				//!< �T���v�����O�{�����[�����ɑ΂���\�z����钸�_��(nx*ny*store)

	int m_iMeshMaxN;				//!< ���b�V�����O���b�h��(���E�������Ƃ������������̕�����)
	int m_iMeshN[3];				//!< ���b�V�����O���b�h��

	Vec3 m_vMeshBoundaryExt;		//!< ���b�V�����E�{�b�N�X�̊e�ӂ̒�����1/2
	Vec3 m_vMeshBoundaryCen;		//!< ���b�V�����E�{�b�N�X�̒��S���W

	rxPolygons m_Poly;				//!< ���b�V��
	//vector<rxPolygons*> m_vSolidPoly;//!< �ő̃��b�V��
	rxMaterialOBJ m_matPoly;

	GLuint m_uVrtVBO;				//!< ���b�V�����_(VBO)
	GLuint m_uTriVBO;				//!< ���b�V���|���S��(VBO)
	GLuint m_uNrmVBO;				//!< ���b�V�����_�@��(VBO)

	//�ǉ��F�ő̗p
	uint m_iNumVrts_solid, m_iNumTris_solid;	//!< �������ꂽ���b�V���̒��_���ƃ��b�V����
	GLuint m_uVrtVBO_solid;			//!< ���b�V�����_(VBO)
	GLuint m_uTriVBO_solid;			//!< ���b�V���|���S��(VBO)
	GLuint m_uNrmVBO_solid;			//!< ���b�V�����_�@��(VBO)

	int m_iDimVBO;

	// ���b�V���o��
	int m_iSaveMeshSpacing;

	// �w�i�摜
	bool m_bUseCubeMap;				//!< �L���[�u�}�b�v�g�p�t���O
	rxCubeMapData m_CubeMap;		//!< �L���[�u�}�b�v

	// ���b�V������
	rxMCMeshCPU *m_pMCMeshCPU;
	rxMCMeshGPU *m_pMCMeshGPU;
	rxSSMeshCPU *m_pSSMCPU;			//!< Screen Space Mesh
	rxSSMeshGPU *m_pSSMGPU;			//!< Screen Space Mesh by CUDA
	int m_iDepthFiltering;			//!< ������(�f�v�X�}�b�v0x01, �֊s0x02)
	int m_iSSDebugOutput;

	double m_fSpacing;				//!< �f�v�X�}�b�v�̃T���v�����O�Ԋu
	double m_fPrtRad;				//!< �p�[�e�B�N���̔��a
	double m_fZmax;					//!< �֊s�ƂȂ�f�v�X����臒l
	int m_iNfilter;					//!< �f�v�X�l�������̃t�B���^�T�C�Y
	int m_iNiters;					//!< �֊s�������̔�����

	double m_fProjectionMatrix[16];	//!< �������e�s��
	double m_fModelviewMatrix[16];	//!< ���f���r���[�s��

	int m_iPickedParticle;			//!< �}�E�X�s�b�N���ꂽ�p�[�e�B�N��


public:
	// �`��t���O
	int m_iDraw;					//!< �`��t���O
	int m_iDrawPS;					//!< �p�[�e�B�N���`����@
	int m_iColorType;				//!< �p�[�e�B�N���`�掞�̐F
	int m_iTriangulationMethod;		//!< �O�p�`���b�V�������@
	int m_iSolidDraw;				//!< �ő̕`��

	// �V�~�����[�V�����ݒ�
	bitset<32> m_bsSimuSetting;		//!< �~�����[�V�����ݒ�t���O
	double m_fVScale;				//!< �x�N�g����`�掞�̃X�P�[��
	double m_fMeshThr;				//!< �A�֐����b�V��������臒l

	// �V�[�����X�g
	vector<string> m_vSceneTitles;	//!< �V�[���t�@�C�����X�g
	int m_iCurrentSceneIdx;			//!< ���݂̃V�[���t�@�C��

	// �摜�o��
	int m_iSaveImageSpacing;		//!< �摜�ۑ��Ԋu(=-1�Ȃ�ۑ����Ȃ�)

public:
	//! �R���X�g���N�^
	rxFlGLWindow(int x, int y, int w, int h, const char* l, void *parent);

	//! �f�X�g���N�^
	~rxFlGLWindow();

public:
	// OpenGL������
	void InitGL(void);

	// OpenGL�`��
	void Projection(void);
	vector<string> SetDrawString(void);
	void ReDisplay(void);

	// GUI�R�[���o�b�N
	void Display(void);
	void Resize(int w, int h);
	void Mouse(int button, int state, int x, int y);
	void Motion(int x, int y);
	void PassiveMotion(int x, int y);
	void Idle(void);
	void Timer(void);
	void Keyboard(int key, int x, int y);
	void SpecialKey(int key, int x, int y);

	// �}�E�X�s�b�N�p
	static void Projection_s(void* x);
	static void DisplayForPick_s(void* x);
	void DisplayForPick(void);
	void PickParticle(int x, int y);

	// ���_
	void InitView(void);

	// �A�j���[�V����
	static void OnTimer_s(void* x);
	static void OnIdle_s(void* x);
	
	// �A�j���[�V�����؂�ւ�
	bool SwitchIdle(int on);

	// �t���X�N���[���؂�ւ�
	void SwitchFullScreen(int win = 1);
	int  IsFullScreen(void);


	// �t�@�C�����o��
	void OpenFile(const string &fn);
	void SaveFile(const string &fn);
	void SaveDisplay(const string &fn);
	void SaveDisplay(const int &stp);
	
	void SaveMesh(const string fn, rxPolygons &polys);
	void SaveMesh(const int &stp, rxPolygons &polys);

	// FTGL�t�H���g�ݒ�
	int SetupFonts(const char* file);


public:
	// ���b�V������
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

	//�ǉ��F�X�p
	void DrawSolidSurface(void);

	//�ǉ��F�������p�X
	void DrawFastPath(RXREAL *prts, RXREAL *vels, int n, int d, double *c0, double *c1, double len = 0.1);

	void SetParticleColorType(int type, int change = 0);

	void RenderSphScene(void);

	//�ǉ��F�F�X
	void InitIceObj(void);
	void StepTimeEvent(void);

	void CountTetraHedra(int tIndx, vector<int>& pList);
	void MakeTetraInfo(int tIndx, int* PtoTNum);
	void MakeTetraInfo(int tIndx, vector<int> pList);

	//�ǉ��F�F���q�x�[�X�F�F�N���X�^
	void MakeCluster(int pIndx);
	void MakeClusterFromNeight();

	void StepIceObj();
	void StepIceStructure();
	void StepHeatTransfer();

	void StepClusterCPU(double dt);
	void StepClusterGPU(double dt);

	//�ǉ��F�F���q�x�[�X�F�F�ő̏��i�N���X�^���j
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

	//�ǉ��F�F�l�ʑ̃x�[�X�F�FSM�@
	void InitSM(rxSPHConfig &sph_scene);
	void InitSM_Layer(rxSPHConfig &sph_scene);
	void StepSM(double dt);

	void AddVertexToCluster(int pIndx, int cIndx);

	//�ǉ��F�F�l�ʑ̃x�[�X�F�F�X�\��
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

	//�ǉ��F�F���̑�����
	void Collision(Vec3 &p, Vec3 &np, Vec3 &v, int obj);	//�Փ˔���N���X
	void ClearPick(void);
	void ClearRect(void);
	void StepParticleColor(void);
	void UpdateInfo();

	//�ǉ��F�F�l�ʑ̍쐬�̂��߂̏���
	void MakeFreezeTetrahedra(vector<int>& pList, vector<int>& tList);
	void MakeFreezeTetrahedra_OnlyFreezeParticle(const vector<int>& pList, vector<int>& tList);
	void AddFreezeTetrahedra(const vector<int>& pList, vector<int>& tList);

	void DumpParticleData();

private:
	//! �`��R�[���o�b�N�֐��̃I�[�o�[���C�h
	void draw(void)
	{
		if(!context_valid()) InitGL();
		if(!valid()) Resize(w(), h());
		Display();    // OpenGL�`��
	}

	//! ���T�C�Y�R�[���o�b�N�֐��̃I�[�o�[���C�h
	void resize(int x_, int y_, int w_, int h_)
	{
		Fl_Gl_Window::resize(x_, y_, w_, h_);
		//Resize(w_, h_);
	}


public:
	// �C�x���g�n���h��
	int handle(int e);	// handle�֐��̃I�[�o�[���C�h

	// ���j���[�C�x���g�n���h��
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
