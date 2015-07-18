/*!
  @file rx_fltk_window.h
	
  @brief FLTK�ɂ��E�B���h�E�N���X
 
  @author Makoto Fujisawa 
  @date   2011-08
*/

#ifndef _RX_FLTK_WINDOW_H_
#define _RX_FLTK_WINDOW_H_


//-----------------------------------------------------------------------------
// �C���N���[�h�t�@�C��
//-----------------------------------------------------------------------------
#include <FL/Fl.H>
#include <FL/Fl_Double_Window.H>
#include <FL/Fl_Menu_Bar.H>
#include <FL/Fl_Button.H>
#include <FL/Fl_Scroll.H>
#include <FL/Fl_Box.H>
#include <FL/filename.H>		// �t�@�C����
#include <FL/Fl_File_Chooser.H>
#include <FL/Fl_File_Icon.H>


#include "rx_fltk_glcanvas.h"

using namespace std;

//-----------------------------------------------------------------------------
// MARK:��`/�萔
//-----------------------------------------------------------------------------
const string RX_PROGRAM_NAME = "rx_fltk_sph";
extern vector<string> g_vDefaultFiles;



//-----------------------------------------------------------------------------
//! rxFlWindow�N���X - fltk�ɂ��E�B���h�E
//-----------------------------------------------------------------------------
class rxFlWindow : public Fl_Double_Window
{
protected:
	// �����o�ϐ�
	int m_iWinW;		//!< �`��E�B���h�E�̕�
	int m_iWinH;		//!< �`��E�B���h�E�̍���
	int m_iWinX;		//!< �`��E�B���h�E�ʒux
	int m_iWinY;		//!< �`��E�B���h�E�ʒuy

	// �E�B�W�b�g
	Fl_Menu_Bar *m_pMenuBar;	//!< ���j���[�o�[

	Fl_Spinner *m_pSpinCoefEt;
	Fl_Spinner *m_pSpinMaxEt;
	Fl_Spinner *m_pSpinSpEt;
	Fl_Spinner *m_pSpinSpEtMesh;
	Fl_Spinner *m_pSpinEtCri;
	Fl_Spinner *m_pSpinWaveletScale;
	Fl_Value_Slider *m_pSliderMeshThr;
	Fl_Value_Slider *m_pSliderVScale;
	Fl_Check_Button *m_pCheckRefraction;
	Fl_Check_Button *m_pCheckMesh;

	//�����@Physical Property
	Fl_Value_Slider *m_pSliderIntrPltn;

	Fl_Check_Button *m_pCheckSpears;
	Fl_Spinner* m_pSpinRadiusSpears;

	Fl_Check_Button *m_pCheckWeight;
	Fl_Spinner* m_pSpinWeight;

	Fl_Check_Button *m_pCheckIteration;
	Fl_Value_Slider *m_pSliderItr;

	Fl_Check_Button *m_pCheckPath;

	Fl_Check_Button *m_pCheckStiff;
	Fl_Spinner* m_pSpinStiff;

	//�����_�����O�p�����[�^
	Fl_Spinner* m_pSpinEtaRatio;
	Fl_Spinner* m_pSpinFresnelBias; 
	Fl_Spinner* m_pSpinFresnelPower;
	Fl_Spinner* m_pSpinFresnelScale;

	//�f�o�b�O���[�h
	Fl_Check_Button *m_pCheckDebug;
	Fl_Check_Button* m_pCheckHeatObj;

	rxFlDndBox *m_pDndBox;		//!< D&D�̈�
	Fl_Box *m_pBoxStatus;		//!< �X�e�[�^�X�o�[

	rxFlGLWindow *m_pGLCanvas;	//!< OpenGL�`��L�����p�X

	// �t�@�C�����
	string m_strFileName;		//!< �ǂݍ��񂾃t�@�C����
	string m_strFullPath;		//!< �ǂݍ��񂾃t�@�C����(�t���p�X)

	char *m_pStatusLabel;		//!< �X�e�[�^�X�o�[������

	bool m_bFullScreen;			//!< �t���X�N���[���\��
	int m_iIdle;				//!< �A�C�h��,�^�C�}�[�̏��

public:
	//! �R���X�g���N�^
	rxFlWindow(int w, int h, const char* title);

	//! �f�X�g���N�^
	~rxFlWindow();

public:
	// �R�[���o�b�N�֐�
	static void OnMenuFile_s(Fl_Widget *widget, void* x);
	inline void OnMenuFileOpen(void);
	inline void OnMenuFileSave(void);
	inline void OnMenuFileSaveFrame(void);
	inline void OnMenuFileQuit(void);
	static void OnMenuStep_s(Fl_Widget *widget, void* x);
	inline void OnMenuStep(string label);
	static void OnMenuWindow_s(Fl_Widget *widget, void* x);
	inline void OnMenuWindow(string label);
	static void OnMenuHelpVersion_s(Fl_Widget *widget, void* x);

	static void OnDnd_s(Fl_Widget *widget, void* x);
	inline void OnDnd(void);

	static void OnButtonStart_s(Fl_Widget *widget, void* x);
	inline void OnButtonStart(void);
	static void OnButtonApply_s(Fl_Widget *widget, void* x);
	inline void OnButtonApply(Fl_Widget *widget);
	static void OnButtonClear_s(Fl_Widget *widget, void* x);
	inline void OnButtonClear(Fl_Widget *widget);


	static void OnSpinSimulation_s(Fl_Widget *widget, void* x);
	inline void OnSpinSimulation(Fl_Widget *widget);
	static void OnSliderDraw_s(Fl_Widget *widget, void* x);
	inline void OnSliderDraw(Fl_Widget *widget);

	static void OnCheckDraw_s(Fl_Widget *widget, void* x);
	inline void OnCheckDraw(Fl_Widget *widget);

	//�ǉ�
	static void OnSliderParam_s(Fl_Widget *widget, void* x);
	inline void OnSliderParam(Fl_Widget *widget);

	static void OnCheckMode_Path_s(Fl_Widget *widget, void* x);
	inline void OnCheckMode_Path(Fl_Widget *widget);

	static void OnRadioMode_Spears_s(Fl_Widget *widget, void* x);
	inline void OnRadioMode_Spears(Fl_Widget *widget);

	static void OnSpinRadius_Spears_s(Fl_Widget *widget, void* x);
	inline void OnSpinRadius_Spears(Fl_Widget *widget);

	static void OnRadioMode_Convolution_s(Fl_Widget *widget, void* x);
	inline void OnRadioMode_Convolution(Fl_Widget *widget);

	static void OnSpinDegree_Weight_s(Fl_Widget *widget, void* x);
	inline void OnSpinDegree_Weight(Fl_Widget *widget);

	static void OnSpinDegree_Shader_s(Fl_Widget *widget, void* x);
	inline void OnSpinDegree_Shader(Fl_Widget *widget);

	//static void OnSpinDegree_DirVec_s(Fl_Widget *widget, void* x);
	//inline void OnSpinDegree_DirVec(Fl_Widget *widget);

	static void OnRadioMode_CalcMethod_s(Fl_Widget *widget, void* x);
	inline void OnRadioMode_CalcMethod(Fl_Widget *widget);

	static void OnRadioMode_SimulationMethod_s(Fl_Widget* widget, void* x);
	inline void OnRadioMode_SimulationMethod(Fl_Widget* widget);

	static void OnSliderItrNum_s(Fl_Widget *widget, void* x);
	inline void OnSliderItrNum(Fl_Widget *widget);

	static void OnButtonDebugDeformation_s(Fl_Widget *widget, void* x);
	inline void OnButtonDebugDeformation(Fl_Widget *widget);

	static void OnCheckMode_HeatObj_s(Fl_Widget *widget, void* x);
	inline void OnCheckMode_HeatObj(Fl_Widget *widget);

	static void OnSpinStiff_s(Fl_Widget *widget, void* x);
	inline void OnSpinStiff(Fl_Widget *widget);

	static void OnCheckMode_Debug_s(Fl_Widget *widget, void* x);
	inline void OnCheckMode_Debug(Fl_Widget *widget);

	// ���j���[�̍X�V
	void UpdateMenuState(void);

	// �X�e�[�^�X�o�[������̕ύX
	void SetStatusLabel(const string &label);

	// �t���X�N���[���ؑ�
	void SwitchFullScreen(void);
	bool IsFullScreen(void) const;

	// �t�@�C�����o��
	void Open(const string &fn);
	void Save(const string &fn);
	void OpenImage(const string &fn);

	// �ݒ�t�@�C��
	void ReadConfig(const string &fn);
	void WriteConfig(const string &fn);

protected:
	int handle(int ev);
};





#endif // #ifdef _RX_FLTK_WINDOW_H_
