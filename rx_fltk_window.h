/*!
  @file rx_fltk_window.h
	
  @brief FLTKによるウィンドウクラス
 
  @author Makoto Fujisawa 
  @date   2011-08
*/

#ifndef _RX_FLTK_WINDOW_H_
#define _RX_FLTK_WINDOW_H_


//-----------------------------------------------------------------------------
// インクルードファイル
//-----------------------------------------------------------------------------
#include <FL/Fl.H>
#include <FL/Fl_Double_Window.H>
#include <FL/Fl_Menu_Bar.H>
#include <FL/Fl_Button.H>
#include <FL/Fl_Scroll.H>
#include <FL/Fl_Box.H>
#include <FL/filename.H>		// ファイル名
#include <FL/Fl_File_Chooser.H>
#include <FL/Fl_File_Icon.H>


#include "rx_fltk_glcanvas.h"

using namespace std;

//-----------------------------------------------------------------------------
// MARK:定義/定数
//-----------------------------------------------------------------------------
const string RX_PROGRAM_NAME = "rx_fltk_sph";
extern vector<string> g_vDefaultFiles;



//-----------------------------------------------------------------------------
//! rxFlWindowクラス - fltkによるウィンドウ
//-----------------------------------------------------------------------------
class rxFlWindow : public Fl_Double_Window
{
protected:
	// メンバ変数
	int m_iWinW;		//!< 描画ウィンドウの幅
	int m_iWinH;		//!< 描画ウィンドウの高さ
	int m_iWinX;		//!< 描画ウィンドウ位置x
	int m_iWinY;		//!< 描画ウィンドウ位置y

	// ウィジット
	Fl_Menu_Bar *m_pMenuBar;	//!< メニューバー

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

	//物性　Physical Property
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

	//レンダリングパラメータ
	Fl_Spinner* m_pSpinEtaRatio;
	Fl_Spinner* m_pSpinFresnelBias; 
	Fl_Spinner* m_pSpinFresnelPower;
	Fl_Spinner* m_pSpinFresnelScale;

	//デバッグモード
	Fl_Check_Button *m_pCheckDebug;
	Fl_Check_Button* m_pCheckHeatObj;

	rxFlDndBox *m_pDndBox;		//!< D&D領域
	Fl_Box *m_pBoxStatus;		//!< ステータスバー

	rxFlGLWindow *m_pGLCanvas;	//!< OpenGL描画キャンパス

	// ファイル情報
	string m_strFileName;		//!< 読み込んだファイル名
	string m_strFullPath;		//!< 読み込んだファイル名(フルパス)

	char *m_pStatusLabel;		//!< ステータスバー文字列

	bool m_bFullScreen;			//!< フルスクリーン表示
	int m_iIdle;				//!< アイドル,タイマーの状態

public:
	//! コンストラクタ
	rxFlWindow(int w, int h, const char* title);

	//! デストラクタ
	~rxFlWindow();

public:
	// コールバック関数
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

	//追加
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

	// メニューの更新
	void UpdateMenuState(void);

	// ステータスバー文字列の変更
	void SetStatusLabel(const string &label);

	// フルスクリーン切替
	void SwitchFullScreen(void);
	bool IsFullScreen(void) const;

	// ファイル入出力
	void Open(const string &fn);
	void Save(const string &fn);
	void OpenImage(const string &fn);

	// 設定ファイル
	void ReadConfig(const string &fn);
	void WriteConfig(const string &fn);

protected:
	int handle(int ev);
};





#endif // #ifdef _RX_FLTK_WINDOW_H_
