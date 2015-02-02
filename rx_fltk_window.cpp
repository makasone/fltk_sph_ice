/*!
  @file rx_fltk_window.cpp
	
  @brief FLTKによるウィンドウクラス
 
  @author Makoto Fujisawa 
  @date   2011-08
*/


//-----------------------------------------------------------------------------
// インクルードファイル
//-----------------------------------------------------------------------------
#include "rx_fltk_window.h"
#include "rx_atom_ini.h"

#include "rx_sph.h"

#include <stdio.h>
#include <math.h>

//-----------------------------------------------------------------------------
// 変数
//-----------------------------------------------------------------------------
// 設定ファイル
rxINI *g_pINI = new rxINI;

// デフォルト読み込みファイル
vector<string> g_vDefaultFiles;

// 乱流パラメータ
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



//-----------------------------------------------------------------------------
// rxFlWindowクラスの実装
//-----------------------------------------------------------------------------

//! コンストラクタ
rxFlWindow::rxFlWindow(int w_, int h_, const char* title)
	: Fl_Double_Window(w_, h_, title), m_iWinX(100), m_iWinY(100), m_iWinW(w_), m_iWinH(h_)
{
	m_pStatusLabel = 0;
	m_bFullScreen = false;

	resizable(this);

	int hs_menu = 20;	// メニューバーの高さ
	int hs_para = 150;	// パラメータ調整用ウィジット配置領域の高さ
	int hs_stat = 24;	// ステータスバーの高さ

	int hor_margin = 5;	// 水平方向マージン
	int ver_margin = 5;	// 垂直方向マージン

	int xs = hor_margin;

	begin();
	{
		// 描画領域
		int ys = hs_menu+ver_margin;
		int ws = w()-hor_margin*2;
		int hs = h()-(hs_menu+hs_para+hs_stat+4*ver_margin);
		Fl_Group *g = new Fl_Group(xs, ys, ws, hs, 0);
		g->box(FL_NO_BOX);
		{
			m_pGLCanvas = new rxFlGLWindow(xs+2, ys+2, ws-4, hs-4, 0, this);
			m_pGLCanvas->box(FL_FLAT_BOX);
			m_pGLCanvas->mode(FL_RGB | FL_ALPHA | FL_DOUBLE | FL_DEPTH | FL_MULTISAMPLE);
			m_pGLCanvas->align(FL_ALIGN_INSIDE | FL_ALIGN_CENTER);
			cout << ws << ", " << hs << endl;

			// D&Dボックス
			m_pDndBox = new rxFlDndBox(xs, ys, ws, hs, 0);
			m_pDndBox->callback(OnDnd_s, this);
		}
		g->end();
		Fl_Group::current()->resizable(g);
	}

	ReadConfig(RX_PROGRAM_NAME);

	{
		// メニューバー
		m_pMenuBar = new Fl_Menu_Bar(0, 0, w(), hs_menu, 0);

		// Fileメニュー
		m_pMenuBar->add("File/Open File", FL_CTRL+'f', OnMenuFile_s, this); 
		m_pMenuBar->add("File/Save As", FL_CTRL+FL_SHIFT+'s', OnMenuFile_s, this, FL_MENU_DIVIDER); 
		m_pMenuBar->add("File/Save FrameBuffer ", FL_CTRL+'s', OnMenuFile_s, this, FL_MENU_DIVIDER); 
		m_pMenuBar->add("File/Quit", FL_CTRL+'q', OnMenuFile_s, this); 

		// Drawメニュー
		int count = 0;
		while(RX_DRAW_STR[2*count] != "-1"){
			string label = "Draw/"+RX_DRAW_STR[2*count]+"  ";
			string shortcut = RX_DRAW_STR[2*count+1];
			m_pMenuBar->add(RX_TO_CHAR(label), RX_TO_CHAR(shortcut), rxFlGLWindow::OnMenuDraw_s, m_pGLCanvas, FL_MENU_TOGGLE);
			count++;
		}

		// Simulationメニュー
		m_pMenuBar->add("Simulation/Reset",						'R', rxFlGLWindow::OnMenuSimulation_s, m_pGLCanvas, FL_MENU_DIVIDER);
#ifdef RX_USE_PBD
		m_pMenuBar->add("Simulation/Artificial Pressure",		"^t", rxFlGLWindow::OnMenuSimulation_s, m_pGLCanvas, FL_MENU_TOGGLE);
#endif
		m_pMenuBar->add("Simulation/Wavelet Turbulence",		'w', rxFlGLWindow::OnMenuSimulation_s, m_pGLCanvas, FL_MENU_TOGGLE);
		m_pMenuBar->add("Simulation/SPS Turbulence",			'W', rxFlGLWindow::OnMenuSimulation_s, m_pGLCanvas, FL_MENU_TOGGLE);
		m_pMenuBar->add("Simulation/Vorticity Confinement  ",	'V', rxFlGLWindow::OnMenuSimulation_s, m_pGLCanvas, FL_MENU_TOGGLE);
		m_pMenuBar->add("Simulation/SPH Only",					0, rxFlGLWindow::OnMenuSimulation_s, m_pGLCanvas, FL_MENU_TOGGLE);	//追加			
		m_pMenuBar->add("Simulation/HeatTransfar",				0, rxFlGLWindow::OnMenuSimulation_s, m_pGLCanvas, FL_MENU_TOGGLE);	//追加	
		m_pMenuBar->add("Simulation/ShapeMatching",				0, rxFlGLWindow::OnMenuSimulation_s, m_pGLCanvas, FL_MENU_TOGGLE);	//追加
		m_pMenuBar->add("Simulation/IceStructure",				0, rxFlGLWindow::OnMenuSimulation_s, m_pGLCanvas, FL_MENU_TOGGLE | FL_MENU_DIVIDER);	//追加
		
		m_pMenuBar->add("Simulation/Particle Data Input",		"^i", rxFlGLWindow::OnMenuSimulation_s, m_pGLCanvas, FL_MENU_TOGGLE);
		m_pMenuBar->add("Simulation/Particle Data Output",		"^o", rxFlGLWindow::OnMenuSimulation_s, m_pGLCanvas, FL_MENU_TOGGLE);
		m_pMenuBar->add("Simulation/Mesh Saving",				"^m", rxFlGLWindow::OnMenuSimulation_s, m_pGLCanvas, FL_MENU_TOGGLE);
		m_pMenuBar->add("Simulation/Image Saving",				"^a", rxFlGLWindow::OnMenuSimulation_s, m_pGLCanvas, FL_MENU_TOGGLE);

		// Particleメニュー
		m_pMenuBar->add("Particle/Anisotoropic Kernel  ",	'a', rxFlGLWindow::OnMenuParticle_s, m_pGLCanvas, FL_MENU_TOGGLE);

		// Particle/Colorメニュー
		m_pMenuBar->add("Particle/Color/Ramp",				0, rxFlGLWindow::OnMenuParticle_s, m_pGLCanvas, FL_MENU_RADIO);
		m_pMenuBar->add("Particle/Color/Constant",			0, rxFlGLWindow::OnMenuParticle_s, m_pGLCanvas, FL_MENU_RADIO);
		m_pMenuBar->add("Particle/Color/Density",			0, rxFlGLWindow::OnMenuParticle_s, m_pGLCanvas, FL_MENU_RADIO);
		m_pMenuBar->add("Particle/Color/Energy Spectrum",	0, rxFlGLWindow::OnMenuParticle_s, m_pGLCanvas, FL_MENU_RADIO);
		m_pMenuBar->add("Particle/Color/Pressure",			0, rxFlGLWindow::OnMenuParticle_s, m_pGLCanvas, FL_MENU_RADIO);
		m_pMenuBar->add("Particle/Color/Surface",			0, rxFlGLWindow::OnMenuParticle_s, m_pGLCanvas, FL_MENU_RADIO);
		m_pMenuBar->add("Particle/Color/Temperature",		0, rxFlGLWindow::OnMenuParticle_s, m_pGLCanvas, FL_MENU_RADIO);	//追加
		m_pMenuBar->add("Particle/Color/Ice_Cnct",	  FL_CTRL+FL_SHIFT+'a', rxFlGLWindow::OnMenuParticle_s, m_pGLCanvas, FL_MENU_RADIO);	//追加
		m_pMenuBar->add("Particle/Color/Ice_Calc",	  FL_CTRL+FL_SHIFT+'q', rxFlGLWindow::OnMenuParticle_s, m_pGLCanvas, FL_MENU_RADIO);	//追加
		m_pMenuBar->add("Particle/Color/FAST_PATH",	  FL_CTRL+FL_SHIFT+'p', rxFlGLWindow::OnMenuParticle_s, m_pGLCanvas, FL_MENU_RADIO);	//追加
		m_pMenuBar->add("Particle/Color/SELECTED",	  FL_CTRL+FL_SHIFT+'w', rxFlGLWindow::OnMenuParticle_s, m_pGLCanvas, FL_MENU_RADIO);	//追加
		m_pMenuBar->add("Particle/Color/DEFORMATION", FL_CTRL+FL_SHIFT+'e', rxFlGLWindow::OnMenuParticle_s, m_pGLCanvas, FL_MENU_RADIO);	//追加
		m_pMenuBar->add("Particle/Color/Edge",				0, rxFlGLWindow::OnMenuParticle_s, m_pGLCanvas, FL_MENU_RADIO);	//追加
		m_pMenuBar->add("Particle/Color/None",				0, rxFlGLWindow::OnMenuParticle_s, m_pGLCanvas, FL_MENU_RADIO);

		// Particle/Drawメニュー
		count = 0;
		while(RX_PARTICLE_DRAW[2*count] != "-1"){
			string label = "Particle/Draw/"+RX_PARTICLE_DRAW[2*count]+"  ";
			string shortcut = RX_PARTICLE_DRAW[2*count+1];
			m_pMenuBar->add(RX_TO_CHAR(label), RX_TO_CHAR(shortcut), rxFlGLWindow::OnMenuParticle_s, m_pGLCanvas, FL_MENU_RADIO);
			count++;
		}

		// Solidメニュー
		count = 0;
		while(RXS_STR[2*count] != "-1"){
			string label = "Solid/"+RXS_STR[2*count]+"  ";
			string shortcut = RXS_STR[2*count+1];
			m_pMenuBar->add(RX_TO_CHAR(label), RX_TO_CHAR(shortcut), rxFlGLWindow::OnMenuSolid_s, m_pGLCanvas, FL_MENU_TOGGLE);
			count++;
		}

		// Meshメニュー
		count = 0;
		while(RX_TRIANGULATION_METHOD[2*count] != "-1"){
			string label = "Mesh/"+RX_TRIANGULATION_METHOD[2*count]+"  ";
			string shortcut = RX_TRIANGULATION_METHOD[2*count+1];
			m_pMenuBar->add(RX_TO_CHAR(label), RX_TO_CHAR(shortcut), rxFlGLWindow::OnMenuTriangulation_s, m_pGLCanvas, FL_MENU_RADIO);
			count++;
		}

		// Sceneメニュー
		count = 0;
		while(!m_pGLCanvas->m_vSceneTitles[count].empty() && count < 12){
			string label = "Scene/"+m_pGLCanvas->m_vSceneTitles[count]+"  ";
			m_pMenuBar->add(RX_TO_CHAR(label), FL_F+count+1, rxFlGLWindow::OnMenuScene_s, m_pGLCanvas, FL_MENU_RADIO);
			count++;
		}
		
		// Stepメニュー
		m_pMenuBar->add("Step/Step", ' ', OnMenuStep_s, this); 
		m_pMenuBar->add("Step/Animation  ", 's', OnMenuStep_s, this); 

		// Windowメニュー
		m_pMenuBar->add("Window/FullScreen (Window)  ", FL_CTRL+'f', OnMenuWindow_s, this, FL_MENU_TOGGLE); 
		m_pMenuBar->add("Window/FullScreen (GLCanvas)  ", FL_ALT+FL_Enter, OnMenuWindow_s, this, FL_MENU_TOGGLE); 
		count = 0;
		while(RX_CANVAS_SIZE_STR[2*count] != "-1"){
			string label = "Window/Canvas Size/"+RX_CANVAS_SIZE_STR[2*count]+"  ";
			string shortcut = RX_CANVAS_SIZE_STR[2*count+1];
			m_pMenuBar->add(RX_TO_CHAR(label), RX_TO_CHAR(shortcut), OnMenuWindow_s, this);
			count++;
		}

		// Helpメニュー
		m_pMenuBar->add("Help/Version  ", 0, OnMenuHelpVersion_s, this); 
	}
	{
		// 下部GUI領域
		xs = hor_margin;
		int ys = h()-(hs_stat+ver_margin+hs_para);
		Fl_Scroll *g = new Fl_Scroll(hor_margin, ys, w()-2*hor_margin, hs_para, 0);
		g->box(FL_FLAT_BOX);

		xs += 7;
		ys += 5;
		int ws = 80;
		int hs = 25;

		Fl_Boxtype boxtype = FL_FLAT_BOX;
		Fl_Button *button;

		// Start/Stopボタン
		button = new Fl_Button(xs, ys, ws, hs, "Start/Stop");
		button->callback(OnButtonStart_s, this);
		button->down_box(boxtype);
		button->clear_visible_focus();

		// Applyボタン
		button = new Fl_Button(xs, ys+1*(hs+5), ws, hs, "Apply");
		button->callback(OnButtonApply_s, this);
		button->down_box(boxtype);
		button->clear_visible_focus();

		// Clearボタン
		button = new Fl_Button(xs, ys+2*(hs+5), ws, hs, "Clear");
		button->callback(OnButtonClear_s, this);
		button->down_box(boxtype);
		button->clear_visible_focus();

		//
		// シミュレーションパラメータ
		//
		xs += ws+7;
		ys = h()-(hs_stat+ver_margin+hs_para)+5;
		ws = 300;
		hs = hs_para-6;

		Fl_Group* sg = new Fl_Group(xs, ys, ws, hs, "Draw");
		{
			sg->box(FL_DOWN_BOX);
			sg->labelsize(12);
			sg->align(FL_ALIGN_TOP | FL_ALIGN_INSIDE);

			int dx = 125;
			ys += 30;

			m_pSliderVScale = new Fl_Value_Slider(xs+dx, ys, 145, 25, "Veloc. Scale ");
			m_pSliderVScale->type(1);
			m_pSliderVScale->callback(OnSliderDraw_s, this);
			m_pSliderVScale->minimum(0.0);
			m_pSliderVScale->maximum(0.1);
			m_pSliderVScale->step(0.001);
			m_pSliderVScale->value(m_pGLCanvas->m_fVScale);
			m_pSliderVScale->align(Fl_Align(FL_ALIGN_LEFT));
			m_pSliderVScale->clear_visible_focus();
			
			m_pSliderMeshThr = new Fl_Value_Slider(xs+dx, ys+30, 145, 25, "Mesh Threshold ");
			m_pSliderMeshThr->type(1);
			m_pSliderMeshThr->callback(OnSliderDraw_s, this);
			m_pSliderMeshThr->minimum(0.0);
			m_pSliderMeshThr->maximum(2000.0);
			m_pSliderMeshThr->step(5.0);
			m_pSliderMeshThr->value(m_pGLCanvas->m_fMeshThr);
			m_pSliderMeshThr->align(Fl_Align(FL_ALIGN_LEFT));
			m_pSliderMeshThr->clear_visible_focus();

			m_pCheckMesh = new Fl_Check_Button(xs+dx, ys+60, 25, 25, "Surface Mesh ");
			m_pCheckMesh->down_box(FL_DOWN_BOX);
			m_pCheckMesh->callback(OnCheckDraw_s, this);
			m_pCheckMesh->align(Fl_Align(FL_ALIGN_LEFT));
			m_pCheckMesh->clear_visible_focus();

			m_pCheckRefraction = new Fl_Check_Button(xs+dx, ys+80, 25, 25, "Refraction ");
			m_pCheckRefraction->down_box(FL_DOWN_BOX);
			m_pCheckRefraction->callback(OnCheckDraw_s, this);
			m_pCheckRefraction->align(Fl_Align(FL_ALIGN_LEFT));
			m_pCheckRefraction->clear_visible_focus();

			sg->resizable(NULL);
			sg->end();
		}

		//追加：：物性のパラメータを制御
		xs += ws+5;
		ys = h()-(hs_stat+ver_margin+hs_para)+5;
		ws = 320;

		sg = new Fl_Group(xs, ys, ws, hs, "Physical Property");
		{
			sg->box(FL_DOWN_BOX);
			sg->labelsize(12);
			sg->align(FL_ALIGN_TOP | FL_ALIGN_INSIDE);

			int dx = 145;
			ys += 30;

			//液体と固体の線形補間係数
				//値を数値ではなく相対量として考えるため，スライダーを採用
			m_pSliderIntrPltn = new Fl_Value_Slider(xs+dx, ys, 145, 25, "SPH+SM InterPolation");
			m_pSliderIntrPltn->type(1);
			m_pSliderIntrPltn->callback(OnSliderParam_s, this);
			m_pSliderIntrPltn->minimum(0.0);
			m_pSliderIntrPltn->maximum(1.0);
			m_pSliderIntrPltn->step(0.05);
			m_pSliderIntrPltn->value(1.0);
			m_pSliderIntrPltn->align(Fl_Align(FL_ALIGN_LEFT));
			m_pSliderIntrPltn->clear_visible_focus();

			sg->resizable(NULL);
			sg->end();
		}

		//追加：：熱計算パラメータ
		xs += ws+5;
		ys = h()-(hs_stat+ver_margin+hs_para)+5;
		ws = 280;

		sg = new Fl_Group(xs, ys, ws, hs, "Temparature");
		{
			sg->box(FL_DOWN_BOX);
			sg->labelsize(12);
			sg->align(FL_ALIGN_TOP | FL_ALIGN_INSIDE);

			int dx = 125;
			ys += 30;

			//空気温度
			m_pSliderVScale = new Fl_Value_Slider(xs+dx, ys, 145, 25, "AirTemp");
			m_pSliderVScale->type(1);
			m_pSliderVScale->callback(OnSliderParam_s, this);
			m_pSliderVScale->minimum(-1000);
			m_pSliderVScale->maximum(1000);
			m_pSliderVScale->step(10);
			m_pSliderVScale->value(0);
			m_pSliderVScale->align(Fl_Align(FL_ALIGN_LEFT));
			m_pSliderVScale->clear_visible_focus();

			//TODO::熱拡散係数


			//TODO::伝熱係数

			sg->resizable(NULL);
			sg->end();
		}

		//追加：：シミュレーションモード
		xs += ws+5;
		ys = h()-(hs_stat+ver_margin+hs_para)+5;
		ws = 500;

		sg = new Fl_Group(xs, ys, ws, hs, "Mode");
		{
			sg->box(FL_DOWN_BOX);
			sg->labelsize(12);
			sg->align(FL_ALIGN_TOP | FL_ALIGN_INSIDE);

			ys += 15;

		//各モード
			//各クラスタからの位置決定
			int weight_xs = xs+5;
			Fl_Group *gr = new Fl_Group(weight_xs, ys, 200, 60);
			{
				gr->box(FL_THIN_UP_FRAME);
				Fl_Round_Button *radio;

				//ノーマル
				radio = new Fl_Round_Button(weight_xs, ys, 20, 20, "Normal");
				radio->type(102);
				radio->down_box(FL_ROUND_DOWN_BOX);
				radio->callback(OnRadioMode_Convolution_s, this);
				radio->value(1);

				//重み付き平均
				radio = new Fl_Round_Button(weight_xs, ys+20, 20, 20, "Weight");
				radio->type(102);
				radio->down_box(FL_ROUND_DOWN_BOX);
				radio->callback(OnRadioMode_Convolution_s, this);

					//カーネル関数の次数
					m_pSpinRadiusSpears = new Fl_Spinner(weight_xs+80, ys+20, 70, 20);
					m_pSpinRadiusSpears->type(FL_FLOAT_INPUT);
					m_pSpinRadiusSpears->callback(OnSpinDegree_Weight_s, this);
					m_pSpinRadiusSpears->minimum(0.5);
					m_pSpinRadiusSpears->maximum(4.0);
					m_pSpinRadiusSpears->step(0.1);
					m_pSpinRadiusSpears->value(1.0);
					m_pSpinRadiusSpears->clear_visible_focus();

				//異方性のある運動　方向ベクトルの指定
				radio = new Fl_Round_Button(weight_xs, ys+40, 20, 20, "Anisotropic");
				radio->type(102);
				radio->down_box(FL_ROUND_DOWN_BOX);
				radio->callback(OnRadioMode_Convolution_s, this);

				gr->end();
				gr->resizable(0);
			}

			//運動計算クラスタ選択
			int spears_xs = xs+210;
			gr = new Fl_Group(spears_xs, ys, 160, 60);
			{
				gr->box(FL_THIN_UP_FRAME);
				Fl_Round_Button *radio;

				//ノーマル
				radio = new Fl_Round_Button(spears_xs, ys, 20, 20, "Normal");
				radio->type(102);
				radio->down_box(FL_ROUND_DOWN_BOX);
				radio->callback(OnRadioMode_Spears_s, this);
				radio->value(1);

				//スパース
				radio = new Fl_Round_Button(spears_xs, ys+20, 20, 20, "Spears");
				radio->type(102);
				radio->down_box(FL_ROUND_DOWN_BOX);
				radio->callback(OnRadioMode_Spears_s, this);

					//影響半径 Poission Disk Sampling
					m_pSpinRadiusSpears = new Fl_Spinner(spears_xs+80, ys+20, 70, 20);
					m_pSpinRadiusSpears->type(FL_FLOAT_INPUT);
					m_pSpinRadiusSpears->callback(OnSpinRadius_Spears_s, this);
					m_pSpinRadiusSpears->minimum(0.1);
					m_pSpinRadiusSpears->maximum(2.5);
					m_pSpinRadiusSpears->step(0.05);
					m_pSpinRadiusSpears->value(0.1);
					m_pSpinRadiusSpears->clear_visible_focus();

				//その他
				radio = new Fl_Round_Button(spears_xs, ys+40, 20, 20, "Other");
				radio->type(102);
				radio->down_box(FL_ROUND_DOWN_BOX);
				radio->callback(OnRadioMode_Spears_s, this);

				gr->end();
				gr->resizable(0);
			}

			//運動計算の反復
			int calc_xs = xs+5;
			int calc_ys = ys+65;
			gr = new Fl_Group(calc_xs, calc_ys, 380, 60);
			{
				gr->box(FL_THIN_UP_FRAME);
				Fl_Round_Button *radio;

				//ノーマル
				radio = new Fl_Round_Button(calc_xs, calc_ys, 20, 20, "Normal");
				radio->type(102);
				radio->down_box(FL_ROUND_DOWN_BOX);
				radio->callback(OnRadioMode_CalcMethod_s, this);
				radio->value(1);

				//反復　回数指定
				radio = new Fl_Round_Button(calc_xs, calc_ys+20, 20, 20, "Itr_Num");
				radio->type(102);
				radio->down_box(FL_ROUND_DOWN_BOX);
				radio->callback(OnRadioMode_CalcMethod_s, this);

					//反復回数
					m_pSliderItr = new Fl_Value_Slider(calc_xs+80, calc_ys+20, 110, 20);
					m_pSliderItr->type(1);
					m_pSliderItr->callback(OnSliderItrNum_s, this);
					m_pSliderItr->minimum(1);
					m_pSliderItr->maximum(50);
					m_pSliderItr->step(1);
					m_pSliderItr->value(1);
					m_pSliderItr->align(Fl_Align(FL_ALIGN_LEFT));
					m_pSliderItr->clear_visible_focus();

				//反復　許容変形量指定
				radio = new Fl_Round_Button(calc_xs, calc_ys+40, 20, 20, "Itr_Stiff");
				radio->type(102);
				radio->down_box(FL_ROUND_DOWN_BOX);
				radio->callback(OnRadioMode_CalcMethod_s, this);

					//許容変形量
					m_pSpinStiff = new Fl_Spinner(calc_xs+80, calc_ys+40, 70, 20);
					m_pSpinStiff->type(FL_FLOAT_INPUT);
					m_pSpinStiff->callback(OnSpinStiff_s, this);
					m_pSpinStiff->minimum(0.01);
					m_pSpinStiff->maximum(100.0);
					m_pSpinStiff->step(0.01);
					m_pSpinStiff->value(20.0);
					m_pSpinStiff->clear_visible_focus();

				calc_xs += 205;

				//反復　クラスタの粒子増加
				radio = new Fl_Round_Button(calc_xs, calc_ys, 20, 20, "Itr_Expand");
				radio->type(102);
				radio->down_box(FL_ROUND_DOWN_BOX);
				radio->callback(OnRadioMode_CalcMethod_s, this);

				gr->end();
				gr->resizable(0);
			}

			//パス
			m_pCheckPath = new Fl_Check_Button(xs+450, ys+85, 25, 25, "Path ");
			m_pCheckPath->down_box(FL_DOWN_BOX);
			m_pCheckPath->callback(OnCheckMode_Path_s, this);
			m_pCheckPath->align(Fl_Align(FL_ALIGN_LEFT));
			m_pCheckPath->clear_visible_focus();

			sg->resizable(NULL);
			sg->end();
		}

		//追加：：デバッグ
		xs += ws+5;
		ys = h()-(hs_stat+ver_margin+hs_para)+5;
		ws = 280;

		sg = new Fl_Group(xs, ys, ws, hs, "Debug");
		{
			sg->box(FL_DOWN_BOX);
			sg->labelsize(12);
			sg->align(FL_ALIGN_TOP | FL_ALIGN_INSIDE);

			int dx = 125;
			ys += 30;

			//デバッグボタン
			m_pCheckDebug = new Fl_Check_Button(xs+dx, ys, 25, 25, "DebugMode");
			m_pCheckDebug->down_box(FL_DOWN_BOX);
			m_pCheckDebug->callback(OnCheckMode_Debug_s, this);
			m_pCheckDebug->align(Fl_Align(FL_ALIGN_LEFT));
			m_pCheckDebug->clear_visible_focus();

			//グラフ作成ボタン
			Fl_Boxtype boxtype = FL_FLAT_BOX;
			Fl_Button *button;

			button = new Fl_Button(xs+dx+50, ys, 80, 25, "MakeGraph");
			button->callback(OnButtonDebugDeformation_s, this);
			button->down_box(boxtype);
			button->clear_visible_focus();

			sg->resizable(NULL);
			sg->end();
		}

		g->resizable(NULL);
		g->end();
	}
	{
		// ステータスバー(Fl_Box)
		int ys = h()-hs_stat;
		m_pBoxStatus = new Fl_Box(0, ys, w(), hs_stat, "status");
		m_pBoxStatus->box(FL_EMBOSSED_BOX);
		m_pBoxStatus->align(FL_ALIGN_INSIDE | FL_ALIGN_WRAP | FL_ALIGN_RIGHT);
		//m_pBoxStatus->color(color());
		m_pBoxStatus->labelfont(FL_HELVETICA_BOLD);
		m_pBoxStatus->clear_visible_focus();
	}
	end();

	resize(m_iWinX, m_iWinY, m_iWinW, m_iWinH);

	UpdateMenuState();

	show();
}

//! デストラクタ
rxFlWindow::~rxFlWindow()
{
	WriteConfig(RX_PROGRAM_NAME);

	if(m_pStatusLabel) delete [] m_pStatusLabel;

	if(m_pMenuBar) delete m_pMenuBar;


	if(m_pSpinCoefEt) delete m_pSpinCoefEt;
	if(m_pSpinMaxEt) delete m_pSpinMaxEt;
	if(m_pSpinSpEt) delete m_pSpinSpEt;
	if(m_pSpinSpEtMesh) delete m_pSpinSpEtMesh;
	if(m_pSpinEtCri) delete m_pSpinEtCri;
	if(m_pSpinWaveletScale) delete m_pSpinWaveletScale;
	if(m_pSliderMeshThr) delete m_pSliderMeshThr;
	if(m_pSliderVScale) delete m_pSliderVScale;
	if(m_pCheckMesh) delete m_pCheckMesh;
	if(m_pCheckRefraction) delete m_pCheckRefraction;

	if(m_pGLCanvas) delete m_pGLCanvas;
	if(m_pDndBox) delete m_pDndBox;
	if(m_pBoxStatus) delete m_pBoxStatus;
}


/*!
 * トグル付メニュー項目の更新
 */
void rxFlWindow::UpdateMenuState(void)
{
	int current, count;
	// Draw - FL_MENU_TOGGLE
	current = m_pGLCanvas->m_iDraw;
	count = 0;
	while(RX_DRAW_STR[2*count] != "-1"){
		string label = "Draw/"+RX_DRAW_STR[2*count]+"  ";
		SetMenuItemState(m_pMenuBar, RX_TO_CHAR(label), ((0x01 << count) & current), 1);
		count++;
	}

	bool turb = m_pGLCanvas->m_bsSimuSetting.at(ID_SPH_PS_TURB);

	// Simulation - FL_MENU_TOGGLE
	SetMenuItemState(m_pMenuBar, "Simulation/Wavelet Turbulence",		(m_pGLCanvas->m_bsSimuSetting.at(ID_SPH_PS_TURB)));
	SetMenuItemState(m_pMenuBar, "Simulation/SPS Turbulence",			(m_pGLCanvas->m_bsSimuSetting.at(ID_SPH_SPS_TURB)));
	SetMenuItemState(m_pMenuBar, "Simulation/SPH Only",					(m_pGLCanvas->m_bsSimuSetting.at(ID_HEAT)));		//追加	
	SetMenuItemState(m_pMenuBar, "Simulation/HeatTransfar",				(m_pGLCanvas->m_bsSimuSetting.at(ID_HEAT)));		//追加	
	SetMenuItemState(m_pMenuBar, "Simulation/ShapeMatching",			(m_pGLCanvas->m_bsSimuSetting.at(ID_SM)));			//追加	
	SetMenuItemState(m_pMenuBar, "Simulation/IceStructure",				(m_pGLCanvas->m_bsSimuSetting.at(ID_ICE)));			//追加	

	SetMenuItemState(m_pMenuBar, "Simulation/Vorticity Confinement  ",	(m_pGLCanvas->m_bsSimuSetting.at(ID_SPH_VC)));
	SetMenuItemState(m_pMenuBar, "Simulation/Particle Data Input",		(m_pGLCanvas->m_bsSimuSetting.at(ID_SPH_INPUT)));
	SetMenuItemState(m_pMenuBar, "Simulation/Particle Data Output",		(m_pGLCanvas->m_bsSimuSetting.at(ID_SPH_OUTPUT)));
	SetMenuItemState(m_pMenuBar, "Simulation/Mesh Saving",				(m_pGLCanvas->m_bsSimuSetting.at(ID_SPH_MESH_OUTPUT)));
	SetMenuItemState(m_pMenuBar, "Simulation/Image Saving",				(m_pGLCanvas->m_iSaveImageSpacing != -1));

	// Particle - FL_MENU_TOGGLE
	SetMenuItemState(m_pMenuBar, "Particle/Anisotoropic Kernel  ",	(m_pGLCanvas->m_bsSimuSetting.at(ID_SPH_ANISOTROPIC)));

	// Particle/Color - FL_MENU_RADIO
	SetMenuItemState(m_pMenuBar, "Particle/Color/Ramp",				(m_pGLCanvas->m_iColorType == rxParticleSystemBase::RX_RAMP));
	SetMenuItemState(m_pMenuBar, "Particle/Color/Density",			(m_pGLCanvas->m_iColorType == rxParticleSystemBase::RX_DENSITY));
	SetMenuItemState(m_pMenuBar, "Particle/Color/Energy Spectrum",	(m_pGLCanvas->m_iColorType == rxParticleSystemBase::RX_ENERGY_SPECTRUM));
	SetMenuItemState(m_pMenuBar, "Particle/Color/Pressure",			(m_pGLCanvas->m_iColorType == rxParticleSystemBase::RX_PRESSURE));
	SetMenuItemState(m_pMenuBar, "Particle/Color/Temperature",		(m_pGLCanvas->m_iColorType == rxParticleSystemBase::RX_TEMP));			//追加
	SetMenuItemState(m_pMenuBar, "Particle/Color/Ice_Connect",		(m_pGLCanvas->m_iColorType == rxParticleSystemBase::RX_ICE_CONNECT));	//追加
	SetMenuItemState(m_pMenuBar, "Particle/Color/Ice_Calc",			(m_pGLCanvas->m_iColorType == rxParticleSystemBase::RX_ICE_CALC));		//追加
	SetMenuItemState(m_pMenuBar, "Particle/Color/Edge",				(m_pGLCanvas->m_iColorType == rxParticleSystemBase::RX_EDGE));			//追加
	SetMenuItemState(m_pMenuBar, "Particle/Color/ICE_FAST_PATH",	(m_pGLCanvas->m_iColorType == rxParticleSystemBase::RX_ICE_FAST_PATH));	//追加
	SetMenuItemState(m_pMenuBar, "Particle/Color/SELECTED",			(m_pGLCanvas->m_iColorType == rxParticleSystemBase::RX_ICE_SELECTED));	//追加	
	SetMenuItemState(m_pMenuBar, "Particle/Color/DEFORMATION",		(m_pGLCanvas->m_iColorType == rxParticleSystemBase::RX_ICE_DEFORMATION));	//追加	
	SetMenuItemState(m_pMenuBar, "Particle/Color/None",				(m_pGLCanvas->m_iColorType == rxParticleSystemBase::RX_NONE));

	// Particle/Color - FL_MENU_RADIO
	current = m_pGLCanvas->m_iDrawPS;
	count = 0;
	while(RX_PARTICLE_DRAW[2*count] != "-1"){
		string label = "Particle/Draw/"+RX_PARTICLE_DRAW[2*count]+"  ";
		SetMenuItemState(m_pMenuBar, RX_TO_CHAR(label), (count == current), 1);
		count++;
	}

	// Solid - FL_MENU_TOGGLE
	current = m_pGLCanvas->m_iSolidDraw;
	count = 0;
	while(RXS_STR[2*count] != "-1"){
		string label = "Solid/"+RXS_STR[2*count]+"  ";
		SetMenuItemState(m_pMenuBar, RX_TO_CHAR(label), ((0x01 << count) & current), 1);
		count++;
	}

	// Mesh - FL_MENU_RADIO
	current = m_pGLCanvas->m_iTriangulationMethod;
	count = 0;
	while(RX_TRIANGULATION_METHOD[2*count] != "-1"){
		string label = "Mesh/"+RX_TRIANGULATION_METHOD[2*count]+"  ";
		SetMenuItemState(m_pMenuBar, RX_TO_CHAR(label), (count == current), 1);
		count++;
	}

	// Scene - FL_MENU_RADIO
	current = m_pGLCanvas->m_iCurrentSceneIdx;
	count = 0;
	while(!m_pGLCanvas->m_vSceneTitles[count].empty() && count < 12){
		string label = "Scene/"+m_pGLCanvas->m_vSceneTitles[count]+"  ";
		SetMenuItemState(m_pMenuBar, RX_TO_CHAR(label), (count == current), 1);
		count++;
	}

	// Windowメニュー
	SetMenuItemState(m_pMenuBar, "Window/FullScreen  ", m_bFullScreen, 1);

	// Refractionチェックボタン
	m_pCheckRefraction->value((m_pGLCanvas->m_iDraw & RXD_REFRAC ? 1 : 0));

	// Meshチェックボタン
	m_pCheckMesh->value((m_pGLCanvas->m_iDraw & RXD_MESH ? 1 : 0));
}

/*!
 * Fl_Buttonのコールバック関数 - Start/Stopボタン
 * @param[in] widget ウィジットの親クラスオブジェクト
 * @param[in] x ユーザ定義変数
 */
void rxFlWindow::OnButtonStart_s(Fl_Widget *widget, void* x)
{
	((rxFlWindow*)x)->OnButtonStart();
}
void rxFlWindow::OnButtonStart(void)
{
	m_pGLCanvas->SwitchIdle(-1);
}

/*!
 * Fl_Buttonのコールバック関数 - Applyボタン
 * @param[in] widget ウィジットの親クラスオブジェクト
 * @param[in] x ユーザ定義変数
 */
void rxFlWindow::OnButtonApply_s(Fl_Widget *widget, void* x)
{
	((rxFlWindow*)x)->OnButtonApply(widget);
}
void rxFlWindow::OnButtonApply(Fl_Widget *widget)
{
	g_fCoefEt		   = m_pSpinCoefEt->value();
	g_fMaxEt		   = m_pSpinMaxEt->value();
	g_fCoefTurb		   = m_pSpinSpEt->value();
	g_fCoefTurbForMesh = m_pSpinSpEtMesh->value();
	g_fEtCri		   = m_pSpinEtCri->value();
	g_fWaveletScale	   = m_pSpinWaveletScale->value();

	m_pGLCanvas->m_fMeshThr = m_pSliderMeshThr->value();
	m_pGLCanvas->m_fVScale = m_pSliderVScale->value();
}

void rxFlWindow::OnButtonClear_s(Fl_Widget *widget, void* x)
{
	((rxFlWindow*)x)->OnButtonClear(widget);
}

void rxFlWindow::OnButtonClear(Fl_Widget *widget)
{
	m_pGLCanvas->ResetSimulation();
}

/*!
 * Fl_Spinnerのコールバック関数 - Simulation
 * @param[in] widget ウィジットの親クラスオブジェクト
 * @param[in] x ユーザ定義変数
 */
void rxFlWindow::OnSpinSimulation_s(Fl_Widget *widget, void* x)
{
	((rxFlWindow*)x)->OnSpinSimulation(widget);
}
void rxFlWindow::OnSpinSimulation(Fl_Widget *widget)
{
	Fl_Spinner *spin = (Fl_Spinner*)widget;
	string label = spin->label();
	double val = spin->value();

	if(label.find("Ecoef") != string::npos){
		g_fCoefEt = m_pSpinCoefEt->value();
	}
	else if(label.find("Emax") != string::npos){
		g_fMaxEt = m_pSpinMaxEt->value();
	}
	else if(label.find("Scale") != string::npos){
		g_fWaveletScale = m_pSpinWaveletScale->value();
	}
	else if(label.find("SPS Scale") != string::npos){
		g_fCoefTurb	= m_pSpinSpEt->value();
	}
	else if(label.find("SPS Mesh Scale") != string::npos){
		g_fCoefTurbForMesh = m_pSpinSpEtMesh->value();
	}
	else if(label.find("SPS Ecri") != string::npos){
		g_fEtCri = m_pSpinEtCri->value();
	}
}

/*!
 * Fl_Value_Sliderのコールバック関数 - Draw
 * @param[in] widget ウィジットの親クラスオブジェクト
 * @param[in] x ユーザ定義変数
 */
void rxFlWindow::OnSliderDraw_s(Fl_Widget *widget, void* x)
{
	((rxFlWindow*)x)->OnSliderDraw(widget);
}
void rxFlWindow::OnSliderDraw(Fl_Widget *widget)
{
	Fl_Value_Slider *slider = (Fl_Value_Slider*)widget;
	string label = slider->label();
	double val = slider->value();

	if(label.find("Veloc. Scale") != string::npos){
		m_pGLCanvas->m_fVScale = val;
	}
	else if(label.find("Mesh Threshold") != string::npos){
		m_pGLCanvas->m_fMeshThr = val;
	}
}

/*!
 * Fl_Value_Sliderのコールバック関数
 * @param[in] widget ウィジットの親クラスオブジェクト
 * @param[in] x ユーザ定義変数
 */
void rxFlWindow::OnSliderParam_s(Fl_Widget *widget, void* x)
{
	((rxFlWindow*)x)->OnSliderParam(widget);
}
void rxFlWindow::OnSliderParam(Fl_Widget *widget)
{
	Fl_Value_Slider *slider = (Fl_Value_Slider*)widget;
	string label = slider->label();
	double val = slider->value();

	if(label.find("AirTemp") != string::npos)
	{
		m_pGLCanvas->m_iceObj->SetAirTemp(val);
	}
	//else if(label.find("Mesh Threshold") != string::npos){
	//	m_pGLCanvas->m_fMeshThr = val;
	//}
	else if(label.find("SPH+SM InterPolation") != string::npos)
	{
		for(int i = 0; i < m_pGLCanvas->m_iIcePrtNum; i++)
		{
			m_pGLCanvas->m_iceObj->SetInterPolationCff(i, val);
		}
	}
}

//モード切替のチェックボックスのコールバック関数
//高速化パス
void rxFlWindow::OnCheckMode_Path_s(Fl_Widget *widget, void* x)
{
	((rxFlWindow*)x)->OnCheckMode_Path(widget);
}
void rxFlWindow::OnCheckMode_Path(Fl_Widget *widget)
{
	Fl_Check_Button* check = (Fl_Check_Button*)widget;
	bool flag = check->value();

	cout << __FUNCTION__ << "flag = " << flag << endl;

	if(flag){	m_pGLCanvas->m_iceObj->ChangeMode_ClusterMove_Path();	}
	else	{	m_pGLCanvas->m_iceObj->ChangeMode_ClusterMove_Normal();	}
}

//スパースな運動計算
void rxFlWindow::OnRadioMode_Spears_s(Fl_Widget *widget, void* x)
{
	((rxFlWindow*)x)->OnRadioMode_Spears(widget);
}

void rxFlWindow::OnRadioMode_Spears(Fl_Widget *widget)
{
	Fl_Round_Button* check = (Fl_Round_Button*)widget;
	int type = check->value();
	string label = check->label();

	cout << __FUNCTION__ << "label  " << label << endl;

	if(label == "Normal")
	{
		m_pGLCanvas->m_iceObj->ChangeMode_JudgeMove_Normal();
		m_pGLCanvas->m_iceObj->ChangeMode_IntrpJudge_Normal();
	}
	else if(label == "Spears")
	{	
		m_pGLCanvas->m_iceObj->ChangeMode_JudgeMove_Spears();
		m_pGLCanvas->m_iceObj->ChangeMode_IntrpJudge_Spears();
	}
}

//影響半径
void rxFlWindow::OnSpinRadius_Spears_s(Fl_Widget *widget, void* x)
{
	((rxFlWindow*)x)->OnSpinRadius_Spears(widget);
}
void rxFlWindow::OnSpinRadius_Spears(Fl_Widget *widget)
{
	Fl_Spinner* check = (Fl_Spinner*)widget;
	double radius = check->value();

	cout << __FUNCTION__ << " radius = " << radius << endl;

	m_pGLCanvas->m_iceObj->SetSelectRadius(radius);
	m_pGLCanvas->m_iceObj->InitSelectCluster();
}

//重み付き平均
void rxFlWindow::OnRadioMode_Convolution_s(Fl_Widget *widget, void* x)
{
	((rxFlWindow*)x)->OnRadioMode_Convolution(widget);
}

void rxFlWindow::OnRadioMode_Convolution(Fl_Widget *widget)
{
	Fl_Round_Button* check = (Fl_Round_Button*)widget;
	int type = check->value();
	string label = check->label();

	cout << __FUNCTION__ << "label  " << label << endl;

	if(label == "Normal")
	{
		m_pGLCanvas->m_iceObj->ChangeMode_Convolution_Normal();
	}
	else if(label == "Weight")
	{	
		m_pGLCanvas->m_iceObj->ChangeMode_Convolution_Weight();		
	}
	else if(label == "Anisotropic")
	{
		m_pGLCanvas->m_iceObj->ChangeMode_Convolution_Anisotropic();	
	}
}

//関数の次数
void rxFlWindow::OnSpinDegree_Weight_s(Fl_Widget *widget, void* x)
{
	((rxFlWindow*)x)->OnSpinDegree_Weight(widget);
}

void rxFlWindow::OnSpinDegree_Weight(Fl_Widget *widget)
{
	Fl_Spinner* check = (Fl_Spinner*)widget;
	double degree = check->value();

	const char* convoWeightName = typeid(Ice_Convolution_Weight).name();
	const char* convoNowName = typeid(*(m_pGLCanvas->m_iceObj->GetConvoObj())).name();

	//現在がWeightモードかの確認　strcmpは文字列が一致すると0になる
	if(strcmp(convoWeightName, convoNowName) != 0) return ;
	cout << __FUNCTION__ << " degree = " << degree << endl;
	cout << convoWeightName << endl;
	cout << convoNowName << endl;

	Ice_Convolution_Weight* convoObj = (Ice_Convolution_Weight*)(m_pGLCanvas->m_iceObj->GetConvoObj());
	convoObj->SetKernelDegree(degree);
}

//反復処理モード
void rxFlWindow::OnRadioMode_CalcMethod_s(Fl_Widget *widget, void* x)
{
	((rxFlWindow*)x)->OnRadioMode_CalcMethod(widget);
}
void rxFlWindow::OnRadioMode_CalcMethod(Fl_Widget *widget)
{
	Fl_Round_Button* check = (Fl_Round_Button*)widget;
	int type = check->value();
	string label = check->label();

	cout << __FUNCTION__ << "label  " << label << endl;

	if(label == "Normal")
	{
		m_pGLCanvas->m_iceObj->ChangeMode_CalcMethod_Normal();
	}
	else if(label == "Itr_Num")
	{	
		m_pGLCanvas->m_iceObj->ChangeMode_CalcMethod_Itr_Num();		
	}
	else if(label == "Itr_Stiff")
	{
		m_pGLCanvas->m_iceObj->ChangeMode_CalcMethod_Itr_Stiff();
	}
	else if(label == "Itr_Expand")
	{
		m_pGLCanvas->m_iceObj->ChangeMode_CalcMethod_Itr_Expand();
	}
}

//反復回数
void rxFlWindow::OnSliderItrNum_s(Fl_Widget *widget, void* x)
{
	((rxFlWindow*)x)->OnSliderItrNum(widget);
}

void rxFlWindow::OnSliderItrNum(Fl_Widget *widget)
{
	Fl_Value_Slider *slider = (Fl_Value_Slider*)widget;
	double val = slider->value();

	Ice_SM::SetIterationNum(val);
}

//許容変形量　硬さ
void rxFlWindow::OnSpinStiff_s(Fl_Widget *widget, void* x)
{
	((rxFlWindow*)x)->OnSpinStiff(widget);
}
void rxFlWindow::OnSpinStiff(Fl_Widget *widget)
{
	Fl_Spinner* check = (Fl_Spinner*)widget;
	double stiffness = check->value();

	cout << __FUNCTION__ << " stiffness = " << stiffness << endl;

	Ice_SM::SetItrStiffness(stiffness);
}

//デバッグモード
void rxFlWindow::OnCheckMode_Debug_s(Fl_Widget *widget, void* x)
{
	((rxFlWindow*)x)->OnCheckMode_Debug(widget);
}

void rxFlWindow::OnCheckMode_Debug(Fl_Widget *widget)
{	cout << __FUNCTION__ << endl;

	Fl_Check_Button* check = (Fl_Check_Button*)widget;
	bool flag = check->value();

	cout << __FUNCTION__ << "flag = " << flag << endl;

	if(flag)
	{
		m_pGLCanvas->m_iceObj->ChangeMode_Debug();		
	}
	else
	{	
		m_pGLCanvas->m_iceObj->ChangeMode_Normal();
	}
}

//変形量データのグラフ作成
void rxFlWindow::OnButtonDebugDeformation_s(Fl_Widget *widget, void* x)
{
	((rxFlWindow*)x)->OnButtonDebugDeformation(widget);
}

void rxFlWindow::OnButtonDebugDeformation(Fl_Widget *widget)
{
	cout << __FUNCTION__ << endl;

	m_pGLCanvas->m_iceObj->DebugMakeGraph();
}

/*!
 * Fl_Check_Buttonのコールバック関数 - Draw
 * @param[in] widget ウィジットの親クラスオブジェクト
 * @param[in] x ユーザ定義変数
 */
void rxFlWindow::OnCheckDraw_s(Fl_Widget *widget, void* x)
{
	((rxFlWindow*)x)->OnCheckDraw(widget);
}
void rxFlWindow::OnCheckDraw(Fl_Widget *widget)
{
	Fl_Check_Button *check = (Fl_Check_Button*)widget;
	string label = check->label();
	double val = check->value();

	if(label.find("Refraction") != string::npos){
		m_pGLCanvas->m_iDraw ^= RXD_REFRAC;
	}
	else if(label.find("Surface Mesh") != string::npos){
		m_pGLCanvas->m_iDraw ^= RXD_MESH;
		m_pGLCanvas->SetMeshCreation();
	}
	//運動計算モード
	else if(label.find("Path") != string::npos)
	{
		cout << "Mode Path" << endl;
		//m_pGLCanvas->m_iceObj->
	}
	else if(label.find("Iteration") != string::npos)
	{
		cout << "Mode Iteration" << endl;
	}
	else if(label.find("Weight") != string::npos)
	{
		cout << "Mode Weight" << endl;
	}
	else if(label.find("Spears") != string::npos)
	{
		cout << "Mode Spears" << endl;
	}

	UpdateMenuState();
}


/*!
 * rxFlDndBoxのコールバック関数
 * @param[in] widget ウィジットの親クラスオブジェクト
 * @param[in] x ユーザ定義変数
 */
void rxFlWindow::OnDnd_s(Fl_Widget *widget, void* x)
{
	((rxFlWindow*)x)->OnDnd();
}
void rxFlWindow::OnDnd(void)
{
	if(m_pDndBox->Event() == FL_PASTE){
		int dnd_text_len = m_pDndBox->EventTextLength();
		string dnd_text = m_pDndBox->EventText();

		// テキストを\nで分割
		vector<string> fns;
		size_t pos0 = 0, pos1 = 0;
		do{
			pos1 = dnd_text.find('\n', pos0);
			string fn = dnd_text.substr(pos0, pos1-pos0);
			if(fn.empty()) break;

			fns.push_back(fn);
			pos0 = pos1+1;
		}while(pos1 != string::npos);

		int n = (int)fns.size();
		for(int i = 0; i < n; ++i){
			cout << "file" << i << " : " << fns[i] << endl;
		}

		if(n == 0) return;

		for(int i = 0; i < n; ++i){
			Open(fns[i]);
		}
	}
}

/*!
 * Fileメニューコールバック関数
 * @param[in] widget ウィジットの親クラスオブジェクト
 * @param[in] x ユーザ定義変数
 */
void rxFlWindow::OnMenuFile_s(Fl_Widget *widget, void* x)
{
	Fl_Menu_Bar *menubar = (Fl_Menu_Bar*)widget;
	char picked[80];
	menubar->item_pathname(picked, sizeof(picked)-1);	// メニュー名

	string label = picked;
	string menu_name = "Draw/";
	label = label.substr(menu_name.size(), string::npos);

	if(label.find("Open File") == 0){
		((rxFlWindow*)x)->OnMenuFileOpen();
	}
	else if(label.find("Save FrameBuffer") == 0){
		((rxFlWindow*)x)->OnMenuFileSaveFrame();
	}
	else if(label.find("Save As") == 0){
		((rxFlWindow*)x)->OnMenuFileSave();
	}
	else if(label.find("Quit") == 0){
		((rxFlWindow*)x)->OnMenuFileQuit();
	}
}

void rxFlWindow::OnMenuFileOpen(void)
{
	string filter = "3D Files (*.obj;*.dxf;*.wrl;*.3ds;*.stl;*.ply)|*.obj;*.dxf;*.wrl;*.3ds;*.stl;*.ply|All Files|*.*";
	vector<string> fns;
	int n = ShowFileDialog(fns, "Open 3D file", filter, false);

	if(n > 0){
		for(int i = 0; i < n; ++i){
			Open(fns[i]);
		}
	}
}
void rxFlWindow::OnMenuFileSave(void)
{
	string filter = "3D Files (*.obj;*.dxf;*.wrl;*.3ds;*.stl;*.ply)|*.obj;*.dxf;*.wrl;*.3ds;*.stl;*.ply|All Files|*.*";
	vector<string> fns;
	int n = ShowFileDialog(fns, "Save 3D file", filter, false);

	if(n > 0){
		for(int i = 0; i < n; ++i){
			Save(fns[i]);
		}
	}
}
void rxFlWindow::OnMenuFileSaveFrame(void)
{
	string filter = "Image Files (*.png;*.bmp)|*.bmp;*.png|All Files|*.*";
	vector<string> fns;
	int n = ShowFileDialog(fns, "Save FrameBuffer", filter, false);

	if(n > 0){
		for(int i = 0; i < n; ++i){
			string ext = GetExtension(fns[i]);
			if(ext.empty()){
				fns[i] += ".png";
			}

			Save(fns[i]);
		}
	}
}
void rxFlWindow::OnMenuFileQuit(void)
{
	WriteConfig(RX_PROGRAM_NAME);
	exit(0);
}


/*!
 * Drawメニューのコールバック関数
 * @param[in] widget ウィジットの親クラスオブジェクト
 * @param[in] x ユーザ定義変数
 */
void rxFlWindow::OnMenuStep_s(Fl_Widget *widget, void* x)
{
	Fl_Menu_Bar *menubar = (Fl_Menu_Bar*)widget;
	char picked[80];
	menubar->item_pathname(picked, sizeof(picked)-1);	// メニュー名

	string label = picked;
	string menu_name = "Step/";
	label = label.substr(menu_name.size(), string::npos);
	
	((rxFlWindow*)x)->OnMenuStep(label);
}
void rxFlWindow::OnMenuStep(string label)
{
	if(label.find("Step") == 0){// アニメーション1ステップだけ実行
		m_pGLCanvas->Idle();
		m_pGLCanvas->m_bFall = !m_pGLCanvas->m_bFall;
	}
	else if(label.find("Animation") == 0){	// アニメーションON/OFF
		m_pGLCanvas->SwitchIdle(-1);
	}
}


/*!
 * Windowメニューのコールバック関数
 * @param[in] widget ウィジットの親クラスオブジェクト
 * @param[in] x ユーザ定義変数
 */
void rxFlWindow::OnMenuWindow_s(Fl_Widget *widget, void* x)
{
	Fl_Menu_Bar *menubar = (Fl_Menu_Bar*)widget;
	char picked[80];
	menubar->item_pathname(picked, sizeof(picked)-1);	// メニュー名

	string label = picked;
	string menu_name = "Window/";
	label = label.substr(menu_name.size(), string::npos);
	
	((rxFlWindow*)x)->OnMenuWindow(label);
}
void rxFlWindow::OnMenuWindow(string label)
{
	if(label.find("FullScreen (Window)") == 0){	// フルスクリーンON/OFF
		SwitchFullScreen();
	}
	else if(label.find("FullScreen (GLCanvas)") == 0){	// フルスクリーンON/OFF
		m_pGLCanvas->SwitchFullScreen(0);
	}
	else if(label.find("Canvas Size/") == 0){	// キャンバスサイズ
		string menu_label = "Canvas Size/";
		string size_str = label.substr(menu_label.size(), string::npos);
		int canvas_w = atoi(size_str.substr(0, size_str.find("x")).c_str());
		int canvas_h = atoi(size_str.substr(size_str.find("x")+1, string::npos).c_str());

		int new_win_w = canvas_w+(w()-m_pGLCanvas->w());
		int new_win_h = canvas_h+(h()-m_pGLCanvas->h());

		resize(x(), y(), new_win_w, new_win_h);
	}
}

/*!
 * メニュー:Help -> Versionのコールバック関数
 * @param[in] widget ウィジットの親クラスオブジェクト
 * @param[in] x ユーザ定義変数
 */
void rxFlWindow::OnMenuHelpVersion_s(Fl_Widget *widget, void* x)
{
	fl_message("OpenGL Application by FLTK\n  version 1.0");
}

/*!
 * ステータスバーに文字列を設定
 * @param[in] label 表示文字列
 */
void rxFlWindow::SetStatusLabel(const string &label)
{
	if(m_pStatusLabel) delete [] m_pStatusLabel;
	m_pStatusLabel = RX_TO_CHAR(label);
	m_pBoxStatus->label(m_pStatusLabel);
}

/*!
 * ファイル読み込み
 * @param[in] fn ファイルパス
 */
void rxFlWindow::Open(const string &fn)
{
	string ext = GetExtension(fn);

	if(ext == "obj" || ext == "dxf" || ext == "wrl" || ext == "3ds" || ext == "stl" || ext == "ply"){
		m_pGLCanvas->OpenFile(fn);
	}
	else if(ext == "bmp" || ext == "jpg" || ext == "png" || ext == "gif" || ext == "tif"){
		return;
	}

	// 読み込んだファイル情報を格納
	m_strFullPath = fn;
	m_strFileName = GetFileName(fn);

	// ファイル名をステータスバーに表示
	SetStatusLabel(m_strFileName);
}

/*!
 * ファイル書き込み
 * @param[in] fn ファイルパス
 */
void rxFlWindow::Save(const string &fn)
{
	string ext = GetExtension(fn);

	if(ext == "obj" || ext == "dxf" || ext == "wrl" || ext == "3ds" || ext == "stl" || ext == "ply"){
		m_pGLCanvas->SaveFile(fn);

		// 読み込んだファイル情報を格納
		m_strFullPath = fn;
		m_strFileName = GetFileName(fn);

		// ファイル名をステータスバーに表示
		SetStatusLabel(m_strFileName);
	}
	else if(ext == "bmp" || ext == "png"){
		m_pGLCanvas->SaveDisplay(fn);
	}
}



/*!
 * 設定ファイル読み込み
 * @param[in] fn 設定ファイル名(拡張子無し)
 */
void rxFlWindow::ReadConfig(const string &fn)
{
	// アプリケーション管理設定ファイル
	if(!g_pINI) g_pINI = new rxINI();

	g_pINI->Set("window", "width",  &m_iWinW, m_iWinW);
	g_pINI->Set("window", "height", &m_iWinH, m_iWinH);
	g_pINI->Set("window", "pos_x",  &m_iWinX, m_iWinX);
	g_pINI->Set("window", "pos_y",  &m_iWinY, m_iWinY);

	if(!(g_pINI->Load(fn+".ini"))){
		cout << "Failed opening the " << fn << ".ini file!" << endl;
	}
}
/*!
 * 設定ファイル書き込み
 * @param[in] fn 設定ファイル名(拡張子無し)
 */
void rxFlWindow::WriteConfig(const string &fn)
{
	m_iWinW = w();
	m_iWinH = h();
	m_iWinX = x();
	m_iWinY = y();
	if(g_pINI->Save(fn+".ini")){
		cout << "save : " << fn << ".ini" << endl;
	}
}

/*!
 * フルスクリーン/ウィンドウ表示の切り替え
 */
void rxFlWindow::SwitchFullScreen(void)
{
	static int pos0[2] = { 0, 0 };
	static int win0[2] = { 500, 500 };
	if(m_bFullScreen){
		fullscreen_off(pos0[0], pos0[1], win0[0], win0[1]);
		m_bFullScreen = false;
	}
	else{
		pos0[0] = x();
		pos0[1] = y();
		win0[0] = w();
		win0[1] = h();
		fullscreen();
		m_bFullScreen = true;
	}

}

/*!
 * フルスクリーン/ウィンドウ表示の状態取得
 */
bool rxFlWindow::IsFullScreen(void) const
{
	return m_bFullScreen;
}

/*!
 * イベントハンドラ
 * @param[in] ev イベントID
 */
int rxFlWindow::handle(int ev)
{
	switch(ev){
	case FL_DND_ENTER:
	case FL_DND_RELEASE:
	case FL_DND_LEAVE:
	case FL_DND_DRAG:
	case FL_PASTE:
		return 1;

	case FL_PUSH:		// マウスボタンダウン
		m_pGLCanvas->Mouse(Fl::event_button(), 1, Fl::event_x(), Fl::event_y());
		break;
	case FL_RELEASE:	// マウスボタンアップ
		m_pGLCanvas->Mouse(Fl::event_button(), 0, Fl::event_x(), Fl::event_y());
		break;
	case FL_DRAG:		// マウスドラッグ
		m_pGLCanvas->Motion(Fl::event_x(), Fl::event_y());
		break;
	case FL_MOVE:		// マウス移動
		m_pGLCanvas->PassiveMotion(Fl::event_x(), Fl::event_y());
		break;

	case FL_KEYBOARD:	// キーダウン
		m_pGLCanvas->Keyboard(Fl::event_key(), Fl::event_x(), Fl::event_y());
		UpdateMenuState();
		break;

	case FL_SHORTCUT:	// グローバルショートカット
		break;

	default:
		break;
	}

	return Fl_Window::handle(ev);
}