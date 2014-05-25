/*!
  @file rx_fltk_window.cpp
	
  @brief FLTK�ɂ��E�B���h�E�N���X
 
  @author Makoto Fujisawa 
  @date   2011-08
*/


//-----------------------------------------------------------------------------
// �C���N���[�h�t�@�C��
//-----------------------------------------------------------------------------
#include "rx_fltk_window.h"
#include "rx_atom_ini.h"

#include "rx_sph.h"

//-----------------------------------------------------------------------------
// �ϐ�
//-----------------------------------------------------------------------------
// �ݒ�t�@�C��
rxINI *g_pINI = new rxINI;

// �f�t�H���g�ǂݍ��݃t�@�C��
vector<string> g_vDefaultFiles;

// �����p�����[�^
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
// rxFlWindow�N���X�̎���
//-----------------------------------------------------------------------------

//! �R���X�g���N�^
rxFlWindow::rxFlWindow(int w_, int h_, const char* title)
	: Fl_Double_Window(w_, h_, title), m_iWinX(100), m_iWinY(100), m_iWinW(w_), m_iWinH(h_)
{
	m_pStatusLabel = 0;
	m_bFullScreen = false;

	resizable(this);

	int hs_menu = 20;	// ���j���[�o�[�̍���
	int hs_para = 150;	// �p�����[�^�����p�E�B�W�b�g�z�u�̈�̍���
	int hs_stat = 24;	// �X�e�[�^�X�o�[�̍���

	int hor_margin = 5;	// ���������}�[�W��
	int ver_margin = 5;	// ���������}�[�W��

	int xs = hor_margin;

	begin();
	{
		// �`��̈�
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

			// D&D�{�b�N�X
			m_pDndBox = new rxFlDndBox(xs, ys, ws, hs, 0);
			m_pDndBox->callback(OnDnd_s, this);
		}
		g->end();
		Fl_Group::current()->resizable(g);
	}

	ReadConfig(RX_PROGRAM_NAME);

	{
		// ���j���[�o�[
		m_pMenuBar = new Fl_Menu_Bar(0, 0, w(), hs_menu, 0);

		// File���j���[
		m_pMenuBar->add("File/Open File", FL_CTRL+'f', OnMenuFile_s, this); 
		m_pMenuBar->add("File/Save As", FL_CTRL+FL_SHIFT+'s', OnMenuFile_s, this, FL_MENU_DIVIDER); 
		m_pMenuBar->add("File/Save FrameBuffer ", FL_CTRL+'s', OnMenuFile_s, this, FL_MENU_DIVIDER); 
		m_pMenuBar->add("File/Quit", FL_CTRL+'q', OnMenuFile_s, this); 

		// Draw���j���[
		int count = 0;
		while(RX_DRAW_STR[2*count] != "-1"){
			string label = "Draw/"+RX_DRAW_STR[2*count]+"  ";
			string shortcut = RX_DRAW_STR[2*count+1];
			m_pMenuBar->add(RX_TO_CHAR(label), RX_TO_CHAR(shortcut), rxFlGLWindow::OnMenuDraw_s, m_pGLCanvas, FL_MENU_TOGGLE);
			count++;
		}

		// Simulation���j���[
		m_pMenuBar->add("Simulation/Reset",						'R', rxFlGLWindow::OnMenuSimulation_s, m_pGLCanvas, FL_MENU_DIVIDER);
#ifdef RX_USE_PBD
		m_pMenuBar->add("Simulation/Artificial Pressure",		"^t", rxFlGLWindow::OnMenuSimulation_s, m_pGLCanvas, FL_MENU_TOGGLE);
#endif
		m_pMenuBar->add("Simulation/Wavelet Turbulence",		'w', rxFlGLWindow::OnMenuSimulation_s, m_pGLCanvas, FL_MENU_TOGGLE);
		m_pMenuBar->add("Simulation/SPS Turbulence",			'W', rxFlGLWindow::OnMenuSimulation_s, m_pGLCanvas, FL_MENU_TOGGLE);
		m_pMenuBar->add("Simulation/Vorticity Confinement  ",	'V', rxFlGLWindow::OnMenuSimulation_s, m_pGLCanvas, FL_MENU_TOGGLE);
		m_pMenuBar->add("Simulation/SPH Only",					0, rxFlGLWindow::OnMenuSimulation_s, m_pGLCanvas, FL_MENU_TOGGLE);	//�ǉ�			
		m_pMenuBar->add("Simulation/HeatTransfar",				0, rxFlGLWindow::OnMenuSimulation_s, m_pGLCanvas, FL_MENU_TOGGLE);	//�ǉ�	
		m_pMenuBar->add("Simulation/ShapeMatching",				0, rxFlGLWindow::OnMenuSimulation_s, m_pGLCanvas, FL_MENU_TOGGLE);	//�ǉ�
		m_pMenuBar->add("Simulation/IceStructure",				0, rxFlGLWindow::OnMenuSimulation_s, m_pGLCanvas, FL_MENU_TOGGLE | FL_MENU_DIVIDER);	//�ǉ�
		
		m_pMenuBar->add("Simulation/Particle Data Input",		"^i", rxFlGLWindow::OnMenuSimulation_s, m_pGLCanvas, FL_MENU_TOGGLE);
		m_pMenuBar->add("Simulation/Particle Data Output",		"^o", rxFlGLWindow::OnMenuSimulation_s, m_pGLCanvas, FL_MENU_TOGGLE);
		m_pMenuBar->add("Simulation/Mesh Saving",				"^m", rxFlGLWindow::OnMenuSimulation_s, m_pGLCanvas, FL_MENU_TOGGLE);
		m_pMenuBar->add("Simulation/Image Saving",				"^a", rxFlGLWindow::OnMenuSimulation_s, m_pGLCanvas, FL_MENU_TOGGLE);

		// Particle���j���[
		m_pMenuBar->add("Particle/Anisotoropic Kernel  ",	'a', rxFlGLWindow::OnMenuParticle_s, m_pGLCanvas, FL_MENU_TOGGLE);

		// Particle/Color���j���[
		m_pMenuBar->add("Particle/Color/Ramp",				0, rxFlGLWindow::OnMenuParticle_s, m_pGLCanvas, FL_MENU_RADIO);
		m_pMenuBar->add("Particle/Color/Constant",			0, rxFlGLWindow::OnMenuParticle_s, m_pGLCanvas, FL_MENU_RADIO);
		m_pMenuBar->add("Particle/Color/Density",			0, rxFlGLWindow::OnMenuParticle_s, m_pGLCanvas, FL_MENU_RADIO);
		m_pMenuBar->add("Particle/Color/Energy Spectrum",	0, rxFlGLWindow::OnMenuParticle_s, m_pGLCanvas, FL_MENU_RADIO);
		m_pMenuBar->add("Particle/Color/Pressure",			0, rxFlGLWindow::OnMenuParticle_s, m_pGLCanvas, FL_MENU_RADIO);
		m_pMenuBar->add("Particle/Color/Surface",			0, rxFlGLWindow::OnMenuParticle_s, m_pGLCanvas, FL_MENU_RADIO);
		m_pMenuBar->add("Particle/Color/Temperature",		0, rxFlGLWindow::OnMenuParticle_s, m_pGLCanvas, FL_MENU_RADIO);	//�ǉ�
		m_pMenuBar->add("Particle/Color/Ice_Cnct",	  FL_CTRL+FL_SHIFT+'a', rxFlGLWindow::OnMenuParticle_s, m_pGLCanvas, FL_MENU_RADIO);	//�ǉ�
		m_pMenuBar->add("Particle/Color/Ice_Calc",	  FL_CTRL+FL_SHIFT+'q', rxFlGLWindow::OnMenuParticle_s, m_pGLCanvas, FL_MENU_RADIO);	//�ǉ�
		m_pMenuBar->add("Particle/Color/FAST_PATH",	  FL_CTRL+FL_SHIFT+'p', rxFlGLWindow::OnMenuParticle_s, m_pGLCanvas, FL_MENU_RADIO);	//�ǉ�
		m_pMenuBar->add("Particle/Color/Edge",				0, rxFlGLWindow::OnMenuParticle_s, m_pGLCanvas, FL_MENU_RADIO);	//�ǉ�
		m_pMenuBar->add("Particle/Color/None",				0, rxFlGLWindow::OnMenuParticle_s, m_pGLCanvas, FL_MENU_RADIO);

		// Particle/Draw���j���[
		count = 0;
		while(RX_PARTICLE_DRAW[2*count] != "-1"){
			string label = "Particle/Draw/"+RX_PARTICLE_DRAW[2*count]+"  ";
			string shortcut = RX_PARTICLE_DRAW[2*count+1];
			m_pMenuBar->add(RX_TO_CHAR(label), RX_TO_CHAR(shortcut), rxFlGLWindow::OnMenuParticle_s, m_pGLCanvas, FL_MENU_RADIO);
			count++;
		}

		// Solid���j���[
		count = 0;
		while(RXS_STR[2*count] != "-1"){
			string label = "Solid/"+RXS_STR[2*count]+"  ";
			string shortcut = RXS_STR[2*count+1];
			m_pMenuBar->add(RX_TO_CHAR(label), RX_TO_CHAR(shortcut), rxFlGLWindow::OnMenuSolid_s, m_pGLCanvas, FL_MENU_TOGGLE);
			count++;
		}

		// Mesh���j���[
		count = 0;
		while(RX_TRIANGULATION_METHOD[2*count] != "-1"){
			string label = "Mesh/"+RX_TRIANGULATION_METHOD[2*count]+"  ";
			string shortcut = RX_TRIANGULATION_METHOD[2*count+1];
			m_pMenuBar->add(RX_TO_CHAR(label), RX_TO_CHAR(shortcut), rxFlGLWindow::OnMenuTriangulation_s, m_pGLCanvas, FL_MENU_RADIO);
			count++;
		}

		// Scene���j���[
		count = 0;
		while(!m_pGLCanvas->m_vSceneTitles[count].empty() && count < 12){
			string label = "Scene/"+m_pGLCanvas->m_vSceneTitles[count]+"  ";
			m_pMenuBar->add(RX_TO_CHAR(label), FL_F+count+1, rxFlGLWindow::OnMenuScene_s, m_pGLCanvas, FL_MENU_RADIO);
			count++;
		}
		
		// Step���j���[
		m_pMenuBar->add("Step/Step", ' ', OnMenuStep_s, this); 
		m_pMenuBar->add("Step/Animation  ", 's', OnMenuStep_s, this); 

		// Window���j���[
		m_pMenuBar->add("Window/FullScreen (Window)  ", FL_CTRL+'f', OnMenuWindow_s, this, FL_MENU_TOGGLE); 
		m_pMenuBar->add("Window/FullScreen (GLCanvas)  ", FL_ALT+FL_Enter, OnMenuWindow_s, this, FL_MENU_TOGGLE); 
		count = 0;
		while(RX_CANVAS_SIZE_STR[2*count] != "-1"){
			string label = "Window/Canvas Size/"+RX_CANVAS_SIZE_STR[2*count]+"  ";
			string shortcut = RX_CANVAS_SIZE_STR[2*count+1];
			m_pMenuBar->add(RX_TO_CHAR(label), RX_TO_CHAR(shortcut), OnMenuWindow_s, this);
			count++;
		}

		// Help���j���[
		m_pMenuBar->add("Help/Version  ", 0, OnMenuHelpVersion_s, this); 
	}
	{
		// ����GUI�̈�
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

		// Start/Stop�{�^��
		button = new Fl_Button(xs, ys, ws, hs, "Start/Stop");
		button->callback(OnButtonStart_s, this);
		button->down_box(boxtype);
		button->clear_visible_focus();

		// Apply�{�^��
		button = new Fl_Button(xs, ys+1*(hs+5), ws, hs, "Apply");
		button->callback(OnButtonApply_s, this);
		button->down_box(boxtype);
		button->clear_visible_focus();

		// Clear�{�^��
		button = new Fl_Button(xs, ys+2*(hs+5), ws, hs, "Clear");
		button->down_box(boxtype);
		button->clear_visible_focus();

		//
		// �V�~�����[�V�����p�����[�^
		//
		xs += ws+7;
		ys = h()-(hs_stat+ver_margin+hs_para)+5;
		ws = 390;
		hs = hs_para-6;

		Fl_Group* sg = new Fl_Group(xs, ys, ws, hs, "Turbulence");
		{
			sg->box(FL_DOWN_BOX);
			sg->labelsize(12);
			sg->align(FL_ALIGN_TOP | FL_ALIGN_INSIDE);

			int dx = 70;
			ys += 30;

			m_pSpinCoefEt = new Fl_Spinner(xs+dx, ys, 87, 25, "Ecoef ");
			m_pSpinCoefEt->type(1);
			m_pSpinCoefEt->callback(OnSpinSimulation_s, this);
			m_pSpinCoefEt->minimum(0);
			m_pSpinCoefEt->maximum(0.005);
			m_pSpinCoefEt->step(1e-005);
			m_pSpinCoefEt->value(g_fCoefEt);
			m_pSpinCoefEt->clear_visible_focus();

			m_pSpinMaxEt = new Fl_Spinner(xs+dx, ys+30, 87, 25, "Emax ");
			m_pSpinMaxEt->type(1);
			m_pSpinMaxEt->callback(OnSpinSimulation_s, this);
			m_pSpinMaxEt->minimum(0);
			m_pSpinMaxEt->maximum(0.2);
			m_pSpinMaxEt->step(0.0002);
			m_pSpinMaxEt->value(g_fMaxEt);
			m_pSpinMaxEt->clear_visible_focus();
		
			m_pSpinWaveletScale = new Fl_Spinner(xs+dx, ys+60, 87, 25, "Scale ");
			m_pSpinWaveletScale->type(1);
			m_pSpinWaveletScale->callback(OnSpinSimulation_s, this);
			m_pSpinWaveletScale->minimum(0);
			m_pSpinWaveletScale->maximum(10);
			m_pSpinWaveletScale->step(0.1);
			m_pSpinWaveletScale->value(g_fWaveletScale);
			m_pSpinWaveletScale->clear_visible_focus();

			dx = 280;

			m_pSpinSpEt = new Fl_Spinner(xs+dx, ys, 87, 25, "SPS Scale ");
			m_pSpinSpEt->type(1);
			m_pSpinSpEt->callback(OnSpinSimulation_s, this);
			m_pSpinSpEt->minimum(0);
			m_pSpinSpEt->maximum(200);
			m_pSpinSpEt->step(0.2);
			m_pSpinSpEt->value(g_fCoefTurb);
			m_pSpinSpEt->clear_visible_focus();
		
			m_pSpinSpEtMesh = new Fl_Spinner(xs+dx, ys+30, 87, 25, "SPS Mesh Scale ");
			m_pSpinSpEtMesh->type(1);
			m_pSpinSpEtMesh->callback(OnSpinSimulation_s, this);
			m_pSpinSpEtMesh->minimum(0);
			m_pSpinSpEtMesh->maximum(20);
			m_pSpinSpEtMesh->step(0.02);
			m_pSpinSpEtMesh->value(g_fCoefTurbForMesh);
			m_pSpinSpEtMesh->clear_visible_focus();
		
			m_pSpinEtCri = new Fl_Spinner(xs+dx, ys+60, 87, 25, "SPS Ecri ");
			m_pSpinEtCri->type(1);
			m_pSpinEtCri->callback(OnSpinSimulation_s, this);
			m_pSpinEtCri->minimum(0);
			m_pSpinEtCri->maximum(1);
			m_pSpinEtCri->step(0.001);
			m_pSpinEtCri->value(g_fEtCri);
			m_pSpinEtCri->clear_visible_focus();
	
			sg->resizable(NULL);
			sg->end();
		}

		xs += ws+5;
		ys = h()-(hs_stat+ver_margin+hs_para)+5;
		ws = 280;

		sg = new Fl_Group(xs, ys, ws, hs, "Draw");
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
		//�ǉ��F�F�M�v�Z�p�����[�^��GUI
		xs += ws+5;
		ys = h()-(hs_stat+ver_margin+hs_para)+5;
		ws = 280;

		sg = new Fl_Group(xs, ys, ws, hs, "TEMP");
		{
			sg->box(FL_DOWN_BOX);
			sg->labelsize(12);
			sg->align(FL_ALIGN_TOP | FL_ALIGN_INSIDE);

			int dx = 125;
			ys += 30;

			m_pSliderVScale = new Fl_Value_Slider(xs+dx, ys, 145, 25, "AirTemp");
			m_pSliderVScale->type(1);
			m_pSliderVScale->callback(OnSliderTemp_s, this);
			m_pSliderVScale->minimum(-1000);
			m_pSliderVScale->maximum(1000);
			m_pSliderVScale->step(10);
//			m_pSliderVScale->value(m_pGLCanvas->m_ht->getAirTemp());	//�������O�ɎQ�Ƃ��Ă��܂��G���[
			m_pSliderVScale->value(0);
			m_pSliderVScale->align(Fl_Align(FL_ALIGN_LEFT));
			m_pSliderVScale->clear_visible_focus();

			sg->resizable(NULL);
			sg->end();
		}

		//�ǉ��F�F�N���X�^�p�����[�^��GUI
		xs += ws+5;
		ys = h()-(hs_stat+ver_margin+hs_para)+5;
		ws = 280;

		sg = new Fl_Group(xs, ys, ws, hs, "ICE");
		{
			sg->box(FL_DOWN_BOX);
			sg->labelsize(12);
			sg->align(FL_ALIGN_TOP | FL_ALIGN_INSIDE);

			int dx = 125;
			ys += 30;

			m_pSliderVScale = new Fl_Value_Slider(xs+dx, ys, 145, 25, "AirTemp");
			m_pSliderVScale->type(1);
			m_pSliderVScale->callback(OnSliderTemp_s, this);
			m_pSliderVScale->minimum(-1000);
			m_pSliderVScale->maximum(1000);
			m_pSliderVScale->step(10);
//			m_pSliderVScale->value(m_pGLCanvas->m_ht->getAirTemp());	//�������O�ɎQ�Ƃ��Ă��܂��G���[
			m_pSliderVScale->value(0);
			m_pSliderVScale->align(Fl_Align(FL_ALIGN_LEFT));
			m_pSliderVScale->clear_visible_focus();

			sg->resizable(NULL);
			sg->end();
		}

		g->resizable(NULL);
		g->end();
	}
	{
		// �X�e�[�^�X�o�[(Fl_Box)
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

//! �f�X�g���N�^
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
 * �g�O���t���j���[���ڂ̍X�V
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
	SetMenuItemState(m_pMenuBar, "Simulation/SPH Only",					(m_pGLCanvas->m_bsSimuSetting.at(ID_HEAT)));		//�ǉ�	
	SetMenuItemState(m_pMenuBar, "Simulation/HeatTransfar",				(m_pGLCanvas->m_bsSimuSetting.at(ID_HEAT)));		//�ǉ�	
	SetMenuItemState(m_pMenuBar, "Simulation/ShapeMatching",			(m_pGLCanvas->m_bsSimuSetting.at(ID_SM)));			//�ǉ�	
	SetMenuItemState(m_pMenuBar, "Simulation/IceStructure",				(m_pGLCanvas->m_bsSimuSetting.at(ID_ICE)));			//�ǉ�	

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
	SetMenuItemState(m_pMenuBar, "Particle/Color/Temperature",		(m_pGLCanvas->m_iColorType == rxParticleSystemBase::RX_TEMP));			//�ǉ�
	SetMenuItemState(m_pMenuBar, "Particle/Color/Ice_Connect",		(m_pGLCanvas->m_iColorType == rxParticleSystemBase::RX_ICE_CONNECT));	//�ǉ�
	SetMenuItemState(m_pMenuBar, "Particle/Color/Ice_Calc",			(m_pGLCanvas->m_iColorType == rxParticleSystemBase::RX_ICE_CALC));		//�ǉ�
	SetMenuItemState(m_pMenuBar, "Particle/Color/Edge",				(m_pGLCanvas->m_iColorType == rxParticleSystemBase::RX_EDGE));			//�ǉ�
	SetMenuItemState(m_pMenuBar, "Particle/Color/ICE_FAST_PATH",	(m_pGLCanvas->m_iColorType == rxParticleSystemBase::RX_ICE_FAST_PATH));	//�ǉ�
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

	// Window���j���[
	SetMenuItemState(m_pMenuBar, "Window/FullScreen  ", m_bFullScreen, 1);

	// Refraction�`�F�b�N�{�^��
	m_pCheckRefraction->value((m_pGLCanvas->m_iDraw & RXD_REFRAC ? 1 : 0));

	// Mesh�`�F�b�N�{�^��
	m_pCheckMesh->value((m_pGLCanvas->m_iDraw & RXD_MESH ? 1 : 0));
}

/*!
 * Fl_Button�̃R�[���o�b�N�֐� - Start/Stop�{�^��
 * @param[in] widget �E�B�W�b�g�̐e�N���X�I�u�W�F�N�g
 * @param[in] x ���[�U��`�ϐ�
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
 * Fl_Button�̃R�[���o�b�N�֐� - Apply�{�^��
 * @param[in] widget �E�B�W�b�g�̐e�N���X�I�u�W�F�N�g
 * @param[in] x ���[�U��`�ϐ�
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

/*!
 * Fl_Spinner�̃R�[���o�b�N�֐� - Simulation
 * @param[in] widget �E�B�W�b�g�̐e�N���X�I�u�W�F�N�g
 * @param[in] x ���[�U��`�ϐ�
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
 * Fl_Value_Slider�̃R�[���o�b�N�֐� - Draw
 * @param[in] widget �E�B�W�b�g�̐e�N���X�I�u�W�F�N�g
 * @param[in] x ���[�U��`�ϐ�
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
 * Fl_Value_Slider�̃R�[���o�b�N�֐� - TEMP
 * @param[in] widget �E�B�W�b�g�̐e�N���X�I�u�W�F�N�g
 * @param[in] x ���[�U��`�ϐ�
 */
void rxFlWindow::OnSliderTemp_s(Fl_Widget *widget, void* x)
{
	((rxFlWindow*)x)->OnSliderTemp(widget);
}
void rxFlWindow::OnSliderTemp(Fl_Widget *widget)
{
	Fl_Value_Slider *slider = (Fl_Value_Slider*)widget;
	string label = slider->label();
	double val = slider->value();

	if(label.find("AirTemp") != string::npos){
		m_pGLCanvas->m_ht->setAirTemp(val);
	}
	//else if(label.find("Mesh Threshold") != string::npos){
	//	m_pGLCanvas->m_fMeshThr = val;
	//}
}


/*!
 * Fl_Check_Button�̃R�[���o�b�N�֐� - Draw
 * @param[in] widget �E�B�W�b�g�̐e�N���X�I�u�W�F�N�g
 * @param[in] x ���[�U��`�ϐ�
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

	UpdateMenuState();
}


/*!
 * rxFlDndBox�̃R�[���o�b�N�֐�
 * @param[in] widget �E�B�W�b�g�̐e�N���X�I�u�W�F�N�g
 * @param[in] x ���[�U��`�ϐ�
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

		// �e�L�X�g��\n�ŕ���
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
 * File���j���[�R�[���o�b�N�֐�
 * @param[in] widget �E�B�W�b�g�̐e�N���X�I�u�W�F�N�g
 * @param[in] x ���[�U��`�ϐ�
 */
void rxFlWindow::OnMenuFile_s(Fl_Widget *widget, void* x)
{
	Fl_Menu_Bar *menubar = (Fl_Menu_Bar*)widget;
	char picked[80];
	menubar->item_pathname(picked, sizeof(picked)-1);	// ���j���[��

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
 * Draw���j���[�̃R�[���o�b�N�֐�
 * @param[in] widget �E�B�W�b�g�̐e�N���X�I�u�W�F�N�g
 * @param[in] x ���[�U��`�ϐ�
 */
void rxFlWindow::OnMenuStep_s(Fl_Widget *widget, void* x)
{
	Fl_Menu_Bar *menubar = (Fl_Menu_Bar*)widget;
	char picked[80];
	menubar->item_pathname(picked, sizeof(picked)-1);	// ���j���[��

	string label = picked;
	string menu_name = "Step/";
	label = label.substr(menu_name.size(), string::npos);
	
	((rxFlWindow*)x)->OnMenuStep(label);
}
void rxFlWindow::OnMenuStep(string label)
{
	if(label.find("Step") == 0){// �A�j���[�V����1�X�e�b�v�������s
		m_pGLCanvas->Idle();
		m_pGLCanvas->m_bFall = !m_pGLCanvas->m_bFall;
	}
	else if(label.find("Animation") == 0){	// �A�j���[�V����ON/OFF
		m_pGLCanvas->SwitchIdle(-1);
	}
}


/*!
 * Window���j���[�̃R�[���o�b�N�֐�
 * @param[in] widget �E�B�W�b�g�̐e�N���X�I�u�W�F�N�g
 * @param[in] x ���[�U��`�ϐ�
 */
void rxFlWindow::OnMenuWindow_s(Fl_Widget *widget, void* x)
{
	Fl_Menu_Bar *menubar = (Fl_Menu_Bar*)widget;
	char picked[80];
	menubar->item_pathname(picked, sizeof(picked)-1);	// ���j���[��

	string label = picked;
	string menu_name = "Window/";
	label = label.substr(menu_name.size(), string::npos);
	
	((rxFlWindow*)x)->OnMenuWindow(label);
}
void rxFlWindow::OnMenuWindow(string label)
{
	if(label.find("FullScreen (Window)") == 0){	// �t���X�N���[��ON/OFF
		SwitchFullScreen();
	}
	else if(label.find("FullScreen (GLCanvas)") == 0){	// �t���X�N���[��ON/OFF
		m_pGLCanvas->SwitchFullScreen(0);
	}
	else if(label.find("Canvas Size/") == 0){	// �L�����o�X�T�C�Y
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
 * ���j���[:Help -> Version�̃R�[���o�b�N�֐�
 * @param[in] widget �E�B�W�b�g�̐e�N���X�I�u�W�F�N�g
 * @param[in] x ���[�U��`�ϐ�
 */
void rxFlWindow::OnMenuHelpVersion_s(Fl_Widget *widget, void* x)
{
	fl_message("OpenGL Application by FLTK\n  version 1.0");
}

/*!
 * �X�e�[�^�X�o�[�ɕ������ݒ�
 * @param[in] label �\��������
 */
void rxFlWindow::SetStatusLabel(const string &label)
{
	if(m_pStatusLabel) delete [] m_pStatusLabel;
	m_pStatusLabel = RX_TO_CHAR(label);
	m_pBoxStatus->label(m_pStatusLabel);
}

/*!
 * �t�@�C���ǂݍ���
 * @param[in] fn �t�@�C���p�X
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

	// �ǂݍ��񂾃t�@�C�������i�[
	m_strFullPath = fn;
	m_strFileName = GetFileName(fn);

	// �t�@�C�������X�e�[�^�X�o�[�ɕ\��
	SetStatusLabel(m_strFileName);
}

/*!
 * �t�@�C����������
 * @param[in] fn �t�@�C���p�X
 */
void rxFlWindow::Save(const string &fn)
{
	string ext = GetExtension(fn);

	if(ext == "obj" || ext == "dxf" || ext == "wrl" || ext == "3ds" || ext == "stl" || ext == "ply"){
		m_pGLCanvas->SaveFile(fn);

		// �ǂݍ��񂾃t�@�C�������i�[
		m_strFullPath = fn;
		m_strFileName = GetFileName(fn);

		// �t�@�C�������X�e�[�^�X�o�[�ɕ\��
		SetStatusLabel(m_strFileName);
	}
	else if(ext == "bmp" || ext == "png"){
		m_pGLCanvas->SaveDisplay(fn);
	}
}



/*!
 * �ݒ�t�@�C���ǂݍ���
 * @param[in] fn �ݒ�t�@�C����(�g���q����)
 */
void rxFlWindow::ReadConfig(const string &fn)
{
	// �A�v���P�[�V�����Ǘ��ݒ�t�@�C��
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
 * �ݒ�t�@�C����������
 * @param[in] fn �ݒ�t�@�C����(�g���q����)
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
 * �t���X�N���[��/�E�B���h�E�\���̐؂�ւ�
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
 * �t���X�N���[��/�E�B���h�E�\���̏�Ԏ擾
 */
bool rxFlWindow::IsFullScreen(void) const
{
	return m_bFullScreen;
}

/*!
 * �C�x���g�n���h��
 * @param[in] ev �C�x���gID
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

	case FL_PUSH:		// �}�E�X�{�^���_�E��
		m_pGLCanvas->Mouse(Fl::event_button(), 1, Fl::event_x(), Fl::event_y());
		break;
	case FL_RELEASE:	// �}�E�X�{�^���A�b�v
		m_pGLCanvas->Mouse(Fl::event_button(), 0, Fl::event_x(), Fl::event_y());
		break;
	case FL_DRAG:		// �}�E�X�h���b�O
		m_pGLCanvas->Motion(Fl::event_x(), Fl::event_y());
		break;
	case FL_MOVE:		// �}�E�X�ړ�
		m_pGLCanvas->PassiveMotion(Fl::event_x(), Fl::event_y());
		break;

	case FL_KEYBOARD:	// �L�[�_�E��
		m_pGLCanvas->Keyboard(Fl::event_key(), Fl::event_x(), Fl::event_y());
		UpdateMenuState();
		break;

	case FL_SHORTCUT:	// �O���[�o���V���[�g�J�b�g
		break;

	default:
		break;
	}

	return Fl_Window::handle(ev);
}