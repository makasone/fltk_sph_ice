/*!
  @file main.cpp
	
  @brief FLTK��OpenGL
 
  @author Makoto Fujisawa
  @date   2011-08
*/

#ifdef _DEBUG
#pragma comment(lib, "fltkd.lib")
//#pragma comment(lib, "fltkformsd.lib")
#pragma comment(lib, "fltkgld.lib")
#pragma comment(lib, "fltkimagesd.lib")
#pragma comment(lib, "fltkjpegd.lib")
#pragma comment(lib, "fltkpngd.lib")
#pragma comment(lib, "fltkzlibd.lib")
#else
#pragma comment(lib, "fltk.lib")
//#pragma comment(lib, "fltkforms.lib")
#pragma comment(lib, "fltkgl.lib")
#pragma comment(lib, "fltkimages.lib")
#pragma comment(lib, "fltkjpeg.lib")
#pragma comment(lib, "fltkpng.lib")
#pragma comment(lib, "fltkzlib.lib")
#endif


#pragma comment(lib, "glew32.lib")
#pragma comment(lib, "cudart.lib")

#pragma comment (lib, "freetype242.lib")
#pragma comment (lib, "ftgl.lib")

#ifdef _DEBUG
#pragma comment(lib, "rx_modeld.lib")
//#pragma comment (lib, "ode_doubled.lib")
//#pragma comment (lib, "libjpegd.lib")
//#pragma comment (lib, "libpngd.lib")
//#pragma comment (lib, "zlibd.lib")
#else
#pragma comment(lib, "rx_model.lib")
//#pragma comment (lib, "ode_double.lib")
//#pragma comment (lib, "libjpeg.lib")
//#pragma comment (lib, "libpng.lib")
//#pragma comment (lib, "zlib.lib")
#endif



//#pragma comment (lib, "opencv_core231.lib")
//#pragma comment (lib, "opencv_imgproc231.lib")
//#pragma comment (lib, "opencv_highgui231.lib")
//#pragma comment (lib, "opencv_video231.lib")



// �R�}���h�v�����v�g���o�������Ȃ��ꍇ�͂������R�����g�A�E�g
//#pragma comment(linker, "/subsystem:\"windows\" /entry:\"mainCRTStartup\"")


//-----------------------------------------------------------------------------
// �C���N���[�h�t�@�C��
//-----------------------------------------------------------------------------
#include "rx_fltk_window.h"

// CUDA
#include "rx_cu_funcs.cuh"

//-----------------------------------------------------------------------------
// ���C���֐�
//-----------------------------------------------------------------------------
/*!
 * ���C�����[�`��
 * @param[in] argc �R�}���h���C�������̐�
 * @param[in] argv �R�}���h���C������
 */
int main(int argc, char *argv[])
{	
	// �R�}���h���C������
	if(argc >= 2){
		for(int i = 1; i < argc; ++i){
			string fn = argv[i];
			g_vDefaultFiles.push_back(fn);
		}
	}

	glutInit(&argc, argv);
	CuInit(argc, argv);

	Fl::visual(FL_DOUBLE | FL_INDEX);
	Fl::get_system_colors();
	fl_register_images();
	Fl::scheme("gtk+");

	rxFlWindow win(480, 480, "opengl application");

	return Fl::run();
}


