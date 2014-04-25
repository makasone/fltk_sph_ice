/*!
  @file rx_fltk_widgets.h
	
  @brief FLTK�ɂ��J�X�^���E�B�W�b�g
 
  @author Makoto Fujisawa 
  @date   2011-08
*/

#ifndef _RX_FLTK_WIDGETS_H_
#define _RX_FLTK_WIDGETS_H_


//-----------------------------------------------------------------------------
// �C���N���[�h�t�@�C��
//-----------------------------------------------------------------------------
#include <sstream>
#include <string>

#include <FL/Fl.H>
#include <FL/Fl_Window.H>
#include <FL/Fl_Input.H>
#include <FL/Fl_Output.H>
#include <FL/Fl_Slider.H>
#include <FL/Fl_Box.H>
#include <FL/Fl_Shared_Image.H> // �摜

#include <GL/glew.h>
#include <GL/glut.h>

#ifdef WIN32
#include <windows.h>
#include <commdlg.h>
#endif



using namespace std;




//-----------------------------------------------------------------------------
// �����񏈗��֐�
//-----------------------------------------------------------------------------



/*!
 * string����char[]�ɕϊ�(\0�I���)
 * @param[in] s string������
 * @return char������
 */
inline char* RX_TO_CHAR(const string &s)
{
	if(s.empty()) return 0;

	int n = (int)s.size();
	char* c = new char[n+1];
	for(int i = 0; i < n; ++i) c[i] = s[i];
	c[n] = '\0';
	return c;
}

/*!
 * string����char[]�ɕϊ�(\0�I���)
 * @param[in] s string������
 * @param[out] c char������
 */
inline void RX_TO_CHAR(const string &s, char c[])
{
	if(s.empty()) return;

	int n = (int)s.size();
	for(int i = 0; i < n; ++i) c[i] = s[i];
	c[n] = '\0';
}



/*!
 * �p�X����t�@�C��������菜�����p�X�𒊏o
 * @param[in] path �p�X
 * @return �t�H���_�p�X
 */
inline string GetFolderPath(const string &path)
{
	size_t pos1;
 
	pos1 = path.rfind('\\');
	if(pos1 != string::npos){
		return path.substr(0, pos1+1);
		
	}
 
	pos1 = path.rfind('/');
	if(pos1 != string::npos){
		return path.substr(0, pos1+1);
	}
 
	return "";
}




//-----------------------------------------------------------------------------
// FLTK�֘A�֐�
//-----------------------------------------------------------------------------
/*!
 * ���j���[�A�C�e��(FL_MENU_TOGGLE)�̏�Ԃ�ύX
 * @param[in] menubar ���j���[�o�[�I�u�W�F�N�g
 * @param[in] name ���j���[��("File/Open"�Ȃ�)
 * @param[in] state �g�O��ON/OFF
 * @param[in] enable �L��/����
 * @return ���j���[�����݂��Ȃ����-1��Ԃ�
 */
inline int SetMenuItemState(Fl_Menu_Bar *menubar, string name, int state, int enable = true)
{
	Fl_Menu_Item *m = (Fl_Menu_Item*)menubar->find_item(RX_TO_CHAR(name));
	if(!m) return -1;

	if(enable){
		m->activate();
	}
	else{
		m->deactivate();
	}

	if(state){
		m->set();
	}
	else{
		m->clear();
	}
	return(0);
}



//-----------------------------------------------------------------------------
//! rxFlDndBox�N���X - �h���b�O�A���h�h���b�v�p�{�b�N�X
//-----------------------------------------------------------------------------
class rxFlDndBox : public Fl_Box
{
	// MARK:rxFlDndBox
protected:
	// �����o�ϐ�
	int m_iEvent;	//!< �C�x���gID

	string m_strEventText;
	int m_iEventTextLen;

public:
	//! �R���X�g���N�^
	rxFlDndBox(int x, int y, int w, int h, const char *l = 0)
		 : Fl_Box(x, y, w, h, l), m_iEvent(FL_NO_EVENT), m_iEventTextLen(0)
	{
		labeltype(FL_NO_LABEL);
		box(FL_NO_BOX);
		clear_visible_focus();
	}

	//! �f�X�g���N�^
	virtual ~rxFlDndBox()
	{
	}

public:
	// �����o�֐�
	static void CallbackS(void *v)
	{
		rxFlDndBox *w = (rxFlDndBox*)v;
		w->do_callback();
	}

	int Event()
	{
		return m_iEvent;
	}

	string EventText()
	{
		return m_strEventText;
	}

	int EventTextLength()
	{
		return m_iEventTextLen;
	}

	int handle(int e)
	{
		switch(e){
			case FL_DND_ENTER:
				cout << "rxFlDndBox::FL_DND_ENTER" << endl;
				m_iEvent = e;
				return 1;
			case FL_DND_RELEASE:
				cout << "rxFlDndBox::FL_DND_RELEASE" << endl;
				m_iEvent = e;
				return 1;
			case FL_DND_LEAVE:
				cout << "rxFlDndBox::FL_DND_LEAVE" << endl;
				m_iEvent = e;
				return 1;
			case FL_DND_DRAG:
				cout << "rxFlDndBox::FL_DND_DRAG" << endl;
				m_iEvent = e;
				return 1;


			case FL_PASTE:
				cout << "rxFlDndBox::FL_PASTE" << endl;
				m_iEvent = e;

				m_iEventTextLen = Fl::event_length();
				m_strEventText = Fl::event_text();

				if(callback() && ((when() & FL_WHEN_RELEASE) || (when() & FL_WHEN_CHANGED)))
					Fl::add_timeout(0.0, rxFlDndBox::CallbackS, (void*)this);
				return 1;
		}

		return Fl_Box::handle(e);
	}
};



//-----------------------------------------------------------------------------
// FLTK�ł̉摜�ǂݍ���
//-----------------------------------------------------------------------------
inline unsigned char* ReadImageByFLTK(const string &fn, int &w, int &h, int &c)
{
	char* cfn = RX_TO_CHAR(fn);
	Fl_Shared_Image *simg = Fl_Shared_Image::get(cfn);

	if(simg->count() != 1){
		return 0;
	}

	w = simg->w();
	h = simg->h();
	c = simg->d();
	unsigned char *dat = new unsigned char[w*h*c];

	// �t�@�C���o��
	const char *buf = simg->data()[0];
	unsigned char r, g, b, a;
	for(int j = 0; j < h; ++j){
		for(int i = 0; i < w; ++i){
			long idx = j*w*c+i*c;
			r = g = b = a = 0;
			switch(c){
			case 1:
				r = g = b = *(buf+idx);
				break;

			case 3:
				r = *(buf+idx+0);
				g = *(buf+idx+1);
				b = *(buf+idx+2);
				break;

			case 4:
				r = *(buf+idx+0);
				g = *(buf+idx+1);
				b = *(buf+idx+2);
				a = *(buf+idx+3);
				break;
							
			default:
				break;
			}

			dat[idx+0] = r;
			dat[idx+1] = g;
			dat[idx+2] = b;
		}
	}

	return dat;
}

/*!
 * OpenGL�e�N�X�`���o�^
 * @param[in] fn �t�@�C����
 * @param[inout] tex_name �e�N�X�`����(0�Ȃ�V���ɐ���)
 * @param[in] mipmap �~�b�v�}�b�v�g�p�t���O
 * @param[in] compress �e�N�X�`�����k�g�p�t���O
 */
static int LoadGLTexture(const string &fn, GLuint &tex_name, bool mipmap, bool compress)
{
	// �摜�ǂݍ���
	unsigned char *buf = 0;
	int w, h, c;
	buf = ReadImageByFLTK(fn, w, h, c);
	if(buf == 0){
		return -1;
	}

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
 
		glTexImage2D(GL_TEXTURE_2D, 0, iformat, w, h, 0, format, GL_UNSIGNED_BYTE, buf);
 
		if(mipmap){
			glGenerateMipmapEXT(GL_TEXTURE_2D);
		}
	}
	else{
		glBindTexture(GL_TEXTURE_2D, tex_name);
		//glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, w, h, format, GL_UNSIGNED_BYTE, pimg);
		glTexImage2D(GL_TEXTURE_2D, 0, iformat, w, h, 0, format, GL_UNSIGNED_BYTE, buf);

		if(mipmap){
			glGenerateMipmapEXT(GL_TEXTURE_2D);
		}
	}

	delete [] buf;
 
	glBindTexture(GL_TEXTURE_2D, 0);

	return 1;
}




//-----------------------------------------------------------------------------
// MARK:�t�@�C���_�C�A���O 
//-----------------------------------------------------------------------------
/*!
 * "|"�ŋ�؂�ꂽ�t�@�C���t�B���^��������Ɗg���q�w�蕔�ɕ���
 *  [��] "Text Files|*.txt;*.log;*.dat|All Files|*.*"
 *       -> descs : "Text Files", "All Files"
 *       -> exts  : "*.txt;*.log;*.dat", "*.*"
 * @param[in] filter "|"�ŋ�؂�ꂽ�t�@�C���t�B���^������
 * @param[out] descs ������
 * @param[out] exts  �g���q�w�蕔
 * @return �t�B���^��( = descs.size() = exts.size() )
 */
inline int ParseFilterString(const string &filter, vector<string> &descs, vector<string> &exts)
{
	int nfilter = 0;
	size_t pos0 = 0, pos1 = 0;
	do{
		// �t�B���^�����̒��o
		pos1 = filter.find('|', pos0);
		if(pos0 < pos1){
			descs.push_back(filter.substr(pos0, pos1-pos0));
			pos0 = pos1+1;

			// �t�B���^�g���q�̒��o
			pos1 = filter.find('|', pos0);
			if(pos0 < pos1){
				exts.push_back(filter.substr(pos0, pos1-pos0));
				nfilter++;
			}
		}

		pos0 = pos1+1;
	}while(pos1 != string::npos);

	return nfilter;
}


/*!
 * fltk��File_Chooser�p�̃t�B���^����
 *  [��] "Text Files (*.{txt,log,dat})\tAll Files (*)"
 * @param[in] descs ������
 * @param[in] exts  �g���q�w�蕔
 * @return �t�B���^������
 */
inline string GenFlFileFilter(const vector<string> &descs, const vector<string> &exts)
{
	int nfilter = (int)descs.size();
	string fl_filter;
	for(int i = 0; i < nfilter; ++i){
		// ����
		fl_filter += descs[i];

		// �g���q
		fl_filter += " (";

		// *.txt;*.log;*.dat �̂悤�ȕ����񂩂�Ctxt,log,dat�𒊏o
		vector<string> es;
		size_t epos0 = 0, epos1 = 0;
		do{
			epos1 = exts[i].find(';', epos0);
			if(epos0 < epos1){
				// "*."��";"�̊Ԃ̕�����𒊏o
				es.push_back(exts[i].substr(epos0+2, epos1-epos0-2));
			}
			epos0 = epos1+1;
		}while(epos1 != string::npos);

		// fltk�̃t�B���^�`���ɕϊ� *.{txt,log,dat}
		if((int)es.size() > 2){
			// �����g���q�w�莞
			fl_filter += "*.{";
			for(int j = 0; j < (int)es.size(); ++j){
				fl_filter += es[j];
				if(j != (int)es.size()-1) fl_filter += ",";
			}
			fl_filter += "}";
		}
		else if(!es.empty() && es[0] != "*"){
			// �P��g���q
			fl_filter += "*."+es[0];
		}
		else{
			// �C�ӊg���q(*)
			fl_filter += "*";
		}

		fl_filter += ")";
		if(i < nfilter-1) fl_filter += "\t";
	}

	return fl_filter;
}

/*!
 * Win32��GetOpenFileName�p�̃t�B���^����(�k��������؂�)
 *  [��] "Text Files\0*.txt;*.log;*.dat\0All Files\0*.*\0\0"
 * @param[in] descs ������
 * @param[in] exts  �g���q�w�蕔
 * @param[out] cfilter �t�B���^������
 * @param[out] n �t�B���^������ő�T�C�Y
 */
inline int GenWin32FileFilter(const vector<string> &descs, const vector<string> &exts, char cfilter[], int n)
{
	int nfilter = (int)descs.size();
	int c = 0, j = 0, k = 0;
	while(c < n){
		if(k%2 == 0){
			if(k != 0 && j == 0){	// ��؂�̃k������
				cfilter[c++] = NULL;
			}

			// ����
			cfilter[c++] = descs[k/2][j++];
			if(j >= (int)descs[k/2].size()){
				j = 0;
				k++;
			}
		}
		else{
			if(k != 0 && j == 0){	// ��؂�̃k������
				cfilter[c++] = NULL;
			}

			// �g���q
			cfilter[c++] = exts[k/2][j++];
			if(j >= (int)exts[k/2].size()){
				j = 0;
				k++;
			}
		}

		if(k >= nfilter*2){	// �Ō�̃k������x2
			cfilter[c++] = NULL;
			cfilter[c++] = NULL;
			break;
		}
	}
	return 1;
}

/*!
 * �t�@�C���I���_�C�A���O�̕\��
 *  Windows�ł̓G�N�X�v���[���`���C���̑��ł�fltk��File_Chooser
 * @param[out] fns �I�����ꂽ�t�@�C��
 * @param[in] title �_�C�A���O�^�C�g��
 * @param[in] filter �t�@�C���t�B���^("|"��؂�)
 * @param[in] multi �����I���̉�
 */
inline int ShowFileDialog(vector<string> &fns, 
						  const string &title, const string &filter, bool multi = false)
{
	vector<string> descs, exts;
	int nfilter = ParseFilterString(filter, descs, exts);

#ifdef WIN32

	OPENFILENAME ofn;
	memset((void*)&ofn, 0, sizeof(OPENFILENAME));

	// �t�@�C���I���_�C�A���O�̐ݒ�
	ofn.lStructSize  = sizeof(OPENFILENAME);
	ofn.Flags       |= OFN_NOVALIDATE;          // �����ȕ������������t�@�C������L���ɂ���(/����͂��܂�t�@�C���p�X��L����)
	ofn.Flags       |= OFN_HIDEREADONLY;        // �ǂݎ���p�`�F�b�N�{�b�N�X���B��
	ofn.Flags       |= OFN_EXPLORER;            // �V�����G�N�X�v���[���E�B���h�E���g�p
	ofn.Flags       |= OFN_ENABLESIZING;        // �_�C�A���O���T�C�Y��
	ofn.Flags       |= OFN_NOCHANGEDIR;         // �ߋ��ɊJ�����t�H���_���f�t�H���g�ɂ���
	if(multi) ofn.Flags |= OFN_ALLOWMULTISELECT;	// �����I��
	ofn.nMaxFile     = 4096-1;
	ofn.lpstrFile    = new char[4096];
	ofn.lpstrFile[0] = 0;
	ofn.hwndOwner    = GetForegroundWindow();
	ofn.lpstrTitle   = title.c_str();

	// �t�@�C���t�B���^
	char cfilter[1024];
	GenWin32FileFilter(descs, exts, cfilter, 1024);
	ofn.lpstrFilter  = cfilter;

	// �t�@�C���I���_�C�A���O���J��
	int err = GetOpenFileName(&ofn);
	if(err == 0){
		err = CommDlgExtendedError();
		if(err == 0) return 0;
		fprintf(stderr, "CommDlgExtendedError() code = %d", err);
		return 0;
	}

	// �����t�@�C���I������ofn.lpstrFile�Ƀf�B���N�g�����ƃt�@�C����(����)���k�������ŋ�؂��Ċi�[����Ă���D
	// �Ō�̃t�@�C�����̌�ɂ�2�̘A�������k���������i�[����Ă���D
	string tmp;
	vector<string> fns0;
	int null_num = 0;
	for(int i = 0; i < (int)ofn.nMaxFile; ++i){
		if(ofn.lpstrFile[i] == NULL){	// �k�������ɂ���؂�
			if(!null_num){
				fns0.push_back(tmp);
				tmp.clear();
			}
			null_num++;
		}
		else{
			tmp.push_back(ofn.lpstrFile[i]);
			null_num = 0;
		}

		// 2�̘A�������k��������������΃��[�v�𔲂���
		if(null_num >= 2) break;
	}

	// �����t�@�C���I������fns�̃T�C�Y��3�ȏ�
	int n = (int)fns0.size();
	//vector<string> fns;
	if(n >= 3){
		fns.resize(n-1);
		for(int i = 0; i < n-1; ++i){
			// �f�B���N�g�����ƃt�@�C����������
			fns[i] = fns0[0]+"/"+fns0[i+1];
		}
		n -= 1;
	}
	else{
		// �P��t�@�C��
		fns.resize(1);
		fns[0] = fns0[0];
		n = 1;
	}

#else
	string fl_filter = GenFlFileFilter(descs, exts);
	//cout << fl_filter << endl;
	//fl_filter = "Movie Files (*.{avi,mp4,flv,mov,mkv})\tAll Files (*)";

	int t = (multi ? Fl_File_Chooser::MULTI : Fl_File_Chooser::SINGLE);
	Fl_File_Chooser *fc = new Fl_File_Chooser(".", fl_filter.c_str(), t, title.c_str());
	fc->show();
	while(fc->visible()){
		Fl::wait();
	}

	int n = fc->count();
	for(int i = 0; i < n; ++i){
		if(!fc->value(i+1)) continue;
		fns.push_back(fc->value(i+1));
	}
	n = (int)fns.size();

	delete fc;

#endif

	return n;
}


#endif // #ifndef _RX_FLTK_WIDGETS_H_