/*! @file
 @brief		gnuplot�N���X
 @author	hecomi (http://d.hatena.ne.jp/hecomi/)
 @date		December 09, 2010. ver. 1.03

 10/07/05 ����
 10/07/12 ���O�v���b�g/reset/PNG�o�͒ǉ�
 10/07/13 �}���`�v���b�g
 10/12/09 �w�b�_�t�@�C���݂̂ɕύX
*/

#pragma once
#pragma warning(disable: 4996)

#include <stdio.h>
#include <stdarg.h>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>

#define __KEYWAIT__ { int i; std::cin >> i; } 

#ifdef USE_BOOST
#include <boost/multi_array.hpp>
#endif

namespace gnuplot {

/*!
@brief gnuplot���������C�u����
	gnuplot�������N���X�Ɋւ��郉�C�u�����ł��D
	pgnuplot�ւƃp�C�v���q���C�֐���f�[�^��gnuplot��p���ăv���b�g����菕�������܂��D
	pgnuplot.exe�����ϐ��ɓo�^����Ă���K�v������܂��D
*/
class CGnuplot
{
private:
/* ----------------------------------------------------------------------
 �v���b�g�̃X�^�C����ύX�������ꍇ�͂�����ύX
---------------------------------------------------------------------- */
	/*!
	@brief �X�^�C��������
	*/
	void InitialStyle()
	{
		Command("set style line  1 linetype  1 linewidth 2");
		Command("set style line  2 linetype  2 linewidth 2");
		Command("set style line  3 linetype  3 linewidth 2");
		Command("set style line  4 linetype  4 linewidth 2");
		Command("set style line  5 linetype  5 linewidth 2");
		Command("set style line  6 linetype  6 linewidth 2");
		Command("set style line  7 linetype  7 linewidth 2");
		Command("set style line  8 linetype  8 linewidth 2");
		Command("set style line  9 linetype  9 linewidth 2");
		Command("set style line 10 linetype 10 linewidth 2");
	}

	/*!
	@brief �X�^�C��������
	*/
	void OutputStyle()
	{
		Command("set style line  1 linetype  1 linewidth 6");
		Command("set style line  2 linetype  2 linewidth 6");
		Command("set style line  3 linetype  3 linewidth 6");
		Command("set style line  4 linetype  4 linewidth 6");
		Command("set style line  5 linetype  5 linewidth 6");
		Command("set style line  6 linetype  6 linewidth 6");
		Command("set style line  7 linetype  7 linewidth 6");
		Command("set style line  8 linetype  8 linewidth 6");
		Command("set style line  9 linetype  9 linewidth 6");
		Command("set style line 10 linetype 10 linewidth 6");
	}

/* ----------------------------------------------------------------------
 �����o�ϐ�
---------------------------------------------------------------------- */
	//! �v���b�g����v�f�̊J�n�ʒu�ƏI���ʒu
	FILE* Fp;

	//! @brief �ꎞ�t�@�C����
	const std::string TempFileName;

	//! �v���b�g����v�f�̊J�n�ʒu�ƏI���ʒu
	template <class Container>
	struct PlotInfo
	{
		Container begin;	//!< �v���b�g�J�n�ʒu
		Container end;		//!< �v���b�g�I���ʒu
	};

	//! �}���`�v���b�g���ǂ���
	bool IfMultiplot;

	//! �v���b�g�t�@�C���p�i���o�[�����l
	const unsigned int DefaultNo;

	//! �v���b�g�t�@�C���p�i���o�[�iTemp(No).dat�j
	unsigned int No;

	//! �v���b�g�t�@�C���p�i���o�[�̍ő�l
	unsigned int MaxNo;

/* ----------------------------------------------------------------------
 �����ň����֐��B
---------------------------------------------------------------------- */
	//! �t�@�C�����擾
	std::string GetPlotFileName()
	{
		std::stringstream ss;
		ss << TempFileName << No << ".dat";
		return ss.str();
	}

	/*!
	@brief ������
	*/
	void Ini()
	{
		//Fp = _popen("pgnuplot", "w");
		Fp = _popen("D:/gnuplot/bin/pgnuplot.exe", "w");
		if (Fp == NULL) {
			printf("pipe error\n");
			exit(EXIT_FAILURE);
		}
		InitialStyle();
	}

	/*!
	@brief �o�b�t�@�t���b�V��
	*/
	void Flush()
	{
		fflush(Fp);
	}

	/*!
	@brief �v���b�g�R�}���h
	*/
	void Plot()
	{
		Command("plot '%s' w %s ls %d notitle", GetPlotFileName().c_str(), GetPlotType().c_str(), No);

		// �}���`�v���b�g�̏ꍇ�͔ԍ����Z�b�g
		if (IfMultiplot) {
			No++;
			if (No > MaxNo) {
				MaxNo = No;
			}
		}
	}

	/*!
	@brief s�v���b�g�R�}���h
	*/
	void SPlot(std::string option = "")
	{
		Command("set hidden3d");
		Command(("splot '%s' notitle " + option).c_str(), GetPlotFileName().c_str(), GetPlotType().c_str(), No);

		// �}���`�v���b�g�̏ꍇ�͔ԍ����Z�b�g
		if (IfMultiplot) {
			No++;
			if (No > MaxNo) {
				MaxNo = No;
			}
		}
	}

	/*!
	@brief 1�����v�f���v���b�g����
	@param[in] x �v���b�g�����i�[����x�f�[�^
	*/
	template <class T>
	void PlotX(PlotInfo<T> x)
	{
		T it = x.begin;
		std::ofstream fout(GetPlotFileName().c_str());
		if (fout.fail()) {
			std::cout << "Error! (@PlotX)" << std::endl;
			return;
		}
		while (it != x.end) {
			fout << *it << std::endl;
			it++;
		}
		Plot();
	}

	/*!
	@brief 2�����v�f���v���b�g����
	@param[in] x �v���b�g�����i�[����x�f�[�^
	@param[in] y �v���b�g�����i�[����y�f�[�^
	*/
	template <class T>
	void PlotXY(PlotInfo<T> x, PlotInfo<T> y)
	{
		T itX = x.begin, itY = y.begin;
		std::ofstream fout(GetPlotFileName().c_str());
		if (fout.fail()) {
			std::cout << "Error! (@PlotXY)" << std::endl;
			return;
		}
		while (itX != x.end && itY != y.end) {
			fout << *itX << " " << *itY << std::endl;
			itX++; itY++;
		}
		Plot();
	}

	/*!
	@brief 3�����v�f���v���b�g����
	@param[in] x �v���b�g�����i�[����x�f�[�^
	@param[in] y �v���b�g�����i�[����y�f�[�^
	@param[in] z �v���b�g�����i�[����z�f�[�^
	*/
	template <class T>
	void PlotXYZ(PlotInfo<T> x, PlotInfo<T> y, PlotInfo<T> z)
	{
		T itX = x.begin, itY = y.begin, itZ = z.begin;
		std::ofstream fout(GetPlotFileName().c_str());
		if (fout.fail()) {
			std::cout << "Error! (@PlotXYZ)" << std::endl;
			return;
		}
		while (itX != x.end && itY != y.end && itZ != z.end) {
			fout << *itX << " " << *itY << " " << *itZ << std::endl;
			itX++; itY++; itZ++;
		}
		SPlot();
	}

	/*!
	@brief �v���b�g�^�C�v�ɉ�������������擾
	*/
	std::string GetPlotType()
	{
		switch (PlotType) {
			case PLOT_TYPE_LINES:			return "lines";			break;
			case PLOT_TYPE_POINTS:			return "points";		break;
			case PLOT_TYPE_LINES_POINTS:	return "linespoints";	break;
			default:						return "lines";
		}
	}


public:
	/*! @brief �R���X�g���N�^ */
	CGnuplot() :
		TempFileName("Temp"),
		PlotType(PLOT_TYPE_LINES),
		IfMultiplot(false),
		DefaultNo(1),
		No(DefaultNo),
		MaxNo(DefaultNo)
	{
		Ini();
	}

	/*!
	@brief �R���X�g���N�^
	������Gnuplot���N�������ۂɁC�����̃v���b�g�f�[�^�t�@�C�����Փ˂���ƁC
	���]�̃v���b�g���ʂ������Ȃ��̂ŁC�v���b�g�f�[�^�t�@�C�������w�肵�Ȃ���΂Ȃ�Ȃ��D
	@param[in] fileName �v���b�g���ɐ�������f�[�^�t�@�C���̖��O
	*/
	CGnuplot(const char* fileName) :
		TempFileName(fileName),
		PlotType(PLOT_TYPE_LINES),
		IfMultiplot(false),
		DefaultNo(1),
		No(DefaultNo),
		MaxNo(DefaultNo)
	{
		Ini();
	}

	/*! @brief �f�X�g���N�^ */
	~CGnuplot()
	{
		Command("exit");

		// �v���b�g�t�@�C���̍폜
		for (No=DefaultNo; No<=MaxNo; No++) {
			remove(GetPlotFileName().c_str());
		}

		_pclose(Fp);
	}

	/*! @brief ����ɋ@�\���Ă��邩�ǂ��� */
	bool Check()
	{
		if (Fp == NULL) {
			return false;
		}
		return true;
	}
	
	/*!
	@brief printf���C�N��gnuplot�̃R�}���h�����s
	@param[in] format printf�ɗp����t�H�[�}�b�g
	@param[in] ... �ϒ�����
	*/
	void Command(const char* format, ...)
	{
		char buf[1024];
		va_list ap;

		va_start(ap, format);
		vsprintf(buf, format, ap);
		va_end(ap);

		fprintf(Fp, "%s\n", buf);
		Flush();
	}

	/*!
	@brief �֐����v���b�g
	Ex. DrawFunc("sin (x)")
	@param[in] format �v���b�g�Ώ̊֐�������
	*/
	void DrawFunc(const char* func)
	{
		Command("plot %s", func);
	}

	/*!
	@brief 1�����v�f���v���b�g����
	@param[in] cont �v���b�g�Ώ̃R���e�i
	*/
	template <class T, template <class A, class Allocator = std::allocator<A> > class Container>
	void Plot(Container<T> cont)
	{
		PlotInfo<Container<T>::iterator> pi = { cont.begin(), cont.end() };
		PlotX(pi);
	}

	/*!
	@brief 1�����v�f���v���b�g����i�z��j
	@param[in] cont �v���b�g�Ώ̔z��
	*/
	template <class T, int N>
	void Plot(T (&cont)[N])
	{
		PlotInfo<T*> pi = { &cont[0], &cont[N-1] };
		PlotX(pi);
	}

	/*!
	@brief 2�����v�f���v���b�g����
	@param[in] contX �v���b�g�Ώ̃R���e�i�ix�j
	@param[in] contY �v���b�g�Ώ̃R���e�i�iy�j
	*/
	template <class T, template <class A, class Allocator = std::allocator<A> > class Container>
	void Plot(Container<T> contX, Container<T> contY)
	{
		PlotInfo<Container<T>::iterator> 
			piX = { contX.begin(), contX.end() },
			piY = { contY.begin(), contY.end() };
		PlotXY(piX, piY);
	}

	/*!
	@brief 2�����v�f���v���b�g����i�z��j
	@param[in] contX �v���b�g�Ώ̔z��ix�j
	@param[in] contY �v���b�g�Ώ̔z��iy�j
	*/
	template <class T, int N, int M>
	void Plot(T (&contX)[N], T (&contY)[M])
	{
		PlotInfo<T*> 
			piX = { &contX[0], &contX[N] },
			piY = { &contY[0], &contY[M] };
		PlotXY(piX, piY);
	}

	/*!
	@brief 3�����v�f���v���b�g����
	@param[in] contX �v���b�g�Ώ̃R���e�i�ix�j
	@param[in] contY �v���b�g�Ώ̃R���e�i�iy�j
	@param[in] contZ �v���b�g�Ώ̃R���e�i�iz�j
	*/
	template <class T, template <class A, class Allocator = std::allocator<A> > class Container>
	void Plot(Container<T> contX, Container<T> contY, Container<T> contZ)
	{
		PlotInfo<Container<T>::iterator> 
			piX = { contX.begin(), contX.end() },
			piY = { contY.begin(), contY.end() },
			piZ = { contZ.begin(), contZ.end() };
		PlotXYZ(piX, piY, piZ);
	}

	/*!
	@brief 3�����v�f���v���b�g����i�z��j
	@param[in] contX �v���b�g�Ώ̔z��ix�j
	@param[in] contY �v���b�g�Ώ̔z��iy�j
	@param[in] contZ �v���b�g�Ώ̔z��iz�j
	*/
	template <class T, int N, int M, int L>
	void Plot(T (&contX)[N], T (&contY)[M], T (&contZ)[L])
	{
		PlotInfo<T*> 
			piX = { &contX[0], &contX[N] },
			piY = { &contY[0], &contY[M] },
			piZ = { &contZ[0], &contZ[L] };
		PlotXYZ(piX, piY, piZ);
	}

	/*!
	@brief 2�����}���`�v���b�g�iSTL�R���e�i����j�ix�����ʁj
	@param[in] x �v���b�g�px���x�N�g��
	@param[in] ys y�x�N�g���̏W��
	*/
	template <class T, template <class A, class Allocator = std::allocator<A> > class Container>
	void Multiplot(Container<T> x, Container<std::pair<std::string, Container<T> > > ys)
	{
		// �O�̂���
		No = DefaultNo;

		// �v���b�g�J�n
		Command("plot \\");

		// �t�@�C���֏����o��
		std::ofstream fout(GetPlotFileName().c_str());
		if (fout.fail()) {
			std::cout << "Error! (@Multiplot)" << std::endl;
			return;
		}

		// ��2�����̓W�J
		Container<std::pair<std::string, Container<T> > >::iterator itYs = ys.begin();
		unsigned int i = 1;
		while (itYs != ys.end()) {
			// �y�A�̓W�J
			std::pair<std::string, Container<T> > pair = *itYs;
			std::string title = pair.first;
			Container<T> y = pair.second;

			// �C�e���[�^�w��
			Container<T>::iterator itX = x.begin(), itY = y.begin();

			while (itX != x.end() && itY != y.end()) {
				fout << *itX << " " << *itY << std::endl;
				itX++; itY++;
			}
			fout << std::endl << std::endl;

			// �v���b�g
			Command("'%s' ind %d w %s ls %d ti '%s'\\", GetPlotFileName().c_str(), i-1, GetPlotType().c_str(), i, title.c_str());
			if (i < ys.size()) {
				Command(",\\");
			} else {
				Command("");
				break;
			}

			itYs++;
			i++;
		}
	}

	/*!
	@brief 2�����}���`�v���b�g�iSTL�R���e�i����j�ix���ʁj
	@param[in] plot �v���b�g�f�[�^
	*/
	template <class T, template <class A, class Allocator = std::allocator<A> > class Container>
	void Multiplot(Container<std::pair<std::string, std::pair<Container<T>, Container<T> > > > plot)
	{
		// �O�̂���
		No = DefaultNo;

		// �v���b�g�J�n
		Command("plot \\");

		// �t�@�C���֏����o��
		std::ofstream fout(GetPlotFileName().c_str());
		if (fout.fail()) {
			std::cout << "Error! (@Multiplot)" << std::endl;
			return;
		}

		// ��2�����̓W�J
		Container<std::pair<std::string, std::pair<Container<T>, Container<T> > > >::iterator it = plot.begin();
		unsigned int i = 1;
		while (it != plot.end()) {
			// �y�A�̓W�J
			std::pair<std::string, std::pair<Container<T>, Container<T> > > pair = *it;
			std::string title = pair.first;
			Container<T> x = pair.second.first;
			Container<T> y = pair.second.second;

			// �C�e���[�^�w��
			Container<T>::iterator itX = x.begin(), itY = y.begin();

			while (itX != x.end() && itY != y.end()) {
				fout << *itX << " " << *itY << std::endl;
				itX++; itY++;
			}
			fout << std::endl << std::endl;

			// �v���b�g
			Command("'%s' ind %d w %s ls %d ti '%s'\\", GetPlotFileName().c_str(), i-1, GetPlotType().c_str(), i, title.c_str());
			if (i < plot.size()) {
				Command(",\\");
			} else {
				Command("");
				break;
			}

			it++;
			i++;
		}
	}

/* ----------------------------------------------------------------------
 boost::multi_array ��3D�v���b�g�ł��܂��D
 �g���ۂ�USE_BOOST��define���Ă��������D
---------------------------------------------------------------------- */
#ifdef USE_BOOST
	/*!
	@brief 3�����v�f���v���b�g����i���ʁj
	@param[in] contX �v���b�g�Ώ̃R���e�i�ix�j
	@param[in] contY �v���b�g�Ώ̃R���e�i�iy�j
	@param[in] contZ �v���b�g�Ώ̃R���e�i�iz�j
	*/
	template <class T, template <class A, class Allocator = std::allocator<A> > class Container>
	void Plot(Container<T> contX, Container<T> contY, boost::multi_array<T, 2> contZ)
	{
		// �v���b�g�T�C�Y
		const int nX = contX.size();
		const int nY = contY.size();

		// �t�@�C���֏����o��
		std::ofstream fout(GetPlotFileName().c_str());
		if (fout.fail()) {
			std::cout << "Error! (@Plot)" << std::endl;
			return;
		}

		for (int i=0; i<nX; i++) {
			for (int j=0; j<nY; j++) {
				fout << contX[i] << " " << contY[j] << " " << contZ[i][j] << std::endl;
			}
			fout << std::endl;
		}

		SPlot("with pm3d");
	}
#endif

	/*!
	@brief 3�����v�f���v���b�g����i���ʁj
	@param[in] contX �v���b�g�Ώ̃R���e�i�ix�j
	@param[in] contY �v���b�g�Ώ̃R���e�i�iy�j
	@param[in] contZ �v���b�g�Ώ̃R���e�i�iz�j
	*/
	template <class T, template <class A, class Allocator = std::allocator<A> > class Container>
	void Plot(Container<T> contX, Container<T> contY, Container<Container<T> > contZ)
	{
		// �v���b�g�T�C�Y
		const int nX = contX.size();
		const int nY = contY.size();

		// �t�@�C���֏����o��
		std::ofstream fout(GetPlotFileName().c_str());
		if (fout.fail()) {
			std::cout << "Error! (@Plot)" << std::endl;
			return;
		}

		for (int i=0; i<nX; i++) {
			for (int j=0; j<nY; j++) {
				fout << contX[i] << " " << contY[j] << " " << contZ[i][j] << std::endl;
			}
			fout << std::endl;
		}

		SPlot("with pm3d");
	}

	/*!
	@brief �v���b�g�^�C�v
		PLOT_TYPE_LINES			���̂�
		PLOT_TYPE_POINTS		�_�̂�
		PLOT_TYPE_LINES_POINTS	���Ɠ_
	*/
	enum PLOT_TYPE {
		PLOT_TYPE_LINES,
		PLOT_TYPE_POINTS,
		PLOT_TYPE_LINES_POINTS
	} PlotType;

	/*!
	@brief �v���b�g�^�C�v�i���C�_�C���Ɠ_�j���Z�b�g
	*/
	void SetPlotType(const PLOT_TYPE pt)
	{
		PlotType = pt;
	}

	/*!
	@brief �}���`�v���b�g�i�����̐}���v���b�g�j�̐؂�ւ�
	@param[in] sw true: �}���`�v���b�g false: ����
	*/
	void SetMultiplot(const bool sw = true)
	{
		if (sw) {
			Command("set multiplot");
			IfMultiplot = true;

		} else {
			Command("unset multiplot");
			IfMultiplot = false;
			No = DefaultNo;
		}
	}

	/*!
	@brief ���x�����Z�b�g
	@param[in] title �^�C�g��
	*/
	void SetTitle(const char* title)
	{
		Command("set title '%s'", title);
	}

	/*!
	@brief ���x�����Z�b�g
	@param[in] labelX X���x����
	@param[in] labelY Y���x����
	*/
	void SetLabel(const char* labelX, const char* labelY)
	{
		Command("set xlabel '%s'", labelX);
		Command("set ylabel '%s'", labelY);
	}

	/*!
	@brief �v���b�g�͈͂��w��
	@param[in] min x���v���b�g�͈͍ŏ��l
	@param[in] min y���v���b�g�͈͍ŏ��l
	*/
	void SetXRange(const double min, const double max)
	{
		Command("set xrange [%f:%f]", min, max);
	}

	/*!
	@brief �v���b�g�͈͂��w��
	@param[in] min y���v���b�g�͈͍ŏ��l
	@param[in] min y���v���b�g�͈͍ŏ��l
	*/
	void SetYRange(const double min, const double max)
	{
		Command("set yrange [%f:%f]", min, max);
	}

	/*!
	@brief ���v���b�g
	*/
	void Replot()
	{
		Command("replot");
	}

	/*!
	@brief �ݒ�S���Z�b�g
	*/
	void Reset()
	{
		Command("reset");
	}

	/*!
	@brief �Бΐ��v���b�g�ɕύX
	*/
	void SetLogPlotY(const bool sw = true)
	{
		if (sw) {
			Command("set logscale y");
			Command("set format y \"10^{%%L}\"");
		} else {
			Command("set nolog y");
			Command("set format y");
		}
	}

	/*!
	@brief �t�Бΐ��v���b�g�ɕύX
	*/
	void SetLogPlotX(const bool sw = true)
	{
		if (sw) {
			Command("set logscale x");
			Command("set format x \"10^{%%L}\"");
		} else {
			Command("set nolog x");
			Command("set format x");
		}
	}

	/*!
	@brief ���ΐ��v���b�g�ɕύX
	*/
	void SetLogPlot(const bool sw = true)
	{
		SetLogPlotY(sw);
		SetLogPlotX(sw);
	}

	/*!
	@brief ���݃v���b�g���Ă���f�[�^���t�@�C���ɏ����o��
	���g��TempFileName�ɕۑ�����Ă���f�[�^���R�s�[���Ă��邾��
	@param[in] fileName �����o����t�@�C����
	*/
	void DumpToFile(const char* fileName)
	{
		std::ofstream fout(fileName, std::ios_base::binary);
		std::ifstream fin(GetPlotFileName().c_str(), std::ios_base::binary);
		if (fin.fail() || fout.fail()) {
			std::cout << "Error! (@DumpToFile)" << std::endl;
			return;
		}
		while (!fin.eof()) {
			const int BUF_SIZE = 4096;
			char buf[BUF_SIZE];
			fin.read(buf, BUF_SIZE);
			fout.write(buf, BUF_SIZE);
		}
	}

	/*!
	@brief ���݃v���b�g���Ă���f�[�^��EPS�ŏ����o��
	@param[in] fileName �����o����t�@�C����
	*/
	void DumpToEps(const char* fileName)
	{
		OutputStyle();
		Command("set term postscript eps enhanced color \"Tahoma\" 35");
		Command("set output '%s.eps'", fileName);
		Command("set size 1.6,1.6");
		Command("set ticscale 3");
		Command("replot");
		Command("set output");
		Command("set terminal window");
		Command("set ticscale");
		Command("set size");
		InitialStyle();
	}

	/*!
	@brief ���݃v���b�g���Ă���f�[�^��PNG�ŏ����o��
	@param[in] fileName �����o����t�@�C����
	*/
	void DumpToPng(const char* fileName)
	{
		Command("set term png");
		Command("set output '%s.png'", fileName);
		Command("set size");
		Command("set ticscale 3");
		Command("replot");
		Command("set output");
		Command("set terminal window");
		Command("set ticscale");
		InitialStyle();
	}

	/*!
	@brief �O���t�@�C������R�}���h��ǂݍ��ݎ��s
	@param[in] fileName �ǂݍ��݌��t�@�C����
	*/
	void CommandFromFile(const char* fileName)
	{
		std::ifstream fin(fileName);
		if (fin.fail()) {
			std::cout << "Error! '" << fileName << "' is not found. (@CommandFromFile)" << std::endl;
			return;
		}
		std::string str;
		while (!fin.eof()) {
			std::getline(fin, str);
			Command(str.c_str());
		}
	}

};


/* ----------------------------------------------------------------------
 �悭�g���Ǝv����^
---------------------------------------------------------------------- */
/*! @brief x�x�N�g���Œ�v���b�g�̗v�f */
typedef std::pair<std::string, std::vector<double> > StrVec;

/*! @brief x�x�N�g���Œ�v���b�g�̈����̌^ */
typedef std::vector<std::pair<std::string, std::vector<double> > > SemiMultiPlot;

/*! @brief x�x�N�g���σv���b�g�̗v�f */
typedef std::pair<std::string, std::pair<std::vector<double>, std::vector<double> > > StrVecVec;

/*! @brief x�x�N�g���σv���b�g */
typedef std::vector<std::pair<std::string, std::pair<std::vector<double>, std::vector<double> > > > MultiPlot;


/* ----------------------------------------------------------------------
 �����̌^��Ԃ��֐�
---------------------------------------------------------------------- */
static StrVecVec svv(const std::string str, const std::vector<double> vecX, const std::vector<double> vecY)
{
	return StrVecVec(str, make_pair(vecX, vecY));
}


}