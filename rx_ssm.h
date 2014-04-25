/*!
  @file rx_ssm.h
	
  @brief Screen Space Mesh�쐬
		- M. Muller et al., Screen space meshes, SCA2007, 2007. 

  @author Makoto Fujisawa
  @date 2011-05
*/
// FILE --rx_ssm.h--

#ifndef _RX_SSM_H_
#define _RX_SSM_H_

#pragma warning (disable: 4244)


//-----------------------------------------------------------------------------
// �C���N���[�h�t�@�C��
//-----------------------------------------------------------------------------
// C�W��
#include <cmath>
#include <cstdlib>

#include <iostream>
#include <fstream>

// STL
#include <vector>
#include <string>

// OpenGL
#include <GL/glew.h>
#include <GL/glut.h>

// �f�[�^�\��
#include "rx_utility.h"
#include "rx_matrix.h"
#include "rx_mesh.h"

// Marching Square�e�[�u��
#include "rx_ssm_tables.h"

// CUDA
#include "rx_cu_common.cuh"

using namespace std;



//-----------------------------------------------------------------------------
// ��`
//-----------------------------------------------------------------------------
#ifndef RXREAL
	#define RXREAL float
#endif
#ifndef RXREAL2
	#define RXREAL2 float2
#endif
#ifndef RXREAL3
	#define RXREAL3 float3
#endif



//! Screen Space�ł̃p�[�e�B�N��
struct rxSSParticle
{
	Vec3 xp;	//!< �X�N���[���X�y�[�X�ł̒��S���W
	Vec3 rp;	//!< �X�N���[���X�y�[�X�ł̔��a
};

//! Screen Space �ł̒��_
struct rxSSVertex
{
	Vec3 pos;		//!< Screen Space���W�l
	RXREAL depth;	//!< �f�v�X�l
	int edge;		//!< �G�b�W���_�ł��邩�ǂ���

	Vec3 avg_pos;	//!< �אڒ��_���ύ��W�l(�֊s�������p)
	RXREAL avg_num;	//!< �אڒ��_��(�֊s�������p)

	rxSSVertex(Vec3 p, int e = 0) : pos(Vec3(p[0], p[1], p[2])), depth(p[2]), edge(e) {}
};

//! �֊s�G�b�W
struct rxSSEdge
{
	Vec3 x0, x1;	//!< �[�_���W
	RXREAL depth;	//!< �G�b�W�f�v�X�l
	RXREAL d0, d1;	//!< �[�_�f�v�X�l
	int xy;			//!< �G�b�W����(x:0, y:1)
	int silhouette;	
	int front_vertex;	//!< �G�b�W���_�̃C���f�b�N�X
	RXREAL dx;			//!< �f�v�X�l���������[�_����G�b�W���_�܂ł̋���
};

//! ���b�V�������p�O���b�h
struct rxSSGrid
{
	int i, j;
	int node_vrts[4];	//!< �m�[�h���_�C���f�b�N�X
	int num_nv;			//!< �m�[�h���_��
	int edge_vrts[4];	//!< �G�b�W���_(front vertex)
	int num_ev;			//!< �G�b�W���_��(front vertex)
	int back_vrts[6];	//!< �G�b�W���_(back vertex, back-2 vertex)
	int num_bv;			//!< �G�b�W���_��(back vertex)

	RXREAL node_depth[4];	//!< �m�[�h�̃f�v�X�l

	int table_index0;	//!< �f�o�b�O�p:���b�V�����̂��߂̃C���f�b�N�X�l
	int table_index1;	//!< �f�o�b�O�p:���b�V�����̂��߂̃C���f�b�N�X�l
	int mesh_num;		//!< �f�o�b�O�p:���b�V����
	int mesh[6];		//!< �f�o�b�O�p:���b�V���C���f�b�N�X
	int back2;			//!< �f�o�b�O�p
	int v[14];			//!< �f�o�b�O�p
};


// �m�[�h���_
// 3 - 2
// |   |
// 0 - 1
	 
// �G�b�W���_
// - 2 -
// 3   1
// - 0 -


//-----------------------------------------------------------------------------
// �֐��v���g�^�C�v�錾
//-----------------------------------------------------------------------------
inline RXREAL GetDepthNearestT(double x, double y, int nx, int ny, double dx, double dy, const vector<RXREAL> &dmap);
inline RXREAL GetDepthInterpT(double x, double y, int nx, int ny, double dx, double dy, const vector<RXREAL> &dmap);


//-----------------------------------------------------------------------------
// Screen Space Mesh - ���N���X
// MARK:rxSSMesh
//-----------------------------------------------------------------------------
class rxSSMesh
{

protected:
	// �����o�ϐ�
	double m_fDmx, m_fDmy;	//!< �f�v�X�}�b�v�̊e�O���b�h��
	int m_iNumNodeVrts;		//!< �m�[�h���_��
	int m_iNumEdgeVrts;		//!< �G�b�W���_��
	int m_iNumMesh;			//!< ���b�V����

	double m_fSpacing;		//!< �f�v�X�}�b�v�̃T���v�����O�Ԋu
	double m_fPrtRad;		//!< �p�[�e�B�N���̔��a
	double m_fSSZmax;		//!< �֊s�ƂȂ�f�v�X����臒l
	int m_iNfilter;			//!< �f�v�X�l�������̃t�B���^�T�C�Y
	int m_iNiters;			//!< �֊s�������̔�����
	int m_iNgx, m_iNgy;		//!< ���b�V�������p�O���b�h�̉𑜓x

protected:
	//! �f�t�H���g�R���X�g���N�^
	rxSSMesh(){}	// �����w��Ȃ��ŃI�u�W�F�N�g���쐬�����̂�h��

public:
	//! �R���X�g���N�^
	rxSSMesh(double zmax, double h, double r, int n_filter, int n_iter)
	{
		m_fSSZmax = zmax;
		m_fSpacing = h;
		m_fPrtRad = r;
		m_iNfilter = n_filter;
		m_iNiters = n_iter;
		m_iNgx = m_iNgy = 24;

		m_iNumNodeVrts = 0;
		m_iNumEdgeVrts = 0;
		m_iNumMesh = 0;
	}

	//! �f�X�g���N�^
	virtual ~rxSSMesh()
	{
	}

public:
	//
	// �A�N�Z�X���\�b�h
	//
	//! �O���b�h��
	int GetNx(void) const { return m_iNgx; }
	int GetNy(void) const { return m_iNgy; }
	void GetDepthMapSize(int &nx, int &ny){ nx = m_iNgx+1; ny = m_iNgy+1; }

	//! �f�v�X�}�b�v�̃T���v�����O�Ԋu
	void   SetSpacing(double h){ m_fSpacing = h; }
	double GetSpacing(void) const { return m_fSpacing; }

	//! �p�[�e�B�N�����a
	void   SetRadius(double r){ m_fPrtRad = r; }
	double GetRadius(void) const { return m_fPrtRad; }

	//! �֊s�ƂȂ�f�v�X����臒l
	void   SetZMax(double r){ m_fSSZmax = r; }
	double GetZMax(void) const { return m_fSSZmax; }

	//! �f�v�X�l�������̃t�B���^�T�C�Y
	void SetFilterRadius(int r){ m_iNfilter = r; }
	int  GetFilterRadius(void) const { return m_iNfilter; }

	//! �֊s�������̔�����
	void SetSmoothIter(int r){ m_iNiters = r; }
	int  GetSmoothIter(void) const { return m_iNiters; }

	//! ���_��
	int GetVertexNum(void) const { return m_iNumNodeVrts+m_iNumEdgeVrts; }

	//! ���b�V����
	int GetMeshNum(void) const { return m_iNumMesh; }

	//! �f�v�X�}�b�v
	virtual RXREAL* GetDepthMap(void) = 0;

	//! ���_���
	virtual rxSSVertex* GetSSVertex(void) = 0;

	//! ���b�V�������p�O���b�h���
	virtual rxSSGrid GetMeshGrid(int idx) = 0;
	virtual rxSSGrid GetMeshGrid(int i, int j) = 0;

public:
	//
	// �����֐�
	//

	/*!
	 * �X�N���[���X�y�[�X���b�V������
	 * @param[in] proj OpenGL�������e�s��
	 * @param[in] modelview OpenGL���f���r���[�s��
	 * @param[in] W,H ��ʉ𑜓x
	 * @param[in] prts �p�[�e�B�N�����W
	 * @param[in] pnum �p�[�e�B�N����
	 * @param[out] vrts ���b�V�����_��
	 * @param[out] polys ���b�V����
	 * @param[in] filtering �t�B���^�����O�t���O(�f�v�X�l0x01, �֊s0x02)
	 * @param[in] debug_output ���ʂ̉�ʏo�̗͂L��
	 */
	virtual void CreateMesh(double *proj, double *modelview, int W, int H, vector<Vec3> &prts, int pnum, 
							vector<Vec3> &vrts, vector< vector<int> > &mesh, int filtering = 0, int debug_output = 0) = 0;

	/*!
	 * �X�N���[���X�y�[�X���b�V������(�@���v�Z�܂�)
	 * @param[in] proj OpenGL�������e�s��
	 * @param[in] modelview OpenGL���f���r���[�s��
	 * @param[in] W,H ��ʉ𑜓x
	 * @param[in] prts �p�[�e�B�N�����W
	 * @param[in] pnum �p�[�e�B�N����
	 * @param[out] vrts ���b�V�����_��
	 * @param[out] nrms ���_�@����
	 * @param[out] polys ���b�V����
	 * @param[in] filtering �t�B���^�����O�t���O(�f�v�X�l0x01, �֊s0x02)
	 * @param[in] debug_output ���ʂ̉�ʏo�̗͂L��
	 */
	virtual void CreateMesh(double *proj, double *modelview, int W, int H, vector<Vec3> &prts, int pnum, 
							vector<Vec3> &vrts, vector<Vec3> &nrms, vector< vector<int> > &mesh, 
							int filtering = 0, int debug_output = 0) = 0;

	/*!
	 * �}�b�v��O���b�h�z��̃T�C�Y��ύX
	 * @param[in] W,H ��ʉ𑜓x
	 * @param[in] spacing �f�v�X�}�b�v�̃T���v�����O�Ԋu
	 */
	virtual void Resize(int W, int H, int spacing = -1) = 0;

	/*!
	 * ���_�@���v�Z
	 * @param[in] vrts ���_���W
	 * @param[in] nvrts ���_��
	 * @param[in] tris �O�p�`�|���S���􉽏��
	 * @param[in] ntris �O�p�`�|���S����
	 * @param[out] nrms �@��
	 */
	virtual void CalVertexNormals(const vector<Vec3> &vrts, unsigned int nvrts, const vector< vector<int> > &tris, unsigned int ntris, 
								  vector<Vec3> &nrms){}

public:
	// OpenGL�`��
	virtual void DrawSilhouetteEdge(void){}
	virtual void DrawSSVertex(Vec3 node_color, Vec3 edge_color){}
	virtual void DrawMeshGrid(int grid, const Vec3 colors[]){}
	virtual void DrawField(double minpos[2], double maxpos[2]){}

	virtual void DrawSSEdge(void){}

public:
	// �f�o�b�O�p
	virtual void OutputGridInfo(int grid){}
	virtual void OutputGridVertexInfo(int grid){}
};


//-----------------------------------------------------------------------------
// Screen Space Mesh - CPU�ł̎���
// MARK:rxSSMeshCPU
//-----------------------------------------------------------------------------
class rxSSMeshCPU : public rxSSMesh
{
protected:
	// �e�[�u���C���f�b�N�X�X�V�֐��|�C���^
	typedef void (rxSSMeshCPU::*FuncTableIndex)(int&, int&, int[], rxSSGrid*);

protected:
	// �����o�ϐ�
	vector<RXREAL> m_vSSDMap;		//!< �f�v�X�}�b�v

	vector<rxSSParticle> m_vSSPrts;		//!< ���K�����W�n�ł̃p�[�e�B�N��

	vector<rxSSEdge> m_vSSEdge;			//!< �֊s�G�b�W
	vector<rxSSVertex> m_vSSEdgeVertex;	//!< �G�b�W���_
	vector<rxSSVertex> m_vSSVertex;		//!< �X�N���[���X�y�[�X���b�V�����_

	vector<rxSSGrid> m_vSSMGrid;		//!< ���b�V�������p�O���b�h
	vector<int> m_vMeshGrid;			//!< �f�o�b�O�p:���b�V����������O���b�h

	vector< vector<double> > m_vFilter;	//!< �������pbinomial�W��

	FuncTableIndex m_FuncTableIndex[25];//!< �e�[�u���C���f�b�N�X�X�V�֐��|�C���^

protected:
	//! �f�t�H���g�R���X�g���N�^
	rxSSMeshCPU(){}	// �����w��Ȃ��ŃI�u�W�F�N�g���쐬�����̂�h��

public:
	//! �R���X�g���N�^
	rxSSMeshCPU(double zmax, double h, double r, int n_filter, int n_iter);

	//! �f�X�g���N�^
	virtual ~rxSSMeshCPU();

public:
	//
	// �A�N�Z�X���\�b�h
	//
	
	//! �f�v�X�}�b�v
	virtual RXREAL* GetDepthMap(void);

	//! ���_���
	virtual rxSSVertex* GetSSVertex(void);

	//! ���b�V�������p�O���b�h���
	virtual rxSSGrid GetMeshGrid(int idx);
	virtual rxSSGrid GetMeshGrid(int i, int j){ return GetMeshGrid(i+j*m_iNgx); }

public:
	//
	// �����֐�
	//

	/*!
	 * �X�N���[���X�y�[�X���b�V������
	 * @param[in] proj OpenGL�������e�s��
	 * @param[in] modelview OpenGL���f���r���[�s��
	 * @param[in] W,H ��ʉ𑜓x
	 * @param[in] prts �p�[�e�B�N�����W
	 * @param[in] pnum �p�[�e�B�N����
	 * @param[out] vrts ���b�V�����_��
	 * @param[out] polys ���b�V����
	 * @param[in] filtering �t�B���^�����O�t���O(�f�v�X�l0x01, �֊s0x02)
	 * @param[in] debug_output ���ʂ̉�ʏo�̗͂L��
	 */
	virtual void CreateMesh(double *proj, double *modelview, int W, int H, vector<Vec3> &prts, int pnum, 
							vector<Vec3> &vrts, vector< vector<int> > &mesh, int filtering = 0, int debug_output = 0);

	/*!
	 * �X�N���[���X�y�[�X���b�V������(�@���v�Z�܂�)
	 * @param[in] proj OpenGL�������e�s��
	 * @param[in] modelview OpenGL���f���r���[�s��
	 * @param[in] W,H ��ʉ𑜓x
	 * @param[in] prts �p�[�e�B�N�����W
	 * @param[in] pnum �p�[�e�B�N����
	 * @param[out] vrts ���b�V�����_��
	 * @param[out] nrms ���_�@����
	 * @param[out] polys ���b�V����
	 * @param[in] filtering �t�B���^�����O�t���O(�f�v�X�l0x01, �֊s0x02)
	 * @param[in] debug_output ���ʂ̉�ʏo�̗͂L��
	 */
	virtual void CreateMesh(double *proj, double *modelview, int W, int H, vector<Vec3> &prts, int pnum, 
							vector<Vec3> &vrts, vector<Vec3> &nrms, vector< vector<int> > &mesh, 
							int filtering = 0, int debug_output = 0);

	/*!
	 * �}�b�v��O���b�h�z��̃T�C�Y��ύX
	 * @param[in] W,H ��ʉ𑜓x
	 * @param[in] spacing �f�v�X�}�b�v�̃T���v�����O�Ԋu
	 */
	virtual void Resize(int W, int H, int spacing = -1);

	/*!
	 * ���_�@���v�Z
	 * @param[in] vrts ���_���W
	 * @param[in] nvrts ���_��
	 * @param[in] tris �O�p�`�|���S���􉽏��
	 * @param[in] ntris �O�p�`�|���S����
	 * @param[out] nrms �@��
	 */
	virtual void CalVertexNormals(const vector<Vec3> &vrts, unsigned int nvrts, const vector< vector<int> > &tris, unsigned int ntris, vector<Vec3> &nrms);

	/*!
	 * �p�[�e�B�N�����W�Ɣ��a�𓧎����e�ϊ����ăf�v�X�}�b�v�𐶐�
	 * @param[in] P �������e�s��
	 * @param[in] MV ���f���r���[�s��
	 * @param[in] W,H ��ʉ𑜓x
	 * @param[in] prts �p�[�e�B�N�����W
	 * @param[in] pnum �p�[�e�B�N����
	 */
	void CalDepthMap(const rxMatrix4 &P, const rxMatrix4 &MV, int W, int H, vector<Vec3> &prts, int pnum);

	/*!
	 * �f�v�X�}�b�v�ɕ��������{��
	 * @param[inout] dmap �f�v�X�}�b�v
	 * @param[in] nx,ny �}�b�v�𑜓x
	 * @param[in] n_filter �t�B���^�����O��
	 */
	void ApplyDepthFilter(vector<RXREAL> &dmap, int nx, int ny, int n_filter);

	/*!
	 * �֊s�ɕ��������{��
	 * @param[inout] ssvrts �X�N���[�����W�n�ł̒��_��
	 * @param[in] polys ���b�V��(�\�����钸�_��)
	 * @param[in] n_iter �t�B���^�����O������
	 */
	void ApplySilhoutteSmoothing(vector<rxSSVertex> &ssvrts, const vector< vector<int> > &polys, int n_iter);

	/*!
	 * �֊s�G�b�W�̌��o��front edge vertex�̌v�Z
	 * @param[in] nx,ny ���b�V�������p�O���b�h�𑜓x
	 * @param[in] dw,dh ���b�V�������p�O���b�h��
	 * @param[in] dgrid ���b�V�������p�f�v�X�}�b�v
	 * @param[in] ssprts Screen Space�ł̃p�[�e�B�N��
	 * @param[in] W,H �X�N���[���𑜓x
	 * @return ���o���ꂽ�֊s�G�b�W��
	 */
	int DetectSilhouetteEdgeVertex(int nx, int ny, double dw, double dh, const vector<RXREAL> &dgrid, 
								   const vector<rxSSParticle> &ssprts, int W, int H);

	/*!
	 * �֊s�G�b�W�̌��o
	 * @param[in] nx,ny ���b�V�������p�O���b�h�𑜓x
	 * @param[in] dw,dh ���b�V�������p�O���b�h��
	 * @param[in] dgrid ���b�V�������p�f�v�X�}�b�v
	 * @return ���o���ꂽ�֊s�G�b�W��
	 */
	int DetectSilhouetteEdge(int nx, int ny, double dw, double dh, const vector<RXREAL> &dgrid);


	/*!
	 * �m�[�h���_����
	 * @param[in] nx,ny ���b�V�������p�O���b�h�𑜓x
	 * @param[in] dw,dh ���b�V�������p�O���b�h��
	 * @param[in] dgrid ���b�V�������p�f�v�X�}�b�v
	 * @return �������ꂽ�m�[�h���_��
	 */
	int CalNodeVertex(int nx, int ny, double dw, double dh, const vector<RXREAL> &dgrid);

	/*!
	 * �G�b�W���_����
	 * @param[in] nx,ny ���b�V�������p�O���b�h�𑜓x
	 * @param[in] dw,dh ���b�V�������p�O���b�h��
	 * @param[in] dgrid ���b�V�������p�f�v�X�}�b�v
	 * @param[in] edges �֊s�G�b�W���X�g
	 * @return �������ꂽ�G�b�W���_��
	 */
	int CalEdgeVertex(int nx, int ny, double dw, double dh, const vector<RXREAL> &dgrid, const vector<rxSSEdge> &edges);

	/*!
	 * �O�p�`���b�V������
	 * @param[in] nx,ny ���b�V�������p�O���b�h�𑜓x
	 * @param[in] dgrid ���b�V���O���b�h
	 * @param[out] polys �O�p�`���b�V��
	 * @param[in] vstart ���_�C���f�b�N�X�̎n�_
	 * @return �������ꂽ�O�p�`���b�V����
	 */
	int CalMesh(int nx, int ny, vector<rxSSGrid> &grid, vector< vector<int> > &polys, int vstart = 0);

protected:
	// �񕪒T���ɂ�钸�_�ʒu�v�Z
	double binarySearchDepth(Vec3 v1, Vec3 v2, Vec3 &vr, double zmax);

	// �f�v�X�l�̎Q��
	RXREAL getDepthValue(Vec3 x)
	{
		return GetDepthInterpT(x[0], x[1], m_iNgx+1, m_iNgy+1, m_fDmx, m_fDmy, m_vSSDMap);
		//return GetDepthNearestT(x[0], x[1], m_iNgx+1, m_iNgy+1, m_fDmx, m_fDmy, m_vSSDMap);
	}

	// ���b�V���e�[�u���C���f�b�N�X�̍X�V
	void updateTableIndexAny(int &table_index, int &vrot, int v[], rxSSGrid *g);
	void updateTableIndexE0N4(int &table_index, int &vrot, int v[], rxSSGrid *g);
	void updateTableIndexE1(int &table_index, int &vrot, int v[], rxSSGrid *g);
	void updateTableIndexE2N4(int &table_index, int &vrot, int v[], rxSSGrid *g);
	void updateTableIndexE3N23(int &table_index, int &vrot, int v[], rxSSGrid *g);
	void updateTableIndexE3N4(int &table_index, int &vrot, int v[], rxSSGrid *g);
	void updateTableIndexE4N2(int &table_index, int &vrot, int v[], rxSSGrid *g);
	void updateTableIndexE4N3(int &table_index, int &vrot, int v[], rxSSGrid *g);
	void updateTableIndexE4N4(int &table_index, int &vrot, int v[], rxSSGrid *g);


public:
	// OpenGL�`��
	virtual void DrawSilhouetteEdge(void);
	virtual void DrawSSVertex(Vec3 node_color, Vec3 edge_color);
	virtual void DrawMeshGrid(int grid, const Vec3 colors[]);
	virtual void DrawField(double minpos[2], double maxpos[2]);

	virtual void DrawSSEdge(void);

public:
	// �f�o�b�O�p
	virtual void OutputGridInfo(int grid);
	virtual void OutputGridVertexInfo(int grid);

	int Mesh2Grid(int mesh_idx){ return m_vMeshGrid[mesh_idx]; }

};


//-----------------------------------------------------------------------------
// Screen Space Mesh - GPU�ł̎���
// MARK:rxSSMeshGPU
//-----------------------------------------------------------------------------
class rxSSMeshGPU : public rxSSMesh
{
protected:
	// �����o�ϐ�
	vector<float> m_vBinomials;			//!< �������pbinomial�W��(1D)

	vector<RXREAL> m_vSSDMap;			//!< �f�v�X�}�b�v
	vector<rxSSVertex> m_vSSVertex;		//!< �X�N���[���X�y�[�X���b�V�����_

	// CUDA�p�ϐ�
	RXREAL *m_dSSPrtsCen;				//!< �X�N���[�����W�n�ł̃p�[�e�B�N�����W
	RXREAL *m_dSSPrtsRad;				//!< �X�N���[�����W�n�ł̃p�[�e�B�N�����a

	RXREAL *m_dSSDMap;					//!< �f�v�X�}�b�v
	rxVPacke m_dSSEdge;					//!< �֊s�G�b�W
	rxVPackf m_dSSNodeVrts;				//!< �m�[�h���_
	rxVPackf m_dSSEdgeFrontVrts;		//!< �G�b�W���_(front)
	rxVPackf m_dSSEdgeBackVrts;			//!< �G�b�W���_(back)
	rxVPackf m_dSSEdgeBack2Vrts;		//!< �G�b�W���_(back-2)
	RXREAL *m_dSSVertex;				//!< �S���_���(�m�[�h���_�C�O�ʃG�b�W���_�C�w�ʃG�b�W���_�C�Ŕw�ʃG�b�W���_)
	rxSSGridG *m_dSSMGrid;				//!< ���b�V�������p�O���b�h
	unsigned int *m_dSSTriNum;					//!< �O���b�h���̃��b�V����
	unsigned int *m_dSSTriNumScan;				//!< �O���b�h���̃��b�V������Scan
	//unsigned int *m_dSSTriArray;				//!< �O�p�`���b�V��

	RXREAL *m_hSSPrtsCen;				//!< �X�N���[�����W�n�ł̃p�[�e�B�N�����W(�f�o�C�X�������Ƃ̌����p)
	RXREAL *m_hSSPrtsRad;				//!< �X�N���[�����W�n�ł̃p�[�e�B�N�����a(�f�o�C�X�������Ƃ̌����p)

	rxSSEdgeG *m_hSSEdge;				//!< �X�N���[�����W�n�ł̃p�[�e�B�N�����W(�f�o�C�X�������Ƃ̌����p)
	RXREAL *m_hSSVertex;				//!< �X�N���[�����W�n�ł̃p�[�e�B�N�����W(�f�o�C�X�������Ƃ̌����p)

	rxSsmParams m_ssm_params;				//!< ���b�V�����p�����[�^


protected:
	//! �f�t�H���g�R���X�g���N�^
	rxSSMeshGPU(){}	// �����w��Ȃ��ŃI�u�W�F�N�g���쐬�����̂�h��

public:
	//! �R���X�g���N�^
	rxSSMeshGPU(double zmax, double h, double r, int n_filter, int n_iter);

	//! �f�X�g���N�^
	virtual ~rxSSMeshGPU();

public:
	//
	// �A�N�Z�X���\�b�h
	//
	
	//! �f�v�X�}�b�v
	virtual RXREAL* GetDepthMap(void);

	//! ���_���
	virtual rxSSVertex* GetSSVertex(void);

	//! ���b�V�������p�O���b�h���
	virtual rxSSGrid GetMeshGrid(int idx);
	virtual rxSSGrid GetMeshGrid(int i, int j){ return GetMeshGrid(i+j*m_iNgx); }

public:
	//
	// �����֐�
	//

	/*!
	 * �X�N���[���X�y�[�X���b�V������
	 * @param[in] proj OpenGL�������e�s��
	 * @param[in] modelview OpenGL���f���r���[�s��
	 * @param[in] W,H ��ʉ𑜓x
	 * @param[in] prts �p�[�e�B�N�����W
	 * @param[in] pnum �p�[�e�B�N����
	 * @param[out] vrts ���b�V�����_��
	 * @param[out] polys ���b�V����
	 * @param[in] filtering �t�B���^�����O�t���O(�f�v�X�l0x01, �֊s0x02)
	 * @param[in] debug_output ���ʂ̉�ʏo�̗͂L��
	 */
	virtual void CreateMesh(double *proj, double *modelview, int W, int H, vector<Vec3> &prts, int pnum, 
							vector<Vec3> &vrts, vector< vector<int> > &mesh, 
							int filtering = 0, int debug_output = 0);

	/*!
	 * �X�N���[���X�y�[�X���b�V������(�@���v�Z�܂�)
	 * @param[in] proj OpenGL�������e�s��
	 * @param[in] modelview OpenGL���f���r���[�s��
	 * @param[in] W,H ��ʉ𑜓x
	 * @param[in] prts �p�[�e�B�N�����W
	 * @param[in] pnum �p�[�e�B�N����
	 * @param[out] vrts ���b�V�����_��
	 * @param[out] nrms ���_�@����
	 * @param[out] polys ���b�V����
	 * @param[in] filtering �t�B���^�����O�t���O(�f�v�X�l0x01, �֊s0x02)
	 * @param[in] debug_output ���ʂ̉�ʏo�̗͂L��
	 */
	virtual void CreateMesh(double *proj, double *modelview, int W, int H, vector<Vec3> &prts, int pnum, 
							vector<Vec3> &vrts, vector<Vec3> &nrms, vector< vector<int> > &mesh, 
							int filtering = 0, int debug_output = 0);

	/*!
	 * �X�N���[���X�y�[�X���b�V������
	 *  - �p�[�e�B�N�����f�o�C�X�������Ɋi�[����Ă���ꍇ
	 *  - �@�����v�Z
	 * @param[in] proj OpenGL�������e�s��
	 * @param[in] modelview OpenGL���f���r���[�s��
	 * @param[in] W,H ��ʉ𑜓x
	 * @param[in] prts �p�[�e�B�N�����W
	 * @param[in] pnum �p�[�e�B�N����
	 * @param[out] dvrt ���b�V�����_��(�f�o�C�X������)
	 * @param[out] dtri ���b�V����(�f�o�C�X������)
	 * @param[in] filtering �t�B���^�����O�t���O(�f�v�X�l0x01, �֊s0x02)
	 * @param[in] debug_output ���ʂ̉�ʏo�̗͂L��
	 */
	void CreateMeshD(double *proj, double *modelview, int W, int H, RXREAL *dppos, RXREAL *dprad, int pnum, int pdim, 
					 RXREAL* &dvrt, unsigned int* &dtri, int filtering = 0, int debug_output = 0);
	/*!
	 * �X�N���[���X�y�[�X���b�V������
	 *  - �p�[�e�B�N�����f�o�C�X�������Ɋi�[����Ă���ꍇ
	 *  - �@�����v�Z
	 * @param[in] proj OpenGL�������e�s��
	 * @param[in] modelview OpenGL���f���r���[�s��
	 * @param[in] W,H ��ʉ𑜓x
	 * @param[in] prts �p�[�e�B�N�����W
	 * @param[in] pnum �p�[�e�B�N����
	 * @param[out] dvrt ���_��(�f�o�C�X������)
	 * @param[out] dnrm ���_�@����(�f�o�C�X������)
	 * @param[out] dtri ���b�V����(�f�o�C�X������)
	 * @param[in] filtering �t�B���^�����O�t���O(�f�v�X�l0x01, �֊s0x02)
	 * @param[in] debug_output ���ʂ̉�ʏo�̗͂L��
	 */
	void CreateMeshD(double *proj, double *modelview, int W, int H, RXREAL *dppos, RXREAL *dprad, int pnum, int pdim, 
					 RXREAL* &dvrt, RXREAL* &dnrm, unsigned int* &dtri, int filtering = 0, int debug_output = 0);

	/*!
	 * �X�N���[���X�y�[�X���b�V������
	 *  - �p�[�e�B�N�����f�o�C�X�������Ɋi�[����Ă���ꍇ
	 *  - VBO�ɏo��
	 * @param[in] proj OpenGL�������e�s��
	 * @param[in] modelview OpenGL���f���r���[�s��
	 * @param[in] W,H ��ʉ𑜓x
	 * @param[in] prts �p�[�e�B�N�����W
	 * @param[in] pnum �p�[�e�B�N����
	 * @param[out] uvrt_vbo ���b�V�����_��VBO
	 * @param[out] unrm_vbo ���b�V�����_��VBO
	 * @param[out] utri_vbo ���b�V�����_��VBO
	 * @param[in] filtering �t�B���^�����O�t���O(�f�v�X�l0x01, �֊s0x02)
	 * @param[in] debug_output ���ʂ̉�ʏo�̗͂L��
	 */
	void CreateMeshVBO(double *proj, double *modelview, int W, int H, RXREAL *dppos, RXREAL *dprad, int pnum, int pdim, 
					   GLuint &uVrtVBO, GLuint &uNrmVBO, GLuint &uTriVBO, int filtering = 0, int debug_output = 0);

	/*!
	 * �}�b�v��O���b�h�z��̃T�C�Y��ύX
	 * @param[in] W,H ��ʉ𑜓x
	 * @param[in] spacing �f�v�X�}�b�v�̃T���v�����O�Ԋu
	 */
	virtual void Resize(int W, int H, int spacing = -1);

	/*!
	 * ���_�@���v�Z
	 * @param[in] vrts ���_���W
	 * @param[in] nvrts ���_��
	 * @param[in] tris �O�p�`�|���S���􉽏��
	 * @param[in] ntris �O�p�`�|���S����
	 * @param[out] nrms �@��
	 */
	virtual void CalVertexNormals(const vector<Vec3> &vrts, unsigned int nvrts, const vector< vector<int> > &tris, unsigned int ntris, 
								  vector<Vec3> &nrms);

protected:
	//! GPU���̃p�����[�^�t�@�C�����X�V
	void updateParams(int W, int H, const rxMatrix4 &P, const rxMatrix4 &PMV);

	void calVertexNormalD(RXREAL* dvrt, unsigned int nvrts, unsigned int* dtri, unsigned int ntris, RXREAL* &nrms);

public:
	// OpenGL�`��
	virtual void DrawSilhouetteEdge(void);
	virtual void DrawSSVertex(Vec3 node_color, Vec3 edge_color);
	virtual void DrawMeshGrid(int grid, const Vec3 colors[]);
	virtual void DrawField(double minpos[2], double maxpos[2]);

	virtual void DrawSSEdge(void);

public:
	// �f�o�b�O�p
	virtual void OutputGridInfo(int grid);
	virtual void OutputGridVertexInfo(int grid);

};



//-----------------------------------------------------------------------------
// �֐�
//-----------------------------------------------------------------------------
/*!
 * �N�����v
 * @param[inout] x �l
 * @param[in] a,b [a,b]�ŃN�����v
 */
template<class T> 
inline void RX_CLAMP2(T &x, const T &a, const T &b){ x = (x < a ? a : (x > b ? b : x)); }

inline void Swap(Vec4 &a, Vec4 &b)
{
	Vec4 tmp = a;
	a = b;
	b = tmp;
}


inline RXREAL GetDepthInterpT(double x, double y, int nx, int ny, double dx, double dy, const vector<RXREAL> &dmap)
{
	int ix0 = (int)(x/dx);
	int ix1 = ix0+1;
	double ddx = x/dx-ix0;

	RX_CLAMP2(ix0, 0, nx-1);
	RX_CLAMP2(ix1, 0, nx-1);

	int iy0 = (int)(y/dx);
	int iy1 = iy0+1;
	double ddy = y/dy-iy0;

	RX_CLAMP2(iy0, 0, ny-1);
	RX_CLAMP2(iy1, 0, ny-1);

	RXREAL d00 = dmap[iy0*nx+ix0];
	RXREAL d01 = dmap[iy0*nx+ix1];
	RXREAL d10 = dmap[iy1*nx+ix0];
	RXREAL d11 = dmap[iy1*nx+ix1];

	return (d00*(1.0-ddx)+d01*ddx)*(1.0-ddy)+(d10*(1.0-ddx)+d11*ddx)*ddy;
}


inline RXREAL GetDepthNearestT(double x, double y, int nx, int ny, double dx, double dy, const vector<RXREAL> &dmap)
{
	int ix = (int)(x/dx);
	RX_CLAMP2(ix, 0, nx-1);

	int iy = (int)(y/dx);
	RX_CLAMP2(iy, 0, ny-1);

	return dmap[iy*nx+ix];
}

// �񕪒T��(1D)
inline double rtbis(double func(const double), const double x1, const double x2, const double xacc)
{
	const int JMAX = 40;
	double dx, f, fmid, xmid, rtb;

	f = func(x1);
	fmid = func(x2);
	if(f*fmid >= 0.0) return 0.0;

	rtb = f < 0.0 ? (dx = x2-x1, x1) : (dx = x1-x2, x2);
	for(int j = 0; j < JMAX; ++j){
		dx *= 0.5;
		xmid = rtb+dx;

		fmid = func(xmid);

		if(fmid <= 0.0){
			rtb = xmid;
		}

		if(fabs(dx) < xacc || fmid == 0.0){
			return rtb;
		}
	}
	
	return 0.0;
}


/*!
 * OpenGL�ϊ��s���rxMatrix4�ɕϊ�
 *  column major �� row major �ɕϊ�����rxMatrix4�Ɋi�[
 * @param[in] m OpenGL�ϊ��s��(	glGetDoublev(GL_PROJECTION_MATRIX, m); �ȂǂŎ擾)
 * @return �ϊ���̍s��
 */
inline rxMatrix4 GetMatrixGL(double *m)
{
	return rxMatrix4(m[0], m[4], m[8],  m[12], 
					 m[1], m[5], m[9],  m[13], 
					 m[2], m[6], m[10], m[14], 
					 m[3], m[7], m[11], m[15]);
}


/*!
 * ������2�i��������ɕϊ�
 * @param[in] x ���̐���
 * @param[in] bit 2�i������
 * @return 2�i��������
 */
inline string GetBitArray(int x, int bit)
{
	string s;
	s.resize(bit, '0');
	for(int i = 0; i < bit; ++i){
		s[bit-i-1] = ((x >> i) & 0x01) ? '1' : '0';
	}
	return s;
}

/*!
 * �f�v�X�l��[0,1]�ɕϊ�
 * @param[in] depth �f�v�X�l
 * @param[in] w division�l
 * @return [0,1]�̒l
 */
inline double RX_DEPTH2COLORf(RXREAL depth)
{
	if(depth == RX_FEQ_INF){
		depth = 1.0;
	}
	else{
		depth *= 0.5;
	}
	return 1.0-RX_CLAMP(depth, (RXREAL)0.0, (RXREAL)1.0);
}
inline unsigned char RX_DEPTH2COLOR(RXREAL d){ return (unsigned char)(RX_DEPTH2COLORf(d)*255); }


/*!
 * �������p��binominal�W�����v�Z
 * @param[in] b ����
 * @return binominal�W����
 */
static vector< vector<double> > CalBinomials(int b)
{
	vector< vector<double> > bs;
	vector<double> f, tmp;
	f.resize(b+1);
	tmp.resize(b+1);

	bs.resize(b+1);

	double a = 1.0;

	for(int i = 0; i < b+1; ++i){
		f[i]   = (i == 0 ? 1 : 0);
		tmp[i] = (i == 0 ? 1 : 0);	
	}

	for(int k = 0; k < b+1; ++k){
		for(int i = 1; i < k+1; ++i){
			tmp[i] = f[i-1]+f[i];
		}

		for(int i = 1; i < k+1; ++i){
			f[i] = tmp[i];
		}

		bs[k].resize(k+1);
		for(int i = 0; i < k+1; ++i){
			bs[k][i] = f[i]*a;
		}

		a *= 0.5;
	}

	return bs;
}

/*!
 * �������p��binominal�W�����v�Z(1�����z���, �t�B���^�ɗp���镨�̂�)
 * @param[in] r �t�B���^��(����b = 2*r+1)
 * @return binominal�W����
 */
static vector<float> CalBinomialsForFilter(int r)
{
	int b = 2*r+1;
	vector<float> bs;
	vector<float> f, tmp;
	f.resize(b+1);
	tmp.resize(b+1);

	bs.resize((r+1)*(r+1));

	float a = 1.0;

	for(int i = 0; i < b+1; ++i){
		f[i]   = (i == 0 ? 1.0f : 0.0f);
		tmp[i] = (i == 0 ? 1.0f : 0.0f);	
	}

	int c = 0;
	for(int k = 0; k < b+1; ++k){
		for(int i = 1; i < k+1; ++i){
			tmp[i] = f[i-1]+f[i];
		}

		for(int i = 1; i < k+1; ++i){
			f[i] = tmp[i];
		}

		if(!(k%2)){
			for(int i = 0; i < k+1; ++i){
				bs[c++] = f[i]*a;
			}
		}

		a *= 0.5;
	}

	return bs;
}

/*!
 * �����Ɖ~�̌�������(2D)
 * @param[in] A,B �����̗��[�_���W
 * @param[in] C �~�̒��S
 * @param[in] r �~�̔��a
 * @param[out] P ��_���W
 * @return ��_��
 */
static int LineCircleIntersection(const Vec2 &A, const Vec2 &B, const Vec2 &C, const double &r, Vec2 P[2], double t[2])
{
	double rr = r*r;
	Vec2 AC = C-A;
	Vec2 BC = C-B;

	Vec2 v = B-A;
	double l = norm(v);
	v /= l;

	double td = dot(v, AC);
	Vec2 D = A+td*v;
	double dd = norm2(D-C);

	if(dd < rr){
		double dt = sqrt(rr-dd);

		double da = rr-norm2(AC);
		double db = rr-norm2(BC);

		int inter = 0;
		double t1 = td-dt;
		double t2 = td+dt;
		if(t1 >= 0 && t1 <= l){
			P[inter] = A+t1*v;
			t[inter] = t1;
			inter++;
		}
		if(t2 >= 0 && t2 <= l){
			P[inter] = A+t2*v;
			t[inter] = t2;
			inter++;
		}

		return inter;
	}
	else{
		return 0;
	}
}




#endif // #ifndef _RX_SSM_H_