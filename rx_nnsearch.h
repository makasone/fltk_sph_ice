/*!
  @file rx_nnsearch.h
	
  @brief ��`�O���b�h�����ɂ��ߖT�T��
 
  @author Makoto Fujisawa
  @date 2012-08
*/
// FILE --rx_nnsearch.h--

#ifndef _RX_NNSEARCH_H_
#define _RX_NNSEARCH_H_


//-----------------------------------------------------------------------------
// �C���N���[�h�t�@�C��
//-----------------------------------------------------------------------------
#include "rx_utility.h"

#include <algorithm>

#include <GL/glut.h>

#include "rx_pcube.h"


using namespace std;

typedef unsigned int uint;

#ifndef RXREAL
	#define RXREAL float
#endif


//-----------------------------------------------------------------------------
//! �n�b�V���l�ɂ��\�[�g�p�̍\����
//-----------------------------------------------------------------------------
struct rxHashSort
{
	uint hash;
	uint value;
};

/*!
 * �n�b�V���l�̔�r�֐�
 * @param[in] left,right ��r����l
 * @return left < right
 */
inline bool LessHash(const rxHashSort &left, const rxHashSort &right)
{
	return left.hash < right.hash;
}


//-----------------------------------------------------------------------------
//! �ߖT�T�����ʍ\����
//-----------------------------------------------------------------------------
struct rxNeigh
{
	int Idx;		//!< �ߖT�p�[�e�B�N���C���f�b�N�X
	RXREAL Dist;	//!< �ߖT�p�[�e�B�N���܂ł̋���
	RXREAL Dist2;	//!< �ߖT�p�[�e�B�N���܂ł�2�拗��
};


//-----------------------------------------------------------------------------
//! rxNNGrid�N���X - �O���b�h�����@�ɂ��ߖT�T��(3D)
//-----------------------------------------------------------------------------
class rxNNGrid
{
public:
	//! �T���p��ԕ����O���b�h�̊e�Z��
	struct rxCell
	{
		uint* hSortedIndex;			//!< �n�b�V���l�Ń\�[�g�����p�[�e�B�N���C���f�b�N�X
		uint* hGridParticleHash;	//!< �e�p�[�e�B�N���̃O���b�h�n�b�V���l
		uint* hCellStart;			//!< �\�[�g���X�g���̊e�Z���̃X�^�[�g�C���f�b�N�X
		uint* hCellEnd;				//!< �\�[�g���X�g���̊e�Z���̃G���h�C���f�b�N�X
		uint  uNumCells;			//!< ���Z����

		uint* hSortedPolyIdx;		//!< �n�b�V���l�Ń\�[�g�����|���S���C���f�b�N�X(�d���L��)
		uint* hGridPolyHash;		//!< �e�|���S���̃O���b�h�n�b�V���l
		uint* hPolyCellStart;		//!< �\�[�g���X�g���̊e�Z���̃X�^�[�g�C���f�b�N�X
		uint* hPolyCellEnd;			//!< �\�[�g���X�g���̊e�Z���̃G���h�C���f�b�N�X	

		uint uNumPolyHash;			//!< �|���S�����܂ރZ���̐�
	};


protected:
	// ��ԕ����i�q
	rxCell m_hCellData;				//!< �T���p��ԕ����O���b�h
	int m_iGridSize[3];				//!< �i�q�̐�
	double m_fCellWidth[3];			//!< �i�q��Ђ̒���

	int m_iDim;						//!< ������(���W�z��̃X�e�b�v��)
	Vec3 m_v3EnvMin;				//!< ���ŏ����W

public:
	//! �f�t�H���g�R���X�g���N�^
	rxNNGrid(int dim) : m_iDim(dim)
	{
		m_hCellData.uNumCells = 0;
		m_hCellData.hSortedIndex = 0;
		m_hCellData.hGridParticleHash = 0;
		m_hCellData.hCellStart = 0;
		m_hCellData.hCellEnd = 0;

		m_hCellData.hSortedPolyIdx = 0;
		m_hCellData.hGridPolyHash = 0;
		m_hCellData.hPolyCellStart = 0;
		m_hCellData.hPolyCellEnd = 0;
		m_hCellData.uNumPolyHash = 0;
	}

	//! �f�X�g���N�^
	~rxNNGrid()
	{
		if(m_hCellData.hSortedIndex) delete [] m_hCellData.hSortedIndex;
		if(m_hCellData.hGridParticleHash) delete [] m_hCellData.hGridParticleHash;
		if(m_hCellData.hCellStart) delete [] m_hCellData.hCellStart;
		if(m_hCellData.hCellEnd) delete [] m_hCellData.hCellEnd;

		if(m_hCellData.hSortedPolyIdx) delete [] m_hCellData.hSortedPolyIdx;
		if(m_hCellData.hGridPolyHash) delete [] m_hCellData.hGridPolyHash;
		if(m_hCellData.hPolyCellStart) delete [] m_hCellData.hPolyCellStart;
		if(m_hCellData.hPolyCellEnd) delete [] m_hCellData.hPolyCellEnd;
	}

public:
	// �����Z���̏����ݒ�
	void Setup(Vec3 vMin, Vec3 vMax, double h, int n);

	// �����Z���փp�[�e�B�N�����i�[
	void SetObjectToCell(RXREAL *p, uint n);

	// �����Z���փ|���S�����i�[
	void SetPolygonsToCell(RXREAL *vrts, int nv, int* tris, int nt);

	// �ߖT�擾
	void GetNN_Direct(Vec3 pos, RXREAL *p, uint n, vector<rxNeigh> &neighs, RXREAL h = -1.0);
	void GetNN(Vec3 pos, RXREAL *p, uint n, vector<rxNeigh> &neighs, RXREAL h = -1.0);

	// �Z�����̃|���S���擾
	int  GetPolygonsInCell(uint grid_hash, vector<int> &polys);
	int  GetPolygonsInCell(int gi, int gj, int gk, vector<int> &polys);
	bool IsPolygonsInCell(int gi, int gj, int gk);

	// OpenGL�`��
	void DrawCell(int i, int j, int k);
	void DrawCells(Vec3 col, Vec3 col2, int sel = 0, RXREAL *p = 0);

	// �O���b�h�n�b�V���̌v�Z
	uint CalGridHash(int x, int y, int z);
	uint CalGridHash(Vec3 pos);

public:
	rxCell& GetCellData(void){ return m_hCellData; }

protected:
	// �����Z������ߖT�p�[�e�B�N�����擾
	void getNeighborsInCell(Vec3 pos, RXREAL *p, int gi, int gj, int gk, vector<rxNeigh> &neighs, RXREAL h);
};


	
/*!
 * ��ԕ����@�̏���
 * @param[in] vMin ���̍ŏ����W
 * @param[in] vMax ���̍ő���W
 * @param[in] h �e�����a
 */
inline void rxNNGrid::Setup(Vec3 vMin, Vec3 vMax, double h, int n)
{
	if(h < RX_EPS) return;

	Vec3 world_size = vMax-vMin;
	Vec3 world_origin = vMin;

	double max_axis = RXFunc::Max3(world_size);

	int d = (int)(log(max_axis/h)/log(2.0)+0.5);
	int m = (int)(pow(2.0, (double)d)+0.5);
	double cell_width = max_axis/m;

	d = (int)(log(world_size[0]/cell_width)/log(2.0)+0.5);
	m_iGridSize[0] = (int)(pow(2.0, (double)d)+0.5);
	d = (int)(log(world_size[1]/cell_width)/log(2.0)+0.5);
	m_iGridSize[1] = (int)(pow(2.0, (double)d)+0.5);;
	d = (int)(log(world_size[2]/cell_width)/log(2.0)+0.5);
	m_iGridSize[2] = (int)(pow(2.0, (double)d)+0.5);;

	m_fCellWidth[0] = cell_width;
	m_fCellWidth[1] = cell_width;
	m_fCellWidth[2] = cell_width;

	m_v3EnvMin = world_origin;

	m_hCellData.uNumCells = m_iGridSize[0]*m_iGridSize[1]*m_iGridSize[2];

	cout << "grid for nn search : " << endl;
	cout << "  size   : " << m_iGridSize[0] << "x" << m_iGridSize[1] << endl;
	cout << "  num    : " << m_hCellData.uNumCells << endl;
	cout << "  origin : " << m_v3EnvMin << endl;
	cout << "  width  : " << m_fCellWidth[0] << "x" << m_fCellWidth[1] << endl;

	// �����O���b�h�\���̂̔z��m��
	m_hCellData.hSortedIndex = new uint[n];
	m_hCellData.hGridParticleHash = new uint[n];
	m_hCellData.hCellStart = new uint[m_hCellData.uNumCells];
	m_hCellData.hCellEnd = new uint[m_hCellData.uNumCells];

	m_hCellData.hPolyCellStart = new uint[m_hCellData.uNumCells];
	m_hCellData.hPolyCellEnd = new uint[m_hCellData.uNumCells];

}

/*!
 * �p�[�e�B�N���𕪊��Z���Ɋi�[
 *  - �p�[�e�B�N���̑�����O���b�h�n�b�V�����v�Z���Ċi�[����
 * @param[in] p �i�[�������S�p�[�e�B�N���̍��W���L�q�����z��
 * @param[in] n �p�[�e�B�N����
 */
inline void rxNNGrid::SetObjectToCell(RXREAL *p, uint n)
{
	int mem_size1 = n*sizeof(uint);
	int mem_size2 = m_hCellData.uNumCells*sizeof(uint);
	memset(m_hCellData.hSortedIndex, 0, mem_size1);
	memset(m_hCellData.hGridParticleHash, 0, mem_size1);
	memset(m_hCellData.hCellStart, 0xffffffff, mem_size2);
	memset(m_hCellData.hCellEnd, 0xffffffff, mem_size2);

	if(n == 0) return;

	// �e�p�[�e�B�N���̃O���b�h�n�b�V���̌v�Z
	for(uint i = 0; i < n; ++i){
		Vec3 pos;
		pos[0] = p[m_iDim*i+0];
		pos[1] = p[m_iDim*i+1];
		pos[2] = p[m_iDim*i+2];

		// �n�b�V���l�v�Z
		uint hash = CalGridHash(pos);

		m_hCellData.hSortedIndex[i] = i;
		m_hCellData.hGridParticleHash[i] = hash;
	}

	// �O���b�h�n�b�V���Ń\�[�g
	vector<rxHashSort> hash_and_value;
	hash_and_value.resize(n);
	for(uint i = 0; i < n; ++i){
		hash_and_value[i].hash = m_hCellData.hGridParticleHash[i];
		hash_and_value[i].value  = m_hCellData.hSortedIndex[i];
	}
	std::sort(hash_and_value.begin(), hash_and_value.end(), LessHash);
	for(uint i = 0; i < n; ++i){
		m_hCellData.hSortedIndex[i] = hash_and_value[i].value;
		m_hCellData.hGridParticleHash[i] = hash_and_value[i].hash;
	}

	// �p�[�e�B�N���z����\�[�g���ꂽ���Ԃɕ��ёւ��C
	// �e�Z���̎n�܂�ƏI���̃C���f�b�N�X������
	for(uint i = 0; i < n; i++){
		int hash = m_hCellData.hGridParticleHash[i];

		if(i == 0){
			m_hCellData.hCellStart[hash] = i;
			m_hCellData.hCellEnd[hash] = i;
		}
		else{
			int prev_hash = m_hCellData.hGridParticleHash[i-1];

			if(i == 0 || hash != prev_hash){
				m_hCellData.hCellStart[hash] = i;
				if(i > 0){
					m_hCellData.hCellEnd[prev_hash] = i;
				}
			}

			if(i == n-1){
				m_hCellData.hCellEnd[hash] = i+1;
			}
		}
	}
}


/*!
 * �|���S���𕪊��Z���Ɋi�[
 * @param[in] vrts �|���S�����_
 * @param[in] nv ���_��
 * @param[in] tris ���b�V��
 * @param[in] nt ���b�V����
 */
inline void rxNNGrid::SetPolygonsToCell(RXREAL *vrts, int nv, int* tris, int nt)
{
	int mem_size2 = m_hCellData.uNumCells*sizeof(uint);
	memset(m_hCellData.hPolyCellStart, 0xffffffff, mem_size2);
	memset(m_hCellData.hPolyCellEnd, 0, mem_size2);

	int num_hash = 0;

	// �e�|���S���̃O���b�h�n�b�V���̌v�Z
	vector<uint> tri_hash, tri_idx;
	vector<Vec3> tri_vrts, tri_vrts_c;
	tri_vrts.resize(3);
	tri_vrts_c.resize(3);
	for(int i = 0; i < nt; i++){
		for(int j = 0; j < 3; ++j){
			Vec3 pos;
			pos[0] = vrts[3*tris[3*i+j]+0];
			pos[1] = vrts[3*tris[3*i+j]+1];
			pos[2] = vrts[3*tris[3*i+j]+2];
			tri_vrts[j] = pos;
		}

		Vec3 nrm = Unit(cross(tri_vrts[1]-tri_vrts[0], tri_vrts[2]-tri_vrts[0]));

		// �|���S����BBox
		Vec3 bmin, bmax;
		bmin = tri_vrts[0];
		bmax = tri_vrts[0];
		for(int j = 1; j < 3; ++j){
			for(int k = 0; k < 3; ++k){
				if(tri_vrts[j][k] < bmin[k]) bmin[k] = tri_vrts[j][k];
				if(tri_vrts[j][k] > bmax[k]) bmax[k] = tri_vrts[j][k];
			}
		}

		// BBox�Əd�Ȃ�Z��
		bmin -= m_v3EnvMin;
		bmax -= m_v3EnvMin;

		// �����Z���C���f�b�N�X�̎Z�o
		int bmin_gidx[3], bmax_gidx[3];
		for(int k = 0; k < 3; ++k){
			bmin_gidx[k] = bmin[k]/m_fCellWidth[k];
			bmax_gidx[k] = bmax[k]/m_fCellWidth[k];

			bmin_gidx[k] = RX_CLAMP(bmin_gidx[k], 0, m_iGridSize[k]-1);
			bmax_gidx[k] = RX_CLAMP(bmax_gidx[k], 0, m_iGridSize[k]-1);
		}

		// �e�Z���Ƀ|���S�����܂܂�邩���`�F�b�N
		Vec3 len = Vec3(m_fCellWidth[0], m_fCellWidth[1], m_fCellWidth[2]);
		Vec3 cen(0.0);
		for(int x = bmin_gidx[0]; x <= bmax_gidx[0]; ++x){
			for(int y = bmin_gidx[1]; y <= bmax_gidx[1]; ++y){
				for(int z = bmin_gidx[2]; z <= bmax_gidx[2]; ++z){
					cen = m_v3EnvMin+Vec3(x+0.5, y+0.5, z+0.5)*len;

					for(int j = 0; j < 3; ++j){
						tri_vrts_c[j] = (tri_vrts[j]-cen)/len;
					}

					if(RXFunc::polygon_intersects_cube(tri_vrts_c, nrm, 0)){
						// �n�b�V���l�v�Z
						uint hash = CalGridHash(x, y, z);

						tri_idx.push_back((uint)i);
						tri_hash.push_back(hash);

						//m_hCellData.hPolyCellStart[hash] = 0;

						num_hash++;
					}
				}
			}
		}
	}

	m_hCellData.uNumPolyHash = (uint)num_hash;

	int mem_size1 = m_hCellData.uNumPolyHash*sizeof(uint);
	m_hCellData.hSortedPolyIdx = new uint[m_hCellData.uNumPolyHash];
	m_hCellData.hGridPolyHash = new uint[m_hCellData.uNumPolyHash];
	memcpy(m_hCellData.hSortedPolyIdx, &tri_idx[0], mem_size1);
	memcpy(m_hCellData.hGridPolyHash, &tri_hash[0], mem_size1);

	// �O���b�h�n�b�V���Ń\�[�g
	vector<rxHashSort> hash_and_value;
	hash_and_value.resize(m_hCellData.uNumPolyHash);
	for(int i = 0; i < (int)m_hCellData.uNumPolyHash; ++i){
		hash_and_value[i].hash  = m_hCellData.hGridPolyHash[i];
		hash_and_value[i].value = m_hCellData.hSortedPolyIdx[i];
	}
	std::sort(hash_and_value.begin(), hash_and_value.end(), LessHash);
	for(int i = 0; i < (int)m_hCellData.uNumPolyHash; ++i){
		m_hCellData.hSortedPolyIdx[i] = hash_and_value[i].value;
		m_hCellData.hGridPolyHash[i]  = hash_and_value[i].hash;
	}

	// �p�[�e�B�N���z����\�[�g���ꂽ���Ԃɕ��ёւ��C
	// �e�Z���̎n�܂�ƏI���̃C���f�b�N�X������
	for(uint i = 0; i < m_hCellData.uNumPolyHash; ++i){
		uint hash = m_hCellData.hGridPolyHash[i];

		if(i == 0){
			m_hCellData.hPolyCellStart[hash] = i;
		}
		else{
			uint prev_hash = m_hCellData.hGridPolyHash[i-1];

			if(i == 0 || hash != prev_hash){
				m_hCellData.hPolyCellStart[hash] = i;
				if(i > 0){
					m_hCellData.hPolyCellEnd[prev_hash] = i;
				}
			}

			if(i == nt-1){
				m_hCellData.hPolyCellEnd[hash] = i+1;
			}
		}
	}
}

/*!
 * �ߖT���q�T��(��������)
 * @param[in] pos �T�����S
 * @param[in] p �p�[�e�B�N���ʒu
 * @param[out] neighs �T�����ʊi�[����ߖT���R���e�i
 * @param[in] h �L�����a
 */
inline void rxNNGrid::GetNN_Direct(Vec3 pos0, RXREAL *p, uint n, vector<rxNeigh> &neighs, RXREAL h)
{
	RXREAL h2 = h*h;

	for(uint i = 0; i < n; i++){
		Vec3 pos1;
		pos1[0] = p[m_iDim*i+0];
		pos1[1] = p[m_iDim*i+1];
		pos1[2] = p[m_iDim*i+2];

		rxNeigh neigh;
		neigh.Dist2 = (RXREAL)norm2(pos0-pos1);

		if(neigh.Dist2 <= h2 && neigh.Dist2 > RX_FEQ_EPS){
			neigh.Idx = i;
			neighs.push_back(neigh);
		}
	}
}

/*!
 * �ߖT���q�T��
 * @param[in] pos �T�����S
 * @param[in] p �p�[�e�B�N���ʒu
 * @param[out] neighs �T�����ʊi�[����ߖT���R���e�i
 * @param[in] h �L�����a
 */
inline void rxNNGrid::GetNN(Vec3 pos, RXREAL *p, uint n, vector<rxNeigh> &neighs, RXREAL h)
{
	// �����Z���C���f�b�N�X�̎Z�o
	int x = (pos[0]-m_v3EnvMin[0])/m_fCellWidth[0];
	int y = (pos[1]-m_v3EnvMin[1])/m_fCellWidth[1];
	int z = (pos[2]-m_v3EnvMin[2])/m_fCellWidth[2];

	int numArdGrid = (int)(h/m_fCellWidth[0])+1;
	for(int k = -numArdGrid; k <= numArdGrid; ++k){
		for(int j = -numArdGrid; j <= numArdGrid; ++j){
			for(int i = -numArdGrid; i <= numArdGrid; ++i){
				int i1 = x+i;
				int j1 = y+j;
				int k1 = z+k;
				if(i1 < 0 || i1 >= m_iGridSize[0] || j1 < 0 || j1 >= m_iGridSize[1] || k1 < 0 || k1 >= m_iGridSize[2]){
					continue;
				}

				getNeighborsInCell(pos, p, i1, j1, k1, neighs, h);
			}
		}
	}
}


/*!
 * �����Z�����̗��q����ߖT�����o
 * @param[in] pos �T�����S
 * @param[in] p �p�[�e�B�N���ʒu
 * @param[in] gi,gj,gk �Ώە����Z��
 * @param[out] neighs �T�����ʊi�[����ߖT���R���e�i
 * @param[in] h �L�����a
 */
inline void rxNNGrid::getNeighborsInCell(Vec3 pos, RXREAL *p, int gi, int gj, int gk, vector<rxNeigh> &neighs, RXREAL h)
{
	RXREAL h2 = h*h;

	uint grid_hash = CalGridHash(gi, gj, gk);

	uint start_index = m_hCellData.hCellStart[grid_hash];
	if(start_index != 0xffffffff){	// �Z������łȂ����̃`�F�b�N
		uint end_index = m_hCellData.hCellEnd[grid_hash];
		for(uint j = start_index; j < end_index; ++j){
			uint idx = m_hCellData.hSortedIndex[j];

			Vec3 xij;
			xij[0] = pos[0]-p[m_iDim*idx+0];
			xij[1] = pos[1]-p[m_iDim*idx+1];
			xij[2] = pos[2]-p[m_iDim*idx+2];

			rxNeigh neigh;
			neigh.Dist2 = norm2(xij);

			if(neigh.Dist2 <= h2){
				neigh.Idx = idx;
				//neigh.Dist = sqrt(neigh.Dist2);

				neighs.push_back(neigh);
			}
		}
	}
}

/*!
 * �����Z���Ɋi�[���ꂽ�|���S�������擾
 * @param[in] gi,gj,gk �Ώە����Z��
 * @param[out] polys �|���S��
 * @return �i�[�|���S����
 */
inline int rxNNGrid::GetPolygonsInCell(uint grid_hash, vector<int> &polys)
{
	uint start_index = m_hCellData.hPolyCellStart[grid_hash];
	if(start_index != 0xffffffff){	// �Z������łȂ����̃`�F�b�N

		int cnt = 0;
		uint end_index = m_hCellData.hPolyCellEnd[grid_hash];
		for(uint j = start_index; j < end_index; ++j){
			uint idx = m_hCellData.hSortedPolyIdx[j];

			polys.push_back(idx);

			cnt++;
		}

		return cnt;
	}

	return 0;
}
inline int rxNNGrid::GetPolygonsInCell(int gi, int gj, int gk, vector<int> &polys)
{
	return GetPolygonsInCell(CalGridHash(gi, gj, gk), polys);
}
/*!
 * �����Z�����̃|���S���̗L���𒲂ׂ�
 * @param[in] gi,gj,gk �Ώە����Z��
 * @return �|���S�����i�[����Ă����true
 */
inline bool rxNNGrid::IsPolygonsInCell(int gi, int gj, int gk)
{
	uint grid_hash = CalGridHash(gi, gj, gk);

	uint start_index = m_hCellData.hPolyCellStart[grid_hash];
	if(start_index != 0xffffffff){	// �Z������łȂ����̃`�F�b�N

		int cnt = 0;
		uint end_index = m_hCellData.hPolyCellEnd[grid_hash];
		for(uint j = start_index; j < end_index; ++j){
			uint idx = m_hCellData.hSortedPolyIdx[j];
			cnt++;
			break;
		}

		return (cnt > 0);
	}

	return false;
}

/*!
 * �O���b�h�n�b�V���l�̌v�Z
 * @param[in] x,y,z �O���b�h�ʒu
 * @return �O���b�h�n�b�V���l
 */
inline uint rxNNGrid::CalGridHash(int x, int y, int z)
{
	x = (x < 0 ? 0 : (x >= m_iGridSize[0] ? m_iGridSize[0]-1 : x));
	y = (y < 0 ? 0 : (y >= m_iGridSize[1] ? m_iGridSize[1]-1 : y));
	z = (z < 0 ? 0 : (z >= m_iGridSize[2] ? m_iGridSize[2]-1 : z));
	return z*m_iGridSize[1]*m_iGridSize[0]+y*m_iGridSize[0]+x;
}
/*!
 * �O���b�h�n�b�V���l�̌v�Z
 * @param[in] pos �p�[�e�B�N�����W
 * @return �O���b�h�n�b�V���l
 */
inline uint rxNNGrid::CalGridHash(Vec3 pos)
{
	pos -= m_v3EnvMin;

	// �����Z���C���f�b�N�X�̎Z�o
	int x = pos[0]/m_fCellWidth[0];
	int y = pos[1]/m_fCellWidth[1];
	int z = pos[2]/m_fCellWidth[2];
	return CalGridHash(x, y, z);
}


//-----------------------------------------------------------------------------
// OpenGL�`��
//-----------------------------------------------------------------------------
/*!
 * �T���p�Z���̕`��
 * @param[in] i,j,k �O���b�h��̃C���f�b�N�X
 */
inline void rxNNGrid::DrawCell(int i, int j, int k)
{
	glPushMatrix();
	glTranslated(m_v3EnvMin[0], m_v3EnvMin[1], m_v3EnvMin[2]);
	glTranslatef((i+0.5)*m_fCellWidth[0], (j+0.5)*m_fCellWidth[1], (k+0.5)*m_fCellWidth[2]);
	glutWireCube(m_fCellWidth[0]);
	glPopMatrix();
}

/*!
 * �T���p�O���b�h�̕`��
 * @param[in] col �p�[�e�B�N�����܂܂��Z���̐F
 * @param[in] col2 �|���S�����܂܂��Z���̐F
 * @param[in] sel �����_���ɑI�����ꂽ�Z���̂ݕ`��(1�ŐV�����Z����I���C2�ł��łɑI������Ă���Z����`��C0�ł��ׂẴZ����`��)
 * @param[in] p �p�[�e�B�N���ʒu
 */
inline void rxNNGrid::DrawCells(Vec3 col, Vec3 col2, int sel, RXREAL *p)
{
	glPushMatrix();

	if(sel){
		// �����_���ɑI�񂾃Z���Ƃ��̒��̃p�[�e�B�N���̂ݕ`��
		static int grid_hash = 0;
		static uint start_index = 0xffffffff;
		if(sel == 1){
			do{
				grid_hash = RXFunc::Nrand(m_hCellData.uNumCells-1);
				start_index = m_hCellData.hCellStart[grid_hash];
			}while(start_index == 0xffffffff);
		}

		uint w = grid_hash%(m_iGridSize[0]*m_iGridSize[1]);
		DrawCell(w%m_iGridSize[0], w/m_iGridSize[0], grid_hash/(m_iGridSize[0]*m_iGridSize[1]));

		if(p){
			glColor3d(1.0, 0.0, 0.0);
			glPointSize(10.0);
			glBegin(GL_POINTS);

			int c = 0;
			uint end_index = m_hCellData.hCellEnd[grid_hash];
			for(uint j = start_index; j < end_index; ++j){
				uint idx = m_hCellData.hSortedIndex[j];
				Vec3 pos;
				pos[0] = p[m_iDim*idx+0];
				pos[1] = p[m_iDim*idx+1];
				pos[2] = p[m_iDim*idx+2];
			
				glVertex3dv(pos);

				c++;
			}
			glEnd();

			cout << "cell(" << grid_hash << ") : " << c << endl;
		}

	}
	else{
		// �p�[�e�B�N�� or �|���S�����܂ޑS�Z���̕`��
		RXFOR3(0, m_iGridSize[0], 0, m_iGridSize[1], 0, m_iGridSize[2]){
			bool disp = false;
			uint grid_hash = CalGridHash(i, j, k);
			uint start_index = m_hCellData.hCellStart[grid_hash];
			uint start_index_poly = 0xffffffff;

			if(m_hCellData.uNumPolyHash) start_index_poly = m_hCellData.hPolyCellStart[grid_hash];
		
			if(start_index != 0xffffffff){
				glColor3dv(col2.data);
				disp = true;
			}
			if(start_index_poly != 0xffffffff){
				glColor3dv(col.data);
				disp = true;
			}

			if(disp){
				DrawCell(i, j, k);
			}
		}
	}

	glPopMatrix();
}

#endif // #ifndef _RX_NNSEARCH_H_