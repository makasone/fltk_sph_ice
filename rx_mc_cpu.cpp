/*!
  @file rx_mc_cpu.cpp
	
  @brief �A�֐��\�ʂ���̃|���S������(MC�@)
	
	http://local.wasp.uwa.edu.au/~pbourke/geometry/polygonise/
 
  @author Raghavendra Chandrashekara (basesd on source code
			provided by Paul Bourke and Cory Gene Bloyd)
  @date   2010-03
*/


//-----------------------------------------------------------------------------
// �C���N���[�h�t�@�C��
//-----------------------------------------------------------------------------
#include <math.h>
#include "rx_mc.h"
#include "rx_mc_tables.h"



//-----------------------------------------------------------------------------
// rxMCMeshCPU�̎���
//-----------------------------------------------------------------------------
/*!
 * �R���X�g���N�^
 */
rxMCMeshCPU::rxMCMeshCPU()
{
	// MARK:�R���X�g���N�^
	m_Grid.fMin = Vec3(0.0);
	m_Grid.fWidth = Vec3(0.0);
	m_Grid.iNum[0] = 0;
	m_Grid.iNum[1] = 0;
	m_Grid.iNum[2] = 0;

	m_nTriangles = 0;
	m_nNormals = 0;
	m_nVertices = 0;

	m_ptScalarField = NULL;
	m_fpScalarFunc = 0;
	m_tIsoLevel = 0;
	m_bValidSurface = false;
}

/*!
 * �f�X�g���N�^
 */
rxMCMeshCPU::~rxMCMeshCPU()
{
	DeleteSurface();
}


/*!
 * �A�֐�����O�p�`���b�V���𐶐�
 * @param[in] func �A�֐��l�擾�p�֐��|�C���^
 * @param[in] min_p �O���b�h�̍ŏ����W
 * @param[in] h �O���b�h�̕�
 * @param[in] n[3] �O���b�h��(x,y,z)
 * @param[in] threshold �������l(�A�֐��l�����̒l�̂Ƃ�������b�V����)
 * @param[in] method ���b�V���������@("mc", "rmt", "bloomenthal")
 * @param[out] vrts ���_���W
 * @param[out] nrms ���_�@��
 * @param[out] tris ���b�V��
 * @retval true  ���b�V����������
 * @retval false ���b�V���������s
 */
bool rxMCMeshCPU::CreateMesh(RXREAL (*func)(double, double, double), Vec3 min_p, double h, int n[3], RXREAL threshold, 
							 vector<Vec3> &vrts, vector<Vec3> &nrms, vector<rxFace> &face)
{
	if(func == NULL) return false;

	RxScalarField sf;
	for(int i = 0; i < 3; ++i){
		sf.iNum[i] = n[i];
		sf.fWidth[i] = h;
		sf.fMin[i] = min_p[i];
	}

	vector<int> tris;
	GenerateSurfaceV(sf, func, threshold, vrts, nrms, tris);

	if(IsSurfaceValid()){
		int nv = (int)GetNumVertices();
		int nm = (int)GetNumTriangles();
		int nn = (int)GetNumNormals();
		//cout << "mesh was created : " << nv << ", " << nm << ", " << nn << endl;

		// �@�����]
		for(int i = 0; i < nn; ++i){
			nrms[i] *= -1.0;
		}

		face.resize(nm);
		for(int i = 0; i < nm; ++i){
			face[i].vert_idx.resize(3);
			for(int j = 0; j < 3; ++j){
				face[i][j] = tris[3*i+(2-j)];
			}
		}

		return true;
	}

	return false;
}




/*!
 * �A�֐�����O�p�`���b�V���𐶐�
 * @param[in] field �T���v���{�����[��
 * @param[in] min_p �O���b�h�̍ŏ����W
 * @param[in] h �O���b�h�̕�
 * @param[in] n[3] �O���b�h��(x,y,z)
 * @param[in] threshold �������l(�A�֐��l�����̒l�̂Ƃ�������b�V����)
 * @param[in] method ���b�V���������@("mc", "rmt", "bloomenthal")
 * @param[out] vrts ���_���W
 * @param[out] nrms ���_�@��
 * @param[out] tris ���b�V��
 * @retval true  ���b�V����������
 * @retval false ���b�V���������s
 */
bool rxMCMeshCPU::CreateMeshV(RXREAL *field, Vec3 min_p, double h, int n[3], RXREAL threshold, 
								  vector<Vec3> &vrts, vector<Vec3> &nrms, vector<rxFace> &face)
{
	if(field == NULL) return false;

	RxScalarField sf;
	for(int i = 0; i < 3; ++i){
		sf.iNum[i] = n[i];
		sf.fWidth[i] = h;
		sf.fMin[i] = min_p[i];
	}

	vector<int> tris;
	GenerateSurface(sf, field, threshold, vrts, nrms, tris);

	if(IsSurfaceValid()){
		int nv = (int)GetNumVertices();
		int nm = (int)GetNumTriangles();
		int nn = (int)GetNumNormals();
		//cout << "mesh was created : " << nv << ", " << nm << ", " << nn << endl;

		// �@�����]
		for(int i = 0; i < nn; ++i){
			nrms[i] *= -1.0;
		}

		face.resize(nm);
		for(int i = 0; i < nm; ++i){
			face[i].vert_idx.resize(3);
			for(int j = 0; j < 3; ++j){
				face[i][j] = tris[3*i+(2-j)];
			}
		}

		return true;
	}

	return false;
}


/*!
 * ���b�V������
 * @param[in] sf �����O���b�h���
 * @param[in] field �T���v���{�����[��
 * @param[in] threshold 臒l
 * @param[out] vrts ���b�V�����_
 * @param[out] nrms ���b�V�����_�@��
 * @param[out] tris ���b�V���􉽏��(���_�ڑ����)
 */
void rxMCMeshCPU::GenerateSurface(const RxScalarField sf, RXREAL *field, RXREAL threshold, 
								  vector<Vec3> &vrts, vector<Vec3> &nrms, vector<int> &tris)
{
	// MARK:GenerateSurface
	if(m_bValidSurface){
		DeleteSurface();
	}

	m_tIsoLevel = threshold;
	m_Grid.iNum[0] = sf.iNum[0];
	m_Grid.iNum[1] = sf.iNum[1];
	m_Grid.iNum[2] = sf.iNum[2];
	m_Grid.fWidth = sf.fWidth;
	m_Grid.fMin = sf.fMin;
	m_ptScalarField = field;

	uint slice0 = (m_Grid.iNum[0] + 1);
	uint slice1 = slice0*(m_Grid.iNum[1] + 1);

	// ���l�ʂ̐���
	for(uint z = 0; z < m_Grid.iNum[2]; ++z){
		for(uint y = 0; y < m_Grid.iNum[1]; ++y){
			for(uint x = 0; x < m_Grid.iNum[0]; ++x){
				// �O���b�h���̒��_�z�u���e�[�u���Q�Ɨp�C���f�b�N�X�̌v�Z
				uint tableIndex = 0;
				if(m_ptScalarField[z*slice1 + y*slice0 + x] < m_tIsoLevel)
					tableIndex |= 1;
				if(m_ptScalarField[z*slice1 + (y+1)*slice0 + x] < m_tIsoLevel)
					tableIndex |= 2;
				if(m_ptScalarField[z*slice1 + (y+1)*slice0 + (x+1)] < m_tIsoLevel)
					tableIndex |= 4;
				if(m_ptScalarField[z*slice1 + y*slice0 + (x+1)] < m_tIsoLevel)
					tableIndex |= 8;
				if(m_ptScalarField[(z+1)*slice1 + y*slice0 + x] < m_tIsoLevel)
					tableIndex |= 16;
				if(m_ptScalarField[(z+1)*slice1 + (y+1)*slice0 + x] < m_tIsoLevel)
					tableIndex |= 32;
				if(m_ptScalarField[(z+1)*slice1 + (y+1)*slice0 + (x+1)] < m_tIsoLevel)
					tableIndex |= 64;
				if(m_ptScalarField[(z+1)*slice1 + y*slice0 + (x+1)] < m_tIsoLevel)
					tableIndex |= 128;

				if(edgeTable[tableIndex] != 0){
					// �G�b�W��̒��_�Z�o
					if(edgeTable[tableIndex] & 8){
						RxVertexID pt = CalculateIntersection(x, y, z, 3);
						uint id = GetEdgeID(x, y, z, 3);
						m_i2pt3idVertices.insert(ID2VertexID::value_type(id, pt));
					}
					if(edgeTable[tableIndex] & 1){
						RxVertexID pt = CalculateIntersection(x, y, z, 0);
						uint id = GetEdgeID(x, y, z, 0);
						m_i2pt3idVertices.insert(ID2VertexID::value_type(id, pt));
					}
					if(edgeTable[tableIndex] & 256){
						RxVertexID pt = CalculateIntersection(x, y, z, 8);
						uint id = GetEdgeID(x, y, z, 8);
						m_i2pt3idVertices.insert(ID2VertexID::value_type(id, pt));
					}
					
					if(x == m_Grid.iNum[0] - 1){
						if(edgeTable[tableIndex] & 4){
							RxVertexID pt = CalculateIntersection(x, y, z, 2);
							uint id = GetEdgeID(x, y, z, 2);
							m_i2pt3idVertices.insert(ID2VertexID::value_type(id, pt));
						}
						if(edgeTable[tableIndex] & 2048){
							RxVertexID pt = CalculateIntersection(x, y, z, 11);
							uint id = GetEdgeID(x, y, z, 11);
							m_i2pt3idVertices.insert(ID2VertexID::value_type(id, pt));
						}
					}
					if(y == m_Grid.iNum[1] - 1){
						if(edgeTable[tableIndex] & 2){
							RxVertexID pt = CalculateIntersection(x, y, z, 1);
							uint id = GetEdgeID(x, y, z, 1);
							m_i2pt3idVertices.insert(ID2VertexID::value_type(id, pt));
						}
						if(edgeTable[tableIndex] & 512){
							RxVertexID pt = CalculateIntersection(x, y, z, 9);
							uint id = GetEdgeID(x, y, z, 9);
							m_i2pt3idVertices.insert(ID2VertexID::value_type(id, pt));
						}
					}
					if(z == m_Grid.iNum[2] - 1){
						if(edgeTable[tableIndex] & 16){
							RxVertexID pt = CalculateIntersection(x, y, z, 4);
							uint id = GetEdgeID(x, y, z, 4);
							m_i2pt3idVertices.insert(ID2VertexID::value_type(id, pt));
						}
						if(edgeTable[tableIndex] & 128){
							RxVertexID pt = CalculateIntersection(x, y, z, 7);
							uint id = GetEdgeID(x, y, z, 7);
							m_i2pt3idVertices.insert(ID2VertexID::value_type(id, pt));
						}
					}
					if((x==m_Grid.iNum[0] - 1) && (y==m_Grid.iNum[1] - 1))
						if(edgeTable[tableIndex] & 1024){
							RxVertexID pt = CalculateIntersection(x, y, z, 10);
							uint id = GetEdgeID(x, y, z, 10);
							m_i2pt3idVertices.insert(ID2VertexID::value_type(id, pt));
						}
					if((x==m_Grid.iNum[0] - 1) && (z==m_Grid.iNum[2] - 1))
						if(edgeTable[tableIndex] & 64){
							RxVertexID pt = CalculateIntersection(x, y, z, 6);
							uint id = GetEdgeID(x, y, z, 6);
							m_i2pt3idVertices.insert(ID2VertexID::value_type(id, pt));
						}
					if((y==m_Grid.iNum[1] - 1) && (z==m_Grid.iNum[2] - 1))
						if(edgeTable[tableIndex] & 32){
							RxVertexID pt = CalculateIntersection(x, y, z, 5);
							uint id = GetEdgeID(x, y, z, 5);
							m_i2pt3idVertices.insert(ID2VertexID::value_type(id, pt));
						}
					
					// �|���S������
					for(uint i = 0; triTable[tableIndex][i] != 255; i += 3){
						RxTriangle triangle;
						uint pointID0, pointID1, pointID2;
						pointID0 = GetEdgeID(x, y, z, triTable[tableIndex][i]);
						pointID1 = GetEdgeID(x, y, z, triTable[tableIndex][i+1]);
						pointID2 = GetEdgeID(x, y, z, triTable[tableIndex][i+2]);
						triangle.pointID[0] = pointID0;
						triangle.pointID[1] = pointID1;
						triangle.pointID[2] = pointID2;
						m_trivecTriangles.push_back(triangle);
					}
				}
			}
		}
	}


	RenameVerticesAndTriangles(vrts, m_nVertices, tris, m_nTriangles);
	CalculateNormals(vrts, m_nVertices, tris, m_nTriangles, nrms, m_nNormals);

	m_bValidSurface = true;
}

/*!
 * ���b�V������(�T���v���{�����[���쐬)
 * @param[in] sf �����O���b�h���
 * @param[in] func �A�֐��l�擾�֐��|�C���^
 * @param[in] threshold 臒l
 * @param[out] vrts ���b�V�����_
 * @param[out] nrms ���b�V�����_�@��
 * @param[out] tris ���b�V���􉽏��(���_�ڑ����)
 */
void rxMCMeshCPU::GenerateSurfaceV(const RxScalarField sf, RXREAL (*func)(double, double, double), RXREAL threshold, 
										  vector<Vec3> &vrts, vector<Vec3> &nrms, vector<int> &tris)
{
	int nx, ny, nz;
	nx = sf.iNum[0]+1;
	ny = sf.iNum[1]+1;
	nz = sf.iNum[2]+1;

	Vec3 minp = sf.fMin;
	Vec3 d = sf.fWidth;

	RXREAL *field = new RXREAL[nx*ny*nz];
	for(int k = 0; k < nz; ++k){
		for(int j = 0; j < ny; ++j){
			for(int i = 0; i < nx; ++i){
				int idx = k*nx*ny+j*nx+i;
				Vec3 pos = minp+Vec3(i, j, k)*d;

				RXREAL val = func(pos[0], pos[1], pos[2]);
				field[idx] = val;
			}
		}
	}

	GenerateSurface(sf, field, threshold, vrts, nrms, tris);

	delete [] field;
}

/*!
 * ���b�V������(�֐�����)
 * @param[in] sf �����O���b�h���
 * @param[in] func �A�֐��l�擾�֐��|�C���^
 * @param[in] threshold 臒l
 * @param[out] vrts ���b�V�����_
 * @param[out] nrms ���b�V�����_�@��
 * @param[out] tris ���b�V���􉽏��(���_�ڑ����)
 */
void rxMCMeshCPU::GenerateSurfaceF(const RxScalarField sf, RXREAL (*func)(double, double, double), RXREAL threshold, 
										  vector<Vec3> &vrts, vector<Vec3> &nrms, vector<int> &tris)
{
	// MARK:GenerateSurfaceF
	if(m_bValidSurface){
		DeleteSurface();
	}

	m_tIsoLevel = threshold;
	m_Grid.iNum[0] = sf.iNum[0];
	m_Grid.iNum[1] = sf.iNum[1];
	m_Grid.iNum[2] = sf.iNum[2];
	m_Grid.fWidth = sf.fWidth;
	m_Grid.fMin = sf.fMin;
	m_fpScalarFunc = func;

	uint slice0 = (m_Grid.iNum[0] + 1);
	uint slice1 = slice0*(m_Grid.iNum[1] + 1);

	double dx = m_Grid.fWidth[0];
	double dy = m_Grid.fWidth[1];
	double dz = m_Grid.fWidth[2];
	
	for(uint k = 0; k < m_Grid.iNum[2]; ++k){
		for(uint j = 0; j < m_Grid.iNum[1]; ++j){
			for(uint i = 0; i < m_Grid.iNum[0]; ++i){
				double x, y, z;
				x = m_Grid.fMin[0]+i*m_Grid.fWidth[0];
				y = m_Grid.fMin[1]+j*m_Grid.fWidth[1];
				z = m_Grid.fMin[2]+k*m_Grid.fWidth[2];

				// �O���b�h���̒��_�z�u���e�[�u���Q�Ɨp�C���f�b�N�X�̌v�Z
				uint tableIndex = 0;
				if(m_fpScalarFunc(x,    y,    z) < m_tIsoLevel) tableIndex |= 1;
				if(m_fpScalarFunc(x,    y+dy, z) < m_tIsoLevel) tableIndex |= 2;
				if(m_fpScalarFunc(x+dx, y+dy, z) < m_tIsoLevel) tableIndex |= 4;
				if(m_fpScalarFunc(x+dx, y,    z) < m_tIsoLevel) tableIndex |= 8;
				if(m_fpScalarFunc(x,    y,    z+dz) < m_tIsoLevel) tableIndex |= 16;
				if(m_fpScalarFunc(x,    y+dy, z+dz) < m_tIsoLevel) tableIndex |= 32;
				if(m_fpScalarFunc(x+dx, y+dy, z+dz) < m_tIsoLevel) tableIndex |= 64;
				if(m_fpScalarFunc(x+dx, y,    z+dz) < m_tIsoLevel) tableIndex |= 128;

				// ���_���C�􉽏�񐶐�
				if(edgeTable[tableIndex] != 0){
					if(edgeTable[tableIndex] & 8){
						RxVertexID pt = CalculateIntersectionF(i, j, k, 3);
						uint id = GetEdgeID(i, j, k, 3);
						m_i2pt3idVertices.insert(ID2VertexID::value_type(id, pt));
					}
					if(edgeTable[tableIndex] & 1){
						RxVertexID pt = CalculateIntersectionF(i, j, k, 0);
						uint id = GetEdgeID(i, j, k, 0);
						m_i2pt3idVertices.insert(ID2VertexID::value_type(id, pt));
					}
					if(edgeTable[tableIndex] & 256){
						RxVertexID pt = CalculateIntersectionF(i, j, k, 8);
						uint id = GetEdgeID(i, j, k, 8);
						m_i2pt3idVertices.insert(ID2VertexID::value_type(id, pt));
					}
					
					if(i == m_Grid.iNum[0] - 1){
						if(edgeTable[tableIndex] & 4){
							RxVertexID pt = CalculateIntersectionF(i, j, k, 2);
							uint id = GetEdgeID(i, j, k, 2);
							m_i2pt3idVertices.insert(ID2VertexID::value_type(id, pt));
						}
						if(edgeTable[tableIndex] & 2048){
							RxVertexID pt = CalculateIntersectionF(i, j, k, 11);
							uint id = GetEdgeID(i, j, k, 11);
							m_i2pt3idVertices.insert(ID2VertexID::value_type(id, pt));
						}
					}
					if(j == m_Grid.iNum[1] - 1){
						if(edgeTable[tableIndex] & 2){
							RxVertexID pt = CalculateIntersectionF(i, j, k, 1);
							uint id = GetEdgeID(i, j, k, 1);
							m_i2pt3idVertices.insert(ID2VertexID::value_type(id, pt));
						}
						if(edgeTable[tableIndex] & 512){
							RxVertexID pt = CalculateIntersectionF(i, j, k, 9);
							uint id = GetEdgeID(i, j, k, 9);
							m_i2pt3idVertices.insert(ID2VertexID::value_type(id, pt));
						}
					}
					if(k == m_Grid.iNum[2] - 1){
						if(edgeTable[tableIndex] & 16){
							RxVertexID pt = CalculateIntersectionF(i, j, k, 4);
							uint id = GetEdgeID(i, j, k, 4);
							m_i2pt3idVertices.insert(ID2VertexID::value_type(id, pt));
						}
						if(edgeTable[tableIndex] & 128){
							RxVertexID pt = CalculateIntersectionF(i, j, k, 7);
							uint id = GetEdgeID(i, j, k, 7);
							m_i2pt3idVertices.insert(ID2VertexID::value_type(id, pt));
						}
					}
					if((i == m_Grid.iNum[0]-1) && (j == m_Grid.iNum[1]-1))
						if(edgeTable[tableIndex] & 1024){
							RxVertexID pt = CalculateIntersectionF(i, j, k, 10);
							uint id = GetEdgeID(i, j, k, 10);
							m_i2pt3idVertices.insert(ID2VertexID::value_type(id, pt));
						}
					if((i == m_Grid.iNum[0]-1) && (k == m_Grid.iNum[2]-1))
						if(edgeTable[tableIndex] & 64){
							RxVertexID pt = CalculateIntersectionF(i, j, k, 6);
							uint id = GetEdgeID(i, j, k, 6);
							m_i2pt3idVertices.insert(ID2VertexID::value_type(id, pt));
						}
					if((j == m_Grid.iNum[1]-1) && (k == m_Grid.iNum[2]-1))
						if(edgeTable[tableIndex] & 32){
							RxVertexID pt = CalculateIntersectionF(i, j, k, 5);
							uint id = GetEdgeID(i, j, k, 5);
							m_i2pt3idVertices.insert(ID2VertexID::value_type(id, pt));
						}
					
					for(uint t = 0; triTable[tableIndex][t] != 255; t += 3){
						RxTriangle triangle;
						uint pointID0, pointID1, pointID2;
						pointID0 = GetEdgeID(i, j, k, triTable[tableIndex][t]);
						pointID1 = GetEdgeID(i, j, k, triTable[tableIndex][t+1]);
						pointID2 = GetEdgeID(i, j, k, triTable[tableIndex][t+2]);
						triangle.pointID[0] = pointID0;
						triangle.pointID[1] = pointID1;
						triangle.pointID[2] = pointID2;
						m_trivecTriangles.push_back(triangle);
					}
				}
			}
		}
	}
	
	RenameVerticesAndTriangles(vrts, m_nVertices, tris, m_nTriangles);
	CalculateNormals(vrts, m_nVertices, tris, m_nTriangles, nrms, m_nNormals);
	m_bValidSurface = true;
}

/*!
 * ���l�ʃ��b�V���̔j��
 */
void rxMCMeshCPU::DeleteSurface()
{
	m_Grid.fWidth[0] = 0;
	m_Grid.fWidth[1] = 0;
	m_Grid.fWidth[2] = 0;
	m_Grid.iNum[0] = 0;
	m_Grid.iNum[1] = 0;
	m_Grid.iNum[2] = 0;

	m_nTriangles = 0;
	m_nNormals = 0;
	m_nVertices = 0;
	
	//m_vVertices.clear();
	//m_vNormals.clear();
	//m_vTriangles.clear();

	m_ptScalarField = NULL;
	m_tIsoLevel = 0;
	m_bValidSurface = false;
}

/*!
 * ���b�V�����ɗp�����O���b�h�̑傫��
 * @param[out] fVolLength* �O���b�h�̑傫��
 * @return ���b�V����������Ă����1, �����łȂ����-1
 */
int rxMCMeshCPU::GetVolumeLengths(double& fVolLengthX, double& fVolLengthY, double& fVolLengthZ)
{
	if(IsSurfaceValid()){
		fVolLengthX = m_Grid.fWidth[0]*m_Grid.iNum[0];
		fVolLengthY = m_Grid.fWidth[1]*m_Grid.iNum[1];
		fVolLengthZ = m_Grid.fWidth[2]*m_Grid.iNum[2];
		return 1;
	}
	else
		return -1;
}

/*!
 * �G�b�WID�̎擾
 * @param[in] nX,nY,nZ �O���b�h�ʒu
 * @param[in] nEdgeNo �G�b�W�ԍ�
 * @return �G�b�WID
 */
uint rxMCMeshCPU::GetEdgeID(uint nX, uint nY, uint nZ, uint nEdgeNo)
{
	switch(nEdgeNo){
	case 0:
		return GetVertexID(nX, nY, nZ) + 1;
	case 1:
		return GetVertexID(nX, nY + 1, nZ);
	case 2:
		return GetVertexID(nX + 1, nY, nZ) + 1;
	case 3:
		return GetVertexID(nX, nY, nZ);
	case 4:
		return GetVertexID(nX, nY, nZ + 1) + 1;
	case 5:
		return GetVertexID(nX, nY + 1, nZ + 1);
	case 6:
		return GetVertexID(nX + 1, nY, nZ + 1) + 1;
	case 7:
		return GetVertexID(nX, nY, nZ + 1);
	case 8:
		return GetVertexID(nX, nY, nZ) + 2;
	case 9:
		return GetVertexID(nX, nY + 1, nZ) + 2;
	case 10:
		return GetVertexID(nX + 1, nY + 1, nZ) + 2;
	case 11:
		return GetVertexID(nX + 1, nY, nZ) + 2;
	default:
		// Invalid edge no.
		return -1;
	}
}

/*!
 * ���_ID�̎擾
 * @param[in] nX,nY,nZ �O���b�h�ʒu
 * @return ���_ID
 */
uint rxMCMeshCPU::GetVertexID(uint nX, uint nY, uint nZ)
{
	return 3*(nZ*(m_Grid.iNum[1] + 1)*(m_Grid.iNum[0] + 1) + nY*(m_Grid.iNum[0] + 1) + nX);
}


/*!
 * ��Ԃɂ��G�b�W��̓��l�_���v�Z(�T���v���{�����[�����)
 * @param[in] nX,nY,nZ �O���b�h�ʒu
 * @param[in] nEdgeNo �G�b�W�ԍ�
 * @return ���b�V�����_���
 */
RxVertexID rxMCMeshCPU::CalculateIntersection(uint nX, uint nY, uint nZ, uint nEdgeNo)
{
	double x1, y1, z1, x2, y2, z2;
	uint v1x = nX, v1y = nY, v1z = nZ;
	uint v2x = nX, v2y = nY, v2z = nZ;
	
	switch(nEdgeNo){
	case 0:
		v2y += 1;
		break;
	case 1:
		v1y += 1;
		v2x += 1;
		v2y += 1;
		break;
	case 2:
		v1x += 1;
		v1y += 1;
		v2x += 1;
		break;
	case 3:
		v1x += 1;
		break;
	case 4:
		v1z += 1;
		v2y += 1;
		v2z += 1;
		break;
	case 5:
		v1y += 1;
		v1z += 1;
		v2x += 1;
		v2y += 1;
		v2z += 1;
		break;
	case 6:
		v1x += 1;
		v1y += 1;
		v1z += 1;
		v2x += 1;
		v2z += 1;
		break;
	case 7:
		v1x += 1;
		v1z += 1;
		v2z += 1;
		break;
	case 8:
		v2z += 1;
		break;
	case 9:
		v1y += 1;
		v2y += 1;
		v2z += 1;
		break;
	case 10:
		v1x += 1;
		v1y += 1;
		v2x += 1;
		v2y += 1;
		v2z += 1;
		break;
	case 11:
		v1x += 1;
		v2x += 1;
		v2z += 1;
		break;
	}

	x1 = m_Grid.fMin[0]+v1x*m_Grid.fWidth[0];
	y1 = m_Grid.fMin[1]+v1y*m_Grid.fWidth[1];
	z1 = m_Grid.fMin[2]+v1z*m_Grid.fWidth[2];
	x2 = m_Grid.fMin[0]+v2x*m_Grid.fWidth[0];
	y2 = m_Grid.fMin[1]+v2y*m_Grid.fWidth[1];
	z2 = m_Grid.fMin[2]+v2z*m_Grid.fWidth[2];

	uint slice0 = (m_Grid.iNum[0] + 1);
	uint slice1 = slice0*(m_Grid.iNum[1] + 1);
	RXREAL val1 = m_ptScalarField[v1z*slice1 + v1y*slice0 + v1x];
	RXREAL val2 = m_ptScalarField[v2z*slice1 + v2y*slice0 + v2x];
	RxVertexID intersection = Interpolate(x1, y1, z1, x2, y2, z2, val1, val2);
	
	return intersection;
}

/*!
 * ��Ԃɂ��G�b�W��̓��l�_���v�Z(�֐����)
 * @param[in] nX,nY,nZ �O���b�h�ʒu
 * @param[in] nEdgeNo �G�b�W�ԍ�
 * @return ���b�V�����_���
 */
RxVertexID rxMCMeshCPU::CalculateIntersectionF(uint nX, uint nY, uint nZ, uint nEdgeNo)
{
	double x1, y1, z1, x2, y2, z2;
	uint v1x = nX, v1y = nY, v1z = nZ;
	uint v2x = nX, v2y = nY, v2z = nZ;
	
	switch(nEdgeNo){
	case 0:
		v2y += 1;
		break;
	case 1:
		v1y += 1;
		v2x += 1;
		v2y += 1;
		break;
	case 2:
		v1x += 1;
		v1y += 1;
		v2x += 1;
		break;
	case 3:
		v1x += 1;
		break;
	case 4:
		v1z += 1;
		v2y += 1;
		v2z += 1;
		break;
	case 5:
		v1y += 1;
		v1z += 1;
		v2x += 1;
		v2y += 1;
		v2z += 1;
		break;
	case 6:
		v1x += 1;
		v1y += 1;
		v1z += 1;
		v2x += 1;
		v2z += 1;
		break;
	case 7:
		v1x += 1;
		v1z += 1;
		v2z += 1;
		break;
	case 8:
		v2z += 1;
		break;
	case 9:
		v1y += 1;
		v2y += 1;
		v2z += 1;
		break;
	case 10:
		v1x += 1;
		v1y += 1;
		v2x += 1;
		v2y += 1;
		v2z += 1;
		break;
	case 11:
		v1x += 1;
		v2x += 1;
		v2z += 1;
		break;
	}

	x1 = m_Grid.fMin[0]+v1x*m_Grid.fWidth[0];
	y1 = m_Grid.fMin[1]+v1y*m_Grid.fWidth[1];
	z1 = m_Grid.fMin[2]+v1z*m_Grid.fWidth[2];
	x2 = m_Grid.fMin[0]+v2x*m_Grid.fWidth[0];
	y2 = m_Grid.fMin[1]+v2y*m_Grid.fWidth[1];
	z2 = m_Grid.fMin[2]+v2z*m_Grid.fWidth[2];

	RXREAL val1 = m_fpScalarFunc(x1, y1, z1);
	RXREAL val2 = m_fpScalarFunc(x2, y2, z2);
	RxVertexID intersection = Interpolate(x1, y1, z1, x2, y2, z2, val1, val2);
	
	return intersection;
}

/*!
 * �O���b�h�G�b�W���[�̉A�֐��l������^��Ԃœ��l�_���v�Z
 * @param[in] fX1,fY1,fZ1 �[�_���W1
 * @param[in] fX2,fY2,fZ2 �[�_���W2
 * @param[in] tVal1 �[�_���W1�ł̃X�J���[�l
 * @param[in] tVal2 �[�_���W2�ł̃X�J���[�l
 * @return ���_���
 */
RxVertexID rxMCMeshCPU::Interpolate(double fX1, double fY1, double fZ1, double fX2, double fY2, double fZ2, RXREAL tVal1, RXREAL tVal2)
{
	RxVertexID interpolation;
	RXREAL mu;

	mu = RXREAL((m_tIsoLevel - tVal1))/(tVal2 - tVal1);
	interpolation.x = fX1 + mu*(fX2 - fX1);
	interpolation.y = fY1 + mu*(fY2 - fY1);
	interpolation.z = fZ1 + mu*(fZ2 - fZ1);

	return interpolation;
}


/*!
 * ���b�V�����_�C�􉽏����o�͌`���Ŋi�[
 * @param[out] vrts ���_���W
 * @param[out] nvrts ���_��
 * @param[out] tris �O�p�`�|���S���􉽏��
 * @param[out] ntris �O�p�`�|���S����
 */
void rxMCMeshCPU::RenameVerticesAndTriangles(vector<Vec3> &vrts, uint &nvrts, vector<int> &tris, uint &ntris)
{
	uint nextID = 0;
	ID2VertexID::iterator mapIterator = m_i2pt3idVertices.begin();
	RxTriangleVector::iterator vecIterator = m_trivecTriangles.begin();

	// Rename vertices.
	while(mapIterator != m_i2pt3idVertices.end()){
		(*mapIterator).second.newID = nextID;
		nextID++;
		mapIterator++;
	}

	// Now rename triangles.
	while(vecIterator != m_trivecTriangles.end()){
		for(uint i = 0; i < 3; i++){
			uint newID = m_i2pt3idVertices[(*vecIterator).pointID[i]].newID;
			(*vecIterator).pointID[i] = newID;
		}
		vecIterator++;
	}

	// Copy all the vertices and triangles into two arrays so that they
	// can be efficiently accessed.
	// Copy vertices.
	mapIterator = m_i2pt3idVertices.begin();
	nvrts = (int)m_i2pt3idVertices.size();
	vrts.resize(nvrts);
	for(uint i = 0; i < nvrts; i++, mapIterator++){
		vrts[i][0] = (*mapIterator).second.x;
		vrts[i][1] = (*mapIterator).second.y;
		vrts[i][2] = (*mapIterator).second.z;
	}
	// Copy vertex indices which make triangles.
	vecIterator = m_trivecTriangles.begin();
	ntris = (int)m_trivecTriangles.size();
	tris.resize(ntris*3);
	for(uint i = 0; i < ntris; i++, vecIterator++){
		tris[3*i+0] = (*vecIterator).pointID[0];
		tris[3*i+1] = (*vecIterator).pointID[1];
		tris[3*i+2] = (*vecIterator).pointID[2];
	}

	m_i2pt3idVertices.clear();
	m_trivecTriangles.clear();
}

/*!
 * ���_�@���v�Z
 * @param[in] vrts ���_���W
 * @param[in] nvrts ���_��
 * @param[in] tris �O�p�`�|���S���􉽏��
 * @param[in] ntris �O�p�`�|���S����
 * @param[out] nrms �@��
 * @param[out] nnrms �@����(=���_��)
 */
void rxMCMeshCPU::CalculateNormals(const vector<Vec3> &vrts, uint nvrts, const vector<int> &tris, uint ntris, 
									   vector<Vec3> &nrms, uint &nnrms)
{
	nnrms = nvrts;
	nrms.resize(nnrms);
	
	// Set all normals to 0.
	for(uint i = 0; i < nnrms; i++){
		nrms[i][0] = 0;
		nrms[i][1] = 0;
		nrms[i][2] = 0;
	}

	// Calculate normals.
	for(uint i = 0; i < ntris; i++){
		Vec3 vec1, vec2, normal;
		uint id0, id1, id2;
		id0 = tris[3*i+0];
		id1 = tris[3*i+1];
		id2 = tris[3*i+2];

		vec1 = vrts[id1]-vrts[id0];
		vec2 = vrts[id2]-vrts[id0];
		normal = cross(vec1, vec2);

		nrms[id0] += normal;
		nrms[id1] += normal;
		nrms[id2] += normal;
	}

	// Normalize normals.
	for(uint i = 0; i < nnrms; i++){
		normalize(nrms[i]);
	}
}

