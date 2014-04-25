/*!
  @file rx_ssm_cpu.cpp
	
  @brief Screen Space Mesh�쐬
		- M. Muller et al., Screen space meshes, SCA2007, 2007. 

  @author Makoto Fujisawa
  @date 2011-05
*/
// FILE --rx_ssm_cpu.cpp--


//-----------------------------------------------------------------------------
// �C���N���[�h�t�@�C��
//-----------------------------------------------------------------------------
#include "rx_ssm.h"


//-----------------------------------------------------------------------------
// �萔�E�ϐ�
//-----------------------------------------------------------------------------
const int BITS[4] = { 1, 2, 4, 8 };


//-----------------------------------------------------------------------------
// rxSSMeshCPU�N���X�̎���
//-----------------------------------------------------------------------------
/*!
 * �R���X�g���N�^
 */
rxSSMeshCPU::rxSSMeshCPU(double zmax, double h, double r, int n_filter, int n_iter)
	: rxSSMesh(zmax, h, r, n_filter, n_iter)
{
	// MARK:�R���X�g���N�^
	m_vFilter = CalBinomials(21);

	// ���b�V���\���e�[�u���p�C���f�b�N�X�X�V�֐��|�C���^ ev*5+nv
	// �G�b�W���_��0
	m_FuncTableIndex[0]  = 0;
	m_FuncTableIndex[1]  = 0;
	m_FuncTableIndex[2]  = 0;
	m_FuncTableIndex[3]  = 0;
	m_FuncTableIndex[4]  = &rxSSMeshCPU::updateTableIndexE0N4;
	// �G�b�W���_��1
	m_FuncTableIndex[5]  = &rxSSMeshCPU::updateTableIndexE1;
	m_FuncTableIndex[6]  = &rxSSMeshCPU::updateTableIndexE1;
	m_FuncTableIndex[7]  = &rxSSMeshCPU::updateTableIndexE1;
	m_FuncTableIndex[8]  = &rxSSMeshCPU::updateTableIndexE1;
	m_FuncTableIndex[9]  = &rxSSMeshCPU::updateTableIndexE1;
	// �G�b�W���_��2
	m_FuncTableIndex[10] = 0;
	m_FuncTableIndex[11] = 0;
	m_FuncTableIndex[12] = 0;
	m_FuncTableIndex[13] = 0;
	m_FuncTableIndex[14] = &rxSSMeshCPU::updateTableIndexE2N4;
	// �G�b�W���_��3
	m_FuncTableIndex[15] = 0;
	m_FuncTableIndex[16] = 0;
	m_FuncTableIndex[17] = &rxSSMeshCPU::updateTableIndexE3N23;
	m_FuncTableIndex[18] = &rxSSMeshCPU::updateTableIndexE3N23;
	m_FuncTableIndex[19] = &rxSSMeshCPU::updateTableIndexE3N4;
	// �G�b�W���_��4
	m_FuncTableIndex[20] = 0;
	m_FuncTableIndex[21] = 0;
	m_FuncTableIndex[22] = &rxSSMeshCPU::updateTableIndexE4N2;
	m_FuncTableIndex[23] = &rxSSMeshCPU::updateTableIndexE4N3;
	m_FuncTableIndex[24] = &rxSSMeshCPU::updateTableIndexE4N4;
}


/*!
 * �f�X�g���N�^
 */
rxSSMeshCPU::~rxSSMeshCPU()
{
}

/*!
 * �}�b�v��O���b�h�z��̃T�C�Y��ύX
 * @param[in] W,H ��ʉ𑜓x
 * @param[in] spacing �f�v�X�}�b�v�̃T���v�����O�Ԋu
 */
void rxSSMeshCPU::Resize(int W, int H, int spacing)
{
	// MARK:Resize
	if(spacing != -1){
		m_fSpacing = spacing;
	}

	m_iNgx = W/m_fSpacing;
	m_iNgy = H/m_fSpacing;

	m_fDmx = (double)W/(double)(m_iNgx);
	m_fDmy = (double)H/(double)(m_iNgy);

	//
	// CPU�������m��
	//
	// �f�v�X�}�b�v
	if((int)m_vSSDMap.size() != (m_iNgx+1)*(m_iNgy+1)){
		m_vSSDMap.clear();
		m_vSSDMap.resize((m_iNgx+1)*(m_iNgy+1));
	}

	// ���b�V�������p�O���b�h
	if((int)m_vSSMGrid.size() != m_iNgx*m_iNgy){
		m_vSSMGrid.clear();
		m_vSSMGrid.resize(m_iNgx*m_iNgy);
	}
}

/*!
 * �X�N���[���X�y�[�X���b�V������(�@���v�Z�܂�)
 *  - M. Muller et al., Screen space meshes, SCA2007, 2007. 
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
void rxSSMeshCPU::CreateMesh(double *proj, double *modelview, int W, int H, vector<Vec3> &prts, int pnum, 
							 vector<Vec3> &vrts, vector<Vec3> &nrms, vector< vector<int> > &polys, 
							 int filtering, int debug_output)
{
	CreateMesh(proj, modelview, W, H, prts, pnum, vrts, polys, filtering, debug_output);

	nrms.clear();
	CalVertexNormals(vrts, vrts.size(), polys, polys.size(), nrms);
}

/*!
 * �X�N���[���X�y�[�X���b�V������
 *  - M. Muller et al., Screen space meshes, SCA2007, 2007. 
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
void rxSSMeshCPU::CreateMesh(double *proj, double *modelview, int W, int H, vector<Vec3> &prts, int pnum, 
							 vector<Vec3> &vrts, vector< vector<int> > &polys, int filtering, int debug_output)
{
	// MARK:CreateMesh
	int nx = W/m_fSpacing;
	int ny = H/m_fSpacing;

	if(nx != m_iNgx || ny != m_iNgy || m_vSSDMap.empty()){
		Resize(W, H, m_fSpacing);
	}

	if(debug_output){
		cout << "window : " << W << " x " << H << endl;
		cout << "depth map : " << m_iNgx+1 << " x " << m_iNgy+1 << endl;
		cout << "grid : " << m_iNgx << " x " << m_iNgy << endl;
	}

	//
	// �R���e�i������
	//

	// Screen Space�ł̃p�[�e�B�N��
	if((int)m_vSSPrts.size() != pnum){
		m_vSSPrts.clear();
		m_vSSPrts.resize(pnum);
	}

	// �֊s�G�b�W
	m_vSSEdge.clear();

	// Screen Space�ł̃��b�V�����_
	m_vSSVertex.clear();

	for(int i = 0; i < m_iNgx; ++i){
		for(int j = 0; j < m_iNgy; ++j){
			int idx = i+j*m_iNgx;
			for(int k = 0; k < 4; ++k){
				m_vSSMGrid[idx].node_vrts[k] = -1;
				m_vSSMGrid[idx].edge_vrts[k] = -1;
				m_vSSMGrid[idx].back_vrts[k] = -1;
			}
			m_vSSMGrid[idx].back_vrts[4] = -1;
			m_vSSMGrid[idx].back_vrts[5] = -1;
			m_vSSMGrid[idx].num_nv = 0;
			m_vSSMGrid[idx].num_ev = 0;
			m_vSSMGrid[idx].num_bv = 0;

			m_vSSMGrid[idx].i = i;
			m_vSSMGrid[idx].j = j;
		}
	}

	// ���b�V��������
	vrts.clear();
	polys.clear();
	m_vMeshGrid.clear();


	
	// �������e�ϊ��s��
	const rxMatrix4 P = GetMatrixGL(proj);

	// ���f���r���[�ϊ��s��
	const rxMatrix4 MV = GetMatrixGL(modelview);


	//
	// �p�[�e�B�N�����W�Ɣ��a�𓧎����e�ϊ����ăf�v�X�}�b�v�𐶐�
	//
	// �f�v�X�}�b�v������	
	for(int i = 0; i < (m_iNgx+1)*(m_iNgy+1); ++i){
		m_vSSDMap[i] = RX_FEQ_INF;
	}

	// �f�v�X�}�b�v�쐬
	CalDepthMap(P, MV, W, H, prts, pnum);



	// 
	// �f�v�X�}�b�v�̕�����
	//
	if(filtering & 0x01){
		ApplyDepthFilter(m_vSSDMap, m_iNgx+1, m_iNgy+1, m_iNfilter);
	}
	



	//
	// �m�[�h���_����
	//
	m_iNumNodeVrts = CalNodeVertex(m_iNgx, m_iNgy, m_fDmx, m_fDmy, m_vSSDMap);
	if(debug_output) cout << "the number of node vertices = " << m_iNumNodeVrts << endl;


	//
	// �֊s�G�b�W�̌��o�ƃG�b�W���_����
	//
	int edge_num = DetectSilhouetteEdgeVertex(m_iNgx, m_iNgy, m_fDmx, m_fDmy, m_vSSDMap, m_vSSPrts, W, H);
	if(debug_output) cout << "the number of silhouette edges = " << edge_num << endl;
	//if(debug_output) cout << "the number of edges vertices 2 = " << m_vSSEdgeVertex.size() << endl;
	
	if(debug_output){
		int evn = 0;
		for(int i = 0; i < (int)m_vSSEdge.size(); ++i){
			if(m_vSSEdge[i].front_vertex != -1){
				evn++;
			}
		}
		cout << "the number of front edge vertices : " << evn << endl;
	}
	
	// �֊s�G�b�W���back vertex
	m_iNumEdgeVrts = CalEdgeVertex(m_iNgx, m_iNgy, m_fDmx, m_fDmy, m_vSSDMap, m_vSSEdge);
	if(debug_output) cout << "the number of edge vertices = " << m_iNumEdgeVrts << endl;



	// 
	// ���b�V������
	// 
	m_iNumMesh = CalMesh(m_iNgx, m_iNgy, m_vSSMGrid, polys, 0);
	if(debug_output) cout << "the number of mesh = " << m_iNumMesh << endl;



	// 
	// �֊s�̕�����
	//
	if(filtering & 0x02){
		ApplySilhoutteSmoothing(m_vSSVertex, polys, m_iNiters);
	}



	//
	// 3D��ԂɃ��b�V����߂�
	//
	vrts.clear();

	// �������e�ϊ��t�s��
	rxMatrix4 Q = P.Inverse();

	// ���f���r���[�ϊ��t�s��
	rxMatrix4 IMV = MV.Inverse();

	rxMatrix4 IMVQ = IMV*Q;

	int nv = (int)m_vSSVertex.size();
	vrts.resize(nv);
	for(int i = 0; i < nv; ++i){
		Vec3 xp = m_vSSVertex[i].pos;
		Vec4 xd;
		xd[0] = -1.0+2.0*xp[0]/W;
		xd[1] = -1.0+2.0*xp[1]/H;
		xd[3] = (1.0-Q(3,2)*xp[2])/(Q(3,0)*xd[0]+Q(3,1)*xd[1]+Q(3,3));

		xd[0] *= xd[3];
		xd[1] *= xd[3];
		xd[2] = xp[2];

		Vec4 x = IMVQ*xd;
	
		vrts[i] = Vec3(x[0], x[1], x[2]);
	}

	if(debug_output) cout << "the number of vertices = " << vrts.size() << endl;

}

/*!
 * �p�[�e�B�N�����W�Ɣ��a�𓧎����e�ϊ����ăf�v�X�}�b�v�𐶐�
 * @param[in] P �������e�s��
 * @param[in] MV ���f���r���[�s��
 * @param[in] W,H ��ʉ𑜓x
 * @param[in] prts �p�[�e�B�N�����W
 * @param[in] pnum �p�[�e�B�N����
 */
void rxSSMeshCPU::CalDepthMap(const rxMatrix4 &P, const rxMatrix4 &MV, int W, int H, vector<Vec3> &prts, int pnum)
{
	// CalDepthMap
	rxMatrix4 PMV = P*MV;

	Vec3 tr;
	tr[0] = m_fPrtRad*W*sqrt(P(0, 0)*P(0, 0)+P(0, 1)*P(0, 1)+P(0, 2)*P(0, 2));
	tr[1] = m_fPrtRad*H*sqrt(P(1, 0)*P(1, 0)+P(1, 1)*P(1, 1)+P(1, 2)*P(1, 2));
	tr[2] = m_fPrtRad*sqrt(P(2, 0)*P(2, 0)+P(2, 1)*P(2, 1)+P(2, 2)*P(2, 2));

	for(int k = 0; k < pnum; ++k){
		Vec4 x = Vec4(prts[k], 1.0);

		// ���e�ϊ�
		Vec4 xd = PMV*x;

		// w�Ŋ��邱�Ƃ�[-1, 1]�̐��K�����W�n�ɕϊ�
		Vec3 xp;
		xp[0] = W*(0.5+0.5*xd[0]/xd[3]);
		xp[1] = H*(0.5+0.5*xd[1]/xd[3]);
		xp[2] = xd[2];		// z�����͌��̒l��p����

		m_vSSPrts[k].xp = xp;

		Vec4 rd = Vec4(m_fPrtRad, m_fPrtRad, m_fPrtRad, 1.0);
		rd = PMV*rd;

		// ���K�����W�n�ł̔��a�l
		Vec3 rp;
		rp[0] = tr[0]/xd[3]/2;
		rp[1] = tr[1]/xd[3]/2;
		rp[2] = tr[2];

		m_vSSPrts[k].rp = rp;

		double rrp = rp[0]*rp[0];

		if(xp[0] < 0 || xp[0] >= W || xp[1] < 0 || xp[1] >= H){
			continue;
		}

		// �f�v�X�}�b�v��ł̃p�[�e�B�N���͈̔�
		int cen[2];	// �p�[�e�B�N�����S
		cen[0] = int(xp[0]/m_fDmx)+1;
		cen[1] = int(xp[1]/m_fDmy)+1;

		int minp[2], maxp[2];
		minp[0] = cen[0]-(rp[0]/m_fDmx+2);
		minp[1] = cen[1]-(rp[1]/m_fDmy+2);
		maxp[0] = cen[0]+(rp[0]/m_fDmx+2);
		maxp[1] = cen[1]+(rp[1]/m_fDmy+2);

		// �͈͂��}�b�v�O�ɂȂ�Ȃ��悤�ɃN�����v
		RX_CLAMP2(minp[0], 0, m_iNgx);
		RX_CLAMP2(minp[1], 0, m_iNgy);
		RX_CLAMP2(maxp[0], 0, m_iNgx);
		RX_CLAMP2(maxp[1], 0, m_iNgy);

		// �p�[�e�B�N���f�v�X�l�X�V
		for(int j = minp[1]; j <= maxp[1]; ++j){
			for(int i = minp[0]; i <= maxp[0]; ++i){
				double rr = (i*m_fDmx-xp[0])*(i*m_fDmx-xp[0])+(j*m_fDmy-xp[1])*(j*m_fDmy-xp[1]);
				if(rr <= rrp){
					double zij = m_vSSDMap[i+j*(m_iNgx+1)];
					double hij = sqrt(1.0-rr/rrp);

					double z = xp[2]-rp[2]*hij;
					if(z >= 0 && z < zij){
						m_vSSDMap[i+j*(m_iNgx+1)] = z;
					}
				}
			}
		}
	}
	//RXTIMER("ssm(depth map1)");
}

/*!
 * �f�v�X�}�b�v�ɕ��������{��
 * @param[inout] dmap �f�v�X�}�b�v
 * @param[in] nx,ny �}�b�v�𑜓x
 * @param[in] n_filter �t�B���^�����O��
 */
void rxSSMeshCPU::ApplyDepthFilter(vector<RXREAL> &dmap, int nx, int ny, int n_filter)
{
	int b = 2*n_filter+1;

	// binomial�W���̌v�Z
	if((int)m_vFilter.size() != b){
		m_vFilter = CalBinomials(b);
	}

	vector<RXREAL> tmp_map;
	tmp_map = dmap;

	vector<RXREAL> d;
	d.resize(b);
	
	// ���������̃t�B���^�[��K�p
	for(int j = 0; j < ny; ++j){
		for(int i = n_filter; i < nx-n_filter; ++i){
			d[0] = tmp_map[i+j*nx];
			if(d[0] < RX_FEQ_INF-1){	// != ���̃s�N�Z���݂̂ɓK�p
				// ����Nfilter�O���b�h�̒l���擾
				int n = 0;
				for(int k = 0; k < n_filter; ++k){
					d[2*n+1] = tmp_map[(i-(k+1))+j*nx];
					d[2*n+2] = tmp_map[(i+(k+1))+j*nx];
					if(fabs(d[2*n+1]-d[0]) > m_fSSZmax || fabs(d[2*n+2]-d[0]) > m_fSSZmax){
						break;
					}
					n++;
				}

				int bn = 2*n;

				RXREAL new_depth = m_vFilter[bn][n]*d[0];
				for(int k = 1; k <= n; ++k){
					new_depth += m_vFilter[bn][n-k]*d[2*k-1];
					new_depth += m_vFilter[bn][n+k]*d[2*k];
				}

				dmap[i+j*nx] = new_depth;
			}
		}
	}

	tmp_map = dmap;

	// ���������̃t�B���^�[��K�p
	for(int i = 0; i < nx; ++i){
		for(int j = n_filter; j < ny-n_filter; ++j){
			d[0] = tmp_map[i+j*nx];
			if(d[0] < RX_FEQ_INF-1){	// != ���̃s�N�Z���݂̂ɓK�p
				// ����Nfilter�O���b�h�̒l���擾
				int n = 0;
				for(int k = 0; k < n_filter; ++k){
					d[2*n+1] = tmp_map[i+(j-(k+1))*nx];
					d[2*n+2] = tmp_map[i+(j+(k+1))*nx];
					if(fabs(d[2*n+1]-d[0]) > m_fSSZmax || fabs(d[2*n+2]-d[0]) > m_fSSZmax){
						break;
					}
					n++;
				}

				int bn = 2*n;

				RXREAL new_depth = m_vFilter[bn][n]*d[0];
				for(int k = 1; k <= n; ++k){
					new_depth += m_vFilter[bn][n-k]*d[2*k-1];
					new_depth += m_vFilter[bn][n+k]*d[2*k];
				}

				dmap[i+j*nx] = new_depth;
			}
		}
	}

}

/*!
 * �֊s�ɕ��������{��
 * @param[inout] ssvrts �X�N���[�����W�n�ł̒��_��
 * @param[in] polys ���b�V��(�\�����钸�_��)
 * @param[in] n_iter �t�B���^�����O������
 */
void rxSSMeshCPU::ApplySilhoutteSmoothing(vector<rxSSVertex> &ssvrts, const vector< vector<int> > &polys, int n_iter)
{
	// ApplySilhoutteSmoothing
	int vn = (int)ssvrts.size();
	int pn = (int)polys.size();

	for(int l = 0; l < n_iter; ++l){
		for(int i = 0; i < vn; ++i){
			ssvrts[i].avg_pos = Vec3(0.0);
			ssvrts[i].avg_num = 0.0;
		}

		int neigh[3][2] = { {1, 2}, {2, 0}, {0, 1} };

		for(int i = 0; i < pn; ++i){
			for(int j = 0; j < 3; ++j){
				int idx = polys[i][j];
				if(ssvrts[idx].edge){
					if(ssvrts[idx].avg_num == 0){
						ssvrts[idx].avg_pos += ssvrts[idx].pos;
						ssvrts[idx].avg_num += 1.0;
					}

					// �אڒ��_
					for(int k = 0; k < 2; ++k){
						int nidx = polys[i][neigh[j][k]];
						if(ssvrts[nidx].edge){
							ssvrts[idx].avg_pos += ssvrts[nidx].pos;
							ssvrts[idx].avg_num += 1.0;
						}
						else{
							ssvrts[idx].avg_pos += 0.5*ssvrts[nidx].pos;
							ssvrts[idx].avg_num += 0.5;
						}
					}
				}
			}
		}

		for(int i = 0; i < vn; ++i){
			if(ssvrts[i].edge && ssvrts[i].avg_num){
				ssvrts[i].avg_pos /= ssvrts[i].avg_num;
				ssvrts[i].pos[0] = ssvrts[i].avg_pos[0];
				ssvrts[i].pos[1] = ssvrts[i].avg_pos[1];
				ssvrts[i].pos[2] = ssvrts[i].avg_pos[2];
			}
		}
	}
}

/*!
 * �֊s�G�b�W�̌��o��front edge vertex�̌v�Z
 * @param[in] nx,ny ���b�V�������p�O���b�h�𑜓x
 * @param[in] dw,dh ���b�V�������p�O���b�h��
 * @param[in] dgrid ���b�V�������p�f�v�X�}�b�v
 * @param[in] ssprts Screen Space�ł̃p�[�e�B�N��
 * @param[in] W,H �X�N���[���𑜓x
 * @return ���o���ꂽ�֊s�G�b�W��
 */
int rxSSMeshCPU::DetectSilhouetteEdgeVertex(int nx, int ny, double dw, double dh, const vector<RXREAL> &dgrid, 
										 const vector<rxSSParticle> &ssprts, int W, int H)
{
	// DetectSilhouetteEdgeVertex
	int en = 0;
	int pnum = (int)ssprts.size();

	m_vSSEdge.resize((nx)*(ny+1)+(nx+1)*(ny));
	int yoffset = (nx)*(ny+1);
	
	// x�����G�b�W������
	for(int j = 0; j <= ny; ++j){
		for(int i = 0; i < nx; ++i){
			int i00 = i+j*(nx+1);
			int i01 = (i+1)+j*(nx+1);

			rxSSEdge *e = &m_vSSEdge[i+j*(nx)];

			e->d0 = dgrid[i00];
			e->d1 = dgrid[i01];
			e->depth = 0.5*(e->d0+e->d1);
			e->x0 = Vec3(dw*i, dh*j, 0.0);
			e->x1 = Vec3(dw*(i+1), dh*j, 0.0);
			e->xy = 0;
			e->silhouette = 0;
			e->front_vertex = -1;
			e->dx = 0.0;

			if(fabs(e->d0-e->d1) > m_fSSZmax){
				e->silhouette = 1;
				en++;
			}
		}
	}

	// y�����G�b�W�̃`�F�b�N
	for(int j = 0; j < ny; ++j){
		for(int i = 0; i <= nx; ++i){
			int i00 = i+j*(nx+1);
			int i10 = i+(j+1)*(nx+1);

			rxSSEdge *e = &m_vSSEdge[yoffset+i+j*(nx+1)];

			e->d0 = dgrid[i00];
			e->d1 = dgrid[i10];
			e->depth = 0.5*(e->d0+e->d1);
			e->x0 = Vec3(dw*i, dh*j, 0.0);
			e->x1 = Vec3(dw*i, dh*(j+1), 0.0);
			e->xy = 1;
			e->silhouette = 0;
			e->front_vertex = -1;
			e->dx = 0.0;

			if(fabs(e->d0-e->d1) > m_fSSZmax){
				e->silhouette = 1;
				en++;
			}
		}
	}

	m_vSSEdgeVertex.clear();
	for(int k = 0; k < pnum; ++k){
		Vec3 xp = ssprts[k].xp;
		Vec3 rp = ssprts[k].rp;

		if(xp[0] < 0 || xp[0] >= W || xp[1] < 0 || xp[1] >= H){
			continue;
		}

		// �f�v�X�}�b�v��ł̃p�[�e�B�N���͈̔�
		int cen[2];	// �p�[�e�B�N�����S
		cen[0] = int(xp[0]/dw)+1;
		cen[1] = int(xp[1]/dh)+1;

		int minp[2], maxp[2];
		minp[0] = cen[0]-(rp[0]/dw+1);
		minp[1] = cen[1]-(rp[1]/dh+1);
		maxp[0] = cen[0]+(rp[0]/dw+1);
		maxp[1] = cen[1]+(rp[1]/dh+1);

		// �͈͂��}�b�v�O�ɂȂ�Ȃ��悤�ɃN�����v
		RX_CLAMP2(minp[0], 0, nx-1);
		RX_CLAMP2(minp[1], 0, ny-1);
		RX_CLAMP2(maxp[0], 0, nx-1);
		RX_CLAMP2(maxp[1], 0, ny-1);

		// �͈͓��̃G�b�W�𒲍�(x����)
		for(int j = minp[1]; j <= maxp[1]+1; ++j){
			for(int i = minp[0]; i <= maxp[0]; ++i){
				rxSSEdge &e = m_vSSEdge[i+j*(nx)];
				if(!e.silhouette) continue;
				if(e.depth < xp[2]) continue;

				// �~�ƃG�b�W�̌�_
				Vec2 A, B, C, P[2];
				double r, t[2];
				if(e.d0 <= e.d1){
					A = Vec2(e.x0[0], e.x0[1]);
					B = Vec2(e.x1[0], e.x1[1]);
				}
				else{
					A = Vec2(e.x1[0], e.x1[1]);
					B = Vec2(e.x0[0], e.x0[1]);
				}
				C = Vec2(xp[0], xp[1]);
				r = rp[0];
				int inter = LineCircleIntersection(A, B, C, r, P, t);
				if(inter == 1){
					if(e.front_vertex == -1 || t[0] > e.dx){
						Vec3 vrt_pos;
						vrt_pos[0] = P[0][0];
						vrt_pos[1] = P[0][1];
						vrt_pos[2] = xp[2];
						m_vSSEdgeVertex.push_back(rxSSVertex(vrt_pos, 1));
						int eidx = (int)m_vSSEdgeVertex.size()-1;

						e.front_vertex = eidx;
						e.dx = t[0];
					}
				}

			}
		}
		// �͈͓��̃G�b�W�𒲍�(x����)
		for(int i = minp[0]; i <= maxp[0]+1; ++i){
			for(int j = minp[1]; j <= maxp[1]; ++j){
				rxSSEdge &e = m_vSSEdge[yoffset+i+j*(nx+1)];
				if(!e.silhouette) continue;
				if(e.depth < xp[2]) continue;

				// �~�ƃG�b�W�̌�_
				Vec2 A, B, C, P[2];
				double r, t[2];
				if(e.d0 <= e.d1){
					A = Vec2(e.x0[0], e.x0[1]);
					B = Vec2(e.x1[0], e.x1[1]);
				}
				else{
					A = Vec2(e.x1[0], e.x1[1]);
					B = Vec2(e.x0[0], e.x0[1]);
				}
				C = Vec2(xp[0], xp[1]);
				r = rp[0];
				int inter = LineCircleIntersection(A, B, C, r, P, t);
				if(inter == 1){
					if(e.front_vertex == -1 || t[0] > e.dx){
						Vec3 vrt_pos;
						vrt_pos[0] = P[0][0];
						vrt_pos[1] = P[0][1];
						vrt_pos[2] = xp[2];
						m_vSSEdgeVertex.push_back(rxSSVertex(vrt_pos, 1));
						int eidx = (int)m_vSSEdgeVertex.size()-1;

						e.front_vertex = eidx;
						e.dx = t[0];
					}
				}

			}
		}
	}

	return en;
}

/*!
 * �֊s�G�b�W�̌��o
 * @param[in] nx,ny ���b�V�������p�O���b�h�𑜓x
 * @param[in] dw,dh ���b�V�������p�O���b�h��
 * @param[in] dgrid ���b�V�������p�f�v�X�}�b�v
 * @return ���o���ꂽ�֊s�G�b�W��
 */
int rxSSMeshCPU::DetectSilhouetteEdge(int nx, int ny, double dw, double dh, const vector<RXREAL> &dgrid)
{
	int en = 0;

	// x�����G�b�W�̃`�F�b�N
	for(int j = 0; j <= ny; ++j){
		for(int i = 0; i < nx; ++i){
			int i00 = i+j*(nx+1);
			int i01 = (i+1)+j*(nx+1);
			double d00 = dgrid[i00];
			double d01 = dgrid[i01];

			if(fabs(d00-d01) > m_fSSZmax){
				rxSSEdge e;
				e.depth = (d00 < d01 ? d00 : d01);
				e.d0 = d00;
				e.d1 = d01;
				e.x0 = Vec3(dw*i, dh*j, 0.0);
				e.x1 = Vec3(dw*(i+1), dh*j, 0.0);
				e.xy = 0;
				e.front_vertex = -1;
				e.dx = 0.0;
				e.silhouette = 1;
				m_vSSEdge.push_back(e);

				en++;
			}
		}
	}

	// y�����G�b�W�̃`�F�b�N
	for(int j = 0; j < ny; ++j){
		for(int i = 0; i <= nx; ++i){
			int i00 = i+j*(nx+1);
			int i10 = i+(j+1)*(nx+1);
			double d00 = dgrid[i00];
			double d10 = dgrid[i10];

			if(fabs(d00-d10) > m_fSSZmax){
				rxSSEdge e;
				e.depth = (d00 < d10 ? d00 : d10);
				e.d0 = d00;
				e.d1 = d10;
				e.x0 = Vec3(dw*i, dh*j, 0.0);
				e.x1 = Vec3(dw*i, dh*(j+1), 0.0);
				e.xy = 1;
				e.front_vertex = -1;
				e.dx = 0.0;
				e.silhouette = 1;
				m_vSSEdge.push_back(e);

				en++;
			}
		}
	}

	return en;
}


/*!
 * �m�[�h���_����
 * @param[in] nx,ny ���b�V�������p�O���b�h�𑜓x
 * @param[in] dw,dh ���b�V�������p�O���b�h��
 * @param[in] dgrid ���b�V�������p�f�v�X�}�b�v
 * @return �������ꂽ�m�[�h���_��
 */
int rxSSMeshCPU::CalNodeVertex(int nx, int ny, double dw, double dh, const vector<RXREAL> &dgrid)
{
	// �O���b�h�m�[�h��̒��_
	int nv_num = 0;
	for(int j = 0; j <= ny; ++j){
		for(int i = 0; i <= nx; ++i){
			double d = dgrid[i+j*(nx+1)];
			if(d < RX_FEQ_INF-1){
				Vec3 vrt_pos = Vec3(dw*i, dh*j, d);

				m_vSSVertex.push_back(rxSSVertex(vrt_pos, 0));
				nv_num++;
				int vidx = (int)m_vSSVertex.size()-1;

				// ���_�����i�[
				if(i != 0){
					if(j != 0){
						m_vSSMGrid[(i-1)+(j-1)*nx].node_vrts[2] = vidx;
						m_vSSMGrid[(i-1)+(j-1)*nx].num_nv++;
					}
					if(j != ny){
						m_vSSMGrid[(i-1)+(j)*nx].node_vrts[1] = vidx;
						m_vSSMGrid[(i-1)+(j)*nx].num_nv++;
					}
				}
				if(i != nx){
					if(j != 0){
						m_vSSMGrid[(i)+(j-1)*nx].node_vrts[3] = vidx;
						m_vSSMGrid[(i)+(j-1)*nx].num_nv++;
					}
					if(j != ny){
						m_vSSMGrid[(i)+(j)*nx].node_vrts[0] = vidx;
						m_vSSMGrid[(i)+(j)*nx].num_nv++;
					}
				}
			}

			// �f�v�X�l�̊i�[
			if(i != 0){
				if(j != 0){
					m_vSSMGrid[(i-1)+(j-1)*nx].node_depth[2] = d;
				}
				if(j != ny){
					m_vSSMGrid[(i-1)+(j)*nx].node_depth[1] = d;
				}
			}
			if(i != nx){
				if(j != 0){
					m_vSSMGrid[(i)+(j-1)*nx].node_depth[3] = d;
				}
				if(j != ny){
					m_vSSMGrid[(i)+(j)*nx].node_depth[0] = d;
				}
			}

		}
	}

	return nv_num;
}

/*!
 * �G�b�W���_����
 * @param[in] nx,ny ���b�V�������p�O���b�h�𑜓x
 * @param[in] dw,dh ���b�V�������p�O���b�h��
 * @param[in] dgrid ���b�V�������p�f�v�X�}�b�v
 * @param[in] edges �֊s�G�b�W���X�g
 * @return �������ꂽ�G�b�W���_��
 */
int rxSSMeshCPU::CalEdgeVertex(int nx, int ny, double dw, double dh, const vector<RXREAL> &dgrid, const vector<rxSSEdge> &edges)
{
	int ev_num = 0;
	int yoffset = (nx)*(ny+1);
	for(int i = 0; i < (int)edges.size(); ++i){
		rxSSEdge e = edges[i];
		if(!e.silhouette) continue;

		int ei0, ei1, ej0, ej1;
		if(i < yoffset){	// x�����s�G�b�W
			ei0 = i%nx;
			ej0 = i/nx;
			ei1 = ei0+1;
			ej1 = ej0;
		}
		else{
			ei0 = (i-yoffset)%(nx+1);
			ej0 = (i-yoffset)/(nx+1);
			ei1 = ei0;
			ej1 = ej0+1;
		}

		Vec3 vrt_pos(0.0);
		double dx = 0.0;

		if(e.front_vertex != -1){
			vrt_pos = m_vSSEdgeVertex[e.front_vertex].pos;
			dx = e.dx;
		}
		else{
			if(e.d0 <= e.d1){
				dx = binarySearchDepth(e.x0, e.x1, vrt_pos, m_fSSZmax);
			}
			else{
				dx = binarySearchDepth(e.x1, e.x0, vrt_pos, m_fSSZmax);
			}

			if(vrt_pos[0] == 0 && vrt_pos[1] == 0){
				cout << "edge vertex error : " << i << endl;
			}
		}

		m_vSSVertex.push_back(rxSSVertex(vrt_pos, 1));
		ev_num++;
		int vidx = (int)m_vSSVertex.size()-1;


		// back vertex�����邩�ǂ������`�F�b�N
		int bidx = -1;
		Vec3 back_vrt;
		if(e.d0 < RX_FEQ_INFM && e.d1 < RX_FEQ_INFM){	// �G�b�W�[�_�������Ƃ� != �� �Ȃ�� back vertex ������
			RXREAL back_node_depth, nn_back_node_depth, back_depth;
			int back_node, nn_back_node;	// back layer�ɑ�����G�b�W�[�_�Ƃ��̗אڃm�[�h
			double l = 1.0;	// �G�b�W�̒���

			// �f�v�X�l���傫������ back layer �ɑ�����
			if(e.d0 > e.d1){
				back_node = ei0+ej0*(nx+1);
				if(i < yoffset){	// x�����s�G�b�W
					nn_back_node = (ei0 == 0 ? ei0 : ei0-1)+ej0*(nx+1);
					l = dw;
				}
				else{			// y�����s�G�b�W
					nn_back_node = ei0+(ej0 == 0 ? ej0 : ej0-1)*(nx+1);
					l = dh;
				}
			}
			else{
				back_node = ei1+ej1*(nx+1);
				if(i < yoffset){	// x�����s�G�b�W
					nn_back_node = (ei1 == nx ? ei1 : ei1+1)+ej1*(nx+1);
					l = dw;
				}
				else{			// y�����s�G�b�W
					nn_back_node = ei1+(ej1 == ny ? ej1 : ej1+1)*(nx+1);
					l = dh;
				}
			}

			// back layer�ɑ�����G�b�W�[�_�Ƃ��̗אڃm�[�h�̃f�v�X�l
			back_node_depth = dgrid[back_node];
			nn_back_node_depth = dgrid[nn_back_node];

			// �f�v�X�l���O�}�ɂ��ߎ�
			back_depth = back_node_depth;
			//back_depth = back_node_depth*((2*l-dx)/l)-nn_back_node_depth*((l-dx)/l);

			// back vertex��ݒ�
			back_vrt[0] = vrt_pos[0];
			back_vrt[1] = vrt_pos[1];
			back_vrt[2] = back_depth;

			m_vSSVertex.push_back(rxSSVertex(back_vrt, 1));
			ev_num++;
			bidx = (int)m_vSSVertex.size()-1;
		}

		//if(bidx != -1 && back_vrt[2] < vrt_pos[2]){
		//	int tmp = vidx;
		//	vidx = bidx;
		//	bidx = tmp;
		//}

		if(e.xy == 0){	// x�����s�G�b�W
			if(ej0 != 0){	// �����E�ł͂Ȃ�
				m_vSSMGrid[ei0+(ej0-1)*nx].edge_vrts[2] = vidx;
				m_vSSMGrid[ei0+(ej0-1)*nx].num_ev++;
				if(bidx != -1){
					m_vSSMGrid[ei0+(ej0-1)*nx].back_vrts[2] = bidx;
					m_vSSMGrid[ei0+(ej0-1)*nx].num_bv++;
				}
			}
			if(ej0 != ny){	// �㋫�E�ł͂Ȃ�
				m_vSSMGrid[ei0+ej0*nx].edge_vrts[0] = vidx;
				m_vSSMGrid[ei0+ej0*nx].num_ev++;
				if(bidx != -1){
					m_vSSMGrid[ei0+ej0*nx].back_vrts[0] = bidx;
					m_vSSMGrid[ei0+ej0*nx].num_bv++;
				}
			}
		}
		else{	// y�����s�G�b�W
			if(ei0 != 0){
				m_vSSMGrid[(ei0-1)+ej0*nx].edge_vrts[1] = vidx;
				m_vSSMGrid[(ei0-1)+ej0*nx].num_ev++;
				if(bidx != -1){
					m_vSSMGrid[(ei0-1)+ej0*nx].back_vrts[1] = bidx;
					m_vSSMGrid[(ei0-1)+ej0*nx].num_bv++;
				}
			}
			if(ei0 != nx){
				m_vSSMGrid[ei0+ej0*nx].edge_vrts[3] = vidx;
				m_vSSMGrid[ei0+ej0*nx].num_ev++;
				if(bidx != -1){
					m_vSSMGrid[ei0+ej0*nx].back_vrts[3] = bidx;
					m_vSSMGrid[ei0+ej0*nx].num_bv++;
				}
			}
		}
	}

	return ev_num;
}

/*!
 * �O�p�`���b�V������
 * @param[in] nx,ny ���b�V�������p�O���b�h�𑜓x
 * @param[in] dgrid ���b�V���O���b�h
 * @param[out] polys �O�p�`���b�V��
 * @param[in] vstart ���_�C���f�b�N�X�̎n�_
 * @return �������ꂽ�O�p�`���b�V����
 */
int rxSSMeshCPU::CalMesh(int nx, int ny, vector<rxSSGrid> &grid, vector< vector<int> > &polys, int vstart)
{
	int v[14];
	// �m�[�h���_
	// 3 - 2
	// |   |
	// 0 - 1
		 
	// �G�b�W���_
	// - 6 -
	// 7   5
	// - 4 -
		 
	// �G�b�W���_(back vertex)
	//  - 10 -
	// 11     9
	//  -  8 -
	// back-2 vertex : 12,13

	int num_mesh = 0;
	for(int j = 0; j < ny; ++j){
		for(int i = 0; i < nx; ++i){
			rxSSGrid *g = &grid[i+j*nx];

			int table_index = 0;
			for(int k = 0; k < 4; ++k){
				v[k]   = g->node_vrts[k];
				v[k+4] = g->edge_vrts[k];
				v[k+8] = g->back_vrts[k];

				table_index |= ((v[k] != -1) ? BITS[k] : 0);		// �m�[�h���_����4�r�b�g
				table_index |= ((v[k+4] != -1) ? BITS[k]*16 : 0);	// �G�b�W���_���4�r�b�g
			}
			v[12] = -1;
			v[13] = -1;

			int rotation = 0;

			g->table_index0 = table_index;	// �f�o�b�O�p

			int update_idx = g->num_ev*5+g->num_nv;
			if(m_FuncTableIndex[update_idx] != 0){
				(this->*m_FuncTableIndex[update_idx])(table_index, rotation, v, g);
			}

			g->table_index1 = table_index;	// �f�o�b�O�p

			if(g_MeshTable[table_index][0] > 0){	// �O���b�h���̃��b�V������0���傫�������烁�b�V������
				vector<int> tri;
				tri.resize(3);

				int m = g_MeshTable[table_index][0];	// �O���b�h���̃��b�V����
				for(int k = 0; k < m; ++k){
					for(int l = 0; l < 3; ++l){
						int idx = g_VrtRotTable[rotation][g_MeshTable[table_index][k*3+l+1]];

						if(v[idx] == -1){
							v[idx] = 0;
							cout << "mesh error : " << i << ", " << j << endl;
						}

						tri[l] = v[idx]+vstart;
					}
					polys.push_back(tri);

					m_vMeshGrid.push_back(i+j*nx);	// �f�o�b�O�p
					g->mesh[k] = num_mesh;				// �f�o�b�O�p

					num_mesh++;
				}

				// �f�o�b�O�p
				g->mesh_num = m;
			}
		}
	}

	return num_mesh;
}



/*!
 * ���_�@���v�Z
 * @param[in] vrts ���_���W
 * @param[in] nvrts ���_��
 * @param[in] tris �O�p�`�|���S���􉽏��
 * @param[in] ntris �O�p�`�|���S����
 * @param[out] nrms �@��
 */
void rxSSMeshCPU::CalVertexNormals(const vector<Vec3> &vrts, unsigned int nvrts, const vector< vector<int> > &tris, unsigned int ntris, vector<Vec3> &nrms)
{
	unsigned int nnrms = nvrts;
	nrms.resize(nnrms);
	
	// �@���z��̏�����
	for(unsigned int i = 0; i < nnrms; i++){
		nrms[i][0] = 0;
		nrms[i][1] = 0;
		nrms[i][2] = 0;
	}

	// ���_�@���̌v�Z
	for(unsigned int i = 0; i < ntris; i++){
		Vec3 vec1, vec2, normal;
		int id0, id1, id2;
		id0 = tris[i][0];
		id1 = tris[i][1];
		id2 = tris[i][2];

		vec1 = vrts[id1]-vrts[id0];
		vec2 = vrts[id2]-vrts[id0];
		normal = cross(vec1, vec2);

		nrms[id0] += normal;
		nrms[id1] += normal;
		nrms[id2] += normal;
	}

	// �@�����K��
	for(unsigned int i = 0; i < nnrms; i++){
		normalize(nrms[i]);
	}
}



/*!
 * �񕪒T��(2�_�� in 2D)
 * @param[in] v1,v2 �T���n�_�C�I�_ (v1�̕����f�v�X�l��������)
 * @param[out] vr ���_�ʒu
 * @param[in] zmax �f�v�X�l��
 * @return v1����̋���
 */
double rxSSMeshCPU::binarySearchDepth(Vec3 v1, Vec3 v2, Vec3 &vr, double zmax)
{
	//vr = 0.5*(v1+v2);

	const int JMAX = 40;
	double dx, f, fmid, xmid, rtb;
	double x1 = 0;
	double x2 = norm(v1-v2);
	Vec3 dv = Unit(v2-v1);

	f = getDepthValue(v1);
	fmid = getDepthValue(v2);
	if(fabs(f-fmid) < zmax){
		vr = v1+2*RX_FEQ_EPS*dv;
		vr[2] = f+2*RX_FEQ_EPS;
		return 0.0;
	}


	rtb = f < fmid ? (dx = x2-x1, x1) : (dx = x1-x2, x2);
	vr = (f < fmid) ? v1 : v2;
	for(int j = 0; j < JMAX; ++j){
		dx *= 0.5;
		xmid = rtb+dx;

		vr = v1+dv*xmid;
		fmid = getDepthValue(vr);
		
		if(fabs(f-fmid) < zmax){
			rtb = xmid;
		}

		if(fabs(dx) < 0.5){
			vr[2] = getDepthValue(vr-dx*dv);
			return xmid;
		}
	}

	return 0.0;
}

/*!
 * �e�[�u���C���f�b�N�X�X�V�Ȃ��p�֐�
 * @param[inout] table_index ���̃e�[�u���C���f�b�N�X(���4�r�b�g���G�b�W���_,����4�r�b�g���m�[�h���_)
 * @param[inout] vrot ���_���[�e�[�V������
 * @param[in] v �O���b�h���̒��_�C���f�b�N�X
 * @param[inout] g ���b�V���𐶐�����O���b�h
 */
void rxSSMeshCPU::updateTableIndexAny(int &table_index, int &vrot, int v[], rxSSGrid *g)
{
}

/*!
 * �p�^�[��0 �������b�V���p�̃e�[�u���C���f�b�N�X�X�V
 * @param[inout] table_index ���̃e�[�u���C���f�b�N�X(���4�r�b�g���G�b�W���_,����4�r�b�g���m�[�h���_)
 * @param[inout] vrot ���_���[�e�[�V������
 * @param[in] v �O���b�h���̒��_�C���f�b�N�X
 * @param[inout] g ���b�V���𐶐�����O���b�h
 */
void rxSSMeshCPU::updateTableIndexE0N4(int &table_index, int &vrot, int v[], rxSSGrid *g)
{
	// �O���b�h�ԍ��ɂ��90�x��]������
	table_index -= ((g->i+g->j) & 0x01) ? 1 : 0;
}

/*!
 * �p�^�[��1,2,4,8 �����֊s(�֊s�̎n�_)�p�̃e�[�u���C���f�b�N�X�X�V
 * @param[inout] table_index ���̃e�[�u���C���f�b�N�X(���4�r�b�g���G�b�W���_,����4�r�b�g���m�[�h���_)
 * @param[inout] vrot ���_���[�e�[�V������
 * @param[in] v �O���b�h���̒��_�C���f�b�N�X
 * @param[inout] g ���b�V���𐶐�����O���b�h
 */
void rxSSMeshCPU::updateTableIndexE1(int &table_index, int &vrot, int v[], rxSSGrid *g)
{
	// �m�[�h���_�� back layer �ɑ����镨��T��
	// ���_������G�b�W�̒[�_�Ńf�v�X�l���傫������ back layer
	int kb = 0;
	for(int k = 0; k < 4; ++k){
		if(v[k+4] != -1){
			int k1 = (k == 3 ? 0 : k+1);
			kb = ((g->node_depth[k] > g->node_depth[k1]) ? BITS[k] : BITS[k1]);
			break;
		}
	}

	kb = (~kb & 0x0F);
	

	// �o�b�N���C���[�ɂ���m�[�h���_�r�b�g��0�ɂ���
	table_index = (table_index & 0xF0)+kb;
}

/*!
 * �p�^�[��3,5,6,9,10,12 �����֊s�p�̃e�[�u���C���f�b�N�X�X�V
 * @param[inout] table_index ���̃e�[�u���C���f�b�N�X(���4�r�b�g���G�b�W���_,����4�r�b�g���m�[�h���_)
 * @param[inout] vrot ���_���[�e�[�V������
 * @param[in] v �O���b�h���̒��_�C���f�b�N�X
 * @param[inout] g ���b�V���𐶐�����O���b�h
 */
void rxSSMeshCPU::updateTableIndexE2N4(int &table_index, int &vrot, int v[], rxSSGrid *g)
{
	// back vertex�̃m�[�h���_�r�b�g��0�Ƃ����r�b�g��𐶐�
	int btable = 0;
	int k0 = 0, k1;
	for(int k = 0; k <= 4; ++k){
		k1 = (k0 == 3 ? 0 : k0+1);
		if(v[k0+4] != -1){	// �m�[�h�ԂɃG�b�W���_����
			btable |= (g->node_depth[k0] > g->node_depth[k1]) ? BITS[k0] : BITS[k1];
		}
		else{
			btable |= ((btable & BITS[k0]) ? BITS[k1] : 0);
		}
		k0++;
		if(k0 == 4) k0 = 0;
	}

	// �O���֊s�Ƌ�ʂ��邽�߂ɁC�����֊s�̏ꍇ�r�b�g���+2����
	btable = (btable+2 & 0x0F);

	table_index = (table_index & 0xF0)+btable;
}

/*!
 * 7,11,13,14 �O��/�����֊s���ݗp�̃e�[�u���C���f�b�N�X�X�V
 * @param[inout] table_index ���̃e�[�u���C���f�b�N�X(���4�r�b�g���G�b�W���_,����4�r�b�g���m�[�h���_)
 * @param[inout] vrot ���_���[�e�[�V������
 * @param[in] v �O���b�h���̒��_�C���f�b�N�X
 * @param[inout] g ���b�V���𐶐�����O���b�h
 */
void rxSSMeshCPU::updateTableIndexE3N23(int &table_index, int &vrot, int v[], rxSSGrid *g)
{
	int btable = 0;	// �m�[�h���_��2,3�̏ꍇ�C1�Ԗڂ̃r�b�g��0(0xxx)
	int pattern = (table_index >> 4);
	int node = g_EdgeTable[pattern][0]-4;
	int R[4];	// ���_���Ȃ��G�b�W���甽���v���Ƀm�[�h���_�ԍ�����ׂ����X�g
	for(int k = 0; k < 4; ++k){
		R[k] = node++;
		if(node == 4) node = 0;
	}

	int ntable = (table_index & 0x0F);	// �m�[�h���_�r�b�g��

	// 2�Ԗڂ̃r�b�g
	btable |= (((ntable >> R[1]) & 1) ? 4 : 0);

	// 3�Ԗڂ̃r�b�g
	btable |= (((ntable >> R[2]) & 1) ? 2 : 0);

	// R�ɏ]��btable����ёւ�
	int btable0 = 0;
	btable0 |= (((ntable >> R[0]) & 1) ? 8 : 0);
	btable0 |= (((ntable >> R[1]) & 1) ? 4 : 0);
	btable0 |= (((ntable >> R[2]) & 1) ? 2 : 0);
	btable0 |= (((ntable >> R[3]) & 1) ? 1 : 0);

	// ���ёւ���btable�̉���2�r�b�g-1
	int n0 = (btable0 & 3)-1;
	int n1 = n0+1;

	if(n0 != 1){
		// �ǉ���back vertex��g_EdgeTable[pattern][1]�̈ʒu�ɒǉ�
		int add_edge = g_EdgeTable[pattern][1]-4;	// �ǉ�����G�b�W�̈ʒu
		Vec3 add_vrt = m_vSSVertex[g->edge_vrts[add_edge]].pos;	// �ǉ����_

		// �ǉ����_�f�v�X�l
		int ref_edge = g_EdgeTable[pattern][0]-4;	// �ǉ����_�̃f�v�X�l���Q�Ƃ���G�b�W���_
		int ref_node = (btable & 4) ? R[1] : R[2];
		if(fabs(m_vSSVertex[g->edge_vrts[ref_edge]].pos[2]-g->node_depth[ref_node]) > m_fSSZmax){
			// edge vertex���g�p
			add_vrt[2] = m_vSSVertex[g->edge_vrts[ref_edge]].pos[2];	// �ǉ����_�̃f�v�X�l
		}
		else{
			// back vertex���g�p
			if(g->back_vrts[ref_edge] == -1){
				ref_edge = g_EdgeTable[pattern][2]-4;
			}
			add_vrt[2] = m_vSSVertex[g->back_vrts[ref_edge]].pos[2];	// �ǉ����_�̃f�v�X�l
		}


		// ���_�����X�g�ɒǉ�
		m_vSSVertex.push_back(rxSSVertex(add_vrt, 1));
		m_iNumEdgeVrts++;
		int vidx = (int)m_vSSVertex.size()-1;

		// �O���b�h�̒��_���X�g�ɒǉ�
		if(m_vSSVertex[vidx].pos[2] < m_vSSVertex[v[add_edge+4]].pos[2]){
			v[add_edge+8] = v[add_edge+4];
			g->back_vrts[add_edge] = v[add_edge+4];
			g->num_bv++;

			v[add_edge+4] = vidx;
			g->edge_vrts[add_edge] = vidx;
		}
		else{
			v[add_edge+8] = vidx;
			g->back_vrts[add_edge] = vidx;
			g->num_bv++;
		}					
			
	}

	// 4�Ԗڂ̃r�b�g : �m�[�h�̃f�v�X�l���r
	btable |= ((g->node_depth[R[n0]] >= g->node_depth[R[n1]]) ? 0 : 1);

	// �e�[�u���C���f�b�N�X�̃G�b�W���_�������X�V
	table_index &= 0xF0;
	table_index |= btable;
}


/*!
 * �p�^�[��7,11,13,14 �����֊s�p�̃e�[�u���C���f�b�N�X�X�V
 *  - ����Ȃ�back vertex���K�v
 * @param[inout] table_index ���̃e�[�u���C���f�b�N�X(���4�r�b�g���G�b�W���_,����4�r�b�g���m�[�h���_)
 * @param[inout] vrot ���_���[�e�[�V������
 * @param[in] v �O���b�h���̒��_�C���f�b�N�X
 * @param[inout] g ���b�V���𐶐�����O���b�h
 */
void rxSSMeshCPU::updateTableIndexE3N4(int &table_index, int &vrot, int v[], rxSSGrid *g)
{
	int btable = 8;	// �m�[�h���_��4�̏ꍇ�C1�Ԗڂ̃r�b�g��1(1xxx)
	int pattern = (table_index >> 4);

	// �m�[�h�̃f�v�X�l���ׂ鏇�Ԃ̌���
	int node = g_EdgeTable[pattern][0]-4;
	int R[4];	// ���_���Ȃ��G�b�W���甽���v���Ƀm�[�h���_�ԍ�����ׂ����X�g
	for(int k = 0; k < 4; ++k){
		R[k] = node++;
		if(node == 4) node = 0;
	}

	// �m�[�h�̃f�v�X�l�̑召�Ńr�b�g���ύX(2-4�Ԗڂ̃r�b�g)
	for(int k = 0; k < 3; ++k){
		// R[k] > R[k+1]�Ȃ�ΑΉ�����r�b�g��1�ɂ���
		if(g->node_depth[R[k]] > g->node_depth[R[k+1]]){
			btable |= BITS[2-k];
		}
	}

	// �e�[�u���C���f�b�N�X�̃G�b�W���_�������X�V
	table_index &= 0xF0;
	table_index |= btable;


	// g_EdgeTable[pattern][1]�̈ʒu�ɒ��_��ǉ�
	int add_edge = g_EdgeTable[pattern][1]-4;	// �ǉ�����G�b�W�̈ʒu
	Vec3 add_vrt = m_vSSVertex[g->edge_vrts[add_edge]].pos;	// �ǉ����_

	// ���b�V��C�̒��_�����ԂŃf�v�X�l�����߂�
	RXREAL ref_depths[4];
	// �^�_
	// 2�[3
	// |�_|
	// 0�[1

	ref_depths[0] = g->node_depth[R[3]];
	ref_depths[1] = g->node_depth[R[0]];

	int e2 = g_EdgeTable[pattern][2]-4;
	int e3 = g_EdgeTable[pattern][0]-4;

	if(fabs(m_vSSVertex[g->edge_vrts[e2]].pos[2]-ref_depths[0]) < m_fSSZmax){
		ref_depths[2] = m_vSSVertex[g->edge_vrts[e2]].pos[2];
	}
	else{
		ref_depths[2] = m_vSSVertex[g->back_vrts[e2]].pos[2];
	}
	if(fabs(m_vSSVertex[g->edge_vrts[e3]].pos[2]-ref_depths[1]) < m_fSSZmax){
		ref_depths[3] = m_vSSVertex[g->edge_vrts[e3]].pos[2];
	}
	else{
		ref_depths[3] = m_vSSVertex[g->back_vrts[e3]].pos[2];
	}

	// �ǉ����_�̃f�v�X�l
	add_vrt[2] = 0.5*(ref_depths[2]+ref_depths[3]);

	// ���_�����X�g�ɒǉ�
	m_vSSVertex.push_back(rxSSVertex(add_vrt, 1));
	m_iNumEdgeVrts++;
	int vidx = (int)m_vSSVertex.size()-1;

	// �O���b�h�̒��_���X�g�ɒǉ�
	//v[12] = vidx;
	//g->back_vrts[4] = vidx;

	// �O���b�h�̒��_���X�g�ɒǉ�
	if(m_vSSVertex[vidx].pos[2] < m_vSSVertex[v[add_edge+8]].pos[2]){
		if(m_vSSVertex[vidx].pos[2] < m_vSSVertex[v[add_edge+4]].pos[2]){
			// front vertex�Ƃ��đ}��
			v[12] = v[add_edge+8];
			g->back_vrts[4] = v[add_edge+8];

			v[add_edge+8] = v[add_edge+4];
			g->back_vrts[add_edge] = v[add_edge+4];

			v[add_edge+4] = vidx;
			g->edge_vrts[add_edge] = vidx;
		}
		else{
			// back vertex�Ƃ��đ}��
			v[12] = v[add_edge+8];
			g->back_vrts[4] = v[add_edge+8];

			v[add_edge+8] = vidx;
			g->back_vrts[add_edge] = vidx;
		}
	}
	else{
		// back-2 vertex�Ƃ��đ}��
		v[12] = vidx;
		g->back_vrts[4] = vidx;
	}					
}

/*!
 * �p�^�[��15 �O���֊s�p�̃e�[�u���C���f�b�N�X�X�V
 * @param[inout] table_index ���̃e�[�u���C���f�b�N�X(���4�r�b�g���G�b�W���_,����4�r�b�g���m�[�h���_)
 * @param[inout] vrot ���_���[�e�[�V������
 * @param[in] v �O���b�h���̒��_�C���f�b�N�X
 * @param[inout] g ���b�V���𐶐�����O���b�h
 */
void rxSSMeshCPU::updateTableIndexE4N2(int &table_index, int &vrot, int v[], rxSSGrid *g)
{
	int ntable = (table_index & 0x0F);	// �m�[�h���_�r�b�g��
	ntable = (ntable == 5 ? 0 : 15);

	// �e�[�u���C���f�b�N�X�̃m�[�h���_�������X�V
	table_index &= 0xF0;
	table_index |= ntable;
}

/*!
 * �p�^�[��15 �O��/�����֊s���ݗp�̃e�[�u���C���f�b�N�X�X�V
 * @param[inout] table_index ���̃e�[�u���C���f�b�N�X(���4�r�b�g���G�b�W���_,����4�r�b�g���m�[�h���_)
 * @param[inout] vrot ���_���[�e�[�V������
 * @param[in] v �O���b�h���̒��_�C���f�b�N�X
 * @param[inout] g ���b�V���𐶐�����O���b�h
 */
void rxSSMeshCPU::updateTableIndexE4N3(int &table_index, int &vrot, int v[], rxSSGrid *g)
{
	int btable = 0;	// �m�[�h���_��3�̏ꍇ�C1�Ԗڂ̃r�b�g��0(0xxx)
	int pattern = (table_index >> 4);
	int ntable = (table_index & 0x0F);	// �m�[�h���_�r�b�g��

	// ���_���Ȃ��m�[�h
	int zero_node = log((double)(~ntable & 0x0F))/log(2.0);

	// �m�[�h�̃f�v�X�l�̑召�Ńr�b�g���ύX(2-4�Ԗڂ̃r�b�g)
	for(int k = 0; k < 3; ++k){
		int k0 = g_NodeTable[zero_node][k];
		int k1 = g_NodeTable[zero_node][k+1];

		// k0 > k1�Ȃ�ΑΉ�����r�b�g��1�ɂ���
		if(g->node_depth[k0] > g->node_depth[k1]){
			btable |= BITS[2-k];
		}
	}

	// ���_���[�e�[�V����
	vrot = zero_node;

	// �e�[�u���C���f�b�N�X�̃m�[�h���_�������X�V
	table_index &= 0xF0;
	table_index |= btable;
}

/*!
 * �p�^�[��15 �����֊s�p�̃e�[�u���C���f�b�N�X�X�V
 *  - �ǉ���back vertex��2�K�v
 * @param[inout] table_index ���̃e�[�u���C���f�b�N�X(���4�r�b�g���G�b�W���_,����4�r�b�g���m�[�h���_)
 * @param[inout] vrot ���_���[�e�[�V������
 * @param[in] v �O���b�h���̒��_�C���f�b�N�X
 * @param[inout] g ���b�V���𐶐�����O���b�h
 */
void rxSSMeshCPU::updateTableIndexE4N4(int &table_index, int &vrot, int v[], rxSSGrid *g)
{
	int btable = 8;	// �m�[�h���_��4�̏ꍇ�C1�Ԗڂ̃r�b�g��1(1xxx)
	int pattern = (table_index >> 4);
	int ntable = (table_index & 0x0F);	// �m�[�h���_�r�b�g��

	// �f�v�X�l�������Ƃ��傫���m�[�h
	int zero_node = 0;
	double max_depth = 0.0;
	for(int k = 1; k < 4; ++k){
		if(g->node_depth[k] > max_depth){
			max_depth = g->node_depth[k];
			zero_node = k;
		}
	}

	// �m�[�h�̃f�v�X�l�̑召�Ńr�b�g���ύX(2-4�Ԗڂ̃r�b�g)
	for(int k = 0; k < 3; ++k){
		int k0 = g_NodeTable[zero_node][k];
		int k1 = g_NodeTable[zero_node][k+1];

		// k0 > k1�Ȃ�ΑΉ�����r�b�g��1�ɂ���
		if(g->node_depth[k0] > g->node_depth[k1]){
			btable |= BITS[2-k];
		}
	}

	// ���_���[�e�[�V����
	vrot = zero_node;

	// �e�[�u���C���f�b�N�X�̃G�b�W���_�������X�V
	table_index &= 0xF0;
	table_index |= btable;

	// 
	// g_NodeTable[zero_node][5,6]�̈ʒu�ɒ��_(back-2 vertex)��ǉ�
	int add_edge1 = g_NodeTable[zero_node][5]-4;	// �ǉ�����G�b�W�̈ʒu
	int add_edge2 = g_NodeTable[zero_node][6]-4;	// �ǉ�����G�b�W�̈ʒu
	Vec3 add_vrt1 = m_vSSVertex[g->edge_vrts[add_edge1]].pos;	// �ǉ����_
	Vec3 add_vrt2 = m_vSSVertex[g->edge_vrts[add_edge2]].pos;	// �ǉ����_

	// g_NodeTable[zero_node][4,7]�̈ʒu��back vertex�̃f�v�X��ݒ�
	int ref_edge1 = g_NodeTable[zero_node][4]-4;	// �Q�Ƃ���G�b�W�̈ʒu
	int ref_edge2 = g_NodeTable[zero_node][7]-4;	// �Q�Ƃ���G�b�W�̈ʒu
	add_vrt1[2] = m_vSSVertex[g->back_vrts[ref_edge1]].pos[2];
	add_vrt2[2] = m_vSSVertex[g->back_vrts[ref_edge2]].pos[2];

	// ���_�����X�g�ɒǉ�
	m_vSSVertex.push_back(rxSSVertex(add_vrt1, 1));
	m_vSSVertex.push_back(rxSSVertex(add_vrt2, 1));
	m_iNumEdgeVrts += 2;
	int vidx1 = (int)m_vSSVertex.size()-2;
	int vidx2 = (int)m_vSSVertex.size()-1;

	// �O���b�h�̒��_���X�g�ɒǉ�
	v[12] = vidx1;
	v[13] = vidx2;
	g->back_vrts[4] = vidx1;
	g->back_vrts[5] = vidx2;
}


/*!
 * �f�v�X�}�b�v�̎擾
 * @return �f�v�X�}�b�v
 */
RXREAL* rxSSMeshCPU::GetDepthMap(void)
{ 
	return &m_vSSDMap[0];
}


/*!
 * ���b�V�������p�O���b�h�̎擾
 * @param[in] idx �O���b�h�C���f�b�N�X
 * @return ���b�V�������p�O���b�h(rxSSGrid)
 */
rxSSGrid rxSSMeshCPU::GetMeshGrid(int idx)
{
	return m_vSSMGrid[idx];
}

/*!
 * �X�N���[����Ԃł̒��_�����擾
 * @return ���_���(rxSSVertex)
 */
rxSSVertex* rxSSMeshCPU::GetSSVertex(void)
{
	return &m_vSSVertex[0];
}


/*!
 * �֊s�G�b�W��OpenGL�`��
 */
void rxSSMeshCPU::DrawSSEdge(void)
{
	glBegin(GL_LINES);
	for(int i = 0; i < (int)m_vSSEdge.size(); ++i){
		rxSSEdge &e = m_vSSEdge[i];
		Vec3 x0 = e.x0;
		Vec3 x1 = e.x1;

		if(e.silhouette){
			glColor3d(1.0, 0.0, 0.0);
		}
		else{
			glColor3d(0.0, 0.0, 1.0);
		}

		glVertex3d(x0[0], x0[1], 0.01);
		glVertex3d(x1[0], x1[1], 0.01);
	}
	glEnd();
}

/*!
 * �֊s�G�b�W��OpenGL�`��
 */
void rxSSMeshCPU::DrawSilhouetteEdge(void)
{
	glBegin(GL_LINES);
	for(int i = 0; i < (int)m_vSSEdge.size(); ++i){
		rxSSEdge &e = m_vSSEdge[i];
		if(e.silhouette){
			Vec3 x0 = e.x0;
			Vec3 x1 = e.x1;
			glVertex3d(x0[0], x0[1], 0.01);
			glVertex3d(x1[0], x1[1], 0.01);
		}
	}
	glEnd();
}

/*!
 * Screen Space�ł̒��_��OpenGL�`��
 * @param[in] node_color �m�[�h���_�̕`��F
 * @param[in] edge_color �G�b�W���_�̕`��F
 */
void rxSSMeshCPU::DrawSSVertex(Vec3 node_color, Vec3 edge_color)
{
	glBegin(GL_POINTS);
	// �m�[�h���_
	glColor3dv(node_color);
	for(int i = 0; i < m_iNumNodeVrts; ++i){
		Vec3 x = m_vSSVertex[i].pos;
		glVertex3d(x[0], x[1], 0.02);
	}
	// �֊s�G�b�W���_
	glColor3dv(edge_color);
	for(int i = m_iNumNodeVrts; i < (int)m_vSSVertex.size(); ++i){
		Vec3 x = m_vSSVertex[i].pos;
		glVertex3d(x[0], x[1], 0.025);
	}
	glEnd();
}

/*!
 * �O���b�h�����_��OpenGL�`��
 * @param[in] grid �O���b�h�C���f�b�N�X
 * @param[in] colors ���_�̐F
 */
void rxSSMeshCPU::DrawMeshGrid(int grid, const Vec3 colors[])
{
	// �O���b�h���̒��_
	rxSSGrid g = GetMeshGrid(grid);
	for(int i = 0; i < 4; ++i){
		glPointSize(5.0);
		// �m�[�h���_
		if(g.node_vrts[i] != -1){
			Vec3 x = m_vSSVertex[g.node_vrts[i]].pos;
			glColor3dv(colors[0]);
			glBegin(GL_POINTS);
			glVertex3d(x[0], x[1], 0.03);
			glEnd();
		}
		// �G�b�W���_
		if(g.edge_vrts[i] != -1){
			Vec3 x = m_vSSVertex[g.edge_vrts[i]].pos;
			glColor3dv(colors[1]);
			glBegin(GL_POINTS);
			glVertex3d(x[0], x[1], 0.04);
			glEnd();
		}
		glPointSize(8.0);
		// �G�b�W���_(back vertex)
		if(g.back_vrts[i] != -1){
			Vec3 x = m_vSSVertex[g.back_vrts[i]].pos;
			glColor3dv(colors[2]);
			glBegin(GL_POINTS);
			glVertex3d(x[0], x[1], 0.039);
			glEnd();
		}
	}
	// �G�b�W���_(back-2 vertex)
	glPointSize(12.0);
	glColor3dv(colors[3]);
	glBegin(GL_POINTS);
	for(int i = 0; i < 2; ++i){
		if(g.back_vrts[i+4] != -1){
			Vec3 x = m_vSSVertex[g.back_vrts[i+4]].pos;
			glVertex3d(x[0], x[1], 0.038);
		}
	}
	glEnd();
	glPointSize(1.0);
}



/*!
 * ���x��̕`��
 * @param[in] minpos[2],maxpos[2] �`��̈�͈̔�
 */
void rxSSMeshCPU::DrawField(double minpos[2], double maxpos[2])
{
	double x, y, h, d00, d01, d10, d11;
	double lx = maxpos[0]-minpos[0];
	double ly = maxpos[1]-minpos[1];

	h = (lx < ly) ? lx/m_iNgx : ly/m_iNgy;

	glBegin(GL_QUADS);

	double a = 1.0, b = 0.0;
	for(int i = 0 ; i < m_iNgx; ++i){
		x = (i)*h;
		for(int j = 0; j < m_iNgy; ++j){
			y = (j)*h;

			d00 = RX_DEPTH2COLORf(m_vSSDMap[i+j*(m_iNgx+1)])*a+b;
			d01 = RX_DEPTH2COLORf(m_vSSDMap[i+(j+1)*(m_iNgx+1)])*a+b;
			d10 = RX_DEPTH2COLORf(m_vSSDMap[(i+1)+j*(m_iNgx+1)])*a+b;
			d11 = RX_DEPTH2COLORf(m_vSSDMap[(i+1)+(j+1)*(m_iNgx+1)])*a+b;

			glColor3d(d00, d00, d00); glVertex3d(x, y, 0.0);
			glColor3d(d10, d10, d10); glVertex3d(x+h, y, 0.0);
			glColor3d(d11, d11, d11); glVertex3d(x+h, y+h, 0.0);
			glColor3d(d01, d01, d01); glVertex3d(x, y+h, 0.0);
		}
	}

	glEnd();
}


/*!
 * �O���b�h�Ɋւ�����̏o��
 * @param[in] grid �O���b�h�C���f�b�N�X(i+j*nx)
 */
void rxSSMeshCPU::OutputGridInfo(int grid)
{
	if(grid < 0){
		return;
	}

	rxSSGrid g = m_vSSMGrid[grid];
	cout << "  node vertex (" << g.num_nv << ") : ";
	for(int i = 0; i < 4; ++i) cout << g.node_vrts[i] << ((i == 3) ? "\n" : ", ");

	cout << "  edge vertex (" << g.num_ev << ") : ";
	for(int i = 0; i < 4; ++i) cout << g.edge_vrts[i] << ((i == 3) ? "\n" : ", ");

	cout << "  back vertex (" << g.num_bv << ") : ";
	for(int i = 0; i < 4; ++i) cout << g.back_vrts[i] << ((i == 3) ? "\n" : ", ");
	
	int b2num = 0;
	if(g.back_vrts[4] != -1) b2num++;
	if(g.back_vrts[5] != -1) b2num++;
	cout << "  back-2 vertex (" << b2num << ") : ";
	cout << g.back_vrts[4] << ", " << g.back_vrts[5] << endl;

	cout << "  mesh pattern : " << (g.table_index1 >> 4) << "(" << GetBitArray((g.table_index1 >> 4), 4) << ")" << endl;
	cout << "  node pattern : " << GetBitArray((g.table_index0 & 0x0F), 4);
	cout << " -> " << GetBitArray((g.table_index1 & 0x0F), 4);
	cout << endl;
}

/*!
 * �O���b�h�Ɋւ�����̏o��(���_���܂�)
 * @param[in] grid �O���b�h�C���f�b�N�X(i+j*nx)
 */
void rxSSMeshCPU::OutputGridVertexInfo(int grid)
{
	if(grid < 0){
		cout << "no grid is selected." << endl;
		return;
	}

	int gx = grid%m_iNgx;
	int gy = grid/m_iNgx;
	cout << "grid : " << gx << ", " << gy << endl;

	rxSSGrid g = m_vSSMGrid[grid];
	cout << "  node vertex (" << g.num_nv << ") : ";
	for(int i = 0; i < 4; ++i){
		cout << ((i == 0) ? "" : "                    ");
		cout << g.node_vrts[i] << " - ";

		if(g.node_vrts[i] != -1){
			Vec3 x0 = m_vSSVertex[g.node_vrts[i]].pos;
			cout << Vec3(x0[0], x0[1], x0[2]);
		}
		else{
			cout << "no vertex";
		}

		cout << " - d = " << g.node_depth[i];

		cout << endl;
	}

	cout << "  edge vertex (" << g.num_ev << ") : ";
	for(int i = 0; i < 4; ++i){
		cout << ((i == 0) ? "" : "                    ");
		cout << g.edge_vrts[i] << " - ";

		if(g.edge_vrts[i] != -1){
			Vec3 x0 = m_vSSVertex[g.edge_vrts[i]].pos;
			cout << Vec3(x0[0], x0[1], x0[2]);
		}
		else{
			cout << "no vertex";
		}

		cout << endl;
	}

	cout << "  back vertex (" << g.num_bv << ") : ";
	for(int i = 0; i < 4; ++i){
		cout << ((i == 0) ? "" : "                    ");
		cout << g.back_vrts[i] << " - ";

		if(g.back_vrts[i] != -1){
			Vec3 x0 = m_vSSVertex[g.back_vrts[i]].pos;
			cout << Vec3(x0[0], x0[1], x0[2]);
		}
		else{
			cout << "no vertex";
		}

		cout << endl;
	}

	int b2num = 0;
	if(g.back_vrts[4] != -1) b2num++;
	if(g.back_vrts[5] != -1) b2num++;
	cout << "  back-2 vertex (" << b2num << ") : ";
	for(int i = 0; i < 2; ++i){
		cout << ((i == 0) ? "" : "                    ");
		cout << g.back_vrts[i+4] << " - ";

		if(g.back_vrts[i+4] != -1){
			Vec3 x0 = m_vSSVertex[g.back_vrts[i+4]].pos;
			cout << Vec3(x0[0], x0[1], x0[2]);
		}
		else{
			cout << "no vertex";
		}

		cout << endl;
	}

	cout << endl;

	cout << "  mesh pattern : " << (g.table_index1 >> 4) << "(" << GetBitArray((g.table_index1 >> 4), 4) << ")" << endl;
	cout << "  node pattern : " << GetBitArray((g.table_index0), 8);
	cout << " -> " << GetBitArray((g.table_index1), 8);
	cout << endl;

	cout << "  back2 : " << GetBitArray(g.back2, 8) << endl;
	cout << "  v : ";
	for(int i = 0; i < 14; ++i) cout << g.v[i] << ((i == 13) ? "\n" : ", ");

	cout << "  the number of mesh = " << g.mesh_num << "  (";
	for(int i = 0; i < g.mesh_num; ++i){
		cout << g.mesh[i] << ((i == g.mesh_num-1) ? ")\n" : ", ");
	}

	cout << endl;


}
