/*!
  @file rx_ssm_cpu.cpp
	
  @brief Screen Space Mesh作成
		- M. Muller et al., Screen space meshes, SCA2007, 2007. 

  @author Makoto Fujisawa
  @date 2011-05
*/
// FILE --rx_ssm_cpu.cpp--


//-----------------------------------------------------------------------------
// インクルードファイル
//-----------------------------------------------------------------------------
#include "rx_ssm.h"


//-----------------------------------------------------------------------------
// 定数・変数
//-----------------------------------------------------------------------------
const int BITS[4] = { 1, 2, 4, 8 };


//-----------------------------------------------------------------------------
// rxSSMeshCPUクラスの実装
//-----------------------------------------------------------------------------
/*!
 * コンストラクタ
 */
rxSSMeshCPU::rxSSMeshCPU(double zmax, double h, double r, int n_filter, int n_iter)
	: rxSSMesh(zmax, h, r, n_filter, n_iter)
{
	// MARK:コンストラクタ
	m_vFilter = CalBinomials(21);

	// メッシュ構成テーブル用インデックス更新関数ポインタ ev*5+nv
	// エッジ頂点数0
	m_FuncTableIndex[0]  = 0;
	m_FuncTableIndex[1]  = 0;
	m_FuncTableIndex[2]  = 0;
	m_FuncTableIndex[3]  = 0;
	m_FuncTableIndex[4]  = &rxSSMeshCPU::updateTableIndexE0N4;
	// エッジ頂点数1
	m_FuncTableIndex[5]  = &rxSSMeshCPU::updateTableIndexE1;
	m_FuncTableIndex[6]  = &rxSSMeshCPU::updateTableIndexE1;
	m_FuncTableIndex[7]  = &rxSSMeshCPU::updateTableIndexE1;
	m_FuncTableIndex[8]  = &rxSSMeshCPU::updateTableIndexE1;
	m_FuncTableIndex[9]  = &rxSSMeshCPU::updateTableIndexE1;
	// エッジ頂点数2
	m_FuncTableIndex[10] = 0;
	m_FuncTableIndex[11] = 0;
	m_FuncTableIndex[12] = 0;
	m_FuncTableIndex[13] = 0;
	m_FuncTableIndex[14] = &rxSSMeshCPU::updateTableIndexE2N4;
	// エッジ頂点数3
	m_FuncTableIndex[15] = 0;
	m_FuncTableIndex[16] = 0;
	m_FuncTableIndex[17] = &rxSSMeshCPU::updateTableIndexE3N23;
	m_FuncTableIndex[18] = &rxSSMeshCPU::updateTableIndexE3N23;
	m_FuncTableIndex[19] = &rxSSMeshCPU::updateTableIndexE3N4;
	// エッジ頂点数4
	m_FuncTableIndex[20] = 0;
	m_FuncTableIndex[21] = 0;
	m_FuncTableIndex[22] = &rxSSMeshCPU::updateTableIndexE4N2;
	m_FuncTableIndex[23] = &rxSSMeshCPU::updateTableIndexE4N3;
	m_FuncTableIndex[24] = &rxSSMeshCPU::updateTableIndexE4N4;
}


/*!
 * デストラクタ
 */
rxSSMeshCPU::~rxSSMeshCPU()
{
}

/*!
 * マップやグリッド配列のサイズを変更
 * @param[in] W,H 画面解像度
 * @param[in] spacing デプスマップのサンプリング間隔
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
	// CPUメモリ確保
	//
	// デプスマップ
	if((int)m_vSSDMap.size() != (m_iNgx+1)*(m_iNgy+1)){
		m_vSSDMap.clear();
		m_vSSDMap.resize((m_iNgx+1)*(m_iNgy+1));
	}

	// メッシュ生成用グリッド
	if((int)m_vSSMGrid.size() != m_iNgx*m_iNgy){
		m_vSSMGrid.clear();
		m_vSSMGrid.resize(m_iNgx*m_iNgy);
	}
}

/*!
 * スクリーンスペースメッシュ生成(法線計算含む)
 *  - M. Muller et al., Screen space meshes, SCA2007, 2007. 
 * @param[in] proj OpenGL透視投影行列
 * @param[in] modelview OpenGLモデルビュー行列
 * @param[in] W,H 画面解像度
 * @param[in] prts パーティクル座標
 * @param[in] pnum パーティクル数
 * @param[out] vrts メッシュ頂点列
 * @param[out] nrms 頂点法線列
 * @param[out] polys メッシュ列
 * @param[in] filtering フィルタリングフラグ(デプス値0x01, 輪郭0x02)
 * @param[in] debug_output 結果の画面出力の有無
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
 * スクリーンスペースメッシュ生成
 *  - M. Muller et al., Screen space meshes, SCA2007, 2007. 
 * @param[in] proj OpenGL透視投影行列
 * @param[in] modelview OpenGLモデルビュー行列
 * @param[in] W,H 画面解像度
 * @param[in] prts パーティクル座標
 * @param[in] pnum パーティクル数
 * @param[out] vrts メッシュ頂点列
 * @param[out] polys メッシュ列
 * @param[in] filtering フィルタリングフラグ(デプス値0x01, 輪郭0x02)
 * @param[in] debug_output 結果の画面出力の有無
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
	// コンテナ初期化
	//

	// Screen Spaceでのパーティクル
	if((int)m_vSSPrts.size() != pnum){
		m_vSSPrts.clear();
		m_vSSPrts.resize(pnum);
	}

	// 輪郭エッジ
	m_vSSEdge.clear();

	// Screen Spaceでのメッシュ頂点
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

	// メッシュ初期化
	vrts.clear();
	polys.clear();
	m_vMeshGrid.clear();


	
	// 透視投影変換行列
	const rxMatrix4 P = GetMatrixGL(proj);

	// モデルビュー変換行列
	const rxMatrix4 MV = GetMatrixGL(modelview);


	//
	// パーティクル座標と半径を透視投影変換してデプスマップを生成
	//
	// デプスマップ初期化	
	for(int i = 0; i < (m_iNgx+1)*(m_iNgy+1); ++i){
		m_vSSDMap[i] = RX_FEQ_INF;
	}

	// デプスマップ作成
	CalDepthMap(P, MV, W, H, prts, pnum);



	// 
	// デプスマップの平滑化
	//
	if(filtering & 0x01){
		ApplyDepthFilter(m_vSSDMap, m_iNgx+1, m_iNgy+1, m_iNfilter);
	}
	



	//
	// ノード頂点生成
	//
	m_iNumNodeVrts = CalNodeVertex(m_iNgx, m_iNgy, m_fDmx, m_fDmy, m_vSSDMap);
	if(debug_output) cout << "the number of node vertices = " << m_iNumNodeVrts << endl;


	//
	// 輪郭エッジの検出とエッジ頂点生成
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
	
	// 輪郭エッジ上のback vertex
	m_iNumEdgeVrts = CalEdgeVertex(m_iNgx, m_iNgy, m_fDmx, m_fDmy, m_vSSDMap, m_vSSEdge);
	if(debug_output) cout << "the number of edge vertices = " << m_iNumEdgeVrts << endl;



	// 
	// メッシュ生成
	// 
	m_iNumMesh = CalMesh(m_iNgx, m_iNgy, m_vSSMGrid, polys, 0);
	if(debug_output) cout << "the number of mesh = " << m_iNumMesh << endl;



	// 
	// 輪郭の平滑化
	//
	if(filtering & 0x02){
		ApplySilhoutteSmoothing(m_vSSVertex, polys, m_iNiters);
	}



	//
	// 3D空間にメッシュを戻す
	//
	vrts.clear();

	// 透視投影変換逆行列
	rxMatrix4 Q = P.Inverse();

	// モデルビュー変換逆行列
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
 * パーティクル座標と半径を透視投影変換してデプスマップを生成
 * @param[in] P 透視投影行列
 * @param[in] MV モデルビュー行列
 * @param[in] W,H 画面解像度
 * @param[in] prts パーティクル座標
 * @param[in] pnum パーティクル数
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

		// 投影変換
		Vec4 xd = PMV*x;

		// wで割ることで[-1, 1]の正規化座標系に変換
		Vec3 xp;
		xp[0] = W*(0.5+0.5*xd[0]/xd[3]);
		xp[1] = H*(0.5+0.5*xd[1]/xd[3]);
		xp[2] = xd[2];		// zだけは元の値を用いる

		m_vSSPrts[k].xp = xp;

		Vec4 rd = Vec4(m_fPrtRad, m_fPrtRad, m_fPrtRad, 1.0);
		rd = PMV*rd;

		// 正規化座標系での半径値
		Vec3 rp;
		rp[0] = tr[0]/xd[3]/2;
		rp[1] = tr[1]/xd[3]/2;
		rp[2] = tr[2];

		m_vSSPrts[k].rp = rp;

		double rrp = rp[0]*rp[0];

		if(xp[0] < 0 || xp[0] >= W || xp[1] < 0 || xp[1] >= H){
			continue;
		}

		// デプスマップ上でのパーティクルの範囲
		int cen[2];	// パーティクル中心
		cen[0] = int(xp[0]/m_fDmx)+1;
		cen[1] = int(xp[1]/m_fDmy)+1;

		int minp[2], maxp[2];
		minp[0] = cen[0]-(rp[0]/m_fDmx+2);
		minp[1] = cen[1]-(rp[1]/m_fDmy+2);
		maxp[0] = cen[0]+(rp[0]/m_fDmx+2);
		maxp[1] = cen[1]+(rp[1]/m_fDmy+2);

		// 範囲がマップ外にならないようにクランプ
		RX_CLAMP2(minp[0], 0, m_iNgx);
		RX_CLAMP2(minp[1], 0, m_iNgy);
		RX_CLAMP2(maxp[0], 0, m_iNgx);
		RX_CLAMP2(maxp[1], 0, m_iNgy);

		// パーティクルデプス値更新
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
 * デプスマップに平滑化を施す
 * @param[inout] dmap デプスマップ
 * @param[in] nx,ny マップ解像度
 * @param[in] n_filter フィルタリング幅
 */
void rxSSMeshCPU::ApplyDepthFilter(vector<RXREAL> &dmap, int nx, int ny, int n_filter)
{
	int b = 2*n_filter+1;

	// binomial係数の計算
	if((int)m_vFilter.size() != b){
		m_vFilter = CalBinomials(b);
	}

	vector<RXREAL> tmp_map;
	tmp_map = dmap;

	vector<RXREAL> d;
	d.resize(b);
	
	// 水平方向のフィルターを適用
	for(int j = 0; j < ny; ++j){
		for(int i = n_filter; i < nx-n_filter; ++i){
			d[0] = tmp_map[i+j*nx];
			if(d[0] < RX_FEQ_INF-1){	// != ∞のピクセルのみに適用
				// 周囲Nfilterグリッドの値を取得
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

	// 垂直方向のフィルターを適用
	for(int i = 0; i < nx; ++i){
		for(int j = n_filter; j < ny-n_filter; ++j){
			d[0] = tmp_map[i+j*nx];
			if(d[0] < RX_FEQ_INF-1){	// != ∞のピクセルのみに適用
				// 周囲Nfilterグリッドの値を取得
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
 * 輪郭に平滑化を施す
 * @param[inout] ssvrts スクリーン座標系での頂点列
 * @param[in] polys メッシュ(構成する頂点列)
 * @param[in] n_iter フィルタリング反復数
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

					// 隣接頂点
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
 * 輪郭エッジの検出とfront edge vertexの計算
 * @param[in] nx,ny メッシュ生成用グリッド解像度
 * @param[in] dw,dh メッシュ生成用グリッド幅
 * @param[in] dgrid メッシュ生成用デプスマップ
 * @param[in] ssprts Screen Spaceでのパーティクル
 * @param[in] W,H スクリーン解像度
 * @return 検出された輪郭エッジ数
 */
int rxSSMeshCPU::DetectSilhouetteEdgeVertex(int nx, int ny, double dw, double dh, const vector<RXREAL> &dgrid, 
										 const vector<rxSSParticle> &ssprts, int W, int H)
{
	// DetectSilhouetteEdgeVertex
	int en = 0;
	int pnum = (int)ssprts.size();

	m_vSSEdge.resize((nx)*(ny+1)+(nx+1)*(ny));
	int yoffset = (nx)*(ny+1);
	
	// x方向エッジ初期化
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

	// y方向エッジのチェック
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

		// デプスマップ上でのパーティクルの範囲
		int cen[2];	// パーティクル中心
		cen[0] = int(xp[0]/dw)+1;
		cen[1] = int(xp[1]/dh)+1;

		int minp[2], maxp[2];
		minp[0] = cen[0]-(rp[0]/dw+1);
		minp[1] = cen[1]-(rp[1]/dh+1);
		maxp[0] = cen[0]+(rp[0]/dw+1);
		maxp[1] = cen[1]+(rp[1]/dh+1);

		// 範囲がマップ外にならないようにクランプ
		RX_CLAMP2(minp[0], 0, nx-1);
		RX_CLAMP2(minp[1], 0, ny-1);
		RX_CLAMP2(maxp[0], 0, nx-1);
		RX_CLAMP2(maxp[1], 0, ny-1);

		// 範囲内のエッジを調査(x方向)
		for(int j = minp[1]; j <= maxp[1]+1; ++j){
			for(int i = minp[0]; i <= maxp[0]; ++i){
				rxSSEdge &e = m_vSSEdge[i+j*(nx)];
				if(!e.silhouette) continue;
				if(e.depth < xp[2]) continue;

				// 円とエッジの交点
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
		// 範囲内のエッジを調査(x方向)
		for(int i = minp[0]; i <= maxp[0]+1; ++i){
			for(int j = minp[1]; j <= maxp[1]; ++j){
				rxSSEdge &e = m_vSSEdge[yoffset+i+j*(nx+1)];
				if(!e.silhouette) continue;
				if(e.depth < xp[2]) continue;

				// 円とエッジの交点
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
 * 輪郭エッジの検出
 * @param[in] nx,ny メッシュ生成用グリッド解像度
 * @param[in] dw,dh メッシュ生成用グリッド幅
 * @param[in] dgrid メッシュ生成用デプスマップ
 * @return 検出された輪郭エッジ数
 */
int rxSSMeshCPU::DetectSilhouetteEdge(int nx, int ny, double dw, double dh, const vector<RXREAL> &dgrid)
{
	int en = 0;

	// x方向エッジのチェック
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

	// y方向エッジのチェック
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
 * ノード頂点生成
 * @param[in] nx,ny メッシュ生成用グリッド解像度
 * @param[in] dw,dh メッシュ生成用グリッド幅
 * @param[in] dgrid メッシュ生成用デプスマップ
 * @return 生成されたノード頂点数
 */
int rxSSMeshCPU::CalNodeVertex(int nx, int ny, double dw, double dh, const vector<RXREAL> &dgrid)
{
	// グリッドノード上の頂点
	int nv_num = 0;
	for(int j = 0; j <= ny; ++j){
		for(int i = 0; i <= nx; ++i){
			double d = dgrid[i+j*(nx+1)];
			if(d < RX_FEQ_INF-1){
				Vec3 vrt_pos = Vec3(dw*i, dh*j, d);

				m_vSSVertex.push_back(rxSSVertex(vrt_pos, 0));
				nv_num++;
				int vidx = (int)m_vSSVertex.size()-1;

				// 頂点情報を格納
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

			// デプス値の格納
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
 * エッジ頂点生成
 * @param[in] nx,ny メッシュ生成用グリッド解像度
 * @param[in] dw,dh メッシュ生成用グリッド幅
 * @param[in] dgrid メッシュ生成用デプスマップ
 * @param[in] edges 輪郭エッジリスト
 * @return 生成されたエッジ頂点数
 */
int rxSSMeshCPU::CalEdgeVertex(int nx, int ny, double dw, double dh, const vector<RXREAL> &dgrid, const vector<rxSSEdge> &edges)
{
	int ev_num = 0;
	int yoffset = (nx)*(ny+1);
	for(int i = 0; i < (int)edges.size(); ++i){
		rxSSEdge e = edges[i];
		if(!e.silhouette) continue;

		int ei0, ei1, ej0, ej1;
		if(i < yoffset){	// x軸平行エッジ
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


		// back vertexがあるかどうかをチェック
		int bidx = -1;
		Vec3 back_vrt;
		if(e.d0 < RX_FEQ_INFM && e.d1 < RX_FEQ_INFM){	// エッジ端点が両方とも != ∞ ならば back vertex が存在
			RXREAL back_node_depth, nn_back_node_depth, back_depth;
			int back_node, nn_back_node;	// back layerに属するエッジ端点とその隣接ノード
			double l = 1.0;	// エッジの長さ

			// デプス値が大きい方が back layer に属する
			if(e.d0 > e.d1){
				back_node = ei0+ej0*(nx+1);
				if(i < yoffset){	// x軸平行エッジ
					nn_back_node = (ei0 == 0 ? ei0 : ei0-1)+ej0*(nx+1);
					l = dw;
				}
				else{			// y軸平行エッジ
					nn_back_node = ei0+(ej0 == 0 ? ej0 : ej0-1)*(nx+1);
					l = dh;
				}
			}
			else{
				back_node = ei1+ej1*(nx+1);
				if(i < yoffset){	// x軸平行エッジ
					nn_back_node = (ei1 == nx ? ei1 : ei1+1)+ej1*(nx+1);
					l = dw;
				}
				else{			// y軸平行エッジ
					nn_back_node = ei1+(ej1 == ny ? ej1 : ej1+1)*(nx+1);
					l = dh;
				}
			}

			// back layerに属するエッジ端点とその隣接ノードのデプス値
			back_node_depth = dgrid[back_node];
			nn_back_node_depth = dgrid[nn_back_node];

			// デプス値を外挿により近似
			back_depth = back_node_depth;
			//back_depth = back_node_depth*((2*l-dx)/l)-nn_back_node_depth*((l-dx)/l);

			// back vertexを設定
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

		if(e.xy == 0){	// x軸平行エッジ
			if(ej0 != 0){	// 下境界ではない
				m_vSSMGrid[ei0+(ej0-1)*nx].edge_vrts[2] = vidx;
				m_vSSMGrid[ei0+(ej0-1)*nx].num_ev++;
				if(bidx != -1){
					m_vSSMGrid[ei0+(ej0-1)*nx].back_vrts[2] = bidx;
					m_vSSMGrid[ei0+(ej0-1)*nx].num_bv++;
				}
			}
			if(ej0 != ny){	// 上境界ではない
				m_vSSMGrid[ei0+ej0*nx].edge_vrts[0] = vidx;
				m_vSSMGrid[ei0+ej0*nx].num_ev++;
				if(bidx != -1){
					m_vSSMGrid[ei0+ej0*nx].back_vrts[0] = bidx;
					m_vSSMGrid[ei0+ej0*nx].num_bv++;
				}
			}
		}
		else{	// y軸平行エッジ
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
 * 三角形メッシュ生成
 * @param[in] nx,ny メッシュ生成用グリッド解像度
 * @param[in] dgrid メッシュグリッド
 * @param[out] polys 三角形メッシュ
 * @param[in] vstart 頂点インデックスの始点
 * @return 生成された三角形メッシュ数
 */
int rxSSMeshCPU::CalMesh(int nx, int ny, vector<rxSSGrid> &grid, vector< vector<int> > &polys, int vstart)
{
	int v[14];
	// ノード頂点
	// 3 - 2
	// |   |
	// 0 - 1
		 
	// エッジ頂点
	// - 6 -
	// 7   5
	// - 4 -
		 
	// エッジ頂点(back vertex)
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

				table_index |= ((v[k] != -1) ? BITS[k] : 0);		// ノード頂点下位4ビット
				table_index |= ((v[k+4] != -1) ? BITS[k]*16 : 0);	// エッジ頂点上位4ビット
			}
			v[12] = -1;
			v[13] = -1;

			int rotation = 0;

			g->table_index0 = table_index;	// デバッグ用

			int update_idx = g->num_ev*5+g->num_nv;
			if(m_FuncTableIndex[update_idx] != 0){
				(this->*m_FuncTableIndex[update_idx])(table_index, rotation, v, g);
			}

			g->table_index1 = table_index;	// デバッグ用

			if(g_MeshTable[table_index][0] > 0){	// グリッド内のメッシュ数が0より大きかったらメッシュ生成
				vector<int> tri;
				tri.resize(3);

				int m = g_MeshTable[table_index][0];	// グリッド内のメッシュ数
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

					m_vMeshGrid.push_back(i+j*nx);	// デバッグ用
					g->mesh[k] = num_mesh;				// デバッグ用

					num_mesh++;
				}

				// デバッグ用
				g->mesh_num = m;
			}
		}
	}

	return num_mesh;
}



/*!
 * 頂点法線計算
 * @param[in] vrts 頂点座標
 * @param[in] nvrts 頂点数
 * @param[in] tris 三角形ポリゴン幾何情報
 * @param[in] ntris 三角形ポリゴン数
 * @param[out] nrms 法線
 */
void rxSSMeshCPU::CalVertexNormals(const vector<Vec3> &vrts, unsigned int nvrts, const vector< vector<int> > &tris, unsigned int ntris, vector<Vec3> &nrms)
{
	unsigned int nnrms = nvrts;
	nrms.resize(nnrms);
	
	// 法線配列の初期化
	for(unsigned int i = 0; i < nnrms; i++){
		nrms[i][0] = 0;
		nrms[i][1] = 0;
		nrms[i][2] = 0;
	}

	// 頂点法線の計算
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

	// 法線正規化
	for(unsigned int i = 0; i < nnrms; i++){
		normalize(nrms[i]);
	}
}



/*!
 * 二分探索(2点間 in 2D)
 * @param[in] v1,v2 探索始点，終点 (v1の方がデプス値が小さい)
 * @param[out] vr 頂点位置
 * @param[in] zmax デプス値差
 * @return v1からの距離
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
 * テーブルインデックス更新なし用関数
 * @param[inout] table_index 元のテーブルインデックス(上位4ビットがエッジ頂点,下位4ビットがノード頂点)
 * @param[inout] vrot 頂点ローテーション数
 * @param[in] v グリッド内の頂点インデックス
 * @param[inout] g メッシュを生成するグリッド
 */
void rxSSMeshCPU::updateTableIndexAny(int &table_index, int &vrot, int v[], rxSSGrid *g)
{
}

/*!
 * パターン0 内部メッシュ用のテーブルインデックス更新
 * @param[inout] table_index 元のテーブルインデックス(上位4ビットがエッジ頂点,下位4ビットがノード頂点)
 * @param[inout] vrot 頂点ローテーション数
 * @param[in] v グリッド内の頂点インデックス
 * @param[inout] g メッシュを生成するグリッド
 */
void rxSSMeshCPU::updateTableIndexE0N4(int &table_index, int &vrot, int v[], rxSSGrid *g)
{
	// グリッド番号により90度回転させる
	table_index -= ((g->i+g->j) & 0x01) ? 1 : 0;
}

/*!
 * パターン1,2,4,8 内部輪郭(輪郭の始点)用のテーブルインデックス更新
 * @param[inout] table_index 元のテーブルインデックス(上位4ビットがエッジ頂点,下位4ビットがノード頂点)
 * @param[inout] vrot 頂点ローテーション数
 * @param[in] v グリッド内の頂点インデックス
 * @param[inout] g メッシュを生成するグリッド
 */
void rxSSMeshCPU::updateTableIndexE1(int &table_index, int &vrot, int v[], rxSSGrid *g)
{
	// ノード頂点で back layer に属する物を探す
	// 頂点があるエッジの端点でデプス値が大きい方が back layer
	int kb = 0;
	for(int k = 0; k < 4; ++k){
		if(v[k+4] != -1){
			int k1 = (k == 3 ? 0 : k+1);
			kb = ((g->node_depth[k] > g->node_depth[k1]) ? BITS[k] : BITS[k1]);
			break;
		}
	}

	kb = (~kb & 0x0F);
	

	// バックレイヤーにあるノード頂点ビットを0にする
	table_index = (table_index & 0xF0)+kb;
}

/*!
 * パターン3,5,6,9,10,12 内部輪郭用のテーブルインデックス更新
 * @param[inout] table_index 元のテーブルインデックス(上位4ビットがエッジ頂点,下位4ビットがノード頂点)
 * @param[inout] vrot 頂点ローテーション数
 * @param[in] v グリッド内の頂点インデックス
 * @param[inout] g メッシュを生成するグリッド
 */
void rxSSMeshCPU::updateTableIndexE2N4(int &table_index, int &vrot, int v[], rxSSGrid *g)
{
	// back vertexのノード頂点ビットを0としたビット列を生成
	int btable = 0;
	int k0 = 0, k1;
	for(int k = 0; k <= 4; ++k){
		k1 = (k0 == 3 ? 0 : k0+1);
		if(v[k0+4] != -1){	// ノード間にエッジ頂点あり
			btable |= (g->node_depth[k0] > g->node_depth[k1]) ? BITS[k0] : BITS[k1];
		}
		else{
			btable |= ((btable & BITS[k0]) ? BITS[k1] : 0);
		}
		k0++;
		if(k0 == 4) k0 = 0;
	}

	// 外部輪郭と区別するために，内部輪郭の場合ビット列に+2する
	btable = (btable+2 & 0x0F);

	table_index = (table_index & 0xF0)+btable;
}

/*!
 * 7,11,13,14 外部/内部輪郭混在用のテーブルインデックス更新
 * @param[inout] table_index 元のテーブルインデックス(上位4ビットがエッジ頂点,下位4ビットがノード頂点)
 * @param[inout] vrot 頂点ローテーション数
 * @param[in] v グリッド内の頂点インデックス
 * @param[inout] g メッシュを生成するグリッド
 */
void rxSSMeshCPU::updateTableIndexE3N23(int &table_index, int &vrot, int v[], rxSSGrid *g)
{
	int btable = 0;	// ノード頂点数2,3の場合，1番目のビットは0(0xxx)
	int pattern = (table_index >> 4);
	int node = g_EdgeTable[pattern][0]-4;
	int R[4];	// 頂点がないエッジから反時計回りにノード頂点番号を並べたリスト
	for(int k = 0; k < 4; ++k){
		R[k] = node++;
		if(node == 4) node = 0;
	}

	int ntable = (table_index & 0x0F);	// ノード頂点ビット列

	// 2番目のビット
	btable |= (((ntable >> R[1]) & 1) ? 4 : 0);

	// 3番目のビット
	btable |= (((ntable >> R[2]) & 1) ? 2 : 0);

	// Rに従いbtableを並び替え
	int btable0 = 0;
	btable0 |= (((ntable >> R[0]) & 1) ? 8 : 0);
	btable0 |= (((ntable >> R[1]) & 1) ? 4 : 0);
	btable0 |= (((ntable >> R[2]) & 1) ? 2 : 0);
	btable0 |= (((ntable >> R[3]) & 1) ? 1 : 0);

	// 並び替えたbtableの下位2ビット-1
	int n0 = (btable0 & 3)-1;
	int n1 = n0+1;

	if(n0 != 1){
		// 追加のback vertexをg_EdgeTable[pattern][1]の位置に追加
		int add_edge = g_EdgeTable[pattern][1]-4;	// 追加するエッジの位置
		Vec3 add_vrt = m_vSSVertex[g->edge_vrts[add_edge]].pos;	// 追加頂点

		// 追加頂点デプス値
		int ref_edge = g_EdgeTable[pattern][0]-4;	// 追加頂点のデプス値を参照するエッジ頂点
		int ref_node = (btable & 4) ? R[1] : R[2];
		if(fabs(m_vSSVertex[g->edge_vrts[ref_edge]].pos[2]-g->node_depth[ref_node]) > m_fSSZmax){
			// edge vertexを使用
			add_vrt[2] = m_vSSVertex[g->edge_vrts[ref_edge]].pos[2];	// 追加頂点のデプス値
		}
		else{
			// back vertexを使用
			if(g->back_vrts[ref_edge] == -1){
				ref_edge = g_EdgeTable[pattern][2]-4;
			}
			add_vrt[2] = m_vSSVertex[g->back_vrts[ref_edge]].pos[2];	// 追加頂点のデプス値
		}


		// 頂点をリストに追加
		m_vSSVertex.push_back(rxSSVertex(add_vrt, 1));
		m_iNumEdgeVrts++;
		int vidx = (int)m_vSSVertex.size()-1;

		// グリッドの頂点リストに追加
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

	// 4番目のビット : ノードのデプス値を比較
	btable |= ((g->node_depth[R[n0]] >= g->node_depth[R[n1]]) ? 0 : 1);

	// テーブルインデックスのエッジ頂点部分を更新
	table_index &= 0xF0;
	table_index |= btable;
}


/*!
 * パターン7,11,13,14 内部輪郭用のテーブルインデックス更新
 *  - さらなるback vertexが必要
 * @param[inout] table_index 元のテーブルインデックス(上位4ビットがエッジ頂点,下位4ビットがノード頂点)
 * @param[inout] vrot 頂点ローテーション数
 * @param[in] v グリッド内の頂点インデックス
 * @param[inout] g メッシュを生成するグリッド
 */
void rxSSMeshCPU::updateTableIndexE3N4(int &table_index, int &vrot, int v[], rxSSGrid *g)
{
	int btable = 8;	// ノード頂点数4の場合，1番目のビットは1(1xxx)
	int pattern = (table_index >> 4);

	// ノードのデプス値を比べる順番の決定
	int node = g_EdgeTable[pattern][0]-4;
	int R[4];	// 頂点がないエッジから反時計回りにノード頂点番号を並べたリスト
	for(int k = 0; k < 4; ++k){
		R[k] = node++;
		if(node == 4) node = 0;
	}

	// ノードのデプス値の大小でビット列を変更(2-4番目のビット)
	for(int k = 0; k < 3; ++k){
		// R[k] > R[k+1]ならば対応するビットを1にする
		if(g->node_depth[R[k]] > g->node_depth[R[k+1]]){
			btable |= BITS[2-k];
		}
	}

	// テーブルインデックスのエッジ頂点部分を更新
	table_index &= 0xF0;
	table_index |= btable;


	// g_EdgeTable[pattern][1]の位置に頂点を追加
	int add_edge = g_EdgeTable[pattern][1]-4;	// 追加するエッジの位置
	Vec3 add_vrt = m_vSSVertex[g->edge_vrts[add_edge]].pos;	// 追加頂点

	// メッシュCの頂点から補間でデプス値を求める
	RXREAL ref_depths[4];
	// ／＼
	// 2ー3
	// |＼|
	// 0ー1

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

	// 追加頂点のデプス値
	add_vrt[2] = 0.5*(ref_depths[2]+ref_depths[3]);

	// 頂点をリストに追加
	m_vSSVertex.push_back(rxSSVertex(add_vrt, 1));
	m_iNumEdgeVrts++;
	int vidx = (int)m_vSSVertex.size()-1;

	// グリッドの頂点リストに追加
	//v[12] = vidx;
	//g->back_vrts[4] = vidx;

	// グリッドの頂点リストに追加
	if(m_vSSVertex[vidx].pos[2] < m_vSSVertex[v[add_edge+8]].pos[2]){
		if(m_vSSVertex[vidx].pos[2] < m_vSSVertex[v[add_edge+4]].pos[2]){
			// front vertexとして挿入
			v[12] = v[add_edge+8];
			g->back_vrts[4] = v[add_edge+8];

			v[add_edge+8] = v[add_edge+4];
			g->back_vrts[add_edge] = v[add_edge+4];

			v[add_edge+4] = vidx;
			g->edge_vrts[add_edge] = vidx;
		}
		else{
			// back vertexとして挿入
			v[12] = v[add_edge+8];
			g->back_vrts[4] = v[add_edge+8];

			v[add_edge+8] = vidx;
			g->back_vrts[add_edge] = vidx;
		}
	}
	else{
		// back-2 vertexとして挿入
		v[12] = vidx;
		g->back_vrts[4] = vidx;
	}					
}

/*!
 * パターン15 外部輪郭用のテーブルインデックス更新
 * @param[inout] table_index 元のテーブルインデックス(上位4ビットがエッジ頂点,下位4ビットがノード頂点)
 * @param[inout] vrot 頂点ローテーション数
 * @param[in] v グリッド内の頂点インデックス
 * @param[inout] g メッシュを生成するグリッド
 */
void rxSSMeshCPU::updateTableIndexE4N2(int &table_index, int &vrot, int v[], rxSSGrid *g)
{
	int ntable = (table_index & 0x0F);	// ノード頂点ビット列
	ntable = (ntable == 5 ? 0 : 15);

	// テーブルインデックスのノード頂点部分を更新
	table_index &= 0xF0;
	table_index |= ntable;
}

/*!
 * パターン15 外部/内部輪郭混在用のテーブルインデックス更新
 * @param[inout] table_index 元のテーブルインデックス(上位4ビットがエッジ頂点,下位4ビットがノード頂点)
 * @param[inout] vrot 頂点ローテーション数
 * @param[in] v グリッド内の頂点インデックス
 * @param[inout] g メッシュを生成するグリッド
 */
void rxSSMeshCPU::updateTableIndexE4N3(int &table_index, int &vrot, int v[], rxSSGrid *g)
{
	int btable = 0;	// ノード頂点数3の場合，1番目のビットは0(0xxx)
	int pattern = (table_index >> 4);
	int ntable = (table_index & 0x0F);	// ノード頂点ビット列

	// 頂点がないノード
	int zero_node = log((double)(~ntable & 0x0F))/log(2.0);

	// ノードのデプス値の大小でビット列を変更(2-4番目のビット)
	for(int k = 0; k < 3; ++k){
		int k0 = g_NodeTable[zero_node][k];
		int k1 = g_NodeTable[zero_node][k+1];

		// k0 > k1ならば対応するビットを1にする
		if(g->node_depth[k0] > g->node_depth[k1]){
			btable |= BITS[2-k];
		}
	}

	// 頂点ローテーション
	vrot = zero_node;

	// テーブルインデックスのノード頂点部分を更新
	table_index &= 0xF0;
	table_index |= btable;
}

/*!
 * パターン15 内部輪郭用のテーブルインデックス更新
 *  - 追加のback vertexが2つ必要
 * @param[inout] table_index 元のテーブルインデックス(上位4ビットがエッジ頂点,下位4ビットがノード頂点)
 * @param[inout] vrot 頂点ローテーション数
 * @param[in] v グリッド内の頂点インデックス
 * @param[inout] g メッシュを生成するグリッド
 */
void rxSSMeshCPU::updateTableIndexE4N4(int &table_index, int &vrot, int v[], rxSSGrid *g)
{
	int btable = 8;	// ノード頂点数4の場合，1番目のビットは1(1xxx)
	int pattern = (table_index >> 4);
	int ntable = (table_index & 0x0F);	// ノード頂点ビット列

	// デプス値がもっとも大きいノード
	int zero_node = 0;
	double max_depth = 0.0;
	for(int k = 1; k < 4; ++k){
		if(g->node_depth[k] > max_depth){
			max_depth = g->node_depth[k];
			zero_node = k;
		}
	}

	// ノードのデプス値の大小でビット列を変更(2-4番目のビット)
	for(int k = 0; k < 3; ++k){
		int k0 = g_NodeTable[zero_node][k];
		int k1 = g_NodeTable[zero_node][k+1];

		// k0 > k1ならば対応するビットを1にする
		if(g->node_depth[k0] > g->node_depth[k1]){
			btable |= BITS[2-k];
		}
	}

	// 頂点ローテーション
	vrot = zero_node;

	// テーブルインデックスのエッジ頂点部分を更新
	table_index &= 0xF0;
	table_index |= btable;

	// 
	// g_NodeTable[zero_node][5,6]の位置に頂点(back-2 vertex)を追加
	int add_edge1 = g_NodeTable[zero_node][5]-4;	// 追加するエッジの位置
	int add_edge2 = g_NodeTable[zero_node][6]-4;	// 追加するエッジの位置
	Vec3 add_vrt1 = m_vSSVertex[g->edge_vrts[add_edge1]].pos;	// 追加頂点
	Vec3 add_vrt2 = m_vSSVertex[g->edge_vrts[add_edge2]].pos;	// 追加頂点

	// g_NodeTable[zero_node][4,7]の位置のback vertexのデプスを設定
	int ref_edge1 = g_NodeTable[zero_node][4]-4;	// 参照するエッジの位置
	int ref_edge2 = g_NodeTable[zero_node][7]-4;	// 参照するエッジの位置
	add_vrt1[2] = m_vSSVertex[g->back_vrts[ref_edge1]].pos[2];
	add_vrt2[2] = m_vSSVertex[g->back_vrts[ref_edge2]].pos[2];

	// 頂点をリストに追加
	m_vSSVertex.push_back(rxSSVertex(add_vrt1, 1));
	m_vSSVertex.push_back(rxSSVertex(add_vrt2, 1));
	m_iNumEdgeVrts += 2;
	int vidx1 = (int)m_vSSVertex.size()-2;
	int vidx2 = (int)m_vSSVertex.size()-1;

	// グリッドの頂点リストに追加
	v[12] = vidx1;
	v[13] = vidx2;
	g->back_vrts[4] = vidx1;
	g->back_vrts[5] = vidx2;
}


/*!
 * デプスマップの取得
 * @return デプスマップ
 */
RXREAL* rxSSMeshCPU::GetDepthMap(void)
{ 
	return &m_vSSDMap[0];
}


/*!
 * メッシュ生成用グリッドの取得
 * @param[in] idx グリッドインデックス
 * @return メッシュ生成用グリッド(rxSSGrid)
 */
rxSSGrid rxSSMeshCPU::GetMeshGrid(int idx)
{
	return m_vSSMGrid[idx];
}

/*!
 * スクリーン空間での頂点情報を取得
 * @return 頂点情報(rxSSVertex)
 */
rxSSVertex* rxSSMeshCPU::GetSSVertex(void)
{
	return &m_vSSVertex[0];
}


/*!
 * 輪郭エッジのOpenGL描画
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
 * 輪郭エッジのOpenGL描画
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
 * Screen Spaceでの頂点のOpenGL描画
 * @param[in] node_color ノード頂点の描画色
 * @param[in] edge_color エッジ頂点の描画色
 */
void rxSSMeshCPU::DrawSSVertex(Vec3 node_color, Vec3 edge_color)
{
	glBegin(GL_POINTS);
	// ノード頂点
	glColor3dv(node_color);
	for(int i = 0; i < m_iNumNodeVrts; ++i){
		Vec3 x = m_vSSVertex[i].pos;
		glVertex3d(x[0], x[1], 0.02);
	}
	// 輪郭エッジ頂点
	glColor3dv(edge_color);
	for(int i = m_iNumNodeVrts; i < (int)m_vSSVertex.size(); ++i){
		Vec3 x = m_vSSVertex[i].pos;
		glVertex3d(x[0], x[1], 0.025);
	}
	glEnd();
}

/*!
 * グリッド内頂点のOpenGL描画
 * @param[in] grid グリッドインデックス
 * @param[in] colors 頂点の色
 */
void rxSSMeshCPU::DrawMeshGrid(int grid, const Vec3 colors[])
{
	// グリッド内の頂点
	rxSSGrid g = GetMeshGrid(grid);
	for(int i = 0; i < 4; ++i){
		glPointSize(5.0);
		// ノード頂点
		if(g.node_vrts[i] != -1){
			Vec3 x = m_vSSVertex[g.node_vrts[i]].pos;
			glColor3dv(colors[0]);
			glBegin(GL_POINTS);
			glVertex3d(x[0], x[1], 0.03);
			glEnd();
		}
		// エッジ頂点
		if(g.edge_vrts[i] != -1){
			Vec3 x = m_vSSVertex[g.edge_vrts[i]].pos;
			glColor3dv(colors[1]);
			glBegin(GL_POINTS);
			glVertex3d(x[0], x[1], 0.04);
			glEnd();
		}
		glPointSize(8.0);
		// エッジ頂点(back vertex)
		if(g.back_vrts[i] != -1){
			Vec3 x = m_vSSVertex[g.back_vrts[i]].pos;
			glColor3dv(colors[2]);
			glBegin(GL_POINTS);
			glVertex3d(x[0], x[1], 0.039);
			glEnd();
		}
	}
	// エッジ頂点(back-2 vertex)
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
 * 密度場の描画
 * @param[in] minpos[2],maxpos[2] 描画領域の範囲
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
 * グリッドに関する情報の出力
 * @param[in] grid グリッドインデックス(i+j*nx)
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
 * グリッドに関する情報の出力(頂点情報含む)
 * @param[in] grid グリッドインデックス(i+j*nx)
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
