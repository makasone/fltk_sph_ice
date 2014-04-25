/*!
  @file rx_ssm_gpu.cpp
	
  @brief Screen Space Mesh作成
		- M. Muller et al., Screen space meshes, SCA2007, 2007. 

  @author Makoto Fujisawa
  @date 2011-05
*/
// FILE --rx_ssm_gpu.cpp--


//-----------------------------------------------------------------------------
// インクルードファイル
//-----------------------------------------------------------------------------
//#include "rx_gltexture.h"
#include <GL/glew.h>
#include <GL/glut.h>

#include "rx_ssm.h"
#include "rx_cu_funcs.cuh"

#include "rx_mc.h"

//-----------------------------------------------------------------------------
// 定数・変数
//-----------------------------------------------------------------------------
const int BITS[4] = { 1, 2, 4, 8 };



//-----------------------------------------------------------------------------
// rxSSMeshGPUクラスの実装
//-----------------------------------------------------------------------------
/*!
 * コンストラクタ
 */
rxSSMeshGPU::rxSSMeshGPU(double zmax, double h, double r, int n_filter, int n_iter)
	: rxSSMesh(zmax, h, r, n_filter, n_iter)
{
	// MARK:コンストラクタ
	m_dSSDMap = 0;
	m_dSSPrtsCen = 0;
	m_dSSPrtsRad = 0;

	m_dSSTriNum = 0;
	m_dSSTriNumScan = 0;

	m_dSSVertex = 0;
	m_dSSMGrid = 0;

	m_hSSPrtsCen = 0;
	m_hSSPrtsRad = 0;

	m_hSSEdge = 0;
	m_hSSVertex = 0;

	m_vBinomials = CalBinomialsForFilter(RX_MAX_FILTER_SIZE);
	CuInitTable();
}


/*!
 * デストラクタ
 */
rxSSMeshGPU::~rxSSMeshGPU()
{
	CuCleanTable();

	if(m_dSSDMap) CuFreeArray(m_dSSDMap);
	if(m_dSSPrtsCen) CuFreeArray(m_dSSPrtsCen);
	if(m_dSSPrtsRad) CuFreeArray(m_dSSPrtsRad);

	CuCleanVPacke(m_dSSEdge);
	CuCleanVPackf(m_dSSNodeVrts);
	CuCleanVPackf(m_dSSEdgeFrontVrts);
	CuCleanVPackf(m_dSSEdgeBackVrts);
	CuCleanVPackf(m_dSSEdgeBack2Vrts);

	if(m_dSSTriNum) CuFreeArray(m_dSSTriNum);
	if(m_dSSTriNumScan) CuFreeArray(m_dSSTriNumScan);

	if(m_dSSVertex) CuFreeArray(m_dSSVertex);
	if(m_dSSMGrid) CuFreeArray(m_dSSMGrid);

	if(m_hSSPrtsCen) delete [] m_hSSPrtsCen;
	if(m_hSSPrtsRad) delete [] m_hSSPrtsRad;

	if(m_hSSEdge) delete [] m_hSSEdge;
	if(m_hSSVertex) delete [] m_hSSVertex;
}

/*!
 * マップやグリッド配列のサイズを変更
 * @param[in] W,H 画面解像度
 * @param[in] spacing デプスマップのサンプリング間隔
 */
void rxSSMeshGPU::Resize(int W, int H, int spacing)
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
	// GPUメモリ確保
	//
	unsigned int size  = (m_iNgx+1)*(m_iNgy+1);
	unsigned int size_e = m_iNgx*(m_iNgy+1)+(m_iNgx+1)*m_iNgy;
	unsigned int size_m = m_iNgx*m_iNgy;

	// デプスマップ
	if(m_dSSDMap) CuFreeArray(m_dSSDMap);
	CuAllocateArray((void**)&m_dSSDMap, size*sizeof(RXREAL));

	// メッシュ生成用グリッド
	if(m_dSSMGrid) CuFreeArray(m_dSSMGrid);
	CuAllocateArray((void**)&m_dSSMGrid, size_m*sizeof(rxSSGridG));

	// ノード頂点
	CuInitVPackf(m_dSSNodeVrts, size);

	// 輪郭エッジ
	CuInitVPacke(m_dSSEdge, size_e);	// 前面
	
	// エッジ頂点
	CuInitVPackf(m_dSSEdgeFrontVrts, size_e);	// 前面
	CuInitVPackf(m_dSSEdgeBackVrts, size_e);		// 背面
	CuInitVPackf(m_dSSEdgeBack2Vrts, size_e);	// 最背面

	// メッシュ数情報
	if(m_dSSTriNum) CuFreeArray(m_dSSTriNum);
	CuAllocateArray((void**)&m_dSSTriNum, size_m*sizeof(unsigned int));
	if(m_dSSTriNumScan) CuFreeArray(m_dSSTriNumScan);
	CuAllocateArray((void**)&m_dSSTriNumScan, size_m*sizeof(unsigned int));
}


/*!
 * GPU側のパラメータファイルを更新
 */
void rxSSMeshGPU::updateParams(int W, int H, const rxMatrix4 &P, const rxMatrix4 &PMV)
{
	m_ssm_params.W = W;
	m_ssm_params.H = H;

	for(int j = 0; j < 4; j++){
		m_ssm_params.PMV[j].x = (float)PMV(0,j);
		m_ssm_params.PMV[j].y = (float)PMV(1,j);
		m_ssm_params.PMV[j].z = (float)PMV(2,j);
		m_ssm_params.PMV[j].w = (float)PMV(3,j);
	}

	m_ssm_params.Tr.x = (float)(m_fPrtRad*W*sqrt(P(0, 0)*P(0, 0)+P(0, 1)*P(0, 1)+P(0, 2)*P(0, 2)));
	m_ssm_params.Tr.y = (float)(m_fPrtRad*H*sqrt(P(1, 0)*P(1, 0)+P(1, 1)*P(1, 1)+P(1, 2)*P(1, 2)));
	m_ssm_params.Tr.z = (float)(m_fPrtRad*sqrt(P(2, 0)*P(2, 0)+P(2, 1)*P(2, 1)+P(2, 2)*P(2, 2)));

	//cout << "tr(gpu) = " << m_ssm_params.Tr.x << ", " << m_ssm_params.Tr.y << ", " << m_ssm_params.Tr.z << endl;


	m_ssm_params.PrtRad = (RXREAL)(m_fPrtRad);
	m_ssm_params.Spacing = (RXREAL)(m_fSpacing);
	m_ssm_params.Zmax = (RXREAL)(m_fSSZmax);
	m_ssm_params.Nfilter = m_iNfilter;
	m_ssm_params.Niters = m_iNiters;
	m_ssm_params.Ngx = m_iNgx;
	m_ssm_params.Ngy = m_iNgy;

	CuSetSSMParameters(&m_ssm_params);
}

/*!
 * デプス値を画像で保存
 * @param[in] fn 保存ファイル名
 * @param[in] nx,ny デプスマップの解像度
 * @param[in] dm デプスマップ
void SaveDepthMapf(string fn, int nx, int ny, RXREAL *dm)
{
	vector<unsigned char> img;
	img.resize(nx*ny*3);

	for(int j = 0; j < ny; ++j){
		for(int i = 0; i < nx; ++i){
			unsigned int ud = RX_DEPTH2COLOR(dm[j*nx+i]);

			img[3*(j*nx+i)+0] = ud;
			img[3*(j*nx+i)+1] = ud;
			img[3*(j*nx+i)+2] = ud;
		}
	}

	SavePngFile(fn, &img[0], nx, ny, 3);
}
 */

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
void rxSSMeshGPU::CreateMesh(double *proj, double *modelview, int W, int H, vector<Vec3> &prts, int pnum, 
							 vector<Vec3> &vrts, vector<Vec3> &nrms, vector< vector<int> > &polys, 
							 int filtering, int debug_output)
{
	// スクリーンスペースでのパーティクル
	unsigned int size_p = pnum*3;
	unsigned int mem_size_p = sizeof(RXREAL)*size_p;

	RXREAL *dppos = 0;
	CuAllocateArray((void**)&dppos, mem_size_p);

	if(m_hSSPrtsCen) delete [] m_hSSPrtsCen;
	m_hSSPrtsCen = new RXREAL[size_p];

	// パーティクルデータをデバイスメモリに転送
	for(int i = 0; i < pnum; ++i){
		m_hSSPrtsCen[3*i+0] = (RXREAL)prts[i][0];
		m_hSSPrtsCen[3*i+1] = (RXREAL)prts[i][1];
		m_hSSPrtsCen[3*i+2] = (RXREAL)prts[i][2];
	}
	CuCopyArrayToDevice(dppos, m_hSSPrtsCen, 0, mem_size_p);

	// メッシュ初期化
	vrts.clear();
	nrms.clear();
	polys.clear();

	//
	// メッシュ生成
	//
	RXREAL *dvrt = 0;
	unsigned int *dtri = 0;
	CreateMeshD(proj, modelview, W, H, dppos, 0, pnum, 3, dvrt, dtri, filtering, debug_output);
	if(dppos) CuFreeArray(dppos);

	RXREAL *dnrm = 0;
	int num_vrts = m_iNumNodeVrts+m_iNumEdgeVrts;
	calVertexNormalD(dvrt, num_vrts, dtri, m_iNumMesh, dnrm);


	//
	// GPU -> CPU のデータ転送
	// 
	// 頂点
	RXREAL *hVrts = new RXREAL[num_vrts*3];
	CuCopyArrayFromDevice(hVrts, dvrt, 0, num_vrts*3*sizeof(RXREAL));

	vrts.resize(num_vrts);
	for(int i = 0; i < num_vrts; ++i){
		vrts[i] = Vec3(hVrts[3*i+0], hVrts[3*i+1], hVrts[3*i+2]);
	}

	delete [] hVrts;

	if(dnrm != 0){
		RXREAL *hNrms = new RXREAL[num_vrts*3];
		CuCopyArrayFromDevice(hNrms, dnrm, 0, num_vrts*3*sizeof(RXREAL));

		nrms.resize(num_vrts);
		for(int i = 0; i < num_vrts; ++i){
			nrms[i] = Vec3(hNrms[3*i+0], hNrms[3*i+1], hNrms[3*i+2]);
		}

		delete [] hNrms;
	}

	// メッシュ
	unsigned int *hTris = new unsigned int[m_iNumMesh*3];
	CuCopyArrayFromDevice(hTris, dtri, 0, m_iNumMesh*3*sizeof(unsigned int));

	polys.resize(m_iNumMesh);
	for(int i = 0; i < m_iNumMesh; ++i){
		polys[i].resize(3);
		for(int j = 0; j < 3; ++j){
			polys[i][j] = (int)hTris[3*i+j];
		}
	}

	delete [] hTris;



	if(dvrt) CuFreeArray(dvrt);
	if(dtri) CuFreeArray(dtri);
	m_vSSVertex.clear();

	if(m_hSSVertex) delete [] m_hSSVertex;
	m_hSSVertex = new RXREAL[num_vrts*3];
	CuCopyArrayFromDevice(m_hSSVertex, m_dSSVertex, 0, num_vrts*3*sizeof(float));
	for(int i = 0; i < m_iNumNodeVrts; ++i){
		Vec3 p;
		p[0] = m_hSSVertex[3*i+0];
		p[1] = m_hSSVertex[3*i+1];
		p[2] = m_hSSVertex[3*i+2];
		m_vSSVertex.push_back(rxSSVertex(p, 0));
	}
	for(int i = m_iNumNodeVrts; i < num_vrts; ++i){
		Vec3 p;
		p[0] = m_hSSVertex[3*i+0];
		p[1] = m_hSSVertex[3*i+1];
		p[2] = m_hSSVertex[3*i+2];
		m_vSSVertex.push_back(rxSSVertex(p, 1));
	}
}

/*!
 * スクリーンスペースメッシュ生成
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
void rxSSMeshGPU::CreateMesh(double *proj, double *modelview, int W, int H, vector<Vec3> &prts, int pnum, 
							 vector<Vec3> &vrts, vector< vector<int> > &polys, int filtering, int debug_output)
{
	unsigned int size_p = pnum*3;
	unsigned int mem_size_p = sizeof(RXREAL)*size_p;

	RXREAL *dppos = 0;
	CuAllocateArray((void**)&dppos, mem_size_p);

	if(m_hSSPrtsCen) delete [] m_hSSPrtsCen;
	m_hSSPrtsCen = new RXREAL[size_p];

	// パーティクルデータをデバイスメモリに転送
	for(int i = 0; i < pnum; ++i){
		m_hSSPrtsCen[3*i+0] = (RXREAL)prts[i][0];
		m_hSSPrtsCen[3*i+1] = (RXREAL)prts[i][1];
		m_hSSPrtsCen[3*i+2] = (RXREAL)prts[i][2];
	}
	CuCopyArrayToDevice(dppos, m_hSSPrtsCen, 0, mem_size_p);

	// メッシュ初期化
	vrts.clear();
	polys.clear();

	//
	// メッシュ生成
	//
	RXREAL *dvrt = 0;
	unsigned int *dtri = 0;
	CreateMeshD(proj, modelview, W, H, dppos, 0, pnum, 3, dvrt, dtri, filtering, debug_output);
	if(dppos) CuFreeArray(dppos);

	//
	// GPU -> CPU のデータ転送
	// 
	// 頂点
	int num_vrts = m_iNumNodeVrts+m_iNumEdgeVrts;
	RXREAL *hVrts = new RXREAL[num_vrts*3];
	CuCopyArrayFromDevice(hVrts, dvrt, 0, num_vrts*3*sizeof(RXREAL));

	vrts.resize(num_vrts);
	for(int i = 0; i < num_vrts; ++i){
		vrts[i] = Vec3(hVrts[3*i+0], hVrts[3*i+1], hVrts[3*i+2]);
	}

	delete [] hVrts;

	// メッシュ
	unsigned int *hTris = new unsigned int[m_iNumMesh*3];
	CuCopyArrayFromDevice(hTris, dtri, 0, m_iNumMesh*3*sizeof(unsigned int));

	polys.resize(m_iNumMesh);
	for(int i = 0; i < m_iNumMesh; ++i){
		polys[i].resize(3);
		for(int j = 0; j < 3; ++j){
			polys[i][j] = (int)hTris[3*i+j];
		}
	}

	delete [] hTris;
	

	if(dvrt) CuFreeArray(dvrt);
	if(dtri) CuFreeArray(dtri);
	m_vSSVertex.clear();

	if(m_hSSVertex) delete [] m_hSSVertex;
	m_hSSVertex = new RXREAL[num_vrts*3];
	CuCopyArrayFromDevice(m_hSSVertex, m_dSSVertex, 0, num_vrts*3*sizeof(float));
	for(int i = 0; i < m_iNumNodeVrts; ++i){
		Vec3 p;
		p[0] = m_hSSVertex[3*i+0];
		p[1] = m_hSSVertex[3*i+1];
		p[2] = m_hSSVertex[3*i+2];
		m_vSSVertex.push_back(rxSSVertex(p, 0));
	}
	for(int i = m_iNumNodeVrts; i < num_vrts; ++i){
		Vec3 p;
		p[0] = m_hSSVertex[3*i+0];
		p[1] = m_hSSVertex[3*i+1];
		p[2] = m_hSSVertex[3*i+2];
		m_vSSVertex.push_back(rxSSVertex(p, 1));
	}

}

/*!
 * スクリーンスペースメッシュ生成
 *  - パーティクルがデバイスメモリに格納されている場合
 *  - VBOに出力
 * @param[in] proj OpenGL透視投影行列
 * @param[in] modelview OpenGLモデルビュー行列
 * @param[in] W,H 画面解像度
 * @param[in] prts パーティクル座標
 * @param[in] pnum パーティクル数
 * @param[out] uvrt_vbo メッシュ頂点列VBO
 * @param[out] unrm_vbo メッシュ頂点列VBO
 * @param[out] utri_vbo メッシュ頂点列VBO
 * @param[in] filtering フィルタリングフラグ(デプス値0x01, 輪郭0x02)
 * @param[in] debug_output 結果の画面出力の有無
 */
void rxSSMeshGPU::CreateMeshVBO(double *proj, double *modelview, int W, int H, RXREAL *dppos, RXREAL *dprad, int pnum, int pdim, 
								GLuint &uvrt_vbo, GLuint &unrm_vbo, GLuint &utri_vbo, int filtering, int debug_output)
{
	RXREAL *dvrt = 0;
	RXREAL *dnrm = 0;
	unsigned int *dtri = 0;

	// メッシュ生成
	CreateMeshD(proj, modelview, W, H, dppos, dprad, pnum, pdim, dvrt, dnrm, dtri, filtering, debug_output);

	int num_vrts = m_iNumNodeVrts+m_iNumEdgeVrts;
	if(num_vrts){
		AssignArrayBuffers(num_vrts, 3, uvrt_vbo, unrm_vbo, utri_vbo);

		// 頂点VBO
		glBindBuffer(GL_ARRAY_BUFFER, uvrt_vbo);
		RXREAL *vrt_ptr = (RXREAL*)glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);

		CuCopyArrayFromDevice(vrt_ptr, dvrt, 0, num_vrts*3*sizeof(RXREAL));

		glUnmapBuffer(GL_ARRAY_BUFFER);
		glBindBuffer(GL_ARRAY_BUFFER, 0);

		// 法線VBO
		glBindBuffer(GL_ARRAY_BUFFER, unrm_vbo);
		RXREAL *nrm_ptr = (RXREAL*)glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);

		CuCopyArrayFromDevice(nrm_ptr, dnrm, 0, num_vrts*3*sizeof(RXREAL));

		glUnmapBuffer(GL_ARRAY_BUFFER);
		glBindBuffer(GL_ARRAY_BUFFER, 0);

		// メッシュVBO
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, utri_vbo);
		unsigned int *tri_ptr = (unsigned int*)glMapBuffer(GL_ELEMENT_ARRAY_BUFFER, GL_WRITE_ONLY);

		CuCopyArrayFromDevice(tri_ptr, dtri, 0, m_iNumMesh*3*sizeof(unsigned int));

		glUnmapBuffer(GL_ELEMENT_ARRAY_BUFFER);
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
	}

	if(dvrt) CuFreeArray(dvrt);
	if(dnrm) CuFreeArray(dnrm);
	if(dtri) CuFreeArray(dtri);
}


/*!
 * スクリーンスペースメッシュ生成
 *  - パーティクルがデバイスメモリに格納されている場合
 *  - 法線も計算
 * @param[in] proj OpenGL透視投影行列
 * @param[in] modelview OpenGLモデルビュー行列
 * @param[in] W,H 画面解像度
 * @param[in] prts パーティクル座標
 * @param[in] pnum パーティクル数
 * @param[out] dvrt 頂点列(デバイスメモリ)
 * @param[out] dnrm 頂点法線列(デバイスメモリ)
 * @param[out] dtri メッシュ列(デバイスメモリ)
 * @param[in] filtering フィルタリングフラグ(デプス値0x01, 輪郭0x02)
 * @param[in] debug_output 結果の画面出力の有無
 */
void rxSSMeshGPU::CreateMeshD(double *proj, double *modelview, int W, int H, RXREAL *dppos, RXREAL *dprad, int pnum, int pdim, 
							  RXREAL* &dvrt, RXREAL* &dnrm, unsigned int* &dtri, int filtering, int debug_output)
{
	// メッシュ生成
	CreateMeshD(proj, modelview, W, H, dppos, dprad, pnum, pdim, dvrt, dtri, filtering, debug_output);

	int num_vrts = m_iNumNodeVrts+m_iNumEdgeVrts;
	if(num_vrts){
		// 頂点法線計算
		calVertexNormalD(dvrt, num_vrts, dtri, m_iNumMesh, dnrm);
	}

}

/*!
 * スクリーンスペースメッシュ生成
 *  - パーティクルがデバイスメモリに格納されている場合
 * @param[in] proj OpenGL透視投影行列
 * @param[in] modelview OpenGLモデルビュー行列
 * @param[in] W,H 画面解像度
 * @param[in] prts パーティクル座標
 * @param[in] pnum パーティクル数
 * @param[out] dvrt メッシュ頂点列(デバイスメモリ)
 * @param[out] dtri メッシュ列(デバイスメモリ)
 * @param[in] filtering フィルタリングフラグ(デプス値0x01, 輪郭0x02)
 * @param[in] debug_output 結果の画面出力の有無
 */
void rxSSMeshGPU::CreateMeshD(double *proj, double *modelview, int W, int H, RXREAL *dppos, RXREAL *dprad, int pnum, int pdim, 
							  RXREAL* &dvrt, unsigned int* &dtri, int filtering, int debug_output)
{
	// MARK:CreateMeshD
	int nx = W/m_fSpacing;
	int ny = H/m_fSpacing;

	if(nx != m_iNgx || ny != m_iNgy || !m_dSSDMap){
		Resize(W, H, m_fSpacing);
	}

	if(debug_output){
		cout << "window : " << W << " x " << H << endl;
		cout << "depth map : " << m_iNgx+1 << " x " << m_iNgy+1 << endl;
		cout << "grid : " << m_iNgx << " x " << m_iNgy << endl;
	}

	//
	// GPUメモリ確保
	//
	unsigned int size_p = pnum*pdim;

	// スクリーンスペースでのパーティクル
	if(m_dSSPrtsCen) CuFreeArray(m_dSSPrtsCen);
	CuAllocateArray((void**)&m_dSSPrtsCen, pnum*pdim*sizeof(RXREAL));
	if(m_dSSPrtsRad) CuFreeArray(m_dSSPrtsRad);
	CuAllocateArray((void**)&m_dSSPrtsRad, pnum*3*sizeof(RXREAL));

	CuCopyArrayD2D(m_dSSPrtsCen, dppos, pnum*pdim*sizeof(RXREAL));
	if(dprad == 0){
		CuInitFloatArray(m_dSSPrtsRad, m_fPrtRad, pnum*3);
	}
	else{
		CuInitRadArray(m_dSSPrtsRad, dprad, pnum);
	}

	
	// 透視投影変換行列
	const rxMatrix4 P = GetMatrixGL(proj);

	// モデルビュー変換行列
	const rxMatrix4 MV = GetMatrixGL(modelview);

	rxMatrix4 PMV = P*MV;

	//updateParams(W, H, P, PMV);


	//
	// パーティクル座標と半径を透視投影変換してデプスマップを生成
	//
	float tr[3];
	tr[0] = (float)(W*sqrt(P(0, 0)*P(0, 0)+P(0, 1)*P(0, 1)+P(0, 2)*P(0, 2)));
	tr[1] = (float)(H*sqrt(P(1, 0)*P(1, 0)+P(1, 1)*P(1, 1)+P(1, 2)*P(1, 2)));
	tr[2] = (float)(sqrt(P(2, 0)*P(2, 0)+P(2, 1)*P(2, 1)+P(2, 2)*P(2, 2)));

	float pmvf[16];
	int c = 0;
	for(int i = 0; i < 4; i++){
		for(int j = 0; j < 4; j++){
			pmvf[c++] = (float)PMV(i, j);
		}
	}

	CuInitDepthMap(m_dSSDMap, m_iNgx, m_iNgy);
	CuCreateDepthMap(m_dSSDMap, m_dSSPrtsCen, m_dSSPrtsRad, pnum, pdim, tr, pmvf, W, H, m_fDmx, m_fDmy, m_iNgx, m_iNgy);


	// 
	// デプスマップの平滑化
	//
	if(filtering & 0x01){
		CuDepthSmoothing(m_dSSDMap, m_iNgx+1, m_iNgy+1, m_iNfilter, &m_vBinomials[0], m_fSSZmax);
	}


	//
	// ノード頂点生成
	//
	m_iNumNodeVrts = 0;
	CuCalNodeVertex(m_dSSDMap, m_iNgx, m_iNgy, m_fDmx, m_fDmy, m_fSSZmax, 
					m_dSSNodeVrts, m_dSSMGrid, m_iNumNodeVrts);


	//
	// 輪郭エッジの検出とエッジ頂点生成
	//
	int num_edge = 0, num_front_vrts = 0, num_back_vrts = 0;
	CuDetectSilhouetteEdgeVertex(m_dSSDMap, m_iNgx, m_iNgy, m_fDmx, m_fDmy, m_fSSZmax, 
								 m_dSSPrtsCen, m_dSSPrtsRad, pnum, pdim, W, H, 
								 m_dSSEdge, m_dSSEdgeFrontVrts, m_dSSEdgeBackVrts, m_dSSMGrid, 
								 m_iNumNodeVrts, num_edge, num_front_vrts, num_back_vrts);

	m_iNumEdgeVrts = num_front_vrts+num_back_vrts;
	int num_vrts = m_iNumNodeVrts+m_iNumEdgeVrts;

	// すべての頂点をまとめる
	if(m_dSSVertex) CuFreeArray(m_dSSVertex);
	CuAllocateArray((void**)&m_dSSVertex, sizeof(float)*3*num_vrts);

	int offset1 = m_iNumNodeVrts*3;
	int offset2 = offset1+num_front_vrts*3;
	CuCopyArrayD2D(m_dSSVertex, m_dSSNodeVrts.dCompactedPos, m_iNumNodeVrts*3*sizeof(float));
	CuCopyArrayD2D(m_dSSVertex+offset1, m_dSSEdgeFrontVrts.dCompactedPos, num_front_vrts*3*sizeof(float));
	if(num_back_vrts){
		CuCopyArrayD2D(m_dSSVertex+offset2, m_dSSEdgeBackVrts.dCompactedPos, num_back_vrts*3*sizeof(float));
	}


	// 
	// メッシュ生成
	// 
	int num_mesh = 0, num_back2 = 0;
	if(num_vrts){
		CuCalMesh(m_dSSMGrid, m_iNgx, m_iNgy, m_fDmx, m_fDmy, m_fSSZmax, m_dSSTriNum, m_dSSTriNumScan, 
				  m_dSSEdgeBack2Vrts, m_dSSVertex, num_vrts, num_back2, dtri, num_mesh);
		m_iNumMesh = num_mesh;

		if(num_back2){
			m_iNumEdgeVrts += num_back2;

			// back-2 vertexを追加
			RXREAL *dTmp;
			CuAllocateArray((void**)&dTmp, sizeof(RXREAL)*3*(num_vrts+num_back2));

			CuCopyArrayD2D(dTmp, m_dSSVertex, num_vrts*3*sizeof(RXREAL));
			CuCopyArrayD2D(dTmp+num_vrts*3, m_dSSEdgeBack2Vrts.dCompactedPos, num_back2*3*sizeof(RXREAL));
		
			CuFreeArray(m_dSSVertex);
			m_dSSVertex = dTmp;
			num_vrts += num_back2;
		}
	}
	else{
		m_iNumMesh = 0;
	}

	if(debug_output){
		cout << "the number of node vertices = " << m_iNumNodeVrts << endl;
		cout << "the number of silhouette edges = " << num_edge << endl;
		cout << "the number of edge vertices = " << num_front_vrts << ", " << num_back_vrts << ", " << num_back2 << endl;
		cout << "the number of vertices = " << num_vrts << endl;
		cout << "the number of mesh = " << num_mesh << endl;
	}


	// 
	// 輪郭の平滑化
	//
	if((filtering & 0x02) && num_vrts && m_iNumMesh){
		CuSilhouetteSmoothing(m_dSSVertex, num_vrts, m_iNumNodeVrts, dtri, num_mesh, 1);
		//ApplySilhoutteSmoothing(m_vSSVertex, polys, m_iNiters);
	}


	//
	// 3D空間にメッシュを戻す
	//
	if(num_vrts){
		// 透視投影変換逆行列
		rxMatrix4 Q = P.Inverse();

		// モデルビュー変換逆行列
		rxMatrix4 IMV = MV.Inverse();

		rxMatrix4 IMVQ = IMV*Q;

		float qf[4];
		float imvqf[16];
		int cq = 0;
		for(int i = 0; i < 4; i++){
			for(int j = 0; j < 4; j++){
				imvqf[cq++] = (float)IMVQ(i, j);
			}
			qf[i] = (float)Q(3, i);
		}

		CuAllocateArray((void**)&dvrt, num_vrts*3*sizeof(RXREAL));
		CuTransformBack(m_dSSVertex, num_vrts, imvqf, qf, W, H, dvrt);
	}
}


/*!
 * 頂点法線計算
 * @param[in] vrts 頂点座標
 * @param[in] nvrts 頂点数
 * @param[in] tris 三角形ポリゴン幾何情報
 * @param[in] ntris 三角形ポリゴン数
 * @param[out] nrms 法線
 */
void rxSSMeshGPU::CalVertexNormals(const vector<Vec3> &vrts, unsigned int nvrts, const vector< vector<int> > &tris, unsigned int ntris, vector<Vec3> &nrms)
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
 * 頂点法線計算
 * @param[in] dvrt 頂点座標
 * @param[in] nvrts 頂点数
 * @param[in] dtri 三角形ポリゴン幾何情報
 * @param[in] ntris 三角形ポリゴン数
 * @param[out] nrms 法線
 */
void rxSSMeshGPU::calVertexNormalD(RXREAL* dvrt, unsigned int nvrts, unsigned int* dtri, unsigned int ntris, RXREAL* &dnrms)
{
	// MARK:calVertexNormalD

	if(dnrms) CuFreeArray(dnrms);
	CuAllocateArray((void**)&dnrms, sizeof(float)*3*nvrts);

	CuCalVertexNormal(dvrt, nvrts, dtri, ntris, dnrms);
}

/*!
 * デプスマップの取得
 * @return デプスマップ
 */
RXREAL* rxSSMeshGPU::GetDepthMap(void)
{ 
	unsigned int size  = (m_iNgx+1)*(m_iNgy+1);

	if((int)m_vSSDMap.size() != size){
		m_vSSDMap.clear();
		m_vSSDMap.resize(size);
	}

	unsigned int mem_size  = sizeof(RXREAL)*size;
	CuCopyArrayFromDevice(&m_vSSDMap[0], m_dSSDMap, 0, mem_size);

	return &m_vSSDMap[0];
}


/*!
 * メッシュ生成用グリッドの取得
 * @param[in] idx グリッドインデックス
 * @return メッシュ生成用グリッド(rxSSGrid)
 */
rxSSGrid rxSSMeshGPU::GetMeshGrid(int idx)
{
	rxSSGridG g0;
	CuCopyArrayFromDevice(&g0, m_dSSMGrid+idx, 0, 1*sizeof(rxSSGridG));

	rxSSGrid g;
	for(int j = 0; j < 4; ++j){
		g.node_vrts[j] = g0.node_vrts[j];
		g.edge_vrts[j] = g0.edge_vrts[j];
		g.back_vrts[j] = g0.back_vrts[j];
		g.node_depth[j] = g0.node_depth[j];
	}
	g.back_vrts[4] = g0.back_vrts[4];
	g.back_vrts[5] = g0.back_vrts[5];
	g.num_nv = g0.num_nv;
	g.num_ev = g0.num_ev;
	g.num_bv = g0.num_bv;

	g.num_bv = g0.num_bv;

	g.i = g0.i;
	g.j = g0.j;

	g.mesh_num = g0.mesh_num;
	for(int j = 0; j < g0.mesh_num; ++j){
		g.mesh[j] = g0.mesh[j];
	}

	g.back2 = g0.back2;
	for(int j = 0; j < 14; ++j){
		g.v[j] = g0.v[j];
	}

	return g;
}

/*!
 * スクリーン空間での頂点情報を取得
 * @return 頂点情報(rxSSVertex)
 */
rxSSVertex* rxSSMeshGPU::GetSSVertex(void)
{
	return &m_vSSVertex[0];
}


/*!
 * 輪郭エッジのOpenGL描画
 */
void rxSSMeshGPU::DrawSSEdge(void)
{
	// DrawSSEdge
	unsigned int size_e = m_iNgx*(m_iNgy+1)+(m_iNgx+1)*m_iNgy;

	if(m_hSSEdge) delete [] m_hSSEdge;
	m_hSSEdge = new rxSSEdgeG[size_e];
	CuCopyArrayFromDevice(m_hSSEdge, m_dSSEdge.dPos, 0, size_e*sizeof(rxSSEdgeG));

	glBegin(GL_LINES);
	for(int i = 0; i < (int)size_e; ++i){
		rxSSEdgeG e = m_hSSEdge[i];
		Vec3 x0 = Vec3(e.x0.x, e.x0.y, e.x0.z);
		Vec3 x1 = Vec3(e.x1.x, e.x1.y, e.x1.z);

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
void rxSSMeshGPU::DrawSilhouetteEdge(void)
{
	unsigned int size_e = m_iNgx*(m_iNgy+1)+(m_iNgx+1)*m_iNgy;

	if(m_hSSEdge) delete [] m_hSSEdge;
	m_hSSEdge = new rxSSEdgeG[size_e];
	CuCopyArrayFromDevice(m_hSSEdge, m_dSSEdge.dPos, 0, size_e*sizeof(rxSSEdgeG));

	glBegin(GL_LINES);
	for(int i = 0; i < (int)size_e; ++i){
		rxSSEdgeG e = m_hSSEdge[i];

		if(e.silhouette){
			Vec3 x0 = Vec3(e.x0.x, e.x0.y, e.x0.z);
			Vec3 x1 = Vec3(e.x1.x, e.x1.y, e.x1.z);

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
void rxSSMeshGPU::DrawSSVertex(Vec3 node_color, Vec3 edge_color)
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
void rxSSMeshGPU::DrawMeshGrid(int grid, const Vec3 colors[])
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
void rxSSMeshGPU::DrawField(double minpos[2], double maxpos[2])
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
void rxSSMeshGPU::OutputGridInfo(int grid)
{
	if(grid < 0){
		return;
	}

	rxSSGridG g;
	CuCopyArrayFromDevice(&g, m_dSSMGrid+grid, 0, 1*sizeof(rxSSGridG));

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
void rxSSMeshGPU::OutputGridVertexInfo(int grid)
{
	if(grid < 0){
		cout << "no grid is selected." << endl;
		return;
	}

	int gx = grid%m_iNgx;
	int gy = grid/m_iNgx;
	cout << "grid : " << gx << ", " << gy << endl;

	rxSSGridG g;
	CuCopyArrayFromDevice(&g, m_dSSMGrid+grid, 0, 1*sizeof(rxSSGridG));

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
