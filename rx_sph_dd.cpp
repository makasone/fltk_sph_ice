/*!
  @file rx_sph_dd.cpp
	
  @brief SPH法(Double Density)の実装
	- S.Clavet et al., "Particle-based Viscoelastic Fluid Simulation", SCA2005, 2005. 
	- http://www.iro.umontreal.ca/labs/infographie/papers/Clavet-2005-PVFS/index.html
 
  @author Makoto Fujisawa
  @date   2008-10,2012-11
*/
// FILE --rx_sph_dd.cpp--


//-----------------------------------------------------------------------------
// インクルードファイル
//-----------------------------------------------------------------------------
#include "rx_sph.h"

#include "rx_cu_funcs.cuh"
#include "rx_cu_common.cuh"

#include "rx_pcube.h"


//-----------------------------------------------------------------------------
// グローバル変数
//-----------------------------------------------------------------------------
double g_fPresK = 2.0;


//-----------------------------------------------------------------------------
// rxSPHクラスの実装
//-----------------------------------------------------------------------------
/*!
 * コンストラクタ
 * @param[in] use_opengl VBO使用フラグ
 */
rxDDSPH::rxDDSPH(bool use_opengl) : 
	rxParticleSystemBase(use_opengl), 
	m_hNrm(0), 
	m_hFrc(0),
	m_hVelOld(0),
	m_hDens(0), 
	m_hPres(0), 
	m_hSurf(0), 
	m_pBoundary(0)
{
	m_v3Gravity = Vec3(0.0, -9.80665*0.005, 0.0);

	m_fBuoyancy = 0.0f;

	// 衝突判定用
	m_fDamping = 0.0;
	m_fRestitution = 0.0;

	// double density
	m_fK = (RXREAL)g_fPresK;
	m_fKnear = 2.0f*m_fK;
	m_fViscC = 2.0f;
	m_fViscBeta = 1.0f;
	m_fElasKspring = 0.3f;
	m_fPlasAlpha = 0.3f;
	m_fPlasGamma = 0.1f;	// [0, 0.2]ぐらい

	// 近傍探索セル
	m_pNNGrid = new rxNNGrid(DIM);

	m_uNumParticles = 0;


	m_iColorType = RX_RAMP;
}

/*!
 * デストラクタ
 */
rxDDSPH::~rxDDSPH()
{
	Finalize();
}


/*!
 * シミュレーションの初期化
 * @param[in] max_particles 最大パーティクル数
 * @param[in] boundary_ext 境界の大きさ(各辺の長さの1/2)
 * @param[in] dens 初期密度
 * @param[in] mass パーティクルの質量
 * @param[in] kernel_particle 有効半径h以内のパーティクル数
 */
void rxDDSPH::Initialize(const rxSPHEnviroment &env)
{
	// MARK:Initialize
	RXCOUT << "[rxDDSPH::Initialize]" << endl;

	m_fInitDens        = env.dens;
	m_fMass            = env.mass;
	m_iKernelParticles = env.kernel_particles;

	RXREAL volume = env.max_particles*m_fMass/m_fInitDens;

	m_fEffectiveRadius = pow(((m_iKernelParticles*volume)/(env.max_particles*RX_PI)), 1.0/2.0);
	m_fParticleRadius  = 0.5f*m_fEffectiveRadius;

	//m_fInitDens = 10.0;
	m_fInitDensNear = 2.0*m_fInitDens;
	m_fInitMaxDens = m_fInitDens;

	//m_fMaxDens = 10.0;
	m_fSurfDens = 0.82;

	
	RXCOUT << "particle : " << endl;
	RXCOUT << " n_max = " << env.max_particles << endl;
	RXCOUT << " h = " << m_fEffectiveRadius << endl;
	RXCOUT << " r = " << m_fParticleRadius << endl;
	RXCOUT << " dens = " << m_fInitDens << endl;
	RXCOUT << " mass = " << m_fMass << endl;
	RXCOUT << " kernel_particles = " << m_iKernelParticles << endl;
	//RXCOUT << " volume = " << volume << endl;

	RXREAL h = m_fEffectiveRadius;
	RXREAL r = m_fParticleRadius;

	//
	// 境界設定
	//
	m_pBoundary = new rxSolidBox(env.boundary_cen-env.boundary_ext, env.boundary_cen+env.boundary_ext, -1);
	//m_pBoundary = new rxSolidSphere(Vec3(0.0, 0.1, 0.0), 0.25, -1);
	
	// シミュレーション環境の大きさ
	m_v3EnvMin = m_pBoundary->GetMin();
	m_v3EnvMax = m_pBoundary->GetMax();
	RXCOUT << "simlation range : " << m_v3EnvMin << " - " << m_v3EnvMax << endl;

	Vec3 world_size = m_v3EnvMax-m_v3EnvMin;
	Vec3 world_origin = m_v3EnvMin;
	
	double expansion = 0.01;
	world_origin -= 0.5*expansion*world_size;
	world_size *= (1.0+expansion); // シミュレーション環境全体を覆うように設定

	m_v3EnvMin = world_origin;
	m_v3EnvMax = world_origin+world_size;
	

	// カーネル関数の定数
	m_fWpoly6	=  315.0/(64.0*RX_PI*pow(h, (RXREAL)9.0));
	m_fGWpoly6	= -945.0/(32.0*RX_PI*pow(h, (RXREAL)9.0));
	m_fLWpoly6	= -945.0/(32.0*RX_PI*pow(h, (RXREAL)9.0));

	m_uNumParticles = 0;
	m_bCalDens = true;

	Allocate(env.max_particles);
}

/*!
 * メモリの確保
 *  - 最大パーティクル数で確保
 * @param[in] max_particles 最大パーティクル数
 */
void rxDDSPH::Allocate(int max_particles)
{
	// MARK:Allocate
	assert(!m_bInitialized);

	//m_uNumParticles = max_particles;
	m_uMaxParticles = max_particles;
	
	unsigned int size  = m_uMaxParticles*DIM;

	unsigned int size1 = m_uMaxParticles*2;
	unsigned int mem_size  = sizeof(RXREAL)*size;
	unsigned int mem_size1 = sizeof(RXREAL)*size1;

	//
	// メモリ確保
	//
	m_hPos = new RXREAL[size];
	m_hVel = new RXREAL[size];
	m_hNrm = new RXREAL[size];
	m_hFrc = new RXREAL[size];
	m_hVelOld = new RXREAL[size];
	memset(m_hPos, 0, mem_size);
	memset(m_hVel, 0, mem_size);
	memset(m_hNrm, 0, mem_size);
	memset(m_hFrc, 0, mem_size);
	memset(m_hVelOld, 0, mem_size);

	m_hDens = new RXREAL[size1];
	m_hPres = new RXREAL[size1];
	memset(m_hDens, 0, mem_size1);
	memset(m_hPres, 0, mem_size1);

	m_hPredictPos = new RXREAL[size];
	memset(m_hPredictPos, 0, mem_size);
	m_hDist = new RXREAL[size1];
	memset(m_hDist, 0, mem_size1);

	m_hSurf = new uint[m_uMaxParticles];
	memset(m_hSurf, 0, sizeof(RXREAL)*m_uMaxParticles);

	m_hTmp = new RXREAL[m_uMaxParticles];
	memset(m_hTmp, 0, sizeof(RXREAL)*m_uMaxParticles);

	if(m_bUseOpenGL){
		m_posVBO = createVBO(mem_size);	
		m_colorVBO = createVBO(m_uMaxParticles*4*sizeof(RXREAL));

		SetColorVBO(RX_RAMP);
	}

	// 分割セル設定
	m_pNNGrid->Setup(m_v3EnvMin, m_v3EnvMax, m_fEffectiveRadius, m_uMaxParticles);
	m_vNeighs.resize(m_uMaxParticles);

	m_bInitialized = true;
}

/*!
 * 確保したメモリの解放
 */
void rxDDSPH::Finalize(void)
{
	assert(m_bInitialized);

	// メモリ解放
	if(m_hPos) delete [] m_hPos;
	if(m_hVel) delete [] m_hVel;
	if(m_hNrm) delete [] m_hNrm;
	if(m_hFrc) delete [] m_hFrc;
	if(m_hVelOld) delete [] m_hVelOld;

	if(m_hDens) delete [] m_hDens;
	if(m_hPres) delete [] m_hPres;

	if(m_hPredictPos) delete [] m_hPredictPos;
	if(m_hDist) delete [] m_hDist;

	if(m_hSurf) delete [] m_hSurf;

	if(m_hTmp) delete [] m_hTmp;

	if(m_bUseOpenGL){
		glDeleteBuffers(1, (const GLuint*)&m_posVBO);
		glDeleteBuffers(1, (const GLuint*)&m_colorVBO);
	}

	if(m_pNNGrid) delete m_pNNGrid;
	m_vNeighs.clear();

	if(m_pBoundary) delete m_pBoundary;

	int num_solid = (int)m_vSolids.size();
	for(int i = 0; i < num_solid; ++i){
		delete m_vSolids[i];
	}
	m_vSolids.clear();	

	m_uNumParticles = 0;
	m_uMaxParticles = 0;
}




/*!
 * パーティクルデータの取得
 * @return パーティクルデータのデバイスメモリポインタ
 */
RXREAL* rxDDSPH::GetParticle(void)
{
	return m_hPos;
}

/*!
 * パーティクルデータの取得
 * @return パーティクルデータのデバイスメモリポインタ
 */
RXREAL* rxDDSPH::GetParticleDevice(void)
{
	return 0;
}

/*!
 * カラー値用VBOの編集
 * @param[in] type 色のベースとする物性値
 */
void rxDDSPH::SetColorVBO(int type)
{
	// MRK:SetColorVBO

	switch(type){
	case RX_DENSITY:
		SetColorVBOFromArray(m_hDens, 2, false, 2.0*m_fInitMaxDens);
		break;

	case RX_PRESSURE:
		SetColorVBOFromArray(m_hPres, 2);
		break;

	case RX_TEST:
		SetColorVBOFromArray(m_hTmp, 1, false, m_fTmpMax);
		break;
		
	case RX_SURFACE:	// 表面パーティクルだけ色を変えて描画
		if(m_bUseOpenGL){
			// カラーバッファに値を設定
			glBindBufferARB(GL_ARRAY_BUFFER, m_colorVBO);
			RXREAL *data = (RXREAL*)glMapBufferARB(GL_ARRAY_BUFFER, GL_WRITE_ONLY);
			RXREAL *ptr = data;
			for(uint i = 0; i < m_uNumParticles; ++i){
				uint s = m_hSurf[i];

				*ptr++ = (RXREAL)s;
				*ptr++ = 0.0f;
				*ptr++ = 0.0f;
				*ptr++ = 1.0f;
			}
			glUnmapBufferARB(GL_ARRAY_BUFFER);
		}
		break;

	case RX_RAMP:
		if(m_bUseOpenGL){
			// カラーバッファに値を設定
			glBindBufferARB(GL_ARRAY_BUFFER, m_colorVBO);
			RXREAL *data = (RXREAL*)glMapBufferARB(GL_ARRAY_BUFFER, GL_WRITE_ONLY);
			RXREAL *ptr = data;
			for(uint i = 0; i < m_uNumParticles; ++i){
				RXREAL t = i/(RXREAL)m_uNumParticles;
#if 0
				*ptr++ = rand()/(RXREAL)RAND_MAX;
				*ptr++ = rand()/(RXREAL)RAND_MAX;
				*ptr++ = rand()/(RXREAL)RAND_MAX;
#else
				RX_COLOR_RAMP(t, ptr);
				ptr += 3;
#endif
				*ptr++ = 1.0f;
			}
			glUnmapBufferARB(GL_ARRAY_BUFFER);
		}
		break;

	default:
		break;
	}

}


/*!
 * SPHを1ステップ進める
 * @param[in] dt 時間ステップ幅
 * @param[in] step 現在のステップ数
 * @retval ture  計算完了
 * @retval false 最大ステップ数を超えています
 */
bool rxDDSPH::Update(RXREAL dt, int step)
{
	// MARK:Update : 時間更新
	assert(m_bInitialized);

	Vec3 v_half, acc_f;

	static bool init = false;

	// 初期配置から初期密度を計算
	if(m_bCalDens){
		// 近傍粒子探索
		SetParticlesToCell();

		calDoubleDensity(m_hPos, m_fInitDens, m_fInitDensNear, m_fInitMaxDens);
		if(m_fInitMaxDens > RX_FEQ_EPS){
			m_bCalDens = false;
		}
	}


	for(uint j = 0; j < m_solverIterations; ++j){
		// 近傍粒子探索
		SetParticlesToCell();

		// 外力項(重力，粘性など)
		calExternalForceToVel(dt);

		// 予測位置更新(prediction-relaxation scheme)
		for(uint i = 0; i < m_uNumParticles; ++i){
			for(uint j = 0; j < DIM; ++j){
				m_hPredictPos[DIM*i+j] = m_hPos[DIM*i+j]+dt*m_hVel[DIM*i+j];
			}
		}

		// 予測位置での密度計算
		calDoubleDensity(m_hPredictPos);

		// 粒子間バネの追加と削除
		adjustSprings(dt);

		// 粒子間バネによるElasticity,Plasticity
		applySpringDisplacements(dt);

		// 非圧縮性(double density relaxation)
		calDoubleDensityRelaxation(m_hPredictPos, dt);

		// 位置，速度更新
		for(uint i = 0; i < m_uNumParticles; ++i){
			Vec3 x0, x1, v;
			for(int k = 0; k < 3; ++k){
				int idx = DIM*i+k;
				x0[k] = m_hPos[idx];
				x1[k] = m_hPredictPos[idx];
				v[k] = m_hVel[idx];
			}

			// 境界との衝突判定
			calCollisionSolid(x0, x1, v, dt);

			// 新しい速度と位置で更新
			for(int k = 0; k < 3; ++k){
				int idx = DIM*i+k;
				m_hPredictPos[idx] = x1[k];
				m_hVel[idx] = (m_hPredictPos[idx]-m_hPos[idx])/dt;
				m_hPos[idx] = m_hPredictPos[idx];
			}
		}

		RXTIMER("sph");

		init = false;
	}

	if(m_bCalNormal){
		CalNormal();
	}


	SetArrayVBO(RX_POSITION, m_hPos, 0, m_uNumParticles);

	SetColorVBO(m_iColorType);

	RXTIMER("color(vbo)");

	return true;
}





//-----------------------------------------------------------------------------
// Double Density
//-----------------------------------------------------------------------------
/*!
 * 密度を計算(初期密度計算用)
 * @param[in] ppos パーティクル座標
 * @param[out] avg_dens 平均密度
 * @param[out] avg_dens_near 平均密度(near)
 * @param[out] max_dens 最大密度
 */
void rxDDSPH::calDoubleDensity(RXREAL *ppos, double &avg_dens, double &avg_dens_near, double &max_dens)
{
	avg_dens = 0.0;
	avg_dens_near = 0.0;
	max_dens = 0.0;
	double min_dens = RX_FEQ_INF;

	RXREAL h = m_fEffectiveRadius;
	RXREAL r0 = m_fInitDens;

	RXREAL k = m_fK;
	RXREAL knear = m_fKnear;

	for(uint i = 0; i < m_uNumParticles; ++i){
		Vec3 pos0 = Vec3(ppos[DIM*i+0], ppos[DIM*i+1], ppos[DIM*i+2]);

		m_hDens[2*i+0] = 0.0;	// Density
		m_hDens[2*i+1] = 0.0;	// Near Density

		// 近傍粒子
		for(vector<rxNeigh>::iterator itr = m_vNeighs[i].begin() ; itr != m_vNeighs[i].end(); ++itr){
			int j = itr->Idx;
			if(j < 0 || i == j) continue;

			Vec3 pos1 = Vec3(ppos[DIM*j+0], ppos[DIM*j+1], ppos[DIM*j+2]);
			Vec3 rij = pos1-pos0;

			RXREAL r = norm(rij);//sqrt(itr->Dist2);
			RXREAL q = r/h;

			if(q <= 1.0 && q >= 0.0){
				// density
				m_hDens[2*i]   += m_fMass*(1.0-q)*(1.0-q);

				// near density
				m_hDens[2*i+1] += m_fMass*(1.0-q)*(1.0-q)*(1.0-q);
			}
		}

		// 平均密度の計算
		avg_dens += m_hDens[2*i];
		avg_dens_near += m_hDens[2*i+1];

		if(m_hDens[2*i] > max_dens) max_dens = m_hDens[2*i];
		if(m_hDens[2*i] < min_dens) min_dens = m_hDens[2*i];
	}

	if(m_uNumParticles){
		avg_dens /= (double)m_uNumParticles;
		avg_dens_near /= (double)m_uNumParticles;
	}

	RXCOUT << "minimum density : " << min_dens << endl;
	RXCOUT << "maximum density : " << max_dens << endl;
}

/*!
 * 密度と圧力を計算
 * @param[in] ppos パーティクル座標
 */
void rxDDSPH::calDoubleDensity(RXREAL *ppos)
{
	// MRK:calDensity
	RXREAL h = m_fEffectiveRadius;
	RXREAL r0 = m_fInitDens;

	RXREAL k = m_fK;
	RXREAL knear = m_fKnear;

	for(uint i = 0; i < m_uNumParticles; ++i){
		Vec3 pos0 = Vec3(ppos[DIM*i+0], ppos[DIM*i+1], ppos[DIM*i+2]);

		m_hDens[2*i+0] = 0.0;	// Density
		m_hDens[2*i+1] = 0.0;	// Near Density

		// 近傍粒子
		for(vector<rxNeigh>::iterator itr = m_vNeighs[i].begin() ; itr != m_vNeighs[i].end(); ++itr){
			int j = itr->Idx;
			if(j < 0 || i == j) continue;

			Vec3 pos1 = Vec3(ppos[DIM*j+0], ppos[DIM*j+1], ppos[DIM*j+2]);
			Vec3 rij = pos1-pos0;

			RXREAL r = norm(rij);//sqrt(itr->Dist2);
			RXREAL q = r/h;

			if(q <= 1.0 && q >= 0.0){
				// density
				m_hDens[2*i]   += m_fMass*(1.0-q)*(1.0-q);

				// near density
				m_hDens[2*i+1] += m_fMass*(1.0-q)*(1.0-q)*(1.0-q);
			}
		}

		m_hPres[2*i]   = k*(m_hDens[2*i]-r0);	// pressure
		m_hPres[2*i+1] = knear*m_hDens[2*i+1];	// near-pressure
	}
}

/*!
 * double density relaxationによる非圧縮性条件の適用
 * @param[in] ppos パーティクル座標
 * @param[in] dt タイムステップ幅
 */
void rxDDSPH::calDoubleDensityRelaxation(RXREAL *ppos, RXREAL dt)
{
	RXREAL h = m_fEffectiveRadius;

	Vec3 avg_dx = 0.0;
	for(uint i = 0; i < m_uNumParticles; ++i){
		Vec3 pos0 = Vec3(ppos[DIM*i+0], ppos[DIM*i+1], ppos[DIM*i+2]);

		RXREAL prsi  = m_hPres[2*i];
		RXREAL prsni = m_hPres[2*i+1];

		Vec3 dx(0.0);
		for(vector<rxNeigh>::iterator itr = m_vNeighs[i].begin() ; itr != m_vNeighs[i].end(); ++itr){
			int j = itr->Idx;
			if(j < 0 || i == j) continue;

			Vec3 pos1 = Vec3(ppos[DIM*j+0], ppos[DIM*j+1], ppos[DIM*j+2]);
			Vec3 rij = pos1-pos0;

			RXREAL r = norm(rij);//sqrt(itr->Dist2);
			RXREAL q = r/h;

			if(q <= 1.0 && q > 0.0){
				rij = Unit(rij);

				RXREAL prsj  = m_hPres[2*j];
				RXREAL prsnj = m_hPres[2*j+1];

				RXREAL prs, prsn;
				prs  = 0.5*(prsi+prsj);
				prsn = 0.5*(prsni+prsnj);

				Vec3 D = m_fMass*(0.5*prs*(1.0-q)+0.5*prsn*(1.0-q)*(1.0-q))*rij/m_hDens[2*j];

				dx -= 0.5*dt*dt*D;
				//da += 0.5*D*m_params.Mass;
			}
		}

		m_hPredictPos[DIM*i]   += dx[0];
		m_hPredictPos[DIM*i+1] += dx[1];
	}

}

/*!
 * 外力項，粘性項の計算
 * @param[in] dt タイムステップ幅
 */
void rxDDSPH::calExternalForceToVel(RXREAL dt)
{
	RXREAL h  = m_fEffectiveRadius;
	RXREAL r0 = m_fInitDens;

	RXREAL c  = m_fViscC;
	RXREAL B  = m_fViscBeta;

	// 速度場の保存
	for(uint i = 0; i < m_uNumParticles; ++i){
		for(uint j = 0; j < DIM; ++j){
			int idx = DIM*i+j;
			m_hVelOld[idx] = m_hVel[idx];
		}
	}

	for(uint i = 0; i < m_uNumParticles; ++i){
		Vec3 pos0, vel0;
		pos0 = Vec3(m_hPos[DIM*i+0], m_hPos[DIM*i+1], m_hPos[DIM*i+2]);
		vel0 = Vec3(m_hVelOld[DIM*i+0], m_hVelOld[DIM*i+1], m_hVelOld[DIM*i+2]);

		// 粘性項の計算
		Vec3 rij, vji;
		Vec3 pfrc(0.0), vfrc(0.0);
		for(vector<rxNeigh>::iterator itr = m_vNeighs[i].begin() ; itr != m_vNeighs[i].end(); ++itr){
			int j = itr->Idx;
			if(j < 0) continue;

			if((int)i < j){
				RXREAL r = sqrt(itr->Dist2);
				RXREAL q = r/h;

				if(q < 1.0 && q > RX_FEQ_EPS){
					Vec3 pos1, vel1;
					pos1 = Vec3(m_hPos[DIM*j+0], m_hPos[DIM*j+1], m_hPos[DIM*j+2]);
					vel1 = Vec3(m_hVelOld[DIM*j+0], m_hVelOld[DIM*j+1], m_hVelOld[DIM*j+2]);

					rij = Unit(pos1-pos0);
					RXREAL u = (RXREAL)dot(vel0-vel1, rij);

					if(u > 0.0){
						// 粘性
						Vec3 I = dt*(1.0-q)*(c*u+B*u*u)*rij;
						m_hVel[DIM*i+0] -= 0.5*I[0];
						m_hVel[DIM*i+1] -= 0.5*I[1];
						m_hVel[DIM*i+2] -= 0.5*I[2];
						m_hVel[DIM*j+0] += 0.5*I[0];
						m_hVel[DIM*j+1] += 0.5*I[1];
						m_hVel[DIM*j+2] += 0.5*I[2];
					}
				}
			}
		}

		// 外力項の計算
		Vec3 g;
		g[0] = m_v3Gravity[0];
		g[1] = m_v3Gravity[1];
		g[2] = m_v3Gravity[2];
		Vec3 v = dt*g;

		m_hVel[DIM*i+0] += v[0];
		m_hVel[DIM*i+1] += v[1];
		m_hVel[DIM*i+2] += v[2];
	}
}


/*!
 * 粒子間バネの追加と削除
 * @param[in] dt タイムステップ幅
 */
void rxDDSPH::adjustSprings(RXREAL dt)
{
	RXREAL h = m_fEffectiveRadius;
	RXREAL r0 = m_fInitDens;

	RXREAL gamma = m_fPlasGamma; // fraction of L : [0, 0.2]
	RXREAL alpha = m_fPlasAlpha; // plasticity constant

	map<int, rxSPHSpring>::iterator itr;
	for(itr = m_mapSpring.begin(); itr != m_mapSpring.end(); ++itr){
		itr->second.enable = false;
	}

	// 粒子間バネの追加と調整
	for(uint i = 0; i < m_uNumParticles; ++i){
		for(vector<rxNeigh>::iterator itr = m_vNeighs[i].begin() ; itr != m_vNeighs[i].end(); ++itr){
			int j = itr->Idx;
			if(j < 0) continue;

			if((int)i < j){
				RXREAL r = sqrt(itr->Dist2);
				RXREAL q = r/h;

				int idx = j*m_uNumParticles+i;
				if(!m_mapSpring[idx].enable){
					// 粒子間バネを新たに追加
					m_mapSpring[idx].pi = i;
					m_mapSpring[idx].pj = j;
					m_mapSpring[idx].L = r;
				}

				m_mapSpring[idx].enable = true;
				m_mapSpring[idx].r = r;

				RXREAL d = gamma*m_mapSpring[idx].L;

				if(r > m_mapSpring[idx].L+d){
					m_mapSpring[idx].L += dt*alpha*(r-m_mapSpring[idx].L-d);
				}
				else if(r < m_mapSpring[idx].L-d){
					m_mapSpring[idx].L -= dt*alpha*(m_mapSpring[idx].L-d-r);
				}
			}
		}
	}

	// 粒子間バネの削除
	map<int, rxSPHSpring>::iterator iter = m_mapSpring.begin();
	while(iter != m_mapSpring.end()){
		if(!iter->second.enable){
			m_mapSpring.erase(iter++);
		}
		else{
			++iter;
		}
	}
}


/*!
 * Elasticity,Plasticity
 * @param[in] dt タイムステップ幅
 */
void rxDDSPH::applySpringDisplacements(RXREAL dt)
{
	RXREAL h = m_fEffectiveRadius;
	RXREAL r0 = m_fInitDens;

	RXREAL ks = m_fElasKspring; // バネ定数

	map<int, rxSPHSpring>::iterator itr;
	for(itr = m_mapSpring.begin(); itr != m_mapSpring.end(); ++itr){
		if(itr->second.enable){
			int i = itr->second.pi;
			int j = itr->second.pj;

			Vec3 pos0, pos1;
			pos0 = Vec3(m_hPos[DIM*i+0], m_hPos[DIM*i+1], m_hPos[DIM*i+2]);
			pos1 = Vec3(m_hPos[DIM*j+0], m_hPos[DIM*j+1], m_hPos[DIM*j+2]);

			Vec3 rij = Unit(pos1-pos0);
			Vec3 D = dt*dt*ks*(1.0-itr->second.L/h)*(itr->second.L-itr->second.r)*rij;

			pos0 -= 0.5*D;
			pos1 += 0.5*D;

			for(int k = 0; k < DIM; ++k){
				m_hPos[DIM*i+k] = pos0[k];
				m_hPos[DIM*j+k] = pos1[k];
			}
		}
	}
}


/*!
 * 固体オブジェクトとの衝突判定，衝突応答
 * @param[in] pos0 前ステップの位置
 * @param[inout] pos1 新しい位置
 * @param[inout] vel 速度
 * @param[in] dt タイムステップ幅
 * @return 衝突オブジェクトの数
 */
int rxDDSPH::calCollisionSolid(Vec3 &pos0, Vec3 &pos1, Vec3 &vel, RXREAL dt)
{
	int c = 0;
	rxCollisionInfo coli;

	// 固体オブジェクトとの衝突処理
	for(vector<rxSolid*>::iterator i = m_vSolids.begin(); i != m_vSolids.end(); ++i){
		if((*i)->GetDistanceR(pos1, m_fParticleRadius, coli)){
			RXREAL res = m_fDamping;
			res = (res > 0) ? (res*fabs(coli.Penetration())/(dt*norm(vel))) : 0.0f;
			vel -= (1+res)*dot(vel, coli.Normal())*coli.Normal();
			pos1 = coli.Contact();
		}
	}

	// シミュレーション空間境界との衝突処理
	if(m_pBoundary->GetDistanceR(pos1, m_fParticleRadius, coli)){
		RXREAL res = m_fDamping;
		res = (res > 0) ? (res*fabs(coli.Penetration())/(dt*norm(vel))) : 0.0f;
		vel -= (1+res)*dot(vel, coli.Normal())*coli.Normal();
		pos1 = coli.Contact()+coli.Normal()*RXFunc::Rand(0.0, 0.001);
	}

	return c;
}



//-----------------------------------------------------------------------------
// 表面パーティクルの検出
//-----------------------------------------------------------------------------
/*!
 * 表面パーティクル検出
 *  - B. Solenthaler, Y. Zhang and R. Pajarola, "Efficient Refinement of Dynamic Point Data",
 *    Proceedings Eurographics/IEEE VGTC Symposium on Point-Based Graphics, 2007.
 *  - 3.1 Surface Particle Detection の 式(2) 
 *  - d_i,cm が閾値以上で表面パーティクルと判定と書かれているが，
 *    d_i,cm は近傍パーティクルの重心へのベクトルで実際にはそのベクトルの長さ |d_i,cm| を使って判定
 */
void rxDDSPH::DetectSurfaceParticles(void)
{
	RXREAL h = m_fEffectiveRadius;
	RXREAL r = m_fParticleRadius;

	//m_hPos = GetArrayVBO(RX_POSITION);

	// 近傍粒子探索
	SetParticlesToCell();

	for(uint i = 0; i < m_uNumParticles; ++i){
		double d = CalDistToNormalizedMassCenter(i);
		int nn_num = (int)m_vNeighs[i].size();

		// デバッグ用
		//m_hTmp[i] = (RXREAL)d;
		//m_fTmpMax = r;

		//m_hSurf[i] = 0;
		if(nn_num <= 3){	// 近傍パーティクル数が3以下ならば表面
			m_hSurf[i] = 1;
		}
		else{				// 3より大きい場合は近傍重心までの距離で判断
			if(m_hSurf[i]){
				// 前ステップで表面だったら小さい閾値で判断
				if(d < g_fSurfThr[0]*r) m_hSurf[i] = 0;
			}
			else{
				if(d > g_fSurfThr[1]*r || d < 0.0) m_hSurf[i] = 1;
			}
		}

		// 密度を使って判断する場合
		//if(m_hDens[2*i] < 0.7*m_fInitMaxDens){
		//	m_hSurf[i] = 1;
		//}
	}
}

/*!
 * 近傍パーティクルの正規化重心までの距離を計算
 * @param[in] i パーティクルインデックス
 */
double rxDDSPH::CalDistToNormalizedMassCenter(const int i)
{
	Vec3 sum_pos(0.0);
	double sum_mass = 0.0;

	Vec3 pos0 = Vec3(m_hPos[DIM*i+0], m_hPos[DIM*i+1], m_hPos[DIM*i+2]);

	for(vector<rxNeigh>::iterator itr = m_vNeighs[i].begin() ; itr != m_vNeighs[i].end(); ++itr){
		int j = itr->Idx;
		if(j < 0 && i == j) continue;

		//RXREAL r = sqrt(itr->Dist2);
		Vec3 pos1 = Vec3(m_hPos[DIM*j+0], m_hPos[DIM*j+1], m_hPos[DIM*j+2]);

		sum_pos  += (pos0-pos1)*m_fMass;
		sum_mass += m_fMass;
	}

	double dis = DBL_MIN;
	if(sum_mass != 0.0){
		sum_pos /= sum_mass;
		dis     = norm(sum_pos);
	}

	return dis;
}

/*!
 * 表面パーティクル情報の取得
 *  パーティクル数と同じ大きさのuint配列で表面パーティクルならば1, そうでなければ0が格納されている
 */
uint* rxDDSPH::GetArraySurf(void)
{
	return m_hSurf;
}

/*!
 * 指定した座標の近傍の表面パーティクル情報を取得する
 * @param[in] pos 探索中心座標
 * @param[in] h 探索半径(0以下の値を設定したら有効半径を用いる)
 * @param[out] sp 表面パーティクル
 * @return 見つかったパーティクル数
 */
int rxDDSPH::GetSurfaceParticles(const Vec3 pos, RXREAL h, vector<rxSurfaceParticle> &sps)
{
	int num = 0;
	sps.clear();

	DetectSurfaceParticles();

	if(pos[0] < m_v3EnvMin[0]) return 0;
	if(pos[0] > m_v3EnvMax[0]) return 0;
	if(pos[1] < m_v3EnvMin[1]) return 0;
	if(pos[1] > m_v3EnvMax[1]) return 0;
	if(pos[2] < m_v3EnvMin[2]) return 0;
	if(pos[2] > m_v3EnvMax[2]) return 0;

	if(h <= 0.0) h = m_fEffectiveRadius;

	vector<rxNeigh> ne;
	m_pNNGrid->GetNN(pos, m_hPos, m_uNumParticles, ne, h);

	// 近傍粒子
	for(vector<rxNeigh>::iterator itr = ne.begin(); itr != ne.end(); ++itr){
		int j = itr->Idx;
		if(j < 0) continue;

		int surf = m_hSurf[j];

		if(surf){
			rxSurfaceParticle sp;
			sp.pos = Vec3(m_hPos[DIM*j+0], m_hPos[DIM*j+1], m_hPos[DIM*j+2]);
			sp.vel = Vec3(m_hVel[DIM*j+0], m_hVel[DIM*j+1], m_hVel[DIM*j+2]);
			sp.nrm = Vec3(m_hNrm[DIM*j+0], m_hNrm[DIM*j+1], m_hNrm[DIM*j+2]);
			sp.idx = j;
			sp.d = sqrt(itr->Dist2);
			sps.push_back(sp);
			num++;
		}
	}

	return num;
}



//-----------------------------------------------------------------------------
// 法線計算
//-----------------------------------------------------------------------------
/*!
 * 法線を計算
 *  - 越川純弥，"粒子法による流体シミュレーションの高品質レンダリング"，静岡大学卒業論文，2007.
 *  - p.24, 3.4.3 法線の修正 の 式(3.8) 
 *  - 式(3.8) では Σw(v-r)(v-r) となっているが Σw(v-r) の間違い
 */
void rxDDSPH::CalNormal(void)
{
	RXREAL h = m_fEffectiveRadius;

	// 近傍粒子探索
	SetParticlesToCell();

	//CalNormalFromDensity();

	for(uint i = 0; i < m_uNumParticles; ++i){
		Vec3 pos0 = Vec3(m_hPos[DIM*i+0], m_hPos[DIM*i+1], m_hPos[DIM*i+2]);
		Vec3 nrm(0.0);

		// 近傍粒子
		for(vector<rxNeigh>::iterator itr = m_vNeighs[i].begin() ; itr != m_vNeighs[i].end(); ++itr){
			int j = itr->Idx;
			if(j < 0 || i == j) continue;

			Vec3 pos1 = Vec3(m_hPos[DIM*j+0], m_hPos[DIM*j+1], m_hPos[DIM*j+2]);
			Vec3 rij = pos0-pos1;

			RXREAL rr = sqrt(itr->Dist2)/h;

			RXREAL w = (rr < 1.0 ? 1.0/rr-1.0 : 0.0);

			nrm += w*rij;
		}

		//normalize(nrm);
		//nrm *= -1;

		m_hNrm[DIM*i+0] = nrm[0];
		m_hNrm[DIM*i+1] = nrm[1];
	}
}

/*!
 * 密度分布から法線を計算
 */
void rxDDSPH::CalNormalFromDensity(void)
{
	RXREAL h = m_fEffectiveRadius;

	for(uint i = 0; i < m_uNumParticles; ++i){
		Vec3 pos0;
		pos0[0] = m_hPos[DIM*i+0];
		pos0[1] = m_hPos[DIM*i+1];
		pos0[2] = m_hPos[DIM*i+2];

		Vec3 nrm(0.0);

		// 近傍粒子
		for(vector<rxNeigh>::iterator itr = m_vNeighs[i].begin() ; itr != m_vNeighs[i].end(); ++itr){
			int j = itr->Idx;
			if(j < 0 || i == j) continue;

			Vec3 pos1;
			pos1[0] = m_hPos[DIM*j+0];
			pos1[1] = m_hPos[DIM*j+1];
			pos1[2] = m_hPos[DIM*j+2];

			Vec3 rij = pos0-pos1;

			RXREAL r = sqrt(itr->Dist2);
			RXREAL q = h*h-r*r;

			nrm += (m_fMass/m_hDens[j])*m_fGWpoly6*q*q*rij;
		}

		//normalize(nrm);
		//nrm *= -1;

		m_hNrm[DIM*i+0] = nrm[0];
		m_hNrm[DIM*i+1] = nrm[1];
	}
}



//-----------------------------------------------------------------------------
// 近傍探索
//-----------------------------------------------------------------------------
/*!
 * 近傍粒子探索
 * @param[in] idx 探索中心パーティクルインデックス
 * @param[in] prts パーティクル位置
 * @param[out] neighs 探索結果格納する近傍情報コンテナ
 * @param[in] h 有効半径
 */
void rxDDSPH::GetNearestNeighbors(int idx, RXREAL *prts, vector<rxNeigh> &neighs, RXREAL h)
{
	if(idx < 0 || idx >= (int)m_uNumParticles) return;

	Vec3 pos;
	pos[0] = prts[DIM*idx+0];
	pos[1] = prts[DIM*idx+1];
	pos[2] = prts[DIM*idx+2];

	if(h < 0.0) h = m_fEffectiveRadius;

	m_pNNGrid->GetNN(pos, prts, m_uNumParticles, neighs, h);
	//m_pNNGrid->GetNN_Direct(pos, prts, m_uNumParticles, neighs, h);	// グリッドを使わない総当たり
}

/*!
 * 近傍粒子探索
 * @param[in] idx 探索中心パーティクルインデックス
 * @param[out] neighs 探索結果格納する近傍情報コンテナ
 * @param[in] h 有効半径
 */
void rxDDSPH::GetNearestNeighbors(Vec3 pos, vector<rxNeigh> &neighs, RXREAL h)
{
	if(h < 0.0) h = m_fEffectiveRadius;

	m_pNNGrid->GetNN(pos, m_hPos, m_uNumParticles, neighs, h);
	//m_pNNGrid->GetNN_Direct(pos, prts, m_uNumParticles, neighs, h);	// グリッドを使わない総当たり
}


/*!
 * 全パーティクルを分割セルに格納
 */
void rxDDSPH::SetParticlesToCell(RXREAL *prts, int n, RXREAL h)
{
	// 分割セルに粒子を登録
	m_pNNGrid->SetObjectToCell(prts, n);

	// 近傍粒子探索
	if(h < 0.0) h = m_fEffectiveRadius;
	for(int i = 0; i < (int)m_uNumParticles; i++){
		m_vNeighs[i].clear();
		GetNearestNeighbors(i, prts, m_vNeighs[i], h);
	}
}
void rxDDSPH::SetParticlesToCell(void)
{
	SetParticlesToCell(m_hPos, m_uNumParticles, m_fEffectiveRadius);
}

/*!
 * 分割セルに格納されたポリゴン情報を取得
 * @param[in] gi,gj,gk 対象分割セル
 * @param[out] polys ポリゴン
 * @return 格納ポリゴン数
 */
int rxDDSPH::GetPolygonsInCell(int gi, int gj, int gk, vector<int> &polys)
{
	return m_pNNGrid->GetPolygonsInCell(gi, gj, gk, polys);
}
/*!
 * 分割セル内のポリゴンの有無を調べる
 * @param[in] gi,gj,gk 対象分割セル
 * @return ポリゴンが格納されていればtrue
 */
bool rxDDSPH::IsPolygonsInCell(int gi, int gj, int gk)
{
	return m_pNNGrid->IsPolygonsInCell(gi, gj, gk);
}

/*!
 * ポリゴンを分割セルに格納
 */
void rxDDSPH::SetPolygonsToCell(void)
{
	m_pNNGrid->SetPolygonsToCell(m_hVrts, m_iNumVrts, m_hTris, m_iNumTris);
}


//-----------------------------------------------------------------------------
// OpenGL描画
//-----------------------------------------------------------------------------
/*!
 * 探索用セルの描画
 * @param[in] i,j,k グリッド上のインデックス
 */
void rxDDSPH::DrawCell(int i, int j, int k)
{
	if(m_pNNGrid) m_pNNGrid->DrawCell(i, j, k);
}

/*!
 * 探索用グリッドの描画
 * @param[in] col パーティクルが含まれるセルの色
 * @param[in] col2 ポリゴンが含まれるセルの色
 * @param[in] sel ランダムに選択されたセルのみ描画(1で新しいセルを選択，2ですでに選択されているセルを描画，0ですべてのセルを描画)
 */
void rxDDSPH::DrawCells(Vec3 col, Vec3 col2, int sel)
{
	if(m_pNNGrid) m_pNNGrid->DrawCells(col, col2, sel, m_hPos);
}


/*!
 * 固体障害物の描画
 */
void rxDDSPH::DrawObstacles(void)
{
	for(vector<rxSolid*>::iterator i = m_vSolids.begin(); i != m_vSolids.end(); ++i){
		(*i)->Draw();
	}
}



//-----------------------------------------------------------------------------
// 固体オブジェクト
//-----------------------------------------------------------------------------
/*!
 * 三角形ポリゴンによる障害物
 * @param[in] vrts 頂点
 * @param[in] tris メッシュ
 */
void rxDDSPH::SetPolygonObstacle(const vector<Vec3> &vrts, const vector<Vec3> &nrms, const vector< vector<int> > &tris, Vec3 vel)
{
	int vn = (int)vrts.size();
	int n = (int)tris.size();

	if(m_hVrts) delete [] m_hVrts;
	if(m_hTris) delete [] m_hTris;
	m_hVrts = new RXREAL[vn*3];
	m_hTris = new int[n*3];

	for(int i = 0; i < vn; ++i){
		for(int j = 0; j < 3; ++j){
			m_hVrts[3*i+j] = vrts[i][j];
		}
	}

	for(int i = 0; i < n; ++i){
		for(int j = 0; j < 3; ++j){
			m_hTris[3*i+j] = tris[i][j];
		}
	}

	m_iNumVrts = vn;
	m_iNumTris = n;
	RXCOUT << "the number of triangles : " << m_iNumTris << endl;

	m_pNNGrid->SetPolygonsToCell(m_hVrts, m_iNumVrts, m_hTris, m_iNumTris);

	//setPolysToCell(m_hVrts, m_iNumVrts, m_hTris, m_iNumTris);
}

/*!
 * ボックス型障害物
 * @param[in] cen ボックス中心座標
 * @param[in] ext ボックスの大きさ(辺の長さの1/2)
 * @param[in] ang ボックスの角度(オイラー角)
 * @param[in] flg 有効/無効フラグ
 */
void rxDDSPH::SetBoxObstacle(Vec3 cen, Vec3 ext, Vec3 ang, Vec3 vel, int flg)
{
	rxSolidBox *box = new rxSolidBox(cen-ext, cen+ext, 1);
	double m[16];
	EulerToMatrix(m, ang[0], ang[1], ang[2]);
	box->SetMatrix(m);

	m_vSolids.push_back(box);

}

/*!
 * 球型障害物
 * @param[in] cen 球体中心座標
 * @param[in] rad 球体の半径
 * @param[in] flg 有効/無効フラグ
 */
void rxDDSPH::SetSphereObstacle(Vec3 cen, double rad, Vec3 vel, int flg)
{
	rxSolidSphere *sphere = new rxSolidSphere(cen, rad, 1);
	m_vSolids.push_back(sphere);
}

/*!
 * 球型障害物を動かす
 * @param[in] b 物体番号
 * @param[in] disp 移動量
 */
void rxDDSPH::MoveSphereObstacle(int b, Vec3 disp)
{
	if(m_vSolids.empty()) return;
	m_vSolids[b]->SetPosition(m_vSolids[b]->GetPosition()+disp);
}

/*!
 * 球型障害物の位置を取得
 * @param[in] b 物体番号
 */
Vec3 rxDDSPH::GetSphereObstaclePos(int b)
{
	if(m_vSolids.empty()) return Vec3(0.0);
	return m_vSolids[b]->GetPosition();
}

/*!
 * VBOからホストメモリへデータを転送，取得
 * @param[in] type データの種類
 * @return ホストメモリ上のデータ
 */
RXREAL* rxDDSPH::GetArrayVBO(rxParticleArray type, bool d2h, int num)
{
	assert(m_bInitialized);
 
	if(num == -1) num = m_uNumParticles;
	RXREAL* hdata = 0;

	unsigned int vbo = 0;

	switch(type){
	default:
	case RX_POSITION:
		hdata = m_hPos;
		vbo = m_posVBO;
		break;

	case RX_VELOCITY:
		hdata = m_hVel;
		break;

	case RX_DENSITY:
		hdata = m_hDens;
		break;

	case RX_PRESSURE:
		hdata = m_hPres;
		break;

	case RX_NORMAL:
		hdata = m_hNrm;
		break;

	case RX_FORCE:
		hdata = m_hFrc;
		break;
	}

	return hdata;
}

/*!
 * ホストメモリからVBOメモリへデータを転送
 * @param[in] type データの種類
 * @param[in] data ホストメモリ上のデータ
 * @param[in] start データの開始インデックス
 * @param[in] count 追加数
 */
void rxDDSPH::SetArrayVBO(rxParticleArray type, const RXREAL* data, int start, int count)
{
	assert(m_bInitialized);
 
	switch(type){
	default:
	case RX_POSITION:
		{
			if(m_bUseOpenGL){
				glBindBuffer(GL_ARRAY_BUFFER, m_posVBO);
				glBufferSubData(GL_ARRAY_BUFFER, start*2*sizeof(RXREAL), count*2*sizeof(RXREAL), data);
				glBindBuffer(GL_ARRAY_BUFFER, 0);
			}
		}
		break;

	case RX_VELOCITY:
		//CuCopyArrayToDevice(m_dVel, data, start*DIM*sizeof(RXREAL), count*DIM*sizeof(RXREAL));
		break;

	case RX_NORMAL:
		//CuCopyArrayToDevice(m_dNrm, data, start*DIM*sizeof(RXREAL), count*DIM*sizeof(RXREAL));
		break;
	}	   
}


//-----------------------------------------------------------------------------
// 陰関数値
//-----------------------------------------------------------------------------
double rxDDSPH::GetImplicit(double x, double y, double z)
{
	return CalColorField(x, y, z);
}

/*!
 * パーティクルからグリッドの陰関数値を計算
 * @param[in] n グリッド数
 * @param[in] minp グリッドの最小座標
 * @param[in] d グリッド幅
 * @param[out] hF 陰関数値(nx×ny×nzの配列)
 */
void rxDDSPH::CalImplicitField(int n[3], Vec3 minp, Vec3 d, RXREAL *hF)
{
	int slice0 = n[0];
	int slice1 = n[0]*n[1];

	for(int k = 0; k < n[2]; ++k){
		for(int j = 0; j < n[1]; ++j){
			for(int i = 0; i < n[0]; ++i){
				int idx = k*slice1+j*slice0+i;
				Vec3 pos = minp+Vec3(i, j, k)*d;
				hF[idx] = GetImplicit(pos[0], pos[1], pos[2]);
			}
		}
	}
}

/*!
 * パーティクルからグリッドの陰関数値を計算
 * @param[in] n グリッド数
 * @param[in] minp グリッドの最小座標
 * @param[in] d グリッド幅
 * @param[out] hF 陰関数値(nx×ny×nzの配列)
 */
void rxDDSPH::CalImplicitFieldDevice(int n[3], Vec3 minp, Vec3 d, RXREAL *dF)
{
	RXREAL *hF = new RXREAL[n[0]*n[1]*n[2]];

	CalImplicitField(n, minp, d, hF);
	//CuCopyArrayToDevice(dF, hF, 0, n[0]*n[1]*n[2]*sizeof(RXREAL));

	delete [] hF;
}


/*!
 * カラーフィールド値計算
 *  - Poly6カーネル
 * @param[in] x,y,z 計算位置
 * @return カラーフィールド値
 */
double rxDDSPH::CalColorField(double x, double y, double z)
{
	// MRK:CalColorField
	RXREAL c = 0.0;
	Vec3 pos(x, y, z);

	if(pos[0] < m_v3EnvMin[0]) return c;
	if(pos[0] > m_v3EnvMax[0]) return c;
	if(pos[1] < m_v3EnvMin[1]) return c;
	if(pos[1] > m_v3EnvMax[1]) return c;
	if(pos[2] < m_v3EnvMin[2]) return c;
	if(pos[2] > m_v3EnvMax[2]) return c;

	RXREAL h = m_fEffectiveRadius;

	vector<rxNeigh> ne;
	m_pNNGrid->GetNN(pos, m_hPos, m_uNumParticles, ne, h);

	// 近傍粒子
	for(vector<rxNeigh>::iterator itr = ne.begin(); itr != ne.end(); ++itr){
		int j = itr->Idx;
		if(j < 0) continue;

		Vec3 pos1;
		pos1[0] = m_hPos[DIM*j+0];
		pos1[1] = m_hPos[DIM*j+1];
		pos1[2] = m_hPos[DIM*j+2];

		RXREAL r = sqrt(itr->Dist2);

		RXREAL q = h*h-r*r;
		c += m_fMass*m_fWpoly6*q*q*q;

		//RXREAL q = r/h;
		//if(q <= 1.0 && q >= 0.0){
		//	c += m_fMass*(1.0-q)*(1.0-q);
		//}
	}

	return c;
}


//-----------------------------------------------------------------------------
// MARK:シミュデータの出力
//-----------------------------------------------------------------------------
/*!
 * シミュレーション設定(パーティクル数，範囲，密度，質量など)
 * @param[in] fn 出力ファイル名
 */
void rxDDSPH::OutputSetting(string fn)
{
	ofstream fout;
	fout.open(fn.c_str());
	if(!fout){
		RXCOUT << fn << " couldn't open." << endl;
		return;
	}

	fout << m_uNumParticles << endl;
	fout << m_v3EnvMin[0] << " " << m_v3EnvMin[1] << " " << m_v3EnvMin[2] << endl;
	fout << m_v3EnvMax[0] << " " << m_v3EnvMax[1] << " " << m_v3EnvMax[2] << endl;
	fout << m_fInitDens << endl;
	fout << m_fMass << endl;
	fout << m_iKernelParticles << endl;

	fout.close();
}


