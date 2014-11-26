/*!
  @file Ice_SM.h
	
  @brief ShapeMatching法(GPU実装)
  @ref M. Muller et al., "Meshless deformations based on shape matching", SIGGRAPH2005. 
 
  @author Ryo Nakasone
  @date 2014-7
*/

#include "Ice_SM.h"

const float* Ice_SM::s_pfPrtPos;
const float* Ice_SM::s_pfPrtVel;

float* Ice_SM::s_pfSldPos;
float* Ice_SM::s_pfSldVel;

//デバイスポインタ
float* Ice_SM::sd_PrtPos;	
float* Ice_SM::sd_PrtVel;

float* Ice_SM::d_OrgPos;
float* Ice_SM::d_CurPos;

float* Ice_SM::d_OrgCm;
float* Ice_SM::d_CurCm;

float* Ice_SM::d_Apq;

float* Ice_SM::d_Mass;
float* Ice_SM::d_Vel;

bool* Ice_SM::d_Fix;

int* Ice_SM::d_PIndxes;
int* Ice_SM::d_IndxSet;					//クラスタのデータの開始添字と終了添字を保存

int Ice_SM::s_vertNum;
int Ice_SM::s_vertSum;					//全クラスタに含まれる粒子の総数

int Ice_SM::s_iIterationNum;


Ice_SM::Ice_SM(int obj) : rxShapeMatching(obj)
{
	m_mtrx3Apq = rxMatrix3(0.0);
	m_mtrxBeforeU.makeIdentity();

	m_iIndxNum = 0;

	m_fDefAmount = 0.0f;

	m_pPrePos = new float[MAXPARTICLE*SM_DIM];

	m_fpAlphas = new float[MAXPARTICLE];
	m_fpBetas =	new float[MAXPARTICLE];

	m_ipLayeres = new int[MAXPARTICLE];
}

/*!
 * デストラクタ
 */
Ice_SM::~Ice_SM()
{
}

void Ice_SM::InitFinalParamPointer(int vrtxNum)
{
	s_pfSldPos = new float[vrtxNum*SM_DIM];
	s_pfSldVel = new float[vrtxNum*SM_DIM];

	for(int i = 0; i < vrtxNum; ++i)
	{
		int pIndx = i*4;
		int cIndx = i*SM_DIM;

		for(int j = 0; j < SM_DIM; j++)
		{
			s_pfSldPos[cIndx+j] = s_pfPrtPos[pIndx+j];
			s_pfSldVel[cIndx+j] = s_pfPrtVel[pIndx+j];
		}
	}
}

//GPU計算のための初期化
void Ice_SM::InitGPU(const vector<Ice_SM*>& ice_sm, float* d_pos, float* d_vel, int prtNum, int MAXCLUSTER)
{
	//デバイスポインタのアドレスを保存
	sd_PrtPos = d_pos;
	sd_PrtVel = d_vel;

	s_vertNum = prtNum;

	//デバイス側のメモリを確保
	//（最大クラスタ数）×（クラスタが保存できる最大粒子数）　でメモリを確保．
	cudaMalloc((void**)&d_OrgPos,	sizeof(float) * MAXCLUSTER * MAXPARTICLE * SM_DIM);
	cudaMalloc((void**)&d_CurPos,	sizeof(float) * MAXCLUSTER * MAXPARTICLE * SM_DIM);
	cudaMalloc((void**)&d_Vel,		sizeof(float) * MAXCLUSTER * MAXPARTICLE * SM_DIM);

	cudaMalloc((void**)&d_OrgCm,	sizeof(float) * MAXCLUSTER * SM_DIM);
	cudaMalloc((void**)&d_CurCm,	sizeof(float) * MAXCLUSTER * SM_DIM);

	cudaMalloc((void**)&d_Apq,		sizeof(float) * MAXCLUSTER * 9);	//TODO: 初期化していないのに注意

	cudaMalloc((void**)&d_Mass,		sizeof(float) * MAXCLUSTER * MAXPARTICLE);
	cudaMalloc((void**)&d_Fix,		sizeof(bool)  * MAXCLUSTER * MAXPARTICLE);
	cudaMalloc((void**)&d_PIndxes,	sizeof(int)   * MAXCLUSTER * MAXPARTICLE);

	cudaMalloc((void**)&d_IndxSet,	sizeof(int)   * MAXCLUSTER * 2);

	//CPUのデータを１次元配列にコピー
	//ホスト側のデータをデバイス側へ転送して初期化
		//一気に大量のメモリは確保できないので、分割してコピーする
	float* oPoses = new float[MAXPARTICLE * SM_DIM];
	float* cPoses = new float[MAXPARTICLE * SM_DIM];
	float* veles  = new float[MAXPARTICLE * SM_DIM];

	float* orgCm = new float[MAXCLUSTER * SM_DIM];
	float* curCm = new float[MAXCLUSTER * SM_DIM];

	float* mass = new float[MAXPARTICLE];
	bool* fix = new bool[MAXPARTICLE];
	int* indx = new int[MAXPARTICLE];

	int* set = new int[2];

	s_vertSum = 0;

	for(int i = 0; i < MAXCLUSTER; i++)
	{
		Ice_SM* sm = ice_sm[i];

		for(int j = 0; j < MAXPARTICLE; j++)
		{
			//現在位置，初期位置，速度
			Vec3 oPos = sm->GetVertexPos(j);
			Vec3 cPos = sm->GetOrgPos(j);
			Vec3 vel  = sm->GetVertexVel(j);

			for(int k = 0; k < SM_DIM; k++)
			{
				oPoses[j*SM_DIM + k] = oPos[k];
				cPoses[j*SM_DIM + k] = cPos[k];
				veles [j*SM_DIM + k] = vel[k];
			}

			//質量，選択フラグ，粒子番号
			mass[j] = sm->GetMass(j);
			fix[j] = false;
			indx[j] = sm->GetParticleIndx(j);
		}

		Vec3 orgCmVec = sm->GetOrgCm();
		Vec3 curCmVec = sm->GetCm();
		
		for(int j = 0; j < SM_DIM; j++)
		{
			orgCm[i*SM_DIM+j] = orgCmVec[j];
			curCm[i*SM_DIM+j] = curCmVec[j];
		}

		//配列のどこからどこまでが埋まっているか
		set[0] = i*MAXPARTICLE;
		set[1] = i*MAXPARTICLE + sm->GetNumVertices()-1;

		//総粒子数
		s_vertSum += sm->GetNumVertices();

		//初期化
		int vecSize = i * MAXPARTICLE * SM_DIM;

		cudaMemcpy(d_OrgPos+vecSize,	oPoses,	sizeof(float) * MAXPARTICLE * SM_DIM, cudaMemcpyHostToDevice);
		cudaMemcpy(d_CurPos+vecSize,	cPoses,	sizeof(float) * MAXPARTICLE * SM_DIM, cudaMemcpyHostToDevice);
		cudaMemcpy(d_Vel   +vecSize,	veles,	sizeof(float) * MAXPARTICLE * SM_DIM, cudaMemcpyHostToDevice);

		int smSize = i * MAXPARTICLE;

		cudaMemcpy(d_Mass   +smSize,	mass,	sizeof(float) * MAXPARTICLE, cudaMemcpyHostToDevice);
		cudaMemcpy(d_Fix    +smSize,	fix,	sizeof(bool)  * MAXPARTICLE, cudaMemcpyHostToDevice);
		cudaMemcpy(d_PIndxes+smSize,	indx,	sizeof(int)   * MAXPARTICLE, cudaMemcpyHostToDevice);

		cudaMemcpy(d_IndxSet+i*2,		set,	sizeof(int) * 2, cudaMemcpyHostToDevice);
	}

	//初期化
	cudaMemcpy(d_OrgCm,	orgCm,	sizeof(float) * MAXCLUSTER * SM_DIM, cudaMemcpyHostToDevice);
	cudaMemcpy(d_CurCm,	curCm,	sizeof(float) * MAXCLUSTER * SM_DIM, cudaMemcpyHostToDevice);

//デバッグ
	////一時配列のリセット
	//for(int i = 0; i < MAXCLUSTER; i++)
	//{
	//	Ice_SM* sm = ice_sm[i];

	//	for(int j = 0; j < MAXPARTICLE; j++)
	//	{
	//		for(int k = 0; k < SM_DIM; k++)
	//		{
	//			oPoses[j*SM_DIM + k] = 0.0f;
	//			cPoses[j*SM_DIM + k] = 0.0f;
	//			veles [j*SM_DIM + k] = 0.0f;
	//		}

	//		mass[j] = 0.0f;
	//		fix [j] = false;
	//		indx[j] = 0;
	//	}
	//}

	////ホスト側のデータを転送した結果をダンプ
	//ofstream ofs( "DtoH_Test.txt" );
	//ofs << "DtoH_Test" << endl;

	////デバイス側のデータをホストへ分割して転送
	//for(int i = 0; i < MAXCLUSTER; i++)
	//{
	//	ofs << "クラスタ " << i << endl;

	//	int vecSize = i * MAXPARTICLE * SM_DIM;
	//	Ice_SM* sm = ice_sm[i];

	//	cudaMemcpy(oPoses,	d_OrgPos+vecSize,	sizeof(float) * MAXPARTICLE * SM_DIM, cudaMemcpyDeviceToHost);
	//	cudaMemcpy(cPoses,	d_CurPos+vecSize,	sizeof(float) * MAXPARTICLE * SM_DIM, cudaMemcpyDeviceToHost);
	//	cudaMemcpy(veles,	d_Vel   +vecSize,	sizeof(float) * MAXPARTICLE * SM_DIM, cudaMemcpyDeviceToHost);

	//	int smSize = i * MAXPARTICLE;

	//	cudaMemcpy(mass,	d_Mass   +smSize,	sizeof(float) * MAXPARTICLE, cudaMemcpyDeviceToHost);
	//	cudaMemcpy(fix,		d_Fix    +smSize,	sizeof(bool)  * MAXPARTICLE, cudaMemcpyDeviceToHost);
	//	cudaMemcpy(indx,	d_PIndxes+smSize,	sizeof(int)   * MAXPARTICLE, cudaMemcpyDeviceToHost);

	//	cudaMemcpy(set,		d_IndxSet+i*2,		sizeof(int)   * 2, cudaMemcpyDeviceToHost);

	//	//転送できているかの確認
	//	for(int j = 0; j < MAXPARTICLE; j++)
	//	{
	//		//ofs << j << " :: "; 
	//		
	//		//初期位置
	//		//ofs << "oPoses = ( "; 
	//		//
	//		//for(int k = 0; k < SM_DIM; k++)
	//		//{
	//		//	ofs << oPoses[j*SM_DIM + k] << " ";
	//		//}
	//		//ofs << "), ";

	//		////mass[j] = j;
	//		////fix[j] = false;
	//		//if(indx[j] == -1){	ofs << "end" << endl; break;		}
	//		//else{				ofs << "indx = " << indx[j] << ",";	}

	//		//ofs << endl;
	//	}

	//	//開始位置と終了位置
	//	ofs << "set = ( " << set[0] << " ," << set[1] << ")" << endl;

	//	//差が１であればよい
	//	//cout << "start = " << set[0] << ", end = " << set[1] << ", Size = " << set[1]-set[0] << ", sm = " << sm->GetNumVertices();
	//	//cout << endl;
	//}

//転送直後のデータを比べてみる
	//for(int i = 0; i < MAXCLUSTER;i++)
	//{
	//	ice_sm[i]->CopyDeviceToInstance(i);
	//}

	delete[] oPoses;
	delete[] cPoses;
	delete[] veles;

	delete[] orgCm;
	delete[] curCm;

	delete[] mass;
	delete[] fix;
	delete[] indx;
	delete[] set;
}

void Ice_SM::AddVertex(const Vec3 &pos, double mass, int pIndx)
{
	rxShapeMatching::AddVertex(pos, mass, pIndx);

	//最大添字番号の更新
	if(m_iNumVertices > m_iIndxNum)
	{
		m_iIndxNum = m_iNumVertices;
	}

	//m_iLinearDeformation.push_back(0);
	//m_iVolumeConservation.push_back(0);

	//重心の更新
	CalcCm();
	m_vec3PreCm = m_vec3OrgCm;

	//変形行列Aqqの更新
	Vec3 p, q;
	for(int i = 0; i < m_iIndxNum; ++i)
	{
		if( CheckHole(i) ){	continue;	}

		q = Vec3(m_pOrgPos[i*SM_DIM+0], m_pOrgPos[i*SM_DIM+1], m_pOrgPos[i*SM_DIM+2])-m_vec3OrgCm;
		double m = m_pMass[i];

		m_mtrx3AqqInv(0,0) += m*q[0]*q[0];
		m_mtrx3AqqInv(0,1) += m*q[0]*q[1];
		m_mtrx3AqqInv(0,2) += m*q[0]*q[2];
		m_mtrx3AqqInv(1,0) += m*q[1]*q[0];
		m_mtrx3AqqInv(1,1) += m*q[1]*q[1];
		m_mtrx3AqqInv(1,2) += m*q[1]*q[2];
		m_mtrx3AqqInv(2,0) += m*q[2]*q[0];
		m_mtrx3AqqInv(2,1) += m*q[2]*q[1];
		m_mtrx3AqqInv(2,2) += m*q[2]*q[2];
	}

	m_mtrx3AqqInv = m_mtrx3AqqInv.Inverse();

	////Qの追加
	//m_vvec3OrgQ.resize(m_iNumVertices);
	//for(int i = 0; i < m_iNumVertices;++i)
	//{
	//	m_vvec3OrgQ[i] = m_vOrgPos[i]-m_vec3OrgCm;
	//}

	//前フレームの位置更新
	for(int i = 0; i < m_iIndxNum; ++i)
	{
		m_pPrePos[i*SM_DIM+0] = m_pCurPos[i*SM_DIM+0];
		m_pPrePos[i*SM_DIM+1] = m_pCurPos[i*SM_DIM+1];
		m_pPrePos[i*SM_DIM+2] = m_pCurPos[i*SM_DIM+2];
	}
}

void Ice_SM::Remove(int indx)
{
	if(m_iNumVertices <= 0) return;
	m_iNumVertices--;

	m_iPIndxes[indx] = MAXINT;

	for(int i = 0; i < SM_DIM; i++)
	{
		m_pOrgPos[indx*SM_DIM+i] = 0.0;
		m_pCurPos[indx*SM_DIM+i] = 0.0;
		m_pVel[indx*SM_DIM+i] = 0.0;
	}

	m_pMass[indx] = 0.0;
	m_pFix[indx] = false;

	m_ipLayeres[indx] = -1;

	m_fpAlphas[indx] = -1;
	m_fpBetas[indx] = -1;

	//m_iLinearDeformation.erase	( m_iLinearDeformation	.begin() + indx);
	//m_iVolumeConservation.erase	( m_iVolumeConservation	.begin() + indx);
}

void Ice_SM::Clear()
{
	rxShapeMatching::Clear();

	m_iIndxNum = 0;

	for(int i = 0; i < MAXPARTICLE; i++)
	{
		m_fpAlphas[i] = 0.0f;
		m_fpBetas[i] = 0.0f;

		m_ipLayeres[i] = -1;
	}

	//m_iLinearDeformation	.clear();
	//m_iVolumeConservation	.clear();
}

bool Ice_SM::CheckIndx(int pIndx)
{//	cout << __FUNCTION__ << endl;

	//vector<int>::iterator check = find( m_iPIndxes.begin(), m_iPIndxes.end(), pIndx);
	//if( check == m_iPIndxes.end() )	return false;
	//return true;

	for(int i = 0; i < m_iIndxNum; i++)
	{
		if(m_iPIndxes[i] == pIndx)
		{
			return true;
		}
	}

	return false;
}

int	Ice_SM::SearchIndx(int pIndx)
{
	////vector findバージョン　めちゃくちゃ重い
	//vector<int>::iterator begin = m_iPIndxes.begin();
	//vector<int>::iterator end = m_iPIndxes.end();
	//vector<int>::iterator check = find(begin, end, pIndx);

	//if(check == end)	return MAXINT;
	//return m_iPIndxes.size() - (end - check);

	////vector binary_searchバージョン
	//ソートできないのでだめ
	//vector<int>::iterator begin = m_iPIndxes.begin();
	//vector<int>::iterator end = m_iPIndxes.end();

	//if(!binary_search(begin, end, pIndx)) return MAXINT;
	//return *( std::lower_bound(begin, end, pIndx) );

	//配列バージョン
	for(int i = 0; i < m_iIndxNum; i++)
	{
		if(m_iPIndxes[i] == pIndx){	return i;	}
	}

	return MAXINT;

	////配列＋algorizm 以下の方法ではできない
	//int* last = m_iPIndxes + (int)m_iIndxNum;
	//int* result = find(m_iPIndxes, last, pIndx);

	//if(last == result) return MAXINT;
	//return *result;
}

/*
 *	添字が穴かどうかのチェック
 */
bool Ice_SM::CheckHole(int oIndx)
{
	return (m_iPIndxes[oIndx] == MAXINT);
}

void Ice_SM::ResetFinalParamPointer(unsigned clusterNum)
{
	for(unsigned cIndx = 0; cIndx < clusterNum; cIndx++)
	{
		for(unsigned dim = 0; dim < SM_DIM; dim++)
		{
			s_pfSldPos[cIndx*SM_DIM+dim] = 0.0f;
			s_pfSldVel[cIndx*SM_DIM+dim] = 0.0f;
		}
	}
}

/*!
 * Shape Matching法
 *  - 目標位置を計算して，m_vNewPosをその位置に近づける
 *  - 各粒子にαとβを持たせたバージョン
 * @param[in] dt タイムステップ幅
 */
void Ice_SM::ShapeMatchingUsePath()
{
	if(m_iNumVertices <= 1) return;

	double mass = 0.0;	// 総質量

	// 重心座標の計算
	for(int i = 0; i < m_iIndxNum; ++i){
		if( CheckHole(i) ){	continue;	}
		//if(m_pFix[i]){	mass += m_pMass[i]*300;	}	// 固定点の質量を大きくする
		mass += m_pMass[i];
	}

	Vec3 cm(1 / mass * m_vec3NowCm);	// 現在の重心
	Vec3 q(0.0);

	//Apqの行列式を求め，反転するかを判定
	//不安定な場合が多いので×
	if( m_mtrx3Apq.Determinant() < 0.0 && m_iNumVertices >= 4)
	{
		//cout << "before det < 0" << endl;
		//１　符号を反転
		m_mtrx3Apq(0,2) = -m_mtrx3Apq(0,2);
		m_mtrx3Apq(1,2) = -m_mtrx3Apq(1,2);
		m_mtrx3Apq(2,2) = -m_mtrx3Apq(2,2);

		//２　a2とa3を交換
		//double tmp;
		//tmp = Apq(0,2);
		//Apq(0,2) = Apq(0,1);
		//Apq(0,1) = tmp;

		//tmp = Apq(1,2);
		//Apq(1,2) = Apq(1,1);
		//Apq(1,1) = tmp;

		//tmp = Apq(2,2);
		//Apq(2,2) = Apq(2,1);
		//Apq(2,1) = tmp;
	}

	rxMatrix3 R, S;
	//PolarDecomposition(m_mtrx3Apq, R, S, m_mtrxBeforeU);
	PolarDecomposition(m_mtrx3Apq, R, S);

	if(m_bLinearDeformation)
	{
		//// Linear Deformations
		//rxMatrix3 A(m_mtrx3Apq * m_mtrx3AqqInv);	// A = Apq*Aqq^-1

		////体積保存のために√(det(A))で割る
		//if(m_bVolumeConservation){
		//	double det = fabs(A.Determinant());
		//	if(det > RX_FEQ_EPS){
		//		det = 1.0/sqrt(det);
		//		if(det > 2.0) det = 2.0;
		//		A *= det;
		//	}
		//}

		//// 目標座標を計算し，現在の頂点座標を移動
		//for(int i = 0; i <  m_iIndxNum; ++i)
		//{
		//	if( CheckHole(i) ){	continue;	}
		//	if(m_pFix[i]) continue;

		//	// 回転行列Rの代わりの行列RL=βA+(1-β)Rを計算
		//	int pIndx = m_iPIndxes[i] * SM_DIM;
		//	int cIndx = i*SM_DIM;

		//	for(int j = 0; j < SM_DIM; j++)
		//	{
		//		q[j] = m_pOrgPos[cIndx+j]-m_vec3OrgCm[j];
		//	}

		//	Vec3 gp(R*q+cm);

		//	for(int j = 0; j < SM_DIM; j++)
		//	{
		//		m_pCurPos[cIndx+j] += (gp[j]-m_pCurPos[cIndx+j]) * m_fpAlphas[i];
		//	}
		//}

		m_fDefAmount = 0.0;

		// 目標座標を計算し，現在の頂点座標を移動
		for(int i = 0; i < m_iIndxNum; ++i)
		{
			if( CheckHole(i) ){	continue;	}
			if(m_pFix[i]) continue;

			int cIndx = i*SM_DIM;

			// 回転行列Rの代わりの行列RL=βA+(1-β)Rを計算
			for(int j = 0; j < SM_DIM; j++)
			{
				q[j] = m_pOrgPos[cIndx+j]-m_vec3OrgCm[j];
			}

			Vec3 gp(R*q+cm);

			for(int j = 0; j < SM_DIM; j++)
			{
				float defAmount = (gp[j]-m_pCurPos[cIndx+j]) * m_fpAlphas[i];
				m_pCurPos[cIndx+j] += defAmount;
				m_fDefAmount += abs(defAmount);
			}
		}
	}

	//境界条件　これがないと現在位置が境界を無視するので，粒子がすり抜ける
	for(int i = 0; i < m_iIndxNum; ++i)
	{
		clamp(m_pCurPos, i*SM_DIM);
	}
}

/*!
 * 外力
 *  - 重力と境界壁からの力の影響
 * @param[in] dt タイムステップ幅
 */
void Ice_SM::calExternalForces()
{
	double dt = m_dDt;

	//最終位置・速度をクラスタの各粒子に反映
	// 重力の影響を付加，速度を反映
	for(int i = 0; i < m_iIndxNum; ++i)
	{
		if( CheckHole(i) ){	continue;	}
		int pIndx = m_iPIndxes[i]*4;
		int cIndx = i*SM_DIM;

		for(int j = 0; j < SM_DIM; j++)
		{
			int jpIndx = pIndx+j;
			int jcIndx = cIndx+j;
			
			m_pCurPos[jcIndx] = s_pfPrtPos[jpIndx] + s_pfPrtVel[jpIndx]*dt;
		}
	}

	// 境界壁の影響
	//処理がかなり重くなるが，安定はするみたい
	double res = 0.9;	// 反発係数
	for(int i = 0; i < m_iIndxNum; ++i)
	{
		if( CheckHole(i) ){	continue;	}
		
		////if(m_pFix[i]) continue;
		//Vec3 &p = m_pCurPos[i];

		//Vec3 &v = m_pVel[i];
		//if(np[0] < m_v3Min[0] || np[0] > m_v3Max[0]){
		//	np[0] = p[0]-v[0]*dt*res;
		//	np[1] = p[1];
		//	np[2] = p[2];
		//}
		//if(np[1] < m_v3Min[1] || np[1] > m_v3Max[1]){
		//	np[1] = p[1]-v[1]*dt*res;
		//	np[0] = p[0] ;
		//	np[2] = p[2];
		//}
		//if(np[2] < m_v3Min[2] || np[2] > m_v3Max[2]){
		//	np[2] = p[2]-v[2]*dt*res;
		//	np[0] = p[0];
		//	np[1] = p[1];
		//}

		clamp(m_pCurPos, i*SM_DIM);
	}
}


//----------------------------------------------ソリッド版---------------------------------------------
/*!
 * Shape Matching法
 *  - 目標位置を計算して，m_vNewPosをその位置に近づける
 *  - 各粒子にαとβを持たせたバージョン
 * @param[in] dt タイムステップ幅
 */
void Ice_SM::ShapeMatchingSolid()
{
	if(m_iNumVertices <= 1) return;

	rxMatrix3 Apq(0.0), Aqq(0.0);
	Vec3 p, q;
	Vec3 cm(0.0), cm_org(0.0);	// 重心
	double mass = 0.0;	// 総質量
	rxMatrix3 R, S;

	//前フレームと現フレームの位置の差の大きさから質量を決定
	//CalcMass();

	// 重心座標の計算
	for(int i = 0; i < m_iIndxNum; ++i)
	{
		if( CheckHole(i) ){	continue;	}

		double m = m_pMass[i];
		if(m_pFix[i]) m *= 300.0;	// 固定点の質量を大きくする

		int cIndx = i*SM_DIM;
		mass += m;

		for(int j = 0; j < SM_DIM; j++)
		{
			cm[j] += m_pCurPos[cIndx+j]*m;
		}
	}

	cm /= mass;
	m_vec3NowCm = cm;
	cm_org = m_vec3OrgCm;

	// Apq = Σmpq^T
	// Aqq = Σmqq^T
	for(int i = 0; i < m_iIndxNum; ++i)
	{
		if( CheckHole(i) ){	continue;	}

		int cIndx = i*SM_DIM;

		for(int j = 0; j < SM_DIM; j++)
		{
			p[j] = m_pCurPos[cIndx+j]-cm[j];
			q[j] = m_pOrgPos[cIndx+j]-cm_org[j];
		}

		double m = m_pMass[i];

		Apq(0,0) += m*p[0]*q[0];
		Apq(0,1) += m*p[0]*q[1];
		Apq(0,2) += m*p[0]*q[2];
		Apq(1,0) += m*p[1]*q[0];
		Apq(1,1) += m*p[1]*q[1];
		Apq(1,2) += m*p[1]*q[2];
		Apq(2,0) += m*p[2]*q[0];
		Apq(2,1) += m*p[2]*q[1];
		Apq(2,2) += m*p[2]*q[2];

		//Aqq(0,0) += m*q[0]*q[0];
		//Aqq(0,1) += m*q[0]*q[1];
		//Aqq(0,2) += m*q[0]*q[2];
		//Aqq(1,0) += m*q[1]*q[0];
		//Aqq(1,1) += m*q[1]*q[1];
		//Aqq(1,2) += m*q[1]*q[2];
		//Aqq(2,0) += m*q[2]*q[0];
		//Aqq(2,1) += m*q[2]*q[1];
		//Aqq(2,2) += m*q[2]*q[2];
	}

	////Apqの行列式を求め，反転するかを判定
	////不安定な場合が多いので×
	if( Apq.Determinant() < 0.0 && m_iNumVertices >= 4)
	{
		//cout << "before det < 0" << endl;
		//１　符号を反転
		Apq(0,2) = -Apq(0,2);
		Apq(1,2) = -Apq(1,2);
		Apq(2,2) = -Apq(2,2);

		//２　a2とa3を交換
		//double tmp;
		//tmp = Apq(0,2);
		//Apq(0,2) = Apq(0,1);
		//Apq(0,1) = tmp;

		//tmp = Apq(1,2);
		//Apq(1,2) = Apq(1,1);
		//Apq(1,1) = tmp;

		//tmp = Apq(2,2);
		//Apq(2,2) = Apq(2,1);
		//Apq(2,1) = tmp;
	}

	//PolarDecomposition(Apq, R, S, m_mtrxBeforeU);
	PolarDecomposition(Apq, R, S);

	if(m_bLinearDeformation)
	{
		//// 体積保存のために√(det(A))で割る
		//if(m_bVolumeConservation){
		//	//beta == 0を前提にしているので，Aを使わないとしている
		//	rxMatrix3 A(Apq*Aqq.Inverse());
		//	double det = fabs(A.Determinant());

		//	if(det > RX_FEQ_EPS){
		//		det = 1.0/sqrt(det);
		//		if(det > 2.0) det = 2.0;
		//		A *= det;
		//	}
		//}

		m_fDefAmount = 0.0f;

		// 目標座標を計算し，現在の頂点座標を移動
		for(int i = 0; i < m_iIndxNum; ++i)
		{
			if( CheckHole(i) ){	continue;	}

			if(m_pFix[i]) continue;

			int cIndx = i*SM_DIM;

			// 回転行列Rの代わりの行列RL=βA+(1-β)Rを計算
			for(int j = 0; j < SM_DIM; j++)
			{
				q[j] = m_pOrgPos[cIndx+j]-cm_org[j];
			}

			Vec3 gp(R*q+cm);

			for(int j = 0; j < SM_DIM; j++)
			{
				float defAmount = (gp[j]-m_pCurPos[cIndx+j]) * m_fpAlphas[i];
				m_pCurPos[cIndx+j] += defAmount;

				m_fDefAmount += abs(defAmount);
			}
		}
	}
}

void Ice_SM::ShapeMatchingSelected(int selected)
{
	if(m_iNumVertices <= 1) return;

	rxMatrix3 Apq(0.0), Aqq(0.0);
	Vec3 p, q;
	Vec3 cm(0.0), cm_org(0.0);	// 重心
	double mass = 0.0;			// 総質量
	rxMatrix3 R, S;

	// 重心座標の計算
	for(int i = 0; i < m_iIndxNum; ++i)
	{
		if( CheckHole(i) ){	continue;	}

		double m = m_pMass[i];
		if(m_pFix[i]) m *= 300.0;	// 固定点の質量を大きくする

		int cIndx = i*SM_DIM;

		mass += m;

		for(int j = 0; j < SM_DIM; j++)
		{
			cm[j] += m_pCurPos[cIndx+j]*m;
		}
	}

	cm /= mass;
	m_vec3NowCm = cm;
	cm_org = m_vec3OrgCm;

	// Apq = Σmpq^T
	// Aqq = Σmqq^T
	for(int i = 0; i < m_iIndxNum; ++i)
	{
		if( CheckHole(i) ){	continue;	}

		int cIndx = i*SM_DIM;

		for(int j = 0; j < SM_DIM; j++)
		{
			p[j] = m_pCurPos[cIndx+j]-cm[j];
			q[j] = m_pOrgPos[cIndx+j]-cm_org[j];
		}

		double m = m_pMass[i];

		Apq(0,0) += m*p[0]*q[0];
		Apq(0,1) += m*p[0]*q[1];
		Apq(0,2) += m*p[0]*q[2];
		Apq(1,0) += m*p[1]*q[0];
		Apq(1,1) += m*p[1]*q[1];
		Apq(1,2) += m*p[1]*q[2];
		Apq(2,0) += m*p[2]*q[0];
		Apq(2,1) += m*p[2]*q[1];
		Apq(2,2) += m*p[2]*q[2];
	}

	PolarDecomposition(Apq, R, S);

	m_fDefAmount = 0.0f;

	if(m_bLinearDeformation)
	{
		// 目標座標を計算し，現在の頂点座標を移動
		for(int i = 0; i < m_iIndxNum; ++i)
		{
			if( CheckHole(i) ){	continue;	}
			if(m_pFix[i]) continue;

			int cIndx = i*SM_DIM;

			// 回転行列Rの代わりの行列RL=βA+(1-β)Rを計算
			for(int j = 0; j < SM_DIM; j++)
			{
				q[j] = m_pOrgPos[cIndx+j]-cm_org[j];
			}

			Vec3 gp(R*q+cm);

			for(int j = 0; j < SM_DIM; j++)
			{
				float defAmount = (gp[j]-m_pCurPos[cIndx+j]) * m_fpAlphas[i];

				m_pCurPos[cIndx+j] += defAmount;
				m_fDefAmount += abs(defAmount);
			}
		}
	}
}

//----------------------------------------------ソリッド版---------------------------------------------

/*!
 * 速度と位置の更新
 *  - 新しい位置と現在の位置座標から速度を算出
 * @param[in] dt タイムステップ幅
 */
void Ice_SM::integrate(double dt)
{
	double dt1 = 1.0/dt;
	for(int i = 0; i < m_iIndxNum; ++i)
	{
		if( CheckHole(i) ){	continue;	}

		int pIndx = m_iPIndxes[i]*4;

		for(int j = 0; j < SM_DIM; j++)
		{
			int cIndx = i*SM_DIM+j;
			m_pVel[cIndx] = (m_pCurPos[cIndx] - s_pfPrtPos[pIndx+j]) * dt1 + m_v3Gravity[j] * dt * 0.1f;	//TODO::0.1fかけるとうまくいって見える
		}
	}

	//前フレームの位置ベクトルを更新
	for(int i = 0; i < m_iIndxNum; ++i)
	{
		UpdatePrePos(i);
	}

	//前フレームの重心ベクトルを更新
	m_vec3PreCm = m_vec3NowCm;
}

void Ice_SM::CalcDisplaceMentVectorCm()
{
	Vec3 preVec(0.0);
	Vec3 nowVec(0.0);
	double mass = 0.0;

	for(int i = 0; i < m_iIndxNum; i++)
	{
		if( CheckHole(i) ){	continue;	}

		int pIndx = m_iPIndxes[i] * 4;
		int cIndx = i*SM_DIM;

		preVec += Vec3(s_pfPrtPos[pIndx+0], s_pfPrtPos[pIndx+1], s_pfPrtPos[pIndx+2]);
		nowVec += Vec3(m_pGoalPos[cIndx+0], m_pGoalPos[cIndx+1], m_pGoalPos[cIndx+2]);
		mass += m_pMass[i];
	}

	preVec /= mass;
	nowVec /= mass;

	m_vec3DisCm = nowVec-preVec;
}

/*!
 * 最終位置を各クラスタに反映
 * @param[in] dt タイムステップ幅
 */
void Ice_SM::calExternalForcesIteration()
{
	// クラスタ内の粒子の位置を，粒子の位置で更新
	for(int i = 0; i < m_iIndxNum; ++i)
	{
		if( CheckHole(i) ){	continue;	}
		int pIndx = m_iPIndxes[i]*SM_DIM;
		int cIndx = i*SM_DIM;

		for(int j = 0; j < SM_DIM; j++)
		{
			int jpIndx = pIndx+j;
			int jcIndx = cIndx+j;

			m_pCurPos[jcIndx] = s_pfSldPos[jpIndx];
		}
	}

	// 境界壁の影響
	//処理がかなり重くなるが，安定する
	double res = 0.9;	// 反発係数
	for(int i = 0; i < m_iIndxNum; ++i)
	{
	//	//if( CheckHole(i) ){	continue;	}
	//	////if(m_pFix[i]) continue;
	//	//Vec3 &p = m_pCurPos[i];

	//	//Vec3 &v = m_pVel[i];
	//	//if(np[0] < m_v3Min[0] || np[0] > m_v3Max[0]){
	//	//	np[0] = p[0]-v[0]*dt*res;
	//	//	np[1] = p[1];
	//	//	np[2] = p[2];
	//	//}
	//	//if(np[1] < m_v3Min[1] || np[1] > m_v3Max[1]){
	//	//	np[1] = p[1]-v[1]*dt*res;
	//	//	np[0] = p[0] ;
	//	//	np[2] = p[2];
	//	//}
	//	//if(np[2] < m_v3Min[2] || np[2] > m_v3Max[2]){
	//	//	np[2] = p[2]-v[2]*dt*res;
	//	//	np[0] = p[0];
	//	//	np[1] = p[1];
	//	//}

		clamp(m_pCurPos, i*SM_DIM);
	}
}


//----------------------------------------------ソリッド版---------------------------------------------
/*!
 * Shape Matching法
 *  - 目標位置を計算して，m_vNewPosをその位置に近づける
 *  - 各粒子にαとβを持たせたバージョン
 * @param[in] dt タイムステップ幅
 */
void Ice_SM::ShapeMatchingIteration()
{
	if(m_iNumVertices <= 1) return;

	rxMatrix3 Apq(0.0), Aqq(0.0);
	Vec3 p, q;
	Vec3 cm(0.0), cm_org(0.0);	// 重心
	double mass = 0.0;	// 総質量
	rxMatrix3 R, S;

	// 重心座標の計算
	for(int i = 0; i < m_iIndxNum; ++i)
	{
		if( CheckHole(i) ){	continue;	}

		double m = m_pMass[i];

		if(m_pFix[i]) m *= 300.0;	// 固定点の質量を大きくする
		mass += m;

		int pIndx = m_iPIndxes[i] * SM_DIM;

		for(int j = 0; j < SM_DIM; j++)
		{
			cm[j] += s_pfSldPos[pIndx+j]*m;
		}
	}

	cm /= mass;
	m_vec3NowCm = cm;
	cm_org = m_vec3OrgCm;

	// Apq = Σmpq^T
	// Aqq = Σmqq^T
	for(int i = 0; i < m_iIndxNum; ++i)
	{
		if( CheckHole(i) ){	continue;	}

		int pIndx = m_iPIndxes[i] * SM_DIM;
		int cIndx = i*SM_DIM;

		for(int j = 0; j < SM_DIM; j++)
		{
			p[j] = s_pfSldPos[pIndx+j]-cm[j];
			q[j] = m_pOrgPos[cIndx+j]-cm_org[j];
		}

		double m = m_pMass[i];

		Apq(0,0) += m*p[0]*q[0];
		Apq(0,1) += m*p[0]*q[1];
		Apq(0,2) += m*p[0]*q[2];
		Apq(1,0) += m*p[1]*q[0];
		Apq(1,1) += m*p[1]*q[1];
		Apq(1,2) += m*p[1]*q[2];
		Apq(2,0) += m*p[2]*q[0];
		Apq(2,1) += m*p[2]*q[1];
		Apq(2,2) += m*p[2]*q[2];

		//Aqq(0,0) += m*q[0]*q[0];
		//Aqq(0,1) += m*q[0]*q[1];
		//Aqq(0,2) += m*q[0]*q[2];
		//Aqq(1,0) += m*q[1]*q[0];
		//Aqq(1,1) += m*q[1]*q[1];
		//Aqq(1,2) += m*q[1]*q[2];
		//Aqq(2,0) += m*q[2]*q[0];
		//Aqq(2,1) += m*q[2]*q[1];
		//Aqq(2,2) += m*q[2]*q[2];
	}

	//Apqの行列式を求め，反転するかを判定
	//不安定な場合が多いので×
	if( Apq.Determinant() < 0.0 && m_iNumVertices >= 4)
	{
		//cout << "before det < 0" << endl;
		//１　符号を反転
		Apq(0,2) = -Apq(0,2);
		Apq(1,2) = -Apq(1,2);
		Apq(2,2) = -Apq(2,2);

		////２　a2とa3を交換
		//double tmp;
		//tmp = Apq(0,2);
		//Apq(0,2) = Apq(0,1);
		//Apq(0,1) = tmp;

		//tmp = Apq(1,2);
		//Apq(1,2) = Apq(1,1);
		//Apq(1,1) = tmp;

		//tmp = Apq(2,2);
		//Apq(2,2) = Apq(2,1);
		//Apq(2,1) = tmp;
	}

	//PolarDecomposition(Apq, R, S, m_mtrxBeforeU);
	PolarDecomposition(Apq, R, S);

	if(m_bLinearDeformation)
	{
		//// Linear Deformations
		//rxMatrix3 A;
		//A = Apq*Aqq.Inverse();	// A = Apq*Aqq^-1

		//// 体積保存のために√(det(A))で割る
		//if(m_bVolumeConservation){
		//	double det = fabs(A.Determinant());
		//	if(det > RX_FEQ_EPS){
		//		det = 1.0/sqrt(det);
		//		if(det > 2.0) det = 2.0;
		//		A *= det;
		//	}
		//}

		//cout << "計測開始2" << endl;
		//qc.Start();

		m_fDefAmount = 0.0f;

		// 目標座標を計算し，現在の頂点座標を移動
		for(int i = 0; i < m_iIndxNum; ++i)
		{
			if( CheckHole(i) ){	continue;	}

			if(m_pFix[i]) continue;

			int cIndx = i*SM_DIM;

			// 回転行列Rの代わりの行列RL=βA+(1-β)Rを計算
			for(int j = 0; j < SM_DIM; j++)
			{
				q[j] = m_pOrgPos[cIndx+j]-cm_org[j];
			}

			Vec3 gp(R*q+cm);

			for(int j = 0; j < SM_DIM; j++)
			{
				float defAmount = (gp[j]-m_pCurPos[cIndx+j]) * m_fpAlphas[i];
				m_pCurPos[cIndx+j] += defAmount;

				m_fDefAmount += abs(defAmount);
			}
		}
	}

	//境界処理
	for(int i = 0; i < m_iIndxNum; ++i)
	{
		clamp(m_pCurPos, i*SM_DIM);
	}
}

/*!
 * 速度と位置の更新
 *  - 新しい位置と現在の位置座標から速度を算出
 * @param[in] dt タイムステップ幅
 */
void Ice_SM::integrateIteration()
{
	double dt1 = 1.0 / m_dDt;
	for(int i = 0; i < m_iIndxNum; ++i)
	{
		if( CheckHole(i) ){	continue;	}
		int pIndx = m_iPIndxes[i]*4;

		for(int j = 0; j < SM_DIM; j++)
		{
			int cIndx = i*SM_DIM+j;
			m_pVel[cIndx] = (m_pCurPos[cIndx] - s_pfPrtPos[pIndx+j]) * dt1 + m_v3Gravity[j] * m_dDt * 0.1f;	//TODO::パラメータ0.1fを使わなくてもいいように
		}
	}

	//前フレームの位置ベクトルを更新
	for(int i = 0; i < m_iIndxNum; ++i)
	{
		UpdatePrePos(i);
	}

	//前フレームの重心ベクトルを更新
	m_vec3PreCm = m_vec3NowCm;
}

//前ステップと現ステップの位置から質量を決定
void Ice_SM::CalcMass()
{
	float distSum = 0.0f;
	for(int i = 0; i < m_iIndxNum; i++)
	{
		float dist =	abs(m_pCurPos[i*SM_DIM+0] - m_pPrePos[i*SM_DIM+0])
					+	abs(m_pCurPos[i*SM_DIM+1] - m_pPrePos[i*SM_DIM+1])
					+	abs(m_pCurPos[i*SM_DIM+2] - m_pPrePos[i*SM_DIM+2]);
		
		distSum += dist*dist;
	}

	float average = 1.0f/(float)m_iIndxNum;
	float massSum = 0.0f;
	
	//総移動量と変形量を閾値とする
	if(m_fDefAmount > 0.1f)
	{
		for(int i = 0; i < m_iIndxNum; i++)
		{
			float dist =	abs(m_pCurPos[i*SM_DIM+0] - m_pPrePos[i*SM_DIM+0])
						+	abs(m_pCurPos[i*SM_DIM+1] - m_pPrePos[i*SM_DIM+1])
						+	abs(m_pCurPos[i*SM_DIM+2] - m_pPrePos[i*SM_DIM+2]);
			
			m_pMass[i] = 1.0f + dist * dist/distSum - average;
			massSum += m_pMass[i];

			//cout << "mass = " << m_pMass[i] << endl;
		}

		//cout << "massSum = " << massSum << " , vertexNum = " << m_iIndxNum << endl;
	}
	else
	{
		for(int i = 0; i < m_iIndxNum; i++)
		{
			m_pMass[i] = 1.0f;
		}
	}
}

//現在の重心を計算
void Ice_SM::CalcCm()
{
	//重心の更新
	double massSum = 0.0;	// 総質量
	m_vec3OrgCm = Vec3(0.0);

	// 重心座標の計算
	for(int i = 0; i < m_iIndxNum;++i)
	{
		if( CheckHole(i) ){	continue;	}

		double m = m_pMass[i];
		int pIndx = i*SM_DIM;

		massSum += m;
		m_vec3OrgCm += Vec3(m_pOrgPos[pIndx+0], m_pOrgPos[pIndx+1], m_pOrgPos[pIndx+2]) * m;
	}

	m_vec3OrgCm /= massSum;
}

void Ice_SM::UpdatePrePos(int pIndx)
{
	pIndx = pIndx*SM_DIM;

	m_pPrePos[pIndx+0] = m_pCurPos[pIndx+0];
	m_pPrePos[pIndx+1] = m_pCurPos[pIndx+1];
	m_pPrePos[pIndx+2] = m_pCurPos[pIndx+2];
}

void Ice_SM::SetPrePos(int pIndx, const Vec3& nowPos)
{
	pIndx = pIndx*SM_DIM;
	m_pPrePos[pIndx+0] = nowPos[0];
	m_pPrePos[pIndx+1] = nowPos[1];
	m_pPrePos[pIndx+2] = nowPos[2];
}

/*!
 *　CPUで運動計算
 */
void Ice_SM::UpdateCPU()
{
	calExternalForces();		//現在位置の計算
	ShapeMatchingSolid();		//現在位置の更新 普通の計算
	integrate(m_dDt);			//速度の計算

	//CalcDisplaceMentVectorCm();	//重心移動量
}

/*!
 *　パスを用いた高速化手法　CPUで運動計算
 */
void Ice_SM::UpdateUsePathCPU()
{
	calExternalForces();			//現在位置の計算
	ShapeMatchingUsePath();			//現在位置の更新 パスを用いた高速化
	integrate(m_dDt);				//速度の計算

	//CalcDisplaceMentVectorCm();	//重心移動量
}

/*
 *	GPUを用いた運動計算
 */
void Ice_SM::UpdateGPU()
{
	LaunchShapeMatchingGPU(s_vertNum, sd_PrtPos, sd_PrtVel, d_OrgPos, d_CurPos, d_OrgCm, d_CurCm, d_Vel, d_PIndxes, d_IndxSet, 0.01);
}

void Ice_SM::UpdateUsePathGPU()
{
	LaunchShapeMatchingUsePathGPU(s_vertNum, sd_PrtPos, sd_PrtVel, d_OrgPos, d_CurPos, d_OrgCm, d_CurCm, d_Apq, d_Vel, d_PIndxes, d_IndxSet, 0.01);
}

/*
 *	GPUを用いた運動計算　反復処理あり
 */
void Ice_SM::UpdateIterationGPU(float* sldPos, float* sldVel)
{
	LaunchShapeMatchingIterationGPU(s_vertNum, sd_PrtPos, sd_PrtVel, sldPos, sldVel, d_OrgPos, d_CurPos, d_OrgCm, d_CurCm, d_Vel, d_PIndxes, d_IndxSet, 0.01);
}

//GPUの計算結果を各インスタンスへコピーする
void Ice_SM::CopyDeviceToInstance(int num)
{	//cout << __FUNCTION__ << " num = " << num << endl;

	//デバイス側のデータをホスト側へ転送
	float* cPoses = new float[MAXPARTICLE * SM_DIM];
	float* veles  = new float[MAXPARTICLE * SM_DIM];

	int vecSize = num * MAXPARTICLE * SM_DIM;

	cudaMemcpy(cPoses,	d_CurPos+vecSize,	sizeof(float) * MAXPARTICLE * SM_DIM, cudaMemcpyDeviceToHost);
	cudaMemcpy(veles,	d_Vel   +vecSize,	sizeof(float) * MAXPARTICLE * SM_DIM, cudaMemcpyDeviceToHost);

	for(int j = 0; j < MAXPARTICLE; j++)
	{
		if(GetParticleIndx(j) == MAXINT){	continue;	}

		//位置，速度
		Vec3 cPos = Vec3(cPoses[j*SM_DIM + 0], cPoses[j*SM_DIM + 1], cPoses[j*SM_DIM + 2]);
		Vec3 vel  = Vec3(veles [j*SM_DIM + 0], veles [j*SM_DIM + 1], veles [j*SM_DIM + 2]);

		SetCurrentPos(j, cPos);
		SetVelocity(j, vel);
	}

	delete[] cPoses;
	delete[] veles;
}

//----------------------------------------デバッグ--------------------------------------
void Ice_SM::DebugIndx()
{	cout << "DebugIndx Indx = ";
	for(int i = 0; i < m_iIndxNum; i++)
	{
		cout << " " << m_iPIndxes[i];
	}
	cout << endl;
}

void Ice_SM::DebugLayer()
{	cout << "DebugLayer layer =";
	for(int i = 0; i < m_iIndxNum; i++)
	{
		cout << " " << m_ipLayeres[i];
	}
	cout << endl;
}