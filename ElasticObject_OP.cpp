/*!
  @file ElasticObject_OP.h
	
  @brief Oriented Particle法による弾性体変形                                                                                                                                                                                                                                                                                                                                                                                                                                                        
  @ref M. Muller et al., "Solid Simulation with Oriented Particles", SIGGRAPH2011. 
  参考にしたページ，プログラム http://yuki-koyama.hatenablog.com/entry/2015/01/30/150419
 
  @author Ryo Nakasone
  @date 2015-4
*/
//数値計算ライブラリ
//TODO: なぜかここに置かないとエラーが出る　他のヘッダをincludeする前に置かないといけないみたい？
#include <Eigen\Core>
#include <Eigen\Dense>
#include <Eigen\Geometry>
using namespace Eigen;

#include "ElasticObject_OP.h"

typedef OrientedParticleBaseElasticObject OrientedPrtObj;

const float* OrientedPrtObj::s_pfPrtPos = 0;
const float* OrientedPrtObj::s_pfPrtVel = 0;

float* OrientedPrtObj::s_pfClstrPos = 0;
float* OrientedPrtObj::s_pfClstrVel = 0;

const unsigned X = 0;
const unsigned Y = 1;
const unsigned Z = 2;

/*!
 * コンストラクタ
 */
OrientedParticleBaseElasticObject::OrientedParticleBaseElasticObject(int obj) : rxShapeMatching(obj)
{
	m_mtrx3Apq = rxMatrix3(0.0);

	m_fDefAmount = 0.0f;

	//AllocateMemory(prtNum);
}

/*!
 * コピーコンストラクタ
 */
OrientedParticleBaseElasticObject::OrientedParticleBaseElasticObject(const OrientedPrtObj& copy) : rxShapeMatching(copy)
{
	Copy(copy);
}

/*!
 * デストラクタ
 */
OrientedParticleBaseElasticObject::~OrientedParticleBaseElasticObject()
{
	ReleaseMemory();
}

//代入演算子でコピー
OrientedPrtObj& OrientedPrtObj::operator=(const OrientedPrtObj& copy)
{
	if(this != &copy){
		rxShapeMatching::ReleaseMemory();
		ReleaseMemory();

		rxShapeMatching::Copy(copy);
		Copy(copy);
	}

	return *this;
}

//メモリ解放
void OrientedPrtObj::ReleaseMemory()
{

}

//データのコピー
void OrientedPrtObj::Copy(const OrientedPrtObj& copy)
{
	m_mtrx3Apq = copy.Apq();
	m_iIndxNum = copy.GetIndxNum();
	m_fDefAmount = copy.DefAmount();
}

void OrientedPrtObj::AllocateStaticMemory(const float* prtPos, const float* prtVel, int vrtxNum)
{
	s_pfPrtPos = prtPos;
	s_pfPrtVel = prtVel;

	s_pfClstrPos = new float[vrtxNum*SM_DIM];
	s_pfClstrVel = new float[vrtxNum*SM_DIM];

	if(s_pfClstrPos == 0 || s_pfClstrVel == 0){
		cout << __FUNCTION__ << " Error!" << endl;
		return;
	}

	for(int i = 0; i < vrtxNum; ++i){
		int pIndx = i*4;
		int cIndx = i*SM_DIM;

		for(int j = 0; j < SM_DIM; j++){
			s_pfClstrPos[cIndx+j] = s_pfPrtPos[pIndx+j];
			s_pfClstrVel[cIndx+j] = s_pfPrtVel[pIndx+j];
		}
	}
}

void OrientedPrtObj::AddParticle(OrientedParticle* orientedPrt)
{
	Vec3 pos = orientedPrt->OrgPos();
	double mass = orientedPrt->Mass();
	int id = orientedPrt->Id();

	m_vOrientedPrtes.push_back(orientedPrt);

	AddParticle(pos, mass, id);
}

void OrientedPrtObj::AddParticle(const Vec3 &pos, double mass, int pIndx)
{
	rxShapeMatching::AddVertex(pos, mass, pIndx);

	//最大添字番号の更新
	if(m_iNumVertices > m_iIndxNum){	m_iIndxNum = m_iNumVertices;	}

	//重心位置の更新
	CalcOrgCm();
}

void OrientedPrtObj::RemoveParticle(int indx)
{
	rxShapeMatching::Remove(indx);
}

void OrientedPrtObj::ClearParticle()
{
	int indx = m_iIndxNum;

	rxShapeMatching::Clear();
}

//主成分分析を行い，クラスタを構成する粒子配置から初期姿勢を決定
void OrientedPrtObj::InitOrientation()
{
	CalcOrgCm();
	Vector3f ave = Vector3f(m_vec3OrgCm[X], m_vec3OrgCm[Y], m_vec3OrgCm[Z]);
	Matrix3f covarianceMtrx = Matrix3f::Zero();

	for(int i = 0; i < m_iIndxNum; i++){
		if(CheckHole(i)){	continue;	}

		int pIndx = i * SM_DIM;
		Vector3f pos = Vector3f(m_pOrgPos[pIndx+X], m_pOrgPos[pIndx+Y], m_pOrgPos[pIndx+Z]);
		Vector3f tmp = pos - ave;

		covarianceMtrx.row(0) += tmp(0) * tmp;
		covarianceMtrx.row(1) += tmp(1) * tmp;
		covarianceMtrx.row(2) += tmp(2) * tmp;
	}

	//EigenSolverだとコンパイルで内部エラーが出る
	//EigenSolver<Matrix3f> es;
	//es.compute(covarianceMtrx, true);

	SelfAdjointEigenSolver<Matrix3f> es;
	es.compute(covarianceMtrx);

	//computeの結果はソートされていないので，固有値の大きい順にソート
	//eigenのVectorはstlと相性が悪いみたい
	//なぜかvector<pair>に入ったVector3dを出力しようとするとコンパイルが通らない
	vector<pair<float, int>> eigenValVec;
	eigenValVec.push_back(pair<float, int>(es.eigenvalues()[0], 0));
	eigenValVec.push_back(pair<float, int>(es.eigenvalues()[1], 1));
	eigenValVec.push_back(pair<float, int>(es.eigenvalues()[2], 2));

	sort(eigenValVec.begin(), eigenValVec.end());

	//固有値，固有ベクトルで楕円の大きさと姿勢を決定
	//ソートした結果から，固有値の小さい順に固有ベクトルを取り出す
	//楕円の各軸の径は固有値から求める

	//本当はeigenのvector3dを使いたいが，OrientedParticleでは先生のVec3を使っているので仕方ない
	//一応，(1.0, 1.0, 1.0)の時に少しでも合わせるために正規化して√3倍してみた
	Vec3 radius = Vec3(0.0);
	Vec3 norm_eigenVal = Vec3(eigenValVec[X].first, eigenValVec[Y].first, eigenValVec[Z].first);
	radius = Unit(norm_eigenVal) * sqrt(3.0f);

	m_vOrientedPrtes[0]->ElipsoidRadius(radius);

	//最大固有値の方向が初期姿勢となる
	int indx = eigenValVec[Z].second;
	Vector3f orientedVec = es.eigenvectors().row(indx);

	//１	複数のベクトルから姿勢を定義する方法がわからない．２つからなら定義できるみたいだけど…　行列に直せる？
	//２	仕方ないので，座標系が回転したと考える．
	//		(1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)の正規直交基底が
	//		固有ベクトルとなるような回転行列を最小二乗法を用いて求める



	////Test
	//Matrix3d testMatrix;
	//testMatrix <<	1.0,0.446,-0.56,

	//				0.446,1.0,-0.239,

	//				-0.56,0.239,1.0;

	//SelfAdjointEigenSolver<Matrix3d> es;
	//es.compute(testMatrix);

	//cout << "Input matrix:" << endl

 //      << testMatrix << endl

 //      << endl

 //      << "Eigenvalues:" << endl

 //      << es.eigenvalues() << endl

 //      << endl

 //      << "Eigenvectors:"<< endl

 //      << es.eigenvectors() << endl;

	//vector<pair<float, int>> eigenValVec;
	//eigenValVec.push_back(pair<float, int>(es.eigenvalues()[1], 1));
	//eigenValVec.push_back(pair<float, int>(es.eigenvalues()[2], 2));
	//eigenValVec.push_back(pair<float, int>(es.eigenvalues()[0], 0));

	//sort(eigenValVec.begin(), eigenValVec.end());

	////テスト
	//Vector3f testVector = Vector3f::Zero();
	//testVector << 1, 2, 3;
	//Matrix3f testMatrix = Matrix3f::Zero();

	//testMatrix.row(0) += testVector(0) * testVector;
	//testMatrix.row(1) += testVector(1) * testVector;
	//testMatrix.row(2) += testVector(2) * testVector;

	//cout << "Test" << endl;
	//cout << testMatrix << endl;
}

void OrientedPrtObj::UpdateCluster_OP()
{
	float dt = m_dDt;

	////CalcForce(dt);
	////CalcVelocity(dt);
	////DampVelocity(dt);
	ProjectConstraint(dt);
	//DistanceConstraint(dt);
	////ApplyEstimatedPosition(dt);
	////CorrectVelocity(dt);

	////従来のシェイプマッチング
	//ShapeMatchingNormal();
}

//うまくいかなかった
void OrientedPrtObj::UpdateCluster_Sampling(const IceStructure* ice_struct)
{
	float dt = m_dDt;

	//更新位置をクラスタの各粒子に反映　各クラスタの各粒子は統一される
	for(int i = 0; i < m_iIndxNum; ++i){
		if( CheckHole(i) ){	continue;	}

		int pIndx = m_iPIndxes[i];
		int cIndx = i * SM_DIM;
		Vec3 prdPos = m_vOrientedPrtes[i]->PrdPos();

		for(int j = 0; j < SM_DIM; j++){
			m_pCurPos[cIndx+j] = prdPos[j];
		}
	}

	// 境界壁の影響
	//処理がかなり重くなるが，安定はするみたい
	double res = 0.9;	// 反発係数
	for(int i = 0; i < m_iIndxNum; ++i){
		if( CheckHole(i) ){	continue;	}
		if(m_pFix[i]) continue;

		int pIndx = m_iPIndxes[i];
		int cIndx = i*SM_DIM;
		Vec3 prdPos = m_vOrientedPrtes[i]->PrdPos();

		if(m_pCurPos[cIndx+0] < m_v3Min[0] || m_pCurPos[cIndx+0] > m_v3Max[0]){
			m_pCurPos[cIndx+0] = prdPos[0] - (s_pfClstrVel[cIndx+0] + (m_v3Gravity[0] * dt))*dt*res;
			m_pCurPos[cIndx+1] = prdPos[1];
			m_pCurPos[cIndx+2] = prdPos[2];
		}

		if(m_pCurPos[cIndx+1] < m_v3Min[1] || m_pCurPos[cIndx+1] > m_v3Max[1]){
			m_pCurPos[cIndx+1] = prdPos[1] - (s_pfClstrVel[cIndx+1] + (m_v3Gravity[1] * dt))*dt*res;
			m_pCurPos[cIndx+0] = prdPos[0] ;
			m_pCurPos[cIndx+2] = prdPos[2];
		}

		if(m_pCurPos[cIndx+2] < m_v3Min[2] || m_pCurPos[cIndx+2] > m_v3Max[2]){
			m_pCurPos[cIndx+2] = prdPos[2] - (s_pfClstrVel[cIndx+2] + (m_v3Gravity[2] * dt))*dt*res;
			m_pCurPos[cIndx+0] = prdPos[0];
			m_pCurPos[cIndx+1] = prdPos[1];
		}

		clamp(m_pCurPos, cIndx);
	}

	//弾性体の制約　シェイプマッチング計算
	CalcNowCm();							//クラスタの重心を更新
	
	////近傍クラスタの回転行列を補間して自身の回転行列を求める
	//rxMatrix3 R_intrp = InterpolateRotation(ice_struct);

	//rxMatrix3 Apq(0.0), Aqq(0.0);
	//CalcClusterMomentMatrix(Apq, Aqq);		//パーティクルの姿勢を考慮したモーメントマトリックス

	//rxMatrix3 R_polor, S;
	//PolarDecomposition(Apq, R_polor, S);

	////if(m_iObjectNo < 5){
	////	cout << m_iObjectNo << endl;
	////	cout << "R_intrp:\n" << R_intrp << endl;
	////	cout << "R_polor:\n" << R_polor << endl;
	////}

	//rxMatrix3 R = R_intrp;

	//近傍クラスタの変形勾配テンソルからApqを近似
	rxMatrix3 Apq = InterpolateApq(ice_struct);
	rxMatrix3 R = Apq;

	//rxMatrix3 R_polor, S;
	//PolarDecomposition(Apq, R_polor, S);
	//rxMatrix3 R = R_polor;

	//補間した回転行列を用いて運動計算
	if(R.Determinant() < 0.0f){
		R = -1.0f * R;
	}

	// 目標座標を計算し，現在の頂点座標を移動
	m_fDefAmount = 0.0f;

	for(int i = 0; i < m_iIndxNum; ++i){
		if( CheckHole(i) ){	continue;	}
		if(m_pFix[i]) continue;

		int cIndx = i*SM_DIM;

		Vec3 q;
		for(int j = 0; j < SM_DIM; j++){
			q[j] = m_pOrgPos[cIndx+j]-m_vec3OrgCm[j];
		}

		Vec3 gp(R*q+m_vec3NowCm);

		for(int j = 0; j < SM_DIM; j++){
			float defAmount = (gp[j]-m_pCurPos[cIndx+j]);
			m_pCurPos[cIndx+j] += defAmount;

			m_fDefAmount += abs(defAmount);
		}
	}

	//姿勢のために回転行列を保存
	m_vOrientedPrtes[0]->Rotation(R);
}

//初期重心位置
void OrientedPrtObj::CalcOrgCm()
{
	//重心の更新
	float massSum = 0.0f;	// 総質量
	m_vec3OrgCm = Vec3(0.0);

	// 重心座標の計算
	for(int i = 0; i < m_iIndxNum;++i){
		if(CheckHole(i)){	continue;	}

		float m = m_pMass[i];
		int pIndx = i*SM_DIM;

		massSum += m;
		m_vec3OrgCm += Vec3(m_pOrgPos[pIndx+0], m_pOrgPos[pIndx+1], m_pOrgPos[pIndx+2]) * m;
	}

	m_vec3OrgCm /= massSum;
}

//現在の重心位置
void OrientedPrtObj::CalcNowCm()
{
	//重心の更新
	float massSum = 0.0f;	// 総質量
	m_vec3NowCm = Vec3(0.0);

	// 重心座標の計算
	for(int i = 0; i < m_iIndxNum;++i){
		if(CheckHole(i)){	continue;	}

		float m = m_pMass[i];
		if(m_pFix[i]) m *= 300.0;	// 固定点の質量を大きくする
		
		int pIndx = i*SM_DIM;

		massSum += m;
		m_vec3NowCm += Vec3(m_pCurPos[pIndx+0], m_pCurPos[pIndx+1], m_pCurPos[pIndx+2]) * m;
	}

	m_vec3NowCm /= massSum;
}

//クラスタに含まれている粒子の姿勢・回転行列から自身の回転行列を補間
//姿勢を補間するだけでよいかもしれないが，一応回転行列で求める
//うまくいかなかった
rxMatrix3 OrientedPrtObj::InterpolateRotation(const IceStructure* ice_struct)
{
	float weight = 0.5f;	//まずは平均
	mk_Quaternion tmp_q = m_vOrientedPrtes[0]->CurOrientation();
	Quaternionf myQuat = Quaternionf(tmp_q.w, tmp_q.x, tmp_q.y, tmp_q.z);
	//Quaternionf myQuat = Quaternionf::Identity();
	rxMatrix3 R = rxMatrix3::Identity();

	//i = 0は自分なので除く
	for(int i = 1; i < m_iIndxNum; i++){
		if(CheckHole(i)){	continue;	}
		int pIndx = m_iPIndxes[i];
		if(ice_struct->GetMotionCalcCluster(pIndx) == 0){	continue;	}

		rxMatrix3 rotateMtrx = m_vOrientedPrtes[i]->Rotation();

		//回転行列→クォータニオン
		Matrix<float, 3, 3, RowMajor> tmp_r1;
		tmp_r1 <<	rotateMtrx(0, 0), rotateMtrx(0, 1), rotateMtrx(0, 2), 
					rotateMtrx(1, 0), rotateMtrx(1, 1), rotateMtrx(1, 2), 
					rotateMtrx(2, 0), rotateMtrx(2, 1), rotateMtrx(2, 2);

		Quaternionf rotateQuat(tmp_r1);

		//クォータニオン同士で球面線形補間
		Quaternionf rotateSlerp = myQuat.slerp(weight, rotateQuat);

		//クォータニオン→回転行列
		Matrix<float, 3, 3, RowMajor> tmp_r2 = rotateSlerp.matrix();
		rxMatrix3 R_Slerp(
			tmp_r2(0, 0), tmp_r2(0, 1), tmp_r2(0, 2), 
			tmp_r2(1, 0), tmp_r2(1, 1), tmp_r2(1, 2), 
			tmp_r2(2, 0), tmp_r2(2, 1), tmp_r2(2, 2)
		);

		//総積計算
		R *= R_Slerp;
	}

	return R;
}

//近傍クラスタの変形勾配テンソルを用いてこのクラスタの変形勾配テンソルを近似
//うまくいかなかった
rxMatrix3 OrientedPrtObj::InterpolateApq(const IceStructure* ice_struct)
{
	float weight = 0.5f;	//てきとうに平均

	//回転成分
	rxMatrix3 R = rxMatrix3::Identity();

	mk_Quaternion tmp_q = m_vOrientedPrtes[0]->CurOrientation();
	//Quaternionf myQuat = Quaternionf(tmp_q.w, tmp_q.x, tmp_q.y, tmp_q.z);

	Quaternionf myQuat = Quaternionf::Identity();

	//i = 0は自分なので除く
	for(int i = 1; i < m_iIndxNum; i++){
		if(CheckHole(i)){	continue;	}
		int pIndx = m_iPIndxes[i];
		if(ice_struct->GetMotionCalcCluster(pIndx) == 0){	continue;	}

		rxMatrix3 rotateMtrx = m_vOrientedPrtes[i]->Rotation();

		//回転行列→クォータニオン
		Matrix<float, 3, 3, RowMajor> tmp_r1;
		tmp_r1 <<	rotateMtrx(0, 0), rotateMtrx(0, 1), rotateMtrx(0, 2), 
					rotateMtrx(1, 0), rotateMtrx(1, 1), rotateMtrx(1, 2), 
					rotateMtrx(2, 0), rotateMtrx(2, 1), rotateMtrx(2, 2);
		
		Quaternionf rotateQuat(tmp_r1);

		//クォータニオン同士で球面線形補間
		Quaternionf rotateSlerp = myQuat.slerp(weight, rotateQuat);

		//クォータニオン→回転行列
		Matrix<float, 3, 3, RowMajor> tmp_r2 = rotateSlerp.matrix();
		rxMatrix3 R_Slerp(
			tmp_r2(0, 0), tmp_r2(0, 1), tmp_r2(0, 2), 
			tmp_r2(1, 0), tmp_r2(1, 1), tmp_r2(1, 2), 
			tmp_r2(2, 0), tmp_r2(2, 1), tmp_r2(2, 2)
		);

		//総積計算
		R *= R_Slerp;
	}

	//回転以外の成分
	Matrix<float, 3, 3, RowMajor> S_temp;
	S_temp <<	0.0f, 0.0f, 0.0f,
				0.0f, 0.0f, 0.0f,
				0.0f, 0.0f, 0.0f;

	for(int i = 1; i < m_iIndxNum; i++){
		if(CheckHole(i)){	continue;	}
		int pIndx = m_iPIndxes[i];
		if(ice_struct->GetMotionCalcCluster(pIndx) == 0){	continue;	}

		rxMatrix3 sym = m_vOrientedPrtes[i]->Symmetric();
		
		Matrix<float, 3, 3, RowMajor> eigen_Sim;
		eigen_Sim <<	sym(0, 0), sym(0, 1), sym(0, 2), 
						sym(1, 0), sym(1, 1), sym(1, 2), 
						sym(2, 0), sym(2, 1), sym(2, 2);

		//指数行列ってこれでいいのか？
		Matrix<float, 3, 3, RowMajor> log_sim = eigen_Sim.array().log();
		S_temp += weight * log_sim;
	}
	//対数行列ってこれでいいのか？
	S_temp = S_temp.array().exp();

	rxMatrix3 S(
			S_temp(0, 0), S_temp(0, 1), S_temp(0, 2), 
			S_temp(1, 0), S_temp(1, 1), S_temp(1, 2), 
			S_temp(2, 0), S_temp(2, 1), S_temp(2, 2)
	);

	rxMatrix3 Apq =R * S;

	return Apq;
}


//シェイプマッチング法
void OrientedPrtObj::UpdateCluster_SM()
{
	double dt = m_dDt;

	//最終位置・速度をクラスタの各粒子に反映　各クラスタの各粒子は統一される
	for(int i = 0; i < m_iIndxNum; ++i)
	{
		if( CheckHole(i) ){	continue;	}
		int pIndx = m_iPIndxes[i]*4;
		int cIndx = i*SM_DIM;
		int indx = m_iPIndxes[i]*3;


		Vec3 prdPos = m_vOrientedPrtes[i]->PrdPos();
		Vec3 vel = m_vOrientedPrtes[i]->Velocity();

		for(int j = 0; j < SM_DIM; j++)
		{
			int jpIndx = pIndx+j;
			int jcIndx = cIndx+j;

			m_pCurPos[jcIndx] = s_pfPrtPos[jpIndx] + (s_pfPrtVel[jpIndx] + (m_v3Gravity[j] * dt)) *dt;
			//m_pCurPos[jcIndx] = s_pfPrtPos[jpIndx] + (s_pfClstrVel[indx+j] + (m_v3Gravity[j] * dt)) *dt;	//かわらない
			//m_pCurPos[jcIndx] = s_pfPrtPos[jpIndx] + (vel[j] + (m_v3Gravity[j] * dt)) *dt;	//なぜかできない
		}
	}

	// 境界壁の影響
	//処理がかなり重くなるが，安定はするみたい
	double res = 0.9;	// 反発係数
	for(int i = 0; i < m_iIndxNum; ++i)
	{
		if( CheckHole(i) ){	continue;	}
		if(m_pFix[i]) continue;

		int pIndx =  m_iPIndxes[i]*4;
		int cIndx = i*SM_DIM;
		Vec3 vel = m_vOrientedPrtes[i]->Velocity();

		if(m_pCurPos[cIndx+0] < m_v3Min[0] || m_pCurPos[cIndx+0] > m_v3Max[0]){
			m_pCurPos[cIndx+0] = s_pfPrtPos[pIndx+0] - (s_pfPrtVel[pIndx+0] + (m_v3Gravity[0] * dt))*dt*res;
			m_pCurPos[cIndx+1] = s_pfPrtPos[pIndx+1];
			m_pCurPos[cIndx+2] = s_pfPrtPos[pIndx+2];
		}

		if(m_pCurPos[cIndx+1] < m_v3Min[1] || m_pCurPos[cIndx+1] > m_v3Max[1]){
			m_pCurPos[cIndx+1] = s_pfPrtPos[pIndx+1] - (s_pfPrtVel[pIndx+1] + (m_v3Gravity[1] * dt))*dt*res;
			m_pCurPos[cIndx+0] = s_pfPrtPos[pIndx+0] ;
			m_pCurPos[cIndx+2] = s_pfPrtPos[pIndx+2];
		}

		if(m_pCurPos[cIndx+2] < m_v3Min[2] || m_pCurPos[cIndx+2] > m_v3Max[2]){
			m_pCurPos[cIndx+2] = s_pfPrtPos[pIndx+2] - (s_pfPrtVel[pIndx+2] + (m_v3Gravity[2] * dt))*dt*res;
			m_pCurPos[cIndx+0] = s_pfPrtPos[pIndx+0];
			m_pCurPos[cIndx+1] = s_pfPrtPos[pIndx+1];
		}

		clamp(m_pCurPos, cIndx);
	}

	if(m_iNumVertices <= 1) return;

	
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
	rxMatrix3 Apq(0.0), Aqq(0.0);
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
		p *= m;

		Apq(0,0) += p[0] * q[0];
		Apq(0,1) += p[0] * q[1];
		Apq(0,2) += p[0] * q[2];
		Apq(1,0) += p[1] * q[0];
		Apq(1,1) += p[1] * q[1];
		Apq(1,2) += p[1] * q[2];
		Apq(2,0) += p[2] * q[0];
		Apq(2,1) += p[2] * q[1];
		Apq(2,2) += p[2] * q[2];
	}


	PolarDecomposition(Apq, R, S);

	if(m_bLinearDeformation){
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
				float defAmount = (gp[j]-m_pCurPos[cIndx+j]);
				m_pCurPos[cIndx+j] += defAmount;

				m_fDefAmount += abs(defAmount);
			}
		}
	}

	m_vOrientedPrtes[0]->Rotation(R);

	double dt1 = 1.0/dt;
	for(int i = 0; i < m_iIndxNum; ++i)
	{
		if( CheckHole(i) ){	continue;	}

		int pIndx = m_iPIndxes[i]*4;

		for(int j = 0; j < SM_DIM; j++)
		{
			int cIndx = i*SM_DIM+j;
			m_pVel[cIndx] = (m_pCurPos[cIndx] - s_pfPrtPos[pIndx+j]) * dt1/* + m_v3Gravity[j] * dt*/;
		}
	}
}

void OrientedPrtObj::UpdateCluster_SM_Itr()
{
	double dt = m_dDt;

	//// 境界壁の影響
	////処理がかなり重くなるが，安定はするみたい
	//double res = 0.9;	// 反発係数
	//for(int i = 0; i < m_iIndxNum; ++i)
	//{
	//	if( CheckHole(i) ){	continue;	}
	//	if(m_pFix[i]) continue;

	//	int pIndx =  m_iPIndxes[i]*4;
	//	int cIndx = i*SM_DIM;

	//	if(m_pCurPos[cIndx+0] < m_v3Min[0] || m_pCurPos[cIndx+0] > m_v3Max[0]){
	//		m_pCurPos[cIndx+0] = s_pfPrtPos[pIndx+0] - (s_pfPrtVel[pIndx+0] + (m_v3Gravity[0] * dt))*dt*res;
	//		m_pCurPos[cIndx+1] = s_pfPrtPos[pIndx+1];
	//		m_pCurPos[cIndx+2] = s_pfPrtPos[pIndx+2];
	//	}

	//	if(m_pCurPos[cIndx+1] < m_v3Min[1] || m_pCurPos[cIndx+1] > m_v3Max[1]){
	//		m_pCurPos[cIndx+1] = s_pfPrtPos[pIndx+1] - (s_pfPrtVel[pIndx+1] + (m_v3Gravity[1] * dt))*dt*res;
	//		m_pCurPos[cIndx+0] = s_pfPrtPos[pIndx+0] ;
	//		m_pCurPos[cIndx+2] = s_pfPrtPos[pIndx+2];
	//	}

	//	if(m_pCurPos[cIndx+2] < m_v3Min[2] || m_pCurPos[cIndx+2] > m_v3Max[2]){
	//		m_pCurPos[cIndx+2] = s_pfPrtPos[pIndx+2] - (s_pfPrtVel[pIndx+2] + (m_v3Gravity[2] * dt))*dt*res;
	//		m_pCurPos[cIndx+0] = s_pfPrtPos[pIndx+0];
	//		m_pCurPos[cIndx+1] = s_pfPrtPos[pIndx+1];
	//	}

	//	clamp(m_pCurPos, cIndx);
	//}

	//更新位置をクラスタの各粒子に反映　各クラスタの各粒子は統一される
	for(int i = 0; i < m_iIndxNum; ++i){
		if( CheckHole(i) ){	continue;	}

		int pIndx = m_iPIndxes[i];
		int cIndx = i * SM_DIM;
		Vec3 prdPos = m_vOrientedPrtes[i]->PrdPos();

		for(int j = 0; j < SM_DIM; j++){
			m_pCurPos[cIndx+j] = prdPos[j];
		}
	}

	if(m_iNumVertices <= 1) return;

	Vec3 p, q;
	Vec3 cm(0.0), cm_org(0.0);	// 重心
	double mass = 0.0;	// 総質量
	rxMatrix3 R, S;

	CalcNowCm();							//クラスタの重心を更新

	if(m_iNumVertices <= 1) return;

	rxMatrix3 Apq(0.0), Aqq(0.0);

	cm = m_vec3NowCm;
	cm_org = m_vec3OrgCm;

	// Apq = Σmpq^T
	// Aqq = Σmqq^T
	for(int i = 0; i < m_iIndxNum; ++i)
	{
		if( CheckHole(i) ){	continue;	}

		int pIndx = m_iPIndxes[i] * SM_DIM;
		int cIndx = i*SM_DIM;

		Vec3 prdPos = m_vOrientedPrtes[i]->PrdPos();

		for(int j = 0; j < SM_DIM; j++)
		{
			p[j] = /*prdPos[j]*/m_pCurPos[cIndx+j]-cm[j];
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
	//TODO:不安定な場合が多いので×
	if(Apq.Determinant() < 0.0 && m_iNumVertices >= 10){
		//１　符号を反転
		Apq(0,2) = -Apq(0,2);
		Apq(1,2) = -Apq(1,2);
		Apq(2,2) = -Apq(2,2);
	}

	//PolarDecomposition(Apq, R, S, m_mtrxBeforeU); //warm start
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
				float defAmount = (gp[j]-m_pCurPos[cIndx+j]);
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

	m_vOrientedPrtes[0]->Rotation(R);
}

//距離制約による運動計算
void OrientedPrtObj::UpdateCluster_DC()
{
	double dt = m_dDt;

	//最終位置・速度をクラスタの各粒子に反映　各クラスタの各粒子は統一される
	for(int i = 0; i < m_iIndxNum; ++i)
	{
		if( CheckHole(i) ){	continue;	}
		int pIndx = m_iPIndxes[i]*4;
		int cIndx = i*SM_DIM;
		int indx = m_iPIndxes[i]*3;


		Vec3 prdPos = m_vOrientedPrtes[i]->PrdPos();
		Vec3 vel = m_vOrientedPrtes[i]->Velocity();

		for(int j = 0; j < SM_DIM; j++)
		{
			int jpIndx = pIndx+j;
			int jcIndx = cIndx+j;

			m_pCurPos[jcIndx] = s_pfPrtPos[jpIndx] + (s_pfPrtVel[jpIndx] + (m_v3Gravity[j] * dt)) *dt;
			//m_pCurPos[jcIndx] = s_pfPrtPos[jpIndx] + (s_pfClstrVel[indx+j] + (m_v3Gravity[j] * dt)) *dt;	//かわらない
			//m_pCurPos[jcIndx] = s_pfPrtPos[jpIndx] + (vel[j] + (m_v3Gravity[j] * dt)) *dt;	//なぜかできない
		}
	}

	// 境界壁の影響
	//処理がかなり重くなるが，安定はするみたい
	double res = 0.9;	// 反発係数
	for(int i = 0; i < m_iIndxNum; ++i)
	{
		if( CheckHole(i) ){	continue;	}
		if(m_pFix[i]) continue;

		int pIndx =  m_iPIndxes[i]*4;
		int cIndx = i*SM_DIM;
		Vec3 vel = m_vOrientedPrtes[i]->Velocity();

		if(m_pCurPos[cIndx+0] < m_v3Min[0] || m_pCurPos[cIndx+0] > m_v3Max[0]){
			m_pCurPos[cIndx+0] = s_pfPrtPos[pIndx+0] - (s_pfPrtVel[pIndx+0] + (m_v3Gravity[0] * dt))*dt*res;
			m_pCurPos[cIndx+1] = s_pfPrtPos[pIndx+1];
			m_pCurPos[cIndx+2] = s_pfPrtPos[pIndx+2];
		}

		if(m_pCurPos[cIndx+1] < m_v3Min[1] || m_pCurPos[cIndx+1] > m_v3Max[1]){
			m_pCurPos[cIndx+1] = s_pfPrtPos[pIndx+1] - (s_pfPrtVel[pIndx+1] + (m_v3Gravity[1] * dt))*dt*res;
			m_pCurPos[cIndx+0] = s_pfPrtPos[pIndx+0] ;
			m_pCurPos[cIndx+2] = s_pfPrtPos[pIndx+2];
		}

		if(m_pCurPos[cIndx+2] < m_v3Min[2] || m_pCurPos[cIndx+2] > m_v3Max[2]){
			m_pCurPos[cIndx+2] = s_pfPrtPos[pIndx+2] - (s_pfPrtVel[pIndx+2] + (m_v3Gravity[2] * dt))*dt*res;
			m_pCurPos[cIndx+0] = s_pfPrtPos[pIndx+0];
			m_pCurPos[cIndx+1] = s_pfPrtPos[pIndx+1];
		}

		clamp(m_pCurPos, cIndx);
	}

	if(m_iNumVertices <= 1) return;

	
	Vec3 cm(0.0), cm_org(0.0);	// 重心
	double mass = 0.0;	// 総質量
	
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

	//距離制約
	Vec3 mainPrtPos_cur = GetVertexPos(0);
	Vec3 mainPrtPos_org = GetOrgPos(0);

	//重心との距離制約　重心位置は固定
	Vec3 nowDx_center = m_vec3NowCm - mainPrtPos_cur;
	float nowDist_center = norm(nowDx_center);

	Vec3 orgDx_center = m_vec3OrgCm - mainPrtPos_org;
	float orgDist_center = norm(orgDx_center);

	Vec3 dirVec_center = nowDist_center > 0.0f ? nowDx_center / nowDist_center : Vec3(0.0f);

	mainPrtPos_cur += dirVec_center * (orgDist_center - nowDist_center);
	SetCurrentPos(0, mainPrtPos_cur);

	m_fDefAmount = 0.0f;

	//距離制約　クラスタに対応する粒子と，その他の粒子
	for(int i = 1; i < m_iIndxNum; i++){
		if(CheckHole(i)){	continue;	}
		if(m_pFix[i]){		continue;	}

		Vec3 vertexPos = GetVertexPos(i);

		Vec3 nowDx = vertexPos - mainPrtPos_cur;
		float nowDist = norm(nowDx);

		Vec3 orgDx = GetOrgPos(i) - mainPrtPos_org;
		float orgDist = norm(orgDx);

		Vec3 dirVec = nowDist > 0.0f ? nowDx / nowDist : Vec3(0.0f);
		Vec3 dx = 1.0f/2.0f * dirVec * (orgDist - nowDist);

		mainPrtPos_cur -= dx;
		SetCurrentPos(i, vertexPos + dx);

		m_fDefAmount += abs(orgDist - nowDist);
	}

	SetCurrentPos(0, mainPrtPos_cur);

	double dt1 = 1.0/dt;
	for(int i = 0; i < m_iIndxNum; ++i)
	{
		if( CheckHole(i) ){	continue;	}

		int pIndx = m_iPIndxes[i]*4;

		for(int j = 0; j < SM_DIM; j++)
		{
			int cIndx = i*SM_DIM+j;
			m_pVel[cIndx] = (m_pCurPos[cIndx] - s_pfPrtPos[pIndx+j]) * dt1/* + m_v3Gravity[j] * dt*/;
		}
	}
}

//反復処理　距離制約による運動計算
void OrientedPrtObj::UpdateCluster_DC_Itr()
{
	double dt = m_dDt;

	//// 境界壁の影響
	////処理がかなり重くなるが，安定はするみたい
	//double res = 0.9;	// 反発係数
	//for(int i = 0; i < m_iIndxNum; ++i)
	//{
	//	if( CheckHole(i) ){	continue;	}
	//	if(m_pFix[i]) continue;

	//	int pIndx =  m_iPIndxes[i]*4;
	//	int cIndx = i*SM_DIM;

	//	if(m_pCurPos[cIndx+0] < m_v3Min[0] || m_pCurPos[cIndx+0] > m_v3Max[0]){
	//		m_pCurPos[cIndx+0] = s_pfPrtPos[pIndx+0] - (s_pfPrtVel[pIndx+0] + (m_v3Gravity[0] * dt))*dt*res;
	//		m_pCurPos[cIndx+1] = s_pfPrtPos[pIndx+1];
	//		m_pCurPos[cIndx+2] = s_pfPrtPos[pIndx+2];
	//	}

	//	if(m_pCurPos[cIndx+1] < m_v3Min[1] || m_pCurPos[cIndx+1] > m_v3Max[1]){
	//		m_pCurPos[cIndx+1] = s_pfPrtPos[pIndx+1] - (s_pfPrtVel[pIndx+1] + (m_v3Gravity[1] * dt))*dt*res;
	//		m_pCurPos[cIndx+0] = s_pfPrtPos[pIndx+0] ;
	//		m_pCurPos[cIndx+2] = s_pfPrtPos[pIndx+2];
	//	}

	//	if(m_pCurPos[cIndx+2] < m_v3Min[2] || m_pCurPos[cIndx+2] > m_v3Max[2]){
	//		m_pCurPos[cIndx+2] = s_pfPrtPos[pIndx+2] - (s_pfPrtVel[pIndx+2] + (m_v3Gravity[2] * dt))*dt*res;
	//		m_pCurPos[cIndx+0] = s_pfPrtPos[pIndx+0];
	//		m_pCurPos[cIndx+1] = s_pfPrtPos[pIndx+1];
	//	}

	//	clamp(m_pCurPos, cIndx);
	//}

	//更新位置をクラスタの各粒子に反映　各クラスタの各粒子は統一される
	for(int i = 0; i < m_iIndxNum; ++i){
		if( CheckHole(i) ){	continue;	}

		int pIndx = m_iPIndxes[i];
		int cIndx = i * SM_DIM;
		Vec3 prdPos = m_vOrientedPrtes[i]->PrdPos();

		for(int j = 0; j < SM_DIM; j++){
			m_pCurPos[cIndx+j] = prdPos[j];
		}
	}

	if(m_iNumVertices <= 1) return;

	Vec3 cm(0.0), cm_org(0.0);	// 重心
	double mass = 0.0;	// 総質量
	
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

	//距離制約
	Vec3 mainPrtPos_cur = GetVertexPos(0);
	Vec3 mainPrtPos_org = GetOrgPos(0);

	//重心との距離制約　重心位置は固定
	Vec3 nowDx_center = m_vec3NowCm - mainPrtPos_cur;
	float nowDist_center = norm(nowDx_center);

	Vec3 orgDx_center = m_vec3OrgCm - mainPrtPos_org;
	float orgDist_center = norm(orgDx_center);

	Vec3 dirVec_center = nowDist_center > 0.0f ? 1.0f / nowDist_center * nowDx_center : Vec3(0.0f);

	Vec3 fixPos = GetVertexPos(0) + dirVec_center * (orgDist_center - nowDist_center);
	SetCurrentPos(0, fixPos);

	mainPrtPos_cur = GetVertexPos(0);
	m_fDefAmount = 0.0f;

	//距離制約　クラスタに対応する粒子と，その他の粒子
	for(int i = 1; i < m_iIndxNum; i++){
		if(CheckHole(i)){	continue;	}
		if(m_pFix[i]){		continue;	}

		Vec3 nowDx = GetVertexPos(i) - mainPrtPos_cur;
		float nowDist = norm(nowDx);

		Vec3 orgDx = GetOrgPos(i) - mainPrtPos_org;
		float orgDist = norm(orgDx);

		Vec3 dirVec = nowDist > 0.0f ? 1.0f / nowDist * nowDx : Vec3(0.0f);
		Vec3 dx = 1.0f/2.0f * dirVec * (orgDist - nowDist);

		SetCurrentPos(0, mainPrtPos_cur - dx);
		SetCurrentPos(i, GetVertexPos(i) + dx);

		mainPrtPos_cur = GetVertexPos(0);

		m_fDefAmount += norm2(dx);
	}

	//境界処理
	for(int i = 0; i < m_iIndxNum; ++i)
	{
		clamp(m_pCurPos, i*SM_DIM);
	}
}

//各粒子に加わる力を更新　主には重力など
void OrientedPrtObj::CalcForce(float dt)
{
	//for(vector<Vec3>::iterator it = s_vvec3Force.begin(); it != s_vvec3Force.end(); it++){
	//	*it = Vec3(0.0);

	//	//重力を適用
	//	(*it) += m_v3Gravity;
	//}
}

void OrientedPrtObj::CalcVelocity(float dt)
{
}

void OrientedPrtObj::DampVelocity(float dt)
{
}

//マウスによるドラッグを反映させるために，無理やり値を更新
void OrientedPrtObj::CopyPrtToClstrPos(unsigned prtNum)
{
	int cIndx = 0;
	int pIndx = 0;
	int num = prtNum;

	//#pragma omp parallel for private(cIndx, pIndx)
	for(int i = 0; i < num; i++){
		cIndx = i * SM_DIM;
		pIndx = i * 4;

		s_pfClstrPos[cIndx+X] = s_pfPrtPos[pIndx+X];
		s_pfClstrVel[cIndx+X] = s_pfPrtVel[pIndx+X];

		s_pfClstrPos[cIndx+Y] = s_pfPrtPos[pIndx+Y];
		s_pfClstrVel[cIndx+Y] = s_pfPrtVel[pIndx+Y];

		s_pfClstrPos[cIndx+Z] = s_pfPrtPos[pIndx+Z];
		s_pfClstrVel[cIndx+Z] = s_pfPrtVel[pIndx+Z];
	}
}

void OrientedPrtObj::CalcClusterMomentMatrix(rxMatrix3& Apq, rxMatrix3& Aqq)
{
	rxMatrix3 Asum = rxMatrix3(0.0f);
	for(int i = 0; i < m_iIndxNum; ++i){
		if( CheckHole(i) ){	continue;	}

		OrientedParticle* prt = m_vOrientedPrtes[i];
		Asum += prt->A_elipsoid();
		Asum += prt->MomentMatrix();
	}

	rxMatrix3 MccT;

	MccT(0,0) = 1.0f * m_vec3NowCm[0] * m_vec3OrgCm[0];
	MccT(0,1) = 1.0f * m_vec3NowCm[0] * m_vec3OrgCm[1];
	MccT(0,2) = 1.0f * m_vec3NowCm[0] * m_vec3OrgCm[2];
	MccT(1,0) = 1.0f * m_vec3NowCm[1] * m_vec3OrgCm[0];
	MccT(1,1) = 1.0f * m_vec3NowCm[1] * m_vec3OrgCm[1];
	MccT(1,2) = 1.0f * m_vec3NowCm[1] * m_vec3OrgCm[2];
	MccT(2,0) = 1.0f * m_vec3NowCm[2] * m_vec3OrgCm[0];
	MccT(2,1) = 1.0f * m_vec3NowCm[2] * m_vec3OrgCm[1];
	MccT(2,2) = 1.0f * m_vec3NowCm[2] * m_vec3OrgCm[2];

	float massSum = (float)m_iIndxNum;

	Apq = Asum - massSum * MccT;

	//Aqq
	Vec3 q(0.0f);

	for(int i = 0; i < m_iIndxNum; ++i){
		if( CheckHole(i) ){	continue;	}
		Vec3 op(m_pOrgPos[i*SM_DIM+0], m_pOrgPos[i*SM_DIM+1], m_pOrgPos[i*SM_DIM+2]);

		q = op - m_vec3OrgCm;
		double m = m_pMass[i];

		Aqq(0,0) += m*q[0]*q[0];
		Aqq(0,1) += m*q[0]*q[1];
		Aqq(0,2) += m*q[0]*q[2];
		Aqq(1,0) += m*q[1]*q[0];
		Aqq(1,1) += m*q[1]*q[1];
		Aqq(1,2) += m*q[1]*q[2];
		Aqq(2,0) += m*q[2]*q[0];
		Aqq(2,1) += m*q[2]*q[1];
		Aqq(2,2) += m*q[2]*q[2];
	}
}

void OrientedPrtObj::ProjectConstraint(float dt)
{
	//更新位置をクラスタの各粒子に反映　各クラスタの各粒子は統一される
	for(int i = 0; i < m_iIndxNum; ++i){
		if( CheckHole(i) ){	continue;	}

		int pIndx = m_iPIndxes[i];
		int cIndx = i * SM_DIM;
		Vec3 prdPos = m_vOrientedPrtes[i]->PrdPos();

		for(int j = 0; j < SM_DIM; j++){
			m_pCurPos[cIndx+j] = prdPos[j];
		}
	}

	// 境界壁の影響
	//処理がかなり重くなるが，安定はするみたい
	double res = 0.9;	// 反発係数
	for(int i = 0; i < m_iIndxNum; ++i){
		if( CheckHole(i) ){	continue;	}
		if(m_pFix[i]) continue;

		int pIndx = m_iPIndxes[i];
		int cIndx = i*SM_DIM;
		Vec3 prdPos = m_vOrientedPrtes[i]->PrdPos();

		if(m_pCurPos[cIndx+0] < m_v3Min[0] || m_pCurPos[cIndx+0] > m_v3Max[0]){
			m_pCurPos[cIndx+0] = prdPos[0] - (s_pfClstrVel[pIndx * SM_DIM+0] + (m_v3Gravity[0] * dt))*dt*res;
			m_pCurPos[cIndx+1] = prdPos[1];
			m_pCurPos[cIndx+2] = prdPos[2];
		}

		if(m_pCurPos[cIndx+1] < m_v3Min[1] || m_pCurPos[cIndx+1] > m_v3Max[1]){
			m_pCurPos[cIndx+1] = prdPos[1] - (s_pfClstrVel[pIndx * SM_DIM+1] + (m_v3Gravity[1] * dt))*dt*res;
			m_pCurPos[cIndx+0] = prdPos[0] ;
			m_pCurPos[cIndx+2] = prdPos[2];
		}

		if(m_pCurPos[cIndx+2] < m_v3Min[2] || m_pCurPos[cIndx+2] > m_v3Max[2]){
			m_pCurPos[cIndx+2] = prdPos[2] - (s_pfClstrVel[pIndx * SM_DIM+2] + (m_v3Gravity[2] * dt))*dt*res;
			m_pCurPos[cIndx+0] = prdPos[0];
			m_pCurPos[cIndx+1] = prdPos[1];
		}

		clamp(m_pCurPos, cIndx);
	}

	//弾性体の制約　シェイプマッチング計算
	CalcNowCm();							//クラスタの重心を更新
	
	rxMatrix3 Apq(0.0), Aqq(0.0);
	CalcClusterMomentMatrix(Apq, Aqq);		//パーティクルの姿勢を考慮したモーメントマトリックス

	rxMatrix3 R, S;
	PolarDecomposition(Apq, R, S);
	if(R.Determinant() < 0.0f){
		R = -1.0f * R;
	}

	rxMatrix3 A(0.0);
	A = Apq*Aqq.Inverse();	// A = Apq*Aqq^-1

	// 目標座標を計算し，現在の頂点座標を移動
	m_fDefAmount = 0.0f;

	float alpha = 0.1f;
	float beta = 0.2f;

	rxMatrix3 RL = beta*A+(1.0f-beta)*R;

	for(int i = 0; i < m_iIndxNum; ++i){
		if( CheckHole(i) ){	continue;	}
		if(m_pFix[i]) continue;

		int cIndx = i*SM_DIM;

		Vec3 q;
		for(int j = 0; j < SM_DIM; j++){
			q[j] = m_pOrgPos[cIndx+j]-m_vec3OrgCm[j];
		}

		Vec3 gp(R*q+m_vec3NowCm);
		//Vec3 gp(RL*q+m_vec3NowCm);

		for(int j = 0; j < SM_DIM; j++){
			//float defAmount = (gp[j]-m_pCurPos[cIndx+j]) * alpha;
			float defAmount = (gp[j]-m_pCurPos[cIndx+j]);
			m_pCurPos[cIndx+j] += defAmount;

			m_fDefAmount += abs(defAmount);
		}
	}

	//姿勢のために回転行列を保存
	m_vOrientedPrtes[0]->Rotation(R);

	//補間のために回転以外の部分を含む行列を保存
	m_vOrientedPrtes[0]->Symmetric(S);
}

//距離制約
void OrientedPrtObj::DistanceConstraint(float dt)
{

}

void OrientedPrtObj::ApplyEstimatedPosition(float dt)
{
}

void OrientedPrtObj::CorrectVelocity(float dt)
{
}

//クラスタ内にある粒子間の内積の総和　クォータニオンでの内積
//全て同じ姿勢なら1，全て反対なら-1
float OrientedPrtObj::QuatInnerDotSum() const
{
	int prtNum = 1;
	float dot = 0.0f;
	mk_Quaternion baseQuat = m_vOrientedPrtes[0]->CurOrientation();

	for(unsigned i = 1; i < m_iIndxNum; i++){
		if( CheckHole(i) ){	continue;	}
		if(m_pFix[i]){		continue;	}

		mk_Quaternion quat = m_vOrientedPrtes[i]->CurOrientation();

		dot += baseQuat.x * quat.x + baseQuat.y * quat.y + baseQuat.z * quat.z + baseQuat.w * quat.w;

		prtNum++;
	}

	return dot/(float)prtNum;
}

float OrientedPrtObj::DirectionalVecInnerDotSum() const
{
	int prtNum = 1;
	float dot = 0.0f;
	Vec3 baseAxis = m_vOrientedPrtes[0]->Axis();

	for(unsigned i = 1; i < m_iIndxNum; i++){
		if( CheckHole(i) ){	continue;	}
		if(m_pFix[i]){		continue;	}

		Vec3 axis = m_vOrientedPrtes[i]->Axis();

		dot += baseAxis[0] * axis[0] + baseAxis[1] * axis[1] + baseAxis[2] * axis[2];

		prtNum++;
	}

	return dot/(float)prtNum;
}

float OrientedPrtObj::FinalDefAmount() const
{
	float defAmount = 0.0f;

	//最終位置の重心
	Vec3 finalCm(0.0f);
	for(unsigned i = 0; i < m_iIndxNum; i++){
		if( CheckHole(i) ){	continue;	}
		if(m_pFix[i]){		continue;	}

		//int pIndx = GetParticleIndx(i) * SM_DIM;
		//finalCm += Vec3(s_pfClstrPos[pIndx + X], s_pfClstrPos[pIndx + Y], s_pfClstrPos[pIndx + Z]);

		int pIndx = GetParticleIndx(i) * 4;
		finalCm += Vec3(s_pfPrtPos[pIndx + X], s_pfPrtPos[pIndx + Y], s_pfPrtPos[pIndx + Z]);
	}

	finalCm /= GetNumVertices();

	for(unsigned i = 1; i < m_iIndxNum; i++){
		if( CheckHole(i) ){	continue;	}
		if(m_pFix[i]){		continue;	}

		//初期位置を局所座標系に直して，回転行列を掛ける
		rxMatrix3 R = m_vOrientedPrtes[i]->Rotation();
		Vec3 orgPos = R * (GetOrgPos(i) - m_vec3OrgCm)/* + finalCm*/;

		//現在位置を局所座標系に変換
		//int pIndx = GetParticleIndx(i) * SM_DIM;
		//Vec3 finalPos = Vec3(s_pfClstrPos[pIndx + X], s_pfClstrPos[pIndx + Y], s_pfClstrPos[pIndx + Z])/* - finalCm*/;

		int pIndx = GetParticleIndx(i) * 4;
		Vec3 finalPos = Vec3(s_pfPrtPos[pIndx + X], s_pfPrtPos[pIndx + Y], s_pfPrtPos[pIndx + Z]) - finalCm;

		//局所座標系における移動量を変形量とする
		defAmount += norm2(orgPos-finalPos);
	}

	return defAmount / GetNumVertices();
}

float OrientedPrtObj::VolumeChange() const
{
	float defAmount = 0.0f;

	//最終位置の重心
	Vec3 finalCm(0.0f);
	for(unsigned i = 0; i < m_iIndxNum; i++){
		if( CheckHole(i) ){	continue;	}
		if(m_pFix[i]){		continue;	}

		int pIndx = GetParticleIndx(i) * 4;
		finalCm += Vec3(s_pfPrtPos[pIndx + X], s_pfPrtPos[pIndx + Y], s_pfPrtPos[pIndx + Z]);
	}

	finalCm /= GetNumVertices();

	rxMatrix3 Apq(0.0);
	Vec3 p(0.0f);
	Vec3 q(0.0f);

	for(int i = 0; i < m_iIndxNum; ++i)
	{
		if( CheckHole(i) ){	continue;	}

		int cIndx = i*SM_DIM;
		int pIndx = GetParticleIndx(i)*4;

		for(int j = 0; j < SM_DIM; j++)
		{
			p[j] = s_pfPrtPos[pIndx+j]-finalCm[j];
			q[j] = m_pOrgPos[cIndx+j]-m_vec3OrgCm[j];
		}

		double m = m_pMass[i];
		p *= m;

		Apq(0,0) += p[0] * q[0];
		Apq(0,1) += p[0] * q[1];
		Apq(0,2) += p[0] * q[2];
		Apq(1,0) += p[1] * q[0];
		Apq(1,1) += p[1] * q[1];
		Apq(1,2) += p[1] * q[2];
		Apq(2,0) += p[2] * q[0];
		Apq(2,1) += p[2] * q[1];
		Apq(2,2) += p[2] * q[2];
	}

	defAmount = Apq.Determinant();

	return defAmount;
}

bool OrientedPrtObj::IsThereOutlier() const
{
	float defAmount = 0.0f;
	int outlierNum = 0;

	//最終位置の重心
	Vec3 finalCm(0.0f);
	for(unsigned i = 0; i < m_iIndxNum; i++){
		if( CheckHole(i) ){	continue;	}
		if(m_pFix[i]){		continue;	}

		int pIndx = GetParticleIndx(i) * 4;
		finalCm += Vec3(s_pfPrtPos[pIndx + X], s_pfPrtPos[pIndx + Y], s_pfPrtPos[pIndx + Z]);
	}

	finalCm /= GetNumVertices();

	for(unsigned i = 1; i < m_iIndxNum; i++){
		if( CheckHole(i) ){	continue;	}
		if(m_pFix[i]){		continue;	}

		//初期位置を局所座標系に直して，回転行列を掛ける
		rxMatrix3 R = m_vOrientedPrtes[i]->Rotation();
		Vec3 orgPos = R * (GetOrgPos(i) - m_vec3OrgCm);

		//現在位置を局所座標系に変換
		int pIndx = GetParticleIndx(i) * 4;
		Vec3 finalPos = Vec3(s_pfPrtPos[pIndx + X], s_pfPrtPos[pIndx + Y], s_pfPrtPos[pIndx + Z]) - finalCm;

		//大きく移動した粒子があるかどうか
		if(norm2(orgPos-finalPos) > 0.007f){
			outlierNum++;
		}
	}

	//return outlierNum > 0;
	return outlierNum > (GetNumVertices()/3);
}