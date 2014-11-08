#include "IceObject.h"

float* IceObject::s_sphPrtPos;
float* IceObject::s_sphPrtVel;

float* IceObject::m_fInterPolationCoefficience;

//デバイスポインタ
float* IceObject::sd_sphPrtPos;
float* IceObject::sd_sphPrtVel;

float* IceObject::sd_sldPrtPos;	
float* IceObject::sd_sldPrtVel;

int IceObject::sm_particleNum;
int IceObject::sm_tetraNum;
int IceObject::sm_clusterNum;
int IceObject::sm_layerNum;
int IceObject::sm_maxParticleNum;

IceObject::IceObject(int pMaxNum, int cMaxNum, int tMaxNum, int prtNum, float* hp, float* hv, float* dp, float* dv, int layer, int maxParticleNum)
{
	InitIceObj(pMaxNum, cMaxNum, tMaxNum);
	SetParticleNum(prtNum);

	SetSPHHostPointer(hp, hv);						//CPU処理で用いるポインタの登録
	SetSPHDevicePointer(dp, dv);					//GPU処理で用いるポインタの登録

	SetSearchLayerNum(layer);						//探索レイヤー数
	SetMaxParticleNum(maxParticleNum);				//最大粒子数
}

IceObject::~IceObject()
{
}

//それぞれのクラス・変数の初期化
void IceObject::InitIceObj(int pMaxNum, int cMaxNum, int tMaxNum)
{	cout << __FUNCTION__ << endl;
	//物体の構造の初期化
	m_iceStrct = new IceStructure(pMaxNum, cMaxNum, tMaxNum, sm_layerNum);

	//運動計算を行うSMクラスタの初期化

	//補間処理のためのパラメータの初期化
	InitInterPolation();

}

void IceObject::InitTetra()
{	cout << __FUNCTION__ << endl;

	IceTetrahedra &tetra = IceTetrahedra::GetInstance();		//シングルトンな四面体管理クラスを用意してみた
	tetra.InitTetra(s_sphPrtPos, sm_particleNum);

cout << __FUNCTION__ << " check1" << endl;

	m_iceStrct->SetTetraNum(tetra.GetTetraListSize());			//現四面体数を登録

cout << __FUNCTION__ << " check2" << endl;

	//各四面体に含まれる粒子数のカウント
	for(unsigned i = 0; i < tetra.GetTetraListSize(); i++)
	{
		m_iceStrct->CountTetrahedra(i, tetra.GetTetraList(i));
	}

cout << __FUNCTION__ << " check3" << endl;

	//メモリ確保
	m_iceStrct->InitTetraInfo();

cout << __FUNCTION__ << " check4" << endl;

	//粒子が所属しているクラスタ数の配列をコピー
	int *PtoTNum = new int[sm_particleNum];

	for(int i = 0; i < sm_particleNum; i++)
	{
		PtoTNum[i] = m_iceStrct->GetPtoTNum(i);
	}

cout << __FUNCTION__ << " check5" << endl;

	//四面体データ登録
	for(unsigned i = 0; i < tetra.GetTetraListSize(); i++)
	{
		m_iceStrct->SetTetraInfo(i, PtoTNum);
	}
	delete[] PtoTNum;

cout << __FUNCTION__ << " check6" << endl;

	//近傍四面体データ登録
	for(unsigned i = 0; i < tetra.GetTetraListSize(); i++)
	{
		m_iceStrct->SetNeighborTetra(i, sm_layerNum);
	}

cout << __FUNCTION__ << " check7" << endl;

	//デバッグ
	//m_iceObj->DebugTetraInfo();
}

//クラスタの初期化
void IceObject::InitCluster(Vec3 boundarySpaceHigh, Vec3 boundarySpaceLow, float timeStep, int itr)
{	cout << __FUNCTION__ << endl;

	//変数の初期化
	for(vector<Ice_SM*>::iterator it = m_iceMove.begin(); it != m_iceMove.end(); ++it)
	{
		if(*it) delete *it;
	}
	
	sm_clusterNum = 0;	//クラスタ数の初期化

	Ice_SM::SetIterationNum(itr);
	Ice_SM::SetPrtPointerPosAndVel(s_sphPrtPos, s_sphPrtVel);

	cout << __FUNCTION__ << ", check1" << endl;

	//四面体リストを元に，粒子毎にクラスタ作成
	for(int i = 0; i < sm_particleNum; ++i)
	{
		//クラスタ初期化
		m_iceMove.push_back(new Ice_SM(sm_clusterNum));
		m_iceMove[sm_clusterNum]->SetSimulationSpace(boundarySpaceLow, boundarySpaceHigh);
		m_iceMove[sm_clusterNum]->SetTimeStep(timeStep);
		m_iceMove[sm_clusterNum]->SetCollisionFunc(0);
		m_iceMove[sm_clusterNum]->SetStiffness(1.0, 1.0);

		//四面体リストを元に，クラスタへ粒子を登録
		SetClusterMoveInfo(i);

		sm_clusterNum++;
	}

	cout << __FUNCTION__ << ", check2" << endl;

	Ice_SM::InitFinalParamPointer(sm_clusterNum);

	////MakeClusterFromNeight();
	////MakeOneCluster();
	
	//TODO::浮力を生むために粒子質量を下げる

//デバッグ
	//クラスタに含まれる粒子
	//for(int i = 0; i < ICENUM; i++)
	//{
	//	cout << "pIndx = " << endl;

	//	for(int j = 0; j < m_sm_cluster[i]->GetNumVertices(); j++)
	//	{
	//		cout << " j = " << m_sm_cluster[i]->GetParticleIndx(j);
	//	}
	//	cout << endl;
	//}

	//DebugClusterInfo();
}

//オブジェクトの構造の初期化
void IceObject::InitStrct()
{	cout << __FUNCTION__ << endl;

	//カウント
	for(int i = 0; i < sm_particleNum; i++)
	{
		int pNum = m_iceMove[i]->GetNumVertices();
		vector<int> pList;

		for(int j = 0; j < pNum; j++)
		{
			pList.push_back(m_iceMove[i]->GetParticleIndx(j));
		}
		m_iceStrct->CountClusterParticle(i, pList, pNum);
	}
	cout << __FUNCTION__ << " カウント" << endl;

	//メモリ確保
	InitClusterInfo();
	cout << __FUNCTION__ << " メモリ確保" << endl;

	//粒子が所属しているクラスタ数の配列をコピー
	int *PtoCNum = new int[sm_particleNum];
	cout << __FUNCTION__ << " 配列をコピー　メモリ確保" << endl;

	for(int i = 0; i < sm_particleNum; i++)
	{
		PtoCNum[i] = GetPtoCNum(i);
	}
	cout << __FUNCTION__ << " 配列をコピー" << endl;

	//クラスタと粒子の関連情報の登録
	for(int i = 0; i < sm_particleNum; i++)
	{
		SetClusterStrctInfo(i, PtoCNum);	//カウントを前で行っているため，こちらを使う
	}
	cout << __FUNCTION__ << " クラスタと粒子の関連情報の登録登録" << endl;

	delete[] PtoCNum;

//デバッグ
	//for(int i = 0; i < m_iClusteresNum; i++)
	//{
	//	m_ice->DebugCtoP(i);
	//}

	//for(int i = 0; i < m_pPS->GetNumParticles(); i++)
	//{
	//	m_ice->DebugPtoC(i);
	//}
}

//熱処理初期化
void IceObject::InitHeatTransfer(float effectiveRadius, float timeStep, float tempMax, float tempMin, float latentHeat, float cffcntHt, float cffcntTd)
{//	cout << __FUNCTION__ << endl;
	m_heatTransfer = new HeatTransfar(sm_particleNum*2);		//最初に最大数を確保しておいて，使うのは作成されたパーティクルまでとする
	m_heatTransfer->setCarnelConstant(effectiveRadius);			//カーネル関数の定数のための処理
	m_heatTransfer->setNumVertices(sm_particleNum);				//パーティクルの数を取得

	m_heatTransfer->setTimeStep(timeStep);
	m_heatTransfer->setTempMax(tempMax);
	m_heatTransfer->setTempMin(tempMin);
	m_heatTransfer->setLatentHeat(latentHeat);
	m_heatTransfer->setCffCntHt(cffcntHt);
	m_heatTransfer->setCffCntTd(cffcntTd);
}

//高速化に用いるパスの初期化
void IceObject::InitPath()
{
	m_SurfSm.InitPath(s_sphPrtPos, s_sphPrtVel, m_iceMove, m_iceStrct, sm_clusterNum, sm_particleNum);
}

//クラスタの運動計算情報の登録
void IceObject::SetClusterMoveInfo(int pIndx)
{
	//ある粒子が含まれる四面体を近傍四面体とし，各近傍四面体に含まれる粒子でクラスタを作成
	for(int i = 0; i < GetPtoTIndx(pIndx); i++)
	{
		int itIndx = GetPtoT(pIndx, i, 0);
		int ioIndx = GetPtoT(pIndx, i, 1);

		if(itIndx == -1 || ioIndx == -1){ continue;	}

		//四面体の全ての粒子をクラスタに登録
		for(int j = 0; j < GetTtoPIndx(itIndx); j++)
		{
			int jpIndx = GetTtoP(itIndx, j);

			if(jpIndx == -1){	continue;	}
			if(m_iceMove[pIndx]->CheckIndx(jpIndx)){	continue;	}

			int pNum = m_iceMove[pIndx]->GetNumVertices();
			float mass = 1.0f;

			m_iceMove[pIndx]->AddVertex( Vec3(s_sphPrtPos[jpIndx*4+0], s_sphPrtPos[jpIndx*4+1], s_sphPrtPos[jpIndx*4+2] ), mass, jpIndx);
			m_iceMove[pIndx]->SetAlphas(pNum, 1.0);
			m_iceMove[pIndx]->SetBetas (pNum, 0.0);
			m_iceMove[pIndx]->SetLayer (pNum, 0);
		}
	}

	//近傍四面体のlayerたどり，粒子を追加していく
	//TODO::不安定になるなら，layerが遠いほどbetaを下げる
	for(int i = 0; i < GetPtoTIndx(pIndx); i++)
	{
		int itIndx = GetPtoT(pIndx, i, 0);
		int ioIndx = GetPtoT(pIndx, i, 1);

		if(itIndx == -1 || ioIndx == -1){	continue;	}

		for(int j = 0; j < GetNTNum(itIndx); j++)
		{
			int jtIndx = GetNeighborTetra(itIndx, j, 0);
			int jlIndx = GetNeighborTetra(itIndx, j, 1);

			if(jtIndx == -1 || jlIndx == -1){	continue;	}

			//四面体の全ての粒子をクラスタに登録
			for(int k = 0; k < GetTtoPIndx(jtIndx); k++)
			{
				int kpIndx = GetTtoP(jtIndx, k);

				if(kpIndx == -1){	continue;	}
				if(m_iceMove[pIndx]->CheckIndx(kpIndx)){	/*cout << "contineue kpIndx = " << kpIndx << endl;*/ continue;	}

				int pNum = m_iceMove[pIndx]->GetNumVertices();
				float mass = 1.0f;

				m_iceMove[pIndx]->AddVertex( Vec3(s_sphPrtPos[kpIndx*4+0], s_sphPrtPos[kpIndx*4+1], s_sphPrtPos[kpIndx*4+2] ), mass, kpIndx);
				m_iceMove[pIndx]->SetAlphas(pNum, 1.0);
				m_iceMove[pIndx]->SetBetas (pNum, 0.0);
				m_iceMove[pIndx]->SetLayer (pNum, jlIndx);
				//cout << "pIndx = " << pIndx << " GetNumVertices = " << m_iceMove[pIndx]->GetNumVertices() << endl;
			}
		}
	}
}

//粒子→クラスタやクラスタ→粒子を登録
void IceObject::SetClusterStrctInfo(int cIndx, int* PtoCNum)
{
	//粒子が属している四面体の番号を登録するための準備
	//pCountListには，cIndx番目のクラスタに含まれる各粒子が，それぞれいくつのクラスタに属するかを求めて保存する
	int* pCountList = new int[m_iceMove[cIndx]->GetNumVertices()];
	int* pLayerList = new int[m_iceMove[cIndx]->GetNumVertices()];			//粒子の所属レイヤー
	
	for(int i = 0; i < m_iceMove[cIndx]->GetNumVertices(); i++)
	{
		int pIndx = m_iceMove[cIndx]->GetParticleIndx(i);

		//TODO::穴あきを想定していないため，誤って上書きしてしまっている
		//TODO::-1を探索して上書きするのに切り替える．
		for(int j = 0; j < GetPtoCMax(); j++)
		{
			//配列上で空いているスペースを探索（-1でない場所を探す）
			if(GetPtoC(pIndx, j, 0) != -1 || GetPtoC(pIndx, j, 1) != -1){	continue;	}
			
			if(j >= GetPtoCIndx(pIndx))
			{
				SetPtoCIndx(pIndx, j+1);		//現在の末尾Indxより大きいなら更新
			}

			pCountList[i] = j;					//配列上での場所を保存
			break;
		}

		pLayerList[i] = m_iceMove[cIndx]->GetLayer(i);					//粒子が何層目の近傍なのかを取得
	}

	//粒子とクラスタの情報登録
	vector<int> pIndxList;

	for(int i = 0; i < GetCtoPNum(cIndx); i++)
	{
		int pIndx = m_iceMove[cIndx]->GetParticleIndx(i);
			
		SetPtoC(pIndx, pCountList[i], cIndx, i, pLayerList[i]);		//粒子が所属しているクラスタを登録

		pIndxList.push_back(pIndx);									//クラスタに含まれる粒子集合
	}

	SetCtoP(cIndx, pIndxList, pLayerList);
	
	delete[] pCountList;
	delete[] pLayerList;
}

//補間処理のためのパラメータの初期化
void IceObject::InitInterPolation()
{	cout << __FUNCTION__ << endl;
	m_fInterPolationCoefficience = new float[sm_particleNum];	//線形補間係数

	for(int i = 0; i < sm_particleNum; ++i)
	{
		m_fInterPolationCoefficience[i] = 1.0f;
	}
}

//GPU処理で用いる変数の初期化
void IceObject::InitGPU()
{
	//固体構造の初期化
	m_iceStrct->InitGPU();

	//クラスタに関するGPUの初期化
	//すべてのクラスタにアクセスできるようにポインタを渡している
	Ice_SM::InitGPU(m_iceMove, sd_sphPrtPos, sd_sphPrtVel, sm_particleNum, sm_maxParticleNum);

#if defined(USE_PATH)
{
	m_SurfSm.InitPathGPU();		//高速化に用いるパスの初期化
}
#endif
	InitIceObjGPU();
}

void IceObject::InitIceObjGPU()
{
	int MAXCLUSTER = GetMaxClusterNum();

	//各粒子の最終的な位置・速度データ
	cudaMalloc((void**)&sd_sldPrtPos,	sizeof(float) * MAXCLUSTER * SM_DIM);
	cudaMalloc((void**)&sd_sldPrtVel,	sizeof(float) * MAXCLUSTER * SM_DIM);

	//最終位置・速度を現在のデータで初期化
	float* fPoses = new float[MAXCLUSTER * SM_DIM];
	float* fVeles = new float[MAXCLUSTER * SM_DIM];

	//s_pfPrtPosなどはデータの中身がDIM=4で作られているので，こうしないといけない
	//TODO::粒子サイズが大きくなると，メモリが確保できないかもしれないのに注意
	int sphDIM = 4;
	for(int i = 0; i < MAXCLUSTER; ++i)
	{
		for(int j = 0; j < SM_DIM; ++j)
		{
			fPoses[i*SM_DIM+j] = s_sphPrtPos[i*sphDIM+j];
			fVeles[i*SM_DIM+j] = s_sphPrtVel[i*sphDIM+j];
		}
	}

	cudaMemcpy(sd_sldPrtPos, fPoses, sizeof(float) * MAXCLUSTER * SM_DIM, cudaMemcpyHostToDevice);
	cudaMemcpy(sd_sldPrtVel, fVeles, sizeof(float) * MAXCLUSTER * SM_DIM, cudaMemcpyHostToDevice);

	////初期化の転送がうまくいったかの確認
	////一時配列のリセット
	//for(int i = 0; i < MAXCLUSTER; ++i)
	//{
	//	for(int j = 0; j < SM_DIM; ++j)
	//	{
	//		fPoses[i*SM_DIM+j] = 0.0f;
	//		fVeles[i*SM_DIM+j] = 0.0f;
	//	}
	//}

	////データを転送
	//cudaMemcpy(fPoses, d_FinalPos, sizeof(float) * MAXCLUSTER * SM_DIM, cudaMemcpyDeviceToHost);
	//cudaMemcpy(fVeles, d_FinalVel, sizeof(float) * MAXCLUSTER * SM_DIM, cudaMemcpyDeviceToHost);

	////ホスト側のデータを転送した結果をダンプ
	//ofstream ofs( "DtoH_Test.txt" );
	//ofs << "DtoH_Test" << endl;
	//
	////デバイス側のデータを転送
	//for(int i = 0; i < MAXCLUSTER; i++)
	//{
	//	ofs << "particle" << i << " pos::(" << fPoses[i*SM_DIM+0] << ", " << fPoses[i*SM_DIM+1] << ", " << fPoses[i*SM_DIM+2] << ")" << endl;
	//	ofs << "sph     " << i << " pos::(" << s_pfPrtPos[i*sphDIM+0] << ", " << s_pfPrtPos[i*sphDIM+1] << ", " << s_pfPrtPos[i*sphDIM+2] << ")" << endl;
	//	ofs << "particle" << i << " vel::(" << fVeles[i*SM_DIM+0] << ", " << fVeles[i*SM_DIM+1] << ", " << fVeles[i*SM_DIM+2] << ")" << endl;
	//	ofs << "sph     " << i << " vel::(" << s_pfPrtVel[i*sphDIM+0] << ", " << s_pfPrtVel[i*sphDIM+1] << ", " << s_pfPrtVel[i*sphDIM+2] << ")" << endl;
	//}

	delete[] fPoses;
	delete[] fVeles;
}

//固体の運動計算
void IceObject::StepObjMoveCPU()
{
	#pragma omp parallel
	{
	#pragma omp for
		for(int i = 0; i < sm_clusterNum; ++i)
		{	
			if(GetPtoCNum(i) == 0){	continue;	}
			m_iceMove[i]->UpdateCPU();				//運動計算
		}
	}//#pragma omp parallel
}

void IceObject::StepObjMoveSelected()
{
	#pragma omp parallel
	{
	#pragma omp for
		for(int i = 0; i < sm_clusterNum; ++i)
		{	
			if(GetPtoCNum(i) == 0){	continue;	}

			if(i%SELECTED != 1){	continue;	}	//テスト
			m_iceMove[i]->calExternalForces();					//現在位置の計算
			m_iceMove[i]->ShapeMatchingSelected(SELECTED);		//現在位置の更新 普通の計算
			m_iceMove[i]->integrate(0.02f);						//速度の計算
		}
	}//#pragma omp parallel
}

//GPU
void IceObject::StepObjMoveGPU()
{
	Ice_SM::SetDevicePosPointer(sd_sphPrtPos);	//VBOを使っているので，これがないとエラーが出るのに注意
	Ice_SM::UpdateGPU();
}

//パスを用いた高速化　固体の運動計算
void IceObject::StepObjMoveUsePath()
{
	//prefixSumの更新
	m_SurfSm.UpdatePrefixSum();

	//クラスタのパラメータ更新
	#pragma omp parallel for
	for(int i = 0; i < sm_clusterNum; ++i)
	{	
		if(GetPtoCNum(i) == 0){	continue;	}
		m_iceMove[i]->SetNowCm(m_SurfSm.CalcCmSum(i));		//重心の更新
		m_iceMove[i]->SetApq(m_SurfSm.CalcApqSum(i));		//変形行列の更新
	}

	//運動計算
	#pragma omp parallel for
	for(int i = 0; i < sm_clusterNum; ++i)
	{	
		if(GetPtoCNum(i) == 0){	continue; }
		m_iceMove[i]->UpdateUsePathCPU();
	}
}

void IceObject::StepObjMoveGPUUsePath()
{
	//prefixSumの更新
	m_SurfSm.SetDevicePointer(sd_sphPrtPos, sd_sphPrtVel);	//VBOを使っているので，これがないとエラーが出るのに注意
	m_SurfSm.UpdatePrefixSumGPU();
	RXTIMER("UpdatePrefixSumGPU");

	//クラスタのパラメータ更新
	LauchUpdateSMFromPath
	(
		sm_clusterNum,
		sd_sphPrtPos,
		sd_sphPrtVel,
		//SM-------------------------------------------------
		Ice_SM::GetOrgPosPointer(),
		Ice_SM::GetDevicePosPointer(),
		Ice_SM::GetOrgCmPointer(),
		Ice_SM::GetCurCmPointer(),
		Ice_SM::GetDeviceApqPointer(),
		Ice_SM::GetDeviceVelPointer(),
		Ice_SM::GetDeviceIndexesPointer(),
		Ice_SM::GetDeviceIndexSetPointer(),
		//Path-----------------------------------------------
		m_SurfSm.GetDevicePRTtoPTHPointer(),
		m_SurfSm.GetDevicePTHandPrfxSetPointer(),
		m_SurfSm.GetDecvicePrfxPos(),
		m_SurfSm.GetDecvicePrfxApq(),
		//IceStruct------------------------------------------
		m_iceStrct->GetDeviceCtoPointer(),
		m_iceStrct->GetDeviceCtoPNumPointer(),
		m_iceStrct->GetCtoPMax(),
		2,
		0.01		//TODO: 0.02では？いろいろばらばらっぽい
	);

	RXTIMER("LauchUpdateSMFromPath");

	//運動計算
	Ice_SM::SetDevicePosPointer(sd_sphPrtPos);				//VBOを使っているので，これがないとエラーが出るのに注意
	Ice_SM::UpdateUsePathGPU();
	RXTIMER("UpdateUsePathGPU");

//デバッグ
	//DebugObjMoveUsePathWithGPU();
}

//反復を用いた固体の運動計算
void IceObject::StepObjMoveIteration()
{
	//運動計算を反復
	for(int itr = 0; itr < Ice_SM::GetIteration(); ++itr)
	{
		//初回
		if(itr == 0)
		{
			#pragma omp parallel for
			for(int i = 0; i < sm_clusterNum; ++i)
			{	
				if(GetPtoCNum(i) == 0){	continue;	}
				GetMoveObj(i)->calExternalForces();				//速度を反映して位置更新
				GetMoveObj(i)->ShapeMatchingSolid();			//速度を反映した状態でSM法
			}
		}
		//反復
		else
		{
			#pragma omp parallel for
			for(int i = 0; i < sm_clusterNum; ++i)
			{	
				if(GetPtoCNum(i) == 0){	continue;	}

				GetMoveObj(i)->ShapeMatchingIteration();		//現在の粒子位置を用いてSM法
			}
		}

		StepInterPolationForCluster();							//各クラスタの総和計算，線形補間をして固体の最終位置を決定
	}

	//速度算出
	#pragma omp parallel for
	for(int i = 0; i < sm_clusterNum; ++i)
	{	
		if(GetPtoCNum(i) == 0){	continue;	}

		GetMoveObj(i)->integrateIteration();
	}
}

//パスを用いた反復運動計算
void IceObject::StepObjMoveIterationUsePath()
{
	//運動計算を反復
	for(int itr = 0; itr < Ice_SM::GetIteration(); ++itr)
	{
		//初回
		if(itr == 0)
		{
			m_SurfSm.UpdatePrefixSum();								//prefixSumの更新
		}
		//反復
		else
		{
			m_SurfSm.UpdatePrefixSumItr();							//prefixSumの更新
		}

		//各クラスタのデータ更新
		#pragma omp parallel for
		for(int i = 0; i < sm_clusterNum; ++i)
		{	
			if(GetPtoCNum(i) == 0){	continue;	}
			GetMoveObj(i)->SetNowCm(m_SurfSm.CalcCmSum(i));		//重心の更新
			GetMoveObj(i)->SetApq(m_SurfSm.CalcApqSum(i));		//変形行列の更新
		}

		//運動計算
		#pragma omp parallel for
		for(int i = 0; i < sm_clusterNum; ++i)
		{	
			if(GetPtoCNum(i) == 0){	continue;	}
			GetMoveObj(i)->ShapeMatchingUsePath();				//現在の位置でSM法
		}

		StepInterPolationForCluster();							//各クラスタの総和計算，線形補間をして固体の最終位置を決定
	}

	//速度算出
	#pragma omp parallel for
	for(int i = 0; i < sm_clusterNum; ++i)
	{	
		if(GetPtoCNum(i) == 0){	continue;	}
		GetMoveObj(i)->integrateIteration();
	}
}

//固体の運動計算，総和計算，補間処理　GPU処理　反復処理あり
void IceObject::StepObjCalcWidhIteration()
{
	for(int itr = 0; itr < g_iterationNum; ++itr)
	{
		if(itr == 0)
		{
			Ice_SM::UpdateGPU();
		}
		else
		{
			Ice_SM::UpdateIterationGPU(sd_sldPrtPos, sd_sldPrtVel);
		}

		StepInterPolationForCluster();			//総和計算，線形補間
	}

	StepInterPolation();
}

//固体の運動計算　重み付け＋反復
void IceObject::StepObjMoveIterationWeighted()
{
	//運動計算を反復
	for(int itr = 0; itr < Ice_SM::GetIteration(); ++itr)
	{
		//初回
		if(itr == 0)
		{
			#pragma omp parallel for
			for(int i = 0; i < sm_clusterNum; ++i)
			{	
				if(GetPtoCNum(i) == 0){	continue;	}
				GetMoveObj(i)->calExternalForces();				//速度を反映して位置更新
				GetMoveObj(i)->ShapeMatchingSolid();			//速度を反映した状態でSM法
			}
		}
		//反復
		else
		{
			#pragma omp parallel for
			for(int i = 0; i < sm_clusterNum; ++i)
			{	
				if(GetPtoCNum(i) == 0){	continue;	}
				GetMoveObj(i)->ShapeMatchingIteration();		//現在の粒子位置を用いてSM法
			}
		}

		StepInterPolationForClusterWeighted();					//各クラスタの総和計算，線形補間をして固体の最終位置を決定
	}

	//速度算出
	#pragma omp parallel for
	for(int i = 0; i < sm_clusterNum; ++i)
	{	
		if(GetPtoCNum(i) == 0){	continue;	}
		GetMoveObj(i)->integrateIteration();
	}
}





//流体と固体の最終的な運動計算
void IceObject::StepInterPolation()
{
#if defined(MK_USE_GPU)

	//固体運動の最終位置計算
	float* smPrtPos = Ice_SM::GetDevicePosPointer();
	float* smPrtVel = Ice_SM::GetDeviceVelPointer();
	int* indxSet = Ice_SM::GetDeviceIndexSetPointer();

	sd_sphPrtPos = Ice_SM::GetDeviceSPHPosPointer();
	sd_sphPrtVel = Ice_SM::GetDeviceSPHVelPointer();

	int* PtoCIndx = IceStructure::GetDevicePtoCIndxPointer();
	int* PtoC = IceStructure::GetDevicePtoCPointer();
	int PNumMax = m_iceStrct->GetPNumMax();
	int PtoCMax = m_iceStrct->GetPtoCMax();
	int PtoCParamSize = 3;

	//固体運動の最終位置計算
	LaunchCalcAverageGPU(sm_particleNum, sd_sldPrtPos, sd_sldPrtVel, sd_sphPrtPos, sd_sphPrtVel, smPrtPos, smPrtVel, indxSet, PtoCIndx, PtoC, PNumMax, PtoCMax, PtoCParamSize);

	//液体と固体の補間
	LaunchInterPolationGPU(sm_particleNum, sd_sldPrtPos, sd_sldPrtVel, sd_sphPrtPos, sd_sphPrtVel);

#else

	Vec3 pos,vel;

	#pragma omp parallel for private(pos, vel)
	for(int i = 0; i < sm_particleNum; ++i)
	{
		if(GetParticleNum() <= i){	continue;	}	//融解のみの実験のときに必要になる．
		if(GetPtoCNum(i) <= 0){		continue;	}
	
		//固体運動の最終位置計算
		CalcAverageCPU(i, pos, vel);

		//液体と固体の補間
		LinerInterPolationCPU(i, pos, vel);			//線形補間
	}

//デバッグ
	//DebugDeformationAmount();
	DebugDeformationAverage();

#endif
}

//粒子を選択的に運動計算させるパターンのテスト
void IceObject::StepInterPolationSelected()
{
	Vec3 pos,vel;

	vector<Vec3> preSphPos;
	preSphPos.resize(sm_particleNum);
	for(int i = 0; i < sm_particleNum; i++)
	{
		Vec3 pos;
		pos[0] = s_sphPrtPos[i*4+0];
		pos[1] = s_sphPrtPos[i*4+1];
		pos[2] = s_sphPrtPos[i*4+2];

		preSphPos.push_back(pos);
	}

	#pragma omp parallel for private(pos, vel)
	for(int i = 0; i < sm_particleNum; ++i)
	{
		if(GetParticleNum() <= i){	continue;	}	//融解のみの実験のときに必要になる．
		if(GetPtoCNum(i) <= 0){		continue;	}
	
		//固体運動の最終位置計算
		CalcAverageSelected(i, pos, vel);

		//液体と固体の補間
		LinerInterPolationCPU(i, pos, vel);			//線形補間
	}

	#pragma omp parallel for private(pos, vel)
	for(int i = 0; i < sm_particleNum; ++i)
	{
		if(GetParticleNum() <= i){	continue;	}	//融解のみの実験のときに必要になる．
		if(GetPtoCNum(i) <= 0){		continue;	}

		//変形量の反映
		UpdateUnSelectedCluster(i, pos, vel, preSphPos);
	}


//デバッグ
	//DebugDeformationAmount();
	DebugDeformationAverage();
}

void IceObject::UpdateUnSelectedCluster(int cIndx, const Vec3& pos, const Vec3& vel, const vector<Vec3>& preSphPos)
{
	if(cIndx%SELECTED == 1)	return;	//運動計算していないクラスタが対象

	float defAmount = 0.0f;

	for(int i = 0; i < m_iceMove[cIndx]->GetVertexNum(); ++i)
	{
		int pIndx = m_iceMove[cIndx]->GetParticleIndx(i);
		
		if(pIndx%SELECTED == 1)	continue;

		Vec3 nowPos = m_iceMove[cIndx]->GetVertexPos(i);
		Vec3 prePos = preSphPos[pIndx];

		defAmount += abs(nowPos[0]-prePos[0]);
		defAmount += abs(nowPos[1]-prePos[1]);
		defAmount += abs(nowPos[2]-prePos[2]);
	}

	m_iceMove[cIndx]->SetDefAmount(defAmount);
}

void IceObject::StepWeightedInterPolation()
{
	Vec3 pos,vel;

	#pragma omp parallel for private(pos, vel)
	for(int i = 0; i < sm_particleNum; ++i)
	{
		if(GetParticleNum() <= i){	continue;	}	//融解のみの実験のときに必要になる．
		if(GetPtoCNum(i) <= 0){		continue;	}
	
		//固体運動の最終位置計算
		CalcWeightedVector(i, pos, vel);

		//液体と固体の補間
		LinerInterPolationCPU(i, pos, vel);			//線形補間
	}

//デバッグ
	//DebugDeformationAmount();
}

void IceObject::StepInterPolationForClusterWeighted()
{
	Vec3 pos,vel;

	#pragma omp parallel for private(pos, vel)
	for(int i = 0; i < sm_particleNum; ++i)
	{
		if(GetParticleNum() <= i){	continue;	}		//融解のみの実験のときに必要になる．
		if(GetPtoCNum(i) <= 0){		continue;	}

		//固体運動の最終位置計算
		//CalcAverageCPU(i, pos, vel);
		CalcWeightedVector(i, pos, vel);

		//液体と固体の補間
		LinerInterPolationForClusterCPU(i, pos, vel);	//線形補間
	}
}

//反復処理を用いた流体と固体の最終的な運動計算
void IceObject::StepInterPolationForCluster()
{
#if defined(MK_USE_GPU)
	//固体運動の最終位置計算
	float* smPrtPos = Ice_SM::GetDevicePosPointer();
	float* smPrtVel = Ice_SM::GetDeviceVelPointer();
	int* indxSet = Ice_SM::GetDeviceIndexSetPointer();

	sd_sphPrtPos = Ice_SM::GetDeviceSPHPosPointer();
	sd_sphPrtVel = Ice_SM::GetDeviceSPHVelPointer();

	int* PtoCIndx = IceStructure::GetDevicePtoCIndxPointer();
	int* PtoC = IceStructure::GetDevicePtoCPointer();
	int PNumMax = m_iceStrct->GetPNumMax();
	int PtoCMax = m_iceStrct->GetPtoCMax();
	int PtoCParamSize = 3;

	//固体運動の最終位置計算
	LaunchCalcAverageGPU(sm_particleNum, sd_sldPrtPos, sd_sldPrtVel, sd_sphPrtPos, sd_sphPrtVel, smPrtPos, smPrtVel, indxSet, PtoCIndx, PtoC, PNumMax, PtoCMax, PtoCParamSize);

#else
	Vec3 pos,vel;

	#pragma omp parallel for private(pos, vel)
	for(int i = 0; i < sm_particleNum; ++i)
	{
		if(GetParticleNum() <= i){	continue;	}		//融解のみの実験のときに必要になる．
		if(GetPtoCNum(i) <= 0){		continue;	}

		//固体運動の最終位置計算
		CalcAverageCPU(i, pos, vel);

		//液体と固体の補間
		LinerInterPolationForClusterCPU(i, pos, vel);	//線形補間
	}

#endif
}

//各クラスタの計算結果の平均を，固体の最終的な運動計算結果とする
void IceObject::CalcAverageCPU(int pIndx, Vec3& pos, Vec3& vel)
{
	//それぞれのベクトルを合成し平均をとる
	pos = Vec3(0.0, 0.0, 0.0);
	vel = Vec3(0.0, 0.0, 0.0);
	double shapeNum = 0.0;		//クラスタの数
	int cIndx, oIndx;

	for(int i = 0; i < GetPtoCIndx(pIndx); ++i)
	{
		cIndx = GetPtoC(pIndx, i, 0);
		oIndx = GetPtoC(pIndx, i, 1);

		if(cIndx == -1 || oIndx == -1){	continue;	}

		pos += m_iceMove[cIndx]->GetVertexPos(oIndx);
		vel += m_iceMove[cIndx]->GetVertexVel(oIndx);

		shapeNum += 1.0;
	}

	//クラスタの数で割る
	if(shapeNum != 0.0)
	{
		pos /= shapeNum;
		vel /= shapeNum;
	}		
	//どのクラスタにも含まれていない場合，運動はSPH法に従う
	else
	{
		int jpIndx = pIndx*4;
		pos = Vec3(s_sphPrtPos[jpIndx+0], s_sphPrtPos[jpIndx+1], s_sphPrtPos[jpIndx+2]);
		vel = Vec3(s_sphPrtVel[jpIndx+0], s_sphPrtVel[jpIndx+1], s_sphPrtVel[jpIndx+2]);
	}
}

void IceObject::CalcAverageSelected(int pIndx, Vec3& pos, Vec3& vel)
{
	//それぞれのベクトルを合成し平均をとる
	pos = Vec3(0.0, 0.0, 0.0);
	vel = Vec3(0.0, 0.0, 0.0);
	double shapeNum = 0.0;		//クラスタの数
	int cIndx, oIndx;

	for(int i = 0; i < GetPtoCIndx(pIndx); ++i)
	{
		cIndx = GetPtoC(pIndx, i, 0);
		oIndx = GetPtoC(pIndx, i, 1);

		if(cIndx == -1 || oIndx == -1){	continue;	}

		//補間粒子クラスタに含まれるデータを使わず，計算粒子クラスタに含まれるデータのみを使う
		if(cIndx%SELECTED != 1){		continue;	}

		pos += m_iceMove[cIndx]->GetVertexPos(oIndx);
		vel += m_iceMove[cIndx]->GetVertexVel(oIndx);

		shapeNum += 1.0;
	}

	//クラスタの数で割る
	if(shapeNum != 0.0)
	{
		pos /= shapeNum;
		vel /= shapeNum;
	}
	//どのクラスタにも含まれていない場合，運動はSPH法に従う
	else
	{
		int jpIndx = pIndx*4;
		pos = Vec3(s_sphPrtPos[jpIndx+0], s_sphPrtPos[jpIndx+1], s_sphPrtPos[jpIndx+2]);
		vel = Vec3(s_sphPrtVel[jpIndx+0], s_sphPrtVel[jpIndx+1], s_sphPrtVel[jpIndx+2]);
	}
}

void IceObject::CalcWeightedVector(int pIndx, Vec3& pos, Vec3& vel)
{
	//それぞれのベクトルを合成し平均をとる
	pos = Vec3(0.0, 0.0, 0.0);
	vel = Vec3(0.0, 0.0, 0.0);

	int cIndx, oIndx;
	int clusterNum = GetPtoCIndx(pIndx);	//テストということで
	int average = 1.0f/(float)clusterNum;

	//位置，速度ベクトル算出
	float deformAmountSum = 0.0f;

	for(int i = 0; i < GetPtoCIndx(pIndx); ++i)
	{
		cIndx = GetPtoC(pIndx, i, 0);
		oIndx = GetPtoC(pIndx, i, 1);

		if(cIndx == -1 || oIndx == -1){	continue;	}

		float defAmount = m_iceMove[cIndx]->GetDefAmount();
		//defAmount *= defAmount;				//二乗版
		//defAmount *= defAmount * defAmount;	//三乗版

		pos += m_iceMove[cIndx]->GetVertexPos(oIndx) * defAmount;
		vel += m_iceMove[cIndx]->GetVertexVel(oIndx) * defAmount;

		deformAmountSum += defAmount;
	}

	pos /= deformAmountSum;
	vel /= deformAmountSum;
	

	//for(int i = 0; i < GetPtoCIndx(pIndx); ++i)
	//{
	//	cIndx = GetPtoC(pIndx, i, 0);
	//	oIndx = GetPtoC(pIndx, i, 1);

	//	if(cIndx == -1 || oIndx == -1){	continue;	}

	//	float defAmount = m_iceMove[cIndx]->GetDefAmount()/deformAmountSum;
	//	//defAmount *= defAmount;		//二乗版

	//	m_iceMove[cIndx]->SetDefPriority(defAmount*10.0f/* - average*/);	//平均とどのくらい差があるかを表す
	//}

	//TODO::逆バージョン　うまく動かない　初期のdeformAmountSumが0だから
	//float deformAmountSum = 0.0f;
	//for(int i = 0; i < GetPtoCIndx(pIndx); ++i)
	//{
	//	cIndx = GetPtoC(pIndx, i, 0);
	//	oIndx = GetPtoC(pIndx, i, 1);

	//	if(cIndx == -1 || oIndx == -1){	continue;	}

	//	deformAmountSum += m_iceMove[cIndx]->GetDefAmount();
	//}

	////位置，速度ベクトル算出
	//for(int i = 0; i < GetPtoCIndx(pIndx); ++i)
	//{
	//	cIndx = GetPtoC(pIndx, i, 0);
	//	oIndx = GetPtoC(pIndx, i, 1);

	//	if(cIndx == -1 || oIndx == -1){	continue;	}

	//	float defAmount = m_iceMove[cIndx]->GetDefAmount();
	//	pos += m_iceMove[cIndx]->GetVertexPos(oIndx) * (deformAmountSum - defAmount)/deformAmountSum;
	//	vel += m_iceMove[cIndx]->GetVertexVel(oIndx) * (deformAmountSum - defAmount)/deformAmountSum;
	//}
}

//SPH法とSM法で求めた速度と位置を線形補間　CPU
void IceObject::LinerInterPolationCPU(int pIndx, const Vec3& pos, const Vec3& vel)
{
	int sphIndx = pIndx*4;
	double intrps = 1.0-m_fInterPolationCoefficience[pIndx];	//補間係数

	s_sphPrtVel[sphIndx+0] = vel[0] * m_fInterPolationCoefficience[pIndx] + s_sphPrtVel[sphIndx+0] * intrps;
	s_sphPrtVel[sphIndx+1] = vel[1] * m_fInterPolationCoefficience[pIndx] + s_sphPrtVel[sphIndx+1] * intrps;
	s_sphPrtVel[sphIndx+2] = vel[2] * m_fInterPolationCoefficience[pIndx] + s_sphPrtVel[sphIndx+2] * intrps;

	s_sphPrtPos[sphIndx+0] = pos[0] * m_fInterPolationCoefficience[pIndx] + s_sphPrtPos[sphIndx+0] * intrps;
	s_sphPrtPos[sphIndx+1] = pos[1] * m_fInterPolationCoefficience[pIndx] + s_sphPrtPos[sphIndx+1] * intrps;
	s_sphPrtPos[sphIndx+2] = pos[2] * m_fInterPolationCoefficience[pIndx] + s_sphPrtPos[sphIndx+2] * intrps;
}

//SPH法とSM法で求めた速度と位置を線形補間　反映先はクラスタ CPU
void IceObject::LinerInterPolationForClusterCPU(int pIndx, const Vec3& pos, const Vec3& vel)
{
	int sphIndx = pIndx*4;
	double intrps = 1.0-m_fInterPolationCoefficience[pIndx];	//補間係数
	float* p = Ice_SM::GetSldPosPointer();
	float* v = Ice_SM::GetSldVelPointer();

	v[pIndx*3+0] = vel[0] * m_fInterPolationCoefficience[pIndx] + s_sphPrtVel[sphIndx+0] * intrps;
	v[pIndx*3+1] = vel[1] * m_fInterPolationCoefficience[pIndx] + s_sphPrtVel[sphIndx+1] * intrps;
	v[pIndx*3+2] = vel[2] * m_fInterPolationCoefficience[pIndx] + s_sphPrtVel[sphIndx+2] * intrps;

	p[pIndx*3+0] = pos[0] * m_fInterPolationCoefficience[pIndx] + s_sphPrtPos[sphIndx+0] * intrps;
	p[pIndx*3+1] = pos[1] * m_fInterPolationCoefficience[pIndx] + s_sphPrtPos[sphIndx+1] * intrps;
	p[pIndx*3+2] = pos[2] * m_fInterPolationCoefficience[pIndx] + s_sphPrtPos[sphIndx+2] * intrps;

}




//熱処理
void IceObject::StepHeatTransfer(const int* surfParticles, const vector<vector<rxNeigh>>& neights, float floor, float effRadius, const float* pos, const float* dens)
{
	vector<int> ids;														//表面粒子の添え字
	vector<float> dists;													//表面粒子の距離

	//初期化
	m_heatTransfer->resetNeighborhoodsId();
	m_heatTransfer->resetNeighborhoodsDist();

	//近傍粒子の添え字と距離の設定
	for( int i = 0; i < sm_particleNum; i++ )
	{
		//表面粒子判定情報　１ならば表面粒子，０なら内部粒子　その際の数値は近傍粒子総数を表す
		//近傍粒子総数で表面積を近似
		if( surfParticles[i] == 1 )
		{
			//床に近く，圧力の高い粒子は，表面ではなく底面とする
			//粒子数によってパラメータを変えないといけない．
			//表面粒子判定に不具合があるので修正が必要．
			//0.75で下２段とれる
			if(pos[i*4+1] < floor+effRadius*0.2)				//1331　下１段
			{
				if(dens[i] < 950.0)
				{
					m_heatTransfer->setSurfaceParticleNums(i, (int)( neights[i].size() ));		//表面扱い
				}
				else
				{
					m_heatTransfer->setSurfaceParticleNums(i, -1);								//底面扱い
				}
			}
			else
			{
				m_heatTransfer->setSurfaceParticleNums(i, (int)( neights[i].size() ));
			}
		}
		else
		{
			m_heatTransfer->setSurfaceParticleNums(i, -1);
		}
		
		//初期化
		ids.clear();
		dists.clear();

		for( unsigned j = 0; j < neights[i].size(); j++)
		{
			if( i == (int)( neights[i][j].Idx ) ) continue;							//自分自身を省く
			ids.push_back( (int)( neights[i][j].Idx ) );
			dists.push_back( (float)( neights[i][j].Dist ) );
		}

		m_heatTransfer->AddNeighborhoodsId( ids );
		m_heatTransfer->AddNeighborhoodsDist( dists );
	}

	//熱処理計算
	m_heatTransfer->heatAirAndParticle(); 		 										//熱処理　空気と粒子
	m_heatTransfer->heatParticleAndParticle(dens, effRadius);	//熱処理　粒子間
	m_heatTransfer->calcTempAndHeat();													//熱量の温度変換，温度の熱量変換

//デバッグ
	//for(int i = 0; i < sm_particleNum; i++)
	//{
	//	cout << m_heatTransfer->getTemps()[i] << endl;
	//}
}

//伝熱処理，相変化処理や固体データの更新
void IceObject::StepIceStructure()
{
	StepMelting();			//融解
	//StepFreezing();		//凝固
}

void IceObject::StepMelting()
{
	vector<unsigned> prtList;			//融解した粒子の集合
	vector<unsigned> clusterList;		//再定義するクラスタの集合
	vector<unsigned> cLayerList;		//再定義するクラスタのレイヤー
	vector<unsigned> tetraList;			//再定義する四面体の集合
	vector<unsigned> tLayerList;		//再定義する四面体のレイヤー

	SearchMeltParticle(prtList);		//融解粒子探索

	m_iceStrct->StepObjMelt(prtList,clusterList, tetraList, cLayerList, tLayerList);	//融解処理
	
	ReConstructCluster(prtList, clusterList);	//TODO::クラスタの再構築
}

void IceObject::StepFreezing()
{
	vector<unsigned> vuFreezePrtList;			//凝固した粒子の集合

	SearchMeltParticle(vuFreezePrtList);		//凝固粒子探索
	m_iceStrct->StepObjFreeze();
}

void IceObject::SearchMeltParticle(vector<unsigned>& pList)
{
	for(int i = 0; i < sm_particleNum; ++i)
	{	
		if( m_heatTransfer->getPhaseChange(i) != 1 )	continue;	//相転移の条件を満たしている
		if( m_heatTransfer->getPhase(i) != 2 )			continue;	//水へと相転移している
		if( m_iceStrct->GetParticleNum() <= i)			continue;	//融解のみの実験のときに必要になる．
		if( m_iceStrct->GetPtoCNum(i) == 0 )			continue;	//クラスタに含まれていない
		if( m_iceStrct->GetPtoTNum(i) == 0 )			continue;	//四面体に含まれていない

		if(pList.size() > 50){	break;	}							//融解粒子数の制限

		m_fInterPolationCoefficience[i] = 0.0f;						//線形補間しない
		m_heatTransfer->setPhaseChange(i, 0);						//相転移し終わったことを伝える
		pList.push_back(i);											//融解粒子の記録
	}

//デバッグ
	//cout << __FUNCTION__ << endl;

	//for(unsigned i = 0; i < pList.size(); i++)
	//{
	//	cout << pList[i] << endl;
	//}

	//cout << "pList.size() = " << pList.size() << endl;
}

void IceObject::SearchFreezeParticle(vector<unsigned>& pList)
{
//	for(int i = 0; i < m_pPS->GetNumParticles(); i++)
//	{	
//		if(m_ht->getPhaseChange(i) != 1)				continue;	//相転移の条件を満たしていない場合は戻る
//		if(m_ht->getPhase(i) != -2)						continue;	//氷へと相転移していない場合は戻る
//		if(m_ice->GetParticleNum() <= i)				continue;	//融解のみの実験のときに必要になる．
//		if(m_ice->GetPtoCNum(i) != 0)					continue;	//クラスタに含まれている
//		if(m_ice->GetPtoTNum(i) != 0)					continue;	//クラスタに含まれている
////		if(pList.size() > 1){	break;	}							//凝固粒子数の制限
//		
//		pList.push_back(i);											//凝固粒子の記録
//	}
}

void IceObject::ReConstructCluster(vector<unsigned>& pList, vector<unsigned>& cList)
{
	int pListSize = pList.size();
	int cListSize = cList.size();

	if(pListSize == 0 || cListSize == 0){	return; }

	vector<int> checkTList;
	int j = 0, k = 0, l = 0;
	int icIndx = 0;
	int jtIndx = 0, joIndx = 0;
	int kpIndx = 0, ktIndx = 0, klIndx = 0;
	int lpIndx = 0;
	int pNum = 0;

	//再定義クラスタの初期化
	#pragma omp parallel
	{
	#pragma omp for
		for(int i = 0; i < cListSize; ++i)
		{
			if(std::find(pList.begin(), pList.end(), cList[i]) != pList.end()){	continue;	}
			GetMoveObj(cList[i])->Clear();
		}
	}//end #pragma omp parallel

	//クラスタの再定義
	#pragma omp parallel
	{
	#pragma omp for private(checkTList, j, k, l, icIndx, jtIndx, joIndx, kpIndx, ktIndx, klIndx, lpIndx, pNum)
		for(int i = 0; i < cListSize; ++i)
		{
			checkTList.clear();
			icIndx = cList[i];
			if(std::find(pList.begin(), pList.end(), icIndx) != pList.end()){	continue;	}

			//クラスタを再定義する際，基本となる粒子が属する四面体から初期粒子を得る．
			//クラスタ番号＝＝粒子番号なのに注意
			//以下を関数にすると，エラーが出てうまくいかない
			for(j = 0; j < GetPtoTIndx(icIndx); j++)
			{
				jtIndx = GetPtoT(icIndx, j, 0);
				joIndx = GetPtoT(icIndx, j, 1);

				if(jtIndx == -1 || joIndx == -1){ continue;	}
				if(std::find(checkTList.begin(), checkTList.end(), jtIndx) != checkTList.end())
				{	continue;	}
				else
				{	checkTList.push_back(jtIndx);	}

				//四面体の全ての粒子をクラスタに登録
				for(k = 0; k < GetTtoPIndx(jtIndx); k++)
				{
					kpIndx = GetTtoP(jtIndx, k);

					if(kpIndx == -1){	continue;	}
					if(GetMoveObj(icIndx)->CheckIndx(kpIndx)){	continue;	}

					pNum = GetMoveObj(icIndx)->GetNumVertices();

					GetMoveObj(icIndx)->AddVertex(Vec3(s_sphPrtPos[kpIndx*4+0], s_sphPrtPos[kpIndx*4+1], s_sphPrtPos[kpIndx*4+2]), 1.0, kpIndx);
					//GetMoveObj(icIndx)->SetVelocity(pNum, Vec3(s_sphPrtVel[kpIndx*4+0], s_sphPrtVel[kpIndx*4+1], s_sphPrtVel[kpIndx*4+2]));
					GetMoveObj(icIndx)->SetAlphas(pNum, 1.0);
					GetMoveObj(icIndx)->SetBetas (pNum, 0.0);
					GetMoveObj(icIndx)->SetLayer (pNum, 0);
				}
			}

			//近傍四面体のlayerたどり，粒子を追加していく
			//TODO::不安定になるなら，layerが遠いほどbetaを下げる
			for(j = 0; j < GetPtoTIndx(icIndx); ++j)
			{
				jtIndx = GetPtoT(icIndx, j, 0);
				joIndx = GetPtoT(icIndx, j, 1);

				if(jtIndx == -1 || joIndx == -1){ continue;	}

				for(k = 0; k < GetNTNum(jtIndx); k++)
				{
					ktIndx = GetNeighborTetra(jtIndx, k, 0);
					klIndx = GetNeighborTetra(jtIndx, k, 1);

					if(ktIndx == -1 || klIndx == -1){	continue;	}
					if(std::find(checkTList.begin(), checkTList.end(), ktIndx) != checkTList.end())
					{	continue;	}
					else
					{	checkTList.push_back(ktIndx);	}

					//四面体の全ての粒子をクラスタに登録
					for(l = 0; l < GetTtoPIndx(ktIndx); l++)
					{
						lpIndx = GetTtoP(ktIndx, l);

						if(lpIndx == -1){	continue;	}
						if(GetMoveObj(icIndx)->CheckIndx(lpIndx)){	/*cout << "contineue kpIndx = " << kpIndx << endl;*/ continue;	}

						pNum = GetMoveObj(icIndx)->GetNumVertices();

						GetMoveObj(icIndx)->AddVertex( Vec3(s_sphPrtPos[lpIndx*4+0], s_sphPrtPos[lpIndx*4+1], s_sphPrtPos[lpIndx*4+2] ), 1.0, lpIndx);
						//GetMoveObj(icIndx)->SetVelocity(pNum, Vec3(s_sphPrtVel[lpIndx*4+0], s_sphPrtVel[lpIndx*4+1], s_sphPrtVel[lpIndx*4+2]));
						GetMoveObj(icIndx)->SetAlphas(pNum, 1.0);
						GetMoveObj(icIndx)->SetBetas (pNum, 0.0);
						GetMoveObj(icIndx)->SetLayer (pNum, klIndx);
					}
				}
			}
		}
	}//end #pragma omp parallel
}



//---------------------------------------------デバッグ------------------------------

//四面体
void IceObject::DebugTetraInfo()
{	cout << __FUNCTION__ << endl;

	unsigned tetraSize = IceTetrahedra::GetInstance().GetTetraListSize();

	//四面体→粒子
	for(unsigned i = 0; i < tetraSize; i++ )
	{
		m_iceStrct->DebugTtoP(i);
	}

	//粒子→四面体
	for(int i = 0; i < sm_particleNum; i++)
	{
		m_iceStrct->DebugPtoT(i);
	}

	//近傍四面体
	for(unsigned i = 0; i < tetraSize; i++ )
	{
		m_iceStrct->DebugNeighborTetra(i);
	}
}

//クラスタ
void IceObject::DebugClusterInfo()
{	cout << __FUNCTION__ << endl;

	//含有粒子数の総数，平均，最小値，最大値
	int sum = 0;
	int max = m_iceMove[0]->GetNumVertices();
	int min = m_iceMove[0]->GetNumVertices();

	for(int i = 0; i < sm_clusterNum; i++)
	{
		int num =  m_iceMove[i]->GetNumVertices();
		
		if(max < num){	max = num;	}

		if(min > num){	min = num;	}

		sum += num;
	}

	cout << "sum = " << sum << " ave = " << (float)(sum/sm_clusterNum) << endl;
	cout << "max = " << max << " min = " << min << endl;
}

void IceObject::DebugDeformationAmount()
{	//cout << __FUNCTION__ << endl;
	
	int deformClusterNum = 0;
	for(int i = 0; i < sm_clusterNum; ++i)
	{	
		if(GetPtoCNum(i) == 0){	continue;	}

		float defAmount = m_iceMove[i]->GetDefAmount();
		if(defAmount < 0.5f){	continue;	}

		//cout << "i = " << i << "::" <<  defAmount << endl;

#ifdef USE_SELECTED
		deformClusterNum += SELECTED;
#else
		deformClusterNum++;
#endif
	}

	//cout << "Deformation Amount deformClusterNum::" << deformClusterNum << endl;

	ofstream ofs( "DeformationAmount.txt", ios::app);
	ofs << deformClusterNum << endl;
}

void IceObject::DebugDeformationAverage()
{
	float deformClusterAverage = 0.0f;
	for(int i = 0; i < sm_clusterNum; ++i)
	{	
		if(GetPtoCNum(i) == 0){	continue;	}

		float defAmount = m_iceMove[i]->GetDefAmount();

		deformClusterAverage += defAmount;
	}

	//cout << "Deformation Average deformClusterNum::" << deformClusterNum << endl;

	ofstream ofs( "DeformationAverage.txt", ios::app);
	ofs << deformClusterAverage/(float)sm_clusterNum << endl;
}

void IceObject::DebugObjMoveUsePathWithGPU()
{
//	QueryCounter counter;
//	QueryCounter counter2;	
//	QueryCounter counter3;
//
//counter.Start();
//
//	//prefixSumの更新
//	m_SurfSm.SetDevicePointer(sd_sphPrtPos, sd_sphPrtVel);	//VBOを使っているので，これがないとエラーが出るのに注意
//	m_SurfSm.UpdatePrefixSumGPU();
//
//	double prefixSumTime = counter.End()/1000.0;
////cout << "UpdatePrefixSumGPU：	" << prefixSumTime << endl;
//
//counter2.Start();
//
//	//クラスタのパラメータ更新
//	//TestUpdateSMFromPath();
//	LauchUpdateSMFromPath
//	(
//		sm_clusterNum,
//		sd_sphPrtPos,
//		sd_sphPrtVel,
//		//SM-------------------------------------------------
//		Ice_SM::GetOrgPosPointer(),
//		Ice_SM::GetDevicePosPointer(),
//		Ice_SM::GetOrgCmPointer(),
//		Ice_SM::GetCurCmPointer(),
//		Ice_SM::GetDeviceApqPointer(),
//		Ice_SM::GetDeviceVelPointer(),
//		Ice_SM::GetDeviceIndexesPointer(),
//		Ice_SM::GetDeviceIndexSetPointer(),
//		//Path-----------------------------------------------
//		m_SurfSm.GetDevicePRTtoPTHPointer(),
//		m_SurfSm.GetDevicePTHandPrfxSetPointer(),
//		m_SurfSm.GetDecvicePrfxPos(),
//		m_SurfSm.GetDecvicePrfxApq(),
//		//IceStruct------------------------------------------
//		m_iceStrct->GetDeviceCtoPointer(),
//		m_iceStrct->GetDeviceCtoPNumPointer(),
//		m_iceStrct->GetCtoPMax(),
//		2,
//		0.01		//TODO: 0.02では？いろいろばらばらっぽい
//	);
//
//	double updatePosApqTime = counter2.End()/1000.0;
////cout << "LauchUpdateSMFromPath：	" << updatePosApqTime << endl;
//
//counter3.Start();
//
//	//運動計算
//	Ice_SM::SetDevicePosPointer(sd_sphPrtPos);				//VBOを使っているので，これがないとエラーが出るのに注意
//	Ice_SM::UpdateUsePathGPU();
//
//	double updateSMTime = counter3.End()/1000.0;
////cout << "UpdateUsePathGPU：	" << updateSMTime << endl;
//
////cout << "合計：		" << prefixSumTime+updatePosApqTime+updateSMTime << endl;
//
//
////ホスト側のデータを転送した結果をダンプ
//	//古いファイルを削除
//	
//
//	//prefixSum
//	ofstream ofs( "LOG_prefixSumTime.txt", ios::app);
//	ofs << prefixSumTime << endl;
//
//	//updatePosApq
//	ofstream ofs1( "LOG_updateTime.txt", ios::app);
//	ofs1 << updatePosApqTime << endl;
//
//	//ShapeMatching
//	ofstream ofs2( "LOG_SMTime.txt", ios::app);
//	ofs2 << updateSMTime << endl;
}

//---------------------------------------------------テスト---------------------------------------------
void IceObject::TestUpdateSMFromPath()
{	cout << __FUNCTION__ << endl;

//GPU側
	LauchUpdateSMFromPath
	(
		sm_clusterNum,
		sd_sphPrtPos,
		sd_sphPrtVel,
		//SM-------------------------------------------------
		Ice_SM::GetOrgPosPointer(),
		Ice_SM::GetDevicePosPointer(),
		Ice_SM::GetOrgCmPointer(),
		Ice_SM::GetCurCmPointer(),
		Ice_SM::GetDeviceApqPointer(),
		Ice_SM::GetDeviceVelPointer(),
		Ice_SM::GetDeviceIndexesPointer(),
		Ice_SM::GetDeviceIndexSetPointer(),
		//Path-----------------------------------------------
		m_SurfSm.GetDevicePRTtoPTHPointer(),
		m_SurfSm.GetDevicePTHandPrfxSetPointer(),
		m_SurfSm.GetDecvicePrfxPos(),
		m_SurfSm.GetDecvicePrfxApq(),
		//IceStruct------------------------------------------
		m_iceStrct->GetDeviceCtoPointer(),
		m_iceStrct->GetDeviceCtoPNumPointer(),
		m_iceStrct->GetCtoPMax(),
		2,
		0.01		//TODO: 0.02では？いろいろばらばらっぽい
	);

//CPU側
	//prefixSumの更新
	m_SurfSm.UpdatePrefixSum();

	//クラスタのパラメータ更新
	#pragma omp parallel for
	for(int i = 0; i < sm_clusterNum; ++i)
	{	
		if(GetPtoCNum(i) == 0){	continue;	}
		m_iceMove[i]->SetNowCm(m_SurfSm.CalcCmSum(i));		//重心の更新
		m_iceMove[i]->SetApq(m_SurfSm.CalcApqSum(i));		//変形行列の更新
	}

	////CurCmの比較
	//float* dCurCm = new float[MAXCLUSTER * SM_DIM];
	//cudaMemcpy(dCurCm, Ice_SM::GetCurCmPointer(), sizeof(float) * MAXCLUSTER * SM_DIM, cudaMemcpyDeviceToHost);

	//for(int iCrstr = 0; iCrstr < MAXCLUSTER; iCrstr++)
	//{
	//	int cIndx = iCrstr * 3;
	//	cout << "GPU:[" << iCrstr << "] = " << endl;
	//	cout <<	"(" << dCurCm[cIndx+0] << ", " << dCurCm[cIndx+1] << ", " << dCurCm[cIndx+2] << ")" << endl;
	//	cout << endl;

	//	Vec3 cm = m_iceMove[iCrstr]->GetCm();
	//	cout << "CPU:[" << iCrstr << "] = " << endl;
	//	cout <<	"(" << cm[0] << ", " << cm[1] << ", " << cm[2] << ")" << endl;
	//	cout << endl;
	//}

	//delete[] dCurCm;

	////CurApqの比較
	//float* dCurApq = new float[MAXCLUSTER * 9];
	//cudaMemcpy(dCurApq, Ice_SM::GetDeviceApqPointer(), sizeof(float) * MAXCLUSTER * 9, cudaMemcpyDeviceToHost);

	//for(int iCrstr = 0; iCrstr < MAXCLUSTER; iCrstr++)
	//{
	//	int cIndx = iCrstr * 9;
	//	cout << "GPU:[" << iCrstr << "] = " << endl;
	//	cout <<	"(" << dCurApq[cIndx+0] << ", " << dCurApq[cIndx+1] << ", " << dCurApq[cIndx+2] << ")" << endl;
	//	cout <<	"(" << dCurApq[cIndx+3] << ", " << dCurApq[cIndx+4] << ", " << dCurApq[cIndx+5] << ")" << endl;
	//	cout <<	"(" << dCurApq[cIndx+6] << ", " << dCurApq[cIndx+7] << ", " << dCurApq[cIndx+8] << ")" << endl;
	//	cout << endl;

	//	rxMatrix3 apq = m_iceMove[iCrstr]->GetApq();
	//	cout << "CPU:[" << iCrstr << "] = " << endl;
	//	cout <<	"(" << apq(0, 0) << ", " << apq(0, 1) << ", " << apq(0, 2) << ")" << endl;
	//	cout <<	"(" << apq(1, 0) << ", " << apq(1, 1) << ", " << apq(1, 2) << ")" << endl;
	//	cout <<	"(" << apq(2, 0) << ", " << apq(2, 1) << ", " << apq(2, 2) << ")" << endl;
	//	cout << endl;
	//}

	//delete[] dCurApq;

}

void IceObject::TestStepInterPolation()
{
	Vec3 pos,vel;

	#pragma omp parallel for private(pos, vel)
	for(int i = 0; i < sm_particleNum; ++i)
	{
		if(GetParticleNum() <= i){	continue;	}	//融解のみの実験のときに必要になる．
		if(GetPtoCNum(i) <= 0){		continue;	}
	
		//固体運動の最終位置計算
		CalcAverageCPU(i, pos, vel);

		//液体と固体の補間
		LinerInterPolationCPU(i, pos, vel);			//線形補間
	}
}