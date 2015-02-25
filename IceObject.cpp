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

	InitTetra();									//四面体の初期化							
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

//影響半径を元に運動計算するクラスタを選択　なるべく疎になるように選ぶ
void IceObject::InitSelectCluster(float radius)
{
	m_iceStrct->SetSelectRadius(radius);
	m_iceStrct->InitSelectCluster(m_iceSM);

	//DebugUpdateSelectCluster();
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
void IceObject::InitCluster(Vec3 boundarySpaceHigh, Vec3 boundarySpaceLow, float timeStep, int itr, const vector<vector<rxNeigh>>& neights)
{	cout << __FUNCTION__ << endl;

	//変数の初期化
	for(vector<Ice_SM*>::iterator it = m_iceSM.begin(); it != m_iceSM.end(); ++it)
	{
		if(*it) delete *it;
	}
	
	sm_clusterNum = 0;	//クラスタ数の初期化

	Ice_SM::SetIterationNum(itr);
	Ice_SM::SetItrStiffness(1.0f);	//TODO::適当
	Ice_SM::SetPrtPointerPosAndVel(s_sphPrtPos, s_sphPrtVel);

cout << __FUNCTION__ << ", check1" << endl;

	//四面体リストを元に，粒子毎にクラスタ作成
	for(int i = 0; i < sm_particleNum; ++i)
	{
		//クラスタ初期化
		m_iceSM.push_back(new Ice_SM(sm_clusterNum));
		m_iceSM[sm_clusterNum]->SetSimulationSpace(boundarySpaceLow, boundarySpaceHigh);
		m_iceSM[sm_clusterNum]->SetTimeStep(timeStep);
		m_iceSM[sm_clusterNum]->SetCollisionFunc(0);
		m_iceSM[sm_clusterNum]->SetStiffness(1.0, 0.0);

		//四面体リストを元に，クラスタへ粒子を登録
#ifdef USE_NEIGHT
		SetClusterMoveInfoFromNeight(i, neights);
#else
		SetClusterMoveInfo(i);
#endif
		sm_clusterNum++;
	}

cout << __FUNCTION__ << ", check2" << endl;

	//質量の修正
	UpdateParticleMass_Normal();		//一定に
	//UpdateParticleMass_Average();		//平均　LSMの方法
	//UpdateParticleMass_Direction();	//方向ベクトルとの類似度で重み付け　Ijiriらの方法

	////MakeOneCluster();
	
	Ice_SM::InitFinalParamPointer(sm_clusterNum);

//デバッグ
	//DebugClusterInfo();
	//DebugNeights(neights);
}

//オブジェクトの構造の初期化
void IceObject::InitStrct()
{	cout << __FUNCTION__ << endl;

	//カウント
	for(int i = 0; i < sm_particleNum; i++)
	{
		int pNum = m_iceSM[i]->GetNumVertices();
		vector<int> pList;

		for(int j = 0; j < pNum; j++)
		{
			pList.push_back(m_iceSM[i]->GetParticleIndx(j));
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

	m_iceStrct->InitSelectClusterFromClusterSet(m_iceSM);

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

//運動系やそれぞれの手法のための初期化
void IceObject::InitMoveMethod()
{	cout << __FUNCTION__ << endl;

	//NULLで初期化
	m_fpStepObjMove = NULL;
	m_iceClsuterMove = NULL;
	m_iceJudeMove = NULL;
	m_iceConvolution = NULL;
	m_iceConvoJudge = NULL;
	m_iceCalcMethod = NULL;

	//全てノーマル
	ChangeMode_Normal();
	ChangeMode_ClusterMove_Normal();
	ChangeMode_JudgeMove_Normal();
	ChangeMode_Convolution_Normal();
	ChangeMode_IntrpJudge_Normal();
	ChangeMode_CalcMethod_Normal();

	m_SurfSm = NULL;

	IceObject::InitInterPolation();					//線形補間の初期化

#ifdef MK_USE_GPU
	InitGPU();										//構造の情報が完成したらGPUを初期化
#endif
}

//高速化に用いるパスの初期化
void IceObject::InitPath()
{
	m_SurfSm = new Surf_SM(s_sphPrtPos, s_sphPrtVel, m_iceSM, m_iceStrct, sm_clusterNum, sm_particleNum);
}

//モード切替
void IceObject::ChangeMode_Debug()
{	FUNCTION_NAME

	m_fpStepObjMove = &IceObject::StepObjMoveCPUDebug;

	//データフォルダの中身を全て削除
	string path = RESULT_DATA_PATH;
	DebugClearDirectry(path);
}

void IceObject::ChangeMode_Normal()
{	FUNCTION_NAME
	
	m_fpStepObjMove = &IceObject::StepObjMoveCPUNormal;
}

//計算手法
void IceObject::ChangeMode_CalcMethod_Normal()
{	FUNCTION_NAME

	if(m_iceCalcMethod != NULL) delete m_iceCalcMethod;
	m_iceCalcMethod = new Ice_CalcMethod_Normal(m_iceClsuterMove);
}

void IceObject::ChangeMode_CalcMethod_Itr_Num()
{	FUNCTION_NAME

	if(m_iceCalcMethod != NULL) delete m_iceCalcMethod;
	m_iceCalcMethod = new Ice_CalcMethod_Iteration(m_iceSM, m_iceClsuterMove, m_iceConvolution);
}

void IceObject::ChangeMode_CalcMethod_Itr_Stiff()
{	FUNCTION_NAME;
	
	if(m_iceCalcMethod != NULL) delete m_iceCalcMethod;
	m_iceCalcMethod = new Ice_CalcMethod_Itr_Stiffness(m_iceSM, m_iceClsuterMove, m_iceConvolution);
}

void IceObject::ChangeMode_CalcMethod_Itr_Expand()
{	FUNCTION_NAME;

	if(m_iceCalcMethod != NULL) delete m_iceCalcMethod;
	m_iceCalcMethod = new Ice_CalcMethod_Itr_Expand(m_iceSM, m_iceClsuterMove, m_iceConvolution);
}

void IceObject::ChangeMode_CalcMethod_Itr_Exp_Stiff()
{	FUNCTION_NAME;

	if(m_iceCalcMethod != NULL) delete m_iceCalcMethod;
	m_iceCalcMethod = new Ice_CalcMethod_Itr_Exp_Stiff(m_iceSM, m_iceClsuterMove, m_iceConvolution);
}


//運動計算クラスタ選択
void IceObject::ChangeMode_JudgeMove_Normal()
{	FUNCTION_NAME

	ResetSelectCluster();		//選択的運動計算の初期化

	if(m_iceJudeMove != NULL) delete m_iceJudeMove;
	m_iceJudeMove = new Ice_JudgeMove_Normal(m_iceSM);

	m_iceClsuterMove->SetJudgeMove(m_iceJudeMove);
}

void IceObject::ChangeMode_JudgeMove_Spears()
{	FUNCTION_NAME

	if(m_iceJudeMove != NULL) delete m_iceJudeMove;
	m_iceJudeMove = new Ice_JudgeMove_Spears(m_iceSM, m_iceStrct);

	m_iceClsuterMove->SetJudgeMove(m_iceJudeMove);
}

//運動計算
void IceObject::ChangeMode_ClusterMove_Normal()
{	FUNCTION_NAME

	if(m_iceClsuterMove != NULL) delete m_iceClsuterMove;
	m_iceClsuterMove = new Ice_ClusterMove_Normal(m_iceSM);

	m_iceClsuterMove->SetJudgeMove(m_iceJudeMove);

	if(m_iceCalcMethod != NULL)	m_iceCalcMethod->SetObjMove(m_iceClsuterMove);
}

void IceObject::ChangeMode_ClusterMove_Path()
{	FUNCTION_NAME

	if(m_SurfSm == NULL){	InitPath();	}	//高速な運動計算のためのパスを準備

	if(m_iceClsuterMove != NULL) delete m_iceClsuterMove;
	m_iceClsuterMove = new Ice_ClusterMove_FastPath(m_iceSM, m_SurfSm);

	m_iceClsuterMove->SetJudgeMove(m_iceJudeMove);

	if(m_iceCalcMethod != NULL)	m_iceCalcMethod->SetObjMove(m_iceClsuterMove);
}

//最終結果補間クラスタ選択
void IceObject::ChangeMode_IntrpJudge_Normal()
{	cout << __FUNCTION__ << endl;

	if(m_iceConvoJudge != NULL) delete m_iceConvoJudge;
	m_iceConvoJudge = new Ice_ConvoJudge_Normal(m_iceSM);

	m_iceConvolution->SetConvoJudge(m_iceConvoJudge);
}

void IceObject::ChangeMode_IntrpJudge_Spears()
{	cout << __FUNCTION__ << endl;

	if(m_iceConvoJudge != NULL) delete m_iceConvoJudge;
	m_iceConvoJudge = new Ice_ConvoJudge_Spears(m_iceSM, m_iceStrct);

	m_iceConvolution->SetConvoJudge(m_iceConvoJudge);
}

//最終結果補間
void IceObject::ChangeMode_Convolution_Normal()
{	cout << __FUNCTION__ << endl;

	if(m_iceConvolution != NULL) delete m_iceConvolution;
	m_iceConvolution = new Ice_Convolution_Normal(m_iceSM, m_iceStrct);

	m_iceConvolution->SetConvoJudge(m_iceConvoJudge);
}

//変形量で重みをつけて補間
void IceObject::ChangeMode_Convolution_Weight()
{	cout << __FUNCTION__ << endl;

	if(m_iceConvolution != NULL) delete m_iceConvolution;
	m_iceConvolution = new Ice_Convolution_Weight(m_iceSM, m_iceStrct);

	m_iceConvolution->SetConvoJudge(m_iceConvoJudge);
}

//方向ベクトルとの類似度で重みを付けて補間
void IceObject::ChangeMode_Convolution_Anisotropic()
{
	FUNCTION_NAME;

	if(m_iceConvolution != NULL) delete m_iceConvolution;
	m_iceConvolution = new Ice_Convolution_Anisotropic(m_iceSM, m_iceStrct);

	m_iceConvolution->SetConvoJudge(m_iceConvoJudge);
}

//クラスタの運動計算情報の登録
void IceObject::SetClusterMoveInfo(int pIndx)
{
	//最初の一つ目は必ずpIndxとする
	int pNum = m_iceSM[pIndx]->GetNumVertices();
	float mass = 1.0f;

	m_iceSM[pIndx]->AddVertex( Vec3(s_sphPrtPos[pIndx*4+0], s_sphPrtPos[pIndx*4+1], s_sphPrtPos[pIndx*4+2] ), mass, pIndx);
	m_iceSM[pIndx]->SetAlphas(pNum, 1.0);
	m_iceSM[pIndx]->SetBetas (pNum, 0.0);
	m_iceSM[pIndx]->SetLayer (pNum, 0);

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
			if(m_iceSM[pIndx]->CheckIndx(jpIndx)){	continue;	}

			int pNum = m_iceSM[pIndx]->GetNumVertices();
			float mass = 1.0f;
			int sphIndx = jpIndx*4;

			m_iceSM[pIndx]->AddVertex( Vec3(s_sphPrtPos[sphIndx+0], s_sphPrtPos[sphIndx+1], s_sphPrtPos[sphIndx+2] ), mass, jpIndx);
			m_iceSM[pIndx]->SetAlphas(pNum, 1.0);
			m_iceSM[pIndx]->SetBetas (pNum, 0.0);
			m_iceSM[pIndx]->SetLayer (pNum, 0);
		}
	}

	if(sm_layerNum == 0){	return;	}

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
				if(m_iceSM[pIndx]->CheckIndx(kpIndx)){	/*cout << "contineue kpIndx = " << kpIndx << endl;*/ continue;	}

				int pNum = m_iceSM[pIndx]->GetNumVertices();
				float mass = 1.0f;

				m_iceSM[pIndx]->AddVertex( Vec3(s_sphPrtPos[kpIndx*4+0], s_sphPrtPos[kpIndx*4+1], s_sphPrtPos[kpIndx*4+2] ), mass, kpIndx);
				m_iceSM[pIndx]->SetAlphas(pNum, 1.0);
				m_iceSM[pIndx]->SetBetas (pNum, 0.0);
				m_iceSM[pIndx]->SetLayer (pNum, jlIndx);
				//cout << "pIndx = " << pIndx << " GetNumVertices = " << m_iceSM[pIndx]->GetNumVertices() << endl;
			}
		}
	}
}

void IceObject::SetClusterMoveInfoFromNeight(int pIndx, const vector<vector<rxNeigh>>& neights)
{
	//最初の一つ目は必ずpIndxとする
	int pNum = m_iceSM[pIndx]->GetNumVertices();
	float mass = 1.0f;

	m_iceSM[pIndx]->AddVertex( Vec3(s_sphPrtPos[pIndx*4+0], s_sphPrtPos[pIndx*4+1], s_sphPrtPos[pIndx*4+2]), mass, pIndx);
	m_iceSM[pIndx]->SetAlphas(pNum, 1.0);
	m_iceSM[pIndx]->SetBetas (pNum, 0.0);
	m_iceSM[pIndx]->SetLayer (pNum, 0);

	//近傍粒子をクラスタに追加
	for(unsigned id_np = 0; id_np < neights[pIndx].size(); id_np++)
	{
		int np_pIndx = neights[pIndx][id_np].Idx;
		int pNum = m_iceSM[pIndx]->GetNumVertices();
		double mass = 1.0;

		if(m_iceSM[pIndx]->CheckIndx(np_pIndx)){	continue;	}

		m_iceSM[pIndx]->AddVertex( Vec3(s_sphPrtPos[np_pIndx*4+0], s_sphPrtPos[np_pIndx*4+1], s_sphPrtPos[np_pIndx*4+2]), mass, np_pIndx);
		m_iceSM[pIndx]->SetAlphas(pNum, 1.0);
		m_iceSM[pIndx]->SetBetas (pNum, 0.0);
		m_iceSM[pIndx]->SetLayer (pNum, 0);
	}
}

//粒子→クラスタやクラスタ→粒子を登録
void IceObject::SetClusterStrctInfo(int cIndx, int* PtoCNum)
{
	//粒子が属している四面体の番号を登録するための準備
	//pCountListには，cIndx番目のクラスタに含まれる各粒子が，それぞれいくつのクラスタに属するかを求めて保存する
	int* pCountList = new int[m_iceSM[cIndx]->GetNumVertices()];
	int* pLayerList = new int[m_iceSM[cIndx]->GetNumVertices()];			//粒子の所属レイヤー
	
	for(int i = 0; i < m_iceSM[cIndx]->GetNumVertices(); i++)
	{
		int pIndx = m_iceSM[cIndx]->GetParticleIndx(i);

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

		pLayerList[i] = m_iceSM[cIndx]->GetLayer(i);					//粒子が何層目の近傍なのかを取得
	}

	//粒子とクラスタの情報登録
	vector<int> pIndxList;

	for(int i = 0; i < GetCtoPNum(cIndx); i++)
	{
		int pIndx = m_iceSM[cIndx]->GetParticleIndx(i);
			
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
	Ice_SM::InitGPU(m_iceSM, sd_sphPrtPos, sd_sphPrtVel, sm_particleNum, sm_maxParticleNum);

#if defined(USE_PATH)
{
	m_SurfSm->InitPathGPU();		//高速化に用いるパスの初期化
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

//デバッグ情報の描画
void IceObject::Display()
{
	//テスト
	//現在のCalcMethodのモード確認
	const char* calcMethodName1 = typeid(Ice_CalcMethod_Itr_Stiffness).name();
	const char* calcMethodName2 = typeid(Ice_CalcMethod_Itr_Exp_Stiff).name();
	const char* nowName = typeid(*(m_iceCalcMethod)).name();

	//現在iteration_stiffnessモードかの確認　strcmpは文字列が一致すると0になる
	if(strcmp(calcMethodName1, nowName) != 0
	&& strcmp(calcMethodName2, nowName) != 0) return ;

	//rigidオブジェクトの位置を描画してみる
	if(strcmp(calcMethodName1, nowName) == 0){
		((Ice_CalcMethod_Itr_Stiffness*)m_iceCalcMethod)->DebugStiffness();
	}
	else if(strcmp(calcMethodName2, nowName) == 0){
		((Ice_CalcMethod_Itr_Exp_Stiff*)m_iceCalcMethod)->DebugStiffness();
	}
}

//固体の運動計算，計算結果の統合
//処理内容は別クラスに記述
void IceObject::StepObjMoveCPU()
{
	//CPUNormalかCPUDebug
	(this->*m_fpStepObjMove)();
}

void IceObject::StepObjMoveCPUNormal()
{
	//固体の運動計算
	m_iceCalcMethod->StepObjMove();			RXTIMER("StepObjMove");

	//固体の最終位置決定
	m_iceConvolution->StepConvolution();	RXTIMER("StepConvolution");

	//液体と固体の線形補間
	StepInterPolation();					RXTIMER("StepInterPolation");
	//StepInterPolationKenjya();
}

void IceObject::StepObjMoveCPUDebug()
{
	//固体の運動計算
	m_iceCalcMethod->StepObjMoveDebug();

	//固体の最終位置決定
	m_iceConvolution->StepConvolutionDebug();

	//液体と固体の線形補間
	StepInterPolation();

	//データ作成
	//TODO::後にクラス化する
	DebugDeformationAmount();
	DebugDeformationAverage();
}

//液体と固体の線形補間
void IceObject::StepInterPolation()
{
	for(int pIndx = 0; pIndx < sm_particleNum; pIndx++)
	{
		if(m_iceSM[pIndx]->GetNumVertices() == 0){			continue;	}

		int sphIndx = pIndx*4;
		int sldIndx = pIndx*3;
		double intrps = 1.0 - m_fInterPolationCoefficience[pIndx];	//補間係数
	
		float* sldVel = Ice_SM::GetSldVelPointer();
		float* sldPos = Ice_SM::GetSldPosPointer();
	
		for(int i = 0; i < 3; i++)
		{
			s_sphPrtVel[sphIndx+i] = sldVel[sldIndx+i] * m_fInterPolationCoefficience[pIndx] + s_sphPrtVel[sphIndx+i] * intrps;
			s_sphPrtPos[sphIndx+i] = sldPos[sldIndx+i] * m_fInterPolationCoefficience[pIndx] + s_sphPrtPos[sphIndx+i] * intrps;
		}
	}
}

void IceObject::StepInterPolationKenjya()
{
	//総変形量
	float defSum = 0.0f;
	for(int pIndx = 0; pIndx < m_iceSM.size(); pIndx++){
		defSum += m_iceSM[pIndx]->GetDefAmount();
	}

	for(int pIndx = 0; pIndx < sm_particleNum; pIndx++)
	{
		if(m_iceSM[pIndx]->GetNumVertices() == 0){			continue;	}

		int sphIndx = pIndx*4;
		int sldIndx = pIndx*3;

		float defAmount = m_iceSM[pIndx]->GetDefAmount() / defSum;
		if(defAmount > 1.0f) defAmount = 1.0f;
		if(defAmount < 0.1f) defAmount = 0.1f;
		double intrps = 1.0 - defAmount;	//補間係数
	
		float* sldVel = Ice_SM::GetSldVelPointer();
		float* sldPos = Ice_SM::GetSldPosPointer();
	
		for(int i = 0; i < 3; i++)
		{
			s_sphPrtVel[sphIndx+i] = sldVel[sldIndx+i] * defAmount + s_sphPrtVel[sphIndx+i] * intrps;
			s_sphPrtPos[sphIndx+i] = sldPos[sldIndx+i] * defAmount + s_sphPrtPos[sphIndx+i] * intrps;
		}
	}
}

//GPU
void IceObject::StepObjMoveGPU()
{
	#ifdef USE_ITR
	//TODO: 未実装
	//StepObjCalcWidhIteration();	RXTIMER("StepObjCalcWidhIteration");
#else

	#ifdef USE_PATH
		StepObjMoveGPUUsePath();	RXTIMER("StepObjMoveGPUUsePath");
	#else
		StepObjMoveGPUNormal();		RXTIMER("StepObjMoveGPU");			//クラスタの運動計算
	#endif
		
	StepInterPolationNormal();		RXTIMER("StepInterpolationGPU");	//総和計算，線形補間

#endif
}

void IceObject::StepObjMoveGPUNormal()
{
	Ice_SM::SetDevicePosPointer(sd_sphPrtPos);	//VBOを使っているので，これがないとエラーが出るのに注意
	Ice_SM::UpdateGPU();
}

void IceObject::StepObjMoveGPUUsePath()
{
	//prefixSumの更新
	m_SurfSm->SetDevicePointer(sd_sphPrtPos, sd_sphPrtVel);	//VBOを使っているので，これがないとエラーが出るのに注意
	m_SurfSm->UpdatePrefixSumGPU();
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
		m_SurfSm->GetDevicePRTtoPTHPointer(),
		m_SurfSm->GetDevicePTHandPrfxSetPointer(),
		m_SurfSm->GetDecvicePrfxPos(),
		m_SurfSm->GetDecvicePrfxApq(),
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

	StepInterPolationNormal();
}

void IceObject::StepInterPolationNormal()
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

#endif
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
#endif
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

	SearchMeltParticle(prtList);		//RXTIMER("SearchMeltParticle");	//融解粒子探索

	//ここが計算時間の殆どを占めている
	//m_iceStrct->StepObjMelt
	//	(	prtList,
	//		clusterList,
	//		tetraList,
	//		cLayerList,
	//		tLayerList);				//融解処理
	//RXTIMER("m_iceStrct->StepObjMelt");

	//m_iceStrct->TestStepObjMelt
	//	(	prtList,
	//		clusterList,
	//		tetraList,
	//		cLayerList,
	//		tLayerList);				//融解処理
	
	//運動計算するクラスタの再定義
	if(prtList.size() == 0) return;

QueryCounter counter1;
QueryCounter counter2;

counter1.Start();
	vector<unsigned> neighborList;
	//CashNeighborList(prtList, neighborList);
	ReConstructCluster(prtList, clusterList);	//RXTIMER("ReConstructCluster");	//クラスタの再構築
double end1 = counter1.End();

counter2.Start();
	m_iceStrct->UpdateSelectCluster(prtList, neighborList, m_iceSM);
double end2 = counter2.End();

	cout << "Time" << endl;
	cout << "ReConstructCluster:	" << end1 << endl;
	cout << "UpdateSelectCLuster:	" << end2 << endl;
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
		if(m_iceSM[i]->GetNumVertices() == 0)			continue;	//クラスタを持っていない
		if( m_heatTransfer->getPhaseChange(i) != 1 )	continue;	//相転移の条件を満たしている
		if( m_heatTransfer->getPhase(i) != 2 )			continue;	//水へと相転移している
		//if( m_iceStrct->GetPtoCNum(i) == 0 )			continue;	//クラスタに含まれていない
		//if( m_iceStrct->GetPtoTNum(i) == 0 )			continue;	//四面体に含まれていない

		//if(pList.size() > 500){	break;	}							//融解粒子数の制限

		m_fInterPolationCoefficience[i] = 0.0f;						//線形補間しない
		m_heatTransfer->setPhaseChange(i, 0);						//相転移し終わったことを伝える
		pList.push_back(i);											//融解粒子の記録
	}

//デバッグ
	if(pList.size() == 0) return;

	cout << __FUNCTION__ << endl;

	//for(unsigned i = 0; i < pList.size(); i++)
	//{
	//	cout << pList[i] << endl;
	//}

	cout << "pList.size() = " << pList.size() << endl;
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

	if(pListSize == 0 /*|| cListSize == 0*/){	return; }

	////あとで二分探索するためにpListをsortする
	//std::sort(pList.begin(), pList.end());

	vector<bool> clusterFlags(sm_particleNum, true);
	for(unsigned i = 0; i < pList.size(); i++)
	{
		clusterFlags[pList[i]] = false;
	}

	//クラスタから粒子を取り除く
	for(int i = 0; i < pListSize; i++)
	{
		int pIndx = pList[i];

		//融解クラスタのlayer0を取得
		vector<int> connectPIndx;
		for(unsigned oIndx = 0; oIndx < GetMoveObj(pIndx)->GetIndxNum(); oIndx++)
		{
			int layer = GetMoveObj(pIndx)->GetLayer(oIndx);

			if(layer == -1) continue;
			if(layer != 0)	break;

			int ipIndx = GetMoveObj(pIndx)->GetParticleIndx(oIndx);

			connectPIndx.push_back(ipIndx);
		}

		//その他のクラスタから融解粒子を取り除く
		int oIndx = 0;
		int ccpIndx = 0;
		int coIndx = 0;
		int ccIndx = 0;
		int connectPIndxSize = connectPIndx.size();
		Ice_SM* moveObj;
		vector<unsigned>::iterator begin = pList.begin();
		vector<unsigned>::iterator end = pList.end();

		#pragma omp parallel
		{
		#pragma omp for private(oIndx, ccpIndx, coIndx, moveObj, ccIndx, begin, end)
			for(int cIndx = 0; cIndx < sm_particleNum; cIndx++)
			{
				//if(std::binary_search(begin, end, (unsigned)cIndx)){	continue;	}	//重い
				if(clusterFlags[cIndx] == false)	continue;

				moveObj = GetMoveObj(cIndx);
				oIndx = moveObj->SearchIndx(pIndx);

				if(oIndx == MAXINT) continue;	//存在しない

				moveObj->Remove(oIndx);

#ifndef USE_NEIGHT
				//融解粒子のlayer0となる粒子を全て取り除く
				//ただし自分のlayer0は残す
				for(ccIndx = 0; ccIndx < connectPIndxSize; ccIndx++)
				{
					ccpIndx = connectPIndx[ccIndx];
					coIndx = moveObj->SearchIndx(ccpIndx);
	
					if(coIndx == MAXINT)				continue;	//存在しない
					if(moveObj->GetLayer(coIndx) == 0)	continue;	//自分のlayer0は残す

					moveObj->Remove(coIndx);
				}
#endif
				moveObj->CalcOrgCm();
				moveObj->ClassifyAllOrgParticle();
			}
		}//end #pragma omp parallel
	}

	#pragma omp parallel
	{
	#pragma omp for
		for(int pIndx = 0; pIndx < pListSize; pIndx++)
		{
			GetMoveObj(pList[pIndx])->Clear();
		}
	}

	//vector<int> checkTList;
	//int j = 0, k = 0, l = 0;
	//int icIndx = 0;
	//int jtIndx = 0, joIndx = 0;
	//int kpIndx = 0, ktIndx = 0, klIndx = 0;
	//int lpIndx = 0;
	//int pNum = 0;

	////再定義クラスタの初期化
	//#pragma omp parallel
	//{
	//#pragma omp for
	//	for(int i = 0; i < cListSize; ++i)
	//	{
	//		if(std::find(pList.begin(), pList.end(), cList[i]) != pList.end()){	continue;	}
	//		GetMoveObj(cList[i])->Clear();
	//	}
	//}//end #pragma omp parallel

	////クラスタの再定義
	//#pragma omp parallel
	//{
	//#pragma omp for private(checkTList, j, k, l, icIndx, jtIndx, joIndx, kpIndx, ktIndx, klIndx, lpIndx, pNum)
	//	for(int i = 0; i < cListSize; ++i)
	//	{
	//		checkTList.clear();
	//		icIndx = cList[i];
	//		if(std::find(pList.begin(), pList.end(), icIndx) != pList.end()){	continue;	}

	//		//クラスタを再定義する際，基本となる粒子が属する四面体から初期粒子を得る．
	//		//クラスタ番号＝＝粒子番号なのに注意
	//		//以下を関数にすると，エラーが出てうまくいかない
	//		for(j = 0; j < GetPtoTIndx(icIndx); j++)
	//		{
	//			jtIndx = GetPtoT(icIndx, j, 0);
	//			joIndx = GetPtoT(icIndx, j, 1);

	//			if(jtIndx == -1 || joIndx == -1){ continue;	}
	//			if(std::find(checkTList.begin(), checkTList.end(), jtIndx) != checkTList.end())
	//			{	continue;	}
	//			else
	//			{	checkTList.push_back(jtIndx);	}

	//			//四面体の全ての粒子をクラスタに登録
	//			for(k = 0; k < GetTtoPIndx(jtIndx); k++)
	//			{
	//				kpIndx = GetTtoP(jtIndx, k);

	//				if(kpIndx == -1){	continue;	}
	//				if(GetMoveObj(icIndx)->CheckIndx(kpIndx)){	continue;	}

	//				pNum = GetMoveObj(icIndx)->GetNumVertices();

	//				GetMoveObj(icIndx)->AddVertex(Vec3(s_sphPrtPos[kpIndx*4+0], s_sphPrtPos[kpIndx*4+1], s_sphPrtPos[kpIndx*4+2]), 1.0, kpIndx);
	//				//GetMoveObj(icIndx)->SetVelocity(pNum, Vec3(s_sphPrtVel[kpIndx*4+0], s_sphPrtVel[kpIndx*4+1], s_sphPrtVel[kpIndx*4+2]));
	//				GetMoveObj(icIndx)->SetAlphas(pNum, 1.0);
	//				GetMoveObj(icIndx)->SetBetas (pNum, 0.0);
	//				GetMoveObj(icIndx)->SetLayer (pNum, 0);
	//			}
	//		}

	//		//近傍四面体のlayerたどり，粒子を追加していく
	//		//TODO::不安定になるなら，layerが遠いほどbetaを下げる
	//		for(j = 0; j < GetPtoTIndx(icIndx); ++j)
	//		{
	//			jtIndx = GetPtoT(icIndx, j, 0);
	//			joIndx = GetPtoT(icIndx, j, 1);

	//			if(jtIndx == -1 || joIndx == -1){ continue;	}

	//			for(k = 0; k < GetNTNum(jtIndx); k++)
	//			{
	//				ktIndx = GetNeighborTetra(jtIndx, k, 0);
	//				klIndx = GetNeighborTetra(jtIndx, k, 1);

	//				if(ktIndx == -1 || klIndx == -1){	continue;	}
	//				if(std::find(checkTList.begin(), checkTList.end(), ktIndx) != checkTList.end())
	//				{	continue;	}
	//				else
	//				{	checkTList.push_back(ktIndx);	}

	//				//四面体の全ての粒子をクラスタに登録
	//				for(l = 0; l < GetTtoPIndx(ktIndx); l++)
	//				{
	//					lpIndx = GetTtoP(ktIndx, l);

	//					if(lpIndx == -1){	continue;	}
	//					if(GetMoveObj(icIndx)->CheckIndx(lpIndx)){	/*cout << "contineue kpIndx = " << kpIndx << endl;*/ continue;	}

	//					pNum = GetMoveObj(icIndx)->GetNumVertices();

	//					GetMoveObj(icIndx)->AddVertex( Vec3(s_sphPrtPos[lpIndx*4+0], s_sphPrtPos[lpIndx*4+1], s_sphPrtPos[lpIndx*4+2] ), 1.0, lpIndx);
	//					//GetMoveObj(icIndx)->SetVelocity(pNum, Vec3(s_sphPrtVel[lpIndx*4+0], s_sphPrtVel[lpIndx*4+1], s_sphPrtVel[lpIndx*4+2]));
	//					GetMoveObj(icIndx)->SetAlphas(pNum, 1.0);
	//					GetMoveObj(icIndx)->SetBetas (pNum, 0.0);
	//					GetMoveObj(icIndx)->SetLayer (pNum, klIndx);
	//				}
	//			}
	//		}
	//	}
	//}//end #pragma omp parallel


	////融解クラスタをリセット
	//#pragma omp parallel
	//{
	//#pragma omp for
	//	for(int i = 0; i < pListSize; ++i)
	//	{
	//		GetMoveObj(pList[i])->Clear();
	//	}
	//}//end #pragma omp parallel


	//////クラスタの固体情報を初期化
	////for(int i = 0; i < (int)cList.size(); ++i)
	////{
	////	m_iceStrct->ClearCtoP(cList[i]);
	////}

	//////クラスタの固体情報を再度登録
	////for(unsigned i = 0; i < cList.size(); ++i)
	////{
	////	if(std::find(pList.begin(), pList.end(), cList[i]) != pList.end())
	////	{
	////		continue;	
	////	}

	////	int cIndx = cList[i];

	////	vector<int> pIndxList;
	////	int pNum = GetMoveObj(cIndx)->GetNumVertices();
	////	int* pLayerList = new int[pNum];									//粒子の所属レイヤー
	////
	////	//layerのためのコピー
	////	for(int i = 0; i < pNum; i++)
	////	{
	////		pLayerList[i] = GetMoveObj(cIndx)->GetLayer(i);		//粒子が何層目の近傍なのかを取得
	////	}
	////
	////	for(int i = 0; i < pNum; i++)
	////	{
	////		int pIndx = GetMoveObj(cIndx)->GetParticleIndx(i);
	////		int freeIndx = m_iceStrct->GetPtoCFreeIndx(pIndx);
	////
	////		m_iceStrct->SetPtoC(pIndx, freeIndx, cIndx, i, pLayerList[i]);		//粒子が所属している四面体を登録
	////		pIndxList.push_back(pIndx);
	////	}
	////
	////	m_iceStrct->SetCtoP(cIndx, pIndxList, pLayerList);						//四面体が含んでいる粒子を登録
	////
	////	delete[] pLayerList;
	////}
}

void IceObject::ResetSelectCluster()
{
	//全てのクラスタを選択する
	for(int cIndx = 0; cIndx < sm_particleNum; cIndx++)
	{
		m_iceStrct->UpdateMotionCalcCluster(cIndx, 1);		
	}
}

void IceObject::CashNeighborList(const vector<unsigned>& prtList, vector<unsigned>& neighborList)
{
	for(unsigned lIndx = 0; lIndx < prtList.size(); lIndx++)
	{
		unsigned cIndx = prtList[lIndx];

		for(int indx = 0; indx < (int)m_iceSM[cIndx]->GetIndxNum(); indx++)
		{
			int pIndx = m_iceSM[cIndx]->GetParticleIndx(indx);
			if(pIndx == MAXINT) continue;

			if(find(neighborList.begin(), neighborList.end(), (unsigned)pIndx) != neighborList.end()) continue;

			neighborList.push_back((unsigned)pIndx);
		}
	}
}

//質量修正　全て1.0f
void IceObject::UpdateParticleMass_Normal()
{
	unsigned pNum = IceObject::GetParticleNum();
	float mass = 1.0f;

	//各粒子の各クラスタでの質量を修正
	for(int cIndx = 0; cIndx < pNum; cIndx++)
	{
		for(int oIndx = 0; oIndx < m_iceSM[cIndx]->GetIndxNum(); oIndx++)
		{
			int pIndx = m_iceSM[cIndx]->GetParticleIndx(oIndx);
			if(pIndx == MAXINT){	continue;	}

			m_iceSM[cIndx]->SetMass(oIndx, mass);
		}
	}

	//各クラスタの重心を更新
	for(int pIndx = 0; pIndx < sm_particleNum; pIndx++)
	{
		m_iceSM[pIndx]->CalcOrgCm();
		m_iceSM[pIndx]->ClassifyAllOrgParticle();
	}
}

//質量修正　クラスタ数で平均
void IceObject::UpdateParticleMass_Average()
{
	unsigned pNum = IceObject::GetParticleNum();
	vector<unsigned> addParticleNum(pNum, 0);

	//各粒子がいくつのクラスタに含まれているか，のカウント
	for(int cIndx = 0; cIndx < pNum; cIndx++)
	{
		for(int oIndx = 0; oIndx < m_iceSM[cIndx]->GetIndxNum(); oIndx++)
		{
			int pIndx = m_iceSM[cIndx]->GetParticleIndx(oIndx);
			if(pIndx == MAXINT){	continue;	}

			//粒子数のカウント
			addParticleNum[pIndx] += 1;
		}
	}

	float mass = 1.0f;

	//各粒子の各クラスタでの質量を修正
	for(int cIndx = 0; cIndx < pNum; cIndx++)
	{
		for(int oIndx = 0; oIndx < m_iceSM[cIndx]->GetIndxNum(); oIndx++)
		{
			int pIndx = m_iceSM[cIndx]->GetParticleIndx(oIndx);
			if(pIndx == MAXINT){	continue;	}

			m_iceSM[cIndx]->SetMass(oIndx, mass/(float)addParticleNum[pIndx]);
		}
	}

	//各クラスタの重心を更新
	for(int pIndx = 0; pIndx < sm_particleNum; pIndx++)
	{
		m_iceSM[pIndx]->CalcOrgCm();
		m_iceSM[pIndx]->ClassifyAllOrgParticle();
	}
}

//質量修正　方向ベクトルとの近似度合
void IceObject::UpdateParticleMass_Direction()
{
	//各クラスタの重心を更新　質量1.0fで重心を計算
	for(int pIndx = 0; pIndx < sm_particleNum; pIndx++){	m_iceSM[pIndx]->CalcOrgCm();	}

	unsigned pNum = IceObject::GetParticleNum();
	vector<unsigned> addParticleNum(pNum, 0);

	//各粒子がいくつのクラスタに含まれているか，のカウント
	for(int cIndx = 0; cIndx < pNum; cIndx++)
	{
		for(int oIndx = 0; oIndx < m_iceSM[cIndx]->GetIndxNum(); oIndx++)
		{
			int pIndx = m_iceSM[cIndx]->GetParticleIndx(oIndx);
			if(pIndx == MAXINT){	continue;	}

			//粒子数のカウント
			addParticleNum[pIndx] += 1;
		}
	}

	Vec3 dirVec = Vec3(0.0f, -1.0f, 1.0f);
	dirVec = Unit(dirVec);

	//各粒子の各クラスタでの質量を修正
	for(int cIndx = 0; cIndx < pNum; cIndx++)
	{
		Vec3 orgCm = m_iceSM[cIndx]->GetOrgCm();

		for(int oIndx = 0; oIndx < m_iceSM[cIndx]->GetIndxNum(); oIndx++)
		{
			int pIndx = m_iceSM[cIndx]->GetParticleIndx(oIndx);
			if(pIndx == MAXINT){	continue;	}

			//ijiriら2012の手法
			float mass = 1.0f;
			Vec3 prtPos = m_iceSM[cIndx]->GetOrgPos(oIndx);
			Vec3 prtVec = Unit(prtPos-orgCm);
			float weight = mass * (dot(dirVec, prtVec) * dot(dirVec, prtVec) + 0.01f);

			m_iceSM[cIndx]->SetMass(oIndx, weight);
		}
	}

	//各クラスタの重心を更新　修正質量で重心を計算
	for(int pIndx = 0; pIndx < sm_particleNum; pIndx++){	m_iceSM[pIndx]->CalcOrgCm();	}
}

//初期位置を保存
void IceObject::SaveInitPos()
{
	for(int pIndx = 0; pIndx < sm_particleNum; pIndx++)
	{
		int sphIndx = pIndx*4;
		m_VecInitPos.push_back(Vec3(s_sphPrtPos[sphIndx+0], s_sphPrtPos[sphIndx+1], s_sphPrtPos[sphIndx+2]));
	}
}

//位置と速度を初期化する
void IceObject::ResetPosAndVel()
{
//粒子
	for(int pIndx = 0; pIndx < sm_particleNum; pIndx++)
	{
		int sphIndx = pIndx*4;

		s_sphPrtPos[sphIndx+0] = m_VecInitPos[pIndx][0];
		s_sphPrtPos[sphIndx+1] = m_VecInitPos[pIndx][1];
		s_sphPrtPos[sphIndx+2] = m_VecInitPos[pIndx][2];

		s_sphPrtVel[sphIndx+0] = 0.0f;
		s_sphPrtVel[sphIndx+1] = 0.0f;
		s_sphPrtVel[sphIndx+2] = 0.0f;
	}

//クラスタ
	for(int pIndx = 0; pIndx < sm_particleNum; pIndx++)
	{
		int smIndx = pIndx*3;
		float* p = Ice_SM::GetSldPosPointer();
		float* v = Ice_SM::GetSldVelPointer();
	
		v[smIndx+0] = 0.0f;
		v[smIndx+1] = 0.0f;
		v[smIndx+2] = 0.0f;
	
		p[smIndx+0] = m_VecInitPos[pIndx][0];
		p[smIndx+1] = m_VecInitPos[pIndx][1];
		p[smIndx+2] = m_VecInitPos[pIndx][2];
	}
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
	int max = m_iceSM[0]->GetNumVertices();
	int min = m_iceSM[0]->GetNumVertices();

	for(int i = 0; i < sm_clusterNum; i++)
	{
		int num =  m_iceSM[i]->GetNumVertices();
		
		if(max < num){	max = num;	}
		if(min > num){	min = num;	}

		sum += num;
	}

	cout << "sum = " << sum << " ave = " << (float)(sum/sm_clusterNum) << endl;
	cout << "max = " << max << " min = " << min << endl;
}

//近傍粒子
void IceObject::DebugNeights(const vector<vector<rxNeigh>>& neights)
{	cout << __FUNCTION__ << endl;
	
	string result = RESULT_DATA_PATH;
	result += "DebugNeights.txt";
	ofstream ofs(result, ios::app);

	for(unsigned pIndx = 0; pIndx < neights.size(); pIndx++){
		if(neights[pIndx].size() == 0){
			continue;
		}
		ofs << "Particle:" << pIndx << endl;
		for(unsigned id_np = 0; id_np < neights[pIndx].size(); id_np++){
			int np_pIndx = neights[pIndx][id_np].Idx;
			ofs << "	" << np_pIndx << endl;
		}
	}
}

//ある閾値を超えた粒子をカウント
void IceObject::DebugDeformationAmount()
{	//cout << __FUNCTION__ << endl;

	int deformClusterNum = 0;
	for(int i = 0; i < sm_clusterNum; ++i)
	{	
		if(m_iceSM[i]->GetNumVertices() == 0){	continue;	}

		float defAmount = m_iceSM[i]->GetDefAmount();
		if(defAmount < 1.0f){	continue;	}

		deformClusterNum++;
	}

	string result = RESULT_DATA_PATH;
	result += "DeformationClusterNum.txt";
	ofstream ofs(result, ios::app);
	ofs << deformClusterNum << endl;
}

//平均変形量
void IceObject::DebugDeformationAverage()
{
	float deformClusterAverage = 0.0f;
	for(int i = 0; i < sm_clusterNum; ++i)
	{	
		if(m_iceSM[i]->GetNumVertices() == 0){	continue;	}

		float defAmount = m_iceSM[i]->GetDefAmount();

		deformClusterAverage += defAmount;
	}

	string result = RESULT_DATA_PATH;
	result += "DeformationAverage.txt";
	ofstream ofs(result, ios::app);
	ofs << deformClusterAverage/(float)sm_clusterNum << endl;
}

//ディレクトリ内をすべて削除
void IceObject::DebugClearDirectry(string path)
{	cout << __FUNCTION__ << " path = " << path << endl;

	const boost::filesystem::path dir(path);
	boost::filesystem::directory_iterator end;

	//イテレータ内で削除の仕方がわからなかったので，全ファイルを取得してから削除している
	vector<boost::filesystem::path> files;

	//ディレクトリ内のファイル取得
	for(boost::filesystem::directory_iterator p(dir); p != end; ++p)
	{
		if(boost::filesystem::is_directory(*p))
		{  
			cout << "	" << p->path() << " is directry" << std::endl;	//ディレクトリの時
		}
		else
		{
			cout << "	" << p->path() << " is file" <<  endl;
		}

		files.push_back(p->path());
	}

	cout << endl;

	//ファイル削除
	vector<boost::filesystem::path>::iterator deleteFile;
	for(deleteFile = files.begin(); deleteFile != files.end(); deleteFile++)
	{
		if(boost::filesystem::is_directory(*deleteFile))
		{
			cout << "	" <<  *deleteFile << " はディレクトリです" << endl;
			continue;
		}

		if(boost::filesystem::remove(*deleteFile))
		{
			cout << "	" << *deleteFile << " を削除しました" << endl;
		}
		else
		{
			cout << "	" << "エラー " <<  *deleteFile << " を削除できません" << endl;
		}
	}
}

//運動計算していないクラスタの変形量を，前フレームとの差分から求める
//→結果が不完全だったので使わない
//→→結局，全てのクラスタを運動計算させてから，結果の算出に使う使わないを判断する，という方法にした
void IceObject::UpdateUnSelectedClusterDefAmount()
{
	float* sldPos = Ice_SM::GetSldPosPointer();

	for(int cIndx = 0; cIndx < sm_particleNum; cIndx++)
	{
		if(m_iceStrct->GetMotionCalcCluster(cIndx) != 0){	continue;	}
	
		float defAmount = 0.0f;
	
		//現在の重心を計算
		Vec3 nowCm(0.0);
		for(int i = 0; i < m_iceSM[cIndx]->GetIndxNum(); ++i)
		{
			int pIndx = m_iceSM[cIndx]->GetParticleIndx(i);			
			if(pIndx == MAXINT){	continue;	}

			nowCm += Vec3(sldPos[pIndx*3+0], sldPos[pIndx*3+1], sldPos[pIndx*3+2]);
		}

		nowCm /= (float)m_iceSM[cIndx]->GetNumVertices();

		//変形量を求める
		for(int i = 0; i < m_iceSM[cIndx]->GetIndxNum(); ++i)
		{
			int pIndx = m_iceSM[cIndx]->GetParticleIndx(i);			
			if(pIndx == MAXINT){	continue;	}
	
			//前フレームの相対座標
			Vec3 prePos = m_iceSM[cIndx]->GetPrePos(i) - m_iceSM[cIndx]->GetPreCm();

			//現在の相対座標
			Vec3 nowPos = Vec3(sldPos[pIndx*3+0], sldPos[pIndx*3+1], sldPos[pIndx*3+2]) - nowCm;

			defAmount += abs(nowPos[0]-prePos[0]);
			defAmount += abs(nowPos[1]-prePos[1]);
			defAmount += abs(nowPos[2]-prePos[2]);

			//前フレームの位置ベクトルを更新
			m_iceSM[cIndx]->SetPrePos(i, Vec3(sldPos[pIndx*3+0], sldPos[pIndx*3+1], sldPos[pIndx*3+2]));
		}
	
		m_iceSM[cIndx]->SetDefAmount(defAmount);
		m_iceSM[cIndx]->SetPreCm(nowCm);

		cout << cIndx << ":" << defAmount << endl;
	}
}

void IceObject::MakeDirectory(string pathName)
{	cout << __FUNCTION__ << " " << pathName << endl;

	const boost::filesystem::path path(pathName);

    boost::system::error_code error;
    const bool result = boost::filesystem::create_directory(path, error);
    if (!result || error) {
        std::cout << "ディレクトリの作成に失敗" << std::endl;
    }
}

//ファイル名に使える文字列を，現在時刻から作成
string IceObject::MakeDirNameTimeStamp()
{	
	FUNCTION_NAME;

	//現在の日時
	struct tm;
	unsigned int size = 26;
	char timeStr[80];
	time_t localtime;
	time(&localtime);
	ctime_s(timeStr, size, &localtime);
	
	//ディレクトリ名を作成
	string dirName(timeStr);
	replace(dirName.begin(), dirName.end(), ' ', '_');	//スペースをアンダースコアに置換 stringの関数でないのに注意
	replace(dirName.begin(), dirName.end(), ':', '_');	//コロンをアンダースコアに置換
	struct pred { bool operator()(char c) const { return /*c == ' ' ||*/ c == '\t' || c == '\r' || c == '\n'; } };
	dirName.erase(std::remove_if(dirName.begin(), dirName.end(), pred()), dirName.end());	//タブ，リターン，改行を削除

	return dirName;
}

void IceObject::DebugMakeGraph()
{	

	FUNCTION_NAME;

	//デバッグモードかの確認
	if(m_fpStepObjMove != &IceObject::StepObjMoveCPUDebug)	return;

	//現在時刻からディレクトリ名に用いる文字列を作成
	string dirName = MakeDirNameTimeStamp();
	dirName.insert(0, RESULT_IMAGE_PATH);			//ディレクトリを作成するパス

	//ディレクトリ作成
	MakeDirectory(dirName);

	//変形クラスタ数
	DebugMakeGraph_DefClusterNum(dirName);

	//平均変形量
	DebugMakeGraph_DefAmountAverage(dirName);

	//反復回数　TODO::きれいな方法じゃない…
	DebugMakeGraph_IterationNum(dirName);

	//現在のモードをテキストで保存
	ofstream ofs(dirName+"/nowMode.txt");

	string nowMode = DebugNowMoveMethod();
	ofs << dirName << endl;
	ofs << nowMode << endl;

////サンプル
//	vector<double> x,y;
//	x.resize(100);
//	y.resize(100);
//	//double x[100], y[100];
//
//	for (int i = 0; i<100; i++) {
//		x[i] = 0.0628*i;
//		y[i] = sin(x[i]);
//	}
//
//	vector<double> x2,y2;
//	x2.resize(100);
//	y2.resize(100);
//	//double x2[100], y2[100];
//
//	for (int i = 0; i<100; i++) {
//		x2[i] = 0.0628*i;
//		y2[i] = -sin(x[i]);
//	}
//
//	CGnuplot gp;
//	gp.SetTitle("Trigonometric Functions");
//	gp.Command("set key left bottom");
//	gp.SetLabel("x", "y");
//	gp.SetPlotType(CGnuplot::PLOT_TYPE_LINES_POINTS);
//	gp.SetXRange(-1.0, 1.0);
//	gp.SetYRange(-1.0, 1.0);
//
//	MultiPlot plotData;
//	plotData.push_back(svv("sin(x)", x, y));
//	plotData.push_back(svv("cos(x)", x2, y2));
//	gp.Multiplot(plotData);
//
//	gp.DumpToPng("result/images/test");
//	//うまくファイルを作れない　最後の行がおかしくなる
//	gp.DumpToFile("result/data/test.txt");

	////すぐに終了するとファイルが作られなかったのでちょっと待つ
	//int sleepTime = 3;	//秒
	//Sleep(sleepTime*1000);
}

//変形クラスタ数
//ファイルを読み込んでデータ作成
//TODO::いずれはファイルを直接渡してグラフを作る
void IceObject::DebugMakeGraph_DefClusterNum(string dirName)
{
	FUNCTION_NAME;

	vector<int> defNum_x, defNum_y;

	//データファイルを開く
	string filePath = RESULT_DATA_PATH;
	filePath += "DeformationClusterNum.txt";
	ifstream ifs(filePath);

	//ファイルの存在確認
	if(ifs.fail()) {	cerr << filePath << " is do not exist.\n";	return;	}

	//変数の用意，初期化
	int count = 0;
	string str;

	//文字列の読み込み
	while( getline(ifs, str) )
	{
		int data = 0;

		sscanf(str.data(), "%d", &data);

		defNum_x.push_back(count++);
		defNum_y.push_back(data);
	}

	ifs.close();

	//使用したデータファイルを移動
	string newFileName = dirName;
	newFileName += "/DeformationClusterNum.txt";
	const boost::filesystem::path path(filePath);
	const boost::filesystem::path dest(newFileName);

    try {
        boost::filesystem::rename(path, dest);
    }
    catch (boost::filesystem::filesystem_error& ex) {
        std::cout << ex.what() << std::endl;
        throw;
    }

	//グラフ作成
	CGnuplot gp;
	gp.SetTitle("Deformation Cluster Num Graph");
	gp.Command("set key left bottom");
	gp.SetLabel("TimeStep", "ClusterNum");
	gp.SetPlotType(CGnuplot::PLOT_TYPE_LINES_POINTS);
	gp.SetXRange(0, (int)defNum_x.size());
	gp.SetYRange(0, sm_particleNum);
	gp.Plot(defNum_x, defNum_y);

	string defNumResultPath = dirName;
	defNumResultPath += "/DeformationClusterNum";
	gp.DumpToPng(defNumResultPath.c_str());

	//すぐに終了するとファイルが作られなかったのでちょっと待つ
	int sleepTime = 2;	//秒
	Sleep(sleepTime*1000);
}

//平均変形量クラスタ数
void IceObject::DebugMakeGraph_DefAmountAverage(string dirName)
{
	FUNCTION_NAME;

	vector<float> defNum_x;
	vector<float> defNum_y;

	//ファイルを読み込み
	string filePath = RESULT_DATA_PATH;
	filePath += "DeformationAverage.txt";
	ifstream ifs(filePath);

	//ファイルの存在確認
	if(ifs.fail()) {	cerr << filePath << " is do not exist.\n";	return;	}

	//変数の用意，初期化
	int count = 0;
	float max = 0;
	string str;

	//文字列の読み込み
	while( getline(ifs, str) )
	{
		float data = 0;

		sscanf(str.data(), "%f", &data);

		if(data > max){	max = data;	}

		defNum_x.push_back(count++);
		defNum_y.push_back(data);
	}

	ifs.close();

	//使用したデータファイルを移動
	string newFileName = dirName;
	newFileName += "/DeformationAverage.txt";
	const boost::filesystem::path path(filePath);
	const boost::filesystem::path dest(newFileName);

    try {
        boost::filesystem::rename(path, dest);
    }
    catch (boost::filesystem::filesystem_error& ex) {
        std::cout << ex.what() << std::endl;
        throw;
    }

	//グラフ作成
	CGnuplot gp;
	gp.SetTitle("Deformation Average Graph");
	gp.Command("set key left bottom");
	gp.SetLabel("TimeStep", "Average Deformation Amount");
	gp.SetPlotType(CGnuplot::PLOT_TYPE_LINES_POINTS);
	gp.SetXRange(0, (int)defNum_x.size());
	//gp.SetYRange(0, max*2.0);	//縦軸は最大値の二倍
	gp.SetYRange(0, 0.5);			//比較しやすいように固定
	gp.Plot(defNum_x, defNum_y);

	string defNumResultPath = dirName;
	defNumResultPath += "/DefotmationAverage";
	gp.DumpToPng(defNumResultPath.c_str());

	//すぐに終了するとファイルが作られなかったのでちょっと待つ
	int sleepTime = 2;	//秒
	Sleep(sleepTime*1000);
}

void IceObject::DebugMakeGraph_IterationNum(string dirName)
{
	//モード確認
	//TODO: あんまり綺麗なやり方じゃない
	const char* calcMethodName1 = typeid(Ice_CalcMethod_Itr_Stiffness).name();
	const char* calcMethodName2 = typeid(Ice_CalcMethod_Itr_Exp_Stiff).name();
	const char* nowName = typeid(*(m_iceCalcMethod)).name();

	//現在iteration_stiffnessモードかの確認　strcmpは文字列が一致すると0になる
	if(strcmp(calcMethodName1, nowName) != 0
	&& strcmp(calcMethodName2, nowName) != 0) return ;

	FUNCTION_NAME;

	vector<float> defNum_x;
	vector<float> defNum_y;

	//ファイルを読み込み
	string filePath = RESULT_DATA_PATH;
	filePath += "Itr_Stiffness_CountNum.txt";
	ifstream ifs(filePath);

	//ファイルの存在確認
	if(ifs.fail()) {	cerr << filePath << " is do not exist.\n";	return;	}

	//変数の用意，初期化
	int count = 0;
	float max = 0;
	string str;

	//文字列の読み込み
	while( getline(ifs, str) )
	{
		float data = 0;

		sscanf(str.data(), "%f", &data);

		if(data > max){	max = data;	}

		defNum_x.push_back(count++);
		defNum_y.push_back(data);
	}

	ifs.close();

	//使用したデータファイルを移動
	string newFileName = dirName;
	newFileName += "/Itr_Stiffness_CountNum.txt";
	const boost::filesystem::path path(filePath);
	const boost::filesystem::path dest(newFileName);

    try {
        boost::filesystem::rename(path, dest);
    }
    catch (boost::filesystem::filesystem_error& ex) {
        std::cout << ex.what() << std::endl;
        throw;
    }

	//グラフ作成
	CGnuplot gp;
	gp.SetTitle("Iteration Num Graph");
	gp.Command("set key left bottom");
	gp.SetLabel("TimeStep", "Iteration Num");
	gp.SetPlotType(CGnuplot::PLOT_TYPE_LINES_POINTS);
	gp.SetXRange(0, (int)defNum_x.size());
	//gp.SetYRange(0, max*2.0);	//縦軸は最大値の二倍
	gp.SetYRange(0, 250);			//比較しやすいように固定
	gp.Plot(defNum_x, defNum_y);

	string defNumResultPath = dirName;
	defNumResultPath += "/Itr_Stiffness_CountNum";
	gp.DumpToPng(defNumResultPath.c_str());

	//すぐに終了するとファイルが作られなかったのでちょっと待つ
	int sleepTime = 2;	//秒
	Sleep(sleepTime*1000);
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

void IceObject::DebugUpdateSelectCluster()
{	cout << __FUNCTION__ << endl;
	for(int cIndx = 0; cIndx < sm_particleNum; cIndx++)
	{
		if(m_iceStrct->GetMotionCalcCluster(cIndx) == 0){	continue;	}
		cout << "cIndx = " << cIndx << " , " << m_iceStrct->GetMotionCalcCluster(cIndx) << endl;
	}
}

//各手法のクラス名を表示
string IceObject::DebugNowMoveMethod()
{	cout << __FUNCTION__ << endl;

	string nowMode;
	string end = "\n";
	nowMode = string("CalcMethod	") +	typeid(*m_iceCalcMethod).name() +	end		//計算方法
			+ string("JudgeMove	") +		typeid(*m_iceJudeMove).name() +		end		//運動計算対象を判定
			+ string("ClsuterMove	") +	typeid(*m_iceClsuterMove).name() +	end		//運動計算方法
			+ string("ConvopoJudge	") +	typeid(*m_iceConvoJudge).name() +	end		//最終統合結果に用いる対象を判定
			+ string("Convolution	" ) +	typeid(*m_iceConvolution).name() +	end;	//最終統合結果を求める

	cout << nowMode << endl;

	return nowMode;
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
		m_SurfSm->GetDevicePRTtoPTHPointer(),
		m_SurfSm->GetDevicePTHandPrfxSetPointer(),
		m_SurfSm->GetDecvicePrfxPos(),
		m_SurfSm->GetDecvicePrfxApq(),
		//IceStruct------------------------------------------
		m_iceStrct->GetDeviceCtoPointer(),
		m_iceStrct->GetDeviceCtoPNumPointer(),
		m_iceStrct->GetCtoPMax(),
		2,
		0.01		//TODO: 0.02では？いろいろばらばらっぽい
	);

//CPU側
	//prefixSumの更新
	m_SurfSm->UpdatePrefixSum();

	//クラスタのパラメータ更新
	#pragma omp parallel for
	for(int i = 0; i < sm_clusterNum; ++i)
	{	
		if(GetPtoCNum(i) == 0){	continue;	}
		m_iceSM[i]->SetNowCm(m_SurfSm->CalcCmSum(i));		//重心の更新
		m_iceSM[i]->SetApq(m_SurfSm->CalcApqSum(i));		//変形行列の更新
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

	//	Vec3 cm = m_iceSM[iCrstr]->GetCm();
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

	//	rxMatrix3 apq = m_iceSM[iCrstr]->GetApq();
	//	cout << "CPU:[" << iCrstr << "] = " << endl;
	//	cout <<	"(" << apq(0, 0) << ", " << apq(0, 1) << ", " << apq(0, 2) << ")" << endl;
	//	cout <<	"(" << apq(1, 0) << ", " << apq(1, 1) << ", " << apq(1, 2) << ")" << endl;
	//	cout <<	"(" << apq(2, 0) << ", " << apq(2, 1) << ", " << apq(2, 2) << ")" << endl;
	//	cout << endl;
	//}

	//delete[] dCurApq;

}

//固体の上部を固定
void IceObject::TestFixUpperPos()
{
	int side = 17;

	//粒子
	for(int pIndx = 0; pIndx < sm_particleNum; pIndx++)
	{
		if(pIndx%(side*side) < side * (side-1)) continue;
		//if( pIndx%(side*side) % side >= 10) continue;

		int sphIndx = pIndx*4;

		s_sphPrtPos[sphIndx+0] = m_VecInitPos[pIndx][0];
		s_sphPrtPos[sphIndx+1] = m_VecInitPos[pIndx][1];
		s_sphPrtPos[sphIndx+2] = m_VecInitPos[pIndx][2];

		s_sphPrtVel[sphIndx+0] = 0.0f;
		s_sphPrtVel[sphIndx+1] = 0.0f;
		s_sphPrtVel[sphIndx+2] = 0.0f;
	}

//クラスタ
	for(int pIndx = 0; pIndx < sm_particleNum; pIndx++)
	{
		if(pIndx%(side*side) < side * (side-1)) continue;
		//if( pIndx%(side*side) % side >= 10) continue;

		int smIndx = pIndx*3;
		float* p = Ice_SM::GetSldPosPointer();
		float* v = Ice_SM::GetSldVelPointer();
	
		v[smIndx+0] = 0.0f;
		v[smIndx+1] = 0.0f;
		v[smIndx+2] = 0.0f;
	
		p[smIndx+0] = m_VecInitPos[pIndx][0];
		p[smIndx+1] = m_VecInitPos[pIndx][1];
		p[smIndx+2] = m_VecInitPos[pIndx][2];
	}
}

//固体の側面を固定
void IceObject::TestFixSidePos()
{
	int side = 17;

	//粒子
	for(int pIndx = 0; pIndx < sm_particleNum/side; pIndx++)
	{
		int sphIndx = pIndx*4;

		s_sphPrtPos[sphIndx+0] = m_VecInitPos[pIndx][0];
		s_sphPrtPos[sphIndx+1] = m_VecInitPos[pIndx][1];
		s_sphPrtPos[sphIndx+2] = m_VecInitPos[pIndx][2];

		s_sphPrtVel[sphIndx+0] = 0.0f;
		s_sphPrtVel[sphIndx+1] = 0.0f;
		s_sphPrtVel[sphIndx+2] = 0.0f;
	}

//クラスタ
	for(int pIndx = 0; pIndx < sm_particleNum/side; pIndx++)
	{
		int smIndx = pIndx*3;
		float* p = Ice_SM::GetSldPosPointer();
		float* v = Ice_SM::GetSldVelPointer();
	
		v[smIndx+0] = 0.0f;
		v[smIndx+1] = 0.0f;
		v[smIndx+2] = 0.0f;
	
		p[smIndx+0] = m_VecInitPos[pIndx][0];
		p[smIndx+1] = m_VecInitPos[pIndx][1];
		p[smIndx+2] = m_VecInitPos[pIndx][2];
	}
}

//ファイルを読み込んだパラメータでシミュレーションを実行
void IceObject::TestSimulationFromFile(string fileName)
{
	cout << __FUNCTION__ << " " << fileName << endl;

	//ファイルオープン
	//データ読み込み
	//データ反映
}