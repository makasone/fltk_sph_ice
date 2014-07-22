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

IceObject::IceObject()
{
}

IceObject::~IceObject()
{
}

//それぞれのクラス・変数の初期化
void IceObject::InitIceObj(int pMaxNum, int cMaxNum, int tMaxNum)
{	cout << __FUNCTION__ << endl;
	//物体の構造の初期化
	m_iceStrct = new IceStructure(pMaxNum, cMaxNum, tMaxNum);

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

//クラスタの初期化	//一時的な実装
void IceObject::InitCluster(Ice_SM* sm)
{
	m_iceMove.push_back(sm);	//ポインタをコピーしているだけ
}

//クラスタの初期化
void IceObject::InitCluster(Vec3 boundarySpaceHigh, Vec3 boundarySpaceLow, float timeStep)
{
	//以下はきちんと移行した場合のコード
	//変数の初期化
	for(vector<Ice_SM*>::iterator it = m_iceMove.begin(); it != m_iceMove.end(); ++it)
	{
		if(*it) delete *it;
	}
	
	sm_clusterNum = 0;	//クラスタ数の初期化

	Ice_SM::SetParticlePosAndVel(s_sphPrtPos, s_sphPrtVel);	//TOD::ポインタをコピーするのではなく上からプッシュしたほうがいい

	////階層構造による剛性の確保
	////実験１
	////MakeClusterHigh();

	//四面体リストを元に，粒子毎にクラスタ作成
	for(int i = 0; i < sm_particleNum; ++i)
	{
		//クラスタ初期化
		m_iceMove.push_back(new Ice_SM(sm_clusterNum));
		m_iceMove[sm_clusterNum]->SetSimulationSpace(boundarySpaceLow, boundarySpaceHigh);
		m_iceMove[sm_clusterNum]->SetTimeStep(timeStep);
		m_iceMove[sm_clusterNum]->SetCollisionFunc(0);
		m_iceMove[sm_clusterNum]->SetStiffness(1.0, 1.0);

		//四面体リストを元に，粒子毎にクラスタ作成
		SetClusterMoveInfo(i);

		sm_clusterNum++;
	}

	////MakeClusterFromNeight();
	////MakeOneCluster();

	//クラスタに関するGPUの初期化
	//すべてのクラスタにアクセスできるようにポインタを渡している
	Ice_SM::InitGPU(m_iceMove, sd_sphPrtPos, sd_sphPrtVel, sm_particleNum);
	Ice_SM::InitFinalParamPointer(sm_clusterNum);

	//クラスタのデータをGPUへ初期化
	for(int i = 0; i < sm_particleNum; i++)
	{
		m_iceMove[i]->InitGPU_Instance();
	}

	//表面粒子に重みを付ける
	//((RXSPH*)m_pPS)->DetectSurfaceParticles();								//表面粒子検出
	//int *surfaceParticles = (int *)( ((RXSPH*)m_pPS)->GetArraySurf() );		//表面粒子
	//float weightedMass = 300.0f;

	////表面粒子の質量だけ
	//for(int pIndx = 0; pIndx < ICENUM; pIndx++)
	//{
	//	if(surfaceParticles[pIndx] != 1){ continue;	}
	//	else{	cout << "surfaceParticles[" << pIndx << "] = " << surfaceParticles[pIndx] << endl;	}

	//	for(int pNum = 0; pNum < m_ice->GetPtoCIndx(pIndx); ++pNum)
	//	{
	//		int cIndx = m_ice->GetPtoT(pIndx, pNum, 0);
	//		int oIndx = m_ice->GetPtoT(pIndx, pNum, 1);

	//		if(cIndx == -1 || oIndx == -1){ continue;	}

	//		m_sm_cluster[cIndx]->SetMass(oIndx, weightedMass);
	//	}
	//}

	//表面粒子のクラスタすべての粒子
	//for(int cIndx = 0; cIndx < ICENUM; ++cIndx)
	//{
	//	if(surfaceParticles[cIndx] != 1){ continue;	}
	//	else{	cout << "surfaceParticles[" << cIndx << "] = " << surfaceParticles[cIndx] << endl;	}

	//	for(int pIndx = 0; pIndx < /*m_ice->GetCtoPIndx(cIndx)*/m_sm_cluster[cIndx]->GetNumVertices(); ++pIndx)
	//	{
	//		m_sm_cluster[cIndx]->SetMass(pIndx, weightedMass);
	//	}

	//	//	for(int pNum = 0; pNum < m_ice->GetPtoCIndx(pIndx); ++pNum)
	//	//	{
	//	//		int cIndx = m_ice->GetPtoT(pIndx, pNum, 0);
	//	//		int oIndx = m_ice->GetPtoT(pIndx, pNum, 1);

	//	//		if(cIndx == -1 || oIndx == -1){ continue;	}

	//	//		m_sm_cluster[cIndx]->SetMass(oIndx, weightedMass);
	//	//	}
	//}

	//全体の半分の粒子
	//なぜかうまくいかない？
	//for(int cIndx = 0; cIndx < ICENUM; cIndx++)
	//{
	//	//for(int pNum = 0; pNum < m_sm_cluster[cIndx]->GetNumVertices(); pNum++)
	//	//{
	//	//	int pIndx = m_sm_cluster[cIndx]->GetParticleIndx(pNum);
	//	//	
	//	//	if(surfaceParticles[pIndx] != 1){ continue;	}
	//	//	else{	cout << "surfaceParticles[" << pIndx << "] = " << surfaceParticles[pIndx] << endl;	}
	//	//
	//	//	m_sm_cluster[cIndx]->SetMass(pNum, weightedMass);
	//	//}
	//	
	//	for(int pNum = 0; pNum < m_sm_cluster[cIndx]->GetNumVertices(); pNum++)
	//	{
	//		int pIndx = m_sm_cluster[cIndx]->GetParticleIndx(pNum);
	//		
	//		//if(surfaceParticles[pIndx] != 1){ continue;	}
	//		//else{	cout << "surfaceParticles[" << pIndx << "] = " << surfaceParticles[pIndx] << endl;	}
	//	
	//		if(pIndx > ICENUM/2 || pIndx < 0){ continue;	}
	//		else{	/*cout << "pInx = " << pIndx << endl;*/	}

	//		m_sm_cluster[cIndx]->SetMass(pNum, weightedMass);
	//	}
	//}
	

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

	DebugClusterInfo();
}

//オブジェクトの構造の初期化
void IceObject::InitStrct()
{
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
	cout << __FUNCTION__ << " check0" << endl;

	//メモリ確保
	InitClusterInfo();
	cout << __FUNCTION__ << " check1" << endl;

	//粒子が所属しているクラスタ数の配列をコピー
	int *PtoCNum = new int[sm_particleNum];
	cout << __FUNCTION__ << " check2" << endl;
	for(int i = 0; i < sm_particleNum; i++)
	{
		PtoCNum[i] = GetPtoCNum(i);
	}
	cout << __FUNCTION__ << " check3" << endl;

	//クラスタと粒子の関連情報の登録
	for(int i = 0; i < sm_particleNum; i++)
	{
		SetClusterStrctInfo(i, PtoCNum);	//カウントを前で行っているため，こちらを使う
	}
	cout << __FUNCTION__ << " check4" << endl;

	delete[] PtoCNum;

	//TODO::クラスタと四面体の関連情報の登録
	cout << __FUNCTION__ << " check5" << endl;
	//m_ice->InitPath(p, v, m_sm_cluster, ICENUM);			//高速化のためのパス作成


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


//クラスタの運動計算情報の登録
void IceObject::SetClusterMoveInfo(int pIndx)
{
	//((RXSPH*)m_pPS)->DetectSurfaceParticles();								//表面粒子検出
	//int *surfaceParticles = (int *)( ((RXSPH*)m_pPS)->GetArraySurf() );		//表面粒子
	//float weightedMass = 300.0f;
	//int layer = 1;

	//ある粒子が含まれる四面体を近傍クラスタとし，各四面体に含まれる粒子でクラスタを作成
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

			//if(jpIndx < ICENUM/2)					//粒子の半分に重み付け
			//if(surfaceParticles[jpIndx] == 1)
			//if( (jpIndx < 13*13*layer) || (ICENUM-13*13*layer < jpIndx) )
			//if( ()									//x層
			//||										//y層
			//)
			//{
			//	//cout << "jpIndx == " << jpIndx << endl;
			//	mass = weightedMass;
			//}

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

				//if(kpIndx < ICENUM/2)					//粒子の半分に重み付け
				//if(surfaceParticles[kpIndx] == 1)
				//if( (kpIndx < 13*13*layer) || (ICENUM-13*13*layer < kpIndx) ) 
				//{
				//	//cout << "kpIndx == " << kpIndx << endl;
				//	mass = weightedMass;
				//}

				m_iceMove[pIndx]->AddVertex( Vec3(s_sphPrtPos[kpIndx*4+0], s_sphPrtPos[kpIndx*4+1], s_sphPrtPos[kpIndx*4+2] ), mass, kpIndx);
				m_iceMove[pIndx]->SetAlphas(pNum, 1.0);
				m_iceMove[pIndx]->SetBetas (pNum, 0.0);
				m_iceMove[pIndx]->SetLayer (pNum, jlIndx);
				//cout << "pIndx = " << pIndx << " vNum = " << m_sm_cluster[pIndx]->GetNumVertices() << endl;
			}
		}
	}
}

//
void IceObject::SetClusterStrctInfo(int cIndx, int* PtoCNum)
{
	//粒子が属している四面体の番号を登録するための準備
	//pCountListには，cIndx番目のクラスタに含まれる各粒子が，それぞれいくつのクラスタに属するかを求めて保存する
	int* pCountList = new int[m_iceMove[cIndx]->GetNumVertices()];
	int* pLayerList = new int[m_iceMove[cIndx]->GetNumVertices()];			//粒子の所属レイヤー
	
	for(int i = 0; i < m_iceMove[cIndx]->GetNumVertices(); i++)
	{
		int pIndx = m_iceMove[cIndx]->GetParticleIndx(i);

		//穴あきを想定していないため，誤って上書きしてしまっている
		//-1を探索して上書きするのに切り替える．
		for(int j = 0; j < GetPtoCMax(); j++)
		{
			if(GetPtoC(pIndx, j, 0) != -1 || GetPtoC(pIndx, j, 1) != -1){	continue;	}
			
			if(j >= GetPtoCIndx(pIndx))
			{
				SetPtoCIndx(pIndx, j+1);		//現在のIndxより大きいなら更新
			}

			pCountList[i] = j;
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
		pIndxList.push_back(pIndx);
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

	//TODO::クラスタのGPU初期化もここに置く

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
void IceObject::StepObjMove()
{

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
	for(int i = 0; i < sm_particleNum; ++i)
	{
		if(GetParticleNum() <= i){	continue;	}	//融解のみの実験のときに必要になる．
		if(GetPtoCNum(i) <= 0){		continue;	}
	
		Vec3 pos,vel;
	
		//固体運動の最終位置計算
		CalcAverageCPU(i, pos, vel);

		//液体と固体の補間
		LinerInterPolationCPU(i, pos, vel);	//線形補間
	}
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
#else
	for(int i = 0; i < sm_particleNum; ++i)
	{
		if(GetParticleNum() <= i){	continue;	}	//融解のみの実験のときに必要になる．
		if(GetPtoCNum(i) <= 0){		continue;	}
	
		Vec3 pos,vel;
	
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

		for(int j = 0; j < GetPtoCIndx(pIndx); ++j)
		{
			int jcIndx = GetPtoC(pIndx, j, 0);
			int joIndx = GetPtoC(pIndx, j, 1);

			if(jcIndx == -1 || joIndx == -1){	continue;	}

			pos += m_iceMove[jcIndx]->GetVertexPos(joIndx);
			vel += m_iceMove[jcIndx]->GetVertexVel(joIndx);

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

//デバッグ

//四面体
void IceObject::DebugTetraInfo()
{
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
{
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