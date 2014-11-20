//-----------------------------------------------------------------------------
// インクルードファイル
//-----------------------------------------------------------------------------
#include "IceStructure.h"

//粒子→
int* IceStructure::sd_piPtoC;
int* IceStructure::sd_piPtoT;

int* IceStructure::sd_piPtoCNum;
int* IceStructure::sd_piPtoTNum;

int* IceStructure::sd_piPtoCIndx;
int* IceStructure::sd_piPtoTIndx;

//クラスタ→
int* IceStructure::sd_piCtoP;

int* IceStructure::sd_piCtoPNum;

int* IceStructure::sd_piCtoPIndx;

//四面体→
int* IceStructure::sd_piTtoP;

int* IceStructure::sd_piTtoPNum;
int* IceStructure::sd_piTtoCNum;

int* IceStructure::sd_piTtoPIndx;
int* IceStructure::sd_piTtoCIndx;

int* IceStructure::sd_piNeighborTetra;

int* IceStructure::sd_piNeighborTetraTNum;

/*!
 * @param[in] pNumMax　最大粒子数
 * @param[in] cNumMax　最大クラスタ数
 * @param[in] tNumMax　最大四面体数
 */
IceStructure::IceStructure()
{
}

IceStructure::IceStructure(int pNumMax, int cNumMax, int tNumMax, int layer)
{	cout << __FUNCTION__ << endl;

	//最大数の登録
	m_iPNumMax = pNumMax;
	m_iCNumMax = cNumMax;
	m_iTNumMax = tNumMax;

	//粒子情報の初期化
	m_iPtoCMax = m_iCNumMax*0.1f;	//1331 layer2 0.4 layer3 0.75
									//2197 layer2 0.4 layer3 0.4 layer4 0.5
									//4913 layer2 0.3
									//4913 layer1 0.1
									//12167 layer1 0.01
									//19683 layer1 0.01

	m_iPtoTMax = m_iTNumMax*0.1f;	//1331 layer2 0.3 layer3 0.5
									//2197 layer2 0.3 layer3 0.3
									//4913 layer2 0.3
									//12167 layer1 0.01
									//19683 layer1 0.01

cout << __FUNCTION__ << " check1" << endl;

	m_piPtoCNum = new int[m_iPNumMax];
	m_piPtoTNum = new int[m_iPNumMax];

	m_piPtoCIndx = new int[m_iPNumMax];
	m_piPtoTIndx = new int[m_iPNumMax];

	for(int i = 0; i < m_iPNumMax; i++)
	{
		m_piPtoCNum[i] = 0;
		m_piPtoTNum[i] = 0;
		
		m_piPtoCIndx[i] = 0;
		m_piPtoTIndx[i] = 0;
	}

	//クラスタ情報の初期化
	//CtoTMaxは粒子数と等しいので定義しない
	//高速化をテストするときはCtoTはコメントに
	m_iCtoPMax = m_iPNumMax*0.1f;	//1331 layer2 0.5 layer3 0.75
									//2197 layer2 0.5 layre3 0.5
									//4913 layer2
									//4913 layer1 0.1
									//12167 layer1 0.01
									//19683 layer1 0.01
	//m_iCtoPMax = m_iPNumMax;		//単一クラスタ

	m_piCtoPNum = new int[m_iCNumMax];
	m_piCtoPIndx = new int[m_iCNumMax];

cout << __FUNCTION__ << " check2" << endl;

	for(int i = 0; i < m_iCNumMax; i++)
	{
		m_piCtoPNum[i] = 0;
		//m_piCtoTNum[i] = 0;
		
		m_piCtoPIndx[i] = 0;
		//m_piCtoTIndx[i] = 0;
	}

	//四面体情報の初期化
	//TtoPMaxは最大４で固定
	//TtoCMaxは必要ない
	m_piTtoPNum = new int[m_iTNumMax];
	m_piTtoCNum = new int[m_iTNumMax];

	m_piTtoPIndx = new int[m_iTNumMax];
	m_piTtoCIndx = new int[m_iTNumMax];

cout << __FUNCTION__ << " check3" << endl;

	for(int i = 0; i < m_iTNumMax; i++)
	{
		m_piTtoPNum[i] = 0;
		m_piTtoCNum[i] = 0;
		
		m_piTtoPIndx[i] = 0;
		m_piTtoCIndx[i] = 0;
	}

	//近傍四面体
	int ntnSize = m_iTNumMax * 1.0f;
	m_piNTNum = new int[ntnSize];

cout << __FUNCTION__ << " check3.5" << endl;

	//m_iNeighborMax = m_iTNumMax*0.1f;		//1331  layer2 0.3 layer3 0.75
											//2197  layer2 0.3 layre3 0.3 layer4 0.4
											//4913  layer2 0.1
											//6859  layer1 0.01
											//9261  layer1 0.005
											//12167 layer1

	m_iNeighborMax = 150;					//layer1なら大丈夫

	cout << "m_iTNumMax = " << m_iTNumMax << ", m_iNeighborMax = " << m_iNeighborMax << endl;

	m_mk3DiNeighborTetra.SetSize(m_iTNumMax, m_iNeighborMax, 2);

cout << __FUNCTION__ << " check4" << endl;

	for(int i = 0; i < m_iTNumMax; i++)
	{
		m_piNTNum[i] = 0;

		for(int j = 0; j < m_iNeighborMax; j++)
		{
			for(int k = 0; k < 2; k++)
			{
				m_mk3DiNeighborTetra(i, j, k) = -1;
			}
		}
	}

	m_iLayer = layer;

	//フラグ
	//m_pbPFlag = new bool[m_iPNumMax];
	//m_pbCFlag = new bool[m_iCNumMax];
	m_pbTFlag = new bool[m_iTNumMax];

	//ResetPFlag(m_iPNumMax);
	//ResetCFlag(m_iCNumMax);
	ResetTFlag(m_iTNumMax);

	//選択的運動計算
	m_psuSelectClusterIndx = new short unsigned[m_iPNumMax];
	
	for(int i = 0; i < m_iPNumMax; i++)
	{
		UpdateMotionCalcCluster(i, 1);
	}

cout << __FUNCTION__ << " check5" << endl;
}

IceStructure::~IceStructure(void)
{
}

//------------------------------------------フラグ管理-------------------------------------------------
/*!
 * 探索用粒子フラグの初期化
 */
void IceStructure::ResetPFlag(int endNum)
{
	for(int i = 0; i < endNum; i++)
	{
		m_pbPFlag[i] = false;
	}
}

/*!
 * 探索用クラスタフラグの初期化
 */
void IceStructure::ResetCFlag(int endNum)
{
	for(int i = 0; i < endNum; i++)
	{
		m_pbCFlag[i] = false;
	}
}

/*!
 * 探索用四面体フラグの初期化
 */
void IceStructure::ResetTFlag(int endNum)
{
	for(int i = 0; i < endNum; i++)
	{
		m_pbTFlag[i] = false;
	}
}
//------------------------------------------フラグ管理-------------------------------------------------

//------------------------------------------初期化-------------------------------------------------
/*!
 * 計算情報　四面体情報の領域確保　粒子ベース処理
 */
void IceStructure::InitTetraInfo()
{
	//粒子→四面体
	m_mk3DiPtoT.SetSize(m_iPNumMax, m_iPtoTMax, 2);

	for(int i = 0; i < m_iPNumMax; i++)
	{		
		for(int j = 0; j < m_iPtoTMax; j++)
		{
			for(int k = 0; k < 2; k++)
			{
				m_mk3DiPtoT(i, j, k) = -1;
			}
		}
	}

	//四面体→粒子
	m_mk2DiTtoP.SetSize(m_iTNumMax, 4);
	
	for(int i = 0; i < m_iTNumMax; i++)
	{
		for(int j = 0; j < 4; j++)
		{
			m_mk2DiTtoP(i, j) = -1;
		}
	}

	//配列の添え字の初期化
	for(int i = 0; i < m_iPNumMax; i++)
	{
		m_piPtoTIndx[i] = m_piPtoTNum[i];
	}

	for(int i = 0; i < m_iTNumMax; i++)
	{
		m_piTtoPIndx[i] = m_piTtoPNum[i];
	}
}

/*!
 * 計算情報　クラスタ情報の領域確保　粒子ベース処理
 */
void IceStructure::InitClusterInfo()
{	cout << __FUNCTION__ << endl;

	//粒子→クラスタ
	m_mk3DiPtoC.SetSize(m_iPNumMax, m_iPtoCMax, 3);
	cout << __FUNCTION__ << ", check1" << endl;

	for(int i = 0; i < m_iPNumMax; i++)
	{
		for(int j = 0; j < m_iPtoCMax; j++)
		{
			//[0] = クラスタ番号，[1] = クラスタ内での番号，[2] = layer番号
			for(int k = 0; k < 3; k++)
			{
				m_mk3DiPtoC(i, j, k) = -1;
			}
		}
	}

	//クラスタ→粒子
	m_mk3DiCtoP.SetSize(m_iCNumMax, m_iCtoPMax, 2);
	cout << __FUNCTION__ << ", check2" << endl;

	for(int i = 0; i < m_iCNumMax; i++)
	{
		for(int j = 0; j < m_iCtoPMax; j++)
		{
			//[0] = 粒子番号，[1] = layer番号
			for(int k = 0; k < 2; k++)
			{
				m_mk3DiCtoP(i, j, k) = -1;
			}
		}
	}

	//配列の添え字Indxの初期化
	for(int i = 0; i < m_iPNumMax; i++)
	{
		m_piPtoCIndx[i] = m_piPtoCNum[i];
	}

	for(int i = 0; i < m_iCNumMax; i++)
	{
		m_piCtoPIndx[i] = m_piCtoPNum[i];
	}

	cout << __FUNCTION__ << ", check3" << endl;
}

//GPU処理で用いるデータの初期化
//今のところはコピーでまかなう
void IceStructure::InitGPU()
{	cout << __FUNCTION__ << endl;

	//デバイス側のメモリを確保
		//粒子→
	cudaMalloc((void**)&sd_piPtoT,		sizeof(int) * m_iPNumMax * m_iPtoTMax * 2);
	cudaMalloc((void**)&sd_piPtoC,		sizeof(int) * m_iPNumMax * m_iPtoCMax * 3);

	cudaMalloc((void**)&sd_piPtoCNum,	sizeof(int) * m_iPNumMax);
	cudaMalloc((void**)&sd_piPtoTNum,	sizeof(int) * m_iPNumMax);

	cudaMalloc((void**)&sd_piPtoCIndx,	sizeof(int) * m_iPNumMax);
	cudaMalloc((void**)&sd_piPtoTIndx,	sizeof(int) * m_iPNumMax);

		//クラスタ→
	cudaMalloc((void**)&sd_piCtoP,		sizeof(int) * m_iCNumMax * m_iCtoPMax * 2);

	cudaMalloc((void**)&sd_piCtoPNum,	sizeof(int) * m_iCNumMax);
	cudaMalloc((void**)&sd_piCtoPIndx,	sizeof(int) * m_iCNumMax);

		//四面体→
	cudaMalloc((void**)&sd_piTtoP,	sizeof(int) * m_iTNumMax * 4);

	cudaMalloc((void**)&sd_piTtoPNum,	sizeof(int) * m_iTNumMax);
	cudaMalloc((void**)&sd_piTtoCNum,	sizeof(int) * m_iTNumMax);

	cudaMalloc((void**)&sd_piTtoPIndx,	sizeof(int) * m_iTNumMax);
	cudaMalloc((void**)&sd_piTtoCIndx,	sizeof(int) * m_iTNumMax);

		//近傍四面体
	cudaMalloc((void**)&sd_piNeighborTetra,		sizeof(int) * m_iTNumMax * m_iNeighborMax * 2);
	cudaMalloc((void**)&sd_piNeighborTetraTNum,	sizeof(int) * m_iTNumMax);

	//初期化
	cudaMemcpy(sd_piPtoCNum, m_piPtoCNum, sizeof(int) * m_iPNumMax, cudaMemcpyHostToDevice);
	cudaMemcpy(sd_piPtoTNum, m_piPtoTNum, sizeof(int) * m_iPNumMax, cudaMemcpyHostToDevice);

	cudaMemcpy(sd_piPtoCIndx, m_piPtoCIndx, sizeof(int) * m_iPNumMax, cudaMemcpyHostToDevice);
	cudaMemcpy(sd_piPtoTIndx, m_piPtoTIndx, sizeof(int) * m_iPNumMax, cudaMemcpyHostToDevice);

	cudaMemcpy(sd_piCtoPNum,	m_piCtoPNum,	sizeof(int) * m_iCNumMax, cudaMemcpyHostToDevice);
	cudaMemcpy(sd_piCtoPIndx,	m_piCtoPIndx,	sizeof(int) * m_iCNumMax, cudaMemcpyHostToDevice);

	cudaMemcpy(sd_piTtoPNum,	m_piTtoPNum,	sizeof(int) * m_iTNumMax, cudaMemcpyHostToDevice);
	cudaMemcpy(sd_piTtoCNum,	m_piTtoCNum,	sizeof(int) * m_iTNumMax, cudaMemcpyHostToDevice);

	cudaMemcpy(sd_piTtoPIndx,	m_piTtoPIndx,	sizeof(int) * m_iTNumMax, cudaMemcpyHostToDevice);
	cudaMemcpy(sd_piTtoCIndx,	m_piTtoCIndx,	sizeof(int) * m_iTNumMax, cudaMemcpyHostToDevice);

		//Vectorで管理しているデータを配列にするためにdata()を使っている
	cudaMemcpy(sd_piNeighborTetra,		m_mk3DiNeighborTetra.Get().data(),	sizeof(int) * m_iTNumMax * m_iNeighborMax * 2,	cudaMemcpyHostToDevice);
	cudaMemcpy(sd_piNeighborTetraTNum,	m_piNTNum,							sizeof(int) * m_iTNumMax,						cudaMemcpyHostToDevice);

	cudaMemcpy(sd_piPtoT,	m_mk3DiPtoT.Get().data(),	sizeof(int) * m_iPNumMax * m_iPtoTMax * 2,	cudaMemcpyHostToDevice);
	cudaMemcpy(sd_piPtoC,	m_mk3DiPtoC.Get().data(),	sizeof(int) * m_iPNumMax * m_iPtoCMax * 3,	cudaMemcpyHostToDevice);

	cudaMemcpy(sd_piTtoP,	m_mk2DiTtoP.Get().data(),	sizeof(int) * m_iTNumMax * 4,				cudaMemcpyHostToDevice);

	cudaMemcpy(sd_piCtoP,	m_mk3DiCtoP.Get().data(),	sizeof(int) * m_iCNumMax * m_iCtoPMax * 2,	cudaMemcpyHostToDevice);

////デバッグ
//	int* testA = new int[m_iPNumMax];
//	int* testB = new int[m_iPNumMax];
//	int* testC = new int[m_iCNumMax];
//	int* testD = new int[m_iCNumMax];
//
//	//デバイスからホストへのコピー
//	cudaMemcpy(testA, sd_piPtoCIndx, sizeof(int) * m_iPNumMax, cudaMemcpyDeviceToHost);
//	cudaMemcpy(testB, sd_piPtoTIndx, sizeof(int) * m_iPNumMax, cudaMemcpyDeviceToHost);
//	cudaMemcpy(testC, sd_piCtoPNum,	 sizeof(int) * m_iCNumMax, cudaMemcpyDeviceToHost);
//	cudaMemcpy(testD, sd_piCtoPIndx, sizeof(int) * m_iCNumMax, cudaMemcpyDeviceToHost);
//
	////ホスト側のデータを転送した結果をダンプ
	//ofstream ofs( "DtoH_Test.txt" );
	//ofs << "DtoH_Test" << endl;

	////for(int i = 0; i < m_iPNumMax; i++)
	//for(int i = 0; i < m_iTNumMax * m_iNeighborMax * 2; i++)
	//{
	//	////等しいなら0になるはず
	//	//ofs << i << "      host-device:: "
	//	//	<< abs(m_piPtoCIndx[i]-testA[i]) << ", "
	//	//	<< abs(m_piPtoTIndx[i]-testB[i]) << ", "
	//	//	<< abs(m_piCtoPNum[i] -testC[i]) << ", "
	//	//	<< abs(m_piCtoPIndx[i]-testD[i]) << ", "
	//	//	<< endl;

	//	//if(m_mk3DiNeighborTetra.Get()[i] == -1 && m_mk3DiNeighborTetra.Get().data()[i]){	continue;	}

	//	//ofs << i << "        data()-Get() = " << m_mk3DiNeighborTetra.Get().data()[i] - m_mk3DiNeighborTetra.Get()[i] << endl;
	//	//ofs << i << " Get() = " << m_mk3DiNeighborTetra.Get()[i] << endl;
	//	//ofs << i << " Get() = " << m_mk3DiNeighborTetra.Get().data()[i] << endl;
	//}
//
//	delete[] testA;
//	delete[] testB;
//	delete[] testC;
//	delete[] testD;
}

//------------------------------------------__初期化-------------------------------------------------

//------------------------------------------＿相変化-------------------------------------------------
void IceStructure::StepObjMelt(
	vector<unsigned>& pList,
	vector<unsigned>& cList,
	vector<unsigned>& tList,
	vector<unsigned>& cLayerList,
	vector<unsigned>& tLayerList)
{
	QueryCounter counter;
	counter.Start();
		double end = counter.End();
	SearchReconstruct_Tetra_Melt(pList, tList, tLayerList);			//再定義四面体の探索
	SearchReconstruct_Cluster_Melt(pList, cList, cLayerList);		//再定義クラスタの探索
	UpdateInfo_Melt_PandT(pList);									//粒子・四面体情報の更新
	UpdateInfo_Melt_PandC(pList, cList);							//粒子・クラスタ情報の更新

	////CheckDeleteCluster();													//同一，包含関係にあるクラスタを削除
	////CheckDeleteTetra(viTetraList, viTLayerList);							//同一，包含関係にある四面体を削除

	SetInfo_Tetra(pList, tList, tLayerList);						//粒子・近傍四面体情報の再定義
	//RXTIMER("SetInfo_Tetra");

//デバッグ
	//DebugStepObjMelt(pList, cList);
}

void IceStructure::StepObjFreeze()
{
	vector<int> viParticleList;														//凝固した粒子集合
	vector<int> viClusterList;														//再定義するクラスタの集合
	vector<int> viCLayerList;														//再定義するクラスタのレイヤー
	vector<int> viTetraList;														//再定義する四面体の集合
	vector<int> viTLayerList;														//再定義する四面体のレイヤー

	//SearchFreezeParticle(viParticleList);											//凝固粒子の探索
	//SetFreezeTetraInfo(viParticleList);												//凝固粒子に関する四面体の作成
	//SetFreezeClusterInfo(viParticleList);											//凝固粒子に関するクラスタの作成
	//SearchReconstructTetra_Freeze(viParticleList, viTetraList, viTLayerList);		//再定義四面体の探索
	//SearchReconstructCluster_Freeze(viParticleList, viClusterList, viCLayerList);	//再定義クラスタの探索

	////CheckDeleteCluster();															//同一，包含関係にあるクラスタを削除
	////CheckDeleteTetra(viTetraList, viTLayerList);									//同一，包含関係にある四面体を削除

	//SetTetraInfo(viParticleList, viTetraList, viTLayerList);						//粒子・近傍四面体情報の再定義
	//SetClusterInfo(viParticleList, viClusterList, viCLayerList);					//粒子・クラスタ情報の再定義

	////デバッグ
	//if(viParticleList.size() == 0 || viClusterList.size() == 0){	return;	}
	//cout << "Debug " << __FUNCTION__ << "viParticleList.size = " << viParticleList.size() << " " << endl;

	//for(unsigned i = 0; i < viParticleList.size(); i++)
	//{
	//	cout << " " << viParticleList[i];
	//}
	//cout << endl;

	//cout << "viClusterList.size =  " << viClusterList.size() << " ";
	//for(unsigned i = 0; i < viClusterList.size(); i++)
	//{
	//	cout << " " << viClusterList[i];
	//}
	//cout << endl;

	//cout << "viCLayerList:: ";
	//for(unsigned i = 0; i < viCLayerList.size(); i++)
	//{
	//	cout << " " << viCLayerList[i];
	//}
	//cout << endl;

	//cout << "viTetraList.size = " << viTetraList.size() << " ";
	//for(unsigned i = 0; i < viTetraList.size(); i++)
	//{
	//	cout << " " << viTetraList[i];
	//}
	//cout << endl;

	//cout << "viTLayerList:: "3
	//for(unsigned i = 0; i < viTLayerList.size(); i++)
	//{
	//	cout << " " << viTLayerList[i];
	//}
	//cout << endl;

	//クラスタ→粒子
	//for(int i = 0; i < m_iClusteresNum; i++){	m_ice->DebugCtoP(i);	}

	//粒子→クラスタ
	//for(int i = 0; i < m_pPS->GetNumParticles(); i++){	m_ice->DebugPtoC(i);	}

	//SMクラスタに含まれる粒子は機能で確認できる
	//for(int i = 0; i < ICENUM; i++)
	//{
	//	cout << "pIndx = " << endl;

	//	for(int j = 0; j < m_sm_cluster[i]->GetNumVertices(); j++)
	//	{
	//		cout << " j = " << m_sm_cluster[i]->GetParticleIndx(j);
	//	}
	//	cout << endl;
	//}

	//四面体→粒子は機能で確認できる

	//粒子→四面体
	//for(int i = 0; i < m_pPS->GetNumParticles(); i++){	m_ice->DebugPtoT(i);	}

	//近傍四面体
	//for(unsigned i = 0; i < m_vviTetraList.size(); i++ ){	m_ice->DebugNeighborTetra(i);	}
}

void IceStructure::SearchReconstruct_Tetra_Melt(const vector<unsigned>& pList, vector<unsigned>& tList, vector<unsigned>& lList)
{
	unsigned pListSize = pList.size();
	if(pListSize == 0){	return;}

	ResetTFlag(m_iTNum);	//四面体探索フラグの初期化

	//１　粒子が含まれていた，いる四面体
	//２　１の四面体の近傍四面体　めちゃくちゃ重い
	//つまりはクラスタを構成した四面体　TODO::覚えられる情報なので，ここの計算コストが高ければ修正可能
	
	//１ 融解粒子が含まれていた四面体
	for(unsigned i = 0; i < pListSize; i++)
	{
		int ipIndx = pList[i];

		for(int j = 0; j < GetPtoTIndx(ipIndx); j++)
		{
			if(GetPtoT(ipIndx, j, 0) == -1 || GetPtoT(ipIndx, j, 1) == -1){	continue;	}

			//四面体探索フラグから既に探索したかどうかを判定
			if( GetTFlag(GetPtoT(ipIndx, j, 0)) )	{	continue;								}
			else									{	SetTFlag(GetPtoT(ipIndx, j, 0), true);	}

			tList.push_back(GetPtoT(ipIndx, j, 0));
			lList.push_back(1);								//0か1かの判断はできないので1に合わせる．
		}
	}

	//return;	//２を処理しないならめちゃくちゃ早くなる

	//２　１の四面体の近傍四面体
	int tetraNum = tList.size();
	for(int i = 0; i < tetraNum; i++)
	{
		int itIndx = tList[i];

		for(int j = 0; j < GetNTNum(itIndx); j++)
		{
			int jtIndx = GetNeighborTetra(itIndx, j, 0);
			int jlIndx = GetNeighborTetra(itIndx, j, 1);

			if(GetTFlag(jtIndx))
			{
				vector<unsigned>::iterator check = std::find(tList.begin(), tList.end(), jtIndx);
				
				//TODO::既に含まれているのなら，layerを比べて小さいほうを優先する
				int layerIndx = check - tList.begin();
				if(lList[layerIndx] > jlIndx)
				{
					lList[layerIndx] = jlIndx;
				}
				continue;
			}
			else
			{
				SetTFlag(jtIndx, true);	
			}

			tList.push_back(jtIndx);
			lList.push_back(jlIndx);
		}
	}
}

void IceStructure::SearchReconstruct_Cluster_Melt(const vector<unsigned>& pList, vector<unsigned>& cList, vector<unsigned>& lList)
{
	unsigned pListSize = pList.size();
	if(pListSize == 0){	return;}	

	//再定義クラスタは，融解粒子が所属していたクラスタ
	for(unsigned i = 0;  i < pListSize; i++)
	{
		unsigned ipIndx = pList[i];

		for(unsigned j = 0, ctopIndx = GetCtoPIndx(ipIndx); j < ctopIndx; j++)
		{
			unsigned jcIndx = GetPtoC(ipIndx, j, 0);
			unsigned joIndx = GetPtoC(ipIndx, j, 1);
			unsigned jlIndx = GetPtoC(ipIndx, j, 2);

			if(jcIndx == -1 || joIndx == -1){	continue;	}

			vector<unsigned>::iterator check = std::find(cList.begin(), cList.end(), jcIndx);

			//TODO::既に再定義クラスタとして取得されているなら，layerを比べて小さいほうを優先する
			if(check != cList.end())
			{
				int layerIndx = check - cList.begin();
				if(lList[layerIndx] > jlIndx)
				{
					lList[layerIndx] = jlIndx;
				}
				continue;
			}

			cList.push_back(jcIndx);
			lList.push_back(jlIndx);
		}
	}
}

void IceStructure::UpdateInfo_Melt_PandT(const vector<unsigned>& pList)
{
	unsigned pListSize = pList.size();
	if(pListSize == 0){	return;}	


	for(unsigned i = 0; i < pListSize; i++)
	{
		int ipIndx = pList[i];

		for(int j = 0; j < GetPtoTIndx(ipIndx); j++)
		{
			int tIndx = GetPtoT(ipIndx, j, 0);
			int oIndx = GetPtoT(ipIndx, j, 1);	//この場合はそのまま添え字
			
			if(tIndx == -1 || oIndx == -1){ continue;	}

			DeleteTtoP(tIndx, oIndx);
		}

		ClearPtoT(ipIndx);
	}
}

void IceStructure::UpdateInfo_Melt_PandC(const vector<unsigned>& pList, const vector<unsigned>& cList)
{
	int pListSize = pList.size();
	int cListSize = cList.size();

	if(pListSize == 0 || cListSize == 0){	return; }

	//並列処理で用いる変数をまとめて定義
	int j= 0, k = 0;
	int icIndx = 0;
	int jpIndx = 0;

	//融解粒子＝クラスタ情報を，融解クラスタに含まれていた粒子から取り除く
	#pragma omp parallel
	{
	#pragma omp for private(j,k,icIndx,jpIndx)
		for(int i = 0; i < pListSize; i++)
		{
			icIndx = pList[i];		//各粒子にクラスタが用意されているのでpIndx->cIndxとしている
	
			for(j = 0; j < GetCtoPIndx(icIndx); j++)
			{
				jpIndx = GetCtoP(icIndx, j, 0);		//融解クラスタに含まれていた粒子
													//この粒子から，融解したクラスタの情報を取り除く
				
				if(jpIndx == -1){	continue;	}

				for(k = 0; k < GetPtoCIndx(jpIndx); k++)
				{					
					if(GetPtoC(jpIndx, k, 0) == -1
					|| GetPtoC(jpIndx, k, 1) == -1
					|| GetPtoC(jpIndx, k, 0) != icIndx)
					{
						continue;
					}
	
					#pragma omp critical (DeletePtoC)	//TODO：：後にカウントしたほうが並列化できてよい
					{
						DeletePtoC(jpIndx, k);
					}

					break;			//同じクラスタに複数所属することは無いので，break
				}
			}
		}
	}//end #pragma omp parallel

	//融解粒子→クラスタ情報を削除
	#pragma omp parallel
	{
	#pragma omp for
		for(int i = 0; i < pListSize; i++)
		{
			ClearPtoC(pList[i]);
		}
	}//end #pragma omp parallel

	//融解クラスタ→粒子情報を削除
	#pragma omp parallel
	{
	#pragma omp for
		for(int i = 0; i < pListSize; i++)
		{
			ClearCtoP(pList[i]);
			//m_iceObj->GetMoveObj(pList[i])->Clear();
		}
	}//end #pragma omp parallel

	//再定義クラスタに含まれる粒子から、再定義クラスタの情報を消去
	#pragma omp parallel
	{
	#pragma omp for private(j,k,icIndx,jpIndx)
		for(int i = 0; i < cListSize; i++)
		{
			icIndx = cList[i];
	
			for(j = 0; j < GetCtoPIndx(icIndx); j++)
			{
				jpIndx = GetCtoP(icIndx, j, 0);
				
				if(jpIndx == -1){	continue;	}

				for(k = 0; k < GetPtoCIndx(jpIndx); k++)
				{
					if(GetPtoC(jpIndx, k, 0) == -1
					|| GetPtoC(jpIndx, k, 1) == -1
					|| GetPtoC(jpIndx, k, 0) != icIndx)
					{
						continue;
					}

					#pragma omp critical (DeletePtoC)	//TODO：：後にカウントしたほうが並列化できてよい
					{
						DeletePtoC(jpIndx, k);
					}

					break;			//同じクラスタに複数所属することは無いので，break
				}
			}
		}
	}//end #pragma omp parallel
}

void IceStructure::SetInfo_Tetra(const vector<unsigned>& pList, const vector<unsigned>& tList, const vector<unsigned>& lList)
{
	int pListSize = pList.size();
	int tListSize = tList.size();

	if(pListSize == 0 || tListSize == 0){	return; }

	int itIndx = 0;
	int ilayer = 0;

	#pragma omp parallel
	{
	#pragma omp for private(itIndx,ilayer)
		//近傍四面体の再定義
		for(int i = 0; i < tListSize; i++)
		{
			itIndx = tList[i];
			ilayer = lList[i];

			ClearNeighborTetraFromLayer(itIndx, ilayer);
			SetNeighborTetraFromLayer(itIndx, m_iLayer, ilayer);	//ここが非常に重い
		}
	}//end #pragma omp parallel
}


//------------------------------------------相変化＿-------------------------------------------------




//-------------------------------------------取得----------------------------------------
/*!
 * 粒子→四面体　配列の空き場所を探して添え字を返す
 */
int IceStructure::GetPtoTFreeIndx(int pIndx)
{
	int freeIndx = -1;

	for(int i = 0; i < m_iPtoTMax; i++)
	{
		if(GetPtoT(pIndx, i, 0) != -1 || GetPtoT(pIndx, i, 1) != -1){	continue;	}
		freeIndx = i;	break;
	}

	if(freeIndx == m_iPtoTMax || freeIndx == -1)
	{	
		cout << __FUNCTION__ << " Error::配列に空きがありません " << endl;
		freeIndx = 0;
	}

	return freeIndx;
}

/*!
 * 粒子→クラスタ　配列の空き場所を探して添え字を返す
 */
int IceStructure::GetPtoCFreeIndx(int pIndx)
{
	int freeIndx = -1;

	for(int i = 0; i < m_iPtoCMax; i++)
	{
		//if(GetPtoC(pIndx, i)[0] != -1 || GetPtoC(pIndx, i)[1] != -1){	continue;	}
		if(GetPtoC(pIndx, i, 0) != -1 || GetPtoC(pIndx, i, 1) != -1 || GetPtoC(pIndx, i, 2) != -1)
		{
			continue;
		}

		freeIndx = i;	break;
	}

	if(freeIndx == m_iPtoCMax || freeIndx == -1)
	{	
		cout << __FUNCTION__ << " Error::配列に空きがありません " << endl;
		freeIndx = 0;
	}

	return freeIndx;
}

//-------------------------------------------取得----------------------------------------

//-------------------------------------------書き込み----------------------------------------
/*!
 * 登録処理　粒子→四面体　TODO::こっちで-1を見つけて埋めていったほうがいい
 * @param[in] pIndx　粒子番号
 * @param[in] lIndx　粒子が属するl番目のクラスタ
 * @param[in] cIndx　クラスタ番号
 * @param[in] oIndx　クラスタ内での順序番号
 */
void IceStructure::SetPtoT(int pIndx, int lIndx, int tIndx, int oIndx)
{
	//エラーチェック
	if(m_iPtoTMax <= lIndx)
	{
		cout << __FUNCTION__ << " Error::粒子が属する四面体情報を格納できなくなりました．メモリが足りません" << endl;
		cout << m_iPtoTMax << "<" << lIndx << endl;
		return;
	}
	else if(lIndx < 0)
	{
		cout << __FUNCTION__ << " Error::添え字が不正な値です．" << endl;
		cout << "lIndx = " << lIndx << endl;
	}

	m_mk3DiPtoT(pIndx, lIndx, 0) = tIndx;
	m_mk3DiPtoT(pIndx, lIndx, 1) = oIndx;
}

/*!
 * 登録処理　四面体→粒子
 * @param[in] tIndx　　　四面体番号
 * @param[in] pIndxList　所属粒子配列
 */
void IceStructure::SetTtoP(int tIndx, vector<int>& pIndxList)
{
	//エラーチェック
	if(4 < pIndxList.size())
	{
		cout << __FUNCTION__ << " Error::四面体が含む粒子の情報を格納できなくなりました．メモリが足りません" << endl;
		cout << 4 << "<" << pIndxList.size() << endl;
		return;
	}

	for(int i = 0; i < pIndxList.size(); i++)
	{
		m_mk2DiTtoP(tIndx, i) = pIndxList[i];
	}
}

/*!
 * 登録処理　粒子→クラスタ
 * @param[in] pIndx　粒子番号
 * @param[in] lIndx　粒子が属するl番目のクラスタ
 * @param[in] cIndx　クラスタ番号
 * @param[in] oIndx　クラスタ内での順序番号
 * @param[in] layer　クラスタ内での層数
 */
void IceStructure::SetPtoC(int pIndx, int lIndx, int cIndx, int oIndx, int layer)
{
	//エラーチェック
	if(m_iPtoCMax <= lIndx)
	{
		cout << __FUNCTION__ << " Error::粒子が属するクラスタ情報を格納できなくなりました．メモリが足りません" << endl;
		cout << m_iPtoCMax << "<" << lIndx << endl;
		return;
	}
	
	m_mk3DiPtoC(pIndx, lIndx, 0) = cIndx;
	m_mk3DiPtoC(pIndx, lIndx, 1) = oIndx;
	m_mk3DiPtoC(pIndx, lIndx, 2) = layer;
}

/*!
 * 登録処理　クラスタ→粒子
 * @param[in] cIndx　　　クラスタ番号
 * @param[in] pIndxList　所属粒子配列
 */
void IceStructure::SetCtoP(int cIndx, const vector<int>& pIndxList, int* pLayerList)
{
	//エラーチェック
	if(m_iCtoPMax <= pIndxList.size())
	{
		cout << __FUNCTION__ << " Error::クラスタが含む粒子の情報を格納できなくなりました．メモリが足りません" << endl;
		cout << m_iCtoPMax << "<" << pIndxList.size() << endl;
		return;
	}

	for(int i = 0; i < pIndxList.size(); i++)
	{
		m_mk3DiCtoP(cIndx, i, 0) = pIndxList[i];
		m_mk3DiCtoP(cIndx, i, 1) = pLayerList[i];
	}
}

/*!
 * 登録処理　初期の近傍四面体
 * @param[in] tIndx　　四面体番号
 * @param[in] PtoTNum　
 */
void IceStructure::SetTetraInfo(int tIndx, int* PtoTNum)
{
	IceTetrahedra &tetra = IceTetrahedra::GetInstance();		//あんまりここでは呼びたくなかった

	//粒子が属している四面体の番号を登録するための準備
	//pCountListには，tIndx番目の四面体に含まれる各粒子が，それぞれいくつの四面体に属するかを求めて保存する
	int* pCountList = new int[tetra.GetTetraList(tIndx).size()];
	//cout << __FUNCTION__ << "::check0" << endl;
	for(int j = 0; j < tetra.GetTetraList(tIndx).size(); j++)
	{
		int pIndx = tetra.GetTetraList(tIndx)[j];
		pCountList[j] = GetPtoTNum(pIndx)-PtoTNum[pIndx];
		PtoTNum[pIndx]--;
	}
	//cout << __FUNCTION__ << "::check1" << endl;
	//粒子と四面体の情報登録
	vector<int>& pIndxList = tetra.GetTetraList(tIndx);

	for(int i = 0; i < GetTtoPNum(tIndx); i++)
	{
		SetPtoT(pIndxList[i], pCountList[i], tIndx, i);
	}
	//cout << __FUNCTION__ << "::check2" << endl;

	SetTtoP(tIndx, pIndxList);
	//cout << __FUNCTION__ << "::check3" << endl;

	delete[] pCountList;
}


/*!
 * 登録処理　初期の近傍四面体
 * @param[in] tIndx　四面体番号
 * @param[in] layer　探索階層
 */
void IceStructure::SetNeighborTetra(int tIndx, int layer)
{
	/*	
		ある四面体Aに含まれている粒子が，他の複数の四面体に含まれている場合，
		その複数の四面体は四面体Aの近傍四面体となる
	*/
	vector<int> pIndxList;

	//layer=1層目
	for(int i = 0; i < GetTtoPIndx(tIndx); i++)
	{
		int ipIndx = GetTtoP(tIndx, i);
		if(ipIndx == -1){	continue;	}

		if(find(pIndxList.begin(), pIndxList.end(), ipIndx) != pIndxList.end())
		{	
			continue;
		}
		pIndxList.push_back(ipIndx);

		//探索
		for(int j = 0; j < GetPtoTIndx(ipIndx); j++)
		{
			if(GetPtoT(ipIndx, j, 0) == -1
			|| GetPtoT(ipIndx, j, 1) == -1
			|| GetPtoT(ipIndx, j, 0) == tIndx)
			{	
				continue;	
			}

			//同じクラスタを既に含んでいないかのチェック
			if(CheckNeighborTetra( tIndx, GetPtoT(ipIndx, j, 0) ) != -1){	continue;	}

			//近傍クラスタの登録＋カウント
			if(GetNTNum(tIndx) >= m_iNeighborMax)
			{
				cout << __FUNCTION__ << "近傍粒子のメモリが足りなくなりました" << endl;
				cout << GetNTNum(tIndx) << " < " << m_iNeighborMax << endl;
			}
			else
			{
				m_mk3DiNeighborTetra(tIndx, GetNTNum(tIndx), 0) = GetPtoT(ipIndx, j, 0);
				m_mk3DiNeighborTetra(tIndx, GetNTNum(tIndx), 1) = 1;
				CountNT(tIndx);
			}
		}
	}

	//layer層目をたどって近傍四面体を取得する
	int nowSize = 0, nowIndx = 0;
	int d_addNum = 0;
	for(int i = 2; i <= layer; i++)
	{
		nowSize = GetNTNum(tIndx);
		d_addNum = GetNTNum(tIndx);

		//探索するクラスタをnowSizeとnowIndxで制限している
		for(int j = nowIndx; j < nowSize; j++)
		{
			int jtIndx = GetNeighborTetra(tIndx, j, 0);			//近傍四面体のひとつ
			
			//近傍四面体に含まれる粒子が他の四面体にも含まれている場合，その四面体を近傍として登録
			for(int k = 0; k < GetTtoPIndx(jtIndx); k++)
			{
				int kpIndx = GetTtoP(jtIndx, k);				//近傍四面体に含まれている粒子のひとつ
				if(kpIndx == -1){ continue;	}

				if(find(pIndxList.begin(), pIndxList.end(), kpIndx) != pIndxList.end())
				{
					continue;
				}
				pIndxList.push_back(kpIndx);

				for(int l = 0; l < GetPtoTIndx(kpIndx); l++)
				{
					if(GetPtoT(kpIndx, l, 0) == -1 
					|| GetPtoT(kpIndx, l, 1) == -1
					|| GetPtoT(kpIndx, l, 0) == tIndx)
					{
						continue;
					}

					//同じ四面体を既に含んでいないかのチェック
					if(CheckNeighborTetra( tIndx, GetPtoT(kpIndx, l, 0) ) != -1)
					{	//cout << "check4" << endl;
						continue;
					}

					//近傍クラスタの登録＋カウント
					if(GetNTNum(tIndx) >= m_iNeighborMax)
					{
						cout << __FUNCTION__ << " 近傍四面体のメモリが足りなくなりました " << GetNTNum(tIndx) << endl;
						cout << GetNTNum(tIndx) << " < " << m_iNeighborMax << endl;
					}
					else
					{
//						if(GetNTNum(tIndx) > 200 ) return;	//近傍四面体数の制限
						m_mk3DiNeighborTetra(tIndx, GetNTNum(tIndx), 0) = GetPtoT(kpIndx, l, 0);
						m_mk3DiNeighborTetra(tIndx, GetNTNum(tIndx), 1) = i;
						CountNT(tIndx);
					}
				}
			}
		}
		nowIndx = nowSize;										//次のループ開始時のスタート番号を更新
		//cout << "tIndx = " << tIndx << " addNum = " << GetNTNum(tIndx)-d_addNum << endl;
	}
}

/*!
 * 登録処理　初期の近傍四面体
 * @param[in] tIndx　四面体番号
 * @param[in] searchLayer　探索終了階層 searchLayer >= 1
 * @param[in] deleteLayer　探索開始階層 deleteLayer >= 1
 */
void IceStructure::SetNeighborTetraFromLayer(int tIndx, int searchLayer, int deleteLayer)
{//	cout << __FUNCTION__ << endl;
	/*	
		ある四面体Aに含まれている粒子が，他の複数の四面体に含まれている場合，
		その複数の四面体は四面体Aの近傍四面体となる
	*/
	vector<int> pIndxList;		//探索済み粒子を保存　これのおかげでかなり早くなる．

	//layer=1層目
	for(int i = 0; i < GetTtoPIndx(tIndx); i++)
	{
		int ipIndx = GetTtoP(tIndx, i);
		if(ipIndx == -1){	continue;	}

		if(find(pIndxList.begin(), pIndxList.end(), ipIndx) != pIndxList.end()){	continue;		}
		pIndxList.push_back(ipIndx);

		if(deleteLayer != 1){ continue;	}		//1層目のみ行う

		//探索
		for(int j = 0; j < GetPtoTIndx(ipIndx); j++)
		{
			if(GetPtoT(ipIndx, j, 0) == -1
			|| GetPtoT(ipIndx, j, 1) == -1
			|| GetPtoT(ipIndx, j, 0) == tIndx)
			{
				continue;
			}

			//同じクラスタを既に含んでいないかのチェック
			if(CheckNeighborTetra(tIndx, GetPtoT(ipIndx, j, 0)) != -1){	continue;	}

			//近傍クラスタの登録＋カウント
			if(GetNTNum(tIndx) >= m_iNeighborMax)
			{
				cout << __FUNCTION__ << "近傍粒子のメモリが足りなくなりました" << endl;
				cout << GetNTNum(tIndx) << " < " << m_iNeighborMax << endl;
			}
			else
			{
				m_mk3DiNeighborTetra(tIndx, GetNTNum(tIndx), 0) = GetPtoT(ipIndx, j, 0);
				m_mk3DiNeighborTetra(tIndx, GetNTNum(tIndx), 1) = 1;
				CountNT(tIndx);
			}
		}
	}

	if(deleteLayer == 1)
	{
		deleteLayer += 1;			//全て初期化する場合は+1する
	}
		
	//layer層目をたどって近傍四面体を取得する
	int nowSize = 0, nowIndx = 0;

	for(int i = 2; i <= searchLayer; i++)
	{
		nowSize = GetNTNum(tIndx);

		//探索するクラスタをnowSizeとnowIndxで制限している
		for(int j = nowIndx; j < nowSize; j++)
		{
			int jtIndx = GetNeighborTetra(tIndx, j, 0);				//近傍四面体のひとつ
			
			//近傍四面体に含まれる粒子が他の四面体にも含まれている場合，その四面体を近傍として登録
			//TODO::四面体→粒子→近傍四面体　ではなく，四面体→近傍四面体　とする
			//TODO::AはBの近傍であるなら，BはAの近傍である　を利用する
			//TODO::融解した粒子が所属していた四面体から分離するかを逆算する．
			//TODO::再構成しない四面体を利用する
			for(int k = 0; k < GetTtoPIndx(jtIndx); k++)
			{
				int kpIndx = GetTtoP(jtIndx, k);					//近傍四面体に含まれている粒子のひとつ
				if(kpIndx == -1){	continue;	}
				if(find(pIndxList.begin(), pIndxList.end(), kpIndx) != pIndxList.end()){	continue;}
				pIndxList.push_back(kpIndx);

				if(i < deleteLayer){	continue;	}				//際探索する必要ないなら，粒子を追加しただけで戻る

				//粒子が所属している近傍四面体を探索
				for(int l = 0; l < GetPtoTIndx(kpIndx); l++)
				{
					if(GetPtoT(kpIndx, l, 0) == -1
					|| GetPtoT(kpIndx, l, 1) == -1
					|| GetPtoT(kpIndx, l, 0) == tIndx)
					{
						continue;
					}

					if(CheckNeighborTetra( tIndx, GetPtoT(kpIndx, l, 0) ) != -1){	continue;	}	//同じ四面体を既に含んでいないかのチェック

					//近傍クラスタの登録＋カウント
					if(GetNTNum(tIndx) >= m_iNeighborMax)
					{
						cout << __FUNCTION__ << " 近傍四面体のメモリが足りなくなりました " << endl;
						cout << GetNTNum(tIndx) << " < " << m_iNeighborMax << endl;
					}
					else
					{
//						if(GetNTNum(tIndx) > 200 ) return;	//近傍四面体数の制限
						m_mk3DiNeighborTetra(tIndx, GetNTNum(tIndx), 0) = GetPtoT(kpIndx, l, 0);
						m_mk3DiNeighborTetra(tIndx, GetNTNum(tIndx), 1) = i;
						CountNT(tIndx);
					}
				}
			}
		}
		nowIndx = nowSize;											//次のループ開始時のスタート番号を更新
	}
}

//-------------------------------------------書き込み----------------------------------------

/*!
 * 取得処理　粒子→四面体
 * @param[in] pIndx　粒子番号
 * @param[in] lIndx　粒子内番号
 * @param[in] oIndx	0->クラスタの番号, 1->クラスタ内での粒子番号, 2->階層
 */

int IceStructure::GetPtoT(int pIndx, int lIndx, int oIndx)
{
	return m_mk3DiPtoT(pIndx, lIndx, oIndx);
}

/*!
 * 取得処理　四面体→粒子
 * @param[in] tIndx　四面体番号
 * @param[in] lIndx　四面体内番号
 */
int IceStructure::GetTtoP(int tIndx, int lIndx)
{
	return m_mk2DiTtoP(tIndx, lIndx);
}

/*!
 * 取得処理　粒子→クラスタ
 * @param[in] pIndx	粒子番号
 * @param[in] lIndx	粒子内番号
 * @param[in] oIndx	0->クラスタの番号, 1->クラスタ内での粒子番号, 2->階層
 */
int IceStructure::GetPtoC(int pIndx, int lIndx, int oIndx)
{
	return m_mk3DiPtoC(pIndx, lIndx, oIndx);
}

/*!
 * 取得処理　クラスタ→粒子
 * @param[in] cIndx クラスタ番号
 * @param[in] lIndx クラスタ内粒子番号
 * @param[in] oIndx 
 */
const int& IceStructure::GetCtoP(const int& cIndx, const int& lIndx, const int& oIndx)
{
	return m_mk3DiCtoP(cIndx, lIndx, oIndx);
}

/*!
 * 取得処理　四面体→近傍四面体
 * @param[in] tIndx　四面体番号
 * @param[in] lIndx　四面体内番号
 * @param[in] oIndx 
 */
int IceStructure::GetNeighborTetra(int tIndx, int lIndx, int oIndx)
{
	return m_mk3DiNeighborTetra(tIndx, lIndx, oIndx);
}

/*!
 * 四面体に含まれている粒子数のカウント，粒子が所属する四面体数のカウント
 * @param[in] tIndx　四面体番号
 * @param[in] pList　四面体に含まれる粒子リスト
 */
void IceStructure::CountTetrahedra(int tIndx, vector<int>& pList)
{
	for(unsigned i = 0; i < pList.size(); i++)
	{
		int pIndx = pList[i];		

		CountPtoT(pIndx);
		CountTtoP(tIndx);

		//Indxの更新
		if(GetPtoTNum(pIndx) >= GetPtoTIndx(pIndx))
		{
			SetPtoTIndx(pIndx, GetPtoTNum(pIndx));
		}
	}
	
	//Indxの更新
	if(GetTtoPNum(tIndx) >= GetTtoPIndx(tIndx))
	{
		SetTtoPIndx(tIndx, GetTtoPNum(tIndx));
	}
}

/*!
 * クラスタに含まれている粒子数のカウント，粒子が所属する四面体数のカウント
 * @param[in] cIndx　粒子番号
 * @param[in] pList　クラスタに含まれる粒子リスト
 */
void IceStructure::CountClusterParticle(int cIndx, vector<int>& pList, int pNum)
{
	for(int j = 0; j < pNum; j++)
	{
		int jpIndx = pList[j];
		CountPtoC(jpIndx);										//粒子が接続クラスタに所属する個数のカウント
		CountCtoP(cIndx);										//接続クラスタが粒子を含む個数のカウント
	
		//Indxの更新
		if(GetPtoCNum(jpIndx) >= GetPtoCIndx(jpIndx))
		{
			SetPtoCIndx(jpIndx, GetPtoCNum(jpIndx));
		}
	}

	//Indxの更新
	if(GetCtoPNum(cIndx) >= GetCtoPIndx(cIndx))
	{
		SetCtoPIndx(cIndx, GetCtoPNum(cIndx));
	}
}

/*!
 * 削除処理　四面体→粒子
 * @param[in] tIndx　四面体番号
 * @param[in] lIndx　四面体内番号
 */
void IceStructure::DeleteTtoP(int tIndx, int lIndx)
{//	cout << __FUNCTION__ << " check1" << endl;
	//今回，lIndxは添え字
	if(4 <= lIndx)
	{
		cout << __FUNCTION__ << " Error:: 4 < lIndx 四面体配列へのアクセスエラー" << endl;
		cout << 4 << "<" << lIndx << endl;
		return;
	}

	if(m_piTtoPNum[tIndx] > 0)
	{
		m_piTtoPNum[tIndx]--;
	}
	else
	{
		cout << __FUNCTION__ << " Error::m_piTtoPNum[" << tIndx << "] < 0" << endl;
		return;
	}

	m_mk2DiTtoP(tIndx, lIndx) = -1;
}

/*!
 * 削除処理　粒子→クラスタ
 * @param[in] pIndx　粒子番号
 * @param[in] lIndx　粒子内番号
 */
void IceStructure::DeletePtoC(int pIndx, int lIndx)
{//	cout << __FUNCTION__ << " check1" << endl;
	//今回，lIndxは添え字
	if(m_iPtoCMax <= lIndx)
	{
		cout << __FUNCTION__ << " Error:: m_iPtoCMax < lIndx" << endl;
		cout << m_iPtoCMax << "<" << lIndx << endl;
		return;
	}

	if(m_piPtoCNum[pIndx] > 0)
	{
		m_piPtoCNum[pIndx]--;
	}
	else
	{
		cout << __FUNCTION__ << " Error::m_piPtoCNum[" << pIndx << "] < 0" << endl;
		return;
	}

	//m_pppiPtoC[pIndx][lIndx][0] = -1;
	//m_pppiPtoC[pIndx][lIndx][1] = -1;
	//m_pppiPtoC[pIndx][lIndx][2] = -1;

	for(int i = 0; i < 3; i++)
	{
		m_mk3DiPtoC(pIndx, lIndx, i) = -1;
	}
}

/*!
 * 削除処理　粒子→四面体
 * @param[in] pIndx　粒子番号
 * @param[in] lIndx　粒子内番号
 */
void IceStructure::DeletePtoT(int pIndx, int lIndx)
{//	cout << __FUNCTION__ << " check1" << endl;
	//今回，lIndxは添え字
	if(m_iPtoTMax <= lIndx)
	{
		cout << __FUNCTION__ << " Error:: m_iPtoTMax < lIndx" << endl;
		cout << m_iPtoTMax << "<" << lIndx << endl;
		return;
	}

	if(m_piPtoTNum[pIndx] > 0)
	{
		m_piPtoTNum[pIndx]--;
	}
	else
	{
		cout << __FUNCTION__ << " Error::m_piPtoTNum[" << pIndx << "] < 0, lIndx = " << lIndx << endl;
		return;
	}

	for(int i = 0; i < 2; i++)
	{
		m_mk3DiPtoT(pIndx, lIndx, i) = -1;
	}
}

/*!
 * 初期化処理　粒子→四面体
 */
void IceStructure::ClearPtoT(int pIndx)
{
	for(int i = 0; i < m_iPtoTMax; i++)
	{
		for(int j = 0; j < 2; j++)
		{
			m_mk3DiPtoT(pIndx, i, j) = -1;
		}
	}

	m_piPtoTIndx[pIndx] = 0;
	m_piPtoTNum[pIndx] = 0;
}

/*!
 * 初期化処理　粒子→クラスタ
 */
void IceStructure::ClearPtoC(int pIndx)
{
	for(int i = 0; i < m_iPtoCMax; i++)
	{
		for(int j = 0; j < 3; j++)
		{
			m_mk3DiPtoC(pIndx, i, j) = -1;
		}
	}

	m_piPtoCIndx[pIndx] = 0;
	m_piPtoCNum[pIndx] = 0;
}

/*!
 * 初期化処理　クラスタ→粒子
 */
void IceStructure::ClearCtoP(int cIndx)
{
	for(int i = 0; i < m_piCtoPIndx[cIndx]; i++)
	{
		m_mk3DiCtoP(cIndx, i, 0) = -1;
		m_mk3DiCtoP(cIndx, i, 1) = -1;
	}

	m_piCtoPIndx[cIndx] = 0;
	m_piCtoPNum[cIndx] = 0;
}

/*!
 * 初期化処理　四面体→粒子
 */
void IceStructure::ClearTtoP(int tIndx)
{
	for(int i = 0; i < m_piTtoPIndx[tIndx]; i++)
	{
		m_mk2DiTtoP(tIndx, i) = -1;
	}

	m_piTtoPIndx[tIndx] = 0;
	m_piTtoPNum[tIndx] = 0;

	ClearNeighborTetra(tIndx);		//近傍四面体も初期化
}

/*!
 * 初期化処理　近傍四面体
 */
void IceStructure::ClearNeighborTetra(int tIndx)
{//	cout << __FUNCTION__ << endl;
	for(int i = 0; i < m_piNTNum[tIndx]; i++)
	{
		m_mk3DiNeighborTetra(tIndx, i, 0) = -1;
		m_mk3DiNeighborTetra(tIndx, i, 1) = -1;
	}

	m_piNTNum[tIndx] = 0;
}

/*!
 * 部分的初期化処理　近傍四面体
 *@paramm[in] layer 削除する層 layer >= 1
 */
void IceStructure::ClearNeighborTetraFromLayer(int tIndx, int layer)
{
	for(int i = m_piNTNum[tIndx]-1; 0 <=i ; i--)
	{
		if(m_mk3DiNeighborTetra(tIndx, i, 1) >= layer)
		{
			m_mk3DiNeighborTetra(tIndx, i, 0) = -1;
			m_mk3DiNeighborTetra(tIndx, i, 1) = -1;
			m_piNTNum[tIndx]--;
		}
		else
		{
			break;
		}
	}
}

/*!
 * 判定処理　ある四面体が近傍四面体として含まれているかのチェック
 * @param[in] tIndx　　　四面体番号
 * @param[in] checkTIndx　確認する四面体番号
 */
int IceStructure::CheckNeighborTetra(int tIndx, int checkTIndx)
{
	int findIndx = -1;
	
	for(int k = 0; k < GetNTNum(tIndx); k++)
	{
		if(checkTIndx == GetNeighborTetra(tIndx, k, 0))
		{
			findIndx = k; break;
		}
	}

	return findIndx;
}

//-------------------------------------更新----------------------------------------

//-------------------------------------選択的運動計算------------------------------------------

void IceStructure::UpdateMotionCalcCluster(unsigned cIndx, short unsigned num)
{
	//cout << __FUNCTION__ << ", cIndx = " << cIndx << ", num = " << num << endl;
	m_psuSelectClusterIndx[cIndx] = num;
}

short unsigned IceStructure::GetMotionCalcCluster(unsigned cIndx)
{
	return m_psuSelectClusterIndx[cIndx];
}


//-------------------------------------デバッグ----------------------------------------
void IceStructure::DebugPtoT(int pIndx)
{	cout << __FUNCTION__ << " pIndx = " << pIndx;
	cout << " num=" << GetPtoTNum(pIndx) << " Indx=" << GetPtoTIndx(pIndx);
	
	for(int i = 0; i < GetPtoTIndx(pIndx); i++)
	{
		cout << " c=" << GetPtoT(pIndx, i, 0) << " o=" << GetPtoT(pIndx, i, 1);
	}
	cout << endl;
}

void IceStructure::DebugPtoC(int pIndx)
{	cout << __FUNCTION__ << " pIndx = " << pIndx;
	cout << " num=" << GetPtoCNum(pIndx) << " Indx=" << GetPtoCIndx(pIndx);
	
	for(int i = 0; i < GetPtoCIndx(pIndx); i++)
	{
		cout << " c=" << GetPtoC(pIndx, i, 0) << " o=" << GetPtoC(pIndx, i, 1);
	}
	cout << endl;

}

void IceStructure::DebugCtoP(int cIndx)
{	cout << __FUNCTION__ << " cIndx = " << cIndx;
	cout << " num=" << GetCtoPNum(cIndx) << " Indx=" << GetCtoPIndx(cIndx);

	for(int i = 0; i < GetCtoPIndx(cIndx); i++)
	{
		cout <<" p=" << GetCtoP(cIndx, i, 0) << " l=" << GetCtoP(cIndx, i, 1);
	}
	cout << endl;
}

void IceStructure::DebugTtoP(int tIndx)
{	cout << __FUNCTION__ << " tIndx=" << tIndx;
	cout << " num=" << GetTtoPNum(tIndx) << " Indx=" << GetTtoPIndx(tIndx);

	for(int i = 0; i < GetTtoPIndx(tIndx); i++)
	{
		cout << " " << GetTtoP(tIndx, i);
	}
	cout << endl;
}

void IceStructure::DebugNeighborTetra(int tIndx)
{	cout << __FUNCTION__ << " tIndx=" << tIndx;
	cout << " num = " << GetNTNum(tIndx);
	for(unsigned j = 0; j < GetNTNum(tIndx); j++)
	{
		cout << " NC=" << GetNeighborTetra(tIndx, j, 0) << " Ly=" << GetNeighborTetra(tIndx, j, 1);
	}
	cout << endl;
}

void IceStructure::DebugStepObjMelt(vector<unsigned>& pList, vector<unsigned>& cList)
{	cout << __FUNCTION__ << endl;
	if(pList.size() == 0){	return;	}
	
	cout << "pList.size = " << pList.size() << " ";
	sort(pList.begin(), pList.end());

	for(unsigned i = 0; i < pList.size(); i++)
	{
		cout << " " << pList[i];
	}
	cout << endl;

	cout << "cList.size =  " << cList.size() << " ";
	sort(cList.begin(), cList.end());

	for(unsigned i = 0; i < cList.size(); i++)
	{
		cout << " " << cList[i];
	}
	cout << endl;

	//cout << "viCLayerList:: ";
	//for(unsigned i = 0; i < viCLayerList.size(); i++)
	//{
	//	cout << " " << viCLayerList[i];
	//}
	//cout << endl;

	//cout << "viTetraList.size = " << viTetraList.size() << " ";
	//for(unsigned i = 0; i < viTetraList.size(); i++)
	//{
	//	cout << " " << viTetraList[i];
	//}
	//cout << endl;

	//cout << "viTLayerList:: ";
	//for(unsigned i = 0; i < viTLayerList.size(); i++)
	//{
	//	cout << " " << viTLayerList[i];
	//}
	//cout << endl;

	////クラスタ→粒子
	//for(int i = 0; i < m_iClusteresNum; i++){	m_ice->DebugCtoP(i);	}

	//粒子→クラスタ
	//for(int i = 0; i < m_pPS->GetNumParticles(); i++){	m_ice->DebugPtoC(i);	}

	//SMクラスタに含まれる粒子は機能で確認できる
	//for(int i = 0; i < ICENUM; i++)
	//{
	//	cout << "pIndx = " << endl;

	//	for(int j = 0; j < m_sm_cluster[i]->GetNumVertices(); j++)
	//	{
	//		cout << " j = " << m_sm_cluster[i]->GetParticleIndx(j);
	//	}
	//	cout << endl;
	//}

	//四面体→粒子は機能で確認できる

	//粒子→四面体
	//for(int i = 0; i < m_pPS->GetNumParticles(); i++){	m_ice->DebugPtoT(i);	}

	//近傍四面体
	//for(unsigned i = 0; i < m_vviTetraList.size(); i++ ){	m_ice->DebugNeighborTetra(i);	}
}

//テスト
void IceStructure::TestStepObjMelt(
	vector<unsigned>& pList,
	vector<unsigned>& cList,
	vector<unsigned>& tList,
	vector<unsigned>& cLayerList,
	vector<unsigned>& tLayerList)
{
	QueryCounter counter1;
	QueryCounter counter2;
	QueryCounter counter3;
	QueryCounter counter4;
	QueryCounter counter5;

counter1.Start();
	SearchReconstruct_Tetra_Melt(pList, tList, tLayerList);			//再定義四面体の探索
double end1 = counter1.End();

counter2.Start();
	SearchReconstruct_Cluster_Melt(pList, cList, cLayerList);		//再定義クラスタの探索
double end2 = counter2.End();

counter3.Start();
	UpdateInfo_Melt_PandT(pList);									//粒子・四面体情報の更新
double end3 = counter3.End();

counter4.Start();
	UpdateInfo_Melt_PandC(pList, cList);							//粒子・クラスタ情報の更新
double end4 = counter4.End();

	////CheckDeleteCluster();													//同一，包含関係にあるクラスタを削除
	////CheckDeleteTetra(viTetraList, viTLayerList);							//同一，包含関係にある四面体を削除

counter5.Start();
	SetInfo_Tetra(pList, tList, tLayerList);						//粒子・近傍四面体情報の再定義
double end5 = counter5.End();

	if(pList.size() == 0) return;
	cout << "SearchReconstruct_Tetra_Melt	:" << end1 << endl;
	cout << "SearchReconstruct_Cluster_Melt	:" << end2 << endl;
	cout << "UpdateInfo_Melt_PandT		:" << end3 << endl;
	cout << "UpdateInfo_Melt_PandC		:" << end4 << endl;
	cout << "SetInfo_Tetra			:" << end5 << endl;

//デバッグ
	//DebugStepObjMelt(pList, cList);
}