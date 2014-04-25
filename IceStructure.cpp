//-----------------------------------------------------------------------------
// インクルードファイル
//-----------------------------------------------------------------------------
#include "IceStructure.h"

/*!
 * @param[in] pNumMax　最大粒子数
 * @param[in] cNumMax　最大クラスタ数
 * @param[in] tNumMax　最大四面体数
 */
IceStructure::IceStructure(int pNumMax, int cNumMax, int tNumMax)
{
	//最大数の登録
	m_iPNumMax = pNumMax;
	m_iCNumMax = cNumMax;
	m_iTNumMax = tNumMax;

	//粒子情報の初期化
	m_iPtoCMax = m_iCNumMax*0.4;	//1331 layer2 0.4 layer3 0.75
									//2197 layer2 0.4 layer3 0.4 layer4 0.5
	m_iPtoTMax = m_iTNumMax*0.4;	//1331 layer2 0.3 layer3 0.5
									//2197 layer2 0.3 layer3 0.3

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
	m_iCtoPMax = m_iPNumMax*0.5;	//1331 layer2 0.5 layer3 0.75
									//2197 layer2 0.5 layre3 0.5
	m_piCtoPNum = new int[m_iCNumMax];
	m_piCtoTNum = new int[m_iCNumMax];

	m_piCtoPIndx = new int[m_iCNumMax];
	m_piCtoTIndx = new int[m_iCNumMax];

	for(int i = 0; i < m_iCNumMax; i++)
	{
		m_piCtoPNum[i] = 0;
		m_piCtoTNum[i] = 0;
		
		m_piCtoPIndx[i] = 0;
		m_piCtoTIndx[i] = 0;
	}

	//四面体情報の初期化
	//TtoPMaxは最大４で固定
	//TtoCMaxは必要ない
	m_piTtoPNum = new int[m_iTNumMax];
	m_piTtoCNum = new int[m_iTNumMax];

	m_piTtoPIndx = new int[m_iTNumMax];
	m_piTtoCIndx = new int[m_iTNumMax];

	for(int i = 0; i < m_iTNumMax; i++)
	{
		m_piTtoPNum[i] = 0;
		m_piTtoCNum[i] = 0;
		
		m_piTtoPIndx[i] = 0;
		m_piTtoCIndx[i] = 0;
	}

	//近傍四面体
	m_pppiNeighborTetra = new int**[m_iTNumMax];
	m_piNTNum = new int[m_iTNumMax];

	m_iNeighborMax = m_iTNumMax*0.3;		//1331 layer2 0.3 layer3 0.75
											//2197 layer2 0.3 layre3 0.3 layer4 0.4
											//3375 layer2
	for(int i = 0; i < m_iTNumMax; i++)
	{
		m_piNTNum[i] = 0;
		m_pppiNeighborTetra[i] = new int*[m_iNeighborMax];

		for(int j = 0; j < m_iNeighborMax; j++)
		{
			m_pppiNeighborTetra[i][j] = new int[2];
			m_pppiNeighborTetra[i][j][0] = -1;
			m_pppiNeighborTetra[i][j][1] = -1;
		}
	}

	//フラグ
	m_pbPFlag = new bool[m_iPNumMax];
	m_pbCFlag = new bool[m_iCNumMax];
	m_pbTFlag = new bool[m_iTNumMax];

	ResetPFlag(m_iPNumMax);
	ResetCFlag(m_iCNumMax);
	ResetTFlag(m_iTNumMax);
}


/*!
 * @param[in] pNum　最大粒子数
 * @param[in] cNum　最大クラスタ数
 */
IceStructure::IceStructure(int pNum, int cNum)
{
	m_iPNumMax = pNum;
	m_iCNumMax = cNum;

	//粒子
//	m_iPtoCMax = m_iCNumMax*0.5;
	m_iPtoCMax = m_iCNumMax*0.65;
	m_piPtoCNum_Connect = new int[m_iPNumMax];
	m_piPtoCNum_Calc	= new int[m_iPNumMax];

	m_piPtoCIndx_Connect	= new int[m_iPNumMax];
	m_piPtoCIndx_Calc		= new int[m_iPNumMax];

	//初期化　粒子がクラスタに所属する個数　配列の添え字
	for(int i = 0; i < m_iPNumMax; i++)
	{
		m_piPtoCNum_Connect[i] = 0;
		m_piPtoCNum_Calc[i] = 0;

		m_piPtoCIndx_Connect[i] = 0;
		m_piPtoCIndx_Calc[i] = 0;
	}

	//クラスタ
	m_iCtoPMax = m_iPNumMax*0.75;
	m_piCtoPNum_Connect = new int[m_iCNumMax];
	m_piCtoPNum_Calc	= new int[m_iCNumMax];

	m_piCtoPIndx_Connect	= new int[m_iCNumMax];
	m_piCtoPIndx_Calc		= new int[m_iCNumMax];

	//初期化　クラスタに含まれる粒子の個数　配列の添え字
	for(int i = 0; i < m_iCNumMax; i++)
	{
		m_piCtoPNum_Connect[i] = 0;
		m_piCtoPNum_Calc[i] = 0;

		m_piCtoPIndx_Connect[i] = 0;
		m_piCtoPIndx_Calc[i] = 0;
	}

	//近傍クラスタ
	m_ppiNeighborCluster = new int**[m_iCNumMax];				//接続クラスタの近傍となるクラスタを保存
	m_ppiNCNum = new int[m_iCNumMax];							//接続クラスタの近傍となるクラスタの個数
	
//	m_iNeighborMax = m_iCNumMax * 0.35;
	m_iNeighborMax = m_iCNumMax * 0.5;							//近傍クラスタの個数

	for(int i = 0; i < m_iCNumMax; i++)
	{
		m_ppiNCNum[i] = 0;
		m_ppiNeighborCluster[i] = new int*[m_iNeighborMax];
		for(int j = 0; j < m_iNeighborMax; j++)
		{
			m_ppiNeighborCluster[i][j] = new int[2];
			m_ppiNeighborCluster[i][j][0] = -1;
			m_ppiNeighborCluster[i][j][1] = -1;
		}
	}

	//フラグ
	m_pbPFlag = new bool[m_iPNumMax];
	m_pbCFlag = new bool[m_iCNumMax];

	ResetPFlag(m_iPNumMax);
	ResetCFlag(m_iCNumMax);
}

IceStructure::~IceStructure(void)
{
}

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

/*!
 * 接続情報　粒子とクラスタの領域確保
 */
void IceStructure::InitStateConnect()
{//	cout << __FUNCTION__ << endl;
 
	//粒子→クラスタ
	m_pppiPtoC_Connect = new int**[m_iPNumMax];

	for(int i = 0; i < m_iPNumMax; i++)
	{
		m_pppiPtoC_Connect[i] = new int*[m_iPtoCMax/4];

		for(int j = 0; j < m_iPtoCMax/4; j++)
		{
			m_pppiPtoC_Connect[i][j] = new int[2];				//[0] = クラスタ番号　[1] = クラスタ内での番号
			m_pppiPtoC_Connect[i][j][0] = -1;
			m_pppiPtoC_Connect[i][j][1] = -1;
		}
	}

	//クラスタ→粒子
	m_ppiCtoP_Connect = new int*[m_iCNumMax];

	for(int i = 0; i < m_iCNumMax; i++)
	{
		m_ppiCtoP_Connect[i] = new int[4];						//４で固定

		for(int j = 0; j < 4; j++)
		{
			m_ppiCtoP_Connect[i][j] = -1;
		}
	}

	//添え字の初期化
	for(int i = 0; i < m_iPNumMax; i++)
	{
		m_piPtoCIndx_Connect[i] = m_piPtoCNum_Connect[i];		//添え字の初期化
	}

	for(int i = 0; i < m_iCNumMax; i++)
	{
		m_piCtoPIndx_Connect[i] = m_piCtoPNum_Connect[i];		//添え字の初期化
	}
}


/*!
 * 計算情報　粒子とクラスタの初期化
 */
void IceStructure::InitStateCalc()
{//	cout << __FUNCTION__ << endl;
	//粒子
	m_pppiPtoC_Calc = new int**[m_iPNumMax];
	int j = 0;

	for(int i = 0; i < m_iPNumMax; i++)
	{
		m_pppiPtoC_Calc[i] = new int*[m_iPtoCMax];

		for(j = 0; j < m_iPtoCMax; j++)
		{
			m_pppiPtoC_Calc[i][j] = new int[3];
			m_pppiPtoC_Calc[i][j][0] = -1;				//[0] = クラスタ番号
			m_pppiPtoC_Calc[i][j][1] = -1;				//[1] = クラスタ内での番号
			m_pppiPtoC_Calc[i][j][2] = -1;				//[2] = layer番号
		}
	}

//	cout << __FUNCTION__ << "Particle end" << endl;

	//クラスタ
	m_ppiCtoP_Calc = new int**[m_iCNumMax];

	for(int i = 0; i < m_iCNumMax; i++)
	{
		m_ppiCtoP_Calc[i] = new int*[m_iCtoPMax];

		for(j = 0; j < m_iCtoPMax; j++)
		{
			m_ppiCtoP_Calc[i][j] = new int[2];
			m_ppiCtoP_Calc[i][j][0] = -1;
			m_ppiCtoP_Calc[i][j][1] = -1;
		}
	}

	//添え字の初期化
	for(int i = 0; i < m_iPNumMax; i++)
	{
		m_piPtoCIndx_Calc[i] = m_piPtoCNum_Calc[i];		//添え字の初期化
	}

	for(int i = 0; i < m_iCNumMax; i++)
	{
		m_piCtoPIndx_Calc[i] = m_piCtoPNum_Calc[i];		//添え字の初期化
	}
}

/*!
 * 計算情報　四面体情報の領域確保　粒子ベース処理
 */
void IceStructure::InitTetraInfo()
{
	//粒子→四面体
	m_pppiPtoT = new int**[m_iPNumMax];

	for(int i = 0; i < m_iPNumMax; i++)
	{
		m_pppiPtoT[i] = new int*[m_iPtoTMax];
		
		for(int j = 0; j < m_iPtoTMax; j++)
		{
			m_pppiPtoT[i][j] = new int[2];
			m_pppiPtoT[i][j][0] = -1;
			m_pppiPtoT[i][j][1] = -1;
		}
	}

	//四面体→粒子
	m_ppiTtoP = new int*[m_iTNumMax];

	for(int i = 0; i < m_iTNumMax; i++)
	{
		m_ppiTtoP[i] = new int[4];
		
		for(int j = 0; j < 4; j++)
		{
			m_ppiTtoP[i][j] = -1;
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
{
	//粒子→クラスタ
	m_pppiPtoC = new int**[m_iPNumMax];

	for(int i = 0; i < m_iPNumMax; i++)
	{
		m_pppiPtoC[i] = new int*[m_iPtoCMax];
		
		for(int j = 0; j < m_iPtoCMax; j++)
		{
			m_pppiPtoC[i][j] = new int[3];
			m_pppiPtoC[i][j][0] = -1;				//[0] = クラスタ番号
			m_pppiPtoC[i][j][1] = -1;				//[1] = クラスタ内での番号
			m_pppiPtoC[i][j][2] = -1;				//[2] = layer番号
		}
	}

	//クラスタ→粒子
	m_pppiCtoP = new int**[m_iCNumMax];

	for(int i = 0; i < m_iCNumMax; i++)
	{
		m_pppiCtoP[i] = new int*[m_iCtoPMax];

		for(int j = 0; j < m_iCtoPMax; j++)
		{
			m_pppiCtoP[i][j] = new int[2];
			m_pppiCtoP[i][j][0] = -1;				//[0] = 粒子番号
			m_pppiCtoP[i][j][1] = -1;				//[1] = layer番号
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
}

//-------------------------------------粒子ベース処理----------------------------------------
/*!
 * 粒子→四面体　配列の空き場所を探して添え字を返す
 */
int IceStructure::GetPtoTFreeIndx(int pIndx)
{
	int freeIndx = -1;

	for(int i = 0; i < m_iPtoTMax; i++)
	{
		if(GetPtoT(pIndx, i)[0] != -1 || GetPtoT(pIndx, i)[1] != -1){	continue;	}
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
		if(GetPtoC(pIndx, i)[0] != -1 || GetPtoC(pIndx, i)[1] != -1){	continue;	}
		freeIndx = i;	break;
	}

	if(freeIndx == m_iPtoCMax || freeIndx == -1)
	{	
		cout << __FUNCTION__ << " Error::配列に空きがありません " << endl;
		freeIndx = 0;
	}

	return freeIndx;
}

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

	m_pppiPtoT[pIndx][lIndx][0] = tIndx;
	m_pppiPtoT[pIndx][lIndx][1] = oIndx;
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
		m_ppiTtoP[tIndx][i] = pIndxList[i];
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

	m_pppiPtoC[pIndx][lIndx][0] = cIndx;
	m_pppiPtoC[pIndx][lIndx][1] = oIndx;
	m_pppiPtoC[pIndx][lIndx][2] = layer;
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
		m_pppiCtoP[cIndx][i][0] = pIndxList[i];
		m_pppiCtoP[cIndx][i][1] = pLayerList[i];
	}
}

/*!
 * 登録処理　初期の近傍四面体
 * @param[in] tIndx　四面体番号
 * @param[in] layer　探索階層
 */
void IceStructure::SetNeighborTetra(int tIndx, int layer)
{	//	cout << __FUNCTION__ << endl;
	
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
			int* jtSet = GetPtoT(ipIndx, j);

			if(jtSet[0] == -1 || jtSet[1] == -1 || jtSet[0] == tIndx){	continue;	}

			//同じクラスタを既に含んでいないかのチェック
			if(CheckNeighborTetra(tIndx, jtSet[0]) != -1){	continue;	}

			//近傍クラスタの登録＋カウント
			if(GetNTNum(tIndx) >= m_iNeighborMax)
			{
				cout << __FUNCTION__ << "近傍粒子のメモリが足りなくなりました" << endl;
				cout << GetNTNum(tIndx) << " < " << m_iNeighborMax << endl;
			}
			else
			{
				m_pppiNeighborTetra[tIndx][GetNTNum(tIndx)][0] = jtSet[0];
				m_pppiNeighborTetra[tIndx][GetNTNum(tIndx)][1] = 1;
				CountNT(tIndx);
			}
		}
	}

	//layer層目をたどって近傍四面体を取得する
	int nowSize = 0, nowIndx = 0;

	for(int i = 2; i <= layer; i++)
	{
		nowSize = GetNTNum(tIndx);

		//探索するクラスタをnowSizeとnowIndxで制限している
		for(int j = nowIndx; j < nowSize; j++)
		{
			int jtIndx = GetNeighborTetra(tIndx, j)[0];			//近傍四面体のひとつ
			
			//近傍四面体に含まれる粒子が他の四面体にも含まれている場合，その四面体を近傍として登録
			for(int k = 0; k < GetTtoPIndx(jtIndx); k++)
			{
				int kpIndx = GetTtoP(jtIndx, k);				//近傍四面体に含まれている粒子のひとつ
				if(kpIndx == -1){	continue;	}

				if(find(pIndxList.begin(), pIndxList.end(), kpIndx) != pIndxList.end())
				{	
					continue;
				}
				pIndxList.push_back(kpIndx);

				for(int l = 0; l < GetPtoTIndx(kpIndx); l++)
				{
					int* ltSet = GetPtoT(kpIndx, l);
					if(ltSet[0] == -1 || ltSet[1] == -1 || ltSet[0] == tIndx){	continue;	}

					//同じ四面体を既に含んでいないかのチェック
					if(CheckNeighborTetra(tIndx, ltSet[0]) != -1){	continue;	}

					//近傍クラスタの登録＋カウント
					if(GetNTNum(tIndx) >= m_iNeighborMax)
					{
						cout << __FUNCTION__ << " 近傍四面体のメモリが足りなくなりました " << GetNTNum(tIndx) << endl;
						cout << GetNTNum(tIndx) << " < " << m_iNeighborMax << endl;
					}
					else
					{
//						if(GetNTNum(tIndx) > 200 ) return;	//近傍四面体数の制限
						m_pppiNeighborTetra[tIndx][GetNTNum(tIndx)][0] = ltSet[0];
						m_pppiNeighborTetra[tIndx][GetNTNum(tIndx)][1] = i;
						CountNT(tIndx);
					}
				}
			}
		}
		nowIndx = nowSize;										//次のループ開始時のスタート番号を更新
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
			int* jtSet = GetPtoT(ipIndx, j);

			if(jtSet[0] == -1 || jtSet[1] == -1 || jtSet[0] == tIndx){	continue;	}

			//同じクラスタを既に含んでいないかのチェック
			if(CheckNeighborTetra(tIndx, jtSet[0]) != -1){	continue;	}

			//近傍クラスタの登録＋カウント
			if(GetNTNum(tIndx) >= m_iNeighborMax)
			{
				cout << __FUNCTION__ << "近傍粒子のメモリが足りなくなりました" << endl;
				cout << GetNTNum(tIndx) << " < " << m_iNeighborMax << endl;
			}
			else
			{
				m_pppiNeighborTetra[tIndx][GetNTNum(tIndx)][0] = jtSet[0];
				m_pppiNeighborTetra[tIndx][GetNTNum(tIndx)][1] = 1;
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
			int jtIndx = GetNeighborTetra(tIndx, j)[0];				//近傍四面体のひとつ
			
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
					int* ltSet = GetPtoT(kpIndx, l);
					if(ltSet[0] == -1 || ltSet[1] == -1 || ltSet[0] == tIndx){	continue;	}
					if(CheckNeighborTetra(tIndx, ltSet[0]) != -1){	continue;	}	//同じ四面体を既に含んでいないかのチェック

					//近傍クラスタの登録＋カウント
					if(GetNTNum(tIndx) >= m_iNeighborMax)
					{
						cout << __FUNCTION__ << " 近傍四面体のメモリが足りなくなりました " << endl;
						cout << GetNTNum(tIndx) << " < " << m_iNeighborMax << endl;
					}
					else
					{
//						if(GetNTNum(tIndx) > 200 ) return;	//近傍四面体数の制限
						m_pppiNeighborTetra[tIndx][GetNTNum(tIndx)][0] = ltSet[0];
						m_pppiNeighborTetra[tIndx][GetNTNum(tIndx)][1] = i;
						CountNT(tIndx);
					}
				}
			}
		}
		nowIndx = nowSize;											//次のループ開始時のスタート番号を更新
	}
}


/*!
 * 取得処理　粒子→四面体
 * @param[in] pIndx　粒子番号
 * @param[in] lIndx　粒子内番号
 */
int* IceStructure::GetPtoT(int pIndx, int lIndx)
{//	cout << __FUNCTION__ << endl;
	return m_pppiPtoT[pIndx][lIndx];
}

/*!
 * 取得処理　四面体→粒子
 * @param[in] tIndx　四面体番号
 * @param[in] lIndx　四面体内番号
 */
int IceStructure::GetTtoP(int tIndx, int lIndx)
{//	cout << __FUNCTION__ << endl;
	return m_ppiTtoP[tIndx][lIndx];
}

/*!
 * 取得処理　粒子→クラスタ
 * @param[in] pIndx　粒子番号
 * @param[in] lIndx　粒子内番号
 */
int* IceStructure::GetPtoC(int pIndx, int lIndx)
{
	return m_pppiPtoC[pIndx][lIndx];
}

/*!
 * 取得処理　クラスタ→粒子
 * @param[in] cIndx　粒子番号
 * @param[in] lIndx　粒子内番号
 */
int* IceStructure::GetCtoP(int cIndx, int lIndx)
{
	return m_pppiCtoP[cIndx][lIndx];
}

/*!
 * 取得処理　四面体→近傍四面体
 * @param[in] tIndx　四面体番号
 * @param[in] lIndx　四面体内番号
 */
int* IceStructure::GetNeighborTetra(int tIndx, int lIndx)
{//	cout << __FUNCTION__ << endl;
	return m_pppiNeighborTetra[tIndx][lIndx];
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

	m_ppiTtoP[tIndx][lIndx] = -1;
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

	m_pppiPtoC[pIndx][lIndx][0] = -1;
	m_pppiPtoC[pIndx][lIndx][1] = -1;
	m_pppiPtoC[pIndx][lIndx][2] = -1;
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

	m_pppiPtoT[pIndx][lIndx][0] = -1;
	m_pppiPtoT[pIndx][lIndx][1] = -1;
}

/*!
 * 初期化処理　粒子→四面体
 */
void IceStructure::ClearPtoT(int pIndx)
{
	for(int i = 0; i < m_iPtoTMax; i++)
	{
		m_pppiPtoT[pIndx][i][0] = -1;
		m_pppiPtoT[pIndx][i][1] = -1;
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
		m_pppiPtoC[pIndx][i][0] = -1;
		m_pppiPtoC[pIndx][i][1] = -1;
		m_pppiPtoC[pIndx][i][2] = -1;
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
		m_pppiCtoP[cIndx][i][0] = -1;
		m_pppiCtoP[cIndx][i][1] = -1;
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
		m_ppiTtoP[tIndx][i] = -1;
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
		m_pppiNeighborTetra[tIndx][i][0] = -1;
		m_pppiNeighborTetra[tIndx][i][1] = -1;
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
		if(m_pppiNeighborTetra[tIndx][i][1] >= layer)
		{
			m_pppiNeighborTetra[tIndx][i][0] = -1;
			m_pppiNeighborTetra[tIndx][i][1] = -1;
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
		if(checkTIndx == GetNeighborTetra(tIndx, k)[0])
		{
			findIndx = k; break;
		}
	}

	return findIndx;
}

//-------------------------------------粒子ベース処理----------------------------------------


//---------------------------------------接続情報-------------------------------------------
/*!
 * 追加処理　粒子→クラスタ
 */
void IceStructure::AddPtoC_Connect()
{

}

void IceStructure::AddPtoC_Connect(int cIndx, int oIndx)
{

}


/*!
 * 追加処理　クラスタ→粒子
 */
void IceStructure::AddCtoP_Connect(vector<int> pIndxes)
{
}

/*!
 * 削除処理　粒子→クラスタ
 */
void IceStructure::DeletePtoC_Connect(int pIndx, int coIndx)
{//	cout << __FUNCTION__ << " check1" << endl;
	if(m_pppiPtoC_Connect[pIndx][coIndx][0] == -1 || m_pppiPtoC_Connect[pIndx][coIndx][1] == -1)
	{	cout << __FUNCTION__ << " Error::-1 " << endl;
		return;
	}

	m_piPtoCNum_Connect[pIndx]--;

	if(m_piPtoCNum_Connect[pIndx] < 0)
	{
		cout << __FUNCTION__ << " Error::minus m_piPtoCNum_Connect[" << pIndx << "] = " << m_piPtoCNum_Connect[pIndx] << endl;
		m_piPtoCNum_Connect[pIndx] = 0;
		return;
	}

	//A　-1巡回方式　　早い　無駄なメモリが増えていく，巡回の場合によっては挙動がおかしくなる
	//とりあえず-1にしただけ．巡回は別の機会に．
	m_pppiPtoC_Connect[pIndx][coIndx][0] = -1;
	m_pppiPtoC_Connect[pIndx][coIndx][1] = -1;
//	cout << __FUNCTION__ << " check2" << endl;

	//B　左詰方式　効率的　コピーが発生し重くなる
}

/*!
 * 削除処理　クラスタ→粒子
 */
void IceStructure::DeleteCtoP_Connect(int cIndx, int oIndx)
{//	cout << __FUNCTION__ << " check1" << endl;
	//oIndxは添え字ではなく個数なので，探索してoIndx個目のやつを書き換える．
	int indx = 0;
	for(int i = 0; i < 4; i++)
	{
		if(m_ppiCtoP_Connect[cIndx][i] == -1){	continue;	}
		if(indx == oIndx){ indx = i;	break;	}				//oIndx個目の場合，Indxを保存してbreak
		indx++;
	}

	if(indx == 4 || m_ppiCtoP_Connect[cIndx][indx] == -1){	return;	}		//breakで抜けなかった場合のチェック


	if(m_piCtoPNum_Connect[cIndx] <= 0)
	{
		cout << __FUNCTION__ << " Error::minus m_piCtoPNum_Connect[cIndx] = " << m_piCtoPNum_Connect[cIndx] << endl;
		return;
	}
	else
	{
		m_piCtoPNum_Connect[cIndx]--;
		if(m_piCtoPNum_Connect[cIndx]==0){		m_piCtoPIndx_Connect[cIndx] = 0;	}
	}

	//A　-1巡回方式　　早い　無駄なメモリが増えていく，巡回の場合によっては挙動がおかしくなる
	m_ppiCtoP_Connect[cIndx][indx] = -1;
}

/*!
 * 変更処理　粒子→クラスタ
 * @param[in] pIndx　粒子番号
 * @param[in] lIndx　粒子が属するl番目のクラスタ
 * @param[in] cIndx　クラスタ番号
 * @param[in] oIndx　クラスタ内での順序番号
 */
void IceStructure::SetPtoC_Connect(int pIndx, int lIndx, int cIndx, int oIndx)
{//	cout << "SetParticleToCluster" << endl;
	//エラーチェック
	if(m_iPtoCMax/4 < lIndx)
	{
		cout << __FUNCTION__ << " Error::粒子が属する接続クラスタ情報を格納できなくなりました．メモリが足りません" << endl;
		cout << m_iPtoCMax/4 << "<" << lIndx << endl;
		return;
	}

	m_pppiPtoC_Connect[pIndx][lIndx][0] = cIndx;
	m_pppiPtoC_Connect[pIndx][lIndx][1] = oIndx;
}

/*!
 * 変更処理　クラスタ→粒子
 */
void IceStructure::SetCtoP_Connect(int cIndx, vector<int> pIndxList)
{
	if(4 < pIndxList.size())
	{
		cout << __FUNCTION__ << " Error::計算クラスタが含む粒子の情報を格納できなくなりました．メモリが足りません" << endl;
		cout << 4 << "<" << pIndxList.size() << endl;
		return;
	}

	for(int i = 0; i < pIndxList.size(); i++)
	{
		m_ppiCtoP_Connect[cIndx][i] = pIndxList[i];
	}
}

void IceStructure::SetPtoCIndx_Connect(int pIndx, int pNum)
{
	//エラーチェック
	if(pNum >= m_iPtoCMax/4)
	{
		cout << __FUNCTION__ << " pIndx = " << pIndx << " pNum = " << pNum << endl;
		cout << "配列に空きがなくなりました．" << endl;

		m_piPtoCIndx_Calc[pIndx] = 0;
		return;
	}
	m_piPtoCIndx_Connect[pIndx] = pNum;
}

void IceStructure::SetCtoPIndx_Connect(int cIndx, int cNum)
{
	//エラーチェック
	if(cNum > 4)
	{
		cout << __FUNCTION__ << " cIndx = " << cIndx << " cNum = " << cNum << endl;
		cout << "接続クラスタが保持できる粒子数が限界です．配列に空きがなくなりました．" << endl;
		return;
	}
	m_piCtoPIndx_Connect[cIndx] = cNum;
}

/*!
 * 取得処理　粒子→クラスタ
 */
int* IceStructure::GetPtoC_Connect(int pIndx, int lIndx)
{
	return m_pppiPtoC_Connect[pIndx][lIndx];
}


/*!
 * 取得処理　クラスタ→粒子
 */
int IceStructure::GetCtoP_Connect(int cIndx, int lIndx)
{//	cout << __FUNCTION__ << endl;
	return m_ppiCtoP_Connect[cIndx][lIndx];
}

/*!
 * 更新処理　粒子→クラスタ
 */
void IceStructure::ResetOrderPToC(int cIndx, int oIndx)
{//	cout << __FUNCTION__ << endl;
	for(int i = 0; i < m_iPNum; i++)
	{
		for(int j = 0; j < m_piPtoCIndx_Connect[i]; j++)
		{
			int* coSet = GetPtoC_Connect(i, j);
			int jcIndx = coSet[0];
			int joIndx = coSet[1];

			if(jcIndx == -1 || joIndx == -1){	continue;	}
			if(jcIndx != cIndx){				continue;	}
			if(joIndx < oIndx){					continue;	}

			coSet[1]--;
		}
	}
//	cout << __FUNCTION__ << "End" << endl;
}

/*!
 * 初期化処理　粒子→クラスタ
 */
void IceStructure::ClearPtoC_Connect(int pIndx)
{	
	for(int i = 0; i < m_iPtoCMax/4; i++)
	{
		m_pppiPtoC_Connect[pIndx][i][0] = -1;
		m_pppiPtoC_Connect[pIndx][i][1] = -1;
	}

	m_piPtoCNum_Connect[pIndx] = 0;
	m_piPtoCIndx_Connect[pIndx] = 0;
}

/*!
 * 初期化処理　クラスタ→粒子
 */
void IceStructure::ClearCtoP_Connect(int cIndx)
{
	for(int i = 0; i < 4; i++)
	{
		m_ppiCtoP_Connect[cIndx][i] = -1;
	}

	m_piCtoPIndx_Connect[cIndx] = 0;
	m_piCtoPNum_Connect[cIndx] = 0;
}

//------------------------------------計算処理クラスタ-------------------------------------
/*!
 * 追加処理　粒子→クラスタ
 */
void IceStructure::AddPtoC_Calc()
{

}

void IceStructure::AddPtoC_Calc(int cIndx, int oIndx)
{

}

/*!
 * 追加処理　クラスタ→粒子
 */
void IceStructure::AddCtoP_Calc(vector<int> pIndxes)
{

}

/*!
 * 削除処理　粒子→クラスタ
 */
void IceStructure::DeletePtoC_Calc(int pIndx, int oIndx)
{
	if(m_pppiPtoC_Calc[pIndx][oIndx][0] == -1 || m_pppiPtoC_Calc[pIndx][oIndx][1] == -1)
	{
		return;
	}

	m_piPtoCNum_Calc[pIndx]--;

	if(m_piPtoCNum_Calc[pIndx] < 0)
	{
		cout << __FUNCTION__ << " Error::minus m_piPtoCNum_Calc[pIndx] = " << m_piPtoCNum_Calc[pIndx] << endl;
		m_piPtoCNum_Calc[pIndx] = 0;
		return;
	}

	m_pppiPtoC_Calc[pIndx][oIndx][0] = -1;
	m_pppiPtoC_Calc[pIndx][oIndx][1] = -1;
	m_pppiPtoC_Calc[pIndx][oIndx][2] = -1;
}

/*!
 * 削除処理　クラスタ→粒子
 */
void IceStructure::DeleteCtoP_Calc(int cIndx, int oIndx)
{
	//oIndxは添え字ではなく個数なので，探索してoIndx個目のやつを書き換える．
	int indx = -1;
	for(int i = 0; i < GetCtoPIndx_Calc(cIndx); i++)
	{
		if(m_ppiCtoP_Calc[cIndx][i][0] == -1 || m_ppiCtoP_Calc[cIndx][i][1] == -1){	continue;	}
		indx++;
		if(indx == oIndx){ indx = i;	break;	}				//oIndx個目の場合，Indxを保存してbreak
	}

	if(m_ppiCtoP_Calc[cIndx][indx][0] == -1 || m_ppiCtoP_Calc[cIndx][indx][1] == -1){	return;	}

	m_piCtoPNum_Calc[cIndx]--;
	if(m_piCtoPNum_Calc[cIndx] < 0)
	{
		cout << __FUNCTION__ << " Error::minus m_piCtoPNum_Calc[cIndx] = " << m_piCtoPNum_Calc[cIndx] << endl;
		return;
	}

	//A　-1巡回方式　　早い　無駄なメモリが増えていく，巡回の場合によっては挙動がおかしくなる
	//この変数は実際は循環はしない
	m_ppiCtoP_Calc[cIndx][indx][0] = -1;
	m_ppiCtoP_Calc[cIndx][indx][1] = -1;
}

/*!
 * 変更処理　粒子→クラスタ
 * @param[in] pIndx　粒子番号
 * @param[in] lIndx　粒子が属するl番目のクラスタ
 * @param[in] cIndx　クラスタ番号
 * @param[in] oIndx　クラスタ内での順序番号
 * @param[in] layer  粒子の属する層番号
 */
void IceStructure::SetPtoC_Calc(int pIndx, int lIndx, int cIndx, int oIndx, int layer)
{//	cout << "SetParticleToCluster" << endl;
	if(m_iPtoCMax < lIndx)
	{	
		cout << __FUNCTION__ << " Error::粒子が属する計算クラスタ情報を格納できなくなりました．メモリが足りません" << endl;
		cout << m_iPtoCMax << "<" << lIndx << endl;
		cout << "pIndx = " << pIndx << " cIndx = " << cIndx << " oIndx = " << oIndx << " layer = " << layer << endl;
		return;
	}

	m_pppiPtoC_Calc[pIndx][lIndx][0] = cIndx;
	m_pppiPtoC_Calc[pIndx][lIndx][1] = oIndx;
	m_pppiPtoC_Calc[pIndx][lIndx][2] = layer;
}

/*!
 * 変更処理　クラスタ→粒子
 */
void IceStructure::SetCtoP_Calc(int cIndx, vector<int> pIndxList, int* pLayerList)
{
	if(m_iCtoPMax < pIndxList.size())
	{
		cout << __FUNCTION__ << " Error::計算クラスタが含む粒子の情報を格納できなくなりました．メモリが足りません" << endl;
		cout << m_iCtoPMax << "<" << pIndxList.size() << endl;
		return;
	}

	for(int i = 0; i < pIndxList.size(); i++)
	{
		m_ppiCtoP_Calc[cIndx][i][0] = pIndxList[i];
		m_ppiCtoP_Calc[cIndx][i][1] = pLayerList[i];
	}
}

void IceStructure::SetPtoCIndx_Calc(int pIndx, int pNum)
{
	//エラーチェック
	if(pNum >= m_iPtoCMax)
	{
		cout << __FUNCTION__ << " pIndx = " << pIndx << " pNum = " << pNum << endl;
		cout << "配列に空きがなくなりました．" << endl;

		m_piPtoCIndx_Calc[pIndx] = 0;
		return;
	}
	m_piPtoCIndx_Calc[pIndx] = pNum;
}

void IceStructure::SetCtoPIndx_Calc(int cIndx, int cNum)
{
	//エラーチェック
	if(cNum >= m_iCtoPMax)
	{
		cout << __FUNCTION__ << " cIndx = " << cIndx << " cNum = " << cNum << endl;
		cout << "計算クラスタが保持できる粒子数が限界です．配列に空きがなくなりました．" << endl;
		return;
	}
	m_piCtoPIndx_Calc[cIndx] = cNum;
}

/*!
 * 初期化処理　粒子→クラスタ
 */
void IceStructure::ClearPtoC_Calc(int pIndx)
{	
	for(int i = 0; i < m_iPtoCMax; i++)
	{
		m_pppiPtoC_Calc[pIndx][i][0] = -1;
		m_pppiPtoC_Calc[pIndx][i][1] = -1;
		m_pppiPtoC_Calc[pIndx][i][1] = -1;
	}

	m_piPtoCNum_Calc[pIndx] = 0;
	m_piPtoCIndx_Calc[pIndx] = 0;
}

/*!
 * 初期化処理　クラスタ→粒子
 */
void IceStructure::ClearCtoP_Calc(int cIndx)
{
	for(int i = 0; i < m_piCtoPIndx_Calc[cIndx]; i++)
	{
		m_ppiCtoP_Calc[cIndx][i][0] = -1;
		m_ppiCtoP_Calc[cIndx][i][1] = -1;
	}

	m_piCtoPIndx_Calc[cIndx] = 0;
	m_piCtoPNum_Calc[cIndx] = 0;
}

/*!
 * 取得処理　粒子→クラスタ
 */
int* IceStructure::GetPtoC_Calc(int pIndx, int lIndx)
{//	cout << __FUNCTION__ << endl;
//	cout << "m_iPtoCMax = " << m_iPtoCMax << " lIndx = " << lIndx << endl;

	if(m_iPtoCMax <= lIndx)
	{	
		cout << __FUNCTION__ << " Error::粒子が属する計算クラスタ情報を格納できなくなりました．メモリが足りません" << endl;
		cout << m_iPtoCMax << "<" << lIndx << " pIndx = " << pIndx << endl;
		return m_pppiPtoC_Calc[pIndx][0];
	}
	return m_pppiPtoC_Calc[pIndx][lIndx];
}

/*!
 * 取得処理　クラスタ→粒子
 */
int* IceStructure::GetCtoP_Calc(int cIndx, int lIndx)
{
	return m_ppiCtoP_Calc[cIndx][lIndx];
}

/*!
 * 更新処理　粒子→クラスタ
 */
void IceStructure::ResetOrderPToCCalc(int pIndx, int cIndx, int oIndx)
{

}

//-----------------------------------近傍クラスタ-------------------------------------------
/*!
 * 追加処理
 */
void IceStructure::AddNeighborCluster(int cIndx, int layer)
{
}

/*!
 * 変更処理　初期化専用　接続クラスタの近傍クラスタを登録
 */
void IceStructure::SetNeighborCluster(int cIndx, int layer)
{//	cout << __FUNCTION__ << endl;
	//初期クラスタに含まれる粒子が他のクラスタにも含まれている場合，そのクラスタを近傍として登録
	int j = 0, k = 0;

	int ipIndx = 0, kpIndx = 0;
	int jcIndx = 0;

	int* jcIndxList;
	bool bFind = false;

	for(int i = 0; i < GetCtoPIndx_Connect(cIndx); i++)
	{
		ipIndx = GetCtoP_Connect(cIndx, i);
		if(ipIndx == -1){	continue;	}

		//探索
		for(j = 0; j < GetPtoCIndx_Connect(ipIndx); j++)
		{
			jcIndxList = GetPtoC_Connect(ipIndx, j);
			if(jcIndxList[0] == -1 || jcIndxList[1] == -1){	continue;	}
			if(jcIndxList[0] == cIndx){	continue;	}

			//同じクラスタを既に含んでいないかのチェック
			bFind = false;
			for(k = 0; k < m_ppiNCNum[cIndx]; k++)
			{
				if(jcIndxList[0] == m_ppiNeighborCluster[cIndx][k][0]){ bFind = true; break;}
			}
			if( bFind == true ){	continue;	}

			//近傍クラスタの登録＋カウント
			if(m_ppiNCNum[cIndx] >= m_iNeighborMax){ cout << "Init 近傍粒子のメモリが足りなくなりました" << endl;	}
			else
			{
				m_ppiNeighborCluster[cIndx][m_ppiNCNum[cIndx]][0] = jcIndxList[0];
				m_ppiNeighborCluster[cIndx][m_ppiNCNum[cIndx]][1] = 0;
				m_ppiNCNum[cIndx]++;
			}
		}
	}

//	cout << __FUNCTION__ << "Check1" << endl;
	//並列化の名残
	int nowSize = 0, nowIndx = 0;
	int l = 0, m = 0;
	int* lcIndxList;
	
	//layer層をたどって近傍クラスタを取得する
	for(int i = 0; i < layer; i++)
	{
		nowSize = m_ppiNCNum[cIndx];

		//探索するクラスタをnowSizeとnowIndxで制限している
		for(j = nowIndx; j < nowSize; j++)
		{
			jcIndx = m_ppiNeighborCluster[cIndx][j][0];				//近傍クラスタのひとつ
			//近傍クラスタに含まれる粒子が他のクラスタにも含まれている場合，そのクラスタを近傍として登録
			for(k = 0; k < GetCtoPIndx_Connect(jcIndx); k++)
			{		
				kpIndx = GetCtoP_Connect(jcIndx, k);				//近傍クラスタに含まれている粒子のひとつ
				if(kpIndx == -1){	continue;	}

				for(l = 0; l < GetPtoCIndx_Connect(kpIndx); l++)
				{
					lcIndxList = GetPtoC_Connect(kpIndx, l);
					if(lcIndxList[0] == -1 || lcIndxList[1] == -1){	continue;	}
					if(lcIndxList[0] == cIndx){	continue;	}

					//同じクラスタを既に含んでいないかのチェック
					bFind = false;
					for(m = 0; m < m_ppiNCNum[cIndx]; m++)
					{
						if(lcIndxList[0] == m_ppiNeighborCluster[cIndx][m][0]){ bFind = true; break;}
					}
					if(bFind == true){	continue;	}

					//近傍クラスタの登録＋カウント
					if(m_ppiNCNum[cIndx] >= m_iNeighborMax){ cout << __FUNCTION__ << " 近傍粒子のメモリが足りなくなりました " << m_ppiNCNum[cIndx] << endl;	}
					else
					{
//						if(m_ppiNCNum[cIndx] > 200 ) return;	//近傍クラスタ数の制限
						m_ppiNeighborCluster[cIndx][m_ppiNCNum[cIndx]][0] = lcIndxList[0];
						m_ppiNeighborCluster[cIndx][m_ppiNCNum[cIndx]][1] = i+1;
						m_ppiNCNum[cIndx]++;
					}
				}
			}
		}
		nowIndx = nowSize;											//次のループ開始時のスタート番号を更新
	}

//	cout << __FUNCTION__ << " m_ppiNCNum[cIndx] = " << m_ppiNCNum[cIndx] << endl;
//	cout << __FUNCTION__ << "Check2" << endl;
}


/*!
 * 変更処理　現在の近傍クラスタの情報を基に，接続クラスタの近傍クラスタを登録
 * @param[in] cIndx　クラスタ番号
 * @param[in] layer　探索終了レイヤー数
 * @param[in] jlayer 探索開始レイヤー数
 */

void IceStructure::SetNeighborClusterFromCluster(int cIndx, int layer, int jlayer)
{//	cout << __FUNCTION__ << endl;

	//近傍クラスタ情報の初期化
	//jlayer以上を初期化，探索開始位置を求める
	//例，jlayer=2なら，初めてlayer=1となる場所がstartNumとなる　また，layer>=2となる場所が初期化される
	int startNum = 0;

	for(int i = m_ppiNCNum[cIndx]-1; 0 <= i; i--)
	{
		if(m_ppiNeighborCluster[cIndx][i][1] < jlayer-1)
		{
			startNum = i+1;
			break;
		}
		else if(m_ppiNeighborCluster[cIndx][i][1] >= jlayer)
		{
			m_ppiNeighborCluster[cIndx][i][0] = -1;
			m_ppiNeighborCluster[cIndx][i][1] = -1;
		}
	}

	m_ppiNCNum[cIndx] = startNum;

	//近傍クラスタの取得
	//jlayer層から開始　jlayer-1層までのクラスタは取得済み
	int j = 0, k = 0, l = 0, m = 0;
	int ipIndx = 0, kpIndx = 0;
	int jcIndx = 0;
	int nowSize = 0, nowIndx = 0;

	int* lcIndxList;

	bool bFind = false;

	//layer層をたどって近傍クラスタを取得する
	for(int i = 0; i < layer; i++)
	{
		nowSize = startNum;

		//探索するクラスタをnowSizeとnowIndxで制限している
		for(j = nowIndx; j < nowSize; j++)
		{
			jcIndx = m_ppiNeighborCluster[cIndx][j][0];				//近傍クラスタのひとつ
			//近傍クラスタに含まれる粒子が他のクラスタにも含まれている場合，そのクラスタを近傍として登録
			for(k = 0; k < GetCtoPIndx_Connect(jcIndx); k++)
			{		
				kpIndx = GetCtoP_Connect(jcIndx, k);				//近傍クラスタに含まれている粒子のひとつ
				if(kpIndx == -1){	continue;	}

				for(l = 0; l < GetPtoCIndx_Connect(kpIndx); l++)
				{
					lcIndxList = GetPtoC_Connect(kpIndx, l);
					if(lcIndxList[0] == -1 || lcIndxList[1] == -1){	continue;	}
					if(lcIndxList[0] == cIndx){	continue;	}

					//同じクラスタを既に含んでいないかのチェック
					bFind = false;
					for(m = 0; m < m_ppiNCNum[cIndx]; m++)
					{
						if(lcIndxList[0] == m_ppiNeighborCluster[cIndx][m][0]){ bFind = true; break;}
					}
					if(bFind == true){	continue;	}

					//近傍クラスタの登録＋カウント
					if(m_ppiNCNum[cIndx] >= m_iNeighborMax){ cout << __FUNCTION__ << " 近傍粒子のメモリが足りなくなりました " << m_ppiNCNum[cIndx] << endl;	}
					else
					{
						m_ppiNeighborCluster[cIndx][m_ppiNCNum[cIndx]][0] = lcIndxList[0];
						m_ppiNeighborCluster[cIndx][m_ppiNCNum[cIndx]][1] = i+1;
//						if(m_ppiNCNum[cIndx] > 200 ) return;	//近傍クラスタ数の制限
						m_ppiNCNum[cIndx]++;
					}
				}
			}
		}
		nowIndx = nowSize;											//次のループ開始時のスタート番号を更新
	}

//	cout << __FUNCTION__ << " m_ppiNCNum[cIndx] = " << m_ppiNCNum[cIndx] << endl;
}

/*!
 * 初期化処理
 */
void IceStructure::ClearNeighborCluster(int cIndx)
{//	cout << __FUNCTION__ << "check1" << endl;
	for(int i = 0; i < m_ppiNCNum[cIndx]; i++)
	{
		m_ppiNeighborCluster[cIndx][i][0] = -1;
		m_ppiNeighborCluster[cIndx][i][1] = -1;
	}

	m_ppiNCNum[cIndx] = 0;
//	cout << __FUNCTION__ << "check1" << endl;
}

/*!
 * 取得処理
 */
int* IceStructure::GetNeighborCluster(int cIndx, int lIndx)
{//	cout << __FUNCTION__ << endl;
	return m_ppiNeighborCluster[cIndx][lIndx];
}

//--------------------------------接続情報と計算処理の関連-------------------------------
/*!
 * 初期化処理
 */
void IceStructure::MakeCalcToConnect()
{
	m_pppiCalcToConnect = new int**[m_iCNumMax];
	int j = 0;

	for(int i = 0; i < m_iCNumMax; i++)
	{
//		cout << "i = " << i << " maxCnum = " << m_iCNumMax << endl;
		m_pppiCalcToConnect[i] = new int*[m_iCtoPMax];

		for(j = 0; j < m_iCtoPMax; j++)
		{
			m_pppiCalcToConnect[i][j] = new int[2];
			m_pppiCalcToConnect[i][j][0] = -1;
			m_pppiCalcToConnect[i][j][1] = -1;
		}
	}
}

/*!
 * 追加処理
 */
void IceStructure::AddCalcToConnect(int caIndx, int coIndx, int oIndx)
{

}

/*!
 * 削除処理
 */
void IceStructure::DeleteCalcToConnect(int caIndx, int poIndx)
{

}

/*!
 * 登録処理
 * @param[in] caIndx　計算クラスタ
 * @param[in] ocaIndx　計算クラスタに所属する粒子の順序番号
 * @param[in] coIndx　粒子が所属している接続番号
 * @param[in] oIndx　　粒子が接続クラスタに所属している順序番号
 */
void IceStructure::SetCalcToConnect(int caIndx, int ocaIndx, int coIndx, int oIndx)
{//	cout << __FUNCTION__ << endl;
	if(m_iCtoPMax <= ocaIndx)
	{
		cout << __FUNCTION__ << " メモリが足りません" << endl;
		cout << m_iCtoPMax << "<" << ocaIndx << endl;
		return;
	}

	m_pppiCalcToConnect[caIndx][ocaIndx][0] = coIndx;
	m_pppiCalcToConnect[caIndx][ocaIndx][1] = oIndx;
}

/*!
 * 取得処理
 */
int* IceStructure::GetCalcToConnect(int caIndx, int ocaIndx)
{//	cout << __FUNCTION__ << " caIndx = " << caIndx << " ocaIndx = " << ocaIndx;
	if(m_iCtoPMax <= ocaIndx)
	{
		cout << __FUNCTION__ << " 不正アクセスです" << endl;
		cout << m_iCtoPMax << "<" << ocaIndx << endl;
		return m_pppiCalcToConnect[0][0];
	}
//	cout << __FUNCTION__ << "End" << endl;
	return m_pppiCalcToConnect[caIndx][ocaIndx];
}

/*!
 * 更新処理
 */
void IceStructure::ResetOrderCaToCo(int caIndx, int poIndx)
{

}

/*!
 * 初期化処理
 */
void IceStructure::ClearCalcToConnect(int caIndx)
{//	cout << __FUNCTION__ << "check1" << endl;
	for(int i = 0; i < m_iCtoPMax; i++)
	{
		m_pppiCalcToConnect[caIndx][i][0] = -1;
		m_pppiCalcToConnect[caIndx][i][1] = -1;
	}
//	cout << __FUNCTION__ << "check2" << endl;
}

//----------------------------------デバッグ機能---------------------------------------------
void IceStructure::DebugPtoC_Connect(int pIndx)
{	cout << "DebugPtoC_Connect pIndx=" << pIndx;
	cout << " num=" << GetPtoCNum_Connect(pIndx) << " Indx=" << GetPtoCIndx_Connect(pIndx);
	
	for(int i = 0; i < GetPtoCIndx_Connect(pIndx); i++)
	{
			int* coSet = GetPtoC_Connect(pIndx, i);
			cout <<" c=" << coSet[0] << " o=" << coSet[1];
	}
	cout << endl;
}

void IceStructure::DebugCtoP_Connect(int cIndx)
{	cout << "DebugCtoP_Connect cIndx=" << cIndx;
	cout << " num=" << GetCtoPNum_Connect(cIndx) << " Indx=" << GetCtoPIndx_Connect(cIndx);

	for(int i = 0; i < GetCtoPIndx_Connect(cIndx); i++)
	{
		cout << " " << GetCtoP_Connect(cIndx, i);
	}
	cout << endl;
}

void IceStructure::DebugPtoC_Calc(int pIndx)
{	cout << "DebugPtoC_Calc	 pIndx=" << pIndx;
	cout << " num=" << GetPtoCNum_Calc(pIndx) << " Indx=" << GetPtoCIndx_Calc(pIndx);

	for(int i = 0; i < GetPtoCIndx_Calc(pIndx); i++)
	{
		int* coSet = GetPtoC_Calc(pIndx, i);
		cout << " c=" << coSet[0] << " o=" << coSet[1] << " L=" << coSet[2];
	}
	cout << endl;
}

void IceStructure::DebugCtoP_Calc(int cIndx)
{	cout << "DebugCtoP_Calc	cIndx=" << cIndx;
	cout << " num=" << GetCtoPNum_Calc(cIndx) << " Indx=" << GetCtoPIndx_Calc(cIndx);

	for(int i = 0; i < GetCtoPIndx_Calc(cIndx); i++)
	{
		int* coSet = GetCtoP_Calc(cIndx, i);
		cout << " p=" << coSet[0] << " L=" << coSet[1];
	}
	cout << endl;
}

void IceStructure::DebugCalcToConnect(int cIndx)
{	cout << "DebugCalcToConnect	 cIndx=" << cIndx;

	for(unsigned j = 0; j < GetCtoPIndx_Calc(cIndx); j++ )
	{
		cout << " co=" << m_pppiCalcToConnect[cIndx][j][0] << " o=" << m_pppiCalcToConnect[cIndx][j][1];
	}
	cout << endl;
}

void IceStructure::DebugNeighborCluster(int cIndx)
{	cout << "DebugNeighborCluster cIndx=" << cIndx;
	cout << " num = " << GetNCNum(cIndx);
	for(unsigned j = 0; j < GetNCNum(cIndx); j++)
	{
		cout << " NC=" << GetNeighborCluster(cIndx, j)[0] << " Ly=" << GetNeighborCluster(cIndx, j)[1];
	}
	cout << endl;
}

//--------------------------------------------粒子ベース-------------------------------------
void IceStructure::DebugPtoT(int pIndx)
{	cout << __FUNCTION__ << " pIndx = " << pIndx;
	cout << " num=" << GetPtoTNum(pIndx) << " Indx=" << GetPtoTIndx(pIndx);
	
	for(int i = 0; i < GetPtoTIndx(pIndx); i++)
	{
			int* coSet = GetPtoT(pIndx, i);
			cout <<" c=" << coSet[0] << " o=" << coSet[1];
	}
	cout << endl;
}

void IceStructure::DebugPtoC(int pIndx)
{	cout << __FUNCTION__ << " pIndx = " << pIndx;
	cout << " num=" << GetPtoCNum(pIndx) << " Indx=" << GetPtoCIndx(pIndx);
	
	for(int i = 0; i < GetPtoCIndx(pIndx); i++)
	{
			int* coSet = GetPtoC(pIndx, i);
			cout <<" c=" << coSet[0] << " o=" << coSet[1];
	}
	cout << endl;

}

void IceStructure::DebugCtoP(int cIndx)
{	cout << __FUNCTION__ << " cIndx = " << cIndx;
	cout << " num=" << GetCtoPNum(cIndx) << " Indx=" << GetCtoPIndx(cIndx);

	for(int i = 0; i < GetCtoPIndx(cIndx); i++)
	{
		int* plSet = GetCtoP(cIndx, i);
		cout <<" p=" << plSet[0] << " l=" << plSet[1];
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
		cout << " NC=" << GetNeighborTetra(tIndx, j)[0] << " Ly=" << GetNeighborTetra(tIndx, j)[1];
	}
	cout << endl;
}