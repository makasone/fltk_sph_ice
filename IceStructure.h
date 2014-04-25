/*!
  @file IceStructure.h
	
  @brief 点，辺，面を用いて氷構造を操作するクラス
 
  @author Ryou Nakasone
  @date 2013-07
*/
//---------------------------------------------------------------------------
//このファイルのクラス・関数が二回以上コンパイルされるのを防ぐための処理（インクルードガード）
#ifndef IceStructure_H
#define IceStructure_H
//---------------------------------------------------------------------------

#include <stdio.h>
#include <iostream>
#include <vector>
#include <map>
#include <algorithm>
#include <omp.h>
#include "rx_matrix.h"


using namespace std;

class IceStructure
{
public:
	IceStructure(int pNumMax, int cNumMax, int tNumMax);				//四面体ベース
	IceStructure(int pNum, int cNum);
	~IceStructure(void);

	void SetParticleNum(int pNum){	m_iPNum = pNum; }		//現在の粒子数
	void SetClusterNum(int cNum){	m_iCNum = cNum; }		//現在のクラスタ数
	void SetTetraNum(int tNum){		m_iTNum = tNum;	}		//現在の四面体数

	int GetParticleNum(void){		return m_iPNum;	}
	int GetClusterNum(void){		return m_iCNum;	}
	int GetTetraNum(void){			return m_iTNum;	}

//-------------------------------------粒子ベース処理----------------------------------------
	void InitTetraInfo();											//四面体情報のメモリ確保
	void InitClusterInfo();											//クラスタ情報のメモリ確保

	int  GetPtoCMax(void){	return m_iPtoCMax;	}

	void CountPtoC(int pIndx){	m_piPtoCNum[pIndx]++;	}
	void CountPtoT(int pIndx){	m_piPtoTNum[pIndx]++;	}
	void CountCtoP(int cIndx){	m_piCtoPNum[cIndx]++;	}
	void CountCtoT(int cIndx){	}
	void CountTtoP(int tIndx){	m_piTtoPNum[tIndx]++;	}
	void CountTtoC(int tIndx){	}

	void CountNT(int tIndx){	m_piNTNum[tIndx]++;		}

	void SetPtoCIndx(int pIndx, int indx){	m_piPtoCIndx[pIndx] = indx;	} 
	void SetCtoPIndx(int cIndx, int indx){	m_piCtoPIndx[cIndx] = indx;	}
	void SetPtoTIndx(int pIndx, int indx){	m_piPtoTIndx[pIndx] = indx; }
	void SetTtoPIndx(int tIndx, int indx){	m_piTtoPIndx[tIndx] = indx; }

	int	 GetPtoTNum(int pIndx){		return m_piPtoTNum[pIndx];	}
	int	 GetPtoCNum(int pIndx){		return m_piPtoCNum[pIndx];	}
	int  GetTtoPNum(int tIndx){		return m_piTtoPNum[tIndx];	}
	int  GetCtoPNum(int cIndx){		return m_piCtoPNum[cIndx];	}

	int  GetNTNum(int tIndx){		return m_piNTNum[tIndx];	}

	int  GetPtoTIndx(int pIndx){	return m_piPtoTIndx[pIndx];	}
	int  GetPtoCIndx(int pIndx){	return m_piPtoCIndx[pIndx];	}
	int  GetTtoPIndx(int tIndx){	return m_piTtoPIndx[tIndx];	}
	int  GetCtoPIndx(int cIndx){	return m_piCtoPIndx[cIndx];	}

	int  GetPtoTFreeIndx(int pIndx);
	int  GetPtoCFreeIndx(int pIndx);

	void SetPtoT(int pIndx, int lIndx, int tIndx, int oIndx);				//粒子が属する四面体の登録　　粒子番号，粒子内順序，四面体番号，四面体内順序
	void SetPtoC(int pIndx, int lIndx, int cIndx, int oIndx, int layer);	//粒子が属するクラスタの登録　粒子番号，粒子内順序，クラスタ番号，クラスタ内順序
	
	void SetCtoP(int cIndx, const vector<int>& pIndxList, int* pLayerList);	//クラスタに属する粒子の登録　クラスタ番号，粒子リスト
	void SetTtoP(int tIndx, vector<int>& pIndxList);						//四面体に属する粒子の登録　　四面体番号，粒子リスト
	
	void SetNeighborTetra(int tIndx, int layer);
	void SetNeighborTetraFromLayer(int tIndx, int searchLayer, int deleteLayer);

	void DeleteTtoP(int tIndx, int lIndx);
	void DeletePtoC(int pIndx, int lIndx);
	void DeletePtoT(int pIndx, int lIndx);

	int* GetPtoT(int pIndx, int lIndx);
	int* GetPtoC(int pIndx, int lIndx);
	int* GetCtoP(int cIndx, int lIndx);
	int  GetTtoP(int tIndx, int lIndx);

	int* GetNeighborTetra(int tIndx, int lIndx);

	void ClearPtoT(int pIndx);
	void ClearPtoC(int pIndx);
	void ClearCtoP(int cIndx);
	void ClearTtoP(int tIndx);

	void ClearNeighborTetra(int tIndx);
	void ClearNeighborTetraFromLayer(int tIndx, int layer);

	int  CheckNeighborTetra(int tIndx, int checkTIndx);

	//デバッグ
	void DebugPtoT(int pIndx);
	void DebugPtoC(int pIndx);
	void DebugCtoP(int cIndx);
	void DebugTtoP(int tIndx);
	void DebugNeighborTetra(int tIndx);

//-------------------------------------粒子ベース処理----------------------------------------


//-------------------------------------四面体ベース処理-------------------------------------
	void InitStateConnect();								//接続情報のサイズ確保
	void InitStateCalc();									//計算情報のサイズ確保

	//接続情報用アクセッサ
	void AddPtoC_Connect();
	void AddPtoC_Connect(int cIndx, int oIndx);
	void AddCtoP_Connect(vector<int> pIndxes);

	void DeletePtoC_Connect(int pIndx, int coIndx);
	void DeleteCtoP_Connect(int cIndx, int oIndx);

	void SetPtoC_Connect(int pIndx, int lIndx, int cIndx, int oIndx);			//粒子番号，クラスタ番号，順序番号　接続情報の作成
	void SetCtoP_Connect(int cIndx, vector<int> pIndxList);						//書き換えるクラスタ番号，粒子番号のリスト

	void SetPtoCNum_Connect(int pIndx, int pNum){	m_piPtoCNum_Connect[pIndx] = pNum;	}
	void SetCtoPNum_Connect(int cIndx, int cNum){	m_piCtoPNum_Connect[cIndx] = cNum;	}

	void SetPtoCIndx_Connect(int pIndx, int pNum);
	void SetCtoPIndx_Connect(int cIndx, int cNum);

	void CountPtoC_Connect(int pIndx){ m_piPtoCNum_Connect[pIndx]++; }			//粒子が所属している接続クラスタのカウンタ
	void CountCtoP_Connect(int cIndx){ m_piCtoPNum_Connect[cIndx]++; }			//接続クラスタに所属している粒子のカウンタ

	int* GetPtoC_Connect(int pIndx, int lIndx);
	int GetCtoP_Connect(int cIndx, int lIndx);

	int GetPtoCNum_Connect(int pIndx){ return m_piPtoCNum_Connect[pIndx]; }		//SM法と連携するためにも重要になる
	int GetCtoPNum_Connect(int cIndx){ return m_piCtoPNum_Connect[cIndx]; }		//SM法と連携するためにも重要になる

	int GetPtoCIndx_Connect(int pIndx){ return m_piPtoCIndx_Connect[pIndx];	}
	int GetCtoPIndx_Connect(int cIndx){ return m_piCtoPIndx_Connect[cIndx];	}

	void ClearPtoC_Connect(int pIndx);
	void ClearCtoP_Connect(int cIndx);

	void ResetOrderPToC(int cIndx, int oIndx);

	//計算処理クラスタ用アクセッサ
	void AddPtoC_Calc();
	void AddPtoC_Calc(int cIndx, int oIndx);
	void AddCtoP_Calc(vector<int> pIndxes);

	void DeletePtoC_Calc(int pIndx, int coIndx);
	void DeleteCtoP_Calc(int cIndx, int oIndx);

	void SetPtoC_Calc(int pIndx, int lIndx, int cIndx, int oIndx, int layer);	//粒子番号，クラスタ番号，順序番号　計算処理クラスタの情報
	void SetCtoP_Calc(int cIndx, vector<int> pIndxList, int* pLayerList);		//書き換えるクラスタ番号，粒子番号のリスト

	void SetPtoCNum_Calc(int pIndx, int pNum){	m_piPtoCNum_Calc[pIndx] = pNum;	}
	void SetCtoPNum_Calc(int cIndx, int cNum){	m_piCtoPNum_Calc[cIndx] = cNum;	}

	void SetPtoCIndx_Calc(int pIndx, int pNum);
	void SetCtoPIndx_Calc(int cIndx, int cNum);

	void CountPtoC_Calc(int pIndx){ m_piPtoCNum_Calc[pIndx]++; }				//粒子が所属している計算クラスタのカウンタ
	void CountCtoP_Calc(int cIndx){ m_piCtoPNum_Calc[cIndx]++; }				//計算クラスタに所属している粒子のカウンタ

	int* GetPtoC_Calc(int pIndx, int lIndx);
	int* GetCtoP_Calc(int cIndx, int lIndx);

	int GetPtoCNum_Calc(int pIndx){	return m_piPtoCNum_Calc[pIndx]; }
	int GetCtoPNum_Calc(int cIndx){ return m_piCtoPNum_Calc[cIndx]; }

	int GetPtoCIndx_Calc(int pIndx){ return m_piPtoCIndx_Calc[pIndx];	}
	int GetCtoPIndx_Calc(int cIndx){ return m_piCtoPIndx_Calc[cIndx];	}

	void ClearPtoC_Calc(int pIndx);
	void ClearCtoP_Calc(int cIndx);

	void ResetOrderPToCCalc(int pIndx, int cIndx, int oIndx);

	//近傍クラスタ用アクセッサ
	void AddNeighborCluster(int cIndx, int layer);							//近傍クラスタの追加　ある近傍クラスタリストのデータを追加
	void SetNeighborCluster(int cIndx, int layer);							//近傍クラスタの修正　ある近傍クラスタリストのデータを書き換え
	void SetNeighborClusterFromCluster(int cIndx, int layer, int jlayer);	//近傍クラスタの修正　近傍クラスタの持つ情報を利用した書き換え
	int* GetNeighborCluster(int cIndx, int lIndx);							//近傍クラスタ番号　クラスタの近傍クラスタのリストを返す
	int GetNCNum(int cIndx){ return m_ppiNCNum[cIndx]; }					//近傍クラスタの数
	void ClearNeighborCluster(int cIndx);

	////計算処理クラスタ→接続情報クラスタ
	void MakeCalcToConnect(void);
	void AddCalcToConnect(int caIndx, int coIndx, int oIndx);
	void DeleteCalcToConnect(int caIndx, int poIndx);
	int* GetCalcToConnect(int cIndx, int lIndx);
	void SetCalcToConnect(int caIndx, int ocaIndx, int coIndx, int oIndx);
	void ResetOrderCaToCo(int cIndx, int poIndx);
	void ClearCalcToConnect(int caIndx);

	//デバッグ
	void DebugPtoC_Connect(int pIndx);
	void DebugCtoP_Connect(int cIndx);

	void DebugPtoC_Calc(int pIndx);
	void DebugCtoP_Calc(int cIndx);

	void DebugCalcToConnect(int cIndx);

	void DebugNeighborCluster(int cIndx);
//-------------------------------------四面体ベース処理-------------------------------------

	//フラグ　未使用
	void SetPFlag(int indx, bool state){	m_pbPFlag[indx] = state;	}
	void SetCFlag(int indx, bool state){	m_pbCFlag[indx] = state;	}
	void SetTFlag(int indx, bool state){	m_pbTFlag[indx] = state;	}

	bool GetPFlag(int indx){	return m_pbPFlag[indx];	}
	bool GetCFlag(int indx){	return m_pbCFlag[indx];	}
	bool GetTFlag(int indx){	return m_pbTFlag[indx];	}

	void ResetPFlag(int endNum);
	void ResetCFlag(int endNum);
	void ResetTFlag(int endNum);

protected:

	int m_iPNumMax;								//最大粒子数
	int m_iPNum;								//現在の粒子数

	int m_iNeighborMax;							//近傍粒子の最大数
	
	int m_iPtoCMax;								//粒子がクラスタに所属する最大数　connectはcalcの半分でいい
	int m_iCtoPMax;								//クラスタが含む粒子の最大数　　　connectは最大４で固定	
//-------------------------------------粒子ベース処理----------------------------------------
	int m_iTNumMax;
	int m_iTNum;

	int m_iPtoTMax;
	int m_iTtoPMax;

	//粒子→
	int*** m_pppiPtoC;							//粒子→クラスタ　何番目のクラスタ内で何番目なのかを判定　0番から始まるのに注意
	int*** m_pppiPtoT;							//粒子→四面体

	int*   m_piPtoCNum;							//粒子→クラスタの個数
	int*   m_piPtoTNum;							//粒子→四面体の個数

	int*   m_piPtoCIndx;
	int*   m_piPtoTIndx;

	//クラスタ→
	int***  m_pppiCtoP;							//クラスタ→粒子　クラスタに所属する粒子を返す　粒子の接続情報
	int*** m_pppiCtoT;							//クラスタ→四面体

	int*   m_piCtoPNum;							//クラスタ→粒子の個数
	int*   m_piCtoTNum;							//クラスタ→四面体の個数

	int*   m_piCtoPIndx;
	int*   m_piCtoTIndx;

	//四面体→
	int**  m_ppiTtoP;							//四面体→粒子
	int*** m_pppiTtoC;							//四面体→クラスタ

	int*   m_piTtoPNum;							//四面体→粒子の個数
	int*   m_piTtoCNum;							//四面体→クラスタの個数

	int*   m_piTtoPIndx;
	int*   m_piTtoCIndx;

	//近傍四面体
	int*** m_pppiNeighborTetra;					//近傍四面体
	int*   m_piNTNum;							//各近傍四面体の個数
//-------------------------------------粒子ベース処理----------------------------------------

//-------------------------------------四面体ベース処理-------------------------------------
	int m_iCNumMax;								//最大クラスタ数
	int m_iCNum;								//現在のクラスタ数

	//接続情報
	int*** m_pppiPtoC_Connect;					//粒子→クラスタ　何番目のクラスタ内で何番目なのかを判定　0番から始まるのに注意
	int**  m_ppiCtoP_Connect;					//クラスタ→粒子　クラスタに所属する粒子を返す　粒子の接続情報

	int* m_piPtoCNum_Connect;					//粒子が接続クラスタに所属している個数
	int* m_piCtoPNum_Connect;					//クラスタが粒子を含んでいる個数

	int* m_piPtoCIndx_Connect;					//粒子が接続クラスタに所属していることを保存する配列の，現在の添え字番号
	int* m_piCtoPIndx_Connect;					//接続クラスタが粒子を含んでいることを保存する配列の，現在の添え字番号

	//計算処理クラスタ
	int*** m_pppiPtoC_Calc;						//粒子→クラスタ　何番目のクラスタ内で何番目なのかを判定　0番から始まるのに注意
	int*** m_ppiCtoP_Calc;						//クラスタ→粒子　クラスタに所属する粒子を返す　粒子の接続情報

	int* m_piPtoCNum_Calc;						//粒子が計算クラスタに所属している個数
	int* m_piCtoPNum_Calc;						//クラスタが粒子を含んでいる個数

	int* m_piPtoCIndx_Calc;						//粒子が計算クラスタに所属していることを保存する配列の，現在の添え字番号
	int* m_piCtoPIndx_Calc;						//計算クラスタが粒子を含んでいることを保存する配列の，現在の添え字番号

	//近傍クラスタ
	int*** m_ppiNeighborCluster;				//近傍クラスタの組　それぞれに粒子のリストを用意
	int*  m_ppiNCNum;							//接続クラスタの近傍となるクラスタの個数

	//計算処理クラスタ→接続情報クラスタ
	int*** m_pppiCalcToConnect;					//計算処理クラスタ→接続情報クラスタ
//-------------------------------------四面体ベース処理-------------------------------------

	//探索用フラグ　未使用
	bool* m_pbPFlag;							//粒子
	bool* m_pbCFlag;							//クラスタ
	bool* m_pbTFlag;							//四面体

	vector<vector<int>>	m_vviEdgeToCluster;					//辺→クラスタ	ある辺がどのクラスタに所属しているか　必ず１個以上
	vector<vector<int>>	m_vviEdgeToParticle;				//辺→粒子　　　ある辺がどの粒子を接続しているか　　　必ず２個

	vector<vector<int>>	m_vviClusterToEdge;					//クラスタ→辺　あるクラスタにどの辺が所属しているか　必ず１個か３個
	vector<vector<int>>	m_vviParticleToEdge;				//粒子→辺　　　ある粒子はどの辺で接続されているか

};
#endif