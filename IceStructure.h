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

#include <iostream>
#include <vector>
#include <algorithm>

//#include "rx_matrix.h"
#include <UtilityScript\mk_Vector2D.h>
#include <UtilityScript\mk_Vector3D.h>

using namespace std;

class IceStructure
{
public:
	IceStructure(int pNum, int cNum, int tNum);
	~IceStructure(void);

	void SetParticleNum(int pNum){	m_iPNum = pNum; }		//現在の粒子数
	void SetClusterNum(int cNum){	m_iCNum = cNum; }		//現在のクラスタ数
	void SetTetraNum(int tNum){		m_iTNum = tNum;	}		//現在の四面体数

	int GetParticleNum(void){		return m_iPNum;	}
	int GetClusterNum(void){		return m_iCNum;	}
	int GetTetraNum(void){			return m_iTNum;	}

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


	int  GetTtoP(int tIndx, int lIndx);

	//修正
	int GetPtoC(int pIndx, int lIndx, int oIndx);
	int GetPtoT(int pIndx, int lIndx, int oIndx);
	int GetCtoP(int cIndx, int lIndx, int oIndx);

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
	
	int m_iPtoCMax;								//粒子がクラスタに所属する最大数　connectはcalcの半分でいい
	int m_iCtoPMax;								//クラスタが含む粒子の最大数　　　connectは最大４で固定	

	int m_iTNumMax;
	int m_iTNum;

	int m_iPtoTMax;
	int m_iTtoPMax;

	int m_iCNumMax;								//最大クラスタ数
	int m_iCNum;								//現在のクラスタ数

	int m_iNeighborMax;							//近傍粒子の最大数

	//とりあえず，ポインタは使わない
	//粒子→
	mk_Vector3D<int> m_mk3DiPtoC;				//粒子→クラスタ
	mk_Vector3D<int> m_mk3DiPtoT;				//粒子→四面体

	int*   m_piPtoCNum;							//粒子→クラスタの個数
	int*   m_piPtoTNum;							//粒子→四面体の個数

	int*   m_piPtoCIndx;
	int*   m_piPtoTIndx;

	//クラスタ→
	mk_Vector3D<int> m_mk3DiCtoP;				//クラスタ→粒子　クラスタに所属する粒子を返す　粒子の接続情報
	mk_Vector3D<int> m_mk3DiCtoT;				//クラスタ→四面体				

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