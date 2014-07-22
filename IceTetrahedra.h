//相変化オブジェクトの四面体に関するクラス
//Tetgenを用いて初期粒子配置から四面体を作成
//TODO::CGALの準備

#ifndef _ICE_TETRAHEDRA_
#define _ICE_TETRAHEDRA_

#include <iostream>
#include <fstream>
#include <string>

#include <vector>

// メッシュ化
//#include "rx_model.h"

#include "tetgen.h"

using namespace std;

#define TETGENCOMMAND "-q10.0a0.5"		//点の追加を許可するコマンド

class IceTetrahedra
{
public:
	static IceTetrahedra& GetInstance();

	void InitTetra(float* pos, int vrtxNum);

	vector<int>& GetTetraList(int listIndx){	return m_vviTetraList[listIndx];	}
	unsigned GetTetraListSize(){	return m_vviTetraList.size();	}

private:
	IceTetrahedra(){};
	IceTetrahedra(const IceTetrahedra &other){};
	IceTetrahedra &operator=(const IceTetrahedra &other){}

	void MakeTetrahedraFromCube(float* pos, int vrtxNum);							//立方体のための処理
	void MakeTetrahedraRectParallele(float* pos, int vrtxNum, int x, int y, int z);	//未完成　直方体のための処理　TODO::立方体と統合したい
	void MakeTetrahedraOnlySurface(float* pos,int vrtxNum);
	//void MakeTetrahedraFromObj(rxPolygons poly,int vrtxNum);

	void Load_obj_File(const string name);
	void Load_ELE_File(const string name);
	void Load_NODE_File(const string name, float* p);

	void Save_POLY_File(const string name, float* pos, int vrtxNum);
	void Save_NODE_File(const string name, float* pos, int vrtxNum);

private:
	vector<vector<int>> m_vviTetraList;		//四面体の組み合わせリスト

};

#endif