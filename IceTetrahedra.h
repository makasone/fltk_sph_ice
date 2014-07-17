//���ω��I�u�W�F�N�g�̎l�ʑ̂Ɋւ���N���X
//Tetgen��p���ď������q�z�u����l�ʑ̂��쐬

#ifndef _ICE_TETRAHEDRA_
#define _ICE_TETRAHEDRA_

#include <iostream>
#include <fstream>
#include <string>

#include <vector>

// ���b�V����
//#include "rx_model.h"

#include "tetgen.h"

using namespace std;

#define TETGENCOMMAND "-q10.0a0.5"		//�_�̒ǉ���������R�}���h

class IceTetrahedra
{
public:
	static IceTetrahedra& GetInstance();

	void InitTetra(float* pos, int vrtxNum);

	vector<int>& GetTetraList(int listIndx){	return m_vviTetraList[listIndx];	}
	unsigned GetTetraListSize(){	return m_vviTetraList.size();	}

	//�{���ɃV���O���g�����̃e�X�g
	void Test1(int num)
	{
		m_num = num;
		cout << __FUNCTION__ << " " << m_num << endl;
	}
	
	void Test2()
	{
		cout << __FUNCTION__ << " " << m_num << endl;
	}
private:
	IceTetrahedra(){};
	IceTetrahedra(const IceTetrahedra &other){};
	IceTetrahedra &operator=(const IceTetrahedra &other){}

	void MakeTetrahedraFromCube(float* pos, int vrtxNum);							//�����̂̂��߂̏���
	void MakeTetrahedraRectParallele(float* pos, int vrtxNum, int x, int y, int z);	//�������@�����̂̂��߂̏����@TODO::�����̂Ɠ���������
	void MakeTetrahedraOnlySurface(float* pos,int vrtxNum);
	//void MakeTetrahedraFromObj(rxPolygons poly,int vrtxNum);

	void Load_ELE_File(const string name);
	void Load_NODE_File(const string name, float* p);

	void Save_POLY_File(const string name, float* pos, int vrtxNum);
	void Save_NODE_File(const string name, float* pos, int vrtxNum);

private:
	int m_num;

	vector<vector<int>> m_vviTetraList;		//�l�ʑ̂̑g�ݍ��킹���X�g
};

#endif