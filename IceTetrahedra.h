//���ω��I�u�W�F�N�g�̎l�ʑ̂Ɋւ���N���X
//Tetgen��p���ď������q�z�u����l�ʑ̂��쐬
//TODO::CGAL�̏���

#ifndef _ICE_TETRAHEDRA_
#define _ICE_TETRAHEDRA_

#include <iostream>
#include <fstream>
#include <string>

#include <vector>

#include "tetgen.h"

//#include <cassert>
//#include <list>
//#include <vector>

//���܂������Ȃ������H
//#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
//#include <CGAL/Triangulation_3.h>
//#include <CGAL/Simple_cartesian.h>
//#include <CGAL/Polyhedron_3.h>
//#include <CGAL/IO/Polyhedron_iostream.h>
//#include <CGAL/IO/Polyhedron_VRML_1_ostream.h>

//typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
//typedef CGAL::Triangulation_3<K> Triangulation;
//typedef CGAL::Simple_cartesian<double> Kernel;
//typedef CGAL::Polyhedron_3<Kernel> Polyhedron;

//typedef Triangulation::Cell_handle Cell_handle;
//typedef Triangulation::Vertex_handle Vertex_handle;
//typedef Triangulation::Locate_type Locate_type;
//typedef Triangulation::Point Point;

using namespace std;

#define TETGENCOMMAND "-q10.0a0.5"		//�_�̒ǉ���������R�}���h
#define OBJ_NAME	"obj/bunny4913.obj"
//#define OBJ_NAME	"obj/slime.obj"
#define ELE_FILE	"obj/bunny1331.ele"
#define NODE_FILE	"obj/voxel.node"

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

	void MakeTetrahedraFromCube(float* pos, int vrtxNum);							//�����̂̂��߂̏���
	void MakeTetrahedraRectParallele(float* pos, int vrtxNum, int x, int y, int z);	//�������@�����̂̂��߂̏����@TODO::�����̂Ɠ���������
	void MakeTetrahedraOnlySurface(float* pos,int vrtxNum);
	//void MakeTetrahedraFromObj(rxPolygons poly,int vrtxNum);

	void Load_obj_File(const string name, float* pos);
	void Load_ELE_File(const string name);
	void Load_NODE_File(const string name, float* p);
	void LoadTest_VoxelFile(const string name, float* pos);

	void Save_POLY_File(const string name, float* pos, int vrtxNum);
	void Save_NODE_File(const string name, float* pos, int vrtxNum);

private:
	vector<vector<int>> m_vviTetraList;		//�l�ʑ̂̑g�ݍ��킹���X�g

};

#endif