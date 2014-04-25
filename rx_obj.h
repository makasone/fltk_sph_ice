/*! @file rx_obj.h
	
	@brief OBJ/MTL File Input/Output
 
	@author Makoto Fujisawa
	@date  
*/


#ifndef _RX_OBJ_H_
#define _RX_OBJ_H_


//-----------------------------------------------------------------------------
// Include Files
//-----------------------------------------------------------------------------
#include "rx_mesh.h"


//-----------------------------------------------------------------------------
// Name Space
//-----------------------------------------------------------------------------
using namespace std;


//-----------------------------------------------------------------------------
// rxOBJ�N���X�̐錾 - OBJ�`���̓ǂݍ���
//-----------------------------------------------------------------------------
class rxOBJ
{
	rxMTL m_mapMaterials;	//!< ���x���ƃf�[�^�̃}�b�v
	string m_strCurrentMat;				//!< ���݂̃f�[�^���������x��

	//vector<rxMaterialOBJ> m_vMaterials;
	//int m_iCurrentMat;

public:
	rxOBJ();
	~rxOBJ();

	bool Read(string file_name, vector<Vec3> &vrts, vector<Vec3> &vnms, vector<rxFace> &plys, rxMTL &mats, bool triangle = true);
	bool Save(string file_name, const vector<Vec3> &vrts, const vector< vector<int> >&idxes);

	rxMTL GetMaterials(void){ return m_mapMaterials; }

private:
	int loadFace(string &buf, vector<int> &vidxs, vector<int> &nidxs, vector<int> &tidxs);
	int loadMTL(const string &mtl_fn);
};




#endif // _RX_VRML_H_
