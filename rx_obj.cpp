/*! @file rx_obj.cpp
	
	@brief OBJ/MTL File Input/Output
 
	@author Makoto Fujisawa
	@date  
*/


//-----------------------------------------------------------------------------
// Include Files
//-----------------------------------------------------------------------------
#include "rx_obj.h"



//-----------------------------------------------------------------------------
// rxVRMLクラスの実装
//-----------------------------------------------------------------------------
/*!
 * コンストラクタ
 */
rxOBJ::rxOBJ(void)
{
	m_strCurrentMat = "";
}

/*!
 * デストラクタ
 */
rxOBJ::~rxOBJ()
{
}

/*!
 * OBJファイル読み込み
 * @param[in] file_name ファイル名(フルパス)
 * @param[out] vrts 頂点座標
 * @param[out] vnms 頂点法線
 * @param[out] txcs 頂点テクスチャ座標
 * @param[out] tris ポリゴン頂点インデックス
 * @param[out] mats ポリゴン材質
 * @param[in] triangle ポリゴンの三角形分割フラグ
 * @return 
 */
bool rxOBJ::Read(string file_name, vector<Vec3> &vrts, vector<Vec3> &vnms, vector<rxFace> &plys, rxMTL &mats, bool triangle)
{
	ifstream file;

	file.open(file_name.c_str());
	if(!file || !file.is_open() || file.bad() || file.fail()){
		cout << "rxOBJ::Read : Invalid file specified" << endl;
		return false;
	}

	vector<Vec3> vnormals;
	vector<Vec2> vtexcoords;

	string buf;
	string::size_type comment_start = 0;
	while(!file.eof()){
		getline(file, buf);

		// '#'以降はコメントとして無視
		if( (comment_start = buf.find('#')) != string::size_type(-1) )
			buf = buf.substr(0, comment_start);

		// 行頭のスペース，タブを削除
		DeleteSpace(buf);

		// 空行は無視
		if(buf.empty())
			continue;

		if(buf[0] == 'v'){
			if(buf[1] == 'n'){		// 頂点法線
				Vec3 n;
				if(StringToVec3s(buf, "vn", n)){
					vnormals.push_back(n);
				}
#if RX_DEBUG_OUT
				cout << "vn " << n << endl;
#endif
			}
			else if(buf[1] == 't'){	// テクスチャ座標
				Vec2 tc;
				if(StringToVec2s(buf, "vt", tc)){
					vtexcoords.push_back(tc);
				}
#if RX_DEBUG_OUT
				cout << "vt " << tc << endl;
#endif
			}
			else{					// 頂点座標
				Vec3 v;
				if(StringToVec3s(buf, "v", v)){
					vrts.push_back(v);
					vnms.push_back(Vec3(0.0));
				}
#if RX_DEBUG_OUT
				cout << "v " << v << endl;
#endif
			}
		}
		else if(buf[0] == 'f'){		// 面
			int num_face = 0;
			vector<int> vidxs, nidxs, tidxs;
			if(!(num_face = loadFace(buf, vidxs, nidxs, tidxs))) continue;

#if RX_DEBUG_OUT
			cout << "f ";
			for(int i = 0; i < (int)vidxs.size(); ++i){
				cout << vidxs[i]+1 << " ";
			}
			cout << endl;
#endif

			if(triangle && num_face >= 4){
				PolyToTri(plys, vidxs, tidxs, vtexcoords, m_strCurrentMat);
			}
			else{
				rxFace face;
				face.vert_idx.resize(num_face);;
				face.material_name = m_strCurrentMat;
				face.texcoords.resize(num_face);
				bool tc = !vtexcoords.empty();
				for(int i = 0; i < num_face; ++i){
					face[i] = vidxs[i];
					if(tc){
						face.texcoords[i] = vtexcoords[tidxs[i]];
					}
					else{
						face.texcoords[i] = Vec2(0.0);
					}
				}
				plys.push_back(face);
			}

			// 頂点法線
			if(!vnormals.empty()){
				for(int i = 0; i < num_face; ++i){
					vnms[vidxs[i]] += vnormals[nidxs[i]];
				}
			}
		}
		else if(buf[0] == 'u'){		// 材質名
			string mat_name;
			if(!StringToString(buf, "usemtl", mat_name)) continue;

			if(mat_name == "(null)" || m_mapMaterials.empty()) continue;

			map<string, rxMaterialOBJ>::iterator i = m_mapMaterials.find(mat_name);

			if(i != m_mapMaterials.end()){
				m_strCurrentMat = mat_name;
			}

			//cout << "usemtl " << mat_name << endl;
		}
		else if(buf[0] == 'm'){		// 材質情報ファイル
			string mat_file_name;
			if(!StringToString(buf, "mtllib", mat_file_name)) continue;

			string dir = ExtractDirPath(file_name);
			if(!dir.empty()){
				mat_file_name = dir+"\\"+mat_file_name;
			}
			//cout << "mtllib " << mat_file_name << endl;

			if(loadMTL(mat_file_name) == 0){
				//cout << "fail to read " << mat_file_name << endl;
				continue;
			}

		}
	}

	if(vnormals.empty()){
		vnms.clear();
	}

	if(!vnms.empty()){
		for(int i = 0; i < (int)vnms.size(); ++i){
			normalize(vnms[i]);
		}
	}

	mats = m_mapMaterials;

	return true;
}

// OBJ形式で保存
bool rxOBJ::Save(string file_name, const vector<Vec3> &vrts, const vector< vector<int> >&idxes)
{
	if((int)idxes.size() == 0) return false;

	return true;
}


int rxOBJ::loadFace(string &buf, vector<int> &vidxs, vector<int> &nidxs, vector<int> &tidxs)
{
	int num_face = CountString(buf, 2, " ")+1;

	if(num_face >= 3){
		vidxs.resize(num_face);
		nidxs.resize(num_face);
		tidxs.resize(num_face);

		int vidx = 0;	// 頂点インデックス
		int tidx = 0;	// テクスチャ座標インデックス
		int nidx = 0;	// 法線インデックス
		int offset = 2;
		for(int i = 0; i < num_face; ++i){
			if(sscanf(&buf[0]+offset, "%d/%d/%d", &vidx, &tidx, &nidx) != 3){
				if(sscanf(&buf[0]+offset, "%d//%d", &vidx, &nidx) != 2){
					if(sscanf(&buf[0]+offset, "%d/%d", &vidx, &tidx) != 2){
						sscanf(&buf[0]+offset, "%d", &vidx);	// 頂点座標のみ
						tidx = 0;
						nidx = 0;
						offset += (int)log10((double)vidx)+1;
					}
					else{	// 頂点座標，テクスチャ座標
						nidx = 0;
						offset += (int)log10((double)vidx)+(int)log10((double)tidx)+3;
					}
				}
				else{	// 頂点座標，法線
					tidx = 0;
					offset += (int)log10((double)vidx)+(int)log10((double)nidx)+4;
				}
			}
			else{ // 頂点座標，テクスチャ座標，法線
				offset += (int)log10((double)vidx)+(int)log10((double)tidx)+1+(int)log10((double)nidx)+4;
			}

			offset++;	// スペース

			if(vidx >= 1){
				vidxs[i]   = vidx-1;
				nidxs[i] = nidx-1;
				tidxs[i]  = tidx-1;
			}
			else{
				return 0;
			}
		}
	}

	return num_face;
}

int rxOBJ::loadMTL(const string &mtl_fn)
{
	ifstream file;

	file.open(mtl_fn.c_str());
	if(!file || !file.is_open() || file.bad() || file.fail()){
		cout << "rxOBJ::loadMTL : Invalid file specified" << endl;
		return false;
	}

	int num_mat = 0;
	rxMaterialOBJ *cur_mat;

	string buf;
	string::size_type comment_start = 0;
	while(!file.eof()){
		getline(file, buf);

		// '#'以降はコメントとして無視
		if( (comment_start = buf.find('#')) != string::size_type(-1) )
			buf = buf.substr(0, comment_start);

		// 行頭のスペース，タブを削除
		DeleteSpace(buf);

		// 空行は無視
		if(buf.empty())
			continue;

		if(buf[0] == 'n'){
			string mat_name;
			if(!StringToString(buf, "newmtl", mat_name)) continue;

			rxMaterialOBJ mat;
			mat.name      = mat_name;
			mat.diffuse   = Vec4(0.0);
			mat.specular  = Vec4(0.0);
			mat.ambient   = Vec4(0.0);
			mat.color     = Vec4(0.0);
			mat.emission  = Vec4(0.0);
			mat.shininess = 0.0;
			mat.illum = 2;
			mat.tex_file = "";
			mat.tex_name = 0;

			m_mapMaterials[mat_name] = mat;
			cur_mat = &m_mapMaterials[mat_name];

			num_mat++;
		}
		else if(buf[0] == 'K'){
			Vec3 col;
			if(buf[1] == 'd'){		// 拡散反射成分
				if(StringToVec3s(buf, "Kd", col)){
					cur_mat->diffuse = Vec4(col, cur_mat->diffuse[3]);
				}
			}
			else if(buf[1] == 's'){	// 鏡面反射成分
				if(StringToVec3s(buf, "Ks", col)){
					cur_mat->specular = Vec4(col, cur_mat->specular[3]);
				}
			}
			else if(buf[1] == 'a'){	// 環境光反射成分
				if(StringToVec3s(buf, "Ka", col)){
					cur_mat->ambient = Vec4(col, cur_mat->ambient[3]);
				}
			}
		}
		else if(buf[0] == 'd'){	// アルファ値(透過率)
			double alpha;
			if(StringToDouble(buf, "d", alpha)){
				cur_mat->diffuse[3]  = alpha;
				cur_mat->specular[3] = alpha;
				cur_mat->ambient[3]  = alpha;
				cur_mat->color[3]    = alpha;
				cur_mat->emission[3] = alpha;
			}
		}
		else if(buf[0] == 'T' && buf[1] == 'r'){	// 発光色
			Vec3 emit;
			if(StringToVec3s(buf, "Tr", emit)){
				cur_mat->emission = Vec4(emit, cur_mat->emission[3]);
			}
		}
		else if(buf[0] == 'N' && buf[1] == 's'){	// Shininess
			double shine = 0.0;
			if(StringToDouble(buf, "Ns", shine)){
				cur_mat->shininess = shine;
			}
		}
		else if(buf[0] == 'i' && buf[1] == 'l' && buf[2] == 'l'){	// 証明モデル(1で鏡面反射無効, 2で有効)
			string illum;
			if(StringToString(buf, "illum", illum)){
				cur_mat->illum = atoi(illum.c_str());
			}
		}
		else if(buf[0] == 'm' && buf[1] == 'a' && buf[2] == 'p'){	// テクスチャ名
			string texfn;
			if(StringToString(buf, "map_Kd", texfn)){
				string dir = ExtractDirPath(mtl_fn);
				if(!dir.empty()){
					texfn = dir+"\\"+texfn;
				}

				cur_mat->tex_file = texfn;
			}
		}
	}

	if(!m_mapMaterials.empty()){
		map<string, rxMaterialOBJ>::iterator iter = m_mapMaterials.begin();
		for(; iter != m_mapMaterials.end(); ++iter){
			cout << "material : " << iter->first << endl;
			cout << "  diffuse   = " << iter->second.diffuse << endl;
			cout << "  specular  = " << iter->second.specular << endl;
			cout << "  ambient   = " << iter->second.ambient << endl;
			cout << "  color     = " << iter->second.color << endl;
			cout << "  emission  = " << iter->second.emission << endl;
			cout << "  shininess = " << iter->second.shininess << endl;
			cout << "  illum     = " << iter->second.illum << endl;
			cout << "  tex_file  = " << iter->second.tex_file << endl;
		}

		m_strCurrentMat = m_mapMaterials.begin()->first;
	}

	return num_mat;
}