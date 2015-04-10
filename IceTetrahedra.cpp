#include "IceTetrahedra.h"
#include "rx_model.h"

IceTetrahedra &IceTetrahedra::GetInstance() {
    static IceTetrahedra instance;
    return instance;
}

/*!
 * 四面体情報の初期化
 */
void IceTetrahedra::InitTetra(float* pos, int vertexNum)
{	cout << __FUNCTION__ << endl;

	//Load_ELE_File(ELE_FILE);					//うまくいかない　eleファイルを読み込みリスト作成
	//Load_obj_File(OBJ_NAME, pos);				//objファイルを読み込みリスト作成
	Load_NODE_File(NODE_FILE, pos);			//バニーモデルを用いる場合に使う

	//MakeTetrahedraFromCube(pos, vertexNum);
}

/*!
 * tetgenを用いて点群の位置情報，メッシュ情報から初期の四面体構造を作成.
 * ダンプしたファイルは src/fltk_sph_turb/bin　に作られる．
 */
void IceTetrahedra::MakeTetrahedraFromCube(float* pos, int vrtxNum)
{	cout << __FUNCTION__ << endl;

	tetgenio in, out;	// 入力メッシュと出力四面体メッシュ

	// ポリゴン頂点インデックスのスタート(0スタートか1スタート)
	in.firstnumber = 1;

	// メッシュ頂点の設定
	in.numberofpoints = vrtxNum;
	in.pointlist = new double[in.numberofpoints*3];
	for(int i = 0; i < in.numberofpoints; ++i){
		for(int j = 0; j < 3; ++j){
			//粒子位置情報の登録
			in.pointlist[3*i+j] = (double)pos[4*i+j];
		}
	}

	// ポリゴンの設定
	int n = pow( in.numberofpoints, 1.0/3.0 ) + 0.5;	//立方体の１辺の頂点数
	if( n == 1 ){	cout << "error::n == 1" << endl;	}

	//各ポリゴンの頂点番号を計算
	//格子状に配置された粒子で作られる面を作っている．
	//格子の内部に存在する面も考慮する必要があるのでややこしくなっている．
	vector<vector<int>> poligonList;
	int poligonNum = 0;

	//上下面
	for( int k = 0; k < n; k++ ){
		for( int i = 1; i < n; i++ ){
			for( int j = 1; j < n; j++ ){
				vector<int> list;
				list.push_back(k*n*n+(i-1)*n+j);
				list.push_back(k*n*n+(i-1)*n+j+1);
				list.push_back(k*n*n+(i-1)*n+j+1+n);
				list.push_back(k*n*n+(i-1)*n+j+n);
				poligonList.push_back(list);
	}}}

	//左右面
	for( int k = 0; k < n; k++ ){
		for( int i = 1; i < n; i++ ){
			for( int j = 1; j < n; j++ ){
				vector<int> list;
				list.push_back(k+(i-1)*n*n+(j-1)*n+1);
				list.push_back(k+(i-1)*n*n+(j-1)*n+1+n);
				list.push_back(k+i*n*n+(j-1)*n+1+n);
				list.push_back(k+i*n*n+(j-1)*n+1);
				poligonList.push_back(list);
	}}}

	//前後面
	for( int k = 0; k < n; k++ ){
		for( int i = 1; i < n; i++ ){
			for( int j = 1; j < n; j++ ){
				vector<int> list;
				list.push_back(k*n+(i-1)*n*n+j);
				list.push_back(k*n+(i-1)*n*n+j+1);
				list.push_back(k*n+(i-1)*n*n+j+1+n*n);
				list.push_back(k*n+(i-1)*n*n+j+n*n);
				poligonList.push_back(list);
	}}}

	in.numberoffacets = (n-1) * (n-1) * n * 3;				//面の総数 (一辺の頂点数-1)の2乗*(一辺の頂点数)*（軸の数）
	in.facetlist = new tetgenio::facet[in.numberoffacets];
	in.facetmarkerlist = new int[in.numberoffacets];

	for(int i = 0; i < in.numberoffacets; ++i){
		tetgenio::facet *f = &in.facetlist[i];

		// ポリゴンリスト
		f->numberofpolygons = 1;
		f->polygonlist = new tetgenio::polygon[f->numberofpolygons];

		// hole(穴)リスト
		f->numberofholes = 0;
		f->holelist = NULL;

		// ポリゴン頂点インデックス
		tetgenio::polygon *p = &f->polygonlist[0];

		p->numberofvertices = 4;
		p->vertexlist = new int[p->numberofvertices];
		for(int j = 0; j < p->numberofvertices; ++j)
		{
			p->vertexlist[j] = poligonList[i][j];
		}

		in.facetmarkerlist[i] =  1;	//??
	}

	// 入力メッシュ情報をファイルにダンプ
	//in.save_poly("test_poly");	// test.poly
	
	// 第一引数は，"p":PLC読み込み，"q":quality mesh generation(qの後にquality boundを数値で指定)，
	// "a":最大体積制限(aの後に体積を数値で指定)
	tetgenbehavior b;
//	if(!b.parse_commandline("-pq1.414a0.1n")){
	if(!b.parse_commandline("-pn")){
		terminatetetgen(0);
	}

	tetrahedralize(&b, &in, &out);				// 四面体メッシュ生成

	//std::cout << " the number of node   : "    << out.numberofpoints << endl;
	//std::cout << " the number of element   : " << out.numberoftetrahedra << endl;
	//std::cout << " the number of face  : "     << out.numberoftrifaces << endl;

	// 出力メッシュ情報をファイルにダンプ
//	out.save_elements("test_out");	// .ele

	// 出力四面体頂点番号を配列に格納
	int nelem = out.numberoftetrahedra;
	int nc = out.numberofcorners;
	for(int i = 0; i < nelem; ++i)
	{
		vector<int> list;
		for(int j = 0; j < nc; ++j)
		{
			int pIndx = out.tetrahedronlist[nc*i+j]-out.firstnumber;
			list.push_back(pIndx);
		}
		m_vviTetraList.push_back( list );
	}

	cout << "m_vviTetraList.size() = " << m_vviTetraList.size() << endl;
//	delete[] in.pointlist;	//必要？
}

/*!
 * tetgenを用いて点群の位置情報，メッシュ情報から初期の四面体構造を作成.
 * 直方体　未完成
 * ダンプしたファイルは src/fltk_sph_turb/bin　に作られる．
 */
void IceTetrahedra::MakeTetrahedraRectParallele(float* pos, int vrtxNum, int x, int y, int z)
{	cout << __FUNCTION__ << endl;
	
	tetgenio in, out;	// 入力メッシュと出力四面体メッシュ

	// ポリゴン頂点インデックスのスタート(0スタートか1スタート)
	in.firstnumber = 1;

	// メッシュ頂点の設定
//	in.numberofpoints = m_pPS->GetNumParticles();
	in.numberofpoints = vrtxNum;
	in.pointlist = new double[in.numberofpoints*3];
	for(int i = 0; i < in.numberofpoints; ++i){
		for(int j = 0; j < 3; ++j){
			//粒子位置情報の登録
			in.pointlist[3*i+j] = (double)pos[4*i+j];
		}
	}

	//各ポリゴンの頂点番号を計算
	//格子状に配置された粒子で作られる面を作っている．
	//格子の内部に存在する面も考慮する必要があるのでややこしくなっている．
	vector<vector<int>> poligonList;
	int poligonNum = 0;

	//前後面
	//Ｘ軸側から見て，X軸に増加→Y軸に増加→Z軸に増加
	//四角形が敷き詰められた層がいくつあるのか，と考えるのでxx=1,yy=1,zz=0
	for(int zz = 0; zz < z; zz++){
		for(int yy = 1; yy < y; yy++){
			for(int xx = 1; xx < x; xx++){
				vector<int> list;
				int Layer = zz*x*y;
				int Heigh = (yy-1)*x;
				int Line  = xx;

				//層＋高さ＋四角
				list.push_back(Layer+Heigh+Line);
				list.push_back(Layer+Heigh+Line+1);
				list.push_back(Layer+Heigh+Line+1+x);
				list.push_back(Layer+Heigh+Line+x);
				poligonList.push_back(list);
	}}}

	//左右面
	//Ｚ軸側から見て，Y軸に増加→Z軸に増加→X軸に増加
	for(int xx = 0; xx < x; xx++){
		for(int zz = 1; zz < z; zz++){
			for(int yy = 1; yy < y; yy++){
				vector<int> list;
				int Layer = xx;
				int Heigh = (zz-1)*x*y;
				int Line  = 1+(yy-1)*x;

				//層＋高さ＋四角 xLineで基点を決めて，四角形を作る
				list.push_back(Layer+Heigh+Line );
				list.push_back(Layer+Heigh+Line +x);
				list.push_back(Layer+Heigh+Line +x*y+x);
				list.push_back(Layer+Heigh+Line +x*y);
				poligonList.push_back(list);
	}}}

	////上下面
	for(int yy = 0; yy < y; yy++){
		for(int zz = 1; zz < z; zz++){
			for(int xx = 1; xx < x; xx++){
				vector<int> list;
				int Layer = yy*x;
				int Heigh = (zz-1)*x*y;
				int Line  = xx;

				//層＋高さ＋四角 xLineで基点を決めて，四角形を作る
				list.push_back(Layer+Heigh+Line );
				list.push_back(Layer+Heigh+Line +1);
				list.push_back(Layer+Heigh+Line +x*y+1);
				list.push_back(Layer+Heigh+Line +x*y);
				poligonList.push_back(list);
	}}}

	in.numberoffacets = (x-1) * (y-1) * z + (y-1) * (z-1) * x + (x-1) * (z-1) * y;	//面の総数
	in.facetlist = new tetgenio::facet[in.numberoffacets];
	in.facetmarkerlist = new int[in.numberoffacets];

	for(int i = 0; i < in.numberoffacets; ++i){
		tetgenio::facet *f = &in.facetlist[i];

		// ポリゴンリスト
		f->numberofpolygons = 1;
		f->polygonlist = new tetgenio::polygon[f->numberofpolygons];

		// hole(穴)リスト
		f->numberofholes = 0;
		f->holelist = NULL;

		// ポリゴン頂点インデックス
		tetgenio::polygon *p = &f->polygonlist[0];

		p->numberofvertices = 4;
		p->vertexlist = new int[p->numberofvertices];
		for(int j = 0; j < p->numberofvertices; ++j)
		{
			p->vertexlist[j] = poligonList[i][j];
		}

		in.facetmarkerlist[i] =  1;	//??
	}

	//入力情報をsファイルにダンプ
	//in.save_poly("Input_RectParallele");
	//in.save_nodes("Input_RectParallele");

	// 第一引数は，"p":PLC読み込み，"q":quality mesh generation(qの後にquality boundを数値で指定)，
	// "a":最大体積制限(aの後に体積を数値で指定)
	tetgenbehavior b;
//	if(!b.parse_commandline("-pq1.414a0.1n")){
	if(!b.parse_commandline("-pn")){
//	if(!b.parse_commandline(TETGENCOMMAND)){
		terminatetetgen(0);
	}

	tetrahedralize(&b, &in, &out);				// 四面体メッシュ生成

	std::cout << " the number of node   : "    << out.numberofpoints << endl;
	std::cout << " the number of element   : " << out.numberoftetrahedra << endl;
	std::cout << " the number of face  : "     << out.numberoftrifaces << endl;

	// 出力メッシュ情報をファイルにダンプ
	//out.save_elements("Output_RectParallele");	// .ele

	// 出力四面体頂点番号を配列に格納
	int nelem = out.numberoftetrahedra;
	int nc = out.numberofcorners;
	for(int i = 0; i < nelem; ++i)
	{
		vector<int> list;
		for(int j = 0; j < nc; ++j)
		{
			int pIndx = out.tetrahedronlist[nc*i+j]-out.firstnumber;
			list.push_back(pIndx);
		}
		m_vviTetraList.push_back( list );
	}
}

/*!
 * tetgenを用いて点群の位置情報，メッシュ情報から初期の四面体構造を作成.
 * ファイルは src/fltk_sph_turb/bin　に作られる．
 */
void IceTetrahedra::MakeTetrahedraOnlySurface(float* pos,int vrtxNum)
{	cout << __FUNCTION__ << endl;
	tetgenio in, out;	// 入力メッシュと出力四面体メッシュ

	// ポリゴン頂点インデックスのスタート(0スタートか1スタート)
	in.firstnumber = 1;

	// メッシュ頂点の設定
	in.numberofpoints = vrtxNum;
	in.pointlist = new double[in.numberofpoints*3];
	for(int i = 0; i < in.numberofpoints; ++i){
		for(int j = 0; j < 3; ++j){
			in.pointlist[3*i+j] = (double)pos[4*i+j];
		}
	}

	// ポリゴンの設定
	int n = pow( in.numberofpoints, 1.0/3.0 ) + 0.5;	//立方体の１辺の頂点数
	if( n == 1 ){	cout << "error::n == 1" << endl;	}

	//各ポリゴンの頂点番号を計算
	vector<vector<int>> poligonList;
	int poligonNum = 0;

	//上下面
	for( int k = 0; k < n; k++ ){
		for( int i = 1; i < n; i++ ){
			for( int j = 1; j < n; j++ ){
				vector<int> list;
				list.push_back(k*n*n+(i-1)*n+j);
				list.push_back(k*n*n+(i-1)*n+j+1);
				list.push_back(k*n*n+(i-1)*n+j+1+n);
				list.push_back(k*n*n+(i-1)*n+j+n);
				poligonList.push_back(list);
	}}}

	//左右面
	for( int k = 0; k < n; k++ ){
		for( int i = 1; i < n; i++ ){
			for( int j = 1; j < n; j++ ){
				vector<int> list;
				list.push_back(k+(i-1)*n*n+(j-1)*n+1);
				list.push_back(k+(i-1)*n*n+(j-1)*n+1+n);
				list.push_back(k+i*n*n+(j-1)*n+1+n);
				list.push_back(k+i*n*n+(j-1)*n+1);
				poligonList.push_back(list);
	}}}

	//前後面
	for( int k = 0; k < n; k++ ){
		for( int i = 1; i < n; i++ ){
			for( int j = 1; j < n; j++ ){
				vector<int> list;
				list.push_back(k*n+(i-1)*n*n+j);
				list.push_back(k*n+(i-1)*n*n+j+1);
				list.push_back(k*n+(i-1)*n*n+j+1+n*n);
				list.push_back(k*n+(i-1)*n*n+j+n*n);
				poligonList.push_back(list);
	}}}

	in.numberoffacets = (n-1) * (n-1) * n * 3;
	in.facetlist = new tetgenio::facet[in.numberoffacets];
	in.facetmarkerlist = new int[in.numberoffacets];

	for(int i = 0; i < in.numberoffacets; ++i){
		tetgenio::facet *f = &in.facetlist[i];

		// ポリゴンリスト
		f->numberofpolygons = 1;
		f->polygonlist = new tetgenio::polygon[f->numberofpolygons];

		// hole(穴)リスト
		f->numberofholes = 0;
		f->holelist = NULL;

		// ポリゴン頂点インデックス
		tetgenio::polygon *p = &f->polygonlist[0];

		p->numberofvertices = 4;
		p->vertexlist = new int[p->numberofvertices];
		for(int j = 0; j < p->numberofvertices; ++j)
		{
			p->vertexlist[j] = poligonList[i][j];
		}

		in.facetmarkerlist[i] =  1;	//??
	}

	// 入力メッシュ情報をファイルにダンプ
	//in.save_poly("test_poly");	// test.poly
	
	// 第一引数は，"p":PLC読み込み，"q":quality mesh generation(qの後にquality boundを数値で指定)，
	// "a":最大体積制限(aの後に体積を数値で指定)
	tetgenbehavior b;
//	if(!b.parse_commandline("-pq1.414a0.1n")){
	if(!b.parse_commandline("-pn")){
		terminatetetgen(0);
	}

	tetrahedralize(&b, &in, &out);				// 四面体メッシュ生成

	//std::cout << " the number of node   : "    << out.numberofpoints << endl;
	//std::cout << " the number of element   : " << out.numberoftetrahedra << endl;
	//std::cout << " the number of face  : "     << out.numberoftrifaces << endl;

	// 出力メッシュ情報をファイルにダンプ
//	out.save_elements("test_out");	// .ele

	// 出力四面体頂点番号を配列に格納
	int nelem = out.numberoftetrahedra;
	int nc = out.numberofcorners;
	for(int i = 0; i < nelem; ++i)
	{
		vector<int> list;
		for(int j = 0; j < nc; ++j)
		{
			int pIndx = out.tetrahedronlist[nc*i+j]-out.firstnumber;
			list.push_back(pIndx);
		}
		m_vviTetraList.push_back( list );
	}
}

///*!
// * tetgenを用いて点群の位置情報，メッシュ情報から初期の四面体構造を作成.
// * ファイルは src/fltk_sph_turb/bin　に作られる．
// */
//void IceTetrahedra::MakeTetrahedraFromObj(rxPolygons poly,int vrtxNum)
//{	cout << __FUNCTION__ << endl;
//
//	tetgenio in, out;	// 入力メッシュと出力四面体メッシュ
//
//	// ポリゴン頂点インデックスのスタート(0スタートか1スタート)
//	in.firstnumber = 0;
//
//	// メッシュ頂点の設定
//	in.numberofpoints = vrtxNum;
//	in.pointlist = new double[in.numberofpoints*3];
//	for(int i = 0; i < in.numberofpoints; ++i){
//		for(int j = 0; j < 3; ++j){
//			//粒子位置情報の登録
//			in.pointlist[i*3+j] = poly.vertices[i][j];
//		}
//	}
//
//	// ポリゴンの設定
//	//ポリゴンの読み込み情報を使う
//	vector<vector<int>> poligonList;
//	int poligonNum = 0;
//
//	for(int i = 0; i < poly.faces.size(); i++)
//	{
//		vector<int> list;
//		for(int j = 0; j < 3; j++)
//		{
//			list.push_back(poly.faces[i][j]);
//		}
//		poligonList.push_back(list);
//	}
//
//	in.numberoffacets = poly.faces.size();				//面の総数 (一辺の頂点数-1)の2乗*(一辺の頂点数)*（軸の数）
//	in.facetlist = new tetgenio::facet[in.numberoffacets];
//	in.facetmarkerlist = new int[in.numberoffacets];
//
//	for(int i = 0; i < in.numberoffacets; ++i)
//	{
//		tetgenio::facet *f = &in.facetlist[i];
//
//		// ポリゴンリスト
//		f->numberofpolygons = 1;
//		f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
//
//		// hole(穴)リスト
//		f->numberofholes = 0;
//		f->holelist = NULL;
//
//		// ポリゴン頂点インデックス
//		tetgenio::polygon *p = &f->polygonlist[0];
//
//		p->numberofvertices = 3;
//		p->vertexlist = new int[p->numberofvertices];
//		for(int j = 0; j < p->numberofvertices; ++j)
//		{
//			p->vertexlist[j] = poligonList[i][j];
//		}
//
//		in.facetmarkerlist[i] =  1;	//??
//	}
//
//	// 入力メッシュ情報をファイルにダンプ
//	in.save_poly("test_poly");	// 不完全らしい
//	in.save_nodes("test_node");
//	//in.save_faces("test_face");
//
//	// 第一引数は，"p":PLC読み込み，"q":quality mesh generation(qの後にquality boundを数値で指定)，
//	// "a":最大体積制限(aの後に体積を数値で指定)
//	tetgenbehavior b;
////	if(!b.parse_commandline("-pq1.414a0.1n")){
////	if(!b.parse_commandline("-pn")){
//	if(!b.parse_commandline(TETGENCOMMAND)){
//		terminatetetgen(0);
//	}
//
//	tetrahedralize(&b, &in, &out);				// 四面体メッシュ生成
//
//	std::cout << " the number of node   : "    << out.numberofpoints << endl;
//	std::cout << " the number of element   : " << out.numberoftetrahedra << endl;
//	std::cout << " the number of face  : "     << out.numberoftrifaces << endl;
//
//	// 出力メッシュ情報をファイルにダンプ
//	out.save_elements("test_out");	// .ele
//
//	// 出力四面体頂点番号を配列に格納
//	int nelem = out.numberoftetrahedra;
//	int nc = out.numberofcorners;
//	for(int i = 0; i < nelem; ++i)
//	{
//		vector<int> list;
//		for(int j = 0; j < nc; ++j)
//		{
//			int pIndx = out.tetrahedronlist[nc*i+j]-out.firstnumber;
//			list.push_back(pIndx);
//		}
//		m_vviTetraList.push_back( list );
//	}
//}

/*!
 * tetgenで得られたファイルを読み込み，初期状態を作成．四面体情報．
 * ファイルは src/fltk_sph_turb/bin　から読み込む．
 */
void IceTetrahedra::Load_ELE_File(const string name)
{	cout << __FUNCTION__ << endl;

	//ファイルを読み込み，四面体となる点の組み合わせをListに入れる．
	ifstream ifs( name );
	string str;

	//ファイルの存在確認
	if(ifs.fail()) 
	{
		cerr << "File do not exist.\n";
		exit(0);
	}

	//変数の用意，初期化
	int a=0, b=0, c=0, d=0, e=0, f=0;
	bool line_1 = false;

	m_vviTetraList.clear();

	//文字列の読み込み
	while( getline(ifs, str) )
	{
		a=0; b=0; c=0; d=0; e=0;
		//無理やりだけどとりあえず
		if( !line_1 )
		{
			line_1 = true;
			sscanf(str.data(), "%d %d %d", &a, &b, &c);
			
			//cout << "a = " << a << "\t";
			//cout << "b = " << b << "\t";
			//cout << "c = " << c << endl;
		}
		else
		{
			if( str[0] == '#' )
			{
				continue;
			}

			sscanf(str.data(), "%d %d %d %d %d", &a, &b, &c, &d, &e);
			//cout << "a = " << a << "\t";
			//cout << "b = " << b << "\t";
			//cout << "c = " << c << "\t";
			//cout << "d = " << d << "\t";
			//cout << "e = " << e << endl;

			vector<int> list;
			list.push_back(b);
			list.push_back(c);
			list.push_back(d);
			list.push_back(e);

			m_vviTetraList.push_back( list );
		}
	}

	//デバッグ
	cout << "m_vviTetraList.size() = " << m_vviTetraList.size() << endl;
}

/*!
 * tetgenで得られたファイルを読み込み，初期状態を作成．四面体情報．
 * ファイルは src/fltk_sph_turb/bin　から読み込む．
 */
void IceTetrahedra::Load_NODE_File(const string name, float* p)
{	cout << __FUNCTION__ << endl;

	//ファイルを読み込み，四面体となる点の組み合わせをListに入れる．
	ifstream ifs( name );
	string str;

	//ファイルの存在確認
	if(ifs.fail()) 
	{
		cerr << "File do not exist.\n";
		exit(0);
	}

	//変数の用意，初期化
	int ia=0, ib=0, ic=0, id=0;
	double da = 0.0, db = 0.0, dc = 0.0, dd = 0.0;
	bool line_1 = false;

	//文字列の読み込み
	while( getline(ifs, str) )
	{
		ia=0; ib=0; ic=0; id=0;
		//無理やりだけどとりあえず
		if( !line_1 )
		{
			line_1 = true;
			sscanf(str.data(), "%d %d %d %d", &ia, &ib, &ic, &id);
			
			//cout << "ia = " << ia << "\t" << "ib = " << ib << "\t" << "ic = " << ic << endl;
		}
		else
		{
			if( str[0] == '#' )
			{
				continue;
			}

			sscanf(str.data(), "%lf %lf %lf %lf", &da, &db, &dc, &dd);
			//cout << "a = " << a << "\t";
			//cout << "b = " << b << "\t";
			//cout << "c = " << c << "\t";
			//cout << "d = " << d << "\t";
			//cout << "e = " << e << endl;

			if(int(da) > 6525)	break;
			int pIndx = (int(da)-1)*4;	//ファイルによっては-1することもある
			
			double radius = 0.75f;
			p[pIndx+0] = db * radius;
			p[pIndx+1] = dc * radius - 0.25;	//TODO:ここで初期位置を調整しないといけない…
			p[pIndx+2] = dd * radius;
		}
	}
}

/*!
 * 初期状態をtetgen用に変換したファイルを作成．頂点情報，面情報．
 * ファイルは src/fltk_sph_turb/bin　に作成される．
 */
void IceTetrahedra::Save_POLY_File(const string name, float* pos, int vrtxNum)
{	cout << __FUNCTION__ << endl;

	//ファイルを作成し，フォーマットにしたがって位置情報と面情報を書き込む　#include <fstream>の位置に注意
	ofstream ofs( name+".poly" );

	//Polyファイル用テスト
//１　頂点の位置情報
	//頂点全体の情報　頂点数，次元数（３で固定），attribute，boundarymark．
	ofs << "# Part 1 - node list" << endl;
	ofs << "       " << vrtxNum << " ";
	ofs << "3 ";
	ofs << "0 ";
	ofs << "0 ";
	ofs << endl;

	//頂点それぞれの情報
	for(int i = 0; i < vrtxNum; i++ )
	{
		ofs << "       " << i+1 << " " << pos[4*i+0] << " " << pos[4*i+1] << " " << pos[4*i+2] << endl;
	}
//２　面の接続情報
	//頂点で作成される面の情報
	ofs << "# Part 2 - facet list" << endl;

	int n = pow( vrtxNum, 1.0/3.0 ) + 0.5;	//立方体の１辺の頂点数
	if( n == 1 )
	{
		cout << "error::n == 1" << endl;
	}

	////表面粒子のみバージョン
//	ofs << "\t" << (n-1) * (n-1) * 6 << "\t";		//初期状態は必ず立方体，とした場合の面の数
	ofs << "\t" << (n-1) * (n-1) * n * 3 << "\t";		//初期状態は必ず立方体，とした場合の面の数
	ofs << "\t" << 1 << "\t";
	ofs << endl;

	//上下面
	for( int k = 0; k < n; k++ )
	{
		for( int i = 1; i < n; i++ )
		{
			for( int j = 1; j < n; j++ )
			{
				ofs << "\t" << 1 << "\t" << 0 << "\t" << 1 << endl;		//パラメータ
				ofs << "\t" << 4					<< "\t"				//面の数
							<< k*n*n+(i-1)*n+j		<< "\t" 
							<< k*n*n+(i-1)*n+j+1	<< "\t"
							<< k*n*n+(i-1)*n+j+1+n	<< "\t"
							<< k*n*n+(i-1)*n+j+n
				<< endl;
			}
		}
	}

	//左右面
	for( int k = 0; k < n; k++ )
	{
		for( int i = 1; i < n; i++ )
		{
			for( int j = 1; j < n; j++ )
			{
				ofs << "\t" << 1 << "\t" << 0 << "\t" << 1 << endl;		//パラメータ
				ofs << "\t" << 4						<< "\t"			//面の数
							<< k+(i-1)*n*n+(j-1)*n+1	<< "\t" 
							<< k+(i-1)*n*n+(j-1)*n+1+n	<< "\t"
							<< k+i*n*n+(j-1)*n+1+n		<< "\t"
							<< k+i*n*n+(j-1)*n+1
				<< endl;
			}
		}
	}

	//前後面
	for( int k = 0; k < n; k++ )
	{
		for( int i = 1; i < n; i++ )
		{
			for( int j = 1; j < n; j++ )
			{
				ofs << "\t" << 1 << "\t" << 0 << "\t" << 1 << endl;		//パラメータ
				ofs << "\t" << 4						<< "\t"			//面の数
							<< k*n+(i-1)*n*n+j			<< "\t" 
							<< k*n+(i-1)*n*n+j+1		<< "\t"
							<< k*n+(i-1)*n*n+j+1+n*n	<< "\t"
							<< k*n+(i-1)*n*n+j+n*n
				<< endl;
			}
		}
	}

//３　
	ofs << "# Part 3 - hole list" << endl;
	ofs << "0";
	ofs << endl;

//４
	ofs << "# Part 4 - region list" << endl;
	ofs << "0";
	ofs << endl;
}

/*!
 * 初期状態をtetgen用に変換したファイルを作成．頂点情報のみ．
 * ファイルは src/fltk_sph_turb/bin　に作成される．
 */
void IceTetrahedra::Save_NODE_File(const string name, float* pos, int vrtxNum)
{	cout << __FUNCTION__ << endl;

	//ファイルを作成し，フォーマットにしたがって位置情報を書き込む　#include <fstream>の位置に注意
	std::ofstream ofs( name+".node" );

	//頂点全体の情報　頂点数，次元数（３で固定），attribute，boundarymark．
	ofs << "       " << vrtxNum << " ";
	ofs << "3 ";
	ofs << "0 ";
	ofs << "0 ";
	ofs << endl;

	//頂点それぞれの情報
	for(int i = 0; i < vrtxNum; i++ )
	{
		ofs << "       " << i+1 << " " << pos[4*i+0] << " " << pos[4*i+1] << " " << pos[4*i+2] << endl;
	}
}

//objファイルを読み込んで粒子位置を初期化
void IceTetrahedra::Load_obj_File(const string objName, float* pos)
{	cout << __FUNCTION__ << endl;

	//objファイル
	rxPolygons m_poly;

	//OpenFile(name);						//モデル読み込み　パスはRelease/binのやつ
	RxModel::Read(objName, m_poly);

	int vertex_count = 0;
	int index_count = 0; 

	vertex_count = (int)m_poly.vertices.size(); // 総頂点数
	index_count = (int)m_poly.faces.size();		// 総ポリゴン数

	cout << "vertex = " << vertex_count << " index_count = " << index_count << endl;


	//3Dモデルから粒子位置を初期化
	for(int i = 0; i < vertex_count; i++)
	{
		for(int j = 0; j < 3; j++)
		{
			pos[i*4+j] = m_poly.vertices[i][j] * 0.2f;
		}
	}

	//tetgenで頂点追加＋四面体分割
	//MakeTetrahedraFromObj();
	//MakeTetrahedraFromCube();

	////node,eleファイルを読み込んで頂点情報と四面体情報を取得
	//Load_ELE_File(ELE_FILE);
	//Load_NODE_File(NODE_FILE, p);

	////CGALを用いて3Dオブジェクト内部に点を追加＋四面体作成
	////test.test();
}

void IceTetrahedra::LoadTest_VoxelFile(const string name, float* pos)
{	cout << __FUNCTION__ << endl;
	
	////ファイルを読み込み，頂点位置を
	//ifstream ifs( name );
	//string str;

	////ファイルの存在確認
	//if(ifs.fail()) 
	//{
	//	cerr << "File do not exist.\n";
	//	exit(0);
	//}

	////変数の用意，初期化
	//int ia=0, ib=0, ic=0, id=0;
	//double da = 0.0, db = 0.0, dc = 0.0, dd = 0.0;
	//bool line_1 = false;

	////文字列の読み込み
	//while( getline(ifs, str) )
	//{
	//	ia=0; ib=0; ic=0; id=0;
	//	//無理やりだけどとりあえず
	//	if( !line_1 )
	//	{
	//		line_1 = true;
	//		sscanf(str.data(), "%d %d %d %d", &ia, &ib, &ic, &id);
	//		
	//		//cout << "ia = " << ia << "\t" << "ib = " << ib << "\t" << "ic = " << ic << endl;
	//	}
	//	else
	//	{
	//		if( str[0] == '#' )
	//		{
	//			continue;
	//		}

	//		sscanf(str.data(), "%lf %lf %lf %lf", &da, &db, &dc, &dd);
	//		//cout << "a = " << a << "\t";
	//		//cout << "b = " << b << "\t";
	//		//cout << "c = " << c << "\t";
	//		//cout << "d = " << d << "\t";
	//		//cout << "e = " << e << endl;

	//		int pIndx = (int(da)/*-1*/)*4;	//ファイルによっては-1することもある
	//		double radius = 0.08;
	//		p[pIndx+0] = db * radius;
	//		p[pIndx+1] = dc * radius;
	//		p[pIndx+2] = dd * radius;
	//	}
	//}
}