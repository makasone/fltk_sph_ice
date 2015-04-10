#include "IceTetrahedra.h"
#include "rx_model.h"

IceTetrahedra &IceTetrahedra::GetInstance() {
    static IceTetrahedra instance;
    return instance;
}

/*!
 * �l�ʑ̏��̏�����
 */
void IceTetrahedra::InitTetra(float* pos, int vertexNum)
{	cout << __FUNCTION__ << endl;

	//Load_ELE_File(ELE_FILE);					//���܂������Ȃ��@ele�t�@�C����ǂݍ��݃��X�g�쐬
	//Load_obj_File(OBJ_NAME, pos);				//obj�t�@�C����ǂݍ��݃��X�g�쐬
	Load_NODE_File(NODE_FILE, pos);			//�o�j�[���f����p����ꍇ�Ɏg��

	//MakeTetrahedraFromCube(pos, vertexNum);
}

/*!
 * tetgen��p���ē_�Q�̈ʒu���C���b�V����񂩂珉���̎l�ʑ̍\�����쐬.
 * �_���v�����t�@�C���� src/fltk_sph_turb/bin�@�ɍ����D
 */
void IceTetrahedra::MakeTetrahedraFromCube(float* pos, int vrtxNum)
{	cout << __FUNCTION__ << endl;

	tetgenio in, out;	// ���̓��b�V���Əo�͎l�ʑ̃��b�V��

	// �|���S�����_�C���f�b�N�X�̃X�^�[�g(0�X�^�[�g��1�X�^�[�g)
	in.firstnumber = 1;

	// ���b�V�����_�̐ݒ�
	in.numberofpoints = vrtxNum;
	in.pointlist = new double[in.numberofpoints*3];
	for(int i = 0; i < in.numberofpoints; ++i){
		for(int j = 0; j < 3; ++j){
			//���q�ʒu���̓o�^
			in.pointlist[3*i+j] = (double)pos[4*i+j];
		}
	}

	// �|���S���̐ݒ�
	int n = pow( in.numberofpoints, 1.0/3.0 ) + 0.5;	//�����̂̂P�ӂ̒��_��
	if( n == 1 ){	cout << "error::n == 1" << endl;	}

	//�e�|���S���̒��_�ԍ����v�Z
	//�i�q��ɔz�u���ꂽ���q�ō����ʂ�����Ă���D
	//�i�q�̓����ɑ��݂���ʂ��l������K�v������̂ł�₱�����Ȃ��Ă���D
	vector<vector<int>> poligonList;
	int poligonNum = 0;

	//�㉺��
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

	//���E��
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

	//�O���
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

	in.numberoffacets = (n-1) * (n-1) * n * 3;				//�ʂ̑��� (��ӂ̒��_��-1)��2��*(��ӂ̒��_��)*�i���̐��j
	in.facetlist = new tetgenio::facet[in.numberoffacets];
	in.facetmarkerlist = new int[in.numberoffacets];

	for(int i = 0; i < in.numberoffacets; ++i){
		tetgenio::facet *f = &in.facetlist[i];

		// �|���S�����X�g
		f->numberofpolygons = 1;
		f->polygonlist = new tetgenio::polygon[f->numberofpolygons];

		// hole(��)���X�g
		f->numberofholes = 0;
		f->holelist = NULL;

		// �|���S�����_�C���f�b�N�X
		tetgenio::polygon *p = &f->polygonlist[0];

		p->numberofvertices = 4;
		p->vertexlist = new int[p->numberofvertices];
		for(int j = 0; j < p->numberofvertices; ++j)
		{
			p->vertexlist[j] = poligonList[i][j];
		}

		in.facetmarkerlist[i] =  1;	//??
	}

	// ���̓��b�V�������t�@�C���Ƀ_���v
	//in.save_poly("test_poly");	// test.poly
	
	// �������́C"p":PLC�ǂݍ��݁C"q":quality mesh generation(q�̌��quality bound�𐔒l�Ŏw��)�C
	// "a":�ő�̐ϐ���(a�̌�ɑ̐ς𐔒l�Ŏw��)
	tetgenbehavior b;
//	if(!b.parse_commandline("-pq1.414a0.1n")){
	if(!b.parse_commandline("-pn")){
		terminatetetgen(0);
	}

	tetrahedralize(&b, &in, &out);				// �l�ʑ̃��b�V������

	//std::cout << " the number of node   : "    << out.numberofpoints << endl;
	//std::cout << " the number of element   : " << out.numberoftetrahedra << endl;
	//std::cout << " the number of face  : "     << out.numberoftrifaces << endl;

	// �o�̓��b�V�������t�@�C���Ƀ_���v
//	out.save_elements("test_out");	// .ele

	// �o�͎l�ʑ̒��_�ԍ���z��Ɋi�[
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
//	delete[] in.pointlist;	//�K�v�H
}

/*!
 * tetgen��p���ē_�Q�̈ʒu���C���b�V����񂩂珉���̎l�ʑ̍\�����쐬.
 * �����́@������
 * �_���v�����t�@�C���� src/fltk_sph_turb/bin�@�ɍ����D
 */
void IceTetrahedra::MakeTetrahedraRectParallele(float* pos, int vrtxNum, int x, int y, int z)
{	cout << __FUNCTION__ << endl;
	
	tetgenio in, out;	// ���̓��b�V���Əo�͎l�ʑ̃��b�V��

	// �|���S�����_�C���f�b�N�X�̃X�^�[�g(0�X�^�[�g��1�X�^�[�g)
	in.firstnumber = 1;

	// ���b�V�����_�̐ݒ�
//	in.numberofpoints = m_pPS->GetNumParticles();
	in.numberofpoints = vrtxNum;
	in.pointlist = new double[in.numberofpoints*3];
	for(int i = 0; i < in.numberofpoints; ++i){
		for(int j = 0; j < 3; ++j){
			//���q�ʒu���̓o�^
			in.pointlist[3*i+j] = (double)pos[4*i+j];
		}
	}

	//�e�|���S���̒��_�ԍ����v�Z
	//�i�q��ɔz�u���ꂽ���q�ō����ʂ�����Ă���D
	//�i�q�̓����ɑ��݂���ʂ��l������K�v������̂ł�₱�����Ȃ��Ă���D
	vector<vector<int>> poligonList;
	int poligonNum = 0;

	//�O���
	//�w�������猩�āCX���ɑ�����Y���ɑ�����Z���ɑ���
	//�l�p�`���~���l�߂�ꂽ�w����������̂��C�ƍl����̂�xx=1,yy=1,zz=0
	for(int zz = 0; zz < z; zz++){
		for(int yy = 1; yy < y; yy++){
			for(int xx = 1; xx < x; xx++){
				vector<int> list;
				int Layer = zz*x*y;
				int Heigh = (yy-1)*x;
				int Line  = xx;

				//�w�{�����{�l�p
				list.push_back(Layer+Heigh+Line);
				list.push_back(Layer+Heigh+Line+1);
				list.push_back(Layer+Heigh+Line+1+x);
				list.push_back(Layer+Heigh+Line+x);
				poligonList.push_back(list);
	}}}

	//���E��
	//�y�������猩�āCY���ɑ�����Z���ɑ�����X���ɑ���
	for(int xx = 0; xx < x; xx++){
		for(int zz = 1; zz < z; zz++){
			for(int yy = 1; yy < y; yy++){
				vector<int> list;
				int Layer = xx;
				int Heigh = (zz-1)*x*y;
				int Line  = 1+(yy-1)*x;

				//�w�{�����{�l�p xLine�Ŋ�_�����߂āC�l�p�`�����
				list.push_back(Layer+Heigh+Line );
				list.push_back(Layer+Heigh+Line +x);
				list.push_back(Layer+Heigh+Line +x*y+x);
				list.push_back(Layer+Heigh+Line +x*y);
				poligonList.push_back(list);
	}}}

	////�㉺��
	for(int yy = 0; yy < y; yy++){
		for(int zz = 1; zz < z; zz++){
			for(int xx = 1; xx < x; xx++){
				vector<int> list;
				int Layer = yy*x;
				int Heigh = (zz-1)*x*y;
				int Line  = xx;

				//�w�{�����{�l�p xLine�Ŋ�_�����߂āC�l�p�`�����
				list.push_back(Layer+Heigh+Line );
				list.push_back(Layer+Heigh+Line +1);
				list.push_back(Layer+Heigh+Line +x*y+1);
				list.push_back(Layer+Heigh+Line +x*y);
				poligonList.push_back(list);
	}}}

	in.numberoffacets = (x-1) * (y-1) * z + (y-1) * (z-1) * x + (x-1) * (z-1) * y;	//�ʂ̑���
	in.facetlist = new tetgenio::facet[in.numberoffacets];
	in.facetmarkerlist = new int[in.numberoffacets];

	for(int i = 0; i < in.numberoffacets; ++i){
		tetgenio::facet *f = &in.facetlist[i];

		// �|���S�����X�g
		f->numberofpolygons = 1;
		f->polygonlist = new tetgenio::polygon[f->numberofpolygons];

		// hole(��)���X�g
		f->numberofholes = 0;
		f->holelist = NULL;

		// �|���S�����_�C���f�b�N�X
		tetgenio::polygon *p = &f->polygonlist[0];

		p->numberofvertices = 4;
		p->vertexlist = new int[p->numberofvertices];
		for(int j = 0; j < p->numberofvertices; ++j)
		{
			p->vertexlist[j] = poligonList[i][j];
		}

		in.facetmarkerlist[i] =  1;	//??
	}

	//���͏���s�t�@�C���Ƀ_���v
	//in.save_poly("Input_RectParallele");
	//in.save_nodes("Input_RectParallele");

	// �������́C"p":PLC�ǂݍ��݁C"q":quality mesh generation(q�̌��quality bound�𐔒l�Ŏw��)�C
	// "a":�ő�̐ϐ���(a�̌�ɑ̐ς𐔒l�Ŏw��)
	tetgenbehavior b;
//	if(!b.parse_commandline("-pq1.414a0.1n")){
	if(!b.parse_commandline("-pn")){
//	if(!b.parse_commandline(TETGENCOMMAND)){
		terminatetetgen(0);
	}

	tetrahedralize(&b, &in, &out);				// �l�ʑ̃��b�V������

	std::cout << " the number of node   : "    << out.numberofpoints << endl;
	std::cout << " the number of element   : " << out.numberoftetrahedra << endl;
	std::cout << " the number of face  : "     << out.numberoftrifaces << endl;

	// �o�̓��b�V�������t�@�C���Ƀ_���v
	//out.save_elements("Output_RectParallele");	// .ele

	// �o�͎l�ʑ̒��_�ԍ���z��Ɋi�[
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
 * tetgen��p���ē_�Q�̈ʒu���C���b�V����񂩂珉���̎l�ʑ̍\�����쐬.
 * �t�@�C���� src/fltk_sph_turb/bin�@�ɍ����D
 */
void IceTetrahedra::MakeTetrahedraOnlySurface(float* pos,int vrtxNum)
{	cout << __FUNCTION__ << endl;
	tetgenio in, out;	// ���̓��b�V���Əo�͎l�ʑ̃��b�V��

	// �|���S�����_�C���f�b�N�X�̃X�^�[�g(0�X�^�[�g��1�X�^�[�g)
	in.firstnumber = 1;

	// ���b�V�����_�̐ݒ�
	in.numberofpoints = vrtxNum;
	in.pointlist = new double[in.numberofpoints*3];
	for(int i = 0; i < in.numberofpoints; ++i){
		for(int j = 0; j < 3; ++j){
			in.pointlist[3*i+j] = (double)pos[4*i+j];
		}
	}

	// �|���S���̐ݒ�
	int n = pow( in.numberofpoints, 1.0/3.0 ) + 0.5;	//�����̂̂P�ӂ̒��_��
	if( n == 1 ){	cout << "error::n == 1" << endl;	}

	//�e�|���S���̒��_�ԍ����v�Z
	vector<vector<int>> poligonList;
	int poligonNum = 0;

	//�㉺��
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

	//���E��
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

	//�O���
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

		// �|���S�����X�g
		f->numberofpolygons = 1;
		f->polygonlist = new tetgenio::polygon[f->numberofpolygons];

		// hole(��)���X�g
		f->numberofholes = 0;
		f->holelist = NULL;

		// �|���S�����_�C���f�b�N�X
		tetgenio::polygon *p = &f->polygonlist[0];

		p->numberofvertices = 4;
		p->vertexlist = new int[p->numberofvertices];
		for(int j = 0; j < p->numberofvertices; ++j)
		{
			p->vertexlist[j] = poligonList[i][j];
		}

		in.facetmarkerlist[i] =  1;	//??
	}

	// ���̓��b�V�������t�@�C���Ƀ_���v
	//in.save_poly("test_poly");	// test.poly
	
	// �������́C"p":PLC�ǂݍ��݁C"q":quality mesh generation(q�̌��quality bound�𐔒l�Ŏw��)�C
	// "a":�ő�̐ϐ���(a�̌�ɑ̐ς𐔒l�Ŏw��)
	tetgenbehavior b;
//	if(!b.parse_commandline("-pq1.414a0.1n")){
	if(!b.parse_commandline("-pn")){
		terminatetetgen(0);
	}

	tetrahedralize(&b, &in, &out);				// �l�ʑ̃��b�V������

	//std::cout << " the number of node   : "    << out.numberofpoints << endl;
	//std::cout << " the number of element   : " << out.numberoftetrahedra << endl;
	//std::cout << " the number of face  : "     << out.numberoftrifaces << endl;

	// �o�̓��b�V�������t�@�C���Ƀ_���v
//	out.save_elements("test_out");	// .ele

	// �o�͎l�ʑ̒��_�ԍ���z��Ɋi�[
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
// * tetgen��p���ē_�Q�̈ʒu���C���b�V����񂩂珉���̎l�ʑ̍\�����쐬.
// * �t�@�C���� src/fltk_sph_turb/bin�@�ɍ����D
// */
//void IceTetrahedra::MakeTetrahedraFromObj(rxPolygons poly,int vrtxNum)
//{	cout << __FUNCTION__ << endl;
//
//	tetgenio in, out;	// ���̓��b�V���Əo�͎l�ʑ̃��b�V��
//
//	// �|���S�����_�C���f�b�N�X�̃X�^�[�g(0�X�^�[�g��1�X�^�[�g)
//	in.firstnumber = 0;
//
//	// ���b�V�����_�̐ݒ�
//	in.numberofpoints = vrtxNum;
//	in.pointlist = new double[in.numberofpoints*3];
//	for(int i = 0; i < in.numberofpoints; ++i){
//		for(int j = 0; j < 3; ++j){
//			//���q�ʒu���̓o�^
//			in.pointlist[i*3+j] = poly.vertices[i][j];
//		}
//	}
//
//	// �|���S���̐ݒ�
//	//�|���S���̓ǂݍ��ݏ����g��
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
//	in.numberoffacets = poly.faces.size();				//�ʂ̑��� (��ӂ̒��_��-1)��2��*(��ӂ̒��_��)*�i���̐��j
//	in.facetlist = new tetgenio::facet[in.numberoffacets];
//	in.facetmarkerlist = new int[in.numberoffacets];
//
//	for(int i = 0; i < in.numberoffacets; ++i)
//	{
//		tetgenio::facet *f = &in.facetlist[i];
//
//		// �|���S�����X�g
//		f->numberofpolygons = 1;
//		f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
//
//		// hole(��)���X�g
//		f->numberofholes = 0;
//		f->holelist = NULL;
//
//		// �|���S�����_�C���f�b�N�X
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
//	// ���̓��b�V�������t�@�C���Ƀ_���v
//	in.save_poly("test_poly");	// �s���S�炵��
//	in.save_nodes("test_node");
//	//in.save_faces("test_face");
//
//	// �������́C"p":PLC�ǂݍ��݁C"q":quality mesh generation(q�̌��quality bound�𐔒l�Ŏw��)�C
//	// "a":�ő�̐ϐ���(a�̌�ɑ̐ς𐔒l�Ŏw��)
//	tetgenbehavior b;
////	if(!b.parse_commandline("-pq1.414a0.1n")){
////	if(!b.parse_commandline("-pn")){
//	if(!b.parse_commandline(TETGENCOMMAND)){
//		terminatetetgen(0);
//	}
//
//	tetrahedralize(&b, &in, &out);				// �l�ʑ̃��b�V������
//
//	std::cout << " the number of node   : "    << out.numberofpoints << endl;
//	std::cout << " the number of element   : " << out.numberoftetrahedra << endl;
//	std::cout << " the number of face  : "     << out.numberoftrifaces << endl;
//
//	// �o�̓��b�V�������t�@�C���Ƀ_���v
//	out.save_elements("test_out");	// .ele
//
//	// �o�͎l�ʑ̒��_�ԍ���z��Ɋi�[
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
 * tetgen�œ���ꂽ�t�@�C����ǂݍ��݁C������Ԃ��쐬�D�l�ʑ̏��D
 * �t�@�C���� src/fltk_sph_turb/bin�@����ǂݍ��ށD
 */
void IceTetrahedra::Load_ELE_File(const string name)
{	cout << __FUNCTION__ << endl;

	//�t�@�C����ǂݍ��݁C�l�ʑ̂ƂȂ�_�̑g�ݍ��킹��List�ɓ����D
	ifstream ifs( name );
	string str;

	//�t�@�C���̑��݊m�F
	if(ifs.fail()) 
	{
		cerr << "File do not exist.\n";
		exit(0);
	}

	//�ϐ��̗p�ӁC������
	int a=0, b=0, c=0, d=0, e=0, f=0;
	bool line_1 = false;

	m_vviTetraList.clear();

	//������̓ǂݍ���
	while( getline(ifs, str) )
	{
		a=0; b=0; c=0; d=0; e=0;
		//������肾���ǂƂ肠����
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

	//�f�o�b�O
	cout << "m_vviTetraList.size() = " << m_vviTetraList.size() << endl;
}

/*!
 * tetgen�œ���ꂽ�t�@�C����ǂݍ��݁C������Ԃ��쐬�D�l�ʑ̏��D
 * �t�@�C���� src/fltk_sph_turb/bin�@����ǂݍ��ށD
 */
void IceTetrahedra::Load_NODE_File(const string name, float* p)
{	cout << __FUNCTION__ << endl;

	//�t�@�C����ǂݍ��݁C�l�ʑ̂ƂȂ�_�̑g�ݍ��킹��List�ɓ����D
	ifstream ifs( name );
	string str;

	//�t�@�C���̑��݊m�F
	if(ifs.fail()) 
	{
		cerr << "File do not exist.\n";
		exit(0);
	}

	//�ϐ��̗p�ӁC������
	int ia=0, ib=0, ic=0, id=0;
	double da = 0.0, db = 0.0, dc = 0.0, dd = 0.0;
	bool line_1 = false;

	//������̓ǂݍ���
	while( getline(ifs, str) )
	{
		ia=0; ib=0; ic=0; id=0;
		//������肾���ǂƂ肠����
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
			int pIndx = (int(da)-1)*4;	//�t�@�C���ɂ���Ă�-1���邱�Ƃ�����
			
			double radius = 0.75f;
			p[pIndx+0] = db * radius;
			p[pIndx+1] = dc * radius - 0.25;	//TODO:�����ŏ����ʒu�𒲐����Ȃ��Ƃ����Ȃ��c
			p[pIndx+2] = dd * radius;
		}
	}
}

/*!
 * ������Ԃ�tetgen�p�ɕϊ������t�@�C�����쐬�D���_���C�ʏ��D
 * �t�@�C���� src/fltk_sph_turb/bin�@�ɍ쐬�����D
 */
void IceTetrahedra::Save_POLY_File(const string name, float* pos, int vrtxNum)
{	cout << __FUNCTION__ << endl;

	//�t�@�C�����쐬���C�t�H�[�}�b�g�ɂ��������Ĉʒu���Ɩʏ����������ށ@#include <fstream>�̈ʒu�ɒ���
	ofstream ofs( name+".poly" );

	//Poly�t�@�C���p�e�X�g
//�P�@���_�̈ʒu���
	//���_�S�̂̏��@���_���C�������i�R�ŌŒ�j�Cattribute�Cboundarymark�D
	ofs << "# Part 1 - node list" << endl;
	ofs << "       " << vrtxNum << " ";
	ofs << "3 ";
	ofs << "0 ";
	ofs << "0 ";
	ofs << endl;

	//���_���ꂼ��̏��
	for(int i = 0; i < vrtxNum; i++ )
	{
		ofs << "       " << i+1 << " " << pos[4*i+0] << " " << pos[4*i+1] << " " << pos[4*i+2] << endl;
	}
//�Q�@�ʂ̐ڑ����
	//���_�ō쐬�����ʂ̏��
	ofs << "# Part 2 - facet list" << endl;

	int n = pow( vrtxNum, 1.0/3.0 ) + 0.5;	//�����̂̂P�ӂ̒��_��
	if( n == 1 )
	{
		cout << "error::n == 1" << endl;
	}

	////�\�ʗ��q�̂݃o�[�W����
//	ofs << "\t" << (n-1) * (n-1) * 6 << "\t";		//������Ԃ͕K�������́C�Ƃ����ꍇ�̖ʂ̐�
	ofs << "\t" << (n-1) * (n-1) * n * 3 << "\t";		//������Ԃ͕K�������́C�Ƃ����ꍇ�̖ʂ̐�
	ofs << "\t" << 1 << "\t";
	ofs << endl;

	//�㉺��
	for( int k = 0; k < n; k++ )
	{
		for( int i = 1; i < n; i++ )
		{
			for( int j = 1; j < n; j++ )
			{
				ofs << "\t" << 1 << "\t" << 0 << "\t" << 1 << endl;		//�p�����[�^
				ofs << "\t" << 4					<< "\t"				//�ʂ̐�
							<< k*n*n+(i-1)*n+j		<< "\t" 
							<< k*n*n+(i-1)*n+j+1	<< "\t"
							<< k*n*n+(i-1)*n+j+1+n	<< "\t"
							<< k*n*n+(i-1)*n+j+n
				<< endl;
			}
		}
	}

	//���E��
	for( int k = 0; k < n; k++ )
	{
		for( int i = 1; i < n; i++ )
		{
			for( int j = 1; j < n; j++ )
			{
				ofs << "\t" << 1 << "\t" << 0 << "\t" << 1 << endl;		//�p�����[�^
				ofs << "\t" << 4						<< "\t"			//�ʂ̐�
							<< k+(i-1)*n*n+(j-1)*n+1	<< "\t" 
							<< k+(i-1)*n*n+(j-1)*n+1+n	<< "\t"
							<< k+i*n*n+(j-1)*n+1+n		<< "\t"
							<< k+i*n*n+(j-1)*n+1
				<< endl;
			}
		}
	}

	//�O���
	for( int k = 0; k < n; k++ )
	{
		for( int i = 1; i < n; i++ )
		{
			for( int j = 1; j < n; j++ )
			{
				ofs << "\t" << 1 << "\t" << 0 << "\t" << 1 << endl;		//�p�����[�^
				ofs << "\t" << 4						<< "\t"			//�ʂ̐�
							<< k*n+(i-1)*n*n+j			<< "\t" 
							<< k*n+(i-1)*n*n+j+1		<< "\t"
							<< k*n+(i-1)*n*n+j+1+n*n	<< "\t"
							<< k*n+(i-1)*n*n+j+n*n
				<< endl;
			}
		}
	}

//�R�@
	ofs << "# Part 3 - hole list" << endl;
	ofs << "0";
	ofs << endl;

//�S
	ofs << "# Part 4 - region list" << endl;
	ofs << "0";
	ofs << endl;
}

/*!
 * ������Ԃ�tetgen�p�ɕϊ������t�@�C�����쐬�D���_���̂݁D
 * �t�@�C���� src/fltk_sph_turb/bin�@�ɍ쐬�����D
 */
void IceTetrahedra::Save_NODE_File(const string name, float* pos, int vrtxNum)
{	cout << __FUNCTION__ << endl;

	//�t�@�C�����쐬���C�t�H�[�}�b�g�ɂ��������Ĉʒu�����������ށ@#include <fstream>�̈ʒu�ɒ���
	std::ofstream ofs( name+".node" );

	//���_�S�̂̏��@���_���C�������i�R�ŌŒ�j�Cattribute�Cboundarymark�D
	ofs << "       " << vrtxNum << " ";
	ofs << "3 ";
	ofs << "0 ";
	ofs << "0 ";
	ofs << endl;

	//���_���ꂼ��̏��
	for(int i = 0; i < vrtxNum; i++ )
	{
		ofs << "       " << i+1 << " " << pos[4*i+0] << " " << pos[4*i+1] << " " << pos[4*i+2] << endl;
	}
}

//obj�t�@�C����ǂݍ���ŗ��q�ʒu��������
void IceTetrahedra::Load_obj_File(const string objName, float* pos)
{	cout << __FUNCTION__ << endl;

	//obj�t�@�C��
	rxPolygons m_poly;

	//OpenFile(name);						//���f���ǂݍ��݁@�p�X��Release/bin�̂��
	RxModel::Read(objName, m_poly);

	int vertex_count = 0;
	int index_count = 0; 

	vertex_count = (int)m_poly.vertices.size(); // �����_��
	index_count = (int)m_poly.faces.size();		// ���|���S����

	cout << "vertex = " << vertex_count << " index_count = " << index_count << endl;


	//3D���f�����痱�q�ʒu��������
	for(int i = 0; i < vertex_count; i++)
	{
		for(int j = 0; j < 3; j++)
		{
			pos[i*4+j] = m_poly.vertices[i][j] * 0.2f;
		}
	}

	//tetgen�Œ��_�ǉ��{�l�ʑ̕���
	//MakeTetrahedraFromObj();
	//MakeTetrahedraFromCube();

	////node,ele�t�@�C����ǂݍ���Œ��_���Ǝl�ʑ̏����擾
	//Load_ELE_File(ELE_FILE);
	//Load_NODE_File(NODE_FILE, p);

	////CGAL��p����3D�I�u�W�F�N�g�����ɓ_��ǉ��{�l�ʑ̍쐬
	////test.test();
}

void IceTetrahedra::LoadTest_VoxelFile(const string name, float* pos)
{	cout << __FUNCTION__ << endl;
	
	////�t�@�C����ǂݍ��݁C���_�ʒu��
	//ifstream ifs( name );
	//string str;

	////�t�@�C���̑��݊m�F
	//if(ifs.fail()) 
	//{
	//	cerr << "File do not exist.\n";
	//	exit(0);
	//}

	////�ϐ��̗p�ӁC������
	//int ia=0, ib=0, ic=0, id=0;
	//double da = 0.0, db = 0.0, dc = 0.0, dd = 0.0;
	//bool line_1 = false;

	////������̓ǂݍ���
	//while( getline(ifs, str) )
	//{
	//	ia=0; ib=0; ic=0; id=0;
	//	//������肾���ǂƂ肠����
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

	//		int pIndx = (int(da)/*-1*/)*4;	//�t�@�C���ɂ���Ă�-1���邱�Ƃ�����
	//		double radius = 0.08;
	//		p[pIndx+0] = db * radius;
	//		p[pIndx+1] = dc * radius;
	//		p[pIndx+2] = dd * radius;
	//	}
	//}
}