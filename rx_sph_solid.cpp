/*!
  @file rx_sph_solid.cpp
	
  @brief SPH�p�ő̒�`
 
  @author Makoto Fujisawa
  @date 2008-12
*/


//-----------------------------------------------------------------------------
// �C���N���[�h�t�@�C��
//-----------------------------------------------------------------------------
#include "rx_sph_solid.h"

#include "rx_particle_on_surf.h"
#include "rx_mc.h"


//-----------------------------------------------------------------------------
// �萔�E�ϐ�
//-----------------------------------------------------------------------------
const Vec3 RX_AABB_NORMALS[6] = { Vec3( 1.0,  0.0,  0.0), 
								  Vec3(-1.0,  0.0,  0.0), 
								  Vec3( 0.0,  1.0,  0.0), 
								  Vec3( 0.0, -1.0,  0.0), 
								  Vec3( 0.0,  0.0,  1.0), 
								  Vec3( 0.0,  0.0, -1.0) };


//-----------------------------------------------------------------------------
// �����v�Z�֐�
//-----------------------------------------------------------------------------
/*!
 * �����̂Ɠ_�̋���
 * @param[in] spos �����̂̒��S�����_�Ƃ������΍��W�l
 * @param[in] r    ���a(���̏ꍇ)
 * @param[in] sgn  �����̂̓��ŋ�������:1,�O�Ő�:-1
 * @param[in] vMin �����̂̍ŏ����W�l(���΍��W)
 * @param[in] vMax �����̂̍ő���W�l(���΍��W)
 * @param[out] d   �����t�����l
 * @param[out] n   �ŋߖT�_�̖@������
 */
bool AABB_point_dist(const Vec3 &spos, const double &r, const int &sgn, 
					 const Vec3 &vMin, const Vec3 &vMax, 
					 double &d, Vec3 &n)
{
	bitset<6> bout;
	bout.reset();
	double d0[6];
	int idx0 = -1;

	// �e�����Ƃɍŏ��ƍő勫�E�O�ɂȂ��Ă��Ȃ������ׂ�
	for(int i = 0; i < 3; ++i){
		int idx = 2*i;
		if((d0[idx] = (spos[i]-r)-vMin[i]) < 0.0){
			bout[idx] = true;
			idx0 = idx;
		}
		idx = 2*i+1;
		if((d0[idx] = vMax[i]-(spos[i]+r)) < 0.0){
			bout[idx] = true;
			idx0 = idx;
		}
	}

	// �����̓�(�S���ŋ��E��)
	if(bout.none()){
		double min_d = 1e10;
		int idx1 = -1;
		for(int i = 0; i < 6; ++i){
			if(d0[i] <= min_d){
				min_d = d0[i];
				idx1 = i;
			}
		}

		d = sgn*min_d;
		n = (idx1 != -1) ? sgn*RX_AABB_NORMALS[idx1] : Vec3(0.0);
		return true;
	}


	Vec3 x(0.0);
	for(int i = 0; i < 3; ++i){
		if(bout[2*i]){
			x[i] = d0[2*i];
		}
		else if(bout[2*i+1]){
			x[i] = -d0[2*i+1];
		}
	}

	// sgn = 1:���C-1:�I�u�W�F�N�g
	int c = (int)bout.count();
	if(c == 1){
		// ���ʋߖT
		d = sgn*d0[idx0];
		n = sgn*RX_AABB_NORMALS[idx0];
	}
	else{
		// �G�b�W/�R�[�i�[�ߖT
		d = -sgn*norm(x);
		n = sgn*(-Unit(x));
	}

	return false;
}


/*!
 * �A�֐��l�̌v�Z
 * @param[in] x,y,z �v�Z�ʒu
 * @return ���z�ƉA�֐��l���i�[����4�����x�N�g��(0�`2:���z,3:�l)
 */
RXREAL rxSolid::GetImplicit_s(void* ptr, double x, double y, double z)
{
	return ((rxSolid*)ptr)->GetImplicit(Vec3(x, y, z));
}
RXREAL rxSolid::GetImplicit(Vec3 pos)
{
	rxCollisionInfo col;
	GetDistance(pos, col);

	return col.Penetration()+m_fOffset;
}

/*!
 * �A�֐��l�Ƃ��̌��z�̌v�Z
 * @param[in] x,y,z �v�Z�ʒu
 * @return ���z�ƉA�֐��l���i�[����4�����x�N�g��(0�`2:���z,3:�l)
 */
Vec4 rxSolid::GetImplicitG_s(void* ptr, double x, double y, double z)
{
	return ((rxSolid*)ptr)->GetImplicitG(Vec3(x, y, z));
}
Vec4 rxSolid::GetImplicitG(Vec3 pos)
{
	rxCollisionInfo col;
	GetDistance(pos, col);

	return Vec4(col.Normal(), col.Penetration()+m_fOffset);
}

/*!
 * �\�ʃp�[�e�B�N������
 */
int rxSolid::GenerateParticlesOnSurf(RXREAL rad, RXREAL **ppos)
{
	int num_particles = 0;

	vector<Vec3> vrts, nrms;
	vector<rxFace> faces;

	Vec3 minp = GetMin()-Vec3(6.0*rad);
	Vec3 maxp = GetMax()+Vec3(6.0*rad);
	Vec3 dim = maxp-minp;

	m_fOffset = rad;

	// ���b�V�������i�q�̌���
	double h = 2.0*rad;
	int n[3];
	for(int i = 0; i < 3; ++i){
		n[i] = (int)(dim[i]/h)+1;
	}

	// �A�֐��l���i�[�����z��̐���
	RxScalarField sf;
	for(int i = 0; i < 3; ++i){
		sf.iNum[i] = n[i];
		sf.fWidth[i] = h;
		sf.fMin[i] = minp[i];
	}
	RXREAL *field = 0;
	GenerateValueArray(&field, rxSolid::GetImplicit_s, this, sf);
	
	// Marching Cubes�œ��l�ʃ��b�V������
	rxMCMeshCPU *mc = new rxMCMeshCPU;
	mc->CreateMeshV(field, minp, h, n, 0.0, vrts, nrms, faces);
	delete mc;
	delete [] field;

	// ���_��
	num_particles = (int)vrts.size();

	if(num_particles){
		cout << num_particles << " particles are generated for the boundary." << endl;

		rxParticleOnSurf sp;

		// �p�[�e�B�N�������z�u
		sp.Initialize(vrts, rad, minp, maxp, rxSolid::GetImplicitG_s, this);

		int iter = 50;
		RXREAL eps = 1.0e-4;

		// �p�[�e�B�N���ʒu�C��
		sp.Update(0.01, iter, eps);

		cout << iter << " iterations and the error value is " << eps << endl;

		num_particles = sp.GetNumParticles();

		if((*ppos)) delete [] (*ppos);
		*ppos = new RXREAL[DIM*num_particles];
		memset(*ppos, 0, DIM*num_particles*sizeof(RXREAL));

		RXREAL *spos = sp.GetPositionArray();
		for(int i = 0; i < DIM*num_particles; ++i){
			(*ppos)[i] = spos[i];
		}
	}

	return num_particles;
}



//-----------------------------------------------------------------------------
// MARK:rxSolidBox�N���X�̎���
//-----------------------------------------------------------------------------
/*!
 * �����l�v�Z
 * @param[in] pos �O���[�o�����W�ł̈ʒu
 * @param[out] d �ő̋��E�ߖT�_�܂ł̕����t����
 * @param[out] n �@������
 * @param[out] v �ߖT�_�ł̑��x
 * @return 
 */
bool rxSolidBox::GetDistance(const Vec3 &pos, rxCollisionInfo &col)
{
	return GetDistanceR(pos, 0.0, col);
}

/*!
 * �����l�v�Z(����)
 * @param[in] pos �O���[�o�����W�ł̈ʒu
 * @param[out] d �ő̋��E�ߖT�_�܂ł̕����t����
 * @param[out] n �@������
 * @param[out] v �ߖT�_�ł̑��x
 * @return 
 */
bool rxSolidBox::GetDistanceR(const Vec3 &pos, const double &r, rxCollisionInfo &col)
{
	Vec3 spos;

	spos = CalLocalCoord(pos);

	col.Velocity() = GetVelocityAtGrobal(pos);//Vec3(0.0);
	int sgn = -m_iSgn;

	double d;
	Vec3 n;
	AABB_point_dist(spos, r, sgn, m_vMin, m_vMax, d, n);

	col.Penetration() = d;
	col.Normal() = n;
	col.Contact() = pos+n*fabs(d);

	return (col.Penetration() <= 0.0);
}

bool rxSolidBox::GetCurvature(const Vec3 &pos, double &k)
{
	//return 0.0;
	return CalCurvature(pos, k, boost::bind(&rxSolidBox::GetDistance, this, _1, _2));
}

void rxSolidBox::SetGLMatrix(void)
{
	glTranslatef(m_vMassCenter[0], m_vMassCenter[1], m_vMassCenter[2]);
	glMultMatrixd(m_matRot.GetValue());
}

void rxSolidBox::Draw(const bool wire)
{
	glPushMatrix();

	SetGLMatrix();

	Vec3 sl = 0.5*(m_vMax-m_vMin);
	sl = RXFunc::Fabs(sl);
	glScalef(2.0*sl[0], 2.0*sl[1], 2.0*sl[2]);

	if(wire){
		glutWireCube(1.0);
	}
	else{
		glutSolidCube(1.0);
	}

	glPopMatrix();
}



//-----------------------------------------------------------------------------
// MARK:rxSolidOpenBox�N���X�̎���
//-----------------------------------------------------------------------------
/*!
 * �����l�v�Z
 * @param[in] pos �O���[�o�����W�ł̈ʒu
 * @param[out] d �ő̋��E�ߖT�_�܂ł̕����t����
 * @param[out] n �@������
 * @param[out] v �ߖT�_�ł̑��x
 * @return 
 */
bool rxSolidOpenBox::GetDistance(const Vec3 &pos, rxCollisionInfo &col)
{
	return GetDistanceR(pos, 0.0, col);
}

/*!
 * �����l�v�Z(����)
 * @param[in] pos �O���[�o�����W�ł̈ʒu
 * @param[out] d �ő̋��E�ߖT�_�܂ł̕����t����
 * @param[out] n �@������
 * @param[out] v �ߖT�_�ł̑��x
 * @return 
 */
bool rxSolidOpenBox::GetDistanceR(const Vec3 &pos, const double &r, rxCollisionInfo &col)
{
	Vec3 spos;

	spos = CalLocalCoord(pos);

	col.Velocity() = GetVelocityAtGrobal(pos);//Vec3(0.0);
	int sgn = -m_iSgn;
	int t = 2;
	double d = RX_FEQ_INF;
	Vec3 n;

	double td;
	Vec3 tn;

	// ��
	Vec3 m0, m1;
	m0 = -m_vSLenOut;
	m1 =  m_vSLenOut;
	m1[t] = m_vSLenIn[t];
	AABB_point_dist(spos, r, sgn, m0, m1, td, tn);
	if(td < d){
		d = td;
		n = tn;
	}

	// ���� -x
	m0 = Vec3(-m_vSLenOut[0], -m_vSLenOut[1], -m_vSLenIn[2]);
	m1 = Vec3(-m_vSLenIn[0],   m_vSLenOut[1],  m_vSLenOut[2]);
	AABB_point_dist(spos, r, sgn, m0, m1, td, tn);
	if(td < d){
		d = td;
		n = tn;
	}

	// ���� +x
	m0 = Vec3( m_vSLenIn[0],  -m_vSLenOut[1], -m_vSLenIn[2]);
	m1 = Vec3( m_vSLenOut[0],  m_vSLenOut[1],  m_vSLenOut[2]);
	AABB_point_dist(spos, r, sgn, m0, m1, td, tn);
	if(td < d){
		d = td;
		n = tn;
	}

	// ���� -y
	m0 = Vec3(-m_vSLenIn[0], -m_vSLenOut[1], -m_vSLenIn[2]);
	m1 = Vec3( m_vSLenIn[0], -m_vSLenIn[1],   m_vSLenOut[2]);
	AABB_point_dist(spos, r, sgn, m0, m1, td, tn);
	if(td < d){
		d = td;
		n = tn;
	}

	// ���� +y
	m0 = Vec3(-m_vSLenIn[0],  m_vSLenIn[1],  -m_vSLenIn[2]);
	m1 = Vec3( m_vSLenIn[0],  m_vSLenOut[1],  m_vSLenOut[2]);
	AABB_point_dist(spos, r, sgn, m0, m1, td, tn);
	if(td < d){
		d = td;
		n = tn;
	}

	col.Penetration() = d;
	col.Normal() = n;
	col.Contact() = pos+n*fabs(d);

	return (col.Penetration() <= 0.0);
}

bool rxSolidOpenBox::GetCurvature(const Vec3 &pos, double &k)
{
	return CalCurvature(pos, k, boost::bind(&rxSolidOpenBox::GetDistance, this, _1, _2));
}

void rxSolidOpenBox::SetGLMatrix(void)
{
	glTranslatef(m_vMassCenter[0], m_vMassCenter[1], m_vMassCenter[2]);
	glMultMatrixd(m_matRot.GetValue());
}

inline void SetVerticesCube(const Vec3 &cn, const Vec3 &sl, Vec3 v[8])
{
	v[0] = cn+Vec3(-sl[0], -sl[1], -sl[2]);
	v[1] = cn+Vec3(-sl[0],  sl[1], -sl[2]);
	v[2] = cn+Vec3(-sl[0],  sl[1],  sl[2]);
	v[3] = cn+Vec3(-sl[0], -sl[1],  sl[2]);

	v[4] = cn+Vec3( sl[0], -sl[1], -sl[2]);
	v[5] = cn+Vec3( sl[0],  sl[1], -sl[2]);
	v[6] = cn+Vec3( sl[0],  sl[1],  sl[2]);
	v[7] = cn+Vec3( sl[0], -sl[1],  sl[2]);
}

inline void CreateBoxPolygon(const Vec3 &sl0, const Vec3 &sl1, const int &d, 
							 vector<Vec3> &vrts, vector< vector<int> > &idxs)
{
	if(d < 0 || d > 2) return;

	double h = sl1[d]-sl0[d];
	Vec3 cn(0.0);
	

	vrts.resize(16);

	// �O���̒��_
	SetVerticesCube(cn, sl1, &vrts[0]);

	// �����̒��_
	cn[d] += h;
	SetVerticesCube(cn, sl0, &vrts[8]);


	int idxs0[5][4] = { {0, 3, 2, 1}, 
						{1, 2, 6, 5}, 
						{5, 6, 7, 4}, 
						{4, 7, 3, 0}, 
						{0, 1, 5, 4} };
	
	int idxs1[4][4] = { {2, 3, 11, 10}, 
						{3, 7, 15, 11}, 
						{7, 6, 14, 15}, 
						{6, 2, 10, 14} };
	
	// �O�p�`�쐬
	idxs.resize(28);
	for(int i = 0; i < 28; ++i) idxs[i].resize(3);

	int c = 0;

	// �O���̔�
	for(int i = 0; i < 5; ++i){
		for(int j = 0; j < 3; ++j){
			idxs[c][j] = idxs0[i][j];
		}
		c++;
		for(int j = 0; j < 3; ++j){
			idxs[c][j] = idxs0[i][((j+2 > 3) ? 0 : j+2)];
		}
		c++;
	}

	// �����̔�
	for(int i = 0; i < 5; ++i){
		for(int j = 0; j < 3; ++j){
			idxs[c][j] = idxs0[i][2-j]+8;
		}
		c++;
		for(int j = 0; j < 3; ++j){
			idxs[c][j] = idxs0[i][(((2-j)+2 > 3) ? 0 : (2-j)+2)]+8;
		}
		c++;
	}

	// �㕔
	for(int i = 0; i < 4; ++i){
		for(int j = 0; j < 3; ++j){
			idxs[c][j] = idxs1[i][j];
		}
		c++;
		for(int j = 0; j < 3; ++j){
			idxs[c][j] = idxs1[i][((j+2 > 3) ? 0 : j+2)];
		}
		c++;
	}

}

void rxSolidOpenBox::Draw(const bool wire)
{
	glPushMatrix();

	SetGLMatrix();

	Vec3 len0 = 2.0*m_vSLenIn;
	Vec3 len1 = 2.0*m_vSLenOut;
	double d = 0.5*(len1[2]-len0[2]);

	glTranslatef(0.0, 0.0, 0.5*len1[2]);
	if(wire){
		glPushMatrix();
		glTranslatef(0.0, 0.0, d);
		glScalef(len0[0], len0[1], len0[2]);
		glutWireCube(1.0);
		glPopMatrix();

		glPushMatrix();
		glScalef(len1[0], len1[1], len1[2]);
		glutWireCube(1.0);
		glPopMatrix();
	}
	else{
		vector<Vec3> vrts;
		vector< vector<int> > idxs;
		CreateBoxPolygon(m_vSLenIn, m_vSLenOut, 2, vrts, idxs);

		// �C���f�b�N�X��1�n�܂��
		int n = (int)idxs.size();
		for(int i = 0; i < n; ++i){
			for(int j = 0; j < 3; ++j){
				idxs[i][j]++;
			}
		}

		// ���C���[�t���[���`��
		glDisable(GL_LIGHTING);
		glPushMatrix();
		glTranslatef(0.0, 0.0, d);
		glScalef(len0[0], len0[1], len0[2]);
		glutWireCube(1.0);
		glPopMatrix();

		glPushMatrix();
		glScalef(len1[0], len1[1], len1[2]);
		glutWireCube(1.0);
		glPopMatrix();

		// �ʕ`��
		glEnable(GL_LIGHTING);
		rxMaterial mat;
		mat.SetGL();
		glColor3f(0.0, 0.0, 1.0);
		for(int i = 0; i < (int)idxs.size(); ++i){
			glBegin(GL_POLYGON);
			for(int j = 0; j < (int)idxs[i].size(); ++j){
				glVertex3dv(vrts[idxs[i][j]-1].data);
			}
			glEnd();
		}
	}


	glPopMatrix();
}



//-----------------------------------------------------------------------------
// MARK:rxSolidSphere�N���X�̎���
//-----------------------------------------------------------------------------
bool rxSolidSphere::GetDistance(const Vec3 &pos, rxCollisionInfo &col)
{
	return GetDistanceR(pos, 0.0, col);
}

bool rxSolidSphere::GetDistanceR(const Vec3 &pos, const double &r, rxCollisionInfo &col)
{
	Vec3 rpos = pos-m_vMassCenter;
	double d = m_iSgn*(norm(rpos)-m_fRadius);
	if(d < r){
		Vec3 n = Unit(rpos);

		col.Penetration() = d-r;
		col.Contact() = m_vMassCenter+n*(m_fRadius+m_iSgn*r);
		col.Normal() = m_iSgn*n;

		col.Velocity() = GetVelocityAtGrobal(pos);
	}
	else{
		return false;
	}

	return (col.Penetration() <= 0.0);
}

void rxSolidSphere::SetGLMatrix(void)
{
	glTranslatef(m_vMassCenter[0], m_vMassCenter[1], m_vMassCenter[2]);
	glMultMatrixd(m_matRot.GetValue());
}

/*!
 * ���_���S�̉~�̃��C���[�t���[���`��
 * @param rad �~�̔��a
 * @param n ������
 */
static void DrawWireCircle(const double &rad, const int &n)
{
	double t = 0.0;
	double dt = 2.0*RX_PI/(double)n;
 
	glBegin(GL_LINE_LOOP);
	do{
		glVertex3f(rad*cos(t), rad*sin(t), 0.0);
		t += dt;
	}while(t < 2.0*RX_PI);
	glEnd();
}
 
/*!
 * ���_���S�̉~�̃��C���[�t���[���`��(XZ����)
 * @param rad �~�̔��a
 * @param n ������
 */
static void DrawWireCircleXZ(const double &rad, const int &n)
{
	double t = 0.0;
	double dt = 2.0*RX_PI/(double)n;
 
	glBegin(GL_LINE_LOOP);
	do{
		glVertex3f(rad*cos(t), 0.0, rad*sin(t));
		t += dt;
	}while(t < 2.0*RX_PI);
	glEnd();
}
 
/*!
 * ���̃��C���[�t���[���`��
 * @param cen ���̒��S
 * @param rad ���̔��a
 * @param col �`��F
 */
void DrawWireSphere(const Vec3 &cen, const float &rad, const Vec3 &col)
{
	glDisable(GL_LIGHTING);
	glPushMatrix();
	glTranslatef(cen[0], cen[1], cen[2]);
	glRotatef(90, 1.0, 0.0, 0.0);
	glColor3f(col[0], col[1], col[2]);
 
	// �ܓx(x-y���ʂɕ��s)
	float z, dz;
	dz = 2.0*rad/8.0f;
	z = -(rad-dz);
	do{
		glPushMatrix();
		glTranslatef(0.0, 0.0, z);
		DrawWireCircle(sqrt(rad*rad-z*z), 32);
		glPopMatrix();
		z += dz;
	}while(z < rad);
 
	// �o�x(z���܂��ɉ�])
	float t, dt;
	t = 0.0f;
	dt = 180.0/8.0;
	do{
		glPushMatrix();
		glRotatef(t,  0.0, 0.0, 1.0);
		DrawWireCircleXZ(rad, 32);
		glPopMatrix();
 
		t += dt;
	}while(t < 180);
 
	//glutWireSphere(rad, 10, 5);
	glPopMatrix();
}

void rxSolidSphere::Draw(const bool wire)
{
	glPushMatrix();
	SetGLMatrix();

	if(wire){
		DrawWireSphere(Vec3(0.0), m_fRadius, Vec3(0.0, 1.0, 0.0));
	}
	else{
		glPushMatrix();
		glRotated(90, 1.0, 0.0, 0.0);
		glutSolidSphere(m_fRadius, 20, 10);
		glPopMatrix();
	}

	glPopMatrix();
}



