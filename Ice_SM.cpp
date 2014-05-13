/*!
  @file Ice_SM.h
	
  @brief ShapeMatching�@����ɂ������ω��V�~�����[�V����
  @ref M. Muller et al., "Meshless deformations based on shape matching", SIGGRAPH2005. 
 
  @author Ryo Nakasone
  @date 2013-10
*/

#include "Ice_SM.h"

Ice_SM::Ice_SM(int obj) : rxShapeMatching(obj)
{
	m_mtrx3Apq = rxMatrix3(0.0);
}

/*!
 * �f�X�g���N�^
 */
Ice_SM::~Ice_SM()
{
}

void Ice_SM::AddVertex(const Vec3 &pos, double mass, int pIndx)
{
	rxShapeMatching::AddVertex(pos, mass);
	
	m_iParticleIndxes.push_back(pIndx);

	m_dAlphas.push_back(1.0);
	m_dBetas.push_back(0.0);

	m_iLayeres.push_back(0);

	m_iLinearDeformation.push_back(0);
	m_iVolumeConservation.push_back(0);

	//�d�S�̍X�V
	double massSum = 0.0;	// ������
	m_vec3OrgCm = Vec3(0.0);

	// �d�S���W�̌v�Z
	for(int i = 0; i < m_iNumVertices;++i){
		double m = m_vMass[i];
		//if(m_vFix[i]) m *= 300.0;	// �Œ�_�̎��ʂ�傫������
		massSum += m;
		m_vec3OrgCm += m_vOrgPos[i]*m;
	}
	m_vec3OrgCm /= massSum;

	//�ό`�s��Aqq�̍X�V
	Vec3 p, q;
	for(int i = 0; i < m_iNumVertices; ++i)
	{
		q = m_vOrgPos[i]-m_vec3OrgCm;
		double m = m_vMass[i];

		m_mtrx3AqqInv(0,0) += m*q[0]*q[0];
		m_mtrx3AqqInv(0,1) += m*q[0]*q[1];
		m_mtrx3AqqInv(0,2) += m*q[0]*q[2];
		m_mtrx3AqqInv(1,0) += m*q[1]*q[0];
		m_mtrx3AqqInv(1,1) += m*q[1]*q[1];
		m_mtrx3AqqInv(1,2) += m*q[1]*q[2];
		m_mtrx3AqqInv(2,0) += m*q[2]*q[0];
		m_mtrx3AqqInv(2,1) += m*q[2]*q[1];
		m_mtrx3AqqInv(2,2) += m*q[2]*q[2];
	}

	m_mtrx3AqqInv = m_mtrx3AqqInv.Inverse();
}

void Ice_SM::Remove(int indx)
{
	if(m_iNumVertices <= 0) return;
	m_iNumVertices--;
	m_vOrgPos.erase	( m_vOrgPos	.begin() + indx );
	m_vCurPos.erase	( m_vCurPos	.begin() + indx );
	m_vNewPos.erase	( m_vNewPos	.begin() + indx );
	m_vGoalPos.erase( m_vGoalPos.begin() + indx );
	m_vMass.erase	( m_vMass	.begin() + indx );
	m_vVel.erase	( m_vVel	.begin() + indx );
	m_vFix.erase	( m_vFix	.begin() + indx );

	m_iParticleIndxes.erase( m_iParticleIndxes.begin() + indx );

	m_iLayeres.erase( m_iLayeres.begin() + indx );

	m_dAlphas.erase	( m_dAlphas	.begin() + indx );
	m_dBetas.erase	( m_dBetas	.begin() + indx );

	m_iLinearDeformation.erase	( m_iLinearDeformation	.begin() + indx);
	m_iVolumeConservation.erase	( m_iVolumeConservation	.begin() + indx);
}

void Ice_SM::Clear()
{
	rxShapeMatching::Clear();

	m_iParticleIndxes.clear();
	
	m_iLayeres.clear();

	m_dAlphas	.clear();
	m_dBetas	.clear();
	
	m_iLinearDeformation	.clear();
	m_iVolumeConservation	.clear();
}

bool Ice_SM::CheckIndx(int pIndx)
{
	vector<int>::iterator check = find( m_iParticleIndxes.begin(), m_iParticleIndxes.end(), pIndx);

	if( check == m_iParticleIndxes.end() )	return false;

	return true;
}

int	Ice_SM::SearchIndx(int pIndx)
{
	vector<int>::iterator check = find( m_iParticleIndxes.begin(), m_iParticleIndxes.end(), pIndx);

	if( check == m_iParticleIndxes.end() )	return -1;

	return m_iParticleIndxes.size() - (m_iParticleIndxes.end() - check);
}

/*!
 * Shape Matching�@
 *  - �ڕW�ʒu���v�Z���āCm_vNewPos�����̈ʒu�ɋ߂Â���
 *  - �e���q�Ƀ��ƃ������������o�[�W����
 * @param[in] dt �^�C���X�e�b�v��
 */
void Ice_SM::ShapeMatching(double dt)
{
	if(m_iNumVertices <= 1) return;

	//clock_t oldTime, newTime;
	//oldTime = clock();
	//cout << "�v���J�n" << endl;

	Vec3 cm(0.0), cm_org(0.0);	// �d�S
	Vec3 p(0.0), q(0.0);
	rxMatrix3 Apq(0.0), Aqq(0.0);
	double mass = 0.0;	// ������

	// �d�S���W�̌v�Z
	for(int i = 0; i < m_iNumVertices;++i){
		double m = m_vMass[i];
		if(m_vFix[i]){	/*m *= 300.0;*/ m *= 1.0f;	}	// �Œ�_�̎��ʂ�傫������
		mass += m;
		//cm += m_vNewPos[i]*m;
	}
	//cm /= mass;

	cm = m_vec3NowCm / mass;
	cm_org = m_vec3OrgCm;
	
	Apq = m_mtrx3Apq;
	Aqq = m_mtrx3AqqInv;

	//Apq�̍s�񎮂����߁C���]���邩�𔻒�
	//�s����ȏꍇ�������̂Ł~
	//if( Apq.Determinant() < 0.0 && m_iNumVertices >= 4)
	//{
		//cout << "before det < 0" << endl;
		//�P�@�����𔽓]
		//Apq(0,2) = -Apq(0,2);
		//Apq(1,2) = -Apq(1,2);
		//Apq(2,2) = -Apq(2,2);

		//�Q�@a2��a3������
		//double tmp;
		//tmp = Apq(0,2);
		//Apq(0,2) = Apq(0,1);
		//Apq(0,1) = tmp;

		//tmp = Apq(1,2);
		//Apq(1,2) = Apq(1,1);
		//Apq(1,1) = tmp;

		//tmp = Apq(2,2);
		//Apq(2,2) = Apq(2,1);
		//Apq(2,1) = tmp;
	//}

	rxMatrix3 R, S;
	PolarDecomposition(Apq, R, S);

	if(m_bLinearDeformation){
		// Linear Deformations
		rxMatrix3 A;
		//A = Apq*Aqq.Inverse();	// A = Apq*Aqq^-1
		A = Apq*Aqq;	// A = Apq*Aqq^-1

		// �̐ϕۑ��̂��߂Ɂ�(det(A))�Ŋ���
		if(m_bVolumeConservation){
			double det = fabs(A.Determinant());
			if(det > RX_FEQ_EPS){
				det = 1.0/sqrt(det);
				if(det > 2.0) det = 2.0;
				A *= det;
			}
		}

		// �ڕW���W���v�Z���C���݂̒��_���W���ړ�
		for(int i = 0; i < m_iNumVertices; ++i){
			if(m_vFix[i]) continue;
			// ��]�s��R�̑���̍s��RL=��A+(1-��)R���v�Z
			rxMatrix3 RL = m_dBetas[i]*A+(1.0-m_dBetas[i])*R;
			q = m_vOrgPos[i]-cm_org;
			m_vGoalPos[i] = RL*q+cm;
			m_vNewPos[i] += (m_vGoalPos[i]-m_vNewPos[i])*m_dAlphas[i];
		}
	}
	else{
		// Quadratic Deformations
		double Atpq[3][9];
		for(int j = 0; j < 9; ++j){
			Atpq[0][j] = 0.0;
			Atpq[1][j] = 0.0;
			Atpq[2][j] = 0.0;
		}
		rxMatrixN<double,9> Atqq;
		Atqq.SetValue(0.0);

		for(int i = 0; i < m_iNumVertices; ++i){
			p = m_vNewPos[i]-cm;
			q = m_vOrgPos[i]-cm_org;

			// q~�̌v�Z
			double qt[9];
			qt[0] = q[0];      qt[1] = q[1];      qt[2] = q[2];
			qt[3] = q[0]*q[0]; qt[4] = q[1]*q[1]; qt[5] = q[2]*q[2];
			qt[6] = q[0]*q[1]; qt[7] = q[1]*q[2]; qt[8] = q[2]*q[0];

			// A~pq = ��mpq~ �̌v�Z
			double m = m_vMass[i];
			for(int j = 0; j < 9; ++j){
				Atpq[0][j] += m*p[0]*qt[j];
				Atpq[1][j] += m*p[1]*qt[j];
				Atpq[2][j] += m*p[2]*qt[j];
			}

			// A~qq = ��mq~q~ �̌v�Z
			for(int j = 0; j < 9; ++j){
				for(int k = 0; k < 9; ++k){
					Atqq(j,k) += m*qt[j]*qt[k];
				}
			}
		}

		// A~qq�̋t�s��Z�o
		Atqq.Invert();

		double At[3][9];
		for(int i = 0; i < 3; ++i){
			for(int j = 0; j < 9; j++){
				At[i][j] = 0.0f;
				for(int k = 0; k < 9; k++){
					At[i][j] += Atpq[i][k]*Atqq(k,j);
				}

				// ��A~+(1-��)R~
				At[i][j] *= m_dBetas[i];
				if(j < 3){
					At[i][j] += (1.0f-m_dBetas[i])*R(i,j);
				}
			}
		}

		 // a00a11a22+a10a21a02+a20a01a12-a00a21a12-a20a11a02-a10a01a22
		double det = At[0][0]*At[1][1]*At[2][2]+At[1][0]*At[2][1]*At[0][2]+At[2][0]*At[0][1]*At[1][2]
					-At[0][0]*At[2][1]*At[1][2]-At[2][0]*At[1][1]*At[0][2]-At[1][0]*At[0][1]*At[2][2];

		// �̐ϕۑ��̂��߂Ɂ�(det(A))�Ŋ���
		if(m_bVolumeConservation){
			if(det != 0.0f){
				det = 1.0f/sqrt(fabs(det));
				if(det > 2.0f) det = 2.0f;
				for(int i = 0; i < 3; ++i){
					for(int j = 0; j < 3; ++j){
						At[i][j] *= det;
					}
				}
			}
		}


		// �ڕW���W���v�Z���C���݂̒��_���W���ړ�
		for(int i = 0; i < m_iNumVertices; ++i){
			if(m_vFix[i]) continue;
			q = m_vOrgPos[i]-cm_org;

			for(int k = 0; k < 3; ++k){
				m_vGoalPos[i][k] = At[k][0]*q[0]+At[k][1]*q[1]+At[k][2]*q[2]
								  +At[k][3]*q[0]*q[0]+At[k][4]*q[1]*q[1]+At[k][5]*q[2]*q[2]+
								  +At[k][6]*q[0]*q[1]+At[k][7]*q[1]*q[2]+At[k][8]*q[2]*q[0];
			}

			m_vGoalPos[i] += cm;
			m_vNewPos[i] += (m_vGoalPos[i]-m_vNewPos[i])*m_dAlphas[i];
		}
	}

	//newTime = clock();
	//cout << "�v���I��::" << (double)(newTime - oldTime)/CLOCKS_PER_SEC << endl;
}


/*!
 * �O��
 *  - �d�͂Ƌ��E�ǂ���̗͂̉e��
 * @param[in] dt �^�C���X�e�b�v��
 */
void Ice_SM::calExternalForces(double dt)
{
	// �d�͂̉e����t���C���x�𔽉f
	for(int i = 0; i < m_iNumVertices; ++i){
		if(m_vFix[i]) continue;
		//m_vVel[i] += m_v3Gravity*dt;
		m_vNewPos[i] = m_vCurPos[i]+m_vVel[i]*dt;
		m_vGoalPos[i] = m_vOrgPos[i];
	}

	// ���E�ǂ̉e��
	double res = 0.5;	// �����W��
	for(int i = 0; i < m_iNumVertices; ++i){
		//if(m_vFix[i]) continue;
		//Vec3 &p = m_vCurPos[i];
		//Vec3 &np = m_vNewPos[i];
		//Vec3 &v = m_vVel[i];
		//if(np[0] < m_v3Min[0] || np[0] > m_v3Max[0]){
		//	np[0] = p[0]-v[0]*dt*res;
		//	np[1] = p[1];
		//	np[2] = p[2];
		//}
		//if(np[1] < m_v3Min[1] || np[1] > m_v3Max[1]){
		//	np[1] = p[1]-v[1]*dt*res;
		//	np[0] = p[0];
		//	np[2] = p[2];
		//}
		//if(np[2] < m_v3Min[2] || np[2] > m_v3Max[2]){
		//	np[2] = p[2]-v[2]*dt*res;
		//	np[0] = p[0];
		//	np[1] = p[1];
		//}
		//clamp(m_vNewPos[i]);
	}
}

/*!
 * ���x�ƈʒu�̍X�V
 *  - �V�����ʒu�ƌ��݂̈ʒu���W���瑬�x���Z�o
 * @param[in] dt �^�C���X�e�b�v��
 */
void Ice_SM::integrate(double dt)
{
	double dt1 = 1.0/dt;
	for(int i = 0; i < m_iNumVertices; ++i){
		m_vVel[i] = (m_vNewPos[i]-m_vCurPos[i])*dt1;
		m_vCurPos[i] = m_vNewPos[i];
	}
}

/*!
 * �V�~�����[�V�����X�e�b�v��i�߂�
 */
void Ice_SM::Update()
{//	cout << "Ice_update" << endl;
	//calCollision(m_dDt);
	calExternalForces(m_dDt);
	ShapeMatching(m_dDt);
	integrate(m_dDt);
}


//----------------------------------------�f�o�b�O--------------------------------------
void Ice_SM::DebugIndx()
{	cout << "DebugIndx Indx = ";
	for(int i = 0; i < m_iNumVertices; i++)
	{
		cout << " " << m_iParticleIndxes[i];
	}
	cout << endl;
}

void Ice_SM::DebugLayer()
{	cout << "DebugLayer layer =";
	for(int i = 0; i < m_iNumVertices; i++)
	{
		cout << " " << m_iLayeres[i];
	}
	cout << endl;
}