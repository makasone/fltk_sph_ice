/*!
  @file rx_shape_matching.h
	
  @brief Shape Matching�@�ɂ��e���ό`
  @ref M. Muller et al., "Meshless deformations based on shape matching", SIGGRAPH2005. 
 
  @author Makoto Fujisawa
  @date 2013-07
*/

#ifndef _RX_SHAPE_MATCHING_H_
#define _RX_SHAPE_MATCHING_H_


//-----------------------------------------------------------------------------
// �C���N���[�h�t�@�C��
//-----------------------------------------------------------------------------
#include <cstdio>
#include <cmath>
#include <cstdlib>

#include <vector>
#include <string>

#include "rx_utility.h"
#include "rx_matrix.h"

#include "rx_nnsearch.h"

#include <cuda_runtime.h>
#include <cuda_gl_interop.h>

using namespace std;

//! �Փ˔���p�֐�
typedef void (*CollisionFunc)(Vec3&, Vec3&, Vec3&, int);


//GPU����
extern void LaunchShapeMatchingGPU(
	/*unsigned int num_particles,
	float (*pos)[2],
	float time,
	float dt*/);

#define MAXPARTICLE 100
#define SM_DIM 3


//-----------------------------------------------------------------------------
// Shape Matching�N���X�̐錾
//-----------------------------------------------------------------------------
class rxShapeMatching
{
protected:
	// �`��f�[�^
	double* m_pOrgPos;				//!< �I���W�i���̒��_�ʒu
	double* m_pCurPos;				//!< ���݂̒��_�ʒu
	double* m_pNewPos;				//!< ���̃X�e�b�v�̒��_�ʒu
	double* m_pGoalPos;				//!< �ڕW���_�ʒu
	double* m_pMass;				//!< ���_����(�ό`���̏d��)
	double* m_pVel;					//!< ���_���x

	bool* m_pFix;					//!< ���_�Œ�t���O

	int m_iNumVertices;

	// �V�~�����[�V�����p�����[�^
	double m_dDt;					//!< �^�C���X�e�b�v��
	Vec3 m_v3Min, m_v3Max;			//!< �V�~�����[�V������Ԃ̑傫��
	Vec3 m_v3Gravity;				//!< �d�͉����x�x�N�g��
	
	//�ǉ�
	rxMatrix3 m_mtrxBeforeU;		//�O�t���[���ł̒����s��@warm start�Ƃ�����������@
	rxMatrix3 m_mtrxBeforeATA;

	double m_dAlpha;				//!< stiffness�p�����[�^[0,1] (���x�v�Z�Ɏg�p)
	double m_dBeta;					//!< deformation�p�����[�^[0,1]

	bool m_bLinearDeformation;		//!< Linear/Quadratic deformation�؂�ւ��t���O
	bool m_bVolumeConservation;		//!< �ό`���̑̐ϕۑ���(��det(A)�Ŋ��邩�ǂ���)

	int m_iObjectNo;				//!< �I�u�W�F�N�g�ԍ�

	CollisionFunc m_fpCollision;

	//GPU
	//�f�o�C�X���ւ̃|�C���^
	double* d_OrgPos;
	double* d_CurPos;
	double* d_NewPos;
	double* d_GoalPos;
	double* d_Mass;
	double* d_Vel;

	bool* d_Fix;

public:
	//! �R���X�g���N�^�ƃf�X�g���N�^
	rxShapeMatching(int obj);
	~rxShapeMatching();

	void InitGPU();

	void Clear();
	void AddVertex(const Vec3 &pos, double mass);

	void Update();

	// �A�N�Z�X���\�b�h
	void SetTimeStep(double dt){ m_dDt = dt; }
	void SetSimulationSpace(Vec3 minp, Vec3 maxp){ m_v3Min = minp; m_v3Max = maxp; }
	void SetStiffness(double alpha, double beta){ m_dAlpha = alpha; m_dBeta = beta; }

	void SetCurrentPos(int oIndx, const Vec3 &pos)
	{ 
		m_pCurPos[oIndx*SM_DIM+0] = pos[0];
		m_pCurPos[oIndx*SM_DIM+1] = pos[1];
		m_pCurPos[oIndx*SM_DIM+2] = pos[2];
	}

	void SetOriginalPos(int oIndx, Vec3 pos)
	{ 
		m_pOrgPos[oIndx*SM_DIM+0] = pos[0];
		m_pOrgPos[oIndx*SM_DIM+1] = pos[1];
		m_pOrgPos[oIndx*SM_DIM+2] = pos[2];
	}

	void SetNewPos(int oIndx, Vec3 pos)
	{
		m_pNewPos[oIndx*SM_DIM+0] = pos[0];
		m_pNewPos[oIndx*SM_DIM+1] = pos[1];
		m_pNewPos[oIndx*SM_DIM+2] = pos[2];
	}

	void SetGoalPos(int oIndx, Vec3 pos)
	{
		m_pGoalPos[oIndx*SM_DIM+0] = pos[0];
		m_pGoalPos[oIndx*SM_DIM+1] = pos[1];
		m_pGoalPos[oIndx*SM_DIM+2] = pos[2];
	}

	void SetVelocity(int oIndx, const Vec3 &vel)
	{
		m_pVel[oIndx*SM_DIM+0] = vel[0];
		m_pVel[oIndx*SM_DIM+1] = vel[1];
		m_pVel[oIndx*SM_DIM+2] = vel[2];
	}

	void SetCollisionFunc(CollisionFunc func){ m_fpCollision = func; }

	int GetNumVertices() const { return m_iNumVertices; }

	const Vec3 GetVertexPos(int i){ return Vec3(m_pCurPos[i*SM_DIM+0], m_pCurPos[i*SM_DIM+1], m_pCurPos[i*SM_DIM+2]); }
	const Vec3 GetNewPos(int i){ return Vec3(m_pNewPos[i*SM_DIM+0], m_pNewPos[i*SM_DIM+1], m_pNewPos[i*SM_DIM+2]); }
	const Vec3 GetOrgPos(int i){ return Vec3(m_pOrgPos[i*SM_DIM+0], m_pOrgPos[i*SM_DIM+1], m_pOrgPos[i*SM_DIM+2]); }
	const Vec3 GetGoalPos(int i) { return Vec3(m_pGoalPos[i*SM_DIM+0], m_pGoalPos[i*SM_DIM+1], m_pGoalPos[i*SM_DIM+2]); }
	const Vec3 GetVertexVel(int i){ return Vec3(m_pVel[i*SM_DIM+0], m_pVel[i*SM_DIM+1], m_pVel[i*SM_DIM+2]); }
	double GetMass(int i){ return m_pMass[i]; }

	void FixVertex(int i, const Vec3 &pos);
	void UnFixVertex(int i);
	bool IsFixed(int i) { return m_pFix[i]; }

protected:
	//! ���_�ʒu�̏�����
	void initialize(void);

	// Shape Matching�@�̌v�Z
	void calExternalForces(double dt);
	void calCollision(double dt);
	void shapeMatching(double dt);
	void integrate(double dt);

	void clamp(Vec3 &pos) const
	{
		if(pos[0] < m_v3Min[0]) pos[0] = m_v3Min[0];
		if(pos[0] > m_v3Max[0]) pos[0] = m_v3Max[0];
		if(pos[1] < m_v3Min[1]) pos[1] = m_v3Min[1];
		if(pos[1] > m_v3Max[1]) pos[1] = m_v3Max[1];
		if(pos[2] < m_v3Min[2]) pos[2] = m_v3Min[2];
		if(pos[2] > m_v3Max[2]) pos[2] = m_v3Max[2];
	}

	void clamp(double* pos, int cIndx) const
	{
		if(pos[cIndx+0] < m_v3Min[0]) pos[cIndx+0] = m_v3Min[0];
		if(pos[cIndx+0] > m_v3Max[0]) pos[cIndx+0] = m_v3Max[0];
		if(pos[cIndx+1] < m_v3Min[1]) pos[cIndx+1] = m_v3Min[1];
		if(pos[cIndx+1] > m_v3Max[1]) pos[cIndx+1] = m_v3Max[1];
		if(pos[cIndx+2] < m_v3Min[2]) pos[cIndx+2] = m_v3Min[2];
		if(pos[cIndx+2] > m_v3Max[2]) pos[cIndx+2] = m_v3Max[2];
	}

};




/*!
 * Jacobi�@�ɂ��ŗL�l�̎Z�o
 * @param[inout] a ���Ώ̍s��D�v�Z��C�Ίp�v�f�ɌŗL�l������
 * @param[out] v �ŗL�x�N�g��(a�Ɠ����T�C�Y)
 * @param[in] n �s��̃T�C�Y(n�~n)
 * @param[in] eps �����덷
 * @param[in] iter_max �ő唽����
 * @return ������
 */
inline int EigenJacobiMethod(double *a, double *v, int n, double eps = 1e-8, int iter_max = 100)
{
	double *bim, *bjm;
	double bii, bij, bjj, bji;
 
	bim = new double[n];
	bjm = new double[n];
 
	//for(int i = 0; i < n; ++i){
	//	for(int j = 0; j < n; ++j){
	//		v[i*n+j] = (i == j) ? 1.0 : 0.0;
	//	}
	//}
 
	int cnt = 0;
	for(;;){
		int i = -1, j = -1;
 
		double x = 0.0;
		for(int ia = 0; ia < n; ++ia){
			for(int ja = 0; ja < n; ++ja){
				int idx = ia*n+ja;
				if(ia != ja && fabs(a[idx]) > x){
					i = ia;
					j = ja;
					x = fabs(a[idx]);
				}
			}
		}

		if(i == -1 || j == -1) return 0;
 
		double aii = a[i*n+i];
		double ajj = a[j*n+j];
		double aij = a[i*n+j];
 
		double m_dAlpha, m_dBeta;
		m_dAlpha = (aii-ajj)/2.0;
		m_dBeta  = sqrt(m_dAlpha*m_dAlpha+aij*aij);
 
		double st, ct;
		ct = sqrt((1.0+fabs(m_dAlpha)/m_dBeta)/2.0);    // sin��
		st = (((aii-ajj) >= 0.0) ? 1.0 : -1.0)*aij/(2.0*m_dBeta*ct);    // cos��
 
		// A = PAP�̌v�Z
		for(int m = 0; m < n; ++m){
			if(m == i || m == j) continue;
 
			double aim = a[i*n+m];
			double ajm = a[j*n+m];
 
			bim[m] =  aim*ct+ajm*st;
			bjm[m] = -aim*st+ajm*ct;
		}
 
		bii = aii*ct*ct+2.0*aij*ct*st+ajj*st*st;
		bij = 0.0;
 
		bjj = aii*st*st-2.0*aij*ct*st+ajj*ct*ct;
		bji = 0.0;
 
		for(int m = 0; m < n; ++m){
			a[i*n+m] = a[m*n+i] = bim[m];
			a[j*n+m] = a[m*n+j] = bjm[m];
		}
		a[i*n+i] = bii;
		a[i*n+j] = bij;
		a[j*n+j] = bjj;
		a[j*n+i] = bji;
 
		// V = PV�̌v�Z
		for(int m = 0; m < n; ++m){
			double vmi = v[m*n+i];
			double vmj = v[m*n+j];
 
			bim[m] =  vmi*ct+vmj*st;
			bjm[m] = -vmi*st+vmj*ct;
		}
		for(int m = 0; m < n; ++m){
			v[m*n+i] = bim[m];
			v[m*n+j] = bjm[m];
		}
 
		double e = 0.0;
		for(int ja = 0; ja < n; ++ja){
			for(int ia = 0; ia < n; ++ia){
				if(ia != ja){
					e += fabs(a[ja*n+ia]);
				}
			}
		}
		if(e < eps) break;
 
		cnt++;
		if(cnt > iter_max) break;
	}
 
	delete [] bim;
	delete [] bjm;
 
	return cnt;
}


/*!
 * �ɕ����ŉ�]�s��ƑΏ̍s��ɕ��� A=RS
 * @param[in] A ���͍s��
 * @param[out] R ��]�s��(�����s�� R^-1 = R^T)
 * @param[out] S �Ώ̍s��
 */
inline void PolarDecomposition(const rxMatrix3 &A, rxMatrix3 &R, rxMatrix3 &S)
{
	// S = (A^T A)^(1/2)�����߂�
	rxMatrix3 ATA;
	// (A^T A)�̌v�Z
	ATA = A.Transpose()*A;

	rxMatrix3 U;
	R.makeIdentity();
	U.makeIdentity();

	// (A^T A)���ŗL�l�������đΊp�s��ƒ����s������߂�
	// M^(1/2) = U^T M' U 
	//  M = (A^T A), M':�Ίp�s��̕����������������, U:�����s��, 
	EigenJacobiMethod(&ATA, &U, 3);

	// �Ίp�s��̕��������Ƃ��āC�t�s��v�Z�̂��߂ɋt���ɂ��Ă���
	real l0 = (ATA(0,0) <= 0.0) ? 0.0 : 1.0/sqrt(ATA(0,0));
	real l1 = (ATA(1,1) <= 0.0) ? 0.0 : 1.0/sqrt(ATA(1,1));
	real l2 = (ATA(2,2) <= 0.0) ? 0.0 : 1.0/sqrt(ATA(2,2));

	// U^T M' U �̋t�s��v�Z
	rxMatrix3 S1;
	S1(0,0) = l0*U(0,0)*U(0,0) + l1*U(0,1)*U(0,1) + l2*U(0,2)*U(0,2);
	S1(0,1) = l0*U(0,0)*U(1,0) + l1*U(0,1)*U(1,1) + l2*U(0,2)*U(1,2);
	S1(0,2) = l0*U(0,0)*U(2,0) + l1*U(0,1)*U(2,1) + l2*U(0,2)*U(2,2);
	S1(1,0) = S1(0,1);
	S1(1,1) = l0*U(1,0)*U(1,0) + l1*U(1,1)*U(1,1) + l2*U(1,2)*U(1,2);
	S1(1,2) = l0*U(1,0)*U(2,0) + l1*U(1,1)*U(2,1) + l2*U(1,2)*U(2,2);
	S1(2,0) = S1(0,2);
	S1(2,1) = S1(1,2);
	S1(2,2) = l0*U(2,0)*U(2,0) + l1*U(2,1)*U(2,1) + l2*U(2,2)*U(2,2);

	R = A*S1;	// R = A S^-1
	S = R.Transpose()*A; // S = R^-1 A = R^T A
}

/*!
 * �ɕ����ŉ�]�s��ƑΏ̍s��ɕ��� A=RS
 * warm start����	1�t���[���ڂ͒P�ʍs��@���ʂ����������͕s��
 * @param[in] A ���͍s��
 * @param[out] R ��]�s��(�����s�� R^-1 = R^T)
 * @param[out] S �Ώ̍s��
 */
inline void PolarDecomposition(const rxMatrix3 &A, rxMatrix3 &R, rxMatrix3 &S, rxMatrix3& bfrU)
{
	// S = (A^T A)^(1/2)�����߂�
	//warm start	1�t���[���ڂ͒P�ʍs��
	// (A^T A)�̌v�Z
	rxMatrix3 ATA(bfrU.Transpose()*A.Transpose()*A*bfrU);

	R.makeIdentity();

	// (A^T A)���ŗL�l�������đΊp�s��ƒ����s������߂�
	// M^(1/2) = U^T M' U 
	//  M = (A^T A), M':�Ίp�s��̕����������������, U:�����s��, 
	EigenJacobiMethod(&ATA, &bfrU, 3);

	// �Ίp�s��̕��������Ƃ��āC�t�s��v�Z�̂��߂ɋt���ɂ��Ă���
	real l0 = (ATA(0,0) <= 0.0) ? 0.0 : 1.0/sqrt(ATA(0,0));
	real l1 = (ATA(1,1) <= 0.0) ? 0.0 : 1.0/sqrt(ATA(1,1));
	real l2 = (ATA(2,2) <= 0.0) ? 0.0 : 1.0/sqrt(ATA(2,2));

	// U^T M' U �̋t�s��v�Z
	rxMatrix3 S1;
	S1(0,0) = l0*bfrU(0,0)*bfrU(0,0) + l1*bfrU(0,1)*bfrU(0,1) + l2*bfrU(0,2)*bfrU(0,2);
	S1(0,1) = l0*bfrU(0,0)*bfrU(1,0) + l1*bfrU(0,1)*bfrU(1,1) + l2*bfrU(0,2)*bfrU(1,2);
	S1(0,2) = l0*bfrU(0,0)*bfrU(2,0) + l1*bfrU(0,1)*bfrU(2,1) + l2*bfrU(0,2)*bfrU(2,2);
	S1(1,0) = S1(0,1);
	S1(1,1) = l0*bfrU(1,0)*bfrU(1,0) + l1*bfrU(1,1)*bfrU(1,1) + l2*bfrU(1,2)*bfrU(1,2);
	S1(1,2) = l0*bfrU(1,0)*bfrU(2,0) + l1*bfrU(1,1)*bfrU(2,1) + l2*bfrU(1,2)*bfrU(2,2);
	S1(2,0) = S1(0,2);
	S1(2,1) = S1(1,2);
	S1(2,2) = l0*bfrU(2,0)*bfrU(2,0) + l1*bfrU(2,1)*bfrU(2,1) + l2*bfrU(2,2)*bfrU(2,2);

	R = A*S1;	// R = A S^-1
	S = R.Transpose()*A; // S = R^-1 A = R^T A
}

#endif
