#include "Ice_CalcStiffData_CompareRigid.h"

typedef Ice_CalcStiffData_CompareRigid StiffDivision;

StiffDivision::Ice_CalcStiffData_CompareRigid(const vector<Ice_SM*>& iceSM, Ice_JudgeMove* judge)
{
	m_iceSM = iceSM;
	m_iceJudge = judge;

	MakeRigidObj();
}

StiffDivision::~Ice_CalcStiffData_CompareRigid()
{
	delete[] m_rigid;
}

void StiffDivision::MakeRigidObj()
{
	//�I�u�W�F�N�g�쐬
	m_rigid = new Ice_SM(MAXINT, m_iceSM.size());
	Vec3 areaMin = m_iceSM[0]->areaMin();
	Vec3 areaMax = m_iceSM[0]->areaMax();
	m_rigid->SetSimulationSpace(areaMin, areaMax);

	float timeStep = m_iceSM[0]->dt();
	m_rigid->SetTimeStep(timeStep);
	m_rigid->SetCollisionFunc(0);
	m_rigid->SetStiffness(1.0, 0.0);

	//���q�ǉ�
	for(vector<Ice_SM*>::iterator it = m_iceSM.begin(); it != m_iceSM.end(); it++){
		//�擪�̗��q�������W�߂�
		int pIndx = (*it)->GetParticleIndx(0);

		Vec3 pos = (*it)->GetOrgPos(0);
		Vec3 vel = (*it)->GetVertexVel(0);

		int pNum = m_rigid->GetNumVertices();
		m_rigid->AddVertex(pos, vel, 1.0f, pIndx);
		m_rigid->SetAlphas(pNum, 1.0);
		m_rigid->SetBetas (pNum, 0.0);
		m_rigid->SetLayer (pNum, 0);
	}

	cout << __FUNCTION__ << " size:" << m_rigid->GetNumVertices() << endl;
}

//rigid�N���X�^�̉^���v�Z
void StiffDivision::StepUpdate()
{
	unsigned pNum = m_iceSM.size();
	const float* sldPos = Ice_SM::GetSldPosPointer();
	const float* sldVel = Ice_SM::GetSldVelPointer();

	//rigid�̌��݈ʒu���ŏI�v�Z���ʂōX�V�@���̃t���[���̓��͂ƂȂ�
	for(unsigned pIndx = 0; pIndx < pNum; pIndx++){
		int sldIndx = pIndx * SM_DIM;

		Vec3 pos(sldPos[sldIndx+0], sldPos[sldIndx+1], sldPos[sldIndx+2]);
		Vec3 vel(sldVel[sldIndx+0], sldVel[sldIndx+1], sldVel[sldIndx+2]);

		//rigid�̌��݈ʒu���ŏI�v�Z���ʂōX�V�@���̃t���[���̓��͂ƂȂ�
		//TODO: ����͒e���݂̂̂̉^���v�Z��z�肵�Ă���@�t�̂Ƃ̕�Ԃ͍l���Ă��Ȃ�
		m_rigid->SetCurrentPos(pIndx, pos);
		m_rigid->SetNewPos(pIndx, pos);

		//�d�͂����f�����̂́CIce_CalcMethod::StepObjMove�̏���̂�
		m_rigid->SetVelocity(pIndx, vel);
	}

	//rigid���^���v�Z���Č��݈ʒu���X�V
	m_rigid->Update();
}

void StiffDivision::StepUpdateItr()
{
	unsigned pNum = m_iceSM.size();
	const float* sldPos = Ice_SM::GetSldPosPointer();
	const float* sldVel = Ice_SM::GetSldVelPointer();

	//rigid�̌��݈ʒu���ŏI�v�Z���ʂōX�V�@���̃t���[���̓��͂ƂȂ�
	for(unsigned pIndx = 0; pIndx < pNum; pIndx++){
		int sldIndx = pIndx * SM_DIM;

		Vec3 pos(sldPos[sldIndx+0], sldPos[sldIndx+1], sldPos[sldIndx+2]);
		Vec3 vel(sldVel[sldIndx+0], sldVel[sldIndx+1], sldVel[sldIndx+2]);

		//rigid�̌��݈ʒu���ŏI�v�Z���ʂōX�V�@���̃t���[���̓��͂ƂȂ�
		//TODO: ����͒e���݂̂̂̉^���v�Z��z�肵�Ă���@�t�̂Ƃ̕�Ԃ͍l���Ă��Ȃ�
		m_rigid->SetCurrentPos(pIndx, pos);
		m_rigid->SetNewPos(pIndx, pos);

		m_rigid->SetVelocity(pIndx, vel);
	}

	//rigid���^���v�Z���Č��݈ʒu���X�V
	double dt = m_rigid->dt();
	m_rigid->rxShapeMatching::calcBoundary();
	m_rigid->rxShapeMatching::shapeMatching(dt);
	m_rigid->rxShapeMatching::integrate(dt);
}

//���̂̉^���v�Z���ʂƂ̍����e���q�ŋ��߁C���̑��a��Ԃ�
float StiffDivision::StepCalcData()
{
	unsigned pNum = m_iceSM.size();
	const float* sldPos = Ice_SM::GetSldPosPointer();
	const float* sldVel = Ice_SM::GetSldVelPointer();

	//rigid�ƃN���X�^�̍ŏI�v�Z���ʂ̍��������߁C���a�v�Z
	float divSum = 0.0f;
	for(unsigned pIndx = 0; pIndx < pNum; pIndx++){
		int sldIndx = pIndx * SM_DIM;

		Vec3 pos(sldPos[sldIndx+0], sldPos[sldIndx+1], sldPos[sldIndx+2]);
		Vec3 vel(sldVel[sldIndx+0], sldVel[sldIndx+1], sldVel[sldIndx+2]);

		Vec3 rpos = m_rigid->GetVertexPos(pIndx);

		divSum += (float)length(pos-rpos);
	}

	//�l��Ԃ�
	return divSum;
}

void StiffDivision::StepCalcDataDebug()
{
	unsigned pNum = m_rigid->GetNumVertices();
	const float* sldPos = Ice_SM::GetSldPosPointer();
	const float* sldVel = Ice_SM::GetSldVelPointer();

	//rigid�̌��݈ʒu��_�ŕ`��
	glColor4f(0.0f,0.0f,1.0f,1.0f);
	glPointSize(10);
	glBegin(GL_POINTS);
	
	for(unsigned pIndx = 0; pIndx < pNum; pIndx++){
		Vec3 pos = m_rigid->GetVertexPos(pIndx);
		glVertex3d(pos[0], pos[1], pos[2]);
	}

	glEnd();

	////�e�X�g�@�{���ɍ��̂ɂȂ邩����
	//for(unsigned pIndx = 0; pIndx < m_iceSM.size(); pIndx++){
	//	for(unsigned i = 0; i < m_iceSM[pIndx]->GetIndxNum(); i++){
	//		int ipIndx = m_iceSM[pIndx]->GetParticleIndx(i);
	//		if(ipIndx == MAXINT){
	//			continue;
	//		}

	//		Vec3 pos = m_rigid->GetVertexPos(ipIndx);
	//		Vec3 vel = m_rigid->GetVertexVel(ipIndx);

	//		m_iceSM[pIndx]->SetCurrentPos(i, pos);
	//		m_iceSM[pIndx]->SetVelocity(i, vel);
	//	}
	//}
}