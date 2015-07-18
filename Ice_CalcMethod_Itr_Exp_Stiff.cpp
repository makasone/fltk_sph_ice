#include "Ice_CalcMethod_Itr_Stiffness.h"

typedef Ice_CalcMethod_Itr_Exp_Stiff CalcIteration;

#define PNUM 6525
#define ADDLIMIT 50

CalcIteration::Ice_CalcMethod_Itr_Exp_Stiff(const vector<Ice_SM*>& iceSM, Ice_SimuMethod* simuMethod, Ice_Convolution* convo)
{
	m_iceSM = iceSM;
	//SetObjMove(clusterMove);
	SetObjMove(simuMethod);
	SetConvolution(convo);

	//臒l�p�f�[�^�Z�o�N���X
	//m_iceCalcStiff = new Ice_CalcStiffData_Summation(iceSM);	//���ό`��
	m_iceCalcStiff = new Ice_CalcStiffData_Average(iceSM);		//���ϕω���
	//m_iceCalcStiff = new Ice_CalcStiffData_StdDevision(iceSM, m_iceMove->GetJudgeMove());		//���q�ʒu�̕��U
	//m_iceCalcStiff = new Ice_CalcStiffData_CompareRigid(iceSM, m_iceMove->GetJudgeMove());	//���̂Ƃ̍���
}

CalcIteration::~Ice_CalcMethod_Itr_Exp_Stiff()
{
}

//void CalcIteration::SetObjMove(Ice_ClusterMove* clusterMove)
void CalcIteration::SetObjMove(Ice_SimuMethod* simuMethod)
{
	m_simuMethod = simuMethod;
}

void CalcIteration::SetConvolution(Ice_Convolution* convo)
{
	m_iceConvo = convo;
}

void CalcIteration::StepObjMove()
{
	//����̉^���v�Z
	//m_iceMove->StepObjMove();
	m_simuMethod->StepObjMove();

	//rigid�̍ŏI���ʁ@�d�͂𔽉f���Ă���o�[�W����
	m_iceCalcStiff->StepUpdate();

	//�^���v�Z
	m_iceConvo->StepConvolution();

	//�e�N���X�^�̃I���W�i�����q���o���Ă���
	//�i���́C�P�Ȃ�g�������̏����Ȃ�I���W�i�����q�̌������ł悢�j
	vector<vector<unsigned>> copyOriginalParticleIndxes(m_iceSM.size(), vector<unsigned>());
	CopyOriginalObject(copyOriginalParticleIndxes);

	//�e�������ɑ����ς݂̓Y���ԍ����L�^���Ă���
	vector<int> searchFinishIndxes(m_iceSM.size(), 0);

	//�������邩���肷�邽�߂̒l
	float threshold = m_iceCalcStiff->StepCalcData();

	//����	
	while(threshold > Ice_SM::GetItrStiffness()){
		//�N���X�^�̉e���͈͂��g�債�čč\��
		ExpandeCluster(searchFinishIndxes);

		//�������̉^���v�Z
		//m_iceMove->StepObjMoveItr();
		m_simuMethod->StepObjMoveItr();

		//�^���v�Z
		m_iceConvo->StepConvolution();

		//�������邩���肷�邽�߂̒l
		threshold = m_iceCalcStiff->StepCalcData();
	}

	//�R�s�[�𗘗p���ăN���X�^�̍\�������ɖ߂�
	ReplaceCluster(copyOriginalParticleIndxes);

	//���x�Z�o
	CalcVel();
}

//���x�Z�o
void CalcIteration::CalcVel()
{
	#pragma omp parallel for
	for(int i = 0; i < IceObject::GetParticleNum(); ++i)
	{	
		//if(m_iceMove->GetJudgeMove()->JudgeMove(i) == false){	continue;	}
		if(m_simuMethod->GetJudgeMove()->JudgeMove(i) == false){	continue;	}
		m_iceSM[i]->integrateIteration();
	}
}

//���̃I�u�W�F�N�g���R�s�[
void CalcIteration::CopyOriginalObject(vector<vector<unsigned>>& copyIndxes)
{
	for(unsigned i = 0; i < m_iceSM.size(); i++){
		for(unsigned j = 0; j < m_iceSM[i]->GetIndxNum(); j++){
			int jpIndx = m_iceSM[i]->GetParticleIndx(j);

			if(jpIndx == MAXINT){
				continue;
			}

			copyIndxes[i].push_back(jpIndx);
		}
	}
}

//�N���X�^�����ɖ߂�
void CalcIteration::ReplaceCluster(const vector<vector<unsigned>>& copyIndxes)
{
	//copyIndxes�ɑ��݂��闱�q�̓I���W�i��
	//���݂��Ȃ����q�̏���S�ď���
	for(unsigned i = 0; i < m_iceSM.size(); i++){
		for(unsigned j = 0; j < m_iceSM[i]->GetIndxNum(); j++){
			//���݂�m_iceSM�Ɋ܂܂�闱�q�ŁCcopyObj�ɑ��݂��Ȃ����̂��폜
			int jpIndx = m_iceSM[i]->GetParticleIndx(j);
			if(jpIndx == MAXINT){
				continue;
			}

			vector<unsigned>::const_iterator check = find(copyIndxes[i].begin(), copyIndxes[i].end(), jpIndx);
			if(check != copyIndxes[i].end()){
				continue;
			}

			//copyIndxes�ɑ��݂��Ȃ��̂Œǉ����ꂽ���q�@�폜����
			int ooIndx = m_iceSM[i]->SearchIndx(jpIndx);
			m_iceSM[i]->Remove(ooIndx);
		}

		m_iceSM[i]->CalcOrgCm();
		m_iceSM[i]->ClassifyAllOrgParticle();
	}
}

//�N���X�^�̉e���͈͂��g�債�čč\��
void CalcIteration::ExpandeCluster(vector<int>& sfFnIndxes)
{
	//�������Ԃ̌v��
QueryCounter counter1;
QueryCounter counter2;

counter1.Start();
	//�ߖT�ɑ��݂��闱�q�Ƃ��̏����ʒu�̃��X�g���쐬
	vector<vector<int>> addPIndxList;
	addPIndxList.resize(m_iceSM.size(), vector<int>(ADDLIMIT, MAXINT));

	for(unsigned i = 0; i < m_iceSM.size(); ++i){
		if(m_iceSM[i]->GetNumVertices() >= MAXPARTICLE){
			continue;
		}

		//���q��ǉ����Ă��邩�ǂ����̃t���O
		bool isAdd[PNUM] = {};
		int limitPrtNum = MAXPARTICLE - m_iceSM[i]->GetNumVertices();

		//���݂̃N���X�^�Ɋ܂܂�Ă��闱�q�͒ǉ����Ȃ�
		for(unsigned j = 0; j < m_iceSM[i]->GetIndxNum(); j++){
			int orgPrt = m_iceSM[i]->GetParticleIndx(j);

			//�ǉ����q���X�g�Ɋ܂܂�Ă��Ȃ���
			if(orgPrt == MAXINT || isAdd[orgPrt]){	
				continue;
			}
			
			isAdd[orgPrt] = true;
		}

		//�ߖT�N���X�^��T��
		//���ɒT���ς݂̃N���X�^�͏���
		unsigned searchStart = sfFnIndxes[i];

		//���݂̗��q�����Y����
		int addIndxNum = 0;

		for(unsigned clsI = searchStart; clsI < m_iceSM[i]->GetIndxNum(); clsI++){
			int nearCluster = m_iceSM[i]->GetParticleIndx(clsI);
			if(nearCluster == MAXINT || i == nearCluster){	
				continue;
			}

			//�N���X�^�Ɋ܂߂��闱�q���͍ő�300��
			if(addIndxNum >= limitPrtNum
			|| addIndxNum >= ADDLIMIT){
				break;
			}

			//jpIndx�̃N���X�^�Ɋ܂܂�闱�q���擾
			for(unsigned prtI = 0; prtI < m_iceSM[nearCluster]->GetIndxNum(); prtI++){
				//���̃N���X�^���Ō�܂ŒT�����I�����
				if(prtI == m_iceSM[nearCluster]->GetIndxNum()-1){
					sfFnIndxes[i]++;
				}

				int nrPrt = m_iceSM[nearCluster]->GetParticleIndx(prtI);

				//�ǉ����q���X�g�Ɋ܂܂�Ă��Ȃ���
				if(nrPrt == MAXINT || isAdd[nrPrt]){	
					continue;
				}

				//�N���X�^�Ɋ܂߂��闱�q���͍ő�300��
				if(addIndxNum >= limitPrtNum
				|| addIndxNum >= ADDLIMIT){
					break;
				}

				isAdd[nrPrt] = true;
				addPIndxList[i][addIndxNum++] = nrPrt;
			}
		}
	}
double end1 = counter1.End();

counter2.Start();
	//�N���X�^�ɗ��q��ǉ�
	const float* sldPos = Ice_SM::GetSldPosPointer();
	const float* sldVel = Ice_SM::GetSldVelPointer();

	for(unsigned i = 0; i < m_iceSM.size(); i++){

		bool isAdd[PNUM] = {};
		//�N���X�^�Ɋ܂܂�闱�q�̃t���O�I��
		for(unsigned j = 0; j < m_iceSM[i]->GetIndxNum(); j++){
			int orgPrtI = m_iceSM[i]->GetParticleIndx(j);
			if(orgPrtI == MAXINT || isAdd[orgPrtI]){
				continue;
			}

			isAdd[orgPrtI] = true;
		}

		for(vector<int>::iterator it = addPIndxList[i].begin(); it != addPIndxList[i].end(); it++){
			int addpIndx = *it;

			if(addpIndx == MAXINT || m_iceSM[i]->GetNumVertices() >= MAXPARTICLE){
				break;
			}

			//�N���X�^�ɗ��q�����݂��邩���m�F
			if(isAdd[addpIndx]){
				continue;
			}

			isAdd[addpIndx] = true;

			int sldIndx = addpIndx * SM_DIM;
			Vec3 orgPos(m_iceSM[addpIndx]->GetOrgPos(0));
			Vec3 pos(sldPos[sldIndx+0], sldPos[sldIndx+1], sldPos[sldIndx+2]);
			Vec3 vel(sldVel[sldIndx+0], sldVel[sldIndx+1], sldVel[sldIndx+2]);

			int check = m_iceSM[i]->AddAnotherClusterVertex(orgPos, pos, vel, 1.0, addpIndx, 1.0, 0.0, 0);
			if(check == -1){
				cout << "over " << addPIndxList[i].size() << endl;
			}
		}

		m_iceSM[i]->CalcOrgCm();
		m_iceSM[i]->ClassifyAllOrgParticle();
	}
double end2 = counter2.End();

	//cout << "Time:" << "Search:" << end1 << "," << "Add:" << end2 << endl;
}

void CalcIteration::StepObjMoveDebug()
{
	//����̉^���v�Z
	//m_iceMove->StepObjMove();
	m_simuMethod->StepObjMove();

	//rigid�̍ŏI���ʁ@�d�͂𔽉f���Ă���o�[�W����
	m_iceCalcStiff->StepUpdate();

	//�^���v�Z
	m_iceConvo->StepConvolution();

	//�e�N���X�^�̃I���W�i�����q���o���Ă���
	//�i���́C�P�Ȃ�g�������̏����Ȃ�I���W�i�����q�̌������ł悢�j
	vector<vector<unsigned>> copyOriginalParticleIndxes(m_iceSM.size(), vector<unsigned>());
	CopyOriginalObject(copyOriginalParticleIndxes);

	//�e�������ɑ����ς݂̓Y���ԍ����L�^���Ă���
	vector<int> searchFinishIndxes(m_iceSM.size(), 0);

	//�������邩���肷�邽�߂̒l
	float threshold = m_iceCalcStiff->StepCalcData();

	int count = 0;

	//����	
	while(threshold > Ice_SM::GetItrStiffness()){
		//�N���X�^�̉e���͈͂��g�債�čč\��
		ExpandeCluster(searchFinishIndxes);

		//�������̉^���v�Z
		//m_iceMove->StepObjMoveItr();
		m_simuMethod->StepObjMoveItr();

		//�^���v�Z
		m_iceConvo->StepConvolution();

		//�������邩���肷�邽�߂̒l
		threshold = m_iceCalcStiff->StepCalcData();

		count++;

		cout << __FUNCTION__ << " thrshold" << threshold << " count:" << count << endl;
	}

	//�R�s�[�𗘗p���ăN���X�^�̍\�������ɖ߂�
	ReplaceCluster(copyOriginalParticleIndxes);

	//���x�Z�o
	CalcVel();

//�f�o�b�O
	cout << __FUNCTION__ << " " << count << endl;

	string result = RESULT_DATA_PATH;
	result += "Itr_Stiffness_CountNum.txt";
	ofstream ofs(result, ios::app);
	ofs << count << endl;
}

//�d����臒l�𑪒肷�鏈���̃f�o�b�O
void CalcIteration::DebugStiffness()
{
	m_iceCalcStiff->StepCalcDataDebug();
}