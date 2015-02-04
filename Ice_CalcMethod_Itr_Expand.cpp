#include "Ice_CalcMethod_Itr_Expand.h"

typedef Ice_CalcMethod_Itr_Expand CalcIteration;

#define PNUM 2197

CalcIteration::Ice_CalcMethod_Itr_Expand(const vector<Ice_SM*>& iceSM, Ice_ClusterMove* clusterMove, Ice_Convolution* convo)
{
	m_iceSM = iceSM;
	SetObjMove(clusterMove);
	SetConvolution(convo);

	m_isAdd.resize(PNUM);
	ResetFlag();
}

CalcIteration::~Ice_CalcMethod_Itr_Expand()
{
}

void CalcIteration::SetObjMove(Ice_ClusterMove* clusterMove)
{
	m_iceMove = clusterMove;
}

void CalcIteration::SetConvolution(Ice_Convolution* convo)
{
	m_iceConvo = convo;
}

void CalcIteration::StepObjMove()
{
QueryCounter counter1;
QueryCounter counter2;
QueryCounter counter3;
QueryCounter counter4;

counter1.Start();
	//����̉^���v�Z
	m_iceMove->StepObjMove();
double end1 = counter1.End();

counter2.Start();
	//�e�N���X�^�̃I���W�i�����q���o���Ă���
	//�i���́C�P�Ȃ�g�������̏����Ȃ�I���W�i�����q�̌������ł悢�j
	vector<vector<unsigned>> copyOriginalParticleIndxes(m_iceSM.size(), vector<unsigned>());
	CopyOriginalObject(copyOriginalParticleIndxes);

	//�e�������ɑ����ς݂̓Y���ԍ����L�^���Ă���
	vector<int> searchFinishIndxes(m_iceSM.size(), 0);
double end2 = counter2.End();

counter3.Start();
	//����
	for(int itr = 1; itr < Ice_SM::GetIteration(); itr++){	
		//�ő̂ɂ�����ŏI�ʒu����
		m_iceConvo->StepConvolution();
	
		//�N���X�^�̉e���͈͂��g�債�čč\��
		ExpandeCluster(searchFinishIndxes);

		//�������̉^���v�Z
		m_iceMove->StepObjMoveItr();
	}
double end3 = counter3.End();

counter4.Start();
	//�R�s�[�𗘗p���ăN���X�^�̍\�������ɖ߂�
	ReplaceCluster(copyOriginalParticleIndxes);
double end4 = counter4.End();

	//���x�Z�o
	CalcVel();

	cout << "Time:"
		<< " first:" << end1
		<< " prepare:" << end2
		<< " itr:" << end3
		<< " replace:" << end4
		<< endl;
}

//���x�Z�o
void CalcIteration::CalcVel()
{
	#pragma omp parallel for
	for(int i = 0; i < IceObject::GetParticleNum(); ++i)
	{	
		if(m_iceMove->GetJudgeMove()->JudgeMove(i) == false){	continue;	}
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
	}
}

//�č\�z����N���X�^�����o
void CalcIteration::GetExpandeCluster()
{
	for(unsigned i = 0; i < m_iceSM.size(); ++i){

	}
}

//�t���O�̏�����
void CalcIteration::ResetFlag()
{
	for(vector<bool>::iterator it = m_isAdd.begin(); it != m_isAdd.end(); it++){
		*it = false;
	}
}

//�N���X�^�̉e���͈͂��g�債�čč\��
void CalcIteration::ExpandeCluster(vector<int>& srFnIndxes)
{
	//ExpandCluster_Test1();
	ExpandCluster_Test2(srFnIndxes);
}

void CalcIteration::ExpandCluster_Test1()
{
	//�������Ԃ̌v��
QueryCounter counter1;
QueryCounter counter2;

counter1.Start();
	////�t�@�C����ǂݍ���
	//string filePath = RESULT_DATA_PATH;
	//filePath += "ExpandParticle.txt";
	//ofstream ofs(filePath);

	////�t�@�C���̑��݊m�F
	//if(ofs.fail()) {	cerr << filePath << " is do not exist.\n";	return;	}

	//�ߖT�ɑ��݂��闱�q�Ƃ��̏����ʒu�̃��X�g���쐬
	vector<vector<pair<int, Vec3>>> addPIndxList;
	addPIndxList.resize(m_iceSM.size(), vector<pair<int, Vec3>>());

	for(unsigned i = 0; i < m_iceSM.size(); ++i){
		//���X�g�����݂̃N���X�^�Ɋ܂܂�Ă��闱�q�ŏ�����
		for(unsigned j = 0; j < m_iceSM[i]->GetIndxNum(); j++){
			int jpIndx = m_iceSM[i]->GetParticleIndx(j);
			if(jpIndx == MAXINT){	
				continue;
			}

			//�N���X�^�Ɋ܂߂��闱�q���͍ő�300��
			if(addPIndxList[i].size() > 299){
				//cout << "init max" << endl;
				return ;
			}

			Vec3 orgPos = m_iceSM[i]->GetOrgPos(j);
			addPIndxList[i].push_back(pair<int, Vec3>(jpIndx, orgPos));
		}

		//�ߖT�N���X�^��T��
		for(unsigned j = 0; j < m_iceSM[i]->GetIndxNum(); j++){
			int jpIndx = m_iceSM[i]->GetParticleIndx(j);
			if(jpIndx == MAXINT){	
				continue;
			}

			//jpIndx�̃N���X�^�Ɋ܂܂�闱�q���擾
			for(unsigned k = 0; k < m_iceSM[jpIndx]->GetIndxNum(); k++){
				int kpIndx = m_iceSM[jpIndx]->GetParticleIndx(k);
				if(kpIndx == MAXINT){	
					continue;
				}

				//�N���X�^�Ɋ܂߂��闱�q���͍ő�300��
				if(addPIndxList[i].size() > 299){
					//cout << "search max" << endl;
					break ;
				}

				//�ǉ����q���X�g�Ɋ܂܂�Ă��Ȃ��Ȃ�ǉ�����
				bool addCheck = false;
				for(vector<pair<int, Vec3>>::iterator it = addPIndxList[i].begin(); it != addPIndxList[i].end(); it++){
					int indx = it->first;

					if(kpIndx == indx){
						addCheck = true;
						break;
					}
				}

				if(addCheck){
					continue;
				}

				Vec3 orgPos = m_iceSM[jpIndx]->GetOrgPos(k);
				addPIndxList[i].push_back(pair<int, Vec3>(kpIndx, orgPos));
			}

			//�N���X�^�Ɋ܂߂��闱�q���͍ő�300��
			if(addPIndxList[i].size() > 299){
				//cout << "search max" << endl;
				break ;
			}
		}

		////�f�o�b�O
		//std::sort(addPIndxList[i].begin(), addPIndxList[i].end());

		//ofs << i << ":";
		//for(vector<int>::iterator it = addPIndxList[i].begin(); it != addPIndxList[i].end(); it++){
		//	ofs << " " << *it ;
		//}
		//ofs << endl;

		//cout << i << ":";
		//for(vector<int>::iterator it = addPIndxList[i].begin(); it != addPIndxList[i].end(); it++){
		//	cout << " " << *it ;
		//}
		//cout << endl;
	}

double end1 = counter1.End();

counter2.Start();
	//�N���X�^�ɗ��q��ǉ�
	const float* sldPos = Ice_SM::GetSldPosPointer();
	const float* sldVel = Ice_SM::GetSldVelPointer();

	for(unsigned i = 0; i < m_iceSM.size(); i++){
		for(vector<pair<int, Vec3>>::iterator it = addPIndxList[i].begin(); it != addPIndxList[i].end(); it++){
			int addpIndx = it->first;

			//�N���X�^�ɗ��q�����݂��邩���m�F
			if(m_iceSM[i]->CheckIndx(addpIndx)){
				continue;
			}

			int sldIndx = addpIndx * SM_DIM;
			Vec3 pos(sldPos[sldIndx+0], sldPos[sldIndx+1], sldPos[sldIndx+2]);
			Vec3 vel(sldVel[sldIndx+0], sldVel[sldIndx+1], sldVel[sldIndx+2]);
			Vec3 orgPos = it->second;
			//TODO: addpIndx��orgPos�́Cm_iceSM[addpIndx]�ɕK�����݂���̂ŁDvec3��ۑ�����K�v�͂Ȃ�

			m_iceSM[i]->AddAnotherClusterVertex(orgPos, pos, vel, 1.0, addpIndx, 1.0, 0.0, 0);
		}

		m_iceSM[i]->CalcOrgCm();
	}
double end2 = counter2.End();

	cout << "Time:" << "Search:" << end1 << "," << "Add:" << end2 << endl;
}

void CalcIteration::ExpandCluster_Test2(vector<int>& sfFnIndxes)
{
	//�������Ԃ̌v��
//QueryCounter counter1;
//QueryCounter counter2;

//counter1.Start();
	//�ߖT�ɑ��݂��闱�q�Ƃ��̏����ʒu�̃��X�g���쐬
	vector<vector<int>> addPIndxList;
	addPIndxList.resize(m_iceSM.size(), vector<int>());

	for(unsigned i = 0; i < m_iceSM.size(); ++i){
		ResetFlag();
		int nowPrtSize = m_iceSM[i]->GetNumVertices();

		//���X�g�����݂̃N���X�^�Ɋ܂܂�Ă��闱�q�ŏ�����
		for(unsigned j = 0; j < m_iceSM[i]->GetIndxNum(); j++){
			int orgPrt = m_iceSM[i]->GetParticleIndx(j);

			//�ǉ����q���X�g�Ɋ܂܂�Ă��Ȃ���
			if(orgPrt == MAXINT || m_isAdd[orgPrt]){	
				continue;
			}
			
			m_isAdd[orgPrt].flip();
		}

		//�ߖT�N���X�^��T��
		//���ɒT���ς݂̃N���X�^�͏���
		unsigned searchStart = sfFnIndxes[i];
		for(unsigned clsI = searchStart; clsI < m_iceSM[clsI]->GetIndxNum(); clsI++){
			int nearCluster = m_iceSM[i]->GetParticleIndx(clsI);
			if(nearCluster == MAXINT){	
				continue;
			}

			//jpIndx�̃N���X�^�Ɋ܂܂�闱�q���擾
			for(unsigned prtI = 0; prtI < m_iceSM[nearCluster]->GetIndxNum(); prtI++){
				int nrPrt = m_iceSM[nearCluster]->GetParticleIndx(prtI);
				//�ǉ����q���X�g�Ɋ܂܂�Ă��Ȃ���
				if(nrPrt == MAXINT || m_isAdd[nrPrt]){	
					continue;
				}

				//�N���X�^�Ɋ܂߂��闱�q���͍ő�300��
				if(addPIndxList[i].size() >= MAXPARTICLE-nowPrtSize){
					break ;
				}

				m_isAdd[nrPrt].flip();
				addPIndxList[i].push_back(nrPrt);
			}

			//�N���X�^�Ɋ܂߂��闱�q���͍ő�300��
			if(addPIndxList[i].size() >= MAXPARTICLE-nowPrtSize){
				break ;
			}
		}
	}

//double end1 = counter1.End();

//counter2.Start();
	//�N���X�^�ɗ��q��ǉ�
	const float* sldPos = Ice_SM::GetSldPosPointer();
	const float* sldVel = Ice_SM::GetSldVelPointer();

	for(unsigned i = 0; i < m_iceSM.size(); i++){
		for(vector<int>::iterator it = addPIndxList[i].begin(); it != addPIndxList[i].end(); it++){
			int addpIndx = *it;

			if(m_iceSM[i]->GetNumVertices() >= MAXPARTICLE){
				break;
			}

			//�N���X�^�ɗ��q�����݂��邩���m�F
			if(m_iceSM[i]->CheckIndx(addpIndx)){
				continue;
			}

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
		sfFnIndxes[i] += m_iceSM[i]->GetNumVertices();
	}
//double end2 = counter2.End();

	//cout << "Time:" << "Search:" << end1 << "," << "Add:" << end2 << endl;
}

//�N���X�^�̉e���͈͂��k�����čč\��
void CalcIteration::ContractCluster()
{
	//�e�N���X�^�̉e���͈͂��L���邽�߂ɋߖT�ɂ���N���X�^�̗��q��ǉ�����
	for(int i = 0; i < IceObject::GetParticleNum(); ++i){

	}
}

//�f�o�b�O
void CalcIteration::StepObjMoveDebug()
{
	//����
	m_iceMove->StepObjMoveDebug();

	//����
	for(int itr = 1; itr < Ice_SM::GetIteration(); itr++)
	{
		m_iceConvo->StepConvolution();
		m_iceMove->StepObjMoveItrDebug();
	}

	//���x�Z�o
	#pragma omp parallel for
	for(int i = 0; i < IceObject::GetParticleNum(); ++i)
	{	
		if(m_iceMove->GetJudgeMove()->JudgeMoveDebug(i) == false){	continue;	}
		m_iceSM[i]->integrateIteration();
	}
}