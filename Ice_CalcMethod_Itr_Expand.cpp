#include "Ice_CalcMethod_Itr_Expand.h"

typedef Ice_CalcMethod_Itr_Expand CalcIteration;

#define PNUM 1331

CalcIteration::Ice_CalcMethod_Itr_Expand(const vector<Ice_SM*>& iceSM, Ice_ClusterMove* clusterMove, Ice_Convolution* convo)
{
	m_iceSM = iceSM;
	SetObjMove(clusterMove);
	SetConvolution(convo);
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
	//����
	m_iceMove->StepObjMove();

	//���̃I�u�W�F�N�g���R�s�[
	Ice_SM copyObj[PNUM];

	//������Z�q�ŃR�s�[
	for(int i = 0; i < IceObject::GetParticleNum(); ++i){
		copyObj[i] = *m_iceSM[i];
	}

	//����
	for(int itr = 1; itr < Ice_SM::GetIteration(); itr++)
	{	
		m_iceConvo->StepConvolution();
	
		//�N���X�^�̉e���͈͂��g�債�čč\��
		ExpandeCluster();

		m_iceMove->StepObjMoveItr();

		//TODO: �ό`�ʂ̑傫���N���X�^�����o
		//GetExpandeCluster(expandObj);
	}

	//ContractCluster(copyObj);

	//�R�s�[�̍폜
	//DeleteCopyObject(copyObj);

	//�R�s�[�𗘗p���ăN���X�^�̍\�������ɖ߂�
	ReplaceCluster(copyObj);

	//���x�Z�o
	#pragma omp parallel for
	for(int i = 0; i < IceObject::GetParticleNum(); ++i)
	{	
		if(m_iceMove->GetJudgeMove()->JudgeMove(i) == false){	continue;	}
		m_iceSM[i]->integrateIteration();
	}
}

//���̃I�u�W�F�N�g���R�s�[
void CalcIteration::CopyOriginalObject(vector<Ice_SM>& copyObj)
{
	//FUNCTION_NAME;
	
}

void CalcIteration::DeleteCopyObject(vector<Ice_SM>& copyObj)
{
	//FUNCTION_NAME;

}

void CalcIteration::ReplaceCluster(Ice_SM copyObj[PNUM])
{
	//copyObj��m_iPindxes�ɑ��݂��闱�q�̓I���W�i��
	//���݂��Ȃ����q�̏���S�ď���

	//TODO: �K�v�ȏ��́C�����̗��q�ԍ�����������
	for(unsigned i = 0; i < m_iceSM.size(); i++){
		for(unsigned j = 0; j < m_iceSM[i]->GetIndxNum(); j++){
			//���݂�m_iceSM�Ɋ܂܂�闱�q�ŁCcopyObj�ɑ��݂��Ȃ����̂��폜
			int jpIndx = m_iceSM[i]->GetParticleIndx(j);

			if(jpIndx == MAXINT || copyObj[i].CheckIndx(jpIndx)){
				continue;
			}

			//copyObj�ɑ��݂��Ȃ��̂Œǉ����ꂽ���q�@�폜����
			//cout << i << " delete " << jpIndx << endl;
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

//�N���X�^�̉e���͈͂��g�債�čč\��
void CalcIteration::ExpandeCluster()
{
	float* sldPos = Ice_SM::GetSldPosPointer();
	float* sldVel = Ice_SM::GetSldVelPointer();

	////�e�X�g�S�Ă̗��q��remove���Ă݂�
	//for(unsigned i = 0; i < m_iceSM.size(); i++){
	//	for(unsigned j = 0; j < m_iceSM[i]->GetIndxNum(); j++){
	//		int jpIndx = m_iceSM[i]->GetParticleIndx(j);
	//		if(jpIndx == MAXINT){
	//			continue;
	//		}
	//	}
	//}

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
	
	//TODO: ���q��ǉ�
	for(unsigned i = 0; i < m_iceSM.size(); i++){
		for(vector<pair<int, Vec3>>::iterator it = addPIndxList[i].begin(); it != addPIndxList[i].end(); it++){
			int addpIndx = it->first;

			//�N���X�^�ɗ��q�����݂��邩���m�F
			if(m_iceSM[i]->CheckIndx(addpIndx)){
				continue;
			}

			int sldIndx = addpIndx * SM_DIM;
			Vec3 pos = Vec3(sldPos[sldIndx+0], sldPos[sldIndx+1], sldPos[sldIndx+2]);
			Vec3 vel = Vec3(sldVel[sldIndx+0], sldVel[sldIndx+1], sldVel[sldIndx+2]);
			Vec3 orgPos = it->second;

			m_iceSM[i]->AddAnotherClusterVertex(orgPos, pos, vel, 1.0, addpIndx, 1.0, 0.0, 0);
		}

		m_iceSM[i]->CalcOrgCm();

		//�f�o�b�O
	}
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