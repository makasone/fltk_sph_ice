#include "Ice_CalcMethod_Itr_Expand.h"

typedef Ice_CalcMethod_Itr_Expand CalcIteration;

#define PNUM 2197
#define ADDLIMIT 50
#define SEARCH_LAYER 5

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
	//StepObjMoveTest();

	//����̉^���v�Z
	m_iceMove->StepObjMove();

	//�ő̂ɂ�����ŏI�ʒu����
	m_iceConvo->StepConvolution();

	//�e�N���X�^�̃I���W�i�����q���o���Ă���
	//�i���́C�P�Ȃ�g�������̏����Ȃ�I���W�i�����q�̌������ł悢�j
	vector<vector<unsigned>> copyOriginalParticleIndxes(m_iceSM.size(), vector<unsigned>());
	CopyOriginalObject(copyOriginalParticleIndxes);

	//�e�������ɑ����ς݂̓Y���ԍ����L�^���Ă���
	vector<int> searchFinishIndxes(m_iceSM.size(), 0);

	////���񂾂��F�����̃N���X�^�ŉe���͈͂��g��
	//ExpandCluster_Far();

	//����
	for(int itr = 1; itr < Ice_SM::GetIteration(); itr++){	
		////�N���X�^�̉e���͈͂��g�債�čč\��
		//ExpandCluster(searchFinishIndxes);

		//�����̃N���X�^�ŉe���͈͂��g��
		ExpandCluster_Far();

		//�������̉^���v�Z
		m_iceMove->StepObjMoveItr();

		//�ő̂ɂ�����ŏI�ʒu����
		m_iceConvo->StepConvolution();
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
		m_iceSM[i]->ClassifyAllOrgParticle();
	}
}

//�č\�z����N���X�^�����o
void CalcIteration::GetExpandeCluster()
{
	for(unsigned i = 0; i < m_iceSM.size(); ++i){

	}
}

//�N���X�^�̉e���͈͂��g�債�čč\��
void CalcIteration::ExpandCluster(vector<int>& searchFinishIndxes)
{
		//�������Ԃ̌v��
QueryCounter counter1;
QueryCounter counter2;

counter1.Start();
	//�ߖT�ɑ��݂��闱�q�Ƃ��̏����ʒu�̃��X�g���쐬
	vector<vector<int>> addPIndxList;
	addPIndxList.resize(m_iceSM.size(), vector<int>(ADDLIMIT, MAXINT));

	SelectAddParticleFromNearestCluster(addPIndxList, searchFinishIndxes);
double end1 = counter1.End();

counter2.Start();
	AddParticleToCluster(addPIndxList);
double end2 = counter2.End();

	cout << "Time:" << "Search:" << end1 << "," << "Add:" << end2 << endl;
}

void CalcIteration::ExpandCluster_Far()
{
	//�������Ԃ̌v��
QueryCounter counter1;
QueryCounter counter2;

counter1.Start();
	//�ߖT�ɑ��݂��闱�q�Ƃ��̏����ʒu�̃��X�g���쐬
	vector<vector<int>> addPIndxList;
	addPIndxList.resize(m_iceSM.size(), vector<int>(ADDLIMIT, MAXINT));

	SelectAddParticleFromFarCluster(addPIndxList);
double end1 = counter1.End();

counter2.Start();
	AddParticleToCluster(addPIndxList);
double end2 = counter2.End();

	cout << "Time:" << "Search:" << end1 << "," << "Add:" << end2 << endl;
}

//�ߖT�T���F�N���X�^�֒ǉ����闱�q�̓Y�������擾
void CalcIteration::SelectAddParticleFromNearestCluster(vector<vector<int>>& addPIndxList, vector<int>& sfFnIndxes)
{
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

		//����N���X�^�Ɋ܂܂�Ă��闱�q���Ώ�
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
}

//�K�w�T���F�N���X�^�֒ǉ����闱�q�̓Y�������擾
void CalcIteration::SelectAddParticleFromFarCluster(vector<vector<int>>& addPIndxList)
{
	//string result = RESULT_DATA_PATH;
	//result += "SelectAddParticleFromFarCluster.txt";
	//ofstream ofs(result);

	int debug_searchCounter = 0;
	int debug_indxSum = 0;
	int debug_crossHold = 0;
	int debug_secondAdd = 0;

	for(unsigned i = 0; i < m_iceSM.size(); ++i){
		if(m_iceSM[i]->GetNumVertices() >= MAXPARTICLE){
			continue;
		}

		//ofs << "Cluster:" << i << endl;

		//���q��ǉ����Ă��邩�ǂ����̃t���O
		bool isAdd[PNUM] = {};

		//���݂̃N���X�^�Ɋ܂܂�Ă��闱�q�͒ǉ����Ȃ�
		for(unsigned j = 0; j < m_iceSM[i]->GetIndxNum(); j++){
			//���q�ԍ��ł���
			int orgPrt = m_iceSM[i]->GetParticleIndx(j);
			if(orgPrt == MAXINT){	
				continue;
			}

			////������O�̃N���X�^�̒T���Œǉ�����Ă���ꍇ
			//if(addPIndxList[i][j] != MAXINT){
			//	int preAddPrt = addPIndxList[i][j];
			//	isAdd[preAddPrt] = true;
			//	continue;
			//}
			
			isAdd[orgPrt] = true;
		}

		//�w�����T��
		////�ǉ����闱�q���F�����̃N���X�^�Ɋ܂܂�Ă��闱�q��
		//unsigned addParticleNum = m_iceSM[i]->GetNumVertices();

//QueryCounter counter;
//counter.Start();
		//�N���X�^�Ɋ܂܂�Ă��闱�q��T���̋N�_�Ƃ���
		for(unsigned clsI = 1; clsI < m_iceSM[i]->GetIndxNum(); clsI++){
			//debug_indxSum++;

			////���ɒǉ��ς݂��ǂ���
			//if(addPIndxList[i][clsI] != MAXINT){
			//	//cout << "CONTINU+E" << endl;
			//	debug_crossHold++;
			//	continue;
			//}

			//�T���N�_���q
			int nearCluster = m_iceSM[i]->GetParticleIndx(clsI);
			
			//���q�ԍ��ł���
			if(nearCluster == MAXINT){
				continue;
			}

			//ofs << "	pIndx:" << nearCluster << endl;

			//�N�_���q�̕����������P�ʃx�N�g��
			const Vec3 dirVec = Unit(m_iceSM[i]->GetOrgPos(clsI) - m_iceSM[i]->GetOrgCm());

			//�N�_���q�̗̈�
			const Ice_SM::EreaData startErea = m_iceSM[i]->erea(clsI);

			//���������x�N�g���������q��H��
			pair<int, int> resultIndxes = SearchSimilarParticle(nearCluster, dirVec, startErea);

			int searchPrtIndx = resultIndxes.first;
			int searchClsIndx = resultIndxes.second;

			//�T�����q��ǉ�
			if(isAdd[searchPrtIndx]){
				//���Ɋ܂܂�Ă���Ȃ�C�ǉ����Ȃ�
				//TODO:�����Ŗ߂�Ƃ������Ƃ́C���ʂȒT��������Ƃ������ƁD�Ȃ�ׂ����ʂɂȂ�Ȃ��悤�ɗ��q��ǉ�����
				continue;;
			}

			debug_searchCounter++;

			//�Ō�ɑI�񂾗ގ����q���N���X�^�ɒǉ�����
			isAdd[searchPrtIndx] = true;
			addPIndxList[i][clsI] = searchPrtIndx;

			////ofs << "		addPIndx:" << searchPrtIndx
			////	<< "	nowIndxNum:" << clsI << endl;

			////TODO:���������֌W
			////�T�����i�܂��ɑł��؂�ꂽ�ꍇ�͎��������֌W�ɂ͂Ȃ�Ȃ�
			//if(searchClsIndx == MAXINT){
			//	continue;
			//}

			////�ϐ��̊֌W���t�ł�₱�����̂ɒ���
			//int oIndx = m_iceSM[searchPrtIndx]->SearchIndx(searchClsIndx);
			////if(oIndx == MAXINT){
			////	//ofs << __FUNCTION__ << " ERROR" << endl;
			////	continue;
			////}

			////���������֌W��ۑ����������C1�Α��֌W�ɂȂ��Ă���
			//if(addPIndxList[searchPrtIndx][oIndx] != MAXINT){
			//	//ofs << "1�Α�" << "searchPrtIndx:" << searchPrtIndx << "," << "searchClsIndx:" << searchClsIndx << endl;
			//	continue;
			//}

			//addPIndxList[searchPrtIndx][oIndx] = i;
			//debug_secondAdd++;
		}//for(unsigned clsI = 0;
//double end = counter.End();
//cout << "LOOP_TIME:" << end << endl;
	}//for(unsigned i = 0;

	cout << "HOLD:" << debug_crossHold << ",";
	cout << "SEARCH_COUNTER:" << debug_searchCounter << ",";
	cout << "SEDONDADD:" << debug_secondAdd << ",";
	cout << endl;
}

//���������x�N�g���������q��H���āC�������q���擾
pair<int, int> CalcIteration::SearchSimilarParticle(int nearCluster, const Vec3& dirVec, const Ice_SM::EreaData& startErea)
{
	//�T���K�w�F�K���Ɏw��
	unsigned searchLayer = SEARCH_LAYER;

	//�T�����q�̎��N���X�^���擾���C�܂܂�Ă��闱�q�𒲍�
	unsigned nowLayer = 1;
	int mostSimPrtIndx = MAXINT;
	double mostSimilarity = 0.0;

	for(unsigned srchI = 0; srchI <= m_iceSM[nearCluster]->GetIndxNum(); srchI++){
		//�I������
		if(srchI == m_iceSM[nearCluster]->GetIndxNum()){
			//�T���𑱂��邩�̔���
			 if(mostSimPrtIndx == MAXINT){
				//�����₪�����炸�T�����i�܂Ȃ������ꍇ�C���݂̃N���X�^��T�����ʂƂ��ĒT���I��
				return pair<int, int>(nearCluster, MAXINT);
			}
			else if(nowLayer == searchLayer){
				//�T���̉񐔐������}�����ꍇ�C�T���I��
				return pair<int, int>(mostSimPrtIndx, nearCluster);
			}

			//�T���𑱂���Ȃ�C�e�ϐ����X�V
			nowLayer++;
			
			srchI = 0;
			nearCluster = mostSimPrtIndx;
			
			mostSimPrtIndx = MAXINT;
			mostSimilarity = 0.0;

			continue;
		}

		int searchPrt = m_iceSM[nearCluster]->GetParticleIndx(srchI);

		//���q�ԍ��ł���C�������g�łȂ��C���O�Ɣ��f����Ă��Ȃ�
		if(searchPrt == MAXINT || searchPrt == nearCluster){
			continue;
		}

		//�T���̑Ώۗ̈�ł��邩�𔻒�
		const Ice_SM::EreaData searchErea = m_iceSM[nearCluster]->erea(srchI);

		//�̈�Ԃ̋����Ŕ���
		if(Ice_SM::CalcEreaDistance(startErea, searchErea) >= 1){
			continue;
		}

		//�T�����q�̕����������P�ʃx�N�g��
		Vec3 srchDirVec = Unit(m_iceSM[nearCluster]->GetOrgPos(srchI) - m_iceSM[nearCluster]->GetOrgCm());

		//�e���q��cos�ގ��x���ׂ�
		double similarity = dot(dirVec, srchDirVec);

		//�ގ��x���������q��I��
		if(mostSimPrtIndx == MAXINT || mostSimilarity <= similarity){
			mostSimPrtIndx = searchPrt;
			mostSimilarity = similarity;
		}
	}

	return pair<int, int>(nearCluster, MAXINT);
}

//���X�g����ɃN���X�^�֗��q��ǉ�
void CalcIteration::AddParticleToCluster(const vector<vector<int>>& addPIndxList)
{
	//�f�o�b�O
	int debug_addNum = 0;

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

		for(vector<int>::const_iterator it = addPIndxList[i].begin(); it != addPIndxList[i].end(); it++){
			int addpIndx = *it;

			if(m_iceSM[i]->GetNumVertices() >= MAXPARTICLE){
				break;
			}

			//far�̏ꍇ��continue, normal��break
			if(addpIndx == MAXINT){
				//break;
				continue;
			}

			//�N���X�^�ɗ��q�����݂��邩���m�F
			if(isAdd[addpIndx]){
				continue;
			}

			isAdd[addpIndx] = true;

			//addpIndx�Ԗڂ̂O�Ԗڂ̗��q�́C���q�ԍ�addpIndx�ƕۏ؂���Ă���
			int sldIndx = addpIndx * SM_DIM;
			Vec3 orgPos(m_iceSM[addpIndx]->GetOrgPos(0));
			Vec3 pos(sldPos[sldIndx+0], sldPos[sldIndx+1], sldPos[sldIndx+2]);
			Vec3 vel(sldVel[sldIndx+0], sldVel[sldIndx+1], sldVel[sldIndx+2]);

			debug_addNum++;

			int check = m_iceSM[i]->AddAnotherClusterVertex(orgPos, pos, vel, 1.0, addpIndx, 1.0, 0.0, 0);
			if(check == -1){
				cout << "over " << addPIndxList[i].size() << endl;
			}
		}

		m_iceSM[i]->CalcOrgCm();
		m_iceSM[i]->ClassifyAllOrgParticle();
	}

	cout << __FUNCTION__ << " addNum:" << debug_addNum << endl;
}

//�N���X�^�̉e���͈͂��k�����čč\��
void CalcIteration::ContractCluster()
{
	//�e�N���X�^�̉e���͈͂��L���邽�߂ɋߖT�ɂ���N���X�^�̗��q��ǉ�����
	for(int i = 0; i < IceObject::GetParticleNum(); ++i){

	}
}

//�f�o�b�O-----------------------------------------------------------------
void CalcIteration::StepObjMoveDebug()
{
QueryCounter counter1;
QueryCounter counter2;
QueryCounter counter3;
QueryCounter counter4;

counter1.Start();
	//����̉^���v�Z
	m_iceMove->StepObjMove();

	//�ő̂ɂ�����ŏI�ʒu����
	m_iceConvo->StepConvolution();
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
		//�N���X�^�̉e���͈͂��g�債�čč\��
		///ExpandCluster(searchFinishIndxes);
		ExpandCluster_Far();

		//�������̉^���v�Z
		m_iceMove->StepObjMoveItr();

		//�ő̂ɂ�����ŏI�ʒu����
		m_iceConvo->StepConvolution();
	}
double end3 = counter3.End();

counter4.Start();
	//�R�s�[�𗘗p���ăN���X�^�̍\�������ɖ߂�
	ReplaceCluster(copyOriginalParticleIndxes);
double end4 = counter4.End();

	//���x�Z�o
	CalcVel();

	//�v�Z���Ԃ̔�r
	cout << "Time:"
		<< " first:" << end1
		<< " prepare:" << end2
		<< " itr:" << end3
		<< " replace:" << end4
		<< endl;
}

//�����N���X�^����C�ǉ����邽�߂̗��q��I��
void CalcIteration::Debug_SelectAddParticleFromFarCluster(vector<vector<int>>& addPIndxList)
{
	string result = RESULT_DATA_PATH;
	result += "SelectAddParticleFromFarCluster.txt";
	ofstream ofs(result, ios::app);

	for(unsigned i = 0; i < m_iceSM.size(); ++i){
		if(m_iceSM[i]->GetNumVertices() >= MAXPARTICLE){
			continue;
		}

		ofs << "Cluster:" << i << endl;

		//���q��ǉ����Ă��邩�ǂ����̃t���O
		bool isAdd[PNUM] = {};

		//���݂̃N���X�^�Ɋ܂܂�Ă��闱�q�͒ǉ����Ȃ�
		for(unsigned j = 0; j < m_iceSM[i]->GetIndxNum(); j++){
			int orgPrt = m_iceSM[i]->GetParticleIndx(j);

			//�ǉ����q���X�g�Ɋ܂܂�Ă��Ȃ���
			if(orgPrt == MAXINT || isAdd[orgPrt]){	
				continue;
			}
			
			isAdd[orgPrt] = true;
		}

		//�w�����T��
		////�ǉ����闱�q���F�����̃N���X�^�Ɋ܂܂�Ă��闱�q��
		//unsigned addParticleNum = m_iceSM[i]->GetNumVertices();

		//�T���K�w�F�K���Ɏw��
		unsigned searchLayer = SEARCH_LAYER;

		//���݂̗��q�����Y����
		int addIndxNum = 0;

		//�N���X�^�Ɋ܂܂�Ă��闱�q��T���̋N�_�Ƃ���
		for(unsigned clsI = 0; clsI < m_iceSM[i]->GetIndxNum(); clsI++){
			int nearCluster = m_iceSM[i]->GetParticleIndx(clsI);
			if(nearCluster == MAXINT || i == nearCluster){	
				continue;
			}

			//�N�_���q�̕����������P�ʃx�N�g��
			Vec3 dirVec = Unit(m_iceSM[i]->GetOrgPos(clsI) - m_iceSM[i]->GetOrgPos(0));

			//ofs << "	StartCluster:" << nearCluster
			//	<< "	dirVec:" << dirVec
			//	<< endl;

			//�T�����q�̎��N���X�^���擾���C�܂܂�Ă��闱�q�𒲍�
			unsigned nowLayer = 1;
			int mostSimPrtIndx = MAXINT;
			double mostSimilarity = 0.0;

			for(unsigned srchI = 0; srchI < m_iceSM[nearCluster]->GetIndxNum(); srchI++){
				int searchPrt = m_iceSM[nearCluster]->GetParticleIndx(srchI);
				if(searchPrt == MAXINT || searchPrt == nearCluster){	
					continue;
				}

				//�T�����q�̕����������P�ʃx�N�g��
				Vec3 srchDirVec = Unit(m_iceSM[nearCluster]->GetOrgPos(srchI) - m_iceSM[nearCluster]->GetOrgPos(0));

				//�e���q��cos�ގ��x���ׂ�
				double similarity = dot(dirVec, srchDirVec);

				//ofs << "		srchI:" << srchI
				//	<< "	searchPrt:" << searchPrt
				//	<< "	similarity:" << similarity
				//	<< "	srchDirVec:" << srchDirVec << endl;

				//�ގ��x���������q��I��
				if(mostSimPrtIndx == MAXINT || mostSimilarity <= similarity){
					mostSimPrtIndx = searchPrt;
					mostSimilarity = similarity;
				}

				//�I������
				if(srchI < m_iceSM[nearCluster]->GetIndxNum() -1){
					continue;
				}

				//ofs << "	End"
				//	<< "	nowLayer:" << nowLayer
				//	<< "	nearClst:" << nearCluster
				//	<< "	mostSimPrtIndx" << mostSimPrtIndx
				//	<< "	mostSimilarity" << mostSimilarity << endl;

				//�T���𑱂��邩�̕���
				if(nowLayer < searchLayer){
					//�T���𑱂���Ȃ�C�e�ϐ����X�V
					nowLayer++;

					srchI = 0;
					nearCluster = mostSimPrtIndx;

					mostSimPrtIndx = MAXINT;
					mostSimilarity = 0.0;
					continue;
				}
				else{
					ofs << "		addPIndx:" << mostSimPrtIndx
						<< "	nowIndxNum:" << addIndxNum << endl;

					//�T���̉񐔐������}����ƁC�T���I��
					//�Ō�ɑI�񂾗ގ����q���N���X�^�ɒǉ�����
					//TODO: ���̂܂܂��ƒǉ�����Ȃ��ꍇ�������Ȃ�̂ŁCstack���g��
					if(!isAdd[mostSimPrtIndx]){
						isAdd[mostSimPrtIndx] = true;
						addPIndxList[i][addIndxNum++] = mostSimPrtIndx;
					}
					break;
				}
			}//for(unsigned srchI = 0;
		}//for(unsigned clsI = 0;
	}//for(unsigned i = 0;
}