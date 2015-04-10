#include "Ice_CalcMethod_Itr_Expand.h"

typedef Ice_CalcMethod_Itr_Expand CalcIteration;

#define PNUM 7000
#define ADDLIMIT 50
#define SEARCH_LAYER 5

CalcIteration::Ice_CalcMethod_Itr_Expand(const vector<Ice_SM*>& iceSM, IceStructure* iceStrct, Ice_ClusterMove* clusterMove, Ice_Convolution* convo)
{
	m_iceSM = iceSM;
	m_iceStrct = iceStrct;
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
	//����̉^���v�Z
	m_iceMove->StepObjMove();

	//�ő̂ɂ�����ŏI�ʒu����
	m_iceConvo->StepConvolution();

	//�e�N���X�^�̃I���W�i�����q���o���Ă���
	//�i���́C�P�Ȃ�g�������̏����Ȃ�I���W�i�����q�̌������ł悢�j
	vector<vector<unsigned>> copyOriginalParticleIndxes(m_iceSM.size(), vector<unsigned>());
	CopyOriginalObject(copyOriginalParticleIndxes);

	////�e�������ɑ����ς݂̓Y���ԍ����L�^���Ă���
	//vector<int> searchFinishIndxes(m_iceSM.size(), 0);

	//�g�������N���X�^�̔ԍ�
	vector<int> expandClusterIndx;

	//���񂾂��F�����̃N���X�^�ŉe���͈͂��g��
	ExchangeCluster_Far(expandClusterIndx);
	//ExpandCluster_Far();

	//����
	for(int itr = 1; itr < Ice_SM::GetIteration(); itr++){	
		////�N���X�^�̉e���͈͂��g�債�čč\��
		//ExpandCluster(searchFinishIndxes);

		////�����̃N���X�^�ŉe���͈͂��g��
		//ExpandCluster_Far();

		////�����̃N���X�^�ŉe���͈͂��g��@�������邽�тɒǉ����q���X�V����
		//ExpandCluster_Far_Step();

		////�����̗��q�ƌ��݂̗��q�����ւ���
		//ExchangeCluster_Far();

		//�������̉^���v�Z
		m_iceMove->StepObjMoveItr();

		//�ő̂ɂ�����ŏI�ʒu����
		m_iceConvo->StepConvolution();
	}

	//�R�s�[�𗘗p���ăN���X�^�̍\�������ɖ߂�
	ResetCluster(copyOriginalParticleIndxes, expandClusterIndx);

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
void CalcIteration::ResetCluster(const vector<vector<unsigned>>& copyIndxes, const vector<int>& exClstrIndxes)
{
	////���q��ǉ������ꍇ�̏����F
	////copyIndxes�ɑ��݂��闱�q�̓I���W�i��
	////���݂��Ȃ����q�̏���S�ď���
	//for(unsigned i = 0; i < m_iceSM.size(); i++){
	//	for(unsigned j = 0; j < m_iceSM[i]->GetIndxNum(); j++){
	//		//���݂�m_iceSM�Ɋ܂܂�闱�q�ŁCcopyObj�ɑ��݂��Ȃ����̂��폜
	//		int jpIndx = m_iceSM[i]->GetParticleIndx(j);
	//		if(jpIndx == MAXINT){
	//			continue;
	//		}

	//		vector<unsigned>::const_iterator check = find(copyIndxes[i].begin(), copyIndxes[i].end(), jpIndx);
	//		if(check != copyIndxes[i].end()){
	//			continue;
	//		}

	//		//copyIndxes�ɑ��݂��Ȃ��̂Œǉ����ꂽ���q�@�폜����
	//		int ooIndx = m_iceSM[i]->SearchIndx(jpIndx);
	//		m_iceSM[i]->Remove(ooIndx);
	//	}

	//	m_iceSM[i]->CalcOrgCm();
	//	m_iceSM[i]->ClassifyAllOrgParticle();
	//}

QueryCounter count1;
QueryCounter count2;

count1.Start();
	//���q�����������ꍇ�̏���
	const float* sldPos = Ice_SM::GetSldPosPointer();
	const float* sldVel = Ice_SM::GetSldVelPointer();

	//�N���X�^�g���������̂������Z�b�g����
	//���q�����ւ����ꍇ�̏����@�����͏d�����ǈ��S
	//�N���X�^�ԍ������q�ԍ��ƂȂ�ŏ��̈��
	//for(unsigned pIndx = 0; pIndx < m_iceSM.size(); pIndx++){
	for(vector<int>::const_iterator it = exClstrIndxes.begin(); it != exClstrIndxes.end(); it++){
		//�ŏ��̈�����͕ۑ����Ă���
		int pIndx = *it;
		int sldIndx = pIndx * SM_DIM;
		Vec3 orgPos(m_iceSM[pIndx]->GetOrgPos(0));
		Vec3 pos(sldPos[sldIndx+0], sldPos[sldIndx+1], sldPos[sldIndx+2]);
		Vec3 vel(sldVel[sldIndx+0], sldVel[sldIndx+1], sldVel[sldIndx+2]);

		//������
		m_iceSM[pIndx]->Clear();

		//�N���X�^�ɗ��q��ǉ�������
		int check = m_iceSM[pIndx]->AddAnotherClusterVertex(orgPos, pos, vel, 1.0, pIndx, 1.0, 0.0, 0);
	}
double end1 = count1.End();

count2.Start();
	//���̑��̗��q���l�ߒ���
	//for(unsigned pIndx = 0; pIndx < m_iceSM.size(); pIndx++){
	for(vector<int>::const_iterator it = exClstrIndxes.begin(); it != exClstrIndxes.end(); it++){
		int pIndx = *it;
		for(unsigned j = 0; j < copyIndxes[pIndx].size(); j++){

			int copyIndx = copyIndxes[pIndx][j];
			if(pIndx == copyIndx){
				continue;
			}

			Vec3 orgPos(m_iceSM[copyIndx]->GetOrgPos(0));
			Vec3 pos(m_iceSM[copyIndx]->GetVertexPos(0));
			Vec3 vel(m_iceSM[copyIndx]->GetVertexVel(0));

			int check = m_iceSM[pIndx]->AddAnotherClusterVertex(orgPos, pos, vel, 1.0, copyIndx, 1.0, 0.0, 0);
		}

		m_iceSM[pIndx]->CalcOrgCm();
		m_iceSM[pIndx]->ClassifyAllOrgParticle();
	}
double end2 = count2.End();

//	//cout << __FUNCTION__ << "," << "end1:" << end1 << "," << "end2:" << end2 << endl;
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

//�����̗��q���g���ăN���X�^�̉e���͈͂��g��
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
	//���q���N���X�^�ɒǉ�
	AddParticleToCluster(addPIndxList);
double end2 = counter2.End();

	cout << "Time:" << "Search:" << end1 << "," << "Add:" << end2 << endl;
}

//�����̗��q���g���ăN���X�^�̗v�f���������C�e���͈͂��g�債�Ă���
void CalcIteration::ExchangeCluster_Far(vector<int>& exClstrIndxes)
{
//�������Ԃ̌v��
QueryCounter counter1;
QueryCounter counter2;
QueryCounter counter3;

counter1.Start();
	//�ߖT�ɑ��݂��闱�q�Ƃ��̏����ʒu�̃��X�g���쐬
	vector<vector<int>> exchangePIndxList;
	exchangePIndxList.resize(m_iceSM.size(), vector<int>(ADDLIMIT, MAXINT));

	SelectExchangeParticleFromFarCluster(exchangePIndxList);
double end1 = counter1.End();

counter2.Start();
	//�T�����ʂ̗��q�ƃN���X�^�Ɋ܂܂�闱�q�����ւ���
	ExchangeParticleToCluster(exchangePIndxList);
double end2 = counter2.End();

counter3. Start();
	//�������s��ꂽ�N���X�^�ԍ����m�F
	ConfirmExClusterIndx(exchangePIndxList, exClstrIndxes);
double end3 = counter3.End();

	//cout << "Time:" << "Search:" << end1 << "," << "Exchange:" << end2 << "," << "Confirm" << end3 << endl;
}

//�����̗��q���g���ăN���X�^�̉e���͈͂��g��@�������邽�тɊg�嗱�q���X�V����
void CalcIteration::ExpandCluster_Far_Step()
{

}

//���q�����������N���X�^�ԍ����m�F
void CalcIteration::ConfirmExClusterIndx(const vector<vector<int>> exIndxList, vector<int>& exClstrIndxes)
{
	for(unsigned i = 0; i < exIndxList.size(); i++){

		////TODO: ���܂������Ȃ��c
		////exIndxList��MAXINT��菬�����l�����݂��邩�ǂ���
		////���݂���Ȃ�������s���Ă���C���݂��Ȃ��Ȃ�������s���Ă��Ȃ��̂Ŗ߂�
		//vector<int> exList = exIndxList[i];
		//if( exList.end() != find(exList.begin(), exList.end(), bind2nd(less<int>(), MAXINT)) ){
		//	continue;
		//}

		//�ǂ��l���Ă��璷
		bool check = false;
		for(vector<int>::const_iterator it = exIndxList[i].begin(); it != exIndxList[i].end(); it++){
			if(*it != MAXINT){	check = true;	break;	}
		}
		if(!check){	continue;	}

		exClstrIndxes.push_back(i);
	}

	//���x��Exchange�����ꍇ�ɔ���������
	////�\�[�g
	//sort(exClstrIndxes.begin(), exClstrIndxes.end());
	//
	////�d���폜
	//exClstrIndxes.erase( 
	//	unique(exClstrIndxes.begin(), exClstrIndxes.end()),
	//	exClstrIndxes.end()
	//);

	//�T�C�Y����
	exClstrIndxes.shrink_to_fit();

	//cout << __FUNCTION__ << "," << "exClstrIndxes.size:" << exClstrIndxes.size() << endl;
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
			if(addIndxNum >= limitPrtNum || addIndxNum >= ADDLIMIT){
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

int debug_searchAll = 0;
int debug_searchCompare = 0;

//�K�w�T���F�N���X�^�֒ǉ����闱�q�̓Y�������擾
void CalcIteration::SelectAddParticleFromFarCluster(vector<vector<int>>& addPIndxList)
{
	debug_searchAll = 0;
	debug_searchCompare = 0;

	//�a�ȃN���X�^�݂̂�I�Ԏ���
	m_iceStrct->InitSelectClusterFromClusterSet(m_iceSM);

	for(unsigned i = 0; i < m_iceSM.size(); ++i){
		if(m_iceSM[i]->GetNumVertices() >= MAXPARTICLE){
			continue;
		}

		//�����F�a�ȃN���X�^�̂ݒT������
		if(m_iceStrct->GetMotionCalcCluster(i) == 0){
			continue;
		}

		//���q��ǉ����Ă��邩�ǂ����̃t���O
		bool isAdd[PNUM] = {};

		//���݂̃N���X�^�Ɋ܂܂�Ă��闱�q�͒ǉ����Ȃ�
		for(unsigned j = 0; j < m_iceSM[i]->GetIndxNum(); j++){
			//���q�ԍ��ł���
			int orgPrt = m_iceSM[i]->GetParticleIndx(j);
			if(orgPrt == MAXINT){	
				continue;
			}

			isAdd[orgPrt] = true;
		}

		//�w�����T��
		////�ǉ����闱�q���F�����̃N���X�^�Ɋ܂܂�Ă��闱�q��
		//unsigned addParticleNum = m_iceSM[i]->GetNumVertices();
	
		//�N���X�^�Ɋ܂܂�Ă��闱�q��T���̋N�_�Ƃ��� ���q�ԍ����N���X�^�ԍ��͏����̂�i=1����X�^�[�g
		for(unsigned clsI = 1; clsI < m_iceSM[i]->GetIndxNum(); clsI++){

			//�T���N�_���q
			int nearCluster = m_iceSM[i]->GetParticleIndx(clsI);
			
			//���q�ԍ��ł���
			if(nearCluster == MAXINT){
				continue;
			}

			//�N�_���q�̕����������P�ʃx�N�g��
			const Vec3 dirVec = Unit(m_iceSM[i]->GetOrgPos(clsI) - m_iceSM[i]->GetOrgCm());

			//�N�_���q�̗̈�
			const Ice_SM::EreaData startErea = m_iceSM[i]->erea(clsI);

			//���������x�N�g���������q��H��
			const pair<int, int> resultIndxes = SearchSimilarParticle(nearCluster, dirVec, startErea);

			int searchPrtIndx = resultIndxes.first;
			int searchClsIndx = resultIndxes.second;

			//�T�����q��ǉ�
			if(isAdd[searchPrtIndx]){
				//���Ɋ܂܂�Ă���Ȃ�C�ǉ����Ȃ�
				//TODO:�����Ŗ߂�Ƃ������Ƃ́C���ʂȒT��������Ƃ������ƁD�Ȃ�ׂ����ʂɂȂ�Ȃ��悤�ɗ��q��ǉ�����
				continue;;
			}

			//�Ō�ɑI�񂾗ގ����q���N���X�^�ɒǉ�����
			isAdd[searchPrtIndx] = true;
			addPIndxList[i][clsI] = searchPrtIndx;
		}//for(unsigned clsI = 0;
	}//for(unsigned i = 0;

	//����
	m_iceStrct->ResetSelectCluster(m_iceSM);

	//cout << "All:" << debug_searchAll << "," << "Compare:" << debug_searchCompare << endl;
}

//���������x�N�g���������q��H���āC�������q���擾
const pair<int, int> CalcIteration::SearchSimilarParticle(int nearCluster, const Vec3& dirVec, const Ice_SM::EreaData& startErea)
{
	//�T���K�w�F�K���Ɏw��
	unsigned searchLayer = SEARCH_LAYER;

	//�T�����q�̎��N���X�^���擾���C�܂܂�Ă��闱�q�𒲍�
	unsigned nowLayer = 1;
	unsigned IndxNum = m_iceSM[nearCluster]->GetIndxNum();
	int mostSimPrtIndx = MAXINT;
	double mostSimilarity = 0.0;

	for(unsigned srchI = 0; srchI <= IndxNum; srchI++){
		debug_searchAll++;

		//�I������
		if(srchI == IndxNum){
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

			IndxNum = m_iceSM[nearCluster]->GetIndxNum();
			
			continue;
		}

		int searchPrt = m_iceSM[nearCluster]->GetParticleIndx(srchI);

		//���q�ԍ��ł���C�������g�łȂ�
		if(m_iceSM[nearCluster]->CheckHole(srchI) || searchPrt == nearCluster){
			continue;
		}

		//�T���̑Ώۗ̈�ł��邩�𔻒�
		const Ice_SM::EreaData searchErea = m_iceSM[nearCluster]->erea(srchI);

		//�̈�Ԃ̋����Ŕ���
		if(Ice_SM::CalcEreaDistance(startErea, searchErea) >= 1){
			continue;
		}

		//�T�����q�̕����������P�ʃx�N�g��
		const Vec3 srchDirVec = Unit(m_iceSM[nearCluster]->GetOrgPos(srchI) - m_iceSM[nearCluster]->GetOrgCm());

		//�T���Ώۂ͋ߖT���q�̂��߁C���K�����Ȃ��ŋߎ��@���ʂ͑債�ĕς��Ȃ��������������Ƒ����Ȃ�
		//const Vec3 srchDirVec(m_iceSM[nearCluster]->GetOrgPos(srchI) - m_iceSM[nearCluster]->GetOrgCm());

		//�e���q��cos�ގ��x���ׂ�
		double similarity = dot(dirVec, srchDirVec);

		//�ގ��x���������q��I��
		if(mostSimPrtIndx == MAXINT || mostSimilarity <= similarity){
			mostSimPrtIndx = searchPrt;
			mostSimilarity = similarity;
		}

		debug_searchCompare++;
	}

	return pair<int, int>(nearCluster, MAXINT);
}

//�����T���F�����p�̗��q�̓Y�������擾
void CalcIteration::SelectExchangeParticleFromFarCluster(vector<vector<int>>& exchangeParticleList)
{
	//�Ƃ肠�����C�]���Ɠ���
	//����ł́C�������g�̗��q�͎��ւ��Ȃ��̂ň��S
	SelectAddParticleFromFarCluster(exchangeParticleList);
}

//���X�g����ɃN���X�^�֗��q��ǉ�
void CalcIteration::AddParticleToCluster(const vector<vector<int>>& addPIndxList)
{
	//�f�o�b�O
	int debug_addNum = 0;

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

			//�T���Ώۂ�������Ȃ������ꍇ
			if(addpIndx == MAXINT){
				continue;
			}

			//�N���X�^�ɗ��q�����݂��邩���m�F
			if(isAdd[addpIndx]){
				continue;
			}

			isAdd[addpIndx] = true;
			
			//�N���X�^�ɗ��q��ǉ�
			//addpIndx�Ԗڂ̂O�Ԗڂ̗��q�́C���q�ԍ�addpIndx�ƕۏ؂���Ă���
			int sldIndx = addpIndx * SM_DIM;
			Vec3 orgPos(m_iceSM[addpIndx]->GetOrgPos(0));
			Vec3 pos(sldPos[sldIndx+0], sldPos[sldIndx+1], sldPos[sldIndx+2]);
			Vec3 vel(sldVel[sldIndx+0], sldVel[sldIndx+1], sldVel[sldIndx+2]);

			int check = m_iceSM[i]->AddAnotherClusterVertex(orgPos, pos, vel, 1.0, addpIndx, 1.0, 0.0, 0);
			if(check == -1){
				cout << "over " << addPIndxList[i].size() << endl;
			}
			
			debug_addNum++;
		}

		m_iceSM[i]->CalcOrgCm();
		m_iceSM[i]->ClassifyAllOrgParticle();
	}

	cout << __FUNCTION__ << " addNum:" << debug_addNum << endl;
}

//���X�g�����ɗ��q�����ւ���@�Ή��֌W���ۏ؂���Ă��邩�ɒ���
void CalcIteration::ExchangeParticleToCluster(const vector<vector<int>>& exchangePIndxList)
{
	//�f�o�b�O
	int debug_replaceNum = 0;

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

		for(vector<int>::const_iterator it = exchangePIndxList[i].begin(); it != exchangePIndxList[i].end(); it++){
			int exchangePIndx = *it;

			if(m_iceSM[i]->GetNumVertices() >= MAXPARTICLE){
				break;
			}

			//�T���Ώۂ�������Ȃ������ꍇ
			if(exchangePIndx == MAXINT){
				continue;
			}

			//�N���X�^�ɗ��q�����݂��邩���m�F
			if(isAdd[exchangePIndx]){
				continue;
			}

			isAdd[exchangePIndx] = true;
			
			//���q������
			int sldIndx = exchangePIndx * SM_DIM;
			Vec3 orgPos(m_iceSM[exchangePIndx]->GetOrgPos(0));
			Vec3 pos(sldPos[sldIndx+0], sldPos[sldIndx+1], sldPos[sldIndx+2]);
			Vec3 vel(sldVel[sldIndx+0], sldVel[sldIndx+1], sldVel[sldIndx+2]);

			//�������闱�q�́C�N���X�^���ł̓Y��
			int oldIndx = it-exchangePIndxList[i].begin();
			//if(oldIndx == 0){
			//	cout << __FUNCTION__ << " ERROR! Base Exchange!" << endl;
			//}

			m_iceSM[i]->SetOriginalPos(oldIndx, orgPos);
			m_iceSM[i]->SetCurrentPos(oldIndx, pos);
			m_iceSM[i]->SetVelocity(oldIndx, vel);
			m_iceSM[i]->SetParticleIndx(oldIndx, exchangePIndx);

			debug_replaceNum++;
		}

		m_iceSM[i]->CalcOrgCm();
		m_iceSM[i]->ClassifyAllOrgParticle();
	}

	//cout << __FUNCTION__ << " replaceNum:" << debug_replaceNum << endl;
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

	//�g�������N���X�^�̔ԍ�
	vector<int> expandClusterIndx;

	//���񂾂��F�����̃N���X�^�ŉe���͈͂��g��
	ExchangeCluster_Far(expandClusterIndx);

	//����
	for(int itr = 1; itr < Ice_SM::GetIteration(); itr++){	
		//�N���X�^�̉e���͈͂��g�債�čč\��
		///ExpandCluster(searchFinishIndxes);
		//ExpandCluster_Far();
		
		//�����̗��q�ƌ��݂̗��q�����ւ���
		//ExchangeCluster_Far();

		//�������̉^���v�Z
		m_iceMove->StepObjMoveItr();

		//�ő̂ɂ�����ŏI�ʒu����
		m_iceConvo->StepConvolution();
	}
double end3 = counter3.End();

counter4.Start();
	//�R�s�[�𗘗p���ăN���X�^�̍\�������ɖ߂�
	ResetCluster(copyOriginalParticleIndxes, expandClusterIndx);
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