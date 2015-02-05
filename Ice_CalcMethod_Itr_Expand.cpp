#include "Ice_CalcMethod_Itr_Expand.h"

typedef Ice_CalcMethod_Itr_Expand CalcIteration;

#define PNUM 2197
#define ADDLIMIT 50

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

	//初回の運動計算
	m_iceMove->StepObjMove();

	//各クラスタのオリジナル粒子を覚えておく
	//（実は，単なる拡張だけの処理ならオリジナル粒子の個数だけでよい）
	vector<vector<unsigned>> copyOriginalParticleIndxes(m_iceSM.size(), vector<unsigned>());
	CopyOriginalObject(copyOriginalParticleIndxes);

	//各反復時に走査済みの添字番号を記録しておく
	vector<int> searchFinishIndxes(m_iceSM.size(), 0);

	//反復
	for(int itr = 1; itr < Ice_SM::GetIteration(); itr++){	
		//固体における最終位置決定
		m_iceConvo->StepConvolution();
	
		//クラスタの影響範囲を拡大して再構成
		ExpandeCluster(searchFinishIndxes);

		//反復時の運動計算
		m_iceMove->StepObjMoveItr();
	}

	//コピーを利用してクラスタの構造を元に戻す
	ReplaceCluster(copyOriginalParticleIndxes);

	//速度算出
	CalcVel();
}

//速度算出
void CalcIteration::CalcVel()
{
	#pragma omp parallel for
	for(int i = 0; i < IceObject::GetParticleNum(); ++i)
	{	
		if(m_iceMove->GetJudgeMove()->JudgeMove(i) == false){	continue;	}
		m_iceSM[i]->integrateIteration();
	}
}

//元のオブジェクトをコピー
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

//クラスタを元に戻す
void CalcIteration::ReplaceCluster(const vector<vector<unsigned>>& copyIndxes)
{
	//copyIndxesに存在する粒子はオリジナル
	//存在しない粒子の情報を全て消す
	for(unsigned i = 0; i < m_iceSM.size(); i++){
		for(unsigned j = 0; j < m_iceSM[i]->GetIndxNum(); j++){
			//現在のm_iceSMに含まれる粒子で，copyObjに存在しないものを削除
			int jpIndx = m_iceSM[i]->GetParticleIndx(j);
			if(jpIndx == MAXINT){
				continue;
			}

			vector<unsigned>::const_iterator check = find(copyIndxes[i].begin(), copyIndxes[i].end(), jpIndx);
			if(check != copyIndxes[i].end()){
				continue;
			}

			//copyIndxesに存在しないので追加された粒子　削除する
			int ooIndx = m_iceSM[i]->SearchIndx(jpIndx);
			m_iceSM[i]->Remove(ooIndx);
		}

		m_iceSM[i]->CalcOrgCm();
	}
}

//再構築するクラスタを検出
void CalcIteration::GetExpandeCluster()
{
	for(unsigned i = 0; i < m_iceSM.size(); ++i){

	}
}

//クラスタの影響範囲を拡大して再構成
void CalcIteration::ExpandeCluster(vector<int>& srFnIndxes)
{
	ExpandeCluster_Test(srFnIndxes);
	//ExpandCluster_Test2(srFnIndxes);
}

void CalcIteration::ExpandeCluster_Test(vector<int>& sfFnIndxes)
{
	//近傍に存在する粒子とその初期位置のリストを作成
	vector<vector<int>> addPIndxList;
	addPIndxList.resize(m_iceSM.size(), vector<int>());

	for(unsigned i = 0; i < m_iceSM.size(); ++i){
		if(m_iceSM[i]->GetNumVertices() >= MAXPARTICLE){
			continue;
		}

		//粒子を追加しているかどうかのフラグ
		bool isAdd[PNUM] = {};
		int limitPrtNum = MAXPARTICLE - m_iceSM[i]->GetNumVertices();

		//現在のクラスタに含まれている粒子は追加しない
		for(unsigned j = 0; j < m_iceSM[i]->GetIndxNum(); j++){
			int orgPrt = m_iceSM[i]->GetParticleIndx(j);

			//追加粒子リストに含まれていないか
			if(orgPrt == MAXINT || isAdd[orgPrt]){	
				continue;
			}
			
			isAdd[orgPrt] = true;
		}

		//近傍クラスタを探索
		//既に探索済みのクラスタは除く
		unsigned searchStart = sfFnIndxes[i];
		for(unsigned clsI = searchStart; clsI < m_iceSM[clsI]->GetIndxNum(); clsI++){
			int nearCluster = m_iceSM[i]->GetParticleIndx(clsI);
			if(nearCluster == MAXINT || i == nearCluster){	
				continue;
			}

			//jpIndxのクラスタに含まれる粒子を取得
			for(unsigned prtI = 0; prtI < m_iceSM[nearCluster]->GetIndxNum(); prtI++){
				int nrPrt = m_iceSM[nearCluster]->GetParticleIndx(prtI);
				//追加粒子リストに含まれていないか
				if(nrPrt == MAXINT || isAdd[nrPrt]){	
					continue;
				}

				//クラスタに含められる粒子数は最大300個
				if(addPIndxList[i].size() >= limitPrtNum
				|| addPIndxList[i].size() >= ADDLIMIT){
					break ;
				}

				isAdd[nrPrt] = true;
				addPIndxList[i].push_back(nrPrt);
			}

			//クラスタに含められる粒子数は最大300個
			if(addPIndxList[i].size() >= limitPrtNum
			|| addPIndxList[i].size() >= ADDLIMIT){
				break ;
			}
		}
	}

	//クラスタに粒子を追加
	const float* sldPos = Ice_SM::GetSldPosPointer();
	const float* sldVel = Ice_SM::GetSldVelPointer();

	for(unsigned i = 0; i < m_iceSM.size(); i++){
		bool isAdd[PNUM] = {};
		//クラスタに含まれる粒子のフラグオン
		for(int j = 0; j < m_iceSM[i]->GetIndxNum(); j++){
			int orgPrtI = m_iceSM[i]->GetParticleIndx(j);
			if(orgPrtI == MAXINT || isAdd[orgPrtI]){
				continue;
			}

			isAdd[orgPrtI] = true;
		}

		for(vector<int>::iterator it = addPIndxList[i].begin(); it != addPIndxList[i].end(); it++){
			int addpIndx = *it;

			if(m_iceSM[i]->GetNumVertices() >= MAXPARTICLE){
				break;
			}

			//クラスタに粒子が存在するかを確認
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
		sfFnIndxes[i] += m_iceSM[i]->GetNumVertices();
	}
}

void CalcIteration::ExpandCluster_Test2(vector<int>& sfFnIndxes)
{
	//処理時間の計測
QueryCounter counter1;
QueryCounter counter2;

counter1.Start();
	//近傍に存在する粒子とその初期位置のリストを作成
	vector<vector<int>> addPIndxList;
	addPIndxList.resize(m_iceSM.size(), vector<int>());

	for(unsigned i = 0; i < m_iceSM.size(); ++i){
		if(m_iceSM[i]->GetNumVertices() >= MAXPARTICLE){
			continue;
		}

		//粒子を追加しているかどうかのフラグ
		bool isAdd[PNUM] = {};
		int limitPrtNum = MAXPARTICLE - m_iceSM[i]->GetNumVertices();

		//現在のクラスタに含まれている粒子は追加しない
		for(unsigned j = 0; j < m_iceSM[i]->GetIndxNum(); j++){
			int orgPrt = m_iceSM[i]->GetParticleIndx(j);

			//追加粒子リストに含まれていないか
			if(orgPrt == MAXINT || isAdd[orgPrt]){	
				continue;
			}
			
			isAdd[orgPrt] = true;
		}

		//近傍クラスタを探索
		//既に探索済みのクラスタは除く
		unsigned searchStart = sfFnIndxes[i];
		for(unsigned clsI = searchStart; clsI < m_iceSM[clsI]->GetIndxNum(); clsI++){
			int nearCluster = m_iceSM[i]->GetParticleIndx(clsI);
			if(nearCluster == MAXINT || i == nearCluster){	
				continue;
			}

			//jpIndxのクラスタに含まれる粒子を取得
			for(unsigned prtI = 0; prtI < m_iceSM[nearCluster]->GetIndxNum(); prtI++){
				int nrPrt = m_iceSM[nearCluster]->GetParticleIndx(prtI);
				//追加粒子リストに含まれていないか
				if(nrPrt == MAXINT || isAdd[nrPrt]){	
					continue;
				}

				//クラスタに含められる粒子数は最大300個
				if(addPIndxList[i].size() >= limitPrtNum
				|| addPIndxList[i].size() >= ADDLIMIT){
					break ;
				}

				isAdd[nrPrt] = true;
				addPIndxList[i].push_back(nrPrt);
			}

			//クラスタに含められる粒子数は最大300個
			if(addPIndxList[i].size() >= limitPrtNum
			|| addPIndxList[i].size() >= ADDLIMIT){
				break ;
			}
		}
	}
double end1 = counter1.End();

counter2.Start();
	//クラスタに粒子を追加
	const float* sldPos = Ice_SM::GetSldPosPointer();
	const float* sldVel = Ice_SM::GetSldVelPointer();

	for(unsigned i = 0; i < m_iceSM.size(); i++){
		bool isAdd[PNUM] = {};
		//クラスタに含まれる粒子のフラグオン
		for(int j = 0; j < m_iceSM[i]->GetIndxNum(); j++){
			int orgPrtI = m_iceSM[i]->GetParticleIndx(j);
			if(orgPrtI == MAXINT || isAdd[orgPrtI]){
				continue;
			}

			isAdd[orgPrtI] = true;
		}

		for(vector<int>::iterator it = addPIndxList[i].begin(); it != addPIndxList[i].end(); it++){
			int addpIndx = *it;

			if(m_iceSM[i]->GetNumVertices() >= MAXPARTICLE){
				break;
			}

			//クラスタに粒子が存在するかを確認
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
		sfFnIndxes[i] += m_iceSM[i]->GetNumVertices();
	}
double end2 = counter2.End();

	cout << "Time:" << "Search:" << end1 << "," << "Add:" << end2 << endl;
}

//クラスタの影響範囲を縮小して再構成
void CalcIteration::ContractCluster()
{
	//各クラスタの影響範囲を広げるために近傍にあるクラスタの粒子を追加する
	for(int i = 0; i < IceObject::GetParticleNum(); ++i){

	}
}

//デバッグ
void CalcIteration::StepObjMoveDebug()
{
	//初回
	m_iceMove->StepObjMoveDebug();

	//反復
	for(int itr = 1; itr < Ice_SM::GetIteration(); itr++)
	{
		m_iceConvo->StepConvolution();
		m_iceMove->StepObjMoveItrDebug();
	}

	//速度算出
	#pragma omp parallel for
	for(int i = 0; i < IceObject::GetParticleNum(); ++i)
	{	
		if(m_iceMove->GetJudgeMove()->JudgeMoveDebug(i) == false){	continue;	}
		m_iceSM[i]->integrateIteration();
	}
}

void CalcIteration::StepObjMoveTest()
{
QueryCounter counter1;
QueryCounter counter2;
QueryCounter counter3;
QueryCounter counter4;

counter1.Start();
	//初回の運動計算
	m_iceMove->StepObjMove();
double end1 = counter1.End();

counter2.Start();
	//各クラスタのオリジナル粒子を覚えておく
	//（実は，単なる拡張だけの処理ならオリジナル粒子の個数だけでよい）
	vector<vector<unsigned>> copyOriginalParticleIndxes(m_iceSM.size(), vector<unsigned>());
	CopyOriginalObject(copyOriginalParticleIndxes);

	//各反復時に走査済みの添字番号を記録しておく
	vector<int> searchFinishIndxes(m_iceSM.size(), 0);
double end2 = counter2.End();

counter3.Start();
	//反復
	for(int itr = 1; itr < Ice_SM::GetIteration(); itr++){	
		//固体における最終位置決定
		m_iceConvo->StepConvolution();
	
		//クラスタの影響範囲を拡大して再構成
		ExpandeCluster(searchFinishIndxes);

		//反復時の運動計算
		m_iceMove->StepObjMoveItr();
	}
double end3 = counter3.End();

counter4.Start();
	//コピーを利用してクラスタの構造を元に戻す
	ReplaceCluster(copyOriginalParticleIndxes);
double end4 = counter4.End();

	//速度算出
	CalcVel();

	//計算時間の比較
	cout << "Time:"
		<< " first:" << end1
		<< " prepare:" << end2
		<< " itr:" << end3
		<< " replace:" << end4
		<< endl;
}