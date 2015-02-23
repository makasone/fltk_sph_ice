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

	//初回の運動計算
	m_iceMove->StepObjMove();

	//固体における最終位置決定
	m_iceConvo->StepConvolution();

	//各クラスタのオリジナル粒子を覚えておく
	//（実は，単なる拡張だけの処理ならオリジナル粒子の個数だけでよい）
	vector<vector<unsigned>> copyOriginalParticleIndxes(m_iceSM.size(), vector<unsigned>());
	CopyOriginalObject(copyOriginalParticleIndxes);

	//各反復時に走査済みの添字番号を記録しておく
	vector<int> searchFinishIndxes(m_iceSM.size(), 0);

	////初回だけ：遠くのクラスタで影響範囲を拡大
	//ExpandCluster_Far();

	//反復
	for(int itr = 1; itr < Ice_SM::GetIteration(); itr++){	
		////クラスタの影響範囲を拡大して再構成
		//ExpandCluster(searchFinishIndxes);

		//遠くのクラスタで影響範囲を拡大
		ExpandCluster_Far();

		//反復時の運動計算
		m_iceMove->StepObjMoveItr();

		//固体における最終位置決定
		m_iceConvo->StepConvolution();
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
		m_iceSM[i]->ClassifyAllOrgParticle();
	}
}

//再構築するクラスタを検出
void CalcIteration::GetExpandeCluster()
{
	for(unsigned i = 0; i < m_iceSM.size(); ++i){

	}
}

//クラスタの影響範囲を拡大して再構成
void CalcIteration::ExpandCluster(vector<int>& searchFinishIndxes)
{
		//処理時間の計測
QueryCounter counter1;
QueryCounter counter2;

counter1.Start();
	//近傍に存在する粒子とその初期位置のリストを作成
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
	//処理時間の計測
QueryCounter counter1;
QueryCounter counter2;

counter1.Start();
	//近傍に存在する粒子とその初期位置のリストを作成
	vector<vector<int>> addPIndxList;
	addPIndxList.resize(m_iceSM.size(), vector<int>(ADDLIMIT, MAXINT));

	SelectAddParticleFromFarCluster(addPIndxList);
double end1 = counter1.End();

counter2.Start();
	AddParticleToCluster(addPIndxList);
double end2 = counter2.End();

	cout << "Time:" << "Search:" << end1 << "," << "Add:" << end2 << endl;
}

//近傍探索：クラスタへ追加する粒子の添え字を取得
void CalcIteration::SelectAddParticleFromNearestCluster(vector<vector<int>>& addPIndxList, vector<int>& sfFnIndxes)
{
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

		//現在の粒子数兼添え字
		int addIndxNum = 0;

		//あるクラスタに含まれている粒子が対象
		for(unsigned clsI = searchStart; clsI < m_iceSM[i]->GetIndxNum(); clsI++){
			int nearCluster = m_iceSM[i]->GetParticleIndx(clsI);
			if(nearCluster == MAXINT || i == nearCluster){	
				continue;
			}

			//クラスタに含められる粒子数は最大300個
			if(addIndxNum >= limitPrtNum
			|| addIndxNum >= ADDLIMIT){
				break;
			}

			//jpIndxのクラスタに含まれる粒子を取得
			for(unsigned prtI = 0; prtI < m_iceSM[nearCluster]->GetIndxNum(); prtI++){
				//このクラスタを最後まで探索し終わった
				if(prtI == m_iceSM[nearCluster]->GetIndxNum()-1){
					sfFnIndxes[i]++;
				}

				int nrPrt = m_iceSM[nearCluster]->GetParticleIndx(prtI);

				//追加粒子リストに含まれていないか
				if(nrPrt == MAXINT || isAdd[nrPrt]){	
					continue;
				}

				//クラスタに含められる粒子数は最大300個
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

//階層探索：クラスタへ追加する粒子の添え字を取得
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

		//粒子を追加しているかどうかのフラグ
		bool isAdd[PNUM] = {};

		//現在のクラスタに含まれている粒子は追加しない
		for(unsigned j = 0; j < m_iceSM[i]->GetIndxNum(); j++){
			//粒子番号である
			int orgPrt = m_iceSM[i]->GetParticleIndx(j);
			if(orgPrt == MAXINT){	
				continue;
			}

			////これより前のクラスタの探索で追加されている場合
			//if(addPIndxList[i][j] != MAXINT){
			//	int preAddPrt = addPIndxList[i][j];
			//	isAdd[preAddPrt] = true;
			//	continue;
			//}
			
			isAdd[orgPrt] = true;
		}

		//指向性探索
		////追加する粒子数：初期のクラスタに含まれている粒子数
		//unsigned addParticleNum = m_iceSM[i]->GetNumVertices();

//QueryCounter counter;
//counter.Start();
		//クラスタに含まれている粒子を探索の起点とする
		for(unsigned clsI = 1; clsI < m_iceSM[i]->GetIndxNum(); clsI++){
			//debug_indxSum++;

			////既に追加済みかどうか
			//if(addPIndxList[i][clsI] != MAXINT){
			//	//cout << "CONTINU+E" << endl;
			//	debug_crossHold++;
			//	continue;
			//}

			//探索起点粒子
			int nearCluster = m_iceSM[i]->GetParticleIndx(clsI);
			
			//粒子番号である
			if(nearCluster == MAXINT){
				continue;
			}

			//ofs << "	pIndx:" << nearCluster << endl;

			//起点粒子の方向を示す単位ベクトル
			const Vec3 dirVec = Unit(m_iceSM[i]->GetOrgPos(clsI) - m_iceSM[i]->GetOrgCm());

			//起点粒子の領域
			const Ice_SM::EreaData startErea = m_iceSM[i]->erea(clsI);

			//似た方向ベクトルを持つ粒子を辿る
			pair<int, int> resultIndxes = SearchSimilarParticle(nearCluster, dirVec, startErea);

			int searchPrtIndx = resultIndxes.first;
			int searchClsIndx = resultIndxes.second;

			//探索粒子を追加
			if(isAdd[searchPrtIndx]){
				//既に含まれているなら，追加しない
				//TODO:ここで戻るということは，無駄な探索があるということ．なるべく無駄にならないように粒子を追加する
				continue;;
			}

			debug_searchCounter++;

			//最後に選んだ類似粒子をクラスタに追加する
			isAdd[searchPrtIndx] = true;
			addPIndxList[i][clsI] = searchPrtIndx;

			////ofs << "		addPIndx:" << searchPrtIndx
			////	<< "	nowIndxNum:" << clsI << endl;

			////TODO:持ち合い関係
			////探索が進まずに打ち切られた場合は持ち合い関係にはならない
			//if(searchClsIndx == MAXINT){
			//	continue;
			//}

			////変数の関係が逆でややこしいのに注意
			//int oIndx = m_iceSM[searchPrtIndx]->SearchIndx(searchClsIndx);
			////if(oIndx == MAXINT){
			////	//ofs << __FUNCTION__ << " ERROR" << endl;
			////	continue;
			////}

			////持ち合い関係を保存したいが，1対多関係になっている
			//if(addPIndxList[searchPrtIndx][oIndx] != MAXINT){
			//	//ofs << "1対多" << "searchPrtIndx:" << searchPrtIndx << "," << "searchClsIndx:" << searchClsIndx << endl;
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

//似た方向ベクトルを持つ粒子を辿って，遠い粒子を取得
pair<int, int> CalcIteration::SearchSimilarParticle(int nearCluster, const Vec3& dirVec, const Ice_SM::EreaData& startErea)
{
	//探索階層：適当に指定
	unsigned searchLayer = SEARCH_LAYER;

	//探索粒子の持つクラスタを取得し，含まれている粒子を調査
	unsigned nowLayer = 1;
	int mostSimPrtIndx = MAXINT;
	double mostSimilarity = 0.0;

	for(unsigned srchI = 0; srchI <= m_iceSM[nearCluster]->GetIndxNum(); srchI++){
		//終了判定
		if(srchI == m_iceSM[nearCluster]->GetIndxNum()){
			//探索を続けるかの判定
			 if(mostSimPrtIndx == MAXINT){
				//一つも候補が見つからず探索が進まなかった場合，現在のクラスタを探索結果として探索終了
				return pair<int, int>(nearCluster, MAXINT);
			}
			else if(nowLayer == searchLayer){
				//探索の回数制限を迎えた場合，探索終了
				return pair<int, int>(mostSimPrtIndx, nearCluster);
			}

			//探索を続けるなら，各変数を更新
			nowLayer++;
			
			srchI = 0;
			nearCluster = mostSimPrtIndx;
			
			mostSimPrtIndx = MAXINT;
			mostSimilarity = 0.0;

			continue;
		}

		int searchPrt = m_iceSM[nearCluster]->GetParticleIndx(srchI);

		//粒子番号である，自分自身でない，除外と判断されていない
		if(searchPrt == MAXINT || searchPrt == nearCluster){
			continue;
		}

		//探索の対象領域であるかを判定
		const Ice_SM::EreaData searchErea = m_iceSM[nearCluster]->erea(srchI);

		//領域間の距離で判定
		if(Ice_SM::CalcEreaDistance(startErea, searchErea) >= 1){
			continue;
		}

		//探索粒子の方向を示す単位ベクトル
		Vec3 srchDirVec = Unit(m_iceSM[nearCluster]->GetOrgPos(srchI) - m_iceSM[nearCluster]->GetOrgCm());

		//各粒子とcos類似度を比べる
		double similarity = dot(dirVec, srchDirVec);

		//類似度が高い粒子を選択
		if(mostSimPrtIndx == MAXINT || mostSimilarity <= similarity){
			mostSimPrtIndx = searchPrt;
			mostSimilarity = similarity;
		}
	}

	return pair<int, int>(nearCluster, MAXINT);
}

//リストを基にクラスタへ粒子を追加
void CalcIteration::AddParticleToCluster(const vector<vector<int>>& addPIndxList)
{
	//デバッグ
	int debug_addNum = 0;

	//クラスタに粒子を追加
	const float* sldPos = Ice_SM::GetSldPosPointer();
	const float* sldVel = Ice_SM::GetSldVelPointer();

	for(unsigned i = 0; i < m_iceSM.size(); i++){

		bool isAdd[PNUM] = {};

		//クラスタに含まれる粒子のフラグオン
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

			//farの場合はcontinue, normalはbreak
			if(addpIndx == MAXINT){
				//break;
				continue;
			}

			//クラスタに粒子が存在するかを確認
			if(isAdd[addpIndx]){
				continue;
			}

			isAdd[addpIndx] = true;

			//addpIndx番目の０番目の粒子は，粒子番号addpIndxと保証されている
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

//クラスタの影響範囲を縮小して再構成
void CalcIteration::ContractCluster()
{
	//各クラスタの影響範囲を広げるために近傍にあるクラスタの粒子を追加する
	for(int i = 0; i < IceObject::GetParticleNum(); ++i){

	}
}

//デバッグ-----------------------------------------------------------------
void CalcIteration::StepObjMoveDebug()
{
QueryCounter counter1;
QueryCounter counter2;
QueryCounter counter3;
QueryCounter counter4;

counter1.Start();
	//初回の運動計算
	m_iceMove->StepObjMove();

	//固体における最終位置決定
	m_iceConvo->StepConvolution();
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
		//クラスタの影響範囲を拡大して再構成
		///ExpandCluster(searchFinishIndxes);
		ExpandCluster_Far();

		//反復時の運動計算
		m_iceMove->StepObjMoveItr();

		//固体における最終位置決定
		m_iceConvo->StepConvolution();
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

//遠いクラスタから，追加するための粒子を選択
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

		//粒子を追加しているかどうかのフラグ
		bool isAdd[PNUM] = {};

		//現在のクラスタに含まれている粒子は追加しない
		for(unsigned j = 0; j < m_iceSM[i]->GetIndxNum(); j++){
			int orgPrt = m_iceSM[i]->GetParticleIndx(j);

			//追加粒子リストに含まれていないか
			if(orgPrt == MAXINT || isAdd[orgPrt]){	
				continue;
			}
			
			isAdd[orgPrt] = true;
		}

		//指向性探索
		////追加する粒子数：初期のクラスタに含まれている粒子数
		//unsigned addParticleNum = m_iceSM[i]->GetNumVertices();

		//探索階層：適当に指定
		unsigned searchLayer = SEARCH_LAYER;

		//現在の粒子数兼添え字
		int addIndxNum = 0;

		//クラスタに含まれている粒子を探索の起点とする
		for(unsigned clsI = 0; clsI < m_iceSM[i]->GetIndxNum(); clsI++){
			int nearCluster = m_iceSM[i]->GetParticleIndx(clsI);
			if(nearCluster == MAXINT || i == nearCluster){	
				continue;
			}

			//起点粒子の方向を示す単位ベクトル
			Vec3 dirVec = Unit(m_iceSM[i]->GetOrgPos(clsI) - m_iceSM[i]->GetOrgPos(0));

			//ofs << "	StartCluster:" << nearCluster
			//	<< "	dirVec:" << dirVec
			//	<< endl;

			//探索粒子の持つクラスタを取得し，含まれている粒子を調査
			unsigned nowLayer = 1;
			int mostSimPrtIndx = MAXINT;
			double mostSimilarity = 0.0;

			for(unsigned srchI = 0; srchI < m_iceSM[nearCluster]->GetIndxNum(); srchI++){
				int searchPrt = m_iceSM[nearCluster]->GetParticleIndx(srchI);
				if(searchPrt == MAXINT || searchPrt == nearCluster){	
					continue;
				}

				//探索粒子の方向を示す単位ベクトル
				Vec3 srchDirVec = Unit(m_iceSM[nearCluster]->GetOrgPos(srchI) - m_iceSM[nearCluster]->GetOrgPos(0));

				//各粒子とcos類似度を比べる
				double similarity = dot(dirVec, srchDirVec);

				//ofs << "		srchI:" << srchI
				//	<< "	searchPrt:" << searchPrt
				//	<< "	similarity:" << similarity
				//	<< "	srchDirVec:" << srchDirVec << endl;

				//類似度が高い粒子を選択
				if(mostSimPrtIndx == MAXINT || mostSimilarity <= similarity){
					mostSimPrtIndx = searchPrt;
					mostSimilarity = similarity;
				}

				//終了判定
				if(srchI < m_iceSM[nearCluster]->GetIndxNum() -1){
					continue;
				}

				//ofs << "	End"
				//	<< "	nowLayer:" << nowLayer
				//	<< "	nearClst:" << nearCluster
				//	<< "	mostSimPrtIndx" << mostSimPrtIndx
				//	<< "	mostSimilarity" << mostSimilarity << endl;

				//探索を続けるかの分岐
				if(nowLayer < searchLayer){
					//探索を続けるなら，各変数を更新
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

					//探索の回数制限を迎えると，探索終了
					//最後に選んだ類似粒子をクラスタに追加する
					//TODO: このままだと追加されない場合が多くなるので，stackを使う
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