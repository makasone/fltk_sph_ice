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
	//初回
	m_iceMove->StepObjMove();

	//元のオブジェクトをコピー
	Ice_SM copyObj[PNUM];

	//代入演算子でコピー
	for(int i = 0; i < IceObject::GetParticleNum(); ++i){
		copyObj[i] = *m_iceSM[i];
	}

	//反復
	for(int itr = 1; itr < Ice_SM::GetIteration(); itr++)
	{	
		m_iceConvo->StepConvolution();
	
		//クラスタの影響範囲を拡大して再構成
		ExpandeCluster();

		m_iceMove->StepObjMoveItr();

		//TODO: 変形量の大きいクラスタを検出
		//GetExpandeCluster(expandObj);
	}

	//ContractCluster(copyObj);

	//コピーの削除
	//DeleteCopyObject(copyObj);

	//コピーを利用してクラスタの構造を元に戻す
	ReplaceCluster(copyObj);

	//速度算出
	#pragma omp parallel for
	for(int i = 0; i < IceObject::GetParticleNum(); ++i)
	{	
		if(m_iceMove->GetJudgeMove()->JudgeMove(i) == false){	continue;	}
		m_iceSM[i]->integrateIteration();
	}
}

//元のオブジェクトをコピー
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
	//copyObjのm_iPindxesに存在する粒子はオリジナル
	//存在しない粒子の情報を全て消す

	//TODO: 必要な情報は，初期の粒子番号だけだった
	for(unsigned i = 0; i < m_iceSM.size(); i++){
		for(unsigned j = 0; j < m_iceSM[i]->GetIndxNum(); j++){
			//現在のm_iceSMに含まれる粒子で，copyObjに存在しないものを削除
			int jpIndx = m_iceSM[i]->GetParticleIndx(j);

			if(jpIndx == MAXINT || copyObj[i].CheckIndx(jpIndx)){
				continue;
			}

			//copyObjに存在しないので追加された粒子　削除する
			//cout << i << " delete " << jpIndx << endl;
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
void CalcIteration::ExpandeCluster()
{
	float* sldPos = Ice_SM::GetSldPosPointer();
	float* sldVel = Ice_SM::GetSldVelPointer();

	////テスト全ての粒子をremoveしてみる
	//for(unsigned i = 0; i < m_iceSM.size(); i++){
	//	for(unsigned j = 0; j < m_iceSM[i]->GetIndxNum(); j++){
	//		int jpIndx = m_iceSM[i]->GetParticleIndx(j);
	//		if(jpIndx == MAXINT){
	//			continue;
	//		}
	//	}
	//}

	////ファイルを読み込み
	//string filePath = RESULT_DATA_PATH;
	//filePath += "ExpandParticle.txt";
	//ofstream ofs(filePath);

	////ファイルの存在確認
	//if(ofs.fail()) {	cerr << filePath << " is do not exist.\n";	return;	}

	//近傍に存在する粒子とその初期位置のリストを作成
	vector<vector<pair<int, Vec3>>> addPIndxList;
	addPIndxList.resize(m_iceSM.size(), vector<pair<int, Vec3>>());

	for(unsigned i = 0; i < m_iceSM.size(); ++i){
		//リストを現在のクラスタに含まれている粒子で初期化
		for(unsigned j = 0; j < m_iceSM[i]->GetIndxNum(); j++){
			int jpIndx = m_iceSM[i]->GetParticleIndx(j);
			if(jpIndx == MAXINT){	
				continue;
			}

			//クラスタに含められる粒子数は最大300個
			if(addPIndxList[i].size() > 299){
				//cout << "init max" << endl;
				return ;
			}

			Vec3 orgPos = m_iceSM[i]->GetOrgPos(j);
			addPIndxList[i].push_back(pair<int, Vec3>(jpIndx, orgPos));
		}

		//近傍クラスタを探索
		for(unsigned j = 0; j < m_iceSM[i]->GetIndxNum(); j++){
			int jpIndx = m_iceSM[i]->GetParticleIndx(j);
			if(jpIndx == MAXINT){	
				continue;
			}

			//jpIndxのクラスタに含まれる粒子を取得
			for(unsigned k = 0; k < m_iceSM[jpIndx]->GetIndxNum(); k++){
				int kpIndx = m_iceSM[jpIndx]->GetParticleIndx(k);
				if(kpIndx == MAXINT){	
					continue;
				}

				//クラスタに含められる粒子数は最大300個
				if(addPIndxList[i].size() > 299){
					//cout << "search max" << endl;
					break ;
				}

				//追加粒子リストに含まれていないなら追加する
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

			//クラスタに含められる粒子数は最大300個
			if(addPIndxList[i].size() > 299){
				//cout << "search max" << endl;
				break ;
			}
		}

		////デバッグ
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
	
	//TODO: 粒子を追加
	for(unsigned i = 0; i < m_iceSM.size(); i++){
		for(vector<pair<int, Vec3>>::iterator it = addPIndxList[i].begin(); it != addPIndxList[i].end(); it++){
			int addpIndx = it->first;

			//クラスタに粒子が存在するかを確認
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

		//デバッグ
	}
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