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

	//閾値用データ算出クラス
	//m_iceCalcStiff = new Ice_CalcStiffData_Summation(iceSM);	//総変形量
	m_iceCalcStiff = new Ice_CalcStiffData_Average(iceSM);		//平均変化量
	//m_iceCalcStiff = new Ice_CalcStiffData_StdDevision(iceSM, m_iceMove->GetJudgeMove());		//粒子位置の分散
	//m_iceCalcStiff = new Ice_CalcStiffData_CompareRigid(iceSM, m_iceMove->GetJudgeMove());	//剛体との差分
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
	//初回の運動計算
	//m_iceMove->StepObjMove();
	m_simuMethod->StepObjMove();

	//rigidの最終結果　重力を反映しているバージョン
	m_iceCalcStiff->StepUpdate();

	//運動計算
	m_iceConvo->StepConvolution();

	//各クラスタのオリジナル粒子を覚えておく
	//（実は，単なる拡張だけの処理ならオリジナル粒子の個数だけでよい）
	vector<vector<unsigned>> copyOriginalParticleIndxes(m_iceSM.size(), vector<unsigned>());
	CopyOriginalObject(copyOriginalParticleIndxes);

	//各反復時に走査済みの添字番号を記録しておく
	vector<int> searchFinishIndxes(m_iceSM.size(), 0);

	//反復するか判定するための値
	float threshold = m_iceCalcStiff->StepCalcData();

	//反復	
	while(threshold > Ice_SM::GetItrStiffness()){
		//クラスタの影響範囲を拡大して再構成
		ExpandeCluster(searchFinishIndxes);

		//反復時の運動計算
		//m_iceMove->StepObjMoveItr();
		m_simuMethod->StepObjMoveItr();

		//運動計算
		m_iceConvo->StepConvolution();

		//反復するか判定するための値
		threshold = m_iceCalcStiff->StepCalcData();
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
		//if(m_iceMove->GetJudgeMove()->JudgeMove(i) == false){	continue;	}
		if(m_simuMethod->GetJudgeMove()->JudgeMove(i) == false){	continue;	}
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

//クラスタの影響範囲を拡大して再構成
void CalcIteration::ExpandeCluster(vector<int>& sfFnIndxes)
{
	//処理時間の計測
QueryCounter counter1;
QueryCounter counter2;

counter1.Start();
	//近傍に存在する粒子とその初期位置のリストを作成
	vector<vector<int>> addPIndxList;
	addPIndxList.resize(m_iceSM.size(), vector<int>(ADDLIMIT, MAXINT));

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
double end1 = counter1.End();

counter2.Start();
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

		for(vector<int>::iterator it = addPIndxList[i].begin(); it != addPIndxList[i].end(); it++){
			int addpIndx = *it;

			if(addpIndx == MAXINT || m_iceSM[i]->GetNumVertices() >= MAXPARTICLE){
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
		m_iceSM[i]->ClassifyAllOrgParticle();
	}
double end2 = counter2.End();

	//cout << "Time:" << "Search:" << end1 << "," << "Add:" << end2 << endl;
}

void CalcIteration::StepObjMoveDebug()
{
	//初回の運動計算
	//m_iceMove->StepObjMove();
	m_simuMethod->StepObjMove();

	//rigidの最終結果　重力を反映しているバージョン
	m_iceCalcStiff->StepUpdate();

	//運動計算
	m_iceConvo->StepConvolution();

	//各クラスタのオリジナル粒子を覚えておく
	//（実は，単なる拡張だけの処理ならオリジナル粒子の個数だけでよい）
	vector<vector<unsigned>> copyOriginalParticleIndxes(m_iceSM.size(), vector<unsigned>());
	CopyOriginalObject(copyOriginalParticleIndxes);

	//各反復時に走査済みの添字番号を記録しておく
	vector<int> searchFinishIndxes(m_iceSM.size(), 0);

	//反復するか判定するための値
	float threshold = m_iceCalcStiff->StepCalcData();

	int count = 0;

	//反復	
	while(threshold > Ice_SM::GetItrStiffness()){
		//クラスタの影響範囲を拡大して再構成
		ExpandeCluster(searchFinishIndxes);

		//反復時の運動計算
		//m_iceMove->StepObjMoveItr();
		m_simuMethod->StepObjMoveItr();

		//運動計算
		m_iceConvo->StepConvolution();

		//反復するか判定するための値
		threshold = m_iceCalcStiff->StepCalcData();

		count++;

		cout << __FUNCTION__ << " thrshold" << threshold << " count:" << count << endl;
	}

	//コピーを利用してクラスタの構造を元に戻す
	ReplaceCluster(copyOriginalParticleIndxes);

	//速度算出
	CalcVel();

//デバッグ
	cout << __FUNCTION__ << " " << count << endl;

	string result = RESULT_DATA_PATH;
	result += "Itr_Stiffness_CountNum.txt";
	ofstream ofs(result, ios::app);
	ofs << count << endl;
}

//硬さの閾値を測定する処理のデバッグ
void CalcIteration::DebugStiffness()
{
	m_iceCalcStiff->StepCalcDataDebug();
}