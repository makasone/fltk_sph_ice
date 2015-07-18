//計算方法を選択する純粋仮想関数

#ifndef _ICE_CALC_METHOD_ITR_EXPAND_
#define _ICE_CALC_METHOD_ITR_EXPAND_

#include "Ice_CalcMethod.h"

#include "Ice_ClusterMove.h"
#include "Ice_JudgeMove.h"
#include "Ice_Convolution.h"

#include "QueryCounter.h"
#include "IceObject.h"
#include "Ice_SM.h"

#include <functional>
#include <algorithm>

using namespace std;

class Ice_CalcMethod_Itr_Expand: public Ice_CalcMethod
{
public:
	//Ice_CalcMethod_Itr_Expand(const vector<Ice_SM*>& iceSM, IceStructure* iceStrct, Ice_ClusterMove* clusterMove, Ice_Convolution* convo);
	Ice_CalcMethod_Itr_Expand(const vector<Ice_SM*>& iceSM, IceStructure* iceStrct, Ice_SimuMethod* simuMethod, Ice_Convolution* convo);
	~Ice_CalcMethod_Itr_Expand();
	
	//void SetObjMove(Ice_ClusterMove* clusterMove);
	void SetObjMove(Ice_SimuMethod* simuMethod);
	void SetConvolution(Ice_Convolution* convo);

	void StepObjMove();
	void StepObjMoveDebug();

private:
	void CalcVel();

	void CopyOriginalObject(vector<vector<unsigned>>& copyIndxes);
	void ResetCluster(const vector<vector<unsigned>>& copyIndxes, const vector<int>& exClstrIndxes);

	void GetExpandeCluster();

	void ExpandCluster(vector<int>& searchFinishIndxes);
	void ExpandCluster_Far();
	void ExpandCluster_Far_Step();

	void ExchangeCluster_Far(vector<int>& exClstrIndxes);
	void ConfirmExClusterIndx(const vector<vector<int>> exIndxList, vector<int>& exClstrIndxes);

	void SelectAddParticleFromNearestCluster(vector<vector<int>>& addParticleList, vector<int>& searchFinishIndxes);
	void SelectAddParticleFromFarCluster(vector<vector<int>>& addParticleList);
	void SelectExchangeParticleFromFarCluster(vector<vector<int>>& exchangeParticleList);

	const pair<int, int> SearchSimilarParticle(int nearCluster, const Vec3& dirVec, const Ice_SM::EreaData& startErea);

	void AddParticleToCluster(const vector<vector<int>>& addParticleList);	
	void ExchangeParticleToCluster(const vector<vector<int>>& replaceParticleList);

	void ContractCluster();

	void Debug_SelectAddParticleFromFarCluster(vector<vector<int>>& addParticleList);

private:
	vector<Ice_SM*> m_iceSM;

	IceStructure* m_iceStrct;

	//運動計算方法を扱うクラス
	//Ice_ClusterMove* m_iceMove;
	Ice_SimuMethod* m_simuMethod;

	//最終統合結果を求めるクラス
	Ice_Convolution* m_iceConvo;
};

#endif