//閾値を指定して反復計算

#ifndef _ICE_CALC_METHOD_ITR_STIFF_
#define _ICE_CALC_METHOD_ITR_STIFF_

#include "Ice_CalcMethod.h"

#include "Ice_ClusterMove.h"
#include "Ice_JudgeMove.h"
#include "Ice_Convolution.h"
#include "Ice_CalcStiffData.h"
#include "Ice_CalcStiffData_Summation.h"
#include "Ice_CalcStiffData_Average.h"
#include "Ice_CalcStiffData_StdDevision.h"
#include "Ice_CalcStiffData_CompareRigid.h"

#include "IceObject.h"
#include "Ice_SM.h"

using namespace std;

class Ice_CalcMethod_Itr_Stiffness : public Ice_CalcMethod
{
public:
	//Ice_CalcMethod_Itr_Stiffness(const vector<Ice_SM*>& iceSM, Ice_ClusterMove* clusterMove, Ice_Convolution* convo);
	Ice_CalcMethod_Itr_Stiffness(const vector<Ice_SM*>& iceSM, Ice_SimuMethod* simuMethod, Ice_Convolution* convo);
	~Ice_CalcMethod_Itr_Stiffness();
	
	//void SetObjMove(Ice_ClusterMove* clusterMove);
	void SetObjMove(Ice_SimuMethod* simuMethod);
	void SetConvolution(Ice_Convolution* convo);

	void StepObjMove();
	void StepObjMoveDebug();

	void DebugStiffness();

private:
	vector<Ice_SM*> m_iceSM;

	//運動計算方法を扱うクラス
	//Ice_ClusterMove* m_iceMove;
	Ice_SimuMethod* m_iceSimu;

	//最終統合結果を求めるクラス
	Ice_Convolution* m_iceConvo;

	//閾値と比較する値を計測するためのクラス
	Ice_CalcStiffData* m_iceCalcStiff;
};

#endif