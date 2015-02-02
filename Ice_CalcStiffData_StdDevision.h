//位置の標準偏差を計算するクラス

#ifndef _ICE_CALCSTIFFDATA_DIVISION_
#define _ICE_CALCSTIFFDATA_DIVISION_

#include "Ice_CalcStiffData.h"
#include "IceObject.h"
#include "Ice_SM.h"

using namespace std;

class Ice_CalcStiffData_StdDevision : public Ice_CalcStiffData
{
public:
	Ice_CalcStiffData_StdDevision(const vector<Ice_SM*>& iceSM, Ice_JudgeMove* judge);
	~Ice_CalcStiffData_StdDevision();

	float StepCalcData();
	void StepCalcDataDebug();


private:
	vector<Ice_SM*> m_iceSM;

	Ice_JudgeMove* m_iceJudge;
};

#endif