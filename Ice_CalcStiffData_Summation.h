//総変形量を計算するクラス

#ifndef _ICE_CALCSTIFFDATA_SUM_
#define _ICE_CALCSTIFFDATA_SUM_

#include "Ice_CalcStiffData.h"
#include "IceObject.h"
#include "Ice_SM.h"

using namespace std;

class Ice_CalcStiffData_Summation : public Ice_CalcStiffData
{
public:
	Ice_CalcStiffData_Summation(const vector<Ice_SM*>& iceSM);
	~Ice_CalcStiffData_Summation();

	void StepUpdate();
	void StepUpdateItr();

	float StepCalcData();
	void StepCalcDataDebug();


private:
	vector<Ice_SM*> m_iceSM;
};

#endif