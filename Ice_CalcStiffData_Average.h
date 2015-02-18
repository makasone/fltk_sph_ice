//総変形量を計算するクラス

#ifndef _ICE_CALCSTIFFDATA_AVE_
#define _ICE_CALCSTIFFDATA_AVE_

#include "Ice_CalcStiffData.h"
#include "IceObject.h"
#include "Ice_SM.h"

using namespace std;

class Ice_CalcStiffData_Average : public Ice_CalcStiffData
{
public:
	Ice_CalcStiffData_Average(const vector<Ice_SM*>& iceSM);
	~Ice_CalcStiffData_Average();

	void StepUpdate();
	void StepUpdateItr();

	float StepCalcData();
	void StepCalcDataDebug();


private:
	vector<Ice_SM*> m_iceSM;
};

#endif