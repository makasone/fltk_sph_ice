//位置の標準偏差を計算するクラス

#ifndef _ICE_CALCSTIFFDATA_COMPARE_RIGID_
#define _ICE_CALCSTIFFDATA_COMPARE_RIGID_

#include "Ice_CalcStiffData.h"
#include "IceObject.h"
#include "Ice_SM.h"

using namespace std;

class Ice_CalcStiffData_CompareRigid : public Ice_CalcStiffData
{
public:
	Ice_CalcStiffData_CompareRigid(const vector<Ice_SM*>& iceSM, Ice_JudgeMove* judge);
	~Ice_CalcStiffData_CompareRigid();

	void StepUpdate();
	void StepUpdateItr();
	
	float StepCalcData();
	void StepCalcDataDebug();

	void MakeRigidObj();
	Ice_SM* GetRigidPointer(){	return m_rigid;	}

private:
	vector<Ice_SM*> m_iceSM;
	Ice_SM* m_rigid;

	Ice_JudgeMove* m_iceJudge;
};

#endif