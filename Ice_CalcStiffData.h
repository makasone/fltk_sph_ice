//閾値と比較するために用いる値を計算するクラス　純粋仮想関数

#ifndef _ICE_CALCSTIFFDATA_
#define _ICE_CALCSTIFFDATA_

class Ice_CalcStiffData
{
public:
	virtual float StepCalcData() = 0;
	virtual void StepCalcDataDebug() = 0;

	//virtual void SetConvoJudge(Ice_ConvoJudge* judge) = 0;
	//virtual Ice_ConvoJudge* GetConvoJudge() = 0;
};

#endif