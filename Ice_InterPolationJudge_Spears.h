//計算結果に反映させるクラスタを選択するクラス
//疎にクラスタを選択

#ifndef _ICE_JUDE_INTERPOLATION_JUDGE_SPEARS_
#define _ICE_JUDE_INTERPOLATION_JUDGE_SPEARS_

#include "Ice_InterPolationJudge.h"

#include "IceStructure.h"
#include "Ice_SM.h"

using namespace std;

class Ice_InterPolationJudge_Spears : public Ice_InterPolationJudge
{
public:
	Ice_InterPolationJudge_Spears(const vector<Ice_SM*>& iceSM, IceStructure* iceStrct);
	~Ice_InterPolationJudge_Spears();

	bool JudgeInterPolation(unsigned indx);

private:
	vector<Ice_SM*> m_iceSM;
	IceStructure* m_iceStrct;
};

#endif