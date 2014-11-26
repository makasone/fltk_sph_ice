//計算結果に反映させるクラスタを選択するクラス
//疎にクラスタを選択

#ifndef _ICE_CONVO_JUDGE_SPEARS_
#define _ICE_CONVO_JUDGE_SPEARS_

#include "Ice_ConvoJudge.h"

#include "IceStructure.h"
#include "Ice_SM.h"

using namespace std;

class Ice_ConvoJudge_Spears : public Ice_ConvoJudge
{
public:
	Ice_ConvoJudge_Spears(const vector<Ice_SM*>& iceSM, IceStructure* iceStrct);
	~Ice_ConvoJudge_Spears();

	bool JudgeConvolution(unsigned indx);

private:
	vector<Ice_SM*> m_iceSM;
	IceStructure* m_iceStrct;
};

#endif