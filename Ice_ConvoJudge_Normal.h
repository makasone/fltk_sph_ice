//計算結果に反映させるクラスタを選択するクラス
//全クラスタを粒子位置に反映

#ifndef _ICE_CONVO_JUDGE_NORMAL_
#define _ICE_CONVO_JUDGE_NORMAL_

#include "Ice_ConvoJudge.h"

#include "IceStructure.h"
#include "Ice_SM.h"

using namespace std;

class Ice_ConvoJudge_Normal : public Ice_ConvoJudge
{
public:
	Ice_ConvoJudge_Normal(const vector<Ice_SM*>& iceSM);
	~Ice_ConvoJudge_Normal();

	bool JudgeConvolution(unsigned indx);

private:
	vector<Ice_SM*> m_iceSM;
};


#endif