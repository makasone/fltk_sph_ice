//計算結果に反映させるクラスタを判定する純粋仮想関数

#ifndef _ICE_CONVO_JUDGE_
#define _ICE_CONVO_JUDGE_

class Ice_ConvoJudge
{
public:
	virtual bool JudgeConvolution(unsigned indx) = 0;
};

#endif