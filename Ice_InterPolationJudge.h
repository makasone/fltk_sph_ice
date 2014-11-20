//計算結果に反映させるクラスタを判定する純粋仮想関数

#ifndef _ICE_JUDE_INTERPOLATION_JUDGE_
#define _ICE_JUDE_INTERPOLATION_JUDGE_

class Ice_InterPolationJudge
{
public:
	virtual bool JudgeInterPolation(unsigned indx) = 0;
};

#endif