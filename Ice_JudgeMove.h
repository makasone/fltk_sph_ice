//運動計算する粒子を選択する純粋仮想関数

#ifndef _ICE_JUDE_MOVE_
#define _ICE_JUDE_MOVE_

class Ice_JudgeMove
{
public:
	virtual bool JudgeMove(unsigned indx) = 0;
};

#endif