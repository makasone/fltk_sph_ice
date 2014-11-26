//計算クラスタ判定クラス
//クラスタに粒子が含まれていたら計算する

#ifndef _ICE_JUDGE_MOVE_NORMAL_
#define _ICE_JUDGE_MOVE_NORMAL_

#include "Ice_JudgeMove.h"

#include "IceObject.h"
#include "IceStructure.h"

using namespace std;

class Ice_JudgeMove_Normal : public Ice_JudgeMove
{
public:
	Ice_JudgeMove_Normal(const vector<Ice_SM*>& iceSM);
	~Ice_JudgeMove_Normal();

	bool JudgeMove(unsigned indx);
	bool JudgeMoveDebug(unsigned indx);

private:
	vector<Ice_SM*> m_iceSM;
	IceStructure m_iceStrct;
};


#endif