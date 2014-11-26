//計算クラスタ判定クラス
//疎に粒子を選択して運動計算を行う

#ifndef _ICE_JUDGE_MOVE_SPEARS_
#define _ICE_JUDGE_MOVE_SPEARS_

#include "Ice_JudgeMove.h"

#include "IceObject.h"
#include "IceStructure.h"

using namespace std;

class Ice_JudgeMove_Spears : public Ice_JudgeMove
{
public:
	Ice_JudgeMove_Spears(const vector<Ice_SM*>& iceSM, IceStructure* strct);
	~Ice_JudgeMove_Spears();
	
	bool JudgeMove(unsigned indx);
	bool JudgeMoveDebug(unsigned indx);

private:
	vector<Ice_SM*> m_iceSM;
	IceStructure* m_iceStrct;
};


#endif