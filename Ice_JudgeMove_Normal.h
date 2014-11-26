//�v�Z�N���X�^����N���X
//�N���X�^�ɗ��q���܂܂�Ă�����v�Z����

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