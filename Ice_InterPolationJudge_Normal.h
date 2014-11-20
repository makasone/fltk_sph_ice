//�v�Z���ʂɔ��f������N���X�^��I������N���X
//�S�N���X�^�𗱎q�ʒu�ɔ��f

#ifndef _ICE_JUDE_INTERPOLATION_JUDGE_NORMAL_
#define _ICE_JUDE_INTERPOLATION_JUDGE_NORMAL_

#include "Ice_InterPolationJudge.h"

#include "IceStructure.h"
#include "Ice_SM.h"

using namespace std;

class Ice_InterPolationJudge_Normal : public Ice_InterPolationJudge
{
public:
	Ice_InterPolationJudge_Normal(const vector<Ice_SM*>& iceSM);
	~Ice_InterPolationJudge_Normal();

	bool JudgeInterPolation(unsigned indx);

private:
	vector<Ice_SM*> m_iceSM;
};


#endif