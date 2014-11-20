//�N���X�^�̉^���v�Z��@�N���X
//�S�N���X�^�Ōʂɉ^���v�Z

#ifndef _ICE_CLUSTER_MOVE_NORMAL_
#define _ICE_CLUSTER_MOVE_NORMAL_

#include "Ice_ClusterMove.h"

#include "IceObject.h"
#include "Ice_SM.h"

using namespace std;

class Ice_ClusterMove_Normal : public Ice_ClusterMove
{
public:
	Ice_ClusterMove_Normal(const vector<Ice_SM*>& iceSM);
	~Ice_ClusterMove_Normal();

	void SetJudgeMove(Ice_JudgeMove* judge);
	Ice_JudgeMove* GetJudgeMove();

	void StepObjMove();
	void StepObjMoveItr();

private:
	vector<Ice_SM*> m_iceSM;

	Ice_JudgeMove* m_iceJudge;	//�^���v�Z����N���X�^�𔻒肷��I�u�W�F�N�g
};


#endif