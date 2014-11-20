//�^���v�Z�N���X
//�e�N���X�^�̃f�[�^����ׂ��p�X�𗘗p���ĉ^���v�Z��������

#ifndef _ICE_CLUSTER_MOVE_FASTPATH_
#define _ICE_CLUSTER_MOVE_FASTPATH_

#include "Ice_ClusterMove.h"

#include "IceObject.h"
#include "Ice_SM.h"
#include "Surf_SM.h"

using namespace std;

class Ice_ClusterMove_FastPath : public Ice_ClusterMove
{
public:
	Ice_ClusterMove_FastPath(const vector<Ice_SM*>& iceSM, Surf_SM* surfSM);
	~Ice_ClusterMove_FastPath();

	void SetJudgeMove(Ice_JudgeMove* judge);
	Ice_JudgeMove* GetJudgeMove();

	void StepObjMove();
	void StepObjMoveItr();

private:
	vector<Ice_SM*> m_iceSM;
	Surf_SM* m_surfSM;

	Ice_JudgeMove* m_iceJudge;	//�^���v�Z����N���X�^�𔻒肷��I�u�W�F�N�g
};


#endif