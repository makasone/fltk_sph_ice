//運動計算手法クラス Shape Matching法

#ifndef _ICE_SIMU_METHOD_OP_
#define _ICE_SIMU_METHOD_OP_

#include "Ice_ClusterMove.h"
#include "Ice_SimuMethod.h"

#include "IceObject.h"
#include "Ice_SM.h"
#include "Ice_OrientedParticle.h"
#include "OrientedParticle.h"

typedef OrientedParticleBaseElasticObject ElasticObj;

using namespace std;

class Ice_SimuMethod_OP : public Ice_SimuMethod
{
public:
	Ice_SimuMethod_OP(const vector<Ice_SM*>& iceSM);
	Ice_SimuMethod_OP(const vector<ElasticObj*>& elasticObj, const vector<OrientedParticle*>& particles, IceStructure* iceStrct);
	~Ice_SimuMethod_OP();

	void SetJudgeMove(Ice_JudgeMove* judge);
	Ice_JudgeMove* GetJudgeMove();

	void StepObjMove();
	void StepObjMoveItr();

	void StepObjUpdate();
	void StepObjUpdateItr();
	void StepObjUpdateItrEnd();

	void StepObjMoveDebug();
	void StepObjMoveItrDebug();

private:
	vector<Ice_SM*> m_iceSM;

	vector<ElasticObj*> m_elasticObj;
	vector<OrientedParticle*> m_vOrientedPrtes;
	
	IceStructure* m_iceStrct;
	Ice_JudgeMove* m_iceJudge;	//運動計算するクラスタを判定するオブジェクト
};


#endif