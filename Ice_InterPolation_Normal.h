//クラスタの計算結果補間クラス
//計算結果の平均で補間

#ifndef _ICE_INTER_POLATION_NORMAL_
#define _ICE_INTER_POLATION_NORMAL_

#include "Ice_InterPolation.h"

#include "IceObject.h"
#include "IceStructure.h"
#include "Ice_SM.h"

using namespace std;

class Ice_InterPolation_Normal : public Ice_InterPolation
{
public:
	Ice_InterPolation_Normal(const vector<Ice_SM*>& iceSM, IceStructure* iceStrct);
	~Ice_InterPolation_Normal();

	void SetIntrpJudge(Ice_InterPolationJudge* judge);
	Ice_InterPolationJudge* GetIntrpJudge();
	
	void StepInterPolation();
	void StepInterPolationItr();

private:
	vector<Ice_SM*> m_iceSM;
	IceStructure* m_iceStrct;

	Ice_InterPolationJudge* m_iceJudge;	//補間に用いるクラスタを判定するオブジェクト
};


#endif