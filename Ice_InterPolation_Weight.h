//クラスタの計算結果補間クラス
//変形量で重みをつけて平均を取る

#ifndef _ICE_INTER_POLATION_WEIGHT_
#define _ICE_INTER_POLATION_WEIGHT_

#include "Ice_InterPolation.h"

#include "IceObject.h"
#include "IceStructure.h"
#include "Ice_SM.h"

using namespace std;

class Ice_InterPolation_Weight : public Ice_InterPolation
{
public:
	Ice_InterPolation_Weight(const vector<Ice_SM*>& iceSM, IceStructure* iceStrct);
	~Ice_InterPolation_Weight();

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