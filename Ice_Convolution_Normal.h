//クラスタの計算結果補間クラス
//計算結果の平均で補間

#ifndef _ICE_CONVOLUTION_NORMAL_
#define _ICE_CONVOLUTION_NORMAL_

#include "Ice_Convolution.h"

#include "IceObject.h"
#include "IceStructure.h"
#include "Ice_SM.h"

#include "ElasticObject_OP.h"
#include "OrientedParticle.h"

using namespace std;

class Ice_Convolution_Normal : public Ice_Convolution
{
public:
	Ice_Convolution_Normal(const vector<Ice_SM*>& iceSM, IceStructure* iceStrct);
	Ice_Convolution_Normal(const vector<ElasticObj*>& elasticObj, const vector<OrientedParticle*>& particles, IceStructure* iceStrct);
	~Ice_Convolution_Normal();

	void SetConvoJudge(Ice_ConvoJudge* judge);
	Ice_ConvoJudge* GetConvoJudge();
	
	void StepConvolution();
	void StepConvolution2();

	void StepConvolutionDebug();

private:
	vector<Ice_SM*> m_iceSM;

	vector<ElasticObj*> m_elasticObj;
	vector<OrientedParticle*> m_vOrientedPrtes;

	IceStructure* m_iceStrct;

	Ice_ConvoJudge* m_iceJudge;		//補間に用いるクラスタを判定するオブジェクト
};


#endif