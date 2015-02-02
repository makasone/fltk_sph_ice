//クラスタの計算結果補間クラス
//指定した方向ベクトルを優先するように重み付け

#ifndef _ICE_CONVOLUTION_ANISOTROPIC_
#define _ICE_CONVOLUTION_ANISOTROPIC_

#include "Ice_Convolution.h"

#include "IceObject.h"
#include "IceStructure.h"
#include "Ice_SM.h"

using namespace std;

class Ice_Convolution_Anisotropic : public Ice_Convolution
{
public:
	Ice_Convolution_Anisotropic(const vector<Ice_SM*>& iceSM, IceStructure* iceStrct);
	~Ice_Convolution_Anisotropic();

	void SetConvoJudge(Ice_ConvoJudge* judge);
	Ice_ConvoJudge* GetConvoJudge();
	
	void StepConvolution();

	void SetDirectionVector(Vec3 dirVec){	m_dirVec = dirVec;	}

	void StepConvolutionDebug();

private:
	vector<Ice_SM*> m_iceSM;
	IceStructure* m_iceStrct;

	Ice_ConvoJudge* m_iceJudge;		//補間に用いるクラスタを判定するオブジェクト

	Vec3 m_dirVec;
};


#endif