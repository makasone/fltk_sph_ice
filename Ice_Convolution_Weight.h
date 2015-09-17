//�N���X�^�̌v�Z���ʕ�ԃN���X
//�ό`�ʂŏd�݂����ĕ��ς����

#ifndef _ICE_CONVOLUTION_WEIGHT_
#define _ICE_CONVOLUTION_WEIGHT_

#include "Ice_Convolution.h"

#include "IceObject.h"
#include "IceStructure.h"
#include "Ice_SM.h"

#include "ElasticObject_OP.h"
#include "OrientedParticle.h"

using namespace std;

class Ice_Convolution_Weight : public Ice_Convolution
{
public:
	Ice_Convolution_Weight(const vector<Ice_SM*>& iceSM, IceStructure* iceStrct);
	Ice_Convolution_Weight(const vector<ElasticObj*>& elasticObj, const vector<OrientedParticle*>& particles, IceStructure* iceStrct);
	~Ice_Convolution_Weight();

	void SetConvoJudge(Ice_ConvoJudge* judge);
	Ice_ConvoJudge* GetConvoJudge();
	
	void StepConvolution();

	void SetKernelDegree(float dim){	m_kernelDegree = dim;	}

	void StepConvolutionDebug();

private:
	vector<Ice_SM*> m_iceSM;
	
	vector<ElasticObj*> m_elasticObj;
	vector<OrientedParticle*> m_vOrientedPrtes;

	IceStructure* m_iceStrct;

	Ice_ConvoJudge* m_iceJudge;		//��Ԃɗp����N���X�^�𔻒肷��I�u�W�F�N�g

	float m_kernelDegree;
};


#endif