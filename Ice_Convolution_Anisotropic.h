//�N���X�^�̌v�Z���ʕ�ԃN���X
//�w�肵�������x�N�g����D�悷��悤�ɏd�ݕt��

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

	Ice_ConvoJudge* m_iceJudge;		//��Ԃɗp����N���X�^�𔻒肷��I�u�W�F�N�g

	Vec3 m_dirVec;
};


#endif