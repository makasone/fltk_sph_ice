//��Ԃŗp�����@��I�����鏃�����z�֐�

#ifndef _ICE_CONVOLUTION_
#define _ICE_CONVOLUTION_

#include "Ice_ConvoJudge.h"

class Ice_Convolution
{
public:
	virtual void StepConvolution() = 0;
	virtual void StepConvolutionDebug() = 0;

	virtual void SetConvoJudge(Ice_ConvoJudge* judge) = 0;
	virtual Ice_ConvoJudge* GetConvoJudge() = 0;
};

#endif