//臒l�Ɣ�r���邽�߂ɗp����l���v�Z����N���X�@�������z�֐�

#ifndef _ICE_CALCSTIFFDATA_
#define _ICE_CALCSTIFFDATA_

class Ice_CalcStiffData
{
public:
	virtual float StepCalcData() = 0;
	virtual void StepCalcDataDebug() = 0;

	//virtual void SetConvoJudge(Ice_ConvoJudge* judge) = 0;
	//virtual Ice_ConvoJudge* GetConvoJudge() = 0;
};

#endif