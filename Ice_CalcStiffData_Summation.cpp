#include "Ice_CalcStiffData_Summation.h"

typedef Ice_CalcStiffData_Summation StiffSummation;

StiffSummation::Ice_CalcStiffData_Summation(const vector<Ice_SM*>& iceSM)
{
	m_iceSM = iceSM;
}

StiffSummation::~Ice_CalcStiffData_Summation()
{

}

void StiffSummation::StepUpdate()
{
}

void StiffSummation::StepUpdateItr()
{
}

float StiffSummation::StepCalcData()
{
	float sum = 0.0f;

	for(int i = 0; i < IceObject::GetParticleNum(); ++i)
	{	
		//全クラスタが対象
		if(m_iceSM[i]->GetNumVertices() == 0){	continue;	}

		float defAmount = m_iceSM[i]->GetDefAmount();

		sum += defAmount;
	}

	return sum;
}

void StiffSummation::StepCalcDataDebug()
{
}