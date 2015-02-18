#include "Ice_CalcStiffData_Average.h"

typedef Ice_CalcStiffData_Average StiffAverage;

StiffAverage::Ice_CalcStiffData_Average(const vector<Ice_SM*>& iceSM)
{
	m_iceSM = iceSM;
}

StiffAverage::~Ice_CalcStiffData_Average()
{

}

void StiffAverage::StepUpdate()
{
}

void StiffAverage::StepUpdateItr()
{
}

float StiffAverage::StepCalcData()
{
	float sum = 0.0f;

	for(int i = 0; i < IceObject::GetParticleNum(); ++i)
	{	
		//全クラスタが対象
		if(m_iceSM[i]->GetNumVertices() == 0){	continue;	}

		float defAmount = m_iceSM[i]->GetDefAmount() / m_iceSM[i]->GetNumVertices();

		sum += defAmount;
	}

	return sum;
}

void StiffAverage::StepCalcDataDebug()
{
}