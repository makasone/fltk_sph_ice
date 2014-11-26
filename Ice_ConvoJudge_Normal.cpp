#include "Ice_ConvoJudge_Normal.h"

typedef Ice_ConvoJudge_Normal JudgeNormal;

JudgeNormal::Ice_ConvoJudge_Normal(const vector<Ice_SM*>& iceSM)
{
	m_iceSM = iceSM;
}

JudgeNormal::~Ice_ConvoJudge_Normal()
{
}

bool JudgeNormal::JudgeConvolution(unsigned indx)
{
	if(m_iceSM[indx]->GetNumVertices() == 0) return false;

	return true;
}