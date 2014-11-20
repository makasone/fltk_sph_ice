#include "Ice_InterPolationJudge_Normal.h"

typedef Ice_InterPolationJudge_Normal JudgeNormal;

JudgeNormal::Ice_InterPolationJudge_Normal(const vector<Ice_SM*>& iceSM)
{
	m_iceSM = iceSM;
}

JudgeNormal::~Ice_InterPolationJudge_Normal()
{
}

bool JudgeNormal::JudgeInterPolation(unsigned indx)
{
	if(m_iceSM[indx]->GetNumVertices() == 0) return false;

	return true;
}