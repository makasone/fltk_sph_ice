#include "Ice_InterPolationJudge_Spears.h"

typedef Ice_InterPolationJudge_Spears JudgeSpears;

JudgeSpears::Ice_InterPolationJudge_Spears(const vector<Ice_SM*>& iceSM, IceStructure* iceStrct)
{
	m_iceSM = iceSM;
	m_iceStrct = iceStrct;
}

JudgeSpears::~Ice_InterPolationJudge_Spears()
{
}

bool JudgeSpears::JudgeInterPolation(unsigned indx)
{
	if(m_iceStrct->GetMotionCalcCluster(indx) == 0)	return false;
	if(m_iceSM[indx]->GetNumVertices() == 0)		return false;

	return true;
}