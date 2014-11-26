#include "Ice_ConvoJudge_Spears.h"

typedef Ice_ConvoJudge_Spears JudgeSpears;

JudgeSpears::Ice_ConvoJudge_Spears(const vector<Ice_SM*>& iceSM, IceStructure* iceStrct)
{
	m_iceSM = iceSM;
	m_iceStrct = iceStrct;
}

JudgeSpears::~Ice_ConvoJudge_Spears()
{
}

bool JudgeSpears::JudgeConvolution(unsigned indx)
{
	if(m_iceStrct->GetMotionCalcCluster(indx) == 0)	return false;
	if(m_iceSM[indx]->GetNumVertices() == 0)		return false;

	return true;
}