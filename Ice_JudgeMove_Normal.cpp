#include "Ice_JudgeMove_Normal.h"

typedef Ice_JudgeMove_Normal JudgeNormal;

JudgeNormal::Ice_JudgeMove_Normal(const vector<Ice_SM*>& iceSM)
{
	m_iceSM = iceSM;
}

JudgeNormal::~Ice_JudgeMove_Normal()
{
}

bool JudgeNormal::JudgeMove(unsigned indx)
{
	if(m_iceSM[indx]->GetNumVertices() == 0) return false;

	return true;
}