#include "Ice_JudgeMove_Spears.h"

typedef Ice_JudgeMove_Spears JudgeSpears;

JudgeSpears::Ice_JudgeMove_Spears(const vector<Ice_SM*>& iceSM, IceStructure* strct)
{
	m_iceSM = iceSM;
	m_iceStrct = strct;
}

JudgeSpears::~Ice_JudgeMove_Spears()
{
}

bool JudgeSpears::JudgeMove(unsigned indx)
{
	if(m_iceStrct->GetMotionCalcCluster(indx) == 0){	return false;	}
	if(m_iceSM[indx]->GetNumVertices() == 0){			return false;	}

	return true;
}

bool JudgeSpears::JudgeMoveDebug(unsigned indx)
{
	if(m_iceSM[indx]->GetNumVertices() == 0){			return false;	}

	return true;
}