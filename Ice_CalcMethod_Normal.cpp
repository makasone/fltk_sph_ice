#include "Ice_CalcMethod_Normal.h"

typedef Ice_CalcMethod_Normal CalcNormal;

CalcNormal::Ice_CalcMethod_Normal(Ice_ClusterMove* clusterMove)
{
	SetObjMove(clusterMove);
}

CalcNormal::~Ice_CalcMethod_Normal()
{
}

void CalcNormal::SetObjMove(Ice_ClusterMove* clusterMove)
{
	m_iceMove = clusterMove;
}

void CalcNormal::StepObjMove()
{
	m_iceMove->StepObjMove();		//‚»‚Ì‚Ü‚ÜŒÄ‚Ô‚¾‚¯
}