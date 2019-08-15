#include "SolverFactory.h"
#include "LBM.h"

SolverFactory::SolverFactory()
{

}

SolverFactory::~SolverFactory()
{

}

SolverFactory* SolverFactory::Create(SolverType type, double u0_, int tsteps_, int m_, int mx)
{
	if (type == LB)
		return new LBM(u0_,tsteps_,m_,mx);
}
