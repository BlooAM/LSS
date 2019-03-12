#include "SolverFactory.h"
#include "LBM.h"

SolverFactory::SolverFactory()
{

}

SolverFactory::~SolverFactory()
{

}

SolverFactory* SolverFactory::Create(SolverType type)
{
	if (type == LB)
		return new LBM();
}
