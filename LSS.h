#pragma once
#include "SolverFactory.h"

class LSS
{
	size_t m, n;
	SolverFactory *ref_trajectory;
	
public:
	LSS();
	~LSS();

	void CreateCase() { ref_trajectory = SolverFactory::Create(SolverFactory::LB); }
	void SolveCase() { ref_trajectory->Exectue(); }

};

