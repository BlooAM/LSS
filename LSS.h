#pragma once
#include "SolverFactory.h"

class LSS
{
	size_t m, n;
	SolverFactory *ref_trajectory;
	
public:
	LSS();
	~LSS();

	void CreateCase() { ref_trajectory = SolverFactory::Create(SolverFactory::LB); } //One solver available
	void SetCase(double s) { ref_trajectory->SetParameter(s); }
	void SolveCase() { ref_trajectory->Exectue(); }
	double GetObjectiveFuntion() { return ref_trajectory->GetObjectiveFunction(); }
};

