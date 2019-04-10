#pragma once
#include "SolverFactory.h"
#include <cmath>
#include <iostream>
#include <iomanip>

class LSS
{
	int m, n, Q, mstep;
	SolverFactory *solver;
	double **** trajectory;
	
public:
	LSS();
	~LSS();

	//Baseline solvers
	void CreateCase() { solver = SolverFactory::Create(SolverFactory::LB); } //One solver available
	void SetCase(double s) { solver->SetParameter(s); }
	void SolveCase() { solver->Exectue(); }
	double GetObjectiveFuntion() { return solver->GetObjectiveFunction(); }
	double**** GetTrajectory() { return solver->GetTrajectory(); }

	//Auxiliary functions
	double Norm(double *r, int n);
	double ScalProduct(int, double*, double*);
	
	//LSS functions
	void AssemblyArray(double*, double*, double*, int);
	void Precond(int, double *, double*, double*);
	void SolveKKT(void(*mult)(double*, double*, double*, int), double*, double*);
};

