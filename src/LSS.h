#pragma once
#include "SolverFactory.h"
#include "LBM.h"
#include <cmath>
#include <iostream>
#include <iomanip>

class LSS
{
	int m, n, Q, mstep;
	SolverFactory *solver;
	double **** trajectory;
	double ****b;
	double ***eq;
	double **rho, **u, **v;
	double *cx, *cy, *w;
	double p, omega, u0;
	
public:
	LSS();
	~LSS();

	//Baseline solvers
	void CreateCase(double u0_, int tsteps_, int m_, int mx) { solver = SolverFactory::Create(SolverFactory::LB,u0_,tsteps_,m_,mx); } //One solver available
	void SetCase(double s, int t, int m, int mx);
	void SolveCase() { solver->Exectue(); }
	double GetObjectiveFuntion() { return solver->GetObjectiveFunction(); }
	double**** GetTrajectory() { return solver->GetTrajectory(); }

	//Auxiliary functions
	double Norm(double *, int);
	double VectorNorm(double ****);
	double ScalProduct(int, double*, double*);
	double ScalarProduct(double****, double****);

	//LSS functions
	void AssemblyArray(double****,double****);
	void CreateRHSVector();
	void GetMacroscopic(int);
	void Precond(int, double *, double*, double*);
	void Preconditioner(double****, double****);
	void Preprocess();
	void Solve(void(*mult)(double*, double*, double*, int),double*,double*);
	void SolveKKT();

	//Export functions
	void Export(double****);
};

