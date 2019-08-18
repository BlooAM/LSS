#pragma once
#include "SolverFactory.h"
#include "LBM.h"
#include <cmath>
#include <iostream>
#include <iomanip>
#include <chrono>


class LSS
{
	int m, n, Q, mstep;
	SolverFactory *solver;
	float **** trajectory;
	float ****b;
	float ***eq;
	float **rho, **u, **v;
	float *cx, *cy, *w;
	float p, omega, u0;
	
public:
	LSS();
	~LSS();

	//Friend functions
	friend void Obj(float p, float u0, int m, int n, float *cx, float *cy, float *w, float **rho, float **u, float **v, float omega, float ***feq, float*** fin, float*** fout, float *J);

	//Baseline solvers
	void CreateCase(float u0_, int tsteps_, int m_, int mx) { solver = SolverFactory::Create(SolverFactory::LB,u0_,tsteps_,m_,mx); } //One solver available
	void SetCase(float s, int t, int m, int mx);
	void SolveCase() { solver->Exectue(); }
	float GetObjectiveFuntion() { return solver->GetObjectiveFunction(); }
	float**** GetTrajectory() { return solver->GetTrajectory(); }

	//Auxiliary functions
	float Norm(float *, int);
	float VectorNorm(float ****);
	float ScalProduct(int, float*, float*);
	float ScalarProduct(float****, float****);

	//LSS functions
	void AssemblyArray(float****,float****);
	float CalculateSensitivity(float****);
	void CreateRHSVector();
	void GetMacroscopic(int);
	void Precond(int, float *, float*, float*);
	void Preconditioner(float****, float****);
	void Preprocess();
	void Solve(void(*mult)(float*, float*, float*, int),float*,float*);
	void SolveKKT();

	//Export functions
	void Export(float****);
};

