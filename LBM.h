#pragma once
#include "SolverFactory.h"
#include <iostream>
#include <fstream>

class LBM : public SolverFactory
{
	size_t m,n; //grid refinement
	size_t mstep; //number of time steps
	double L, H, p; //domain parameters
	size_t m0, n0; //bump parameters
	int tstep;	//current time step
	double w[9], cx[9], cy[9]; //weights and lattice velcities
	double u0, rho0, alfa, omega;
	double *tSpan, **rho, **u, **v; //macroscopic quantities vectors
	double **feq[9];	//equlibrium distribution vector
	double*** u_ref[9]; //trajectory vector
	
public:

	LBM();
	virtual ~LBM();
	
	//Set functions
	void SetGridRefinementLevel(size_t m_) { m = m_, n = 25 * m_; }
	void SetNoTimeSteps(size_t mstep_) { mstep = mstep_; }
	void SetDomainParameters(double L_, double H_, double p_) { L = L_; H = H_; p = p_; }

	//Main loop methods
	void ApplyBC();
	void CalculateMacroscopic();
	void CollisionStep(); //Involves applying BC
	virtual void Exectue();
	void PostProcess();
	void StreamingStep();
	
	



};

