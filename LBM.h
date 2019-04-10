#pragma once
#include "SolverFactory.h"
#include "Derivatives.h"
#include <iostream>
#include <fstream>
#include "LSS.h"

class LBM : public SolverFactory
{
	int m,n; //grid refinement
	int mstep; //number of time steps
	double L, H, p; //domain parameters
	int m0, n0; //bump parameters
	int tstep;	//current time step
	double J, Jbar; //Temporary and long-time average of objective function
	double w[9], cx[9], cy[9]; //weights and lattice velcities
	double u0, rho0, alfa, omega; //additional parameters
	double *tSpan, **rho, **u, **v; //macroscopic quantities vectors
	double **feq[9];	//equlibrium distribution vector
	double*** u_ref[9]; //trajectory vector
	double**** f; //trajectory vectory in another representation
	
public:

	LBM();
	virtual ~LBM();
	
	//Set functions
	void SetGridRefinementLevel(size_t m_) { m = m_, n = 25 * m_; }
	void SetNoTimeSteps(size_t mstep_) { mstep = mstep_; }
	void SetDomainParameters(double L_, double H_) { L = L_; H = H_; }
	void SetParameter(double u0_) { u0 = u0_; }

	//Get functions
	virtual double GetObjectiveFunction() { return Jbar; }
	virtual double**** GetTrajectory() { return f; }
	virtual double GetNoTimeSteps() { return mstep; }
	virtual double GetGridRefinementX() { return n; }
	virtual double GetGridRefinementY() { return m; }
	virtual int GetLatticeNo() { return 9; }

	//Main loop methods
	void ApplyBC();
	void CalculateMacroscopic();
	void CollisionStep(); 
	void PostProcess();
	void StreamingStep();

	//Virtual methods
	virtual void Exectue(); //Perform simulation

	//Friend functions
	friend void SolveTimeStep(double p, double u0, int m, int n, double *cx, double *cy,
		double *w, double **rho, double **u, double **v, double omega, double ***feq,
		double*** fin, double*** fout);
	friend void dfds_b(double p, double u0, double *u0b, int m, int n, double *cx,
		double *cy, double *w, double **rho, double **u, double **v, double
		omega, double ***feq, double ***feqb, double ***fin, double ***fout,
		double ***foutb);
	friend void dfds_d(double p, double u0, double u0d, int m, int n, double *cx,
		double *cy, double *w, double **rho, double **u, double **v, double
		omega, double ***feq, double ***feqd, double ***fin, double ***fout,
		double ***foutd);
	friend void dfdu_b(double p, double u0, int m, int n, double *cx, double *cy, double *w,
		double **rho, double **u, double **v, double omega, double ***feq, double ***feqb,
		double ***fin, double ***finb, double ***fout, double ***foutb);
	friend void dfdu_d(double p, double u0, int m, int n, double *cx, double *cy, double *w,
		double **rho, double **u, double **v, double omega, double ***feq, double ***feqd,
		double ***fin, double ***find, double ***fout, double ***foutd);
	friend void LSS::AssemblyArray(double *x, double *y, double *d, int N);


};

