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
	float L, H, p; //domain parameters
	int m0, n0; //bump parameters
	int tstep;	//current time step
	float J, Jbar; //Temporary and long-time average of objective function
	float w[9], cx[9], cy[9]; //weights and lattice velcities
	float u0, rho0, alfa, omega; //additional parameters
	float *tSpan, **rho, **u, **v; //macroscopic quantities vectors
	float **feq[9];	//equlibrium distribution vector
	float*** u_ref[9]; //trajectory vector
	float**** f; //trajectory vectory in another representation
	
public:

	LBM();
	LBM(float,int,int,int);
	virtual ~LBM();
	
	//Set functions
	void SetGridRefinementLevel(int m_, int mult) { m = m_, n = mult * m_; }
	void SetNoTimeSteps(int mstep_) { mstep = mstep_; }
	void SetDomainParameters(float L_, float H_) { L = L_; H = H_; }
	void SetParameter(float u0_) { u0 = u0_; }

	//Get functions
	virtual float GetObjectiveFunction() { return Jbar; }
	virtual float**** GetTrajectory() { return f; }
	virtual float GetNoTimeSteps() { return mstep; }
	virtual float GetGridRefinementX() { return n; }
	virtual float GetGridRefinementY() { return m; }
	virtual int GetLatticeNo() { return 9; }
	virtual float* GetLatticeVeliocityX() { return cx; }
	virtual float* GetLatticeVeliocityY() { return cy; }
	virtual float* GetLatticeWeights() { return w; }
	virtual float GetOmega() { return omega; }
	virtual float GetStepParameter() { return p; }
	virtual float GetInletVelocity() { return u0; }
	virtual float** GetRho() { return rho; }
	virtual float** GetVelocityX() { return u; }
	virtual float** GetVelocityY() { return v; }
	virtual void GetEqulibrium(int, float**, float**, float**, float***);


	//Main loop methods
	void ApplyBC();
	void CalculateMacroscopic();
	void CollisionStep(); 
	void PostProcess();
	void StreamingStep();

	//Virtual methods
	virtual void Exectue(); //Perform simulation

	//Friend functions
	friend void SolveTimeStep(float p, float u0, int m, int n, float *cx, float *cy,
		float *w, float **rho, float **u, float **v, float omega, float ***feq,
		float*** fin, float*** fout);
	friend void dfds_b(float p, float u0, float *u0b, int m, int n, float *cx,
		float *cy, float *w, float **rho, float **u, float **v, float
		omega, float ***feq, float ***feqb, float ***fin, float ***fout,
		float ***foutb);
	friend void dfds_d(float p, float u0, float u0d, int m, int n, float *cx,
		float *cy, float *w, float **rho, float **u, float **v, float
		omega, float ***feq, float ***feqd, float ***fin, float ***fout,
		float ***foutd);
	friend void dfdu_b(float p, float u0, int m, int n, float *cx, float *cy, float *w,
		float **rho, float **u, float **v, float omega, float ***feq, float ***feqb,
		float ***fin, float ***finb, float ***fout, float ***foutb);
	friend void dfdu_d(float p, float u0, int m, int n, float *cx, float *cy, float *w,
		float **rho, float **u, float **v, float omega, float ***feq, float ***feqd,
		float ***fin, float ***find, float ***fout, float ***foutd);
	friend class LSS;


};

