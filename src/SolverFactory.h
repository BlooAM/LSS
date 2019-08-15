#pragma once

class SolverFactory
{
	friend class LSS;
protected:
	enum SolverType { LB };
public:
	SolverFactory();
	virtual ~SolverFactory();

	virtual void Exectue() = 0; //Perform simulation
	virtual void SetParameter(double) = 0;
	virtual void SetGridRefinementLevel(int, int) = 0;
	virtual void SetNoTimeSteps(int) = 0;
	virtual double GetObjectiveFunction() = 0;
	virtual double**** GetTrajectory() = 0;
	virtual double GetNoTimeSteps() = 0;
	virtual double GetGridRefinementX() = 0;
	virtual double GetGridRefinementY() = 0;
	virtual int GetLatticeNo() = 0;
	virtual double* GetLatticeVeliocityX() = 0;
	virtual double* GetLatticeVeliocityY() = 0;
	virtual double* GetLatticeWeights() = 0;
	virtual double GetOmega() = 0;
	virtual double GetStepParameter() = 0;
	virtual double GetInletVelocity() = 0;
	virtual double** GetRho() = 0;
	virtual double** GetVelocityX() = 0;
	virtual double** GetVelocityY() = 0;
	virtual void GetEqulibrium(int, double**, double**, double**, double***) = 0;

	static SolverFactory* Create(SolverType,double,int,int,int);
};

