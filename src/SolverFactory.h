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
	virtual void SetParameter(float) = 0;
	virtual void SetGridRefinementLevel(int, int) = 0;
	virtual void SetNoTimeSteps(int) = 0;
	virtual float GetObjectiveFunction() = 0;
	virtual float**** GetTrajectory() = 0;
	virtual float GetNoTimeSteps() = 0;
	virtual float GetGridRefinementX() = 0;
	virtual float GetGridRefinementY() = 0;
	virtual int GetLatticeNo() = 0;
	virtual float* GetLatticeVeliocityX() = 0;
	virtual float* GetLatticeVeliocityY() = 0;
	virtual float* GetLatticeWeights() = 0;
	virtual float GetOmega() = 0;
	virtual float GetStepParameter() = 0;
	virtual float GetInletVelocity() = 0;
	virtual float** GetRho() = 0;
	virtual float** GetVelocityX() = 0;
	virtual float** GetVelocityY() = 0;
	virtual void GetEqulibrium(int, float**, float**, float**, float***) = 0;

	static SolverFactory* Create(SolverType,float,int,int,int);
};

