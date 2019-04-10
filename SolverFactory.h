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
	virtual double GetObjectiveFunction() = 0;
	virtual double**** GetTrajectory() = 0;
	virtual double GetNoTimeSteps() = 0;
	virtual double GetGridRefinementX() = 0;
	virtual double GetGridRefinementY() = 0;
	virtual int GetLatticeNo() = 0;
	static SolverFactory* Create(SolverType);
};

