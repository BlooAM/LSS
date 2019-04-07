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
	static SolverFactory* Create(SolverType);
};

