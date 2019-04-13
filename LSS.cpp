#include "LSS.h"


//Constructor
LSS::LSS()
{
	m = n = Q = mstep = 0;
	solver = nullptr;
	trajectory = nullptr;
	x = nullptr;
	b = nullptr;
}

//Destructor
LSS::~LSS()
{
	delete solver, trajectory, b, x;
}

//Calculate scalar product of two vectors
double LSS::ScalProduct(int N, double*a, double*b)
{
	double result = 0;
	for (int i = 0; i < N; i++)	result += a[i] * b[i];
	return result;
}
double LSS::ScalarProduct(double**a, double**b)
{
	double result = 0;
	for (int i = 0; i < mstep; i++)	
		for(int j = 0; j < m*n*Q; j++)
			result += a[i][j] * b[i][j];
	return result;
}

//Calculate vector norm
double LSS::Norm(double *r, int N)
{
	double res = 0;
	for (int i = 0; i < N; i++)
	{
		res += r[i] * r[i];
	}
	return sqrt(res);
}
double LSS::VectorNorm(double **r)
{
	double res = 0;
	for (int i = 0; i < mstep; i++)
		for(int j=0; j<m*n*Q; j++)
			res += r[i][j] * r[i][j];
	return sqrt(res);
}

//Preconditioner(diagonal)
void LSS::Precond(int N, double *diagA, double*r, double*p)
{
	double temp;
	for (int i = 0; i < N; i++)
	{
		temp = 0;
		//Gauss-Seidel
		//for (int j = 0; j < i; j++) temp += A[i][j] * p[j];
		p[i] = (r[i] - temp) / diagA[i];
	}
}

//Dummy preconditioner
void LSS::Preconditioner(double**r, double**p)
{
	double temp;
	for (int i = 0; i < mstep; i++)
		for(int j=0; j<m*n*Q; j++)
			p[i][j] = r[i][j];
}

//Function for assembling array
void LSS::AssemblyArray(double *x, double *y, double *d, int N)
{
	if (N != mstep * m * n * Q)
	{
		std::cerr << "System size is invalid";
		exit(-1);
	}
	for (int tsteps = 0; tsteps < mstep; tsteps++)
	{
		//First block-row
		if (tsteps == 0)
		{
			
		}
		//Medium block-rows
		if (tsteps < mstep-1)
		{

		}
		//Last block-row
		else
		{

		}
	}
}

//Implementation of preconditioned conjugate gradient method
void LSS::Solve(void(*mult)(double*, double*, double*, int), double*b, double*x)
{
	//Set trajectory
	trajectory = (*this).GetTrajectory();

	//Set auxiliary parameters
	mstep = solver->GetNoTimeSteps();
	m = solver->GetGridRefinementY();
	n = solver->GetGridRefinementX();
	Q = solver->GetLatticeNo();

	//Averging trajectory
	for (int tsteps = 0; tsteps < mstep; tsteps++)
	{
		for (int i=0; i<n; i++)
		{
			for (int j=0; j<m; j++)
			{
				for (int k = 0; k < Q; k++)
				{
					if(tsteps < mstep-1)	trajectory[tsteps][i][j][k] = (trajectory[tsteps + 1][i][j][k] + trajectory[tsteps][i][j][k]) / 2;
					else trajectory[tsteps][i][j][k] = 0;
				}
			}
		}
	}

	//CGM parameters
	int N = mstep * m * n * Q; //System size
	int d = 6, D = pow(10.0, double(d / 2));
	int maxIter = pow(10.0, double(d));
	double eps = 1e-6, alfa = 0.5, beta;
	double *temp = new double[N], *p = new double[N], *diagA = new double[N];
	double *res_b = new double[N], *res_b_n = new double[N], *res = new double[N], *res_n = new double[N];
	
	mult(x, res, diagA, N); //r = Ax
	for (int i = 0; i < N; i++) res[i] = b[i] - res[i]; //r = b - r
	Precond(N, diagA, res, res_b); //calculate res_b
	for (int i = 0; i < N; i++) p[i] = res_b[i];
	int iter;
	for (iter = 0; iter < maxIter; iter++)
	{
		mult(p, temp, diagA, N); //calculate temp = Ap
		alfa = ScalProduct(N, res, res_b) / ScalProduct(N, p, temp);
		for (int i = 0; i < n; i++)
		{
			x[i] = x[i] + alfa * p[i]; //updated value of x
			res_n[i] = res[i] - alfa * temp[i];
		}
		if (Norm(res_n, N) < eps)	break;
		Precond(N, diagA, res_n, res_b_n);
		beta = ScalProduct(N, res_b_n, res_n) / ScalProduct(N, res_b, res);
		for (int i = 0; i < N; i++)
		{
			p[i] = res_b_n[i] + beta * p[i];
			res_b[i] = res_b_n[i];	res[i] = res_n[i];	//update res and z to next iter
		}
		if (iter%D == 0)	std::cout << "Residuum in " << iter + 1 << " iter = " << Norm(res, N)
			<< "\t alfa = " << alfa << "\t beta = " << beta << std::endl;
	}
	std::cout << "Total number of iterations: " << iter << std::endl;
	delete[] temp, p, diagA, res, res_n, res_b, res_b_n;
}

//Implementation of preconditioned conjugate gradient method for KKT system with Shur compliment
void LSS::SolveKKT(void(*mult)(double**, double**, int))
{
	//CGM parameters
	int N = mstep * m * n * Q; //System size
	int d = 6, D = pow(10.0, double(d / 2));
	int maxIter = pow(10.0, double(d));
	double eps = 1e-6, alfa = 0.5, beta;
	double **res_b, **res_b_n, **res, **res_n, **p, **temp;

	//Set trajectory
	trajectory = (*this).GetTrajectory();

	//Set auxiliary parameters
	mstep = solver->GetNoTimeSteps();
	m = solver->GetGridRefinementY();
	n = solver->GetGridRefinementX();
	Q = solver->GetLatticeNo();

	//Averging trajectory
	for (int tsteps = 0; tsteps < mstep; tsteps++)
	{
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < m; j++)
			{
				for (int k = 0; k < Q; k++)
				{
					if (tsteps < mstep - 1)	trajectory[tsteps][i][j][k] = (trajectory[tsteps + 1][i][j][k] + trajectory[tsteps][i][j][k]) / 2;
					else trajectory[tsteps][i][j][k] = 0;
				}
			}
		}
	}

	//Create solution and RHS vectors
	b = new double *[mstep];	x = new double *[mstep];
	res = new double *[mstep];	res_b = new double *[mstep];
	res_n = new double *[mstep];	res_b_n = new double *[mstep];
	p = new double *[mstep];	temp = new double *[mstep];
	for (int i = 0; i < mstep; i++)
	{
		b[i] = new double[m*n*Q];	x[i] = new double[m*n*Q];
		res[i] = new double[m*n*Q];	res_b[i] = new double[m*n*Q];
		res_n[i] = new double[m*n*Q];	res_b_n[i] = new double[m*n*Q];
		p[i] = new double[m*n*Q];	temp[i] = new double[m*n*Q];
		for (int j = 0; j < m*n*Q; j++)
		{
			b[i][j] = 0;	x[i][j] = 0;
			res[i][j] = 0;	res_b[i][j] = 0;
			res_n[i][j] = 0;	res_b_n[i][j] = 0;
			p[i][j] = 0;	temp[i][j] = 0;
		}
	}


	mult(x, res, N); //r = Ax
	for (int i=0; i< mstep; i++)
		for (int j = 0; j < m*n*Q; j++) 
			res[i][j] = b[i][j] - res[i][j]; //r = b - r
	Preconditioner(res, res_b); //calculate res_b
	for(int i = 0; i<mstep; i++)
		for (int j = 0; j < m*n*Q; j++) 
			p[i][j] = res_b[i][j];
	int iter;
	for (iter = 0; iter < maxIter; iter++)
	{
		mult(p, temp, N); //calculate temp = Ap
		alfa = ScalarProduct(res, res_b) / ScalarProduct(p, temp);
		for(int i = 0; i<mstep; i++)
			for (int j = 0; j < m*n*Q; j++)
			{
				x[i][j] = x[i][j] + alfa * p[i][j]; //updated value of x
				res_n[i][j] = res[i][j] - alfa * temp[i][j];
			}
		if (VectorNorm(res_n) < eps)	break;
		Preconditioner(res_n, res_b_n);
		beta = ScalarProduct(res_b_n, res_n) / ScalarProduct(res_b, res);
		for (int i = 0; i < mstep; i++)
			for(int j = 0; j < m*n*Q; j++)
			{
				p[i][j] = res_b_n[i][j] + beta * p[i][j];
				res_b[i][j] = res_b_n[i][j];	res[i][j] = res_n[i][j];	//update res and z to next iter
			}
		if (iter%D == 0)	std::cout << "Residuum in " << iter + 1 << " iter = " << VectorNorm(res)
			<< "\t alfa = " << alfa << "\t beta = " << beta << std::endl;
	}
	std::cout << "Total number of iterations: " << iter << std::endl;
	delete[] temp, p, res, res_n, res_b, res_b_n;
}



