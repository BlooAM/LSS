#include "LSS.h"


//Constructor
LSS::LSS()
{
	m = n = Q = mstep = 0;
	solver = nullptr;	trajectory = nullptr;
	x = nullptr;	b = nullptr;
	cx = nullptr;	cy = nullptr;	w = nullptr;
	rho = nullptr;	u = nullptr;	v = nullptr;
	eq = nullptr;
}

//Destructor
LSS::~LSS()
{
	delete solver, trajectory, b, x, rho, u, v, cx, cy, w, eq;
}

//Calculate scalar product of two vectors
double LSS::ScalProduct(int N, double*a, double*b)
{
	double result = 0;
	for (int i = 0; i < N; i++)	result += a[i] * b[i];
	return result;
}
double LSS::ScalarProduct(double****a, double****b)
{
	double result = 0;
	for (int i = 0; i < mstep; i++)	
		for(int j = 0; j < n; j++)
			for (int k = 0; k < m; k++)
				for (int l = 0; l < Q; l++)
					result += a[i][j][k][l] * b[i][j][k][l];
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
double LSS::VectorNorm(double ****r)
{
	double res = 0;
	for (int i = 0; i < mstep; i++)
		for(int j=0; j<n; j++)
			for (int k = 0; k < m; k++)
				for (int l = 0; l < Q; l++)
					res += r[i][j][k][l] * r[i][j][k][l];
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
void LSS::Preconditioner(double****r, double****p)
{
	double temp;
	for (int i = 0; i < mstep; i++)
		for(int j=0; j<n; j++)
			for (int k = 0; k < m; k++)
				for (int l = 0; l < Q; l++)
					p[i][j][k][l] = r[i][j][k][l];
}

//Function for assembling array
void LSS::AssemblyArray(double ****x, double ****y)
{
	double ***dummy = nullptr, ***temp = new double **[m];
	for (int i = 0; i < n; i++)
		for (int j = 0; j < m; j++)
			for (int k = 0; k < Q; k++)
				temp[i][j][k] = 0;

	for (int i = 0; i < mstep-1; i++)
	{
		//First block-row
		if (i == 0)
		{
			
		}
		//Medium block-rows
		if (i < mstep-2)
		{
			solver->GetMacroscopic(i, rho, u, v); //macroscopic values for (i+0.5) time step
			dfdu_b(p, u0, m, n, cx, cy, w, rho, u, v, omega, eq, dummy, trajectory[i], temp, trajectory[i + 1], x[i]);
		}
		//Last block-row
		else
		{

		}
	}
	delete dummy, temp;
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

void LSS::Preprocess()
{
	//Set trajectory
	trajectory = (*this).GetTrajectory();

	//Set auxiliary parameters
	mstep = solver->GetNoTimeSteps();
	m = solver->GetGridRefinementY();
	n = solver->GetGridRefinementX();
	Q = solver->GetLatticeNo();
	p = solver->GetStepParameter();
	omega = solver->GetOmega();
	u0 = solver->GetInletVelocity();
	cx = solver->GetLatticeVeliocityX();
	cy = solver->GetLatticeVeliocityY();
	w = solver->GetLatticeWeights();


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

	//Set macroscopic values
	rho = new double *[n];	u = new double *[n];	v = new double *[n];
	for (int j = 0; j < n; ++j)
	{
		rho[j] = new double[m];	u[j] = new double[m];	v[j] = new double[m];
		for (int k = 0; k < m; ++k)
		{
			rho[j][k] = 0;	u[j][k] = 0;	v[j][k] = 0;
		}
	}
}

//Implementation of preconditioned conjugate gradient method for KKT system with Shur compliment
void LSS::SolveKKT(void(*mult)(double****, double****))
{
	//CGM parameters
	int N = mstep * m * n * Q; //System size
	int d = 6, D = pow(10.0, double(d / 2));
	int maxIter = pow(10.0, double(d));
	double eps = 1e-6, alfa = 0.5, beta;
	double ****res_b, ****res_b_n, ****res, ****res_n, ****p, ****temp;

	//Create solution and RHS vectors
	b = new double ***[mstep];	x = new double ***[mstep];
	res = new double ***[mstep];	res_b = new double ***[mstep];
	res_n = new double ***[mstep];	res_b_n = new double ***[mstep];
	p = new double ***[mstep];	temp = new double ***[mstep];
	for (int i = 0; i < mstep; i++)
	{
		b[i] = new double **[n];	x[i] = new double **[n];
		res[i] = new double **[n];	res_b[i] = new double **[n];
		res_n[i] = new double **[n];	res_b_n[i] = new double **[n];
		p[i] = new double **[n];	temp[i] = new double **[n];
		for (int j = 0; j < n; j++)
		{
			b[i][j] = new double *[m];	x[i][j] = new double *[m];
			res[i][j] = new double *[m];	res_b[i][j] = new double *[m];
			res_n[i][j] = new double *[m];	res_b_n[i][j] = new double *[m];
			p[i][j] = new double *[m];	temp[i][j] = new double *[m];
			for (int k = 0; k < m; k++)
			{
				b[i][j][k] = new double [Q];	x[i][j][k] = new double [Q];
				res[i][j][k] = new double [Q];	res_b[i][j][k] = new double [Q];
				res_n[i][j][k] = new double [Q];	res_b_n[i][j][k] = new double [Q];
				p[i][j][k] = new double [Q];	temp[i][j][k] = new double [Q];
				for (int l = 0; l < Q; l++)
				{
					b[i][j][k][l] = 0;	x[i][j][k][l] = 0;
					res[i][j][k][l] = 0;	res_b[i][j][k][l] = 0;
					res_n[i][j][k][l] = 0;	res_b_n[i][j][k][l] = 0;
					p[i][j][k][l] = 0;	temp[i][j][k][l] = 0;
				}
			}
		}
	}


	mult(x, res); //r = Ax
	for (int i=0; i< mstep; i++)
		for (int j = 0; j < n; j++)
			for (int k = 0; k < m; k++)
				for (int l = 0; l < Q; l++)
					res[i][j][k][l] = b[i][j][k][l] - res[i][j][k][l]; //r = b - r
	Preconditioner(res, res_b); //calculate res_b
	for(int i = 0; i<mstep; i++)
		for (int j = 0; j < n; j++) 
			for (int k = 0; k < m; k++)
				for (int l = 0; l < Q; l++)
					p[i][j][k][l] = res_b[i][j][k][l];
	int iter;
	for (iter = 0; iter < maxIter; iter++)
	{
		mult(p, temp); //calculate temp = Ap
		alfa = ScalarProduct(res, res_b) / ScalarProduct(p, temp);
		for(int i = 0; i<mstep; i++)
			for (int j = 0; j < n; j++)
				for (int k = 0; k < m; k++)
					for (int l = 0; l < Q; l++)
					{
						x[i][j][k][l] = x[i][j][k][l] + alfa * p[i][j][k][l]; //updated value of x
						res_n[i][j][k][l] = res[i][j][k][l] - alfa * temp[i][j][k][l];
					}
		if (VectorNorm(res_n) < eps)	break;
		Preconditioner(res_n, res_b_n);
		beta = ScalarProduct(res_b_n, res_n) / ScalarProduct(res_b, res);
		for (int i = 0; i < mstep; i++)
			for(int j = 0; j < n; j++)
				for (int k = 0; k < m; k++)
					for (int l = 0; l < Q; l++)
					{
						p[i][j][k][l] = res_b_n[i][j][k][l] + beta * p[i][j][k][l];
						res_b[i][j][k][l] = res_b_n[i][j][k][l];	res[i][j][k][l] = res_n[i][j][k][l];	//update res and res_b to next iter
					}
		if (iter%D == 0)	std::cout << "Residuum in " << iter + 1 << " iter = " << VectorNorm(res)
			<< "\t alfa = " << alfa << "\t beta = " << beta << std::endl;
	}
	std::cout << "Total number of iterations: " << iter << std::endl;
	delete[] temp, p, res, res_n, res_b, res_b_n;
}



