#include "LSS.h"


//Constructor
LSS::LSS()
{
	m = n = Q = mstep = 0;
	solver = nullptr;	trajectory = nullptr;
	eq = nullptr;	b = nullptr;
	cx = nullptr;	cy = nullptr;	w = nullptr;
	rho = nullptr;	u = nullptr;	v = nullptr;
}

//Destructor
LSS::~LSS()
{
	delete solver, trajectory, b, rho, u, v, cx, cy, w, eq;
}

//Set case
void LSS::SetCase(double s, int t, int m, int mx)
{
	solver->SetParameter(s);
	solver->SetNoTimeSteps(t);
	solver->SetGridRefinementLevel(m, mx);
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
	for (int i = 0; i < mstep-1; i++)
		for (int j = 0; j < n; j++)
			for (int k = 0; k < m; k++)
				for (int l = 0; l < Q; l++)
				{
					result += a[i][j][k][l] * b[i][j][k][l];
				}
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
	for (int i = 0; i < mstep - 1; i++)
		for (int j = 0; j < n; j++)
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
	for (int i = 0; i < mstep - 1; i++)
		for (int j = 0; j < n; j++)
			for (int k = 0; k < m; k++)
				for (int l = 0; l < Q; l++)
					p[i][j][k][l] = r[i][j][k][l];
}

//Function for assembling array
void LSS::AssemblyArray(double ****x, double ****y)
{
	double ***dummy = nullptr, ***temp1 = new double **[n], ***temp2 = new double **[n], ***temp3 = new double **[n], ***temp4 = new double **[n], temp = 0;
	double AlfaParameter = 10;
	for (int i = 0; i < n; i++)
	{
		temp1[i] = new double *[m];	temp2[i] = new double *[m];
		temp3[i] = new double *[m];	temp4[i] = new double *[m];
		for (int j = 0; j < m; j++)
		{
			temp1[i][j] = new double[Q]; temp2[i][j] = new double[Q];
			temp3[i][j] = new double[Q]; temp4[i][j] = new double[Q];
			for (int k = 0; k < Q; k++)
			{
				temp1[i][j][k] = 0;
				temp2[i][j][k] = 0;
				temp3[i][j][k] = 0;
				temp4[i][j][k] = 0;
			}
		}
	}
	std::cout << "\nASSEMBLING ARRAY\n";
	for (int iter = 1; iter < mstep - 1; iter++)
	{
		std::cout << "#";
		//First block-row
		if (iter == 1)
		{
			//Above diagonal
			(*this).GetMacroscopic(iter); //macroscopic values for (i+0.5) time step
			dfdu_b(p, u0, m, n, cx, cy, w, rho, u, v, omega, eq, dummy, trajectory[iter], temp4, trajectory[iter + 1], x[iter]);
			dfdu_d(p, u0, m, n, cx, cy, w, rho, u, v, omega, eq, dummy, trajectory[iter], temp4, trajectory[iter + 1], temp4);

			//Calculate diagonal block
			(*this).GetMacroscopic(iter); //macroscopic values for (i+0.5) time step
			dfdu_b(p, u0, m, n, cx, cy, w, rho, u, v, omega, eq, dummy, trajectory[iter], temp1, trajectory[iter + 1], x[iter - 1]);
			dfdu_d(p, u0, m, n, cx, cy, w, rho, u, v, omega, eq, dummy, trajectory[iter], temp1, trajectory[iter + 1], temp1);
			(*this).GetMacroscopic(iter - 1); //macroscopic values for (i-0.5) time step
			dfdu_b(p, u0, m, n, cx, cy, w, rho, u, v, omega, eq, dummy, trajectory[iter - 1], temp2, trajectory[iter], x[iter - 1]);
			dfdu_d(p, u0, m, n, cx, cy, w, rho, u, v, omega, eq, dummy, trajectory[iter - 1], temp2, trajectory[iter], temp2);
			for (int i = 0; i < n; i++)
				for (int j = 0; j < m; j++)
					for (int k = 0; k < Q; k++)
						temp += (trajectory[iter][i][j][k] - trajectory[iter - 1][i][j][k])*x[iter - 1][i][j][k];
			//Calculate output
			for (int i = 0; i < n; i++)
				for (int j = 0; j < m; j++)
					for (int k = 0; k < Q; k++)
						y[iter - 1][i][j][k] = (trajectory[iter][i][j][k] - trajectory[iter - 1][i][j][k])*temp/(AlfaParameter*AlfaParameter) + temp1[i][j][k] + temp2[i][j][k] + temp4[i][j][k];

		}
		//Medium block-rows
		else if (iter < mstep - 1)
		{
			//Above diagonal
			(*this).GetMacroscopic(iter); //macroscopic values for (i+0.5) time step
			dfdu_b(p, u0, m, n, cx, cy, w, rho, u, v, omega, eq, dummy, trajectory[iter], temp4, trajectory[iter + 1], x[iter]);
			dfdu_d(p, u0, m, n, cx, cy, w, rho, u, v, omega, eq, dummy, trajectory[iter], temp4, trajectory[iter + 1], temp4);

			//Calculate diagonal block
			(*this).GetMacroscopic(iter); //macroscopic values for (i+0.5) time step
			dfdu_b(p, u0, m, n, cx, cy, w, rho, u, v, omega, eq, dummy, trajectory[iter], temp1, trajectory[iter + 1], x[iter - 1]);
			dfdu_d(p, u0, m, n, cx, cy, w, rho, u, v, omega, eq, dummy, trajectory[iter], temp1, trajectory[iter + 1], temp1);
			(*this).GetMacroscopic(iter - 1); //macroscopic values for (i-0.5) time step
			dfdu_b(p, u0, m, n, cx, cy, w, rho, u, v, omega, eq, dummy, trajectory[iter - 1], temp2, trajectory[iter], x[iter - 1]);
			dfdu_d(p, u0, m, n, cx, cy, w, rho, u, v, omega, eq, dummy, trajectory[iter - 1], temp2, trajectory[iter], temp2);
			for (int i = 0; i < n; i++)
				for (int j = 0; j < m; j++)
					for (int k = 0; k < Q; k++)
						temp += (trajectory[iter][i][j][k] - trajectory[iter - 1][i][j][k])*x[iter - 1][i][j][k];

			//Below diagonal
			(*this).GetMacroscopic(iter - 1); //macroscopic values for (i+0.5) time step
			dfdu_b(p, u0, m, n, cx, cy, w, rho, u, v, omega, eq, dummy, trajectory[iter - 1], temp3, trajectory[iter], x[iter - 2]);
			dfdu_d(p, u0, m, n, cx, cy, w, rho, u, v, omega, eq, dummy, trajectory[iter - 1], temp3, trajectory[iter], temp3);

			//Calculate output
			for (int i = 0; i < n; i++)
				for (int j = 0; j < m; j++)
					for (int k = 0; k < Q; k++)
						y[iter - 1][i][j][k] = (trajectory[iter][i][j][k] - trajectory[iter - 1][i][j][k])*temp / (AlfaParameter*AlfaParameter) + temp1[i][j][k] + temp2[i][j][k] + temp3[i][j][k] + temp4[i][j][k];


		}
		//Last block-row
		else
		{
			//Calculate diagonal block
			(*this).GetMacroscopic(iter); //macroscopic values for (i+0.5) time step
			dfdu_b(p, u0, m, n, cx, cy, w, rho, u, v, omega, eq, dummy, trajectory[iter], temp1, trajectory[iter + 1], x[iter - 1]);
			dfdu_d(p, u0, m, n, cx, cy, w, rho, u, v, omega, eq, dummy, trajectory[iter], temp1, trajectory[iter + 1], temp1);
			(*this).GetMacroscopic(iter - 1); //macroscopic values for (i-0.5) time step
			dfdu_b(p, u0, m, n, cx, cy, w, rho, u, v, omega, eq, dummy, trajectory[iter - 1], temp2, trajectory[iter], x[iter - 1]);
			dfdu_d(p, u0, m, n, cx, cy, w, rho, u, v, omega, eq, dummy, trajectory[iter - 1], temp2, trajectory[iter], temp2);
			for (int i = 0; i < n; i++)
				for (int j = 0; j < m; j++)
					for (int k = 0; k < Q; k++)
						temp += (trajectory[iter][i][j][k] - trajectory[iter - 1][i][j][k])*x[iter - 1][i][j][k];

			//Below diagonal
			(*this).GetMacroscopic(iter - 1); //macroscopic values for (i+0.5) time step
			dfdu_b(p, u0, m, n, cx, cy, w, rho, u, v, omega, eq, dummy, trajectory[iter - 1], temp3, trajectory[iter], x[iter - 2]);
			dfdu_d(p, u0, m, n, cx, cy, w, rho, u, v, omega, eq, dummy, trajectory[iter - 1], temp3, trajectory[iter], temp3);

			//Calculate output
			for (int i = 0; i < n; i++)
				for (int j = 0; j < m; j++)
					for (int k = 0; k < Q; k++)
						y[iter - 1][i][j][k] = (trajectory[iter][i][j][k] - trajectory[iter - 1][i][j][k])*temp / (AlfaParameter*AlfaParameter) + temp1[i][j][k] + temp2[i][j][k] + temp3[i][j][k];

		}
		temp = 0;
		for (int i = 0; i < n; i++)
			for (int j = 0; j < m; j++)
				for (int k = 0; k < Q; k++)
				{
					temp1[i][j][k] = 0;	temp2[i][j][k] = 0;
					temp3[i][j][k] = 0;	temp4[i][j][k] = 0;
				}
	}
	std::cout << "\nASSEMBLING DONE\n";

	delete dummy, temp1, temp2, temp3, temp4;
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
	std::cout << "\nPREPROCESSING\n";
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
	//Set equlibrium
	eq = new double **[Q];
	for (int i = 0; i < Q; i++)
	{
		eq[i] = new double *[n];
		for (int j = 0; j < n; j++)
		{
			eq[i][j] = new double[m];
			for (int k = 0; k < m; k++)
				eq[i][j][k] = 0;
		}
	}
	std::cout << "\nPREPROCESSING DONE\n";

}

void LSS::GetMacroscopic(int tstep)
{
	double ssum, usum, vsum;
	for (int j = 0; j < m; j++)
	{
		for (int i = 0; i < n; i++)
		{
			ssum = 0;
			for (int k = 0; k < 9; k++)
			{
				ssum = ssum + trajectory[tstep][i][j][k];
			}
			rho[i][j] = ssum;
		}
	}

	for (int i = 0; i < n; i++)
	{
		rho[i][m - 1] = trajectory[tstep][i][m - 1][0] + trajectory[tstep][i][m - 1][1] + trajectory[tstep][i][m - 1][3] + 2 * (trajectory[tstep][i][m - 1][2] + trajectory[tstep][i][m - 1][6] + trajectory[tstep][i][m - 1][5]);
	}

	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m - 1; j++)
		{
			usum = 0;	vsum = 0;
			for (int k = 0; k < 9; k++)
			{
				usum = usum + trajectory[tstep][i][j][k] * cx[k];
				vsum = vsum + trajectory[tstep][i][j][k] * cy[k];
			}
			u[i][j] = usum / rho[i][j];
			v[i][j] = vsum / rho[i][j];
		}
	}

	for (int j = 0; j < m; j++)
	{
		v[n - 1][j] = 0;
	}

	for (int j = 0; j < (int)(p*m); j++)
	{
		for (int i = 0; i < m; i++)
		{
			u[i][j] = 0;
			v[i][j] = 0;
		}
	}
}

void LSS::CreateRHSVector()
{
	bool flag = 0;
	std::cout << "\nCONSTRUCTING RHS VECTOR\n";
	double u0d = 1, ***dummy = nullptr, ***temp = new double **[n];
	for (int j = 0; j < n; j++)
	{
		temp[j] = new double *[m];
		for (int k = 0; k < m; k++)
		{
			temp[j][k] = new double[Q];
			for (int l = 0; l < Q; l++)
				temp[j][k][l] = 0;
		}
	}
	for (int i = 0; i < mstep - 2; i++)
	{
		(*this).GetMacroscopic(i); //macroscopic values for (i+0.5) time step
		dfds_d(p, u0, u0d, m, n, cx, cy, w, rho, u, v, omega, eq, dummy, trajectory[i], trajectory[i + 1], temp);
		(*this).GetMacroscopic(i + 1); //macroscopic values for (i+1+0.5) time step
		dfds_d(p, u0, u0d, m, n, cx, cy, w, rho, u, v, omega, eq, dummy, trajectory[i], trajectory[i + 1], b[i]);
		for (int j = 0; j < n; j++)
			for (int k = 0; k < m; k++)
				for (int l = 0; l < Q; l++)
				{
					b[i][j][k][l] = (b[i][j][k][l] + temp[j][k][l]) / 2;
					//Check if nan
					if (b[i][j][k][l] != b[i][j][k][l])
					{
						std::cout << i << "\t" << j << "\t" << k << "\t" << l << "\t"<<std::endl;
						std::cout << b[i][j][k][l] << std::endl;
						system("pause");
					}
				}

	}
	std::cout << "\RHS VECTOR CREATED\n";

	delete dummy, temp;
}

//Implementation of preconditioned conjugate gradient method for KKT system with Shur compliment
void LSS::SolveKKT()
{
	//CGM parameters
	bool ExportFlag = 0;
	std::ofstream CGM("CGM.txt");
	int N = mstep * m * n * Q; //System size
	int d = 6, D = 1;// pow(10.0, double(d / 2));
	int maxIter = pow(10.0, double(d));
	double eps = 1e-6, alfa = 0.5, beta;
	double ****res_b, ****res_b_n, ****res, ****res_n, ****p, ****temp, ****x;

	//Create solution vectors
	b = new double ***[mstep - 1];	x = new double ***[mstep - 1];
	res = new double ***[mstep - 1];	res_b = new double ***[mstep - 1];
	res_n = new double ***[mstep - 1];	res_b_n = new double ***[mstep - 1];
	p = new double ***[mstep - 1];	temp = new double ***[mstep - 1];
	for (int i = 0; i < mstep - 1; i++)
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
				b[i][j][k] = new double[Q];	x[i][j][k] = new double[Q];
				res[i][j][k] = new double[Q];	res_b[i][j][k] = new double[Q];
				res_n[i][j][k] = new double[Q];	res_b_n[i][j][k] = new double[Q];
				p[i][j][k] = new double[Q];	temp[i][j][k] = new double[Q];
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

	//Construct RHS vector
	(*this).CreateRHSVector();

	std::cout << "\nINITIALIZING RESIDUAL VECTOR\n";
	//(*this).AssemblyArray(x, res); //r = Ax
	for (int i = 0; i < mstep - 1; i++)
		for (int j = 0; j < n; j++)
			for (int k = 0; k < m; k++)
				for (int l = 0; l < Q; l++)
					res[i][j][k][l] = b[i][j][k][l] - res[i][j][k][l]; //r = b - r
	std::cout << "\RESIDUAL VECTOR CREATED\n";
	Preconditioner(res, res_b); //calculate res_b
	std::cout << "\nINITIALIZING FIRST CONJUGATE VECTOR\n";
	for (int i = 0; i < mstep - 1; i++)
		for (int j = 0; j < n; j++)
			for (int k = 0; k < m; k++)
				for (int l = 0; l < Q; l++)
				{
					p[i][j][k][l] = res_b[i][j][k][l];
				}
	std::cout << "\FIRST CONJUGATE VECTOR CREATED\n";
	int iter;
	if (ExportFlag)
	{
		CGM << "ITER\tRESIDUUM\tALFA\tBETA\n";
	}
	for (iter = 0; iter < maxIter; iter++)
	{
		(*this).AssemblyArray(p, temp); //calculate temp = Ap
		alfa = ScalarProduct(res, res_b)/ ScalarProduct(p, temp);
		
		//std::cout << "alfa = " << alfa << std::endl;
		//Export(res);
		
		for (int i = 0; i < mstep - 1; i++)
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
		for (int i = 0; i < mstep - 1; i++)
			for (int j = 0; j < n; j++)
				for (int k = 0; k < m; k++)
					for (int l = 0; l < Q; l++)
					{
						p[i][j][k][l] = res_b_n[i][j][k][l] + beta * p[i][j][k][l];
						res_b[i][j][k][l] = res_b_n[i][j][k][l];	res[i][j][k][l] = res_n[i][j][k][l];	//update res and res_b to next iter
					}
		if (iter%D == 0)	std::cout << "Residuum in " << iter + 1 << " iter = " << VectorNorm(res)
			<< "\t alfa = " << alfa << "\t beta = " << beta << std::endl;
		if (ExportFlag)
		{
			CGM << iter + 1 << "\t" << VectorNorm(res) << "\t" << alfa << "\t" << beta << "\n";
		}
	}
	std::cout << "Total number of iterations: " << iter << std::endl;
	if (ExportFlag) CGM.close();
	delete[] temp, p, res, res_n, res_b, res_b_n;
}

//Export results
void LSS::Export(double ****w)
{
	std::ofstream file("Vector.txt");
	file << "u0\tdJds\n";
	for (int i = 0; i < mstep - 1; i++)
	{
		file << "Time step no: " << i;
		for (int j = 0; j < n; j++)
		{
			for (int k = 0; k < m; k++)
			{
				for (int l = 0; l < Q; l++)
				{
					file << "\t" << w[i][j][k][l];
				}
			}
		}
		file << "\t";
	}

	file.close();
}

