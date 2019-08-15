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

void Obj(double p, double u0, int m, int n, double *cx, double *cy, double *w, double **rho, double **u, double **v, double omega, double ***feq, double*** fin, double*** fout, double *J)
{
	//Buff for fout
	double*** foutBuffer = new double**[n];
	for (int i = 0; i < n; ++i)
	{
		foutBuffer[i] = new double*[m];

		for (int j = 0; j < m; ++j)
		{
			foutBuffer[i][j] = new double[9];

			for (int k = 0; k < 9; ++k)
			{
				foutBuffer[i][j][k] = fout[i][j][k];
			}
		}
	}

	double pin = 0, pout = 0;

	//Time step
	int m0 = p * m, n0 = m;

	//Collision step
	double temp1, temp2, rhow, ssum, usum, vsum;
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++)
		{
			temp1 = u[i][j] * u[i][j] + v[i][j] * v[i][j];
			for (int k = 0; k < 9; k++)
			{
				temp2 = u[i][j] * cx[k] + v[i][j] * cy[k];
				feq[k][i][j] = rho[i][j] * w[k] * (1 + 3 * temp2 + 4.5*temp2*temp2 - 1.5*temp1);
				fout[i][j][k] = omega * feq[k][i][j] + (1 - omega)*fin[i][j][k];
			}
		}
	}

	//Apply BC
	//Inlet
	for (int j = m0; j < m; j++)
	{
		rhow = (fout[0][j][0] + fout[0][j][2] + fout[0][j][4] + 2 * (fout[0][j][3] + fout[0][j][6] + fout[0][j][7])) / (1 - u0);
		fout[0][j][1] = fout[0][j][3] + 2 * rhow*u0 / 3;
		fout[0][j][5] = fout[0][j][7] + rhow * u0 / 6;
		fout[0][j][8] = fout[0][j][6] + rhow * u0 / 6;
	}
	//South
	for (int i = 0; i < n; i++)
	{
		fout[i][0][2] = fout[i][0][4];
		fout[i][0][5] = fout[i][0][7];
		fout[i][0][6] = fout[i][0][8];
	}
	//North
	for (int i = 0; i < n; i++)
	{
		fout[i][m - 1][4] = fout[i][m - 1][2];
		fout[i][m - 1][8] = fout[i][m - 1][6];
		fout[i][m - 1][7] = fout[i][m - 1][5];
	}
	//Outlet - extrapolation
	for (int j = 0; j < m; j++)
	{
		fout[n - 1][j][1] = 2 * fout[n - 2][j][1] - fout[n - 3][j][1];
		fout[n - 1][j][5] = 2 * fout[n - 2][j][5] - fout[n - 3][j][5];
		fout[n - 1][j][8] = 2 * fout[n - 2][j][8] - fout[n - 3][j][8];
	}
	//Back-facing step
	for (int i = 0; i < n0; i++)
	{
		fout[i][m0 - 1][2] = fout[i][m0 - 1][4];
		fout[i][m0 - 1][5] = fout[i][m0 - 1][7];
		fout[i][m0 - 1][6] = fout[i][m0 - 1][8];
	}
	for (int j = 0; j < m0; j++)
	{
		fout[n0 - 1][j][1] = fout[n0 - 1][j][3];
		fout[n0 - 1][j][5] = fout[n0 - 1][j][7];
		fout[n0 - 1][j][8] = fout[n0 - 1][j][6]; //???
	}

	//Streaming step
	for (int j = 0; j < m; j++)
	{
		for (int i = n - 1; i > 0; i--)
		{
			fout[i][j][1] = fout[i - 1][j][1];
		}

		for (int i = 0; i < n - 1; i++)
		{
			fout[i][j][3] = fout[i + 1][j][3];
		}
	}

	for (int j = m - 1; j > 0; j--)
	{
		for (int i = 0; i < n; i++)
		{
			fout[i][j][2] = fout[i][j - 1][2];
		}
		for (int i = n - 1; i > 0; i--)
		{
			fout[i][j][5] = fout[i - 1][j - 1][5];
		}
		for (int i = 0; i < n - 1; i++)
		{
			fout[i][j][6] = fout[i + 1][j - 1][6];
		}
	}
	for (int j = 0; j < m - 1; j++)
	{
		for (int i = 0; i < n; i++)
		{
			fout[i][j][4] = fout[i][j + 1][4];
		}
		for (int i = 0; i < n - 1; i++)
		{
			fout[i][j][7] = fout[i + 1][j + 1][7];
		}
		for (int i = n - 1; i > 0; i--)
		{
			fout[i][j][8] = fout[i - 1][j + 1][8];
		}
	}

	//Calculate macroscopic
	for (int j = 0; j < m; j++)
	{
		for (int i = 0; i < n; i++)
		{
			ssum = 0;
			for (int k = 0; k < 9; k++)
			{
				ssum = ssum + fout[i][j][k];
			}
			rho[i][j] = ssum;
		}
	}

	for (int i = 0; i < n; i++)
	{
		rho[i][m - 1] = fout[i][m - 1][0] + fout[i][m - 1][1] + fout[i][m - 1][3] + 2 * (fout[i][m - 1][2] + fout[i][m - 1][6] + fout[i][m - 1][5]);
	}

	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m - 1; j++)
		{
			usum = 0;	vsum = 0;
			for (int k = 0; k < 9; k++)
			{
				usum = usum + fout[i][j][k] * cx[k];
				vsum = vsum + fout[i][j][k] * cy[k];
			}
			u[i][j] = usum / rho[i][j];
			v[i][j] = vsum / rho[i][j];
		}
	}

	for (int j = 0; j < m; j++)
	{
		v[n - 1][j] = 0;
	}

	for (int j = 0; j < m0; j++)
	{
		for (int i = 0; i < n0; i++)
		{
			u[i][j] = 0;
			v[i][j] = 0;
		}
	}

	//Pressure drop
	for (int j = m0; j < m; j++)	pin += rho[0][j] / 3;
	for (int j = 0; j < m; j++)	pout += rho[n - 1][j] / 3;
	pin /= (m - m0);
	pout /= m;
	(*J) = (pin - pout) / 2;


	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < m; ++j)
		{
			for (int k = 0; k < 9; ++k)
			{
				fout[i][j][k] = foutBuffer[i][j][k];
			}
		}
	}

	for (int iter = 0; iter < n; ++iter)
	{
		for (int j = 0; j < m; ++j)
		{
			delete[] foutBuffer[iter][j];
		}
		delete[] foutBuffer[iter];
	}
	delete[] foutBuffer;
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
	for (int iter = 0; iter < mstep - 1; iter++)
	{
		std::cout << "#";
		//First block-row
		if (iter == 0)
		{
			//Above diagonal
			(*this).GetMacroscopic(iter + 1); //macroscopic values for (i+0.5) time step
			dfdu_b(p, u0, m, n, cx, cy, w, rho, u, v, omega, eq, dummy, trajectory[iter+1], temp4, trajectory[iter + 2], x[iter + 1]);
			(*this).GetMacroscopic(iter);
			dfdu_d(p, u0, m, n, cx, cy, w, rho, u, v, omega, eq, dummy, trajectory[iter], temp4, trajectory[iter + 1], temp4);

			//Calculate diagonal block
			(*this).GetMacroscopic(iter); //macroscopic values for (i+0.5) time step
			dfdu_b(p, u0, m, n, cx, cy, w, rho, u, v, omega, eq, dummy, trajectory[iter], temp1, trajectory[iter + 1], x[iter]);
			dfdu_d(p, u0, m, n, cx, cy, w, rho, u, v, omega, eq, dummy, trajectory[iter], temp1, trajectory[iter + 1], temp1);
			(*this).GetMacroscopic(iter); //macroscopic values for (i-0.5) time step
			dfdu_b(p, u0, m, n, cx, cy, w, rho, u, v, omega, eq, dummy, trajectory[iter], temp2, trajectory[iter + 1], x[iter]);
			dfdu_d(p, u0, m, n, cx, cy, w, rho, u, v, omega, eq, dummy, trajectory[iter], temp2, trajectory[iter + 1], temp2);
			for (int i = 0; i < n; i++)
				for (int j = 0; j < m; j++)
					for (int k = 0; k < Q; k++)
						temp += (trajectory[iter + 1][i][j][k] - trajectory[iter][i][j][k])*x[iter][i][j][k];
			//Calculate output
			for (int i = 0; i < n; i++)
				for (int j = 0; j < m; j++)
					for (int k = 0; k < Q; k++)
						y[iter][i][j][k] = (trajectory[iter + 1][i][j][k] - trajectory[iter][i][j][k])*temp/(AlfaParameter*AlfaParameter) + temp1[i][j][k] + temp2[i][j][k] + temp4[i][j][k];

		}
		//Medium block-rows
		else if (iter < mstep - 2)
		{
			//Above diagonal
			(*this).GetMacroscopic(iter + 1); //macroscopic values for (i+1.5) time step
			dfdu_b(p, u0, m, n, cx, cy, w, rho, u, v, omega, eq, dummy, trajectory[iter + 1], temp4, trajectory[iter + 2], x[iter + 1]);
			(*this).GetMacroscopic(iter);
			dfdu_d(p, u0, m, n, cx, cy, w, rho, u, v, omega, eq, dummy, trajectory[iter], temp4, trajectory[iter + 1], temp4);

			//Calculate diagonal block
			(*this).GetMacroscopic(iter); //macroscopic values for (i+0.5) time step
			dfdu_b(p, u0, m, n, cx, cy, w, rho, u, v, omega, eq, dummy, trajectory[iter], temp1, trajectory[iter + 1], x[iter]);
			dfdu_d(p, u0, m, n, cx, cy, w, rho, u, v, omega, eq, dummy, trajectory[iter], temp1, trajectory[iter + 1], temp1);
			(*this).GetMacroscopic(iter); //macroscopic values for (i-0.5) time step
			dfdu_b(p, u0, m, n, cx, cy, w, rho, u, v, omega, eq, dummy, trajectory[iter], temp2, trajectory[iter + 1], x[iter]);
			dfdu_d(p, u0, m, n, cx, cy, w, rho, u, v, omega, eq, dummy, trajectory[iter], temp2, trajectory[iter + 1], temp2);
			for (int i = 0; i < n; i++)
				for (int j = 0; j < m; j++)
					for (int k = 0; k < Q; k++)
						temp += (trajectory[iter][i][j][k] - trajectory[iter - 1][i][j][k])*x[iter][i][j][k];

			//Below diagonal
			(*this).GetMacroscopic(iter - 1); //macroscopic values for (i+0.5) time step
			dfdu_b(p, u0, m, n, cx, cy, w, rho, u, v, omega, eq, dummy, trajectory[iter - 1], temp3, trajectory[iter], x[iter - 1]);
			(*this).GetMacroscopic(iter); //macroscopic values for (i+0.5) time step
			dfdu_d(p, u0, m, n, cx, cy, w, rho, u, v, omega, eq, dummy, trajectory[iter], temp3, trajectory[iter + 1], temp3);

			//Calculate output
			for (int i = 0; i < n; i++)
				for (int j = 0; j < m; j++)
					for (int k = 0; k < Q; k++)
						y[iter][i][j][k] = (trajectory[iter + 1][i][j][k] - trajectory[iter][i][j][k])*temp / (AlfaParameter*AlfaParameter) + temp1[i][j][k] + temp2[i][j][k] + temp3[i][j][k] + temp4[i][j][k];


		}
		//Last block-row
		else
		{
			//Calculate diagonal block
			(*this).GetMacroscopic(iter); //macroscopic values for (i+0.5) time step
			dfdu_b(p, u0, m, n, cx, cy, w, rho, u, v, omega, eq, dummy, trajectory[iter], temp1, trajectory[iter + 1], x[iter]);
			dfdu_d(p, u0, m, n, cx, cy, w, rho, u, v, omega, eq, dummy, trajectory[iter], temp1, trajectory[iter + 1], temp1);
			(*this).GetMacroscopic(iter); //macroscopic values for (i-0.5) time step
			dfdu_b(p, u0, m, n, cx, cy, w, rho, u, v, omega, eq, dummy, trajectory[iter], temp2, trajectory[iter + 1], x[iter]);
			dfdu_d(p, u0, m, n, cx, cy, w, rho, u, v, omega, eq, dummy, trajectory[iter], temp2, trajectory[iter + 1], temp2);
			for (int i = 0; i < n; i++)
				for (int j = 0; j < m; j++)
					for (int k = 0; k < Q; k++)
						temp += (trajectory[iter][i][j][k] - trajectory[iter - 1][i][j][k])*x[iter][i][j][k];

			//Below diagonal
			(*this).GetMacroscopic(iter - 1); //macroscopic values for (i-0.5) time step
			dfdu_b(p, u0, m, n, cx, cy, w, rho, u, v, omega, eq, dummy, trajectory[iter - 1], temp3, trajectory[iter], x[iter - 1]);
			(*this).GetMacroscopic(iter); //macroscopic values for (i+0.5) time step
			dfdu_d(p, u0, m, n, cx, cy, w, rho, u, v, omega, eq, dummy, trajectory[iter], temp3, trajectory[iter + 1], temp3);

			//Calculate output
			for (int i = 0; i < n; i++)
				for (int j = 0; j < m; j++)
					for (int k = 0; k < Q; k++)
						y[iter][i][j][k] = (trajectory[iter + 1][i][j][k] - trajectory[iter][i][j][k])*temp / (AlfaParameter*AlfaParameter) + temp1[i][j][k] + temp2[i][j][k] + temp3[i][j][k];

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
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++)
		{
			delete[] temp1[i][j]; delete[] temp2[i][j];
			delete[] temp3[i][j]; delete[] temp4[i][j];
		}
		delete[] temp1[i];	delete[] temp2[i];
		delete[] temp3[i];	delete[] temp4[i];
	}
	delete[] temp1, temp2, temp3, temp4;
	delete dummy;
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
					if (trajectory[tsteps][i][j][k] != trajectory[tsteps][i][j][k]) std::cerr << "NaN trajectory element";
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
	bool ExportFlag = 1;
	std::ofstream CGM("CGM.txt");
	int N = mstep * m * n * Q; //System size
	int d = 3, D = 1;// pow(10.0, double(d / 2));
	int maxIter = 3*pow(10.0, double(d));
	double eps = 1e-5, alfa = 0.5, beta, buff = 0;
	double ****res_b, ****res, ****p, ****temp, ****x;

	//Create solution vectors
	b = new double ***[mstep - 1];	x = new double ***[mstep - 1];
	res = new double ***[mstep - 1];	res_b = new double ***[mstep - 1];
	p = new double ***[mstep - 1];	temp = new double ***[mstep - 1];
	for (int i = 0; i < mstep - 1; i++)
	{
		b[i] = new double **[n];	x[i] = new double **[n];
		res[i] = new double **[n];	res_b[i] = new double **[n];
		p[i] = new double **[n];	temp[i] = new double **[n];
		for (int j = 0; j < n; j++)
		{
			b[i][j] = new double *[m];	x[i][j] = new double *[m];
			res[i][j] = new double *[m];	res_b[i][j] = new double *[m];
			p[i][j] = new double *[m];	temp[i][j] = new double *[m];
			for (int k = 0; k < m; k++)
			{
				b[i][j][k] = new double[Q];	x[i][j][k] = new double[Q];
				res[i][j][k] = new double[Q];	res_b[i][j][k] = new double[Q];
				p[i][j][k] = new double[Q];	temp[i][j][k] = new double[Q];
				for (int l = 0; l < Q; l++)
				{
					b[i][j][k][l] = 0;	x[i][j][k][l] = 0;
					res[i][j][k][l] = 0;	res_b[i][j][k][l] = 0;
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
		auto start = std::chrono::steady_clock::now();
		(*this).AssemblyArray(p, temp); //calculate temp = Ap
		alfa = ScalarProduct(res, res_b)/ ScalarProduct(p, temp);
		
		//std::cout << "alfa = " << alfa << std::endl;
		//Export(res);
		buff = ScalarProduct(res_b, res);
		for (int i = 0; i < mstep - 1; i++)
			for (int j = 0; j < n; j++)
				for (int k = 0; k < m; k++)
					for (int l = 0; l < Q; l++)
					{
						x[i][j][k][l] = x[i][j][k][l] + alfa * p[i][j][k][l]; //updated value of x
						res[i][j][k][l] = res[i][j][k][l] - alfa * temp[i][j][k][l];
					}
		if (VectorNorm(res) < eps)	break;
		Preconditioner(res, res_b);
		beta = ScalarProduct(res_b, res) / buff;
		for (int i = 0; i < mstep - 1; i++)
			for (int j = 0; j < n; j++)
				for (int k = 0; k < m; k++)
					for (int l = 0; l < Q; l++)
					{
						p[i][j][k][l] = res_b[i][j][k][l] + beta * p[i][j][k][l];
					}
		auto end = std::chrono::steady_clock::now();
		if (iter%D == 0)	std::cout << "Residuum in " << iter + 1 << " iter = " << VectorNorm(res)
			<< "\t alfa = " << alfa << "\t beta = " << beta << "\tElapsed time(sec) : " 
			<< std::chrono::duration_cast<std::chrono::seconds>(end - start).count() << std::endl;
		if (ExportFlag)
		{
			CGM << iter + 1 << "\t" << VectorNorm(res) << "\t" << alfa << "\t" << beta << "\n";
		}
	}
	std::cout << "Total number of iterations: " << iter << std::endl;
	if (ExportFlag) CGM.close();
	Export(x);
	delete[] temp, p, res, res_b;
}

//Export results
void LSS::Export(double ****w)
{
	std::ofstream file("KKTSolution.txt");
	//First 5 lines are the parameters: u0,mstep,n,m,Q
	file << u0 << "\n" << mstep << "\n" << n << "\n" << m << "\n" << Q << "\n";
	for (int i = 0; i < mstep - 1; i++)
	{
		for (int j = 0; j < n; j++)
		{
			for (int k = 0; k < m; k++)
			{
				for (int l = 0; l < Q; l++)
				{
					file << w[i][j][k][l] ;
					if (~((j == n - 1) && (k == m - 1) && (l == Q - 1))) file << "\n";
				}
			}
		}
		file << "\n";
	}

	file.close();
}

double LSS::CalculateSensitivity(double ****KKT)
{
	double ****sol = KKT;
	double *eta = new double[mstep - 1];
	double ***temp, sum;
	double ****ve = new double ***[mstep-1];
	double J = 0, dJdu = 0, dJds = 0, Jtemp = 0, grad1 = 0, grad2 = 0, Jbar = 0, ***dummy = nullptr, ***temp1 = new double **[n], ***temp2 = new double **[n], ***temp3 = new double **[n], Jbuff = 0, dJdsbuff = 0, dJdubuff = 0;
	double **rhod = new double *[n], ***foutd = new double **[n];
	for (int i = 0; i < n; i++)
	{
		rhod[i] = new double[m];
		foutd[i] = new double *[m];
			for (int j = 0; j < m; j++)
			{
				rhod[i][j] = 0;
				foutd[i][j] = new double[Q];
				for (int k = 0; k < Q; k++)
					foutd[i][j][k] = 0;
			}
	}
	for (int i = 0; i < mstep - 1; i++)
	{
		ve[i] = new double **[n];
		for (int j = 0; j < n; j++)
		{
			ve[i][j] = new double *[m];
			for (int k = 0; k < m; k++)
			{
				ve[i][j][k] = new double[Q];
				for (int l = 0; l < Q; l++)
				{
					ve[i][j][k][l] = 0;
				}
			}
		}
	}
	temp = new double **[n];
	for (int j = 0; j < n; j++)
	{
		temp[j] = new double *[m];
		for (int k = 0; k < m; k++)
		{
			temp[j][k] = new double[Q];
			for (int l = 0; l < Q; l++)
			{
				temp[j][k][l] = 0;
			}
		}
	}
	for (int i = 0; i < n; i++)
	{
		temp1[i] = new double *[m];	
		temp2[i] = new double *[m];
		temp3[i] = new double *[m];
		for (int j = 0; j < m; j++)
		{
			temp1[i][j] = new double[Q]; 
			temp2[i][j] = new double[Q];
			temp3[i][j] = new double[Q];
			for (int k = 0; k < Q; k++)
			{
				temp1[i][j][k] = 0;
				temp2[i][j][k] = 0;
				temp3[i][j][k] = 0;
			}
		}
	}
	std::cout << "\nCALCULATE FIRST CONTRIBUTION\n";
	for (int iter = 0; iter < mstep - 1; iter++)
	{
		sum = 0;
		for (int j = 0; j < n; j++)
			for (int k = 0; k < m; k++)
				for (int l = 0; l < Q; l++)
				{
					temp[j][k][l] = (trajectory[iter + 1][j][k][l] - trajectory[iter][j][k][l]);
					sum += temp[j][k][l] * sol[iter][j][k][l];
				}
		eta[iter] = -sum / 1; //alfa parameter squared instead of 1 here
		//std::cout << "Eta in " << iter << " time step: " << eta[iter] << std::endl;
		(*this).GetMacroscopic(iter);
		dObjds_d(p,u0,1,m,n,cx,cy,w,rho,rhod,u,v,omega,eq,dummy, trajectory[iter],trajectory[iter + 1],foutd,&J,&dJds);
		if (iter < mstep - 2)
		{
			(*this).GetMacroscopic(iter); //macroscopic values for (i+0.5) time step
			dfdu_b(p, u0, m, n, cx, cy, w, rho, u, v, omega, eq, dummy, trajectory[iter], temp1, trajectory[iter + 1], sol[iter]);
			(*this).GetMacroscopic(iter); //macroscopic values for (i+0.5) time step
			dfdu_b(p, u0, m, n, cx, cy, w, rho, u, v, omega, eq, dummy, trajectory[iter], temp2, trajectory[iter + 1], sol[iter+1]);
			for (int j = 0; j < n; j++)
			{
				for (int k = 0; k < m; k++)
					for (int l = 0; l < Q; l++)
					{
						ve[iter][j][k][l] = (temp1[j][k][l] + temp2[j][k][l]);
					}

			}

		}
		else
		{
			(*this).GetMacroscopic(iter); //macroscopic values for (i+0.5) time step
			dfdu_b(p, u0, m, n, cx, cy, w, rho, u, v, omega, eq, dummy, trajectory[iter], temp1, trajectory[iter + 1], sol[iter]);
			for (int j = 0; j < n; j++)
				for (int k = 0; k < m; k++)
					for (int l = 0; l < Q; l++)
						ve[iter][j][k][l] = (temp1[j][k][l]);
		}
		(*this).GetMacroscopic(iter);
		dObjdu_d(p,u0,m,n,cx,cy,w,rho,rhod,u,v,omega,eq,dummy, trajectory[iter],ve[iter],trajectory[iter + 1],temp3,&J,&dJdu);
		
		grad1 += dJds - dJdu;
		dJdsbuff += dJds; dJdubuff += dJdu;
		(*this).GetMacroscopic(iter);
		Obj(p, u0, m, n, cx, cy, w, rho, u, v, omega, eq, trajectory[iter], trajectory[iter+1], &J);
		Jbuff += J;
		std::cout << "#";
	}
	std::cout << "\nFIRST CONTRIBUTION DONE\n";
	grad1 = grad1 / (mstep - 1);
	std::cout << "FIRST CONTRIBUTION VALUE: " << grad1;
	std::cout << "\ndJds CONTRIBUTION VALUE: " << dJdsbuff << "\ndJdu CONTRIBUTION VALUE: " << dJdubuff << std::endl;
	Jbuff = Jbuff / (mstep - 1);
	std::cout << "\nCALCULATE SECOND CONTRIBUTION\n";
	for (int iter = 1; iter < mstep - 1; iter++)
	{
		(*this).GetMacroscopic(iter-1);
		Obj(p, u0, m, n, cx, cy, w, rho, u, v, omega, eq, trajectory[iter - 1], trajectory[iter], &Jtemp);
		(*this).GetMacroscopic(iter);
		Obj(p, u0, m, n, cx, cy, w, rho, u, v, omega, eq, trajectory[iter], trajectory[iter + 1], &J);
		J = (J + Jtemp) / 2;
		J = J - Jbuff;
		grad2 += eta[iter] * J;
		std::cout << "#";
	}
	std::cout << "\SECOND CONTRIBUTION DONE\n";
	grad2 = grad2 / (mstep - 2);
	std::cout << "SECOND CONTRIBUTION VALUE: " << grad2;
	Jbar = grad1 + grad2;
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++)
		{
			delete[] temp1[i][j];
			delete[] temp2[i][j];
			delete[] temp3[i][j];
		}
		delete[] temp1[i];
		delete[] temp2[i];
		delete[] temp3[i];
	}
	delete[] temp1, temp2, temp3;
	delete dummy;
	
	std::cout << "\n\nLSS SENTIVITY: " << Jbar;
	return Jbar;
}