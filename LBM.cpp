#include "LBM.h"


LBM::LBM()
{
	//Case parameters
	tstep = 0;
	L = 50; H = 2;	p = 0.5;
	m = 20; n = 10 * m; //m = 40, n = 25m
	mstep = 3000; //400
	m0 = p * m; n0 = m;
	u0 = 0.1, alfa = 0.01;

	//Grid constants
	w[0] = 4 / 9;
	for (int i = 1; i < 9; i++)
	{
		w[i] = 1 / 9;
		w[i + 4] = 1 / 36;
	}
	cx[0] = 0; cy[0] = 0;
	cx[1] = 1; cy[1] = 0;
	cx[2] = 0; cy[2] = 1;
	cx[3] = -1; cy[3] = 0;
	cx[4] = 0; cy[4] = -1;
	cx[5] = 1; cy[5] = 1;
	cx[6] = -1; cy[6] = 1;
	cx[7] = -1; cy[7] = -1;
	cx[8] = 1; cy[8] = -1;

	//Pointers
	tSpan = nullptr;
	for (int iter = 0; iter < 9; ++iter)
	{
		u_ref[iter] = new double **[mstep];
		for (int i = 0; i < mstep; ++i)
		{
			u_ref[iter][i] = new double*[n];
			for (int j = 0; j < n; ++j)
			{
				u_ref[iter][i][j] = new double[m];
				for (int k = 0; k < m; ++k)
				{
					u_ref[iter][i][j][k] = 0;
				}
			}
		}
	}

	for (int iter = 0; iter < 9; ++iter)
	{
		feq[iter] = new double*[n];

		for (int j = 0; j < n; ++j)
		{
			feq[iter][j] = new double[m];

			for (int k = 0; k < m; ++k)
			{
				feq[iter][j][k] = 0;
			}
		}
	}

	rho = new double *[n];
	u = new double *[n];
	v = new double *[n];
	for (int j = 0; j < n; ++j)
	{
		rho[j] = new double[m];
		u[j] = new double[m];
		v[j] = new double[m];

		for (int k = 0; k < m; ++k)
		{
			rho[j][k] = 0;
			if (j == 0)
				u[j][k] = u0;
			else
				u[j][k] = 0;
			v[j][k] = 0;
		}
	}
}

LBM::~LBM()
{
	for (int j = 0; j < n; ++j)
	{
		delete[] rho[j], u[j], v[j];
	}
	delete[] rho, u, v;

	for (int iter = 0; iter < 9; ++iter)
	{
		for (int j = 0; j < n; ++j)
		{
			delete[] feq[iter][j];
		}
		delete[] feq[iter];
	}

	for (int iter = 0; iter < 9; ++iter)
	{
		for (int i = 0; i < mstep; ++i)
		{
			for (int j = 0; j < n; ++j)
			{
				delete[] u_ref[iter][i][j];
			}
			delete[] u_ref[iter][i];
		}
		delete[] u_ref[iter];
	}

}

void LBM::CalculateEqulibrium()
{

}
void LBM::CalculateMacroscopic()
{

}
void LBM::CollisionStep()
{

}
void LBM::Exectue()
{
	//Parameters for back-facing step case
	double omega, sumvel = 0, rho0 = 5;
	omega = 1 / (3 * alfa + 0.5);
	std::cout << "Reynolds number: " << u0 * m / alfa << std::endl;

	for (int i = 0; i < mstep; i++)
	{
		StreamingStep();

		tstep++;
	}
}

void LBM::StreamingStep()
{
	for (int j = 0; j < m; j++)
	{
		for (int i = n - 1; i > 0; i--)
		{
			u_ref[1][tstep + 1][i][j] = u_ref[1][tstep][i - 1][j];
		}

		for (int i = 0; i < n - 1; i++)
		{
			u_ref[3][tstep + 1][i][j] = u_ref[3][tstep][i + 1][j];
		}
	}

	for (int j = m - 1; j > 0; j--)
	{
		for (int i = 0; i < n; i++)
		{
			u_ref[2][tstep + 1][i][j] = u_ref[2][tstep][i][j - 1];
		}
		for (int i = n - 1; i > 0; i--)
		{
			u_ref[5][tstep + 1][i][j] = u_ref[5][tstep][i - 1][j - 1];
		}
		for (int i = 0; i < n - 1; i++)
		{
			u_ref[6][tstep + 1][i][j] = u_ref[6][tstep][i + 1][j - 1];
		}
	}
	for (int j = 0; j < m - 1; j++)
	{
		for (int i = 0; i < n; i++)
		{
			u_ref[4][tstep + 1][i][j] = u_ref[4][tstep][i][j + 1];
		}
		for (int i = 0; i < n - 1; i++)
		{
			u_ref[7][tstep + 1][i][j] = u_ref[7][tstep][i + 1][j + 1];
		}
		for (int i = n - 1; i > 0; i--)
		{
			u_ref[8][tstep + 1][i][j] = u_ref[8][tstep][i - 1][j + 1];
		}

	}
}

