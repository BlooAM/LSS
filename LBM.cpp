#include "LBM.h"


LBM::LBM()
{
	//Case parameters
	tstep = 0;
	L = 50; H = 2;	p = 0.5;
	m = 20; n = 10 * m; //m = 40, n = 25m
	mstep = 3000; //400
	m0 = p * m; n0 = m;
	u0 = 0.1, rho0 = 5, alfa = 0.01;
	omega = 0;

	//Grid constants
	w[0] = 4./9.;
	for (int i = 1; i < 5; i++)
	{
		w[i] = 1./9.;
	}
	for (int i = 5; i < 9; i++)
	{
		w[i] = 1. / 36.;
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
			rho[j][k] = rho0;
			u[j][k] = 0;
			v[j][k] = 0;
		}
	}
	for (int i = 1; i < m - 1; i++)
	{
		u[0][i] = u0;
		v[0][i] = 0;
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

void LBM::ApplyBC()
{
	double rhow;

	//Inlet
	for (int j = m0; j < m; j++)
	{
		rhow = (u_ref[0][tstep][0][j]+ u_ref[2][tstep][0][j]+ u_ref[4][tstep][0][j]+2*(u_ref[3][tstep][0][j]+ u_ref[6][tstep][0][j]+ u_ref[7][tstep][0][j])) / (1 - u0);
		u_ref[1][tstep][0][j] = u_ref[3][tstep][0][j] + 2*rhow*u0 / 3;
		u_ref[5][tstep][0][j] = u_ref[7][tstep][0][j] +  rhow*u0 / 6;
		u_ref[8][tstep][0][j] = u_ref[6][tstep][0][j] + rhow*u0 / 6;
	}
	//South
	for (int i = 0; i < n; i++)
	{
		u_ref[2][tstep][i][0] = u_ref[4][tstep][i][0];
		u_ref[5][tstep][i][0] = u_ref[7][tstep][i][0];
		u_ref[6][tstep][i][0] = u_ref[8][tstep][i][0];
	}
	//North
	for (int i = 0; i < n; i++)
	{
		u_ref[4][tstep][i][m - 1] = u_ref[2][tstep][i][m - 1];
		u_ref[8][tstep][i][m - 1] = u_ref[6][tstep][i][m - 1];
		u_ref[7][tstep][i][m - 1] = u_ref[5][tstep][i][m - 1];
	}
	//Outlet - extrapolation
	for (int j = 0; j < m; j++)
	{
		u_ref[1][tstep][n - 1][j] = 2 * u_ref[1][tstep][n - 2][j] - u_ref[1][tstep][n - 3][j];
		u_ref[5][tstep][n - 1][j] = 2 * u_ref[5][tstep][n - 2][j] - u_ref[5][tstep][n - 3][j];
		u_ref[8][tstep][n - 1][j] = 2 * u_ref[8][tstep][n - 2][j] - u_ref[8][tstep][n - 3][j];
	}
	//Back-facing step
	for (int i = 0; i < n0; i++)
	{
		u_ref[2][tstep][i][m0 - 1] = u_ref[4][tstep][i][m0 - 1];
		u_ref[5][tstep][i][m0 - 1] = u_ref[7][tstep][i][m0 - 1];
		u_ref[6][tstep][i][m0 - 1] = u_ref[8][tstep][i][m0 - 1];
	}
	for (int j = 0; j < m0; j++)
	{
		u_ref[1][tstep][n0 - 1][j] = u_ref[3][tstep][n0 - 1][j];
		u_ref[5][tstep][n0 - 1][j] = u_ref[7][tstep][n0 - 1][j];
		u_ref[8][tstep][n0 - 1][j] = u_ref[6][tstep][n0 - 1][j]; //???
	}
}

void LBM::CalculateMacroscopic()
{
	double ssum,usum,vsum;
	for (int j = 0; j < m; j++)
	{
		for (int i = 0; i < n; i++)
		{
			ssum = 0;
			for (int k = 0; k < 9; k++)
			{
				ssum = ssum + u_ref[k][tstep][i][j];
			}
			rho[i][j] = ssum;
		}
	}

	for (int i = 0; i < n; i++)
	{
		rho[i][m - 1] = u_ref[0][tstep][i][m - 1] + u_ref[1][tstep][i][m - 1] + u_ref[3][tstep][i][m - 1] + 2 * (u_ref[2][tstep][i][m - 1] + u_ref[6][tstep][i][m - 1] + u_ref[5][tstep][i][m - 1]);
	}

	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m - 1; j++)
		{
			usum = 0;	vsum = 0;
			for (int k = 0; k < 9; k++)
			{
				usum = usum + u_ref[k][tstep][i][j] * cx[k];
				vsum = vsum + u_ref[k][tstep][i][j] * cy[k];
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
}
void LBM::CollisionStep()
{
	double temp1,temp2;
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++)
		{
			temp1 = u[i][j] * u[i][j] + v[i][j] * v[i][j];
			for (int k = 0; k < 9; k++)
			{
				temp2 = u[i][j] * cx[k] + v[i][j] * cy[k];
				feq[k][i][j] = rho[i][j] * w[k] * (1 + 3*temp2 + 4.5*temp2*temp2 - 1.5*temp1);
				u_ref[k][tstep][i][j] = omega * feq[k][i][j] + (1 - omega)*u_ref[k][tstep][i][j];
			}
		}
	}
}
void LBM::Exectue()
{
	//Parameters for back-facing step case
	double sumvel = 0;
	omega = 1 / (3 * alfa + 0.5);
	std::cout << "Reynolds number: " << u0 * m / alfa << std::endl;
	
	for (int i = 0; i < mstep-1; i++) //mstep-1
	{
		std::cout << tstep << std::endl;
		CollisionStep();
		//tstep++;
		StreamingStep();
		ApplyBC();
		CalculateMacroscopic();
	}
	PostProcess();
}

void LBM::StreamingStep()
{
	for (int j = 0; j < m; j++)
	{
		for (int i = n - 1; i > 0; i--)
		{
			u_ref[1][tstep][i][j] = u_ref[1][tstep][i - 1][j];
		}

		for (int i = 0; i < n - 1; i++)
		{
			u_ref[3][tstep][i][j] = u_ref[3][tstep][i + 1][j];
		}
	}

	for (int j = m - 1; j > 0; j--)
	{
		for (int i = 0; i < n; i++)
		{
			u_ref[2][tstep][i][j] = u_ref[2][tstep][i][j - 1];
		}
		for (int i = n - 1; i > 0; i--)
		{
			u_ref[5][tstep][i][j] = u_ref[5][tstep][i - 1][j - 1];
		}
		for (int i = 0; i < n - 1; i++)
		{
			u_ref[6][tstep][i][j] = u_ref[6][tstep][i + 1][j - 1];
		}
	}
	for (int j = 0; j < m - 1; j++)
	{
		for (int i = 0; i < n; i++)
		{
			u_ref[4][tstep][i][j] = u_ref[4][tstep][i][j + 1];
		}
		for (int i = 0; i < n - 1; i++)
		{
			u_ref[7][tstep][i][j] = u_ref[7][tstep][i + 1][j + 1];
		}
		for (int i = n - 1; i > 0; i--)
		{
			u_ref[8][tstep][i][j] = u_ref[8][tstep][i - 1][j + 1];
		}
	}
}
void LBM::PostProcess()
{
	std::ofstream file("result.txt");

	file << "Rho\n";
	for (int j = 0; j < m; j++)
	{
		for (int i = 0; i < n; i++)
		{
			file << "\t" << rho[i][j];
		}
		file << "\n";
	}

	file << "U\n";
	for (int j = 0; j < m; j++)
	{
		for (int i = 0; i < n; i++)
		{
			file << "\t" << u[i][j];
		}
		file << "\n";
	}

	file << "V\n";
	for (int j = 0; j < m; j++)
	{
		for (int i = 0; i < n; i++)
		{
			file << "\t" << v[i][j];
		}
		file << "\n";
	}
	file.close();
}
