#include "LBM.h"

//Friend functions first
void SolveTimeStep(float p, float u0, int m, int n, float *cx, float *cy, float *w, float **rho, float **u, float **v, float omega, float ***feq, float*** fin, float*** fout) //Friend function
{
	int m0 = p * m, n0 = m;
	//Collision step
	float temp1, temp2, rhow, ssum, usum, vsum;
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

	/*for (int i = 0; i < n; i++)
	{
		rho[i][m - 1] = fout[i][m - 1][0] + fout[i][m - 1][1] + fout[i][m - 1][3] + 2 * (fout[i][m - 1][2] + fout[i][m - 1][6] + fout[i][m - 1][5]);
	}*/

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
}

LBM::LBM()
{
	//Case parameters
	tstep = 0;
	L = 50; H = 2;	p = 0.5;
	m = 20; n = 10 * m; //m = 40, n = 25m
	mstep = 50; //4000
	m0 = p * m; n0 = p * m;
	u0 = 0.1, rho0 = 5, alfa = 0.01;
	omega = 0;
	J = 0; Jbar = 0;

	//Grid constants
	w[0] = 4. / 9.;
	for (int i = 1; i < 5; i++)
	{
		w[i] = 1. / 9.;
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
	/*for (int iter = 0; iter < 9; ++iter)
	{
		u_ref[iter] = new float **[mstep];
		for (int i = 0; i < mstep; ++i)
		{
			u_ref[iter][i] = new float*[n];
			for (int j = 0; j < n; ++j)
			{
				u_ref[iter][i][j] = new float[m];
				for (int k = 0; k < m; ++k)
				{
					u_ref[iter][i][j][k] = 0;
				}
			}
		}
	}*/
	f = new float ***[mstep];
	for (int iter = 0; iter < mstep; ++iter)
	{
		f[iter] = new float **[n];
		for (int i = 0; i < n; ++i)
		{
			f[iter][i] = new float*[m];
			for (int j = 0; j < m; ++j)
			{
				f[iter][i][j] = new float[9];
				for (int k = 0; k < 9; ++k)
				{
					f[iter][i][j][k] = 0;
				}
			}
		}
	}

	for (int iter = 0; iter < 9; ++iter)
	{
		feq[iter] = new float*[n];

		for (int j = 0; j < n; ++j)
		{
			feq[iter][j] = new float[m];

			for (int k = 0; k < m; ++k)
			{
				feq[iter][j][k] = 0;
			}
		}
	}

	rho = new float *[n];
	u = new float *[n];
	v = new float *[n];
	for (int j = 0; j < n; ++j)
	{
		rho[j] = new float[m];
		u[j] = new float[m];
		v[j] = new float[m];
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

LBM::LBM(float u0_, int mstep_, int m_, int mx)
{
	//Case parameters
	tstep = 0;
	L = 50; H = 2;	p = 0.5;
	u0 = u0_; mstep = mstep_;
	m = m_; n = mx * m;
	m0 = p * m; n0 = m;
	rho0 = 5, alfa = 0.01;
	omega = 0;
	J = 0; Jbar = 0;

	//Grid constants
	w[0] = 4. / 9.;
	for (int i = 1; i < 5; i++)
	{
		w[i] = 1. / 9.;
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
	f = new float ***[mstep];
	for (int iter = 0; iter < mstep; ++iter)
	{
		f[iter] = new float **[n];
		for (int i = 0; i < n; ++i)
		{
			f[iter][i] = new float*[m];
			for (int j = 0; j < m; ++j)
			{
				f[iter][i][j] = new float[9];
				for (int k = 0; k < 9; ++k)
				{
					f[iter][i][j][k] = 0;
				}
			}
		}
	}

	for (int iter = 0; iter < 9; ++iter)
	{
		feq[iter] = new float*[n];

		for (int j = 0; j < n; ++j)
		{
			feq[iter][j] = new float[m];

			for (int k = 0; k < m; ++k)
			{
				feq[iter][j][k] = 0;
			}
		}
	}

	rho = new float *[n];
	u = new float *[n];
	v = new float *[n];
	for (int j = 0; j < n; ++j)
	{
		rho[j] = new float[m];
		u[j] = new float[m];
		v[j] = new float[m];
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

	/*for (int iter = 0; iter < 9; ++iter)
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
	}*/

	for (int iter = 0; iter < mstep; ++iter)
	{
		for (int i = 0; i < n; ++i)
		{
			for (int j = 0; j < m; ++j)
			{
				delete[] f[iter][i][j];
			}
			delete[] f[iter][i];
		}
		delete[] f[iter];
	}
}

void LBM::ApplyBC()
{
	float rhow;

	//Inlet
	for (int j = m0; j < m; j++)
	{
		rhow = (u_ref[0][tstep][0][j] + u_ref[2][tstep][0][j] + u_ref[4][tstep][0][j] + 2 * (u_ref[3][tstep][0][j] + u_ref[6][tstep][0][j] + u_ref[7][tstep][0][j])) / (1 - u0);
		u_ref[1][tstep][0][j] = u_ref[3][tstep][0][j] + 2 * rhow*u0 / 3;
		u_ref[5][tstep][0][j] = u_ref[7][tstep][0][j] + rhow * u0 / 6;
		u_ref[8][tstep][0][j] = u_ref[6][tstep][0][j] + rhow * u0 / 6;
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
	float ssum, usum, vsum;
	for (int j = 0; j < m; j++)
	{
		for (int i = 0; i < n; i++)
		{
			ssum = 0;
			for (int k = 0; k < 9; k++)
			{
				ssum = ssum + u_ref[k][tstep + 1][i][j];
			}
			rho[i][j] = ssum;
		}
	}

	for (int i = 0; i < n; i++)
	{
		rho[i][m - 1] = u_ref[0][tstep + 1][i][m - 1] + u_ref[1][tstep + 1][i][m - 1] + u_ref[3][tstep + 1][i][m - 1] + 2 * (u_ref[2][tstep + 1][i][m - 1] + u_ref[6][tstep + 1][i][m - 1] + u_ref[5][tstep + 1][i][m - 1]);
	}

	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m - 1; j++)
		{
			usum = 0;	vsum = 0;
			for (int k = 0; k < 9; k++)
			{
				usum = usum + u_ref[k][tstep + 1][i][j] * cx[k];
				vsum = vsum + u_ref[k][tstep + 1][i][j] * cy[k];
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
	float temp1, temp2;
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++)
		{
			temp1 = u[i][j] * u[i][j] + v[i][j] * v[i][j];
			for (int k = 0; k < 9; k++)
			{
				temp2 = u[i][j] * cx[k] + v[i][j] * cy[k];
				feq[k][i][j] = rho[i][j] * w[k] * (1 + 3 * temp2 + 4.5*temp2*temp2 - 1.5*temp1);
				u_ref[k][tstep + 1][i][j] = omega * feq[k][i][j] + (1 - omega)*u_ref[k][tstep][i][j];
			}
		}
	}

	//Apply BC
	float rhow;

	//Inlet
	for (int j = m0; j < m; j++)
	{
		rhow = (u_ref[0][tstep + 1][0][j] + u_ref[2][tstep + 1][0][j] + u_ref[4][tstep + 1][0][j] + 2 * (u_ref[3][tstep + 1][0][j] + u_ref[6][tstep + 1][0][j] + u_ref[7][tstep + 1][0][j])) / (1 - u0);
		u_ref[1][tstep + 1][0][j] = u_ref[3][tstep + 1][0][j] + 2 * rhow*u0 / 3;
		u_ref[5][tstep + 1][0][j] = u_ref[7][tstep + 1][0][j] + rhow * u0 / 6;
		u_ref[8][tstep + 1][0][j] = u_ref[6][tstep + 1][0][j] + rhow * u0 / 6;
	}
	//South
	for (int i = 0; i < n; i++)
	{
		u_ref[2][tstep + 1][i][0] = u_ref[4][tstep + 1][i][0];
		u_ref[5][tstep + 1][i][0] = u_ref[7][tstep + 1][i][0];
		u_ref[6][tstep + 1][i][0] = u_ref[8][tstep + 1][i][0];
	}
	//North
	for (int i = 0; i < n; i++)
	{
		u_ref[4][tstep + 1][i][m - 1] = u_ref[2][tstep + 1][i][m - 1];
		u_ref[8][tstep + 1][i][m - 1] = u_ref[6][tstep + 1][i][m - 1];
		u_ref[7][tstep + 1][i][m - 1] = u_ref[5][tstep + 1][i][m - 1];
	}
	//Outlet - extrapolation
	for (int j = 0; j < m; j++)
	{
		u_ref[1][tstep + 1][n - 1][j] = 2 * u_ref[1][tstep + 1][n - 2][j] - u_ref[1][tstep + 1][n - 3][j];
		u_ref[5][tstep + 1][n - 1][j] = 2 * u_ref[5][tstep + 1][n - 2][j] - u_ref[5][tstep + 1][n - 3][j];
		u_ref[8][tstep + 1][n - 1][j] = 2 * u_ref[8][tstep + 1][n - 2][j] - u_ref[8][tstep + 1][n - 3][j];
	}
	//Back-facing step
	for (int i = 0; i < n0; i++)
	{
		u_ref[2][tstep + 1][i][m0 - 1] = u_ref[4][tstep + 1][i][m0 - 1];
		u_ref[5][tstep + 1][i][m0 - 1] = u_ref[7][tstep + 1][i][m0 - 1];
		u_ref[6][tstep + 1][i][m0 - 1] = u_ref[8][tstep + 1][i][m0 - 1];
	}
	for (int j = 0; j < m0; j++)
	{
		u_ref[1][tstep + 1][n0 - 1][j] = u_ref[3][tstep + 1][n0 - 1][j];
		u_ref[5][tstep + 1][n0 - 1][j] = u_ref[7][tstep + 1][n0 - 1][j];
		u_ref[8][tstep + 1][n0 - 1][j] = u_ref[6][tstep + 1][n0 - 1][j]; //???
	}
}
void LBM::Exectue()
{
	//Parameters for back-facing step case
	float sumvel = 0, pin, pout, temp;
	omega = 1 / (3 * alfa + 0.5);
	std::cout << "Reynolds number: " << u0 * m / alfa << std::endl;
	J = 0; Jbar = 0; temp = 0;
	std::cout << "\nSOLVING TRAJECTORY\n";
	for (int i = 0; i < mstep - 1; i++) //mstep-1
	{
		pin = 0; pout = 0;
		//std::cout << i << std::endl;
		//CollisionStep();
		//StreamingStep();
		//CalculateMacroscopic();
		SolveTimeStep(p, u0, m, n, cx, cy, w, rho, u, v, omega, feq, f[i], f[i + 1]); //Friend function for TAPENADE
		for (int j = m0; j < m; j++)	pin += rho[0][j] / 3;
		for (int j = 0; j < m; j++)	pout += rho[n - 1][j] / 3;
		pin /= (m - m0); pout /= m;
		if (i > 0) J += (pin - pout + temp) / 2;
		temp = pin - pout;
	}
	std::cout << "DONE\n";
	Jbar = J / (mstep - 1);
	std::cout << "Long-time average pressure drop: " << Jbar << "\n\n";
	PostProcess();
}


void LBM::GetEqulibrium(int i, float** rho, float** u, float** v, float ***eq)
{
	float temp1, temp2;
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++)
		{
			temp1 = u[i][j] * u[i][j] + v[i][j] * v[i][j];
			for (int k = 0; k < 9; k++)
			{
				temp2 = u[i][j] * cx[k] + v[i][j] * cy[k];
				eq[k][i][j] = rho[i][j] * w[k] * (1 + 3 * temp2 + 4.5*temp2*temp2 - 1.5*temp1);
			}
		}
	}
}

void LBM::StreamingStep()
{
	for (int j = 0; j < m; j++)
	{
		for (int i = n - 1; i > 0; i--)
		{
			u_ref[1][tstep + 1][i][j] = u_ref[1][tstep + 1][i - 1][j];
		}

		for (int i = 0; i < n - 1; i++)
		{
			u_ref[3][tstep + 1][i][j] = u_ref[3][tstep + 1][i + 1][j];
		}
	}

	for (int j = m - 1; j > 0; j--)
	{
		for (int i = 0; i < n; i++)
		{
			u_ref[2][tstep + 1][i][j] = u_ref[2][tstep + 1][i][j - 1];
		}
		for (int i = n - 1; i > 0; i--)
		{
			u_ref[5][tstep + 1][i][j] = u_ref[5][tstep + 1][i - 1][j - 1];
		}
		for (int i = 0; i < n - 1; i++)
		{
			u_ref[6][tstep + 1][i][j] = u_ref[6][tstep + 1][i + 1][j - 1];
		}
	}
	for (int j = 0; j < m - 1; j++)
	{
		for (int i = 0; i < n; i++)
		{
			u_ref[4][tstep + 1][i][j] = u_ref[4][tstep + 1][i][j + 1];
		}
		for (int i = 0; i < n - 1; i++)
		{
			u_ref[7][tstep + 1][i][j] = u_ref[7][tstep + 1][i + 1][j + 1];
		}
		for (int i = n - 1; i > 0; i--)
		{
			u_ref[8][tstep + 1][i][j] = u_ref[8][tstep + 1][i - 1][j + 1];
		}
	}
}
void LBM::PostProcess()
{
	std::ofstream file("result.txt");

	file << "Rho\n";
	file << "U\n";
	file << "V\n";

	for (int j = 0; j < m; j++)
	{
		for (int i = 0; i < n; i++)
		{
			file << "\t" << rho[i][j];
			if (rho[i][j] != rho[i][j]) std::cerr << "\nNaN rho value\n";
		}
		file << "\n";
	}
	file << "\n";
	for (int j = 0; j < m; j++)
	{
		for (int i = 0; i < n; i++)
		{
			file << "\t" << u[i][j];
			if (u[i][j] != u[i][j]) std::cerr << "\nNaN u value\n";
		}
		file << "\n";
	}
	file << "\n";
	for (int j = 0; j < m; j++)
	{
		for (int i = 0; i < n; i++)
		{
			file << "\t" << v[i][j];
			if (v[i][j] != v[i][j]) std::cerr << "\nNaN v value\n";
		}
		file << "\n";
	}
	file.close();
}
