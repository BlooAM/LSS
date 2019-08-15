#include "Derivatives.h"

void dfds_b(double p, double u0, double *u0b, int m, int n, double *cx,
	double *cy, double *w, double **rho, double **u, double **v, double
	omega, double ***feq, double ***feqb, double ***fin, double ***fout,
	double ***foutb) {
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
	//Friend function
	int m0, n0;
	m0 = p * m;
	n0 = m;
	//Collision step
	double temp1, temp2, rhow = 0, ssum, usum, vsum;
	double rhowb;
	for (int i = 0; i < n; ++i)
		for (int j = 0; j < m; ++j) {
			temp1 = u[i][j] * u[i][j] + v[i][j] * v[i][j];
			for (int k = 0; k < 9; ++k) {
				temp2 = u[i][j] * cx[k] + v[i][j] * cy[k];
				feq[k][i][j] = rho[i][j] * w[k] * (1 + 3 * temp2 + 4.5*temp2*temp2 - 1.5*temp1);
				fout[i][j][k] = omega * feq[k][i][j] + (1 - omega)*fin[i][j][k];
			}
		}
	//Apply BC
	//Inlet
	for (int j = m0; j < m; ++j) {
		pushreal8_(rhow);
		rhow = (fout[0][j][0] + fout[0][j][2] + fout[0][j][4] +
			2 * (fout[0][j][3] + fout[0][j][6] + fout[0][j][7])) / (1 - u0);
		pushreal8_(fout[0][j][1]);
		fout[0][j][1] = fout[0][j][3] + 2 * rhow*u0 / 3;
		pushreal8_(fout[0][j][5]);
		fout[0][j][5] = fout[0][j][7] + rhow * u0 / 6;
		pushreal8_(fout[0][j][8]);
		fout[0][j][8] = fout[0][j][6] + rhow * u0 / 6;
	}
	{
		double tmp;
		double tmp0;
		double tmp1;
	}
	//Streaming step
	for (int j = 0; j < m; ++j) {
		{
			double tmp;
		}
		{
			double tmp;
		}
	}
	for (int j = m - 1; j > 0; --j) {
		{
			double tmp;
		}
		{
			double tmp;
		}
		{
			double tmp;
		}
	}
	for (int j = 0; j < m - 1; ++j) {
		{
			double tmp;
		}
		{
			double tmp;
		}
		{
			double tmp;
		}
	}
	pushreal8_(rhow);
	pushinteger4_(n0);
	pushinteger4_(m0);
	popinteger4_(&m0);
	popinteger4_(&n0);
	popreal8_(&rhow);
	for (int j = m - 2; j > -1; --j) {
		{
			double tmpb;
			for (int i = 1; i < n; ++i) {
				tmpb = foutb[i][j][8];
				foutb[i][j][8] = 0.0;
				foutb[i - 1][j + 1][8] = foutb[i - 1][j + 1][8] + tmpb;
			}
		}
		{
			double tmpb;
			for (int i = n - 2; i > -1; --i) {
				tmpb = foutb[i][j][7];
				foutb[i][j][7] = 0.0;
				foutb[i + 1][j + 1][7] = foutb[i + 1][j + 1][7] + tmpb;
			}
		}
		{
			double tmpb;
			for (int i = n - 1; i > -1; --i) {
				tmpb = foutb[i][j][4];
				foutb[i][j][4] = 0.0;
				foutb[i][j + 1][4] = foutb[i][j + 1][4] + tmpb;
			}
		}
	}
	for (int j = 1; j < m; ++j) {
		{
			double tmpb;
			for (int i = n - 2; i > -1; --i) {
				tmpb = foutb[i][j][6];
				foutb[i][j][6] = 0.0;
				foutb[i + 1][j - 1][6] = foutb[i + 1][j - 1][6] + tmpb;
			}
		}
		{
			double tmpb;
			for (int i = 1; i < n; ++i) {
				tmpb = foutb[i][j][5];
				foutb[i][j][5] = 0.0;
				foutb[i - 1][j - 1][5] = foutb[i - 1][j - 1][5] + tmpb;
			}
		}
		{
			double tmpb;
			for (int i = n - 1; i > -1; --i) {
				tmpb = foutb[i][j][2];
				foutb[i][j][2] = 0.0;
				foutb[i][j - 1][2] = foutb[i][j - 1][2] + tmpb;
			}
		}
	}
	for (int j = m - 1; j > -1; --j) {
		{
			double tmpb;
			for (int i = n - 2; i > -1; --i) {
				tmpb = foutb[i][j][3];
				foutb[i][j][3] = 0.0;
				foutb[i + 1][j][3] = foutb[i + 1][j][3] + tmpb;
			}
		}
		{
			double tmpb;
			for (int i = 1; i < n; ++i) {
				tmpb = foutb[i][j][1];
				foutb[i][j][1] = 0.0;
				foutb[i - 1][j][1] = foutb[i - 1][j][1] + tmpb;
			}
		}
	}
	for (int j = m0 - 1; j > -1; --j) {
		foutb[n0 - 1][j][6] = foutb[n0 - 1][j][6] + foutb[n0 - 1][j][8];
		foutb[n0 - 1][j][8] = 0.0;
		foutb[n0 - 1][j][7] = foutb[n0 - 1][j][7] + foutb[n0 - 1][j][5];
		foutb[n0 - 1][j][5] = 0.0;
		foutb[n0 - 1][j][3] = foutb[n0 - 1][j][3] + foutb[n0 - 1][j][1];
		foutb[n0 - 1][j][1] = 0.0;
	}
	for (int i = n0 - 1; i > -1; --i) {
		foutb[i][m0 - 1][8] = foutb[i][m0 - 1][8] + foutb[i][m0 - 1][6];
		foutb[i][m0 - 1][6] = 0.0;
		foutb[i][m0 - 1][7] = foutb[i][m0 - 1][7] + foutb[i][m0 - 1][5];
		foutb[i][m0 - 1][5] = 0.0;
		foutb[i][m0 - 1][4] = foutb[i][m0 - 1][4] + foutb[i][m0 - 1][2];
		foutb[i][m0 - 1][2] = 0.0;
	}
	{
		double tmpb;
		double tmpb0;
		double tmpb1;
		for (int j = m - 1; j > -1; --j) {
			tmpb = foutb[n - 1][j][8];
			foutb[n - 1][j][8] = 0.0;
			foutb[n - 2][j][8] = foutb[n - 2][j][8] + 2 * tmpb;
			foutb[n - 3][j][8] = foutb[n - 3][j][8] - tmpb;
			tmpb0 = foutb[n - 1][j][5];
			foutb[n - 1][j][5] = 0.0;
			foutb[n - 2][j][5] = foutb[n - 2][j][5] + 2 * tmpb0;
			foutb[n - 3][j][5] = foutb[n - 3][j][5] - tmpb0;
			tmpb1 = foutb[n - 1][j][1];
			foutb[n - 1][j][1] = 0.0;
			foutb[n - 2][j][1] = foutb[n - 2][j][1] + 2 * tmpb1;
			foutb[n - 3][j][1] = foutb[n - 3][j][1] - tmpb1;
		}
	}
	for (int i = n - 1; i > -1; --i) {
		foutb[i][m - 1][5] = foutb[i][m - 1][5] + foutb[i][m - 1][7];
		foutb[i][m - 1][7] = 0.0;
		foutb[i][m - 1][6] = foutb[i][m - 1][6] + foutb[i][m - 1][8];
		foutb[i][m - 1][8] = 0.0;
		foutb[i][m - 1][2] = foutb[i][m - 1][2] + foutb[i][m - 1][4];
		foutb[i][m - 1][4] = 0.0;
	}
	for (int i = n - 1; i > -1; --i) {
		foutb[i][0][8] = foutb[i][0][8] + foutb[i][0][6];
		foutb[i][0][6] = 0.0;
		foutb[i][0][7] = foutb[i][0][7] + foutb[i][0][5];
		foutb[i][0][5] = 0.0;
		foutb[i][0][4] = foutb[i][0][4] + foutb[i][0][2];
		foutb[i][0][2] = 0.0;
	}
	{
		double tempb;
		*u0b = 0.0;
		for (int j = m - 1; j > m0 - 1; --j) {
			popreal8_(&(fout[0][j][8]));
			foutb[0][j][6] = foutb[0][j][6] + foutb[0][j][8];
			rhowb = u0 * foutb[0][j][8] / 6;
			*u0b = *u0b + rhow * foutb[0][j][8] / 6;
			foutb[0][j][8] = 0.0;
			popreal8_(&(fout[0][j][5]));
			foutb[0][j][7] = foutb[0][j][7] + foutb[0][j][5];
			rhowb = rhowb + u0 * foutb[0][j][5] / 6;
			*u0b = *u0b + rhow * foutb[0][j][5] / 6;
			foutb[0][j][5] = 0.0;
			popreal8_(&(fout[0][j][1]));
			foutb[0][j][3] = foutb[0][j][3] + foutb[0][j][1];
			rhowb = rhowb + u0 * 2 * foutb[0][j][1] / 3;
			*u0b = *u0b + 2 * rhow*foutb[0][j][1] / 3;
			foutb[0][j][1] = 0.0;
			popreal8_(&rhow);
			tempb = rhowb / (1 - u0);
			foutb[0][j][0] = foutb[0][j][0] + tempb;
			foutb[0][j][2] = foutb[0][j][2] + tempb;
			foutb[0][j][4] = foutb[0][j][4] + tempb;
			foutb[0][j][3] = foutb[0][j][3] + 2 * tempb;
			foutb[0][j][6] = foutb[0][j][6] + 2 * tempb;
			foutb[0][j][7] = foutb[0][j][7] + 2 * tempb;
			*u0b = *u0b + (fout[0][j][0] + fout[0][j][2] + fout[0][j][4] +
				2 * fout[0][j][3] + 2 * fout[0][j][6] + 2 * fout[0][j][7])*tempb / (1 - u0);
		}
	}
	for (int i = n - 1; i > -1; --i)
		for (int j = m - 1; j > -1; --j)
			for (int k = 8; k > -1; --k)
				foutb[i][j][k] = 0.0;
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

///******************************************************************
///******************************************************************
///******************************************************************

void dfds_d(double p, double u0, double u0d, int m, int n, double *cx,
	double *cy, double *w, double **rho, double **u, double **v, double
	omega, double ***feq, double ***feqd, double ***fin, double ***fout,
	double ***foutd) {
	//Friend function
	int m0, n0;
	m0 = p * m;
	n0 = m;
	//Collision step
	double temp1, temp2, rhow, ssum, usum, vsum;
	double rhowd;
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
	for (int i = 0; i < n; ++i)
		for (int j = 0; j < m; ++j) {
			temp1 = u[i][j] * u[i][j] + v[i][j] * v[i][j];
			for (int k = 0; k < 9; ++k) {
				temp2 = u[i][j] * cx[k] + v[i][j] * cy[k];
				feq[k][i][j] = rho[i][j] * w[k] * (1 + 3 * temp2 + 4.5*temp2*temp2 - 1.5*temp1);
				foutd[i][j][k] = 0.0;
				fout[i][j][k] = omega * feq[k][i][j] + (1 - omega)*fin[i][j][k];
			}
		}
	//***foutd = 0.0;
	//Apply BC
	//Inlet
	for (int j = m0; j < m; ++j) {
		rhowd = ((foutd[0][j][0] + foutd[0][j][2] + foutd[0][j][4] + 2 * foutd[0][j][3
			] + 2 * foutd[0][j][6] + 2 * foutd[0][j][7])*(1 - u0) + (fout[0][j][0] + fout[0]
				[j][2] + fout[0][j][4] + 2 * (fout[0][j][3] + fout[0][j][6] + fout[0][j][7])
				)*u0d) / ((1 - u0)*(1 - u0));
		rhow = (fout[0][j][0] + fout[0][j][2] + fout[0][j][4] + 2 * (fout[0][j][3] +
			fout[0][j][6] + fout[0][j][7])) / (1 - u0);
		foutd[0][j][1] = foutd[0][j][3] + 2 * (rhowd*u0 + rhow * u0d) / 3;
		fout[0][j][1] = fout[0][j][3] + 2 * rhow*u0 / 3;
		foutd[0][j][5] = foutd[0][j][7] + (rhowd*u0 + rhow * u0d) / 6;
		fout[0][j][5] = fout[0][j][7] + rhow * u0 / 6;
		foutd[0][j][8] = foutd[0][j][6] + (rhowd*u0 + rhow * u0d) / 6;
		fout[0][j][8] = fout[0][j][6] + rhow * u0 / 6;
	}
	//South
	for (int i = 0; i < n; ++i) {
		foutd[i][0][2] = foutd[i][0][4];
		fout[i][0][2] = fout[i][0][4];
		foutd[i][0][5] = foutd[i][0][7];
		fout[i][0][5] = fout[i][0][7];
		foutd[i][0][6] = foutd[i][0][8];
		fout[i][0][6] = fout[i][0][8];
	}
	//North
	for (int i = 0; i < n; ++i) {
		foutd[i][m - 1][4] = foutd[i][m - 1][2];
		fout[i][m - 1][4] = fout[i][m - 1][2];
		foutd[i][m - 1][8] = foutd[i][m - 1][6];
		fout[i][m - 1][8] = fout[i][m - 1][6];
		foutd[i][m - 1][7] = foutd[i][m - 1][5];
		fout[i][m - 1][7] = fout[i][m - 1][5];
	}
	//Outlet - extrapolation
	for (int j = 0; j < m; ++j) {
		foutd[n - 1][j][1] = 2 * foutd[n - 2][j][1] - foutd[n - 3][j][1];
		fout[n - 1][j][1] = 2 * fout[n - 2][j][1] - fout[n - 3][j][1];
		foutd[n - 1][j][5] = 2 * foutd[n - 2][j][5] - foutd[n - 3][j][5];
		fout[n - 1][j][5] = 2 * fout[n - 2][j][5] - fout[n - 3][j][5];
		foutd[n - 1][j][8] = 2 * foutd[n - 2][j][8] - foutd[n - 3][j][8];
		fout[n - 1][j][8] = 2 * fout[n - 2][j][8] - fout[n - 3][j][8];
	}
	//Back-facing step
	for (int i = 0; i < n0; ++i) {
		foutd[i][m0 - 1][2] = foutd[i][m0 - 1][4];
		fout[i][m0 - 1][2] = fout[i][m0 - 1][4];
		foutd[i][m0 - 1][5] = foutd[i][m0 - 1][7];
		fout[i][m0 - 1][5] = fout[i][m0 - 1][7];
		foutd[i][m0 - 1][6] = foutd[i][m0 - 1][8];
		fout[i][m0 - 1][6] = fout[i][m0 - 1][8];
	}
	for (int j = 0; j < m0; ++j) {
		foutd[n0 - 1][j][1] = foutd[n0 - 1][j][3];
		fout[n0 - 1][j][1] = fout[n0 - 1][j][3];
		foutd[n0 - 1][j][5] = foutd[n0 - 1][j][7];
		fout[n0 - 1][j][5] = fout[n0 - 1][j][7];
		foutd[n0 - 1][j][8] = foutd[n0 - 1][j][6];
		fout[n0 - 1][j][8] = fout[n0 - 1][j][6];
		//???
	}
	//Streaming step
	for (int j = 0; j < m; ++j) {
		for (int i = n - 1; i > 0; --i) {
			foutd[i][j][1] = foutd[i - 1][j][1];
			fout[i][j][1] = fout[i - 1][j][1];
		}
		for (int i = 0; i < n - 1; ++i) {
			foutd[i][j][3] = foutd[i + 1][j][3];
			fout[i][j][3] = fout[i + 1][j][3];
		}
	}
	for (int j = m - 1; j > 0; --j) {
		for (int i = 0; i < n; ++i) {
			foutd[i][j][2] = foutd[i][j - 1][2];
			fout[i][j][2] = fout[i][j - 1][2];
		}
		for (int i = n - 1; i > 0; --i) {
			foutd[i][j][5] = foutd[i - 1][j - 1][5];
			fout[i][j][5] = fout[i - 1][j - 1][5];
		}
		for (int i = 0; i < n - 1; ++i) {
			foutd[i][j][6] = foutd[i + 1][j - 1][6];
			fout[i][j][6] = fout[i + 1][j - 1][6];
		}
	}
	for (int j = 0; j < m - 1; ++j) {
		for (int i = 0; i < n; ++i) {
			foutd[i][j][4] = foutd[i][j + 1][4];
			fout[i][j][4] = fout[i][j + 1][4];
		}
		for (int i = 0; i < n - 1; ++i) {
			foutd[i][j][7] = foutd[i + 1][j + 1][7];
			fout[i][j][7] = fout[i + 1][j + 1][7];
		}
		for (int i = n - 1; i > 0; --i) {
			foutd[i][j][8] = foutd[i - 1][j + 1][8];
			fout[i][j][8] = fout[i - 1][j + 1][8];
		}
	}
	//Calculate macroscopic
	for (int j = 0; j < m; ++j)
		for (int i = 0; i < n; ++i) {
			ssum = 0;
			for (int k = 0; k < 9; ++k)
				ssum = ssum + fout[i][j][k];
			rho[i][j] = ssum;
		}
	for (int i = 0; i < n; ++i)
		rho[i][m - 1] = fout[i][m - 1][0] + fout[i][m - 1][1] + fout[i][m - 1]
		[3] + 2 * (fout[i][m - 1][2] + fout[i][m - 1][6] + fout[i][m - 1][5]);
	for (int i = 0; i < n; ++i)
		for (int j = 0; j < m - 1; ++j) {
			usum = 0;
			vsum = 0;
			for (int k = 0; k < 9; ++k) {
				usum = usum + fout[i][j][k] * cx[k];
				vsum = vsum + fout[i][j][k] * cy[k];
			}
			u[i][j] = usum / rho[i][j];
			v[i][j] = vsum / rho[i][j];
		}
	for (int j = 0; j < m; ++j)
		v[n - 1][j] = 0;
	for (int j = 0; j < m0; ++j)
		for (int i = 0; i < n0; ++i) {
			u[i][j] = 0;
			v[i][j] = 0;
		}
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

///******************************************************************
///******************************************************************
///******************************************************************

void dfdu_b(double p, double u0, int m, int n, double *cx, double *cy
	, double *w, double **rho, double **u, double **v, double omega,
	double ***feq, double ***feqb, double ***fin, double ***finb, double ***fout, double ***foutb) {
	//Friend function
	int m0, n0;
	m0 = p * m;
	n0 = m;
	//Collision step
	double temp1, temp2, rhow, ssum, usum, vsum;
	double rhowb;
	double*** foutbBuffer = new double**[n];
	for (int i = 0; i < n; ++i)
	{
		foutbBuffer[i] = new double*[m];

		for (int j = 0; j < m; ++j)
		{
			foutbBuffer[i][j] = new double[9];

			for (int k = 0; k < 9; ++k)
			{
				foutbBuffer[i][j][k] = foutb[i][j][k];
			}
		}
	}
	{
		double tmp;
		double tmp0;
		double tmp1;
	}
	//Streaming step
	for (int j = 0; j < m; ++j) {
		{
			double tmp;
		}
		{
			double tmp;
		}
	}
	for (int j = m - 1; j > 0; --j) {
		{
			double tmp;
		}
		{
			double tmp;
		}
		{
			double tmp;
		}
	}
	for (int j = 0; j < m - 1; ++j) {
		{
			double tmp;
		}
		{
			double tmp;
		}
		{
			double tmp;
		}
	}
	pushinteger4_(n0);
	pushinteger4_(m0);
	popinteger4_(&m0);
	popinteger4_(&n0);
	for (int j = m - 2; j > -1; --j) {
		{
			double tmpb;
			for (int i = 1; i < n; ++i) {
				tmpb = foutb[i][j][8];
				foutb[i][j][8] = 0.0;
				foutb[i - 1][j + 1][8] = foutb[i - 1][j + 1][8] + tmpb;
			}
		}
		{
			double tmpb;
			for (int i = n - 2; i > -1; --i) {
				tmpb = foutb[i][j][7];
				foutb[i][j][7] = 0.0;
				foutb[i + 1][j + 1][7] = foutb[i + 1][j + 1][7] + tmpb;
			}
		}
		{
			double tmpb;
			for (int i = n - 1; i > -1; --i) {
				tmpb = foutb[i][j][4];
				foutb[i][j][4] = 0.0;
				foutb[i][j + 1][4] = foutb[i][j + 1][4] + tmpb;
			}
		}
	}
	for (int j = 1; j < m; ++j) {
		{
			double tmpb;
			for (int i = n - 2; i > -1; --i) {
				tmpb = foutb[i][j][6];
				foutb[i][j][6] = 0.0;
				foutb[i + 1][j - 1][6] = foutb[i + 1][j - 1][6] + tmpb;
			}
		}
		{
			double tmpb;
			for (int i = 1; i < n; ++i) {
				tmpb = foutb[i][j][5];
				foutb[i][j][5] = 0.0;
				foutb[i - 1][j - 1][5] = foutb[i - 1][j - 1][5] + tmpb;
			}
		}
		{
			double tmpb;
			for (int i = n - 1; i > -1; --i) {
				tmpb = foutb[i][j][2];
				foutb[i][j][2] = 0.0;
				foutb[i][j - 1][2] = foutb[i][j - 1][2] + tmpb;
			}
		}
	}
	for (int j = m - 1; j > -1; --j) {
		{
			double tmpb;
			for (int i = n - 2; i > -1; --i) {
				tmpb = foutb[i][j][3];
				foutb[i][j][3] = 0.0;
				foutb[i + 1][j][3] = foutb[i + 1][j][3] + tmpb;
			}
		}
		{
			double tmpb;
			for (int i = 1; i < n; ++i) {
				tmpb = foutb[i][j][1];
				foutb[i][j][1] = 0.0;
				foutb[i - 1][j][1] = foutb[i - 1][j][1] + tmpb;
			}
		}
	}
	for (int j = m0 - 1; j > -1; --j) {
		foutb[n0 - 1][j][6] = foutb[n0 - 1][j][6] + foutb[n0 - 1][j][8];
		foutb[n0 - 1][j][8] = 0.0;
		foutb[n0 - 1][j][7] = foutb[n0 - 1][j][7] + foutb[n0 - 1][j][5];
		foutb[n0 - 1][j][5] = 0.0;
		foutb[n0 - 1][j][3] = foutb[n0 - 1][j][3] + foutb[n0 - 1][j][1];
		foutb[n0 - 1][j][1] = 0.0;
	}
	for (int i = n0 - 1; i > -1; --i) {
		foutb[i][m0 - 1][8] = foutb[i][m0 - 1][8] + foutb[i][m0 - 1][6];
		foutb[i][m0 - 1][6] = 0.0;
		foutb[i][m0 - 1][7] = foutb[i][m0 - 1][7] + foutb[i][m0 - 1][5];
		foutb[i][m0 - 1][5] = 0.0;
		foutb[i][m0 - 1][4] = foutb[i][m0 - 1][4] + foutb[i][m0 - 1][2];
		foutb[i][m0 - 1][2] = 0.0;
	}
	{
		double tmpb;
		double tmpb0;
		double tmpb1;
		for (int j = m - 1; j > -1; --j) {
			tmpb = foutb[n - 1][j][8];
			foutb[n - 1][j][8] = 0.0;
			foutb[n - 2][j][8] = foutb[n - 2][j][8] + 2 * tmpb;
			foutb[n - 3][j][8] = foutb[n - 3][j][8] - tmpb;
			tmpb0 = foutb[n - 1][j][5];
			foutb[n - 1][j][5] = 0.0;
			foutb[n - 2][j][5] = foutb[n - 2][j][5] + 2 * tmpb0;
			foutb[n - 3][j][5] = foutb[n - 3][j][5] - tmpb0;
			tmpb1 = foutb[n - 1][j][1];
			foutb[n - 1][j][1] = 0.0;
			foutb[n - 2][j][1] = foutb[n - 2][j][1] + 2 * tmpb1;
			foutb[n - 3][j][1] = foutb[n - 3][j][1] - tmpb1;
		}
	}
	for (int i = n - 1; i > -1; --i) {
		foutb[i][m - 1][5] = foutb[i][m - 1][5] + foutb[i][m - 1][7];
		foutb[i][m - 1][7] = 0.0;
		foutb[i][m - 1][6] = foutb[i][m - 1][6] + foutb[i][m - 1][8];
		foutb[i][m - 1][8] = 0.0;
		foutb[i][m - 1][2] = foutb[i][m - 1][2] + foutb[i][m - 1][4];
		foutb[i][m - 1][4] = 0.0;
	}
	for (int i = n - 1; i > -1; --i) {
		foutb[i][0][8] = foutb[i][0][8] + foutb[i][0][6];
		foutb[i][0][6] = 0.0;
		foutb[i][0][7] = foutb[i][0][7] + foutb[i][0][5];
		foutb[i][0][5] = 0.0;
		foutb[i][0][4] = foutb[i][0][4] + foutb[i][0][2];
		foutb[i][0][2] = 0.0;
	}
	{
		double tempb;
		for (int j = m - 1; j > m0 - 1; --j) {
			foutb[0][j][6] = foutb[0][j][6] + foutb[0][j][8];
			rhowb = u0 * foutb[0][j][8] / 6;
			foutb[0][j][8] = 0.0;
			foutb[0][j][7] = foutb[0][j][7] + foutb[0][j][5];
			rhowb = rhowb + u0 * foutb[0][j][5] / 6;
			foutb[0][j][5] = 0.0;
			foutb[0][j][3] = foutb[0][j][3] + foutb[0][j][1];
			rhowb = rhowb + u0 * 2 * foutb[0][j][1] / 3;
			foutb[0][j][1] = 0.0;
			tempb = rhowb / (1 - u0);
			foutb[0][j][0] = foutb[0][j][0] + tempb;
			foutb[0][j][2] = foutb[0][j][2] + tempb;
			foutb[0][j][4] = foutb[0][j][4] + tempb;
			foutb[0][j][3] = foutb[0][j][3] + 2 * tempb;
			foutb[0][j][6] = foutb[0][j][6] + 2 * tempb;
			foutb[0][j][7] = foutb[0][j][7] + 2 * tempb;
		}
	}
	for (int i = n - 1; i > -1; --i)
		for (int j = m - 1; j > -1; --j)
			for (int k = 8; k > -1; --k) {
				finb[i][j][k] = 0;
				finb[i][j][k] = finb[i][j][k] + (1 - omega)*foutb[i][j][k];
				foutb[i][j][k] = 0.0;
			}
	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < m; ++j)
		{
			for (int k = 0; k < 9; ++k)
			{
				foutb[i][j][k] = foutbBuffer[i][j][k];
			}
		}
	}
	for (int iter = 0; iter < n; ++iter)
	{
		for (int j = 0; j < m; ++j)
		{
			delete[] foutbBuffer[iter][j];
		}
		delete[] foutbBuffer[iter];
	}
	delete[] foutbBuffer;
}

///******************************************************************
///******************************************************************
///******************************************************************

void dfdu_d(double p, double u0, int m, int n, double *cx, double *cy
	, double *w, double **rho, double **u, double **v, double omega,
	double ***feq, double ***feqd, double ***fin, double ***find, double ***fout, double ***foutd) {
	//Friend function
	int m0, n0;
	m0 = p * m;
	n0 = m;
	//Collision step
	double temp1, temp2, rhow, ssum, usum, vsum;
	double rhowd;
	//***foutd = 0.0;
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
	for (int i = 0; i < n; ++i)
		for (int j = 0; j < m; ++j) {
			temp1 = u[i][j] * u[i][j] + v[i][j] * v[i][j];
			for (int k = 0; k < 9; ++k) {
				temp2 = u[i][j] * cx[k] + v[i][j] * cy[k];
				feq[k][i][j] = rho[i][j] * w[k] * (1 + 3 * temp2 + 4.5*temp2*temp2 - 1.5*temp1);
				foutd[i][j][k] = (1 - omega)*find[i][j][k];
				fout[i][j][k] = omega * feq[k][i][j] + (1 - omega)*fin[i][j][k];
			}
		}
	//Apply BC
	//Inlet
	for (int j = m0; j < m; ++j) {
		rhowd = (foutd[0][j][0] + foutd[0][j][2] + foutd[0][j][4] +
			2 * foutd[0][j][3] + 2 * foutd[0][j][6] + 2 * foutd[0][j][7]) / (1 - u0);
		rhow = (fout[0][j][0] + fout[0][j][2] + fout[0][j][4] +
			2 * (fout[0][j][3] + fout[0][j][6] + fout[0][j][7])) / (1 - u0);
		foutd[0][j][1] = foutd[0][j][3] + 2 * u0*rhowd / 3;
		fout[0][j][1] = fout[0][j][3] + 2 * rhow*u0 / 3;
		foutd[0][j][5] = foutd[0][j][7] + u0 * rhowd / 6;
		fout[0][j][5] = fout[0][j][7] + rhow * u0 / 6;
		foutd[0][j][8] = foutd[0][j][6] + u0 * rhowd / 6;
		fout[0][j][8] = fout[0][j][6] + rhow * u0 / 6;
	}
	//South
	for (int i = 0; i < n; ++i) {
		foutd[i][0][2] = foutd[i][0][4];
		fout[i][0][2] = fout[i][0][4];
		foutd[i][0][5] = foutd[i][0][7];
		fout[i][0][5] = fout[i][0][7];
		foutd[i][0][6] = foutd[i][0][8];
		fout[i][0][6] = fout[i][0][8];
	}
	//North
	for (int i = 0; i < n; ++i) {
		foutd[i][m - 1][4] = foutd[i][m - 1][2];
		fout[i][m - 1][4] = fout[i][m - 1][2];
		foutd[i][m - 1][8] = foutd[i][m - 1][6];
		fout[i][m - 1][8] = fout[i][m - 1][6];
		foutd[i][m - 1][7] = foutd[i][m - 1][5];
		fout[i][m - 1][7] = fout[i][m - 1][5];
	}
	//Outlet - extrapolation
	for (int j = 0; j < m; ++j) {
		foutd[n - 1][j][1] = 2 * foutd[n - 2][j][1] - foutd[n - 3][j][1];
		fout[n - 1][j][1] = 2 * fout[n - 2][j][1] - fout[n - 3][j][1];
		foutd[n - 1][j][5] = 2 * foutd[n - 2][j][5] - foutd[n - 3][j][5];
		fout[n - 1][j][5] = 2 * fout[n - 2][j][5] - fout[n - 3][j][5];
		foutd[n - 1][j][8] = 2 * foutd[n - 2][j][8] - foutd[n - 3][j][8];
		fout[n - 1][j][8] = 2 * fout[n - 2][j][8] - fout[n - 3][j][8];
	}
	//Back-facing step
	for (int i = 0; i < n0; ++i) {
		foutd[i][m0 - 1][2] = foutd[i][m0 - 1][4];
		fout[i][m0 - 1][2] = fout[i][m0 - 1][4];
		foutd[i][m0 - 1][5] = foutd[i][m0 - 1][7];
		fout[i][m0 - 1][5] = fout[i][m0 - 1][7];
		foutd[i][m0 - 1][6] = foutd[i][m0 - 1][8];
		fout[i][m0 - 1][6] = fout[i][m0 - 1][8];
	}
	for (int j = 0; j < m0; ++j) {
		foutd[n0 - 1][j][1] = foutd[n0 - 1][j][3];
		fout[n0 - 1][j][1] = fout[n0 - 1][j][3];
		foutd[n0 - 1][j][5] = foutd[n0 - 1][j][7];
		fout[n0 - 1][j][5] = fout[n0 - 1][j][7];
		foutd[n0 - 1][j][8] = foutd[n0 - 1][j][6];
		fout[n0 - 1][j][8] = fout[n0 - 1][j][6];
		//???
	}
	//Streaming step
	for (int j = 0; j < m; ++j) {
		for (int i = n - 1; i > 0; --i) {
			foutd[i][j][1] = foutd[i - 1][j][1];
			fout[i][j][1] = fout[i - 1][j][1];
		}
		for (int i = 0; i < n - 1; ++i) {
			foutd[i][j][3] = foutd[i + 1][j][3];
			fout[i][j][3] = fout[i + 1][j][3];
		}
	}
	for (int j = m - 1; j > 0; --j) {
		for (int i = 0; i < n; ++i) {
			foutd[i][j][2] = foutd[i][j - 1][2];
			fout[i][j][2] = fout[i][j - 1][2];
		}
		for (int i = n - 1; i > 0; --i) {
			foutd[i][j][5] = foutd[i - 1][j - 1][5];
			fout[i][j][5] = fout[i - 1][j - 1][5];
		}
		for (int i = 0; i < n - 1; ++i) {
			foutd[i][j][6] = foutd[i + 1][j - 1][6];
			fout[i][j][6] = fout[i + 1][j - 1][6];
		}
	}
	for (int j = 0; j < m - 1; ++j) {
		for (int i = 0; i < n; ++i) {
			foutd[i][j][4] = foutd[i][j + 1][4];
			fout[i][j][4] = fout[i][j + 1][4];
		}
		for (int i = 0; i < n - 1; ++i) {
			foutd[i][j][7] = foutd[i + 1][j + 1][7];
			fout[i][j][7] = fout[i + 1][j + 1][7];
		}
		for (int i = n - 1; i > 0; --i) {
			foutd[i][j][8] = foutd[i - 1][j + 1][8];
			fout[i][j][8] = fout[i - 1][j + 1][8];
		}
	}
	//Calculate macroscopic
	for (int j = 0; j < m; ++j)
		for (int i = 0; i < n; ++i) {
			ssum = 0;
			for (int k = 0; k < 9; ++k)
				ssum = ssum + fout[i][j][k];
			rho[i][j] = ssum;
		}
	for (int i = 0; i < n; ++i)
		rho[i][m - 1] = fout[i][m - 1][0] + fout[i][m - 1][1] + fout[i][m - 1][3] +
		2 * (fout[i][m - 1][2] + fout[i][m - 1][6] + fout[i][m - 1][5]);
	for (int i = 0; i < n; ++i)
		for (int j = 0; j < m - 1; ++j) {
			usum = 0;
			vsum = 0;
			for (int k = 0; k < 9; ++k) {
				usum = usum + fout[i][j][k] * cx[k];
				vsum = vsum + fout[i][j][k] * cy[k];
			}
			u[i][j] = usum / rho[i][j];
			v[i][j] = vsum / rho[i][j];
		}
	for (int j = 0; j < m; ++j)
		v[n - 1][j] = 0;
	for (int j = 0; j < m0; ++j)
		for (int i = 0; i < n0; ++i) {
			u[i][j] = 0;
			v[i][j] = 0;
		}
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

///******************************************************************
///******************************************************************
///******************************************************************

void dObjdu_d(double p, double u0, int m, int n, double *cx, double *cy, double *w, double **rho, double **rhod, double **u, 
	double **v, double omega,double ***feq, double ***feqd, double ***fin, double ***find,double ***fout, double ***foutd, 
	double *J, double *Jd) {
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
	double pin = 0;
	double pind;
	double pout = 0;
	double poutd;
	//Time step
	int m0 = p * m;
	int n0 = m;
	//Collision step
	double temp1, temp2, rhow, ssum, usum, vsum;
	double rhowd, ssumd;
	for (int i = 0; i < n; ++i)
		for (int j = 0; j < m; ++j) {
			temp1 = u[i][j] * u[i][j] + v[i][j] * v[i][j];
			for (int k = 0; k < 9; ++k) {
				foutd[i][j][k] = 0;
				temp2 = u[i][j] * cx[k] + v[i][j] * cy[k];
				feq[k][i][j] = rho[i][j] * w[k] * (1 + 3 * temp2 + 4.5*temp2*temp2 - 1.5*temp1);
				foutd[i][j][k] = (1 - omega)*find[i][j][k];
				fout[i][j][k] = omega * feq[k][i][j] + (1 - omega)*fin[i][j][k];
			}
		}
	//Apply BC
	//Inlet
	for (int j = m0; j < m; ++j) {
		rhowd = (foutd[0][j][0] + foutd[0][j][2] + foutd[0][j][4] + 2 * foutd[0][j][3]
			+ 2 * foutd[0][j][6] + 2 * foutd[0][j][7]) / (1 - u0);
		rhow = (fout[0][j][0] + fout[0][j][2] + fout[0][j][4] + 2 * (fout[0][j][3] +
			fout[0][j][6] + fout[0][j][7])) / (1 - u0);
		foutd[0][j][1] = foutd[0][j][3] + 2 * u0*rhowd / 3;
		fout[0][j][1] = fout[0][j][3] + 2 * rhow*u0 / 3;
		foutd[0][j][5] = foutd[0][j][7] + u0 * rhowd / 6;
		fout[0][j][5] = fout[0][j][7] + rhow * u0 / 6;
		foutd[0][j][8] = foutd[0][j][6] + u0 * rhowd / 6;
		fout[0][j][8] = fout[0][j][6] + rhow * u0 / 6;
	}
	//South
	for (int i = 0; i < n; ++i) {
		foutd[i][0][2] = foutd[i][0][4];
		fout[i][0][2] = fout[i][0][4];
		foutd[i][0][5] = foutd[i][0][7];
		fout[i][0][5] = fout[i][0][7];
		foutd[i][0][6] = foutd[i][0][8];
		fout[i][0][6] = fout[i][0][8];
	}
	//North
	for (int i = 0; i < n; ++i) {
		foutd[i][m - 1][4] = foutd[i][m - 1][2];
		fout[i][m - 1][4] = fout[i][m - 1][2];
		foutd[i][m - 1][8] = foutd[i][m - 1][6];
		fout[i][m - 1][8] = fout[i][m - 1][6];
		foutd[i][m - 1][7] = foutd[i][m - 1][5];
		fout[i][m - 1][7] = fout[i][m - 1][5];
	}
	//Outlet - extrapolation
	for (int j = 0; j < m; ++j) {
		foutd[n - 1][j][1] = 2 * foutd[n - 2][j][1] - foutd[n - 3][j][1];
		fout[n - 1][j][1] = 2 * fout[n - 2][j][1] - fout[n - 3][j][1];
		foutd[n - 1][j][5] = 2 * foutd[n - 2][j][5] - foutd[n - 3][j][5];
		fout[n - 1][j][5] = 2 * fout[n - 2][j][5] - fout[n - 3][j][5];
		foutd[n - 1][j][8] = 2 * foutd[n - 2][j][8] - foutd[n - 3][j][8];
		fout[n - 1][j][8] = 2 * fout[n - 2][j][8] - fout[n - 3][j][8];
	}
	//Back-facing step
	for (int i = 0; i < n0; ++i) {
		foutd[i][m0 - 1][2] = foutd[i][m0 - 1][4];
		fout[i][m0 - 1][2] = fout[i][m0 - 1][4];
		foutd[i][m0 - 1][5] = foutd[i][m0 - 1][7];
		fout[i][m0 - 1][5] = fout[i][m0 - 1][7];
		foutd[i][m0 - 1][6] = foutd[i][m0 - 1][8];
		fout[i][m0 - 1][6] = fout[i][m0 - 1][8];
	}
	for (int j = 0; j < m0; ++j) {
		foutd[n0 - 1][j][1] = foutd[n0 - 1][j][3];
		fout[n0 - 1][j][1] = fout[n0 - 1][j][3];
		foutd[n0 - 1][j][5] = foutd[n0 - 1][j][7];
		fout[n0 - 1][j][5] = fout[n0 - 1][j][7];
		foutd[n0 - 1][j][8] = foutd[n0 - 1][j][6];
		fout[n0 - 1][j][8] = fout[n0 - 1][j][6];
		//???
	}
	//Streaming step
	for (int j = 0; j < m; ++j) {
		for (int i = n - 1; i > 0; --i) {
			foutd[i][j][1] = foutd[i - 1][j][1];
			fout[i][j][1] = fout[i - 1][j][1];
		}
		for (int i = 0; i < n - 1; ++i) {
			foutd[i][j][3] = foutd[i + 1][j][3];
			fout[i][j][3] = fout[i + 1][j][3];
		}
	}
	for (int j = m - 1; j > 0; --j) {
		for (int i = 0; i < n; ++i) {
			foutd[i][j][2] = foutd[i][j - 1][2];
			fout[i][j][2] = fout[i][j - 1][2];
		}
		for (int i = n - 1; i > 0; --i) {
			foutd[i][j][5] = foutd[i - 1][j - 1][5];
			fout[i][j][5] = fout[i - 1][j - 1][5];
		}
		for (int i = 0; i < n - 1; ++i) {
			foutd[i][j][6] = foutd[i + 1][j - 1][6];
			fout[i][j][6] = fout[i + 1][j - 1][6];
		}
	}
	for (int j = 0; j < m - 1; ++j) {
		for (int i = 0; i < n; ++i) {
			foutd[i][j][4] = foutd[i][j + 1][4];
			fout[i][j][4] = fout[i][j + 1][4];
		}
		for (int i = 0; i < n - 1; ++i) {
			foutd[i][j][7] = foutd[i + 1][j + 1][7];
			fout[i][j][7] = fout[i + 1][j + 1][7];
		}
		for (int i = n - 1; i > 0; --i) {
			foutd[i][j][8] = foutd[i - 1][j + 1][8];
			fout[i][j][8] = fout[i - 1][j + 1][8];
		}
	}
	//Calculate macroscopic
	for (int j = 0; j < m; ++j)
		for (int i = 0; i < n; ++i) {
			ssum = 0;
			ssumd = 0.0;
			rhod[i][j] = 0;
			for (int k = 0; k < 9; ++k) {
				ssumd = ssumd + foutd[i][j][k];
				ssum = ssum + fout[i][j][k];
			}
			rhod[i][j] = ssumd;
			rho[i][j] = ssum;
		}
	for (int i = 0; i < n; ++i) {
		rhod[i][m - 1] = foutd[i][m - 1][0] + foutd[i][m - 1][1] + foutd[i][m
			- 1][3] + 2 * foutd[i][m - 1][2] + 2 * foutd[i][m - 1][6] + 2 * foutd[i][m - 1][5];
		rho[i][m - 1] = fout[i][m - 1][0] + fout[i][m - 1][1] + fout[i][m - 1]
			[3] + 2 * (fout[i][m - 1][2] + fout[i][m - 1][6] + fout[i][m - 1][5]);
	}
	for (int i = 0; i < n; ++i)
		for (int j = 0; j < m - 1; ++j) {
			usum = 0;
			vsum = 0;
			for (int k = 0; k < 9; ++k) {
				usum = usum + fout[i][j][k] * cx[k];
				vsum = vsum + fout[i][j][k] * cy[k];
			}
			u[i][j] = usum / rho[i][j];
			v[i][j] = vsum / rho[i][j];
		}
	for (int j = 0; j < m; ++j)
		v[n - 1][j] = 0;
	for (int j = 0; j < m0; ++j)
		for (int i = 0; i < n0; ++i) {
			u[i][j] = 0;
			v[i][j] = 0;
		}
	pind = 0.0;
	//Pressure drop
	for (int j = m0; j < m; ++j) {
		pind = pind + rhod[0][j] / 3;
		pin += rho[0][j] / 3;
	}
	poutd = 0.0;
	for (int j = 0; j < m; ++j) {
		poutd = poutd + rhod[n - 1][j] / 3;
		pout += rho[n - 1][j] / 3;
	}
	pind = pind / (m - m0);
	pin /= m - m0;
	poutd = poutd / m;
	pout /= m;
	*Jd = (pind - poutd) / 2;
	*J = (pin - pout) / 2;

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

///******************************************************************
///******************************************************************
///******************************************************************

void dObjds_d(double p, double u0, double u0d, int m, int n, double *cx, double *cy, double *w, double **rho, double **rhod, double **u, double **v,
	double omega, double ***feq, double ***feqd, double ***fin, double ***fout, double ***foutd, double *J, double *Jd) {
	double pin = 0;
	double pind;
	double pout = 0;
	double poutd;
	//Time step
	int m0 = p * m;
	int n0 = m;
	//Collision step
	double temp1, temp2, rhow, ssum, usum, vsum;
	double rhowd, ssumd;
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
	for (int i = 0; i < n; ++i)
		for (int j = 0; j < m; ++j) {
			temp1 = u[i][j] * u[i][j] + v[i][j] * v[i][j];
			for (int k = 0; k < 9; ++k) {
				temp2 = u[i][j] * cx[k] + v[i][j] * cy[k];
				feq[k][i][j] = rho[i][j] * w[k] * (1 + 3 * temp2 + 4.5*temp2*temp2 - 1.5*temp1);
				foutd[i][j][k] = 0.0;
				fout[i][j][k] = omega * feq[k][i][j] + (1 - omega)*fin[i][j][k];
			}
		}
	//Apply BC
	//Inlet
	for (int j = m0; j < m; ++j) {
		rhowd = ((foutd[0][j][0] + foutd[0][j][2] + foutd[0][j][4] + 2 * foutd[0][j][3] + 2 * foutd[0][j][6] + 
			2 * foutd[0][j][7])*(1 - u0) + (fout[0][j][0] + fout[0][j][2] + fout[0][j][4] + 
				2 * (fout[0][j][3] + fout[0][j][6] + fout[0][j][7]))*u0d) / ((1 - u0)*(1 - u0));
		rhow = (fout[0][j][0] + fout[0][j][2] + fout[0][j][4] + 
			2 * (fout[0][j][3] +fout[0][j][6] + fout[0][j][7])) / (1 - u0);
		foutd[0][j][1] = foutd[0][j][3] + 2 * (rhowd*u0 + rhow * u0d) / 3;
		fout[0][j][1] = fout[0][j][3] + 2 * rhow*u0 / 3;
		foutd[0][j][5] = foutd[0][j][7] + (rhowd*u0 + rhow * u0d) / 6;
		fout[0][j][5] = fout[0][j][7] + rhow * u0 / 6;
		foutd[0][j][8] = foutd[0][j][6] + (rhowd*u0 + rhow * u0d) / 6;
		fout[0][j][8] = fout[0][j][6] + rhow * u0 / 6;
	}
	//South
	for (int i = 0; i < n; ++i) {
		foutd[i][0][2] = foutd[i][0][4];
		fout[i][0][2] = fout[i][0][4];
		foutd[i][0][5] = foutd[i][0][7];
		fout[i][0][5] = fout[i][0][7];
		foutd[i][0][6] = foutd[i][0][8];
		fout[i][0][6] = fout[i][0][8];
	}
	//North
	for (int i = 0; i < n; ++i) {
		foutd[i][m - 1][4] = foutd[i][m - 1][2];
		fout[i][m - 1][4] = fout[i][m - 1][2];
		foutd[i][m - 1][8] = foutd[i][m - 1][6];
		fout[i][m - 1][8] = fout[i][m - 1][6];
		foutd[i][m - 1][7] = foutd[i][m - 1][5];
		fout[i][m - 1][7] = fout[i][m - 1][5];
	}
	//Outlet - extrapolation
	for (int j = 0; j < m; ++j) {
		foutd[n - 1][j][1] = 2 * foutd[n - 2][j][1] - foutd[n - 3][j][1];
		fout[n - 1][j][1] = 2 * fout[n - 2][j][1] - fout[n - 3][j][1];
		foutd[n - 1][j][5] = 2 * foutd[n - 2][j][5] - foutd[n - 3][j][5];
		fout[n - 1][j][5] = 2 * fout[n - 2][j][5] - fout[n - 3][j][5];
		foutd[n - 1][j][8] = 2 * foutd[n - 2][j][8] - foutd[n - 3][j][8];
		fout[n - 1][j][8] = 2 * fout[n - 2][j][8] - fout[n - 3][j][8];
	}
	//Back-facing step
	for (int i = 0; i < n0; ++i) {
		foutd[i][m0 - 1][2] = foutd[i][m0 - 1][4];
		fout[i][m0 - 1][2] = fout[i][m0 - 1][4];
		foutd[i][m0 - 1][5] = foutd[i][m0 - 1][7];
		fout[i][m0 - 1][5] = fout[i][m0 - 1][7];
		foutd[i][m0 - 1][6] = foutd[i][m0 - 1][8];
		fout[i][m0 - 1][6] = fout[i][m0 - 1][8];
	}
	for (int j = 0; j < m0; ++j) {
		foutd[n0 - 1][j][1] = foutd[n0 - 1][j][3];
		fout[n0 - 1][j][1] = fout[n0 - 1][j][3];
		foutd[n0 - 1][j][5] = foutd[n0 - 1][j][7];
		fout[n0 - 1][j][5] = fout[n0 - 1][j][7];
		foutd[n0 - 1][j][8] = foutd[n0 - 1][j][6];
		fout[n0 - 1][j][8] = fout[n0 - 1][j][6];
		//???
	}
	//Streaming step
	for (int j = 0; j < m; ++j) {
		for (int i = n - 1; i > 0; --i) {
			foutd[i][j][1] = foutd[i - 1][j][1];
			fout[i][j][1] = fout[i - 1][j][1];
		}
		for (int i = 0; i < n - 1; ++i) {
			foutd[i][j][3] = foutd[i + 1][j][3];
			fout[i][j][3] = fout[i + 1][j][3];
		}
	}
	for (int j = m - 1; j > 0; --j) {
		for (int i = 0; i < n; ++i) {
			foutd[i][j][2] = foutd[i][j - 1][2];
			fout[i][j][2] = fout[i][j - 1][2];
		}
		for (int i = n - 1; i > 0; --i) {
			foutd[i][j][5] = foutd[i - 1][j - 1][5];
			fout[i][j][5] = fout[i - 1][j - 1][5];
		}
		for (int i = 0; i < n - 1; ++i) {
			foutd[i][j][6] = foutd[i + 1][j - 1][6];
			fout[i][j][6] = fout[i + 1][j - 1][6];
		}
	}
	for (int j = 0; j < m - 1; ++j) {
		for (int i = 0; i < n; ++i) {
			foutd[i][j][4] = foutd[i][j + 1][4];
			fout[i][j][4] = fout[i][j + 1][4];
		}
		for (int i = 0; i < n - 1; ++i) {
			foutd[i][j][7] = foutd[i + 1][j + 1][7];
			fout[i][j][7] = fout[i + 1][j + 1][7];
		}
		for (int i = n - 1; i > 0; --i) {
			foutd[i][j][8] = foutd[i - 1][j + 1][8];
			fout[i][j][8] = fout[i - 1][j + 1][8];
		}
	}
	for (int j = 0; j < m; ++j)
		for (int i = 0; i < n; ++i)
			rhod[i][j] = 0;
	//Calculate macroscopic
	for (int j = 0; j < m; ++j)
		for (int i = 0; i < n; ++i) {
			ssum = 0;
			ssumd = 0.0;
			for (int k = 0; k < 9; ++k) {
				ssumd = ssumd + foutd[i][j][k];
				ssum = ssum + fout[i][j][k];
			}
			rhod[i][j] = ssumd;
			rho[i][j] = ssum;
		}
	for (int i = 0; i < n; ++i) {
		rhod[i][m - 1] = foutd[i][m - 1][0] + foutd[i][m - 1][1] + foutd[i][m- 1][3] + 
			2 * foutd[i][m - 1][2] + 2 * foutd[i][m - 1][6] + 2 * foutd[i][m - 1][5];
		rho[i][m - 1] = fout[i][m - 1][0] + fout[i][m - 1][1] + fout[i][m - 1][3] + 
			2 * (fout[i][m - 1][2] + fout[i][m - 1][6] + fout[i][m - 1][5]);
	}
	for (int i = 0; i < n; ++i)
		for (int j = 0; j < m - 1; ++j) {
			usum = 0;
			vsum = 0;
			for (int k = 0; k < 9; ++k) {
				usum = usum + fout[i][j][k] * cx[k];
				vsum = vsum + fout[i][j][k] * cy[k];
			}
			u[i][j] = usum / rho[i][j];
			v[i][j] = vsum / rho[i][j];
		}
	for (int j = 0; j < m; ++j)
		v[n - 1][j] = 0;
	for (int j = 0; j < m0; ++j)
		for (int i = 0; i < n0; ++i) {
			u[i][j] = 0;
			v[i][j] = 0;
		}
	pind = 0.0;
	//Pressure drop
	for (int j = m0; j < m; ++j) {
		pind = pind + rhod[0][j] / 3;
		pin += rho[0][j] / 3;
	}
	poutd = 0.0;
	for (int j = 0; j < m; ++j) {
		poutd = poutd + rhod[n - 1][j] / 3;
		pout += rho[n - 1][j] / 3;
	}
	pind = pind / (m - m0);
	pin /= m - m0;
	poutd = poutd / m;
	pout /= m;
	*Jd = (pind - poutd) / 2;
	*J = (pin - pout) / 2;

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