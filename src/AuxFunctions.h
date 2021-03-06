#pragma once
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <iomanip>

extern double *xs; //for tests
extern const int N;
extern double **C;

namespace AuxFun
{
	//Array dynamic allocation
	double ** new_random_matrix(int n, int m) {
		double ** tab;
		srand(time(NULL));
		tab = (double **)malloc(n * sizeof(double*));
		tab[0] = (double *)malloc(n*m * sizeof(double));
		for (int i = 0; i < n; i++) {
			tab[i] = &(tab[0][i*m]);
			for (int j = 0; j < m; j++) {
				tab[i][j] = (double)rand() / RAND_MAX;;
			}
		}
		return tab;
	}

	//Gauss elimination
	void Gauss(int n, double **M, double *f, double *x)
	{
		int i, j, k;
		double s;
		double wsp;

		double **a;
		double *b;

		b = new double[n];
		a = new double*[n];
		for (int i = 0; i < n; ++i)
			a[i] = new double[n];

		for (int i = 0; i < n; ++i)
		{
			for (int j = 0; j < n; ++j)
			{
				a[i][j] = M[i][j];
			}
			b[i] = f[i];
		}

		for (i = 0; i < n - 1; i++)
		{
			for (k = i + 1; k < n; k++)
			{
				wsp = a[k][i] / a[i][i];
				for (j = i + 1; j < n; j++)
					a[k][j] -= a[i][j] * wsp;
				b[k] -= b[i] * wsp;
			}
		}
		x[n - 1] = b[n - 1] / a[n - 1][n - 1];

		for (i = n - 1; i >= 0; i--)
		{
			s = 0.0;
			for (j = i + 1; j < n; j++)
				s += a[i][j] * x[j];

			x[i] = (b[i] - s) / a[i][i];
		}

		for (int i = 0; i < n; ++i)
			delete[] a[i];
		delete[] a;
		delete[] b;
	}
	
	//Vector dynamic allocation
	double * new_vector(int n) {
		double * tab;
		tab = (double *)malloc(n * sizeof(double));
		for (int j = 0; j < n; j++) {
			tab[j] = 0.0;
		}
		return tab;
	}

	void free_vector(double * tab) {
		free(tab);
	}

	//Array dynamic allocation
	double ** new_matrix(int n, int m) {
		double ** tab;
		tab = (double **)malloc(n * sizeof(double*));
		tab[0] = (double *)malloc(n*m * sizeof(double));
		for (int i = 0; i < n; i++) {
			tab[i] = &(tab[0][i*m]);
			for (int j = 0; j < m; j++) {
				tab[i][j] = 0.0;
			}
		}
		return tab;
	}

	void free_matrix(double ** tab) {
		free(tab[0]);
		free(tab);
	}

	//Multiply matrix by vector
	void MatMult(int n, double **A, double *x, double *y)
	{
		for (int i = 0; i < n; i++)
		{
			y[i] = 0.0;
			for (int j = 0; j < n; j++)
			{
				y[i] += A[i][j] * x[j];
			}
		}

	}

	//Display linear system
	void dispSystem(int n, double **A, double*b)
	{
		std::cout << "A = " << std::endl;
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++)
			{
				std::cout << A[i][j] << "\t";
				if (j == n - 1)	std::cout << std::endl;
			}
		}
		std::cout << std::endl << "b = " << std::endl;;
		for (int i = 0; i < n; i++)
		{
			std::cout << b[i] << std::endl;
		}
	}

	//Calculate scalar product
	double skal(int n, double*a, double*b)
	{
		double result = 0;
		for (int i = 0; i < n; i++)	result += a[i] * b[i];
		return result;
	}

	//Calculate vector norm
	double norm(double *r, int n)
	{
		double res = 0;
		for (int i = 0; i < n; i++)
		{
			res += r[i] * r[i];
		}
		return sqrt(res);
	}

	//Preconditioner funtions
	void Precond(int n, double **A, double*r, double*p)
	{
		double temp;
		for (int i = 0; i < n; i++)
		{
			temp = 0;
			//Gauss-Seidel
			//for (int j = 0; j < i; j++) temp += A[i][j] * p[j];
			p[i] = (r[i] - temp) / A[i][i];
		}
	}
	void FastPrecond(int n, double *diagA, double*r, double*p)
	{
		double temp;
		for (int i = 0; i < n; i++)
		{
			temp = 0;
			//Gauss-Seidel
			//for (int j = 0; j < i; j++) temp += A[i][j] * p[j];
			p[i] = (r[i] - temp) / diagA[i];
		}
	}


	void FastSolve(int n, double**A, double*b, double*x)
	{
		int d = 6, D = pow(10.0, double(d / 2));
		int maxIter = pow(10.0, double(d));
		double eps = 1e-6, alfa = 0.5, beta;
		double *temp = new_vector(n), *z = new_vector(n), *z_n = new_vector(n);
		double *p = new_vector(n), *res = new_vector(n), *res_n = new_vector(n);
		std::cout << std::endl << "ITERATIVE SECTION" << std::endl;
		MatMult(n, A, x, res); //r = Ax
		for (int i = 0; i < n; i++) res[i] = b[i] - res[i]; //r = b - r
		Precond(n, A, res, z); //calculate p
		for (int i = 0; i < n; i++) p[i] = z[i];

		for (int iter = 0; iter < maxIter; iter++)
		{
			MatMult(n, A, p, temp); //calculate Ap
			alfa = skal(n, res, z) / skal(n, p, temp);
			for (int i = 0; i < n; i++)
			{
				x[i] = x[i] + alfa * p[i]; //updated value of x
				res_n[i] = res[i] - alfa * temp[i];
			}
			if (norm(res_n, n) < eps)	break;
			Precond(n, A, res_n, z_n);
			beta = skal(n, z_n, res_n) / skal(n, z, res);
			for (int i = 0; i < n; i++)
			{
				p[i] = z_n[i] + beta * p[i];
				z[i] = z_n[i];	res[i] = res_n[i];	//update res and z to next iter
			}
			if (iter%D == 0)	std::cout << "Residuum in " << iter + 1 << " iter = " << norm(res, n)
				<< "\t alfa = " << alfa << "\t beta = " << beta << std::endl;
		}
		std::cout << std::setprecision(10);
		std::cout << "Iterative result: " << "\t Gauss result: " << "\t Realtive error" << std::endl;;
		for (int i = 0; i < n; i++)
		{
			std::cout << x[i] << "     " << xs[i] << "     " << fabs((xs[i] - x[i]) / xs[i]) << std::endl;
		}
		free_vector(p);	free_vector(temp);
		free_vector(res);	free_vector(res_n);
		free_vector(z);	free_vector(z_n);
	}

	//Solve using main array as a function
	void Solve(int n, void (*mult)(double*,double*,double*,int), double*b, double*x)
	{
		int d = 6, D = pow(10.0, double(d / 2));
		int maxIter = pow(10.0, double(d));
		double eps = 1e-6, alfa = 0.5, beta;
		double *temp = new_vector(n), *res_b = new_vector(n), *res_b_n = new_vector(n);
		double *p = new_vector(n), *res = new_vector(n), *res_n = new_vector(n);
		double *diagA = new_vector(n);
		std::cout << std::endl << "ITERATIVE SECTION" << std::endl;
		mult(x, res, diagA, n); //r = Ax
		for (int i = 0; i < n; i++) res[i] = b[i] - res[i]; //r = b - r
		FastPrecond(n, diagA, res, res_b); //calculate zs
		for (int i = 0; i < n; i++) p[i] = res_b[i];
		int iter;
		for (iter = 0; iter < maxIter; iter++)
		{
			mult(p, temp, diagA, n); //calculate temp = Ap
			alfa = skal(n, res, res_b) / skal(n, p, temp);
			for (int i = 0; i < n; i++)
			{
				x[i] = x[i] + alfa * p[i]; //updated value of x
				res_n[i] = res[i] - alfa * temp[i];
			}
			if (norm(res_n, n) < eps)	break;
			FastPrecond(n, diagA, res_n, res_b_n);
			beta = skal(n, res_b_n, res_n) / skal(n, res_b, res);
			for (int i = 0; i < n; i++)
			{
				p[i] = res_b_n[i] + beta * p[i];
				res_b[i] = res_b_n[i];	res[i] = res_n[i];	//update res and z to next iter
			}
			if (iter%D == 0)	std::cout << "Residuum in " << iter + 1 << " iter = " << norm(res, n)
				<< "\t alfa = " << alfa << "\t beta = " << beta << std::endl;
		}
		std::cout << "Total number of iterations: " << iter << std::endl;
		std::cout << std::setprecision(10);
		std::cout << "Iterative result: " << "\t Gauss result: " << "\t Realtive error" << std::endl;;
		for (int i = 0; i < n; i++)
		{
			std::cout << x[i] << "     " << xs[i] << "     " << fabs((xs[i] - x[i]) / xs[i]) << std::endl;
		}
		free_vector(p);	free_vector(temp);
		free_vector(res);	free_vector(res_n);
		free_vector(res_b);	free_vector(res_b_n);
		free_vector(diagA);
	}

	//Function for assembling array
	void AssemblyArray(double *x, double *y, double *d, int n)
	{
		double temp = 0;
		for (int i = 0; i < n; i++)
		{
			y[i] = 0;
			for (int j = 0; j < n; j++)
			{
				temp = 0;
				for (int k = 0; k < n; k++)	temp += C[i][k] * C[j][k];
				temp = temp / 10;
				if (i == j)
				{
					temp *= pow(10.0, 2.0);
					d[i] = temp;
				}
				y[i] += temp * x[j];
			}
		}
	}

	//Main function
	void Calculate()
	{

		int n = N;
		double **A, **B, *b, *x, *y, *res, *temp;
		double alfa = 1.1;

		//Construct matrix
		A = new_matrix(n, n);	B = new_matrix(n, n);	b = new_vector(n);
		x = new_vector(n);	xs = new_vector(n);	 y = new_vector(n);
		temp = new_vector(n);	res = new_vector(n);
		//srand(time(NULL));

		/*for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++)
			{
				B[i][j] = (double)rand() / RAND_MAX;
			}
		}*/
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++)
			{
				A[i][j] = 0;
				for (int k = 0; k < n; k++)	A[i][j] += C[i][k] * C[j][k];
				A[i][j] = A[i][j] / 10;
			}
		}
		for (int i = 0; i < n; i++)
		{
			A[i][i] *= pow(10.0, 2.0);
			b[i] = i;
			x[i] = 0;	y[i] = 0;
			res[i] = 0;	temp[i] = 0;
		}
		if (n < 8)	dispSystem(n, A, b);


		std::cout << "GAUSS SECTION" << std::endl;
		MatMult(n, A, x, y);
		for (int i = 0; i < n; i++) res[i] = b[i] - y[i];
		std::cout << std::endl << "Norm of res (before)= " << norm(res, n);
		Gauss(n, A, b, x);
		MatMult(n, A, x, y);
		for (int i = 0; i < n; i++) res[i] = b[i] - y[i];
		std::cout << std::endl << "Norm of res (after)= " << norm(res, n) << std::endl;
		std::cout << "Result: " << std::endl;
		for (int i = 0; i < n; i++)
		{
			xs[i] = x[i];
			x[i] = temp[i];
		}
		system("pause");
		//FastSolve(n, A, b, x);
		Solve(n, AssemblyArray, b, x);
		free_matrix(A); free_matrix(B); free_vector(b);
		free_vector(x);	free_vector(xs);	free_vector(y);
		free_vector(res);	free_vector(temp);
		free_matrix(C);
	}
	
} //AuxFun namespace