#include "winbgi2.h"
#include "LSS.h"
#include "AuxFunctions.h"
#include <fstream>
#include <string>
#include <cstdlib>

double *xs; //for tests
const int N = 7;
double **C = AuxFun::new_random_matrix(N, N); //Global random matrix as a base for main matrix

int main(int argc, char* argv[])
{
	double s = 0.1;
	int tsteps = 100, m = 20, mx = 10;
	if (std::string(argv[1]) == std::string("test"))
	{
		AuxFun::Calculate();
	}
	else if (std::string(argv[1]) == std::string("SolveTrajectory"))
	{
		LSS obj;
		if (argc > 2)
		{
			s = atof(argv[2]);	tsteps = atoi(argv[3]);
			m = atoi(argv[4]);	mx = atoi(argv[5]);
			std::cout << "Number of passed parameters: " << argc << "\n";
			std::cout << "Passed parameters: u0 = " << s << "\ttsteps = " << tsteps << "\t m = " << m << "\tn/m = " << mx << "\n";
		}
		obj.CreateCase(s, tsteps, m, mx);
		obj.SolveCase();
	}
	else if (std::string(argv[1]) == std::string("ObjectiveFunction"))
	{
		std::ofstream Mfile("PressureDrop.txt");
		Mfile << "u0\tdeltaP\n";
		while (s <= 0.18)
		{
			LSS obj;
			obj.CreateCase(s, tsteps, m, mx);
			obj.SolveCase();
			Mfile << s << "\t" << obj.GetObjectiveFuntion() << "\n";
			s += 0.001;
		}
		Mfile.close();
	}
	else if (std::string(argv[1]) == std::string("ObjectiveFunctionSensitivity"))
	{
		double s = 0.01, deltas = 0.005, temp = 0;
		if (argc > 2)
		{
			tsteps = atoi(argv[2]);
			std::cout << "Number of passed parameters: " << argc << "\n";
			std::cout << "Passed parameters: tsteps = " << tsteps << "\n";
		}
		std::ofstream Mfile("Derivative.txt");
		Mfile << "u0\tdp\tdp_\tdJds\n";
		while (s <= 0.18)
		{
			temp = 0;
			for (int iter = 0; iter < 2; iter++)
			{
				LSS obj;
				obj.CreateCase(s, tsteps, m, mx);
				obj.SolveCase();
				if (iter > 0) Mfile << s << "\t" << temp <<"\t" << obj.GetObjectiveFuntion() << "\t" << (obj.GetObjectiveFuntion() - temp) / (deltas / 100) << "\n";
				temp = obj.GetObjectiveFuntion();
				s = s + deltas / 100;
			}
			s = s - 2 * deltas / 100;
			s += deltas;
		}
		Mfile.close();
	}
	else if (std::string(argv[1]) == std::string("LSS"))
	{
		LSS obj;
		if (argc > 2)
		{
			s = atof(argv[2]);	tsteps = atoi(argv[3]);
			m = atoi(argv[4]);	mx = atoi(argv[5]);
			std::cout << "Number of passed parameters: " << argc << "\n";
			std::cout << "Passed parameters: u0 = " << s << "\ttsteps = " << tsteps << "\t m = " << m << "\tn/m = " << mx << "\n";
		}
		//Solve LB case to obtain base trajectory
		obj.CreateCase(s, tsteps, m, mx);
		obj.SolveCase();
		//Solve KKT system with Shur compliment
		obj.Preprocess();
		obj.SolveKKT();
	}
	else if (std::string(argv[1]) == std::string("SensitivityLSS"))
	{
		LSS obj;
		double ****w, u0;
		int mstep, n, m, mx, Q;
		//Read KKT result for this case
		std::string line;
		std::ifstream file("KKTSolution.txt", std::ifstream::in);
		//try
		{
			for (int iter = 1; iter < 6; iter++)
			{
				getline(file, line);
				switch (iter)
				{
				case 1:
				{
					u0 = strtod(line.c_str(), nullptr);
					break;
				}
				case 2:
				{
					mstep = atoi(line.c_str());
					break;
				}
				case 3:
				{
					n = atoi(line.c_str());
					break;
				}
				case 4:
				{
					m = atoi(line.c_str());
					mx = n / m;
					break;
				}
				case 5:
				{
					Q = atoi(line.c_str());
					std::cout << "Case parameters: \nu0 = " << u0 << "\nmstep = " << mstep << "\nm = " << m << "\nn = " << n << "\nQ = " << Q << std::endl;
					break;
				}
				}
			}
			//Create solution vectors
			std::cout << "CRETING VECTOR FOR KKT SOLUTION\n";
			auto start1 = std::chrono::steady_clock::now();
			w = new double ***[mstep - 1];
			for (int i = 0; i < mstep - 1; i++)
			{
				w[i] = new double **[n];
				for (int j = 0; j < n; j++)
				{
					w[i][j] = new double *[m];
					for (int k = 0; k < m; k++)
					{
						w[i][j][k] = new double[Q];
						for (int l = 0; l < Q; l++)
							w[i][j][k][l] = 0;
					}
				}
				std::cout << "#";
			}
			std::cout << "\nCRETING VECTOR FOR KKT SOLUTION DONE\n";
			auto end1 = std::chrono::steady_clock::now();
			std::cout << "ELAPSED TIME: " << std::chrono::duration_cast<std::chrono::seconds>(end1 - start1).count() << std::endl;

			std::cout << "\nREADING KKT SOLUTION\n";
			auto start2 = std::chrono::steady_clock::now();
			for (int i = 0; i < mstep - 1; i++)
			{
				for (int j = 0; j < n; j++)
				{
					for (int k = 0; k < m; k++)
					{
						for (int l = 0; l < Q; l++)
						{
							getline(file, line);
							w[i][j][k][l] = strtod(line.c_str(), nullptr);
							if (w[i][j][k][l] != w[i][j][k][l]) std::cerr << "NaN in KKT solution.";
						}

					}
				}
				std::cout << "#";
			}
			file.close();
			auto end2 = std::chrono::steady_clock::now();
			std::cout << "\nREADING KKT SOLUTION DONE\n";
			std::cout << "ELAPSED TIME: " << std::chrono::duration_cast<std::chrono::seconds>(end2 - start2).count() << std::endl;
		}
		//catch (...)
		//{
			//std::cerr << "Couldn't open the file";
		//}
		//Solve LB case to obtain base trajectory
		obj.CreateCase(u0, mstep, m, mx);
		obj.SolveCase();
		obj.Preprocess();
		obj.CalculateSensitivity(w);
	}
	return 0;
}
