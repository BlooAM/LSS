#include "winbgi2.h"
#include "LSS.h"
#include "AuxFunctions.h"
#include <fstream>

double *xs; //for tests
const int N = 7;
double **C = AuxFun::new_random_matrix(N, N); //Global random matrix as a base for main matrix

int main(int argc, char* argv[])
{
	double s = 0.1;
	int tsteps = 50, m = 20, mx = 10;
	if (std::string(argv[1]) == std::string("test"))
	{
		AuxFun::Calculate();
	}
	else if (std::string(argv[1]) == std::string("SolveTrajectory"))
	{
		LSS obj;
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
		double s = 0.05, deltas = 0.005, temp = 0;
		std::ofstream Mfile("Derivative.txt");
		Mfile << "u0\tdJds\n";
		while (s <= 0.18)
		{
			temp = 0;
			for (int iter = 0; iter < 2; iter++)
			{
				LSS obj;
				obj.CreateCase(s, tsteps, m, mx);
				obj.SolveCase();
				if (iter > 0) Mfile << s << "\t" << (obj.GetObjectiveFuntion() - temp) / (deltas / 100) << "\n";
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
		//Solve LB case to obtain base trajectory
		obj.CreateCase(s,tsteps,m,mx);
		obj.SolveCase();
		//Solve KKT system with Shur compliment
		obj.Preprocess();
		obj.SolveKKT();
	}

	return 0;
}
