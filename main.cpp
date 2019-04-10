#include "winbgi2.h"
#include "LSS.h"
#include "AuxFunctions.h"
#include <fstream>

double *xs; //for tests
const int N = 7;
double **C = AuxFun::new_random_matrix(N, N); //Global random matrix as a base for main matrix

int main(int argc, char* argv[])
{
	if (std::string(argv[1]) == std::string("test"))
	{
		AuxFun::Calculate();
	}
	else
	if (std::string(argv[1]) == std::string("LSS"))
	{
		double s = 0.1;
		LSS obj;
		obj.CreateCase();
		obj.SetCase(s);
		obj.SolveCase();
	}
	else
	if (std::string(argv[1]) == std::string("LSS_loop"))
	{
		double s = 0.01;
		std::ofstream Mfile("PressureDrop.txt");
		Mfile << "u0\tdeltaP\n";
		while (s <= 0.18)
		{
			LSS obj;
			obj.CreateCase();
			obj.SetCase(s);
			obj.SolveCase();
			Mfile << s << "\t" << obj.GetObjectiveFuntion() << "\n";
			s += 0.001;
		}
		Mfile.close();
	}
	else
	if (std::string(argv[1]) == std::string("LSS_derivative"))
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
				obj.CreateCase();
				obj.SetCase(s);
				obj.SolveCase();
				if (iter > 0) Mfile << s << "\t" << (obj.GetObjectiveFuntion() - temp) / (deltas/100) << "\n";
				temp = obj.GetObjectiveFuntion();
				s = s + deltas / 100;
			}
			s = s - 2*deltas / 100;
			s += deltas;
		}
		Mfile.close();
	}

	return 0;
}
