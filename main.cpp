#include "winbgi2.h"
#include "LSS.h"
#include "AuxFunctions.h"

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
		double s = 0.5;
		LSS obj;
		obj.CreateCase();
		obj.SetCase(s);
		obj.SolveCase();

	}

	
	return 0;
}
