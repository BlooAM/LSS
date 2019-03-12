#include "winbgi2.h"
#include "LSS.h"
#include "AuxFunctions.h"

double *xs; //for tests

int main(int argc, char* argv[])
{
	if (std::string(argv[1]) == std::string("test"))
	{
		AuxFun::Calculate();
	}
	else
	if (std::string(argv[1]) == std::string("LSS"))
	{
		LSS obj;
		obj.CreateCase();
		//obj.SolveCase();

	}

	
	return 0;
}
