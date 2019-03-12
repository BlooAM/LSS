#include "LSS.h"



LSS::LSS()
{
	m = n = 0;
	ref_trajectory = nullptr;
}


LSS::~LSS()
{
	delete ref_trajectory;
}

