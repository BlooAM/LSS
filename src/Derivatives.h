#pragma once
#include "LBM.h"
#include "adBuffer.h"

//Dynamic derivatives
void dfds_b(float p, float u0, float *u0b, int m, int n, float *cx,
	float *cy, float *w, float **rho, float **u, float **v, float
	omega, float ***feq, float ***feqb, float ***fin, float ***fout,
	float ***foutb);

void dfds_d(float p, float u0, float u0d, int m, int n, float *cx,
	float *cy, float *w, float **rho, float **u, float **v, float
	omega, float ***feq, float ***feqd, float ***fin, float ***fout,
	float ***foutd);

void dfdu_b(float p, float u0, int m, int n, float *cx, float *cy
	, float *w, float **rho, float **u, float **v, float omega,
	float ***feq, float ***feqb, float ***fin, float ***finb, float ***fout, float ***foutb);

void dfdu_d(float p, float u0, int m, int n, float *cx, float *cy
	, float *w, float **rho, float **u, float **v, float omega,
	float ***feq, float ***feqd, float ***fin, float ***find, float ***fout, float ***foutd);

//Objective function derivatives
void dObjdu_d(float p, float u0, int m, int n, float *cx, float *cy, float *w, float **rho, float **rhod, float **u, 
	float **v, float omega,float ***feq, float ***feqd, float ***fin, float ***find, float ***fout, float ***foutd, float *J, 
	float *Jd);

void dObjds_d(float p, float u0, float u0d, int m, int n, float *cx, float *cy, float *w, float **rho, float **rhod, 
	float **u, float **v,float omega, float ***feq, float ***feqd, float ***fin, float ***fout, float ***foutd, 
	float *J, float *Jd);
