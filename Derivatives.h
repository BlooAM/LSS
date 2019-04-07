#pragma once
#include "LBM.h"
#include "adBuffer.h"
void dfds_b(double p, double u0, double *u0b, int m, int n, double *cx,
	double *cy, double *w, double **rho, double **u, double **v, double
	omega, double ***feq, double ***feqb, double ***fin, double ***fout,
	double ***foutb);

void dfds_d(double p, double u0, double u0d, int m, int n, double *cx,
	double *cy, double *w, double **rho, double **u, double **v, double
	omega, double ***feq, double ***feqd, double ***fin, double ***fout,
	double ***foutd);

void dfdu_b(double p, double u0, int m, int n, double *cx, double *cy
	, double *w, double **rho, double **u, double **v, double omega,
	double ***feq, double ***feqb, double ***fin, double ***finb, double ***fout, double ***foutb);

void dfdu_d(double p, double u0, int m, int n, double *cx, double *cy
	, double *w, double **rho, double **u, double **v, double omega,
	double ***feq, double ***feqd, double ***fin, double ***find, double ***fout, double ***foutd);
