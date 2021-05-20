#pragma  once

void Jacobi(double *A, int size, double *b, double *x_initial);

void Gauss_Seidel(double *A, int size, double *b, double *x_initial);

void SOR(double *A, int size, double *b, double relaxation_factor, double *x_initial);