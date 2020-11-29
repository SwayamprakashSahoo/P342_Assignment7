#pragma once
void LUdecompose(int var, float mat[4][4]);
void LUsolution(int var, float mat[4][4], float constant[4], float sol[4]);
void LUpivot(float mat[4][4], float permute[4][4], int var);
void LUinverse(float mat[4][4], float permute[4][4], float inv[4][4], int var);
float LUdeterminant(float mat[4][4], int n);
float determinant(float mat[4][4], int n);

void pivot(float mat[][20], int var);
void GaussJordan(float mat[][20], int var);
void Matrix_mult(float mat[][20], float inv[][20], int var);

double func(double n, int x);
double bracket(double arr[2], int x);
double bisection(double a, double b, int x);
double regulaFalsi(double a, double b, int x);
void newtonRaphson(double a, int x);

double* derivative(double* array, int n);
double polVal(double* coeff, double x, int n);
double laguerre(double* coeff, double x, int n);
double* synDiv(double* array, double x, int n);

double midpoint(double ll, double ul, int count, int x);
double trapezoidal(double ll, double ul, int count, int x);
double simpsons(double ll, double ul, int count, int x);
double monteCarlo(int n);

void euler(double x, double y, double h, double t, int n);
void RK4(double u[1], double v[1], double x, double t, double h);
double shoot(double u, double v, double x, double t, double h, int a);