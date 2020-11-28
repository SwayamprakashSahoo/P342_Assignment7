#include "mylib.h"
#include <stdlib.h>
#include <math.h>
#include <stdio.h>

//function to perform LU decomposition on a matrix and store it in the same matrix
void LUdecompose(int var, float mat[4][4])
{
	for (int k = 0; k < var; k++)
	{
		for (int j = k; j < var; j++)
		{
			float sum = 0;
			for (int p = 0; p < k; p++)
			{
				sum += mat[k][p] * mat[p][j];
			}
			mat[k][j] = (mat[k][j] - sum);
		}
		for (int i = k + 1; i < var; i++)
		{
			float sum = 0;
			for (int p = 0; p < k; p++)
			{
				sum += mat[i][p] * mat[p][k];
			}
			mat[i][k] = (mat[i][k] - sum) / mat[k][k];
		}
	}
}

//function to find the solution from the LU decomposition
void LUsolution(int var, float mat[4][4], float constant[4], float sol[4])
{
	float y[4] = { 0 };
	for (int i = 0; i < var; i++)
	{
		float sum = 0;
		for (int k = 0; k < i; k++)
		{
			sum += mat[i][k] * y[k];
		}
		y[i] = (constant[i] - sum);
	}
	for (int i = var - 1; i >= 0; i--)
	{
		float sum = 0;
		for (int k = i + 1; k < var; k++)
		{
			sum += mat[i][k] * sol[k];
		}
		sol[i] = (y[i] - sum) / mat[i][i];
	}
}

//Function to perform row-wise pivoting of the matrix
void LUpivot(float mat[4][4],float permute[4][4], int var)
{
	for (int i = 0; i < var - 1; i++)
	{
		if (mat[i][i] == 0)
		{
			for (int j = i + 1; j < var; j++)
			{
				if (abs(mat[j][i]) > mat[i][i])
				{
					for (int k = 0; k < var; k++)
					{
						float temp1 = mat[i][k];
						mat[i][k] = mat[j][k];
						mat[j][k] = temp1;
						float temp2 = permute[i][k];
						permute[i][k] = permute[j][k];
						permute[j][k] = temp2;
						continue;
					}
				}
			}
		}
	}
}

//function to find determinant of a matrix
float determinant(float mat[4][4], int n)
{
	float Minor[4][4];
	int c1, c2;
	float det;
	float c[4];
	float O = 1;
	det = 0;
	if (n == 2)
	{
		det = mat[0][0] * mat[1][1] - mat[0][1] * mat[1][0];
		return det;
	}
	else
	{
		for (int i = 0; i < n; i++)
		{
			c1 = 0, c2 = 0;
			for (int j = 0; j < n; j++)
			{
				for (int k = 0; k < n; k++)
				{
					if (j != 0 && k != i)
					{
						Minor[c1][c2] = mat[j][k];
						c2++;
						if (c2 > n - 2)
						{
							c1++;
							c2 = 0;
						}
					}
				}
			}
			det = det + O * (mat[0][i] * determinant(Minor, n - 1));
			O = -1 * O;
		}
	}
	return det;
}

//function to find determinant using LU-decomposition
float LUdeterminant(float mat[4][4], int n)
{
	float det = 1;
	for (int i = 0; i < n; i++)
	{
		det = det * mat[i][i];
	}
	return det;
}

//function to find inverse using LU decomposition
void LUinverse(float mat[4][4], float permute[4][4], float inv[4][4], int var)
{
	for (int i = 0; i < var; i++)
	{
		float temp[4] = { 0 };
		for (int j = 0; j < var; j++)
			temp[j] = permute[j][i];
		float y[4];
		for (int j = 0; j < var; j++)
		{
			float x = 0.0;
			for (int k = 0; k <= j - 1; k++)
			{
				x = x + mat[j][k] * y[k];
			}
			y[j] = (temp[j] - x);
		}

		for (int j = var - 1; j >= 0; j--)
		{
			float x = 0.0;
			for (int k = j + 1; k < var; k++)
			{
				x = x + mat[j][k] * inv[k][i];
			}
			inv[j][i] = (y[j] - x) / mat[j][j];
		}
	}
}

//function for pivoting for Gauss-Jordan elimination
void pivot(float mat[][20], int var)
{
	for (int i = 0; i < var - 1; i++)
	{
		if (mat[i][i] == 0)
		{
			for (int j = i + 1; j < var; j++)
			{
				if (abs(mat[j][i]) > mat[i][i])
				{
					for (int k = i; k < var + 1; k++)
					{
						float temp = mat[i][k];
						mat[i][k] = mat[j][k];
						mat[j][k] = temp;
						continue;
					}
				}
			}
		}
	}
}

//function for solve for solution using Gauss-Jordan elimination and print it
void GaussJordan(float mat[][20], int var)
{
	for (int i = 0; i < var; i++)
	{
		pivot(mat, var);
	}

	float c, x[20];
	for (int j = 0; j < var; j++)
	{
		for (int i = 0; i < var; i++)
		{
			if (i != j)
			{
				c = mat[i][j] / mat[j][j];
				for (int k = 0; k < var + 1; k++)
				{
					mat[i][k] = mat[i][k] - c * mat[j][k];
				}
			}
		}
	}
	printf("\nThe solution is:\n");
	for (int i = 0; i < var; i++)
	{
		x[i] = mat[i][var] / mat[i][i];
		printf("\n x%d=%f\n", i + 1, x[i]);
	}
}

//function to multiply two matrices
void Matrix_mult(float mat[][20], float inv[][20], int var)
{
	float product[20][20];
	printf("\nThe product of matrix A and its inverse is:\n");
	for (int i = 0; i < var; i++)
	{
		for (int j = 0; j < var; j++)
		{
			product[i][j] = 0;
			for (int k = 0; k < var; k++)
				product[i][j] = product[i][j] + mat[i][k] * inv[k][j];
			printf("%d ", (int)(product[i][j]));
		}
		printf("\n");
	}
}



//essential functions for assignment 5


//Input function for solving 
double func(double n, int x)
{
	if (x == 1)
		return (log(n) - sin(n));
	else
		return (-n - cos(n));
}
//Function derivative
double g(double n, int x)
{
	if (x == 1)
		return (1 / n - cos(n));
	else
		return (-1 + sin(n));
}

//function to bracket the root of the function in a given interval
double bracket(double arr[2], int x)
{
	double product = func(arr[0], x) * func(arr[1], x);
	if (product < 0)
	{
		printf("Root already exists in the guessed interval");
	}
	else
	{
		while (product > 0)
		{
			if (fabs(func(arr[0], x)) < fabs(func(arr[1], x)))
			{
				arr[0] = arr[0] + 1.5 * (arr[1] - arr[0]);
				product = func(arr[0], x) * func(arr[1], x);
			}
			else
			{
				arr[1] = arr[1] + 1.5 * (arr[1] - arr[0]);
				product = func(arr[0], x) * func(arr[1], x);
			}
		}
		printf("The root exists in the interval: [%lf, %lf]", arr[0], arr[1]);
	}
}

//function to use bisection method to solve an equation
double bisection(double a, double b, int x)
{
	FILE* fptr;
	if(x == 1)
		fptr = fopen("C:/Users/sahoo/Desktop/P342/bisection(a).txt", "w");
	else
		fptr = fopen("C:/Users/sahoo/Desktop/P342/bisection(b).txt", "w");
	if (fptr == NULL)
	{
		// File not created hence exit
		printf("\nUnable to create file.\n");
		exit(EXIT_FAILURE);
	}

	int count = 1;
	double eps = 0.000001;
	double c = a;
	double temp, err;
	printf("\nIt.No.\t\ta\t\t\tb\t\t\tc\t\t\tf(c)\t\t\tError\n");
	do
	{ 
		c = (a + b) / 2;
		temp = c;
		double f = func(c, x);
		printf(" %d\t\t%lf\t\t%lf\t\t%lf\t\t%lf", count, a, b, temp, f);

		if (func(c, x) == 0.0)
			continue;
		else if (func(c, x) * func(a, x) < 0)
			b = c;
		else
			a = c;

		c = (a + b) / 2;
		err = fabs(temp - c);
		fprintf(fptr, "%d\t%lf\n", count, err);
		printf("\t\t%lf\n",err);
		count++;
	} while (fabs(b - a) >= eps && count <= 200 && fabs(func(c, x)) > eps);
	
	printf("\nThe value of root is : %f", c);
}


//function to use Regula Falsi method to solve an equation
double regulaFalsi(double a, double b, int x)
{
	FILE* fptr;
	if(x == 1)
		fptr = fopen("C:/Users/sahoo/Desktop/P342/regulaPhalsi(a).txt", "w");
	else 
		fptr = fopen("C:/Users/sahoo/Desktop/P342/regulaPhalsi(b).txt", "w");
	
	if (fptr == NULL)
	{
		// File not created hence exit
		printf("\nUnable to create file.\n");
		exit(EXIT_FAILURE);
	}

	int count = 1;
	double c, eps = 0.000001;
	double f, temp, err;
	printf("\nIt.No.\t\ta\t\t\tb\t\t\tc\t\t\tf(c)\t\t\tError\n");
	do
	{
		c = (a*func(b, x)-b*func(a, x))/(func(b, x)-func(a, x));
		temp = c;
		f = func(c, x);
		printf(" %d\t\t%lf\t\t%lf\t\t%lf\t\t%lf", count, a, b, temp, f);

		if (func(c, x) == 0.0)
			continue;
		else if (func(c, x) * func(a, x) < 0)
			b = c;
		else
			a = c;
		
		c = (a * func(b, x) - b * func(a, x)) / (func(b, x) - func(a, x));
		err = fabs(temp - c);
		fprintf(fptr, "%d\t%lf\n", count, err);
		printf("\t\t%lf\n", err);
		count++;

	} while (fabs(f) > eps && count <= 200 && f != 0);

	printf("\nThe root, using Regula Falsi method, is: %lf", c);
}

//Function showcasing the use of Newton-Raphson method 
void newtonRaphson(double a, int x)
{
	FILE* fptr;
	if(x == 1)
		fptr = fopen("C:/Users/sahoo/Desktop/P342/newtonRaphson(a).txt", "w");
	else
		fptr = fopen("C:/Users/sahoo/Desktop/P342/newtonRaphson(b).txt", "w");
	if (fptr == NULL)
	{
		// File not created hence exit
		printf("\nUnable to create file.\n");
		exit(EXIT_FAILURE);
	} 

	double g0, f0, f1, c;
	int count = 1;
	double eps = 0.000001;
	double temp, err;
	printf("\nIter.No.\ta\t\tf(a)\t\tc\t\tError\n");
	do
	{
		g0 = g(a, x);
		f0 = func(a, x);
		if (g0 == 0.0)
		{
			printf("Mathematical Error.");
			exit(0);
		}
		c = a - f0 / g0;
		temp = c;
		printf("%d\t\t%lf\t%lf\t%lf", count, a, f0, c);
		a = c;
		g0 = g(a, x);
		f0 = func(a, x);
		c = a - f0 / g0;
		err = fabs(temp - c);
		fprintf(fptr, "%d\t%lf\n", count, err);
		printf("\t%lf\n", err);
		count++;

		if (count > 199)
		{
			printf("The series is not convergent.");
			exit(0);
		}
		f1 = func(c, x);
	} while (fabs(f1) > eps);

	printf("\nRoot is: %f", c);
}

//finding derivative of a polynomial
double* derivative(double* array, int n)
{
	int m = n - 1;
	double* der = malloc(n * sizeof(double));
	for (int i = 0; i <= m; i++)
	{
		der[i] = array[i]*(n-i);
	}
	return der;
}

//value of a polynomial function
double polVal(double* coeff, double x, int n)
{
	double val = 0;
	for (int i = 0; i <= n; i++)
	{
		double prod = 1;
		for (int j = i; j < n; j++)
			prod = prod * x;
		val = val + (coeff[i] * prod);
	}
	return val;
}

//laguerre
double laguerre(double* coeff, double x, int n)
{
	int count = 0;
	double f, g, h, a;
	double eps = 0.000001;
	double* der1 = derivative(&coeff[0], n);
	double* der2 = derivative(&der1[0], n-1);

	do
	{
		f = polVal(&coeff[0], x, n);
		g = polVal(&der1[0], x, n - 1);
		h = polVal(&der2[0], x, n - 2);
		double G = g / f;
		double H = (G * G) - (h / f);
		double F = sqrt((n - 1) * (n * H - G * G));
		if (fabs(g + f) > fabs(g - f))
			a = n / (G + F);
		else
			a = n / (G - F);
		x = x - a;
			
		count++;
	} while (f > eps && count <= 200 && a > eps);
	return x;
}

//synthetic division
double* synDiv(double* array, double x, int n)
{
	double* rem = malloc(n * sizeof(double));
	rem[0] = array[0];
	for (int i = 1; i <= n; i++) 
	{
		rem[i] = (rem[i - 1] * x) + array[i];
	}
	return &rem[0];
}


//essential functions for assignment 6

//input functions for integration
double f(double n, int x)
{
	if (x == 1)
		return n / (n + 1);
	else if (x == 2)
		return exp(-(n * n));
	else
		return (4 / (1 + (n * n)));
}

//Integration using Midpoint rule
double midpoint(double ll, double ul, int count, int x)
{
	double step = (ul - ll) / count;
	double integral = 0;
	double lowerbound, upperbound, midpt;
	lowerbound = ll;
	
	for (int i = 1; i <= count; i++)
	{
		upperbound = lowerbound + step;
		midpt = (lowerbound + upperbound)/ 2;
		integral += step * f(midpt, x);
		
		lowerbound = upperbound;
	}
	return integral;
}

//Integration using Trapezoidal-rule
double trapezoidal(double ll, double ul, int count, int x)
{
	double step = (ul - ll) / count;
	double integral = 0;
	double lowerbound, upperbound;
	lowerbound = ll;
	
	for (int i = 1; i <= count; i++)
	{
		upperbound = lowerbound + step;
		double a = f(lowerbound, x);
		double b = f(upperbound, x);
		integral += step * ((a + b)/2);
		lowerbound = upperbound;
	}
	return integral;
}

//Integration using Simpson's rule
double simpsons(double ll, double ul, int count, int x)
{
	double step = (ul - ll) / count;
	double sum = 0;
	double integral, lowerbound;

	for (int i = 1; i < count; i++) 
	{
		lowerbound = ll + i * step;
		if (i % 2 == 0) 
		{
			sum = sum + 2 * f(lowerbound, x);
		}
		else 
		{
			sum = sum + 4 * f(lowerbound, x);
		}
	}
	integral = (step / 3) * (f(ll, x) + f(ul, x) + sum);
	return integral;
}




//Finding the value of pi using Monte-Carlo method
double monteCarlo(int n)
{
	double sum = 0;
	for (int i = 0; i < n; i++)
	{
		double random = rand() / ((double)RAND_MAX);
		double res = f(random, 3);
		sum += res;
	}
	double avg = sum / n;
	return avg;
}


//essential functions for assignment 7




double funcEu(double x, double y, int n)
{
	if (n == 1)
		return ((y * log(y)) / x);
	if (n == 2)
		return (6 - ((2 * y) / x));
}

double funcRK4(double u, double v, double x, int n, int m)
{
	if (m == 1)
	{
		if (n == 1)
			return v;
		if (n == 2)
			return (1 - x - v);
	}
	else
	{
		if (n == 1)
			return v;
		else
			return (1 + v);
	}
}

//function using explicit Euler's method
void euler(double x, double y,double t, int n)
{
	FILE* fptr;
	if (n == 1)
		fptr = fopen("C:/Users/sahoo/Desktop/P342/euler1.txt", "w");
	else
		fptr = fopen("C:/Users/sahoo/Desktop/P342/euler2.txt", "w");
	if (fptr == NULL)
	{
		// File not created hence exit
		printf("\nUnable to create file.\n");
		exit(EXIT_FAILURE);
	}

	double k;
	double h = (t-x)/1000;
	printf("\n  x\t  y\n");
	while (x <= t)
	{
		k = h * funcEu(x, y, n);
		fprintf(fptr, "%lf\t%lf\n", x, y);
		printf("%0.3lf\t%0.3lf\n", x, y);
		y = y + k;
		x = x + h;
	}
}

//Function to perform Renge-Kutta 4th order
void RK4(double u[1], double v[1], double x, double t, double h, int p)
{
	double k11, k12, k13, k14;
	double k21, k22, k23, k24;
	k11 = h * funcRK4(u[0], v[0], x, 1, p);
	k21 = h * funcRK4(u[0], v[0], x, 2, p);
	k12 = h * funcRK4(u[0] + 0.5 * k11, v[0] + 0.5 * k21, x + 0.5 * h, 1, p);
	k22 = h * funcRK4(u[0] + 0.5 * k11, v[0] + 0.5 * k21, x + 0.5 * h, 2, p);
	k13 = h * funcRK4(u[0] + 0.5 * k12, v[0] + 0.5 * k22, x + 0.5 * h, 1, p);
	k23 = h * funcRK4(u[0] + 0.5 * k12, v[0] + 0.5 * k22, x + 0.5 * h, 2, p);
	k14 = h * funcRK4(u[0] + k13, v[0] + k23, x + h, 1, p);
	k24 = h * funcRK4(u[0] + k13, v[0] + k23, x + h, 2, p);
	u[0] += (k11 + 2 * k12 + 2 * k13 + k14) / 6;
	v[0] += (k21 + 2 * k22 + 2 * k23 + k24) / 6;
}

double shoot(double u, double v, double x, double t, double h, int a)
{
	double y[1], guess[1];
	y[0] = u;
	guess[0] = v;

	while (fabs(x) < fabs(t))
	{
		RK4(y, guess, x, t, h, 2);
		u += y[0];
		v += guess[0];
		x += h;
		if (a == 1)
		{
			printf("\n%lf\t%lf", x, u);
		}
	}
	return u;
}

