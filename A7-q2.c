//Program to solve an ODE using 4th Order Runge Kutta method

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mylib.h"

void main()
{
	double u[1], v[1], x, t;

	printf("Given ODE is :\n\ty'' + y' = 1 - x\n");
	printf("Initial boundary value:\n\ty(0) = 2\n\ty'(0) = 1\n\n");

	FILE* fptr1;
	fptr1 = fopen("C:/Users/sahoo/Desktop/P342/RK4a.txt", "w");

	FILE* fptr2;
	fptr2 = fopen("C:/Users/sahoo/Desktop/P342/RK4b.txt", "w");

	if (fptr1 == NULL || fptr2 == NULL)
	{
		// File not created hence exit
		printf("\nUnable to create file.\n");
		exit(EXIT_FAILURE);
	}

	//initiating for the range [0,5]
	x = 0;
	u[0] = 2;
	v[0] = 1;
	t = 5;
	double h = (t - x) / 1000;
	while (fabs(x) < fabs(t))
	{
		RK4(u, v, x, t, h, 1);
		fprintf(fptr1, "%lf\t%lf\n", x, u[0]);
		x += h;
	}
	fclose(fptr1);

	//initiating for the range [-5,0]
	x = 0;
	u[0] = 2;
	v[0] = 1;
	t = -5;
	h = (t - x) / 1000;
	while (fabs(x) < fabs(t))
	{
		RK4(u, v, x, t, h, 1);
		fprintf(fptr2, "%lf\t%lf\n", x, u[0]);
		x += h;
	}
	fclose(fptr2);

	printf("ODE solved using Runge Kutta method-4th order.\n");
	printf("Data storage locations:\n\tC:/Users/sahoo/Desktop/P342/RK4a.txt\n\tC:/Users/sahoo/Desktop/P342/RK4b.txt");
}






/*********************************OUTPUT************************************************
Given ODE is :
        y'' + y' = 1 - x
Initial boundary value:
        y(0) = 2
        y'(0) = 1

ODE solved using Runge Kutta method-4th order.
Data storage locations:
        C:/Users/sahoo/Desktop/P342/RK4a.txt
        C:/Users/sahoo/Desktop/P342/RK4b.txt
***************************************************************************************/
