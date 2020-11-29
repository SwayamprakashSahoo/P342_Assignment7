#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mylib.h"

void cases(double h);

void main()
{
	printf("Given ODE is :\n\ty'' + y' = 1 - x\n");
	printf("Initial boundary value:\n\ty(0) = 2\n\ty'(0) = 1\n\n");

	cases(0.5);
	cases(0.2);
	cases(0.05);
	cases(0.02);
	
	printf("\nFile path to Data storage locations:\n\tC:/Users/sahoo/Desktop/P342/");
}

void cases(double h)
{
	double u[1], v[1], x, t;
	FILE* fptr1;
	FILE* fptr2;
	if (h == 0.5)
	{
		fptr1 = fopen("C:/Users/sahoo/Desktop/P342/RK41a.txt", "w");
		fptr2 = fopen("C:/Users/sahoo/Desktop/P342/RK41b.txt", "w");
	}
	else if(h==0.2)
	{
		fptr1 = fopen("C:/Users/sahoo/Desktop/P342/RK42a.txt", "w");
		fptr2 = fopen("C:/Users/sahoo/Desktop/P342/RK42b.txt", "w");
	}
	else if(h==0.05)
	{
		fptr1 = fopen("C:/Users/sahoo/Desktop/P342/RK43a.txt", "w");
		fptr2 = fopen("C:/Users/sahoo/Desktop/P342/RK43b.txt", "w");
	}
	else
	{
		fptr1 = fopen("C:/Users/sahoo/Desktop/P342/RK44a.txt", "w");
		fptr2 = fopen("C:/Users/sahoo/Desktop/P342/RK44b.txt", "w");
	}

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
	double p = -h;
	while (fabs(x) < fabs(t))
	{
		RK4(u, v, x, t, p, 1);
		fprintf(fptr2, "%lf\t%lf\n", x, u[0]);
		x += p;
	}
	fclose(fptr2);
	printf("ODE solved using Runge Kutta method-4th order with h = %lf.\n", h);
}