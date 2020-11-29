//Program to solve a first order ode using explicit euler's method

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mylib.h"

void main()
{
	//initial values for function 1
	double x = 2;
	double y = exp(1);
	double t1, t2;

	printf("\nGiven ODE:\n\t(dy/dx) = (y ln y)/x\n");
	printf("Initial value : y(2) = e\n");
	printf("Please enter the final value of x :\n ");
	scanf_s("%lf", &t1);

	printf("\nFor h = 0.5\n");
	euler(x, y, t1, 0.5, 1);
	printf("\nFor h = 0.2\n");
	euler(x, y, t1, 0.2, 1);
	printf("\nFor h = 0.05\n");
	euler(x, y, t1, 0.05, 1);
	printf("\nFor h = 0.02\n");
	euler(x, y, t1, 0.02, 1);


	//initial values for function 2
	x = 3;
	y = 1;

	printf("\nGiven ODE:\n\t(dy/dx) = 6 - (2y)/x\n");
	printf("Initial value : y(3) = 1\n");
	printf("Please enter the final value of x :\n ");
	scanf_s("%lf", &t2);

	printf("\nFor h = 0.5\n");
	euler(x, y, t2, 0.5, 2);
	printf("\nFor h = 0.2\n");
	euler(x, y, t2, 0.2, 2);
	printf("\nFor h = 0.05\n");
	euler(x, y, t2, 0.05, 2);
	printf("\nFor h = 0.02\n");
	euler(x, y, t2, 0.02, 2);
}