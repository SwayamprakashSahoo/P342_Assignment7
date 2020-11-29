#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mylib.h"

void cases(double h);

void main()
{
	printf("Given ODE :\n\ty'' = y' + 1\n");
    printf("Given boundary values :\n\ty(0) = 1\n\ty(1) = 2(e-1)\n");
	
    cases(0.5);
    cases(0.2);
    cases(0.05);
    cases(0.02);
}

void cases(double h)
{
    printf("\n\nFor value of h = %lf\n", h);
    double x, u, t, y, v, e, m1, m2, m3, b, b1, b2, b3;
    int a;
    x = 0;
    u = 1;
    t = 1;
    y = 2 * (exp(1) - 1);
    printf("\nEnter the 1st guess M1: ");
    scanf("%lf", &m1);
    v = m1;
    b = y;
    b1 = shoot(u, v, x, t, h, 1);
    printf("\nB1 is %lf", b1);
    if (fabs(b1 - b) < 0.00005)
    {
        printf("\n  The value of x and respective z are:\n");
        e = shoot(u, v, x, t, h, 1);
        return(0);
    }
    else
    {
        printf("\nEnter the 2nd guess M2:");
        scanf("%lf", &m2);
        v = m2;
        b2 = shoot(u, v, x, t, h, 1);
        printf("\nB2 is %lf", b2);
    }
    if (fabs(b2 - b) < 0.00005)
    {
        printf("\n  The value of x and respective z are\n");
        e = shoot(u, v, x, t, h, 1);
        return(0);
    }
    else
    {
        printf("\n\nM2=%lf\tM1=%lf", m2, m1);
        m3 = m2 + (((m2 - m1) * (b - b2)) / (b2 - b1));
        if (b1 - b2 == 0)
            return;

        printf("\nExact value of M =%lf", m3);
        v = m3;
        b3 = shoot(u, v, x, t, h, 0);
    }
    if (fabs(b3 - b) < 0.00005)
    {
        printf("\nFor which solution is :\n");
        e = shoot(u, v, x, t, h, 1);
        return;
    }
    do
    {
        m1 = m2;
        m2 = m3;
        b1 = b2;
        b2 = b3;
        m3 = m2 + (((m2 - m1) * (b - b2)) / (b2 - b1));
        v = m3;
        b3 = shoot(u, v, x, t, h, 0);

    } while (fabs(b3 - b) < 0.00005);
    v = m3;
    e = shoot(u, v, x, t, h, 1);
}