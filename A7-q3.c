#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mylib.h"

void main()
{
	printf("Given ODE :\n\ty'' = y' + 1\n");
	double x, u, t, y, v, e, h, m1, m2, m3, b, b1, b2, b3;
	int a;
	x = 0;
	u = 1;
	t = 1;
	y = 2 * (exp(1) - 1);
    h = (t - x) / 10;
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
        printf("\nM2=%lf\tM1=%lf", m2, m1);
        m3 = m2 + (((m2 - m1) * (b - b2)) / (b2 - b1));
        if (b1 - b2 == 0)
            exit(0);

        printf("\nExact value of M =%lf", m3);
        v = m3;
        b3 = shoot(u, v, x, t, h, 0);
    }
    if (fabs(b3 - b) < 0.000005)
    {
        printf("\nFor which solution is :\n");
        e = shoot(u, v, x, t, h, 1);
        exit(0);
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

    } while (fabs(b3 - b) < 0.0005);
    v = m3;
    e = shoot(u, v, x, t, h, 1);
}