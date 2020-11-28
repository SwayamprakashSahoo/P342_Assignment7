//Program to solve an ODE using shooting method

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





/*************************************OUTPUT*************************************
Given ODE :
        y'' = y' + 1
Given boundary values :
        y(0) = 1
        y(1) = 2(e-1)

Enter the 1st guess M1: 0

0.100000        2.005171
0.200000        3.026573
0.300000        4.076432
0.400000        5.168256
0.500000        6.316977
0.600000        7.539095
0.700000        8.852846
0.800000        10.278386
0.900000        11.837987
1.000000        13.556267
1.100000        15.460431
B1 is 15.460431
Enter the 2nd guess M2:-2

0.100000        1.794829
0.200000        2.373427
0.300000        2.723568
0.400000        2.831744
0.500000        2.683023
0.600000        2.260905
0.700000        1.547154
0.800000        0.521614
0.900000        -0.837987
1.000000        -2.556267
1.100000        -4.660431
B2 is -4.660431
M2=-2.000000    M1=0.000000
Exact value of M =-1.195164
For which solution is :

0.100000        1.879474
0.200000        2.636265
0.300000        3.267985
0.400000        3.771998
0.500000        4.145391
0.600000        4.384943
0.700000        4.487095
0.800000        4.447913
0.900000        4.263051
1.000000        3.927705
1.100000        3.436564
*****************************************************************************/
