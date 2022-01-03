#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "linearalg.h"
#include "poisson.h"

mtrx poisson(mtrx f, double dx, double dy, int itmax, double tol)
{
    int i, j, k, nx, ny;
    double u_temp, error, w, rho, a, b, c, d, u1max, u2max;
    mtrx u, u0;
    u = initm(f.m, f.n);
    u0 = initm(f.m, f.n);
    nx = f.m;
    ny = f.n;

    rho = 0.5 * (cos(PI / nx) + cos(PI / ny));
    w = 2 / (1 + sqrt(1 - rho * rho)) - 1;
    a = dy * dy;
    b = dx * dx;
    c = 2 * (a + b);
    d = a * b;

    for (k = 0; k < itmax; k++)
    {
        u1max = maxel(u);
        for (i = 1; i < nx - 1; i++)
        {
            for (j = 1; j < ny - 1; j++)
            {
                u_temp = (a * (u.M[i - 1][j] + u.M[i + 1][j]) + b * (u.M[i][j - 1] + u.M[i][j + 1])) / c;
                u.M[i][j] = (1 - w) * u0.M[i][j] + w * (u_temp - d * f.M[i][j] / c);
            }
        }
        u2max = maxel(u);
        error = abs((u2max - u1max) / u1max);
        printf("Abs error: %lf\n", error);
        if (error < tol)
        {
            printf("Poisson equation solved with %d iterations.\n", k);
            return u;
        }
        else
        {
            u0 = u;
        }
    }
    printf("Error: maximum number of iterations achieved for Poisson equation.\n");
    exit(1);
}