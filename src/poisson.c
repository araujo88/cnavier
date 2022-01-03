#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "linearalg.h"
#include "poisson.h"

mtrx poisson(mtrx f, double dx, double dy, int itmax, double tol)
{
    int i, j, k, nx, ny;
    double error, u1max, u2max;
    mtrx u;
    u = initm(f.m, f.n);
    nx = f.m;
    ny = f.n;

    for (k = 0; k < itmax; k++)
    {
        u1max = maxel(u);
        // printf("u1max: %lf\n", u1max);
        for (i = 1; i < nx - 1; i++)
        {
            for (j = 1; j < ny - 1; j++)
            {
                u.M[i][j] = (dy * dy * (u.M[i + 1][j] + u.M[i - 1][j]) + dx * dx * (u.M[i][j + 1] + u.M[i][j - 1]) - dx * dx * dy * dy * f.M[i][j]) / (2 * (dx * dx + dy * dy));
            }
        }
        u2max = maxel(u);
        //printf("u2max: %lf\n", u2max);
        error = fabs((u2max - u1max));
        //printf("Abs error: %lf\n", error);
        if (error < tol)
        {
            printf("Poisson equation solved with %d iterations - abs error: %E\n", k, error);
            return u;
        }
    }
    printf("Error: maximum number of iterations achieved for Poisson equation.\n");
    exit(1);
}