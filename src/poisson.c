#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "linearalg.h"
#include "poisson.h"

double error(mtrx u1, mtrx u2)
{
    double e = 0;
    int i, j;

    for (i = 0; i < u1.m; i++)
    {
        for (j = 0; j < u1.n; j++)
        {
            e += sqrt(pow(u2.M[i][j] - u1.M[i][j], 2));
        }
    }
    return e;
}

mtrx poisson(mtrx f, double dx, double dy, int itmax, double tol)
{
    int i, j, k, nx, ny;
    double e;
    //double u1max, u2max;
    mtrx u, u0;
    u = initm(f.m, f.n);
    u0 = initm(f.m, f.n);
    nx = f.m;
    ny = f.n;

    for (k = 0; k < itmax; k++)
    {
        mtrxcpy(u0, u);
        for (i = 1; i < nx - 1; i++)
        {
            for (j = 1; j < ny - 1; j++)
            {
                u.M[i][j] = (dy * dy * (u.M[i + 1][j] + u.M[i - 1][j]) + dx * dx * (u.M[i][j + 1] + u.M[i][j - 1]) - dx * dx * dy * dy * f.M[i][j]) / (2 * (dx * dx + dy * dy));
            }
        }
        e = error(u, u0);
        if (e < tol)
        {
            printf("Poisson equation solved with %d iterations - root-sum-of-squares error: %E\n", k, e);
            u0.M = freem(u0);
            return u;
        }
    }
    printf("Error: maximum number of iterations achieved for Poisson equation.\n");

    u.M = freem(u);
    u0.M = freem(u0);
    exit(1);
}

mtrx poisson_SOR(mtrx f, double dx, double dy, int itmax, double tol, double beta)
{
    int i, j, k, nx, ny;
    double e;
    //double u1max, u2max;
    mtrx u, u0;
    u = initm(f.m, f.n);
    u0 = initm(f.m, f.n);
    nx = f.m;
    ny = f.n;

    for (k = 0; k < itmax; k++)
    {
        mtrxcpy(u0, u);
        for (i = 1; i < nx - 1; i++)
        {
            for (j = 1; j < ny - 1; j++)
            {
                u.M[i][j] = beta * (dy * dy * (u.M[i + 1][j] + u.M[i - 1][j]) + dx * dx * (u.M[i][j + 1] + u.M[i][j - 1]) - dx * dx * dy * dy * f.M[i][j]) / (2 * (dx * dx + dy * dy)) + (1 - beta) * u0.M[i][j];
            }
        }
        e = error(u, u0);
        if (e < tol)
        {
            printf("Poisson equation solved with %d iterations - root-sum-of-squares error: %E\n", k, e);
            u0.M = freem(u0);
            return u;
        }
    }
    printf("Error: maximum number of iterations achieved for Poisson equation.\n");

    u.M = freem(u);
    u0.M = freem(u0);
    exit(1);
}