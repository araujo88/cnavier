#include <stdio.h>
#include <stdlib.h>
#include "fluiddyn.h"
#include "linearalg.h"

mtrx euler(mtrx w, mtrx dwdx, mtrx dwdy, mtrx d2wdx2, mtrx d2wdy2, mtrx u, mtrx v, int Re, double dt)
{
    int i, j;
    mtrx temp;
    temp = initm(w.m, w.n);

    for (i = 0; i < w.m; i++)
    {
        for (j = 0; j < w.n; j++)
        {
            temp.M[i][j] = (-u.M[i][j] * dwdx.M[i][j] - v.M[i][j] * dwdy.M[i][j] + (double)(1 / Re) * (d2wdx2.M[i][j] + d2wdy2.M[i][j])) * dt + w.M[i][j];
        }
    }

    return temp;
}

mtrx continuity(mtrx dudx, mtrx dvdy)
{
    int i, j;
    mtrx temp;
    temp = initm(dudx.m, dudx.n);

    for (i = 0; i < temp.m; i++)
    {
        for (j = 0; j < temp.n; j++)
        {
            temp.M[i][j] = dudx.M[i][j] + dvdy.M[i][j];
        }
    }
    return temp;
}