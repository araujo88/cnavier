#include <stdio.h>
#include <stdlib.h>
#include "fluiddyn.h"
#include "linearalg.h"

void euler(mtrx w, mtrx dwdx, mtrx dwdy, mtrx d2wdx2, mtrx d2wdy2, mtrx u, mtrx v, double Re, double dt)
{
    int i, j;

    for (i = 0; i < w.m; i++)
    {
        for (j = 0; j < w.n; j++)
        {
            w.M[i][j] = (-u.M[i][j] * dwdx.M[i][j] - v.M[i][j] * dwdy.M[i][j] + (1. / Re) * (d2wdx2.M[i][j] + d2wdy2.M[i][j])) * dt + w.M[i][j];
        }
    }
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

mtrx vorticity(mtrx dudy, mtrx dvdx)
{
    int i, j;
    mtrx temp;
    temp = initm(dudy.m, dudy.n);

    for (i = 0; i < temp.m; i++)
    {
        for (j = 0; j < temp.n; j++)
        {
            temp.M[i][j] = dvdx.M[i][j] - dudy.M[i][j];
        }
    }
    return temp;
}