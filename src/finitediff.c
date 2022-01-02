#include <stdio.h>
#include <stdlib.h>
#include "linearalg.h"
#include "finitediff.h"

// Computes finite-difference matrices for the first derivative
mtrx Diff1(int n, int o)
{
    int i;
    mtrx D;
    D.M = allocm(n, n);
    D.m = n;
    D.n = n;

    if (o == 2) // second order
    {
        D.M[0][0] = -1;
        D.M[0][1] = 1;
        for (i = 1; i < (n - 1); i++)
        {
            D.M[i][i - 1] = -0.5;
            D.M[i][i] = 0;
            D.M[i][i + 1] = 0.5;
        }
        D.M[n - 1][n - 1] = D.M[0][1];
        D.M[n - 1][n - 2] = D.M[0][0];
        return D;
    }
    else if (o == 4) // Fourth-order
    {
        D.M[0][0] = -1;
        D.M[0][1] = 1;
        D.M[1][0] = -0.5;
        D.M[1][1] = 0;
        D.M[1][2] = 0.5;
        for (i = 2; i < (n - 2); i++)
        {
            D.M[i][i - 2] = (double)1 / 12;
            D.M[i][i - 1] = (double)-2 / 3;
            D.M[i][i] = 0;
            D.M[i][i + 1] = (double)2 / 3;
            D.M[i][i + 2] = (double)-1 / 12;
        }
        D.M[n - 1][n - 1] = D.M[0][1];
        D.M[n - 1][n - 2] = D.M[0][0];
        D.M[n - 2][n - 1] = D.M[1][2];
        D.M[n - 2][n - 2] = D.M[1][1];
        D.M[n - 2][n - 3] = D.M[1][0];
        return D;
    }
    else if (o == 6) // Sixth-order
    {
        D.M[0][0] = -1;
        D.M[0][1] = 1;
        D.M[1][0] = -0.5;
        D.M[1][1] = 0;
        D.M[1][2] = 0.5;
        D.M[2][0] = (double)1 / 12;
        D.M[2][1] = (double)-2 / 3;
        D.M[2][2] = 0;
        D.M[2][3] = (double)2 / 3;
        D.M[2][4] = (double)-1 / 12;
        for (i = 3; i < (n - 3); i++)
        {
            D.M[i][i - 3] = (double)-1 / 60;
            D.M[i][i - 2] = (double)3 / 20;
            D.M[i][i - 1] = (double)-3 / 4;
            D.M[i][i] = 0;
            D.M[i][i + 1] = (double)3 / 4;
            D.M[i][i + 2] = (double)-3 / 20;
            D.M[i][i + 3] = (double)1 / 60;
        }
        D.M[n - 1][n - 1] = D.M[0][1];
        D.M[n - 1][n - 2] = D.M[0][0];
        D.M[n - 2][n - 1] = D.M[1][2];
        D.M[n - 2][n - 2] = D.M[1][1];
        D.M[n - 2][n - 3] = D.M[1][0];
        D.M[n - 3][n - 1] = D.M[2][4];
        D.M[n - 3][n - 2] = D.M[2][3];
        D.M[n - 3][n - 3] = D.M[2][2];
        D.M[n - 3][n - 4] = D.M[2][1];
        D.M[n - 3][n - 5] = D.M[2][0];
        return D;
    }
    else
    {
        printf("** Error: valid orders are 2, 4 or 6 **\n");
        exit(1);
    }
}

// Computes finite - difference matrices for the second derivative
mtrx Diff2(int n, int o)
{
    int i;
    mtrx D;
    D.M = allocm(n, n);
    D.m = n;
    D.n = n;
    zerosm(D);

    if (o == 2) // Second-order
    {
        D.M[0][0] = 2; // Forward scheme (second-order)
        D.M[0][1] = -5;
        D.M[0][2] = 4;
        D.M[0][3] = -1;
        for (i = 1; i < (n - 1); i++)
        {
            D.M[i][i - 1] = 1;
            D.M[i][i] = -2;
            D.M[i][i + 1] = 1;
        }
        D.M[n - 1][n - 1] = D.M[0][0];
        D.M[n - 1][n - 2] = D.M[0][1];
        D.M[n - 1][n - 3] = D.M[0][2];
        D.M[n - 1][n - 4] = D.M[0][3];
        return D;
    }
    else if (o == 4) // Fourth-order
    {
        D.M[0][0] = 2; // ForwarD.M scheme (second-order)
        D.M[0][1] = -5;
        D.M[0][2] = 4;
        D.M[0][3] = -1;
        D.M[1][0] = 1; // Central scheme (second-order)
        D.M[1][1] = -2;
        D.M[1][2] = 1;
        for (i = 2; i < (n - 2); i++)
        {
            D.M[i][i - 2] = (double)-1 / 12;
            D.M[i][i - 1] = (double)4 / 3;
            D.M[i][i] = (double)-5 / 2;
            D.M[i][i + 1] = (double)4 / 3;
            D.M[i][i + 2] = (double)-1 / 12;
        }
        D.M[n - 1][n - 1] = D.M[0][0];
        D.M[n - 1][n - 2] = D.M[0][1];
        D.M[n - 1][n - 3] = D.M[0][2];
        D.M[n - 1][n - 4] = D.M[0][3];
        D.M[n - 2][n - 1] = D.M[1][0];
        D.M[n - 2][n - 2] = D.M[1][1];
        D.M[n - 2][n - 3] = D.M[1][2];
        return D;
    }
    else if (o == 6) // Sixth-order
    {
        D.M[0][0] = 2; // Forward-scheme (second-order)
        D.M[0][1] = -5;
        D.M[0][2] = 4;
        D.M[0][3] = -1;
        D.M[1][0] = 1; // Central-scheme (second-order)
        D.M[1][1] = -2;
        D.M[1][2] = 1;
        D.M[2][0] = (double)-1 / 12; // Central-scheme (fourth-order)
        D.M[2][1] = (double)4 / 3;
        D.M[2][2] = (double)-5 / 2;
        D.M[2][3] = (double)4 / 3;
        D.M[2][4] = (double)-1 / 12;
        for (i = 3; i < (n - 3); i++)
        {
            D.M[i][i - 3] = (double)1 / 90;
            D.M[i][i - 2] = (double)-3 / 20;
            D.M[i][i - 1] = (double)3 / 2;
            D.M[i][i] = (double)-49 / 18;
            D.M[i][i + 1] = (double)3 / 2;
            D.M[i][i + 2] = (double)-3 / 20;
            D.M[i][i + 3] = (double)1 / 90;
        }
        D.M[n - 1][n - 1] = D.M[0][0];
        D.M[n - 1][n - 2] = D.M[0][1];
        D.M[n - 1][n - 3] = D.M[0][2];
        D.M[n - 1][n - 4] = D.M[0][3];
        D.M[n - 2][n - 1] = D.M[1][0];
        D.M[n - 2][n - 2] = D.M[1][1];
        D.M[n - 2][n - 3] = D.M[1][2];
        D.M[n - 3][n - 1] = D.M[2][0];
        D.M[n - 3][n - 2] = D.M[2][1];
        D.M[n - 3][n - 3] = D.M[2][2];
        D.M[n - 3][n - 4] = D.M[2][3];
        D.M[n - 3][n - 5] = D.M[2][4];
        return D;
    }
    else
    {
        printf("** Error: valid orders are 2, 4 or 6 **\n");
        exit(1);
    }
}