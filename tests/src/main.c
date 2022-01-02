#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../../include/linearalg.h"

void ReshapeTest()
{
    printf("*** Reshape test ***\n\n");

    mtrx A, B;
    A.M = readm("reshape/matrix.txt", &A.m, &A.n); // read matrix
    printf("Input matrix (%d x %d):\n", A.m, A.n);
    printm(A);       // print matrix
    B.m = A.m * A.n; // number of rows of reshaped matrix
    B.n = 1;         // number of columns of reshaped matrix
    printf("Reshaped matrix (%d x %d):\n", B.m, B.n);
    B = reshape(A, B.m, B.n); // reshape matrix
    printm(B);                // print matrix
    A.M = freem(A);           // Free memory matrix A
    B.M = freem(B);           // Free memory matrix B

    printf("*** End of reshape test ***\n\n");
}

void MatrixMulTest()
{
    printf("*** Matrix multiplication test ***\n\n");

    mtrx A, B, C;
    A.M = readm("matrix_mul/matrix2.txt", &A.m, &A.n); // read matrix 1
    B.M = readm("matrix_mul/matrix1.txt", &B.m, &B.n); // read matrix 1
    printf("Input matrix A (%d x %d):\n", A.m, A.n);
    printm(A); // print matrix 1
    printf("Input matrix B (%d x %d):\n", B.m, B.n);
    printm(B); // print matrix 1
    printf("Matrix multiplication (A x B):\n");
    C = mtrxmul(A, B, &C.m, &C.n); // matrix multiplifaction
    printm(C);                     // print matrix multiplication result
    A.M = freem(A);                // Free memory matrix A
    B.M = freem(B);                // Free memory matrix B
    C.M = freem(C);                // Free memory matrix C

    printf("*** End of matrix multiplication test ***\n\n");
}

void GaussianTest()
{
    printf("*** Gaussian ellimination test ***\n\n");

    mtrx A;
    vec b, x;

    A.M = readm("gaussian/matrix1.txt", &A.m, &A.n); // read matrix
    b.v = readv("gaussian/vector1.txt", &b.n);       // read vector
    printf("Input matrix:\n");
    printm(A); // print matrix
    printf("Input vector:\n");
    printv(b); // print vector
    x = gaussian(A, b);
    printf("Solution:\n");
    printv(x);

    if ((floor(x.v[0]) == 2) && (floor(x.v[1]) == 3) && (floor(x.v[2]) == -1))
    {
        printf("Test passed!\n");
    }
    else
    {
        printf("Test failed!\n");
    }

    A.M = freem(A); // Free memory matrix
    b.v = freev(b); // Free memory vector
    x.v = freev(x); // Free memory solution

    printf("*** End of gaussian ellimination test ***\n\n");
}

void KroneckerTest()
{
    printf("*** Kronecker product test ***\n\n");

    mtrx M1, M2, K;
    M1.M = readm("kronecker/matrix1.txt", &M1.m, &M1.n); // read matrix 1
    M2.M = readm("kronecker/matrix2.txt", &M2.m, &M2.n); // read matrix 2
    printm(M1);                                          // print matrix 1
    printm(M2);                                          // print matrix 2
    K = kronecker(M1, M2);                               // kronecker product
    printm(K);
    M1.M = freem(M1); // Free memory matrix 1
    M2.M = freem(M2); // Free memory matrix 2
    K.M = freem(K);   // Free memory matrix 3

    printf("*** End of kronecker product test ***\n\n");
}

int main(int argc, char *argv[])
{
    GaussianTest();
    KroneckerTest();
    MatrixMulTest();
    ReshapeTest();

    mtrx I = eye(10);
    printm(I);

    return 0;
}