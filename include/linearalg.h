// Linear algebra library

#ifndef LINEARALG_H_INCLUDED
#define LINEARALG_H_INCLUDED

typedef struct matrix
{
    double **M; // matrix double pointer
    int m;      // number of rows
    int n;      // number of columns
} mtrx;

typedef struct vector
{
    double *v; // vector pointer
    int n;     // number of elements
} vec;

void zerosm(mtrx A);                            // Initialize matrix with zeros
double **allocm(int m, int n);                  // Allocate matrix with size m x n
double **freem(mtrx A);                         // Free memory for matrix
double **readm(char *filename, int *m, int *n); // Read matrix
void printm(mtrx A);                            // Print matrix
void zerosv(vec v);                             // Initialize vector with zeros
double *allocv(int n);                          // Allocate vector with size n
double *freev(vec v);                           // Free memory for vector
double *readv(char *filename, int *n);          // Read vector
void printv(vec v);                             // Print vector;
mtrx mtrxmul(mtrx A, mtrx B, int *m, int *n);   // Matrix multiplication
vec gaussian(mtrx A, vec b);                    // Gaussian elimination
mtrx kronecker(mtrx A, mtrx B);                 // Kronecker matrix product
mtrx reshape(mtrx A, int m, int n);             // reshape matrix
mtrx eye(int n);                                // generates identiy matrix

#endif // LINEARALG_H_INCLUDED