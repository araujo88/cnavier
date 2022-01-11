// Poisson solver library

#ifndef POISSON_H_INCLUDED
#define POISSON_H_INCLUDED

#include "linearalg.h"

#define PI 3.14159265359

double error(mtrx u1, mtrx u2);
mtrx poisson(mtrx f, double dx, double dy, int itmax, double tol);
mtrx poisson_SOR(mtrx f, double dx, double dy, int itmax, double tol, double beta);

#endif // POISSON_H_INCLUDED