// Poisson solver library

#ifndef POISSON_H_INCLUDED
#define POISSON_H_INCLUDED

#include "linearalg.h"

#define PI 3.14159265359

mtrx poisson(mtrx f, double dx, double dy, int itmax, double tol);

#endif // POISSON_H_INCLUDED