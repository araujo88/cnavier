// Finite difference library

#ifndef FINITEDIFF_H_INCLUDED
#define FINITEDIFF_H_INCLUDED

#include "linearalg.h"

mtrx Diff1(int n, int o, double dx); // Computes finite-difference matrices for the first derivative
mtrx Diff2(int n, int o, double dx); // Computes finite-difference matrices for the second derivative

#endif // FINITEDIFF_H_INCLUDED