// Fluid dynamics library

#ifndef FLUIDDYN_H_INCLUDED
#define FLUIDDYN_H_INCLUDED

#include "linearalg.h"

mtrx euler(mtrx w, mtrx dwdx, mtrx dwdy, mtrx d2wdx2, mtrx d2wdy2, mtrx u, mtrx v, int Re, double dt); // Euler time-advancement
mtrx continuity(mtrx dudx, mtrx dvdy);                                                                 // computes continuity equation
mtrx pressure(void);                                                                                   // computes pressure

#endif // FLUIDDYN_H_INCLUDED