#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include <time.h>
#include "linearalg.h"
#include "finitediff.h"
#include "utils.h"
#include "poisson.h"

int main(int argc, char *argv[])
{
    srand(time(NULL));

    mtrx D1 = Diff1(10, 4);
    printm(D1);
    D1.M = freem(D1);

    return 0;
}