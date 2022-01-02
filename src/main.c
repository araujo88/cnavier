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
    // srand(time(NULL));
    int i, j, t;

    // Physical parameters
    int Re = 100; // Reynolds number
    int Lx = 1;   // length
    int Ly = 1;   // width

    // Numerical parameters
    int nx = 20;       // number of points in x direction
    int ny = 20;       // number of points in y direction
    double dt = 0.02;  // time step
    double tf = 10;    // final time
    double max_co = 1; // max Courant number
    int order = 2;     // finite difference order for spatial derivatives

    // Boundary conditions (Dirichlet)
    double u0 = 0; // internal field for u
    double v0 = 0; // internal field for v

    double u1 = 0; // bottom boundary condition
    double u2 = 1; // top boundary condition
    double u3 = 0; // right boundary condition
    double u4 = 0; // left boundary condition

    double v1 = 0;
    double v2 = 0;
    double v3 = 0;
    double v4 = 0;

    // Computes cell sizes
    double dx = Lx / nx;
    double dy = Ly / ny;

    // Generates derivatives operators
    mtrx d_x = Diff1(nx, order, dx);
    mtrx d_y = Diff1(ny, order, dy);
    mtrx d_x2 = Diff2(nx, order, dx);
    mtrx d_y2 = Diff2(ny, order, dy);

    mtrx Ix = eye(nx);              // identity matrix
    mtrx Iy = eye(ny);              // identity matrix
    mtrx DX = kronecker(d_x, Ix);   // kronecker product for x first derivative
    mtrx DY = kronecker(Iy, d_y);   // kronecker product for y first derivative
    mtrx DX2 = kronecker(d_x2, Ix); // kronecker product for x second derivative
    mtrx DY2 = kronecker(Iy, d_y2); // kronecker product for y second derivative

    // Maximum number of iterations
    int it_max = (int)((tf / dt) - 1);

    // Courant numbers
    double r1 = u1 * dt / (dx);
    double r2 = u1 * dt / (dy);

    if ((r1 > max_co) || (r2 > max_co))
    {
        printf("Unstable Solution!\n");
        printf("r1: %lf\n", r1);
        printf("r2: %lf\n", r2);
        exit(1);
    }

    // Initialize variables
    mtrx u = initm(nx, ny);   // x-velocity
    mtrx v = initm(nx, ny);   // y-velocity
    mtrx w = initm(nx, ny);   // vorticity
    mtrx psi = initm(nx, ny); // stream-function
    mtrx p = initm(nx, ny);   // pressure

    mtrx dwdx = initm(nx, ny);   // vorticity x-derivative
    mtrx dwdy = initm(nx, ny);   // vorticity y-derivative
    mtrx d2wdx2 = initm(nx, ny); // vorticity x-derivative (2nd)
    mtrx d2wdy2 = initm(nx, ny); // vorticity y-derivative (2nd)

    mtrx dpsidx = initm(nx, ny); // stream-function x-derivative
    mtrx dpsidy = initm(nx, ny); // stream-function y-derivative

    mtrx dudx = initm(nx, ny); // x-velocity x-derivative
    mtrx dudy = initm(nx, ny); // x-velocity y-derivative
    mtrx dvdx = initm(nx, ny); // y-velocity x-derivative
    mtrx dvdy = initm(nx, ny); // y-velocity y-derivative

    mtrx continuity = initm(nx, ny);

    // Initial condition
    for (i = 0; i < nx - 1; i++)
    {
        for (j = 1; j < ny - 1; j++)
        {
            u.M[i][j] = u0;
            v.M[i][j] = v0;
        }
    }

    // Main time - loop
    for (t = 0; t < it_max; t++)
    {
        // Boundary conditions
        for (j = 0; j < ny; j++)
        {
            u.M[0][j] = u3;
            u.M[nx - 1][j] = u4;
            v.M[0][j] = v3;
            v.M[nx - 1][j] = v4;
        }
        for (i = 0; i < nx; i++)
        {
            u.M[i][0] = u1;
            u.M[i][ny - 1] = u2;
            v.M[i][0] = v1;
            v.M[i][ny - 1] = v2;
        }

        dudy = reshape(mtrxmul(DY, reshape(u, (nx * ny, 1)), (nx, ny)));
        dvdx = reshape(mtrxmul(DX, reshape(v, (nx * ny, 1)), (nx, ny)));

        for (j = 0; j < ny; j++)
        {
            w[0, j] = dvdx[0, j] - dudy[0, j];
            w[nx - 1, j] = dvdx[nx - 1, j] - dudy[nx - 1, j];
        }
        for (i = 0; i < nx; i++)
        {
            w[i, 0] = dvdx[i, 0] - dudy[i, 0];
            w[i, ny - 1] = dvdx[i, ny - 1] - dudy[i, ny - 1];
        }
        psi = poisson(-w, dx);

        // Computes derivatives
        dwdx = np.reshape(DX @np.reshape(w, (nx * ny, 1)), (nx, ny));
        dwdy = np.reshape(DY @np.reshape(w, (nx * ny, 1)), (nx, ny));
        d2wdx2 = np.reshape(DX2 @np.reshape(w, (nx * ny, 1)), (nx, ny));
        d2wdy2 = np.reshape(DY2 @np.reshape(w, (nx * ny, 1)), (nx, ny));

        // Time - advancement(Euler)
        w = (-u * dwdx - v * dwdy + (1 / Re) * (d2wdx2 + d2wdy2)) * dt + w;

        // Solves poisson equation for stream function
        psi = poisson(-w, dx);

        // Computes velocities
        dpsidx = np.reshape(DX @np.reshape(psi, (nx * ny, 1)), (nx, ny));
        dpsidy = np.reshape(DY @np.reshape(psi, (nx * ny, 1)), (nx, ny));
        u = dpsidy;
        v = -dpsidx;

        // Checks continuity equation
        dudx = np.reshape(DX @np.reshape(u, (nx * ny, 1)), (nx, ny));
        dvdy = np.reshape(DY @np.reshape(v, (nx * ny, 1)), (nx, ny));
        continuity = dudx + dvdy;
        print('Iteration: ' + str(t));
        print('Continuity max: ' + str(continuity.max()) + ' Continuity min: ' + str(continuity.min()));

        // Computes pressure
        //	dudx = np.reshape(DX @ np.reshape(u,(nx*ny,1)),(nx,ny))
        //	dudy = np.reshape(DY @ np.reshape(u,(nx*ny,1)),(nx,ny))
        //	dvdx = np.reshape(DX @ np.reshape(v,(nx*ny,1)),(nx,ny))
        //	dvdy = np.reshape(DY @ np.reshape(v,(nx*ny,1)),(nx,ny))
        //	f = dudx**2+dvdy**2+2*dudy*dvdx
        //	p = fft_poisson(-f,dx)
    }

    return 0;
}