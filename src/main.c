#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "linearalg.h"
#include "finitediff.h"
#include "utils.h"
#include "poisson.h"
#include "fluiddyn.h"

int main(int argc, char *argv[])
{
    // srand(time(NULL));
    int i, j, t;

    // Physical parameters
    double Re = 100.; // Reynolds number
    int Lx = 1;       // length
    int Ly = 1;       // width

    // Numerical parameters
    int nx = 32;                // number of points in x direction
    int ny = 32;                // number of points in y direction
    double dt = 0.005;          // time step
    double tf = 10;             // final time
    double max_co = 1.;         // max Courant number
    int order = 6;              // finite difference order for spatial derivatives
    int poisson_max_it = 10000; // Poisson equation max number of iterations
    double poisson_tol = 1E-6;  // Poisson equation criterion for convergence

    // Boundary conditions (Dirichlet)
    double u0 = 0.; // internal field for u
    double v0 = 0.; // internal field for v

    double u1 = 0.; // right boundary condition
    double u2 = 0.; // left boundary condition
    double u3 = 0.; // bottom boundary condition
    double u4 = 1.; // top boundary condition

    double v1 = 0.;
    double v2 = 0.;
    double v3 = 0.;
    double v4 = 0.;

    // Computes cell sizes
    double dx = (double)Lx / nx;
    double dy = (double)Ly / ny;

    // Generates derivatives operators
    mtrx d_x = Diff1(nx, order, dx);
    mtrx d_y = Diff1(ny, order, dy);
    mtrx d_x2 = Diff2(nx, order, dx);
    mtrx d_y2 = Diff2(ny, order, dy);

    mtrx Ix = eye(nx);              // identity matrix
    mtrx Iy = eye(ny);              // identity matrix
    mtrx DX = kronecker(Ix, d_x);   // kronecker product for x first derivative
    mtrx DY = kronecker(d_y, Iy);   // kronecker product for y first derivative
    mtrx DX2 = kronecker(Ix, d_x2); // kronecker product for x second derivative
    mtrx DY2 = kronecker(d_y2, Iy); // kronecker product for y second derivative

    // Maximum number of iterations
    int it_max = (int)((tf / dt) - 1);
    int output_interval = 10;

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

    mtrx check_continuity = initm(nx, ny);

    // Initial condition
    for (i = 1; i < nx - 1; i++)
    {
        for (j = 1; j < ny - 1; j++)
        {
            u.M[i][j] = u0;
            v.M[i][j] = v0;
        }
    }

    // Main time loop
    for (t = 0; t <= it_max; t++)
    {
        // Boundary conditions
        for (j = 0; j < ny; j++)
        {
            v.M[0][j] = v3;
            v.M[nx - 1][j] = v4;
            u.M[0][j] = u3;
            u.M[nx - 1][j] = u4;
        }
        for (i = 0; i < nx; i++)
        {
            v.M[i][0] = v1;
            v.M[i][ny - 1] = v2;
            u.M[i][0] = u1;
            u.M[i][ny - 1] = u2;
        }

        dudy = reshape(mtrxmul(DY, reshape(u, nx * ny, 1)), nx, ny);
        dvdx = reshape(mtrxmul(DX, reshape(v, nx * ny, 1)), nx, ny);

        //w = vorticity(dudy, dvdx);

        for (j = 0; j < ny; j++)
        {
            w.M[0][j] = dvdx.M[0][j] - dudy.M[0][j];
            w.M[nx - 1][j] = dvdx.M[nx - 1][j] - dudy.M[nx - 1][j];
        }
        for (i = 0; i < nx; i++)
        {
            w.M[i][0] = dvdx.M[i][0] - dudy.M[i][0];
            w.M[i][ny - 1] = dvdx.M[i][ny - 1] - dudy.M[i][ny - 1];
        }

        // Computes derivatives
        mtrx temp;
        temp = initm(nx * ny, 1);
        temp = reshape(w, nx * ny, 1);

        dwdx = mtrxmul(DX, temp);
        dwdx = reshape(dwdx, nx, ny);

        dwdy = mtrxmul(DY, temp);
        dwdy = reshape(dwdy, nx, ny);

        d2wdx2 = mtrxmul(DX2, temp);
        d2wdx2 = reshape(d2wdx2, nx, ny);

        d2wdy2 = mtrxmul(DY2, temp);
        d2wdy2 = reshape(d2wdy2, nx, ny);
        freem(temp);

        // Time - advancement(Euler)
        euler(w, dwdx, dwdy, d2wdx2, d2wdy2, u, v, Re, dt);

        // Solves poisson equation for stream function
        psi = poisson(invsig(w), dx, dy, poisson_max_it, poisson_tol);

        // Computes velocities
        mtrx temp2;
        temp2 = initm(nx, ny);
        temp2 = reshape(psi, nx * ny, 1);

        dpsidx = mtrxmul(DX, temp2);
        dpsidx = reshape(dpsidx, nx, ny);
        dpsidy = mtrxmul(DY, temp2);
        dpsidy = reshape(dpsidy, nx, ny);

        freem(temp2);

        mtrxcpy(u, dpsidy);
        mtrxcpy(v, invsig(dpsidx));

        // Checks continuity equation
        dudx = reshape(mtrxmul(DX, reshape(u, nx * ny, 1)), nx, ny);
        dvdy = reshape(mtrxmul(DY, reshape(v, nx * ny, 1)), nx, ny);
        check_continuity = continuity(dudx, dvdy);
        printf("Iteration: %d\n", t);
        printf("Continuity max: %E\n", maxel(check_continuity));
        printf("Continuity min: %E\n", minel(check_continuity));

        // Computes pressure
        //	dudx = np.reshape(DX @ np.reshape(u,(nx*ny,1)),(nx,ny))
        //	dudy = np.reshape(DY @ np.reshape(u,(nx*ny,1)),(nx,ny))
        //	dvdx = np.reshape(DX @ np.reshape(v,(nx*ny,1)),(nx,ny))
        //	dvdy = np.reshape(DY @ np.reshape(v,(nx*ny,1)),(nx,ny))
        //	f = dudx**2+dvdy**2+2*dudy*dvdx
        //	p = fft_poisson(-f,dx)

        if (t % output_interval == 0)
        {
            printvtk(psi, "stream-function");
            // printvtk(w, "vorticity");
            // printvtk(u, "x-velocity");
            // printvtk(v, "y-velocity");
            // printvtk(p, "pressure");
        }
    }

    // Free memory
    freem(d_x);
    freem(d_y);
    freem(d_x2);
    freem(d_y2);
    freem(Ix);
    freem(Iy);
    freem(DX);
    freem(DY);
    freem(DX2);
    freem(DY2);
    freem(u);
    freem(v);
    freem(w);
    freem(psi);
    freem(p);
    freem(dwdx);
    freem(dwdy);
    freem(d2wdx2);
    freem(d2wdy2);
    freem(dpsidx);
    freem(dpsidy);
    freem(dudx);
    freem(dudy);
    freem(dvdx);
    freem(dvdy);
    freem(check_continuity);

    return 0;
}