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
    double Re = 1000.; // Reynolds number
    int Lx = 1;        // length
    int Ly = 1;        // width

    // Numerical parameters
    int nx = 64;                // number of points in x direction
    int ny = 64;                // number of points in y direction
    double dt = 0.004;          // time step
    double tf = 30;             // final time
    double max_co = 1.;         // max Courant number
    int order = 6;              // finite difference order for spatial derivatives
    int poisson_max_it = 10000; // Poisson equation max number of iterations
    double poisson_tol = 1E-9;  // Poisson equation criterion for convergence

    // Boundary conditions (Dirichlet)
    double ui = 0.; // internal field for u
    double vi = 0.; // internal field for v

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

    d_x.M = freem(d_x);
    d_y.M = freem(d_y);
    d_x2.M = freem(d_x2);
    d_y2.M = freem(d_y2);
    Ix.M = freem(Ix);
    Iy.M = freem(Iy);

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

    // Variables
    // Initialize velocities
    mtrx u = initm(nx, ny);   // x-velocity
    mtrx v = initm(nx, ny);   // y-velocity
    mtrx w = initm(nx, ny);   // vorticity
    mtrx psi = initm(nx, ny); // stream-function

    // Derivatives
    // mtrx p;   // pressure
    mtrx dwdx;             // vorticity x-derivative
    mtrx dwdy;             // vorticity y-derivative
    mtrx d2wdx2;           // vorticity x-derivative (2nd)
    mtrx d2wdy2;           // vorticity y-derivative (2nd)
    mtrx dpsidx;           // stream-function x-derivative
    mtrx dpsidy;           // stream-function y-derivative
    mtrx dudx;             // x-velocity x-derivative
    mtrx dudy;             // x-velocity y-derivative
    mtrx dvdx;             // y-velocity x-derivative
    mtrx dvdy;             // y-velocity y-derivative
    mtrx check_continuity; // continuity equation

    // Auxiliary variables
    mtrx w0;
    mtrx dwdx0;
    mtrx dwdy0;
    mtrx d2wdx20;
    mtrx d2wdy20;
    mtrx psi0;
    mtrx dpsidx0;
    mtrx dpsidy0;
    mtrx u0;
    mtrx v0;
    mtrx dudx0;
    mtrx dudy0;
    mtrx dvdx0;
    mtrx dvdy0;

    // Initial condition
    for (i = 1; i < nx - 1; i++)
    {
        for (j = 1; j < ny - 1; j++)
        {
            u.M[i][j] = ui;
            v.M[i][j] = vi;
        }
    }

    // Main time loop
    for (t = 0; t <= it_max; t++)
    {
        // Initialize variables

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

        u0 = reshape(u, nx * ny, 1);
        v0 = reshape(v, nx * ny, 1);
        dudy0 = mtrxmul(DY, u0);
        dvdx0 = mtrxmul(DX, v0);

        dudy = reshape(dudy0, nx, ny);
        dvdx = reshape(dvdx0, nx, ny);

        u0.M = freem(u0);
        v0.M = freem(v0);
        dudy0.M = freem(dudy0);
        dvdx0.M = freem(dvdx0);

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
        w0 = reshape(w, nx * ny, 1);
        dwdx0 = mtrxmul(DX, w0);
        dwdy0 = mtrxmul(DY, w0);

        dwdx = reshape(dwdx0, nx, ny);
        dwdy = reshape(dwdy0, nx, ny);

        dwdx0.M = freem(dwdx0);
        dwdy0.M = freem(dwdy0);

        d2wdx20 = mtrxmul(DX2, w0);
        d2wdy20 = mtrxmul(DY2, w0);

        d2wdx2 = reshape(d2wdx20, nx, ny);
        d2wdy2 = reshape(d2wdy20, nx, ny);

        d2wdx20.M = freem(d2wdx20);
        d2wdy20.M = freem(d2wdy20);
        w0.M = freem(w0);

        // Time - advancement(Euler)
        euler(w, dwdx, dwdy, d2wdx2, d2wdy2, u, v, Re, dt);

        // Solves poisson equation for stream function
        psi.M = freem(psi);
        invsig(w);
        psi = poisson(w, dx, dy, poisson_max_it, poisson_tol);
        invsig(w);

        // Computes velocities
        psi0 = reshape(psi, nx * ny, 1);
        dpsidx0 = mtrxmul(DX, psi0);
        dpsidy0 = mtrxmul(DY, psi0);

        dpsidx = reshape(dpsidx0, nx, ny);
        dpsidy = reshape(dpsidy0, nx, ny);

        psi0.M = freem(psi0);
        dpsidx0.M = freem(dpsidx0);
        dpsidy0.M = freem(dpsidy0);

        u.M = freem(u);
        v.M = freem(v);
        u = initm(nx, ny);
        v = initm(nx, ny);
        mtrxcpy(u, dpsidy);
        invsig(dpsidx);
        mtrxcpy(v, dpsidx);

        // Checks continuity equation

        u0 = reshape(u, nx * ny, 1);
        v0 = reshape(v, nx * ny, 1);

        dudx0 = mtrxmul(DX, u0);
        dvdy0 = mtrxmul(DY, v0);

        dudx = reshape(dudx0, nx, ny);
        dvdy = reshape(dvdy0, nx, ny);
        check_continuity = continuity(dudx, dvdy);
        printf("Iteration: %d | ", t);
        printf("Time: %lf | ", (double)t * dt);
        printf("Progress: %.2lf%%\n", (double)100 * t / it_max);
        printf("Continuity max: %E | ", maxel(check_continuity));
        printf("Continuity min: %E\n", minel(check_continuity));

        u0.M = freem(u0);
        v0.M = freem(v0);
        dudx0.M = freem(dudx0);
        dvdy0.M = freem(dvdy0);

        // Computes pressure
        //	dudx = np.reshape(DX @ np.reshape(u,(nx*ny,1)),(nx,ny))
        //	dudy = np.reshape(DY @ np.reshape(u,(nx*ny,1)),(nx,ny))
        //	dvdx = np.reshape(DX @ np.reshape(v,(nx*ny,1)),(nx,ny))
        //	dvdy = np.reshape(DY @ np.reshape(v,(nx*ny,1)),(nx,ny))
        //	f = dudx**2+dvdy**2+2*dudy*dvdx
        //	p = fft_poisson(-f,dx)

        if (t % output_interval == 0)
        {
            //printvtk(psi, "stream-function");
            printvtk(w, "vorticity");
            // printvtk(u, "x-velocity");
            // printvtk(v, "y-velocity");
            // printvtk(p, "pressure");
        }

        // Free memory
        // freem(p);
        dwdx.M = freem(dwdx);
        dwdy.M = freem(dwdy);
        d2wdx2.M = freem(d2wdx2);
        d2wdy2.M = freem(d2wdy2);
        dpsidx.M = freem(dpsidx);
        dpsidy.M = freem(dpsidy);
        dudx.M = freem(dudx);
        dudy.M = freem(dudy);
        dvdx.M = freem(dvdx);
        dvdy.M = freem(dvdy);
        check_continuity.M = freem(check_continuity);
    }
    // Free memory
    u.M = freem(u);
    v.M = freem(v);
    w.M = freem(w);
    psi.M = freem(psi);
    DX.M = freem(DX);
    DY.M = freem(DY);
    DX2.M = freem(DX2);
    DY2.M = freem(DY2);

    return 0;
}