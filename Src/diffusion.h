#ifndef DIFFUSION_H
#define DIFFUSION_H

#include "./Math/largeVectorBasis.h"
#include "./Math/matrixElement.h"
#include "./Math/specialmath.h"

//dp/dt
double lossFunction(double p, double B, double q, double m);

double evaluateDiffusionCoefficient(double B, double p, double q, double m);
void solveDiffusion(double**** F, double*** B, double*** Vx, double*** Vy, double*** Vz, double* xgrid, double* ygrid, double* zgrid, double* pgrid, int Nx, int Ny, int Nz, int Nmomentum, double* Fsource, double sourcePower, double time, double** sourceCoords, int Nsource);
double evaluateDiffusionTimeStep(double*** B, double*** Vx, double*** Vy, double*** Vz, double* xgrid, double* ygrid, double* zgrid, double* pgrid, int Nx, int Ny, int Nz, int Nmomentum);
void evaluateDiffusionMatrix(std::vector<MatrixElement>**** matrix, double*** B, double*** Vx, double*** Vy, double*** Vz, double* xgrid, double* ygrid, double* zgrid, double* pgrid, int Nx, int Ny, int Nz, int Nmomentum, double dt);
void diffusionStep(std::vector<MatrixElement>**** matrix, LargeVectorBasis* basis, double**** F, double**** rightPart, double*** B, double*** Vx, double*** Vy, double*** Vz, double* xgrid, double* ygrid, double* zgrid, double* pgrid, int Nx, int Ny, int Nz, int Nmomentum, double dt);

#endif