#include "stdio.h"
#include "math.h"

#include "util.h"
#include "constants.h"

#include "diffusion.h"

double lossFunction(double p, double B, double q, double m)
{
	double T1 = 2.725;
	double relictEnergyDensity = 0.25 * 1.6E-12;
	double T2 = 140;
	double photonIRenergyDensity = 1.29E-12;
	double magneticEnergyDensity = B * B / (8 * pi);
	double sigmaT = (8.0 * pi / 3.0) * sqr(q * q / (m * speed_of_light2));
	double coef = (4.0 / 3.0) * sigmaT / (m * m * speed_of_light * speed_of_light * speed_of_light);

	double dpdt = speed_of_light * coef * p * p *m *m *speed_of_light*speed_of_light * (magneticEnergyDensity + relictEnergyDensity / pow(1.0 + 4 * 2.8 * kBoltzman * T1 * p / (m * m * speed_of_light * speed_of_light * speed_of_light), 1.5) + photonIRenergyDensity / pow(1.0 + 4 * 2.8 * kBoltzman * T2 * p / (m * m * speed_of_light * speed_of_light * speed_of_light), 1.5));
	return dpdt/(massElectron*speed_of_light);
}

double evaluateDiffusionCoefficient(double B, double p, double q, double m)
{
	//todo
	if (B < 5E-7) {
		B = 5E-7;
	}
	return p*m*speed_of_light*speed_of_light2/(3*q*B);
}

void solveDiffusion(double**** F, double*** B, double*** Vx, double*** Vy, double*** Vz, double* xgrid, double* ygrid, double* zgrid, double* pgrid, int Nx, int Ny, int Nz, int Nmomentum, double* Fsource, double sourcePower, double time, double** sourceCoords, int Nsource)
{
	double**** rightPart = create4dArray(Nx, Ny, Nz, Nmomentum);
	std::vector<MatrixElement>**** matrix = new std::vector<MatrixElement>***[Nx];
	for (int i = 0; i < Nx; ++i) {
		matrix[i] = new std::vector<MatrixElement>**[Ny];
		for (int j = 0; j < Ny; ++j) {
			matrix[i][j] = new std::vector<MatrixElement>*[Nz];
			for (int k = 0; k < Nz; ++k) {
				matrix[i][j][k] = new std::vector<MatrixElement>[Nmomentum];
			}
		}
	}
	LargeVectorBasis* basis = new LargeVectorBasis(10, Nx, Ny, Nz, Nmomentum);
	double dt = evaluateDiffusionTimeStep(B, Vx, Vy, Vz, xgrid, ygrid, zgrid, pgrid, Nx, Ny, Nz, Nmomentum);

	int** sourceIndex = new int* [Nsource];
	for (int i = 0; i < Nsource; ++i) {
		sourceIndex[i] = new int[3];
		sourceIndex[i][0] = 0;
		sourceIndex[i][1] = 0;
		sourceIndex[i][2] = 0;
	}

	//todo middpoint
	for (int l = 0; l < Nsource; ++l) {
		for (int i = 0; i < Nx - 1; ++i) {
			if ((xgrid[i] <= sourceCoords[l][0]) && (xgrid[i + 1] > sourceCoords[l][0])) {
				sourceIndex[l][0] = i;
			}
		}
		for (int j = 0; j < Ny - 1; ++j) {
			if ((ygrid[j] <= sourceCoords[l][1]) && (ygrid[j + 1] > sourceCoords[l][1])) {
				sourceIndex[l][1] = j;
			}
		}
		for (int k = 0; k < Nz - 1; ++k) {
			if ((zgrid[k] <= sourceCoords[l][2]) && (zgrid[k + 1] > sourceCoords[l][2])){ 
				sourceIndex[l][2] = k;
			}
		}
	}

	int Niterations = floor(time / dt) + 1;
	dt = time / Niterations;
	printf("number of diffusion iterations = %d\n", Niterations);
	printLog("number of diffusion iterations = %d\n", Niterations);

	evaluateDiffusionMatrix(matrix, B, Vx, Vy, Vz, xgrid, ygrid, zgrid, pgrid, Nx, Ny, Nz, Nmomentum, dt);

	for (int i = 0; i < Niterations; ++i) {
		printf("diffusion solver iteration %d\n", i);
		printLog("diffusion solver iteration %d\n", i);
		for (int l = 0; l < Nsource; ++l) {
			double dx = 0;
			if (sourceIndex[l][0] == 0) {
				dx = xgrid[1] - xgrid[0];
			}
			else {
				dx = xgrid[sourceIndex[l][0]] - xgrid[sourceIndex[l][0] - 1];
			}
			double dy = 0;
			if (sourceIndex[l][1] == 0) {
				dy = ygrid[1] - ygrid[0];
			}
			else {
				dy = ygrid[sourceIndex[l][1]] - ygrid[sourceIndex[l][1] - 1];
			}
			double dz = 0;
			if (sourceIndex[l][2] == 0) {
				dz = zgrid[1] - zgrid[0];
			}
			else {
				dz = zgrid[sourceIndex[l][2]] - zgrid[sourceIndex[l][2] - 1];
			}
			for (int m = 0; m < Nmomentum; ++m) {
				F[sourceIndex[l][0]][sourceIndex[l][1]][sourceIndex[l][2]][m] += dt * sourcePower * Fsource[m]/(dx*dy*dz);
			}
		}
		diffusionStep(matrix, basis, F, rightPart, B, Vx, Vy, Vz, xgrid, ygrid, zgrid, pgrid, Nx, Ny, Nz, Nmomentum, dt);
	}

	delete4dArray(rightPart, Nx, Ny, Nz, Nmomentum);
	for (int i = 0; i < Nx; ++i) {
		for (int j = 0; j < Ny; ++j) {
			for (int k = 0; k < Nz; ++k) {
				delete[] matrix[i][j][k];
			}
			delete[] matrix[i][j];
		}
		delete[] matrix[i];
	}
	delete[] matrix;

	delete basis;
}

double evaluateDiffusionTimeStep(double*** B, double*** Vx, double*** Vy, double*** Vz, double* xgrid, double* ygrid, double* zgrid, double* pgrid, int Nx, int Ny, int Nz, int Nmomentum)
{
	double inv_dt = 1E-7;
	for (int i = 1; i < Nx-1; ++i) {
		double dx = 0;
		if (i == 0) {
			dx = xgrid[1] - xgrid[0];
		}
		else {
			dx = xgrid[i] - xgrid[i - 1];
		}
		for (int j = 1; j < Ny-1; ++j) {
			double dy = 0;
			if (j == 0) {
				dy = ygrid[1] - ygrid[0];
			}
			else {
				dy = ygrid[j] - ygrid[j - 1];
			}
			for (int k = 1; k < Nz-1; ++k) {
				double dz = 0;
				if (k == 0) {
					dz = zgrid[1] - zgrid[0];
				}
				else {
					dz = zgrid[k] - zgrid[k - 1];
				}

				double divu = 0;
				divu += (Vx[i+1][j][k] - Vx[i-1][j][k]) / (xgrid[i + 1] - xgrid[i - 1]);
				divu += (Vy[i][j+1][k] - Vy[i][j-1][k]) / (ygrid[j + 1] - ygrid[j - 1]);
				divu += (Vz[i][j][k+1] - Vz[i][j][k-1]) / (zgrid[k + 1] - zgrid[k - 1]);


				double inv_dt_new = fabs(4 * Vx[i][j][k] / dx);
				inv_dt = max(inv_dt, inv_dt_new);
				inv_dt_new = fabs(4 * Vy[i][j][k] / dy);
				inv_dt = max(inv_dt, inv_dt_new);
				inv_dt_new = fabs(4 * Vz[i][j][k] / dz);
				inv_dt = max(inv_dt, inv_dt_new);

				for (int l = 0; l < Nmomentum; ++l) {
					double dp = 0;
					if (l == 0) {
						dp = pgrid[1] - pgrid[0];
					}
					else {
						dp = pgrid[l] - pgrid[l - 1];
					}

					inv_dt_new = fabs(100 * divu * pgrid[l] / dp);
					inv_dt = max(inv_dt, inv_dt_new);

					double D = evaluateDiffusionCoefficient(B[i][j][k], pgrid[l], electron_charge, massElectron);

					/*inv_dt_new = fabs(4 * 4 * D / (dx * dx));
					inv_dt = max(inv_dt, inv_dt_new);
					inv_dt_new = fabs(4 * 4 * D / (dy * dy));
					inv_dt = max(inv_dt, inv_dt_new);
					inv_dt_new = fabs(4 * 4 * D / (dz * dz));
					inv_dt = max(inv_dt, inv_dt_new);*/

					/*double lossRate = lossFunction(pgrid[l], B[i][j][k], electron_charge, massElectron);
					inv_dt_new = fabs(4 * lossRate / dp);
					inv_dt = max(inv_dt, inv_dt_new);*/
				}
			}
		}
	}
	return 1.0/inv_dt;
}

void evaluateDiffusionMatrix(std::vector<MatrixElement>**** matrix, double*** B, double*** Vx, double*** Vy, double*** Vz, double* xgrid, double* ygrid, double* zgrid, double* pgrid, int Nx, int Ny, int Nz, int Nmomentum, double dt)
{
	for (int i = 0; i < Nx; ++i) {
		for (int j = 0; j < Ny; ++j) {
			for (int k = 0; k < Nz; ++k) {
				for (int l = 0; l < Nmomentum; ++l) {
					matrix[i][j][k][l].clear();
					matrix[i][j][k][l].push_back(MatrixElement(1.0, i, j, k, l));
				}
			}
		}
	}

	//zero boundary conditions
	printf("writing matrix\n");
	printLog("writing matrix\n");
	int i;
	//#pragma omp parallel for private(i) shared(matrix, rightPart, F, B, Vx, Vy, Vz, xgrid, ygrid, zgrid, pgrid, Nx, Ny, Nz, Nmomentum)
	for (int i = 1; i < Nx - 1; ++i) {
		printf("matrix i %d\n", i);
		printLog("matrix i %d\n", i);
		for (int j = 1; j < Ny - 1; ++j) {
			for (int k = 1; k < Nz - 1; ++k) {
				double divu = 0;
				divu += (Vx[i + 1][j][k] - Vx[i - 1][j][k]) / (xgrid[i + 1] - xgrid[i - 1]);
				divu += (Vy[i][j + 1][k] - Vy[i][j - 1][k]) / (ygrid[j + 1] - ygrid[j - 1]);
				divu += (Vz[i][j][k + 1] - Vz[i][j][k - 1]) / (zgrid[k + 1] - zgrid[k - 1]);
				for (int l = 0; l < Nmomentum - 1; ++l) {
					//printf("matrix %d %d %d %d\n", i, j, k, l);
					//printLog("matrix %d %d %d %d\n", i, j, k, l);

					double D = evaluateDiffusionCoefficient(B[i][j][k], pgrid[l], electron_charge, massElectron);
					double Dleft = evaluateDiffusionCoefficient(B[i - 1][j][k], pgrid[l], electron_charge, massElectron);
					double Dright = evaluateDiffusionCoefficient(B[i + 1][j][k], pgrid[l], electron_charge, massElectron);

					double Dbottom = evaluateDiffusionCoefficient(B[i][j - 1][k], pgrid[l], electron_charge, massElectron);
					double Dtop = evaluateDiffusionCoefficient(B[i][j + 1][k], pgrid[l], electron_charge, massElectron);

					double Dback = evaluateDiffusionCoefficient(B[i][j][k - 1], pgrid[l], electron_charge, massElectron);
					double Dfront = evaluateDiffusionCoefficient(B[i][j][k + 1], pgrid[l], electron_charge, massElectron);

					double dp = 0;
					if (divu > 0) {
						if (l == Nmomentum - 1) {
							dp = pgrid[Nmomentum - 1] - pgrid[Nmomentum - 2];
						}
						else {
							dp = pgrid[l + 1] - pgrid[l];
						}
					}
					else {
						if (l == 0) {
							dp = pgrid[1] - pgrid[0];
						}
						else {
							dp = pgrid[l] - pgrid[l - 1];
						}
					}

					double value = 0;

					if (divu >= 0) {
						value = dt * divu * pgrid[l] / (3.0 * dp);
						matrix[i][j][k][l].push_back(MatrixElement(value, i, j, k, l));
						if ((value != value) || (value * 0 != value * 0)) {
							printf("1 value = NaN for (%d %d %d %d) (%d %d %d %d)\n", i, j, k, l, i, j, k, l);
							exit(0);
						}
						if (l < Nmomentum - 1) {
							value = -dt * divu * pgrid[l] / (3.0 * dp);
							matrix[i][j][k][l].push_back(MatrixElement(value, i, j, k, l + 1));
							if ((value != value) || (value * 0 != value * 0)) {
								printf("2 value = NaN for (%d %d %d %d) (%d %d %d %d)\n", i, j, k, l, i, j, k, l + 1);
								exit(0);
							}
						}
					}
					else {
						value = -dt * divu * pgrid[l] / (3.0 * dp);
						matrix[i][j][k][l].push_back(MatrixElement(value, i, j, k, l));
						if ((value != value) || (value * 0 != value * 0)) {
							printf("3 value = NaN for (%d %d %d %d) (%d %d %d %d)\n", k, j, i, l, i, j, k, l);
							exit(0);
						}
						if (l > 0) {
							value = dt * divu * pgrid[l] / (3.0 * dp);
							matrix[i][j][k][l].push_back(MatrixElement(value, i, j, k, l - 1));
							if ((value != value) || (value * 0 != value * 0)) {
								printf("4 value = NaN for (%d %d %d %d) (%d %d %d %d)\n", i, j, k, l, i, j, k, l - 1);
								exit(0);
							}
						}
					}

					double lossRate = lossFunction(pgrid[l], B[i][j][k], electron_charge, massElectron);
					if (lossRate > 0) {
						if (l < Nmomentum - 1) {
							dp = pgrid[l + 1] - pgrid[l];
						}
						else {
							dp = pgrid[Nmomentum - 1] - pgrid[Nmomentum - 2];
						}
					}
					else {
						if (l > 0) {
							dp = pgrid[l] - pgrid[l - 1];
						}
						else {
							dp = pgrid[1] - pgrid[0];
						}

					}

					if (lossRate > 0) {
						value = dt * lossRate / dp;
						matrix[i][j][k][l].push_back(MatrixElement(value, i, j, k, l));
						if ((value != value) || (value * 0 != value * 0)) {
							printf("1 value = NaN for (%d %d %d %d) (%d %d %d %d)\n", i, j, k, l, i, j, k, l);
							exit(0);
						}
						if (l < Nmomentum - 1) {
							value = -dt * lossRate * pgrid[l] / (pgrid[l + 1] * dp);
							matrix[i][j][k][l].push_back(MatrixElement(value, i, j, k, l + 1));
							if ((value != value) || (value * 0 != value * 0)) {
								printf("1 value = NaN for (%d %d %d %d) (%d %d %d %d)\n", i, j, k, l, i, j, k, l + 1);
								exit(0);
							}
						}
					}
					else {
						value = -dt * lossRate / dp;
						matrix[i][j][k][l].push_back(MatrixElement(value, i, j, k, l));
						if ((value != value) || (value * 0 != value * 0)) {
							printf("1 value = NaN for (%d %d %d %d) (%d %d %d %d)\n", i, j, k, l, i, j, k, l);
							exit(0);
						}
						if (l > 0) {
							value = dt * lossRate * pgrid[l] / (pgrid[l - 1] * dp);
							matrix[i][j][k][l].push_back(MatrixElement(value, i, j, k, l - 1));
							if ((value != value) || (value * 0 != value * 0)) {
								printf("1 value = NaN for (%d %d %d %d) (%d %d %d %d)\n", i, j, k, l, i, j, k, l - 1);
								exit(0);
							}
						}
					}

					value = dt * (((Dright + D) / (xgrid[i + 1] - xgrid[i])) + ((D + Dleft) / (xgrid[i] - xgrid[i - 1]))) / (xgrid[i + 1] - xgrid[i - 1]);
					matrix[i][j][k][l].push_back(MatrixElement(value, i, j, k, l));
					if ((value != value) || (value * 0 != value * 0)) {
						printf("5 value = NaN for (%d %d %d %d) (%d %d %d %d)\n", i, j, k, l, i, j, k, l);
						exit(0);
					}


					value = -dt * ((Dright + D) / (xgrid[i + 1] - xgrid[i])) / (xgrid[i + 1] - xgrid[i - 1]);
					matrix[i][j][k][l].push_back(MatrixElement(value, i + 1, j, k, l));
					if ((value != value) || (value * 0 != value * 0)) {
						printf("6 value = NaN for (%d %d %d %d) (%d %d %d %d)\n", i, j, k, l, i + 1, j, k, l);
						exit(0);
					}


					value = -dt * ((D + Dleft) / (xgrid[i] - xgrid[i - 1])) / (xgrid[i + 1] - xgrid[i - 1]);
					matrix[i][j][k][l].push_back(MatrixElement(value, i - 1, j, k, l));
					if ((value != value) || (value * 0 != value * 0)) {
						printf("7 value = NaN for (%d %d %d %d) (%d %d %d %d)\n", i, j, k, l, i - 1, j, k, l);
						exit(0);
					}


					if (Vx[i][j][k] > 0) {
						value = dt * Vx[i][j][k] / (xgrid[i] - xgrid[i - 1]);
						matrix[i][j][k][l].push_back(MatrixElement(value, i, j, k, l));
						if ((value != value) || (value * 0 != value * 0)) {
							printf("8 value = NaN for (%d %d %d %d) (%d %d %d %d)\n", i, j, k, l, i, j, k, l);
							exit(0);
						}

						value = -dt * Vx[i - 1][j][k] / (xgrid[i] - xgrid[i - 1]);
						matrix[i][j][k][l].push_back(MatrixElement(value, i - 1, j, k, l));
						if ((value != value) || (value * 0 != value * 0)) {
							printf("9 value = NaN for (%d %d %d %d) (%d %d %d %d)\n", i, j, k, l, i - 1, j, k, l);
							exit(0);
						}

					}
					else {
						value = -dt * Vx[i][j][k] / (xgrid[i + 1] - xgrid[i]);
						matrix[i][j][k][l].push_back(MatrixElement(value, i, j, k, l));
						if ((value != value) || (value * 0 != value * 0)) {
							printf("10 value = NaN for (%d %d %d %d) (%d %d %d %d)\n", i, j, k, l, i, j, k, l);
							exit(0);
						}

						value = dt * Vx[i + 1][j][k] / (xgrid[i + 1] - xgrid[i]);
						matrix[i][j][k][l].push_back(MatrixElement(value, i + 1, j, k, l));
						if ((value != value) || (value * 0 != value * 0)) {
							printf("11 value = NaN for (%d %d %d %d) (%d %d %d %d)\n", i, j, k, l, i + 1, j, k, l);
							exit(0);
						}

					}

					value = dt * (((Dtop + D) / (ygrid[j + 1] - ygrid[j])) + ((D + Dbottom) / (ygrid[j] - ygrid[j - 1]))) / (ygrid[j + 1] - ygrid[j - 1]);
					matrix[i][j][k][l].push_back(MatrixElement(value, i, j, k, l));
					if ((value != value) || (value * 0 != value * 0)) {
						printf("12 value = NaN for (%d %d %d %d) (%d %d %d %d)\n", i, j, j, l, i, j, k, l);
						exit(0);
					}


					value = -dt * ((Dtop + D) / (ygrid[j + 1] - ygrid[j])) / (ygrid[j + 1] - ygrid[j - 1]);
					matrix[i][j][k][l].push_back(MatrixElement(value, i, j + 1, k, l));
					if ((value != value) || (value * 0 != value * 0)) {
						printf("13 value = NaN for (%d %d %d %d) (%d %d %d %d)\n", i, j, k, l, i, j + 1, k, l);
						exit(0);
					}


					value = -dt * ((D + Dbottom) / (ygrid[j] - ygrid[j - 1])) / (ygrid[j + 1] - ygrid[j - 1]);
					matrix[i][j][k][l].push_back(MatrixElement(value, i, j - 1, k, l));
					if ((value != value) || (value * 0 != value * 0)) {
						printf("14 value = NaN for (%d %d %d %d) (%d %d %d %d)\n", i, j, k, l, i, j - 1, k, l);
						exit(0);
					}


					if (Vy[i][j][k] > 0) {
						value = dt * Vy[i][j][k] / (ygrid[j] - ygrid[j - 1]);
						matrix[i][j][k][l].push_back(MatrixElement(value, i, j, j, l));
						if ((value != value) || (value * 0 != value * 0)) {
							printf("15 value = NaN for (%d %d %d %d) (%d %d %d %d)\n", i, j, k, l, i, j, k, l);
							exit(0);
						}

						value = -dt * Vy[i][j - 1][k] / (ygrid[j] - ygrid[j - 1]);
						matrix[i][j][k][l].push_back(MatrixElement(value, i, j - 1, k, l));
						if ((value != value) || (value * 0 != value * 0)) {
							printf("16 value = NaN for (%d %d %d %d) (%d %d %d %d)\n", i, j, k, l, i, j - 1, k, l);
							exit(0);
						}

					}
					else {
						value = -dt * Vy[i][j][k] / (ygrid[j + 1] - ygrid[j]);
						matrix[i][j][k][l].push_back(MatrixElement(value, i, j, k, l));
						if ((value != value) || (value * 0 != value * 0)) {
							printf("17 value = NaN for (%d %d %d %d) (%d %d %d %d)\n", i, j, k, l, i, j, k, l);
							exit(0);
						}

						value = dt * Vy[i][j + 1][k] / (ygrid[j + 1] - ygrid[j]);
						matrix[i][j][k][l].push_back(MatrixElement(value, i, j + 1, k, l));
						if ((value != value) || (value * 0 != value * 0)) {
							printf("18 value = NaN for (%d %d %d %d) (%d %d %d %d)\n", i, j, k, l, i, j + 1, k, l);
							exit(0);
						}

					}
					value = dt * (((Dfront + D) / (zgrid[k + 1] - zgrid[k])) + ((D + Dback) / (zgrid[k] - zgrid[k - 1]))) / (zgrid[k + 1] - zgrid[k - 1]);
					matrix[i][j][k][l].push_back(MatrixElement(value, i, j, k, l));
					if ((value != value) || (value * 0 != value * 0)) {
						printf("19 value = NaN for (%d %d %d %d) (%d %d %d %d)\n", k, j, i, l, i, j, k, l);
						exit(0);
					}


					value = -dt * ((Dfront + D) / (zgrid[k + 1] - zgrid[k])) / (zgrid[k + 1] - zgrid[k - 1]);
					matrix[i][j][k][l].push_back(MatrixElement(value, i, j, k + 1, l));
					if ((value != value) || (value * 0 != value * 0)) {
						printf("20 value = NaN for (%d %d %d %d) (%d %d %d %d)\n", i, j, k, l, i, j, k + 1, l);
						exit(0);
					}


					value = -dt * ((D + Dback) / (zgrid[k] - zgrid[k - 1])) / (zgrid[k + 1] - zgrid[k - 1]);
					matrix[i][j][k][l].push_back(MatrixElement(value, i, j, k - 1, l));
					if ((value != value) || (value * 0 != value * 0)) {
						printf("21 value = NaN for (%d %d %d %d) (%d %d %d %d)\n", i, j, k, l, i, j, k - 1, l);
						exit(0);
					}


					if (Vz[i][j][k] > 0) {
						value = dt * Vz[i][j][k] / (zgrid[k] - zgrid[k - 1]);
						matrix[i][j][k][l].push_back(MatrixElement(value, i, j, k, l));
						if ((value != value) || (value * 0 != value * 0)) {
							printf("22 value = NaN for (%d %d %d %d) (%d %d %d %d)\n", i, j, k, l, i, j, k, l);
							exit(0);
						}

						value = -dt * Vz[i][j][k - 1] / (zgrid[k] - zgrid[k - 1]);
						matrix[i][j][k][l].push_back(MatrixElement(value, i, j, k - 1, l));
						if ((value != value) || (value * 0 != value * 0)) {
							printf("23 value = NaN for (%d %d %d %d) (%d %d %d %d)\n", i, j, k, l, i, j, k - 1, l);
							exit(0);
						}

					}
					else {
						value = -dt * Vz[i][j][k] / (zgrid[k + 1] - zgrid[k]);
						matrix[i][j][k][l].push_back(MatrixElement(value, i, j, k, l));
						if ((value != value) || (value * 0 != value * 0)) {
							printf("24 value = NaN for (%d %d %d %d) (%d %d %d %d)\n", i, j, k, l, i, j, k, l);
							exit(0);
						}

						value = dt * Vz[i][j][k + 1] / (zgrid[k + 1] - zgrid[k]);
						matrix[i][j][k][l].push_back(MatrixElement(value, i, j, k + 1, l));
						if ((value != value) || (value * 0 != value * 0)) {
							printf("25 value = NaN for (%d %d %d %d) (%d %d %d %d)\n", i, j, k, l, i, j, k + 1, l);
							exit(0);
						}

					}
					/*double diagonal = 0;
					double notdiagonal = 0;

					for(int ii = 0; ii < matrix[i][j][k][l].size(); ++ii){
						MatrixElement element = matrix[i][j][k][l][ii];
						if ((element.i == i) && (element.j == j) && (element.k == k) && (element.l == l)) {
							diagonal = diagonal + element.value;
						}
						else {
							notdiagonal = notdiagonal + fabs(element.value);
						}
					}

					if (diagonal < notdiagonal) {
						printf("diagonal = %g < notdiagonal = %g, %d %d %d %d\n", diagonal, notdiagonal, i, j, k, l);
						//exit(0);
					}*/
				}
			}
		}
	}
}

void diffusionStep(std::vector<MatrixElement>**** matrix, LargeVectorBasis* basis, double**** F, double**** rightPart, double*** B, double*** Vx, double*** Vy, double*** Vz, double* xgrid, double* ygrid, double* zgrid, double* pgrid, int Nx, int Ny, int Nz, int Nmomentum, double dt)
{
	for(int i = 0; i < Nx; ++i) {
		for (int j = 0; j < Ny; ++j) {
			for (int k = 0; k < Nz; ++k) {
				for (int l = 0; l < Nmomentum; ++l) {
					rightPart[i][j][k][l] = 0.0;
				}
			}
		}
	}

	//zero boundary conditions
	printf("writing right part\n");
	printLog("writing right part\n");
	int i;
//#pragma omp parallel for private(i) shared(matrix, rightPart, F, B, Vx, Vy, Vz, xgrid, ygrid, zgrid, pgrid, Nx, Ny, Nz, Nmomentum)
	for (int i = 1; i < Nx - 1; ++i) {
		for (int j = 1; j < Ny - 1; ++j) {
			for (int k = 1; k < Nz - 1; ++k) {
				for (int l = 0; l < Nmomentum - 1; ++l) {
					//printf("matrix %d %d %d %d\n", i, j, k, l);
					//printLog("matrix %d %d %d %d\n", i, j, k, l);
					rightPart[i][j][k][l] = F[i][j][k][l];
				}
			}
		}
	}

	printf("solving matrix\n");
	printLog("solving matrix\n");

	double precision = 0.1 / pgrid[Nmomentum - 1];
	precision = 1E-5;
	generalizedMinimalResidualMethod(matrix, rightPart, F, Nx, Ny, Nz, Nmomentum, precision, 50, 1, basis);

	for (i = 0; i < Nx; ++i) {
		for (int j = 0; j < Ny; ++j) {
			for (int k = 0; k < Nz; ++k) {
				for (int l = 0; l < Nmomentum; ++l) {
					if (F[i][j][k][l] < 0) {
						F[i][j][k][l] = 0.0;
					}
				}
			}
		}
	}
}
