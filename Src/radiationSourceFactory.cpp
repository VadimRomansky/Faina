#include "stdio.h"
#include "math.h"

#include "constants.h"
#include "radiationSource.h"
#include "util.h"
#include "massiveParticleDistribution.h"

#include "radiationSourceFactory.h"
#include <queue>
#include <tuple>

double RadiationSourceFactory::evaluateTurbulenceAmplitude(const double& k, const double& turbulenceKoef, const double& index, const double& L0)
{
	double k0 = 2 * pi / L0;
	return turbulenceKoef * pow(k0 / k, index);
}

double RadiationSourceFactory::evaluateAnisotropicTurbulenceAmplitude(const double& kx, const double& ky, const double& kz, const double& turbulenceKoef, const double& index, const double& L0, const double& anisotropy)
{
	//if(anisotropy > 1 -> Bz < Bx and By)
	double k0 = 2 * pi / L0;
	double k = sqrt(kz * kz + anisotropy * anisotropy * (kx * kx + ky * ky));
	return turbulenceKoef * pow(k0 / k, index);
}

void RadiationSourceFactory::normalizeTurbulenceKoef(double& turbulenceKoef, const double& index, const double& L0, const int Nmodes, const double& fraction, const double& B0)
{
	double energy = 0;
	for (int i = 0; i < Nmodes; ++i) {
		double kx = i * 2 * pi / L0;
		for (int j = 0; j < Nmodes; ++j) {
			double ky = j * 2 * pi / L0;
			for (int k = 0; k < Nmodes; ++k) {
				double kz = k * 2 * pi / L0;
				//if (i + j + k > 0) {
				if (i + j + k > minModeNumber) {
					double kt = sqrt(kx * kx + ky * ky + kz * kz);
					energy += sqr(evaluateTurbulenceAmplitude(kt, turbulenceKoef, index, L0));
				}
			}
		}
	}

	turbulenceKoef = turbulenceKoef * sqrt(B0 * B0 * fraction / energy);
}

void RadiationSourceFactory::normalizeAnisotropicTurbulenceKoef(double& turbulenceKoef, const double& index, const double& L0, const int Nmodes, const double& fraction, const double& B0, const double& anisotropy) {
	double energy = 0;
	for (int i = 0; i < Nmodes; ++i) {
		double kx = i * 2 * pi / L0;
		for (int j = 0; j < Nmodes; ++j) {
			double ky = j * 2 * pi / L0;
			for (int k = 0; k < Nmodes; ++k) {
				double kz = k * 2 * pi / L0;
				//if (i + j + k > 0) {
				if (i + j + k > minModeNumber) {
					energy += sqr(evaluateAnisotropicTurbulenceAmplitude(kx, ky, kz, turbulenceKoef, index, L0, anisotropy));
				}
			}
		}
	}

	turbulenceKoef = turbulenceKoef * sqrt(B0 * B0 * fraction / energy);
}

void RadiationSourceFactory::sumFields(double*** B, double*** theta, double*** phi, double*** B1, double*** theta1, double*** phi1, double*** B2, double*** theta2, double*** phi2, int Nrho, int Nz, int Nphi) {
	for (int irho = 0; irho < Nrho; ++irho) {
		//double rho = (irho + 0.5) * R / Nrho;
		for (int iz = 0; iz < Nz; ++iz) {
			//double z = 2 * (iz + 0.5) * R / Nz - R;
			//double r = sqrt(z * z + rho * rho);
			//double cosAlpha = z / r;
			//double sinAlpha = sqrt(1.0 - cosAlpha * cosAlpha);
			for (int iphi = 0; iphi < Nphi; ++iphi) {
				//double hi = 2 * pi * (iphi + 0.5) / Nphi;

				double cosTheta1 = cos(theta1[irho][iz][iphi]);
				double cosTheta2 = cos(theta2[irho][iz][iphi]);

				double Bx = B1[irho][iz][iphi] * sin(theta1[irho][iz][iphi]) * cos(phi1[irho][iz][iphi]) +
					B2[irho][iz][iphi] * sin(theta2[irho][iz][iphi]) * cos(phi2[irho][iz][iphi]);
				double By = B1[irho][iz][iphi] * sin(theta1[irho][iz][iphi]) * sin(phi1[irho][iz][iphi]) +
					B2[irho][iz][iphi] * sin(theta2[irho][iz][iphi]) * sin(phi2[irho][iz][iphi]);
				double Bz = B1[irho][iz][iphi] * cosTheta1 + B2[irho][iz][iphi] * cosTheta2;


				double Blocal = sqrt(Bx * Bx + By * By + Bz * Bz);
				double Bxy = sqrt(Bx * Bx + By * By);

				B[irho][iz][iphi] = Blocal;
				theta[irho][iz][iphi] = acos(Bz / Blocal);
				phi[irho][iz][iphi] = atan2(By, Bx);
			}
		}
	}
}

void RadiationSourceFactory::initializeTurbulentField(double*** B, double*** theta, double*** phi, int Nrho, int Nz, int Nphi, const double& B0, const double& theta0, const double& phi0, const double& fraction, const double& index, const double& L0, const int Nmodes, const double& R)
{
	double cosTheta0 = cos(theta0);
	double sinTheta0 = sin(theta0);

	double*** Bx = create3dArray(Nrho, Nz, Nphi, B0 * sinTheta0 * cos(phi0));
	double*** By = create3dArray(Nrho, Nz, Nphi, B0 * sinTheta0 * sin(phi0));
	double*** Bz = create3dArray(Nrho, Nz, Nphi, B0 * cosTheta0);

	double*** phases1 = create3dArray(Nmodes, Nmodes, Nmodes);
	double*** phases2 = create3dArray(Nmodes, Nmodes, Nmodes);
	//srand(1234);
	for (int i = 0; i < Nmodes; ++i) {
		for (int j = 0; j < Nmodes; ++j) {
			for (int k = 0; k < Nmodes; ++k) {
				phases1[i][j][k] = 2 * pi * uniformDistribution();
				phases2[i][j][k] = 2 * pi * uniformDistribution();
			}
		}
	}

	double turbulenceKoef = 1.0;

	normalizeTurbulenceKoef(turbulenceKoef, index, L0, Nmodes, fraction, B0);


	for (int i = 0; i < Nmodes; ++i) {
		double kx = i * 2 * pi / L0;
		for (int j = 0; j < Nmodes; ++j) {
			double ky = j * 2 * pi / L0;
			for (int k = 0; k < Nmodes; ++k) {
				double kz = k * 2 * pi / L0;
				//if (i + j + k > 0) {
				if (i + j + k > minModeNumber) {
					double kt = sqrt(kx * kx + ky * ky + kz * kz);
					double cosThetat = kz / kt;
					double sinThetat = sqrt(1.0 - cosThetat * cosThetat);
					double phit = 0;
					if (i + j > 0) {
						phit = atan2(ky, kx);
					}
					double B = evaluateTurbulenceAmplitude(kt, turbulenceKoef, index, L0);
					for (int irho = 0; irho < Nrho; ++irho) {
						double rho = (irho + 0.5) * R / Nrho;
						for (int iz = 0; iz < Nz; ++iz) {
							double z = 2 * (iz + 0.5) * R / Nz - R;

							for (int iphi = 0; iphi < Nphi; ++iphi) {
								double hi = 2 * pi * (iphi + 0.5) / Nphi;

								double x = rho * cos(hi);
								double y = rho * sin(hi);

								double B1 = B * cos(kx * x + ky * y + kz * z + phases1[i][j][k]);
								double B2 = B * cos(kx * x + ky * y + kz * z + phases2[i][j][k]);

								Bx[irho][iz][iphi] += -B1 * cosThetat * cos(phit) - B2 * sin(phit);
								By[irho][iz][iphi] += -B1 * cosThetat * sin(phit) + B2 * cos(phit);
								Bz[irho][iz][iphi] += B1 * sinThetat;
							}
						}
					}
				}
			}
		}
	}

	for (int irho = 0; irho < Nrho; ++irho) {
		double rho = (irho + 0.5) * R / Nrho;
		for (int iz = 0; iz < Nz; ++iz) {
			double z = 2 * (iz + 0.5) * R / Nz - R;

			for (int iphi = 0; iphi < Nphi; ++iphi) {
				double hi = 2 * pi * (iphi + 0.5) / Nphi;

				double Blocal = sqrt(Bx[irho][iz][iphi] * Bx[irho][iz][iphi] + By[irho][iz][iphi] * By[irho][iz][iphi] + Bz[irho][iz][iphi] * Bz[irho][iz][iphi]);
				double Bxy = sqrt(Bx[irho][iz][iphi] * Bx[irho][iz][iphi] + By[irho][iz][iphi] * By[irho][iz][iphi]);

				B[irho][iz][iphi] = Blocal;
				theta[irho][iz][iphi] = acos(Bz[irho][iz][iphi] / Blocal);
				phi[irho][iz][iphi] = atan2(By[irho][iz][iphi], Bx[irho][iz][iphi]);
			}
		}
	}

	delete3dArray(Bx, Nrho, Nz, Nphi);
	delete3dArray(By, Nrho, Nz, Nphi);
	delete3dArray(Bz, Nrho, Nz, Nphi);

	delete3dArray(phases1, Nmodes, Nmodes, Nmodes);
	delete3dArray(phases2, Nmodes, Nmodes, Nmodes);
}

void RadiationSourceFactory::initializeAnisotropicLocalTurbulentFieldInSphericalSource(double*** B, double*** theta, double*** phi, int Nrho, int Nz, int Nphi, const double& B0, const double& theta0, const double& phi0, const double& fraction, const double& index, const double& L0, int Nmodes, const double& R, const double& anisotropy)
{
	double cosTheta0 = cos(theta0);
	double sinTheta0 = sin(theta0);

	double*** Bx = create3dArray(Nrho, Nz, Nphi, B0 * sinTheta0 * cos(phi0));
	double*** By = create3dArray(Nrho, Nz, Nphi, B0 * sinTheta0 * sin(phi0));
	double*** Bz = create3dArray(Nrho, Nz, Nphi, B0 * cosTheta0);

	double*** phases1 = create3dArray(Nmodes, Nmodes, Nmodes);
	double*** phases2 = create3dArray(Nmodes, Nmodes, Nmodes);
	//srand(1234);
	for (int i = 0; i < Nmodes; ++i) {
		for (int j = 0; j < Nmodes; ++j) {
			for (int k = 0; k < Nmodes; ++k) {
				phases1[i][j][k] = 2 * pi * uniformDistribution();
				phases2[i][j][k] = 2 * pi * uniformDistribution();
			}
		}
	}

	double A0 = 280;
	double A1 = 8.8E15;
	double A2 = 0.27;
	double A3 = -2E16;

	int irho;

#pragma omp parallel for private(irho) shared(Nrho, Nz, Nphi, Bx, By, Bz, A0, A1, A2, A3, phases1, phases2, index, L0, Nmodes, fraction, B0, anisotropy)
	for (int irho = 0; irho < Nrho; ++irho) {
		printf("irho = %d\n", irho);
		double rho = (irho + 0.5) * R / Nrho;
		for (int iz = 0; iz < Nz; ++iz) {
			double z = 2 * (iz + 0.5) * R / Nz - R;
			double r = sqrt(rho * rho + z * z);
			for (int iphi = 0; iphi < Nphi; ++iphi) {
				double turbulenceKoef = 1.0;
				double correction = sqrt(A0 / pow(1 + sqr((R - (r - A3)) / A1), A2));

				normalizeAnisotropicTurbulenceKoef(turbulenceKoef, index, L0, Nmodes, fraction, correction*B0, anisotropy);

				double hi = 2 * pi * (iphi + 0.5) / Nphi;

				double x = rho * cos(hi);
				double y = rho * sin(hi);

				double Bxprimed = 0;
				double Byprimed = 0;
				double Bzprimed = 0;

				double cosThetar = z / sqrt(rho*rho + z*z);
				double sinThetar = sqrt(1.0 - cosThetar * cosThetar);
				for (int i = 0; i < Nmodes; ++i) {
					double kx = i * 2 * pi / L0;
					for (int j = 0; j < Nmodes; ++j) {
						double ky = j * 2 * pi / L0;
						for (int k = 0; k < Nmodes; ++k) {
							double kz = k * 2 * pi / L0;
							//if (i + j + k > 0) {
							//phases1[i][j][k] = 2 * pi * uniformDistribution();
							//phases2[i][j][k] = 2 * pi * uniformDistribution();

							if (i + j + k > minModeNumber) {
								double kt = sqrt(kx * kx + ky * ky + kz * kz);
								double cosThetat = kz / kt;
								double sinThetat = sqrt(1.0 - cosThetat * cosThetat);
								double phit = 0;
								if (i + j > 0) {
									phit = atan2(ky, kx);
								}
								double B = evaluateAnisotropicTurbulenceAmplitude(kx, ky, kz, turbulenceKoef, index, L0, anisotropy);
								double phase = kz * r + kx*rho*hi + ky*r*acos(cosThetar);

								double B1 = B * cos(phase + phases1[i][j][k]);
								double B2 = B * cos(phase + phases2[i][j][k]);

								Bxprimed += -B1 * cosThetat * cos(phit) - B2 * sin(phit);
								Byprimed += -B1 * cosThetat * sin(phit) + B2 * cos(phit);
								Bzprimed += B1 * sinThetat;
							}
						}
					}
				}
				Bx[irho][iz][iphi] = Bxprimed*cos(hi) - sin(hi)*(Bzprimed*sinThetar - Byprimed*cosThetar);
				By[irho][iz][iphi] = (Byprimed*cosThetar + Bzprimed*sinThetar)*cos(hi) + sin(hi)*Bxprimed;
				Bz[irho][iz][iphi] = Bzprimed*cosThetar - Byprimed*sinThetar;
			}
		}
	}

	for (int irho = 0; irho < Nrho; ++irho) {
		double rho = (irho + 0.5) * R / Nrho;
		for (int iz = 0; iz < Nz; ++iz) {
			double z = 2 * (iz + 0.5) * R / Nz - R;

			for (int iphi = 0; iphi < Nphi; ++iphi) {
				double hi = 2 * pi * (iphi + 0.5) / Nphi;

				double Blocal = sqrt(Bx[irho][iz][iphi] * Bx[irho][iz][iphi] + By[irho][iz][iphi] * By[irho][iz][iphi] + Bz[irho][iz][iphi] * Bz[irho][iz][iphi]);
				double Bxy = sqrt(Bx[irho][iz][iphi] * Bx[irho][iz][iphi] + By[irho][iz][iphi] * By[irho][iz][iphi]);

				B[irho][iz][iphi] = Blocal;
				theta[irho][iz][iphi] = acos(Bz[irho][iz][iphi] / Blocal);
				phi[irho][iz][iphi] = atan2(By[irho][iz][iphi], Bx[irho][iz][iphi]);
			}
		}
	}

	delete3dArray(Bx, Nrho, Nz, Nphi);
	delete3dArray(By, Nrho, Nz, Nphi);
	delete3dArray(Bz, Nrho, Nz, Nphi);

	delete3dArray(phases1, Nmodes, Nmodes, Nmodes);
	delete3dArray(phases2, Nmodes, Nmodes, Nmodes);
}

void RadiationSourceFactory::initializeAnisotropicLocalTurbulentFieldInSectoralSphericalSource(double*** B, double*** theta, double*** phi, int Nrho, int Nz, int Nphi, const double& B0, const double& theta0, const double& phi0, const double& fraction, const double& index, const double& L0, int Nmodes, const double& R, const double& Rmin, const double& phiR, const double& anisotropy)
{
	double cosTheta0 = cos(theta0);
	double sinTheta0 = sin(theta0);

	double*** Bx = create3dArray(Nrho, Nz, Nphi, B0 * sinTheta0 * cos(phi0));
	double*** By = create3dArray(Nrho, Nz, Nphi, B0 * sinTheta0 * sin(phi0));
	double*** Bz = create3dArray(Nrho, Nz, Nphi, B0 * cosTheta0);

	double*** phases1 = create3dArray(Nmodes, Nmodes, Nmodes);
	double*** phases2 = create3dArray(Nmodes, Nmodes, Nmodes);
	//srand(1234);
	for (int i = 0; i < Nmodes; ++i) {
		for (int j = 0; j < Nmodes; ++j) {
			for (int k = 0; k < Nmodes; ++k) {
				phases1[i][j][k] = 2 * pi * uniformDistribution();
				phases2[i][j][k] = 2 * pi * uniformDistribution();
			}
		}
	}

	double A0 = 280;
	double A1 = 8.8E15;
	double A2 = 0.27;
	double A3 = -2E16;

	int irho;

	double Z = sqrt(R * R - Rmin * Rmin);

#pragma omp parallel for private(irho) shared(Nrho, Nz, Nphi, Bx, By, Bz, A0, A1, A2, A3, phases1, phases2, index, L0, Nmodes, fraction, B0, anisotropy)
	for (int irho = 0; irho < Nrho; ++irho) {
		printf("irho = %d\n", irho);
		double rho = Rmin + (irho + 0.5) * (R - Rmin) / Nrho;
		for (int iz = 0; iz < Nz; ++iz) {
			double z = 2 * (iz + 0.5) * Z / Nz - Z;
			double r = sqrt(rho * rho + z * z);
			for (int iphi = 0; iphi < Nphi; ++iphi) {
				double turbulenceKoef = 1.0;
				//double correction = sqrt(A0 / pow(1 + sqr((R - (r - A3)) / A1), A2));
				double correction = sqrt(1.0 / pow(1 + sqr((R - (r - A3)) / A1), A2));

				normalizeAnisotropicTurbulenceKoef(turbulenceKoef, index, L0, Nmodes, fraction, correction * B0, anisotropy);

				double hi = phiR * (iphi + 0.5) / Nphi;

				double x = rho * cos(hi);
				double y = rho * sin(hi);

				double Bxprimed = 0;
				double Byprimed = 0;
				double Bzprimed = 0;

				double cosThetar = z / sqrt(rho * rho + z * z);
				double sinThetar = sqrt(1.0 - cosThetar * cosThetar);
				for (int i = 0; i < Nmodes; ++i) {
					double kx = i * 2 * pi / L0;
					for (int j = 0; j < Nmodes; ++j) {
						double ky = j * 2 * pi / L0;
						for (int k = 0; k < Nmodes; ++k) {
							double kz = k * 2 * pi / L0;
							//if (i + j + k > 0) {
							//phases1[i][j][k] = 2 * pi * uniformDistribution();
							//phases2[i][j][k] = 2 * pi * uniformDistribution();

							if (i + j + k > minModeNumber) {
								double kt = sqrt(kx * kx + ky * ky + kz * kz);
								double cosThetat = kz / kt;
								double sinThetat = sqrt(1.0 - cosThetat * cosThetat);
								double phit = 0;
								if (i + j > 0) {
									phit = atan2(ky, kx);
								}
								double B = evaluateAnisotropicTurbulenceAmplitude(kx, ky, kz, turbulenceKoef, index, L0, anisotropy);
								double phase = kz * r + kx * rho * hi + ky * r * acos(cosThetar);

								double B1 = B * cos(phase + phases1[i][j][k]);
								double B2 = B * cos(phase + phases2[i][j][k]);

								Bxprimed += -B1 * cosThetat * cos(phit) - B2 * sin(phit);
								Byprimed += -B1 * cosThetat * sin(phit) + B2 * cos(phit);
								Bzprimed += B1 * sinThetat;
							}
						}
					}
				}
				Bx[irho][iz][iphi] += Bxprimed * cos(hi) - sin(hi) * (Bzprimed * sinThetar - Byprimed * cosThetar);
				By[irho][iz][iphi] += (Byprimed * cosThetar + Bzprimed * sinThetar) * cos(hi) + sin(hi) * Bxprimed;
				Bz[irho][iz][iphi] += Bzprimed * cosThetar - Byprimed * sinThetar;
			}
		}
	}

	for (int irho = 0; irho < Nrho; ++irho) {
		double rho = Rmin + (irho + 0.5) * (R-Rmin) / Nrho;
		for (int iz = 0; iz < Nz; ++iz) {
			double z = 2 * (iz + 0.5) * Z / Nz - Z;

			for (int iphi = 0; iphi < Nphi; ++iphi) {
				double hi = phiR * (iphi + 0.5) / Nphi;

				double Blocal = sqrt(Bx[irho][iz][iphi] * Bx[irho][iz][iphi] + By[irho][iz][iphi] * By[irho][iz][iphi] + Bz[irho][iz][iphi] * Bz[irho][iz][iphi]);
				double Bxy = sqrt(Bx[irho][iz][iphi] * Bx[irho][iz][iphi] + By[irho][iz][iphi] * By[irho][iz][iphi]);

				B[irho][iz][iphi] = Blocal;
				theta[irho][iz][iphi] = acos(Bz[irho][iz][iphi] / Blocal);
				phi[irho][iz][iphi] = atan2(By[irho][iz][iphi], Bx[irho][iz][iphi]);
			}
		}
	}

	delete3dArray(Bx, Nrho, Nz, Nphi);
	delete3dArray(By, Nrho, Nz, Nphi);
	delete3dArray(Bz, Nrho, Nz, Nphi);

	delete3dArray(phases1, Nmodes, Nmodes, Nmodes);
	delete3dArray(phases2, Nmodes, Nmodes, Nmodes);
}

void RadiationSourceFactory::initializeAnisotropicLocalTurbulentFieldInDiskSource(double*** B, double*** theta, double*** phi, int Nrho, int Nz, int Nphi, const double& B0, const double& theta0, const double& phi0, const double& fraction, const double& index, const double& L0, int Nmodes, const double& R, const double& anisotropy)
{
	double cosTheta0 = cos(theta0);
	double sinTheta0 = sin(theta0);

	double*** Bx = create3dArray(Nrho, Nz, Nphi, B0 * sinTheta0 * cos(phi0));
	double*** By = create3dArray(Nrho, Nz, Nphi, B0 * sinTheta0 * sin(phi0));
	double*** Bz = create3dArray(Nrho, Nz, Nphi, B0 * cosTheta0);

	double*** phases1 = create3dArray(Nmodes, Nmodes, Nmodes);
	double*** phases2 = create3dArray(Nmodes, Nmodes, Nmodes);
	//srand(1234);
	for (int i = 0; i < Nmodes; ++i) {
		for (int j = 0; j < Nmodes; ++j) {
			for (int k = 0; k < Nmodes; ++k) {
				phases1[i][j][k] = 2 * pi * uniformDistribution();
				phases2[i][j][k] = 2 * pi * uniformDistribution();
			}
		}
	}

	double turbulenceKoef = 1.0;

	normalizeAnisotropicTurbulenceKoef(turbulenceKoef, index, L0, Nmodes, fraction, B0, anisotropy);


	for (int i = 0; i < Nmodes; ++i) {
		printf("ix = %d\n", i);
		double kx = i * 2 * pi / L0;
		for (int j = 0; j < Nmodes; ++j) {
			printf("iy = %d\n", j);
			double ky = j * 2 * pi / L0;
			for (int k = 0; k < Nmodes; ++k) {
				double kz = k * 2 * pi / L0;
				//if (i + j + k > 0) {
				if (i + j + k > minModeNumber) {
					double kt = sqrt(kx * kx + ky * ky + kz * kz);
					double cosThetat = kz / kt;
					double sinThetat = sqrt(1.0 - cosThetat * cosThetat);
					double phit = 0;
					if (i + j > 0) {
						phit = atan2(ky, kx);
					}
					double B = evaluateAnisotropicTurbulenceAmplitude(kx, ky, kz, turbulenceKoef, index, L0, anisotropy);
					for (int irho = 0; irho < Nrho; ++irho) {
						double rho = (irho + 0.5) * R / Nrho;
						for (int iz = 0; iz < Nz; ++iz) {
							double z = 2 * (iz + 0.5) * R / Nz - R;

							for (int iphi = 0; iphi < Nphi; ++iphi) {
								double hi = 2 * pi * (iphi + 0.5) / Nphi;

								double x = rho * cos(hi);
								double y = rho * sin(hi);

								double B1 = B * cos(kx * x + ky * y + kz * z + phases1[i][j][k]);
								double B2 = B * cos(kx * x + ky * y + kz * z + phases2[i][j][k]);

								Bx[irho][iz][iphi] += -B1 * cosThetat * cos(phit) - B2 * sin(phit);
								By[irho][iz][iphi] += -B1 * cosThetat * sin(phit) + B2 * cos(phit);
								Bz[irho][iz][iphi] += B1 * sinThetat;
							}
						}
					}
				}
			}
		}
	}

	for (int irho = 0; irho < Nrho; ++irho) {
		double rho = (irho + 0.5) * R / Nrho;
		for (int iz = 0; iz < Nz; ++iz) {
			double z = 2 * (iz + 0.5) * R / Nz - R;

			for (int iphi = 0; iphi < Nphi; ++iphi) {
				double hi = 2 * pi * (iphi + 0.5) / Nphi;

				double Blocal = sqrt(Bx[irho][iz][iphi] * Bx[irho][iz][iphi] + By[irho][iz][iphi] * By[irho][iz][iphi] + Bz[irho][iz][iphi] * Bz[irho][iz][iphi]);
				double Bxy = sqrt(Bx[irho][iz][iphi] * Bx[irho][iz][iphi] + By[irho][iz][iphi] * By[irho][iz][iphi]);

				B[irho][iz][iphi] = Blocal;
				theta[irho][iz][iphi] = acos(Bz[irho][iz][iphi] / Blocal);
				phi[irho][iz][iphi] = atan2(By[irho][iz][iphi], Bx[irho][iz][iphi]);
			}
		}
	}

	//write3dArrayToFile(Bx, Nrho, Nz, Nphi, "Bx.dat");
	//write3dArrayToFile(By, Nrho, Nz, Nphi, "By.dat");
	//write3dArrayToFile(Bz, Nrho, Nz, Nphi, "Bz.dat");

	delete3dArray(Bx, Nrho, Nz, Nphi);
	delete3dArray(By, Nrho, Nz, Nphi);
	delete3dArray(Bz, Nrho, Nz, Nphi);

	delete3dArray(phases1, Nmodes, Nmodes, Nmodes);
	delete3dArray(phases2, Nmodes, Nmodes, Nmodes);
}

void RadiationSourceFactory::initializeParkerField(double*** B, double*** theta, double*** phi, double*** concentration, int Nrho, int Nz, int Nphi, const double& B0, const double& n0, const double& v, const double& d, const double& omega, const double& R) {
	for (int irho = 0; irho < Nrho; ++irho) {
		double rho = (irho + 0.5) * R / Nrho;
		for (int iz = 0; iz < Nz; ++iz) {
			double z = 2 * (iz + 0.5) * R / Nz - R;
			double r = sqrt(z * z + rho * rho);
			double cosAlpha = z / r;
			double sinAlpha = sqrt(1.0 - cosAlpha * cosAlpha);
			for (int iphi = 0; iphi < Nphi; ++iphi) {
				double hi = 2 * pi * (iphi + 0.5) / Nphi;

				double Br = B0 * cosAlpha * sqr(d / r);
				double Bphi = B0 * cosAlpha * sinAlpha * (omega / v) * (r - d) * sqr(d / r);

				double Bx = Br * sinAlpha * cos(hi) - Bphi * sin(hi);
				double By = Br * sinAlpha * sin(hi) + Bphi * cos(hi);

				B[irho][iz][iphi] = sqrt(Br * Br + Bphi * Bphi);
				theta[irho][iz][iphi] = acos(Br * cosAlpha / B[irho][iz][iphi]);
				phi[irho][iz][iphi] = atan2(By, Bx);

				concentration[irho][iz][iphi] = n0 * sqr(R / r);
			}
		}
	}
}

void RadiationSourceFactory::initializeParkerFieldWithRotation(double*** B, double*** theta, double*** phi, double*** concentration, int Nrho, int Nz, int Nphi, const double& B0, const double& n0, const double& v, const double& d, const double& omega, const double& R, const double& thetaRot) {
	double cosThetaRot = cos(thetaRot);
	double sinThetaRot = sin(thetaRot);
	for (int irho = 0; irho < Nrho; ++irho) {
		double rho = (irho + 0.5) * R / Nrho;
		for (int iz = 0; iz < Nz; ++iz) {
			double z = 2 * (iz + 0.5) * R / Nz - R;
			double r = sqrt(z * z + rho * rho);
			for (int iphi = 0; iphi < Nphi; ++iphi) {
				double hi = 2 * pi * (iphi + 0.5) / Nphi;

				double x = rho * cos(hi);
				double y1 = rho * sin(hi);
				double z1 = z * cosThetaRot - x * sinThetaRot;
				double x1 = x * cosThetaRot + z * sinThetaRot;
				double r1 = r;
				double cosAlpha1 = z1 / r;
				double sinAlpha1 = sqrt(1.0 - cosAlpha1 * cosAlpha1);
				double hi1 = atan2(y1, x1);


				double Br1 = B0 * cosAlpha1 * sqr(d / r1);
				double Bphi1 = B0 * cosAlpha1 * sinAlpha1 * (omega / v) * (r1 - d) * sqr(d / r1);

				double Bx1 = Br1 * sinAlpha1 * cos(hi1) - Bphi1 * sin(hi1);
				double By1 = Br1 * sinAlpha1 * sin(hi1) + Bphi1 * cos(hi1);
				double Bz1 = Br1 * cosAlpha1;

				B[irho][iz][iphi] = sqrt(Br1 * Br1 + Bphi1 * Bphi1);

				double Bx = Bx1 * cosThetaRot + Bz1 * sinThetaRot;
				double By = By1;
				double Bz = Bz1 * cosThetaRot - Bx1 * sinThetaRot;
				theta[irho][iz][iphi] = acos(Bz / B[irho][iz][iphi]);
				phi[irho][iz][iphi] = atan2(By, Bx);

				concentration[irho][iz][iphi] = n0 * sqr(R / r);
			}
		}
	}
}

void RadiationSourceFactory::initializeAngularMask(bool** mask, int Nrho, int Nphi, const double& angle) {
	if (angle < 0) {
		printf("angle in mask < 0\n");
		printLog("angle in mask < 0\n");
		exit(0);
	}
	if (angle > pi / 2) {
		printf("angle in mask greater than pi/2\n");
		printLog("angle in mask greater than pi/2\n");
		exit(0);
	}
	for (int iphi = 0; iphi < Nphi; ++iphi) {
		double phi = 2 * pi * (iphi + 0.5) / Nphi;
		bool value = true;
		if (((phi > angle) && (phi < pi - angle)) || ((phi > pi + angle) && (phi < 2 * pi - angle))) {
			value = false;
		}
		for (int irho = 0; irho < Nrho; ++irho) {
			mask[irho][iphi] = value;
		}
	}
}

void RadiationSourceFactory::initialize3dAngularMask(bool** mask, int Nrho, int Nphi, const double& angle) {
	if (angle < 0) {
		printf("angle in mask < 0\n");
		printLog("angle in mask < 0\n");
		exit(0);
	}
	if (angle > pi / 2) {
		printf("angle in mask greater than pi/2\n");
		printLog("angle in mask greater than pi/2\n");
		exit(0);
	}
	for (int iphi = 0; iphi < Nphi; ++iphi) {
		double phi = 2 * pi * (iphi + 0.5) / Nphi;
		bool value = true;
		if (((phi > angle) && (phi < pi - angle)) || ((phi > pi + angle) && (phi < 2 * pi - angle))) {
			value = false;
		}
		for (int irho = 0; irho < Nrho; ++irho) {
			double rho = (irho + 0.5) / Nrho;
			if (rho < cos(angle)) {
				mask[irho][iphi] = false;
			}
			else {
				mask[irho][iphi] = value;
			}
		}
	}
}

void RadiationSourceFactory::initializeRhoMask(bool** mask, int Nrho, int Nphi, const double& fraction) {
	for (int irho = 0; irho < Nrho; ++irho) {
		bool value = ((irho + 0.5) / Nrho > fraction);
		for (int iphi = 0; iphi < Nphi; ++iphi) {
			mask[irho][iphi] = value;
		}
	}
}

/*AngleDependentElectronsSphericalSource* RadiationSourceFactory::createAngularDependentSource(MassiveParticleDistribution** electronDistributions, int Ntheta, int Nrho, int Nz, int Nphi, const double& B0, const double& theta0, const double& phi0, const double n0, const double& rho, const double& rhoin, const double& distance)
{
	double*** B = create3dArray(Nrho, Nz, Nphi, B0);
	double*** theta = create3dArray(Nrho, Nz, Nphi, theta0);
	double*** phi = create3dArray(Nrho, Nz, Nphi, phi0);

	double*** concentration = create3dArray(Nrho, Nz, Nphi, n0);

	AngleDependentElectronsSphericalSource* source = new AngleDependentElectronsSphericalSource(Nrho, Nz, Nphi, Ntheta, electronDistributions, B, theta, phi, concentration, rho, rhoin, distance);

	delete3dArray(B, Nrho, Nz, Nphi);
	delete3dArray(theta, Nrho, Nz, Nphi);
	delete3dArray(phi, Nrho, Nz, Nphi);
	delete3dArray(concentration, Nrho, Nz, Nphi);
	return source;
}*/

AngleDependentElectronsSphericalSource* RadiationSourceFactory::createSourceWithTurbulentField(MassiveParticleDistribution** electronDistributions, int Ntheta, int Nrho, int Nz, int Nphi, const double& B0, const double& theta0, const double& phi0, const double n0, const double& turbulentFraction, const double& index, const double& L0, int Nmodes, const double& rho, const double& rhoin, const double& distance)
{
	double*** B = create3dArray(Nrho, Nz, Nphi, B0);
	double*** theta = create3dArray(Nrho, Nz, Nphi, theta0);
	double*** phi = create3dArray(Nrho, Nz, Nphi, phi0);

	double*** concentration = create3dArray(Nrho, Nz, Nphi, n0);

	initializeTurbulentField(B, theta, phi, Nrho, Nz, Nphi, B0, theta0, phi0, turbulentFraction, index, L0, Nmodes, rho);

	AngleDependentElectronsSphericalSource* source = new AngleDependentElectronsSphericalSource(Nrho, Nz, Nphi, Ntheta, electronDistributions, B, theta, phi, concentration, rho, rhoin, distance);

	delete3dArray(B, Nrho, Nz, Nphi);
	delete3dArray(theta, Nrho, Nz, Nphi);
	delete3dArray(phi, Nrho, Nz, Nphi);
	delete3dArray(concentration, Nrho, Nz, Nphi);
	return source;
}

AngleDependentElectronsSphericalSource* RadiationSourceFactory::createSourceWithParkerField(MassiveParticleDistribution** electronDistributions, int Ntheta, int Nrho, int Nz, int Nphi, const double& B0, const double& n0, const double& v, const double& d, const double& omega, const double& rho, const double& rhoin, const double& distance)
{
	double*** B = new double** [Nrho];
	double*** theta = new double** [Nrho];
	double*** phi = new double** [Nrho];
	double*** concentration = new double** [Nrho];
	for (int irho = 0; irho < Nrho; ++irho) {
		B[irho] = new double* [Nz];
		theta[irho] = new double* [Nz];
		phi[irho] = new double* [Nz];
		concentration[irho] = new double* [Nz];
		for (int iz = 0; iz < Nz; ++iz) {
			B[irho][iz] = new double[Nphi];
			theta[irho][iz] = new double[Nphi];
			phi[irho][iz] = new double[Nphi];
			concentration[irho][iz] = new double[Nphi];
		}
	}

	initializeParkerField(B, theta, phi, concentration, Nrho, Nz, Nphi, B0, n0, v, d, omega, rho);

	for (int irho = 0; irho < Nrho; ++irho) {
		for (int iz = 0; iz < Nz; ++iz) {
			delete[] B[irho][iz];
			delete[] theta[irho][iz];
			delete[] phi[irho][iz];
			delete[] concentration[irho][iz];
		}
		delete[] B[irho];
		delete[] theta[irho];
		delete[] phi[irho];
		delete[] concentration[irho];
	}
	delete[] B;
	delete[] theta;
	delete[] phi;
	delete[] concentration;

	return new AngleDependentElectronsSphericalSource(Nrho, Nz, Nphi, Ntheta, electronDistributions, B, theta, phi, concentration, rho, rhoin, distance);
}

AngleDependentElectronsSphericalSource* RadiationSourceFactory::createSourceWithParkerFieldWithRotation(MassiveParticleDistribution** electronDistributions, int Ntheta, int Nrho, int Nz, int Nphi, const double& B0, const double& n0, const double& v, const double& d, const double& omega, const double& thetaRot, const double& rho, const double& rhoin, const double& distance)
{
	double*** B = new double** [Nrho];
	double*** theta = new double** [Nrho];
	double*** phi = new double** [Nrho];
	double*** concentration = new double** [Nrho];
	for (int irho = 0; irho < Nrho; ++irho) {
		B[irho] = new double* [Nz];
		theta[irho] = new double* [Nz];
		phi[irho] = new double* [Nz];
		concentration[irho] = new double* [Nz];
		for (int iz = 0; iz < Nz; ++iz) {
			B[irho][iz] = new double[Nphi];
			theta[irho][iz] = new double[Nphi];
			phi[irho][iz] = new double[Nphi];
			concentration[irho][iz] = new double[Nphi];
		}
	}

	initializeParkerFieldWithRotation(B, theta, phi, concentration, Nrho, Nz, Nphi, B0, n0, v, d, omega, rho, thetaRot);

	for (int irho = 0; irho < Nrho; ++irho) {
		for (int iz = 0; iz < Nz; ++iz) {
			delete[] B[irho][iz];
			delete[] theta[irho][iz];
			delete[] phi[irho][iz];
			delete[] concentration[irho][iz];
		}
		delete[] B[irho];
		delete[] theta[irho];
		delete[] phi[irho];
		delete[] concentration[irho];
	}
	delete[] B;
	delete[] theta;
	delete[] phi;
	delete[] concentration;

	return new AngleDependentElectronsSphericalSource(Nrho, Nz, Nphi, Ntheta, electronDistributions, B, theta, phi, concentration, rho, rhoin, distance);
}

TabulatedDiskSource* RadiationSourceFactory::readDiskSourceFromFile(MassiveParticleDistribution* electronDistribution, const double& rho, const double& zmin, const double zmax, const int Nrho, const int Nz, const int Nphi, const double& distance, SourceInputGeometry geometry, const char* BFileName, const char* concentrationFileName, const double& thetar, const double& phir, const double& psir)
{
	FILE* Bfile = fopen(BFileName, "r");
	FILE* concentrationFile = fopen(concentrationFileName, "r");
	int Nxb, Nyb, Nzb;
	int Nxc, Nyc, Nzc;

	fscanf(Bfile, "%d %d %d", &Nxb, &Nyb, &Nzb);
	fscanf(concentrationFile, "%d %d %d", &Nxc, &Nyc, &Nzc);

	double x1, y1, z1;
	double x2, y2, z2;

	fscanf(Bfile, "%lf %lf %lf", &x1, &y1, &z1);
	fscanf(concentrationFile, "%lf %lf %lf", &x1, &y1, &z1);

	fscanf(Bfile, "%lf %lf %lf", &x2, &y2, &z2);
	fscanf(concentrationFile, "%lf %lf %lf", &x2, &y2, &z2);

	if ((Nxb != Nxc) || (Nyb != Nyc) || (Nzb != Nzc)) {
		printf("Dimensions are different in B and concentrationFile\n");
		printLog("Dimensions are different in B and concentrationFile\n");
		exit(0);
	}

	double*** concentration2 = new double** [Nrho];
	double*** B = new double** [Nrho];
	double*** Btheta = new double** [Nrho];
	double*** Bphi = new double** [Nrho];

	for (int i = 0; i < Nrho; ++i) {
		concentration2[i] = new double* [Nz];
		B[i] = new double* [Nz];
		Btheta[i] = new double* [Nz];
		Bphi[i] = new double* [Nz];
		for (int j = 0; j < Nz; ++j) {
			concentration2[i][j] = new double[Nphi];
			B[i][j] = new double[Nphi];
			Btheta[i][j] = new double[Nphi];
			Bphi[i][j] = new double[Nphi];
		}
	}

	double*** concentration = new double** [Nxb];

	double*** B1 = new double** [Nxb];
	double*** B2 = new double** [Nxb];
	double*** B3 = new double** [Nxb];
	for (int i = 0; i < Nxb; ++i) {
		concentration[i] = new double* [Nyb];
		B1[i] = new double* [Nyb];
		B2[i] = new double* [Nyb];
		B3[i] = new double* [Nyb];
		for (int j = 0; j < Nyb; ++j) {
			concentration[i][j] = new double[Nzb];
			B1[i][j] = new double[Nzb];
			B2[i][j] = new double[Nzb];
			B3[i][j] = new double[Nzb];
			for (int k = 0; k < Nzb; ++k) {
				fscanf(concentrationFile, "%lf", &concentration[i][j][k]);
				fscanf(Bfile, "%lf %lf %lf", &B1[i][j][k], &B2[i][j][k], &B3[i][j][k]);		
				if (checkNaNorInfinity(concentration[i][j][k])) {
					printf("concentration is Nan or Infinity %lf\n", concentration[i][j][k]);
					printLog("concentration is Nan or Infinity %lf\n", concentration[i][j][k]);
					exit(0);
				}
				if (checkNaNorInfinity(B1[i][j][k])) {
					printf("B1 is Nan or Infinity %lf\n", B1[i][j][k]);
					printLog("B1 is Nan or Infinity %lf\n", B1[i][j][k]);
					exit(0);
				}

				if (checkNaNorInfinity(B2[i][j][k])) {
					printf("B2 is Nan or Infinity %lf\n", B2[i][j][k]);
					printLog("B2 is Nan or Infinity %lf\n", B2[i][j][k]);
					exit(0);
				}

				if (checkNaNorInfinity(B3[i][j][k])) {
					printf("B3 is Nan or Infinity %lf\n", B3[i][j][k]);
					printLog("B3 is Nan or Infinity %lf\n", B3[i][j][k]);
					exit(0);
				}
			}
		}
	}

	correctSourceConcentration(concentration, Nxb, Nyb, Nzb, Nxb-1, Nyb-1, Nzb-1, 0.002);
	correctSourceConcentration(concentration, Nxb, Nyb, Nzb, Nxb/2, Nyb/2, Nzb - 1, 0.002);

	transformScalarArrayToCylindrical(concentration, Nxb, Nyb, Nzb, x1, x2, y1, y2, z1, z2, geometry, concentration2, Nrho, Nz, Nphi, 0, rho, zmin, zmax, 0, 2 * pi, thetar, phir, psir);
	transformVectorArraysToCylindrical(B1, B2, B3, Nxb, Nyb, Nzb, x1, x2, y1, y2, z1, z2, geometry, B, Btheta, Bphi, Nrho, Nz, Nphi, 0, rho, zmin, zmax, 0, 2 * pi, thetar, phir, psir);

	TabulatedDiskSource* source = new TabulatedDiskSource(Nrho, Nz, Nphi, electronDistribution, B, Btheta, Bphi, concentration2, rho, zmax - zmin, distance);

	fclose(Bfile);
	fclose(concentrationFile);
	for (int i = 0; i < Nxb; ++i) {
		for (int j = 0; j < Nyb; ++j) {
			delete[] concentration[i][j];
			delete[] B1[i][j];
			delete[] B2[i][j];
			delete[] B3[i][j];
		}
		delete[] concentration[i];
		delete[] B1[i];
		delete[] B2[i];
		delete[] B3[i];
	}
	delete[] concentration;
	delete[] B1;
	delete[] B2;
	delete[] B3;

	for (int i = 0; i < Nrho; ++i) {
		for (int j = 0; j < Nz; ++j) {
			delete[] concentration2[i][j];
			delete[] B[i][j];
			delete[] Btheta[i][j];
			delete[] Bphi[i][j];
		}
		delete[] concentration2[i];
		delete[] B[i];
		delete[] Btheta[i];
		delete[] Bphi[i];
	}
	delete[] concentration2;
	delete[] B;
	delete[] Btheta;
	delete[] Bphi;

	return source;
}

void RadiationSourceFactory::readRectangularSourceArraysFromFile(double***& B, double***& Btheta, double***& Bphi, double***& concentration, const double& minX, const double& maxX, const double& minZ, const double& maxZ, const double& minY, const double& maxY, const int Nx, const int Nz, const int Ny, SourceInputGeometry geometry, const char* BFileName, const char* concentrationFileName, const double& thetar, const double& phir, const double& psir)
{
	FILE* Bfile = fopen(BFileName, "r");
	FILE* concentrationFile = fopen(concentrationFileName, "r");
	int Nxb, Nyb, Nzb;
	int Nxc, Nyc, Nzc;

	fscanf(Bfile, "%d %d %d", &Nxb, &Nyb, &Nzb);
	fscanf(concentrationFile, "%d %d %d", &Nxc, &Nyc, &Nzc);

	double x1, y1, z1;
	double x2, y2, z2;

	fscanf(Bfile, "%lf %lf %lf", &x1, &y1, &z1);
	fscanf(concentrationFile, "%lf %lf %lf", &x1, &y1, &z1);

	fscanf(Bfile, "%lf %lf %lf", &x2, &y2, &z2);
	fscanf(concentrationFile, "%lf %lf %lf", &x2, &y2, &z2);

	if ((Nxb != Nxc) || (Nyb != Nyc) || (Nzb != Nzc)) {
		printf("Dimensions are different in B and concentrationFile\n");
		printLog("Dimensions are different in B and concentrationFile\n");
		exit(0);
	}

	concentration = new double** [Nx];
	B = new double** [Nx];
	Btheta = new double** [Nx];
	Bphi = new double** [Nx];

	for (int i = 0; i < Nx; ++i) {
		concentration[i] = new double* [Nz];
		B[i] = new double* [Nz];
		Btheta[i] = new double* [Nz];
		Bphi[i] = new double* [Nz];
		for (int j = 0; j < Nz; ++j) {
			concentration[i][j] = new double[Ny];
			B[i][j] = new double[Ny];
			Btheta[i][j] = new double[Ny];
			Bphi[i][j] = new double[Ny];
		}
	}

	double*** concentration1 = new double** [Nxb];

	double*** B1 = new double** [Nxb];
	double*** B2 = new double** [Nxb];
	double*** B3 = new double** [Nxb];
	for (int i = 0; i < Nxb; ++i) {
		concentration1[i] = new double* [Nyb];
		B1[i] = new double* [Nyb];
		B2[i] = new double* [Nyb];
		B3[i] = new double* [Nyb];
		for (int j = 0; j < Nyb; ++j) {
			concentration1[i][j] = new double[Nzb];
			B1[i][j] = new double[Nzb];
			B2[i][j] = new double[Nzb];
			B3[i][j] = new double[Nzb];
			for (int k = 0; k < Nzb; ++k) {
				fscanf(concentrationFile, "%lf", &concentration1[i][j][k]);
				fscanf(Bfile, "%lf %lf %lf", &B1[i][j][k], &B2[i][j][k], &B3[i][j][k]);
				if (checkNaNorInfinity(concentration1[i][j][k])) {
					printf("concentration is Nan or Infinity %lf\n", concentration1[i][j][k]);
					printLog("concentration is Nan or Infinity %lf\n", concentration1[i][j][k]);
					exit(0);
				}
				if (checkNaNorInfinity(B1[i][j][k])) {
					printf("B1 is Nan or Infinity %lf\n", B1[i][j][k]);
					printLog("B1 is Nan or Infinity %lf\n", B1[i][j][k]);
					exit(0);
				}

				if (checkNaNorInfinity(B2[i][j][k])) {
					printf("B2 is Nan or Infinity %lf\n", B2[i][j][k]);
					printLog("B2 is Nan or Infinity %lf\n", B2[i][j][k]);
					exit(0);
				}

				if (checkNaNorInfinity(B3[i][j][k])) {
					printf("B3 is Nan or Infinity %lf\n", B3[i][j][k]);
					printLog("B3 is Nan or Infinity %lf\n", B3[i][j][k]);
					exit(0);
				}
			}
		}
	}

	//correctSourceConcentration(concentration, Nxb, Nyb, Nzb, Nxb - 1, Nyb - 1, Nzb - 1, 0.002);
	//correctSourceConcentration(concentration, Nxb, Nyb, Nzb, Nxb / 2, Nyb / 2, Nzb - 1, 0.002);

	transformScalarArrayToCartesian(concentration1, Nxb, Nyb, Nzb, x1, x2, y1, y2, z1, z2, geometry, concentration, Nx, Nz, Ny, minX, maxX, minZ, maxZ, minY, maxY, thetar, phir, psir);
	transformVectorArraysToCartesian(B1, B2, B3, Nxb, Nyb, Nzb, x1, x2, y1, y2, z1, z2, geometry, B, Btheta, Bphi, Nx, Nz, Ny, minX, maxX, minZ, maxZ, minY, maxY, thetar, phir, psir);

	fclose(Bfile);
	fclose(concentrationFile);
	for (int i = 0; i < Nxb; ++i) {
		for (int j = 0; j < Nyb; ++j) {
			delete[] concentration1[i][j];
			delete[] B1[i][j];
			delete[] B2[i][j];
			delete[] B3[i][j];
		}
		delete[] concentration1[i];
		delete[] B1[i];
		delete[] B2[i];
		delete[] B3[i];
	}
	delete[] concentration1;
	delete[] B1;
	delete[] B2;
	delete[] B3;
}

RectangularSource* RadiationSourceFactory::readRectangularSourceFromFile(MassiveParticleDistribution* electronDistribution, const double& minX, const double& maxX, const double& minZ, const double& maxZ, const double& minY, const double& maxY, const int Nx, const int Nz, const int Ny, const double& distance, SourceInputGeometry geometry, const char* BFileName, const char* concentrationFileName, const double& thetar, const double& phir, const double& psir)
{

	double*** concentration;
	double*** B;
	double*** Btheta;
	double*** Bphi;

	RadiationSourceFactory::readRectangularSourceArraysFromFile(B, Btheta, Bphi, concentration, minX, maxX, minZ, maxZ, minY, maxY, Nx, Nz, Ny, geometry, BFileName, concentrationFileName, thetar, phir, psir);

	RectangularSource* source = new RectangularSource(Nx, Ny, Nz, electronDistribution, B, Btheta, Bphi, concentration, minX, maxX, minY, maxY, minZ, maxZ, distance);


	for (int i = 0; i < Nx; ++i) {
		for (int j = 0; j < Nz; ++j) {
			delete[] concentration[i][j];
			delete[] B[i][j];
			delete[] Btheta[i][j];
			delete[] Bphi[i][j];
		}
		delete[] concentration[i];
		delete[] B[i];
		delete[] Btheta[i];
		delete[] Bphi[i];
	}
	delete[] concentration;
	delete[] B;
	delete[] Btheta;
	delete[] Bphi;

	return source;
}

ThermalRectangularSource* RadiationSourceFactory::readThermalRectangularSourceFromFile(const double& minX, const double& maxX, const double& minZ, const double& maxZ, const double& minY, const double& maxY, const int Nx, const int Nz, const int Ny, const double& distance, SourceInputGeometry geometry, const char* BFileName, const char* concentrationFileName, const char* temperatureFileName, const double& thetar, const double& phir, const double& psir)
{
	FILE* Bfile = fopen(BFileName, "r");
	FILE* concentrationFile = fopen(concentrationFileName, "r");
	FILE* temperatureFile = fopen(temperatureFileName, "r");
	int Nxb, Nyb, Nzb;
	int Nxc, Nyc, Nzc;
	int Nxt, Nyt, Nzt;

	fscanf(Bfile, "%d %d %d", &Nxb, &Nyb, &Nzb);
	fscanf(concentrationFile, "%d %d %d", &Nxc, &Nyc, &Nzc);
	fscanf(temperatureFile, "%d %d %d", &Nxt, &Nyt, &Nzt);

	double x1, y1, z1;
	double x2, y2, z2;

	fscanf(Bfile, "%lf %lf %lf", &x1, &y1, &z1);
	fscanf(concentrationFile, "%lf %lf %lf", &x1, &y1, &z1);
	fscanf(temperatureFile, "%lf %lf %lf", &x1, &y1, &z1);

	fscanf(Bfile, "%lf %lf %lf", &x2, &y2, &z2);
	fscanf(concentrationFile, "%lf %lf %lf", &x2, &y2, &z2);
	fscanf(temperatureFile, "%lf %lf %lf", &x2, &y2, &z2);

	if ((Nxb != Nxc) || (Nyb != Nyc) || (Nzb != Nzc) || (Nxb != Nxt) || (Nyb != Nyt) || (Nzb != Nzt)) {
		printf("Dimensions are different in B and concentrationFile\n");
		printLog("Dimensions are different in B and concentrationFile\n");
		exit(0);
	}

	double*** concentration2 = new double** [Nx];
	double*** B = new double** [Nx];
	double*** Btheta = new double** [Nx];
	double*** Bphi = new double** [Nx];
	double*** temperature2 = new double** [Nx];

	for (int i = 0; i < Nx; ++i) {
		concentration2[i] = new double* [Nz];
		B[i] = new double* [Nz];
		Btheta[i] = new double* [Nz];
		Bphi[i] = new double* [Nz];
		temperature2[i] = new double* [Nz];
		for (int j = 0; j < Nz; ++j) {
			concentration2[i][j] = new double[Ny];
			B[i][j] = new double[Ny];
			Btheta[i][j] = new double[Ny];
			Bphi[i][j] = new double[Ny];
			temperature2[i][j] = new double[Ny];
		}
	}

	double*** concentration = new double** [Nxb];

	double*** B1 = new double** [Nxb];
	double*** B2 = new double** [Nxb];
	double*** B3 = new double** [Nxb];
	double*** temperature = new double** [Nxb];
	for (int i = 0; i < Nxb; ++i) {
		concentration[i] = new double* [Nyb];
		B1[i] = new double* [Nyb];
		B2[i] = new double* [Nyb];
		B3[i] = new double* [Nyb];
		temperature[i] = new double* [Nyb];
		for (int j = 0; j < Nyb; ++j) {
			concentration[i][j] = new double[Nzb];
			B1[i][j] = new double[Nzb];
			B2[i][j] = new double[Nzb];
			B3[i][j] = new double[Nzb];
			temperature[i][j] = new double[Nzb];
			for (int k = 0; k < Nzb; ++k) {
				fscanf(concentrationFile, "%lf", &concentration[i][j][k]);
				fscanf(Bfile, "%lf %lf %lf", &B1[i][j][k], &B2[i][j][k], &B3[i][j][k]);
				fscanf(temperatureFile, "%lf", &temperature[i][j][k]);
				if (checkNaNorInfinity(concentration[i][j][k])) {
					printf("concentration is Nan or Infinity %lf\n", concentration[i][j][k]);
					printLog("concentration is Nan or Infinity %lf\n", concentration[i][j][k]);
					exit(0);
				}
				if (checkNaNorInfinity(B1[i][j][k])) {
					printf("B1 is Nan or Infinity %lf\n", B1[i][j][k]);
					printLog("B1 is Nan or Infinity %lf\n", B1[i][j][k]);
					exit(0);
				}

				if (checkNaNorInfinity(B2[i][j][k])) {
					printf("B2 is Nan or Infinity %lf\n", B2[i][j][k]);
					printLog("B2 is Nan or Infinity %lf\n", B2[i][j][k]);
					exit(0);
				}

				if (checkNaNorInfinity(B3[i][j][k])) {
					printf("B3 is Nan or Infinity %lf\n", B3[i][j][k]);
					printLog("B3 is Nan or Infinity %lf\n", B3[i][j][k]);
					exit(0);
				}
				if (checkNaNorInfinity(temperature[i][j][k])) {
					printf("temperature is Nan or Infinity %lf\n", temperature[i][j][k]);
					printLog("temperature is Nan or Infinity %lf\n", temperature[i][j][k]);
					exit(0);
				}
			}
		}
	}

	//correctSourceConcentration(concentration, Nxb, Nyb, Nzb, Nxb - 1, Nyb - 1, Nzb - 1, 0.002);
	//correctSourceConcentration(concentration, Nxb, Nyb, Nzb, Nxb / 2, Nyb / 2, Nzb - 1, 0.002);

	transformScalarArrayToCartesian(concentration, Nxb, Nyb, Nzb, x1, x2, y1, y2, z1, z2, geometry, concentration2, Nx, Nz, Ny, minX, maxX, minZ, maxZ, minY, maxY, thetar, phir, psir);
	transformVectorArraysToCartesian(B1, B2, B3, Nxb, Nyb, Nzb, x1, x2, y1, y2, z1, z2, geometry, B, Btheta, Bphi, Nx, Nz, Ny, minX, maxX, minZ, maxZ, minY, maxY, thetar, phir, psir);
	transformScalarArrayToCartesian(temperature, Nxb, Nyb, Nzb, x1, x2, y1, y2, z1, z2, geometry, temperature2, Nx, Nz, Ny, minX, maxX, minZ, maxZ, minY, maxY, thetar, phir, psir);

	correctSourceTemperature(temperature2, Nx, Ny, Nz, 1.0);

	ThermalRectangularSource* source = new ThermalRectangularSource(Nx, Ny, Nz, massElectron, B, Btheta, Bphi, concentration2, temperature2, minX, maxX, minY, maxY, minZ, maxZ, distance);

	fclose(Bfile);
	fclose(concentrationFile);
	fclose(temperatureFile);
	for (int i = 0; i < Nxb; ++i) {
		for (int j = 0; j < Nyb; ++j) {
			delete[] concentration[i][j];
			delete[] B1[i][j];
			delete[] B2[i][j];
			delete[] B3[i][j];
			delete[] temperature[i][j];
		}
		delete[] concentration[i];
		delete[] B1[i];
		delete[] B2[i];
		delete[] B3[i];
		delete[] temperature[i];
	}
	delete[] concentration;
	delete[] B1;
	delete[] B2;
	delete[] B3;
	delete[] temperature;

	for (int i = 0; i < Nx; ++i) {
		for (int j = 0; j < Nz; ++j) {
			delete[] concentration2[i][j];
			delete[] B[i][j];
			delete[] Btheta[i][j];
			delete[] Bphi[i][j];
			delete[] temperature2[i][j];
		}
		delete[] concentration2[i];
		delete[] B[i];
		delete[] Btheta[i];
		delete[] Bphi[i];
		delete[] temperature2[i];
	}
	delete[] concentration2;
	delete[] B;
	delete[] Btheta;
	delete[] Bphi;
	delete[] temperature2;

	return source;
}

void RadiationSourceFactory::correctSourceConcentration(double*** concentration, const int Nx, const int Ny, const int Nz, const int start1, const int start2, const int start3, const double& threshold)
{
	bool*** sourceMark = breadthFirstSearchOfRegion(concentration, Nx, Ny, Nz, start1, start2, start3, threshold);

	FILE* markFile = fopen("marker.dat", "w");
	for (int i = 0; i < Nx; ++i) {
		for (int j = 0; j < Ny; ++j) {
			for (int k = 0; k < Nz; ++k) {
				if (sourceMark[i][j][k]) {
					concentration[i][j][k] = 0;
					fprintf(markFile, "%lf\n", 0.0);
				}
				else {
					fprintf(markFile, "%lf\n", 1.0);
				}
			}
		}
	}

	fclose(markFile);

	for (int i = 0; i < Nx; ++i) {
		for (int j = 0; j < Ny; ++j) {
			delete[] sourceMark[i][j];
		}
		delete[] sourceMark[i];
	}
	delete[] sourceMark;
}

void RadiationSourceFactory::correctSourceTemperature(double*** temperature, const int Nx, const int Ny, const int Nz, const double& value)
{
	for (int i = 0; i < Nx; ++i) {
		for (int j = 0; j < Nz; ++j) {
			for (int k = 0; k < Ny; ++k) {
				if (temperature[i][j][k] <= 0.0) {
					temperature[i][j][k] = value;
				}
			}
		}
	}
}

//transforms array from arbitrary pluto coordinates to rotated cartesian coordinates
void RadiationSourceFactory::transformScalarArrayToCartesian(double*** inputArray, int N1, int N2, int N3, const double& xmin1, const double& xmax1, const double& xmin2, const double& xmax2, const double& xmin3, const double& xmax3, SourceInputGeometry geometry, double*** outputArray, int Nx, int Nz, int Ny, const double& minX, const double& maxX, const double& minZ, const double& maxZ, const double& minY, const double& maxY, const double& thetar, const double& phir, const double& psir)
{
	if (maxX <= minX) {
		printf("rhomax = %g <= rhomin = %g in transformScalarArray\n", maxX, minX);
		printLog("rhomax = %g <= rhomin = %g in transformScalarArray\n", maxX, minX);
		exit(0);
	}
	if (maxZ <= minZ) {
		printf("zmax = %g <= zmin = %g in transformScalarArray\n", maxZ, minZ);
		printLog("zmax = %g <= zmin = %g in transformScalarArray\n", maxZ, minZ);
		exit(0);
	}
	if (maxY <= minY) {
		printf("phimax = %g <= phimin = %g in transformScalarArray\n", maxY, minY);
		printLog("phimax = %g <= phimin = %g in transformScalarArray\n", maxY, minY);
		exit(0);
	}
	if (xmax1 <= xmin1) {
		printf("xmax1 = %g <= xmin1 = %g in transformScalarArray\n", xmax1, xmin1);
		printLog("xmax1 = %g <= xmin1 = %g in transformScalarArray\n", xmax1, xmin1);
		exit(0);
	}
	if (xmax2 <= xmin2) {
		printf("xmax2 = %g <= xmin2 = %g in transformScalarArray\n", xmax2, xmin2);
		printLog("xmax2 = %g <= xmin2 = %g in transformScalarArray\n", xmax2, xmin2);
		exit(0);
	}
	if (xmax3 <= xmin3) {
		printf("xmax3 = %g <= xmin3 = %g in transformScalarArray\n", xmax3, xmin3);
		printLog("xmax3 = %g <= xmin3 = %g in transformScalarArray\n", xmax3, xmin3);
		exit(0);
	}

	double dx = (maxX - minX) / Nx;
	double dz = (maxZ - minZ) / Nz;
	double dy = (maxY - minY) / Ny;

	double dx1 = (xmax1 - xmin1) / N1;
	double dx2 = (xmax2 - xmin2) / N2;
	double dx3 = (xmax3 - xmin3) / N3;

	for (int i = 0; i < Nx; ++i) {
		double x = minX + dx * (i + 0.5);
		for (int j = 0; j < Nz; ++j) {
			double z = minZ + dz * (j + 0.5);
			for (int k = 0; k < Ny; ++k) {
				double y = minY + dy * (k + 0.5);

				double z1 = z;
				double x1 = x * cos(psir) - y*sin(psir);
				double y1 = y * cos(psir) + x*sin(psir);

				double z2 = z1 * cos(thetar) - x1 * sin(thetar);
				double x2 = (x1 * cos(thetar) + z1 * sin(thetar)) * cos(phir) - y1 * sin(phir);
				double y2 = y1 * cos(phir) + (x1 * cos(thetar) + z1 * sin(thetar)) * sin(phir);

				if (geometry == SourceInputGeometry::CARTESIAN) {
					if ((N1 > 1) && ((x2 < xmin1) || (x2 > xmax1))) {
						outputArray[i][j][k] = 0;
						//printf("warning: transformed array out of boundaries\n");
						//printLog("warning: transformed array out of boundaries\n");
						continue;
					}
					if ((N2 > 1) && ((y2 < xmin2) || (y2 > xmax2))) {
						outputArray[i][j][k] = 0;
						//printf("warning: transformed array out of boundaries\n");
						//printLog("warning: transformed array out of boundaries\n");
						continue;
					}
					if ((N3 > 1) && ((z2 < xmin3) || (z2 > xmax3))) {
						outputArray[i][j][k] = 0;
						//printf("warning: transformed array out of boundaries\n");
						//printLog("warning: transformed array out of boundaries\n");
						continue;
					}

					int i1 = floor((x2 - xmin1) / dx1);
					if (N1 == 1) {
						i1 = 0;
					}
					int i2 = floor((y2 - xmin2) / dx2);
					if (N2 == 1) {
						i2 = 0;
					}
					int i3 = floor((z2 - xmin3) / dx3);
					if (N3 == 1) {
						i3 = 0;
					}

					if ((i2 < 0) || (i2 >= N2)) {
						i2 = 0;
					}
					if ((i3 < 0) || (i3 >= N3)) {
						i3 = 0;
					}

					outputArray[i][j][k] = inputArray[i1][i2][i3];
				}
				else if (geometry == SourceInputGeometry::CYLINDRICAL) {
					double rho2 = sqrt(x2 * x2 + y2 * y2);
					double phi2 = atan2(y2, x2);
					if (phi2 < 0) {
						phi2 += 2 * pi;
					}

					if ((rho2 < xmin1) || (rho2 > xmax1)) {
						outputArray[i][j][k] = 0;
						//printf("warning: transformed array out of boundaries\n");
						//printLog("warning: transformed array out of boundaries\n");
						continue;
					}
					if ((N2 > 1) && ((z2 < xmin2) || (z2 > xmax2))) {
						outputArray[i][j][k] = 0;
						//printf("warning: transformed array out of boundaries\n");
						//printLog("warning: transformed array out of boundaries\n");
						continue;
					}
					if ((N3 > 1) && ((phi2 < xmin3) || (phi2 > xmax3))) {
						outputArray[i][j][k] = 0;
						//printf("warning: transformed array out of boundaries\n");
						//printLog("warning: transformed array out of boundaries\n");
						continue;
					}

					int i1 = floor((rho2 - xmin1) / dx1);
					int i2 = floor((z2 - xmin2) / dx2);
					int i3 = floor((phi2 - xmin3) / dx3);

					if ((i2 < 0) || (i2 >= N2)) {
						i2 = 0;
					}
					if ((i3 < 0) || (i3 >= N3)) {
						i3 = 0;
					}

					outputArray[i][j][k] = inputArray[i1][i2][i3];
				}
				else if (geometry == SourceInputGeometry::SPHERICAL) {
					double r2 = sqrt(x2 * x2 + y2 * y2 + z2 * z2);
					double theta2 = acos(z2 / r2);
					if (r2 == 0) {
						theta2 = 0;
					}
					double phi2 = atan2(y2, x2);
					if (phi2 < 0) {
						phi2 += 2 * pi;
					}

					if ((r2 < xmin1) || (r2 > xmax1)) {
						outputArray[i][j][k] = 0;
						//printf("warning: transformed array out of boundaries\n");
						//printLog("warning: transformed array out of boundaries\n");
						continue;
					}
					if ((N2 > 1) && ((theta2 < xmin2) || (theta2 > xmax2))) {
						outputArray[i][j][k] = 0;
						//printf("warning: transformed array out of boundaries\n");
						//printLog("warning: transformed array out of boundaries\n");
						continue;
					}
					if ((N3 > 1) && ((phi2 < xmin3) || (phi2 > xmax3))) {
						outputArray[i][j][k] = 0;
						//printf("warning: transformed array out of boundaries\n");
						//printLog("warning: transformed array out of boundaries\n");
						continue;
					}

					int i1 = floor((r2 - xmin1) / dx1);
					int i2 = floor((theta2 - xmin2) / dx2);
					int i3 = floor((phi2 - xmin3) / dx3);

					if ((i2 < 0) || (i2 >= N2)) {
						i2 = 0;
					}
					if ((i3 < 0) || (i3 >= N3)) {
						i3 = 0;
					}

					outputArray[i][j][k] = inputArray[i1][i2][i3];
				}

				if (checkNaNorInfinity(outputArray[i][j][k])) {
					printf("outputArray is Nan or Infinity %lf\n", outputArray[i][j][k]);
					printLog("outputArray is Nan or Infinity %lf\n", outputArray[i][j][k]);
					exit(0);
				}
			}
		}
	}
}

void RadiationSourceFactory::transformVectorArraysToCartesian(double*** inputArray1, double*** inputArray2, double*** inputArray3, int N1, int N2, int N3, const double& xmin1, const double& xmax1, const double& xmin2, const double& xmax2, const double& xmin3, const double& xmax3, SourceInputGeometry geometry, double*** outputArray, double*** outputArrayTheta, double*** outputArrayPhi, int Nx, int Nz, int Ny, const double& minX, const double& maxX, const double& minZ, const double& maxZ, const double& minY, const double& maxY, const double& thetar, const double& phir, const double& psir)
{
	if (maxX <= minX) {
		printf("rhomax = %g <= rhomin = %g in transformScalarArray\n", maxX, minX);
		printLog("rhomax = %g <= rhomin = %g in transformScalarArray\n", maxX, minX);
		exit(0);
	}
	if (maxZ <= minZ) {
		printf("zmax = %g <= zmin = %g in transformScalarArray\n", maxZ, minZ);
		printLog("zmax = %g <= zmin = %g in transformScalarArray\n", maxZ, minZ);
		exit(0);
	}
	if (maxY <= minY) {
		printf("phimax = %g <= phimin = %g in transformScalarArray\n", maxY, minY);
		printLog("phimax = %g <= phimin = %g in transformScalarArray\n", maxY, minY);
		exit(0);
	}
	if (xmax1 <= xmin1) {
		printf("xmax1 = %g <= xmin1 = %g in transformScalarArray\n", xmax1, xmin1);
		printLog("xmax1 = %g <= xmin1 = %g in transformScalarArray\n", xmax1, xmin1);
		exit(0);
	}
	if (xmax2 <= xmin2) {
		printf("xmax2 = %g <= xmin2 = %g in transformScalarArray\n", xmax2, xmin2);
		printLog("xmax2 = %g <= xmin2 = %g in transformScalarArray\n", xmax2, xmin2);
		exit(0);
	}
	if (xmax3 <= xmin3) {
		printf("xmax3 = %g <= xmin3 = %g in transformScalarArray\n", xmax3, xmin3);
		printLog("xmax3 = %g <= xmin3 = %g in transformScalarArray\n", xmax3, xmin3);
		exit(0);
	}

	double dx = (maxX - minX) / Nx;
	double dz = (maxZ - minZ) / Nz;
	double dy = (maxY - minY) / Ny;

	double dx1 = (xmax1 - xmin1) / N1;
	double dx2 = (xmax2 - xmin2) / N2;
	double dx3 = (xmax3 - xmin3) / N3;

	double sinthetar = sin(thetar);
	double costhetar = cos(thetar);
	double sinphir = sin(phir);
	double cosphir = cos(phir);
	double sinpsir = sin(psir);
	double cospsir = cos(psir);

	for (int i = 0; i < Nx; ++i) {
		double x = minX + dx * (i + 0.5);
		for (int j = 0; j < Nz; ++j) {
			double z = minZ + dz * (j + 0.5);
			for (int k = 0; k < Ny; ++k) {
				double y = minY + dy * (k + 0.5);

				double z1 = z;
				double x1 = x * cospsir - y*sinpsir;
				double y1 = y*cospsir + x*sinpsir;

				double z2 = z1 * cos(thetar) - x1 * sin(thetar);
				double x2 = (x1 * cos(thetar) + z1 * sin(thetar)) * cos(phir) - y1 * sin(phir);
				double y2 = y1 * cos(phir) + (x1 * cos(thetar) + z1 * sin(thetar)) * sin(phir);

				double ax, ay, az;

				if (geometry == SourceInputGeometry::CARTESIAN) {
					if ((x2 < xmin1) || (x2 > xmax1)) {
						outputArray[i][j][k] = 0;
						outputArrayTheta[i][j][k] = pi / 2;
						outputArrayTheta[i][j][k] = 0;
						//printf("warning: transformed array out of boundaries\n");
						//printLog("warning: transformed array out of boundaries\n");
						continue;
					}
					if ((N2 > 1) && ((y2 < xmin2) || (y2 > xmax2))) {
						outputArray[i][j][k] = 0;
						outputArrayTheta[i][j][k] = pi / 2;
						outputArrayTheta[i][j][k] = 0;
						//printf("warning: transformed array out of boundaries\n");
						//printLog("warning: transformed array out of boundaries\n");
						continue;
					}
					if ((N3 > 1) && ((z2 < xmin3) || (z2 > xmax3))) {
						outputArray[i][j][k] = 0;
						outputArrayTheta[i][j][k] = pi / 2;
						outputArrayTheta[i][j][k] = 0;
						//printf("warning: transformed array out of boundaries\n");
						//printLog("warning: transformed array out of boundaries\n");
						continue;
					}

					int i1 = floor((x2 - xmin1) / dx1);
					int i2 = floor((y2 - xmin2) / dx2);
					int i3 = floor((z2 - xmin3) / dx3);

					if ((i2 < 0) || (i2 >= N2)) {
						i2 = 0;
					}
					if ((i3 < 0) || (i3 >= N3)) {
						i3 = 0;
					}

					ax = inputArray1[i1][i2][i3];
					ay = inputArray2[i1][i2][i3];
					az = inputArray3[i1][i2][i3];

					outputArray[i][j][k] = sqrt(ax * ax + ay * ay + az * az);
					if (checkNaNorInfinity(outputArray[i][j][k])) {
						printf("outputArray is Nan or Infinity %lf\n", outputArray[i][j][k]);
						printLog("outputArray is Nan or Infinity %lf\n", outputArray[i][j][k]);
						exit(0);
					}
				}
				else if (geometry == SourceInputGeometry::CYLINDRICAL) {
					double rho2 = sqrt(x2 * x2 + y2 * y2);
					double phi2 = atan2(y2, x2);
					if (phi2 < 0) {
						phi2 += 2 * pi;
					}

					if ((rho2 < xmin1) || (rho2 > xmax1)) {
						outputArray[i][j][k] = 0;
						outputArrayTheta[i][j][k] = pi / 2;
						outputArrayTheta[i][j][k] = 0;
						//printf("warning: transformed array out of boundaries\n");
						//printLog("warning: transformed array out of boundaries\n");
						continue;
					}
					if ((N2 > 1) && ((z2 < xmin2) || (z2 > xmax2))) {
						outputArray[i][j][k] = 0;
						outputArrayTheta[i][j][k] = pi / 2;
						outputArrayTheta[i][j][k] = 0;
						//printf("warning: transformed array out of boundaries\n");
						//printLog("warning: transformed array out of boundaries\n");
						continue;
					}
					if ((N3 > 1) && ((phi2 < xmin3) || (phi2 > xmax3))) {
						outputArray[i][j][k] = 0;
						outputArrayTheta[i][j][k] = pi / 2;
						outputArrayTheta[i][j][k] = 0;
						//printf("warning: transformed array out of boundaries\n");
						//printLog("warning: transformed array out of boundaries\n");
						continue;
					}

					int i1 = floor((rho2 - xmin1) / dx1);
					int i2 = floor((z2 - xmin2) / dx2);
					int i3 = floor((phi2 - xmin3) / dx3);

					if ((i2 < 0) || (i2 >= N2)) {
						i2 = 0;
					}
					if ((i3 < 0) || (i3 >= N3)) {
						i3 = 0;
					}

					double arho = inputArray1[i1][i2][i3];
					az = inputArray2[i1][i2][i3];
					double aphi = inputArray3[i1][i2][i3];

					ax = arho * cos(phi2) - aphi * sin(phi2);
					ay = aphi * cos(phi2) + arho * sin(phi2);

					outputArray[i][j][k] = sqrt(arho * arho + az * az + aphi * aphi);
					if (checkNaNorInfinity(outputArray[i][j][k])) {
						printf("outputArray is Nan or Infinity %lf\n", outputArray[i][j][k]);
						printLog("outputArray is Nan or Infinity %lf\n", outputArray[i][j][k]);
						exit(0);
					}
				}
				else if (geometry == SourceInputGeometry::SPHERICAL) {
					double r2 = sqrt(x2 * x2 + y2 * y2 + z2 * z2);
					double theta2 = acos(z2 / r2);
					if (r2 == 0) {
						theta2 = 0;
					}
					double phi2 = atan2(y2, x2);
					if (phi2 < 0) {
						phi2 += 2 * pi;
					}

					if ((r2 < xmin1) || (r2 > xmax1)) {
						outputArray[i][j][k] = 0;
						outputArrayTheta[i][j][k] = pi / 2;
						outputArrayTheta[i][j][k] = 0;
						//printf("warning: transformed array out of boundaries\n");
						//printLog("warning: transformed array out of boundaries\n");
						continue;
					}
					if ((N2 > 1) && ((theta2 < xmin2) || (theta2 > xmax2))) {
						outputArray[i][j][k] = 0;
						outputArrayTheta[i][j][k] = pi / 2;
						outputArrayTheta[i][j][k] = 0;
						//printf("warning: transformed array out of boundaries\n");
						//printLog("warning: transformed array out of boundaries\n");
						continue;
					}
					if ((N3 > 1) && ((phi2 < xmin3) || (phi2 > xmax3))) {
						outputArray[i][j][k] = 0;
						outputArrayTheta[i][j][k] = pi / 2;
						outputArrayTheta[i][j][k] = 0;
						//printf("warning: transformed array out of boundaries\n");
						//printLog("warning: transformed array out of boundaries\n");
						continue;
					}

					int i1 = floor((r2 - xmin1) / dx1);
					int i2 = floor((theta2 - xmin2) / dx2);
					int i3 = floor((phi2 - xmin3) / dx3);

					if ((i2 < 0) || (i2 >= N2)) {
						i2 = 0;
					}
					if ((i3 < 0) || (i3 >= N3)) {
						i3 = 0;
					}

					double ar = inputArray1[i1][i2][i3];
					double atheta = inputArray2[i1][i2][i3];
					double aphi = inputArray3[i1][i2][i3];

					az = ar * cos(theta2) + atheta * sin(theta2);
					ax = ar * sin(theta2) * cos(phi2) - atheta * cos(theta2) * cos(phi2) - aphi * sin(phi2);
					ay = ar * sin(theta2) * sin(phi2) - atheta * cos(theta2) * sin(phi2) + aphi * cos(phi2);

					outputArray[i][j][k] = sqrt(ar * ar + atheta * atheta + aphi * aphi);
					if (checkNaNorInfinity(outputArray[i][j][k])) {
						printf("outputArray is Nan or Infinity %lf\n", outputArray[i][j][k]);
						printLog("outputArray is Nan or Infinity %lf\n", outputArray[i][j][k]);
						exit(0);
					}
				}


				if (outputArray[i][j][k] == 0) {
					outputArrayTheta[i][j][k] = pi / 2;
					outputArrayPhi[i][j][k] = 0;
					continue;
				}

				double amultz = (ax * sinthetar * cosphir + ay * sinthetar * sinphir + az * costhetar);
				double amultx = (ax * (costhetar * cosphir * cospsir - sinpsir * sinphir) + ay * (sinpsir * cosphir + costhetar * cospsir * sinphir) - az * sinthetar * cospsir);
				double amulty = (ax * (-sinpsir * costhetar * cosphir - cospsir * sinphir) + ay * (cospsir * cosphir - sinpsir * costhetar * sinphir) + az * sinthetar * sinpsir);

				outputArrayTheta[i][j][k] = acos(amultz / outputArray[i][j][k]);
				outputArrayPhi[i][j][k] = atan2(amulty, amultx);
				if (outputArrayPhi[i][j][k] < 0) {
					outputArrayPhi[i][j][k] += 2 * pi;
				}
			}
		}
	}
}

//transforms array from arbitrary pluto coordinates to rotated cylindrical coordinates
void RadiationSourceFactory::transformScalarArrayToCylindrical(double*** inputArray, int N1, int N2, int N3, const double& xmin1, const double& xmax1, const double& xmin2, const double& xmax2, const double& xmin3, const double& xmax3, SourceInputGeometry geometry, double*** outputArray, int Nrho, int Nz, int Nphi, const double& rhomin, const double& rhomax, const double& zmin, const double& zmax, const double& phimin, const double& phimax, const double& thetar, const double& phir, const double& psir) {
	if (rhomax <= rhomin) {
		printf("rhomax = %g <= rhomin = %g in transformScalarArray\n", rhomax, rhomin);
		printLog("rhomax = %g <= rhomin = %g in transformScalarArray\n", rhomax, rhomin);
		exit(0);
	}
	if (zmax <= zmin) {
		printf("zmax = %g <= zmin = %g in transformScalarArray\n", zmax, zmin);
		printLog("zmax = %g <= zmin = %g in transformScalarArray\n", zmax, zmin);
		exit(0);
	}
	if (phimax <= phimin) {
		printf("phimax = %g <= phimin = %g in transformScalarArray\n", phimax, phimin);
		printLog("phimax = %g <= phimin = %g in transformScalarArray\n", phimax, phimin);
		exit(0);
	}
	if (xmax1 <= xmin1) {
		printf("xmax1 = %g <= xmin1 = %g in transformScalarArray\n", xmax1, xmin1);
		printLog("xmax1 = %g <= xmin1 = %g in transformScalarArray\n", xmax1, xmin1);
		exit(0);
	}
	if (xmax2 <= xmin2) {
		printf("xmax2 = %g <= xmin2 = %g in transformScalarArray\n", xmax2, xmin2);
		printLog("xmax2 = %g <= xmin2 = %g in transformScalarArray\n", xmax2, xmin2);
		exit(0);
	}
	if (xmax3 <= xmin3) {
		printf("xmax3 = %g <= xmin3 = %g in transformScalarArray\n", xmax3, xmin3);
		printLog("xmax3 = %g <= xmin3 = %g in transformScalarArray\n", xmax3, xmin3);
		exit(0);
	}

	double drho = (rhomax - rhomin) / Nrho;
	double dz = (zmax - zmin) / Nz;
	double dphi = (phimax - phimin) / Nphi;

	double dx1 = (xmax1 - xmin1) / N1;
	double dx2 = (xmax2 - xmin2) / N2;
	double dx3 = (xmax3 - xmin3) / N3;

	for (int i = 0; i < Nrho; ++i) {
		double rho = rhomin + drho * (i + 0.5);
		for (int j = 0; j < Nz; ++j) {
			double z = zmin + dz * (j + 0.5);
			for (int k = 0; k < Nphi; ++k) {
				double phi = phimin + dphi * (k + 0.5);

				double z1 = z;
				double x1 = rho * cos(phi + psir);
				double y1 = rho * sin(phi + psir);

				double z2 = z1 * cos(thetar) - x1 * sin(thetar);
				double x2 = (x1 * cos(thetar) + z1 * sin(thetar)) * cos(phir) - y1 * sin(phir);
				double y2 = y1 * cos(phir) + (x1 * cos(thetar) + z1 * sin(thetar)) * sin(phir);

				if (geometry == SourceInputGeometry::CARTESIAN) {
					if ((x2 < xmin1) || (x2 > xmax1)) {
						outputArray[i][j][k] = 0;
						//printf("warning: transformed array out of boundaries\n");
						//printLog("warning: transformed array out of boundaries\n");
						continue;
					}
					if ((N2 > 1) && ((y2 < xmin2) || (y2 > xmax2))) {
						outputArray[i][j][k] = 0;
						//printf("warning: transformed array out of boundaries\n");
						//printLog("warning: transformed array out of boundaries\n");
						continue;
					}
					if ((N3 > 1) && ((z2 < xmin3) || (z2 > xmax3))) {
						outputArray[i][j][k] = 0;
						//printf("warning: transformed array out of boundaries\n");
						//printLog("warning: transformed array out of boundaries\n");
						continue;
					}

					int i1 = floor((x2 - xmin1) / dx1);
					int i2 = floor((y2 - xmin2) / dx2);
					int i3 = floor((z2 - xmin3) / dx3);

					if ((i2 < 0) || (i2 >= N2)) {
						i2 = 0;
					}
					if ((i3 < 0) || (i3 >= N3)) {
						i3 = 0;
					}

					outputArray[i][j][k] = inputArray[i1][i2][i3];
				}
				else if (geometry == SourceInputGeometry::CYLINDRICAL) {
					double rho2 = sqrt(x2 * x2 + y2 * y2);
					double phi2 = atan2(y2, x2);
					if (phi2 < 0) {
						phi2 += 2 * pi;
					}

					if ((rho2 < xmin1) || (rho2 > xmax1)) {
						outputArray[i][j][k] = 0;
						//printf("warning: transformed array out of boundaries\n");
						//printLog("warning: transformed array out of boundaries\n");
						continue;
					}
					if ((N2 > 1) && ((z2 < xmin2) || (z2 > xmax2))) {
						outputArray[i][j][k] = 0;
						//printf("warning: transformed array out of boundaries\n");
						//printLog("warning: transformed array out of boundaries\n");
						continue;
					}
					if ((N3 > 1) && ((phi2 < xmin3) || (phi2 > xmax3))) {
						outputArray[i][j][k] = 0;
						//printf("warning: transformed array out of boundaries\n");
						//printLog("warning: transformed array out of boundaries\n");
						continue;
					}

					int i1 = floor((rho2 - xmin1) / dx1);
					int i2 = floor((z2 - xmin2) / dx2);
					int i3 = floor((phi2 - xmin3) / dx3);

					if ((i2 < 0) || (i2 >= N2)) {
						i2 = 0;
					}
					if ((i3 < 0) || (i3 >= N3)) {
						i3 = 0;
					}

					outputArray[i][j][k] = inputArray[i1][i2][i3];
				}
				else if (geometry == SourceInputGeometry::SPHERICAL) {
					double r2 = sqrt(x2 * x2 + y2 * y2 + z2 * z2);
					double theta2 = acos(z2 / r2);
					if (r2 == 0) {
						theta2 = 0;
					}
					double phi2 = atan2(y2, x2);
					if (phi2 < 0) {
						phi2 += 2 * pi;
					}

					if ((r2 < xmin1) || (r2 > xmax1)) {
						outputArray[i][j][k] = 0;
						//printf("warning: transformed array out of boundaries\n");
						//printLog("warning: transformed array out of boundaries\n");
						continue;
					}
					if ((N2 > 1) && ((theta2 < xmin2) || (theta2 > xmax2))) {
						outputArray[i][j][k] = 0;
						//printf("warning: transformed array out of boundaries\n");
						//printLog("warning: transformed array out of boundaries\n");
						continue;
					}
					if ((N3 > 1) && ((phi2 < xmin3) || (phi2 > xmax3))) {
						outputArray[i][j][k] = 0;
						//printf("warning: transformed array out of boundaries\n");
						//printLog("warning: transformed array out of boundaries\n");
						continue;
					}

					int i1 = floor((r2 - xmin1) / dx1);
					int i2 = floor((theta2 - xmin2) / dx2);
					int i3 = floor((phi2 - xmin3) / dx3);

					if ((i2 < 0) || (i2 >= N2)) {
						i2 = 0;
					}
					if ((i3 < 0) || (i3 >= N3)) {
						i3 = 0;
					}

					outputArray[i][j][k] = inputArray[i1][i2][i3];
				}

				if (checkNaNorInfinity(outputArray[i][j][k])) {
					printf("outputArray is Nan or Infinity %lf\n", outputArray[i][j][k]);
					printLog("outputArray is Nan or Infinity %lf\n", outputArray[i][j][k]);
					exit(0);
				}
			}
		}
	}
}

void RadiationSourceFactory::transformVectorArraysToCylindrical(double*** inputArray1, double*** inputArray2, double*** inputArray3, int N1, int N2, int N3, const double& xmin1, const double& xmax1, const double& xmin2, const double& xmax2, const double& xmin3, const double& xmax3, SourceInputGeometry geometry, double*** outputArray, double*** outputArrayTheta, double*** outputArrayPhi, int Nrho, int Nz, int Nphi, const double& rhomin, const double& rhomax, const double& zmin, const double& zmax, const double& phimin, const double& phimax, const double& thetar, const double& phir, const double& psir){
	if (rhomax <= rhomin) {
		printf("rhomax = %g <= rhomin = %g in transformScalarArray\n", rhomax, rhomin);
		printLog("rhomax = %g <= rhomin = %g in transformScalarArray\n", rhomax, rhomin);
		exit(0);
	}
	if (zmax <= zmin) {
		printf("zmax = %g <= zmin = %g in transformScalarArray\n", zmax, zmin);
		printLog("zmax = %g <= zmin = %g in transformScalarArray\n", zmax, zmin);
		exit(0);
	}
	if (phimax <= phimin) {
		printf("phimax = %g <= phimin = %g in transformScalarArray\n", phimax, phimin);
		printLog("phimax = %g <= phimin = %g in transformScalarArray\n", phimax, phimin);
		exit(0);
	}
	if (xmax1 <= xmin1) {
		printf("xmax1 = %g <= xmin1 = %g in transformScalarArray\n", xmax1, xmin1);
		printLog("xmax1 = %g <= xmin1 = %g in transformScalarArray\n", xmax1, xmin1);
		exit(0);
	}
	if (xmax2 <= xmin2) {
		printf("xmax2 = %g <= xmin2 = %g in transformScalarArray\n", xmax2, xmin2);
		printLog("xmax2 = %g <= xmin2 = %g in transformScalarArray\n", xmax2, xmin2);
		exit(0);
	}
	if (xmax3 <= xmin3) {
		printf("xmax3 = %g <= xmin3 = %g in transformScalarArray\n", xmax3, xmin3);
		printLog("xmax3 = %g <= xmin3 = %g in transformScalarArray\n", xmax3, xmin3);
		exit(0);
	}

	double drho = (rhomax - rhomin) / Nrho;
	double dz = (zmax - zmin) / Nz;
	double dphi = (phimax - phimin) / Nphi;

	double dx1 = (xmax1 - xmin1) / N1;
	double dx2 = (xmax2 - xmin2) / N2;
	double dx3 = (xmax3 - xmin3) / N3;

	double sinthetar = sin(thetar);
	double costhetar = cos(thetar);
	double sinphir = sin(phir);
	double cosphir = cos(phir);
	double sinpsir = sin(psir);
	double cospsir = cos(psir);

	for (int i = 0; i < Nrho; ++i) {
		double rho = rhomin + drho * (i + 0.5);
		for (int j = 0; j < Nz; ++j) {
			double z = zmin + dz * (j + 0.5);
			for (int k = 0; k < Nphi; ++k) {
				double phi = phimin + dphi * (k + 0.5);

				double z1 = z;
				double x1 = rho * cos(phi + psir);
				double y1 = rho * sin(phi + psir);

				double z2 = z1 * cos(thetar) - x1 * sin(thetar);
				double x2 = (x1 * cos(thetar) + z1 * sin(thetar)) * cos(phir) - y1 * sin(phir);
				double y2 = y1 * cos(phir) + (x1 * cos(thetar) + z1 * sin(thetar)) * sin(phir);

				double ax, ay, az;

				if (geometry == SourceInputGeometry::CARTESIAN) {
					if ((x2 < xmin1) || (x2 > xmax1)) {
						outputArray[i][j][k] = 0;
						outputArrayTheta[i][j][k] = pi / 2;
						outputArrayTheta[i][j][k] = 0;
						//printf("warning: transformed array out of boundaries\n");
						//printLog("warning: transformed array out of boundaries\n");
						continue;
					}
					if ((N2 > 1) && ((y2 < xmin2) || (y2 > xmax2))) {
						outputArray[i][j][k] = 0;
						outputArrayTheta[i][j][k] = pi / 2;
						outputArrayTheta[i][j][k] = 0;
						//printf("warning: transformed array out of boundaries\n");
						//printLog("warning: transformed array out of boundaries\n");
						continue;
					}
					if ((N3 > 1) && ((z2 < xmin3) || (z2 > xmax3))) {
						outputArray[i][j][k] = 0;
						outputArrayTheta[i][j][k] = pi / 2;
						outputArrayTheta[i][j][k] = 0;
						//printf("warning: transformed array out of boundaries\n");
						//printLog("warning: transformed array out of boundaries\n");
						continue;
					}

					int i1 = floor((x2 - xmin1) / dx1);
					int i2 = floor((y2 - xmin2) / dx2);
					int i3 = floor((z2 - xmin3) / dx3);

					if ((i2 < 0) || (i2 >= N2)) {
						i2 = 0;
					}
					if ((i3 < 0) || (i3 >= N3)) {
						i3 = 0;
					}

					ax = inputArray1[i1][i2][i3];
					ay = inputArray2[i1][i2][i3];
					az = inputArray3[i1][i2][i3];

					outputArray[i][j][k] = sqrt(ax*ax + ay*ay + az*az);
					if (checkNaNorInfinity(outputArray[i][j][k])) {
						printf("outputArray is Nan or Infinity %lf\n", outputArray[i][j][k]);
						printLog("outputArray is Nan or Infinity %lf\n", outputArray[i][j][k]);
						exit(0);
					}
				}
				else if (geometry == SourceInputGeometry::CYLINDRICAL) {
					double rho2 = sqrt(x2 * x2 + y2 * y2);
					double phi2 = atan2(y2, x2);
					if (phi2 < 0) {
						phi2 += 2 * pi;
					}

					if ((rho2 < xmin1) || (rho2 > xmax1)) {
						outputArray[i][j][k] = 0;
						outputArrayTheta[i][j][k] = pi / 2;
						outputArrayTheta[i][j][k] = 0;
						//printf("warning: transformed array out of boundaries\n");
						//printLog("warning: transformed array out of boundaries\n");
						continue;
					}
					if ((N2 > 1) && ((z2 < xmin2) || (z2 > xmax2))) {
						outputArray[i][j][k] = 0;
						outputArrayTheta[i][j][k] = pi / 2;
						outputArrayTheta[i][j][k] = 0;
						//printf("warning: transformed array out of boundaries\n");
						//printLog("warning: transformed array out of boundaries\n");
						continue;
					}
					if ((N3 > 1) && ((phi2 < xmin3) || (phi2 > xmax3))) {
						outputArray[i][j][k] = 0;
						outputArrayTheta[i][j][k] = pi / 2;
						outputArrayTheta[i][j][k] = 0;
						//printf("warning: transformed array out of boundaries\n");
						//printLog("warning: transformed array out of boundaries\n");
						continue;
					}

					int i1 = floor((rho2 - xmin1) / dx1);
					int i2 = floor((z2 - xmin2) / dx2);
					int i3 = floor((phi2 - xmin3) / dx3);

					if ((i2 < 0) || (i2 >= N2)) {
						i2 = 0;
					}
					if ((i3 < 0) || (i3 >= N3)) {
						i3 = 0;
					}

					double arho = inputArray1[i1][i2][i3];
					az = inputArray2[i1][i2][i3];
					double aphi = inputArray3[i1][i2][i3];

					ax = arho * cos(phi2) - aphi * sin(phi2);
					ay = aphi * cos(phi2) + arho * sin(phi2);

					outputArray[i][j][k] = sqrt(arho * arho + az * az + aphi * aphi);
					if (checkNaNorInfinity(outputArray[i][j][k])) {
						printf("outputArray is Nan or Infinity %lf\n", outputArray[i][j][k]);
						printLog("outputArray is Nan or Infinity %lf\n", outputArray[i][j][k]);
						exit(0);
					}
				}
				else if (geometry == SourceInputGeometry::SPHERICAL) {
					double r2 = sqrt(x2 * x2 + y2 * y2 + z2 * z2);
					double theta2 = acos(z2 / r2);
					if (r2 == 0) {
						theta2 = 0;
					}
					double phi2 = atan2(y2, x2);
					if (phi2 < 0) {
						phi2 += 2 * pi;
					}

					if ((r2 < xmin1) || (r2 > xmax1)) {
						outputArray[i][j][k] = 0;
						outputArrayTheta[i][j][k] = pi / 2;
						outputArrayTheta[i][j][k] = 0;
						//printf("warning: transformed array out of boundaries\n");
						//printLog("warning: transformed array out of boundaries\n");
						continue;
					}
					if ((N2 > 1) && ((theta2 < xmin2) || (theta2 > xmax2))) {
						outputArray[i][j][k] = 0;
						outputArrayTheta[i][j][k] = pi / 2;
						outputArrayTheta[i][j][k] = 0;
						//printf("warning: transformed array out of boundaries\n");
						//printLog("warning: transformed array out of boundaries\n");
						continue;
					}
					if ((N3 > 1) && ((phi2 < xmin3) || (phi2 > xmax3))) {
						outputArray[i][j][k] = 0;
						outputArrayTheta[i][j][k] = pi / 2;
						outputArrayTheta[i][j][k] = 0;
						//printf("warning: transformed array out of boundaries\n");
						//printLog("warning: transformed array out of boundaries\n");
						continue;
					}

					int i1 = floor((r2 - xmin1) / dx1);
					int i2 = floor((theta2 - xmin2) / dx2);
					int i3 = floor((phi2 - xmin3) / dx3);

					if ((i2 < 0) || (i2 >= N2)) {
						i2 = 0;
					}
					if ((i3 < 0) || (i3 >= N3)) {
						i3 = 0;
					}

					double ar = inputArray1[i1][i2][i3];
					double atheta = inputArray2[i1][i2][i3];
					double aphi = inputArray3[i1][i2][i3];

					az = ar * cos(theta2) + atheta * sin(theta2);
					ax = ar * sin(theta2) * cos(phi2) - atheta * cos(theta2) * cos(phi2) - aphi * sin(phi2);
					ay = ar * sin(theta2) * sin(phi2) - atheta * cos(theta2) * sin(phi2) + aphi * cos(phi2);

					outputArray[i][j][k] = sqrt(ar * ar + atheta * atheta + aphi * aphi);
					if (checkNaNorInfinity(outputArray[i][j][k])) {
						printf("outputArray is Nan or Infinity %lf\n", outputArray[i][j][k]);
						printLog("outputArray is Nan or Infinity %lf\n", outputArray[i][j][k]);
						exit(0);
					}
				}


				if (outputArray[i][j][k] == 0) {
					outputArrayTheta[i][j][k] = pi / 2;
					outputArrayPhi[i][j][k] = 0;
					continue;
				}

				double amultz = (ax * sinthetar * cosphir + ay * sinthetar * sinphir + az * costhetar);
				double amultx = (ax * (costhetar * cosphir * cospsir - sinpsir * sinphir) + ay * (sinpsir * cosphir + costhetar * cospsir * sinphir) - az * sinthetar * cospsir);
				double amulty = (ax * (-sinpsir * costhetar * cosphir - cospsir * sinphir) + ay * (cospsir * cosphir - sinpsir * costhetar * sinphir) + az * sinthetar * sinpsir);

				outputArrayTheta[i][j][k] = acos(amultz / outputArray[i][j][k]);
				outputArrayPhi[i][j][k] = atan2(amulty, amultx);
				if (outputArrayPhi[i][j][k] < 0) {
					outputArrayPhi[i][j][k] += 2 * pi;
				}
			}
		}
	}
}

bool*** RadiationSourceFactory::breadthFirstSearchOfRegion(double*** inputArray, const int N1, const int N2, const int N3, const int start1, const int start2, const int start3, const double& threshold) {
	
	bool*** result = new bool** [N1];
	for (int i = 0; i < N1; ++i) {
		result[i] = new bool* [N2];
		for (int j = 0; j < N2; ++j) {
			result[i][j] = new bool[N3];
			for (int k = 0; k < N3; ++k) {
				result[i][j][k] = false;
			}
		}
	}
	
	std::queue<std::tuple<int, int, int> > queue;

	queue.push(std::tuple<int, int, int>(start1, start2, start3));
	result[start1][start2][start3] = true;

	while (queue.size() > 0) {
		std::tuple<int, int, int> tuple = queue.front();
		queue.pop();
		int n1 = std::get<0>(tuple);
		int n2 = std::get<1>(tuple);
		int n3 = std::get<2>(tuple);

		//result[n1][n2][n3] = true;
		double current = inputArray[n1][n2][n3];

		if (n1 > 0) {
			double local = inputArray[n1 - 1][n2][n3];
			if ((!result[n1-1][n2][n3]) && (local > (1.0 - threshold)* current) && (local < (1.0 + threshold) * current)) {
				queue.push(std::tuple<int, int, int>(n1 - 1, n2, n3));
				result[n1 - 1][n2][n3] = true;
			}
		}
		if (n1 < N1 - 1) {
			double local = inputArray[n1 + 1][n2][n3];
			if ((!result[n1+1][n2][n3]) && (local >(1.0 - threshold)* current) && (local < (1.0 + threshold) * current)) {
				queue.push(std::tuple<int, int, int>(n1 + 1, n2, n3));
				result[n1 + 1][n2][n3] = true;
			}
		}

		if (n2 > 0) {
			double local = inputArray[n1][n2-1][n3];
			if ((!result[n1][n2-1][n3]) && (local > (1.0 - threshold) * current) && (local < (1.0 + threshold) * current)) {
				queue.push(std::tuple<int, int, int>(n1, n2-1, n3));
				result[n1][n2 - 1][n3] = true;
			}
		}
		if (n2 < N2 - 1) {
			double local = inputArray[n1][n2 + 1][n3];
			if ((!result[n1][n2 + 1][n3]) && (local > (1.0 - threshold) * current) && (local < (1.0 + threshold) * current)) {
				queue.push(std::tuple<int, int, int>(n1, n2 + 1, n3));
				result[n1][n2 + 1][n3] = true;
			}
		}
		if (n3 > 0) {
			double local = inputArray[n1][n2][n3-1];
			if ((!result[n1][n2][n3-1]) && (local > (1.0 - threshold) * current) && (local < (1.0 + threshold) * current)) {
				queue.push(std::tuple<int, int, int>(n1, n2, n3-1));
				result[n1][n2][n3 - 1] = true;
			}
		}
		if (n3 < N3 - 1) {
			double local = inputArray[n1][n2][n3+1];
			if ((!result[n1][n2][n3 + 1]) && (local > (1.0 - threshold) * current) && (local < (1.0 + threshold) * current)) {
				queue.push(std::tuple<int, int, int>(n1, n2, n3+1));
				result[n1][n2][n3 + 1] = true;
			}
		}
	}

	return result;
}

void RadiationSourceFactory::evaluateDistributionAfterDiffusion(double* energy, double* distribution, const int Ne, double* D, const double& x, const double& y, const double& z, const double& time, const int Nt, double* sourceDistribution, void(*sourceCoordinates)(const double& time, double& x1, double& x2, double& x3), double (*sourcePower)(const double& time))
{
	double dt = time / Nt;

	for (int j = 0; j < Ne; ++j) {
		distribution[j] = 0;
	}

	for (int i = 0; i < Nt; ++i) {
		double currentTime = i * dt;
		double power = (*sourcePower)(currentTime);
		double sourceX, sourceY, sourceZ;
		(*sourceCoordinates)(currentTime, sourceX, sourceY, sourceZ);
		double deltaT = time - currentTime;
		double d2 = sqr(x - sourceX) + sqr(y - sourceY) + sqr(z - sourceZ);

		for (int j = 0; j < Ne; ++j) {
			distribution[j] += dt*power * sourceDistribution[j] * exp(-d2 / (D[j] * deltaT)) / cube(sqrt(4 * pi * D[j] * deltaT));
		}
	}
}

RectangularSourceInhomogenousDistribution* RadiationSourceFactory::createRectangularSourceFromDiffusion(const double& mass, double* energy, double* sourceDistribution, const int Ne, double* D, const int Nx, const int Ny, const int Nz, const double& minX, const double& maxX, const double& minY, const double& maxY, const double& minZ, const double& maxZ, const double& time, const int Nt, void (*sourceCoordinates)(const double& time, double& x1, double& x2, double& x3), double (*sourcePower)(const double& time), const double& B, const double& theta, const double& phi, const double& distance, const double& velocity, const double& redShift)
{

	MassiveParticleDistribution**** distributions1;
	double*** concentration;

	RadiationSourceFactory::createRectangularSourceArraysFromDiffusion(mass, energy, sourceDistribution, Ne, D, Nx, Ny, Nz, minX, maxX, minY, maxY, minZ, maxZ, time, Nt, sourceCoordinates, sourcePower, distributions1, concentration);

	RectangularSourceInhomogenousDistribution* source = new RectangularSourceInhomogenousDistribution(Nx, Ny, Nz, distributions1, B, theta, phi, concentration, minX, maxX, minY, maxY, minZ, maxZ, distance, velocity, redShift);

	for (int i = 0; i < Nx; ++i) {
		for (int j = 0; j < Nz; ++j) {
			delete[] distributions1[i][j];
			delete[] concentration[i][j];
		}
		delete[] distributions1[i];
		delete[] concentration[i];
	}
	delete[] distributions1;
	delete[] concentration;

	return source;
}

void RadiationSourceFactory::createRectangularSourceArraysFromDiffusion(const double& mass, double* energy, double* sourceDistribution, const int Ne, double* D, const int Nx, const int Ny, const int Nz, const double& minX, const double& maxX, const double& minY, const double& maxY, const double& minZ, const double& maxZ, const double& time, const int Nt, void(*sourceCoordinates)(const double& time, double& x1, double& x2, double& x3), double(*sourcePower)(const double& time), MassiveParticleDistribution****& distributions1, double***& concentration)
{
	double dx = (maxX - minX) / Nx;
	double dy = (maxY - minY) / Ny;
	double dz = (maxZ - minZ) / Nz;
	double**** distributions = new double*** [Nx];
	distributions1 = new MassiveParticleDistribution ***[Nx];
	concentration = new double** [Nx];
	for (int i = 0; i < Nx; ++i) {
		printf("evaluzting diffusion for i = %d\n", i);
		distributions[i] = new double** [Nz];
		distributions1[i] = new MassiveParticleDistribution **[Nz];
		concentration[i] = new double* [Nz];
		double x = minX + (i + 0.5) * dx;
		for (int j = 0; j < Nz; ++j) {
			distributions[i][j] = new double* [Ny];
			distributions1[i][j] = new MassiveParticleDistribution * [Ny];
			concentration[i][j] = new double[Ny];
			double z = minZ + (j + 0.5) * dz;
			for (int k = 0; k < Ny; ++k) {
				distributions[i][j][k] = new double[Ne];
				double y = minY + (k + 0.5) * dy;
				evaluateDistributionAfterDiffusion(energy, distributions[i][j][k], Ne, D, x, y, z, time, Nt, sourceDistribution, sourceCoordinates, sourcePower);
				concentration[i][j][k] = MassiveParticleDistributionFactory::evaluateNorm(energy, distributions[i][j][k], Ne);
				distributions1[i][j][k] = new MassiveParticleTabulatedIsotropicDistribution(mass, energy, distributions[i][j][k], Ne, DistributionInputType::ENERGY_FE);
			}
		}
	}

	for (int i = 0; i < Nx; ++i) {
		for (int j = 0; j < Nz; ++j) {
			for (int k = 0; k < Ny; ++k) {
				delete[] distributions[i][j][k];
			}
			delete[] distributions[i][j];
		}
		delete[] distributions[i];
	}
	delete[] distributions;
}
