#include "stdio.h"
#include "math.h"
#include <omp.h>

#include "constants.h"
#include "massiveParticleDistribution.h"
#include "photonDistribution.h"
#include "radiationSource.h"
#include "radiation.h"

RadiationEvaluator::RadiationEvaluator(int Ne, const double& Emin, const double& Emax) {
    my_Ne = Ne;
    my_Emin = Emin;
    my_Emax = Emax;
    my_Ee = new double[my_Ne];
    double factor = pow(my_Emax / my_Emin, 1.0 / (my_Ne - 1));

    my_Ee[0] = my_Emin;
    for (int i = 1; i < my_Ne; ++i) {
        my_Ee[i] = my_Ee[i - 1] * factor;
    }
}
RadiationEvaluator::~RadiationEvaluator() {
    delete[] my_Ee;
}

double RadiationEvaluator::evaluateFluxFromSource(const double& photonFinalEnergy, RadiationSource* source) {
    int Nrho = source->getNrho();
    int Nz = source->getNz();
    int Nphi = source->getNphi();

    double result = 0;
    int irho = 0;

    omp_init_lock(&my_lock);

#pragma omp parallel for private(irho) shared(photonFinalEnergy, source, Nrho, Nz, Nphi) reduction(+:result)

    for (irho = 0; irho < Nrho; ++irho) {
        for (int iphi = 0; iphi < Nphi; ++iphi) {
            /*for (int iz = 0; iz < Nz; ++iz) {
                result += evaluateFluxFromIsotropicFunction(photonFinalEnergy, source->getParticleDistribution(irho, iz, iphi), source->getVolume(irho, iz, iphi), source->getDistance());
            }*/
            if (source->isSource(irho, iphi)) {
                result += evaluateFluxFromSourceAtPoint(photonFinalEnergy, source, irho, iphi);
            }
        }
    }

    omp_destroy_lock(&my_lock);

    return result;
}

double RadiationEvaluator::evaluateTotalFluxInEnergyRange(const double& Ephmin, const double& Ephmax, int Nph, RadiationSource* source) {
    double factor = pow(Ephmax / Ephmin, 1.0 / (Nph - 1));
    double currentE = Ephmin;
    double flux = 0;
    for (int i = 0; i < Nph; ++i) {
        printf("%d\n", i);
        double dE = currentE * (factor - 1.0);
        flux += evaluateFluxFromSource(currentE, source) * dE;
        currentE = currentE * factor;
    }
    return flux;
}

void RadiationEvaluator::writeFluxFromSourceToFile(const char* fileName, RadiationSource* source, const double& Ephmin, const double& Ephmax, const int Nph) {
    double factor = pow(Ephmax / Ephmin, 1.0 / (Nph - 1));
    double currentE = Ephmin;
    FILE* outFile = fopen(fileName, "w");
    for (int i = 0; i < Nph; ++i) {
        printf("%d\n", i);
        double flux = evaluateFluxFromSource(currentE, source);
        fprintf(outFile, "%g %g\n", currentE, flux);
        currentE = currentE * factor;
    }
    fclose(outFile);
}

void RadiationEvaluator::writeImageFromSourceToFile(const char* fileName, RadiationSource* source, const double& Ephmin, const double& Ephmax, const int Nph) {
    int Nrho = source->getNrho();
    int Nz = source->getNz();
    int Nphi = source->getNphi();

    double** image = new double*[Nrho];
    for (int irho = 0; irho < Nrho; ++irho) {
        image[irho] = new double[Nphi];
        for (int iphi = 0; iphi < Nphi; ++iphi) {
            image[irho][iphi] = 0;
        }
    }

    omp_init_lock(&my_lock);

    int irho;
#pragma omp parallel for private(irho) shared(Ephmin, Ephmax, source, Nrho, Nz, Nphi, image)
    for (irho = 0; irho < Nrho; ++irho) {
        printf("evaluating image irho = %d\n", irho);
        for (int iphi = 0; iphi < Nphi; ++iphi) {
            double factor = pow(Ephmax / Ephmin, 1.0 / (Nph - 1));
            double currentE = Ephmin;
            double localFlux = 0;
            double rho0 = irho * source->getMaxRho() / Nrho;
            double rho1 = (irho + 1) * source->getMaxRho() / Nrho;
            double s = 2 * pi * (rho1 * rho1 - rho0 * rho0) / Nphi;
            for (int ie = 0; ie < Nph; ++ie) {
                double dE = currentE * (factor - 1.0);
                localFlux += evaluateFluxFromSourceAtPoint(currentE, source, irho, iphi)*dE/s;
                currentE = currentE * factor;
            }
            image[irho][iphi] = localFlux;
        }
    }

    omp_destroy_lock(&my_lock);

    FILE* outFile = fopen(fileName, "w");
    for (int irho = 0; irho < Nrho; ++irho) {
        for (int iphi = 0; iphi < Nphi; ++iphi) {
            fprintf(outFile, "%g ", image[irho][iphi]);
        }
        fprintf(outFile, "\n");
    }
    fclose(outFile);
}

void RadiationEvaluator::writeImageFromSourceAtEToFile(const double& photonFinalEnergy, const char* fileName, RadiationSource* source) {
    int Nrho = source->getNrho();
    int Nz = source->getNz();
    int Nphi = source->getNphi();

    double** image = new double* [Nrho];
    for (int irho = 0; irho < Nrho; ++irho) {
        image[irho] = new double[Nphi];
        for (int iphi = 0; iphi < Nphi; ++iphi) {
            image[irho][iphi] = 0;
        }
    }

    omp_init_lock(&my_lock);

    int irho;
#pragma omp parallel for private(irho) shared(photonFinalEnergy, source, Nrho, Nz, Nphi, image)
    for (irho = 0; irho < Nrho; ++irho) {
        for (int iphi = 0; iphi < Nphi; ++iphi) {
            double rho0 = irho * source->getMaxRho() / Nrho;
            double rho1 = (irho + 1) * source->getMaxRho() / Nrho;
            double s = 2 * pi * (rho1 * rho1 - rho0 * rho0) / Nphi;
            double localFlux =  evaluateFluxFromSourceAtPoint(photonFinalEnergy, source, irho, iphi)/s;
            image[irho][iphi] = localFlux;
        }
    }

    omp_destroy_lock(&my_lock);

    FILE* outFile = fopen(fileName, "w");
    for (int irho = 0; irho < Nrho; ++irho) {
        for (int iphi = 0; iphi < Nphi; ++iphi) {
            fprintf(outFile, "%g ", image[irho][iphi]);
        }
        fprintf(outFile, "\n");
    }
    fclose(outFile);
}