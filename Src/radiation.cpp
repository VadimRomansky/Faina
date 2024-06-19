#include "stdio.h"
#include "math.h"
#include <omp.h>

#include "constants.h"
#include "massiveParticleDistribution.h"
#include "photonDistribution.h"
#include "radiationSource.h"
#include "radiation.h"

RadiationEvaluator::RadiationEvaluator(int Ne, const double& Emin, const double& Emax, bool absorption, bool doppler) {
    my_Ne = Ne;
    my_Emin = Emin;
    my_Emax = Emax;
    my_Ee = new double[my_Ne];
    double factor = pow(my_Emax / my_Emin, 1.0 / (my_Ne - 1));

    my_Ee[0] = my_Emin;
    for (int i = 1; i < my_Ne; ++i) {
        my_Ee[i] = my_Ee[i - 1] * factor;
    }

    my_absorption = absorption;
    my_doppler = doppler;
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

double RadiationEvaluator::evaluateFluxFromSourceAtPoint(const double& photonFinalEnergy, RadiationSource* source, int irho, int iphi) {
    int Nrho = source->getNrho();
    int Nz = source->getNz();
    int Nphi = source->getNphi();
    double photonFinalFrequency = photonFinalEnergy / hplank;
    double localI = 0;
    double prevArea = 0;
    for (int iz = 0; iz < Nz; ++iz) {
        double area = source->getArea(irho, iz, iphi);
        double A = 0;
        double I = 0;
        double v;
        double theta;
        double phi;
        source->getVelocity(irho, iz, iphi, v, theta, phi);
        if (!my_doppler) {
            v = 0;
        }
        double beta = v / speed_of_light;
        double gamma = 1.0 / sqrt(1 - beta * beta);
        double mu = cos(theta);

        double D = gamma * (1.0 - beta * mu);
        double photonFinalFrequencyPrimed = photonFinalFrequency * D;

        I = evaluateEmissivity(photonFinalEnergy, irho, iz, iphi, source);
        A = evaluateAbsorbtion(photonFinalEnergy, irho, iz, iphi, source);

        double length = source->getLength(irho, iz, iphi);
        if (length > 0) {
            if (my_absorption) {
                double I0 = localI * D * D;
                double tempI01 = I0;
                double tempI02 = 0;
                if (area < prevArea) {
                    tempI01 = I0 * area / prevArea;
                    tempI02 = I0 - tempI01;
                }
                prevArea = area;
                double Q = I * area;
                //todo lorentz length
                double lnorm = fabs(length * sin(theta));
                double lpar = fabs(length * cos(theta));
                double lengthPrimed = sqrt(lnorm * lnorm + lpar * lpar * gamma * gamma);
                double tau = A * lengthPrimed;
                double S = 0;
                if (A > 0) {
                    S = Q / A;
                }
                if (fabs(tau) < 1E-10) {
                    localI = (tempI01 * (1.0 - tau) + S * tau + tempI02) / (D * D);
                }
                else {
                    localI = (S + (tempI01 - S) * exp(-tau) + tempI02) / (D * D);
                }
            }
            else {
                localI = localI + I * area * length / (D * D);
            }
        }
    }
    double distance = source->getDistance();
    return localI / (distance * distance);
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
            double s = source->getCrossSectionArea(irho, iphi);
            double d = source->getDistance();
            for (int ie = 0; ie < Nph; ++ie) {
                double dE = currentE * (factor - 1.0);
                localFlux += evaluateFluxFromSourceAtPoint(currentE, source, irho, iphi)*dE*d*d/s;
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

    for (int irho = 0; irho < Nrho; ++irho) {
        delete[] image[irho];
    }
    delete[] image;
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
        //printf("i = %d\n", irho);
        for (int iphi = 0; iphi < Nphi; ++iphi) {
            double s = source->getCrossSectionArea(irho, iphi);
            //todo? surface emissivity?
            double d = source->getDistance();
            double localFlux =  evaluateFluxFromSourceAtPoint(photonFinalEnergy, source, irho, iphi)*d*d/s;
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

    for (int irho = 0; irho < Nrho; ++irho) {
        delete[] image[irho];
    }
    delete[] image;
}

void RadiationEvaluator::writeImageFromSourceInRangeToFile(const double& photonEmin, const double& photonEmax, int N, const char* fileName, RadiationSource* source) {
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

    double factor = pow(photonEmax / photonEmin, 1.0 / (N - 1));
    double currentE = photonEmin;

    FILE* outFile = fopen(fileName, "w");
    for (int k = 0; k < N; ++k) {
        omp_init_lock(&my_lock);
        int irho;
#pragma omp parallel for private(irho) shared(k, currentE, source, Nrho, Nz, Nphi, image)
        for (irho = 0; irho < Nrho; ++irho) {
            //printf("i = %d\n", irho);
            for (int iphi = 0; iphi < Nphi; ++iphi) {
                double s = source->getCrossSectionArea(irho, iphi);
                double d = source->getDistance();
                double localFlux = evaluateFluxFromSourceAtPoint(currentE, source, irho, iphi)*d*d / s;
                image[irho][iphi] = localFlux;
            }
        }
        currentE = currentE * factor;

        omp_destroy_lock(&my_lock);

        for (int irho = 0; irho < Nrho; ++irho) {
            for (int iphi = 0; iphi < Nphi; ++iphi) {
                fprintf(outFile, "%g ", image[irho][iphi]);
            }
            fprintf(outFile, "\n");
        }
    }
    fclose(outFile);

    for(int irho = 0; irho < Nrho; ++irho) {
        delete[] image[irho];
    }
    delete[] image;
}

RadiationSumEvaluator::RadiationSumEvaluator(int Ne, const double& Emin, const double& Emax, RadiationEvaluator* evaluator1, RadiationEvaluator* evaluator2, bool absorption, bool doppler) : RadiationEvaluator(Ne, Emin, Emax, absorption, doppler) {
    my_Nevaluators = 2;
    my_Evaluators = new RadiationEvaluator*[my_Nevaluators];
    my_Evaluators[0] = evaluator1;
    my_Evaluators[1] = evaluator2;
}

RadiationSumEvaluator::RadiationSumEvaluator(int Ne, const double& Emin, const double& Emax, int Nev, RadiationEvaluator** evaluators, bool absorption, bool doppler) : RadiationEvaluator(Ne, Emin, Emax, absorption, doppler) {
    my_Nevaluators = Nev;
    my_Evaluators = new RadiationEvaluator*[my_Nevaluators];
    for (int i = 0; i < my_Nevaluators; ++i) {
        my_Evaluators[i] = evaluators[i];
    }
}

RadiationSumEvaluator::~RadiationSumEvaluator()
{
    delete[] my_Evaluators;
}

void RadiationSumEvaluator::resetParameters(const double* parameters, const double* normalizationUnits)
{
    for (int i = 0; i < my_Nevaluators; ++i) {
        my_Evaluators[i]->resetParameters(parameters, normalizationUnits);
    }
}

double RadiationSumEvaluator::evaluateFluxFromSource(const double& photonFinalEnergy, RadiationSource* source)
{
    double result = 0;
    for (int i = 0; i < my_Nevaluators; ++i) {
        result = result + my_Evaluators[i]->evaluateFluxFromSource(photonFinalEnergy, source);
    }
    return result;
}

double RadiationSumEvaluator::evaluateFluxFromSourceAtPoint(const double& photonFinalEnergy, RadiationSource* source, int rhoi, int phi)
{
    double result = 0;
    for (int i = 0; i < my_Nevaluators; ++i) {
        result = result + my_Evaluators[i]->evaluateFluxFromSourceAtPoint(photonFinalEnergy, source, rhoi, phi);
    }
    return result;
}

double RadiationSumEvaluator::evaluateEmissivity(const double& photonFinalEnergy, int irho, int iz, int iphi, RadiationSource* source)
{
    double result = 0;
    for (int i = 0; i < my_Nevaluators; ++i) {
        result += my_Evaluators[i]->evaluateEmissivity(photonFinalEnergy, irho, iz, iphi, source);
    }
    return result;
}

double RadiationSumEvaluator::evaluateAbsorbtion(const double& photonFinalEnergy, int irho, int iz, int iphi, RadiationSource* source)
{
    double result = 0;
    for (int i = 0; i < my_Nevaluators; ++i) {
        result += my_Evaluators[i]->evaluateAbsorbtion(photonFinalEnergy, irho, iz, iphi, source);
    }
    return result;
}
