#include "stdio.h"
#include "math.h"
#include <omp.h>

#include "constants.h"
#include "massiveParticleDistribution.h"
#include "photonDistribution.h"
#include "util.h"
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
    int Nx1 = source->getNx1();
    int Nz = source->getNz();
    int Nx2 = source->getNx2();

    double result = 0;
    int irho = 0;

    double redShift = source->getRedShift();

    double photonFinalEnergy1 = photonFinalEnergy * (1 + redShift);

    omp_init_lock(&my_lock);

#pragma omp parallel for private(irho) shared(photonFinalEnergy1, source, Nx1, Nz, Nx2) reduction(+:result)

    for (irho = 0; irho < Nx1; ++irho) {
        //omp_set_lock(&my_lock);
        //printLog("irho = %d\n", irho);
        //omp_unset_lock(&my_lock);
        for (int iphi = 0; iphi < Nx2; ++iphi) {
            //omp_set_lock(&my_lock);
            //printLog("iphi = %d\n", iphi);
            //omp_unset_lock(&my_lock);
            /*for (int iz = 0; iz < Nz; ++iz) {
                result += evaluateFluxFromIsotropicFunction(photonFinalEnergy1, source->getParticleDistribution(irho, iz, iphi), source->getVolume(irho, iz, iphi), source->getDistance());
            }*/
            if (source->isSource(irho, iphi)) {
                result += evaluateFluxFromSourceAtPoint(photonFinalEnergy1, source, irho, iphi);
                if ((result != result) || (0 * result != 0 * result)) {
                    omp_set_lock(&my_lock);
                    printf("flux from source = NaN or Infinity at irho = %d iphi = %d\n", irho, iphi);
                    printLog("flux from source = NaN or Infinity at irho = %d iphi = %d\n", irho, iphi);
                    omp_unset_lock(&my_lock);
                    exit(0);
                }
            }
        }
    }

    omp_destroy_lock(&my_lock);

    return result*(1+redShift)*(1+redShift);
}

double RadiationEvaluator::evaluateFluxFromSourceAtPoint(const double& photonFinalEnergy, RadiationSource* source, int ii, int jj) {
    int Nx1 = source->getNx1();
    int Nz = source->getNz();
    int Nx2 = source->getNx2();
    double localI = 0;
    double prevArea = 0;
    for (int iz = 0; iz < Nz; ++iz) {
        double area = source->getArea(ii, iz, jj);
        double A = 0;
        double I = 0;
        double v;
        double theta;
        double phi;
        /*source->getVelocity(irho, iz, iphi, v, theta, phi);
        if (!my_doppler) {
            v = 0;
        }
        double beta = v / speed_of_light;
        double gamma = 1.0 / sqrt(1 - beta * beta);
        double mu = cos(theta);

        double D = gamma * (1.0 - beta * mu);
        double photonFinalEnergyPrimed = photonFinalEnergy * D;*/

        evaluateEmissivityAndAbsorption(photonFinalEnergy, ii, iz, jj, source, I, A);
        /*evaluateEmissivityAndAbsorption(photonFinalEnergyPrimed, irho, iz, iphi, source, I, A);

        if (my_doppler) {
            I = I / (D * D);
            A = A * D;
        }*/

        double length = source->getLength(ii, iz, jj);
        if (length > 0) {
            /*if (my_absorption) {
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
            }*/
            if (my_absorption) {
                double I0 = localI;
                double tempI01 = I0;
                double tempI02 = 0;
                if (area < prevArea) {
                    tempI01 = I0 * area / prevArea;
                    tempI02 = I0 - tempI01;
                }
                prevArea = area;
                double Q = I * area;

                double tau = A * length;
                double S = 0;
                if (A > 0) {
                    S = Q / A;
                }
                if (fabs(tau) < 1E-10) {
                    localI = (tempI01 * (1.0 - tau) + S * tau + tempI02);
                }
                else {
                    localI = (S + (tempI01 - S) * exp(-tau) + tempI02);
                }
            }
            else {
                localI = localI + I * area * length;
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

void RadiationEvaluator::evaluateEmissivityAndAbsorption(const double& photonFinalEnergy, int ix1, int iz, int ix2, RadiationSource* source, double& I, double& A) {
    I = evaluateEmissivity(photonFinalEnergy, ix1, iz, ix2, source);
    if (my_absorption) {
        A = evaluateAbsorption(photonFinalEnergy, ix1, iz, ix2, source);
    }
    else {
        A = 0;
    }
}

void RadiationEvaluator::writeFluxFromSourceToFile(const char* fileName, RadiationSource* source, const double& Ephmin, const double& Ephmax, const int Nph) {
    double factor = pow(Ephmax / Ephmin, 1.0 / (Nph - 1));
    double currentE = Ephmin;
    FILE* outFile = fopen(fileName, "w");
    for (int i = 0; i < Nph; ++i) {
        //omp_set_lock(&my_lock);
        printf("writeFluxFromSourceToFile iph = %d\n", i);
        printLog("writeFluxFromSourceToFile iph = %d\n", i);
        //omp_unset_lock(&my_lock);
        double flux = evaluateFluxFromSource(currentE, source);
        fprintf(outFile, "%g %g\n", currentE, flux);
        currentE = currentE * factor;
    }
    fclose(outFile);
}

void RadiationEvaluator::writeEFEFromSourceToFile(const char* fileName, RadiationSource* source, const double& Ephmin, const double& Ephmax, const int Nph)
{
    double factor = pow(Ephmax / Ephmin, 1.0 / (Nph - 1));
    double currentE = Ephmin;
    FILE* outFile = fopen(fileName, "w");
    for (int i = 0; i < Nph; ++i) {
        //omp_set_lock(&my_lock);
        printf("writeEFEFromSourceToFile iph = %d\n", i);
        printLog("writeEFEFromSourceToFile iph = %d\n", i);
        //omp_unset_lock(&my_lock);
        double flux = evaluateFluxFromSource(currentE, source);
        fprintf(outFile, "%g %g\n", currentE/1.6E-12, currentE*flux);
        currentE = currentE * factor;
    }
    fclose(outFile);
}

void RadiationEvaluator::writeImageFromSourceToFile(const char* fileName, RadiationSource* source, const double& Ephmin, const double& Ephmax, const int Nph) {
    int Nx1 = source->getNx1();
    int Nz = source->getNz();
    int Nx2 = source->getNx2();

    double** image = new double*[Nx1];
    for (int irho = 0; irho < Nx1; ++irho) {
        image[irho] = new double[Nx2];
        for (int iphi = 0; iphi < Nx2; ++iphi) {
            image[irho][iphi] = 0;
        }
    }

    omp_init_lock(&my_lock);

    int irho;
#pragma omp parallel for private(irho) shared(Ephmin, Ephmax, source, Nx1, Nz, Nx2, image)
    for (irho = 0; irho < Nx1; ++irho) {
        omp_set_lock(&my_lock);
        printf("evaluating image irho = %d\n", irho);
        printLog("evaluating image irho = %d\n", irho);
        omp_unset_lock(&my_lock);
        for (int iphi = 0; iphi < Nx2; ++iphi) {
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
    for (int irho = 0; irho < Nx1; ++irho) {
        for (int iphi = 0; iphi < Nx2; ++iphi) {
            fprintf(outFile, "%g ", image[irho][iphi]);
        }
        fprintf(outFile, "\n");
    }
    fclose(outFile);

    for (int irho = 0; irho < Nx1; ++irho) {
        delete[] image[irho];
    }
    delete[] image;
}

void RadiationEvaluator::writeImageFromSourceAtEToFile(const double& photonFinalEnergy, const char* fileName, RadiationSource* source) {
    int Nx1 = source->getNx1();
    int Nz = source->getNz();
    int Nx2 = source->getNx2();

    double** image = new double* [Nx1];
    for (int irho = 0; irho < Nx1; ++irho) {
        image[irho] = new double[Nx2];
        for (int iphi = 0; iphi < Nx2; ++iphi) {
            image[irho][iphi] = 0;
        }
    }

    omp_init_lock(&my_lock);

    int irho;
#pragma omp parallel for private(irho) shared(photonFinalEnergy, source, Nx1, Nz, Nx2, image)
    for (irho = 0; irho < Nx1; ++irho) {
        //printf("i = %d\n", irho);
        for (int iphi = 0; iphi < Nx2; ++iphi) {
            double s = source->getCrossSectionArea(irho, iphi);
            //todo? surface emissivity?
            double d = source->getDistance();
            double localFlux =  evaluateFluxFromSourceAtPoint(photonFinalEnergy, source, irho, iphi)*d*d/s;
            image[irho][iphi] = localFlux;
        }
    }

    omp_destroy_lock(&my_lock);

    FILE* outFile = fopen(fileName, "w");
    for (int irho = 0; irho < Nx1; ++irho) {
        for (int iphi = 0; iphi < Nx2; ++iphi) {
            fprintf(outFile, "%g ", image[irho][iphi]);
        }
        fprintf(outFile, "\n");
    }
    fclose(outFile);

    for (int irho = 0; irho < Nx1; ++irho) {
        delete[] image[irho];
    }
    delete[] image;
}

void RadiationEvaluator::writeImageFromSourceInRangeToFile(const double& photonEmin, const double& photonEmax, int N, const char* fileName, RadiationSource* source) {
    int Nx1 = source->getNx1();
    int Nz = source->getNz();
    int Nx2 = source->getNx2();

    double** image = new double* [Nx1];
    for (int irho = 0; irho < Nx1; ++irho) {
        image[irho] = new double[Nx2];
        for (int iphi = 0; iphi < Nx2; ++iphi) {
            image[irho][iphi] = 0;
        }
    }

    double factor = pow(photonEmax / photonEmin, 1.0 / (N - 1));
    double currentE = photonEmin;

    FILE* outFile = fopen(fileName, "w");
    for (int k = 0; k < N; ++k) {
        omp_init_lock(&my_lock);
        int irho;
        double dE = currentE * (factor - 1);
#pragma omp parallel for private(irho) shared(k, currentE, source, Nx1, Nz, Nx2, image)
        for (irho = 0; irho < Nx1; ++irho) {
            //printf("i = %d\n", irho);
            for (int iphi = 0; iphi < Nx2; ++iphi) {
                double s = source->getCrossSectionArea(irho, iphi);
                double d = source->getDistance();
                double localFlux = evaluateFluxFromSourceAtPoint(currentE, source, irho, iphi)*d*d / s;
                image[irho][iphi] = localFlux;
            }
        }
        for (int irho = 0; irho < Nx1; ++irho) {
            for (int iphi = 0; iphi < Nx2; ++iphi) {
                fprintf(outFile, "%g ", image[irho][iphi]);
            }
            fprintf(outFile, "\n");
        }
        currentE = currentE * factor;

        omp_destroy_lock(&my_lock);

    }
    fclose(outFile);

    for(int irho = 0; irho < Nx1; ++irho) {
        delete[] image[irho];
    }
    delete[] image;
}

void RadiationEvaluator::initLock()
{
    omp_init_lock(&my_lock);
}

void RadiationEvaluator::destroyLock()
{
    omp_destroy_lock(&my_lock);
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

double RadiationSumEvaluator::evaluateFluxFromSourceAtPoint(const double& photonFinalEnergy, RadiationSource* source, int ii, int jj)
{
    double result = 0;
    for (int i = 0; i < my_Nevaluators; ++i) {
        result = result + my_Evaluators[i]->evaluateFluxFromSourceAtPoint(photonFinalEnergy, source, ii, jj);
    }
    return result;
}

double RadiationSumEvaluator::evaluateEmissivity(const double& photonFinalEnergy, int ix1, int iz, int ix2, RadiationSource* source)
{
    double result = 0;
    for (int i = 0; i < my_Nevaluators; ++i) {
        result += my_Evaluators[i]->evaluateEmissivity(photonFinalEnergy, ix1, iz, ix2, source);
    }
    return result;
}

double RadiationSumEvaluator::evaluateAbsorption(const double& photonFinalEnergy, int ix1, int iz, int ix2, RadiationSource* source)
{
    double result = 0;
    for (int i = 0; i < my_Nevaluators; ++i) {
        result += my_Evaluators[i]->evaluateAbsorption(photonFinalEnergy, ix1, iz, ix2, source);
    }
    return result;
}

void RadiationSumEvaluator::initLock()
{
    RadiationEvaluator::initLock();
    for (int i = 0; i < my_Nevaluators; ++i) {
        my_Evaluators[i]->initLock();
    }
}

void RadiationSumEvaluator::destroyLock()
{
    RadiationEvaluator::destroyLock();
    for (int i = 0; i < my_Nevaluators; ++i) {
        my_Evaluators[i]->destroyLock();
    }
}
