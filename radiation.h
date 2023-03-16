#ifndef RADIATION_H
#define RADIATION_H

#include "math.h"

#include "massiveParticleDistribution.h"
#include "photonDistribution.h"
#include "radiationSource.h"

class RadiationEvaluator{
protected:
    int my_Ne;
    double my_Emin;
    double my_Emax;

    double* my_Ee;
public:
    RadiationEvaluator(int Ne, const double& Emin, const double& Emax){
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
    virtual ~RadiationEvaluator(){
        delete[] my_Ee;
    }

    virtual void resetParameters(const double* parameters, const double* normalizationUnits) = 0;

    //returns fluxes in cm^-2 s^-1 = d energy flux/ d energy
    virtual double evaluateFluxFromIsotropicFunction(const double& photonFinalEnergy, MassiveParticleIsotropicDistribution* electronDistribution, const double& volume, const double& distance) = 0;
    virtual double evaluateFluxFromSource(const double& photonFinalEnergy, RadiationSource* source) = 0;

    double evaluateTotalFluxInEnergyRange(const double& Ephmin, const double& Ephmax, RadiationSource* source) {
        int Nph = 100;
        double factor = pow(Ephmax / Ephmin, 1.0 / (Nph - 1));
        double currentE = Ephmin;
        double flux = 0;
        for (int i = 0; i < Nph; ++i) {
            printf("%d\n", i);
            double dE = currentE * (factor - 1.0);
            flux += evaluateFluxFromSource(currentE, source)*dE;
            currentE = currentE * factor;
        }
        return flux;
    }

    void writeFluxFromSourceToFile(const char* fileName, RadiationSource* source, const double& Ephmin, const double& Ephmax, const int Nph) {
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
};

#endif // RADIATION_H
