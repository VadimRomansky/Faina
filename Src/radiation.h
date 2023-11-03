#ifndef RADIATION_H
#define RADIATION_H

#include "math.h"
#include <omp.h>

#include "massiveParticleDistribution.h"
#include "photonDistribution.h"
#include "radiationSource.h"

class RadiationEvaluator{
protected:
    omp_lock_t my_lock;

    int my_Ne;
    double my_Emin;
    double my_Emax;

    double* my_Ee;
public:
    RadiationEvaluator(int Ne, const double& Emin, const double& Emax);
    virtual ~RadiationEvaluator();

    virtual void resetParameters(const double* parameters, const double* normalizationUnits) = 0;

    //returns fluxes in cm^-2 s^-1 = d energy flux/ d energy
    //virtual double evaluateFluxFromIsotropicFunction(const double& photonFinalEnergy, MassiveParticleIsotropicDistribution* electronDistribution, const double& volume, const double& distance) = 0;

    virtual double evaluateFluxFromSource(const double& photonFinalEnergy, RadiationSource* source);

    virtual double evaluateFluxFromSourceAtPoint(const double& photonFinalEnergy, RadiationSource* source, int rhoi, int phi) = 0;

    double evaluateTotalFluxInEnergyRange(const double& Ephmin, const double& Ephmax, int Nph, RadiationSource* source);

    void writeFluxFromSourceToFile(const char* fileName, RadiationSource* source, const double& Ephmin, const double& Ephmax, const int Nph);

    void writeImageFromSourceToFile(const char* fileName, RadiationSource* source, const double& Ephmin, const double& Ephmax, const int Nph);

    void writeImageFromSourceAtEToFile(const double& photonFinalEnergy, const char* fileName, RadiationSource* source);
};

class RadiationSumEvaluator : public RadiationEvaluator {
protected:
    int my_Nevaluators;
    RadiationEvaluator** my_Evaluators;
public:
    RadiationSumEvaluator(int Ne, const double& Emin, const double& Emax, RadiationEvaluator* evaluator1, RadiationEvaluator* evaluator2);
    RadiationSumEvaluator(int Ne, const double& Emin, const double& Emax, int Nev, RadiationEvaluator** evaluators);
    ~RadiationSumEvaluator();

    virtual void resetParameters(const double* parameters, const double* normalizationUnits);
    virtual double evaluateFluxFromSource(const double& photonFinalEnergy, RadiationSource* source);
    virtual double evaluateFluxFromSourceAtPoint(const double& photonFinalEnergy, RadiationSource* source, int rhoi, int phi);
};

#endif // RADIATION_H
