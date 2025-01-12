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

    bool my_absorption;
    bool my_doppler;
public:
    RadiationEvaluator(int Ne, const double& Emin, const double& Emax, bool absorption, bool doppler);
    virtual ~RadiationEvaluator();

    virtual void resetParameters(const double* parameters, const double* normalizationUnits) = 0;

    //returns fluxes in cm^-2 s^-1 = d energy flux/ d energy
    //virtual double evaluateFluxFromIsotropicFunction(const double& photonFinalEnergy, MassiveParticleIsotropicDistribution* electronDistribution, const double& volume, const double& distance) = 0;

    virtual double evaluateFluxFromSource(const double& photonFinalEnergy, RadiationSource* source);

    virtual double evaluateFluxFromSourceAtPoint(const double& photonFinalEnergy, RadiationSource* source, int ii, int jj);

    double evaluateTotalFluxInEnergyRange(const double& Ephmin, const double& Ephmax, int Nph, RadiationSource* source);

    virtual void evaluateEmissivityAndAbsorption(const double& photonFinalEnergy, int ix1, int iz, int ix2, RadiationSource* source, double& I, double& A);

    virtual double evaluateEmissivity(const double& photonFinalEnergy, int ix1, int iz, int ix2, RadiationSource* source) = 0;

    virtual double evaluateAbsorption(const double& photonFinalEnergy, int ix1, int iz, int ix2, RadiationSource* source) = 0;

    void writeFluxFromSourceToFile(const char* fileName, RadiationSource* source, const double& Ephmin, const double& Ephmax, const int Nph);

    void writeEFEFromSourceToFile(const char* fileName, RadiationSource* source, const double& Ephmin, const double& Ephmax, const int Nph);

    void writeImageFromSourceToFile(const char* fileName, RadiationSource* source, const double& Ephmin, const double& Ephmax, const int Nph);

    void writeImageFromSourceAtEToFile(const double& photonFinalEnergy, const char* fileName, RadiationSource* source);

    void writeImageFromSourceInRangeToFile(const double& photonEmin, const double& photonEmax, int N, const char* fileName, RadiationSource* source);
};

class RadiationSumEvaluator : public RadiationEvaluator {
protected:
    int my_Nevaluators;
    RadiationEvaluator** my_Evaluators;
public:
    RadiationSumEvaluator(int Ne, const double& Emin, const double& Emax, RadiationEvaluator* evaluator1, RadiationEvaluator* evaluator2, bool absorption = true, bool doppler = false);
    RadiationSumEvaluator(int Ne, const double& Emin, const double& Emax, int Nev, RadiationEvaluator** evaluators, bool absorption = true, bool doppler = false);
    ~RadiationSumEvaluator();

    virtual void resetParameters(const double* parameters, const double* normalizationUnits);
    virtual double evaluateFluxFromSource(const double& photonFinalEnergy, RadiationSource* source);
    virtual double evaluateFluxFromSourceAtPoint(const double& photonFinalEnergy, RadiationSource* source, int ii, int jj);
    virtual double evaluateEmissivity(const double& photonFinalEnergy, int ix1, int iz, int ix2, RadiationSource* source);
    virtual double evaluateAbsorption(const double& photonFinalEnergy, int ix1, int iz, int ix2, RadiationSource* source);
};

#endif // RADIATION_H
