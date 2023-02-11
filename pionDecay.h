#ifndef PION_DECAY_H
#define PION_DECAY_H

#include "massiveParticleDistribution.h"
#include "radiationSource.h"
#include "radiation.h"

class PionDecayEvaluatorBase : public RadiationEvaluator{
protected:
    double my_ambientConcentration;
public:
    PionDecayEvaluatorBase(int Ne, double Emin, double Emax, const double& ambientConcentration);
    virtual ~PionDecayEvaluatorBase();
	double sigmaInelastic(const double& energy);
};

class PionDecayEvaluator : public PionDecayEvaluatorBase{
protected:
    void getCoefs(double& alpha, double& beta, double& gamma, double& lambda, const double& protonEnergy);
    void getBcoefs(double& b1, double& b2, double& b3, const double& protonEnergy);
public:
    PionDecayEvaluator(int Ne, double Emin, double Emax, const double& ambientConcentration);
    ~PionDecayEvaluator();
    double sigmaPion(const double& energy);
    double sigma2Pion(const double& energy);
    double sigmaGamma(const double& photonEnergy, const double& protonEnergy);
    void resetParameters(const double *parameters, const double *normalizationUnits);
    double evaluateFluxFromIsotropicFunction(const double& photonFinalEnergy, MassiveParticleIsotropicDistribution* protonDistribution, const double& volume, const double& distance);
    double evaluateFluxFromSource(const double& photonFinalEnergy, RadiationSource* source);
};

class PionDecayEvaluatorKelner : public PionDecayEvaluatorBase{
protected:
public:
    PionDecayEvaluatorKelner(int Ne, double Emin, double Emax, const double& ambientConcentration);
    ~PionDecayEvaluatorKelner();
    double functionKelner(const double& x, const double& protonEnergy);
    void resetParameters(const double *parameters, const double *normalizationUnits);
    double evaluateFluxFromIsotropicFunction(const double& photonFinalEnergy, MassiveParticleIsotropicDistribution* protonDistribution, const double& volume, const double& distance);
    double evaluateFluxFromSource(const double& photonFinalEnergy, RadiationSource* source);
};

#endif
