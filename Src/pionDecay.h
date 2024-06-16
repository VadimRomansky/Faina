#ifndef PION_DECAY_H
#define PION_DECAY_H

#include "massiveParticleDistribution.h"
#include "radiationSource.h"
#include "radiation.h"

class PionDecayEvaluatorBase : public RadiationEvaluator{
protected:
    double my_ambientConcentration;
public:
    PionDecayEvaluatorBase(int Ne, double Emin, double Emax, const double& ambientConcentration, bool absorption = false, bool doppler = false);
    virtual ~PionDecayEvaluatorBase();
	double sigmaInelastic(const double& energy);
};

class PionDecayEvaluator : public PionDecayEvaluatorBase{
protected:
    void getCoefs(double& alpha, double& beta, double& gamma, double& lambda, const double& protonEnergy);
    void getBcoefs(double& b1, double& b2, double& b3, const double& protonEnergy);
    double sigmaPion(const double& energy);
    double sigma2Pion(const double& energy);
public:
    PionDecayEvaluator(int Ne, double Emin, double Emax, const double& ambientConcentration, bool absorption = false, bool doppler = false);
    ~PionDecayEvaluator();
    double sigmaGamma(const double& photonEnergy, const double& protonEnergy);
    void resetParameters(const double *parameters, const double *normalizationUnits);
    double evaluateFluxFromIsotropicFunction(const double& photonFinalEnergy, MassiveParticleIsotropicDistribution* protonDistribution, const double& volume, const double& distance);
    virtual double evaluateEmissivity(const double& photonFinalEnergy, int irho, int iz, int iphi, RadiationSource* source);
    virtual double evaluateAbsorbtion(const double& photonFinalEnergy, int irho, int iz, int iphi, RadiationSource* source);
    //virtual double evaluateFluxFromSourceAtPoint(const double& photonFinalEnergy, RadiationSource* source, int irho, int iphi);
};

class PionDecayEvaluatorKelner : public PionDecayEvaluatorBase{
protected:
    double functionKelner(const double& x, const double& protonEnergy);
public:
    PionDecayEvaluatorKelner(int Ne, double Emin, double Emax, const double& ambientConcentration, bool absorption = false, bool doppler = false);
    ~PionDecayEvaluatorKelner();
    void resetParameters(const double *parameters, const double *normalizationUnits);
    double evaluateFluxFromIsotropicFunction(const double& photonFinalEnergy, MassiveParticleIsotropicDistribution* protonDistribution, const double& volume, const double& distance);
    virtual double evaluateEmissivity(const double& photonFinalEnergy, int irho, int iz, int iphi, RadiationSource* source);
    virtual double evaluateAbsorbtion(const double& photonFinalEnergy, int irho, int iz, int iphi, RadiationSource* source);
    //virtual double evaluateFluxFromSourceAtPoint(const double& photonFinalEnergy, RadiationSource* source, int irho, int iphi);
};

#endif
