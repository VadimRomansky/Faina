#ifndef PION_DECAY_H
#define PION_DECAY_H

#include "massiveParticleDistribution.h"
#include "radiationSource.h"
#include "radiation.h"

class PionDecayEvaluator : public RadiationEvaluator{
protected:
    double my_ambientConcentration;

	void getCoefs(double& alpha, double& beta, double& gamma, double& lambda, const double& protonEnergy);
	void getBcoefs(double& b1, double& b2, double& b3, const double& protonEnergy);
public:
    PionDecayEvaluator(int Ne, double Emin, double Emax, const double& ambientConcentration);
	~PionDecayEvaluator();
	double sigmaInelastic(const double& energy);
	double sigmaPion(const double& energy);
	double sigma2Pion(const double& energy);
	double sigmaGamma(const double& photonEnergy, const double& protonEnergy);
	double functionKelner(const double& x, const double& protonEnergy);
    double evaluateFluxFromIsotropicFunction(const double& photonFinalEnergy, MassiveParticleIsotropicDistribution* protonDistribution, const double& volume, const double& distance);
    double evaluatePionDecayKelnerLuminocityIsotropicFunction(const double& photonFinalEnergy, MassiveParticleIsotropicDistribution* protonDistribution, const double& volume, const double& distance);
    double evaluateFluxFromSource(const double& photonFinalEnergy, RadiationSource* source);
	double evaluatePionDecayKelnerIsotropicFluxFromSource(const double& photonFinalEnergy, RadiationSource* source);
};

#endif
