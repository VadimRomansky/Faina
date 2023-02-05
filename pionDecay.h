#ifndef PION_DECAY_H
#define PION_DECAY_H

#include "massiveParticleDistribution.h"
#include "radiationSource.h"

class PionDecayEvaluator {
protected:
	int my_Ne;
	double my_Emin;
	double my_Emax;

	double* my_Ee;
	void getCoefs(double& alpha, double& beta, double& gamma, double& lambda, const double& protonEnergy);
	void getBcoefs(double& b1, double& b2, double& b3, const double& protonEnergy);
public:
	PionDecayEvaluator(int Ne, double Emin, double Emax);
	~PionDecayEvaluator();
	double sigmaInelastic(const double& energy);
	double sigmaPion(const double& energy);
	double sigma2Pion(const double& energy);
	double sigmaGamma(const double& photonEnergy, const double& protonEnergy);
	double evaluatePionDecayLuminocityIsotropicFunction(const double& photonFinalEnergy, MassiveParticleIsotropicDistribution* protonDistribution, const double& ambientConcentration, const double& volume, const double& distance);
	double evaluatePionDecayIsotropicFluxFromSource(const double& photonFinalEnergy, RadiationSource* source);
};

#endif
