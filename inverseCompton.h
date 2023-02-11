#ifndef inverse_compton_h
#define inverse_compton_h

#include "massiveParticleDistribution.h"
#include "photonDistribution.h"
#include "radiationSource.h"
#include "radiation.h"

class InverseComptonEvaluator : public RadiationEvaluator{
protected:
	int my_Nmu;
	int my_Nphi;

	double* my_cosTheta;
	double* my_cosThetaLeft;
	double* my_dcosTheta;

    PhotonIsotropicDistribution* my_photonDistribution;
public:
    InverseComptonEvaluator(int Ne, int Nmu, int Nphi, double Emin, double Emax, PhotonIsotropicDistribution* photonDistribution);
	~InverseComptonEvaluator();
    void resetParameters(const double *parameters, const double *normalizationUnits);
	double evaluateComptonLuminocity(const double& photonFinalEnergy, const double& photonFinalTheta, const double& photonFinalPhi, PhotonDistribution* photonDistribution, MassiveParticleDistribution* electronDistribution, const double& volume, const double& distance);
    double evaluateFluxFromIsotropicFunction(const double& photonFinalEnergy, MassiveParticleIsotropicDistribution* electronDistribution, const double& volume, const double& distance);
    double evaluateFluxFromSource(const double& photonFinalEnergy, RadiationSource* source);
};
#endif
