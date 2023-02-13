#ifndef inverse_compton_h
#define inverse_compton_h

#include "massiveParticleDistribution.h"
#include "photonDistribution.h"
#include "radiationSource.h"
#include "radiation.h"

enum COMPTON_SOLVER_TYPE {ISOTROPIC_THOMPSON, ISOTROPIC_KANG_JONES, KLEIN_NISINA};

class InverseComptonEvaluator : public RadiationEvaluator{
protected:
	COMPTON_SOLVER_TYPE my_solverType;

	int my_Nmu;
	int my_Nphi;

	double* my_cosTheta;
	double* my_cosThetaLeft;
	double* my_dcosTheta;

    PhotonIsotropicDistribution* my_photonDistribution;
public:
    InverseComptonEvaluator(int Ne, int Nmu, int Nphi, double Emin, double Emax, PhotonIsotropicDistribution* photonDistribution, COMPTON_SOLVER_TYPE solverType);
	~InverseComptonEvaluator();
    void resetParameters(const double *parameters, const double *normalizationUnits);

	double evaluateComptonLuminocityThompsonIsotropic(const double& photonFinalEnergy, PhotonIsotropicDistribution* photonDistribution, MassiveParticleIsotropicDistribution* electronDistribution, const double& volume, const double& distance);
	double evaluateComptonLuminocityKangJonesIsotropic(const double& photonFinalEnergy, PhotonIsotropicDistribution* photonDistribution, MassiveParticleIsotropicDistribution* electronDistribution, const double& volume, const double& distance);
	double evaluateComptonLuminocityKleinNisinaIsotropic(const double& photonFinalEnergyi, PhotonIsotropicDistribution* photonDistribution, MassiveParticleIsotropicDistribution* electronDistribution, const double& volume, const double& distance);
	double evaluateComptonLuminocityKleinNisinaAnisotropic(const double& photonFinalEnergy, const double& photonFinalTheta, const double& photonFinalPhi, PhotonDistribution* photonDistribution, MassiveParticleDistribution* electronDistribution, const double& volume, const double& distance);

	double evaluateComptonLuminocity(const double& photonFinalEnergy, const double& photonFinalTheta, const double& photonFinalPhi, PhotonDistribution* photonDistribution, MassiveParticleDistribution* electronDistribution, const double& volume, const double& distance);
    double evaluateFluxFromIsotropicFunction(const double& photonFinalEnergy, MassiveParticleIsotropicDistribution* electronDistribution, const double& volume, const double& distance);
    double evaluateFluxFromSource(const double& photonFinalEnergy, RadiationSource* source);
};

#endif
